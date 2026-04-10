import copy
import os
import shutil
import subprocess
from pathlib import Path

from dqm_pipeline.core import build_output_root, emit_log, expand_path_patterns, sanitize


def _collect_input_files(era_cfg, config_dir):
    files = []
    if "file" in era_cfg:
        files.extend(expand_path_patterns(era_cfg["file"], config_dir))
    if "files" in era_cfg:
        files.extend(expand_path_patterns(era_cfg["files"], config_dir))
    if "file_glob" in era_cfg:
        files.extend(expand_path_patterns(era_cfg["file_glob"], config_dir))
    if "run_files" in era_cfg:
        for file_list in era_cfg["run_files"].values():
            values = file_list if isinstance(file_list, list) else [file_list]
            files.extend(values)

    out = []
    seen = set()
    for item in files:
        if item not in seen:
            out.append(str(item))
            seen.add(item)
    return out


def _workflow_path(cfg, era, era_cfg, preprocess_cfg):
    template = preprocess_cfg.get("workflow_template")
    if template:
        return str(template).format(
            tag=sanitize(str(cfg.get("tag", "DQMIO2ROOT"))),
            resonance=sanitize(str(cfg.get("resonance", "resonance"))),
            era=sanitize(str(era)),
            label=sanitize(str(era_cfg.get("label", era))),
        )
    return f"/{sanitize(str(cfg.get('tag', 'DQMIO2ROOT')))}/{sanitize(str(era_cfg.get('label', era)))}/DQMIO2ROOT"


def _write_cmsrun_cfg(cfg_path, input_files, workflow, force_run_number, convention):
    input_literal = ",\n        ".join(repr(str(x)) for x in input_files)
    cfg_text = f"""import FWCore.ParameterSet.Config as cms

process = cms.Process("DQMCONVERT")
process.load("DQMServices.Core.DQMStore_cfi")

process.source = cms.Source(
    "DQMRootSource",
    fileNames=cms.untracked.vstring(
        {input_literal}
    ),
)

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(-1))

process.dqmSaver = cms.EDAnalyzer(
    "DQMFileSaver",
    convention=cms.untracked.string({convention!r}),
    workflow=cms.untracked.string({workflow!r}),
    saveByRun=cms.untracked.int32(1),
    saveAtJobEnd=cms.untracked.bool(True),
    forceRunNumber=cms.untracked.int32({int(force_run_number)}),
)

process.p = cms.Path(process.dqmSaver)
"""
    cfg_path.write_text(cfg_text, encoding="utf-8")


def _find_converted_root(workdir, force_run_number):
    candidates = sorted(workdir.glob("DQM_V*.root"))
    if not candidates:
        return None

    run_token = f"R{int(force_run_number):09d}"
    matching = [x for x in candidates if run_token in x.name]
    if matching:
        return matching[-1]
    return candidates[-1]


def _default_cmssw_env_command():
    cmssw_base = os.environ.get("CMSSW_BASE")
    if not cmssw_base:
        return None
    return (
        "source /cvmfs/cms.cern.ch/cmsset_default.sh && "
        f"cd {cmssw_base} && "
        "eval `scram runtime -sh`"
    )


def _build_cmssw_subprocess_env(preprocess_cfg):
    if not bool(preprocess_cfg.get("clean_env", True)):
        return None

    env = dict(os.environ)
    unset_defaults = {
        "PYTHONHOME",
        "PYTHONPATH",
        "PYTHONSTARTUP",
        "PYTHONUSERBASE",
        "PYTHONNOUSERSITE",
        "VIRTUAL_ENV",
        "CONDA_PREFIX",
        "CONDA_DEFAULT_ENV",
        "__PYVENV_LAUNCHER__",
        "ROOTSYS",
        "LD_PRELOAD",
    }
    for key in unset_defaults:
        env.pop(key, None)

    # LCG/env.sh style runtime often poisons cmsRun bootstrap; let scram rebuild these.
    for key in ("LD_LIBRARY_PATH", "LIBRARY_PATH", "CPLUS_INCLUDE_PATH", "C_INCLUDE_PATH"):
        env.pop(key, None)

    for key in preprocess_cfg.get("unset_env_vars", []):
        env.pop(str(key), None)

    cmssw_base = preprocess_cfg.get("cmssw_base", os.environ.get("CMSSW_BASE"))
    if cmssw_base:
        env["CMSSW_BASE"] = str(cmssw_base)

    if not env.get("PATH"):
        env["PATH"] = "/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"
    return env


def _run_cmsrun(cfg_path, cmsrun_bin, era_dir, log_path, env_command=None, env=None):
    with log_path.open("w", encoding="utf-8") as log_handle:
        if env_command:
            cmd = f"{env_command} && {cmsrun_bin} {cfg_path.name}"
            proc = subprocess.run(
                ["bash", "--noprofile", "--norc", "-lc", cmd],
                cwd=str(era_dir),
                stdout=log_handle,
                stderr=subprocess.STDOUT,
                text=True,
                env=env,
            )
        else:
            proc = subprocess.run(
                [cmsrun_bin, str(cfg_path)],
                cwd=str(era_dir),
                stdout=log_handle,
                stderr=subprocess.STDOUT,
                text=True,
                env=env,
            )
    return proc.returncode


def preprocess_dqmio_inputs(cfg, config_dir, strict=False, progress=None):
    preprocess_cfg = cfg.get("preprocess_dqmio")
    if not preprocess_cfg or not preprocess_cfg.get("enabled", False):
        return cfg

    cfg = copy.deepcopy(cfg)
    out_root = build_output_root(cfg)
    preprocess_root = out_root / str(preprocess_cfg.get("output_subdir", "preprocessed_dqmio"))
    preprocess_root.mkdir(parents=True, exist_ok=True)

    eras = cfg.get("eras", {})
    selected_eras = [
        era for era, era_cfg in eras.items()
        if str(era_cfg.get("input_format", "")).lower() == "dqmio" or bool(era_cfg.get("preprocess_dqmio", False))
    ]
    if not selected_eras:
        emit_log(progress, "[preprocess_dqmio] no eras requested preprocessing", style="bright_black")
        return cfg

    emit_log(progress, f"[preprocess_dqmio] start eras={len(selected_eras)} out_dir={preprocess_root}", style="bold cyan")
    task_id = None
    if progress is not None:
        task_id = progress.add_task("Preprocess DQMIO", total=len(selected_eras))

    cmsrun_bin = str(preprocess_cfg.get("cmsrun_bin", "cmsRun"))
    cmsrun_env_command = preprocess_cfg.get("env_command") or _default_cmssw_env_command()
    cmsrun_env = _build_cmssw_subprocess_env(preprocess_cfg)
    convention = str(preprocess_cfg.get("convention", "Offline"))
    skip_existing = bool(preprocess_cfg.get("skip_existing", True))

    for era in selected_eras:
        era_cfg = cfg["eras"][era]
        emit_log(progress, f"[preprocess_dqmio] era start: {era}", style="cyan")
        input_files = _collect_input_files(era_cfg, config_dir)
        if not input_files:
            msg = f"[preprocess_dqmio] era={era}: no input files found for DQMIO preprocessing."
            if strict:
                raise RuntimeError(msg)
            emit_log(progress, f"{msg} Skipping.", style="yellow")
            if progress is not None and task_id is not None:
                progress.update(task_id, advance=1)
            continue

        force_run_number = era_cfg.get("run_number", preprocess_cfg.get("force_run_number", cfg.get("run_number")))
        if force_run_number is None:
            raise RuntimeError(
                f"[preprocess_dqmio] era={era}: run_number is required for DQMIO preprocessing when input filenames do not encode a run."
            )

        era_dir = preprocess_root / sanitize(str(era))
        era_dir.mkdir(parents=True, exist_ok=True)
        output_file = era_dir / f"{sanitize(str(era))}_legacy.root"
        if skip_existing and output_file.exists():
            emit_log(progress, f"[preprocess_dqmio] era={era}: reuse existing {output_file}", style="bright_black")
        else:
            if cmsrun_env_command and "CMSSW_BASE" in cmsrun_env_command and not (cmsrun_env or {}).get("CMSSW_BASE"):
                raise RuntimeError(
                    "[preprocess_dqmio] CMSSW_BASE is not set. Set preprocess_dqmio.cmssw_base "
                    "or run from a CMSSW area."
                )
            if cmsrun_env_command:
                emit_log(
                    progress,
                    f"[preprocess_dqmio] era={era}: using CMSSW env command (clean_env={bool(preprocess_cfg.get('clean_env', True))})",
                    style="bright_black",
                )
            elif shutil.which(cmsrun_bin) is None:
                raise RuntimeError(
                    f"[preprocess_dqmio] cmsRun binary '{cmsrun_bin}' not found in PATH. "
                    "Set preprocess_dqmio.env_command or prepare a CMSSW environment."
                )

            workflow = _workflow_path(cfg, era, era_cfg, preprocess_cfg)
            cfg_path = era_dir / "convert_dqmio_cfg.py"
            log_path = era_dir / "cmsRun.log"
            _write_cmsrun_cfg(
                cfg_path=cfg_path,
                input_files=input_files,
                workflow=workflow,
                force_run_number=force_run_number,
                convention=convention,
            )

            emit_log(progress, f"[preprocess_dqmio] era={era}: cmsRun convert {len(input_files)} file(s)", style="bright_black")
            returncode = _run_cmsrun(
                cfg_path=cfg_path,
                cmsrun_bin=cmsrun_bin,
                era_dir=era_dir,
                log_path=log_path,
                env_command=cmsrun_env_command,
                env=cmsrun_env,
            )
            if returncode != 0:
                raise RuntimeError(
                    f"[preprocess_dqmio] era={era}: cmsRun failed (rc={returncode}). See {log_path}"
                )

            produced = _find_converted_root(era_dir, force_run_number)
            if produced is None:
                raise RuntimeError(
                    f"[preprocess_dqmio] era={era}: cmsRun succeeded but no DQM_V*.root output was found in {era_dir}"
                )
            if produced.resolve() != output_file.resolve():
                if output_file.exists():
                    output_file.unlink()
                shutil.move(str(produced), str(output_file))

        for key in ["file", "files", "file_glob", "run_files", "das", "DAS", "das_instance"]:
            era_cfg.pop(key, None)
        era_cfg["file"] = str(output_file)
        era_cfg["run_number"] = int(force_run_number)
        era_cfg["input_format"] = "legacy_root"
        era_cfg["preprocessed_from_dqmio"] = True
        emit_log(progress, f"[preprocess_dqmio] era done: {era} -> {output_file}", style="cyan")

        if progress is not None and task_id is not None:
            progress.update(task_id, advance=1)

    return cfg
