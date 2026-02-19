import json
import hashlib
import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime
import time
from array import array
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from urllib.parse import urlparse
from urllib.request import urlopen

import ROOT


DEFAULT_XROOTD_REDIRECTOR = "root://cms-xrd-global.cern.ch"


def sanitize(name):
    safe = name
    for bad in ["/", " ", "[", "]", "(", ")", ":", "="]:
        safe = safe.replace(bad, "_")
    while "__" in safe:
        safe = safe.replace("__", "_")
    return safe.strip("_")


def resolve_path(path_str, base_dir):
    path = Path(path_str)
    if path.is_absolute():
        return str(path)
    return str((base_dir / path).resolve())


def is_http_url(value):
    parsed = urlparse(str(value))
    return parsed.scheme in ("http", "https")


def as_root_uri(path_str, redirector):
    if path_str.startswith("root://"):
        return path_str
    if path_str.startswith("/store/"):
        return f"{redirector.rstrip('/')}//{path_str.lstrip('/')}"
    return path_str


def run_command(cmd, progress=None, description=None, env=None, timeout=300):
    task_id = None
    started = time.time()
    if description:
        emit_log(progress, f"[cmd] start: {description}", style="bright_black")
    if progress is not None and description:
        task_id = progress.add_task(f"[yellow]{description}", total=None)

    # Use temporary files instead of PIPE to avoid deadlocks on large command output.
    with tempfile.TemporaryFile() as stdout_tmp, tempfile.TemporaryFile() as stderr_tmp:
        proc = subprocess.Popen(
            cmd,
            stdout=stdout_tmp,
            stderr=stderr_tmp,
            env=env,
        )

        heartbeat_sec = 20
        last_heartbeat = time.time()
        while True:
            rc = proc.poll()
            if rc is not None:
                break

            now = time.time()
            if now - started > timeout:
                proc.kill()
                proc.wait()
                stdout_tmp.seek(0)
                stderr_tmp.seek(0)
                stdout = stdout_tmp.read().decode("utf-8", errors="replace")
                stderr = stderr_tmp.read().decode("utf-8", errors="replace")
                if progress is not None and task_id is not None:
                    progress.remove_task(task_id)
                raise subprocess.TimeoutExpired(cmd=cmd, timeout=timeout, output=stdout, stderr=stderr)

            if description and (now - last_heartbeat) >= heartbeat_sec:
                elapsed = int(now - started)
                emit_log(progress, f"[cmd] running: {description} elapsed={elapsed}s", style="bright_black")
                if progress is not None and task_id is not None:
                    progress.update(task_id, description=f"[yellow]{description} (elapsed {elapsed}s)")
                last_heartbeat = now
            time.sleep(0.2)

        proc.wait()
        stdout_tmp.seek(0)
        stderr_tmp.seek(0)
        stdout = stdout_tmp.read().decode("utf-8", errors="replace")
        stderr = stderr_tmp.read().decode("utf-8", errors="replace")

    if progress is not None and task_id is not None:
        progress.remove_task(task_id)
    elapsed = time.time() - started
    if description:
        emit_log(progress, f"[cmd] done: {description} ({elapsed:.1f}s, rc={proc.returncode})", style="bright_black")
    if proc.returncode != 0:
        raise subprocess.CalledProcessError(proc.returncode, cmd, output=stdout, stderr=stderr)
    return stdout, stderr


def query_das_files(dataset, instance=None, progress=None):
    if shutil.which("dasgoclient") is None:
        raise RuntimeError("dasgoclient not found in PATH. Needed for DAS dataset resolution.")

    query = f"file dataset={dataset}"
    if instance:
        query += f" instance={instance}"
    emit_log(progress, f"[resolve] DAS query: {query}", style="bright_black")

    stdout, _ = run_command(
        ["dasgoclient", "--query", query],
        progress=progress,
        description=f"dasgoclient file query ({dataset})",
        timeout=180,
    )

    files = [line.strip() for line in stdout.splitlines() if line.strip()]
    if not files:
        raise RuntimeError(f"No files found from DAS query: {query}")
    emit_log(progress, f"[resolve] DAS returned {len(files)} files", style="bright_black")
    return files


def query_run_for_file(file_pfn, instance=None):
    if shutil.which("dasgoclient") is None:
        return None

    query = f"run file={file_pfn}"
    if instance:
        query += f" instance={instance}"

    try:
        proc = subprocess.run(
            ["dasgoclient", "--query", query],
            capture_output=True,
            text=True,
            check=True,
            timeout=60,
        )
    except Exception:
        return None

    stdout = proc.stdout
    vals = [line.strip() for line in stdout.splitlines() if line.strip()]
    if not vals:
        return None

    try:
        return int(vals[0])
    except ValueError:
        return None


def extract_run_from_name(path_or_name):
    name = Path(path_or_name).name
    patterns = [r"_R0*([0-9]{5,})__", r"Run0*([0-9]{5,})", r"run0*([0-9]{5,})"]
    for pat in patterns:
        match = re.search(pat, name)
        if match:
            return int(match.group(1))
    return None


def count_lumisections(ls_ranges):
    total = 0
    for pair in ls_ranges:
        if not isinstance(pair, list) or len(pair) != 2:
            continue
        lo = int(pair[0])
        hi = int(pair[1])
        if hi >= lo:
            total += hi - lo + 1
    return total


def load_golden_json(golden_json_source, cache, progress=None):
    if golden_json_source in cache:
        cached = cache[golden_json_source]
        return cached["run_to_ls"], cached["local_json_path"]

    if is_http_url(golden_json_source):
        emit_log(progress, f"[resolve] Downloading golden JSON from URL: {golden_json_source}", style="bright_black")
        with urlopen(golden_json_source, timeout=120) as resp:
            raw = resp.read().decode("utf-8")
        data = json.loads(raw)

        digest = hashlib.sha256(golden_json_source.encode("utf-8")).hexdigest()[:16]
        local_json_path = str(Path(tempfile.gettempdir()) / f"dqm_golden_{digest}.json")
        with open(local_json_path, "w", encoding="utf-8") as out:
            json.dump(data, out)
        emit_log(progress, f"[resolve] Cached golden JSON to: {local_json_path}", style="bright_black")
    else:
        emit_log(progress, f"[resolve] Loading golden JSON from local file: {golden_json_source}", style="bright_black")
        with open(golden_json_source, "r", encoding="utf-8") as f:
            data = json.load(f)
        local_json_path = golden_json_source

    run_to_ls = {}
    for run_str, ls_ranges in data.items():
        try:
            run = int(run_str)
        except ValueError:
            continue
        run_to_ls[run] = count_lumisections(ls_ranges)

    cache[golden_json_source] = {
        "run_to_ls": run_to_ls,
        "local_json_path": local_json_path,
    }
    emit_log(progress, f"[resolve] Golden JSON contains {len(run_to_ls)} runs", style="bright_black")
    return run_to_ls, local_json_path


def parse_brilcalc_csv(stdout_text, selected_runs):
    selected = set(int(x) for x in selected_runs)
    recorded_idx = None
    total = 0.0

    for raw in stdout_text.splitlines():
        line = raw.strip()
        if not line:
            continue

        if line.startswith("#"):
            line = line.lstrip("#").strip()
        if not line:
            continue

        fields = [x.strip() for x in line.split(",")]
        if not fields:
            continue

        if "run:fill" in fields[0].lower():
            for idx, name in enumerate(fields):
                if name.lower().startswith("recorded("):
                    recorded_idx = idx
                    break
            continue

        if ":" not in fields[0]:
            continue

        try:
            run = int(fields[0].split(":")[0])
        except ValueError:
            continue

        if run not in selected:
            continue

        if recorded_idx is None:
            for idx in range(len(fields) - 1, -1, -1):
                try:
                    float(fields[idx])
                    recorded_idx = idx
                    break
                except ValueError:
                    continue

        if recorded_idx is None or recorded_idx >= len(fields):
            continue

        try:
            total += float(fields[recorded_idx])
        except ValueError:
            continue

    return total


def convert_lumi_to_fb(value, unit):
    scale_to_fb = {
        "/fb": 1.0,
        "/pb": 1.0e-3,
        "/nb": 1.0e-6,
        "/ub": 1.0e-9,
        "/mb": 1.0e-12,
    }
    if unit not in scale_to_fb:
        raise RuntimeError(f"Unsupported brilcalc unit '{unit}'. Use one of: {sorted(scale_to_fb.keys())}")
    return value * scale_to_fb[unit]


def build_brilcalc_env(lumi_cfg):
    # Run brilcalc in a cleaned subprocess env to avoid clashes with LCG/venv python vars.
    if not bool(lumi_cfg.get("clean_env", True)):
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
    }
    for key in unset_defaults:
        env.pop(key, None)

    for key in lumi_cfg.get("unset_env_vars", []):
        env.pop(str(key), None)

    if not env.get("PATH"):
        env["PATH"] = "/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"
    return env


def estimate_lumi_fb_with_brilcalc(era, selected_runs, golden_json_path, lumi_cfg, strict=False, progress=None):
    if not selected_runs:
        return None

    brilcalc_bin = lumi_cfg.get("brilcalc_bin", "brilcalc")
    brilcalc_env = lumi_cfg.get("brilcalc_env")
    brilcalc_timeout = int(lumi_cfg.get("timeout_sec", 600))
    unit = lumi_cfg.get("unit", "/fb")
    calibration = lumi_cfg.get("calibration", "web")
    normtag = lumi_cfg.get("normtag")
    run_env = build_brilcalc_env(lumi_cfg)

    begin_run = min(selected_runs)
    end_run = max(selected_runs)

    cmd = [
        brilcalc_bin,
        "lumi",
        "-u",
        unit,
        "-c",
        calibration,
        "--begin",
        str(begin_run),
        "--end",
        str(end_run),
        "-i",
        golden_json_path,
        "--output-style",
        "csv",
    ]
    if normtag:
        cmd.extend(["--normtag", normtag])
    cmd_text = " ".join(shlex.quote(x) for x in cmd)
    print(f"[lumi] era={era} running brilcalc: {cmd_text}")

    proc = None
    stdout = ""
    stderr = ""
    source_attempt_error = None

    # Attempt 1: run with optional brilws env sourcing from config.
    try:
        if brilcalc_env:
            cmd_text = " ".join(shlex.quote(x) for x in cmd)
            shell_cmd = (
                "set -e\n"
                "shopt -s expand_aliases\n"
                f"source {shlex.quote(brilcalc_env)}\n"
                f"eval {shlex.quote(cmd_text)}\n"
            )
            stdout, stderr = run_command(
                ["bash", "--noprofile", "--norc", "-lc", shell_cmd],
                progress=progress,
                description=f"brilcalc (era={era}, sourced env)",
                timeout=brilcalc_timeout,
                env=run_env,
            )
            proc = True
        else:
            stdout, stderr = run_command(
                cmd,
                progress=progress,
                description=f"brilcalc (era={era})",
                timeout=brilcalc_timeout,
                env=run_env,
            )
            proc = True
    except subprocess.CalledProcessError as exc:
        source_attempt_error = exc
    except Exception as exc:
        source_attempt_error = exc

    # Attempt 2: if source-based attempt failed with "command not found", retry current shell env.
    if proc is None and isinstance(source_attempt_error, subprocess.CalledProcessError):
        if source_attempt_error.returncode == 127 and brilcalc_env:
            try:
                stdout, stderr = run_command(
                    cmd,
                    progress=progress,
                    description=f"brilcalc retry without sourcing (era={era})",
                    timeout=brilcalc_timeout,
                    env=run_env,
                )
                proc = True
            except Exception:
                pass

    if proc is None:
        err_text = ""
        if isinstance(source_attempt_error, subprocess.CalledProcessError):
            stderr_txt = (source_attempt_error.stderr or "").strip()
            if stderr_txt:
                err_text = f" stderr={stderr_txt}"
            if source_attempt_error.returncode == 127:
                err_text += " hint=brilcalc not found in shell PATH (or after sourcing brilws-env)."
        msg = f"[lumi][WARN] era={era}: brilcalc failed ({source_attempt_error}){err_text}"
        if strict:
            raise RuntimeError(msg) from source_attempt_error
        print(msg)
        return None

    recorded_in_unit = parse_brilcalc_csv(stdout, selected_runs)
    if recorded_in_unit <= 0.0:
        msg = f"[lumi][WARN] era={era}: brilcalc returned zero recorded lumi in selected runs."
        if strict:
            raise RuntimeError(msg)
        print(msg)
        return None

    lumi_fb = convert_lumi_to_fb(recorded_in_unit, unit)
    print(f"[lumi] era={era} recorded={recorded_in_unit:.6g} {unit} => {lumi_fb:.6g} /fb")
    return lumi_fb


def run_passes_requirement(run, req):
    if not req:
        return True

    if "before" in req and run >= int(req["before"]):
        return False
    if "after" in req and run <= int(req["after"]):
        return False
    if "min" in req and run < int(req["min"]):
        return False
    if "max" in req and run > int(req["max"]):
        return False

    include = req.get("include")
    if include is not None:
        include_set = {int(x) for x in include}
        if run not in include_set:
            return False

    exclude = req.get("exclude")
    if exclude is not None:
        exclude_set = {int(x) for x in exclude}
        if run in exclude_set:
            return False

    return True


def _ansi_wrap(text, style):
    if os.environ.get("NO_COLOR"):
        return text
    if not sys.stdout.isatty():
        return text

    style_map = {
        "yellow": "\033[33m",
        "cyan": "\033[36m",
        "green": "\033[32m",
        "blue": "\033[34m",
        "magenta": "\033[35m",
        "bright_black": "\033[90m",
        "bold blue": "\033[1;34m",
        "bold green": "\033[1;32m",
        "bold cyan": "\033[1;36m",
        "bold magenta": "\033[1;35m",
    }
    prefix = style_map.get(style)
    if not prefix:
        return text
    return f"{prefix}{text}\033[0m"


def emit_log(progress, message, style=None):
    stamp = datetime.now().strftime("%H:%M:%S")
    full_message = f"[{stamp}] {message}"
    if progress is not None and hasattr(progress, "console"):
        progress.console.print(full_message, style=style)
    else:
        print(_ansi_wrap(full_message, style))


def _emit(progress, message, style=None):
    emit_log(progress, message, style=style)


def _run_range_text(run_files):
    if not run_files:
        return "n/a"
    runs = sorted(run_files.keys())
    return f"{runs[0]}-{runs[-1]}"


def resolve_era_source(era, era_cfg, cfg, config_dir, golden_cache, strict=False, progress=None):
    redirector = era_cfg.get("xrootd_redirector", cfg.get("xrootd_redirector", DEFAULT_XROOTD_REDIRECTOR))
    run_req = era_cfg.get("run_requirement", era_cfg.get("run-requirement"))
    # Optional: if missing/null, no run-based filtering is applied.
    if run_req is None:
        run_req = {}

    run_files = defaultdict(list)
    stage_task = None
    if progress is not None:
        stage_task = progress.add_task(f"[white]{era}: resolve inputs", total=4)

    # Legacy single-file mode (kept for backward compatibility).
    if "file" in era_cfg:
        file_path = era_cfg["file"]
        run = extract_run_from_name(file_path)
        if run is None:
            legacy_run = cfg.get("run_number")
            if legacy_run is None:
                raise RuntimeError(
                    f"Era {era}: cannot infer run from file name and no top-level run_number provided."
                )
            run = int(legacy_run)
        run_files[run].append(as_root_uri(file_path, redirector))

    if "run_files" in era_cfg:
        for run_key, file_list in era_cfg["run_files"].items():
            run = int(run_key)
            values = file_list if isinstance(file_list, list) else [file_list]
            for file_path in values:
                run_files[run].append(as_root_uri(file_path, redirector))

    if "files" in era_cfg:
        for file_path in era_cfg["files"]:
            run = extract_run_from_name(file_path)
            if run is None:
                msg = f"Era {era}: cannot infer run from file '{file_path}'."
                if strict:
                    raise RuntimeError(msg)
                _emit(progress, f"[resolve][WARN] {msg} Skipping file.", style="yellow")
                continue
            run_files[run].append(as_root_uri(file_path, redirector))

    dataset = era_cfg.get("das", era_cfg.get("DAS"))
    das_instance = era_cfg.get("das_instance")
    if dataset:
        if progress is not None and stage_task is not None:
            progress.update(stage_task, description=f"[white]{era}: querying DAS")
        das_files = query_das_files(dataset=dataset, instance=das_instance, progress=progress)
        known_mapped = 0
        unresolved_pfns = []

        for pfn in das_files:
            run = extract_run_from_name(pfn)
            if run is not None:
                run_files[run].append(as_root_uri(pfn, redirector))
                known_mapped += 1
            else:
                unresolved_pfns.append(pfn)

        _emit(
            progress,
            f"[resolve] {era}: direct run parsing mapped {known_mapped}/{len(das_files)} files; querying DAS for {len(unresolved_pfns)} unresolved files",
            style="bright_black",
        )

        unresolved_count = 0
        if unresolved_pfns:
            workers = int(era_cfg.get("das_run_query_workers", cfg.get("das_run_query_workers", 8)))
            workers = max(1, workers)
            map_task = None
            if progress is not None:
                map_task = progress.add_task(
                    f"[yellow]{era}: map files -> runs (0/{len(unresolved_pfns)})",
                    total=len(unresolved_pfns),
                )
            last_report = time.time()
            processed = 0

            with ThreadPoolExecutor(max_workers=workers) as pool:
                future_map = {
                    pool.submit(query_run_for_file, pfn, das_instance): pfn
                    for pfn in unresolved_pfns
                }
                for fut in as_completed(future_map):
                    pfn = future_map[fut]
                    run = None
                    try:
                        run = fut.result()
                    except Exception:
                        run = None

                    processed += 1
                    if run is None:
                        unresolved_count += 1
                        msg = f"Era {era}: cannot resolve run for DAS file '{pfn}'."
                        if strict:
                            raise RuntimeError(msg)
                    else:
                        run_files[run].append(as_root_uri(pfn, redirector))

                    if progress is not None and map_task is not None:
                        if processed == 1 or processed == len(unresolved_pfns) or processed % 10 == 0:
                            progress.update(
                                map_task,
                                description=f"[yellow]{era}: map files -> runs ({processed}/{len(unresolved_pfns)})",
                            )
                        progress.update(map_task, advance=1)
                    else:
                        now = time.time()
                        if (
                            processed == 1
                            or processed == len(unresolved_pfns)
                            or (now - last_report) >= 5.0
                        ):
                            pct = 100.0 * processed / len(unresolved_pfns)
                            _emit(
                                progress,
                                f"[resolve] {era}: mapping unresolved files {processed}/{len(unresolved_pfns)} ({pct:.1f}%)",
                                style="bright_black",
                            )
                            last_report = now

            if unresolved_count > 0:
                _emit(
                    progress,
                    f"[resolve][WARN] {era}: {unresolved_count} DAS files still unresolved and were skipped",
                    style="yellow",
                )

        _emit(
            progress,
            f"[resolve] {era}: mapped {len(das_files) - unresolved_count}/{len(das_files)} DAS files to runs",
            style="bright_black",
        )

    if not run_files:
        raise RuntimeError(f"Era {era}: no inputs found. Provide one of: file / files / run_files / DAS(das).")
    if progress is not None and stage_task is not None:
        progress.update(stage_task, advance=1, description=f"[white]{era}: applying filters")

    golden_json_source = era_cfg.get("golden_json", cfg.get("golden_json"))
    golden_runs = None
    run_to_ls = {}
    golden_json_for_brilcalc = None
    if golden_json_source:
        if not is_http_url(golden_json_source):
            golden_json_source = resolve_path(golden_json_source, config_dir)
        run_to_ls, golden_json_for_brilcalc = load_golden_json(golden_json_source, golden_cache, progress=progress)
        golden_runs = set(run_to_ls.keys())

    selected_run_files = {}
    excluded_by_run_requirement = 0
    excluded_by_golden = 0
    for run in sorted(run_files.keys()):
        if not run_passes_requirement(run, run_req):
            excluded_by_run_requirement += 1
            continue
        if golden_runs is not None and run not in golden_runs:
            excluded_by_golden += 1
            continue
        selected_run_files[run] = sorted(set(run_files[run]))

    if not selected_run_files:
        raise RuntimeError(f"Era {era}: no runs left after run_requirement/golden_json filters.")
    if progress is not None and stage_task is not None:
        progress.update(stage_task, advance=1, description=f"[white]{era}: lumi estimation")

    selected_lumisections = 0
    if run_to_ls:
        selected_lumisections = sum(run_to_ls.get(run, 0) for run in selected_run_files)

    lumi_fb = era_cfg.get("lumi_fb")
    lumi_source = "manual" if lumi_fb is not None else "none"

    lumi_cfg = cfg.get("lumi", {})
    use_brilcalc = bool(lumi_cfg.get("use_brilcalc", False) or era_cfg.get("lumi_from_brilcalc", False))
    if lumi_fb is None and use_brilcalc:
        if not golden_json_for_brilcalc:
            msg = f"[lumi][WARN] era={era}: use_brilcalc=true but no golden_json provided."
            if strict:
                raise RuntimeError(msg)
            _emit(progress, msg, style="yellow")
        else:
            lumi_fb = estimate_lumi_fb_with_brilcalc(
                era=era,
                selected_runs=sorted(selected_run_files.keys()),
                golden_json_path=golden_json_for_brilcalc,
                lumi_cfg=lumi_cfg,
                strict=strict,
                progress=progress,
            )
            if lumi_fb is not None:
                lumi_source = "brilcalc"
    if progress is not None and stage_task is not None:
        progress.update(stage_task, advance=1, description=f"[white]{era}: finalizing")

    out = {
        "era": era,
        "run_files": selected_run_files,
        "n_runs": len(selected_run_files),
        "n_runs_discovered": len(run_files),
        "n_runs_excluded_by_run_requirement": excluded_by_run_requirement,
        "n_runs_excluded_by_golden": excluded_by_golden,
        "n_files_discovered": sum(len(v) for v in run_files.values()),
        "n_files_selected": sum(len(v) for v in selected_run_files.values()),
        "run_range_discovered": _run_range_text(run_files),
        "run_range_selected": _run_range_text(selected_run_files),
        "selected_lumisections": selected_lumisections,
        "lumi_fb": lumi_fb,
        "lumi_source": lumi_source,
        "dataset": dataset,
        "golden_json_source": golden_json_source,
    }
    if progress is not None and stage_task is not None:
        progress.update(stage_task, advance=1, description=f"[white]{era}: done")
    return out


def prepare_era_sources(cfg, config_dir, strict=False, progress=None):
    golden_cache = {}
    out = {}
    era_items = list(cfg["eras"].items())

    task_id = None
    if progress is not None:
        task_id = progress.add_task("[cyan]Resolving eras", total=len(era_items))

    total_runs_selected = 0
    total_files_selected = 0

    for era, era_cfg in era_items:
        if progress is not None and task_id is not None:
            progress.update(task_id, description=f"[cyan]Resolving eras ({era})")

        source = resolve_era_source(
            era=era,
            era_cfg=era_cfg,
            cfg=cfg,
            config_dir=config_dir,
            golden_cache=golden_cache,
            strict=strict,
            progress=progress,
        )
        out[era] = source
        total_runs_selected += int(source["n_runs"])
        total_files_selected += int(source["n_files_selected"])

        if source["lumi_fb"] is not None:
            norm_text = f"lumi_fb={float(source['lumi_fb']):.3f} ({source.get('lumi_source', 'unknown')})"
        elif source["selected_lumisections"] > 0:
            norm_text = f"golden_lumisections={source['selected_lumisections']}"
        else:
            norm_text = "no lumi normalization"

        _emit(
            progress,
            (
                f"[resolve] {era}: runs {source['n_runs_discovered']} -> {source['n_runs']} "
                f"(run-req -{source['n_runs_excluded_by_run_requirement']}, golden -{source['n_runs_excluded_by_golden']}), "
                f"files {source['n_files_discovered']} -> {source['n_files_selected']}"
            ),
            style="cyan",
        )
        _emit(
            progress,
            (
                f"          run range: {source['run_range_selected']} "
                f"(discovered {source['run_range_discovered']}) | "
                f"dataset: {source.get('dataset') or 'local files'} | {norm_text}"
            ),
            style="green",
        )
        if source.get("golden_json_source"):
            _emit(
                progress,
                f"          golden_json: {source['golden_json_source']}",
                style="bright_black",
            )

        if progress is not None and task_id is not None:
            progress.update(task_id, advance=1)

    _emit(
        progress,
        f"[resolve][summary] eras={len(era_items)}, selected runs={total_runs_selected}, selected files={total_files_selected}",
        style="bold blue",
    )

    return out


def load_histogram(file_path, hist_path):
    root_file = ROOT.TFile.Open(file_path)
    if not root_file or root_file.IsZombie():
        raise RuntimeError(f"Cannot open file: {file_path}")

    hist = root_file.Get(hist_path)
    if not hist:
        root_file.Close()
        raise RuntimeError(f"Missing histogram '{hist_path}' in {file_path}")

    out = hist.Clone()
    out.SetDirectory(0)
    root_file.Close()
    return out


def aggregate_histogram_for_era(era, source, hist_path_template, fmt_args, strict=False):
    merged = None
    used_runs = []
    runs = sorted(source["run_files"].keys())
    files_total = sum(len(v) for v in source["run_files"].values())
    print(
        f"[aggregate] era={era} start runs={len(runs)} files={files_total} template='{hist_path_template}' args={fmt_args}"
    )

    for i_run, run in enumerate(runs, start=1):
        hist_path = hist_path_template.format(run=run, **fmt_args)
        run_has_hist = False

        for file_path in source["run_files"][run]:
            try:
                hist = load_histogram(file_path, hist_path)
            except Exception as exc:
                print(f"[aggregate][WARN] era={era} run={run}: {exc}")
                if strict:
                    raise
                continue

            if merged is None:
                merged = hist.Clone(f"agg_{sanitize(era)}_{sanitize(hist.GetName())}")
                merged.SetDirectory(0)
            else:
                merged.Add(hist)
            run_has_hist = True

        if run_has_hist:
            used_runs.append(run)
        if i_run == 1 or i_run % 25 == 0 or i_run == len(runs):
            print(f"[aggregate] era={era} progress runs={i_run}/{len(runs)} used_runs={len(used_runs)}")
    print(f"[aggregate] era={era} done used_runs={len(used_runs)}/{len(runs)}")
    return merged, used_runs


def rebin_histogram(hist, bins=None, rebin_factor=1, name_hint="hist"):
    rebinned = hist.Clone(f"{name_hint}_clone")
    rebinned.SetDirectory(0)

    if bins:
        edges = array("d", [float(x) for x in bins])
        custom = rebinned.Rebin(len(edges) - 1, f"{name_hint}_custom", edges)
        custom.SetDirectory(0)
        return custom

    if int(rebin_factor) > 1:
        rebinned.Rebin(int(rebin_factor))
    return rebinned


def build_output_root(cfg):
    return Path(cfg["output_dir"]) / cfg.get("tag", cfg["resonance"])
