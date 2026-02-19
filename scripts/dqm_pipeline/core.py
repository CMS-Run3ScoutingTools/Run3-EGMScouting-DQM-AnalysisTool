import json
import hashlib
import re
import shlex
import shutil
import subprocess
import tempfile
from array import array
from collections import defaultdict
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


def query_das_files(dataset, instance=None):
    if shutil.which("dasgoclient") is None:
        raise RuntimeError("dasgoclient not found in PATH. Needed for DAS dataset resolution.")

    query = f"file dataset={dataset}"
    if instance:
        query += f" instance={instance}"

    proc = subprocess.run(
        ["dasgoclient", "--query", query],
        capture_output=True,
        text=True,
        check=True,
        timeout=180,
    )

    files = [line.strip() for line in proc.stdout.splitlines() if line.strip()]
    if not files:
        raise RuntimeError(f"No files found from DAS query: {query}")
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

    vals = [line.strip() for line in proc.stdout.splitlines() if line.strip()]
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


def load_golden_json(golden_json_source, cache):
    if golden_json_source in cache:
        cached = cache[golden_json_source]
        return cached["run_to_ls"], cached["local_json_path"]

    if is_http_url(golden_json_source):
        with urlopen(golden_json_source, timeout=120) as resp:
            raw = resp.read().decode("utf-8")
        data = json.loads(raw)

        digest = hashlib.sha256(golden_json_source.encode("utf-8")).hexdigest()[:16]
        local_json_path = str(Path(tempfile.gettempdir()) / f"dqm_golden_{digest}.json")
        with open(local_json_path, "w", encoding="utf-8") as out:
            json.dump(data, out)
    else:
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


def estimate_lumi_fb_with_brilcalc(era, selected_runs, golden_json_path, lumi_cfg, strict=False):
    if not selected_runs:
        return None

    brilcalc_bin = lumi_cfg.get("brilcalc_bin", "brilcalc")
    brilcalc_env = lumi_cfg.get("brilcalc_env")
    unit = lumi_cfg.get("unit", "/fb")
    calibration = lumi_cfg.get("calibration", "web")
    normtag = lumi_cfg.get("normtag")

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

    try:
        if brilcalc_env:
            shell_cmd = f"source {shlex.quote(brilcalc_env)} && " + " ".join(shlex.quote(x) for x in cmd)
            proc = subprocess.run(
                ["bash", "-lc", shell_cmd],
                capture_output=True,
                text=True,
                check=True,
                timeout=300,
            )
        else:
            proc = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True,
                timeout=300,
            )
    except Exception as exc:
        msg = f"[lumi][WARN] era={era}: brilcalc failed ({exc})"
        if strict:
            raise RuntimeError(msg) from exc
        print(msg)
        return None

    recorded_in_unit = parse_brilcalc_csv(proc.stdout, selected_runs)
    if recorded_in_unit <= 0.0:
        msg = f"[lumi][WARN] era={era}: brilcalc returned zero recorded lumi in selected runs."
        if strict:
            raise RuntimeError(msg)
        print(msg)
        return None

    return convert_lumi_to_fb(recorded_in_unit, unit)


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


def resolve_era_source(era, era_cfg, cfg, config_dir, golden_cache, strict=False):
    redirector = era_cfg.get("xrootd_redirector", cfg.get("xrootd_redirector", DEFAULT_XROOTD_REDIRECTOR))
    run_req = era_cfg.get("run_requirement", era_cfg.get("run-requirement"))
    # Optional: if missing/null, no run-based filtering is applied.
    if run_req is None:
        run_req = {}

    run_files = defaultdict(list)

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
                print(f"[resolve][WARN] {msg} Skipping file.")
                continue
            run_files[run].append(as_root_uri(file_path, redirector))

    dataset = era_cfg.get("das", era_cfg.get("DAS"))
    das_instance = era_cfg.get("das_instance")
    if dataset:
        das_files = query_das_files(dataset=dataset, instance=das_instance)
        for pfn in das_files:
            run = extract_run_from_name(pfn)
            if run is None:
                run = query_run_for_file(pfn, instance=das_instance)
            if run is None:
                msg = f"Era {era}: cannot resolve run for DAS file '{pfn}'."
                if strict:
                    raise RuntimeError(msg)
                print(f"[resolve][WARN] {msg} Skipping file.")
                continue
            run_files[run].append(as_root_uri(pfn, redirector))

    if not run_files:
        raise RuntimeError(f"Era {era}: no inputs found. Provide one of: file / files / run_files / DAS(das).")

    golden_json_source = era_cfg.get("golden_json", cfg.get("golden_json"))
    golden_runs = None
    run_to_ls = {}
    golden_json_for_brilcalc = None
    if golden_json_source:
        if not is_http_url(golden_json_source):
            golden_json_source = resolve_path(golden_json_source, config_dir)
        run_to_ls, golden_json_for_brilcalc = load_golden_json(golden_json_source, golden_cache)
        golden_runs = set(run_to_ls.keys())

    selected_run_files = {}
    for run in sorted(run_files.keys()):
        if not run_passes_requirement(run, run_req):
            continue
        if golden_runs is not None and run not in golden_runs:
            continue
        selected_run_files[run] = sorted(set(run_files[run]))

    if not selected_run_files:
        raise RuntimeError(f"Era {era}: no runs left after run_requirement/golden_json filters.")

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
            print(msg)
        else:
            lumi_fb = estimate_lumi_fb_with_brilcalc(
                era=era,
                selected_runs=sorted(selected_run_files.keys()),
                golden_json_path=golden_json_for_brilcalc,
                lumi_cfg=lumi_cfg,
                strict=strict,
            )
            if lumi_fb is not None:
                lumi_source = "brilcalc"

    return {
        "run_files": selected_run_files,
        "n_runs": len(selected_run_files),
        "selected_lumisections": selected_lumisections,
        "lumi_fb": lumi_fb,
        "lumi_source": lumi_source,
        "dataset": dataset,
    }


def prepare_era_sources(cfg, config_dir, strict=False, progress=None):
    golden_cache = {}
    out = {}
    era_items = list(cfg["eras"].items())

    task_id = None
    if progress is not None:
        task_id = progress.add_task("[cyan]Resolving eras", total=len(era_items))

    for era, era_cfg in era_items:
        source = resolve_era_source(
            era=era,
            era_cfg=era_cfg,
            cfg=cfg,
            config_dir=config_dir,
            golden_cache=golden_cache,
            strict=strict,
        )
        out[era] = source

        if source["lumi_fb"] is not None:
            norm_text = f"lumi_fb={float(source['lumi_fb']):.3f} ({source.get('lumi_source', 'unknown')})"
        elif source["selected_lumisections"] > 0:
            norm_text = f"golden_lumisections={source['selected_lumisections']}"
        else:
            norm_text = "no lumi normalization"

        print(f"[resolve] {era}: runs={source['n_runs']}, {norm_text}")
        if progress is not None and task_id is not None:
            progress.update(task_id, advance=1)

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

    for run in sorted(source["run_files"].keys()):
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
