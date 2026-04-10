from itertools import product
from pathlib import Path

import ROOT
import cmsstyle as CMS
import yaml

from dqm_pipeline.core import aggregate_histogram_for_era, emit_log, sanitize


MODULE_NAME = "quantity_per_era"


COLOR_TEMPLATE = [
    ROOT.TColor.GetColor("#0072B2"),
    ROOT.TColor.GetColor("#D55E00"),
    ROOT.TColor.GetColor("#009E73"),
    ROOT.TColor.GetColor("#CC79A7"),
    ROOT.TColor.GetColor("#F0E442"),
    ROOT.TColor.GetColor("#56B4E9"),
    ROOT.TColor.GetColor("#E69F00"),
    ROOT.TColor.GetColor("#999999"),
]


def source_display_label(source):
    return str(source.get("display_label", source.get("era", "era")))


def register_canvas_input(canvas, key, obj):
    if not hasattr(canvas, "_input_cache"):
        canvas._input_cache = {}
    canvas._input_cache[str(key)] = obj
    return obj


def build_quantity_hist_name(resonance, job):
    tagging_type = str(job.get("tagging_type", "pat"))
    probe_type = str(job.get("probe_type", "sct"))
    quantity = str(job["quantity"])
    region = str(job["region"])
    return f"{resonance}_Tag_{tagging_type}_Probe_{probe_type}Electron_{quantity}_{region}"


def expand_jobs(section):
    jobs = section.get("jobs")
    if jobs:
        return jobs

    quantities = list(section.get("quantities", []))
    probe_types = list(section.get("probe_types", []))
    regions = list(section.get("regions", []))
    tagging_type = str(section.get("tagging_type", "pat"))

    expanded = []
    for quantity, probe_type, region in product(quantities, probe_types, regions):
        expanded.append(
            {
                "name": f"{quantity}_{probe_type}_{region}",
                "quantity": quantity,
                "probe_type": probe_type,
                "region": region,
                "tagging_type": tagging_type,
            }
        )
    return expanded


def plot_quantity_per_era(cfg, job, era_hists, era_sources, out_base):
    plotting = cfg.get("plotting", {})
    CMS.SetExtraText(str(plotting.get("cms_extra_text", "Preliminary")))
    CMS.SetEnergy(str(plotting.get("energy_tev", cfg.get("energy_tev", 13.6))))
    CMS.SetLumi(str(plotting.get("lumi_text", "Run3 2026")))

    quantity = str(job["quantity"])
    x_title = str(job.get("x_title", quantity))
    y_title = str(job.get("y_title", "A.U."))
    y_min = float(job.get("ymin", 0.01))

    x_low = None
    x_high = None
    ymax = 0.0

    for hist in era_hists.values():
        axis = hist.GetXaxis()
        low = axis.GetBinLowEdge(1)
        high = axis.GetBinLowEdge(axis.GetNbins() + 1)
        x_low = low if x_low is None else min(x_low, low)
        x_high = high if x_high is None else max(x_high, high)
        ymax = max(ymax, hist.GetMaximum())

    x_range = job.get("x_range")
    if x_range and len(x_range) == 2:
        x_low = float(x_range[0])
        x_high = float(x_range[1])

    canvas = CMS.cmsCanvas(
        "",
        float(x_low),
        float(x_high),
        y_min,
        max(ymax * float(job.get("ymax_scale", 5.0)), y_min * 10.0),
        x_title,
        y_title,
        square=CMS.kSquare,
        extraSpace=0.0,
        iPos=0,
    )
    register_canvas_input(canvas, "canvas", canvas)

    legend = CMS.cmsLeg(0.75, 0.70, 0.95, 0.90, textSize=float(job.get("legend_text_size", 0.02)))
    register_canvas_input(canvas, "legend", legend)

    for idx, (era_key, hist) in enumerate(era_hists.items()):
        color = COLOR_TEMPLATE[idx % len(COLOR_TEMPLATE)]
        register_canvas_input(canvas, f"hist::{era_key}", hist)
        CMS.cmsDraw(hist, "HIST", lcolor=color, msize=0, lwidth=3, fstyle=0)
        legend.AddEntry(hist, source_display_label(era_sources[era_key]), "L")

    if bool(job.get("logy", True)):
        canvas.SetLogy()

    CMS.SaveCanvas(canvas, str(out_base.with_suffix(".png")), close=False)
    CMS.SaveCanvas(canvas, str(out_base.with_suffix(".pdf")))


def run_module(cfg, era_sources, out_root, strict=False, progress=None):
    section = cfg.get(MODULE_NAME)
    if not section or not section.get("enabled", True):
        return {}

    resonance = str(cfg["resonance"])
    path_template = section.get(
        "hist_path_template",
        "DQMData/Run {run}/HLT/Run summary/ScoutingOffline/EGamma/TnP/{target_dir}/{hist}",
    )
    target_dir = str(section.get("target_dir", "Tag_PatElectron"))
    output_subdir = str(section.get("output_subdir", "quantityPlot"))
    out_dir = Path(out_root) / output_subdir
    out_dir.mkdir(parents=True, exist_ok=True)

    default_rebin = int(section.get("rebin", 2))
    default_scale_to = float(section.get("scale_to", 100000.0))
    quantity_overrides = section.get("quantity_overrides", {})

    jobs = expand_jobs(section)
    emit_log(progress, f"[{MODULE_NAME}] start jobs={len(jobs)} eras={len(era_sources)} out_dir={out_dir}", style="blue")

    summary = {}
    job_task = None
    if progress is not None:
        job_task = progress.add_task(f"{MODULE_NAME}: jobs", total=len(jobs))

    for job in jobs:
        job = dict(job)
        job.update(quantity_overrides.get(str(job.get("quantity")), {}))
        job_name = str(job.get("name", f"{job['quantity']}_{job['probe_type']}_{job['region']}"))
        emit_log(progress, f"[{MODULE_NAME}] job start: {job_name}", style="blue")
        hist_name = build_quantity_hist_name(resonance, job)
        era_hists = {}
        summary[job_name] = {}

        era_task = None
        era_items = list(era_sources.items())
        if progress is not None:
            era_task = progress.add_task(f"{MODULE_NAME}: {job_name} eras", total=len(era_items))

        for era_key, source in era_items:
            try:
                hist, used_runs = aggregate_histogram_for_era(
                    era=era_key,
                    source=source,
                    hist_path_template=path_template,
                    fmt_args={"target_dir": target_dir, "hist": hist_name},
                    strict=strict,
                )
                if hist is None:
                    summary[job_name][era_key] = {"status": "empty", "used_runs": 0}
                    continue

                integral = hist.Integral()
                if integral <= 0:
                    summary[job_name][era_key] = {"status": "empty", "used_runs": used_runs, "integral": float(integral)}
                    continue

                hist = hist.Clone(f"{sanitize(job_name)}_{sanitize(era_key)}_{sanitize(hist.GetName())}")
                hist.SetDirectory(0)
                hist.Scale(default_scale_to / integral)

                x_range = job.get("x_range")
                if x_range and len(x_range) == 2:
                    hist.GetXaxis().SetRangeUser(float(x_range[0]), float(x_range[1]))

                rebin_factor = int(job.get("rebin", default_rebin))
                if rebin_factor > 1:
                    hist = hist.Rebin(rebin_factor, f"{hist.GetName()}_rebin")
                    hist.SetDirectory(0)

                era_hists[era_key] = hist
                summary[job_name][era_key] = {
                    "status": "ok",
                    "used_runs": used_runs,
                    "integral_before_scale": float(integral),
                }
            except Exception as exc:
                summary[job_name][era_key] = {"status": "error", "message": str(exc)}
                print(f"[{MODULE_NAME}][WARN] job={job_name} era={era_key}: {exc}")
                if strict:
                    raise
            finally:
                if progress is not None and era_task is not None:
                    progress.update(era_task, advance=1)

        if era_hists:
            ordered_hists = {era_key: era_hists[era_key] for era_key in era_sources if era_key in era_hists}
            out_base = out_dir / sanitize(hist_name)
            plot_quantity_per_era(cfg, job, ordered_hists, era_sources, out_base)
            usable = len(ordered_hists)
            skipped = len(era_sources) - usable
            emit_log(progress, f"[{MODULE_NAME}] job done: {job_name} usable_eras={usable} skipped_eras={skipped}", style="blue")
        else:
            emit_log(progress, f"[{MODULE_NAME}][WARN] job done: {job_name} usable_eras=0 skipped_eras={len(era_sources)}", style="yellow")

        if progress is not None and job_task is not None:
            progress.update(job_task, advance=1)

    summary_file = out_dir / f"{MODULE_NAME}_summary.yaml"
    with open(summary_file, "w", encoding="utf-8") as handle:
        yaml.safe_dump(summary, handle, sort_keys=False)

    emit_log(progress, f"[{MODULE_NAME}] done. outputs={out_dir}", style="blue")
    return summary
