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


def apply_poisson_errors_from_counts(hist, count_hist, scale_factor=1.0):
    for idx in range(1, hist.GetNbinsX() + 1):
        raw_value = float(count_hist.GetBinContent(idx))
        error = (raw_value ** 0.5) * float(scale_factor) if raw_value > 0 else 0.0
        hist.SetBinError(idx, error)


def make_error_band(hist, name, color):
    band = hist.Clone(name)
    band.SetDirectory(0)
    band.SetLineColor(color)
    band.SetFillStyle(0)
    band.SetMarkerSize(0)
    band.SetLineWidth(2)
    band.SetLineStyle(2)
    return band


def find_populated_xrange(era_hists, margin_bins=1):
    first_edge = None
    last_edge = None
    for hist in era_hists.values():
        axis = hist.GetXaxis()
        first_bin = None
        last_bin = None
        for idx in range(1, hist.GetNbinsX() + 1):
            if hist.GetBinContent(idx) > 0:
                if first_bin is None:
                    first_bin = idx
                last_bin = idx
        if first_bin is None or last_bin is None:
            continue
        first_bin = max(1, int(first_bin) - int(margin_bins))
        last_bin = min(hist.GetNbinsX(), int(last_bin) + int(margin_bins))
        low = axis.GetBinLowEdge(first_bin)
        high = axis.GetBinUpEdge(last_bin)
        first_edge = low if first_edge is None else min(first_edge, low)
        last_edge = high if last_edge is None else max(last_edge, high)
    return first_edge, last_edge


def find_visible_yextrema(era_hists, x_low, x_high):
    y_max = 0.0
    min_positive = None
    for hist in era_hists.values():
        axis = hist.GetXaxis()
        for idx in range(1, hist.GetNbinsX() + 1):
            center = axis.GetBinCenter(idx)
            if center < x_low or center > x_high:
                continue
            value = float(hist.GetBinContent(idx))
            if value <= 0:
                continue
            y_max = max(y_max, value)
            min_positive = value if min_positive is None else min(min_positive, value)
    return min_positive, y_max


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


def plot_quantity_per_era(cfg, job, era_hists, era_sources, out_base, logy=False):
    plotting = cfg.get("plotting", {})
    CMS.SetExtraText(str(plotting.get("cms_extra_text", "Preliminary")))
    CMS.SetEnergy(str(plotting.get("energy_tev", cfg.get("energy_tev", 13.6))))
    CMS.SetLumi(str(plotting.get("lumi_text", "Run3 2026")))

    quantity = str(job["quantity"])
    x_title = str(job.get("x_title", quantity))
    y_title = str(job.get("y_title", "A.U."))
    ymin_floor = float(job.get("ymin_floor", 1e-4))
    y_min = float(job.get("ymin", 0.01 if logy else 0.0))
    ratio_ymin = float(job.get("ratio_ymin", 0.5))
    ratio_ymax = float(job.get("ratio_ymax", 1.5))

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
    else:
        zoom_low, zoom_high = find_populated_xrange(era_hists, margin_bins=int(job.get("zoom_margin_bins", 1)))
        if zoom_low is not None and zoom_high is not None and zoom_high > zoom_low:
            x_low = float(zoom_low)
            x_high = float(zoom_high)

    visible_ymin, visible_ymax = find_visible_yextrema(era_hists, float(x_low), float(x_high))
    if visible_ymax > 0:
        ymax = visible_ymax
    if visible_ymin is not None:
        if logy:
            y_min = max(ymin_floor, visible_ymin / float(job.get("ymin_divisor", 5.0)))
        else:
            y_min = float(job.get("ymin_linear", 0.0))
    if logy and y_min <= 0:
        if visible_ymin is not None and visible_ymin > 0:
            y_min = max(ymin_floor, visible_ymin / float(job.get("ymin_divisor", 5.0)))
        else:
            y_min = ymin_floor
    ymax_scale_default = 5.0 if logy else 1.18
    min_span_default = 10.0 if logy else 1.05
    y_max = max(
        ymax * float(job.get("ymax_scale_log" if logy else "ymax_scale_linear", job.get("ymax_scale", ymax_scale_default))),
        y_min * float(job.get("min_y_span_factor_log" if logy else "min_y_span_factor_linear", job.get("min_y_span_factor", min_span_default))),
    )

    ordered_eras = [era_key for era_key in era_sources if era_key in era_hists]
    reference_era = ordered_eras[0]
    reference_hist = era_hists[reference_era]
    reference_label = source_display_label(era_sources[reference_era])
    ratio_title = str(job.get("ratio_title", f"other / {reference_label}"))

    canvas = CMS.cmsDiCanvas(
        "",
        float(x_low),
        float(x_high),
        y_min,
        y_max,
        ratio_ymin,
        ratio_ymax,
        x_title,
        y_title,
        ratio_title,
        square=CMS.kSquare,
        extraSpace=0.0,
        iPos=0,
    )
    register_canvas_input(canvas, "canvas", canvas)

    legend = CMS.cmsLeg(0.72, 0.68, 0.95, 0.90, textSize=float(job.get("legend_text_size", 0.03)))
    register_canvas_input(canvas, "legend", legend)

    canvas.cd(1)
    for idx, era_key in enumerate(ordered_eras):
        hist = era_hists[era_key]
        color = COLOR_TEMPLATE[idx % len(COLOR_TEMPLATE)]
        band = make_error_band(hist, f"{hist.GetName()}_band", color)
        register_canvas_input(canvas, f"band::{era_key}", band)
        CMS.cmsDraw(band, "E1 X0", lcolor=color, fstyle=0, lwidth=2, lstyle=2, msize=0)
        register_canvas_input(canvas, f"hist::{era_key}", hist)
        CMS.cmsDraw(hist, "HIST SAME", lcolor=color, msize=0, lwidth=3, fstyle=0)
        legend.AddEntry(hist, source_display_label(era_sources[era_key]), "L")

    if logy:
        canvas.cd(1)
        ROOT.gPad.SetLogy()

    canvas.cd(2)
    first_ratio = True
    for idx, era_key in enumerate(ordered_eras):
        if era_key == reference_era:
            continue
        hist = era_hists[era_key]
        ratio = hist.Clone(f"{hist.GetName()}_ratio")
        ratio.SetDirectory(0)
        ratio.Divide(reference_hist)
        color = COLOR_TEMPLATE[idx % len(COLOR_TEMPLATE)]
        ratio_band = make_error_band(ratio, f"{ratio.GetName()}_band", color)
        register_canvas_input(canvas, f"ratio_band::{era_key}", ratio_band)
        CMS.cmsDraw(
            ratio_band,
            "E1 X0" if first_ratio else "E1 X0 SAME",
            lcolor=color,
            fstyle=0,
            lwidth=2,
            lstyle=2,
            msize=0,
        )
        register_canvas_input(canvas, f"ratio::{era_key}", ratio)
        draw_opt = "HIST SAME"
        CMS.cmsDraw(ratio, draw_opt, lcolor=color, msize=0, lwidth=3, fstyle=0)
        first_ratio = False

    ratio_unity = reference_hist.Clone(f"{reference_hist.GetName()}_unity")
    ratio_unity.SetDirectory(0)
    for idx in range(1, ratio_unity.GetNbinsX() + 1):
        ratio_unity.SetBinContent(idx, 1.0)
        ratio_unity.SetBinError(idx, 0.0)
    register_canvas_input(canvas, "ratio_unity", ratio_unity)
    CMS.cmsDraw(ratio_unity, "HIST" if first_ratio else "HIST SAME", lcolor=ROOT.kBlack, lwidth=2, lstyle=2, fstyle=0)

    suffix = "_log" if logy else "_linear"
    CMS.SaveCanvas(canvas, str(out_base.with_name(f"{out_base.name}{suffix}").with_suffix(".png")), close=False)
    CMS.SaveCanvas(canvas, str(out_base.with_name(f"{out_base.name}{suffix}").with_suffix(".pdf")))


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

                rebin_factor = int(job.get("rebin", default_rebin))
                hist = hist.Clone(f"{sanitize(job_name)}_{sanitize(era_key)}_{sanitize(hist.GetName())}")
                hist.SetDirectory(0)
                if rebin_factor > 1:
                    hist = hist.Rebin(rebin_factor, f"{hist.GetName()}_rebin")
                    hist.SetDirectory(0)

                count_hist = hist.Clone(f"{hist.GetName()}_counts")
                count_hist.SetDirectory(0)
                scale_factor = default_scale_to / integral
                hist.Scale(scale_factor)
                apply_poisson_errors_from_counts(hist, count_hist, scale_factor=scale_factor)

                x_range = job.get("x_range")
                if x_range and len(x_range) == 2:
                    hist.GetXaxis().SetRangeUser(float(x_range[0]), float(x_range[1]))

                era_hists[era_key] = hist
                summary[job_name][era_key] = {
                    "status": "ok",
                    "used_runs": used_runs,
                    "integral_before_scale": float(integral),
                    "scale_factor": float(scale_factor),
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
            save_linear = bool(job.get("save_linear", section.get("save_linear", True)))
            save_log = bool(job.get("save_log", section.get("save_log", True)))
            if save_linear:
                plot_quantity_per_era(cfg, job, ordered_hists, era_sources, out_base, logy=False)
            if save_log:
                plot_quantity_per_era(cfg, job, ordered_hists, era_sources, out_base, logy=True)
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
