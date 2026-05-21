from pathlib import Path

import ROOT
import cmsstyle as CMS
import yaml

from dqm_pipeline.core import aggregate_histogram_for_era, emit_log, sanitize


MODULE_NAME = "dqm_monitoring"


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


def source_lumi_label(cfg, source):
    plotting = cfg.get("plotting", {})
    energy = plotting.get("energy_tev", cfg.get("energy_tev", 13.6))
    label = source_display_label(source)
    if source.get("lumi_fb") is not None:
        return f"{label} {float(source['lumi_fb']):.2f} fb^{{-1}} ({energy} TeV)"
    return f"{label} ({energy} TeV)"


def global_lumi_label(cfg):
    plotting = cfg.get("plotting", {})
    if plotting.get("lumi_text"):
        return str(plotting["lumi_text"])
    energy = plotting.get("energy_tev", cfg.get("energy_tev", 13.6))
    campaign = plotting.get("campaign_label", cfg.get("campaign_label", "Run 3"))
    return f"{campaign} ({energy} TeV)"


def register_canvas_input(canvas, key, obj):
    if not hasattr(canvas, "_input_cache"):
        canvas._input_cache = {}
    canvas._input_cache[str(key)] = obj
    return obj


def normalize_hist(hist, mode):
    if mode in (None, "none"):
        return hist

    integral = hist.Integral()
    if integral <= 0:
        return hist

    if mode in ("unit", "area", "integral"):
        hist.Scale(1.0 / integral)
    return hist


def draw_hist_overlay(cfg, hist_name, era_hists, era_sources, out_base, job):
    plotting = cfg.get("plotting", {})
    CMS.SetExtraText(str(plotting.get("cms_extra_text", "Preliminary")))
    CMS.SetEnergy(str(plotting.get("energy_tev", cfg.get("energy_tev", 13.6))))
    CMS.SetLumi(global_lumi_label(cfg))

    x_low = job.get("x_min")
    x_high = job.get("x_max")
    y_max = 0.0
    min_positive = None

    for hist in era_hists.values():
        axis = hist.GetXaxis()
        x_low = axis.GetXmin() if x_low is None else min(float(x_low), axis.GetXmin())
        x_high = axis.GetXmax() if x_high is None else max(float(x_high), axis.GetXmax())
        y_max = max(y_max, hist.GetMaximum())
        for idx in range(1, hist.GetNbinsX() + 1):
            value = hist.GetBinContent(idx)
            if value > 0:
                min_positive = value if min_positive is None else min(min_positive, value)

    logy = bool(job.get("logy", False))
    y_min = float(job.get("y_min", max(1e-4, min_positive / 5.0) if logy and min_positive else 0.0))
    y_max = float(job.get("y_max", y_max * (5.0 if logy else 1.25) if y_max > 0 else 1.0))
    x_title = str(job.get("x_title", hist_name))
    y_title = str(job.get("y_title", "Events"))

    canvas = CMS.cmsCanvas(
        "",
        float(x_low),
        float(x_high),
        y_min,
        y_max,
        x_title,
        y_title,
        square=CMS.kSquare,
        extraSpace=0.0,
        iPos=0,
    )
    register_canvas_input(canvas, "canvas", canvas)
    if logy:
        ROOT.gPad.SetLogy()

    legend = CMS.cmsLeg(0.58, 0.70, 0.90, 0.90, textSize=float(job.get("legend_text_size", 0.03)))
    register_canvas_input(canvas, "legend", legend)

    for idx, era_key in enumerate(era_sources):
        if era_key not in era_hists:
            continue
        hist = era_hists[era_key]
        color = COLOR_TEMPLATE[idx % len(COLOR_TEMPLATE)]
        hist.SetLineColor(color)
        hist.SetMarkerColor(color)
        hist.SetLineWidth(int(job.get("line_width", 3)))
        register_canvas_input(canvas, f"hist::{era_key}", hist)
        CMS.cmsDraw(hist, "HIST SAME", lcolor=color, msize=0, lwidth=int(job.get("line_width", 3)), fstyle=0)
        legend.AddEntry(hist, source_display_label(era_sources[era_key]), "L")

    CMS.SaveCanvas(canvas, str(out_base.with_suffix(".png")), close=False)
    CMS.SaveCanvas(canvas, str(out_base.with_suffix(".pdf")))


def build_comparison_hist_name(resonance, probe_type, quantity, region):
    return f"{resonance}_Probe_{probe_type}Electron_{quantity}_{region}"


def draw_pat_sct_comparison(cfg, era_key, source, quantity, region, pat_hist, sct_hist, out_base, job):
    plotting = cfg.get("plotting", {})
    CMS.SetExtraText(str(plotting.get("cms_extra_text", "Preliminary")))
    CMS.SetEnergy(str(plotting.get("energy_tev", cfg.get("energy_tev", 13.6))))
    CMS.SetLumi(str(job.get("lumi_text", source_lumi_label(cfg, source))))

    pat = pat_hist.Clone(f"{pat_hist.GetName()}_{sanitize(era_key)}_pat_norm")
    sct = sct_hist.Clone(f"{sct_hist.GetName()}_{sanitize(era_key)}_sct_norm")
    pat.SetDirectory(0)
    sct.SetDirectory(0)

    normalize_hist(pat, job.get("normalize", "unit"))
    normalize_hist(sct, job.get("normalize", "unit"))

    x_axis = pat.GetXaxis()
    x_min = float(job.get("x_min", x_axis.GetXmin()))
    x_max = float(job.get("x_max", x_axis.GetXmax()))
    y_max = max(pat.GetMaximum(), sct.GetMaximum()) * float(job.get("ymax_scale", 1.25))
    if y_max <= 0:
        y_max = 1.0

    canvas = CMS.cmsDiCanvas(
        "",
        x_min,
        x_max,
        float(job.get("y_min", 0.0)),
        y_max,
        float(job.get("ratio_ymin", 0.0)),
        float(job.get("ratio_ymax", 2.0)),
        str(job.get("x_title", quantity)),
        str(job.get("y_title", "A.U.")),
        str(job.get("ratio_title", "Scouting / PAT")),
        square=CMS.kSquare,
        extraSpace=0.0,
        iPos=0,
    )
    register_canvas_input(canvas, "canvas", canvas)

    canvas.cd(1)
    legend = CMS.cmsLeg(0.58, 0.68, 0.90, 0.88, textSize=float(job.get("legend_text_size", 0.035)))
    register_canvas_input(canvas, "legend", legend)
    register_canvas_input(canvas, "pat", pat)
    register_canvas_input(canvas, "sct", sct)
    CMS.cmsDraw(pat, "HIST SAME", lcolor=ROOT.kBlue + 1, msize=0, lwidth=3, fstyle=0)
    CMS.cmsDraw(sct, "HIST SAME", lcolor=ROOT.kRed + 1, msize=0, lwidth=3, fstyle=0)
    legend.AddEntry(pat, "PAT electron", "L")
    legend.AddEntry(sct, "Scouting electron", "L")

    canvas.cd(2)
    ratio = sct.Clone(f"{sct.GetName()}_ratio")
    ratio.SetDirectory(0)
    ratio.Divide(pat)
    register_canvas_input(canvas, "ratio", ratio)
    CMS.cmsDraw(ratio, "HIST", lcolor=ROOT.kRed + 1, msize=0, lwidth=3, fstyle=0)

    CMS.SaveCanvas(canvas, str(out_base.with_suffix(".png")), close=False)
    CMS.SaveCanvas(canvas, str(out_base.with_suffix(".pdf")))


def expand_comparison_jobs(section):
    jobs = section.get("comparison_jobs")
    if jobs:
        return jobs

    quantities = section.get("quantities", [])
    regions = section.get("regions", [])
    return [
        {
            "name": f"{quantity}_{region}",
            "quantity": quantity,
            "region": region,
        }
        for quantity in quantities
        for region in regions
    ]


def run_module(cfg, era_sources, out_root, strict=False, progress=None):
    section = cfg.get(MODULE_NAME)
    if not section or not section.get("enabled", True):
        return {}

    resonance = str(cfg["resonance"])
    path_template = section.get("hist_path_template", "DQMData/Run {run}/ScoutingMonitoring/Run summary/{hist}")
    out_dir = Path(out_root) / section.get("output_subdir", MODULE_NAME)
    out_dir.mkdir(parents=True, exist_ok=True)

    emit_log(progress, f"[{MODULE_NAME}] start eras={len(era_sources)} out_dir={out_dir}", style="blue")

    summary = {
        "histograms": {},
        "comparisons": {},
    }

    single_histograms = section.get("histograms", [])
    for raw_job in single_histograms:
        job = dict(raw_job) if isinstance(raw_job, dict) else {"name": str(raw_job), "hist": str(raw_job)}
        hist_name = str(job["hist"])
        job_name = str(job.get("name", hist_name))
        era_hists = {}
        summary["histograms"][job_name] = {}

        for era_key, source in era_sources.items():
            try:
                hist, used_runs = aggregate_histogram_for_era(
                    era=era_key,
                    source=source,
                    hist_path_template=path_template,
                    fmt_args={"hist": hist_name},
                    strict=strict,
                )
                if hist is None:
                    summary["histograms"][job_name][era_key] = {"status": "empty"}
                    continue
                hist = hist.Clone(f"{sanitize(job_name)}_{sanitize(era_key)}")
                hist.SetDirectory(0)
                rebin = int(job.get("rebin", section.get("rebin", 1)))
                if rebin > 1:
                    hist = hist.Rebin(rebin, f"{hist.GetName()}_rebin")
                    hist.SetDirectory(0)
                normalize_hist(hist, job.get("normalize", section.get("normalize", "none")))
                era_hists[era_key] = hist
                summary["histograms"][job_name][era_key] = {"status": "ok", "used_runs": used_runs}
            except Exception as exc:
                summary["histograms"][job_name][era_key] = {"status": "error", "message": str(exc)}
                print(f"[{MODULE_NAME}][WARN] histogram={hist_name} era={era_key}: {exc}")
                if strict:
                    raise

        if era_hists:
            draw_hist_overlay(
                cfg=cfg,
                hist_name=hist_name,
                era_hists=era_hists,
                era_sources=era_sources,
                out_base=out_dir / sanitize(job_name),
                job=job,
            )

    comparison_jobs = expand_comparison_jobs(section)
    for raw_job in comparison_jobs:
        job = dict(raw_job)
        quantity = str(job["quantity"])
        region = str(job["region"])
        job_name = str(job.get("name", f"{quantity}_{region}"))
        pat_name = str(job.get("pat_hist", build_comparison_hist_name(resonance, "pat", quantity, region)))
        sct_name = str(job.get("sct_hist", build_comparison_hist_name(resonance, "sct", quantity, region)))
        summary["comparisons"][job_name] = {}

        for era_key, source in era_sources.items():
            try:
                pat_hist, pat_runs = aggregate_histogram_for_era(
                    era=era_key,
                    source=source,
                    hist_path_template=path_template,
                    fmt_args={"hist": pat_name},
                    strict=strict,
                )
                sct_hist, sct_runs = aggregate_histogram_for_era(
                    era=era_key,
                    source=source,
                    hist_path_template=path_template,
                    fmt_args={"hist": sct_name},
                    strict=strict,
                )
                if pat_hist is None or sct_hist is None:
                    summary["comparisons"][job_name][era_key] = {"status": "empty"}
                    continue

                rebin = int(job.get("rebin", section.get("comparison_rebin", section.get("rebin", 1))))
                if rebin > 1:
                    pat_hist = pat_hist.Rebin(rebin, f"{sanitize(job_name)}_{sanitize(era_key)}_pat_rebin")
                    sct_hist = sct_hist.Rebin(rebin, f"{sanitize(job_name)}_{sanitize(era_key)}_sct_rebin")
                    pat_hist.SetDirectory(0)
                    sct_hist.SetDirectory(0)

                draw_pat_sct_comparison(
                    cfg=cfg,
                    era_key=era_key,
                    source=source,
                    quantity=quantity,
                    region=region,
                    pat_hist=pat_hist,
                    sct_hist=sct_hist,
                    out_base=out_dir / f"{sanitize(job_name)}_{sanitize(era_key)}",
                    job=job,
                )
                summary["comparisons"][job_name][era_key] = {
                    "status": "ok",
                    "pat_hist": pat_name,
                    "sct_hist": sct_name,
                    "pat_used_runs": pat_runs,
                    "sct_used_runs": sct_runs,
                }
            except Exception as exc:
                summary["comparisons"][job_name][era_key] = {"status": "error", "message": str(exc)}
                print(f"[{MODULE_NAME}][WARN] comparison={job_name} era={era_key}: {exc}")
                if strict:
                    raise

    summary_file = out_dir / f"{MODULE_NAME}_summary.yaml"
    with open(summary_file, "w", encoding="utf-8") as handle:
        yaml.safe_dump(summary, handle, sort_keys=False)

    emit_log(progress, f"[{MODULE_NAME}] done. outputs={out_dir}", style="blue")
    return summary
