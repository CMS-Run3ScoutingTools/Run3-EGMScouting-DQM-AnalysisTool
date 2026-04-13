from pathlib import Path

import ROOT
import cmsstyle as CMS
import yaml
import numpy as np

from dqm_pipeline.core import (
    aggregate_histogram_for_era,
    emit_log,
    load_histogram,
    rebin_histogram,
    sanitize,
)


MODULE_NAME = "tnp_efficiency"


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

MARKER_STYLES = [20, 24, 25, 26, 27, 28]


def source_display_label(source):
    return str(source.get("display_label", source.get("era", "era")))


def global_plot_lumi_text(cfg):
    plotting = cfg.get("plotting", {})
    if plotting.get("lumi_text"):
        return str(plotting["lumi_text"])
    return "Run3"


def register_canvas_input(canvas, key, obj):
    if not hasattr(canvas, "_input_cache"):
        canvas._input_cache = {}
    canvas._input_cache[str(key)] = obj
    return obj


def set_eff_style(eff, color, marker, marker_size=1.0, line_style=1):
    eff.SetLineColor(color)
    eff.SetMarkerColor(color)
    eff.SetMarkerStyle(marker)
    eff.SetMarkerSize(marker_size)
    eff.SetLineWidth(2)
    eff.SetLineStyle(line_style)


def make_legend(x1, y1, x2, y2, text_size=0.02, ncols=1):
    legend = CMS.cmsLeg(x1, y1, x2, y2, textSize=text_size)
    try:
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)
        legend.SetNColumns(ncols)
    except Exception:
        pass
    return legend


def default_bins_for_original_behavior(axis, numerator_name):
    axis = str(axis).lower()
    if axis == "eta":
        return [-2.5 + i * 0.5 for i in range(11)]
    if "Single" in numerator_name:
        return [0, 5, 10, 15, 20, 30, 40, 50, 80, 120, 200]
    return [0, 10, 20, 30, 40, 50, 120]


def build_tnp_hist_names(resonance, job, tag):
    tagging_type = job.get("tagging_type", "pat")
    probe_type = job.get("probe_type", "pat")
    pt_order = job.get("pt_order", "leading")
    axis = job.get("axis", "pt").lower()
    region = job.get("region", "Barrel")

    if axis not in ["pt", "eta"]:
        raise ValueError(f"Unsupported axis '{axis}'. Use 'pt' or 'eta'.")

    if tagging_type == "pat":
        if axis == "pt":
            numerator = f"{resonance}_Tag_{tagging_type}_Probe_{probe_type}Electron_{pt_order}_Pt_{region}_pass{tag}"
            denominator = f"{resonance}_Tag_{tagging_type}_Probe_{probe_type}Electron_{pt_order}_Pt_{region}_passBaseDST"
        else:
            numerator = f"{resonance}_Tag_{tagging_type}_Probe_{probe_type}Electron_{pt_order}_Eta_pass{tag}"
            denominator = f"{resonance}_Tag_{tagging_type}_Probe_{probe_type}Electron_{pt_order}_Eta_passBaseDST"
    else:
        if axis == "pt":
            numerator = f"{resonance}_{pt_order}_Pt_{region}_pass{tag}"
            denominator = f"{resonance}_{pt_order}_Pt_{region}_passBaseDST"
        else:
            numerator = f"{resonance}_{pt_order}_Eta_pass{tag}"
            denominator = f"{resonance}_{pt_order}_Eta_passBaseDST"

    target_dir = job.get("target_dir")
    if not target_dir:
        target_dir = "Tag_PatElectron" if tagging_type == "pat" else "Tag_ScoutingElectron"

    return target_dir, numerator, denominator


def compute_efficiency_points(num_hist, den_hist):
    if not ROOT.TEfficiency.CheckConsistency(num_hist, den_hist):
        raise RuntimeError("TEfficiency consistency check failed (num/den bins mismatch or num>den).")

    teff = ROOT.TEfficiency(num_hist, den_hist)
    n_bins = num_hist.GetNbinsX()

    x, xerr, y, yerr_low, yerr_up = [], [], [], [], []
    for i_bin in range(1, n_bins + 1):
        x.append(num_hist.GetXaxis().GetBinCenter(i_bin))
        xerr.append(0.5 * num_hist.GetXaxis().GetBinWidth(i_bin))
        y.append(teff.GetEfficiency(i_bin))
        yerr_low.append(teff.GetEfficiencyErrorLow(i_bin))
        yerr_up.append(teff.GetEfficiencyErrorUp(i_bin))

    return {
        "x": np.asarray(x),
        "xerr": np.asarray(xerr),
        "y": np.asarray(y),
        "yerr_low": np.asarray(yerr_low),
        "yerr_up": np.asarray(yerr_up),
    }


def build_central_eff_hist(num_hist, den_hist, name):
    ratio = num_hist.Clone(name)
    ratio.SetDirectory(0)
    ratio.Divide(den_hist)
    return ratio


def plot_tnp_efficiency(cfg, job, tag, per_era_data, era_sources, out_base, mc_data=None):
    plotting = cfg.get("plotting", {})
    energy = str(plotting.get("energy_tev", cfg.get("energy_tev", 13.6)))
    CMS.SetExtraText(str(plotting.get("cms_extra_text", "Preliminary")))
    CMS.SetEnergy(energy)
    CMS.SetLumi(str(plotting.get("lumi_text", global_plot_lumi_text(cfg))))

    axis = job.get("axis", "pt").lower()
    pt_order = job.get("pt_order", "leading")
    x_label = job.get("x_label", f"{pt_order} p_{{T}} [GeV]" if axis == "pt" else "#eta")
    ratio_label = str(job.get("ratio_label", "w. ref"))
    eff_ymin = float(job.get("ymin", 0.0))
    eff_ymax = float(job.get("ymax", 1.08))
    ratio_ymin = float(job.get("ratio_ymin", 0.9))
    ratio_ymax = float(job.get("ratio_ymax", 1.1))

    ordered_eras = [era for era in era_sources if era in per_era_data]
    if not ordered_eras:
        raise RuntimeError("No era data to plot.")

    sample_eff = per_era_data[ordered_eras[0]]["eff"]
    x_axis = sample_eff.GetPassedHistogram().GetXaxis()
    x_min = x_axis.GetXmin()
    x_max = x_axis.GetXmax()

    canvas = CMS.cmsDiCanvas(
        "",
        x_min,
        x_max,
        eff_ymin,
        eff_ymax,
        ratio_ymin,
        ratio_ymax,
        x_label,
        "efficiency",
        ratio_label,
        square=CMS.kSquare,
        extraSpace=0.0,
        iPos=0,
    )
    register_canvas_input(canvas, "canvas", canvas)

    canvas.cd(1)
    legend = make_legend(0.40, 0.20, 0.80, 0.40, text_size=float(job.get("legend_text_size", 0.02)), ncols=1)
    register_canvas_input(canvas, "legend", legend)

    reference_era = ordered_eras[0]
    reference_ratio = per_era_data[reference_era]["ratio_hist"]

    for idx, era in enumerate(ordered_eras):
        entry = per_era_data[era]
        color = COLOR_TEMPLATE[idx % len(COLOR_TEMPLATE)]
        marker = MARKER_STYLES[idx % len(MARKER_STYLES)]
        eff = entry["eff"]
        set_eff_style(eff, color, marker, marker_size=1.0, line_style=1)
        register_canvas_input(canvas, f"eff::{era}", eff)
        CMS.cmsDraw(eff, "P E", lcolor=color, mcolor=color, msize=1.0, lwidth=2, fstyle=0)
        legend.AddEntry(eff, source_display_label(era_sources[era]), "PL")

    if mc_data is not None:
        mc_eff = mc_data["eff"]
        set_eff_style(mc_eff, ROOT.kBlack, 34, marker_size=1.0, line_style=1)
        register_canvas_input(canvas, "eff::MC", mc_eff)
        CMS.cmsDraw(mc_eff, "P E", lcolor=ROOT.kBlack, mcolor=ROOT.kBlack, msize=1.0, lwidth=2, fstyle=0)
        legend.AddEntry(mc_eff, "MC", "PL")

    canvas.cd(2)
    unity = reference_ratio.Clone(f"{reference_ratio.GetName()}_unity")
    unity.SetDirectory(0)
    for i_bin in range(1, unity.GetNbinsX() + 1):
        unity.SetBinContent(i_bin, 1.0)
        unity.SetBinError(i_bin, 0.0)
    register_canvas_input(canvas, "unity", unity)
    CMS.cmsDraw(unity, "HIST", lcolor=ROOT.kBlack, lwidth=2, lstyle=2, fstyle=0)

    first_ratio = True
    for idx, era in enumerate(ordered_eras):
        if era == reference_era:
            continue
        entry = per_era_data[era]
        color = COLOR_TEMPLATE[idx % len(COLOR_TEMPLATE)]
        ratio_hist = entry["ratio_hist"].Clone(f"{entry['ratio_hist'].GetName()}_vs_ref")
        ratio_hist.SetDirectory(0)
        ratio_hist.Divide(reference_ratio)
        register_canvas_input(canvas, f"ratio::{era}", ratio_hist)
        CMS.cmsDraw(
            ratio_hist,
            "HIST SAME" if not first_ratio else "HIST SAME",
            lcolor=color,
            lwidth=3,
            fstyle=0,
        )
        first_ratio = False

    if mc_data is not None:
        mc_ratio = mc_data["ratio_hist"].Clone(f"{mc_data['ratio_hist'].GetName()}_vs_ref")
        mc_ratio.SetDirectory(0)
        mc_ratio.Divide(reference_ratio)
        register_canvas_input(canvas, "ratio::MC", mc_ratio)
        CMS.cmsDraw(mc_ratio, "HIST SAME", lcolor=ROOT.kBlack, lwidth=3, lstyle=1, fstyle=0)

    CMS.SaveCanvas(canvas, str(out_base.with_suffix(".png")), close=False)
    CMS.SaveCanvas(canvas, str(out_base.with_suffix(".pdf")))


def run_module(cfg, era_sources, out_root, strict=False, progress=None):
    section = cfg.get(MODULE_NAME)
    if not section or not section.get("enabled", True):
        return {}

    resonance = cfg["resonance"]
    path_template = section.get(
        "hist_path_template",
        "DQMData/Run {run}/HLT/Run summary/ScoutingOffline/EGamma/TnP/{target_dir}/{hist}",
    )

    bins_cfg = section.get("bins", {})
    default_rebin = int(section.get("rebin", 1))
    jobs = section.get("jobs", [])

    out_dir = Path(out_root) / section.get("output_subdir", MODULE_NAME)
    out_dir.mkdir(parents=True, exist_ok=True)
    emit_log(
        progress,
        f"[{MODULE_NAME}] start jobs={len(jobs)} eras={len(era_sources)} out_dir={out_dir}",
        style="bold blue",
    )

    mc_file = section.get("mc_file")
    mc_run_number = section.get("mc_run_number")
    plot_lumi_text = global_plot_lumi_text(cfg)

    all_results = {}
    total_successful_eras = 0
    total_failed_eras = 0
    job_task = None
    if progress is not None:
        job_task = progress.add_task(f"[blue]{MODULE_NAME}: jobs", total=len(jobs))

    for job in jobs:
        job_name = job["name"]
        emit_log(progress, f"[{MODULE_NAME}] job start: {job_name}", style="blue")
        tags = job.get("tags", [])
        axis = job.get("axis", "pt").lower()
        default_bins = bins_cfg.get(axis)
        rebin_factor = int(job.get("rebin", default_rebin))

        all_results[job_name] = {}
        job_successful_eras = 0
        job_failed_eras = 0
        tag_task = None
        if progress is not None:
            tag_task = progress.add_task(f"[blue]{MODULE_NAME}: {job_name} tags", total=len(tags))

        for tag in tags:
            emit_log(progress, f"[{MODULE_NAME}] tag start: {job_name}/{tag}", style="blue")
            target_dir, num_name, den_name = build_tnp_hist_names(resonance, job, tag)
            bins = job.get("bins", default_bins)
            if bins is None:
                bins = default_bins_for_original_behavior(axis=axis, numerator_name=num_name)
            all_results[job_name][tag] = {}
            tag_successful_eras = 0
            tag_failed_eras = 0
            per_era_data = {}
            era_task = None
            era_items = list(era_sources.items())
            if progress is not None:
                era_task = progress.add_task(
                    f"[blue]{MODULE_NAME}: {job_name}/{tag} eras",
                    total=len(era_items),
                )

            for era, source in era_items:
                try:
                    num_hist, used_runs_num = aggregate_histogram_for_era(
                        era=era,
                        source=source,
                        hist_path_template=path_template,
                        fmt_args={"target_dir": target_dir, "hist": num_name},
                        strict=strict,
                    )
                    den_hist, used_runs_den = aggregate_histogram_for_era(
                        era=era,
                        source=source,
                        hist_path_template=path_template,
                        fmt_args={"target_dir": target_dir, "hist": den_name},
                        strict=strict,
                    )

                    if num_hist is None or den_hist is None:
                        raise RuntimeError("No numerator/denominator histogram found in selected runs.")

                    num_hist = rebin_histogram(
                        num_hist,
                        bins=bins,
                        rebin_factor=rebin_factor,
                        name_hint=f"num_{sanitize(job_name)}_{sanitize(tag)}_{sanitize(era)}",
                    )
                    den_hist = rebin_histogram(
                        den_hist,
                        bins=bins,
                        rebin_factor=rebin_factor,
                        name_hint=f"den_{sanitize(job_name)}_{sanitize(tag)}_{sanitize(era)}",
                    )

                    eff = ROOT.TEfficiency(num_hist, den_hist)
                    ratio_hist = build_central_eff_hist(
                        num_hist,
                        den_hist,
                        f"ratio_{sanitize(job_name)}_{sanitize(tag)}_{sanitize(era)}",
                    )
                    points = compute_efficiency_points(num_hist, den_hist)
                    per_era_data[era] = {
                        "eff": eff,
                        "num_hist": num_hist,
                        "den_hist": den_hist,
                        "ratio_hist": ratio_hist,
                        "points": points,
                    }

                    all_results[job_name][tag][era] = {
                        "numerator": num_name,
                        "denominator": den_name,
                        "used_runs_num": used_runs_num,
                        "used_runs_den": used_runs_den,
                        "mean_efficiency": float(np.mean(points["y"])),
                    }
                    tag_successful_eras += 1
                except Exception as exc:
                    tag_failed_eras += 1
                    print(f"[{MODULE_NAME}][WARN] job={job_name} tag={tag} era={era}: {exc}")
                    if strict:
                        raise
                finally:
                    if progress is not None and era_task is not None:
                        progress.update(era_task, advance=1)

            mc_points = None
            mc_data = None
            include_mc_for_this_tag = bool(job.get("include_mc", False)) and (
                len(tags) == 1 or bool(job.get("include_mc_when_multiple_tags", False))
            )
            if include_mc_for_this_tag and mc_file:
                try:
                    if "{run}" in path_template:
                        if mc_run_number is None:
                            raise RuntimeError(
                                "tnp_efficiency.mc_run_number is required when include_mc=true and hist_path_template uses {run}."
                            )
                        num_path_mc = path_template.format(run=int(mc_run_number), target_dir=target_dir, hist=num_name)
                        den_path_mc = path_template.format(run=int(mc_run_number), target_dir=target_dir, hist=den_name)
                    else:
                        num_path_mc = path_template.format(target_dir=target_dir, hist=num_name)
                        den_path_mc = path_template.format(target_dir=target_dir, hist=den_name)

                    num_hist_mc = load_histogram(mc_file, num_path_mc)
                    den_hist_mc = load_histogram(mc_file, den_path_mc)

                    num_hist_mc = rebin_histogram(
                        num_hist_mc,
                        bins=bins,
                        rebin_factor=rebin_factor,
                        name_hint=f"num_mc_{sanitize(job_name)}_{sanitize(tag)}",
                    )
                    den_hist_mc = rebin_histogram(
                        den_hist_mc,
                        bins=bins,
                        rebin_factor=rebin_factor,
                        name_hint=f"den_mc_{sanitize(job_name)}_{sanitize(tag)}",
                    )

                    mc_eff = ROOT.TEfficiency(num_hist_mc, den_hist_mc)
                    mc_ratio_hist = build_central_eff_hist(
                        num_hist_mc,
                        den_hist_mc,
                        f"ratio_mc_{sanitize(job_name)}_{sanitize(tag)}",
                    )
                    mc_points = compute_efficiency_points(num_hist_mc, den_hist_mc)
                    mc_data = {
                        "eff": mc_eff,
                        "num_hist": num_hist_mc,
                        "den_hist": den_hist_mc,
                        "ratio_hist": mc_ratio_hist,
                        "points": mc_points,
                    }
                    all_results[job_name][tag]["MC"] = {
                        "numerator": num_name,
                        "denominator": den_name,
                        "mean_efficiency": float(np.mean(mc_points["y"])),
                    }
                except Exception as exc:
                    print(f"[{MODULE_NAME}][WARN] job={job_name} tag={tag} mc: {exc}")
                    if strict:
                        raise

            if per_era_data or mc_data is not None:
                out_base = out_dir / f"{sanitize(job_name)}_{sanitize(tag)}"
                plot_tnp_efficiency(
                    cfg=cfg,
                    job=job,
                    tag=tag,
                    per_era_data=per_era_data,
                    era_sources=era_sources,
                    out_base=out_base,
                    mc_data=mc_data,
                )
            all_results[job_name][tag]["_status"] = {
                "status": "ok" if tag_successful_eras > 0 else "empty",
                "n_successful_eras": tag_successful_eras,
                "n_failed_eras": tag_failed_eras,
            }
            job_successful_eras += tag_successful_eras
            job_failed_eras += tag_failed_eras
            total_successful_eras += tag_successful_eras
            total_failed_eras += tag_failed_eras
            if tag_successful_eras > 0:
                emit_log(
                    progress,
                    f"[{MODULE_NAME}] tag done: {job_name}/{tag} usable_eras={tag_successful_eras} skipped_eras={tag_failed_eras}",
                    style="blue",
                )
            else:
                emit_log(
                    progress,
                    f"[{MODULE_NAME}][WARN] tag done: {job_name}/{tag} usable_eras=0 skipped_eras={tag_failed_eras}",
                    style="yellow",
                )
            if progress is not None and tag_task is not None:
                progress.update(tag_task, advance=1)
        all_results[job_name]["_status"] = {
            "status": "ok" if job_successful_eras > 0 else "empty",
            "n_successful_eras": job_successful_eras,
            "n_failed_eras": job_failed_eras,
        }
        if job_successful_eras > 0:
            emit_log(
                progress,
                f"[{MODULE_NAME}] job done: {job_name} usable_eras={job_successful_eras} skipped_eras={job_failed_eras}",
                style="blue",
            )
        else:
            emit_log(
                progress,
                f"[{MODULE_NAME}][WARN] job done: {job_name} usable_eras=0 skipped_eras={job_failed_eras}",
                style="yellow",
            )
        if progress is not None and job_task is not None:
            progress.update(job_task, advance=1)

    summary_file = out_dir / "tnp_efficiency_summary.yaml"
    with open(summary_file, "w", encoding="utf-8") as f:
        yaml.safe_dump(all_results, f, sort_keys=False)

    final_style = "bold blue" if total_successful_eras > 0 else "yellow"
    final_tag = "done" if total_successful_eras > 0 else "[WARN] done"
    emit_log(
        progress,
        f"[{MODULE_NAME}] {final_tag}. usable_eras={total_successful_eras} skipped_eras={total_failed_eras} outputs={out_dir}",
        style=final_style,
    )
    return all_results
