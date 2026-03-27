from pathlib import Path
import re
from array import array

import ROOT
import cmsstyle as CMS
import yaml

from dqm_pipeline.core import aggregate_histogram_for_era, emit_log, load_histogram, rebin_histogram, sanitize


MODULE_NAME = "filter_monitoring"


DEFAULT_TRIGGER_SETTING = {
    "DST_PFScouting_DoubleEG_v": [
        "hltEGL1DoubleIsoEG11Filter",
        "hltDoubleEG11CaloIdLEt11Filter",
        "hltDoubleEG11CaloIdLClusterShapeFilter",
        "hltDoubleEG11CaloIdLHEFilter",
    ],
    "DST_PFScouting_SinglePhotonEB_v": [
        "hltEGL1SingleEGOrFilter",
        "hltEG30EBTightIDTightIsoEtFilter",
        "hltEG30EBTightIDTightIsoClusterShapeFilter",
        "hltEG30EBTightIDTightIsoHEFilter",
        "hltEG30EBTightIDTightIsoR9Filter",
        "hltEG30EBTightIDTightIsotEcalIsoFilter",
        "hltEG30EBTightIDTightIsoHcalIsoFilter",
        "hltEG30EBTightIDTightIsoTrackIsoFilter",
    ],
}


def safe_tag(value):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value))


def compact_label(value):
    return str(value)


def source_display_label(source):
    return str(source.get("display_label", source.get("era", "era")))


def source_lumi_label(cfg, source):
    plotting = cfg.get("plotting", {})
    energy = plotting.get("energy_tev", cfg.get("energy_tev", 13.6))
    label = source_display_label(source)
    if source.get("lumi_fb") is not None:
        return f"{label} {float(source['lumi_fb']):.2f} fb^{{-1}}"
    return f"{label} ({energy} TeV)"


def default_bins(axis):
    if str(axis).lower() == "eta":
        return [-2.5 + i * 0.5 for i in range(11)]
    return [0, 5, 10, 15, 20, 30, 40, 50, 80, 120, 200]


def build_hist_names(resonance, job, tag):
    tagging_type = job.get("tagging_type", "pat")
    probe_type = job.get("probe_type", "pat")
    pt_order = job.get("pt_order", "leading")
    axis = str(job.get("axis", "pt")).lower()
    region = job.get("region", "Barrel")

    if tagging_type == "pat":
        if axis == "pt":
            numerator = f"{resonance}_Tag_{tagging_type}_Probe_{probe_type}Electron_{pt_order}_Pt_{region}_pass{tag}"
            denominator = f"{resonance}_Tag_{tagging_type}_Probe_{probe_type}Electron_{pt_order}_Pt_{region}_passBaseDST"
            hist_pass_prefix = f"{resonance}_Tag_{tagging_type}_Probe_{probe_type}Electron_{pt_order}_Pt_{region}_pass"
        else:
            numerator = f"{resonance}_Tag_{tagging_type}_Probe_{probe_type}Electron_{pt_order}_Eta_pass{tag}"
            denominator = f"{resonance}_Tag_{tagging_type}_Probe_{probe_type}Electron_{pt_order}_Eta_passBaseDST"
            hist_pass_prefix = f"{resonance}_Tag_{tagging_type}_Probe_{probe_type}Electron_{pt_order}_Eta_pass"
    else:
        if axis == "pt":
            numerator = f"{resonance}_{pt_order}_Pt_{region}_pass{tag}"
            denominator = f"{resonance}_{pt_order}_Pt_{region}_passBaseDST"
            hist_pass_prefix = f"{resonance}_{pt_order}_Pt_{region}_pass"
        else:
            numerator = f"{resonance}_{pt_order}_Eta_pass{tag}"
            denominator = f"{resonance}_{pt_order}_Eta_passBaseDST"
            hist_pass_prefix = f"{resonance}_{pt_order}_Eta_pass"

    target_dir = job.get("target_dir")
    if not target_dir:
        target_dir = "Tag_PatElectron" if tagging_type == "pat" else "Tag_ScoutingElectron"
    return target_dir, numerator, denominator, hist_pass_prefix


def discover_filter_tags(root_file, dirpath, hist_pass_prefix, base_tag, suffix="Filter"):
    directory = root_file.Get(dirpath)
    if not directory:
        return []

    strict_prefix = hist_pass_prefix + base_tag + "_"
    out = []
    seen = set()
    for key in directory.GetListOfKeys():
        name = key.GetName()
        if not name.startswith(strict_prefix):
            continue
        if suffix and not name.endswith(suffix):
            continue
        tag = name[len(hist_pass_prefix):]
        if tag not in seen:
            out.append(tag)
            seen.add(tag)
    return out


def set_eff_style(eff, color, marker_style, marker_size=1.0, line_width=2, line_style=1):
    eff.SetLineColor(color)
    eff.SetMarkerColor(color)
    eff.SetMarkerStyle(marker_style)
    eff.SetMarkerSize(marker_size)
    eff.SetLineWidth(line_width)
    eff.SetLineStyle(line_style)


def make_compact_legend(x1, y1, x2, y2, text_size=0.016, ncols=2):
    legend = CMS.cmsLeg(x1, y1, x2, y2, textSize=text_size)
    try:
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)
        legend.SetNColumns(ncols)
        legend.SetColumnSeparation(0.15)
        legend.SetEntrySeparation(0.12)
        legend.SetMargin(0.18)
    except Exception:
        pass
    return legend


def aggregate_histogram_for_run(source, run, hist_path_template, fmt_args, strict=False):
    merged = None
    hist_path = hist_path_template.format(run=run, **fmt_args)

    for file_path in source["run_files"].get(run, []):
        try:
            hist = load_histogram(
                file_path,
                hist_path,
                xrootd_redirectors=source.get("xrootd_redirectors"),
                suppress_root_errors=bool(source.get("suppress_root_xrd_errors", True)),
            )
        except Exception:
            if strict:
                raise
            continue

        if merged is None:
            merged = hist.Clone(f"run_{sanitize(str(run))}_{sanitize(hist.GetName())}")
            merged.SetDirectory(0)
        else:
            merged.Add(hist)
    return merged


def mean_efficiency(num_hist, den_hist):
    total = 0.0
    n_valid = 0
    for i_bin in range(1, num_hist.GetNbinsX() + 1):
        den = den_hist.GetBinContent(i_bin)
        if den <= 0:
            continue
        num = num_hist.GetBinContent(i_bin)
        total += num / den
        n_valid += 1
    return total / n_valid if n_valid > 0 else None


def draw_filter_canvas(
    cfg,
    source,
    era_key,
    job,
    base_tag,
    base_eff,
    abs_effs,
    step_effs,
    out_base,
    ratio_ymin,
    ratio_ymax,
    legend_cols,
):
    plotting = cfg.get("plotting", {})
    energy = str(plotting.get("energy_tev", cfg.get("energy_tev", 13.6)))
    CMS.SetExtraText(str(plotting.get("cms_extra_text", "Preliminary")))
    CMS.SetEnergy(energy)
    CMS.SetLumi(source_lumi_label(cfg, source))

    axis = str(job.get("axis", "pt")).lower()
    pt_order = job.get("pt_order", "leading")
    x_title = f"{pt_order} p_{{T}} [GeV]" if axis == "pt" else f"{pt_order} #eta"

    # Leave room above the frame so CMS text sits outside the plot area.
    canvas = CMS.cmsDiCanvas(
        "",
        base_eff.GetPassedHistogram().GetXaxis().GetXmin(),
        base_eff.GetPassedHistogram().GetXaxis().GetXmax(),
        0.0,
        1.08,
        ratio_ymin,
        ratio_ymax,
        x_title,
        "efficiency",
        "w. prev",
        square=CMS.kSquare,
        extraSpace=0.08,
        iPos=11,
    )

    palette_hex = [
        "#3f90da", "#ffa90e", "#bd1f01", "#94a4a2", "#832db6",
        "#a96b59", "#e76300", "#b9ac70", "#717581", "#92dadd",
        "#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e",
    ]
    palette = [ROOT.TColor.GetColor(x) for x in palette_hex]
    marker_styles = [20, 21, 22, 23, 33, 34, 29, 30]

    canvas.cd(1)
    legend = make_compact_legend(0.43, 0.58, 0.93, 0.90, text_size=0.018, ncols=legend_cols)

    set_eff_style(base_eff, ROOT.kBlack, 20, marker_size=1.05, line_style=1)
    CMS.cmsDraw(base_eff, "P E", lcolor=ROOT.kBlack, mcolor=ROOT.kBlack, msize=1.0, lwidth=2, fstyle=0)
    legend.AddEntry(base_eff, compact_label(base_tag), "PL")

    for idx, (label, eff_abs) in enumerate(abs_effs.items()):
        color = palette[idx % len(palette)]
        marker = marker_styles[idx % len(marker_styles)]
        set_eff_style(eff_abs, color, marker, marker_size=0.95, line_style=2)
        CMS.cmsDraw(eff_abs, "P E", lcolor=color, mcolor=color, msize=0.95, lwidth=2, fstyle=0)
        legend.AddEntry(eff_abs, compact_label(label), "PL")

    canvas.cd(2)
    for idx, (label, eff_step) in enumerate(step_effs.items()):
        color = palette[idx % len(palette)]
        marker = marker_styles[idx % len(marker_styles)]
        set_eff_style(eff_step, color, marker, marker_size=0.95, line_style=1)
        CMS.cmsDraw(eff_step, "P E", lcolor=color, mcolor=color, msize=0.95, lwidth=2, fstyle=0)

    CMS.SaveCanvas(canvas, str(out_base.with_suffix(".png")), close=False)
    CMS.SaveCanvas(canvas, str(out_base.with_suffix(".pdf")))


def draw_filter_era_comparison_canvas(
    cfg,
    job,
    base_tag,
    filter_name,
    abs_by_era,
    step_by_era,
    out_base,
    legend_cols,
):
    plotting = cfg.get("plotting", {})
    energy = str(plotting.get("energy_tev", cfg.get("energy_tev", 13.6)))
    CMS.SetExtraText(str(plotting.get("cms_extra_text", "Preliminary")))
    CMS.SetEnergy(energy)
    CMS.SetLumi(str(plotting.get("lumi_text", f"Run 3 ({energy} TeV)")))

    sample_eff = next(iter(abs_by_era.values()))
    axis = str(job.get("axis", "pt")).lower()
    pt_order = job.get("pt_order", "leading")
    x_title = f"{pt_order} p_{{T}} [GeV]" if axis == "pt" else f"{pt_order} #eta"

    canvas = CMS.cmsDiCanvas(
        "",
        sample_eff.GetPassedHistogram().GetXaxis().GetXmin(),
        sample_eff.GetPassedHistogram().GetXaxis().GetXmax(),
        0.0,
        1.08,
        0.0,
        1.08,
        x_title,
        "efficiency",
        "step eff.",
        square=CMS.kSquare,
        extraSpace=0.08,
        iPos=11,
    )

    palette_hex = [
        "#3f90da", "#ffa90e", "#bd1f01", "#94a4a2", "#832db6",
        "#a96b59", "#e76300", "#b9ac70", "#717581", "#92dadd",
    ]
    palette = [ROOT.TColor.GetColor(x) for x in palette_hex]
    marker_styles = [20, 21, 22, 23, 33, 34]

    canvas.cd(1)
    legend = make_compact_legend(0.42, 0.62, 0.93, 0.90, text_size=0.020, ncols=legend_cols)
    for idx, (era_label, eff) in enumerate(abs_by_era.items()):
        color = palette[idx % len(palette)]
        marker = marker_styles[idx % len(marker_styles)]
        set_eff_style(eff, color, marker, marker_size=1.0, line_style=1)
        CMS.cmsDraw(eff, "P E", lcolor=color, mcolor=color, msize=1.0, lwidth=2, fstyle=0)
        legend.AddEntry(eff, compact_label(era_label), "PL")

    canvas.cd(2)
    for idx, (era_label, eff) in enumerate(step_by_era.items()):
        color = palette[idx % len(palette)]
        marker = marker_styles[idx % len(marker_styles)]
        set_eff_style(eff, color, marker, marker_size=1.0, line_style=1)
        CMS.cmsDraw(eff, "P E", lcolor=color, mcolor=color, msize=1.0, lwidth=2, fstyle=0)

    CMS.SaveCanvas(canvas, str(out_base.with_suffix(".png")), close=False)
    CMS.SaveCanvas(canvas, str(out_base.with_suffix(".pdf")))


def draw_run_trend_canvas(
    cfg,
    source,
    base_tag,
    title_suffix,
    trend_data,
    out_base,
    legend_cols,
):
    plotting = cfg.get("plotting", {})
    energy = str(plotting.get("energy_tev", cfg.get("energy_tev", 13.6)))
    CMS.SetExtraText(str(plotting.get("cms_extra_text", "Preliminary")))
    CMS.SetEnergy(energy)
    CMS.SetLumi(source_lumi_label(cfg, source))

    runs = sorted({run for points in trend_data.values() for run, _ in points})
    if not runs:
        return

    canvas = CMS.cmsCanvas(
        "",
        float(min(runs)) - 1.0,
        float(max(runs)) + 1.0,
        0.0,
        1.08,
        "run",
        "mean efficiency",
        square=CMS.kSquare,
        extraSpace=0.08,
        iPos=11,
    )
    legend = make_compact_legend(0.40, 0.60, 0.93, 0.90, text_size=0.018, ncols=legend_cols)

    palette_hex = [
        "#3f90da", "#ffa90e", "#bd1f01", "#94a4a2", "#832db6",
        "#a96b59", "#e76300", "#b9ac70", "#717581", "#92dadd",
    ]
    palette = [ROOT.TColor.GetColor(x) for x in palette_hex]
    marker_styles = [20, 21, 22, 23, 33, 34]

    for idx, (filter_name, points) in enumerate(trend_data.items()):
        if not points:
            continue
        xs = array("d", [float(run) for run, _ in points])
        ys = array("d", [float(val) for _, val in points])
        graph = ROOT.TGraph(len(xs), xs, ys)
        color = palette[idx % len(palette)]
        marker = marker_styles[idx % len(marker_styles)]
        graph.SetLineColor(color)
        graph.SetMarkerColor(color)
        graph.SetMarkerStyle(marker)
        graph.SetMarkerSize(0.9)
        graph.SetLineWidth(2)
        draw_opt = "ALP" if idx == 0 else "LP SAME"
        graph.Draw(draw_opt)
        legend.AddEntry(graph, compact_label(filter_name), "PL")

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
    output_subdir = section.get("output_subdir", MODULE_NAME)
    out_dir = Path(out_root) / output_subdir
    out_dir.mkdir(parents=True, exist_ok=True)

    trigger_setting = dict(DEFAULT_TRIGGER_SETTING)
    trigger_setting.update(section.get("trigger_setting", {}))

    jobs = section.get("jobs", [])
    ratio_ymin = float(section.get("ratio_ymin", 0.9))
    ratio_ymax = float(section.get("ratio_ymax", 1.1))
    legend_cols = int(section.get("legend_cols", 2))
    filter_suffix = str(section.get("filter_suffix", "Filter"))
    save_era_comparison = bool(section.get("save_era_comparison", True))
    save_run_trends = bool(section.get("save_run_trends", False))

    emit_log(progress, f"[{MODULE_NAME}] start jobs={len(jobs)} eras={len(era_sources)} out_dir={out_dir}", style="blue")

    summary = {}
    job_task = None
    if progress is not None:
        job_task = progress.add_task(f"{MODULE_NAME}: jobs", total=len(jobs))

    for job in jobs:
        job_name = str(job["name"])
        summary[job_name] = {}
        era_comparison_store = {}
        run_trend_store = {}
        emit_log(progress, f"[{MODULE_NAME}] job start: {job_name}", style="blue")

        base_tags = list(job.get("base_tags", []))
        if not base_tags:
            raise RuntimeError(f"{MODULE_NAME} job '{job_name}' requires non-empty 'base_tags'.")

        axis = str(job.get("axis", "pt")).lower()
        bins = job.get("bins", section.get("bins", {}).get(axis), default_bins(axis))
        rebin_factor = int(job.get("rebin", section.get("rebin", 1)))

        era_task = None
        era_items = list(era_sources.items())
        if progress is not None:
            era_task = progress.add_task(f"{MODULE_NAME}: {job_name} eras", total=len(era_items))

        for era_key, source in era_items:
            era_result = {"status": "empty", "base_tags": {}}
            try:
                for base_tag in base_tags:
                    target_dir, numerator_name, denominator_name, hist_pass_prefix = build_hist_names(resonance, job, base_tag)

                    base_num, base_num_runs = aggregate_histogram_for_era(
                        era=era_key,
                        source=source,
                        hist_path_template=path_template,
                        fmt_args={"target_dir": target_dir, "hist": numerator_name},
                        strict=strict,
                    )
                    base_den, base_den_runs = aggregate_histogram_for_era(
                        era=era_key,
                        source=source,
                        hist_path_template=path_template,
                        fmt_args={"target_dir": target_dir, "hist": denominator_name},
                        strict=strict,
                    )

                    if base_num is None or base_den is None:
                        raise RuntimeError(f"No numerator/denominator histogram found for base tag '{base_tag}'.")

                    base_num = rebin_histogram(
                        base_num,
                        bins=bins,
                        rebin_factor=rebin_factor,
                        name_hint=f"{sanitize(job_name)}_{sanitize(era_key)}_{sanitize(base_tag)}_base_num",
                    )
                    base_den = rebin_histogram(
                        base_den,
                        bins=bins,
                        rebin_factor=rebin_factor,
                        name_hint=f"{sanitize(job_name)}_{sanitize(era_key)}_{sanitize(base_tag)}_base_den",
                    )

                    if not ROOT.TEfficiency.CheckConsistency(base_num, base_den):
                        raise RuntimeError(f"TEfficiency consistency failed for base tag '{base_tag}'.")
                    base_eff = ROOT.TEfficiency(base_num, base_den)

                    ordered_filters = list(trigger_setting.get(base_tag, []))
                    if not ordered_filters:
                        raise RuntimeError(f"No configured filter chain found for base tag '{base_tag}'.")

                    prev_num = None
                    abs_effs = {}
                    step_effs = {}
                    filter_summary = []

                    for filter_name in ordered_filters:
                        full_tag = f"{base_tag}_{filter_name}"
                        _, filter_num_name, _, _ = build_hist_names(resonance, job, full_tag)
                        filter_num, filter_runs = aggregate_histogram_for_era(
                            era=era_key,
                            source=source,
                            hist_path_template=path_template,
                            fmt_args={"target_dir": target_dir, "hist": filter_num_name},
                            strict=strict,
                        )
                        if filter_num is None:
                            continue

                        filter_num = rebin_histogram(
                            filter_num,
                            bins=bins,
                            rebin_factor=rebin_factor,
                            name_hint=f"{sanitize(job_name)}_{sanitize(era_key)}_{sanitize(full_tag)}_num",
                        )

                        if not ROOT.TEfficiency.CheckConsistency(filter_num, base_den):
                            continue
                        eff_abs = ROOT.TEfficiency(filter_num, base_den)
                        abs_effs[filter_name] = eff_abs

                        step_den = base_den if prev_num is None else prev_num
                        if ROOT.TEfficiency.CheckConsistency(filter_num, step_den):
                            step_effs[filter_name] = ROOT.TEfficiency(filter_num, step_den)

                        prev_num = filter_num
                        filter_summary.append(
                            {
                                "filter": filter_name,
                                "full_tag": full_tag,
                                "used_runs": filter_runs,
                            }
                        )

                        if save_era_comparison:
                            era_comparison_store.setdefault(base_tag, {}).setdefault(filter_name, {"abs": {}, "step": {}})
                            era_comparison_store[base_tag][filter_name]["abs"][source_display_label(source)] = eff_abs
                            if filter_name in step_effs:
                                era_comparison_store[base_tag][filter_name]["step"][source_display_label(source)] = step_effs[filter_name]

                        if save_run_trends:
                            run_trend_store.setdefault(era_key, {}).setdefault(base_tag, {"abs": {}, "step": {}})
                            abs_points = []
                            step_points = []
                            for run in sorted(source["run_files"].keys()):
                                run_base_den = aggregate_histogram_for_run(
                                    source=source,
                                    run=run,
                                    hist_path_template=path_template,
                                    fmt_args={"target_dir": target_dir, "hist": denominator_name},
                                    strict=False,
                                )
                                run_filter_num = aggregate_histogram_for_run(
                                    source=source,
                                    run=run,
                                    hist_path_template=path_template,
                                    fmt_args={"target_dir": target_dir, "hist": filter_num_name},
                                    strict=False,
                                )
                                if run_base_den is None or run_filter_num is None:
                                    continue
                                run_base_den = rebin_histogram(
                                    run_base_den,
                                    bins=bins,
                                    rebin_factor=rebin_factor,
                                    name_hint=f"run_{sanitize(str(run))}_{sanitize(full_tag)}_base_den",
                                )
                                run_filter_num = rebin_histogram(
                                    run_filter_num,
                                    bins=bins,
                                    rebin_factor=rebin_factor,
                                    name_hint=f"run_{sanitize(str(run))}_{sanitize(full_tag)}_num",
                                )
                                run_abs_mean = mean_efficiency(run_filter_num, run_base_den)
                                if run_abs_mean is not None:
                                    abs_points.append((run, run_abs_mean))

                                step_run_den = run_base_den
                                prev_full_tag = None
                                prev_idx = ordered_filters.index(filter_name) - 1
                                if prev_idx >= 0:
                                    prev_full_tag = f"{base_tag}_{ordered_filters[prev_idx]}"
                                    _, prev_num_name, _, _ = build_hist_names(resonance, job, prev_full_tag)
                                    prev_run_num = aggregate_histogram_for_run(
                                        source=source,
                                        run=run,
                                        hist_path_template=path_template,
                                        fmt_args={"target_dir": target_dir, "hist": prev_num_name},
                                        strict=False,
                                    )
                                    if prev_run_num is not None:
                                        step_run_den = rebin_histogram(
                                            prev_run_num,
                                            bins=bins,
                                            rebin_factor=rebin_factor,
                                            name_hint=f"run_{sanitize(str(run))}_{sanitize(prev_full_tag)}_num",
                                        )
                                run_step_mean = mean_efficiency(run_filter_num, step_run_den)
                                if run_step_mean is not None:
                                    step_points.append((run, run_step_mean))

                            run_trend_store[era_key][base_tag]["abs"][filter_name] = abs_points
                            run_trend_store[era_key][base_tag]["step"][filter_name] = step_points

                    if not abs_effs:
                        raise RuntimeError(f"No usable filter histograms found for base tag '{base_tag}'.")

                    out_base = out_dir / f"{sanitize(job_name)}_{sanitize(era_key)}_{sanitize(base_tag)}"
                    draw_filter_canvas(
                        cfg=cfg,
                        source=source,
                        era_key=era_key,
                        job=job,
                        base_tag=base_tag,
                        base_eff=base_eff,
                        abs_effs=abs_effs,
                        step_effs=step_effs,
                        out_base=out_base,
                        ratio_ymin=ratio_ymin,
                        ratio_ymax=ratio_ymax,
                        legend_cols=legend_cols,
                    )

                    era_result["status"] = "ok"
                    era_result["base_tags"][base_tag] = {
                        "used_runs_num": base_num_runs,
                        "used_runs_den": base_den_runs,
                        "filters": filter_summary,
                        "n_filters_drawn": len(abs_effs),
                    }
            except Exception as exc:
                print(f"[{MODULE_NAME}][WARN] job={job_name} era={era_key}: {exc}")
                if strict:
                    raise
            finally:
                summary[job_name][era_key] = era_result
                if progress is not None and era_task is not None:
                    progress.update(era_task, advance=1)

        usable = sum(1 for era_result in summary[job_name].values() if era_result.get("status") == "ok")
        skipped = len(summary[job_name]) - usable
        if usable > 0:
            emit_log(progress, f"[{MODULE_NAME}] job done: {job_name} usable_eras={usable} skipped_eras={skipped}", style="blue")
        else:
            emit_log(progress, f"[{MODULE_NAME}][WARN] job done: {job_name} usable_eras=0 skipped_eras={skipped}", style="yellow")

        if save_era_comparison:
            for base_tag, filter_map in era_comparison_store.items():
                for filter_name, payload in filter_map.items():
                    if payload["abs"]:
                        out_base = out_dir / f"{sanitize(job_name)}_{sanitize(base_tag)}_{sanitize(filter_name)}_era_compare"
                        draw_filter_era_comparison_canvas(
                            cfg=cfg,
                            job=job,
                            base_tag=base_tag,
                            filter_name=filter_name,
                            abs_by_era=payload["abs"],
                            step_by_era=payload["step"],
                            out_base=out_base,
                            legend_cols=legend_cols,
                        )

        if save_run_trends:
            for era_key, era_map in run_trend_store.items():
                source = era_sources[era_key]
                for base_tag, payload in era_map.items():
                    if payload["abs"]:
                        out_base_abs = out_dir / f"{sanitize(job_name)}_{sanitize(era_key)}_{sanitize(base_tag)}_run_trend_abs"
                        draw_run_trend_canvas(
                            cfg=cfg,
                            source=source,
                            base_tag=base_tag,
                            title_suffix="abs",
                            trend_data=payload["abs"],
                            out_base=out_base_abs,
                            legend_cols=legend_cols,
                        )
                    if payload["step"]:
                        out_base_step = out_dir / f"{sanitize(job_name)}_{sanitize(era_key)}_{sanitize(base_tag)}_run_trend_step"
                        draw_run_trend_canvas(
                            cfg=cfg,
                            source=source,
                            base_tag=base_tag,
                            title_suffix="step",
                            trend_data=payload["step"],
                            out_base=out_base_step,
                            legend_cols=legend_cols,
                        )

        if progress is not None and job_task is not None:
            progress.update(job_task, advance=1)

    summary_file = out_dir / f"{MODULE_NAME}_summary.yaml"
    with open(summary_file, "w", encoding="utf-8") as handle:
        yaml.safe_dump(summary, handle, sort_keys=False)

    emit_log(progress, f"[{MODULE_NAME}] done. outputs={out_dir}", style="blue")
    return summary
