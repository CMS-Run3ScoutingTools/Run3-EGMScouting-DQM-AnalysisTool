from pathlib import Path

import ROOT
import yaml
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep

from dqm_pipeline.core import (
    aggregate_histogram_for_era,
    emit_log,
    load_histogram,
    rebin_histogram,
    sanitize,
)


MODULE_NAME = "tnp_efficiency"


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


def plot_tnp_efficiency(job, tag, per_era_points, out_png, mc_points=None):
    fig, ax = plt.subplots(figsize=(8.5, 5.2))

    for era, pts in per_era_points.items():
        ax.errorbar(
            pts["x"],
            pts["y"],
            xerr=pts["xerr"],
            yerr=np.vstack([pts["yerr_low"], pts["yerr_up"]]),
            fmt="o-",
            markersize=4,
            linewidth=1.4,
            capsize=2,
            label=era,
        )

    if mc_points is not None:
        ax.errorbar(
            mc_points["x"],
            mc_points["y"],
            xerr=mc_points["xerr"],
            yerr=np.vstack([mc_points["yerr_low"], mc_points["yerr_up"]]),
            fmt="s--",
            markersize=4,
            linewidth=1.4,
            capsize=2,
            color="black",
            label="MC",
        )

    axis = job.get("axis", "pt").lower()
    pt_order = job.get("pt_order", "leading")
    x_label = job.get("x_label", f"{pt_order} pT [GeV]" if axis == "pt" else "eta")

    ax.set_xlabel(x_label)
    ax.set_ylabel("Efficiency")
    ax.set_ylim(job.get("ymin", 0.0), job.get("ymax", 1.05))
    ax.grid(True, alpha=0.3)

    title = job.get("title", job.get("name", "tnp"))
    ax.set_title(f"{title} | pass{tag}")
    ax.legend(fontsize=9)

    hep.cms.text("Preliminary", loc=2, ax=ax, fontsize=12)
    hep.cms.lumitext("2025 (13.6 TeV)", ax=ax)

    fig.tight_layout()
    fig.savefig(out_png, dpi=140)
    plt.close(fig)


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

    all_results = {}
    job_task = None
    if progress is not None:
        job_task = progress.add_task(f"[blue]{MODULE_NAME}: jobs", total=len(jobs))

    for job in jobs:
        job_name = job["name"]
        emit_log(progress, f"[{MODULE_NAME}] job start: {job_name}", style="blue")
        tags = job.get("tags", [])
        axis = job.get("axis", "pt").lower()
        bins = job.get("bins", bins_cfg.get(axis))
        rebin_factor = int(job.get("rebin", default_rebin))

        all_results[job_name] = {}
        tag_task = None
        if progress is not None:
            tag_task = progress.add_task(f"[blue]{MODULE_NAME}: {job_name} tags", total=len(tags))

        for tag in tags:
            emit_log(progress, f"[{MODULE_NAME}] tag start: {job_name}/{tag}", style="blue")
            target_dir, num_name, den_name = build_tnp_hist_names(resonance, job, tag)
            all_results[job_name][tag] = {}
            per_era_points = {}
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

                    points = compute_efficiency_points(num_hist, den_hist)
                    per_era_points[era] = points

                    all_results[job_name][tag][era] = {
                        "numerator": num_name,
                        "denominator": den_name,
                        "used_runs_num": used_runs_num,
                        "used_runs_den": used_runs_den,
                        "mean_efficiency": float(np.mean(points["y"])),
                    }
                except Exception as exc:
                    print(f"[{MODULE_NAME}][WARN] job={job_name} tag={tag} era={era}: {exc}")
                    if strict:
                        raise
                finally:
                    if progress is not None and era_task is not None:
                        progress.update(era_task, advance=1)

            mc_points = None
            if job.get("include_mc", False) and mc_file:
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

                    mc_points = compute_efficiency_points(num_hist_mc, den_hist_mc)
                    all_results[job_name][tag]["MC"] = {
                        "numerator": num_name,
                        "denominator": den_name,
                        "mean_efficiency": float(np.mean(mc_points["y"])),
                    }
                except Exception as exc:
                    print(f"[{MODULE_NAME}][WARN] job={job_name} tag={tag} mc: {exc}")
                    if strict:
                        raise

            if per_era_points or mc_points is not None:
                out_png = out_dir / f"{sanitize(job_name)}_{sanitize(tag)}.png"
                plot_tnp_efficiency(
                    job=job,
                    tag=tag,
                    per_era_points=per_era_points,
                    out_png=str(out_png),
                    mc_points=mc_points,
                )
            emit_log(progress, f"[{MODULE_NAME}] tag done: {job_name}/{tag}", style="blue")
            if progress is not None and tag_task is not None:
                progress.update(tag_task, advance=1)
        emit_log(progress, f"[{MODULE_NAME}] job done: {job_name}", style="blue")
        if progress is not None and job_task is not None:
            progress.update(job_task, advance=1)

    summary_file = out_dir / "tnp_efficiency_summary.yaml"
    with open(summary_file, "w", encoding="utf-8") as f:
        yaml.safe_dump(all_results, f, sort_keys=False)

    emit_log(progress, f"[{MODULE_NAME}] done. Outputs in {out_dir}", style="bold blue")
    return all_results
