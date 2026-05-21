from pathlib import Path

import yaml
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
from scipy.optimize import curve_fit

from dqm_pipeline.core import aggregate_histogram_for_era, sanitize, emit_log


MODULE_NAME = "mass_fit"


def source_display_label(source):
    return str(source.get("display_label", source.get("era", "era")))


def global_plot_lumi_text(cfg):
    plotting = cfg.get("plotting", {})
    if plotting.get("lumi_text"):
        return str(plotting["lumi_text"])
    energy = plotting.get("energy_tev", cfg.get("energy_tev", 13.6))
    campaign = plotting.get("campaign_label", cfg.get("campaign_label", "Run 3"))
    return f"{campaign} ({energy} TeV)"


def source_plot_lumi_text(cfg, source):
    if source.get("plot_lumi_text"):
        return str(source["plot_lumi_text"])
    energy = cfg.get("plotting", {}).get("energy_tev", cfg.get("energy_tev", 13.6))
    label = source_display_label(source)
    if source.get("lumi_fb") is not None:
        return f"{label} {float(source['lumi_fb']):.2f} $\\mathrm{{fb}}^{{-1}}$ ({energy} TeV)"
    return f"{label} ({energy} TeV)"


def double_crystal_ball(x, mean, sigma, alpha_l, n_l, alpha_r, n_r, amp):
    x = np.asarray(x)
    t = (x - mean) / sigma

    alpha_l = abs(alpha_l) if abs(alpha_l) > 1e-6 else 1e-6
    alpha_r = abs(alpha_r) if abs(alpha_r) > 1e-6 else 1e-6
    n_l = max(n_l, 1.01)
    n_r = max(n_r, 1.01)

    left_a = (n_l / alpha_l) ** n_l * np.exp(-0.5 * alpha_l**2)
    left_b = n_l / alpha_l - alpha_l
    right_a = (n_r / alpha_r) ** n_r * np.exp(-0.5 * alpha_r**2)
    right_b = n_r / alpha_r - alpha_r

    out = np.exp(-0.5 * t**2)
    left_mask = t < -alpha_l
    right_mask = t > alpha_r
    out[left_mask] = left_a * np.power(left_b - t[left_mask], -n_l)
    out[right_mask] = right_a * np.power(right_b + t[right_mask], -n_r)
    return amp * out


def background_curve(x, params, bkg_type):
    x = np.asarray(x, dtype=float)
    if bkg_type == "none":
        return np.zeros_like(x, dtype=float)
    if bkg_type == "linear":
        offset, slope = params
        return np.clip(offset + slope * x, 0.0, None)
    if bkg_type == "constant":
        return np.full_like(x, params[0], dtype=float)
    if bkg_type == "exp":
        amp, slope = params
        return amp * np.exp(slope * x)
    if bkg_type == "quintic":
        return sum(coeff * np.power(x, order) for order, coeff in enumerate(params))
    raise RuntimeError(f"Unsupported mass_fit background type '{bkg_type}'.")


def n_background_params(bkg_type):
    return {"none": 0, "constant": 1, "linear": 2, "exp": 2, "quintic": 6}[bkg_type]


def model_mass_dcb(x, *params, bkg_type="exp"):
    n_bkg = n_background_params(bkg_type)
    bkg_params = params[:n_bkg]
    mean, sigma, alpha_l, n_l, alpha_r, n_r, amp = params[n_bkg:]
    return background_curve(x, bkg_params, bkg_type) + double_crystal_ball(
        x,
        mean=mean,
        sigma=sigma,
        alpha_l=alpha_l,
        n_l=n_l,
        alpha_r=alpha_r,
        n_r=n_r,
        amp=amp,
    )


def initial_dcb_parameters(x_vals, y_vals, xmin, xmax, bkg_type):
    peak_idx = int(np.argmax(y_vals))
    peak_x = float(x_vals[peak_idx])
    peak_y = float(y_vals[peak_idx])
    bkg_guess = max(float(np.percentile(y_vals, 10)), 1.0)
    sigma_guess = max(0.03 * (xmax - xmin), 0.05)

    if bkg_type == "none":
        bkg_p0 = []
        bkg_lo = []
        bkg_hi = []
    elif bkg_type == "constant":
        bkg_p0 = [bkg_guess]
        bkg_lo = [0.0]
        bkg_hi = [np.inf]
    elif bkg_type == "linear":
        bkg_p0 = [bkg_guess, 0.0]
        bkg_lo = [0.0, -np.inf]
        bkg_hi = [np.inf, np.inf]
    elif bkg_type == "exp":
        bkg_p0 = [bkg_guess, 0.0]
        bkg_lo = [0.0, -1.0]
        bkg_hi = [np.inf, 1.0]
    elif bkg_type == "quintic":
        bkg_p0 = [bkg_guess, 0.0, 0.0, 0.0, 0.0, 0.0]
        bkg_lo = [-np.inf] * 6
        bkg_hi = [np.inf] * 6
    else:
        raise RuntimeError(f"Unsupported mass_fit background type '{bkg_type}'.")

    signal_p0 = [
        min(max(peak_x, xmin), xmax),
        sigma_guess,
        1.5,
        3.0,
        1.5,
        3.0,
        max(peak_y - bkg_guess, peak_y, 1.0),
    ]
    signal_lo = [xmin, 0.01, 0.05, 1.01, 0.05, 1.01, 0.0]
    signal_hi = [xmax, max(xmax - xmin, 0.1), 10.0, 200.0, 10.0, 200.0, np.inf]
    return bkg_p0 + signal_p0, bkg_lo + signal_lo, bkg_hi + signal_hi


def extract_points_from_hist(hist, xmin, xmax, error_model="sqrt_y"):
    x_vals, y_vals, y_errs = [], [], []
    for i_bin in range(1, hist.GetNbinsX() + 1):
        x = hist.GetBinCenter(i_bin)
        y = hist.GetBinContent(i_bin)
        if y > 0 and xmin < x < xmax:
            if error_model == "hist":
                err = hist.GetBinError(i_bin)
                if err <= 0:
                    err = np.sqrt(max(y, 0.0))
            else:
                # Original script behavior.
                err = np.sqrt(max(y, 0.0))
            x_vals.append(x)
            y_vals.append(y)
            y_errs.append(err)

    return np.asarray(x_vals), np.asarray(y_vals), np.asarray(y_errs)


def extract_points_with_target_nbins(hist, xmin, xmax, target_nbins, error_model="sqrt_y"):
    if int(target_nbins) <= 0:
        raise RuntimeError(f"target_nbins must be > 0, got {target_nbins}")

    n_bins = int(target_nbins)
    edges = np.linspace(float(xmin), float(xmax), n_bins + 1)
    sums = np.zeros(n_bins, dtype=float)
    errs2 = np.zeros(n_bins, dtype=float)

    for i_bin in range(1, hist.GetNbinsX() + 1):
        x = hist.GetBinCenter(i_bin)
        y = hist.GetBinContent(i_bin)
        if y <= 0 or x <= xmin or x >= xmax:
            continue
        idx = np.searchsorted(edges, x, side="right") - 1
        if idx < 0 or idx >= n_bins:
            continue
        if error_model == "hist":
            err = hist.GetBinError(i_bin)
            if err <= 0:
                err = np.sqrt(max(y, 0.0))
        else:
            err = np.sqrt(max(y, 0.0))
        sums[idx] += y
        errs2[idx] += err * err

    centers = 0.5 * (edges[:-1] + edges[1:])
    errs = np.sqrt(errs2)
    mask = sums > 0
    return centers[mask], sums[mask], errs[mask]


def fit_histogram(
    hist,
    xmin,
    xmax,
    era,
    era_label,
    out_png,
    rebin_factor=1,
    target_nbins=None,
    error_model="sqrt_y",
    plot_lumi_text=None,
    bkg_type="exp",
):
    if target_nbins is not None:
        x_vals, y_vals, y_errs = extract_points_with_target_nbins(
            hist=hist,
            xmin=xmin,
            xmax=xmax,
            target_nbins=int(target_nbins),
            error_model=error_model,
        )
    else:
        fit_hist = hist
        if int(rebin_factor) > 1:
            fit_hist = hist.Clone(
                f"{hist.GetName()}_{sanitize(era)}_fit_rebin_{int(rebin_factor)}_{sanitize(str(xmin))}_{sanitize(str(xmax))}"
            )
            fit_hist.SetDirectory(0)
            fit_hist.Rebin(int(rebin_factor))
        x_vals, y_vals, y_errs = extract_points_from_hist(
            hist=fit_hist,
            xmin=xmin,
            xmax=xmax,
            error_model=error_model,
        )

    if len(x_vals) < 7:
        raise RuntimeError(f"Not enough points to fit in [{xmin}, {xmax}] for era {era}.")

    bkg_type = str(bkg_type).lower()
    p0, lower_bounds, upper_bounds = initial_dcb_parameters(x_vals, y_vals, xmin, xmax, bkg_type)

    popt, pcov = curve_fit(
        lambda x, *params: model_mass_dcb(x, *params, bkg_type=bkg_type),
        x_vals,
        y_vals,
        p0=p0,
        sigma=y_errs,
        absolute_sigma=True,
        maxfev=10000,
        bounds=(lower_bounds, upper_bounds),
    )
    perr = np.sqrt(np.diag(pcov))

    x_fit = np.linspace(x_vals.min(), x_vals.max(), 1000)
    n_bkg = n_background_params(bkg_type)
    bkg_curve = background_curve(x_fit, popt[:n_bkg], bkg_type)
    mean, sigma, alpha_l, n_l, alpha_r, n_r, amp = popt[n_bkg:]
    sig_curve = double_crystal_ball(x_fit, mean, sigma, alpha_l, n_l, alpha_r, n_r, amp)
    tot_curve = bkg_curve + sig_curve

    bkg_integral = float(np.trapz(bkg_curve, x_fit))
    sig_integral = float(np.trapz(sig_curve, x_fit))

    fig, ax = plt.subplots()
    ax.errorbar(x_vals, y_vals, yerr=y_errs, fmt="ko", label="Data", markersize=3)
    ax.plot(x_fit, tot_curve, "m-", label="Signal + Background")
    ax.plot(x_fit, bkg_curve, color="brown", label="Background")
    ax.plot(x_fit, sig_curve, color="cyan", label="Signal(Double Crystal Ball)")
    ax.set_xlabel("Dielectron mass [GeV]")
    ax.set_ylabel("Events / 1 GeV")
    ax.set_xlim(float(xmin), float(xmax))
    ax.grid()
    ax.legend(loc="upper right", frameon=False)

    hep.cms.text("Preliminary", loc=2, ax=ax, fontsize=12)
    hep.cms.lumitext(plot_lumi_text or f"{era_label}", ax=ax)

    plt.text(
        0.05,
        0.80,
        f"background integral: {bkg_integral:.0f}",
        transform=plt.gca().transAxes,
        fontsize=10,
        verticalalignment="top",
        horizontalalignment="left",
        color="black",
    )
    plt.text(
        0.05,
        0.75,
        f"signal integral: {sig_integral:.0f}",
        transform=plt.gca().transAxes,
        fontsize=10,
        verticalalignment="top",
        horizontalalignment="left",
        color="black",
    )
    plt.text(
        0.05,
        0.70,
        f"mean: {mean:.1f} GeV",
        transform=plt.gca().transAxes,
        fontsize=10,
        verticalalignment="top",
        horizontalalignment="left",
        color="black",
    )
    plt.text(
        0.05,
        0.65,
        f"$\\sigma_{{dcb}}$: {sigma:.1f} GeV",
        transform=plt.gca().transAxes,
        fontsize=10,
        verticalalignment="top",
        horizontalalignment="left",
        color="black",
    )
    plt.text(
        0.05,
        0.60,
        f"relative width: {100.0 * sigma / mean:.1f} %",
        transform=plt.gca().transAxes,
        fontsize=10,
        verticalalignment="top",
        horizontalalignment="left",
        color="black",
    )

    fig.tight_layout()
    fig.savefig(out_png)
    plt.close(fig)

    bkg_names = {
        "none": [],
        "constant": ["bkg_amp"],
        "linear": ["bkg_offset", "bkg_slope"],
        "exp": ["bkg_amp", "bkg_slope"],
        "quintic": ["bkg_c0", "bkg_c1", "bkg_c2", "bkg_c3", "bkg_c4", "bkg_c5"],
    }[bkg_type]
    names = bkg_names + ["dcb_mean", "dcb_sigma", "dcb_alpha_l", "dcb_n_l", "dcb_alpha_r", "dcb_n_r", "dcb_amp"]
    out = {}
    for name, val, err in zip(names, popt, perr):
        out[name] = {"value": float(val), "err": float(err)}
    out["cb_mean"] = dict(out["dcb_mean"])
    out["cb_sigma"] = dict(out["dcb_sigma"])
    out["fit_model"] = "double_crystal_ball"
    out["background_type"] = bkg_type
    out["background_integral"] = bkg_integral
    out["signal_integral"] = sig_integral
    return out


def era_label_for_plots(source):
    if source.get("lumi_fb") is not None:
        return f"{float(source['lumi_fb']):.2f} fb$^{{-1}}$"
    if source.get("selected_lumisections", 0) > 0:
        return f"{source['n_runs']} runs, {source['selected_lumisections']} golden LS"
    return f"{source['n_runs']} runs"


def overlay_scale_value(source, mode):
    if mode == "lumi_fb":
        val = source.get("lumi_fb")
        return float(val) if val is not None else None
    if mode == "golden_lumisections":
        val = source.get("selected_lumisections", 0)
        return float(val) if val > 0 else None
    return None


def plot_mass_overlay(variable, era_hists, era_sources, out_png, scale_mode="none", plot_lumi_text="Run 3 (13.6 TeV)"):
    fig, ax = plt.subplots(figsize=(9.2, 6.0))

    for era, hist in era_hists.items():
        source = era_sources[era]
        shown = hist.Clone(f"{hist.GetName()}_{sanitize(era)}_overlay")
        scale = overlay_scale_value(source, scale_mode)
        if scale is not None and scale > 0:
            shown.Scale(1.0 / scale)

        x = np.array([shown.GetBinCenter(i) for i in range(1, shown.GetNbinsX() + 1)])
        y = np.array([shown.GetBinContent(i) for i in range(1, shown.GetNbinsX() + 1)])
        mask = (x > 0) & (y > 0)

        ax.step(
            x[mask],
            y[mask],
            where="mid",
            linewidth=1.6,
            label=f"{source_display_label(source)} [{era_label_for_plots(source)}]",
        )

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("M$_{ee}$ [GeV]", fontsize=22)

    if scale_mode == "lumi_fb":
        ax.set_ylabel("Events / fb", fontsize=22)
    elif scale_mode == "golden_lumisections":
        ax.set_ylabel("Events / golden LS", fontsize=22)
    else:
        ax.set_ylabel("Events", fontsize=22)

    ax.tick_params(axis="both", which="major", labelsize=16, direction="in", top=True, right=True, length=10)
    ax.tick_params(axis="both", which="minor", direction="in", top=True, right=True, length=5)
    ax.grid(True, which="both", alpha=0.25, linestyle=":")

    ax.text(
        0.00,
        1.03,
        "CMS",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=20,
        fontweight="bold",
    )
    ax.text(
        0.10,
        1.03,
        "Preliminary",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=17,
        style="italic",
    )
    ax.text(
        1.00,
        1.03,
        plot_lumi_text,
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=14,
    )

    ax.legend(
        loc="upper right",
        frameon=False,
        fontsize=11,
        handlelength=2.2,
        labelspacing=0.6,
    )
    fig.subplots_adjust(left=0.14, right=0.99, bottom=0.13, top=0.88)
    fig.savefig(out_png, dpi=140)
    plt.close(fig)


def plot_mass_by_era(results_by_era, out_png, ymin, ymax, title, plot_lumi_text="Run 3 (13.6 TeV)"):
    eras = list(results_by_era.keys())
    labels = [results_by_era[e].get("display_label", e) for e in eras]
    masses = [results_by_era[e]["mass"]["value"] for e in eras]
    # Keep original plotting behavior: "fitted error" came from width fit error.
    mass_errs = [results_by_era[e]["width"]["err"] for e in eras]
    widths = [results_by_era[e]["width"]["value"] for e in eras]

    fig, ax = plt.subplots(figsize=(8, 4.8))
    ax.errorbar(labels, masses, yerr=widths, color="tab:blue", fmt="s", capsize=4, label="Fitted mass (w width)")
    ax.errorbar(labels, masses, yerr=mass_errs, color="tab:red", fmt="s", capsize=4, label="Fitted mass (w fitted error)")
    ax.set_ylabel("Mass [GeV]")
    ax.set_ylim(ymin, ymax)
    # ax.set_title(title)
    ax.grid(True, which="both", axis="x", linestyle="--", alpha=0.5)
    ax.legend(loc="upper right", frameon=False)
    hep.cms.text("Preliminary", loc=2, ax=ax, fontsize=12)
    hep.cms.lumitext(plot_lumi_text, ax=ax)
    fig.tight_layout()
    fig.savefig(out_png, dpi=140)
    plt.close(fig)


def run_module(cfg, era_sources, out_root, strict=False, progress=None):
    section = cfg.get(MODULE_NAME)
    if not section or not section.get("enabled", True):
        return {}

    fit_windows = section["fit_windows"]
    variables = section["variables"]
    hist_tpl = section["hist_path_template"]
    default_rebin = int(section.get("rebin", 1))
    error_model = section.get("error_model", "sqrt_y")  # sqrt_y | hist
    default_bkg_type = str(section.get("bkg_type", section.get("background_type", "exp"))).lower()
    overlay_rebin = int(section.get("overlay_rebin", default_rebin))
    overlay_scale_mode = section.get("overlay_scale", "none")
    summary_lumi_text = global_plot_lumi_text(cfg)

    out_dir = Path(out_root) / section.get("output_subdir", MODULE_NAME)
    out_dir.mkdir(parents=True, exist_ok=True)
    emit_log(
        progress,
        f"[{MODULE_NAME}] start variables={len(variables)} eras={len(era_sources)} out_dir={out_dir}",
        style="bold magenta",
    )

    all_results = {}
    total_successful_eras = 0
    total_failed_eras = 0
    var_task = None
    if progress is not None:
        var_task = progress.add_task(f"[magenta]{MODULE_NAME}: variables", total=len(variables))

    for variable in variables:
        emit_log(progress, f"[{MODULE_NAME}] variable start: {variable}", style="magenta")
        era_hists = {}
        variable_results = {w: {} for w in fit_windows}
        variable_successful_eras = 0
        variable_failed_eras = 0
        era_task = None
        era_items = list(era_sources.items())
        if progress is not None:
            era_task = progress.add_task(f"[magenta]{MODULE_NAME}: {variable} eras", total=len(era_items))

        for era, source in era_items:
            try:
                hist, used_runs = aggregate_histogram_for_era(
                    era=era,
                    source=source,
                    hist_path_template=hist_tpl,
                    fmt_args={"var": variable},
                    strict=strict,
                )
                if hist is None:
                    raise RuntimeError("No histogram found in selected runs.")

                hist_for_overlay = hist
                if overlay_rebin > 1:
                    hist_for_overlay = hist.Clone(
                        f"{hist.GetName()}_{sanitize(era)}_overlay_rebin_{overlay_rebin}"
                    )
                    hist_for_overlay.SetDirectory(0)
                    hist_for_overlay.Rebin(overlay_rebin)

                era_hists[era] = hist_for_overlay
                era_text = era_label_for_plots(source)
                individual_lumi_text = source_plot_lumi_text(cfg, source)

                for fit_name, win in fit_windows.items():
                    win_rebin = int(win.get("rebin", default_rebin))
                    win_nbins = win.get("nbins")
                    if win_nbins is not None:
                        win_nbins = int(win_nbins)
                    win_bkg_type = str(win.get("bkg_type", win.get("background_type", default_bkg_type))).lower()

                    fit_png = out_dir / f"{era}_{sanitize(variable)}_{fit_name}.png"
                    fit_out = fit_histogram(
                        hist=hist,
                        xmin=float(win["xmin"]),
                        xmax=float(win["xmax"]),
                        era=era,
                        era_label=era_text,
                        out_png=str(fit_png),
                        rebin_factor=win_rebin,
                        target_nbins=win_nbins,
                        error_model=error_model,
                        plot_lumi_text=individual_lumi_text,
                        bkg_type=win_bkg_type,
                    )
                    fit_binning = (
                        {"mode": "nbins", "value": win_nbins}
                        if win_nbins is not None
                        else {"mode": "rebin", "value": win_rebin}
                    )
                    variable_results[fit_name][era] = {
                        "display_label": source_display_label(source),
                        "mass": fit_out["cb_mean"],
                        "width": fit_out["cb_sigma"],
                        "fit_binning": fit_binning,
                        "used_runs": used_runs,
                        "full": fit_out,
                    }
                variable_successful_eras += 1
            except Exception as exc:
                variable_failed_eras += 1
                print(f"[{MODULE_NAME}][WARN] era={era} variable={variable}: {exc}")
                if strict:
                    raise
            finally:
                if progress is not None and era_task is not None:
                    progress.update(era_task, advance=1)

        if era_hists:
            overlay_png = out_dir / f"overlay_{sanitize(variable)}.png"
            plot_mass_overlay(
                variable=variable,
                era_hists=era_hists,
                era_sources=era_sources,
                out_png=str(overlay_png),
                scale_mode=overlay_scale_mode,
                plot_lumi_text=summary_lumi_text,
            )
        else:
            emit_log(
                progress,
                f"[{MODULE_NAME}][WARN] variable={variable}: no usable eras, overlay skipped",
                style="yellow",
            )

        for fit_name, era_result in variable_results.items():
            if not era_result:
                continue
            win = fit_windows[fit_name]
            mass_png = out_dir / f"mass_by_era_{sanitize(variable)}_{fit_name}.png"
            plot_mass_by_era(
                results_by_era=era_result,
                out_png=str(mass_png),
                ymin=float(win["ymin"]),
                ymax=float(win["ymax"]),
                title=f"{variable} [{fit_name}]",
                plot_lumi_text=summary_lumi_text,
            )

        all_results[variable] = {
            "status": "ok" if variable_successful_eras > 0 else "empty",
            "n_successful_eras": variable_successful_eras,
            "n_failed_eras": variable_failed_eras,
            "fits": variable_results,
        }
        total_successful_eras += variable_successful_eras
        total_failed_eras += variable_failed_eras
        if variable_successful_eras > 0:
            emit_log(
                progress,
                f"[{MODULE_NAME}] variable done: {variable} usable_eras={variable_successful_eras} skipped_eras={variable_failed_eras}",
                style="magenta",
            )
        else:
            emit_log(
                progress,
                f"[{MODULE_NAME}][WARN] variable done: {variable} usable_eras=0 skipped_eras={variable_failed_eras}",
                style="yellow",
            )
        if progress is not None and var_task is not None:
            progress.update(var_task, advance=1)

    summary_file = out_dir / "mass_fit_summary.yaml"
    with open(summary_file, "w", encoding="utf-8") as f:
        yaml.safe_dump(all_results, f, sort_keys=False)

    final_style = "bold magenta" if total_successful_eras > 0 else "yellow"
    final_tag = "done" if total_successful_eras > 0 else "[WARN] done"
    emit_log(
        progress,
        f"[{MODULE_NAME}] {final_tag}. usable_eras={total_successful_eras} skipped_eras={total_failed_eras} outputs={out_dir}",
        style=final_style,
    )
    return all_results
