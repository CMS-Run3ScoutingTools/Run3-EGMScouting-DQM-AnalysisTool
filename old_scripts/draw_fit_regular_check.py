#!/usr/bin/env python3

import ROOT
import cmsstyle as CMS
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.special import erf
import mplhep as hep


def hex_to_rgb(hex_color):
    hex_color = hex_color.lstrip("#")
    r = int(hex_color[0:2], 16) / 255.0
    g = int(hex_color[2:4], 16) / 255.0
    b = int(hex_color[4:6], 16) / 255.0
    return r, g, b


def fit_hist(hist, xmin, xmax, fig_name, era, lumi):
    x_vals = []
    y_vals = []
    y_errs = []

    for i in range(1, hist.GetNbinsX() + 1):
        x = hist.GetBinCenter(i)
        y = hist.GetBinContent(i)
        if (y > 0) and (x > xmin) and (x < xmax):
            x_vals.append(x)
            y_vals.append(y)
            y_errs.append(np.sqrt(y))

    x_vals = np.array(x_vals)
    y_vals = np.array(y_vals)
    y_errs = np.array(y_errs)

    def crystal_ball(x, alpha, n, mean, sigma, amp):
        x = np.array(x)
        t = (x - mean) / sigma
        abs_alpha = np.abs(alpha)
        if abs_alpha < 1e-3:
            abs_alpha = 1e-3
        if n < 1.01:
            n = 1.01

        a = (n / abs_alpha) ** n * np.exp(-0.5 * abs_alpha**2)
        b = n / abs_alpha - abs_alpha
        c = n / abs_alpha / (n - 1) * np.exp(-0.5 * abs_alpha**2)
        d = np.sqrt(np.pi / 2) * (1 + erf(abs_alpha / np.sqrt(2)))
        norm = 1.0 / (sigma * (c + d))

        result = np.zeros_like(t)
        mask = t > -abs_alpha
        result[mask] = np.exp(-0.5 * t[mask] ** 2)
        result[~mask] = a * np.power(b - t[~mask], -n)
        return amp * norm * result

    def model(x, bkg_amp, bkg_slope, cb_amp, cb_alpha, cb_n, cb_mean, cb_sigma):
        background = bkg_amp * np.exp(bkg_slope * x)
        cb = crystal_ball(x, cb_alpha, cb_n, cb_mean, cb_sigma, cb_amp)
        return background + cb

    p0 = [1e6, 0.0, 2e5, 1.5, 3.0, (xmax + xmin) / 2, 1.0]
    bounds = (
        [0, -1, 0, 0.5, 1.0, xmin, 0.01],
        [np.inf, 1, np.inf, 5.0, 20.0, xmax, 10.0],
    )

    popt, pcov = curve_fit(
        model,
        x_vals,
        y_vals,
        p0=p0,
        sigma=y_errs,
        absolute_sigma=True,
        maxfev=10000,
        bounds=bounds,
    )
    perr = np.sqrt(np.diag(pcov))

    x_fit = np.linspace(min(x_vals), max(x_vals), 1000)
    bkg_curve = popt[0] * np.exp(popt[1] * x_fit)
    sig_curve = crystal_ball(x_fit, *popt[3:7], popt[2])
    y_fit = bkg_curve + sig_curve

    fig, ax = plt.subplots()
    ax.errorbar(x_vals, y_vals, yerr=y_errs, fmt="ko", label="Data", markersize=3)
    ax.plot(x_fit, y_fit, "m-", label="Signal + Background")
    ax.plot(x_fit, bkg_curve, "brown", label="Background")
    ax.plot(x_fit, sig_curve, "cyan", label="Signal(Crystal Ball)")

    background_integral = np.trapz(bkg_curve, x_fit)
    signal_integral = np.trapz(sig_curve, x_fit)

    ax.set_xlabel("Dielectron mass [GeV]")
    ax.set_ylabel("Events / 1 GeV")
    ax.legend()
    ax.grid()

    hep.cms.text("Preliminary", loc=2, ax=ax, fontsize=12)
    hep.cms.lumitext(f"Run3 {era} {lumi} $\\mathrm{{fb}}^{{-1}}$ (13.6 TeV)", ax=ax)

    plt.text(0.05, 0.8, f"background integral: {background_integral:.0f}", transform=plt.gca().transAxes, fontsize=10)
    plt.text(0.05, 0.75, f"signal integral: {signal_integral:.0f}", transform=plt.gca().transAxes, fontsize=10)
    plt.text(0.05, 0.7, f"mean: {popt[5]:.1f} GeV", transform=plt.gca().transAxes, fontsize=10)
    plt.text(0.05, 0.65, f"$\\sigma_{{cb}}$: {popt[6]:.1f} GeV", transform=plt.gca().transAxes, fontsize=10)
    plt.text(0.05, 0.6, f"relative width: {100 * popt[6] / popt[5]:.1f} %", transform=plt.gca().transAxes, fontsize=10)

    plt.tight_layout()
    plt.savefig(f"{fig_name}.png")
    plt.close()

    param_names = ["bkg_amp", "bkg_slope", "cb_amp", "cb_alpha", "cb_n", "cb_mean", "cb_sigma"]
    output_dict = {}
    print("Fit results:")
    for name, val, err in zip(param_names, popt, perr):
        print(f"{name:>10} = {val:.6f} +/- {err:.6f}")
        output_dict[name] = {"value": val, "err": err}
    return output_dict


CMS.SetExtraText("Preliminary")
CMS.SetEnergy("2025, 13.6")
CMS.SetLumi("")

resonance = "resonanceZ"
plotdir = f"plot_DQM2026_check/{resonance}"
os.system(f"mkdir -p {plotdir}")
TAG = "EGM_Scouting_DQM_2026_v2"

lumi_dict = {
    "eraG": 21.7
}

ERA = list(lumi_dict.keys())[0]
file = ROOT.TFile.Open(f"/eos/user/t/tihsu/database/{TAG}/ntuple.root")


def Draw(lumi_dict_, name):
    c = None
    hist_dict = {}
    hex_colors = [
        "#3f90da", "#ffa90e", "#bd1f01", "#94a4a2", "#832db6",
        "#a96b59", "#e76300", "#b9ac70", "#717581", "#92dadd",
    ]
    default_colors = [ROOT.TColor.GetColor(hx) for hx in hex_colors]

    mass_width_result_Z = {}
    mass_width_result_JPsi = {}

    for idx, era in enumerate(lumi_dict_):
        lumi = lumi_dict_[era]
        file_ = ROOT.TFile.Open(f"/eos/user/t/tihsu/database/{TAG}/{era}/data.root")
        hist = file_.Get("DQMData/Run 333334/HLT/Run summary/ScoutingOffline/EGamma/Collection/" + name)
        if not hist:
            print(f"[WARN] missing hist for era={era} name={name}")
            file_.Close()
            continue

        hist_dict[era] = hist.Clone()
        hist_dict[era].SetDirectory(0)
        print(era)

        fit_results = fit_hist(
            hist_dict[era],
            xmin=70,
            xmax=110,
            fig_name=os.path.join(plotdir, f"{era}_{name}_Z"),
            era=era,
            lumi=lumi,
        )
        mass_width_result_Z[era] = {"mass": fit_results["cb_mean"], "width": fit_results["cb_sigma"]}

        fit_results = fit_hist(
            hist_dict[era],
            xmin=1,
            xmax=5,
            fig_name=os.path.join(plotdir, f"{era}_{name}_J"),
            era=era,
            lumi=lumi,
        )
        mass_width_result_JPsi[era] = {"mass": fit_results["cb_mean"], "width": fit_results["cb_sigma"]}

        hist_dict[era].Scale(1.0 / lumi)
        hist_dict[era] = hist_dict[era].Rebin(2)

        if c is None:
            c = CMS.cmsCanvas(
                "",
                2,
                140,
                1,
                hist_dict[era].GetMaximum() * 20,
                "M_{e,e} [GeV]",
                "events / int. lumi. [fb]",
                square=CMS.kRectangular,
                extraSpace=0.01,
                iPos=11,
                yTitOffset=0.8,
            )
            legend = CMS.cmsLeg(0.3, 0.65, 0.96, 0.92, textSize=0.033, columns=3)

        CMS.cmsDraw(
            hist_dict[era],
            "HIST",
            mcolor=default_colors[idx % len(default_colors)],
            fstyle=0,
            lwidth=2,
            lcolor=default_colors[idx % len(default_colors)],
            fcolor=0,
            msize=0,
        )
        legend.AddEntry(hist_dict[era], f"{era} [{lumi} fb^{{-1}} ]")
        file_.Close()

    if c is None:
        return

    c.SetLogx()
    c.SetLogy()
    hdf = CMS.GetcmsCanvasHist(c)
    hdf.GetYaxis().SetLabelSize(0.04)
    hdf.GetXaxis().SetLabelSize(0.04)
    hdf.SetMaximum(1e5)
    hdf.SetMinimum(1e-1)

    CMS.SaveCanvas(c, os.path.join(plotdir, f"{name}.png"), close=False)
    CMS.SaveCanvas(c, os.path.join(plotdir, f"{name}.pdf"))

    plot_mass_and_width_by_era(
        mass_width_result_Z,
        fig_name=os.path.join(plotdir, f"Z_mass_fitting_{name}"),
        ymin=80,
        ymax=100,
    )
    plot_mass_and_width_by_era(
        mass_width_result_JPsi,
        fig_name=os.path.join(plotdir, f"J_mass_fitting_{name}"),
        ymin=1,
        ymax=5,
    )


def plot_mass_and_width_by_era(mass_width_result, fig_name, ymin, ymax):
    eras = list(mass_width_result.keys())
    widths = [mass_width_result[era]["width"]["value"] for era in eras]
    masses = [mass_width_result[era]["mass"]["value"] for era in eras]
    mass_errs = [mass_width_result[era]["width"]["err"] for era in eras]

    fig, ax1 = plt.subplots()
    ax1.grid()
    ax2 = ax1.twinx()

    ax2.set_ylabel("Mass [GeV]", color="tab:red")
    ax2.errorbar(eras, masses, yerr=widths, color="tab:blue", fmt="s", label="fitted Mass (w width)", capsize=4)
    ax2.errorbar(eras, masses, yerr=mass_errs, color="tab:red", fmt="s", label="fitted Mass (w fitted error)", capsize=4)
    ax2.tick_params(axis="y", labelcolor="tab:red")
    ax2.set_ylim(ymin, ymax)

    hep.cms.text("Preliminary", loc=2, fontsize=12)
    hep.cms.lumitext("2025 (13.6 TeV)")
    fig.tight_layout()
    plt.grid(True, which="both", axis="x", linestyle="--", alpha=0.5)
    plt.savefig(f"{fig_name}.png")
    plt.close()


for var_ in [
    "sctElectron_EBEB_appliedID_invMass",
    "sctElectron_EBEB_appliedID_invMass_pass_DST_PFScouting_DoubleEG_v",
]:
    Draw(lumi_dict, var_)

file.Close()
