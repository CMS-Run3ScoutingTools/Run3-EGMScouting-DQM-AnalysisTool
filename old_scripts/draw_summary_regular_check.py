#!/usr/bin/env python3

import ROOT
import cmsstyle as CMS
import os
from array import array


CMS.SetExtraText("Preliminary")
CMS.SetEnergy("13.6")
CMS.SetLumi("Run3")

ERA = "eraF"
resonance = "resonanceZ"
plotdir = f"plot_eff_2025/{resonance}"
os.system(f"mkdir -p {plotdir}")

file_Data = ROOT.TFile.Open("/eos/user/t/tihsu/database/EGM_DQM_v4/Data.root")
file_MC = ROOT.TFile.Open("/eos/user/t/tihsu/database/EGM_DQM_TnP/DY/DY.root")

lumi_dict = {
    "eraB": 0.26,
    "eraC": 12.46,
    "eraD": 5.38,
    "eraE": 8.3,
    "eraF_preHV": 7.2,
    "eraF_postHV": 12.6,
    "eraG_digiOff": 8.1,
    "eraG_digiOn": 13.6,
    "eraF": 27.8,
    "eraG": 37.7,
}

file_F = ROOT.TFile.Open("/eos/user/t/tihsu/database/EGM_DQM_PFMonitor/DataEraF_v4/Data.root")
file_G = ROOT.TFile.Open("/eos/user/t/tihsu/database/EGM_DQM_PFMonitor/DataEraG/Data.root")
file_H = ROOT.TFile.Open("/eos/user/t/tihsu/database/EGM_DQM_PFMonitor/DataEraH_l1fixed/Data.root")


DoubleEGL1 = [
    "L1_DoubleEG_LooseIso16_LooseIso12_er1p5",
    "L1_DoubleEG_LooseIso18_LooseIso12_er1p5",
    "L1_DoubleEG_LooseIso20_LooseIso12_er1p5",
    "L1_DoubleEG_LooseIso22_LooseIso12_er1p5",
    "L1_DoubleEG11_er1p2_dR_Max0p6",
]

SinglePhotonL1_LooseIso = [
    "L1_SingleLooseIsoEG26er2p5",
    "L1_SingleLooseIsoEG26er1p5",
    "L1_SingleLooseIsoEG28er2p5",
    "L1_SingleLooseIsoEG28er2p1",
    "L1_SingleLooseIsoEG28er1p5",
    "L1_SingleLooseIsoEG30er2p5",
    "L1_SingleLooseIsoEG30er1p5",
]

SinglePhotonL1_Standard = [
    "L1_SingleEG26er2p5",
    "L1_SingleEG38er2p5",
    "L1_SingleEG40er2p5",
    "L1_SingleEG42er2p5",
    "L1_SingleEG45er2p5",
    "L1_SingleEG60",
    "L1_SingleEG34er2p5",
    "L1_SingleEG36er2p5",
]

SinglePhotonL1_Iso_er2p1 = [
    "L1_SingleIsoEG24er2p1",
    "L1_SingleIsoEG26er2p1",
    "L1_SingleIsoEG28er2p1",
    "L1_SingleIsoEG30er2p1",
    "L1_SingleIsoEG32er2p1",
]

SinglePhotonL1_Iso_er2p5 = [
    "L1_SingleIsoEG26er2p5",
    "L1_SingleIsoEG28er2p5",
    "L1_SingleIsoEG30er2p5",
    "L1_SingleIsoEG32er2p5",
    "L1_SingleIsoEG34er2p5",
]

L1_dict = {
    "double_eg": DoubleEGL1,
    "single_eg_standard": SinglePhotonL1_Standard,
    "single_eg_iso_er2p1": SinglePhotonL1_Iso_er2p1,
    "single_eg_iso_er2p5": SinglePhotonL1_Iso_er2p5,
    "single_eg_loose_iso": SinglePhotonL1_LooseIso,
}

L1Seed = DoubleEGL1 + SinglePhotonL1_LooseIso + SinglePhotonL1_Standard + SinglePhotonL1_Iso_er2p1 + SinglePhotonL1_Iso_er2p5
DST = ["DST_PFScouting_DoubleEG_v", "DST_PFScouting_SinglePhotonEB_v"]


def _get_clone(file_handle, path):
    h = file_handle.Get(path)
    if not h:
        raise RuntimeError(f"Missing histogram: {path}")
    out = h.Clone()
    out.SetDirectory(0)
    return out


def Draw(file_F_, name, era):
    lumi = lumi_dict[era]
    CMS.SetLumi("Run3 " + era + str(lumi))
    hist = _get_clone(file_F_, "DQMData/Run 333334/ScoutingMonitoring/Run summary/" + name)
    x_axis = hist.GetXaxis()
    x_binnings = [x_axis.GetBinLowEdge(bin_ + 1) for bin_ in range(x_axis.GetNbins() + 1)]

    c = CMS.cmsCanvas("", 1, max(x_binnings), 0, hist.GetMaximum() * 1.3, "", "", square=CMS.kSquare, extraSpace=0.0, iPos=0)
    CMS.cmsDraw(hist, "HIST", mcolor=ROOT.kBlack)

    if "invMass" in name:
        c.SetLogx()
    CMS.SaveCanvas(c, os.path.join(plotdir, f"{era}_{name}.png"), close=False)
    CMS.SaveCanvas(c, os.path.join(plotdir, f"{era}_{name}.pdf"))


def DrawCompareOffline(file_F_, name, era, region):
    lumi = lumi_dict[era]
    CMS.SetLumi("Run3 " + era + " " + str(lumi))
    patname = f"{resonance}_Probe_patElectron_{name}_{region}"
    sctname = f"{resonance}_Probe_sctElectron_{name}_{region}"
    hist_pat = _get_clone(file_F_, "DQMData/Run 333334/ScoutingMonitoring/Run summary/" + patname).Rebin(2)
    hist_sct = _get_clone(file_F_, "DQMData/Run 333334/ScoutingMonitoring/Run summary/" + sctname).Rebin(2)

    if hist_pat.Integral() > 0:
        hist_pat.Scale(1.0 / hist_pat.Integral())
    if hist_sct.Integral() > 0:
        hist_sct.Scale(1.0 / hist_sct.Integral())

    x_axis = hist_pat.GetXaxis()
    x_binnings = [x_axis.GetBinLowEdge(bin_ + 1) for bin_ in range(x_axis.GetNbins() + 1)]
    c = CMS.cmsDiCanvas(
        "c",
        min(x_binnings),
        max(x_binnings),
        0,
        hist_pat.GetMaximum() * 1.5,
        0.0,
        2.0,
        name,
        "norm",
        "sct/offline",
        square=CMS.kSquare,
        extraSpace=0.0,
        iPos=0,
    )
    c.cd(1)
    legend = CMS.cmsLeg(0.55, 0.6, 0.8, 0.8, textSize=0.04)
    CMS.cmsDraw(hist_pat, "HIST", lcolor=ROOT.kBlue, msize=0, lwidth=3, fstyle=0)
    CMS.cmsDraw(hist_sct, "HIST SAME", lcolor=ROOT.kRed, msize=0, lwidth=3, fstyle=0)
    legend.AddEntry(hist_pat, "Offline(Pat) Electron", "L")
    legend.AddEntry(hist_sct, "Scouting Electron", "L")

    h_ratio = hist_sct.Clone()
    h_ratio.Divide(hist_pat)
    c.cd(2)
    c.SetGridx()
    c.SetGridy()
    CMS.cmsDraw(h_ratio, "HIST", lcolor=ROOT.kRed, lwidth=3, fstyle=0)
    CMS.SaveCanvas(c, os.path.join(plotdir, f"pat_sct_compare_{name}_{region}_{era}.png"), close=False)
    CMS.SaveCanvas(c, os.path.join(plotdir, f"pat_sct_compare_{name}_{region}_{era}.pdf"))


def DrawEfficiency(fileData, fileMC, era, region, tags, type_, pt_order="leading", eta=False, name="", tagging_type="pat"):
    lumi = lumi_dict[era]
    CMS.SetLumi("Run3 " + era + " " + str(lumi))
    color_template = [ROOT.kBlue + 2, ROOT.kViolet - 2, ROOT.kRed - 7, ROOT.kOrange + 7, ROOT.kGreen - 2]

    histograms = {}
    for tag in tags:
        if tagging_type == "pat":
            if not eta:
                numerator_name = f"{resonance}_Tag_{tagging_type}_Probe_{type_}Electron_{pt_order}_Pt_{region}_pass{tag}"
                denominator_name = f"{resonance}_Tag_{tagging_type}_Probe_{type_}Electron_{pt_order}_Pt_{region}_passBaseDST"
            else:
                numerator_name = f"{resonance}_Tag_{tagging_type}_Probe_{type_}Electron_{pt_order}_Eta_pass{tag}"
                denominator_name = f"{resonance}_Tag_{tagging_type}_Probe_{type_}Electron_{pt_order}_Eta_passBaseDST"
        else:
            if not eta:
                numerator_name = f"{resonance}_{pt_order}_Pt_{region}_pass{tag}"
                denominator_name = f"{resonance}_{pt_order}_Pt_{region}_passBaseDST"
            else:
                numerator_name = f"{resonance}_{pt_order}_Eta_pass{tag}"
                denominator_name = f"{resonance}_{pt_order}_Eta_passBaseDST"

        binning = array("d", [0, 5, 10, 15, 20, 30, 40, 50, 80, 120, 200]) if not eta else array("d", [-2.5 + i * 0.5 for i in range(11)])
        targetdir = "Tag_PatElectron" if tagging_type == "pat" else "Tag_ScoutingElectron"

        hist_num_Data = _get_clone(fileData, f"DQMData/Run 333334/HLT/Run summary/ScoutingOffline/EGamma/TnP/{targetdir}/{numerator_name}").Rebin(len(binning) - 1, "num", binning)
        hist_den_Data = _get_clone(fileData, f"DQMData/Run 333334/HLT/Run summary/ScoutingOffline/EGamma/TnP/{targetdir}/{denominator_name}").Rebin(len(binning) - 1, "den", binning)
        hist_num_MC = _get_clone(fileMC, f"DQMData/Run 333334/HLT/Run summary/ScoutingOffline/EGamma/TnP/{targetdir}/{numerator_name}").Rebin(len(binning) - 1, "num", binning)
        hist_den_MC = _get_clone(fileMC, f"DQMData/Run 333334/HLT/Run summary/ScoutingOffline/EGamma/TnP/{targetdir}/{denominator_name}").Rebin(len(binning) - 1, "den", binning)

        hist_eff_Data = ROOT.TEfficiency(hist_num_Data, hist_den_Data)
        hist_eff_MC = ROOT.TEfficiency(hist_num_MC, hist_den_MC)

        eff_ratio = hist_num_Data.Clone()
        eff_ratio.Divide(hist_den_Data)
        eff_ratio_MC = hist_num_MC.Clone()
        eff_ratio_MC.Divide(hist_den_MC)
        eff_ratio.Divide(eff_ratio_MC)

        histograms[tag] = {"Data": hist_eff_Data.Clone(), "MC": hist_eff_MC.Clone(), "ratio": eff_ratio.Clone()}

    c = CMS.cmsCanvas("", 0, 200, 0, 1.2, f"{pt_order} p_{{T}} [GeV]", "efficiency", square=CMS.kSquare, extraSpace=0.0, iPos=0)
    legend = CMS.cmsLeg(0.4, 0.2, 0.8, 0.4, textSize=0.02)
    iColor = 0
    for tag in histograms:
        iColor += 1
        CMS.cmsDraw(
            histograms[tag]["Data"],
            "P E",
            lcolor=color_template[iColor % len(color_template)],
            mcolor=color_template[iColor % len(color_template)],
            msize=1,
            lwidth=2,
            fstyle=0,
        )
        legend.AddEntry(histograms[tag]["Data"], f"{tag}[Data]", "PL")

    CMS.SaveCanvas(c, os.path.join(plotdir, f"{name}_tagging_{tagging_type}_{pt_order}_{type_}Electron_Efficiency_{region}_comparison.png"), close=False)
    CMS.SaveCanvas(c, os.path.join(plotdir, f"{name}_tagging_{tagging_type}_{pt_order}_{type_}Electron_Efficiency_{region}_comparison.pdf"))


def DrawMass_perEra(file_dict):
    CMS.SetLumi("Run3 ")
    histos = {}
    ymax = 0
    for era, file_ in file_dict.items():
        lumi = lumi_dict[era]
        name = "sctElectron_EBEB_appliedID_invMass_pass_DST_PFScouting_DoubleEG_v"
        hist = _get_clone(file_, "DQMData/Run 333334/HLT/Run summary/ScoutingOffline/EGamma/Collection/" + name)
        hist.Scale(1 / lumi)
        hist.GetXaxis().SetRangeUser(0, 130)
        hist.Rebin(2)
        histos[era] = hist.Clone()
        ymax = max(ymax, hist.GetMaximum() * 2)

    c = CMS.cmsCanvas("", 0.5, 130, 1, ymax, "inv mass", "nEntries / Lumi", square=CMS.kSquare, extraSpace=0.0, iPos=0)
    legend = CMS.cmsLeg(0.55, 0.7, 0.85, 0.9, textSize=0.02)
    for icolor, era_ in enumerate(histos):
        CMS.cmsDraw(histos[era_], "HIST", lcolor=ROOT.kBlue + icolor * 2, msize=0, lwidth=3, lstyle=icolor, fstyle=0)
        legend.AddEntry(histos[era_], f"Run2024 ({era_})", "L")
    c.SetLogy()
    os.makedirs("eraPlot", exist_ok=True)
    CMS.SaveCanvas(c, os.path.join("eraPlot", "compare_mass.png"), close=False)
    CMS.SaveCanvas(c, os.path.join("eraPlot", "compare_mass.pdf"))


def DrawOffline_Comparison(fileData, var, region, era):
    lumi = lumi_dict[era]
    CMS.SetLumi("Run3 " + era + " " + str(lumi))
    hist_Scouting = _get_clone(fileData, "DQMData/Run 333334/HLT/Run summary/ScoutingOffline/EGamma/TnP/Tag_PatElectron/" + f"resonanceZ_Tag_pat_Probe_patElectron_{var}_{region}")
    hist_Offline = _get_clone(fileData, "DQMData/Run 333334/HLT/Run summary/ScoutingOffline/EGamma/TnP/Tag_PatElectron/" + f"resonanceZ_Tag_pat_Probe_sctElectron_{var}_{region}")

    hist_Scouting.Rebin(5)
    hist_Offline.Rebin(5)
    if hist_Scouting.Integral() > 0:
        hist_Scouting.Scale(1.0 / hist_Scouting.Integral())
    if hist_Offline.Integral() > 0:
        hist_Offline.Scale(1.0 / hist_Offline.Integral())

    c = CMS.cmsDiCanvas("", 0, 200, 0, max(hist_Scouting.GetMaximum(), hist_Offline.GetMaximum()) * 1.2, 0.0, 2.0, var, "A.U.", "data/pred", square=CMS.kSquare, extraSpace=0.0, iPos=0)
    legend = CMS.cmsLeg(0.6, 0.6, 0.8, 0.8, textSize=0.04)
    c.SetLogy()
    c.cd(1)
    CMS.cmsDraw(hist_Offline, "HIST", lcolor=ROOT.kOrange, mcolor=ROOT.kOrange, msize=1, lwidth=2, fstyle=1001, fcolor=ROOT.kOrange)
    CMS.cmsDraw(hist_Scouting, "P E", lcolor=ROOT.kBlack, mcolor=ROOT.kBlack, msize=1, lwidth=2, fstyle=0)
    legend.AddEntry(hist_Offline, "Offline Electron", "F")
    legend.AddEntry(hist_Scouting, "Scouting Electron", "PL")
    c.cd(2)
    hist_ratio = hist_Scouting.Clone()
    hist_ratio.Divide(hist_Offline)
    CMS.cmsDraw(hist_ratio, "P E", lcolor=ROOT.kBlack, mcolor=ROOT.kBlack, msize=1, lwidth=2, fstyle=0)
    os.makedirs("offline_comparison_plot", exist_ok=True)
    CMS.SaveCanvas(c, os.path.join("offline_comparison_plot", f"Electron_{var}_{region}.png"), close=False)


def DrawMC_Comparison(fileData, fileMC_, type_, var, region, era):
    lumi = lumi_dict[era]
    CMS.SetLumi("Run3 " + era + " " + str(lumi))
    hist_Data = _get_clone(fileData, "DQMData/Run 333334/HLT/Run summary/ScoutingOffline/EGamma/TnP/Tag_PatElectron/" + f"resonanceZ_Tag_pat_Probe_{type_}Electron_{var}_{region}")
    hist_MC = _get_clone(fileMC_, "DQMData/Run 333334/HLT/Run summary/ScoutingOffline/EGamma/TnP/Tag_PatElectron/" + f"resonanceZ_Tag_pat_Probe_{type_}Electron_{var}_{region}")
    total_entry = fileMC_.Get("DQMData/Run 333334/HLT/Run summary/ScoutingOffline/EGamma/TnP/Tag_PatElectron/total_entry").GetEntries()

    hist_Data.Rebin(5)
    hist_MC.Rebin(5)
    hist_MC.Scale(lumi * 1000 / total_entry * 6077.22)
    if hist_Data.Integral() > 0:
        hist_Data.Scale(1.0 / hist_Data.Integral())
    if hist_MC.Integral() > 0:
        hist_MC.Scale(1.0 / hist_MC.Integral())

    c = CMS.cmsDiCanvas("", 0, 200, 0, max(hist_Data.GetMaximum(), hist_MC.GetMaximum()) * 1.2, 0.0, 2.0, var, "A.U.", "data/pred", square=CMS.kSquare, extraSpace=0.0, iPos=0)
    legend = CMS.cmsLeg(0.6, 0.6, 0.8, 0.8, textSize=0.04)
    c.SetLogy()
    c.cd(1)
    CMS.cmsDraw(hist_MC, "HIST", lcolor=ROOT.kOrange, mcolor=ROOT.kOrange, msize=1, lwidth=2, fstyle=1001, fcolor=ROOT.kOrange)
    CMS.cmsDraw(hist_Data, "P E", lcolor=ROOT.kBlack, mcolor=ROOT.kBlack, msize=1, lwidth=2, fstyle=0)
    legend.AddEntry(hist_MC, "Z#rightarrow ee Summer 24 MC", "F")
    legend.AddEntry(hist_Data, f"PFMonitoring {era}", "PL")
    c.cd(2)
    hist_ratio = hist_Data.Clone()
    hist_ratio.Divide(hist_MC)
    CMS.cmsDraw(hist_ratio, "P E", lcolor=ROOT.kBlack, mcolor=ROOT.kBlack, msize=1, lwidth=2, fstyle=0)
    os.makedirs("mc_plot", exist_ok=True)
    CMS.SaveCanvas(c, os.path.join("mc_plot", f"{type_}Electron_{var}_{region}.png"), close=False)
    CMS.SaveCanvas(c, os.path.join("mc_plot", f"{type_}Electron_{var}_{region}.pdf"))


for variable in ["Pt", "HoverE", "OoEMOoP", "MissingHits", "RelEcalIsolation", "RelHcalIsolation", "RelTrackIsolation", "SigmaIetaIeta", "Trackfbrem", "dEtaIn", "dPhiIn"]:
    for type_ in ["pat", "sct"]:
        for region in ["Barrel", "Endcap"]:
            DrawMC_Comparison(file_Data, file_MC, type_, variable, region, ERA)
            DrawOffline_Comparison(file_Data, variable, region, ERA)

DrawMass_perEra({"eraF": file_F, "eraG": file_G})

file_Dict = {
    "eraC": ROOT.TFile.Open("/eos/user/t/tihsu/database/DPNote/ScoutingPFMonitorData2025/eraC/data.root"),
    "eraD": ROOT.TFile.Open("/eos/user/t/tihsu/database/DPNote/ScoutingPFMonitorData2025/eraD/data.root"),
    "eraE": ROOT.TFile.Open("/eos/user/t/tihsu/database/EGM_DQM_PFMonitor/eraE/data.root"),
    "eraF_preHV": ROOT.TFile.Open("/eos/user/t/tihsu/database/EGM_DQM_PFMonitor/eraF_preHV/data.root"),
    "eraF_postHV": ROOT.TFile.Open("/eos/user/t/tihsu/database/EGM_DQM_PFMonitor/eraF_postHV/data.root"),
    "eraG_digiOff": ROOT.TFile.Open("/eos/user/t/tihsu/database/EGM_DQM_PFMonitor/2025_eraG_before_run_398288_v2/data.root"),
    "eraG_digiOn": ROOT.TFile.Open("/eos/user/t/tihsu/database/EGM_DQM_PFMonitor/2025_eraG_after_run_398288_v2/data.root"),
}

for pt_order in ["leading", "subleading"]:
    for DST_ in DST:
        fireObj = "_fireTrigObj"
        DrawEfficiency(file_Data, file_Data, ERA, "Barrel", [DST_ + fireObj], "sct", pt_order, name=DST_ + fireObj, tagging_type="sct")
        DrawEfficiency(file_Data, file_Data, ERA, "Endcap", [DST_ + fireObj], "sct", pt_order, name=DST_ + fireObj, tagging_type="sct")
        DrawEfficiency(file_Data, file_Data, ERA, "Full", [DST_ + fireObj], "sct", pt_order, eta=True, name=DST_ + fireObj, tagging_type="sct")

    for name in L1_dict:
        DrawEfficiency(file_Data, file_MC, ERA, "Barrel", L1_dict[name], "pat", pt_order, name=name)
        DrawEfficiency(file_Data, file_MC, ERA, "Barrel", L1_dict[name], "sct", pt_order, name=name)
        DrawEfficiency(file_Data, file_MC, ERA, "Endcap", L1_dict[name], "pat", pt_order, name=name)
        DrawEfficiency(file_Data, file_MC, ERA, "Endcap", L1_dict[name], "sct", pt_order, name=name)
