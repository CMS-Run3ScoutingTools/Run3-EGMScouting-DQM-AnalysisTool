import ROOT
import cmsstyle as CMS
import os, sys
from array import array

CMS.SetExtraText("Preliminary")
CMS.SetEnergy("13.6")
CMS.SetLumi("Run3")

ERA = 'eraF'
resonance = 'resonanceZ'
plotdir = f"plot_eff_2025/{resonance}"
os.system("mkdir -p {}".format(plotdir))

# Open the ROOT file
file_Data = ROOT.TFile.Open("/eos/user/t/tihsu/database/EGM_DQM_v4/Data.root")
file_MC = ROOT.TFile.Open("/eos/user/t/tihsu/database/EGM_DQM_TnP/DY/DY.root")
lumi_dict = {'eraF': 27.8, 'eraG': 37.7, 'eraH': 5.4, 'eraI': 1.0}

lumi_dict = {
  'eraB': 0.26,
  'eraC': 12.46,
  'eraD': 5.38,
  'eraE': 8.3,
  'eraF_preHV': 7.2,
  'eraF_postHV': 12.6,
  'eraG_digiOff': 8.1,
  'eraG_digiOn': 13.6,
'eraF': 27.8, 'eraG': 37.7,
}

file_F = ROOT.TFile.Open("/eos/user/t/tihsu/database/EGM_DQM_PFMonitor/DataEraF_v4/Data.root")
file_G = ROOT.TFile.Open("/eos/user/t/tihsu/database/EGM_DQM_PFMonitor/DataEraG/Data.root")
file_H = ROOT.TFile.Open("/eos/user/t/tihsu/database/EGM_DQM_PFMonitor/DataEraH_l1fixed//Data.root")


# L1Seed
DoubleEGL1 = [
  "L1_DoubleEG_LooseIso16_LooseIso12_er1p5",
  "L1_DoubleEG_LooseIso18_LooseIso12_er1p5",
  "L1_DoubleEG_LooseIso20_LooseIso12_er1p5",
  "L1_DoubleEG_LooseIso22_LooseIso12_er1p5",
  "L1_DoubleEG11_er1p2_dR_Max0p6"
]

SinglePhotonL1 = [
  'L1_SingleLooseIsoEG26er2p5',
  'L1_SingleLooseIsoEG26er1p5',
  'L1_SingleLooseIsoEG28er2p5',
  'L1_SingleLooseIsoEG28er2p1',
  'L1_SingleLooseIsoEG28er1p5',
  'L1_SingleLooseIsoEG30er2p5',
  'L1_SingleLooseIsoEG30er1p5',
  'L1_SingleEG26er2p5',
  'L1_SingleEG38er2p5',
  'L1_SingleEG40er2p5',
  'L1_SingleEG42er2p5',
  'L1_SingleEG45er2p5',
  'L1_SingleEG60',
  'L1_SingleEG34er2p5',
  'L1_SingleEG36er2p5',
  'L1_SingleIsoEG24er2p1',
  'L1_SingleIsoEG26er2p1',
  'L1_SingleIsoEG28er2p1',
  'L1_SingleIsoEG30er2p1',
  'L1_SingleIsoEG32er2p1',
  'L1_SingleIsoEG26er2p5',
  'L1_SingleIsoEG28er2p5',
  'L1_SingleIsoEG30er2p5',
  'L1_SingleIsoEG32er2p5',
  'L1_SingleIsoEG34er2p5'
]

DoubleEGL1 = [
  "L1_DoubleEG_LooseIso16_LooseIso12_er1p5",
  "L1_DoubleEG_LooseIso18_LooseIso12_er1p5",
  "L1_DoubleEG_LooseIso20_LooseIso12_er1p5",
  "L1_DoubleEG_LooseIso22_LooseIso12_er1p5",
  "L1_DoubleEG11_er1p2_dR_Max0p6"
]

SinglePhotonL1_LooseIso = [
  'L1_SingleLooseIsoEG26er2p5',
  'L1_SingleLooseIsoEG26er1p5',
  'L1_SingleLooseIsoEG28er2p5',
  'L1_SingleLooseIsoEG28er2p1',
  'L1_SingleLooseIsoEG28er1p5',
  'L1_SingleLooseIsoEG30er2p5',
  'L1_SingleLooseIsoEG30er1p5'
]

SinglePhotonL1_Standard = [
  'L1_SingleEG26er2p5',
  'L1_SingleEG38er2p5',
  'L1_SingleEG40er2p5',
  'L1_SingleEG42er2p5',
  'L1_SingleEG45er2p5',
  'L1_SingleEG60',
  'L1_SingleEG34er2p5',
  'L1_SingleEG36er2p5'
]
SinglePhotonL1_Iso_er2p1 = [
  'L1_SingleIsoEG24er2p1',
  'L1_SingleIsoEG26er2p1',
  'L1_SingleIsoEG28er2p1',
  'L1_SingleIsoEG30er2p1',
  'L1_SingleIsoEG32er2p1',
]

SinglePhotonL1_Iso_er2p5 = [
  'L1_SingleIsoEG26er2p5',
  'L1_SingleIsoEG28er2p5',
  'L1_SingleIsoEG30er2p5',
  'L1_SingleIsoEG32er2p5',
  'L1_SingleIsoEG34er2p5'
]

L1_dict = {
    'double_eg': DoubleEGL1,
    'single_eg_standard': SinglePhotonL1_Standard,
    'single_eg_iso_er2p1': SinglePhotonL1_Iso_er2p1,
    'single_eg_iso_er2p5': SinglePhotonL1_Iso_er2p5,
    'single_eg_loose_iso': SinglePhotonL1_LooseIso
}

L1Seed = DoubleEGL1 + SinglePhotonL1
DST = ["DST_PFScouting_DoubleEG_v", "DST_PFScouting_SinglePhotonEB_v"]


def Draw(file_F, name, era):
    lumi = lumi_dict[era]
    CMS.SetLumi("Run3 " + era + str(lumi))
    hist = file_F.Get("DQMData/Run 333334/ScoutingMonitoring/Run summary/" + name).Clone()
    x_axis = hist.GetXaxis()
    y_axis = hist.GetYaxis()
    nbinX  = x_axis.GetNbins()
    nbinY  = y_axis.GetNbins()
    x_binnings = [x_axis.GetBinLowEdge(bin_+1) for bin_ in range(nbinX+1)]
    y_binnings = [y_axis.GetBinLowEdge(bin_+1) for bin_ in range(nbinY+1)]

    c = CMS.cmsCanvas('', 1, max(x_binnings), 0, hist.GetMaximum()*1.3, '', '', square = CMS.kSquare, extraSpace=0.0, iPos=0)
    legend = CMS.cmsLeg(0.65, 0.2, 0.85, 0.4, textSize=0.04)

    CMS.cmsDraw(hist, 'HIST', mcolor = ROOT.kBlack)
    

    if 'invMass' in name: 
      c.SetLogx()
    CMS.SaveCanvas(c, os.path.join(plotdir, '{}_{}.png'.format(era,name)), close = False)
    CMS.SaveCanvas(c, os.path.join(plotdir, '{}_{}.pdf'.format(era,name)))



def Draw_2D(file_F, name, era):
    lumi = lumi_dict[era]
    CMS.SetLumi("Run3 " + era + " " + str(lumi))
    # Get the 2D histogram from the file
    hist = file_F.Get("DQMData/Run 333334/ScoutingMonitoring/Run summary/" + name).Clone()
    profile_x = hist.ProfileX()
 
    # Check if the histogram exists
    if not hist:
        print(f"Histogram {name} not found in file.")
        return
    
    # Set up axes and binning for 2D
    x_axis = hist.GetXaxis()
    y_axis = hist.GetYaxis()
    nbinX  = x_axis.GetNbins()
    nbinY  = y_axis.GetNbins()
    x_binnings = [x_axis.GetBinLowEdge(bin_+1) for bin_ in range(nbinX+1)]
    y_binnings = [y_axis.GetBinLowEdge(bin_+1) for bin_ in range(nbinY+1)]

    # Create the canvas for 2D plot
    c = CMS.cmsCanvas('', 1, max(x_binnings), 0, hist.GetMaximum()*1.3, '', '', square=CMS.kSquare, extraSpace=0.0, iPos=0)
    
    # Apply color palette for better 2D visualization
    ROOT.gStyle.SetPalette(ROOT.kBird)
    hist.SetContour(100)  # for smooth gradient
    
    # Draw the 2D histogram with color gradient
    hist.Draw("COLZ")

    # Initialize arrays to hold graph points
    n_bins = hist.GetNbinsX()
    x_vals = []
    x_errs = []
    y_vals = []
    y_errs = []

    # Loop through x-bins to collect mean and width
    for i in range(1, n_bins + 1):
        proj_y = hist.ProjectionY("_py", i, i)  # Projection for bin i
        if proj_y.GetEntries() > 0:
            x_vals.append(hist.GetXaxis().GetBinCenter(i))  # Bin center
            x_errs.append(0.0)  # No x-errors
            y_vals.append(proj_y.GetMean())  # Mean of y
            y_errs.append(proj_y.GetStdDev())  # Std deviation of y
        else:
            x_vals.append(hist.GetXaxis().GetBinCenter(i))
            x_errs.append(0.0)
            y_vals.append(0.0)
            y_errs.append(0.0)

    # Create a TGraphErrors for the profile with custom uncertainties
    graph = ROOT.TGraphErrors(len(x_vals), 
                              array('d', x_vals),
                              array('d', y_vals),
                              array('d', x_errs),
                              array('d', y_errs))
    
    graph.Draw("SAME")


    # Set log scale for x-axis if 'invMass' in name
    if 'invMass' in name:
        c.SetLogx()


    # Save canvas as PNG and PDF
    CMS.SaveCanvas(c, os.path.join(plotdir, f"{era}_{name}.png"), close=False)
    CMS.SaveCanvas(c, os.path.join(plotdir, f"{era}_{name}.pdf"))


def DrawCompareOffline(file_F, name, era, region):
    lumi = lumi_dict[era]
    CMS.SetLumi("Run3 " + era + " " + str(lumi))
    patname = "{}_Probe_patElectron_{}_{}".format(resonance, name, region)
    sctname = "{}_Probe_sctElectron_{}_{}".format(resonance, name, region)
    #patname = "{}_Tag_pat_Probe_patElectron_{}_{}".format(resonance, name, region)
    #sctname = "{}_Tag_pat_Probe_sctElectron_{}_{}".format(resonance, name, region)
    hist_pat = file_F.Get("DQMData/Run 333334/ScoutingMonitoring/Run summary/" + patname).Clone().Rebin(2)
    hist_sct = file_F.Get("DQMData/Run 333334/ScoutingMonitoring/Run summary/" + sctname).Clone().Rebin(2)
    hist_pat.Scale(1./hist_pat.Integral())
    hist_sct.Scale(1./hist_sct.Integral())
    hist_pat.GetXaxis().SetTitle(name)
    hist_pat.GetYaxis().SetTitle("norm.")

    x_axis = hist_pat.GetXaxis()
    y_axis = hist_pat.GetYaxis()
    nbinX  = x_axis.GetNbins()
    nbinY  = y_axis.GetNbins()
    x_binnings = [x_axis.GetBinLowEdge(bin_+1) for bin_ in range(nbinX+1)]
    y_binnings = [y_axis.GetBinLowEdge(bin_+1) for bin_ in range(nbinY+1)]

    c = CMS.cmsDiCanvas('c', min(x_binnings), max(x_binnings),  0, hist_pat.GetMaximum()*1.5, 0.0, 2.0, name, 'norm', 'sct/offline', square = CMS.kSquare, extraSpace=0.0, iPos=0)
    c.cd(1)
    legend = CMS.cmsLeg(0.55, 0.6, 0.8, 0.8, textSize=0.04)

    CMS.cmsDraw(hist_pat, 'HIST ', lcolor = ROOT.kBlue, msize=0, lwidth = 3, fstyle = 0)
    CMS.cmsDraw(hist_sct, 'HIST SAME', lcolor = ROOT.kRed, msize=0, lwidth = 3, fstyle = 0)

    legend.AddEntry(hist_pat, 'Offline(Pat) Electron', 'L')
    legend.AddEntry(hist_sct, 'Scouting Electron', 'L')

    h_ratio = hist_sct.Clone()
    h_ratio.Divide(hist_pat)
    c.cd(2)
    c.SetGridx()
    c.SetGridy()
    CMS.cmsDraw(h_ratio, 'HIST', lcolor = ROOT.kRed, lwidth = 3, fstyle=0)

    #c.SetLogx()
    CMS.SaveCanvas(c, os.path.join(plotdir, 'pat_sct_compare_{}_{}_{}.png'.format(name, region, era)), close = False)
    CMS.SaveCanvas(c, os.path.join(plotdir, 'pat_sct_compare_{}_{}_{}.pdf'.format(name, region, era)))

def DrawEfficiency(fileData, fileMC, era, region, tags, type_, pt_order = "leading", eta=False, name = "", tagging_type = "pat"):
    lumi = lumi_dict[era]
    CMS.SetLumi("Run3 " + era + " " + str(lumi))

    color_template = [ROOT.kBlue+2, ROOT.kViolet-2, ROOT.kRed-7, ROOT.kOrange+7, ROOT.kGreen - 2, ROOT.kYellow+2, ROOT.kCyan+1, ROOT.kSpring -7, ROOT.kPink -1]

    Histograms = dict()
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
      if not eta:
           binning = array('d', [0, 5, 10, 15, 20, 30, 40, 50, 80, 120, 200])
      
      else:
         binning = array('d', [-2.5 + i*0.5 for i in range(11)])

      if tagging_type == "pat":
        targetdir = "Tag_PatElectron"
      else:
        targetdir = "Tag_ScoutingElectron"


      print(f"✅ {name} numerator ⚠️  num: DQMData/Run 333334/HLT/Run summary/ScoutingOffline/EGamma/TnP/{targetdir}/{numerator_name}  ⚠️  den: DQMData/Run 333334/HLT/Run summary/ScoutingOffline/EGamma/TnP/{targetdir}/{denominator_name}")
      print(f"DQMData/Run 333334/HLT/Run summary/ScoutingOffline/EGamma/TnP/{targetdir}/" + numerator_name)
      print(f"DQMData/Run 333334/HLT/Run summary/ScoutingOffline/EGamma/TnP/{targetdir}/" + denominator_name)
      hist_num_Data = fileData.Get(f"DQMData/Run 333334/HLT/Run summary/ScoutingOffline/EGamma/TnP/{targetdir}/" + numerator_name).Clone().Rebin(len(binning)-1, "num", binning)
      hist_den_Data = fileData.Get(f"DQMData/Run 333334/HLT/Run summary/ScoutingOffline/EGamma/TnP/{targetdir}/" + denominator_name).Clone().Rebin(len(binning)-1, "den", binning)
      hist_num_MC = fileMC.Get(f"DQMData/Run 333334/HLT/Run summary/ScoutingOffline/EGamma/TnP/{targetdir}/" + numerator_name).Clone().Rebin(len(binning)-1, "num", binning)
      hist_den_MC = fileMC.Get(f"DQMData/Run 333334/HLT/Run summary/ScoutingOffline/EGamma/TnP/{targetdir}/" + denominator_name).Clone().Rebin(len(binning)-1, "den", binning)

      hist_num_Data.SetDirectory(0)
      hist_den_Data.SetDirectory(0)
      hist_num_MC.SetDirectory(0)
      hist_den_MC.SetDirectory(0)

      hist_eff_Data   = ROOT.TEfficiency(hist_num_Data, hist_den_Data)
      hist_eff_MC   = ROOT.TEfficiency(hist_num_MC, hist_den_MC)

      

      eff_ratio = hist_num_Data.Clone()
      eff_ratio.Divide(hist_den_Data)
      eff_ratio_MC = hist_num_MC.Clone()
      eff_ratio_MC.Divide(hist_den_MC)
      eff_ratio.Divide(eff_ratio_MC)

      print(hist_num_MC.Integral())

      x_axis = hist_num_Data.GetXaxis()
      y_axis = hist_num_Data.GetYaxis()
      nbinX  = x_axis.GetNbins()
      nbinY  = y_axis.GetNbins()
      x_binnings = [x_axis.GetBinLowEdge(bin_+1) for bin_ in range(nbinX+1)]
      y_binnings = [y_axis.GetBinLowEdge(bin_+1) for bin_ in range(nbinY+1)]

      Histograms[tag] = {
        "Data": hist_eff_Data.Clone(),
        "MC": hist_eff_MC.Clone(),
        "ratio": eff_ratio.Clone()
      }

    #c = CMS.cmsDiCanvas('', min(x_binnings), max(x_binnings), 0, max(y_binnings) * 1.2, 0.0, 2.0, f"{pt_order} "+"p_{T} [GeV]", 'efficiency', 'data/pred', square = CMS.kSquare, extraSpace=0.0, iPos=0)
    c = CMS.cmsCanvas('', min(x_binnings), max(x_binnings), 0, max(y_binnings) * 1.2, f"{pt_order} "+"p_{T} [GeV]", 'efficiency', square = CMS.kSquare, extraSpace=0.0, iPos=0)
    legend = CMS.cmsLeg(0.4, 0.2, 0.8, 0.4, textSize=0.02)
 
    iColor = 0
    for tag in Histograms:
#      c.cd(1)
      iColor += 1
      CMS.cmsDraw(Histograms[tag]["Data"], 'P E', lcolor = color_template[iColor], mcolor = color_template[iColor], msize = 1, lwidth = 2, fstyle = 0)
      legend.AddEntry(Histograms[tag]["Data"], f"{tag}[Data]", 'PL')

      if len(Histograms) > 1:
        continue
#      CMS.cmsDraw(Histograms[tag]["MC"], 'P E', lcolor = ROOT.kRed + iColor, mcolor = ROOT.kRed + iColor, msize = 1, lwidth = 2, fstyle = 0, marker = 4)
#      legend.AddEntry(Histograms[tag]["MC"], f"{tag}[MC]", 'PL')

#      c.cd(2)
#      CMS.cmsDraw(Histograms[tag]["ratio"], 'P', lcolor = iColor, mcolor = iColor, msize = 1, lwidth = 0, fstyle = 0)

    CMS.SaveCanvas(c, os.path.join(plotdir, f'{name}_tagging_{tagging_type}_{pt_order}_{type_}Electron_Efficiency_{region}_comparison.png'), close = False)
    CMS.SaveCanvas(c, os.path.join(plotdir, f'{name}_tagging_{tagging_type}_{pt_order}_{type_}Electron_Efficiency_{region}_comparison.pdf'))

def DrawMass_perEra(file_dict):

    CMS.SetLumi("Run3 ")
    histos = dict()
    ymax = 0
    for era, file in file_dict.items():
      lumi = lumi_dict[era]
      name = f"sctElectron_EBEB_appliedID_invMass_pass_DST_PFScouting_DoubleEG_v"
      hist = file.Get("DQMData/Run 333334/HLT/Run summary/ScoutingOffline/EGamma/Collection/" + name).Clone()
      hist.Scale(1/lumi)

      hist.GetXaxis().SetRangeUser(0, 130)

      x_axis = hist.GetXaxis()
      y_axis = hist.GetYaxis()
      nbinX  = x_axis.GetNbins()
      nbinY  = y_axis.GetNbins()
      x_binnings = [x_axis.GetBinLowEdge(bin_+1) for bin_ in range(nbinX+1)]
      y_binnings = [y_axis.GetBinLowEdge(bin_+1) for bin_ in range(nbinY+1)]

      hist.Rebin(2)
      histos[era] = hist.Clone()
      if ymax < hist.GetMaximum()*2:
           ymax = hist.GetMaximum()*2

    c = CMS.cmsCanvas('', 0.5, max(x_binnings), 1, ymax, 'inv mass', 'nEntries / Lumi', square = CMS.kSquare, extraSpace=0.0, iPos=0)
    legend = CMS.cmsLeg(0.55, 0.7, 0.85, 0.9, textSize=0.02)


    for icolor, era_ in enumerate(histos):
        CMS.cmsDraw(histos[era_], 'HIST ', lcolor = ROOT.kBlue + icolor * 2, msize=0, lwidth = 3, lstyle = icolor, fstyle = 0)
        legend.AddEntry(histos[era_], f"Run2024 ({era_})", 'L')

    c.SetLogy()

    CMS.SaveCanvas(c, os.path.join("eraPlot", f'compare_mass.png'), close = False)
    CMS.SaveCanvas(c, os.path.join("eraPlot", f'compare_mass.pdf'))


def DrawEfficiency_perEra(file_dict, era, region, tags, type_, pt_order = "leading", eta=False, name_ = "", tagging_type = "pat"):
    lumi = lumi_dict[era]
    CMS.SetLumi("Run3 2025")

    color_template = [
    ROOT.TColor.GetColor("#0072B2"),  # Blue
    ROOT.TColor.GetColor("#D55E00"),  # Vermillion
    ROOT.TColor.GetColor("#009E73"),  # Bluish green
    ROOT.TColor.GetColor("#CC79A7"),  # Reddish purple
    ROOT.TColor.GetColor("#F0E442"),  # Yellow
    ROOT.TColor.GetColor("#56B4E9"),  # Sky blue
    ROOT.TColor.GetColor("#E69F00"),  # Orange
    ROOT.TColor.GetColor("#999999"),  # Gray
    ]
    Histograms = dict()
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
      if not eta:
           binning = array('d', [0, 5, 10, 15, 20, 30, 40, 50, 80, 120, 200]) if "Single" in numerator_name else  array('d', [0, 10, 20, 30, 40, 50, 120])
      
      else:
         binning = array('d', [-2.5 + i*0.5 for i in range(11)])

      if tagging_type == "pat":
        targetdir = "Tag_PatElectron"
      else:
        targetdir = "Tag_ScoutingElectron"

      Histograms = dict()
      for name, file_ in file_dict.items():
        print(f"DQMData/Run 333334/HLT/Run summary/ScoutingOffline/EGamma/TnP/{targetdir}/" + numerator_name)
        print(f"DQMData/Run 333334/HLT/Run summary/ScoutingOffline/EGamma/TnP/{targetdir}/" + denominator_name)
        hist_num_Data = file_.Get(f"DQMData/Run 333334/HLT/Run summary/ScoutingOffline/EGamma/TnP/{targetdir}/" + numerator_name).Clone().Rebin(len(binning)-1, "num", binning)
        hist_den_Data = file_.Get(f"DQMData/Run 333334/HLT/Run summary/ScoutingOffline/EGamma/TnP/{targetdir}/" + denominator_name).Clone().Rebin(len(binning)-1, "den", binning)

        hist_num_Data.SetDirectory(0)
        hist_den_Data.SetDirectory(0)

        hist_eff_Data   = ROOT.TEfficiency(hist_num_Data, hist_den_Data)



        x_axis = hist_num_Data.GetXaxis()
        y_axis = hist_num_Data.GetYaxis()
        nbinX  = x_axis.GetNbins()
        nbinY  = y_axis.GetNbins()
        x_binnings = [x_axis.GetBinLowEdge(bin_+1) for bin_ in range(nbinX+1)]
        y_binnings = [y_axis.GetBinLowEdge(bin_+1) for bin_ in range(nbinY+1)]

        Histograms[name] = hist_eff_Data


    xtitle = f"{pt_order} "+"p_{T} [GeV]" if not region == "Full" else "#eta"
    c = CMS.cmsCanvas('', min(x_binnings), max(x_binnings), 0, max(y_binnings) * 1.5, xtitle, 'efficiency', square = CMS.kSquare, extraSpace=0.0, iPos=0)
    legend = CMS.cmsLeg(0.15, 0.7, 0.9, 0.9, textSize=0.02)
    legend.SetNColumns(3)

 
    iColor = 0
    for name in Histograms:
      c.cd(1)
      CMS.cmsDraw(Histograms[name], 'P E', lcolor = color_template[iColor], mcolor = color_template[iColor], msize = 1, lwidth = 2, fstyle = 0)
      iColor += 1
    
      legend.AddEntry(Histograms[name], f"{name}", 'PL')

    legend.SetHeader("#bf{%s}"%name_)

    plotdir = "plot_eff_perEra"
    os.makedirs(plotdir, exist_ok = True)
    CMS.SaveCanvas(c, os.path.join(plotdir, f'{name_}_tagging_{tagging_type}_{pt_order}_{type_}Electron_Efficiency_{region}_{era}.png'), close = False)
    CMS.SaveCanvas(c, os.path.join(plotdir, f'{name_}_tagging_{tagging_type}_{pt_order}_{type_}Electron_Efficiency_{region}_{era}.pdf'))




def DrawOffline_Comparison(fileData, var, region, era):

    lumi = lumi_dict[era]
    CMS.SetLumi("Run3 " + era + " " + str(lumi))

    Histograms = dict()
    hist_Scouting = fileData.Get("DQMData/Run 333334/HLT/Run summary/ScoutingOffline/EGamma/TnP/Tag_PatElectron/" + f"resonanceZ_Tag_pat_Probe_patElectron_{var}_{region}").Clone()
    hist_Offline = fileData.Get("DQMData/Run 333334/HLT/Run summary/ScoutingOffline/EGamma/TnP/Tag_PatElectron/" +  f"resonanceZ_Tag_pat_Probe_sctElectron_{var}_{region}").Clone()

    hist_Scouting.Rebin(5)
    hist_Offline.Rebin(5)

    hist_Scouting.Scale(1./hist_Scouting.Integral())
    hist_Offline.Scale(1./hist_Offline.Integral())


    x_axis = hist_Scouting.GetXaxis()
    y_axis = hist_Scouting.GetYaxis()
    nbinX  = x_axis.GetNbins()
    nbinY  = y_axis.GetNbins()
    x_binnings = [x_axis.GetBinLowEdge(bin_+1) for bin_ in range(nbinX+1)]
    y_binnings = [y_axis.GetBinLowEdge(bin_+1) for bin_ in range(nbinY+1)]

    c = CMS.cmsDiCanvas('', min(x_binnings), max(x_binnings), 0, max(hist_Scouting.GetMaximum(), hist_Offline.GetMaximum()) * 1.2, 0.0, 2.0, var, 'A.U.', 'data/pred', square = CMS.kSquare, extraSpace=0.0, iPos=0)
    legend = CMS.cmsLeg(0.6, 0.6, 0.8, 0.8, textSize=0.04)

    c.SetLogy()

    c.cd(1)
    CMS.cmsDraw(hist_Offline,   'HIST', lcolor = ROOT.kOrange, mcolor = ROOT.kOrange, msize = 1, lwidth = 2, fstyle = 1001, fcolor = ROOT.kOrange)
    CMS.cmsDraw(hist_Scouting, 'P E', lcolor = ROOT.kBlack, mcolor = ROOT.kBlack, msize = 1, lwidth = 2, fstyle = 0)

    legend.AddEntry(hist_Offline, f"Offline Electron", 'F')
    legend.AddEntry(hist_Scouting, f"Scouting Electron", 'PL')


    c.cd(2)
    hist_ratio = hist_Scouting.Clone()
    hist_ratio.Divide(hist_Offline)
    CMS.cmsDraw(hist_ratio, 'P E', lcolor = ROOT.kBlack, mcolor = ROOT.kBlack, msize = 1, lwidth = 2, fstyle = 0)

    plotdir = "offline_comparison_plot"
    os.makedirs(plotdir, exist_ok = True)
    CMS.SaveCanvas(c, os.path.join(plotdir, f'Electron_{var}_{region}.png'), close = False)




def DrawMC_Comparison(fileData, fileMC, type_, var, region, era):

    lumi = lumi_dict[era]
    CMS.SetLumi("Run3 " + era + " " + str(lumi))

    Histograms = dict()
    print("DQMData/Run 333334/HLT/Run summary/ScoutingOffline/EGamma/TnP/Tag_PatElectron/" + f"resonanceZ_Tag_pat_Probe_{type_}Electron_{var}_{region}")
    hist_Data = fileData.Get("DQMData/Run 333334/HLT/Run summary/ScoutingOffline/EGamma/TnP/Tag_PatElectron/" + f"resonanceZ_Tag_pat_Probe_{type_}Electron_{var}_{region}").Clone()
    hist_MC = fileMC.Get("DQMData/Run 333334/HLT/Run summary/ScoutingOffline/EGamma/TnP/Tag_PatElectron/" +  f"resonanceZ_Tag_pat_Probe_{type_}Electron_{var}_{region}").Clone()
    total_entry = fileMC.Get("DQMData/Run 333334/HLT/Run summary/ScoutingOffline/EGamma/TnP/Tag_PatElectron/total_entry").GetEntries()

    hist_Data.Rebin(5)
    hist_MC.Rebin(5)

    hist_MC.Scale(lumi * 1000  / total_entry * 6077.22)

    hist_Data.Scale(1./hist_Data.Integral())
    hist_MC.Scale(1./hist_MC.Integral())


    x_axis = hist_Data.GetXaxis()
    y_axis = hist_Data.GetYaxis()
    nbinX  = x_axis.GetNbins()
    nbinY  = y_axis.GetNbins()
    x_binnings = [x_axis.GetBinLowEdge(bin_+1) for bin_ in range(nbinX+1)]
    y_binnings = [y_axis.GetBinLowEdge(bin_+1) for bin_ in range(nbinY+1)]

    c = CMS.cmsDiCanvas('', min(x_binnings), max(x_binnings), 0, max(hist_Data.GetMaximum(), hist_MC.GetMaximum()) * 1.2, 0.0, 2.0, var, 'A.U.', 'data/pred', square = CMS.kSquare, extraSpace=0.0, iPos=0)
    legend = CMS.cmsLeg(0.6, 0.6, 0.8, 0.8, textSize=0.04)

    c.SetLogy()

    c.cd(1)
    CMS.cmsDraw(hist_MC,   'HIST', lcolor = ROOT.kOrange, mcolor = ROOT.kOrange, msize = 1, lwidth = 2, fstyle = 1001, fcolor = ROOT.kOrange)
    CMS.cmsDraw(hist_Data, 'P E', lcolor = ROOT.kBlack, mcolor = ROOT.kBlack, msize = 1, lwidth = 2, fstyle = 0)

    legend.AddEntry(hist_MC, f"Z#rightarrow ee Summer 24 MC", 'F')
    legend.AddEntry(hist_Data, f"PFMonitoring {era}", 'PL')
     

    c.cd(2)
    hist_ratio = hist_Data.Clone()
    hist_ratio.Divide(hist_MC)
    CMS.cmsDraw(hist_ratio, 'P E', lcolor = ROOT.kBlack, mcolor = ROOT.kBlack, msize = 1, lwidth = 2, fstyle = 0)

    plotdir = "mc_plot"
    os.makedirs(plotdir, exist_ok = True)
    CMS.SaveCanvas(c, os.path.join(plotdir, f'{type_}Electron_{var}_{region}.png'), close = False)
    CMS.SaveCanvas(c, os.path.join(plotdir, f'{type_}Electron_{var}_{region}.pdf'))


for variable in ["Pt", "HoverE", "OoEMOoP", "MissingHits", "RelEcalIsolation", "RelHcalIsolation", "RelTrackIsolation", "SigmaIetaIeta", "Trackfbrem", "dEtaIn", "dPhiIn"]:
  for type_ in ["pat", "sct"]:
    for region in ["Barrel", "Endcap"]:
      DrawMC_Comparison(file_Data, file_MC, type_, variable, region, ERA)
      DrawOffline_Comparison(file_Data, variable, region, ERA)

DrawMass_perEra({"eraF": file_F, "eraG": file_G})
#file_Dict = {"eraF": file_F, "eraG": file_G, "eraH": file_H}
file_Dict = {
#   "eraB": ROOT.TFile.Open("/eos/user/t/tihsu/database/DPNote/ScoutingPFMonitorData2025/eraB/data.root"),
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

    for fireObj in ["", "_fireTrigObj"]:
      DrawEfficiency_perEra(file_Dict, ERA, 'Barrel', [DST_ + fireObj], "pat", pt_order, name_ = DST_ + fireObj)
      DrawEfficiency_perEra(file_Dict, ERA, 'Barrel', [DST_ + fireObj], "sct", pt_order, name_ = DST_ + fireObj)
      DrawEfficiency_perEra(file_Dict, ERA, 'Endcap', [DST_ + fireObj], "pat", pt_order, name_ = DST_ + fireObj)
      DrawEfficiency_perEra(file_Dict, ERA, 'Endcap', [DST_ + fireObj], "sct", pt_order, name_ = DST_ + fireObj)
      DrawEfficiency_perEra(file_Dict, ERA, 'Full', [DST_ + fireObj], "sct", pt_order, eta = True, name_ = DST_ + fireObj)
      DrawEfficiency_perEra(file_Dict, ERA, 'Full', [DST_ + fireObj], "pat", pt_order, eta = True, name_ = DST_ + fireObj)



for pt_order in ["leading", "subleading"]:
  for DST_ in DST:
    fireObj = "_fireTrigObj"
    DrawEfficiency(file_Data, file_Data, ERA, 'Barrel', [DST_ + fireObj], "sct", pt_order, name = DST_ + fireObj, tagging_type = "sct")
    DrawEfficiency(file_Data, file_Data, ERA, 'Endcap', [DST_ + fireObj], "sct", pt_order, name = DST_ + fireObj, tagging_type = "sct")
    DrawEfficiency(file_Data, file_Data, ERA, 'Full', [DST_ + fireObj], "sct", pt_order, eta = True, name = DST_ + fireObj, tagging_type = "sct")

#    for fireObj in ["", "_fireTrigObj"]:
#      DrawEfficiency(file_Data, file_MC, ERA, 'Barrel', [DST_ + fireObj], "pat", pt_order, name = DST_ + fireObj)
#      DrawEfficiency(file_Data, file_MC, ERA, 'Barrel', [DST_ + fireObj], "sct", pt_order, name = DST_ + fireObj)
#      DrawEfficiency(file_Data, file_MC, ERA, 'Endcap', [DST_ + fireObj], "pat", pt_order, name = DST_ + fireObj)
#      DrawEfficiency(file_Data, file_MC, ERA, 'Endcap', [DST_ + fireObj], "sct", pt_order, name = DST_ + fireObj)
#      DrawEfficiency(file_Data, file_MC, ERA, 'Full', [DST_ + fireObj], "sct", pt_order, eta = True, name = DST_ + fireObj)
#      DrawEfficiency(file_Data, file_MC, ERA, 'Full', [DST_ + fireObj], "pat", pt_order, eta = True, name = DST_ + fireObj)

  for name in L1_dict:
    DrawEfficiency(file_Data, file_MC, ERA, 'Barrel', L1_dict[name], "pat", pt_order, name = name)
    DrawEfficiency(file_Data, file_MC, ERA, 'Barrel', L1_dict[name], "sct", pt_order, name = name)
    DrawEfficiency(file_Data, file_MC, ERA, 'Endcap', L1_dict[name], "pat", pt_order, name = name)
    DrawEfficiency(file_Data, file_MC, ERA, 'Endcap', L1_dict[name], "sct", pt_order, name = name)
