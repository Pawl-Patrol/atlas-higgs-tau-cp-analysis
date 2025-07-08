import ROOT

file_even_hadhad = ROOT.TFile.Open("../run/cp-even-hadhad/hist-dataset.root")
file_odd_hadhad = ROOT.TFile.Open("../run/cp-odd-hadhad/hist-dataset.root")
file_even_hadlep = ROOT.TFile.Open("../run/cp-even-hadlep/hist-dataset.root")
file_odd_hadlep = ROOT.TFile.Open("../run/cp-odd-hadlep/hist-dataset.root")

hist_even_phi_CP_hadhad = file_even_hadhad.Get("phi_CP_pion")
hist_even_phi_CP_reco_hadhad = file_even_hadhad.Get("phi_CP_pion_jet_reco")
hist_odd_phi_CP_hadhad = file_odd_hadhad.Get("phi_CP_pion")
hist_odd_phi_CP_reco_hadhad = file_odd_hadhad.Get("phi_CP_pion_jet_reco")
hist_even_phi_CP_hadlep = file_even_hadlep.Get("phi_CP_pion")
hist_even_phi_CP_reco_hadlep = file_even_hadlep.Get("phi_CP_pion_jet_reco")
hist_odd_phi_CP_hadlep = file_odd_hadlep.Get("phi_CP_pion")
hist_odd_phi_CP_reco_hadlep = file_odd_hadlep.Get("phi_CP_pion_jet_reco")

hist_even_sum = hist_even_phi_CP_hadhad.Clone("phi_CP_sum_even")
hist_even_sum.Add(hist_even_phi_CP_hadlep)
hist_odd_sum = hist_odd_phi_CP_hadhad.Clone("phi_CP_sum_odd")
hist_odd_sum.Add(hist_odd_phi_CP_hadlep)

hist_even_sum_reco = hist_even_phi_CP_reco_hadhad.Clone("phi_CP_sum_reco_even")
hist_even_sum_reco.Add(hist_even_phi_CP_reco_hadlep)
hist_odd_sum_reco = hist_odd_phi_CP_reco_hadhad.Clone("phi_CP_sum_reco_odd")
hist_odd_sum_reco.Add(hist_odd_phi_CP_reco_hadlep)

histograms = {
    # "#phi_{CP} (tau-pion, even)": file_even.Get("phi_CP_tau_pi"),
    # "#phi_{CP} (tau-pion, odd)": file_odd.Get("phi_CP_tau_pi"),
    # "#phi_{CP} (neutrino-pion, even)": file_even.Get("phi_CP_neutrino_pi"),
    # "#phi_{CP} (neutrino-pion, odd)": file_odd.Get("phi_CP_neutrino_pi"),
    # -------------------
    "#phi_{CP} (hadhad, even, truth)": file_even_hadhad.Get("phi_CP_pion"),
    "#phi_{CP} (hadhad, odd, truth)": file_odd_hadhad.Get("phi_CP_pion"),
    # "#phi_{CP} (hadlep, even, truth)": file_even_hadlep.Get("phi_CP_pion"),
    # "#phi_{CP} (hadlep, odd, truth)": file_odd_hadlep.Get("phi_CP_pion"),
    # -------------------
    "#phi_{CP} (hadhad, even, reco)": file_even_hadhad.Get("phi_CP_pion_jet_reco"),
    "#phi_{CP} (hadhad, odd, reco)": file_odd_hadhad.Get("phi_CP_pion_jet_reco"),
    # "#phi_{CP} (hadlep, even, reco)": file_even_hadlep.Get("phi_CP_pion_jet_reco"),
    # "#phi_{CP} (hadlep, odd, reco)": file_odd_hadlep.Get("phi_CP_pion_jet_reco"),
    # -------------------
    # "#phi_{CP} (sum, truth, even)": hist_even_sum,
    # "#phi_{CP} (sum, truth, odd)": hist_odd_sum,
    # "#phi_{CP} (sum, reco, even)": hist_even_sum_reco,
    # "#phi_{CP} (sum, reco, odd)": hist_odd_sum_reco,
}

BINS = 50
DPHI = 2 * ROOT.TMath.Pi() / BINS

for name, hist in histograms.items():
    integral = hist.Integral()
    if integral == 0:
        print(f"Warning: Histogram '{name}' has zero integral, skipping normalization.")
        continue
    else:
        print(f"Normalizing histogram '{name}' with integral {integral:.0f}.")
    hist.Scale(1 / (integral * DPHI))  # Normalization
    hist.SetLineWidth(1)
    hist.SetStats(0)

canvas = ROOT.TCanvas()
canvas.SetLeftMargin(0.15)

global_max_y = max(map(lambda h: h.GetMaximum(), histograms.values()))
global_min_y = min(map(lambda h: h.GetMinimum(), histograms.values()))

print(f"Global max y: {global_max_y}, Global min y: {global_min_y}")

COLORS = [
    ROOT.kBlue,
    ROOT.kBlack,
    ROOT.kRed,
    ROOT.kGreen,
    ROOT.kMagenta,
    ROOT.kCyan,
    ROOT.kOrange,
    ROOT.kViolet,
    ROOT.kPink,
    ROOT.kGray,
]

for i, hist in enumerate(histograms.values()):
    if i == 0:
        hist.SetMaximum(global_max_y * 1.05)
        hist.SetMinimum(global_min_y * 0.95)
        hist.SetTitle("Phi CP Histograms")
        hist.GetXaxis().SetTitle("#phi_{CP}")
        hist.GetYaxis().SetTitle("1/#sigma_{tot} d#sigma/d#phi_{CP}")
    hist.SetLineColor(COLORS[i])
    hist.Draw("HIST" if i == 0 else "HIST SAME")

legend = ROOT.TLegend(0.7, 0.8, 0.9, 0.9)
for name, hist in histograms.items():
    legend.AddEntry(hist, name, "l")
legend.Draw()

canvas.Update()

canvas.SetCanvasSize(1920, 1080)
canvas.SaveAs("histogram.png")

file_even_hadhad.Close()
file_odd_hadhad.Close()
file_even_hadlep.Close()
file_odd_hadlep.Close()
