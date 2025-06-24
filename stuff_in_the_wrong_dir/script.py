import ROOT

file_even = ROOT.TFile.Open("../run/cp-even/hist-dataset.root")
file_odd = ROOT.TFile.Open("../run/cp-odd/hist-dataset.root")

histograms = {
    # "#phi_{CP} (tau-pion, even)": file_even.Get("phi_CP_tau_pi"),
    # "#phi_{CP} (tau-pion, odd)": file_odd.Get("phi_CP_tau_pi"),
    # "#phi_{CP} (neutrino-pion, even)": file_even.Get("phi_CP_neutrino_pi"),
    # "#phi_{CP} (neutrino-pion, odd)": file_odd.Get("phi_CP_neutrino_pi"),
    # "#phi_{CP} (ip-pion, even)": file_even.Get("phi_CP_pion"),
    # "#phi_{CP} (ip-pion, odd)": file_odd.Get("phi_CP_pion"),
    "#phi_{CP} (ip-pion-jets, even)": file_even.Get("phi_CP_pion_jet"),
    "#phi_{CP} (ip-pion-jets, odd)": file_odd.Get("phi_CP_pion_jet"),
    # "#phi_{CP} (ip-pion-jets-vertex, even)": file_even.Get("phi_CP_pion_jet_reco"),
    # "#phi_{CP} (ip-pion-jets-vertex, odd)": file_odd.Get("phi_CP_pion_jet_reco"),
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

canvas.SaveAs("histogram.png")

file_even.Close()
file_odd.Close()
