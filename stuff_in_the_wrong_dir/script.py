import ROOT

file_even = ROOT.TFile.Open("../run/cp-even/hist-dataset.root")
file_odd = ROOT.TFile.Open("../run/cp-odd/hist-dataset.root")

histograms = {
    "#phi_{CP} (tau, even)": file_even.Get("phi_CP_tau_pi"),
    "#phi_{CP} (tau, odd)": file_odd.Get("phi_CP_tau_pi"),
    "#phi_{CP} (neutrino, even)": file_even.Get("phi_CP_neutrino_pi"),
    "#phi_{CP} (neutrino, odd)": file_odd.Get("phi_CP_neutrino_pi"),
    # "#phi_{CP} (pion, even)": file_even.Get("phi_CP_pion"),
    # "#phi_{CP} (pion, odd)": file_odd.Get("phi_CP_pion"),
}

for hist in histograms.values():
    integral = hist.Integral()
    hist.Scale(1.0 / integral)  # Normalization
    hist.SetLineWidth(1)
    hist.SetStats(0)

canvas = ROOT.TCanvas()
canvas.SetLeftMargin(0.15)

global_max_y = max(map(lambda h: h.GetMaximum(), histograms.values()))
global_min_y = min(map(lambda h: h.GetMinimum(), histograms.values()))

COLORS = [ROOT.kBlue, ROOT.kRed, ROOT.kBlack, ROOT.kGreen, ROOT.kMagenta, ROOT.kCyan]

for i, hist in enumerate(histograms.values()):
    if i == 0:
        hist.SetMaximum(global_max_y * 1.05)
        hist.SetMinimum(global_min_y * 0.95)
        hist.SetTitle("Phi CP Histograms")
        hist.GetXaxis().SetTitle("#phi_{CP}")
        hist.GetYaxis().SetTitle("~ 1/#sigma_{tot} d#sigma/d#phi_{CP}")
    hist.SetLineColor(COLORS[i])
    hist.Draw("HIST" if i == 0 else "HIST SAME")

legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
for name, hist in histograms.items():
    legend.AddEntry(hist, name, "l")
legend.Draw()

canvas.Update()

canvas.SaveAs("histogram.png")

file_even.Close()
file_odd.Close()
