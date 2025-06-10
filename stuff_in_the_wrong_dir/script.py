import ROOT

file_even = ROOT.TFile.Open("../run/submitDir-2025-06-10-1240-5ba5/hist-dataset.root")
file_odd = ROOT.TFile.Open("../run/submitDir-2025-06-10-1241-01af/hist-dataset.root")

hist_phi_cp_tau_even = file_even.Get("phi_CP_tau_pi")
hist_phi_cp_neutrino_even = file_even.Get("phi_CP_neutrino_pi")

hist_phi_cp_tau_odd = file_odd.Get("phi_CP_tau_pi")
hist_phi_cp_neutrino_odd = file_odd.Get("phi_CP_neutrino_pi")

for hist in [hist_phi_cp_tau_even, hist_phi_cp_neutrino_even, hist_phi_cp_tau_odd, hist_phi_cp_neutrino_odd]:
    integral = hist.Integral()
    hist.Scale(1.0 / integral)  # Normalization
    hist.SetLineWidth(1)
    hist.SetStats(0)

canvas = ROOT.TCanvas()
canvas.SetLeftMargin(0.15)

global_max_y = max(
    hist_phi_cp_tau_even.GetMaximum(),
    hist_phi_cp_neutrino_even.GetMaximum(),
    hist_phi_cp_tau_odd.GetMaximum(),
    hist_phi_cp_neutrino_odd.GetMaximum()
)

global_min_y = min(
    hist_phi_cp_tau_even.GetMinimum(),
    hist_phi_cp_neutrino_even.GetMinimum(),
    hist_phi_cp_tau_odd.GetMinimum(),
    hist_phi_cp_neutrino_odd.GetMinimum()
)

hist_phi_cp_tau_even.SetMaximum(global_max_y * 1.05)
hist_phi_cp_tau_even.SetMinimum(global_min_y * 0.95)


hist_phi_cp_tau_even.SetLineColor(ROOT.kBlue)
hist_phi_cp_tau_even.SetTitle("Phi CP Histograms")
hist_phi_cp_tau_even.GetXaxis().SetTitle("#phi_{CP}")
hist_phi_cp_tau_even.GetYaxis().SetTitle("~ 1/#sigma_{tot} d#sigma/d#phi_{CP}")
hist_phi_cp_tau_even.Draw("HIST")

hist_phi_cp_neutrino_even.SetLineColor(ROOT.kRed)
hist_phi_cp_neutrino_even.Draw("HIST SAME")

hist_phi_cp_tau_odd.SetLineColor(ROOT.kGreen)
hist_phi_cp_tau_odd.Draw("HIST SAME")

hist_phi_cp_neutrino_odd.SetLineColor(ROOT.kMagenta)
hist_phi_cp_neutrino_odd.Draw("HIST SAME")

legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
legend.AddEntry(hist_phi_cp_tau_even, "Phi CP (tau, even)", "l")
legend.AddEntry(hist_phi_cp_neutrino_even, "Phi CP (neutrino, even)", "l")
legend.AddEntry(hist_phi_cp_tau_odd, "Phi CP (tau, odd)", "l")
legend.AddEntry(hist_phi_cp_neutrino_odd, "Phi CP (neutrino, odd)", "l")
legend.Draw()

canvas.Update()

canvas.SaveAs("histogram.png")

file_even.Close()
file_odd.Close()
