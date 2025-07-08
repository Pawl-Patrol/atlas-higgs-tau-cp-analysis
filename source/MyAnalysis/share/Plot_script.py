#!/usr/bin/env python3

import inquirer
from Run_script import SAMPLES


samples = inquirer.checkbox(
    "Which samples do you want to load?",
    choices=list(SAMPLES.keys()),
)

import ROOT

files = [
    ROOT.TFile.Open(f"/srv/run/{SAMPLES[sample]}/data-ANALYSIS/dataset.root")
    for sample in samples
]
trees = [file.Get("tau_analysis") for file in files]
branch_names = [branch.GetName() for branch in trees[0].GetListOfBranches()]


branches = inquirer.checkbox(
    "Which branches do you want to plot?",
    choices=branch_names,
)

BINS = 50
LIMIT = 2 * ROOT.TMath.Pi()
DPHI = LIMIT / BINS

histograms = {}

for tree, sample in zip(trees, samples):
    for branch in branches:
        hist_name = f"#phi_{{CP}} ({SAMPLES[sample]}/{branch.replace('phiCP_', '')})"
        hist = ROOT.TH1F(hist_name, hist_name, BINS, 0, LIMIT)

        for entry in tree:
            value = getattr(entry, branch)
            hist.Fill(value)

        integral = hist.Integral()
        if integral == 0:
            print(
                f"Warning: Histogram '{hist_name}' has zero integral, skipping normalization."
            )
        else:
            print(f"Normalizing histogram '{hist_name}' with integral {integral:.0f}.")
            hist.Scale(1 / (integral * DPHI))  # Normalization
        hist.SetLineWidth(1)
        hist.SetStats(0)
        histograms[hist_name] = hist

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

for file in files:
    file.Close()
