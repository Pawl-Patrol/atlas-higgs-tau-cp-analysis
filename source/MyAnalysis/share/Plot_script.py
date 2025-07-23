#!/usr/bin/env python3

import inquirer
from Run_script import SAMPLES
import math

samples = inquirer.checkbox(
    "Which samples do you want to load?",
    choices=list(SAMPLES.keys()),
)

lower = 0.0  # float(inquirer.text("What's the lower limit?", default="0.000"))
upper = (
    2  # * math.pi # float(inquirer.text("What's the upper limit?", default="6.283"))
)

x_label = "#phi_{CP}"  # inquirer.text("What's the x label?", default="#phi_{{CP}}")
y_label = "1/#sigma_{tot} d#sigma/d#phi_{CP}"  # inquirer.text("What's the y label?", default="1/#sigma_{{tot}} d#sigma/d#phi_{{CP}}")

BINS = 100
DPHI = (upper - lower) / BINS

scaling_factor = 1.0 / DPHI

print("Loading ntuples...")

import ROOT

files = [
    ROOT.TFile.Open(f"/srv/run/{SAMPLES[sample]}/data-ANALYSIS/dataset.root")
    for sample in samples
]
trees = [file.Get("tau_analysis") for file in files]


branch_names = set()
for tree in trees:
    for branch in tree.GetListOfBranches():
        branch_name = branch.GetName()
        valid = False
        for entry in tree:
            value = getattr(entry, branch_name)
            if value is None:
                continue
            if value < lower:
                continue
            if value > upper:
                continue
            valid = True
            break
        if not valid:
            continue
        branch_names.add(branch_name)

branch_names = sorted(branch_names)

branches = inquirer.checkbox(
    "Which branches do you want to plot?",
    choices=list(branch_names),
)

cuts = {}
while True:
    branch = inquirer.list_input(
        "Select branches to apply cuts on.",
        choices=["No, I'm good!"] + list(branch_names),
        default=None,
    )
    if branch == "No, I'm good!":
        break
    lower_cut = inquirer.text(f"Apply lower limit {branch}?", default="")
    upper_cut = inquirer.text(f"Apply upper limit {branch}?", default="")
    cuts[branch] = (
        float(lower_cut) if lower_cut else -math.inf,
        float(upper_cut) if upper_cut else math.inf,
    )


histograms = {}

for tree, sample in zip(trees, samples):
    for branch in branches:
        hist_name = (
            SAMPLES[sample]
            .replace("-hadhad", "")
            .replace("-hadlep", "")
            .replace("-lephad", "")
        )
        hist_title = branch.replace("phiCP_", "")
        hist = ROOT.TH1F(hist_name, hist_name, BINS, lower, upper)
        hist.SetTitle(hist_title)

        for entry in tree:
            for cut_branch, (cut_lower, cut_upper) in cuts.items():
                cut_value = abs(getattr(entry, cut_branch))
                if cut_lower and cut_value < cut_lower:
                    break
                if cut_upper and cut_value > cut_upper:
                    break
            else:  # Only proceed if all cuts are satisfied
                value = getattr(entry, branch)
                hist.Fill(value)

        integral = hist.Integral()
        if integral == 0:
            print(
                f"Warning: Histogram '{hist_name}' has zero integral, skipping normalization."
            )
        else:
            print(f"Normalizing histogram '{hist_name}' with integral {integral:.0f}.")
            hist.Scale(scaling_factor / integral)  # Normalization
        hist.SetLineWidth(1)
        hist.SetStats(0)
        histograms[hist_name] = hist

canvas = ROOT.TCanvas()
canvas.SetLeftMargin(0.15)

global_max_y = 0.1  # max(map(lambda h: h.GetMaximum(), histograms.values()))
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
        hist.GetXaxis().SetTitle(x_label)
        hist.GetYaxis().SetTitle(y_label)
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
