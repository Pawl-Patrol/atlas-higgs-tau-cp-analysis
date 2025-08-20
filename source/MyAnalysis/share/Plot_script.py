#!/usr/bin/env python3

import inquirer
from Run_script import SAMPLES
import math
import re

# Define mode constants
HISTOGRAM_MODE = "Histogram Mode"
HISTOGRAM_FIT_MODE = "Histogram + Fit Mode"
HISTOGRAM_FIT_UNCERTAINTY_MODE = "Histogram + Fit + Uncertainties Mode"
SCATTER_MODE = "X-Y Scatter Plot Mode"

# Choose the plot mode
plot_mode = inquirer.list_input(
    "Select plot mode:",
    choices=[
        HISTOGRAM_MODE,
        HISTOGRAM_FIT_MODE,
        HISTOGRAM_FIT_UNCERTAINTY_MODE,
        SCATTER_MODE,
    ],
    default=HISTOGRAM_MODE,
)

if plot_mode in [HISTOGRAM_MODE, HISTOGRAM_FIT_MODE, HISTOGRAM_FIT_UNCERTAINTY_MODE]:
    # Current histogram mode
    samples = inquirer.checkbox(
        "Which samples do you want to load?",
        choices=list(SAMPLES.keys()),
    )

    lower = 0.0  # float(inquirer.text("What's the lower limit?", default="0.000"))
    upper = (
        2 * math.pi
    )  # float(inquirer.text("What's the upper limit?", default="6.283"))

    x_label = "#phi_{CP}"  # inquirer.text("What's the x label?", default="#phi_{{CP}}")
    y_label = "1/#sigma_{tot} d#sigma/d#phi_{CP}"  # inquirer.text("What's the y label?", default="1/#sigma_{{tot}} d#sigma/d#phi_{{CP}}")

    BINS = 50
    DPHI = (upper - lower) / BINS

    scaling_factor = 1.0 / DPHI

else:
    # X-Y scatter plot mode
    sample = inquirer.list_input(
        "Select one sample for X-Y plotting:", choices=list(SAMPLES.keys())
    )

    # Add binning parameters for 2D heatmap
    BINS_X = 50
    BINS_Y = 50

    samples = [sample]  # Convert to list for compatibility

print("Loading ntuples...")

import ROOT
import os

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
            if "lower" in locals() and value < lower:
                continue
            if "upper" in locals() and value > upper:
                continue
            valid = True
            break
        if not valid:
            continue
        branch_names.add(branch_name)

branch_names = sorted(branch_names)

if plot_mode in [HISTOGRAM_MODE, HISTOGRAM_FIT_MODE, HISTOGRAM_FIT_UNCERTAINTY_MODE]:
    branches = inquirer.checkbox(
        "Which branches do you want to plot?",
        choices=list(branch_names),
    )
else:
    # X-Y mode: select two branches
    x_branch = inquirer.list_input(
        "Select branch for X-axis:", choices=list(branch_names)
    )
    y_branch = inquirer.list_input(
        "Select branch for Y-axis:", choices=list(branch_names)
    )

    x_label = inquirer.text(f"X-axis label for {x_branch}", default="#phi_{{CP}}")
    y_label = inquirer.text(f"Y-axis label for {y_branch}", default="#phi_{{CP}}")

    branches = [x_branch, y_branch]  # For compatibility with existing code

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
heatmaps = {}


if plot_mode in [HISTOGRAM_MODE, HISTOGRAM_FIT_MODE, HISTOGRAM_FIT_UNCERTAINTY_MODE]:
    hist_title_template = inquirer.text(
        "Enter histogram title template",
        default="",
        # default="{{mode1}}-{{mode2}}",
    )

    hist_name_template = inquirer.text(
        "Enter histogram name template", default="CP {{cp_nature}}"
    )

    # Original histogram mode
    for tree, sample in zip(trees, samples):
        cp_nature, decay_mode, mass = re.match(
            r"^cp-(.+?)-(.+?)(?:-H(.+?))?$", SAMPLES[sample]
        ).groups()
        for branch in branches:
            match = re.match(r"^phiCP_(.+?)_(.+?)_(.+?)$", branch)
            args = dict(
                cp_nature=cp_nature,
                decay_mode=decay_mode,
                mass=mass or "125",
            )
            if match:
                mode1, mode2, info_type = match.groups()
                args.update(mode1=mode1, mode2=mode2, info_type=info_type)

            hist_title = hist_title_template.format(**args)
            hist_name = hist_name_template.format(**args)

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
                print(
                    f"Normalizing histogram '{hist_name}' with integral {integral:.0f}."
                )
                hist.Scale(scaling_factor / integral)  # Normalization
            hist.SetLineWidth(1)
            hist.SetStats(0)
            histograms[hist_name] = hist

else:
    # X-Y heatmap mode (2D binning)
    for tree, sample in zip(trees, samples):
        x_values = []
        y_values = []

        for entry in tree:
            # Apply cuts
            cuts_satisfied = True
            for cut_branch, (cut_lower, cut_upper) in cuts.items():
                cut_value = abs(getattr(entry, cut_branch))
                if cut_lower and cut_value < cut_lower:
                    cuts_satisfied = False
                    break
                if cut_upper and cut_value > cut_upper:
                    cuts_satisfied = False
                    break

            if cuts_satisfied:
                x_val = getattr(entry, x_branch)
                y_val = getattr(entry, y_branch)
                if x_val is not None and y_val is not None:
                    x_values.append(x_val)
                    y_values.append(y_val)

        if len(x_values) > 0:
            # Use fixed ranges for both axes (0 to 2π)
            x_min, x_max = 0.0, 2 * math.pi
            y_min, y_max = 0.0, 2 * math.pi

            heatmap_name = f"heatmap_{sample}"
            heatmap_title = ""
            # heatmap_title = (
            # SAMPLES[sample]
            # .replace("-hadhad", "")
            # .replace("-hadlep", "")
            # .replace("-lephad", "")
            # )

            # Create 2D histogram
            heatmap = ROOT.TH2F(
                heatmap_name, heatmap_title, BINS_X, x_min, x_max, BINS_Y, y_min, y_max
            )

            out_of_range_counter = 0
            # Fill the 2D histogram
            for x_val, y_val in zip(x_values, y_values):
                heatmap.Fill(x_val, y_val)

            heatmap.SetStats(0)  # Disable statistics box
            heatmaps[heatmap_name] = heatmap

            integral = heatmap.Integral()
            print(f"2D heatmap with {BINS_X}x{BINS_Y} with integral {integral:.0f}")

# Define the fit function: y(x) = A*cos(B*x + C) + D (for histogram + fit modes)
if plot_mode in [HISTOGRAM_FIT_MODE, HISTOGRAM_FIT_UNCERTAINTY_MODE]:
    fit_function = ROOT.TF1("cosine_fit", "[0]*cos([1]*x + [2]) + [3]", lower, upper)
    fit_function.SetParNames("A", "B", "C", "D")
    # Set initial parameter estimates (you may need to adjust these based on your data)
    fit_function.SetParameter(0, 1.0)  # A (amplitude)
    fit_function.SetParameter(1, 1.0)  # B (frequency)
    fit_function.SetParameter(2, 0.0)  # C (phase)
    fit_function.SetParameter(3, 0.0)  # D (offset)

# Store fit results for each histogram
fit_results = {}

# Create canvas with ratio plot for uncertainty mode
if plot_mode == HISTOGRAM_FIT_UNCERTAINTY_MODE:
    canvas = ROOT.TCanvas("canvas", "canvas", 800, 800)
    canvas.Divide(1, 2)

    # Main plot pad (top)
    pad1 = canvas.cd(1)
    pad1.SetPad(0, 0.3, 1, 1)
    pad1.SetBottomMargin(0.02)
    pad1.SetLeftMargin(0.18)  # Match ratio pad margin

    # Ratio plot pad (bottom)
    pad2 = canvas.cd(2)
    pad2.SetPad(0, 0, 1, 0.3)
    pad2.SetTopMargin(0.02)
    pad2.SetBottomMargin(0.4)
    pad2.SetLeftMargin(0.18)  # Increased left margin for y-axis labels

    canvas.cd(1)  # Start with main pad
else:
    canvas = ROOT.TCanvas()
    canvas.SetLeftMargin(0.15)

if plot_mode in [HISTOGRAM_MODE, HISTOGRAM_FIT_MODE, HISTOGRAM_FIT_UNCERTAINTY_MODE]:
    global_max_y = max(map(lambda h: h.GetMaximum(), histograms.values()))
    global_min_y = min(map(lambda h: h.GetMinimum(), histograms.values()))
else:
    # For scatter plots, we don't need global max/min for histograms
    pass

if plot_mode in [HISTOGRAM_MODE, HISTOGRAM_FIT_MODE, HISTOGRAM_FIT_UNCERTAINTY_MODE]:
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

if plot_mode in [HISTOGRAM_MODE, HISTOGRAM_FIT_MODE, HISTOGRAM_FIT_UNCERTAINTY_MODE]:
    # Store ratio histograms for uncertainty mode
    ratio_histograms = {}

    # Original histogram plotting and fitting code
    for i, (name, hist) in enumerate(histograms.items()):
        if i == 0:
            hist.SetMaximum(global_max_y * 1.05)
            hist.SetMinimum(global_min_y * 0.95)
            hist.GetXaxis().SetTitle(x_label)
            hist.GetYaxis().SetTitle(y_label)

            # For uncertainty mode, hide x-axis labels on main plot
            if plot_mode == HISTOGRAM_FIT_UNCERTAINTY_MODE:
                hist.GetXaxis().SetLabelSize(0)
                hist.GetXaxis().SetTitleSize(0)

        hist.SetLineColor(COLORS[i])

        # For uncertainty mode, draw with error bars
        if plot_mode == HISTOGRAM_FIT_UNCERTAINTY_MODE:
            hist.Draw("E1" if i == 0 else "E1 SAME")
        else:
            hist.Draw("HIST" if i == 0 else "HIST SAME")

        # Perform the fit for both fit modes
        if plot_mode in [HISTOGRAM_FIT_MODE, HISTOGRAM_FIT_UNCERTAINTY_MODE]:
            print(f"Fitting histogram '{name}' with cosine function...")

            # Create a unique fit function for each histogram
            fit_func = ROOT.TF1(
                f"fit_{name}", "[0]*cos([1]*x + [2]) + [3]", lower, upper
            )
            fit_func.SetParNames("A", "B", "C", "D")

            # Set initial parameter estimates based on histogram properties
            max_val = hist.GetMaximum()
            min_val = hist.GetMinimum()
            amplitude_guess = (max_val - min_val) * 0.5
            offset_guess = (max_val + min_val) * 0.5

            fit_func.SetParameter(0, amplitude_guess)  # A (amplitude)
            fit_func.SetParameter(1, 1.0)  # B (frequency)
            fit_func.SetParameter(2, math.pi * 0.5)  # C (phase)
            fit_func.SetParameter(3, offset_guess)  # D (vertical offset)

            # Set reasonable parameter limits if needed
            fit_func.SetParLimits(
                0, amplitude_guess * 0.5, amplitude_guess * 2
            )  # A should be positive and reasonable
            fit_func.SetParLimits(1, 0.5, 2.0)  # B
            fit_func.SetParLimits(2, 0, math.pi * 2)  # C
            fit_func.SetParLimits(3, min_val, max_val)  # D

            # Perform the fit
            fit_result = hist.Fit(fit_func, "S")  # "S" option returns fit result

            # Set the fit function color to match the histogram
            fit_func.SetLineColor(COLORS[i])
            fit_func.SetLineStyle(2)  # Dashed line for fits

            # Draw the fit function
            fit_func.Draw("SAME")

            # Store fit results
            if fit_result and fit_result.IsValid():
                fit_results[name] = {
                    "A": fit_func.GetParameter(0),
                    "A_err": fit_func.GetParError(0),
                    "B": fit_func.GetParameter(1),
                    "B_err": fit_func.GetParError(1),
                    "C": fit_func.GetParameter(2),
                    "C_err": fit_func.GetParError(2),
                    "D": fit_func.GetParameter(3),
                    "D_err": fit_func.GetParError(3),
                    "chi2_ndf": (
                        fit_result.Chi2() / fit_result.Ndf()
                        if fit_result.Ndf() > 0
                        else -1
                    ),
                    "fit_func": fit_func,  # Store function for ratio calculation
                }

            # Create ratio histogram for uncertainty mode
            if (
                plot_mode == HISTOGRAM_FIT_UNCERTAINTY_MODE
                and fit_result
                and fit_result.IsValid()
            ):
                ratio_hist = hist.Clone(f"ratio_{name}")
                ratio_hist.SetTitle("")
                ratio_hist.Reset()

                # Calculate data/fit ratio for each bin
                for bin_idx in range(1, hist.GetNbinsX() + 1):
                    data_value = hist.GetBinContent(bin_idx)
                    data_error = hist.GetBinError(bin_idx)
                    bin_center = hist.GetBinCenter(bin_idx)
                    fit_value = fit_func.Eval(bin_center)

                    if fit_value > 0 and data_value > 0:
                        ratio = data_value / fit_value
                        ratio_error = data_error / fit_value if fit_value > 0 else 0
                        ratio_hist.SetBinContent(bin_idx, ratio)
                        ratio_hist.SetBinError(bin_idx, ratio_error)

                ratio_hist.SetLineColor(COLORS[i])
                ratio_hist.SetMarkerColor(COLORS[i])
                ratio_hist.SetMarkerStyle(20)
                ratio_hist.SetMarkerSize(0.6)  # Smaller markers
                ratio_histograms[name] = ratio_hist

    # Draw ratio plots for uncertainty mode
    if plot_mode == HISTOGRAM_FIT_UNCERTAINTY_MODE and ratio_histograms:
        canvas.cd(2)  # Switch to ratio pad

        first_ratio = True
        for i, (name, ratio_hist) in enumerate(ratio_histograms.items()):
            if first_ratio:
                ratio_hist.SetMaximum(1.5)
                ratio_hist.SetMinimum(0.5)
                ratio_hist.GetXaxis().SetTitle(x_label)
                ratio_hist.GetYaxis().SetTitle("Data/Fit")
                ratio_hist.GetYaxis().SetTitleSize(0.12)
                ratio_hist.GetYaxis().SetLabelSize(0.1)
                ratio_hist.GetXaxis().SetTitleSize(0.12)
                ratio_hist.GetXaxis().SetLabelSize(0.1)
                ratio_hist.GetYaxis().SetTitleOffset(
                    0.7
                )  # Increased for better spacing
                ratio_hist.GetXaxis().SetTitleOffset(1.0)
                ratio_hist.GetYaxis().SetNdivisions(505)  # Better tick spacing
                ratio_hist.Draw("E1")
                first_ratio = False
            else:
                ratio_hist.Draw("E1 SAME")

        # Draw horizontal line at y=1
        line = ROOT.TLine(lower, 1.0, upper, 1.0)
        line.SetLineStyle(2)
        line.SetLineColor(ROOT.kBlack)
        line.Draw("SAME")

        canvas.cd(1)  # Switch back to main pad

elif plot_mode == SCATTER_MODE:
    # X-Y heatmap mode
    for i, (name, heatmap) in enumerate(heatmaps.items()):
        if i == 0:
            heatmap.GetXaxis().SetTitle(x_label)
            heatmap.GetYaxis().SetTitle(y_label)
            heatmap.GetZaxis().SetTitle("Entries")

            # Create custom color palette with white for low values
            import array

            # Define the color stops (position along the palette)
            stops = array.array("d", [0.00, 1.00])

            # Define the RGB values for each stop (white -> blue)
            red = array.array("d", [1.00, 0.00])
            green = array.array("d", [1.00, 0.00])
            blue = array.array("d", [1.00, 1.00])

            # Create the palette
            ROOT.TColor.CreateGradientColorTable(
                len(stops), stops, red, green, blue, 50
            )
            ROOT.gStyle.SetOptStat(0)  # Disable statistics box

            # Set minimum to 0 to ensure zero bins show as white
            heatmap.SetMinimum(0.0)

        # Use COLZ option for color heatmap
        heatmap.Draw("COLZ" if i == 0 else "COLZ SAME")

if plot_mode in [HISTOGRAM_MODE, HISTOGRAM_FIT_MODE, HISTOGRAM_FIT_UNCERTAINTY_MODE]:
    if plot_mode == HISTOGRAM_FIT_UNCERTAINTY_MODE:
        canvas.cd(1)  # Make sure we're on the main pad for the legend

    legend = ROOT.TLegend(0.8, 0.8, 0.9, 0.9)
    for name, hist in histograms.items():
        legend.AddEntry(hist, name, "l")
    legend.Draw()
    canvas.Update()

base = inquirer.text("Enter output filename for canvas.", default="result")

output_filename = f"{base}.png"

# Ensure the output filename is unique
counter = 1
while os.path.exists(output_filename):
    output_filename = f"{base}_{counter:02d}.png"
    counter += 1

# canvas.SetCanvasSize(1920, 1080)
canvas.SaveAs(output_filename)

for file in files:
    file.Close()


# Save fit results to a text file (for both fit modes)
if plot_mode in [HISTOGRAM_FIT_MODE, HISTOGRAM_FIT_UNCERTAINTY_MODE] and fit_results:
    fit_results_file = output_filename.replace(".png", "_fit_results.txt")
    with open(fit_results_file, "w") as f:
        f.write("Cosine Fit Results: y(x) = A*cos(B*x + C) + D\n")
        f.write("=" * 50 + "\n\n")
        for name, result in fit_results.items():
            f.write(f"Histogram: {name}\n")
            f.write(f"  A = {result['A']:.6f} ± {result['A_err']:.6f}\n")
            f.write(f"  B = {result['B']:.6f} ± {result['B_err']:.6f}\n")
            f.write(f"  C = {result['C']:.6f} ± {result['C_err']:.6f}\n")
            f.write(f"  D = {result['D']:.6f} ± {result['D_err']:.6f}\n")
            f.write(f"  χ²/ndf = {result['chi2_ndf']:.6f}\n")
            f.write("\n")
    print(f"Fit results saved to '{fit_results_file}'")
