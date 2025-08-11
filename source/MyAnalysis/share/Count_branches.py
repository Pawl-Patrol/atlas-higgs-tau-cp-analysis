#!/usr/bin/env python3

from Run_script import SAMPLES
import ROOT
import os

print("Counting events per branch in all processed samples...")
print("=" * 80)

# Open all processed sample files and count events per branch
for sample_name, sample_dir in SAMPLES.items():
    file_path = f"/srv/run/{sample_dir}/data-ANALYSIS/dataset.root"

    # Check if file exists
    if not os.path.exists(file_path):
        print(f"{sample_name:<25}: File not found at {file_path}")
        continue

    try:
        # Open the ROOT file
        root_file = ROOT.TFile.Open(file_path)

        if not root_file or root_file.IsZombie():
            print(f"{sample_name:<25}: Failed to open file")
            continue

        # Get the tree
        tree = root_file.Get("tau_analysis")

        if not tree:
            print(f"{sample_name:<25}: Tree 'tau_analysis' not found")
            root_file.Close()
            continue

        print(f"\n{sample_name}:")
        print("-" * 50)

        # Get all branches and count valid entries for each
        branches = tree.GetListOfBranches()
        branch_counts = {}

        for branch in branches:
            branch_name = branch.GetName()
            valid_count = 0

            # Count entries where the branch value is not -99.0 (the default invalid value)
            for entry in tree:
                value = getattr(entry, branch_name)
                if value is not None and value != -99.0:
                    valid_count += 1

            branch_counts[branch_name] = valid_count

        # Sort branches by name for consistent output
        for branch_name in sorted(branch_counts.keys()):
            count = branch_counts[branch_name]
            total_entries = tree.GetEntries()
            percentage = (count / total_entries * 100) if total_entries > 0 else 0
            print(
                f"  {branch_name:<30}: {count:>8,} / {total_entries:>8,} ({percentage:>5.1f}%)"
            )

        # Close the file
        root_file.Close()

    except Exception as e:
        print(f"{sample_name:<25}: Error - {str(e)}")

print("\n" + "=" * 80)
print("Branch counting completed!")
print(
    "\nNote: Counts show valid entries (not equal to -99.0) vs total entries per branch."
)
