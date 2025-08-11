#!/usr/bin/env python3

from Run_script import SAMPLES
import ROOT
import os

print("Counting events in all samples...")
print("=" * 50)

# Open all sample files and count events
for sample_name, sample_dir in SAMPLES.items():
    sample_path = f"/samples/{sample_dir}"

    # Check if sample directory exists
    if not os.path.exists(sample_path):
        print(f"{sample_name:<25}: Directory not found at {sample_path}")
        continue

    # Find all ROOT files in the sample directory
    root_files = [f for f in os.listdir(sample_path) if f.endswith(".root.1")]

    if not root_files:
        print(f"{sample_name:<25}: No ROOT files found in {sample_path}")
        continue

    total_events = 0
    file_count = 0

    for root_filename in root_files:
        file_path = os.path.join(sample_path, root_filename)

        try:
            # Open the ROOT file
            root_file = ROOT.TFile.Open(file_path)

            if not root_file or root_file.IsZombie():
                print(f"  Warning: Failed to open {root_filename}")
                continue

            # Get the tree (try common tree names for ATLAS data)
            tree = None
            tree_names = ["CollectionTree", "tau_analysis", "nominal"]

            for tree_name in tree_names:
                tree = root_file.Get(tree_name)
                if tree:
                    break

            if not tree:
                print(f"  Warning: No recognized tree found in {root_filename}")
                root_file.Close()
                continue

            # Count events in this file
            events_in_file = tree.GetEntries()
            total_events += events_in_file
            file_count += 1

            # Close the file
            root_file.Close()

        except Exception as e:
            print(f"  Error processing {root_filename}: {str(e)}")

    if file_count > 0:
        print(f"{sample_name:<25}: {total_events:>10,} events ({file_count} files)")
    else:
        print(f"{sample_name:<25}: No valid files processed")

print("=" * 50)
print("Event counting completed!")
