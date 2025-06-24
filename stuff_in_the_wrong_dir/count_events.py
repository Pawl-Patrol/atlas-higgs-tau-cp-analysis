import ROOT
import os

# Directory containing ROOT files
root_dir = "../samples/cp-odd/"

# Name of the TTree to count events in
tree_name = "CollectionTree"

total_events = 0

for filename in os.listdir(root_dir):
    if filename.endswith(".root.1"):
        file_path = os.path.join(root_dir, filename)
        f = ROOT.TFile.Open(file_path)
        if not f or f.IsZombie():
            print(f"Could not open {file_path}")
            continue
        tree = f.Get(tree_name)
        if not tree:
            print(f"No tree named '{tree_name}' in {file_path}")
            f.Close()
            continue
        n_events = tree.GetEntries()
        print(f"{filename}: {n_events} events")
        total_events += n_events
        f.Close()

print(f"Total events in directory '{root_dir}': {total_events}")
