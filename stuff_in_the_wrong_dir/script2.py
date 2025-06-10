import ROOT

file = ROOT.TFile.Open("../run/submitDir/data-ANALYSIS/dataset.root")
file.ls()

tree = file.Get("truth_tau_analysis")
branches = [branch.GetName() for branch in tree.GetListOfBranches()]
print("Branches:", branches)

for i in range(tree.GetEntries()):
    tree.GetEntry(i)
    print(
        tree.run_number,
        tree.event_number,
        tree.phi_CP_tau_pi,
        tree.phi_CP_neutrino_pi,
    )


file.Close()