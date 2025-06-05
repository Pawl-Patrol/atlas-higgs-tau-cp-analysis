import ROOT

file = ROOT.TFile.Open("./run/submitDir/hist-dataset.root")
hist = file.Get("eventNumber")

canvas = ROOT.TCanvas()
hist.Draw()
canvas.SaveAs("histogram.png")

for i in range(1, hist.GetNbinsX() + 1):
    print(f"Bin {i}: {hist.GetBinContent(i)}")


file.Close()