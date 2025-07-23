# Overview

Contains code and scripts used for creating plots of CP-sensitive observables in Higgs -> tau tau decays at ATLAS. Based on [this](https://atlassoftwaredocs.web.cern.ch/analysis-software/AnalysisSWTutorial/).

# Usage

Commands I use regulary.

## Commands on startup

1. `setupATLAS -c el9 --mount /path/to/samples:/samples`
2. `asetup --restore`
3. `source x86_64-el9-gcc11-opt/setup.sh`

As a single command:

`setupATLAS -c el9 --mount /path/to/samples:/samples --postsetup="cd build && asetup --restore && source x86_64-el9-gcc11-opt/setup.sh && cd .."`

## Compilation

1. `cd ./build`
2. `cmake -DCMAKE_EXPORT_COMPILE_COMMANDS=1 ../source`
3. `make`

## Running
- `Run_script.py` - Run algorithm on samples
- `Plot_script.py` - Plot histograms on ntuples

## Other useful commands
- `checkxAOD.py ./sample.root` - List information about ROOT file
- `root -l ./sample.root` - Start interactive session for ROOT file