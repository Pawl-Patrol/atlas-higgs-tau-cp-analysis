#!/usr/bin/env python3

# Import os to be able to read environment variables
import os

# Read the submission directory as a command line argument. You can
# extend the list of arguments with your private ones later on.
import argparse

parser = argparse.ArgumentParser()
parser.add_argument(
    "-c",
    "--config-path",
    dest="configPath",
    action="store",
    type=str,
    required=True,
    help="Path to the directory containing samples to be processed.",
)
parser.add_argument(
    "-s",
    "--submission-dir",
    dest="submitDir",
    action="store",
    type=str,
    default="submitDir",
    help="Submission directory for EventLoop",
)
parser.add_argument(
    "-e",
    "--event-limit",
    dest="eventLimit",
    action="store",
    type=int,
    default=500,
    help="Maximum number of events to process. Use -1 for no limit.",
)
parser.add_argument(
    "-d",
    "--debug",
    dest="debug",
    action="store_true",
    default=False,
    help="Enable debug output.",
)
options = parser.parse_args()

# Set up (Py)ROOT.
import ROOT

ROOT.xAOD.Init().ignore()

# Set up the SampleHandler object to handle the input files
sh = ROOT.SH.SampleHandler()

# Set the name of the tree in our files in the xAOD the TTree
# containing the EDM containers is "CollectionTree"
sh.setMetaString("nc_tree", "CollectionTree")

# Use SampleHandler to get the sample from the defined location
sample = ROOT.SH.SampleLocal("dataset")
for filename in os.listdir(options.configPath):
    sample.add(os.path.join(options.configPath, filename))
sh.add(sample)

# Print information about the sample
sh.printContent()

# Create an EventLoop job.
job = ROOT.EL.Job()
job.sampleHandler(sh)
job.options().setDouble(ROOT.EL.Job.optMaxEvents, options.eventLimit)
job.options().setString(ROOT.EL.Job.optSubmitDirMode, "unique-link")

# Create the algorithm's configuration.
from AnaAlgorithm.DualUseConfig import createAlgorithm

alg = createAlgorithm("TruthLevelAnalysis", "AnalysisAlg")

# Set the algorithm output level.
if options.debug:
    alg.OutputLevel = ROOT.MSG.DEBUG

# Add our algorithm to the job
job.algsAdd(alg)

# Later on we'll add some configuration options for our job here

# Add output stream
job.outputAdd(ROOT.EL.OutputStream("ANALYSIS"))

# Run the job using the direct driver.
driver = ROOT.EL.DirectDriver()
driver.submit(job, options.submitDir)

# Do not add anything here!!!
