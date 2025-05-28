# See: https://twiki.cern.ch/twiki/bin/viewauth/AtlasComputing/SoftwareTutorialxAODAnalysisInCMake for more details about anything here

# User options, which can be set from command line after a "-" character
# athena ATestRun_jobOptions.py - --myOption ...
from AthenaCommon.AthArgumentParser import AthArgumentParser
athArgsParser = AthArgumentParser()
athArgsParser.add_argument( '-c', '--config-path', dest='config_path',
        action='store', default='', 
        help='path to the directory containing samples to be processed.' )
athArgs = athArgsParser.parse_args()

# make sure file exists
configPath = athArgs.config_path
if not os.path.isfile(configPath):
    raise FileNotFoundError(f'{configPath} is not a file.')

# data type: 'data', 'mc'
dataType = 'mc'

# Select the sample associated with the data type used
if dataType not in ["data", "mc"] :
    raise Exception (f"invalid data type: {dataType}")
if dataType == 'mc':
    testFile = os.getenv ('ALRB_TutorialData')+'/mc20_13TeV.312276.aMcAtNloPy8EG_A14N30NLO_LQd_mu_ld_0p3_beta_0p5_2ndG_M1000.deriv.DAOD_PHYS.e7587_e7400_a907_r14861_r14919_p6026/DAOD_PHYS.37773721._000001.pool.root.1'
else:
    testFile = os.getenv('ASG_TEST_FILE_DATA')

# Override next line on command line with: --filesInput=XXX
jps.AthenaCommonFlags.FilesInput = [testFile] 

# get flags used to configure algorithms
from AthenaConfiguration.AllConfigFlags import initConfigFlags
flags = initConfigFlags()
flags.Input.Files = [testFile]
flags.lock()

# Specify AccessMode (read mode) ... ClassAccess is good default for xAOD
jps.AthenaCommonFlags.AccessMode = "ClassAccess" 

# Create output stream
jps.AthenaCommonFlags.HistOutputs = ["ANALYSIS:MyxAODAnalysis.outputs.root"]
svcMgr.THistSvc.MaxFileSize=-1 #speeds up jobs that output lots of histograms

# Create an algSequence form the YAML configuration file


# We need to explicitly instantiate CutFlowSvc in Athena


# Create the algorithm's configuration.


# Add our algorithm to the main alg sequence


# Add all algorithms from the sequence to the job.


# Limit the number of events (for testing purposes)
theApp.EvtMax = 500

# Optional include for reducing printout from athena
include("AthAnalysisBaseComps/SuppressLogging.py")
