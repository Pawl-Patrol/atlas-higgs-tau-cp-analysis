from AnaAlgorithm.AlgSequence import AlgSequence
from AnalysisAlgorithmsConfig.ConfigAccumulator import ConfigAccumulator
from AnalysisAlgorithmsConfig.ConfigText import TextConfig
from Campaigns.Utils import Campaign
from AthenaConfiguration.AllConfigFlags import initConfigFlags
# Import os to be able to read environment variables
import os

# geometry: RUN{1..4}
#   https://gitlab.cern.ch/atlas/athena/-/blob/main/Control/AthenaConfiguration/python/Enums.py
# campaign enum
#   https://gitlab.cern.ch/atlas/athena/-/blob/main/Tools/Campaigns/python/Utils.py

def makeSequence (configPath, dataType='mc', isPhyslite=False, geometry='RUN2',
    campaign=Campaign.MC20e, autoconfigFromFlags=None):
    """
    Read in text configuration and produce configured algSequence
    """

    # make sure file exists
    if not os.path.isfile(configPath):
        raise FileNotFoundError(f'{configPath} cannot be found.')

    # initialize and configure configuration sequence with YAML file
    configSeq = TextConfig(configPath).configure()

    # compile
    algSeq = AlgSequence()
    configAccumulator = ConfigAccumulator(algSeq, dataType=dataType,
            isPhyslite=isPhyslite,
            geometry=geometry,
            campaign=campaign,
            autoconfigFromFlags=autoconfigFromFlags)
    configSeq.fullConfigure(configAccumulator)

    return algSeq

