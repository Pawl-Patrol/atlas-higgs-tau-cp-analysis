#ifndef MyAnalysis_DetectorLevelAnalysis_H
#define MyAnalysis_DetectorLevelAnalysis_H

#include "AsgMessaging/StatusCode.h"
#include <AnaAlgorithm/AnaAlgorithm.h>
#include <AsgMessaging/MessageCheck.h>
#include <Math/GenVector/LorentzVector.h>
#include <TH1.h>
#include <TTree.h>
#include <TruthUtils/AtlasPID.h>
#include <cmath>
#include <xAODEventInfo/EventInfo.h>
#include <xAODTruth/TruthEventContainer.h>
#include <xAODTruth/TruthParticle.h>
#include <xAODTruth/TruthParticleContainer.h>

class DetectorLevelAnalysis : public EL::AnaAlgorithm {
private:
  double m_phiCP = 0.;

public:
  // This is a standard algorithm constructor
  DetectorLevelAnalysis(const std::string &name, ISvcLocator *pSvcLocator);

  // These are the functions inherited from Algorithm
  virtual StatusCode initialize() override;
  virtual StatusCode execute() override;
  virtual StatusCode finalize() override;

private:
  // Configuration, and any other types of variables go here.
};

#endif
