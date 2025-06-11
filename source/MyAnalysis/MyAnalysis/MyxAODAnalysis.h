#ifndef MyAnalysis_MyxAODAnalysis_H
#define MyAnalysis_MyxAODAnalysis_H

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

class MyxAODAnalysis : public EL::AnaAlgorithm {
private:
  unsigned int m_runNumber = 0;
  unsigned long long m_eventNumber = 0;
  double m_phiCP = 0.;
  double m_phiCPNeutri = 0.;
  double m_phiCPPion = 0.;

public:
  // This is a standard algorithm constructor
  MyxAODAnalysis(const std::string &name, ISvcLocator *pSvcLocator);

  // These are the functions inherited from Algorithm
  virtual StatusCode initialize() override;
  virtual StatusCode execute() override;
  virtual StatusCode finalize() override;

private:
  // Configuration, and any other types of variables go here.
};

#endif
