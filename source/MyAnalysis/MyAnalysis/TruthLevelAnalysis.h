#ifndef MyAnalysis_TruthLevelAnalysis_H
#define MyAnalysis_TruthLevelAnalysis_H

#include "AsgMessaging/StatusCode.h"
#include <AnaAlgorithm/AnaAlgorithm.h>
#include <AsgMessaging/MessageCheck.h>
#include <Math/GenVector/LorentzVector.h>
#include <TH1.h>
#include <TLorentzVector.h>
#include <TTree.h>
#include <TruthUtils/AtlasPID.h>
#include <cmath>
#include <xAODEventInfo/EventInfo.h>

class TruthLevelAnalysis : public EL::AnaAlgorithm {
private:
  double m_phiCP = 0.;
  double m_phiCPNeutri = 0.;
  double m_phiCPPion = 0.;
  double m_phiCPPionJet = 0.;
  double m_phiCPPionJetReco = 0.;

public:
  // This is a standard algorithm constructor
  TruthLevelAnalysis(const std::string &name, ISvcLocator *pSvcLocator);

  // These are the functions inherited from Algorithm
  virtual StatusCode initialize() override;
  virtual StatusCode execute() override;
  virtual StatusCode finalize() override;

private:
  // Configuration, and any other types of variables go here.
};

#endif
