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
  double m_phiCP_lept_1p0n_truth = 0.;
  double m_phiCP_lept_1p0n_recon = 0.;
  double m_phiCP_lept_1p1n_truth = 0.;
  double m_phiCP_lept_1p1n_recon = 0.;
  double m_phiCP_lept_1pXn_truth = 0.;
  double m_phiCP_lept_1pXn_recon = 0.;

  double m_phiCP_1p0n_1p0n_truth = 0.;
  double m_phiCP_1p0n_1p0n_recon = 0.;

  double m_phiCP_1p0n_1p0n_prod_vtx_diff_x = 0.;
  double m_phiCP_1p0n_1p0n_prod_vtx_diff_y = 0.;
  double m_phiCP_1p0n_1p0n_prod_vtx_diff_z = 0.;
  double m_phiCP_1p0n_1p0n_prod_vtx_diff_abs = 0.;

  double m_phiCP_1p0n_1p0n_prim_vtx_diff_x = 0.;
  double m_phiCP_1p0n_1p0n_prim_vtx_diff_y = 0.;
  double m_phiCP_1p0n_1p0n_prim_vtx_diff_z = 0.;
  double m_phiCP_1p0n_1p0n_prim_vtx_diff_abs = 0.;

  double m_phiCP_1p1n_1p1n_truth = 0.;
  double m_phiCP_1p1n_1p1n_recon = 0.;
  double m_phiCP_1p1n_1pXn_truth = 0.;
  double m_phiCP_1p1n_1pXn_recon = 0.;

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
