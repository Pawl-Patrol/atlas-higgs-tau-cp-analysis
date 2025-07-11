#ifndef MyAnalysis_TruthLevelAnalysis_H
#define MyAnalysis_TruthLevelAnalysis_H

#include <AnaAlgorithm/AnaAlgorithm.h>

class TruthLevelAnalysis : public EL::AnaAlgorithm {
public:
  TruthLevelAnalysis(const std::string &name, ISvcLocator *pSvcLocator);

  virtual StatusCode initialize() override;
  virtual StatusCode execute() override;
  virtual StatusCode finalize() override;

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
};

#endif
