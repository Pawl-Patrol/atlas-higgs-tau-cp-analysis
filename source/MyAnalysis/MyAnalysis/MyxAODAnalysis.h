#ifndef MyAnalysis_MyxAODAnalysis_H
#define MyAnalysis_MyxAODAnalysis_H

#include <AnaAlgorithm/AnaAlgorithm.h>
#include <TH1.h>

class MyxAODAnalysis : public EL::AnaAlgorithm {
private:
  unsigned int m_runNumber = 0;
  unsigned long long m_eventNumber = 0;
  double m_phi_CP = 0.;
  double m_phi_CP_neutri = 0.;

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
