#ifndef MyAnalysis_Observables_H
#define MyAnalysis_Observables_H

#include <TLorentzVector.h>

double phiCP_Pion_Tau(TLorentzVector higgsP4, TLorentzVector tauPosP4,
                      TLorentzVector tauNegP4, TLorentzVector piPosP4,
                      TLorentzVector piNegP4);

double phiCP_Pion_Neutrino(TLorentzVector higgsP4, TLorentzVector antiNeutriP4,
                           TLorentzVector neutriP4, TLorentzVector piPosP4,
                           TLorentzVector piNegP4);

double phiCP_Pion_ImpactParameter(TLorentzVector pionPosProdVtx,
                                  TLorentzVector pionNegProdVtx,
                                  TLorentzVector tauNegProdVtx,
                                  TLorentzVector pionPosP4,
                                  TLorentzVector pionNegP4);

#endif