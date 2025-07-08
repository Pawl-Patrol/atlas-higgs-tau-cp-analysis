#ifndef MyAnalysis_Observables_H
#define MyAnalysis_Observables_H

#include <TLorentzVector.h>

double phiCP_Pion_Tau(TLorentzVector higgsP4, TLorentzVector tauPosP4,
                      TLorentzVector tauNegP4, TLorentzVector piPosP4,
                      TLorentzVector piNegP4);

double phiCP_Pion_Neutrino(TLorentzVector higgsP4, TLorentzVector antiNeutriP4,
                           TLorentzVector neutriP4, TLorentzVector piPosP4,
                           TLorentzVector piNegP4);

double phiCP_ImpactParameter(TVector3 pionPosImpactParam,
                             TVector3 pionNegImpactParam,
                             TLorentzVector pionPosP4, TLorentzVector pionNegP4,
                             TLorentzVector referenceFrame);

double phiCP_Pion_RhoDecayPlane(TLorentzVector pionPosP4,
                                TLorentzVector pionNeuPosP4,
                                TLorentzVector pionNegP4,
                                TLorentzVector pionNeuNegP4,
                                TLorentzVector referenceFrame);

#endif