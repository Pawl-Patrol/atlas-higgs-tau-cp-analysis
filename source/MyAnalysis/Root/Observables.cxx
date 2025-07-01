#include "MyAnalysis/Observables.h"
#include <TLorentzVector.h>
#include <TVector3.h>

TVector3 getTransverseComponent(const TVector3 &vec1, const TVector3 &vec2) {
  TVector3 unit_vec2 = vec2.Unit();
  return vec1.Dot(unit_vec2) * unit_vec2;
}

TVector3 getPerpendicularComponent(const TVector3 &vec1, const TVector3 &vec2) {
  return vec1 - getTransverseComponent(vec1, vec2);
}

double phiCP_Pion_Tau(TLorentzVector higgsP4, TLorentzVector tauPosP4,
                      TLorentzVector tauNegP4, TLorentzVector piPosP4,
                      TLorentzVector piNegP4) {

  TVector3 higgsBoostVector = higgsP4.BoostVector();
  tauPosP4.Boost(-higgsBoostVector);
  tauNegP4.Boost(-higgsBoostVector);
  piPosP4.Boost(-higgsBoostVector);
  piNegP4.Boost(-higgsBoostVector);

  TVector3 planePos = (tauPosP4.Vect().Cross(piPosP4.Vect())).Unit();
  TVector3 planeNeg = (tauNegP4.Vect().Cross(piNegP4.Vect())).Unit();

  double angleO = piNegP4.Vect().Unit() * (planePos.Cross(planeNeg));
  double phi = acos(planePos * planeNeg);

  return angleO >= 0 ? phi : 2 * M_PI - phi;
}

double phiCP_Pion_Neutrino(TLorentzVector higgsP4, TLorentzVector antiNeutriP4,
                           TLorentzVector neutriP4, TLorentzVector piPosP4,
                           TLorentzVector piNegP4) {
  TVector3 higgsBoostVector = higgsP4.BoostVector();
  antiNeutriP4.Boost(-higgsBoostVector);
  neutriP4.Boost(-higgsBoostVector);
  piPosP4.Boost(-higgsBoostVector);
  piNegP4.Boost(-higgsBoostVector);

  TVector3 planePos = (antiNeutriP4.Vect().Cross(piPosP4.Vect())).Unit();
  TVector3 planeNeg = (neutriP4.Vect().Cross(piNegP4.Vect())).Unit();

  double angleO = piNegP4.Vect().Unit() * (planePos.Cross(planeNeg));
  double phi = acos(planePos * planeNeg);

  return angleO >= 0 ? phi : 2 * M_PI - phi;
}

double
phiCP_Pion_ImpactParameter(TVector3 pionPosProdVtx, TVector3 pionNegProdVtx,
                           TVector3 tauNegProdVtx, TVector3 tauPosProdVtx,
                           TLorentzVector pionPosP4, TLorentzVector pionNegP4) {

  // Using pion impact parameter/momentum planes
  TLorentzVector impactParamPos =
      TLorentzVector(getPerpendicularComponent(pionPosProdVtx - tauPosProdVtx,
                                               pionPosP4.Vect())
                         .Unit(),
                     0.);
  TLorentzVector impactParamNeg =
      TLorentzVector(getPerpendicularComponent(pionNegProdVtx - tauNegProdVtx,
                                               pionNegP4.Vect())
                         .Unit(),
                     0.);

  // Boost into CMF of the pions
  TVector3 cmfBoostVector = (pionPosP4 + pionNegP4).BoostVector();
  impactParamPos.Boost(-cmfBoostVector);
  impactParamNeg.Boost(-cmfBoostVector);
  pionPosP4.Boost(-cmfBoostVector);
  pionNegP4.Boost(-cmfBoostVector);

  // Get the impact parameter component perpendicular to the momentum
  TVector3 planePos =
      getPerpendicularComponent(impactParamPos.Vect(), pionPosP4.Vect()).Unit();
  TVector3 planeNeg =
      getPerpendicularComponent(impactParamNeg.Vect(), pionNegP4.Vect()).Unit();

  double angleO = pionNegP4.Vect().Unit() * (planePos.Cross(planeNeg));
  double phi = acos(planePos * planeNeg);

  return angleO >= 0 ? phi : 2 * M_PI - phi;
}

double phiCP_Pion_RhoDecayPlane(TLorentzVector pionPosP4,
                                TLorentzVector pionNeuPosP4,
                                TLorentzVector pionNegP4,
                                TLorentzVector pionNeuNegP4) {
  TVector3 cmfBoostVecotr =
      (pionPosP4 + pionNeuPosP4 + pionNegP4 + pionNeuNegP4).BoostVector();
  pionPosP4.Boost(-cmfBoostVecotr);
  pionNeuPosP4.Boost(-cmfBoostVecotr);
  pionNegP4.Boost(-cmfBoostVecotr);
  pionNeuNegP4.Boost(-cmfBoostVecotr);

  TVector3 planePos =
      getTransverseComponent(pionNeuPosP4.Vect(), pionPosP4.Vect()).Unit();
  TVector3 planeNeg =
      getTransverseComponent(pionNeuNegP4.Vect(), pionNegP4.Vect()).Unit();

  double angleO = pionNegP4.Vect().Unit() * (planePos.Cross(planeNeg));
  double phiStar = acos(planePos * planeNeg);

  double phiStarPrime = angleO >= 0 ? phiStar : 2 * M_PI - phiStar;

  double yPos =
      (pionPosP4.E() - pionNeuPosP4.E()) / (pionPosP4.E() + pionNeuPosP4.E());
  double yNeg =
      (pionNegP4.E() - pionNeuNegP4.E()) / (pionNegP4.E() + pionNeuNegP4.E());

  return yPos * yNeg >= 0 ? phiStarPrime : phiStarPrime + M_PI;
}
