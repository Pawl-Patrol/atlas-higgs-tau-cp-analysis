#include "MyAnalysis/Observables.h"

TVector3 getPerpendicularComponent(const TVector3 &vec1, const TVector3 &vec2) {
  TVector3 unit_vec2 = vec2.Unit();
  return vec1 - (vec1.Dot(unit_vec2)) * unit_vec2;
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

double phiCP_Pion_ImpactParameter(TLorentzVector pionPosProdVtx,
                                  TLorentzVector pionNegProdVtx,
                                  TLorentzVector tauNegProdVtx,
                                  TLorentzVector pionPosP4,
                                  TLorentzVector pionNegP4) {

  // Using pion impact parameter/momentum planes
  TLorentzVector impactParamPos = TLorentzVector(
      getPerpendicularComponent(pionPosProdVtx.Vect() - tauNegProdVtx.Vect(),
                                pionPosP4.Vect())
          .Unit(),
      0.);
  TLorentzVector impactParamNeg = TLorentzVector(
      getPerpendicularComponent(pionNegProdVtx.Vect() - tauNegProdVtx.Vect(),
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
      getPerpendicularComponent(impactParamPos.Vect(), pionPosP4.Vect());
  TVector3 planeNeg =
      getPerpendicularComponent(impactParamNeg.Vect(), pionNegP4.Vect());

  double angleO = pionNegP4.Vect().Unit() * (planePos.Cross(planeNeg).Unit());
  double phi = acos(planePos.Unit() * planeNeg.Unit());

  return angleO >= 0 ? phi : 2 * M_PI - phi;
}
