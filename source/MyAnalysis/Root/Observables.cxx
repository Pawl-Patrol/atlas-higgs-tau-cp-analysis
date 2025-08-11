#include "MyAnalysis/Observables.h"
#include "MyAnalysis/Utils.h"
#include <TLorentzVector.h>
#include <TVector3.h>

/* IP-method */
double phiCP_ImpactParameter(TVector3 pionPosImpactParam,
                             TVector3 pionNegImpactParam,
                             TLorentzVector pionPosP4, TLorentzVector pionNegP4,
                             TLorentzVector referenceFrame) {
  TLorentzVector pionPosImpactParamP4(pionPosImpactParam.Unit(), 0.0);
  TLorentzVector pionNegImpactParamP4(pionNegImpactParam.Unit(), 0.0);

  // Boost into reference frame
  TVector3 boostVector = referenceFrame.BoostVector();
  pionPosImpactParamP4.Boost(-boostVector);
  pionNegImpactParamP4.Boost(-boostVector);
  pionPosP4.Boost(-boostVector);
  pionNegP4.Boost(-boostVector);

  // Get the impact parameter component perpendicular to the momentum
  TVector3 planePos =
      getPerpendicularComponent(pionPosImpactParamP4.Vect(), pionPosP4.Vect())
          .Unit();
  TVector3 planeNeg =
      getPerpendicularComponent(pionNegImpactParamP4.Vect(), pionNegP4.Vect())
          .Unit();

  double angleO = pionNegP4.Vect().Unit().Dot(planePos.Cross(planeNeg));
  double phi = acos(planePos * planeNeg);

  return angleO >= 0 ? phi : 2 * M_PI - phi;
}

/* ρ-method */
double phiCP_Pion_RhoDecayPlane(TLorentzVector pionPosP4,
                                TLorentzVector pionNeuPosP4,
                                TLorentzVector pionNegP4,
                                TLorentzVector pionNeuNegP4,
                                TLorentzVector referenceFrame) {
  // Calculate y+ and y- in the laboratory frame
  double yPos = upsilon(pionPosP4.E(), pionNeuPosP4.E());
  double yNeg = upsilon(pionNegP4.E(), pionNeuNegP4.E());

  // Boost into reference frame
  TVector3 boostVector = referenceFrame.BoostVector();
  pionPosP4.Boost(-boostVector);
  pionNeuPosP4.Boost(-boostVector);
  pionNegP4.Boost(-boostVector);
  pionNeuNegP4.Boost(-boostVector);

  // Get the neutral p4 component perpendicular to the momentum
  TVector3 planePos =
      getPerpendicularComponent(pionNeuPosP4.Vect(), pionPosP4.Vect()).Unit();
  TVector3 planeNeg =
      getPerpendicularComponent(pionNeuNegP4.Vect(), pionNegP4.Vect()).Unit();

  double angleO = pionNegP4.Vect().Unit().Dot(planePos.Cross(planeNeg));
  double phiStar = acos(planePos * planeNeg);

  double phiStarPrime = angleO >= 0 ? phiStar : 2 * M_PI - phiStar;

  if (yPos * yNeg < 0) {
    if (phiStarPrime < M_PI) {
      return phiStarPrime + M_PI;
    } else {
      return phiStarPrime - M_PI;
    }
  }

  return phiStarPrime;
}

/* IP-ρ-method */
double phiCP_IP_Rho(TVector3 pionImpactParam, TLorentzVector pionP4,
                    TLorentzVector rhoChargedP4, TLorentzVector rhoNeutralP4,
                    TLorentzVector referenceFrame, bool rhoIsPositive) {
  TLorentzVector pionImpactParamP4(pionImpactParam.Unit(), 0.0);

  // Calculate y in the laboratory frame
  double y = (rhoChargedP4.E() - rhoNeutralP4.E()) /
             (rhoChargedP4.E() + rhoNeutralP4.E());

  // Boost into reference frame
  TVector3 cmfBoostVector = referenceFrame.BoostVector();
  pionImpactParamP4.Boost(-cmfBoostVector);
  pionP4.Boost(-cmfBoostVector);
  rhoChargedP4.Boost(-cmfBoostVector);
  rhoNeutralP4.Boost(-cmfBoostVector);

  TVector3 planeIP =
      getPerpendicularComponent(pionImpactParamP4.Vect(), pionP4.Vect()).Unit();
  TVector3 planeRho =
      getPerpendicularComponent(rhoNeutralP4.Vect(), rhoChargedP4.Vect())
          .Unit();

  double angleO;
  if (rhoIsPositive) {
    angleO = pionP4.Vect().Unit().Dot(planeRho.Cross(planeIP));
  } else {
    angleO = rhoChargedP4.Vect().Unit().Dot(planeIP.Cross(planeRho));
  }

  double phiStar = acos(planeRho * planeIP);

  double phiStarPrime = angleO >= 0 ? phiStar : 2 * M_PI - phiStar;

  if (y < 0) {
    if (phiStarPrime < M_PI) {
      return phiStarPrime + M_PI;
    } else {
      return phiStarPrime - M_PI;
    }
  }

  return phiStarPrime;
}