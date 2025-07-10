#include "MyAnalysis/Observables.h"
#include "MyAnalysis/Utils.h"
#include <TLorentzVector.h>
#include <TVector3.h>

double phiCP_ImpactParameter(TVector3 pionPosImpactParam,
                             TVector3 pionNegImpactParam,
                             TLorentzVector pionPosP4, TLorentzVector pionNegP4,
                             TLorentzVector referenceFrame) {

  // Using pion impact parameter/momentum planes
  TLorentzVector impactParamPos(pionPosImpactParam.Unit(), 0.0);
  TLorentzVector impactParamNeg(pionNegImpactParam.Unit(), 0.0);

  // Boost into CMF of the pions
  TVector3 cmfBoostVector = referenceFrame.BoostVector();
  impactParamPos.Boost(-cmfBoostVector);
  impactParamNeg.Boost(-cmfBoostVector);
  pionPosP4.Boost(-cmfBoostVector);
  pionNegP4.Boost(-cmfBoostVector);

  // Get the impact parameter component perpendicular to the momentum
  TVector3 planePos =
      getPerpendicularComponent(impactParamPos.Vect(), pionPosP4.Vect()).Unit();
  TVector3 planeNeg =
      getPerpendicularComponent(impactParamNeg.Vect(), pionNegP4.Vect()).Unit();

  double angleO = pionNegP4.Vect().Unit().Dot(planePos.Cross(planeNeg));
  double phi = acos(planePos * planeNeg);

  return angleO >= 0 ? phi : 2 * M_PI - phi;
}

double phiCP_Pion_RhoDecayPlane(TLorentzVector pionPosP4,
                                TLorentzVector pionNeuPosP4,
                                TLorentzVector pionNegP4,
                                TLorentzVector pionNeuNegP4,
                                TLorentzVector referenceFrame) {
  if (!pionPosP4.Mag2() || !pionNegP4.Mag2() || !pionNeuPosP4.Mag2() ||
      !pionNeuNegP4.Mag2()) {
    return -99; // Invalid input
  }

  double yPos =
      (pionPosP4.E() - pionNeuPosP4.E()) / (pionPosP4.E() + pionNeuPosP4.E());
  double yNeg =
      (pionNegP4.E() - pionNeuNegP4.E()) / (pionNegP4.E() + pionNeuNegP4.E());

  TVector3 cmfBoostVecotr = referenceFrame.BoostVector();
  pionPosP4.Boost(-cmfBoostVecotr);
  pionNeuPosP4.Boost(-cmfBoostVecotr);
  pionNegP4.Boost(-cmfBoostVecotr);
  pionNeuNegP4.Boost(-cmfBoostVecotr);

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

double phiCP_IP_Rho(TVector3 pionImpactParam, TLorentzVector pionP4,
                    TLorentzVector rhoChargedP4, TLorentzVector rhoNeutralP4,
                    TLorentzVector referenceFrame, bool rhoIsPositive) {
  if (!pionP4.Mag2() || !rhoNeutralP4.Mag2() || !rhoNeutralP4.Mag2()) {
    return -99; // Invalid input
  }

  double y = (rhoChargedP4.E() - rhoNeutralP4.E()) /
             (rhoChargedP4.E() + rhoNeutralP4.E());

  TLorentzVector impactParam(pionImpactParam.Unit(), 0.0);

  TVector3 cmfBoostVector = referenceFrame.BoostVector();
  impactParam.Boost(-cmfBoostVector);
  pionP4.Boost(-cmfBoostVector);
  rhoChargedP4.Boost(-cmfBoostVector);
  rhoNeutralP4.Boost(-cmfBoostVector);

  TVector3 planeIP =
      getPerpendicularComponent(impactParam.Vect(), pionP4.Vect()).Unit();
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