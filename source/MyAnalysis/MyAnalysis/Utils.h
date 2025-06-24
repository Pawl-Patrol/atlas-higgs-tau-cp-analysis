#ifndef MyAnalysis_Utils_H
#define MyAnalysis_Utils_H

#include "MyAnalysis/Observables.h"
#include "xAODTau/TauTrack.h"
#include "xAODTracking/TrackParticle.h"
#include <TLorentzVector.h>
#include <TVector3.h>
#include <cmath>

TVector3 constructImpactParameterFromTauTrack(const xAOD::TrackParticle *track,
                                              TVector3 primaryVtx) {
  TVector3 supportVtx = TVector3(track->d0() * cos(track->phi()),
                                 track->d0() * sin(track->phi()), track->z0());
  TVector3 direction = track->p4().Vect();
  TVector3 result =
      getPerpendicularComponent(supportVtx - primaryVtx, direction);
  return supportVtx - primaryVtx;
}

#endif