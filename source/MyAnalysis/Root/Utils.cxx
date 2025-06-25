#include <MyAnalysis/Utils.h>

TVector3 calculatePionTrackVertex(const xAOD::TrackParticle *track,
                                  const TVector3 &beamSpot) {
  return beamSpot + TVector3(-track->d0() * sin(track->phi0()),
                             track->d0() * cos(track->phi0()), track->z0());
}