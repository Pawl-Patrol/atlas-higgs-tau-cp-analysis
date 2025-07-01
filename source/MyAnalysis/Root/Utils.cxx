#include <MyAnalysis/Utils.h>

TVector3 calculatePionTrackVertex(const xAOD::TrackParticle *track,
                                  const TVector3 &beamSpot) {
  return beamSpot + TVector3(-(track->d0()) * sin(track->phi0()),
                             (track->d0()) * cos(track->phi0()), track->z0());
}

TauDecayMode inferTauDecayMode(int nLepton, int nPionCharged, int nPionZero,
                               int nNeutrino) {
  if (nLepton == 1 && nPionCharged == 0 && nPionZero == 0 && nNeutrino == 2) {
    return TauDecayMode::LEPTONIC;
  }

  if (nLepton == 0 && nPionCharged == 1 && nPionZero == 0 && nNeutrino == 1) {
    return TauDecayMode::HADRONIC_1P0N;
  }

  if (nLepton == 0 && nPionCharged == 1 && nPionZero == 1 && nNeutrino == 1) {
    return TauDecayMode::HADRONIC_1P1N;
  }

  if (nLepton == 0 && nPionCharged == 1 && nPionZero >= 2 && nNeutrino == 1) {
    return TauDecayMode::HADRONIC_1PXN;
  }

  if (nLepton == 0 && nPionCharged == 3 && nPionZero == 0 && nNeutrino == 0) {
    return TauDecayMode::HADRONIC_3P0N;
  }

  return TauDecayMode::UNKNOWN;
}

TVector3 GetVertexVector(const xAOD::Vertex *vertex) {
  return TVector3(vertex->x(), vertex->y(), vertex->z());
}