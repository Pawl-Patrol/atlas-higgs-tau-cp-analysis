#include "xAODEgamma/ElectronContainer.h"
#include "xAODTau/TauJet.h"
#include "xAODTau/TauJetContainer.h"
#include <MyAnalysis/Utils.h>
#include <TVector3.h>
#include <cstdlib>

TVector3 getParallelComponent(const TVector3 &vec1, const TVector3 &vec2) {
  return vec1.Dot(vec2) / vec2.Mag2() * vec2;
}

TVector3 getPerpendicularComponent(const TVector3 &vec1, const TVector3 &vec2) {
  return vec1 - getParallelComponent(vec1, vec2);
}

TVector3 calculateImpactParameter(const TVector3 &trackVtx,
                                  const TVector3 &trackDirection,
                                  const TVector3 &primaryVertex) {
  return getPerpendicularComponent(trackVtx - primaryVertex, trackDirection);
}

TVector3 calculateTrackImpactParameter(const xAOD::TrackParticle *track,
                                       const TVector3 &primaryVertex) {
  TVector3 trackPos = TVector3(-track->d0() * sin(track->phi0()),
                               track->d0() * cos(track->phi0()), track->z0());
  return calculateImpactParameter(trackPos, track->p4().Vect(), primaryVertex);
}

TauDecayMode inferTauDecayMode(int nLepton, int nPionCharged, int nPionZero,
                               int nNeutrino) {
  if (nLepton == 1 && nPionCharged == 0 && nPionZero == 0 && nNeutrino == 2) {
    return LEPTONIC;
  }

  if (nLepton == 0 && nPionCharged == 1 && nPionZero == 0 && nNeutrino == 1) {
    return HADRONIC_1P0N;
  }

  if (nLepton == 0 && nPionCharged == 1 && nPionZero == 1 && nNeutrino == 1) {
    return HADRONIC_1P1N;
  }

  if (nLepton == 0 && nPionCharged == 1 && nPionZero >= 2 && nNeutrino == 1) {
    return HADRONIC_1PXN;
  }

  if (nLepton == 0 && nPionCharged == 3 && nPionZero == 0 && nNeutrino == 0) {
    return HADRONIC_3P0N;
  }

  return UNKNOWN;
}

TVector3 GetVertexVector(const xAOD::Vertex *vertex) {
  return TVector3(vertex->x(), vertex->y(), vertex->z());
}

const xAOD::TauJet *GetLeadingJet(const xAOD::TauJetContainer *jets,
                                  bool positive) {
  const xAOD::TauJet *leadingJet = nullptr;
  for (const xAOD::TauJet *jet : *jets) {
    if ((jet->charge() <= 0 && positive) || (jet->charge() >= 0 && !positive)) {
      continue;
    }

    // Require 1 or 3 tracks
    if (jet->nTracks() != 1 && jet->nTracks() != 3) {
      continue;
    }

    // Require leading track
    if (jet->track(0) == nullptr || jet->track(0)->track() == nullptr) {
      continue;
    }

    if (jet->pt() < 20000.0) {
      continue; // Skip jets below 25 GeV
    }

    if (abs(jet->eta()) > 2.47 ||
        (abs(jet->eta()) > 1.37 && abs(jet->eta()) < 1.52)) {
      continue; // Skip jets outside the acceptance range
    }

    if (leadingJet == nullptr || jet->pt() > leadingJet->pt()) {
      leadingJet = jet;
    }
  }

  return leadingJet;
}

const xAOD::Electron *
GetLeadingElectron(const xAOD::ElectronContainer *electrons, bool positive) {
  const xAOD::Electron *leading = nullptr;
  for (const xAOD::Electron *electron : *electrons) {
    if ((electron->charge() <= 0 && positive) ||
        (electron->charge() >= 0 && !positive)) {
      continue;
    }

    if (electron->nTrackParticles() == 0) {
      continue; // Skip electrons without tracks
    }

    if (leading == nullptr || electron->pt() > leading->pt()) {
      leading = electron;
    }
  }

  return leading;
}