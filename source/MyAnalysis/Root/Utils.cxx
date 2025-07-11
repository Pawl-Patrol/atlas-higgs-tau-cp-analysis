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
  TVector3 pointOfClosestApproach =
      TVector3(-track->d0() * sin(track->phi0()),
               track->d0() * cos(track->phi0()), track->z0());
  return calculateImpactParameter(pointOfClosestApproach, track->p4().Vect(),
                                  primaryVertex);
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

/**
 * From the TauJet container, find the best candidate jet for the given charge.
 * Applies some very basic selection criteria and find the jet with the highest
 * pT.
 */
const xAOD::TauJet *GetLeadingJet(const xAOD::TauJetContainer *jets,
                                  bool positive) {
  const xAOD::TauJet *leadingJet = nullptr;
  for (const xAOD::TauJet *jet : *jets) {
    // Require the jet to have the right charge
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

    // Skip jets below 20 GeV
    if (jet->pt() < 20000.0) {
      continue;
    }

    // Skip jets outside the barrel and endcap regions
    if (abs(jet->eta()) > 2.47 ||
        (abs(jet->eta()) > 1.37 && abs(jet->eta()) < 1.52)) {
      continue;
    }

    if (leadingJet == nullptr || jet->pt() > leadingJet->pt()) {
      leadingJet = jet;
    }
  }

  return leadingJet;
}

/**
 * From the Electron container, find the leading electron with the given charge.
 * Applies some very basic selection criteria and finds the electron with the
 * highest pT.
 */
const xAOD::Electron *
GetLeadingElectron(const xAOD::ElectronContainer *electrons, bool positive) {
  const xAOD::Electron *leadingElectron = nullptr;
  for (const xAOD::Electron *electron : *electrons) {
    // Require the electron to have the right charge
    if ((electron->charge() <= 0 && positive) ||
        (electron->charge() >= 0 && !positive)) {
      continue;
    }

    // Skip electrons without tracks
    if (electron->nTrackParticles() == 0) {
      continue;
    }

    if (leadingElectron == nullptr || electron->pt() > leadingElectron->pt()) {
      leadingElectron = electron;
    }
  }

  return leadingElectron;
}