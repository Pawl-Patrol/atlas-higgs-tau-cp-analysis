#ifndef MyAnalysis_Utils_H
#define MyAnalysis_Utils_H

#include "xAODEgamma/ElectronContainer.h"
#include "xAODTau/TauJetContainer.h"
#include "xAODTracking/TrackParticle.h"
#include "xAODTracking/Vertex.h"
#include <TVector3.h>

TVector3 getTransverseComponent(const TVector3 &vec1, const TVector3 &vec2);
TVector3 getPerpendicularComponent(const TVector3 &vec1, const TVector3 &vec2);

TVector3 calculateImpactParameter(const TVector3 &trackVtx,
                                  const TVector3 &trackDirection,
                                  const TVector3 &primaryVertex);
TVector3 calculateTrackImpactParameter(const xAOD::TrackParticle *track,
                                       const TVector3 &primaryVertex,
                                       const TVector3 &beamSpot);

enum TauDecayMode {
  LEPTONIC,
  HADRONIC_1P0N,
  HADRONIC_1P1N,
  HADRONIC_1PXN,
  HADRONIC_3P0N,
  UNKNOWN
};

TauDecayMode inferTauDecayMode(int leptonCount, int pionChargedCount,
                               int pionZeroCount, int neutrinoCount);

TVector3 GetVertexVector(const xAOD::Vertex *vertex);

const xAOD::TauJet *GetLeadingJet(const xAOD::TauJetContainer *jets,
                                  bool positive);
const xAOD::Electron *
GetLeadingElectron(const xAOD::ElectronContainer *electrons, bool positive);

#endif