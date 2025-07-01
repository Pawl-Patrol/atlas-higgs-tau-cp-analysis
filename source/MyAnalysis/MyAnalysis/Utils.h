#ifndef MyAnalysis_Utils_H
#define MyAnalysis_Utils_H

#include "xAODTracking/TrackParticle.h"
#include "xAODTracking/Vertex.h"
#include <TVector3.h>

TVector3 calculatePionTrackVertex(const xAOD::TrackParticle *track,
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

#endif