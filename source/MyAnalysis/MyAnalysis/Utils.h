#ifndef MyAnalysis_Utils_H
#define MyAnalysis_Utils_H

#include "xAODTracking/TrackParticle.h"
#include <TLorentzVector.h>
#include <TVector3.h>

TVector3 calculatePionTrackVertex(const xAOD::TrackParticle *track,
                                  const TVector3 &beamSpot);

#endif