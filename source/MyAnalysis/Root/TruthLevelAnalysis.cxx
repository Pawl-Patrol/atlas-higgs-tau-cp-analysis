#include "AsgMessaging/MessageCheck.h"
#include "MyAnalysis/Observables.h"
#include "TruthUtils/AtlasPID.h"
#include "xAODTau/TauJet.h"
#include "xAODTau/TauTrack.h"
#include "xAODTracking/TrackParticle.h"
#include "xAODTracking/TrackingPrimitives.h"
#include "xAODTracking/VertexContainer.h"
#include "xAODTruth/TruthParticleContainer.h"
#include "xAODTruth/TruthParticleFwd.h"
#include "xAODTruth/TruthVertexContainerFwd.h"
#include "xAODTruth/TruthVertexFwd.h"
#include "xAODTruth/versions/TruthVertex_v1.h"
#include <MyAnalysis/TruthLevelAnalysis.h>
#include <MyAnalysis/Utils.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <boost/date_time/constrained_value.hpp>
#include <cmath>
#include <xAODTau/TauJetContainer.h>
#include <xAODTruth/xAODTruthHelpers.h>

// static const int RHOPLUS = 213;
// static const int RHOMINUS = -213;
// static const int RHO0 = 113;

TruthLevelAnalysis::TruthLevelAnalysis(const std::string &name,
                                       ISvcLocator *pSvcLocator)
    : EL::AnaAlgorithm(name, pSvcLocator) {}

StatusCode TruthLevelAnalysis::initialize() {
  const int BINS = 50;

  // Setup histograms
  ANA_CHECK(book(TH1F("phi_CP_pion", "phi_CP_pion", BINS, 0, 2 * M_PI)));
  ANA_CHECK(
      book(TH1F("phi_CP_pion_jet", "phi_CP_pion_jet", BINS, 0, 2 * M_PI)));
  ANA_CHECK(book(
      TH1F("phi_CP_pion_jet_reco", "phi_CP_pion_jet_reco", BINS, 0, 2 * M_PI)));

  return StatusCode::SUCCESS;
}

StatusCode TruthLevelAnalysis::execute() {

  m_phiCPPion = NAN;
  m_phiCPPionJet = NAN;
  m_phiCPPionJetReco = NAN;

  // Event information
  const xAOD::EventInfo *eventInfo = nullptr;
  // Containers
  const xAOD::TruthParticleContainer *truthHiggsWithDecayParticles = nullptr;
  const xAOD::TruthParticleContainer *truthTausWithDecayParticles = nullptr;
  const xAOD::TauJetContainer *tauJets = nullptr;
  const xAOD::TruthVertexContainer *truthVertices = nullptr;
  const xAOD::VertexContainer *vertices = nullptr;
  // Pointers to truth particles
  const xAOD::TruthParticle *pHiggs = nullptr;
  const xAOD::TruthParticle *pTauPos = nullptr;
  const xAOD::TruthParticle *pTauNeg = nullptr;
  const xAOD::TruthParticle *pPionPos = nullptr;
  const xAOD::TruthParticle *pPionNeg = nullptr;
  const xAOD::TruthParticle *pLepton = nullptr;
  const xAOD::TruthParticle *pPosition = nullptr;
  // Pointers to tau jets
  const xAOD::TauJet *tauPosJet = nullptr;
  const xAOD::TauJet *tauNegJet = nullptr;

  ANA_CHECK(evtStore()->retrieve(eventInfo, "EventInfo"));

  StatusCode result = evtStore()->retrieve(truthHiggsWithDecayParticles,
                                           "TruthBSMWithDecayParticles");

  if (result.isFailure() || truthHiggsWithDecayParticles->empty()) {
    ANA_CHECK(evtStore()->retrieve(truthHiggsWithDecayParticles,
                                   "TruthBosonsWithDecayParticles"));
  }

  ANA_CHECK(evtStore()->retrieve(truthTausWithDecayParticles,
                                 "TruthTausWithDecayParticles"));
  ANA_CHECK(evtStore()->retrieve(tauJets, "TauJets"));
  ANA_CHECK(evtStore()->retrieve(truthVertices, "TruthPrimaryVertices"));
  ANA_CHECK(evtStore()->retrieve(vertices, "PrimaryVertices"));

  TVector3 beamSpot(eventInfo->beamPosX(), eventInfo->beamPosY(),
                    eventInfo->beamPosZ());

  int tauCount = 0;

  for (const xAOD::TruthParticle *particle : *truthHiggsWithDecayParticles) {
    if (particle->nParents() == 0) {
      continue;
    }

    if (particle->pdgId() == TAU) {
      pTauNeg = particle;
      pHiggs = particle->parent(0);
      tauCount++;
    } else if (particle->pdgId() == -TAU) {
      pTauPos = particle;
      tauCount++;
    }
  }

  // Exclusion rules
  if (pHiggs == nullptr) {
    ANA_MSG_VERBOSE("Higgs boson not found. Excluding event.");
    return StatusCode::SUCCESS;
  }

  if (pTauNeg == nullptr || pTauPos == nullptr || tauCount != 2) {
    ANA_MSG_VERBOSE(
        "Could not find tau+ tau- decay products. Excluding event.");
    return StatusCode::SUCCESS;
  }

  int tauPosLeptonCount = 0;
  int tauNegLeptonCount = 0;
  int tauPosPionZeroCount = 0;
  int tauNegPionZeroCount = 0;
  int tauPosPionChargedCount = 0;
  int tauNegPionChargedCount = 0;
  int tauPosNeutrinoCount = 0;
  int tauNegNeutrinoCount = 0;

  for (const xAOD::TruthParticle *particle : *truthTausWithDecayParticles) {
    // Skip tauons
    if (particle->nParents() == 0) {
      continue;
    }

    // Tau-
    if (particle->parent(0)->pdgId() == TAU) {
      switch (particle->pdgId()) {
      case PIPLUS:
        pPionPos = particle;
        tauNegPionChargedCount++;
        break;
      case PIMINUS:
        pPionNeg = particle;
        tauNegPionChargedCount++;
        break;
      case PI0:
        tauNegPionZeroCount++;
        break;
      case ELECTRON:
      case MUON:
        pLepton = particle;
        tauNegLeptonCount++;
        break;
      case NU_TAU:
      case -NU_TAU:
        tauNegNeutrinoCount++;
        break;
      }
    }

    // Tau+
    else if (particle->parent(0)->pdgId() == -TAU) {
      switch (particle->pdgId()) {
      case PIPLUS:
        pPionPos = particle;
        tauPosPionChargedCount++;
        break;
      case PIMINUS:
        pPionNeg = particle;
        tauPosPionChargedCount++;
        break;
      case PI0:
        tauPosPionZeroCount++;
        break;
      case POSITRON:
      case -MUON:
        pPosition = particle;
        tauPosLeptonCount++;
        break;
      case NU_TAU:
      case -NU_TAU:
        tauPosNeutrinoCount++;
        break;
      }
    }
  }

  TauDecayMode tauNegDecayMode =
      inferTauDecayMode(tauNegLeptonCount, tauNegPionChargedCount,
                        tauNegPionZeroCount, tauNegNeutrinoCount);
  TauDecayMode tauPosDecayMode =
      inferTauDecayMode(tauPosLeptonCount, tauPosPionChargedCount,
                        tauPosPionZeroCount, tauPosNeutrinoCount);

  if (tauNegDecayMode != TauDecayMode::HADRONIC_1P0N ||
      tauPosDecayMode != TauDecayMode::HADRONIC_1P0N) {
    ANA_MSG_VERBOSE(
        "Could not identify tau+ tau- decay products. Excluding event.");
    return StatusCode::SUCCESS;
  }

  for (const xAOD::TauJet *tauJet : *tauJets) {
    if (tauJet->charge() > 0) {
      tauPosJet = tauJet;
    } else if (tauJet->charge() < 0) {
      tauNegJet = tauJet;
    }
  }

  if (tauPosJet == nullptr || tauNegJet == nullptr) {
    ANA_MSG_VERBOSE("Could not find tau+ or tau- jets. Excluding event.");
    return StatusCode::SUCCESS;
  }

  if (tauPosJet->vertex() == nullptr || tauNegJet->vertex() == nullptr) {
    ANA_MSG_VERBOSE("Could not find tau+ or tau- vertex. Excluding event.");
    return StatusCode::SUCCESS;
  }

  if (tauPosJet->nTracks() != 1 || tauNegJet->nTracks() != 1) {
    ANA_MSG_VERBOSE("Tau jets do not have exactly one track. Excluding event.");
    return StatusCode::SUCCESS;
  }

  ANA_MSG_DEBUG("Found higgs -> tau+ tau- -> pi+ pi- decay.");

  // Truth level phi_CP
  m_phiCPPion = phiCP_Pion_ImpactParameter(
      pPionPos->prodVtx()->v4().Vect(), pPionNeg->prodVtx()->v4().Vect(),
      pTauPos->prodVtx()->v4().Vect(), pTauNeg->prodVtx()->v4().Vect(),
      pPionPos->p4(), pPionNeg->p4());

  // Get the tracks
  const xAOD::TrackParticle *tauPosTrack = tauPosJet->track(0)->track();
  const xAOD::TrackParticle *tauNegTrack = tauNegJet->track(0)->track();

  TVector3 pionVtxPos = calculatePionTrackVertex(tauPosTrack, beamSpot);
  TVector3 pionVtxNeg = calculatePionTrackVertex(tauNegTrack, beamSpot);

  // phi_CP using tracks and truth/reconstructed vertices
  m_phiCPPionJet = phiCP_Pion_ImpactParameter(
      pionVtxPos, pionVtxNeg, pTauPos->prodVtx()->v4().Vect(),
      pTauNeg->prodVtx()->v4().Vect(), tauPosTrack->p4(), tauNegTrack->p4());

  m_phiCPPionJetReco = phiCP_Pion_ImpactParameter(
      pionVtxPos, pionVtxNeg, GetVertexVector(tauPosJet->vertex()),
      GetVertexVector(tauNegJet->vertex()), tauPosTrack->p4(),
      tauNegTrack->p4());

  hist("phi_CP_pion")->Fill(m_phiCPPion);
  hist("phi_CP_pion_jet")->Fill(m_phiCPPionJet);
  hist("phi_CP_pion_jet_reco")->Fill(m_phiCPPionJetReco);

  return StatusCode::SUCCESS;
}

StatusCode TruthLevelAnalysis ::finalize() { return StatusCode::SUCCESS; }