#include "AsgMessaging/MessageCheck.h"
#include "AsgMessaging/StatusCode.h"
#include "MyAnalysis/Observables.h"
#include "TruthUtils/AtlasPID.h"
#include "xAODEgamma/ElectronContainer.h"
#include "xAODTau/TauJet.h"
#include "xAODTau/TauTrack.h"
#include "xAODTracking/TrackParticle.h"
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
  const xAOD::ElectronContainer *electrons = nullptr;
  const xAOD::TruthVertexContainer *truthVertices = nullptr;
  const xAOD::VertexContainer *vertices = nullptr;
  // Pointers to truth particles
  const xAOD::TruthParticle *pHiggs = nullptr;
  const xAOD::TruthParticle *pTauPos = nullptr;
  const xAOD::TruthParticle *pTauNeg = nullptr;
  const xAOD::TruthParticle *pPionPos = nullptr;
  const xAOD::TruthParticle *pPionNeg = nullptr;
  const xAOD::TruthParticle *pLeptonNeg = nullptr;
  const xAOD::TruthParticle *pLeptonPos = nullptr;

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
  ANA_CHECK(evtStore()->retrieve(electrons, "Electrons"));
  ANA_CHECK(evtStore()->retrieve(truthVertices, "TruthPrimaryVertices"));
  ANA_CHECK(evtStore()->retrieve(vertices, "PrimaryVertices"));

  TVector3 beamSpot(eventInfo->beamPosX(), eventInfo->beamPosY(),
                    eventInfo->beamPosZ());

  int tauCount = 0;

  for (const xAOD::TruthParticle *particle : *truthHiggsWithDecayParticles) {
    if (particle->nParents() == 0) {
      continue;
    }

    if (particle->parent(0)->pdgId() != HIGGSBOSON) {
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

    const xAOD::TruthParticle *parent = particle->parent(0);

    // Resolve tau chains which may appear
    // TODO: is this still needed?
    while (abs(parent->pdgId()) == TAU && parent->nParents() > 0 &&
           abs(parent->parent(0)->pdgId()) == TAU) {
      parent = parent->parent(0);
    }

    // Tau-
    if (parent->barcode() == pTauNeg->barcode()) {
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
      case MUON:
      case ELECTRON:
        pLeptonNeg = particle;
        tauNegLeptonCount++;
        break;
      case NU_TAU:
      case -NU_TAU:
      case NU_E:
      case -NU_E:
      case NU_MU:
      case -NU_MU:
        tauNegNeutrinoCount++;
        break;
      }
    }

    // Tau+
    else if (parent->barcode() == pTauPos->barcode()) {
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
      case -MUON:
      case POSITRON:
        pLeptonPos = particle;
        tauPosLeptonCount++;
        break;
      case NU_TAU:
      case -NU_TAU:
      case NU_E:
      case -NU_E:
      case NU_MU:
      case -NU_MU:
        tauPosNeutrinoCount++;
        break;
      }
    }
  }

  TVector3 primaryVertex(0, 0, 0);
  bool foundPrimaryVertex = false;

  for (const xAOD::Vertex *vertex : *vertices) {
    if (vertex->vertexType() == xAOD::VxType::PriVtx) {
      primaryVertex = GetVertexVector(vertex);
      foundPrimaryVertex = true;
    }
  }

  if (!foundPrimaryVertex) {
    ANA_MSG_VERBOSE("No primary vertex found. Excluding event.");
    return StatusCode::SUCCESS;
  }

  TauDecayMode tauNegDecayMode =
      inferTauDecayMode(tauNegLeptonCount, tauNegPionChargedCount,
                        tauNegPionZeroCount, tauNegNeutrinoCount);
  TauDecayMode tauPosDecayMode =
      inferTauDecayMode(tauPosLeptonCount, tauPosPionChargedCount,
                        tauPosPionZeroCount, tauPosNeutrinoCount);

  if (tauNegDecayMode == TauDecayMode::HADRONIC_1P0N &&
      tauPosDecayMode == TauDecayMode::HADRONIC_1P0N) {
    // if (tauJets->size() != 2) {
    //   ANA_MSG_VERBOSE("TauJets size is not 2, but " << tauJets->size()
    //                                                 << ". Excluding event.");
    //   return StatusCode::SUCCESS;
    // }

    const xAOD::TauJet *tauPosJet = nullptr;
    const xAOD::TauJet *tauNegJet = nullptr;

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

    ANA_MSG_DEBUG("Found higgs -> tau+ tau- -> pion+ pion- decay");

    TVector3 imParamPos = calculateImpactParameter(
        pPionPos->prodVtx()->v4().Vect(), pPionPos->p4().Vect(),
        pTauPos->prodVtx()->v4().Vect());
    TVector3 imParamNeg = calculateImpactParameter(
        pPionNeg->prodVtx()->v4().Vect(), pPionNeg->p4().Vect(),
        pTauNeg->prodVtx()->v4().Vect());
    m_phiCPPion = phiCP_ImpactParameter(imParamPos, imParamNeg, pPionPos->p4(),
                                        pPionNeg->p4());

    const xAOD::TrackParticle *tauPosTrack = tauPosJet->track(0)->track();
    const xAOD::TrackParticle *tauNegTrack = tauNegJet->track(0)->track();

    TVector3 pionPosImParamTruth = calculateTrackImpactParameter(
        tauPosTrack, pTauPos->prodVtx()->v4().Vect(), beamSpot);
    TVector3 pionNegImParamTruth = calculateTrackImpactParameter(
        tauNegTrack, pTauNeg->prodVtx()->v4().Vect(), beamSpot);
    m_phiCPPionJet =
        phiCP_ImpactParameter(pionPosImParamTruth, pionNegImParamTruth,
                              tauPosTrack->p4(), tauNegTrack->p4());

    TVector3 pionPosImParam = calculateTrackImpactParameter(
        tauPosTrack, GetVertexVector(tauPosJet->vertex()), beamSpot);
    TVector3 pionNegImParam = calculateTrackImpactParameter(
        tauNegTrack, GetVertexVector(tauNegJet->vertex()), beamSpot);
    m_phiCPPionJetReco = phiCP_ImpactParameter(
        pionPosImParam, pionNegImParam, tauPosTrack->p4(), tauNegTrack->p4());
  } else if (tauNegDecayMode == TauDecayMode::LEPTONIC &&
             tauPosDecayMode == TauDecayMode::HADRONIC_1P0N) {
    const xAOD::TauJet *tauPosJet = GetLeadingJet(tauJets, true);

    if (tauPosJet == nullptr) {
      ANA_MSG_VERBOSE("Could not find tau+ jet. Excluding event.");
      return StatusCode::SUCCESS;
    }

    const xAOD::Electron *electron = GetLeadingElectron(electrons, false);

    if (electron == nullptr) {
      ANA_MSG_VERBOSE("Could not find tau- lepton. Excluding event.");
      return StatusCode::SUCCESS;
    }

    ANA_MSG_DEBUG("Found higgs -> tau+ tau- -> pion+ lepton- decay");

    TVector3 imParamPos = calculateImpactParameter(
        pPionPos->prodVtx()->v4().Vect(), pPionPos->p4().Vect(),
        pTauPos->prodVtx()->v4().Vect());
    TVector3 imParamNeg = calculateImpactParameter(
        pLeptonNeg->prodVtx()->v4().Vect(), pLeptonNeg->p4().Vect(),
        pTauNeg->prodVtx()->v4().Vect());
    m_phiCPPion = phiCP_ImpactParameter(imParamPos, imParamNeg, pPionPos->p4(),
                                        pLeptonNeg->p4());

    const xAOD::TrackParticle *tauPosTrack = tauPosJet->track(0)->track();
    const xAOD::TrackParticle *tauNegTrack = electron->trackParticle(0);

    TVector3 pionPosImParamTruth = calculateTrackImpactParameter(
        tauPosTrack, pTauPos->prodVtx()->v4().Vect(), beamSpot);
    TVector3 pionNegImParamTruth = calculateTrackImpactParameter(
        tauNegTrack, pTauNeg->prodVtx()->v4().Vect(), beamSpot);
    m_phiCPPionJet =
        phiCP_ImpactParameter(pionPosImParamTruth, pionNegImParamTruth,
                              tauPosTrack->p4(), tauNegTrack->p4());

    TVector3 pionPosImParam =
        calculateTrackImpactParameter(tauPosTrack, primaryVertex, beamSpot);
    TVector3 pionNegImParam =
        calculateTrackImpactParameter(tauNegTrack, primaryVertex, beamSpot);
    m_phiCPPionJetReco = phiCP_ImpactParameter(
        pionPosImParam, pionNegImParam, tauPosTrack->p4(), tauNegTrack->p4());
  } else if (tauNegDecayMode == TauDecayMode::HADRONIC_1P0N &&
             tauPosDecayMode == TauDecayMode::LEPTONIC) {
    ANA_MSG_DEBUG("Found higgs -> tau+ tau- -> lepton+ pion- decay");
    // TVector3 imParamPos = calculateImpactParameter(
    //     pLeptonPos->prodVtx()->v4().Vect(), pLeptonPos->p4().Vect(),
    //     pTauPos->prodVtx()->v4().Vect());
    // TVector3 imParamNeg = calculateImpactParameter(
    //     pPionNeg->prodVtx()->v4().Vect(), pPionNeg->p4().Vect(),
    //     pTauNeg->prodVtx()->v4().Vect());
    // m_phiCPPion = phiCP_ImpactParameter(imParamPos, imParamNeg,
    //                                     pLeptonPos->p4(), pPionNeg->p4());

    // const xAOD::TrackParticle *tauPosTrack = tauPosJet->track(0)->track();
    // const xAOD::TrackParticle *tauNegTrack = tauNegJet->track(0)->track();

    // TVector3 pionPosImParamTruth = calculateTrackImpactParameter(
    //     tauPosTrack, pTauPos->prodVtx()->v4().Vect(), beamSpot);
    // TVector3 pionNegImParamTruth = calculateTrackImpactParameter(
    //     tauNegTrack, pTauNeg->prodVtx()->v4().Vect(), beamSpot);
    // m_phiCPPionJet =
    //     phiCP_ImpactParameter(pionPosImParamTruth, pionNegImParamTruth,
    //                           tauPosTrack->p4(), tauNegTrack->p4());

    // TVector3 pionPosImParam = calculateTrackImpactParameter(
    //     tauPosTrack, GetVertexVector(tauPosJet->vertex()), beamSpot);
    // TVector3 pionNegImParam = calculateTrackImpactParameter(
    //     tauNegTrack, GetVertexVector(tauNegJet->vertex()), beamSpot);
    // m_phiCPPionJetReco = phiCP_ImpactParameter(
    //     pionPosImParam, pionNegImParam, tauPosTrack->p4(),
    //     tauNegTrack->p4());
  } else {
    ANA_MSG_VERBOSE("Unknown tau+ tau- decay mode. Excluding event.");
    return StatusCode::SUCCESS;
  }

  hist("phi_CP_pion")->Fill(m_phiCPPion);
  hist("phi_CP_pion_jet")->Fill(m_phiCPPionJet);
  hist("phi_CP_pion_jet_reco")->Fill(m_phiCPPionJetReco);

  return StatusCode::SUCCESS;
}

StatusCode TruthLevelAnalysis ::finalize() { return StatusCode::SUCCESS; }