#include "AsgMessaging/MessageCheck.h"
#include "MyAnalysis/Observables.h"
#include "xAODTau/TauJet.h"
#include "xAODTau/TauTrack.h"
#include "xAODTracking/TrackParticle.h"
#include "xAODTracking/TrackingPrimitives.h"
#include "xAODTracking/VertexContainer.h"
#include "xAODTruth/TruthParticleContainer.h"
#include "xAODTruth/TruthVertexContainerFwd.h"
#include "xAODTruth/TruthVertexFwd.h"
#include "xAODTruth/versions/TruthVertex_v1.h"
#include <MyAnalysis/TruthLevelAnalysis.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <boost/date_time/constrained_value.hpp>
#include <xAODTau/TauJetContainer.h>
#include <xAODTruth/xAODTruthHelpers.h>

TruthLevelAnalysis::TruthLevelAnalysis(const std::string &name,
                                       ISvcLocator *pSvcLocator)
    : EL::AnaAlgorithm(name, pSvcLocator) {}

StatusCode TruthLevelAnalysis::initialize() {
  const int BINS = 50;

  // Setup histograms
  ANA_CHECK(book(TH1F("phi_CP_tau_pi", "phi_CP_tau_pi", BINS, 0, 2 * M_PI)));
  ANA_CHECK(book(
      TH1F("phi_CP_neutrino_pi", "phi_CP_neutrino_pi", BINS, 0, 2 * M_PI)));
  ANA_CHECK(book(TH1F("phi_CP_pion", "phi_CP_pion", BINS, 0, 2 * M_PI)));
  ANA_CHECK(
      book(TH1F("phi_CP_pion_jet", "phi_CP_pion_jet", BINS, 0, 2 * M_PI)));
  ANA_CHECK(book(
      TH1F("phi_CP_pion_jet_reco", "phi_CP_pion_jet_reco", BINS, 0, 2 * M_PI)));

  return StatusCode::SUCCESS;
}

StatusCode TruthLevelAnalysis::execute() {
  m_phiCP = 0.;
  m_phiCPNeutri = 0.;
  m_phiCPPion = 0.;
  m_phiCPPionJet = 0.;
  m_phiCPPionJetReco = 0.;

  const xAOD::TruthParticleContainer *truthHiggsWithDecayParticles = nullptr;
  const xAOD::TruthParticleContainer *truthTausWithDecayParticles = nullptr;
  const xAOD::TauJetContainer *tauJets = nullptr;
  const xAOD::TruthVertexContainer *truthVertices = nullptr;
  const xAOD::VertexContainer *vertices = nullptr;

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

  const xAOD::TruthParticle *pHiggs = nullptr;
  const xAOD::TruthParticle *pTauPos = nullptr;
  const xAOD::TruthParticle *pTauNeg = nullptr;
  const xAOD::TruthParticle *pPionPos = nullptr;
  const xAOD::TruthParticle *pPionNeg = nullptr;
  const xAOD::TruthParticle *pNeutrino = nullptr;
  const xAOD::TruthParticle *pAntiNeutrino = nullptr;
  const xAOD::TauJet *tauPosJet = nullptr;
  const xAOD::TauJet *tauNegJet = nullptr;
  const xAOD::TrackParticle *tauPosTrack = nullptr;
  const xAOD::TrackParticle *tauNegTrack = nullptr;

  TLorentzVector truthPrimaryVertex = TLorentzVector(0., 0., 0., 0.);
  TLorentzVector primaryVertex = TLorentzVector(0., 0., 0., 0.);

  int tauCount = 0;
  int pionCount = 0;
  int neutrinoCount = 0;

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

  for (const xAOD::TruthParticle *particle : *truthTausWithDecayParticles) {
    if (particle->nParents() == 0) {
      continue;
    }

    // Pions
    if (particle->pdgId() == PI0) {
      pionCount++;
    } else if (particle->pdgId() == PIPLUS &&
               particle->parent(0)->pdgId() == -TAU) {
      pPionPos = particle;
      pionCount++;
    } else if (particle->pdgId() == PIMINUS &&
               particle->parent(0)->pdgId() == TAU) {
      pPionNeg = particle;
      pionCount++;
    }

    // Neutrinos
    else if (particle->pdgId() == -NU_TAU &&
             particle->parent(0)->pdgId() == -TAU) {
      pAntiNeutrino = particle;
      neutrinoCount++;
    } else if (particle->pdgId() == NU_TAU &&
               particle->parent(0)->pdgId() == TAU) {
      pNeutrino = particle;
      neutrinoCount++;
    }
  }

  for (const xAOD::TauJet *tauJet : *tauJets) {
    if (tauJet->charge() > 0) {
      tauPosJet = tauJet;
    } else if (tauJet->charge() < 0) {
      tauNegJet = tauJet;
    }
  }

  for (const xAOD::TruthVertex *vertex : *truthVertices) {
    truthPrimaryVertex = vertex->v4();
  }

  for (const xAOD::Vertex *vertex : *vertices) {
    if (vertex->vertexType() == xAOD::VxType::PriVtx) {
      primaryVertex = TLorentzVector(vertex->x(), vertex->y(), vertex->z(), 0);
    }
  }

  // Excludsion rules
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

  if (pHiggs == nullptr) {
    ANA_MSG_VERBOSE("Higgs boson not found. Excluding event.");
    return StatusCode::SUCCESS;
  }

  if (pTauNeg == nullptr || pTauPos == nullptr || tauCount != 2) {
    ANA_MSG_VERBOSE(
        "Could not find tau+ tau- decay products. Excluding event.");
    return StatusCode::SUCCESS;
  }

  if (pPionPos == nullptr || pPionNeg == nullptr || pionCount > 2) {
    ANA_MSG_VERBOSE("Could not find pi+ pi- decay products. Excluding event.");
    return StatusCode::SUCCESS;
  }

  if (pNeutrino == nullptr || pAntiNeutrino == nullptr || neutrinoCount > 2) {
    ANA_MSG_VERBOSE("Could not find neutrino decay products. Excluding event.");
    return StatusCode::SUCCESS;
  }

  ANA_MSG_DEBUG("Found higgs -> tau+ tau- -> pi+ pi- decay.");

  // Get the tracks
  tauPosTrack = tauPosJet->track(0)->track();
  tauNegTrack = tauNegJet->track(0)->track();

  if (tauPosJet->vertex() != nullptr) {
    primaryVertex.SetXYZT(tauPosJet->vertex()->x(), tauPosJet->vertex()->y(),
                          tauPosJet->vertex()->z(), 0);
  }

  // log out primary vertices
  ANA_MSG_DEBUG("Truth primary vertex: " << truthPrimaryVertex.X() << " "
                                         << truthPrimaryVertex.Y() << " "
                                         << truthPrimaryVertex.Z() << " "
                                         << truthPrimaryVertex.T());
  ANA_MSG_DEBUG("Primary vertex: "
                << primaryVertex.X() << " " << primaryVertex.Y() << " "
                << primaryVertex.Z() << " " << primaryVertex.T());

  // log out pion and tau momenta
  ANA_MSG_DEBUG("Tau+ P4: " << pTauPos->p4().X() << " " << pTauPos->p4().Y()
                            << " " << pTauPos->p4().Z() << " "
                            << pTauPos->p4().T());
  ANA_MSG_DEBUG("Tau- P4: " << pTauNeg->p4().X() << " " << pTauNeg->p4().Y()
                            << " " << pTauNeg->p4().Z() << " "
                            << pTauNeg->p4().T());
  ANA_MSG_DEBUG("Pion+ P4: " << pPionPos->p4().X() << " " << pPionPos->p4().Y()
                             << " " << pPionPos->p4().Z() << " "
                             << pPionPos->p4().T());
  ANA_MSG_DEBUG("Pion- P4: " << pPionNeg->p4().X() << " " << pPionNeg->p4().Y()
                             << " " << pPionNeg->p4().Z() << " "
                             << pPionNeg->p4().T());

  // log out pion and tau production vertices
  ANA_MSG_DEBUG("Tau+ vertex: " << pTauPos->prodVtx()->v4().X() << " "
                                << pTauPos->prodVtx()->v4().Y() << " "
                                << pTauPos->prodVtx()->v4().Z() << " "
                                << pTauPos->prodVtx()->v4().T());
  ANA_MSG_DEBUG("Tau- vertex: " << pTauNeg->prodVtx()->v4().X() << " "
                                << pTauNeg->prodVtx()->v4().Y() << " "
                                << pTauNeg->prodVtx()->v4().Z() << " "
                                << pTauNeg->prodVtx()->v4().T());

  ANA_MSG_DEBUG("Tau+ Jet: "
                << tauPosJet->vertex()->x() << " " << tauPosJet->vertex()->y()
                << " " << tauPosJet->vertex()->z() << " "
                << tauPosJet->p4().Px() << " " << tauPosJet->p4().Py() << " "
                << tauPosJet->p4().Pz() << " " << tauPosJet->p4().E());

  ANA_MSG_DEBUG("Tau- Jet: "
                << tauNegJet->vertex()->x() << " " << tauNegJet->vertex()->y()
                << " " << tauNegJet->vertex()->z() << " "
                << tauNegJet->p4().Px() << " " << tauNegJet->p4().Py() << " "
                << tauNegJet->p4().Pz() << " " << tauNegJet->p4().E());

  ANA_MSG_DEBUG("Tau+ track: "
                << tauPosTrack->d0() << " " << tauPosTrack->z0() << " "
                << tauPosTrack->phi() << " p4 " << tauPosTrack->p4().Px() << " "
                << tauPosTrack->p4().Py() << " " << tauPosTrack->p4().Pz()
                << " " << tauPosTrack->p4().E());
  ANA_MSG_DEBUG("Tau- track: "
                << tauNegTrack->d0() << " " << tauNegTrack->z0() << " "
                << tauNegTrack->phi() << " p4 " << tauNegTrack->p4().Px() << " "
                << tauNegTrack->p4().Py() << " " << tauNegTrack->p4().Pz()
                << " " << tauNegTrack->p4().E());

  ANA_MSG_DEBUG("Pion+ vertex: " << pPionPos->prodVtx()->v4().X() << " "
                                 << pPionPos->prodVtx()->v4().Y() << " "
                                 << pPionPos->prodVtx()->v4().Z() << " "
                                 << pPionPos->prodVtx()->v4().T());
  ANA_MSG_DEBUG("Pion- vertex: " << pPionNeg->prodVtx()->v4().X() << " "
                                 << pPionNeg->prodVtx()->v4().Y() << " "
                                 << pPionNeg->prodVtx()->v4().Z() << " "
                                 << pPionNeg->prodVtx()->v4().T());

  // Calculate observables
  m_phiCP = phiCP_Pion_Tau(pHiggs->p4(), pTauPos->p4(), pTauNeg->p4(),
                           pPionPos->p4(), pPionNeg->p4());

  m_phiCPNeutri =
      phiCP_Pion_Neutrino(pHiggs->p4(), pAntiNeutrino->p4(), pNeutrino->p4(),
                          pPionPos->p4(), pPionNeg->p4());

  // phi_CP using only truth
  m_phiCPPion = phiCP_Pion_ImpactParameter(
      pPionPos->prodVtx()->v4(), pPionNeg->prodVtx()->v4(), truthPrimaryVertex,
      pPionPos->p4(), pPionNeg->p4());

  // phi_CP using truth but pion track momenta
  m_phiCPPionJet = phiCP_Pion_ImpactParameter(
      pPionPos->prodVtx()->v4(), pPionNeg->prodVtx()->v4(), primaryVertex,
      tauPosTrack->p4(), tauNegTrack->p4());

  TLorentzVector pionVtxPos(-tauPosTrack->d0() * sin(tauPosTrack->phi0()),
                            tauPosTrack->d0() * cos(tauPosTrack->phi0()),
                            tauPosTrack->z0(), 0);
  TLorentzVector pionVtxNeg(-tauNegTrack->d0() * sin(tauNegTrack->phi0()),
                            tauNegTrack->d0() * cos(tauNegTrack->phi0()),
                            tauNegTrack->z0(), 0);

  TLorentzVector beamSpot(-0.5, -0.5, 0.0, 0.0);
  pionVtxPos += beamSpot;
  pionVtxNeg += beamSpot;

  // phi_CP using pion track impact parameter and momenta
  m_phiCPPionJetReco =
      phiCP_Pion_ImpactParameter(pionVtxPos, pionVtxNeg, primaryVertex,
                                 tauPosTrack->p4(), tauNegTrack->p4());

  hist("phi_CP_tau_pi")->Fill(m_phiCP);
  hist("phi_CP_neutrino_pi")->Fill(m_phiCPNeutri);
  hist("phi_CP_pion")->Fill(m_phiCPPion);
  hist("phi_CP_pion_jet")->Fill(m_phiCPPionJet);
  hist("phi_CP_pion_jet_reco")->Fill(m_phiCPPionJetReco);

  return StatusCode::SUCCESS;
}

StatusCode TruthLevelAnalysis ::finalize() { return StatusCode::SUCCESS; }