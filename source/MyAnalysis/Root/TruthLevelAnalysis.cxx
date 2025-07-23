#include "AsgMessaging/MessageCheck.h"
#include "MyAnalysis/Observables.h"
#include "xAODEgamma/Electron.h"
#include "xAODTracking/TrackParticlexAODHelpers.h"
#include "xAODTruth/TruthParticleContainer.h"
#include "xAODTruth/versions/TruthVertex_v1.h"
#include <MyAnalysis/TruthLevelAnalysis.h>
#include <MyAnalysis/Utils.h>
#include <TTree.h>
#include <TruthUtils/AtlasPID.h>
#include <xAODEventInfo/EventInfo.h>

TruthLevelAnalysis::TruthLevelAnalysis(const std::string &name,
                                       ISvcLocator *pSvcLocator)
    : EL::AnaAlgorithm(name, pSvcLocator) {}

StatusCode TruthLevelAnalysis::initialize() {
  ANA_CHECK(book(TTree("tau_analysis", "tau analysis")));
  TTree *myTree = tree("tau_analysis");

  // Hadronic observables
  myTree->Branch("phiCP_1p0n_1p0n_truth", &m_phiCP_1p0n_1p0n_truth);
  myTree->Branch("phiCP_1p0n_1p0n_recon", &m_phiCP_1p0n_1p0n_recon);
  myTree->Branch("phiCP_1p1n_1p1n_truth", &m_phiCP_1p1n_1p1n_truth);
  myTree->Branch("phiCP_1p1n_1p1n_recon", &m_phiCP_1p1n_1p1n_recon);
  myTree->Branch("phiCP_1p1n_1pXn_truth", &m_phiCP_1p1n_1pXn_truth);
  myTree->Branch("phiCP_1p1n_1pXn_recon", &m_phiCP_1p1n_1pXn_recon);

  // Leptonic observables
  myTree->Branch("phiCP_lept_1p0n_truth", &m_phiCP_lept_1p0n_truth);
  myTree->Branch("phiCP_lept_1p0n_recon", &m_phiCP_lept_1p0n_recon);
  myTree->Branch("phiCP_lept_1p1n_truth", &m_phiCP_lept_1p1n_truth);
  myTree->Branch("phiCP_lept_1p1n_recon", &m_phiCP_lept_1p1n_recon);
  myTree->Branch("phiCP_lept_1pXn_truth", &m_phiCP_lept_1pXn_truth);
  myTree->Branch("phiCP_lept_1pXn_recon", &m_phiCP_lept_1pXn_recon);

  // For debugging purposes
  myTree->Branch("d0_sig_tau_pos_track", &m_d0_sig_tau_pos_track);
  myTree->Branch("d0_sig_tau_neg_track", &m_d0_sig_tau_neg_track);
  myTree->Branch("tau_jets_vtx_diff", &m_tau_jets_vtx_diff);

  return StatusCode::SUCCESS;
}

StatusCode TruthLevelAnalysis::execute() {
  m_phiCP_lept_1p0n_truth = -99.0;
  m_phiCP_lept_1p0n_recon = -99.0;
  m_phiCP_1p0n_1p0n_truth = -99.0;
  m_phiCP_1p0n_1p0n_recon = -99.0;
  m_phiCP_1p1n_1p1n_truth = -99.0;
  m_phiCP_1p1n_1p1n_recon = -99.0;
  m_phiCP_1p1n_1pXn_truth = -99.0;
  m_phiCP_1p1n_1pXn_recon = -99.0;

  m_d0_sig_tau_pos_track = 0.0;
  m_d0_sig_tau_neg_track = 0.0;
  m_tau_jets_vtx_diff = 0.0;

  // Containers
  const xAOD::EventInfo *eventInfo = nullptr;
  const xAOD::TruthParticleContainer *truthHiggsWithDecayParticles = nullptr;
  const xAOD::TruthParticleContainer *truthTausWithDecayParticles = nullptr;
  const xAOD::TauJetContainer *tauJets = nullptr;
  const xAOD::ElectronContainer *electrons = nullptr;
  const xAOD::VertexContainer *vertices = nullptr;

  // Truth particles
  const xAOD::TruthParticle *pHiggs = nullptr;
  const xAOD::TruthParticle *pTauPos = nullptr;
  const xAOD::TruthParticle *pTauNeg = nullptr;
  const xAOD::TruthParticle *pPionPos = nullptr;
  const xAOD::TruthParticle *pPionNeg = nullptr;
  std::vector<const xAOD::TruthParticle *> pPionZerosOfTauPos;
  std::vector<const xAOD::TruthParticle *> pPionZerosOfTauNeg;
  const xAOD::TruthParticle *pLeptonNeg = nullptr;
  const xAOD::TruthParticle *pLeptonPos = nullptr;

  // Retrieve containers
  ANA_CHECK(evtStore()->retrieve(eventInfo, "EventInfo"));
  ANA_CHECK(evtStore()->retrieve(truthTausWithDecayParticles,
                                 "TruthTausWithDecayParticles"));
  ANA_CHECK(evtStore()->retrieve(tauJets, "TauJets"));
  ANA_CHECK(evtStore()->retrieve(electrons, "Electrons"));
  ANA_CHECK(evtStore()->retrieve(vertices, "PrimaryVertices"));

  StatusCode result = evtStore()->retrieve(truthHiggsWithDecayParticles,
                                           "TruthBSMWithDecayParticles");
  if (result.isFailure() || truthHiggsWithDecayParticles->empty()) {
    ANA_CHECK(evtStore()->retrieve(truthHiggsWithDecayParticles,
                                   "TruthBosonsWithDecayParticles"));
    ANA_MSG_DEBUG("No BSM Higgs found, falling back to standard Higgs.");
  } else {
    ANA_MSG_DEBUG("Found BSM Higgs with decay products.");
  }

  // Retrieve beamspot and primary vertex
  TVector3 beamSpot(eventInfo->beamPosX(), eventInfo->beamPosY(),
                    eventInfo->beamPosZ());

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

  // Retrieve Higgs decay products
  int tauCount = 0;

  for (const xAOD::TruthParticle *particle : *truthHiggsWithDecayParticles) {
    // Skip higgs bosons
    if (particle->nParents() == 0) {
      continue;
    }

    // Skip particles which are not directly decay products of the Higgs boson
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

  if (pHiggs == nullptr) {
    ANA_MSG_VERBOSE("Higgs boson not found. Excluding event.");
    return StatusCode::SUCCESS;
  }

  if (pTauNeg == nullptr || pTauPos == nullptr || tauCount != 2) {
    ANA_MSG_VERBOSE(
        "Could not find tau+ tau- decay products. Excluding event.");
    return StatusCode::SUCCESS;
  }

  // Retrieve tau decay products
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
        pPionZerosOfTauNeg.push_back(particle);
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
        pPionZerosOfTauPos.push_back(particle);
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

  // Identify tau decay mode and construct phiCP observables
  TauDecayMode tauNegDecayMode =
      inferTauDecayMode(tauNegLeptonCount, tauNegPionChargedCount,
                        tauNegPionZeroCount, tauNegNeutrinoCount);
  TauDecayMode tauPosDecayMode =
      inferTauDecayMode(tauPosLeptonCount, tauPosPionChargedCount,
                        tauPosPionZeroCount, tauPosNeutrinoCount);

  if (tauNegDecayMode == TauDecayMode::HADRONIC_1P0N &&
      tauPosDecayMode == TauDecayMode::HADRONIC_1P0N) {
    if (tauJets == nullptr || tauJets->size() < 2) {
      ANA_MSG_VERBOSE("Not enough tau jets found. Excluding event.");
      return StatusCode::SUCCESS;
    }

    const xAOD::TauJet *tauPosJet = GetLeadingJet(tauJets, true);
    const xAOD::TauJet *tauNegJet = GetLeadingJet(tauJets, false);

    if (tauPosJet == nullptr || tauNegJet == nullptr) {
      ANA_MSG_VERBOSE("Could not find tau+ or tau- jets. Excluding event.");
      return StatusCode::SUCCESS;
    }

    if (tauPosJet->vertex() == nullptr || tauNegJet->vertex() == nullptr) {
      ANA_MSG_VERBOSE("Could not find tau+ or tau- vertex. Excluding event.");
      return StatusCode::SUCCESS;
    }

    ANA_MSG_DEBUG("Found higgs -> tau+ tau- -> pion+ pion- decay");

    m_tau_jets_vtx_diff = (GetVertexVector(tauPosJet->vertex()) -
                           GetVertexVector(tauNegJet->vertex()))
                              .Mag();

    TVector3 imParamPos = calculateImpactParameter(
        pPionPos->prodVtx()->v4().Vect(), pPionPos->p4().Vect(),
        pTauPos->prodVtx()->v4().Vect());
    TVector3 imParamNeg = calculateImpactParameter(
        pPionNeg->prodVtx()->v4().Vect(), pPionNeg->p4().Vect(),
        pTauNeg->prodVtx()->v4().Vect());
    m_phiCP_1p0n_1p0n_truth =
        phiCP_ImpactParameter(imParamPos, imParamNeg, pPionPos->p4(),
                              pPionNeg->p4(), pPionPos->p4() + pPionNeg->p4());

    const xAOD::TrackParticle *tauPosTrack = tauPosJet->track(0)->track();
    const xAOD::TrackParticle *tauNegTrack = tauNegJet->track(0)->track();

    m_d0_sig_tau_pos_track = xAOD::TrackingHelpers::d0significance(
        tauPosTrack, eventInfo->beamPosSigmaX(), eventInfo->beamPosSigmaY(),
        eventInfo->beamPosSigmaXY());
    m_d0_sig_tau_neg_track = xAOD::TrackingHelpers::d0significance(
        tauNegTrack, eventInfo->beamPosSigmaX(), eventInfo->beamPosSigmaY(),
        eventInfo->beamPosSigmaXY());

    TVector3 pionPosImParamJetVertex = calculateTrackImpactParameter(
        tauPosTrack, GetVertexVector(tauPosJet->vertex()) - beamSpot);
    TVector3 pionNegImParamJetVertex = calculateTrackImpactParameter(
        tauNegTrack, GetVertexVector(tauNegJet->vertex()) - beamSpot);

    m_phiCP_1p0n_1p0n_recon = phiCP_ImpactParameter(
        pionPosImParamJetVertex, pionNegImParamJetVertex, tauPosTrack->p4(),
        tauNegTrack->p4(), tauPosTrack->p4() + tauNegTrack->p4());
  } else if (tauNegDecayMode == TauDecayMode::LEPTONIC &&
             tauPosDecayMode == TauDecayMode::HADRONIC_1P0N) {
    const xAOD::TauJet *tauPosJet = GetLeadingJet(tauJets, true);
    const xAOD::Electron *electron = GetLeadingElectron(electrons, false);

    if (tauPosJet == nullptr) {
      ANA_MSG_VERBOSE("Could not find tau+ jet. Excluding event.");
      return StatusCode::SUCCESS;
    }

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
    m_phiCP_lept_1p0n_truth = phiCP_ImpactParameter(
        imParamPos, imParamNeg, pPionPos->p4(), pLeptonNeg->p4(),
        pPionPos->p4() + pLeptonNeg->p4());

    const xAOD::TrackParticle *tauPosTrack = tauPosJet->track(0)->track();
    const xAOD::TrackParticle *tauNegTrack = electron->trackParticle(0);

    m_d0_sig_tau_pos_track = xAOD::TrackingHelpers::d0significance(
        tauPosTrack, eventInfo->beamPosSigmaX(), eventInfo->beamPosSigmaY(),
        eventInfo->beamPosSigmaXY());
    m_d0_sig_tau_neg_track = xAOD::TrackingHelpers::d0significance(
        tauNegTrack, eventInfo->beamPosSigmaX(), eventInfo->beamPosSigmaY(),
        eventInfo->beamPosSigmaXY());

    TVector3 pionPosImParam = calculateTrackImpactParameter(
        tauPosTrack, GetVertexVector(tauPosJet->vertex()) - beamSpot);
    TVector3 pionNegImParam = calculateTrackImpactParameter(
        tauNegTrack, GetVertexVector(tauPosJet->vertex()) - beamSpot);
    m_phiCP_lept_1p0n_recon = phiCP_ImpactParameter(
        pionPosImParam, pionNegImParam, tauPosTrack->p4(), tauNegTrack->p4(),
        tauPosJet->p4() + electron->p4());

    // leptonic correction:
    m_phiCP_lept_1p0n_truth = m_phiCP_lept_1p0n_truth < M_PI
                                  ? m_phiCP_lept_1p0n_truth + M_PI
                                  : m_phiCP_lept_1p0n_truth - M_PI;
    m_phiCP_lept_1p0n_recon = m_phiCP_lept_1p0n_recon < M_PI
                                  ? m_phiCP_lept_1p0n_recon + M_PI
                                  : m_phiCP_lept_1p0n_recon - M_PI;
  } else if (tauNegDecayMode == TauDecayMode::HADRONIC_1P0N &&
             tauPosDecayMode == TauDecayMode::LEPTONIC) {
    const xAOD::TauJet *tauNegJet = GetLeadingJet(tauJets, false);
    const xAOD::Electron *positron = GetLeadingElectron(electrons, true);

    if (tauNegJet == nullptr) {
      ANA_MSG_VERBOSE("Could not find tau+ jet. Excluding event.");
      return StatusCode::SUCCESS;
    }

    if (positron == nullptr) {
      ANA_MSG_VERBOSE("Could not find tau- lepton. Excluding event.");
      return StatusCode::SUCCESS;
    }

    ANA_MSG_DEBUG("Found higgs -> tau+ tau- -> lepton+ pion- decay");

    TVector3 imParamPos = calculateImpactParameter(
        pLeptonPos->prodVtx()->v4().Vect(), pLeptonPos->p4().Vect(),
        pTauPos->prodVtx()->v4().Vect());
    TVector3 imParamNeg = calculateImpactParameter(
        pPionNeg->prodVtx()->v4().Vect(), pPionNeg->p4().Vect(),
        pTauNeg->prodVtx()->v4().Vect());
    m_phiCP_lept_1p0n_truth = phiCP_ImpactParameter(
        imParamPos, imParamNeg, pLeptonPos->p4(), pPionNeg->p4(),
        pLeptonPos->p4() + pPionNeg->p4());

    const xAOD::TrackParticle *tauNegTrack = tauNegJet->track(0)->track();
    const xAOD::TrackParticle *tauPosTrack = positron->trackParticle(0);

    m_d0_sig_tau_pos_track = xAOD::TrackingHelpers::d0significance(
        tauPosTrack, eventInfo->beamPosSigmaX(), eventInfo->beamPosSigmaY(),
        eventInfo->beamPosSigmaXY());
    m_d0_sig_tau_neg_track = xAOD::TrackingHelpers::d0significance(
        tauNegTrack, eventInfo->beamPosSigmaX(), eventInfo->beamPosSigmaY(),
        eventInfo->beamPosSigmaXY());

    TVector3 pionPosImParam = calculateTrackImpactParameter(
        tauPosTrack, GetVertexVector(tauNegJet->vertex()) - beamSpot);
    TVector3 pionNegImParam = calculateTrackImpactParameter(
        tauNegTrack, GetVertexVector(tauNegJet->vertex()) - beamSpot);
    m_phiCP_lept_1p0n_recon = phiCP_ImpactParameter(
        pionPosImParam, pionNegImParam, tauPosTrack->p4(), tauNegTrack->p4(),
        tauNegJet->p4() + positron->p4());

    // leptonic correction:
    m_phiCP_lept_1p0n_truth = m_phiCP_lept_1p0n_truth < M_PI
                                  ? m_phiCP_lept_1p0n_truth + M_PI
                                  : m_phiCP_lept_1p0n_truth - M_PI;
    m_phiCP_lept_1p0n_recon = m_phiCP_lept_1p0n_recon < M_PI
                                  ? m_phiCP_lept_1p0n_recon + M_PI
                                  : m_phiCP_lept_1p0n_recon - M_PI;
  } else if ((tauNegDecayMode == TauDecayMode::HADRONIC_1P1N &&
              tauPosDecayMode == TauDecayMode::HADRONIC_1P1N) ||
             (tauNegDecayMode == TauDecayMode::HADRONIC_1P1N &&
              tauPosDecayMode == TauDecayMode::HADRONIC_1PXN) ||
             (tauNegDecayMode == TauDecayMode::HADRONIC_1PXN &&
              tauPosDecayMode == TauDecayMode::HADRONIC_1P1N)) {
    if (tauJets == nullptr || tauJets->size() < 2) {
      ANA_MSG_VERBOSE("Not enough tau jets found. Excluding event.");
      return StatusCode::SUCCESS;
    }

    const xAOD::TauJet *tauPosJet = GetLeadingJet(tauJets, true);
    const xAOD::TauJet *tauNegJet = GetLeadingJet(tauJets, false);

    if (tauPosJet == nullptr || tauNegJet == nullptr) {
      ANA_MSG_VERBOSE("Could not find tau+ or tau- jets. Excluding event.");
      return StatusCode::SUCCESS;
    }

    if (tauPosJet->vertex() == nullptr || tauNegJet->vertex() == nullptr) {
      ANA_MSG_VERBOSE("Could not find tau+ or tau- vertex. Excluding event.");
      return StatusCode::SUCCESS;
    }

    ANA_MSG_DEBUG("Found higgs -> tau+ tau- -> pion+ pion- pion0 decay");

    m_tau_jets_vtx_diff = (GetVertexVector(tauPosJet->vertex()) -
                           GetVertexVector(tauNegJet->vertex()))
                              .Mag();

    // sum over neutral pions
    TLorentzVector chargedP4Pos = pPionPos->p4();
    TLorentzVector neutralP4Pos(0.0, 0.0, 0.0, 0.0);
    for (const xAOD::TruthParticle *pionZero : pPionZerosOfTauPos) {
      neutralP4Pos += pionZero->p4();
    }

    TLorentzVector chargedP4Neg = pPionNeg->p4();
    TLorentzVector neutralP4Neg(0.0, 0.0, 0.0, 0.0);
    for (const xAOD::TruthParticle *pionZero : pPionZerosOfTauNeg) {
      neutralP4Neg += pionZero->p4();
    }

    double phiCP_truth =
        phiCP_Pion_RhoDecayPlane(chargedP4Pos, neutralP4Pos, chargedP4Neg,
                                 neutralP4Neg, pTauPos->p4() + pTauNeg->p4());

    chargedP4Pos.SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
    for (auto track : tauPosJet->tracks()) {
      chargedP4Pos += track->track()->p4();
    }

    neutralP4Pos.SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
    for (size_t i = 0; i < tauPosJet->nNeutralPFOs(); ++i) {
      neutralP4Pos += tauPosJet->neutralPFO(i)->p4();
    }

    chargedP4Neg.SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
    for (auto track : tauNegJet->tracks()) {
      chargedP4Neg += track->track()->p4();
    }

    neutralP4Neg.SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
    for (size_t i = 0; i < tauNegJet->nNeutralPFOs(); ++i) {
      neutralP4Neg += tauNegJet->neutralPFO(i)->p4();
    }

    double phiCP_recon =
        phiCP_Pion_RhoDecayPlane(chargedP4Pos, neutralP4Pos, chargedP4Neg,
                                 neutralP4Neg, chargedP4Pos + chargedP4Neg);

    if (tauNegDecayMode == TauDecayMode::HADRONIC_1PXN ||
        tauPosDecayMode == TauDecayMode::HADRONIC_1PXN) {
      m_phiCP_1p1n_1pXn_truth = phiCP_truth;
      m_phiCP_1p1n_1pXn_recon = phiCP_recon;
    } else {
      m_phiCP_1p1n_1p1n_truth = phiCP_truth;
      m_phiCP_1p1n_1p1n_recon = phiCP_recon;
    }
  } else if (tauNegDecayMode == TauDecayMode::LEPTONIC &&
             (tauPosDecayMode == TauDecayMode::HADRONIC_1P1N ||
              tauPosDecayMode == TauDecayMode::HADRONIC_1PXN)) {
    const xAOD::TauJet *tauPosJet = GetLeadingJet(tauJets, true);
    const xAOD::Electron *electron = GetLeadingElectron(electrons, false);

    if (tauPosJet == nullptr) {
      ANA_MSG_VERBOSE("Could not find tau+ jet. Excluding event.");
      return StatusCode::SUCCESS;
    }

    if (electron == nullptr) {
      ANA_MSG_VERBOSE("Could not find tau- lepton. Excluding event.");
      return StatusCode::SUCCESS;
    }

    ANA_MSG_DEBUG("Found higgs -> tau+ tau- -> pion+ pion0 lepton- decay");

    // sum over neutral pions
    TLorentzVector chargedP4Pos = pPionPos->p4();
    TLorentzVector neutralP4Pos(0.0, 0.0, 0.0, 0.0);
    for (const xAOD::TruthParticle *pionZero : pPionZerosOfTauPos) {
      ANA_MSG_DEBUG("Found neutral pion in tau+ decay: "
                    << pionZero->pdgId() << " " << pionZero->barcode());
      neutralP4Pos += pionZero->p4();
    }

    TVector3 imParamNeg = calculateImpactParameter(
        pLeptonNeg->prodVtx()->v4().Vect(), pLeptonNeg->p4().Vect(),
        pTauNeg->prodVtx()->v4().Vect());
    double phiCP_truth =
        phiCP_IP_Rho(imParamNeg, pLeptonNeg->p4(), chargedP4Pos, neutralP4Pos,
                     pTauPos->p4() + pTauNeg->p4(), true);

    const xAOD::TrackParticle *tauNegTrack = electron->trackParticle(0);

    chargedP4Pos.SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
    for (auto track : tauPosJet->tracks()) {
      chargedP4Pos += track->track()->p4();
    }

    neutralP4Pos.SetPxPyPzE(0.0, 0.0, 0.0, 0.0);
    for (size_t i = 0; i < tauPosJet->nNeutralPFOs(); ++i) {
      neutralP4Pos += tauPosJet->neutralPFO(i)->p4();
    }

    TVector3 pionNegImParam = calculateTrackImpactParameter(
        tauNegTrack, GetVertexVector(tauPosJet->vertex()) - beamSpot);
    double phiCP_recon =
        phiCP_IP_Rho(pionNegImParam, tauNegTrack->p4(), chargedP4Pos,
                     neutralP4Pos, tauPosJet->p4() + electron->p4(), true);

    // leptonic correction:
    phiCP_truth = phiCP_truth < M_PI ? phiCP_truth + M_PI : phiCP_truth - M_PI;
    phiCP_recon = phiCP_recon < M_PI ? phiCP_recon + M_PI : phiCP_recon - M_PI;

    if (tauPosDecayMode == TauDecayMode::HADRONIC_1PXN) {
      m_phiCP_lept_1pXn_truth = phiCP_truth;
      m_phiCP_lept_1pXn_recon = phiCP_recon;
    } else {
      m_phiCP_lept_1p1n_truth = phiCP_truth;
      m_phiCP_lept_1p1n_recon = phiCP_recon;
    }
  } else {
    ANA_MSG_VERBOSE("Unknown tau+ tau- decay mode. Excluding event.");
    return StatusCode::SUCCESS;
  }

  tree("tau_analysis")->Fill();

  return StatusCode::SUCCESS;
}

StatusCode TruthLevelAnalysis ::finalize() { return StatusCode::SUCCESS; }
