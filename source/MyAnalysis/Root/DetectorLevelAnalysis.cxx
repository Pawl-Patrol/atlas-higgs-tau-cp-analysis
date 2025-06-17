#include "AsgMessaging/MessageCheck.h"
#include "MyAnalysis/Observables.h"
#include "xAODBase/IParticle.h"
#include "xAODJet/JetContainer.h"
#include "xAODTracking/TrackParticleFwd.h"
#include "xAODTracking/TrackingPrimitives.h"
#include "xAODTracking/VertexContainer.h"
#include "xAODTruth/TruthParticleContainer.h"
#include "xAODTruth/TruthVertex.h"
#include "xAODTruth/TruthVertexContainer.h"
#include "xAODTruth/versions/TruthVertex_v1.h"
#include <MyAnalysis/DetectorLevelAnalysis.h>
#include <TLorentzVector.h>
#include <xAODTracking/TrackParticle.h>

DetectorLevelAnalysis::DetectorLevelAnalysis(const std::string &name,
                                             ISvcLocator *pSvcLocator)
    : EL::AnaAlgorithm(name, pSvcLocator) {}

StatusCode DetectorLevelAnalysis::initialize() {
  const int BINS = 50;

  // Setup histograms
  ANA_CHECK(book(TH1F("phi_CP_pion", "phi_CP_pion", BINS, 0, 2 * M_PI)));

  return StatusCode::SUCCESS;
}

StatusCode DetectorLevelAnalysis::execute() {
  m_phiCP = 0.;
  const xAOD::TruthParticleContainer *truthHiggsWithDecayParticles = nullptr;
  StatusCode result = evtStore()->retrieve(truthHiggsWithDecayParticles,
                                           "TruthBSMWithDecayParticles");

  if (result.isFailure() || truthHiggsWithDecayParticles->empty()) {
    ANA_CHECK(evtStore()->retrieve(truthHiggsWithDecayParticles,
                                   "TruthBosonsWithDecayParticles"));
  }

  const xAOD::TruthParticleContainer *truthTausWithDecayParticles = nullptr;
  ANA_CHECK(evtStore()->retrieve(truthTausWithDecayParticles,
                                 "TruthTausWithDecayParticles"));

  const xAOD::TruthParticle *pHiggs = nullptr;
  const xAOD::TruthParticle *pTauPos = nullptr;
  const xAOD::TruthParticle *pTauNeg = nullptr;
  const xAOD::TruthParticle *pPionPos = nullptr;
  const xAOD::TruthParticle *pPionNeg = nullptr;
  const xAOD::TruthParticle *pNeutrino = nullptr;
  const xAOD::TruthParticle *pAntiNeutrino = nullptr;

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

  if (pHiggs == nullptr) {
    ANA_MSG_DEBUG("Higgs boson not found. Excluding event.");
    return StatusCode::SUCCESS;
  }

  if (pTauNeg == nullptr || pTauPos == nullptr || tauCount != 2) {
    ANA_MSG_DEBUG("Could not find tau+ tau- decay products. Excluding event.");
    return StatusCode::SUCCESS;
  }

  if (pPionPos == nullptr || pPionNeg == nullptr || pionCount > 2) {
    ANA_MSG_DEBUG("Could not find pi+ pi- decay products. Excluding event.");
    return StatusCode::SUCCESS;
  }

  if (pNeutrino == nullptr || pAntiNeutrino == nullptr || neutrinoCount > 2) {
    ANA_MSG_DEBUG("Could not find neutrino decay products. Excluding event.");
    return StatusCode::SUCCESS;
  }

  ANA_MSG_DEBUG("Found higgs -> tau+ tau- -> pi+ pi- decay.");

  ANA_MSG_DEBUG("TauNeg prodVtx: " << pTauNeg->prodVtx()->v4().X() << " "
                                   << pTauNeg->prodVtx()->v4().Y() << " "
                                   << pTauNeg->prodVtx()->v4().Z());

  const xAOD::TruthVertexContainer *vertices = nullptr;
  ANA_CHECK(evtStore()->retrieve(vertices, "TruthPrimaryVertices"));

  TLorentzVector primaryVertex = TLorentzVector(0., 0., 0., 0.);
  for (const xAOD::TruthVertex *vertex : *vertices) {
    ANA_MSG_DEBUG("Vertex: " << vertex->v4().X() << " " << vertex->v4().Y()
                             << " " << vertex->v4().Z());
  }

  m_phiCP = phiCP_Pion_ImpactParameter(
      pPionPos->prodVtx()->v4(), pPionNeg->prodVtx()->v4(),
      pTauNeg->prodVtx()->v4(), pPionPos->p4(), pPionNeg->p4());

  hist("phi_CP_pion")->Fill(m_phiCP);

  return StatusCode::SUCCESS;
}

StatusCode DetectorLevelAnalysis ::finalize() { return StatusCode::SUCCESS; }