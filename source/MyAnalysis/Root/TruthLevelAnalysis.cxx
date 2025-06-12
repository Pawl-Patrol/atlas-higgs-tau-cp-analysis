#include "AsgMessaging/MessageCheck.h"
#include "MyAnalysis/Observables.h"
#include <MyAnalysis/TruthLevelAnalysis.h>
#include <TLorentzVector.h>

TruthLevelAnalysis::TruthLevelAnalysis(const std::string &name,
                                       ISvcLocator *pSvcLocator)
    : EL::AnaAlgorithm(name, pSvcLocator) {}

StatusCode TruthLevelAnalysis::initialize() {
  ANA_MSG_DEBUG("Initializing");

  const int BINS = 50;

  // Setup histograms
  ANA_CHECK(book(TH1F("phi_CP_tau_pi", "phi_CP_tau_pi", BINS, 0, 2 * M_PI)));
  ANA_CHECK(book(
      TH1F("phi_CP_neutrino_pi", "phi_CP_neutrino_pi", BINS, 0, 2 * M_PI)));
  ANA_CHECK(book(TH1F("phi_CP_pion", "phi_CP_pion", BINS, 0, 2 * M_PI)));

  return StatusCode::SUCCESS;
}

StatusCode TruthLevelAnalysis::execute() {
  // Here you do everything that needs to be done on every single
  // events, e.g. read input variables, apply cuts, and fill
  // histograms and trees. This is where most of your actual analysis
  // code will go.
  const xAOD::EventInfo *eventInfo = nullptr;
  ANA_CHECK(evtStore()->retrieve(eventInfo, "EventInfo"));

  m_runNumber = eventInfo->runNumber();
  m_eventNumber = eventInfo->eventNumber();
  m_phiCP = 0.;
  m_phiCPNeutri = 0.;
  m_phiCPPion = 0.;

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
    // ignore the original higgs boson
    TLorentzVector m = particle->p4();
    ANA_MSG_DEBUG("Processing particle with PDG ID: "
                  << particle->pdgId() << " And id: " << particle->id());

    if (particle->nParents() == 0) {
      continue;
    }

    ANA_MSG_DEBUG("Parent ID of particle: " << particle->parent(0)->id());

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
    ANA_MSG_DEBUG("Processing particle with PDG ID: "
                  << particle->pdgId() << " And id: " << particle->id());
    if (particle->nParents() == 0) {
      continue;
    }
    ANA_MSG_DEBUG("Parent ID of particle: " << particle->parent(0)->id());

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

  // Calculate observables
  m_phiCP = phiCP_Pion_Tau(pHiggs->p4(), pTauPos->p4(), pTauNeg->p4(),
                           pPionPos->p4(), pPionNeg->p4());

  m_phiCPNeutri =
      phiCP_Pion_Neutrino(pHiggs->p4(), pAntiNeutrino->p4(), pNeutrino->p4(),
                          pPionPos->p4(), pPionNeg->p4());

  m_phiCPPion = phiCP_Pion_ImpactParameter(
      pPionPos->prodVtx()->v4(), pPionNeg->prodVtx()->v4(),
      pTauNeg->prodVtx()->v4(), pPionPos->p4(), pPionNeg->p4());

  hist("phi_CP_tau_pi")->Fill(m_phiCP);
  hist("phi_CP_neutrino_pi")->Fill(m_phiCPNeutri);
  hist("phi_CP_pion")->Fill(m_phiCPPion);

  return StatusCode::SUCCESS;
}

StatusCode TruthLevelAnalysis ::finalize() {
  // This method is the mirror image of initialize(), meaning it gets
  // called after the last event has been processed on the worker node
  // and allows you to finish up any objects you created in
  // initialize() before they are written to disk. This is actually
  // fairly rare, since this happens separately for each worker node.
  // Most of the time you want to do your post-processing on the
  // submission node after all your histogram outputs have been
  // merged.
  return StatusCode::SUCCESS;
}