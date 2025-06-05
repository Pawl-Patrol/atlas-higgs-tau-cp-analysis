#include <AsgMessaging/MessageCheck.h>
#include <Math/GenVector/CoordinateSystemTags.h>
#include <Math/GenVector/LorentzVector.h>
#include <MyAnalysis/MyxAODAnalysis.h>
#include <TTree.h>
#include <TruthUtils/AtlasPID.h>
#include <cmath>
#include <xAODEventInfo/EventInfo.h>
#include <xAODTruth/TruthEventContainer.h>
#include <xAODTruth/TruthParticle.h>
#include <xAODTruth/TruthParticleContainer.h>

MyxAODAnalysis ::MyxAODAnalysis(const std::string &name,
                                ISvcLocator *pSvcLocator)
    : EL::AnaAlgorithm(name, pSvcLocator) {
  // Here you put any code for the base initialization of variables,
  // e.g. initialize all pointers to 0. This is also where you
  // declare all properties for your algorithm. Note that things like
  // resetting statistics variables or booking histograms should
  // rather go into the initialize() function.
}

StatusCode MyxAODAnalysis ::initialize() {
  // Here you do everything that needs to be done at the very
  // beginning on each worker node, e.g. create histograms and output
  // trees. This method gets called before any input files are
  // connected.
  ANA_MSG_INFO("Initializing");

  // Setup tree
  ANA_CHECK(book(TTree("truth_tau_analysis", "Zee analysis ntuple")));
  TTree *myTree = tree("truth_tau_analysis");
  myTree->Branch("run_number", &m_runNumber);
  myTree->Branch("event_number", &m_eventNumber);
  myTree->Branch("phi_CP_tau_pi", &m_phi_CP);
  myTree->Branch("phi_CP_neutrino_pi", &m_phi_CP_neutri);

  return StatusCode::SUCCESS;
}

StatusCode MyxAODAnalysis ::execute() {
  // Here you do everything that needs to be done on every single
  // events, e.g. read input variables, apply cuts, and fill
  // histograms and trees. This is where most of your actual analysis
  // code will go.
  const xAOD::EventInfo *eventInfo = nullptr;
  ANA_CHECK(evtStore()->retrieve(eventInfo, "EventInfo"));

  m_runNumber = eventInfo->runNumber();
  m_eventNumber = eventInfo->eventNumber();

  bool switch_BSM = false;
  const xAOD::TruthParticleContainer *truthHiggs = nullptr;
  ANA_CHECK(evtStore()->retrieve(truthHiggs, "TruthBoson"));
  if (evtStore()->contains<xAOD::TruthParticleContainer>("TruthBSM")) {
    if (evtStore()->retrieve(truthHiggs, "TruthBSM")) {
      if (!truthHiggs->empty()) {
        switch_BSM = true;
      }
    }
  }

  if (switch_BSM == true) {
    ANA_CHECK(evtStore()->retrieve(truthHiggs, "TruthBSM"));
  } else {
    ANA_CHECK(evtStore()->retrieve(truthHiggs, "TruthBoson"));
  }

  const xAOD::TruthParticleContainer *truthTauParticles = nullptr;
  ANA_CHECK(evtStore()->retrieve(truthTauParticles, "TruthTaus"));

  const xAOD::TruthParticleContainer *truthTausWithDecayParticles = nullptr;
  ANA_CHECK(evtStore()->retrieve(truthTausWithDecayParticles,
                                 "TruthTausWithDecayParticles"));

  // Get Higgs particle four-momentum
  TLorentzVector higgs_p4;
  for (const xAOD::TruthParticle *higgs : *truthHiggs) {
    higgs_p4 = higgs->p4();
  }

  // Get tau particles and their four-momenta
  TLorentzVector tau_pos_p4 = TLorentzVector(0., 0., 0., 0.);
  TLorentzVector tau_neg_p4 = TLorentzVector(0., 0., 0., 0.);

  int tauCount = 0;
  for (const xAOD::TruthParticle *tau : *truthTauParticles) {
    // Skip taus without parents
    if (tau->nParents() == 0) {
      continue;
    }

    // Check if the parent is a Higgs boson
    if (!tau->parent(0)->isHiggs()) {
      continue;
    }

    tauCount += 1;
    if (tau->pdgId() == -TAU) {
      tau_pos_p4 = tau->p4();
    }
    if (tau->pdgId() == TAU) {
      tau_neg_p4 = tau->p4();
    }
  }

  if (tauCount != 2) {
    ANA_MSG_INFO("Tau count is not 2. Excluding event.");
    tau_pos_p4 = TLorentzVector(0., 0., 0., 0.);
    tau_neg_p4 = TLorentzVector(0., 0., 0., 0.);
  }

  return StatusCode::SUCCESS;
}

StatusCode MyxAODAnalysis ::finalize() {
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