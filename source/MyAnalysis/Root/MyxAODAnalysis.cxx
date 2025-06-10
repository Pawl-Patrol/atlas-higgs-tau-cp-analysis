#include <AsgMessaging/MessageCheck.h>
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

  // Setup histograms
  ANA_CHECK(book(TH1F("phi_CP_tau_pi", "phi_CP_tau_pi", 50, 0, 2 * M_PI)));
  ANA_CHECK(
      book(TH1F("phi_CP_neutrino_pi", "phi_CP_neutrino_pi", 50, 0, 2 * M_PI)));
  ANA_CHECK(book(TH1F("phi_CP_pion", "phi_CP_pion", 50, 0, 2 * M_PI)));

  // Setup tree
  // ANA_CHECK(book(TTree("truth_tau_analysis", "Zee analysis ntuple")));
  // TTree *myTree = tree("truth_tau_analysis");
  // myTree->Branch("run_number", &m_runNumber);
  // myTree->Branch("event_number", &m_eventNumber);
  // myTree->Branch("phi_CP_tau_pi", &m_phiCP);
  // myTree->Branch("phi_CP_neutrino_pi", &m_phiCPNeutri);

  return StatusCode::SUCCESS;
}

TVector3 getPerpendicularComponent(const TVector3 &vec1, const TVector3 &vec2) {
  TVector3 unit_vec2 = vec2.Unit();
  return vec1 - (vec1.Dot(unit_vec2)) * unit_vec2;
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
  m_phiCP = -99999.;
  m_phiCPNeutri = -99999.;

  bool switchBSM = false;
  const xAOD::TruthParticleContainer *truthHiggs = nullptr;
  ANA_CHECK(evtStore()->retrieve(truthHiggs, "TruthBoson"));
  if (evtStore()->contains<xAOD::TruthParticleContainer>("TruthBSM")) {
    if (evtStore()->retrieve(truthHiggs, "TruthBSM")) {
      if (!truthHiggs->empty()) {
        switchBSM = true;
      }
    }
  }

  if (switchBSM == true) {
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
  TLorentzVector higgsP4;
  TLorentzVector higgsProdVtx = TLorentzVector(0., 0., 0., 0.);

  for (const xAOD::TruthParticle *higgs : *truthHiggs) {
    higgsP4 = higgs->p4();
    if (higgs->hasProdVtx()) {
      higgsProdVtx = higgs->prodVtx()->v4();
    }
  }

  // Get tau particles and their four-momenta
  TLorentzVector tauPosP4 = TLorentzVector(0., 0., 0., 0.);
  TLorentzVector tauNegP4 = TLorentzVector(0., 0., 0., 0.);

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
      tauPosP4 = tau->p4();
    }
    if (tau->pdgId() == TAU) {
      tauNegP4 = tau->p4();
    }
  }

  // technically tauCount == 1 should be impossible, but we check for it anyways
  if (tauCount != 2) {
    ANA_MSG_INFO("Tau count is not 2. Excluding event.");
    // tree("truth_tau_analysis")->Fill();
    return StatusCode::SUCCESS;
  }

  // Check for charged pions in the decay products of taus
  TLorentzVector piPosP4 = TLorentzVector(0., 0., 0., 0.);
  TLorentzVector piNegP4 = TLorentzVector(0., 0., 0., 0.);
  TLorentzVector piPosProdVtx = TLorentzVector(0., 0., 0., 0.);
  TLorentzVector piNegProdVtx = TLorentzVector(0., 0., 0., 0.);

  int chPionCount = 0;
  for (const xAOD::TruthParticle *part : *truthTausWithDecayParticles) {
    if (abs(part->pdgId()) == PIPLUS) {
      chPionCount += 1;
    }
  }

  if (chPionCount != 2) {
    ANA_MSG_INFO("Charged pion count is not 2. Excluding event.");
    // tree("truth_tau_analysis")->Fill();
    return StatusCode::SUCCESS;
  }

  for (const xAOD::TruthParticle *particle : *truthTausWithDecayParticles) {
    // Skip particles without parents
    if (particle->nParents() == 0) {
      continue;
    }

    if (particle->pdgId() == PIPLUS && particle->parent(0)->pdgId() == -TAU) {
      piPosP4 = particle->p4();
      piPosProdVtx = particle->prodVtx()->v4();
    }
    if (particle->pdgId() == PIMINUS && particle->parent(0)->pdgId() == TAU) {
      piNegP4 = particle->p4();
      piNegProdVtx = particle->prodVtx()->v4();
    }
  }

  // Check for neutrinos in the decay products of taus
  TLorentzVector neutriP4 = TLorentzVector(0., 0., 0., 0.);
  TLorentzVector antiNeutriP4 = TLorentzVector(0., 0., 0., 0.);

  int chNuTauCount = 0;
  for (auto part : *truthTausWithDecayParticles) {
    if (abs(part->pdgId()) == NU_TAU) {
      chNuTauCount += 1;
    }
  }

  if (chNuTauCount != 2) {
    ANA_MSG_INFO("Charged neutrino count is not 2. Excluding event.");
    // tree("truth_tau_analysis")->Fill();
    return StatusCode::SUCCESS;
  }

  for (const xAOD::TruthParticle *particle : *truthTausWithDecayParticles) {
    // Skip particles without parents
    if (particle->nParents() == 0) {
      continue;
    }

    if (particle->pdgId() == -NU_TAU && particle->parent(0)->pdgId() == -TAU) {
      antiNeutriP4 = particle->p4();
    }
    if (particle->pdgId() == NU_TAU && particle->parent(0)->pdgId() == TAU) {
      neutriP4 = particle->p4();
    }
  }

  ANA_MSG_INFO("Found pi+ pi- decay");

  // Using pion impact parameter/momentum planes
  TLorentzVector dca_pos_pion = TLorentzVector(
      getPerpendicularComponent(piPosProdVtx.Vect() - higgsProdVtx.Vect(),
                                piPosP4.Vect()),
      0.);
  TLorentzVector dca_neg_pion = TLorentzVector(
      getPerpendicularComponent(piNegProdVtx.Vect() - higgsProdVtx.Vect(),
                                piNegP4.Vect()),
      0.);
  TLorentzVector pion_neg_p4 = TLorentzVector(piNegP4);

  // Boost into CMF of the pions
  TVector3 cmfBoostVector = (piPosP4 + piNegP4).BoostVector();
  dca_pos_pion.Boost(-cmfBoostVector);
  dca_neg_pion.Boost(-cmfBoostVector);
  pion_neg_p4.Boost(-cmfBoostVector);

  // Get the impact parameter component perpendicular to the momentum
  TVector3 dca_pos_perp =
      getPerpendicularComponent(dca_pos_pion.Vect(), piPosP4.Vect());
  TVector3 dca_neg_perp =
      getPerpendicularComponent(dca_neg_pion.Vect(), piNegP4.Vect());

  double angleO_pion =
      pion_neg_p4.Vect().Unit() * (dca_pos_perp.Cross(dca_neg_perp).Unit());

  if (angleO_pion >= 0) {
    m_phiCPPion = acos(dca_pos_perp.Unit() * dca_neg_perp.Unit());
  } else {
    m_phiCPPion = 2 * M_PI - acos(dca_pos_perp.Unit() * dca_neg_perp.Unit());
  }

  // Boost all four-momenta to the Higgs boson rest frame = taus CMF
  TVector3 higgsBoostVector = higgsP4.BoostVector();
  tauPosP4.Boost(-higgsBoostVector);
  tauNegP4.Boost(-higgsBoostVector);
  piPosP4.Boost(-higgsBoostVector);
  piNegP4.Boost(-higgsBoostVector);
  neutriP4.Boost(-higgsBoostVector);
  antiNeutriP4.Boost(-higgsBoostVector);

  // Tau/Pion decay planes
  TVector3 norm_vec_pos = (tauPosP4.Vect().Cross(piPosP4.Vect())).Unit();
  TVector3 norm_vec_neg = (tauNegP4.Vect().Cross(piNegP4.Vect())).Unit();

  double angleO = piNegP4.Vect().Unit() * (norm_vec_pos.Cross(norm_vec_neg));

  if (angleO >= 0) {
    m_phiCP = acos(norm_vec_pos * norm_vec_neg);
  } else {
    m_phiCP = 2 * M_PI - acos(norm_vec_pos * norm_vec_neg);
  }

  // Neutrino/Pion decay planes
  TVector3 norm_vec_pos_aneutri =
      (antiNeutriP4.Vect().Cross(piPosP4.Vect())).Unit();
  TVector3 norm_vec_neg_neutri = (neutriP4.Vect().Cross(piNegP4.Vect())).Unit();

  double angleO_neutri =
      piNegP4.Vect().Unit() * (norm_vec_pos_aneutri.Cross(norm_vec_neg_neutri));

  if (angleO_neutri >= 0) {
    m_phiCPNeutri = acos(norm_vec_pos_aneutri * norm_vec_neg_neutri);
  } else {
    m_phiCPNeutri = 2 * M_PI - acos(norm_vec_pos_aneutri * norm_vec_neg_neutri);
  }

  hist("phi_CP_tau_pi")->Fill(m_phiCP);
  hist("phi_CP_neutrino_pi")->Fill(m_phiCPNeutri);
  hist("phi_CP_pion")->Fill(m_phiCPPion);
  // tree("truth_tau_analysis")->Fill();

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