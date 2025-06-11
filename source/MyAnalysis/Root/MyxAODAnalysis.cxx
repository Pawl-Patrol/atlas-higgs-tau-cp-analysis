#include "AsgMessaging/MessageCheck.h"
#include <MyAnalysis/MyxAODAnalysis.h>
#include <TLorentzVector.h>

MyxAODAnalysis ::MyxAODAnalysis(const std::string &name,
                                ISvcLocator *pSvcLocator)
    : EL::AnaAlgorithm(name, pSvcLocator) {}

StatusCode MyxAODAnalysis ::initialize() {
  ANA_MSG_DEBUG("Initializing");

  // Setup histograms
  ANA_CHECK(book(TH1F("phi_CP_tau_pi", "phi_CP_tau_pi", 50, 0, 2 * M_PI)));
  ANA_CHECK(
      book(TH1F("phi_CP_neutrino_pi", "phi_CP_neutrino_pi", 50, 0, 2 * M_PI)));
  ANA_CHECK(book(TH1F("phi_CP_pion", "phi_CP_pion", 50, 0, 2 * M_PI)));

  return StatusCode::SUCCESS;
}

TVector3 getPerpendicularComponent(const TVector3 &vec1, const TVector3 &vec2) {
  TVector3 unit_vec2 = vec2.Unit();
  return vec1 - (vec1.Dot(unit_vec2)) * unit_vec2;
}

StatusCode MyxAODAnalysis::execute() {
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

  // If not found, try to get TruthBosons
  if (!result.isSuccess() || truthHiggsWithDecayParticles->empty()) {
    ANA_CHECK(evtStore()->retrieve(truthHiggsWithDecayParticles,
                                   "TruthBosonsWithDecayParticles"));
  }

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

  const xAOD::TruthParticleContainer *truthTausWithDecayParticles = nullptr;
  ANA_CHECK(evtStore()->retrieve(truthTausWithDecayParticles,
                                 "TruthTausWithDecayParticles"));

  for (const xAOD::TruthParticle *particle : *truthTausWithDecayParticles) {
    if (particle->nParents() == 0) {
      continue;
    }

    if (particle->pdgId() == PIPLUS && particle->parent(0)->pdgId() == -TAU) {
      pPionPos = particle;
      pionCount++;
    } else if (particle->pdgId() == PIMINUS &&
               particle->parent(0)->pdgId() == TAU) {
      pPionNeg = particle;
      pionCount++;
    }

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

  // Using pion impact parameter/momentum planes
  TLorentzVector impact_param_pos_pion = TLorentzVector(
      getPerpendicularComponent(pPionPos->prodVtx()->v4().Vect() -
                                    pTauNeg->prodVtx()->v4().Vect(),
                                pPionPos->p4().Vect())
          .Unit(),
      0.);
  TLorentzVector impact_param_neg_pion = TLorentzVector(
      getPerpendicularComponent(pPionNeg->prodVtx()->v4().Vect() -
                                    pTauNeg->prodVtx()->v4().Vect(),
                                pPionNeg->p4().Vect())
          .Unit(),
      0.);

  TLorentzVector pion_pos_p4 = TLorentzVector(pPionPos->p4());
  TLorentzVector pion_neg_p4 = TLorentzVector(pPionNeg->p4());

  // Boost into CMF of the pions
  TVector3 cmfBoostVector = (pion_pos_p4 + pion_neg_p4).BoostVector();
  impact_param_pos_pion.Boost(-cmfBoostVector);
  impact_param_neg_pion.Boost(-cmfBoostVector);
  pion_pos_p4.Boost(-cmfBoostVector);
  pion_neg_p4.Boost(-cmfBoostVector);

  // Get the impact parameter component perpendicular to the momentum
  TVector3 impact_param_pos_perp = getPerpendicularComponent(
      impact_param_pos_pion.Vect(), pion_pos_p4.Vect());
  TVector3 impact_param_neg_perp = getPerpendicularComponent(
      impact_param_neg_pion.Vect(), pion_neg_p4.Vect());

  double angleO_pion =
      pion_neg_p4.Vect().Unit() *
      (impact_param_pos_perp.Cross(impact_param_neg_perp).Unit());

  if (angleO_pion >= 0) {
    m_phiCPPion =
        acos(impact_param_pos_perp.Unit() * impact_param_neg_perp.Unit());
  } else {
    m_phiCPPion = 2 * M_PI - acos(impact_param_pos_perp.Unit() *
                                  impact_param_neg_perp.Unit());
  }

  // Boost all four-momenta to the Higgs boson rest frame = taus CMF
  TVector3 higgsBoostVector = pHiggs->p4().BoostVector();
  TLorentzVector tauPosP4 = pTauPos->p4();
  TLorentzVector tauNegP4 = pTauNeg->p4();
  TLorentzVector piPosP4 = pPionPos->p4();
  TLorentzVector piNegP4 = pPionNeg->p4();
  TLorentzVector neutriP4 = pNeutrino->p4();
  TLorentzVector antiNeutriP4 = pAntiNeutrino->p4();

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
