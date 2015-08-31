#include <vector>
#include <map>
#include <math.h>

// event loop
#include <EventLoop/Job.h>
#include <EventLoop/Worker.h>
#include "EventLoop/OutputStream.h"
#include "xAODEventInfo/EventInfo.h"

#include <MyAnalysis/WriteTree.h>

#include <xAODJet/FastJetLinkBase.h>
#include "fastjet/JetDefinition.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/Selector.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"
#include <fastjet/tools/Subtractor.hh>
#include "xAODCaloEvent/CaloClusterContainer.h"

#include "xAODCore/ShallowAuxContainer.h"
#include "xAODJet/JetConstituentVector.h"

// EDM
#include "xAODCaloEvent/CaloCluster.h"
#include "xAODTruth/TruthParticleContainer.h"
#include "xAODTruth/TruthEventContainer.h"
#include "xAODTruth/TruthEvent.h"

// Infrastructure include(s):
#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TEvent.h"

// xAH includes
#include "xAODAnaHelpers/HelperFunctions.h"
#include "xAODAnaHelpers/tools/ReturnCheck.h"

namespace HF = HelperFunctions;

#define ARRAY_INIT {}
// this is needed to distribute the algorithm to the workers
ClassImp(WriteTree)

WriteTree :: WriteTree () :
  m_tree(new TTree("oTree", "output tree")),
  m_eventNumber(-999.0),
  m_NPV(-99),
  m_mu(-99),

  m_jvoro0pt{ARRAY_INIT},
  m_jvoro0eta{ARRAY_INIT},
  m_jvoro0phi{ARRAY_INIT},
  m_jvoro0mass{ARRAY_INIT},
  m_jvoro0width{ARRAY_INIT},
  m_jvoro0mindr{ARRAY_INIT},
  m_tjvoro0pt{ARRAY_INIT},
  m_tjvoro0eta{ARRAY_INIT},
  m_tjvoro0phi{ARRAY_INIT},
  m_tjvoro0mass{ARRAY_INIT},
  m_tjvoro0width{ARRAY_INIT},
  m_tjvoro0mindr{ARRAY_INIT}
{}

EL::StatusCode WriteTree :: setupJob (EL::Job& job)
{
  // let's initialize the algorithm to use the xAODRootAccess package
  job.useXAOD ();
  xAOD::Init("VoronoiWeights").ignore(); // call before opening first file

  /* https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/EventLoop#Creating_output_n_tuples */
  // output file for tree
  EL::OutputStream out ("outputTree");
  job.outputAdd (out);

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode WriteTree :: histInitialize () {
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode WriteTree :: fileExecute () {return EL::StatusCode::SUCCESS;}

EL::StatusCode WriteTree :: changeInput (bool /*firstFile*/) {return EL::StatusCode::SUCCESS;}

EL::StatusCode WriteTree :: initialize ()
{
  m_event = wk()->xaodEvent();
  m_store = wk()->xaodStore();
  // as a check, let's see the number of events in our xAOD
  Info("initialize()", "Number of events = %lli", m_event->getEntries() ); // print long long int

  TFile *file = wk()->getOutputFile ("outputTree");
  m_tree->SetDirectory (file);

  m_tree->Branch ("event_number",              &m_eventNumber, "event_number/I");
  m_tree->Branch ("NPV",              &m_NPV, "NPV/I");
  m_tree->Branch ("mu",              &m_mu, "mu/I");

  m_tree->Branch("jvoro0pt","std::vector<float>", &m_jvoro0pt);
  m_tree->Branch("jvoro0eta","std::vector<float>", &m_jvoro0eta);
  m_tree->Branch("jvoro0phi","std::vector<float>", &m_jvoro0phi);
  m_tree->Branch("jvoro0mass","std::vector<float>", &m_jvoro0mass);
  m_tree->Branch("jvoro0width","std::vector<float>", &m_jvoro0width);
  m_tree->Branch("jvoro0mindr","std::vector<float>", &m_jvoro0mindr);
  m_tree->Branch("tjvoro0pt","std::vector<float>", &m_tjvoro0pt);
  m_tree->Branch("tjvoro0eta","std::vector<float>", &m_tjvoro0eta);
  m_tree->Branch("tjvoro0phi","std::vector<float>", &m_tjvoro0phi);
  m_tree->Branch("tjvoro0m","std::vector<float>", &m_tjvoro0mass);
  m_tree->Branch("tjvoro0width","std::vector<float>", &m_tjvoro0width);
  m_tree->Branch("tjvoro0mindr","std::vector<float>", &m_tjvoro0mindr);

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode WriteTree :: execute ()
{
  const char* APP_NAME = "WriteTree::execute()";
  const xAOD::EventInfo*                        eventInfo     (nullptr);
  const xAOD::JetContainer*                     in_jets       (nullptr);
  const xAOD::JetContainer*                     truth_jets    (nullptr);
  const xAOD::JetContainer*                     voronoi_jets  (nullptr);
  const xAOD::VertexContainer*                  vertices      (nullptr);

  // start grabbing all the containers that we can
  RETURN_CHECK("VoronoiWeights::execute()", HF::retrieve(eventInfo,    m_eventInfo,        m_event, m_store, m_debug), "Could not get the EventInfo container.");
  if(!m_jets.empty()) RETURN_CHECK("VoronoiWeights::execute()", HF::retrieve(in_jets,     m_jets,       m_event, m_store, m_debug), "Could not get the jets container.");
  if(!m_truth_jets.empty()) RETURN_CHECK("VoronoiWeights::execute()", HF::retrieve(truth_jets,    m_truth_jets,       m_event, m_store, m_debug), "Could not get the truth jets container.");
  if(!m_voronoi_jets.empty()) RETURN_CHECK("VoronoiWeights::execute()", HF::retrieve(voronoi_jets,    m_voronoi_jets,       m_event, m_store, m_debug), "Could not get the voronoi jets container.");
  if(!m_vertices.empty()) RETURN_CHECK("VoronoiWeights::execute()", HF::retrieve(vertices,    m_vertices,       m_event, m_store, m_debug), "Could not get the vertices container.");

  m_eventNumber = eventInfo->eventNumber();

  if(FillJetVars(voronoi_jets,truth_jets,m_jvoro0pt,m_jvoro0eta,m_jvoro0phi,m_jvoro0mass,m_jvoro0width,m_jvoro0mindr,m_tjvoro0pt,m_tjvoro0eta,m_tjvoro0phi,m_tjvoro0mass,m_tjvoro0width,m_tjvoro0mindr) != EL::StatusCode::SUCCESS)
    Error(APP_NAME,"Error in FillJetVars");

  m_mu = eventInfo->averageInteractionsPerCrossing();
  m_NPV = 0;
  for ( auto *ivert : *vertices ){
    if ( (ivert)->nTrackParticles() >= 2 ) ++m_NPV;
  }
  // fill in all variables
  m_tree->Fill();

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode WriteTree :: postExecute () {return EL::StatusCode::SUCCESS;}

EL::StatusCode WriteTree :: finalize ()
{
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode WriteTree :: histFinalize () {
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode WriteTree :: FillJetVars(const DataVector<xAOD::Jet_v1>* jets,
                                        const DataVector<xAOD::Jet_v1>* tjets,
                                        std::vector<float>& jetpt,
                                        std::vector<float>& jeteta,
                                        std::vector<float>& jetphi,
                                        std::vector<float>& jetmass,
                                        std::vector<float>& jetwidth,
                                        std::vector<float>& jetmindr,
                                        std::vector<float>& tjetpt,
                                        std::vector<float>& tjeteta,
                                        std::vector<float>& tjetphi,
                                        std::vector<float>& tjetmass,
                                        std::vector<float>& tjetwidth,
                                        std::vector<float>& tjetmindr){

  jetpt.clear();
  jeteta.clear();
  jetphi.clear();
  jetmass.clear();
  jetwidth.clear();
  jetmindr.clear();
  tjetpt.clear();
  tjeteta.clear();
  tjetphi.clear();
  tjetmass.clear();
  tjetwidth.clear();
  tjetmindr.clear();

  DataVector<xAOD::Jet_v1> sorted_truth_jets = HF::sort_container_pt(tjets);
  static SG::AuxElement::ConstAccessor< int > truth_match_i("truth_match_i");
  static SG::AuxElement::ConstAccessor< float > minDR("minDR");

  for(const auto jet: *jets){
    jetpt.push_back(jet->pt());
    jeteta.push_back(jet->eta());
    jetphi.push_back(jet->phi());
    jetmass.push_back(jet->m());
    jetwidth.push_back(-99);
    jetmindr.push_back(-99);

    int truthmatch = truth_match_i(*jet);
    if(truthmatch > -1){
      auto tjet = sorted_truth_jets.at(truthmatch);
      tjetpt.push_back(tjet->pt());
      tjeteta.push_back(tjet->eta());
      tjetphi.push_back(tjet->phi());
      tjetmass.push_back(tjet->m());
      tjetwidth.push_back(-99);
      tjetmindr.push_back(minDR(*tjet));
    }
    else{
      tjetpt.push_back(-99);
      tjeteta.push_back(-99);
      tjetphi.push_back(-99);
      tjetmass.push_back(-99);
      tjetwidth.push_back(-99);
      tjetmindr.push_back(-99);
    }
  }
  return EL::StatusCode::SUCCESS;
}
