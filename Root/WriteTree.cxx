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
//#include "xAODAnaHelpers/tools/ReturnCheck.h"

namespace HF = HelperFunctions;

#define ARRAY_INIT {}
// this is needed to distribute the algorithm to the workers
ClassImp(WriteTree)

WriteTree :: WriteTree () :
  m_tree(new TTree("oTree", "output tree")),
  m_eventNumber(-999.0),
  m_eventWeight(-999.0),
  m_NPV(-99),
  m_mu(-99),
  m_rho(-99),

  m_tjpt{ARRAY_INIT},
  m_tjeta{ARRAY_INIT},
  m_tjphi{ARRAY_INIT},
  m_tjmass{ARRAY_INIT},
  m_tjwidth{ARRAY_INIT},
  m_tjmindr{ARRAY_INIT},

  m_jvoro0pt{ARRAY_INIT},
  m_jvoro0eta{ARRAY_INIT},
  m_jvoro0phi{ARRAY_INIT},
  m_jvoro0mass{ARRAY_INIT},
  m_jvoro0width{ARRAY_INIT},
  m_jvoro0mindr{ARRAY_INIT},
  m_jvoro0isPU{ARRAY_INIT},
  m_jvoro0cl0pt{ARRAY_INIT},
  m_tjvoro0pt{ARRAY_INIT},
  m_tjvoro0eta{ARRAY_INIT},
  m_tjvoro0phi{ARRAY_INIT},
  m_tjvoro0mass{ARRAY_INIT},
  m_tjvoro0width{ARRAY_INIT},
  m_tjvoro0mindr{ARRAY_INIT},

  m_jvoro1pt{ARRAY_INIT},
  m_jvoro1eta{ARRAY_INIT},
  m_jvoro1phi{ARRAY_INIT},
  m_jvoro1mass{ARRAY_INIT},
  m_jvoro1width{ARRAY_INIT},
  m_jvoro1mindr{ARRAY_INIT},
  m_jvoro1isPU{ARRAY_INIT},
  m_jvoro1cl0pt{ARRAY_INIT},
  m_tjvoro1pt{ARRAY_INIT},
  m_tjvoro1eta{ARRAY_INIT},
  m_tjvoro1phi{ARRAY_INIT},
  m_tjvoro1mass{ARRAY_INIT},
  m_tjvoro1width{ARRAY_INIT},
  m_tjvoro1mindr{ARRAY_INIT},

  m_jvorospt{ARRAY_INIT},
  m_jvoroseta{ARRAY_INIT},
  m_jvorosphi{ARRAY_INIT},
  m_jvorosmass{ARRAY_INIT},
  m_jvoroswidth{ARRAY_INIT},
  m_jvorosmindr{ARRAY_INIT},
  m_jvorosisPU{ARRAY_INIT},
  m_jvoroscl0pt{ARRAY_INIT},
  m_tjvorospt{ARRAY_INIT},
  m_tjvoroseta{ARRAY_INIT},
  m_tjvorosphi{ARRAY_INIT},
  m_tjvorosmass{ARRAY_INIT},
  m_tjvoroswidth{ARRAY_INIT},
  m_tjvorosmindr{ARRAY_INIT},

  m_jnoarea0pt{ARRAY_INIT},
  m_jnoarea0eta{ARRAY_INIT},
  m_jnoarea0phi{ARRAY_INIT},
  m_jnoarea0mass{ARRAY_INIT},
  m_jnoarea0width{ARRAY_INIT},
  m_jnoarea0mindr{ARRAY_INIT},
  m_jnoarea0isPU{ARRAY_INIT},
  m_jnoarea0cl0pt{ARRAY_INIT},
  m_tjnoarea0pt{ARRAY_INIT},
  m_tjnoarea0eta{ARRAY_INIT},
  m_tjnoarea0phi{ARRAY_INIT},
  m_tjnoarea0mass{ARRAY_INIT},
  m_tjnoarea0width{ARRAY_INIT},
  m_tjnoarea0mindr{ARRAY_INIT},

  m_j0pt{ARRAY_INIT},
  m_j0eta{ARRAY_INIT},
  m_j0phi{ARRAY_INIT},
  m_j0mass{ARRAY_INIT},
  m_j0width{ARRAY_INIT},
  m_j0mindr{ARRAY_INIT},
  m_j0isPU{ARRAY_INIT},
  m_j0cl0pt{ARRAY_INIT},
  m_tj0pt{ARRAY_INIT},
  m_tj0eta{ARRAY_INIT},
  m_tj0phi{ARRAY_INIT},
  m_tj0mass{ARRAY_INIT},
  m_tj0width{ARRAY_INIT},
  m_tj0mindr{ARRAY_INIT}
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
  //if(m_doLC) m_jets = "AntiKt4LCTopoJets";
  //else m_jets = "AntiKt4EMTopoJets";

  // as a check, let's see the number of events in our xAOD
  Info("initialize()", "Number of events = %lli", m_event->getEntries() ); // print long long int

  TFile *file = wk()->getOutputFile ("outputTree");
  m_tree->SetDirectory (file);

  m_tree->Branch ("event_number",              &m_eventNumber, "event_number/I");
  m_tree->Branch ("event_weight",              &m_eventWeight, "event_weight/F");
  m_tree->Branch ("NPV",              &m_NPV, "NPV/I");
  m_tree->Branch ("mu",              &m_mu, "mu/I");
  m_tree->Branch ("rho",              &m_rho, "rho/F");

  m_tree->Branch("tjpt","std::vector<float>", &m_tjpt);
  m_tree->Branch("tjeta","std::vector<float>", &m_tjeta);
  m_tree->Branch("tjphi","std::vector<float>", &m_tjphi);
  m_tree->Branch("tjm","std::vector<float>", &m_tjmass);
  m_tree->Branch("tjwidth","std::vector<float>", &m_tjwidth);
  m_tree->Branch("tjmindr","std::vector<float>", &m_tjmindr);

  m_tree->Branch("jvoro0pt","std::vector<float>", &m_jvoro0pt);
  m_tree->Branch("jvoro0eta","std::vector<float>", &m_jvoro0eta);
  m_tree->Branch("jvoro0phi","std::vector<float>", &m_jvoro0phi);
  m_tree->Branch("jvoro0mass","std::vector<float>", &m_jvoro0mass);
  m_tree->Branch("jvoro0width","std::vector<float>", &m_jvoro0width);
  m_tree->Branch("jvoro0mindr","std::vector<float>", &m_jvoro0mindr);
  m_tree->Branch("jvoro0isPU","std::vector<bool>", &m_jvoro0isPU);
  m_tree->Branch("jvoro0cl0pt","std::vector<float>", &m_jvoro0cl0pt);
  m_tree->Branch("tjvoro0pt","std::vector<float>", &m_tjvoro0pt);
  m_tree->Branch("tjvoro0eta","std::vector<float>", &m_tjvoro0eta);
  m_tree->Branch("tjvoro0phi","std::vector<float>", &m_tjvoro0phi);
  m_tree->Branch("tjvoro0m","std::vector<float>", &m_tjvoro0mass);
  m_tree->Branch("tjvoro0width","std::vector<float>", &m_tjvoro0width);
  m_tree->Branch("tjvoro0mindr","std::vector<float>", &m_tjvoro0mindr);

  m_tree->Branch("jvoro1pt","std::vector<float>", &m_jvoro1pt);
  m_tree->Branch("jvoro1eta","std::vector<float>", &m_jvoro1eta);
  m_tree->Branch("jvoro1phi","std::vector<float>", &m_jvoro1phi);
  m_tree->Branch("jvoro1mass","std::vector<float>", &m_jvoro1mass);
  m_tree->Branch("jvoro1width","std::vector<float>", &m_jvoro1width);
  m_tree->Branch("jvoro1mindr","std::vector<float>", &m_jvoro1mindr);
  m_tree->Branch("jvoro1isPU","std::vector<bool>", &m_jvoro1isPU);
  m_tree->Branch("jvoro1cl0pt","std::vector<float>", &m_jvoro1cl0pt);
  m_tree->Branch("tjvoro1pt","std::vector<float>", &m_tjvoro1pt);
  m_tree->Branch("tjvoro1eta","std::vector<float>", &m_tjvoro1eta);
  m_tree->Branch("tjvoro1phi","std::vector<float>", &m_tjvoro1phi);
  m_tree->Branch("tjvoro1m","std::vector<float>", &m_tjvoro1mass);
  m_tree->Branch("tjvoro1width","std::vector<float>", &m_tjvoro1width);
  m_tree->Branch("tjvoro1mindr","std::vector<float>", &m_tjvoro1mindr);

  m_tree->Branch("jvorospt","std::vector<float>", &m_jvorospt);
  m_tree->Branch("jvoroseta","std::vector<float>", &m_jvoroseta);
  m_tree->Branch("jvorosphi","std::vector<float>", &m_jvorosphi);
  m_tree->Branch("jvorosmass","std::vector<float>", &m_jvorosmass);
  m_tree->Branch("jvoroswidth","std::vector<float>", &m_jvoroswidth);
  m_tree->Branch("jvorosmindr","std::vector<float>", &m_jvorosmindr);
  m_tree->Branch("jvorosisPU","std::vector<bool>", &m_jvorosisPU);
  m_tree->Branch("jvoroscl0pt","std::vector<float>", &m_jvoroscl0pt);
  m_tree->Branch("tjvorospt","std::vector<float>", &m_tjvorospt);
  m_tree->Branch("tjvoroseta","std::vector<float>", &m_tjvoroseta);
  m_tree->Branch("tjvorosphi","std::vector<float>", &m_tjvorosphi);
  m_tree->Branch("tjvorosm","std::vector<float>", &m_tjvorosmass);
  m_tree->Branch("tjvoroswidth","std::vector<float>", &m_tjvoroswidth);
  m_tree->Branch("tjvorosmindr","std::vector<float>", &m_tjvorosmindr);

  m_tree->Branch("jnoarea0pt","std::vector<float>", &m_jnoarea0pt);
  m_tree->Branch("jnoarea0eta","std::vector<float>", &m_jnoarea0eta);
  m_tree->Branch("jnoarea0phi","std::vector<float>", &m_jnoarea0phi);
  m_tree->Branch("jnoarea0mass","std::vector<float>", &m_jnoarea0mass);
  m_tree->Branch("jnoarea0width","std::vector<float>", &m_jnoarea0width);
  m_tree->Branch("jnoarea0mindr","std::vector<float>", &m_jnoarea0mindr);
  m_tree->Branch("jnoarea0isPU","std::vector<bool>", &m_jnoarea0isPU);
  m_tree->Branch("jnoarea0cl0pt","std::vector<float>", &m_jnoarea0cl0pt);
  m_tree->Branch("tjnoarea0pt","std::vector<float>", &m_tjnoarea0pt);
  m_tree->Branch("tjnoarea0eta","std::vector<float>", &m_tjnoarea0eta);
  m_tree->Branch("tjnoarea0phi","std::vector<float>", &m_tjnoarea0phi);
  m_tree->Branch("tjnoarea0m","std::vector<float>", &m_tjnoarea0mass);
  m_tree->Branch("tjnoarea0width","std::vector<float>", &m_tjnoarea0width);
  m_tree->Branch("tjnoarea0mindr","std::vector<float>", &m_tjnoarea0mindr);

  m_tree->Branch("j0pt","std::vector<float>", &m_j0pt);
  m_tree->Branch("j0eta","std::vector<float>", &m_j0eta);
  m_tree->Branch("j0phi","std::vector<float>", &m_j0phi);
  m_tree->Branch("j0mass","std::vector<float>", &m_j0mass);
  m_tree->Branch("j0width","std::vector<float>", &m_j0width);
  m_tree->Branch("j0mindr","std::vector<float>", &m_j0mindr);
  m_tree->Branch("j0isPU","std::vector<bool>", &m_j0isPU);
  m_tree->Branch("j0cl0pt","std::vector<float>", &m_j0cl0pt);
  m_tree->Branch("tj0pt","std::vector<float>", &m_tj0pt);
  m_tree->Branch("tj0eta","std::vector<float>", &m_tj0eta);
  m_tree->Branch("tj0phi","std::vector<float>", &m_tj0phi);
  m_tree->Branch("tj0m","std::vector<float>", &m_tj0mass);
  m_tree->Branch("tj0width","std::vector<float>", &m_tj0width);
  m_tree->Branch("tj0mindr","std::vector<float>", &m_tj0mindr);

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode WriteTree :: execute ()
{
  const char* APP_NAME = "WriteTree::execute()";
  const xAOD::EventInfo*                        eventInfo     (nullptr);
  const xAOD::JetContainer*                     in_jets       (nullptr);
  const xAOD::JetContainer*                     truth_jets    (nullptr);
  const xAOD::JetContainer*                     voronoi0_jets  (nullptr);
  const xAOD::JetContainer*                     voronoi1_jets  (nullptr);
  const xAOD::JetContainer*                     voronois_jets  (nullptr);
  const xAOD::VertexContainer*                  vertices      (nullptr);

  // start grabbing all the containers that we can
  ANA_CHECK(HF::retrieve(eventInfo,    m_eventInfo,        m_event, m_store, m_debug));
  //m_jets = "AntiKt4EMTopoJets";
  if(!m_jets.empty()) ANA_CHECK(HF::retrieve(in_jets,     m_jets,       m_event, m_store, m_debug));
  if(!m_truth_jets.empty()) ANA_CHECK(HF::retrieve(truth_jets,    m_truth_jets,       m_event, m_store, m_debug));
  if(!m_voronoi0_jets.empty()) ANA_CHECK(HF::retrieve(voronoi0_jets,    m_voronoi0_jets,       m_event, m_store, m_debug));
  if(!m_voronoi1_jets.empty()) ANA_CHECK(HF::retrieve(voronoi1_jets,    m_voronoi1_jets,       m_event, m_store, m_debug));
  if(!m_voronois_jets.empty()) ANA_CHECK(HF::retrieve(voronois_jets,    m_voronois_jets,       m_event, m_store, m_debug));
  if(!m_vertices.empty()) ANA_CHECK(HF::retrieve(vertices,    m_vertices,       m_event, m_store, m_debug));

  if(m_debug){
    std::cout << "Voronoi0 Jets" << std::endl;
    for(const auto jet:*voronoi0_jets){
      std::cout << jet->pt() << std::endl;
    }
  }
  m_eventNumber = eventInfo->eventNumber();
  m_eventWeight = eventInfo->mcEventWeight();

  if(FillTJetVars(truth_jets,m_tjpt,m_tjeta,m_tjphi,m_tjmass,m_tjwidth,m_tjmindr) != EL::StatusCode::SUCCESS)
    Error(APP_NAME,"Error in FillTJetVars");

  if(FillJetVars(voronoi0_jets,truth_jets,m_jvoro0pt,m_jvoro0eta,m_jvoro0phi,m_jvoro0mass,m_jvoro0width,m_jvoro0mindr,m_jvoro0isPU,m_jvoro0cl0pt,m_tjvoro0pt,m_tjvoro0eta,m_tjvoro0phi,m_tjvoro0mass,m_tjvoro0width,m_tjvoro0mindr) != EL::StatusCode::SUCCESS)
    Error(APP_NAME,"Error in FillJetVars");
  if(FillJetVars(voronoi1_jets,truth_jets,m_jvoro1pt,m_jvoro1eta,m_jvoro1phi,m_jvoro1mass,m_jvoro1width,m_jvoro1mindr,m_jvoro1isPU,m_jvoro1cl0pt,m_tjvoro1pt,m_tjvoro1eta,m_tjvoro1phi,m_tjvoro1mass,m_tjvoro1width,m_tjvoro1mindr) != EL::StatusCode::SUCCESS)
    Error(APP_NAME,"Error in FillJetVars");
  if(FillJetVars(voronois_jets,truth_jets,m_jvorospt,m_jvoroseta,m_jvorosphi,m_jvorosmass,m_jvoroswidth,m_jvorosmindr,m_jvorosisPU,m_jvoroscl0pt,m_tjvorospt,m_tjvoroseta,m_tjvorosphi,m_tjvorosmass,m_tjvoroswidth,m_tjvorosmindr) != EL::StatusCode::SUCCESS)
    Error(APP_NAME,"Error in FillJetVars");

  if(FillJetVars(in_jets,truth_jets,m_jnoarea0pt,m_jnoarea0eta,m_jnoarea0phi,m_jnoarea0mass,m_jnoarea0width,m_jnoarea0mindr,m_jnoarea0isPU,m_jnoarea0cl0pt,m_tjnoarea0pt,m_tjnoarea0eta,m_tjnoarea0phi,m_tjnoarea0mass,m_tjnoarea0width,m_tjnoarea0mindr,false) != EL::StatusCode::SUCCESS)
    Error(APP_NAME,"Error in FillJetVars");
  if(FillJetVars(in_jets,truth_jets,m_j0pt,m_j0eta,m_j0phi,m_j0mass,m_j0width,m_j0mindr,m_j0isPU,m_j0cl0pt,m_tj0pt,m_tj0eta,m_tj0phi,m_tj0mass,m_tj0width,m_tj0mindr,true) != EL::StatusCode::SUCCESS)
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

EL::StatusCode WriteTree :: FillTJetVars(const DataVector<xAOD::Jet_v1>* tjets,
                                         std::vector<float>& tjetpt,
                                         std::vector<float>& tjeteta,
                                         std::vector<float>& tjetphi,
                                         std::vector<float>& tjetmass,
                                         std::vector<float>& tjetwidth,
                                         std::vector<float>& tjetmindr){
  tjetpt.clear();
  tjeteta.clear();
  tjetphi.clear();
  tjetmass.clear();
  tjetwidth.clear();
  tjetmindr.clear();

  static SG::AuxElement::ConstAccessor< float > minDR("minDR");
  for(const auto jet: *tjets){
    tjetpt.push_back(jet->pt()/1000.);
    tjeteta.push_back(jet->eta());
    tjetphi.push_back(jet->phi());
    tjetmass.push_back(jet->m()/1000.);
    tjetwidth.push_back(-99);
    tjetmindr.push_back(minDR(*jet));
  }
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
                                        std::vector<bool>& jetisPU,
                                        std::vector<float>& cl0pt,
                                        std::vector<float>& tjetpt,
                                        std::vector<float>& tjeteta,
                                        std::vector<float>& tjetphi,
                                        std::vector<float>& tjetmass,
                                        std::vector<float>& tjetwidth,
                                        std::vector<float>& tjetmindr,
                                        bool doareasub){

  jetpt.clear();
  jeteta.clear();
  jetphi.clear();
  jetmass.clear();
  jetwidth.clear();
  jetmindr.clear();
  jetisPU.clear();
  cl0pt.clear();
  tjetpt.clear();
  tjeteta.clear();
  tjetphi.clear();
  tjetmass.clear();
  tjetwidth.clear();
  tjetmindr.clear();

  DataVector<xAOD::Jet_v1> sorted_truth_jets = HF::sort_container_pt(tjets);
  static SG::AuxElement::ConstAccessor< int > truth_match_i("truth_match_i");
  static SG::AuxElement::ConstAccessor< float > minDR("minDR");
  static SG::AuxElement::ConstAccessor< bool > isPU("isPU");
  static SG::AuxElement::ConstAccessor< float > rho("rho");

  for(const auto jet: *jets){
    if(doareasub){
       float jpt = jet->pt();
       float activeArea(-99.0);
       jet->getAttribute("ActiveArea", activeArea);
       //float m_rho = rho(*jet);
       m_rho = rho(*jet);
       jpt = jpt-activeArea*m_rho;
       jetpt.push_back(jpt/1000.);
    }
    else{
      jetpt.push_back(jet->pt()/1000.);
    }
    jeteta.push_back(jet->eta());
    jetphi.push_back(jet->phi());
    jetmass.push_back(jet->m()/1000.);
    jetwidth.push_back(-99);
    jetmindr.push_back(minDR(*jet));
    jetisPU.push_back(isPU(*jet));

    std::vector<float> jetclpt;
    for(auto constit: jet->getConstituents()){
      jetclpt.push_back(constit->pt());
    }
    if(jetclpt.size()>0){
      auto biggest = std::max_element(std::begin(jetclpt),std::end(jetclpt));
      cl0pt.push_back(*biggest/1000.);
    }
    else cl0pt.push_back(-99);
  

    int truthmatch = truth_match_i(*jet);
    if(truthmatch > -1){
      auto tjet = sorted_truth_jets.at(truthmatch);
      tjetpt.push_back(tjet->pt()/1000.);
      tjeteta.push_back(tjet->eta());
      tjetphi.push_back(tjet->phi());
      tjetmass.push_back(tjet->m()/1000.);
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
