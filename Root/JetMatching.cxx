#include <vector>
#include <map>
#include <math.h>

// event loop
#include <EventLoop/Job.h>
#include <EventLoop/Worker.h>
#include "xAODEventInfo/EventInfo.h"

#include <MyAnalysis/JetMatching.h>

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

// this is needed to distribute the algorithm to the workers
ClassImp(JetMatching)

JetMatching :: JetMatching () {}

EL::StatusCode JetMatching :: setupJob (EL::Job& job)
{
  // let's initialize the algorithm to use the xAODRootAccess package
  job.useXAOD ();
  xAOD::Init("VoronoiWeights").ignore(); // call before opening first file

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode JetMatching :: histInitialize () {
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode JetMatching :: fileExecute () {return EL::StatusCode::SUCCESS;}

EL::StatusCode JetMatching :: changeInput (bool /*firstFile*/) {return EL::StatusCode::SUCCESS;}

EL::StatusCode JetMatching :: initialize ()
{
  m_event = wk()->xaodEvent();
  m_store = wk()->xaodStore();

  // as a check, let's see the number of events in our xAOD
  Info("initialize()", "Number of events = %lli", m_event->getEntries() ); // print long long int

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode JetMatching :: execute ()
{
  const char* APP_NAME = "JetMatching::execute()";
  const xAOD::EventInfo*                        eventInfo     (nullptr);
  const xAOD::JetContainer*                     in_jets       (nullptr);
  const xAOD::JetContainer*                     truth_jets       (nullptr);
  const xAOD::JetContainer*                     voronoi0_jets       (nullptr);
  const xAOD::JetContainer*                     voronoi1_jets       (nullptr);
  const xAOD::JetContainer*                     voronoispread_jets       (nullptr);
  
  // start grabbing all the containers that we can
  RETURN_CHECK("VoronoiWeights::execute()", HF::retrieve(eventInfo,    m_eventInfo,        m_event, m_store, m_debug), "Could not get the EventInfo container.");
  //m_jets = "AntiKt4EMTopoJets";
  if(!m_jets.empty()) RETURN_CHECK("VoronoiWeights::execute()", HF::retrieve(in_jets,     m_jets,       m_event, m_store, m_debug), "Could not get the jets container.");
  if(!m_truth_jets.empty()) RETURN_CHECK("VoronoiWeights::execute()", HF::retrieve(truth_jets,    m_truth_jets,       m_event, m_store, m_debug), "Could not get the truth jets container.");
  if(!m_voronoi0_jets.empty()) RETURN_CHECK("VoronoiWeights::execute()", HF::retrieve(voronoi0_jets,    m_voronoi0_jets,       m_event, m_store, m_debug), "Could not get the voronoi jets container.");
  if(!m_voronoi1_jets.empty()) RETURN_CHECK("VoronoiWeights::execute()", HF::retrieve(voronoi1_jets,    m_voronoi1_jets,       m_event, m_store, m_debug), "Could not get the voronoi jets container.");
  if(!m_voronoispread_jets.empty()) RETURN_CHECK("VoronoiWeights::execute()", HF::retrieve(voronoispread_jets,    m_voronoispread_jets,       m_event, m_store, m_debug), "Could not get the voronoi jets container.");

  if(FindTruthMatch(HF::sort_container_pt(voronoi0_jets), HF::sort_container_pt(truth_jets)) != EL::StatusCode::SUCCESS)
    Error(APP_NAME,"Error in FindTruthMatch");
  if(FindTruthMatch(HF::sort_container_pt(voronoi1_jets), HF::sort_container_pt(truth_jets)) != EL::StatusCode::SUCCESS)
    Error(APP_NAME,"Error in FindTruthMatch");
  if(FindTruthMatch(HF::sort_container_pt(voronoispread_jets), HF::sort_container_pt(truth_jets)) != EL::StatusCode::SUCCESS)
    Error(APP_NAME,"Error in FindTruthMatch");
  if(FindTruthMatch(HF::sort_container_pt(in_jets), HF::sort_container_pt(truth_jets)) != EL::StatusCode::SUCCESS)
    Error(APP_NAME,"Error in FindTruthMatch");

  if(SetMinDR(HF::sort_container_pt(truth_jets)) != EL::StatusCode::SUCCESS)
    Error(APP_NAME,"Error in SetMinDR (only one jet in event?)");
  if(SetMinDR(HF::sort_container_pt(in_jets),2000) != EL::StatusCode::SUCCESS)
    Error(APP_NAME,"Error in SetMinDR (only one jet in event?)");
  if(SetMinDR(HF::sort_container_pt(voronoi0_jets),2000) != EL::StatusCode::SUCCESS)
    Error(APP_NAME,"Error in SetMinDR (only one jet in event?)");
  if(SetMinDR(HF::sort_container_pt(voronoi1_jets),2000) != EL::StatusCode::SUCCESS)
    Error(APP_NAME,"Error in SetMinDR (only one jet in event?)");
  if(SetMinDR(HF::sort_container_pt(voronoispread_jets),2000) != EL::StatusCode::SUCCESS)
    Error(APP_NAME,"Error in SetMinDR (only one jet in event?)");

  //To check truth matches:
  if(m_debug){
    int event_number = eventInfo->eventNumber(); //added
    if(event_number == 455939){
      DataVector<xAOD::Jet_v1> sorted_truth_jets = HF::sort_container_pt(truth_jets);
      static SG::AuxElement::ConstAccessor< int > truth_match_i("truth_match_i");
      for(auto jet: *voronoispread_jets){
        std::cout << "Reco jet: " << jet->pt() << "\t" << jet->eta() << "\t" << jet->phi() << "\t" << jet->m() << std::endl;
        int truthmatch = truth_match_i(*jet);
        if(truthmatch>-1){
          auto tjet =  sorted_truth_jets.at(truthmatch);
          std::cout << "True jet: " << tjet->pt() << "\t" << tjet->eta() << "\t" << tjet->phi() << "\t" << tjet->m() << std::endl;
        }
      }
    }

    char filename[50];
    sprintf(filename,"truth_jets_%i",event_number);
    std::ofstream fout(filename);
    for(auto jet: *truth_jets){
      fout << jet->pt() << "\t" << jet->eta() << "\t" << jet->phi() << "\t" << jet->m() << std::endl;
    }
  }
  //To check minDR:
  /*static SG::AuxElement::ConstAccessor< float > minDR("minDR");
  for(auto jet: *truth_jets){
    std::cout << minDR(*jet) << std::endl;
  }*/
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode JetMatching :: postExecute () {return EL::StatusCode::SUCCESS;}

EL::StatusCode JetMatching :: finalize ()
{
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode JetMatching :: histFinalize () {
  return EL::StatusCode::SUCCESS;
}

float deltaR(DataVector<xAOD::Jet_v1>::ElementProxy jet1, DataVector<xAOD::Jet_v1>::ElementProxy jet2){
  fastjet::PseudoJet pjet1 = fastjet::PseudoJet(jet1->p4());
  fastjet::PseudoJet pjet2 = fastjet::PseudoJet(jet2->p4());
  float deltaphi = pjet1.delta_phi_to(pjet2);
  float deltaeta = jet1->eta() - jet2->eta();
  float dR =  sqrt(pow(deltaphi,2) + pow(deltaeta,2));
  return dR;
}

EL::StatusCode JetMatching :: FindTruthMatch(DataVector<xAOD::Jet_v1> jets, DataVector<xAOD::Jet_v1> tjets){
  std::vector<int> matched;
  float MAXJETTRUTHMATCHDR = 0.3;
  float MINPUJETTRUTHMATCHDR = 0.6;
  static SG::AuxElement::Decorator< int > truth_match_i("truth_match_i");
  static SG::AuxElement::Decorator< bool > isPU("isPU");
  for(int iJ=0; iJ<jets.size(); iJ++) {
    DataVector<xAOD::Jet_v1>::ElementProxy jet = jets.at(iJ);
    float mindR= 999.99;
    float maxPt=-999.99;
    //int minDRindex =-1;
    int maxPtIndex =-1;
    for(int iTrueJ=0; iTrueJ<tjets.size(); iTrueJ++){
      DataVector<xAOD::Jet_v1>::ElementProxy tjet = tjets.at(iTrueJ);

      float dR = deltaR(jet,tjet);
      if(dR<mindR){ mindR = dR;}
      if(std::find(matched.begin(), matched.end(), iTrueJ) != matched.end())
        continue; //if true jet has already been matched, skip it -> truth jet is matched to at most one reco jet 
      if(dR < MAXJETTRUTHMATCHDR && maxPt < tjet->pt())
        { maxPt = tjet->pt(); maxPtIndex = iTrueJ;} //match to highest pT truth jet within MAXDR, not closest
    }//true jets

    if(maxPtIndex != -1){
      matched.push_back(maxPtIndex);
      /*if (!(jet(maxPtIndex, TruthJetType).Exists(JetType+"_match")))
        jet(maxPtIndex, TruthJetType).Add(JetType+"_match", thejet, true);*/ //add link from truth jet to jet as well eventually
    }
    truth_match_i(*jet) = maxPtIndex; //if no match truth_match_i == -1
    isPU(*jet) = mindR>MINPUJETTRUTHMATCHDR; //mindR to any truth jet
  }//jet loop
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode JetMatching :: SetMinDR(DataVector<xAOD::Jet_v1> jets,float threshold){
  static SG::AuxElement::Decorator< float > minDR("minDR");
  for(int iJ=0; iJ<jets.size(); iJ++) {
    DataVector<xAOD::Jet_v1>::ElementProxy jet1 = jets.at(iJ);
    float mindR= 999.99;
    for(int jJ=0; jJ<jets.size(); jJ++){
      if(iJ==jJ) continue;
      DataVector<xAOD::Jet_v1>::ElementProxy jet2 = jets.at(jJ);
      if(jet2->pt() < threshold) continue;

      float dR = deltaR(jet1,jet2);

      if(dR<mindR){ mindR = dR;}
    } //inner jet loop

    //if(mindR < 999.99){
    minDR(*jet1) = mindR; 
    //}
    //else{
    //  return EL::StatusCode::FAILURE;
    //}
  }//jet loop
  return EL::StatusCode::SUCCESS;
}
