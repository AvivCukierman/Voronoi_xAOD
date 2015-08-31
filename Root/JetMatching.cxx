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
  const xAOD::CaloClusterContainer*             in_clusters   (nullptr);
  const xAOD::JetContainer*                     in_jets       (nullptr);
  const xAOD::JetContainer*                     truth_jets       (nullptr);
  const xAOD::JetContainer*                     voronoi_jets       (nullptr);

  // start grabbing all the containers that we can
  RETURN_CHECK("VoronoiWeights::execute()", HF::retrieve(eventInfo,    m_eventInfo,        m_event, m_store, m_debug), "Could not get the EventInfo container.");
  if(!m_jets.empty()) RETURN_CHECK("VoronoiWeights::execute()", HF::retrieve(in_jets,     m_jets,       m_event, m_store, m_debug), "Could not get the jets container.");
  if(!m_truth_jets.empty()) RETURN_CHECK("VoronoiWeights::execute()", HF::retrieve(truth_jets,    m_truth_jets,       m_event, m_store, m_debug), "Could not get the truth jets container.");
  if(!m_voronoi_jets.empty()) RETURN_CHECK("VoronoiWeights::execute()", HF::retrieve(voronoi_jets,    m_voronoi_jets,       m_event, m_store, m_debug), "Could not get the voronoi jets container.");

  FindTruthMatch(HF::sort_container_pt(voronoi_jets), HF::sort_container_pt(truth_jets));

  SetMinDR(HF::sort_container_pt(truth_jets));

  //To check truth matches:
  //DataVector<xAOD::Jet_v1> sorted_truth_jets = HF::sort_container_pt(truth_jets);
  /*static SG::AuxElement::ConstAccessor< int > truth_match_i("truth_match_i");
  std::cout << "Voronoi" << std::endl;
  for(auto jet: *voronoi_jets){
    std::cout << jet->pt() << std:: endl;
    int truthmatch = truth_match_i(*jet);
    std::cout << truthmatch  << std::endl;
    std::cout << sorted_truth_jets.at(truthmatch)->pt() << std::endl;
  }*/

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
  float deltaphi = fmod(fabs(jet1->phi() - jet2->phi()),2*3.14159265);
  float deltaeta = jet1->eta() - jet2->eta();
  return sqrt(pow(deltaphi,2) + pow(deltaeta,2));
}

EL::StatusCode JetMatching :: FindTruthMatch(DataVector<xAOD::Jet_v1> jets, DataVector<xAOD::Jet_v1> tjets){
  std::vector<int> matched;
  float MAXJETTRUTHMATCHDR = 0.3;
  static SG::AuxElement::Decorator< int > truth_match_i("truth_match_i");
  for(int iJ=0; iJ<jets.size(); iJ++) {
    DataVector<xAOD::Jet_v1>::ElementProxy jet = jets.at(iJ);
    //float mindR= 999.99;
    float maxPt=-999.99;
    //int minDRindex =-1;
    int maxPtIndex =-1;
    for(int iTrueJ=0; iTrueJ<tjets.size(); iTrueJ++){
      if(std::find(matched.begin(), matched.end(), iTrueJ) != matched.end())
        continue; //if true jet has already been matched, skip it -> bijective matching
      DataVector<xAOD::Jet_v1>::ElementProxy tjet = tjets.at(iTrueJ);

      float dR = deltaR(jet,tjet);

      //if(dR<mindR){ mindR = dR; minDRindex = iTrueJ;}
      if(dR < MAXJETTRUTHMATCHDR && maxPt < tjet->pt())
        { maxPt = tjet->pt(); maxPtIndex = iTrueJ;} //match to highest pT truth jet within MAXDR, not closest
    }//true jets

    if(maxPtIndex != -1){
      matched.push_back(maxPtIndex);
      /*if (!(jet(maxPtIndex, TruthJetType).Exists(JetType+"_match")))
        jet(maxPtIndex, TruthJetType).Add(JetType+"_match", thejet, true);*/ //add link from truth jet to jet as well eventually
    }
    truth_match_i(*jet) = maxPtIndex; //if no match truth_match_i == -1
  }//jet loop
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode JetMatching :: SetMinDR(DataVector<xAOD::Jet_v1> jets){
  static SG::AuxElement::Decorator< float > minDR("minDR");
  for(int iJ=0; iJ<jets.size(); iJ++) {
    DataVector<xAOD::Jet_v1>::ElementProxy jet1 = jets.at(iJ);
    float mindR= 999.99;
    for(int jJ=0; jJ<jets.size(); jJ++){
      if(iJ==jJ) continue;
      DataVector<xAOD::Jet_v1>::ElementProxy jet2 = jets.at(jJ);
      if(jet2->pt() < 5000) continue;

      float dR = deltaR(jet1,jet2);

      if(dR<mindR){ mindR = dR;}
    } //inner jet loop

    if(mindR < 999.99){
      minDR(*jet1) = mindR; 
    }
    else{
      return EL::StatusCode::FAILURE;
    }
  }//jet loop
  return EL::StatusCode::SUCCESS;
}