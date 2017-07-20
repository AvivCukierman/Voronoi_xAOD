// ROOT generic
#include <TSystem.h>
#include <TTree.h>
#include <TFile.h>

#include <vector>

// event loop
#include <EventLoop/Job.h>
#include <EventLoop/Worker.h>
#include "xAODEventInfo/EventInfo.h"

#include <MyAnalysis/MyxAODAnalysis.h>
#include <VoronoiWeightTool/VoronoiWeightTool.h>
#include <JetReclustering/JetReclusteringTool.h>

#include <xAODJet/FastJetLinkBase.h>
#include "JetInterface/IPseudoJetGetter.h"
#include "fastjet/JetDefinition.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/Selector.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"
#include <fastjet/tools/Subtractor.hh>
#include "fastjet/PseudoJet.hh"
#include "xAODCaloEvent/CaloClusterContainer.h"

#include "xAODCore/ShallowAuxContainer.h"
#include "xAODJet/JetConstituentVector.h"

// EDM
#include "PATInterfaces/SystematicVariation.h"
#include "xAODBase/IParticleHelpers.h"
#include "xAODCaloEvent/CaloCluster.h"
#include "xAODCaloEvent/CaloClusterChangeSignalState.h"
#include "xAODTruth/TruthParticleContainer.h"
#include "xAODTruth/TruthEventContainer.h"
#include "xAODTruth/TruthEvent.h"

// utils
#include "CxxUtils/fpcompare.h"

// Infrastructure include(s):
#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TEvent.h"

// xAH includes
#include "xAODAnaHelpers/HelperFunctions.h"
//#include "xAODAnaHelpers/tools/ReturnCheck.h"

// this is needed to distribute the algorithm to the workers
ClassImp(MyxAODAnalysis)

namespace HF = HelperFunctions;

MyxAODAnalysis :: MyxAODAnalysis ():
   m_eventInfo(0),
   m_voronoiWeightTools{{nullptr, nullptr, nullptr}},
   m_jetReclusteringTools{{nullptr, nullptr, nullptr}}
{}

EL::StatusCode MyxAODAnalysis :: setupJob (EL::Job& job)
{
  // let's initialize the algorithm to use the xAODRootAccess package
  job.useXAOD ();
  xAOD::Init(); // call before opening first file

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode MyxAODAnalysis :: histInitialize () {
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode MyxAODAnalysis :: fileExecute () {return EL::StatusCode::SUCCESS;}

EL::StatusCode MyxAODAnalysis :: changeInput (bool /*firstFile*/) {return EL::StatusCode::SUCCESS;}

EL::StatusCode MyxAODAnalysis :: initialize ()
{
  m_event = wk()->xaodEvent();
  m_store = wk()->xaodStore();

  m_eventCounter = 0;
  // as a check, let's see the number of events in our xAOD
  Info("initialize()", "Number of events = %lli", m_event->getEntries() ); // print long long int

  std::string cluster_containers[4] = {"voro0clusters","voro1clusters","vorosclusters",m_inputContainer};
  for(int i=0; i<3; i++){
    m_voronoiWeightTools[i] = new VoronoiWeightTool(m_outputContainer+std::to_string(std::rand()));
    m_voronoiWeightTools[i]->setProperty("InputContainer", m_inputContainer);
    m_voronoiWeightTools[i]->setProperty("doLCWeights", m_doLC);
    switch(i){
      case 0:
        m_voronoiWeightTools[i]->setProperty("doSpread", false);
        m_voronoiWeightTools[i]->setProperty("nSigma", 0);
        break;
      case 1: 
        m_voronoiWeightTools[i]->setProperty("doSpread", false);
        m_voronoiWeightTools[i]->setProperty("nSigma", 1);
        break;
      case 2: 
        m_voronoiWeightTools[i]->setProperty("doSpread", true);
        m_voronoiWeightTools[i]->setProperty("nSigma", 0);
        break;
    }
    m_voronoiWeightTools[i]->setProperty("OutputContainer", cluster_containers[i]);
    m_voronoiWeightTools[i]->initialize();
  }

  std::string outputContainer;
  for(int i=0; i<4; i++){
    switch(i){
      case 0:
        outputContainer = "AntiKt4Voronoi0Jets";
        break;
      case 1: 
        outputContainer = "AntiKt4Voronoi1Jets";
        break;
      case 2: 
        outputContainer = "AntiKt4VoronoiSpreadJets";
        break;
      case 3:
        outputContainer = "AntiKt4NoAreaJets";
        break;
    }
    m_jetReclusteringTools[i] = new JetReclusteringTool(outputContainer+std::to_string(std::rand()));
    ANA_CHECK(m_jetReclusteringTools[i]->setProperty("InputJetContainer",  cluster_containers[i]));
    ANA_CHECK(m_jetReclusteringTools[i]->setProperty("OutputJetContainer", outputContainer));
    ANA_CHECK(m_jetReclusteringTools[i]->setProperty("ReclusterRadius",    0.4));
    //ANA_CHECK(m_jetReclusteringTools[i]->setProperty("ReclusterAlgorithm", fastjet::antikt_algorithm));
    ANA_CHECK(m_jetReclusteringTools[i]->setProperty("ReclusterAlgorithm", "AntiKt"));
    ANA_CHECK(m_jetReclusteringTools[i]->setProperty("InputJetPtMin",      0));
    ANA_CHECK(m_jetReclusteringTools[i]->setProperty("RCJetPtMin",         0.1));
    ANA_CHECK(m_jetReclusteringTools[i]->setProperty("RCJetPtFrac",        0));
    ANA_CHECK(m_jetReclusteringTools[i]->setProperty("DoArea",        true));
    ANA_CHECK(m_jetReclusteringTools[i]->initialize());
  }


  return EL::StatusCode::SUCCESS;
}

EL::StatusCode MyxAODAnalysis :: execute ()
{
  const char* APP_NAME = "MyxAODAnalysis::execute()";

  const xAOD::CaloClusterContainer*             in_clusters   (nullptr);
  const xAOD::EventInfo*                        eventInfo     (nullptr);
  
  // start grabbing all the containers that we can
  std::string m_clust = "CaloCalTopoClusters";
  if(!m_clust.empty())
    ANA_CHECK(HF::retrieve(in_clusters,     m_clust,       m_event, m_store, msg()));
  std::string m_eventInfo      = "EventInfo";
  ANA_CHECK(HF::retrieve(eventInfo,    m_eventInfo,        m_event, m_store, msg()));
   //int event_number = eventInfo->eventNumber(); //added
   //std::cout << event_number << std::endl;

  /*std::cout << "CalCalTopoClusters" << std::endl;
  for(const auto clust: *in_clusters){
    if(clust->e() >= 0) std::cout << clust->e() << std::endl;
  }*/

  for(int i=0; i<3; i++) m_voronoiWeightTools[i]->execute();
 

  // Have to access the clusters in the state they were set in order to get voronoi-subtracted clusters
  CaloClusterChangeSignalStateList stateHelperList;
  std::string cluster_containers[3] = {"voro0clusters","voro1clusters","vorosclusters"};
  const xAOD::CaloClusterContainer*             voronoi_clusters   (nullptr);
  if(!cluster_containers[0].empty())
    ANA_CHECK(HF::retrieve(voronoi_clusters,     cluster_containers[0],       m_event, m_store, msg()));
  const xAOD::CaloClusterContainer*             voronoi1_clusters   (nullptr);
  if(!cluster_containers[1].empty())
    ANA_CHECK(HF::retrieve(voronoi1_clusters,     cluster_containers[1],       m_event, m_store, msg()));
  const xAOD::CaloClusterContainer*             voronoispread_clusters   (nullptr);
  if(!cluster_containers[2].empty())
    ANA_CHECK(HF::retrieve(voronoispread_clusters,     cluster_containers[2],       m_event, m_store, msg()));
  for(const auto clust: *in_clusters){
    if(m_doLC) stateHelperList.add(clust,xAOD::CaloCluster::State(1)); //default is calibrated but we can make it explicit anyway
    else stateHelperList.add(clust,xAOD::CaloCluster::State(0));
    //if(clust->e() >= 100000) std::cout << "EM clust: " << clust->e() << ";" << clust->eta()<< ";" << clust->phi()<< std::endl;
  }
  for(const auto clust: *voronoi_clusters){
    if(m_doLC) stateHelperList.add(clust,xAOD::CaloCluster::State(1)); //default is calibrated but we can make it explicit anyway
    else stateHelperList.add(clust,xAOD::CaloCluster::State(0));
    //if(clust->e() >= 100000) std::cout << "V0 clust: " << clust->e() << ";" << clust->eta()<< ";" << clust->phi()<< std::endl;
  }
  for(const auto clust: *voronoi1_clusters){
    if(m_doLC) stateHelperList.add(clust,xAOD::CaloCluster::State(1)); //default is calibrated but we can make it explicit anyway
    else stateHelperList.add(clust,xAOD::CaloCluster::State(0));
    //if(clust->e() >= 100000) std::cout << "V1 clust: " << clust->e() << ";" << clust->eta()<< ";" << clust->phi()<< std::endl;
  }
  for(const auto clust: *voronoispread_clusters){
    if(m_doLC) stateHelperList.add(clust,xAOD::CaloCluster::State(1)); //default is calibrated but we can make it explicit anyway
    else stateHelperList.add(clust,xAOD::CaloCluster::State(0));
    //if(clust->e() >= 100000) std::cout << "VS clust: " << clust->e() << ";" << clust->eta()<< ";" << clust->phi()<< std::endl;
  }
  

  //jets
  for(int i=0; i<4; i++) m_jetReclusteringTools[i]->execute();

  float event_rho,event_sigma;
  FindRho(in_clusters,event_rho,event_sigma);
  /*std::cout << "event_rho: " << event_rho << std::endl;
  m_store->record(event_rho,std::string("rho")); // This works but I have no idea how to retrieve it.
  float* new_rho;
  std::cout << *new_rho << std::endl;
  HF::retrieve(new_rho, std::string("rho"), m_event, m_store, msg());
  std::cout << "new_rho: " << *new_rho << std::endl;*/

  const xAOD::JetContainer*             test_jets   (nullptr);
  std::string outputContainer = "AntiKt4NoAreaJets";
//  std::cout << outputContainer << std::endl;
    ANA_CHECK(HF::retrieve(test_jets,     outputContainer,       m_event, m_store, msg()));
  SetRho(HF::sort_container_pt(test_jets),event_rho);

  /*static SG::AuxElement::ConstAccessor< float > rho("rho");
  for(auto jet: *test_jets){
    float activeArea(-99.0);
    jet->getAttribute("ActiveArea", activeArea);
    if(jet->pt() > 5000) std::cout << "Custom: " << jet->pt()-rho(*jet)*activeArea << std::endl;
  }

  const xAOD::JetContainer * oldjets (nullptr);
  std::string m_oldjets = "AntiKt4EMTopoJets";
  ANA_CHECK(HF::retrieve(oldjets,     m_oldjets,       m_event, m_store, msg()));
  for(auto jet: *oldjets){
    std::cout << "Old: " << jet->pt() << std::endl;
  }*/
  
  /*const xAOD::JetContainer*             voronoi_jets   (nullptr);
  outputContainer = "AntiKt4Voronoi0Jets";
    ANA_CHECK(HF::retrieve(voronoi_jets,     outputContainer,       m_event, m_store, msg()));
  for(auto jet: *voronoi_jets){
    std::cout << "V0 Jet: " << jet->pt() << "\t" << jet->eta() << "\t" << jet->phi() << "\t" << jet->m() << std::endl;
    for(auto constit: jet->getConstituents()){
      std::cout<< "Constit pT: " << constit->pt() << std::endl; 
    }
  }*/
  /*const xAOD::JetContainer*             voronoispread_jets   (nullptr);
  outputContainer = "AntiKt4VoronoiSpreadJets";
  ANA_CHECK(HF::retrieve(voronoispread_jets,     outputContainer,       m_event, m_store, msg()));
  for(auto jet: *voronoispread_jets){
    std::cout << "VS Jet: " << jet->pt() << "\t" << jet->eta() << "\t" << jet->phi() << "\t" << jet->m() << std::endl;
  }*/

  /*  if(event_number == 455939){
    }
  */
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode MyxAODAnalysis :: postExecute () {return EL::StatusCode::SUCCESS;}

EL::StatusCode MyxAODAnalysis :: finalize ()
{
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode MyxAODAnalysis :: histFinalize () {
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode MyxAODAnalysis::SetRho(DataVector<xAOD::Jet_v1> jets,float event_rho){
  static SG::AuxElement::Decorator< float > rho("rho");
  for(int iJ=0; iJ<jets.size(); iJ++) {
    DataVector<xAOD::Jet_v1>::ElementProxy jet = jets.at(iJ);
    rho(*jet) = event_rho;
  }
}

EL::StatusCode MyxAODAnalysis::FindRho(const xAOD::CaloClusterContainer* in_clusters,float& rho, float& sigma){
  //CaloClusterChangeSignalStateList stateHelperList;
  std::vector<fastjet::PseudoJet> inputConst;
  for(const auto clust: *in_clusters){
    //if(m_doLC) stateHelperList.add(clust,xAOD::CaloCluster::State(1)); //default is calibrated but we can make it explicit anyway
    //else stateHelperList.add(clust,xAOD::CaloCluster::State(0));
    fastjet::PseudoJet test;
    test = fastjet::PseudoJet(clust->p4());
    if(clust->e() >= 0) inputConst.push_back(test);
  }

  fastjet::Selector jselector = fastjet::SelectorAbsRapRange(0.0,2.1);
  fastjet::JetAlgorithm algo = fastjet::kt_algorithm;
  float jetR = 0.4;
  fastjet::JetDefinition jetDef(algo, jetR,fastjet::E_scheme, fastjet::Best);
  fastjet::AreaDefinition area_def(fastjet::voronoi_area, fastjet::VoronoiAreaSpec(0.9));

  fastjet::JetMedianBackgroundEstimator bge(jselector,jetDef,area_def);
  bge.set_particles(inputConst);

  rho = bge.rho();
  sigma = bge.sigma();
}
