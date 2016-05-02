#include <vector>
#include <map>
#include <math.h>

// event loop
#include <EventLoop/Job.h>
#include <EventLoop/Worker.h>
#include "xAODEventInfo/EventInfo.h"

#include <MyAnalysis/VoronoiJets.h>

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
#include "xAODCaloEvent/CaloClusterChangeSignalState.h"
#include "xAODBase/IParticleHelpers.h"

// Infrastructure include(s):
#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TEvent.h"

// xAH includes
#include "xAODAnaHelpers/HelperFunctions.h"
#include "xAODAnaHelpers/tools/ReturnCheck.h"

//Reclustering
#include <xAODJetReclustering/JetReclusteringTool.h>

namespace HF = HelperFunctions;

// this is needed to distribute the algorithm to the workers
ClassImp(VoronoiJets)

VoronoiJets :: VoronoiJets () :
//m_jetReclusteringTool{nullptr}
m_jetReclusteringTools{{nullptr, nullptr, nullptr}}
{}

EL::StatusCode VoronoiJets :: setupJob (EL::Job& job)
{
  // let's initialize the algorithm to use the xAODRootAccess package
  job.useXAOD ();
  xAOD::Init("VoronoiWeights").ignore(); // call before opening first file

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode VoronoiJets :: histInitialize () {
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode VoronoiJets :: fileExecute () {return EL::StatusCode::SUCCESS;}

EL::StatusCode VoronoiJets :: changeInput (bool /*firstFile*/) {return EL::StatusCode::SUCCESS;}

EL::StatusCode VoronoiJets :: initialize ()
{
  m_debug = false;
  m_event = wk()->xaodEvent();
  m_store = wk()->xaodStore();
  // as a check, let's see the number of events in our xAOD
  Info("initialize()", "Number of events = %lli", m_event->getEntries() ); // print long long int

  for(int i=0; i<3; i++){
    //m_jetReclusteringTool = new JetReclusteringTool("JetReclusteringTool");
    std::string outputContainer,inputContainer;
    switch(i){
      case 0:
        outputContainer = "AntiKt4Voronoi0Jets";
        inputContainer = "Voronoi0ClustersCDV";
        break;
      case 1: 
        outputContainer = "AntiKt4Voronoi1Jets";
        inputContainer = "Voronoi1ClustersCDV";
        break;
      case 2: 
        outputContainer = "AntiKt4VoronoiSpreadJets";
        inputContainer = "VoronoiSpreadClustersCDV";
        break;
    }
    m_jetReclusteringTools[i] = new JetReclusteringTool(outputContainer+std::to_string(std::rand()));
    RETURN_CHECK("VoronoiWeights::execute()",m_jetReclusteringTools[i]->setProperty("InputJetContainer",  inputContainer),"Problem with jetReclusteringTool initialization");
    RETURN_CHECK("VoronoiWeights::execute()",m_jetReclusteringTools[i]->setProperty("OutputJetContainer", outputContainer),"Problem with jetReclusteringTool initialization");
    RETURN_CHECK("VoronoiWeights::execute()",m_jetReclusteringTools[i]->setProperty("ReclusterRadius",    0.4),"Problem with jetReclusteringTool initialization");
    RETURN_CHECK("VoronoiWeights::execute()",m_jetReclusteringTools[i]->setProperty("ReclusterAlgorithm", fastjet::antikt_algorithm),"Problem with jetReclusteringTool initialization");
    RETURN_CHECK("VoronoiWeights::execute()",m_jetReclusteringTools[i]->setProperty("InputJetPtMin",      0),"Problem with jetReclusteringTool initialization");
    RETURN_CHECK("VoronoiWeights::execute()",m_jetReclusteringTools[i]->setProperty("RCJetPtMin",         7),"Problem with jetReclusteringTool initialization");
    RETURN_CHECK("VoronoiWeights::execute()",m_jetReclusteringTools[i]->setProperty("RCJetPtFrac",        0),"Problem with jetReclusteringTool initialization");
    RETURN_CHECK("VoronoiWeights::execute()",m_jetReclusteringTools[i]->initialize(),"Problem with jetReclusteringTool initialization");
  }

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode VoronoiJets :: execute ()
{
  const char* APP_NAME = "VoronoiJets::execute()";
  const xAOD::EventInfo*                        eventInfo     (nullptr);
  const xAOD::CaloClusterContainer*             in_clusters   (nullptr);

  // start grabbing all the containers that we can
  RETURN_CHECK("VoronoiWeights::execute()", HF::retrieve(eventInfo,    m_eventInfo,        m_event, m_store, m_debug), "Could not get the EventInfo container.");
  if(!m_clust.empty()) RETURN_CHECK("VoronoiWeights::execute()", HF::retrieve(in_clusters,     m_clust,       m_event, m_store, m_debug), "Could not get the clusters container.");
  //int event_number = eventInfo->eventNumber(); //added

  static SG::AuxElement::ConstAccessor< float > voro0Pt("voro0Pt");
  static SG::AuxElement::ConstAccessor< float > voro1Pt("voro1Pt");
  static SG::AuxElement::ConstAccessor< float > spreadPt("spreadPt");
  typedef xAOD::CaloClusterContainer ccc;

  /*CaloClusterChangeSignalStateList stateHelperList;
  for(auto cluster: *(in_clusters)){
    if(m_doLC) stateHelperList.add(cluster,xAOD::CaloCluster::State(1));
    else stateHelperList.add(cluster,xAOD::CaloCluster::State(0));
  }*/

  for(int i=0; i<3; i++){
    std::pair< ccc*, xAOD::ShallowAuxContainer* > clustersSC = xAOD::shallowCopyContainer( *in_clusters );
    ConstDataVector<ccc>* subset = new ConstDataVector<ccc>(SG::VIEW_ELEMENTS);
    for(auto cluster: *(clustersSC.first)){
      float correctedPt_f;
      if(i==0) correctedPt_f = voro0Pt(*cluster);
      if(i==1) correctedPt_f = voro1Pt(*cluster);
      if(i==2) correctedPt_f = spreadPt(*cluster);
      if(correctedPt_f <= 0) continue;
      if(setClusterP4(cluster,xAOD::JetFourMom_t(correctedPt_f, cluster->eta(), cluster->phi(), cluster->m())) != EL::StatusCode::SUCCESS)
        Error(APP_NAME,"Error in setClusterP4");
      subset->push_back(cluster);
    }

    switch(i){
      case 0:
        m_store->record(clustersSC.first, "CorrectedCaloClusters0");
        m_store->record(clustersSC.second, "CorrectedCaloClusters0Aux.");
        m_store->record(subset, "Voronoi0ClustersCDV");
        break;
      case 1:
        m_store->record(clustersSC.first, "CorrectedCaloClusters1");
        m_store->record(clustersSC.second, "CorrectedCaloClusters1Aux.");
        m_store->record(subset, "Voronoi1ClustersCDV");
        break;
      case 2:
        m_store->record(clustersSC.first, "CorrectedCaloClusters2");
        m_store->record(clustersSC.second, "CorrectedCaloClusters2Aux.");
        m_store->record(subset, "VoronoiSpreadClustersCDV");
        break;
    }
  }
  if(m_debug){
    const xAOD::CaloClusterContainer* voronoi_clusters(nullptr);
    m_store->retrieve(voronoi_clusters, "VoronoiSpreadClustersCDV");
    std::cout << "Clusters" << std::endl;
    for(auto cluster:*voronoi_clusters){
      std::cout << cluster->pt() << ";" << cluster->eta() << ";" << cluster->phi() << ";" << cluster->m() << std::endl;
    }
  }

  for(int i=0; i<3; i++) m_jetReclusteringTools[i]->execute();
  
  if(m_debug){
    const xAOD::JetContainer*                     out_jets       (nullptr);
    RETURN_CHECK("VoronoiWeights::execute()", HF::retrieve(out_jets,     "AntiKt4VoronoiSpreadJets",       m_event, m_store, m_debug), "Could not get the voronoi jets container.");
    for(auto jet: *out_jets){
      std::cout << "Jet: " << jet->pt() << ";" << jet->eta() << ";" << jet->phi() << ";" << jet->m() << std::endl;
      for(auto constit: jet->getConstituents()){
        std::cout<< "Constit pT: " << constit->pt() << std::endl; 
      }
    }
  }
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode VoronoiJets :: postExecute () {return EL::StatusCode::SUCCESS;}

EL::StatusCode VoronoiJets :: finalize ()
{
  for(int i=0; i<3; i++){
    if(m_jetReclusteringTools[i]) delete m_jetReclusteringTools[i];
  }
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode VoronoiJets :: histFinalize () {
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode VoronoiJets::setClusterP4(xAOD::CaloCluster* cluster, const xAOD::JetFourMom_t &p4){
  float pt = p4.pt();
  //float e = p4.e();
  float eta = p4.eta();
  float phi = p4.phi();
  cluster->setE(pt*cosh(eta));
  cluster->setEta(eta);
  cluster->setPhi(phi);
  return EL::StatusCode::SUCCESS;
}
