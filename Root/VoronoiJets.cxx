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
m_jetReclusteringTool{nullptr}
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

EL::StatusCode VoronoiJets :: changeInput (bool firstFile) {return EL::StatusCode::SUCCESS;}

EL::StatusCode VoronoiJets :: initialize ()
{
  m_event = wk()->xaodEvent();
  m_store = wk()->xaodStore();
  // as a check, let's see the number of events in our xAOD
  Info("initialize()", "Number of events = %lli", m_event->getEntries() ); // print long long int

  m_jetReclusteringTool = new JetReclusteringTool("JetReclusteringTool");
  RETURN_CHECK("VoronoiWeights::execute()",m_jetReclusteringTool->setProperty("InputJetContainer",  "VoronoiClustersCDV"),"Problem with jetReclusteringTool initialization");
  RETURN_CHECK("VoronoiWeights::execute()",m_jetReclusteringTool->setProperty("OutputJetContainer", "voronoi_jets"),"Problem with jetReclusteringTool initialization");
  RETURN_CHECK("VoronoiWeights::execute()",m_jetReclusteringTool->setProperty("ReclusterRadius",    0.4),"Problem with jetReclusteringTool initialization");
  RETURN_CHECK("VoronoiWeights::execute()",m_jetReclusteringTool->setProperty("ReclusterAlgorithm", fastjet::antikt_algorithm),"Problem with jetReclusteringTool initialization");
  RETURN_CHECK("VoronoiWeights::execute()",m_jetReclusteringTool->setProperty("InputJetPtMin",      0),"Problem with jetReclusteringTool initialization");
  RETURN_CHECK("VoronoiWeights::execute()",m_jetReclusteringTool->setProperty("RCJetPtMin",         0),"Problem with jetReclusteringTool initialization");
  RETURN_CHECK("VoronoiWeights::execute()",m_jetReclusteringTool->setProperty("RCJetPtFrac",        0),"Problem with jetReclusteringTool initialization");
  RETURN_CHECK("VoronoiWeights::execute()",m_jetReclusteringTool->initialize(),"Problem with jetReclusteringTool initialization");

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode VoronoiJets :: execute ()
{
  const char* APP_NAME = "VoronoiJets::execute()";
  const xAOD::EventInfo*                        eventInfo     (nullptr);
  const xAOD::CaloClusterContainer*             in_clusters   (nullptr);
  const xAOD::JetContainer*                     in_jets       (nullptr);

  // start grabbing all the containers that we can
  RETURN_CHECK("VoronoiWeights::execute()", HF::retrieve(eventInfo,    m_eventInfo,        m_event, m_store, m_debug), "Could not get the EventInfo container.");
  if(!m_clust.empty()) RETURN_CHECK("VoronoiWeights::execute()", HF::retrieve(in_clusters,     m_clust,       m_event, m_store, m_debug), "Could not get the clusters container.");
  if(!m_jets.empty()) RETURN_CHECK("VoronoiWeights::execute()", HF::retrieve(in_jets,     m_jets,       m_event, m_store, m_debug), "Could not get the jets container.");

  static SG::AuxElement::ConstAccessor< float > correctedPt("correctedPt");
  typedef xAOD::CaloClusterContainer ccc;
  ConstDataVector<ccc>* subset = new ConstDataVector<ccc>(SG::VIEW_ELEMENTS);
  std::pair< ccc*, xAOD::ShallowAuxContainer* > clustersSC = xAOD::shallowCopyContainer( *in_clusters );
  for(auto cluster: *(clustersSC.first)){
    float correctedPt_f = correctedPt(*cluster);
    if(correctedPt_f <= 0) continue;
    if(setClusterP4(cluster,xAOD::JetFourMom_t(correctedPt_f, cluster->eta(), cluster->phi(), cluster->m())) != EL::StatusCode::SUCCESS);
    subset->push_back(cluster);
  }

  m_store->record(clustersSC.first, "CorrectedCaloClusters");
  m_store->record(clustersSC.second, "CorrectedCaloClustersAux.");
  m_store->record(subset, "VoronoiClustersCDV");

/*  const xAOD::CaloClusterContainer* voronoi_clusters(nullptr);
  m_store->retrieve(voronoi_clusters, "VoronoiClustersCDV");*/

  m_jetReclusteringTool->execute();

/*  const xAOD::JetContainer*                     out_jets       (nullptr);
  RETURN_CHECK("VoronoiWeights::execute()", HF::retrieve(out_jets,     "voronoi_jets",       m_event, m_store, m_debug), "Could not get the voronoi jets container.");*/

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode VoronoiJets :: postExecute () {return EL::StatusCode::SUCCESS;}

EL::StatusCode VoronoiJets :: finalize ()
{
  if(m_jetReclusteringTool) delete m_jetReclusteringTool;
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode VoronoiJets :: histFinalize () {
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode VoronoiJets::setClusterP4(xAOD::CaloCluster* cluster, const xAOD::JetFourMom_t &p4){
  float pt = p4.pt();
  float e = p4.e();
  float eta = p4.eta();
  float phi = p4.phi();
  cluster->setE(pt*cosh(eta));
  cluster->setEta(eta);
  cluster->setPhi(phi);
  return EL::StatusCode::SUCCESS;
}
