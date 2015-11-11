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
#include "xAODAnaHelpers/tools/ReturnCheck.h"

// this is needed to distribute the algorithm to the workers
ClassImp(MyxAODAnalysis)

namespace HF = HelperFunctions;

MyxAODAnalysis :: MyxAODAnalysis ():
   m_eventInfo(0),
   m_VoronoiTool(0)
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

  m_VoronoiTool = new VoronoiWeightTool("VoronoiTest");
  m_VoronoiTool->setProperty("InputContainer", "CaloCalTopoClusters");
  m_VoronoiTool->setProperty("OutputContainer", "VoronoiClusters");
  m_VoronoiTool->initialize();

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode MyxAODAnalysis :: execute ()
{
  const char* APP_NAME = "MyxAODAnalysis::execute()";

  const xAOD::CaloClusterContainer*             in_clusters   (nullptr);
  
  // start grabbing all the containers that we can
  std::string m_clust = "CaloCalTopoClusters";
  if(!m_clust.empty())
    RETURN_CHECK("VoronoiWeights::execute()", HF::retrieve(in_clusters,     m_clust,       m_event, m_store, m_debug), "Could not get the clusters container.");

  /*for(const auto clust: *in_clusters){
    if(clust->e() >= 0) std::cout << clust->e() << std::endl;
  }*/

  if(m_VoronoiTool){
    m_VoronoiTool->execute();
  }

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
