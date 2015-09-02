#include <vector>
#include <map>
#include <math.h>

// event loop
#include <EventLoop/Job.h>
#include <EventLoop/Worker.h>
#include "xAODEventInfo/EventInfo.h"

#include <MyAnalysis/VoronoiWeights.h>

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
#include "xAODCaloEvent/CaloClusterChangeSignalState.h"
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
ClassImp(VoronoiWeights)

VoronoiWeights :: VoronoiWeights () {}

EL::StatusCode VoronoiWeights :: setupJob (EL::Job& job)
{
  // let's initialize the algorithm to use the xAODRootAccess package
  job.useXAOD ();
  xAOD::Init("VoronoiWeights").ignore(); // call before opening first file

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode VoronoiWeights :: histInitialize () {
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode VoronoiWeights :: fileExecute () {return EL::StatusCode::SUCCESS;}

EL::StatusCode VoronoiWeights :: changeInput (bool /*firstFile*/) {return EL::StatusCode::SUCCESS;}

EL::StatusCode VoronoiWeights :: initialize ()
{
  m_event = wk()->xaodEvent();
  m_store = wk()->xaodStore();

  // as a check, let's see the number of events in our xAOD
  Info("initialize()", "Number of events = %lli", m_event->getEntries() ); // print long long int

  return EL::StatusCode::SUCCESS;
}

//Have to define custom comparator for PseudoJets in order to have a map from PJs to anything
//Comparison is fuzzy to account for rounding errors
struct VoronoiWeights :: PJcomp {
  bool operator() (const std::pair<fastjet::PseudoJet, float>& lhsp, const std::pair<fastjet::PseudoJet, float>& rhsp)
  {
    fastjet::PseudoJet lhs = lhsp.first;
    fastjet::PseudoJet rhs = rhsp.first;
    return lhs.pt()>rhs.pt();
    //The comparator must be a strict weak ordering. 
  }
};

EL::StatusCode VoronoiWeights :: execute ()
{
  const char* APP_NAME = "VoronoiWeights::execute()";

  const xAOD::EventInfo*                        eventInfo     (nullptr);
  const xAOD::CaloClusterContainer*             in_clusters   (nullptr);
  
  // start grabbing all the containers that we can
  RETURN_CHECK("VoronoiWeights::execute()", HF::retrieve(eventInfo,    m_eventInfo,        m_event, m_store, m_debug), "Could not get the EventInfo container.");
  if(!m_clust.empty())
    RETURN_CHECK("VoronoiWeights::execute()", HF::retrieve(in_clusters,     m_clust,       m_event, m_store, m_debug), "Could not get the clusters container.");

  clusters.clear();

  CaloClusterChangeSignalStateList stateHelperList;
  for(const auto clust: *in_clusters){
    //read in clusters as PseudoJets
    if(m_doLC) stateHelperList.add(clust,xAOD::CaloCluster::State(1)); //default is calibrated but we can make it explicit anyway
    else stateHelperList.add(clust,xAOD::CaloCluster::State(0));
    fastjet::PseudoJet test;
    test = fastjet::PseudoJet(clust->p4());
    if(clust->e() >= 0) clusters.push_back(test);
  }

  std::vector< std::pair< fastjet::PseudoJet,float > > ptvec; //vector of pairs of PJs and their corrected pTs
  if(MakeVoronoiClusters(ptvec) != EL::StatusCode::SUCCESS) Error(APP_NAME,"Error in MakeVoronoiClusters");
  std::sort(ptvec.begin(), ptvec.end(), PJcomp());

  int i=0;
  static SG::AuxElement::Decorator< float > correctedPt("correctedPt");
  for(const auto clust: HF::sort_container_pt(in_clusters)){
    correctedPt(*clust) = 0;
    if(m_debug){
      std::cout << "CDV Pt: " << clust->pt() << "; E: " << clust->e()<< std::endl;
      std::cout << "PT Vec Pt: " << ptvec[i].first.pt() << "; E: " << ptvec[i].first.e()<< std::endl;
    }

    //There should be the same number of positive E Clusters in the CDV as clusters in the ptvec
    bool endCDV = clust->e()<=0;
    bool endvec = i==ptvec.size();
    if(endCDV && endvec) continue;
    else if(endCDV || endvec){
      Error(APP_NAME,"Clusters don't have same number of elements.");
      return EL::StatusCode::FAILURE;
    }
    else{
      correctedPt(*clust) = ptvec[i].second;
      i++;
    }

    //And the clusters should match
    float CDVpt = clust->pt();
    float PJpt = ptvec[i-1].first.pt();
    if (fabs(CDVpt-PJpt) > 0.1){
      if(m_debug) std::cout << fabs(CDVpt-PJpt) << std::endl;
      Error(APP_NAME,"Clusters don't match.");
      return EL::StatusCode::FAILURE;
    }
  }

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode VoronoiWeights :: postExecute () {return EL::StatusCode::SUCCESS;}

EL::StatusCode VoronoiWeights :: finalize ()
{
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode VoronoiWeights :: histFinalize () {
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode VoronoiWeights::MakeVoronoiClusters(std::vector< std::pair< fastjet::PseudoJet,float > >& correctedptvec){
  std::vector<fastjet::PseudoJet> inputConst = clusters;
  fastjet::Selector jselector = fastjet::SelectorAbsRapRange(0.0,2.1);
  fastjet::JetAlgorithm algo = fastjet::kt_algorithm;
  float jetR = 0.4;
  fastjet::JetDefinition jetDef(algo, jetR,fastjet::E_scheme, fastjet::Best);
  fastjet::AreaDefinition area_def(fastjet::voronoi_area, fastjet::VoronoiAreaSpec(0.9));

  fastjet::ClusterSequenceArea clustSeq(inputConst, jetDef, area_def);
  fastjet::JetMedianBackgroundEstimator bge(jselector,jetDef,area_def);

  bge.set_particles(inputConst);
  std::vector<fastjet::PseudoJet> inclusiveJets = sorted_by_pt(clustSeq.inclusive_jets(0));

  std::map<int,float> result;
  for(unsigned int iJet = 0 ; iJet < inclusiveJets.size() ; iJet++){
    fastjet::PseudoJet jet = inclusiveJets[iJet];
    std::vector<fastjet::PseudoJet> constituents = jet.constituents();
    for(auto cons : constituents){
      float pt = cons.pt();
      float area = cons.area();
      float rho = bge.rho();
      float correctedPt = pt-rho*area;
      //std::cout << "Area: " << area << "; Rho: " << bge.rho() << "; pt: " << constituents[iCons].pt() << "; corrected: " << correctedPt << std::endl;
      //std::cout << "Pt: " << cons.pt() << "; Eta: " << cons.eta() <<"; Phi: " << cons.phi() << std::endl;
      //fastjet::PseudoJet constituentP;
      /*if(correctedPt<0.) continue;
      constituentP.reset_PtYPhiM(correctedPt, constituents[iCons].rap(), constituents[iCons].phi(), constituents[iCons].m());
      clusters_voronoi.push_back(constituentP);*/
      //correctedptmap[cons] = correctedPt;
      std::pair <fastjet::PseudoJet,float> pjcptpair (cons,correctedPt);
      correctedptvec.push_back(pjcptpair);
    } // end loop over cons
  } // end loop over jets
  //std::cout << "Size: " << correctedptmap.size() << std::endl;

return EL::StatusCode::SUCCESS;
}

