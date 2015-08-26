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

namespace HF = HelperFunctions;

// this is needed to distribute the algorithm to the workers
ClassImp(VoronoiJets)

VoronoiJets :: VoronoiJets () {}

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

  clusters_voronoi.clear();
  for(const auto clust: *in_clusters){
    //std::cout << "Corrected: " << correctedPt(*clust) << std::endl;
    //std::cout << "Original: " << clust->pt() << std::endl;
    fastjet::PseudoJet pjclust;
    if(correctedPt(*clust)>0 && clust->e()>0){
      pjclust.reset_PtYPhiM(correctedPt(*clust), clust->rapidity(), clust->phi(), clust->m());
      clusters_voronoi.push_back(pjclust);
    }
  }
  MakeJetsWArea(fastjet::antikt_algorithm, 0.4, clusters_voronoi, jets_voronoi, false);

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode VoronoiJets :: postExecute () {return EL::StatusCode::SUCCESS;}

EL::StatusCode VoronoiJets :: finalize ()
{
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode VoronoiJets :: histFinalize () {
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode VoronoiJets::MakeJetsWArea(const fastjet::JetAlgorithm algo, const double jetR, const std::vector<fastjet::PseudoJet>& inputConst, std::vector<fastjet::PseudoJet>& out_jets,bool doareasub){

  fastjet::JetDefinition jetDef(algo, jetR,fastjet::E_scheme, fastjet::Best);
  fastjet::AreaDefinition active_area = fastjet::AreaDefinition(fastjet::active_area);
  fastjet::ClusterSequenceArea clustSeq(inputConst, jetDef, active_area);
  std::vector<fastjet::PseudoJet> inclusiveJets = sorted_by_pt(clustSeq.inclusive_jets(5.));
  std::vector<fastjet::PseudoJet> subtractedJets = inclusiveJets;

  if(doareasub){
    fastjet::Selector jselector = fastjet::SelectorAbsRapRange(0.0,2.1);
    fastjet::JetDefinition jetDefkt(fastjet::kt_algorithm,0.4);
    fastjet::AreaDefinition voronoi_area = fastjet::AreaDefinition(fastjet::voronoi_area,fastjet::VoronoiAreaSpec(0.9));
    fastjet::JetMedianBackgroundEstimator bge(jselector,jetDefkt,voronoi_area);
    fastjet::Subtractor subtractor(&bge);
    bge.set_particles(inputConst);
    double rho = bge.rho();
    subtractedJets = subtractor(inclusiveJets);
  }

  for(unsigned int iJet = 0 ; iJet < subtractedJets.size() ; iJet++){
    fastjet::PseudoJet jet = subtractedJets[iJet];
    if(jet.perp2()<5000.*5000.) continue;
    out_jets.push_back(jet);
    //vector<fastjet::PseudoJet> constituents = jet.constituents();
    //jetP->Set("area", jet.area());
    /*for(unsigned int iCons = 0; iCons < constituents.size(); iCons++){
      const PJ_Info* info = &(constituents[iCons].user_info<PJ_Info>());
      jetP->Add(ConsKey, info->Pointer);
      double ptcons = constituents[iCons].pt();
      double drcons = constituents[iCons].delta_R(jet);
      numwidth += drcons*ptcons;
      denwidth += ptcons;
    } // end loop over cons
    jetP->Set("width",numwidth/denwidth);
    jetP->Set("constscale_pt",jet.pt());*/
  }// end loop over jets

return EL::StatusCode::SUCCESS;
}

void VoronoiJets::print(fastjet::PseudoJet p){
  std::cout << "Eta: " << p.eta() << "; Phi: " << p.phi() << "; Pt: " << p.pt() << std::endl;
}

