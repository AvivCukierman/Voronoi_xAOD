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

EL::StatusCode VoronoiWeights :: changeInput (bool firstFile) {return EL::StatusCode::SUCCESS;}

EL::StatusCode VoronoiWeights :: initialize ()
{
  m_event = wk()->xaodEvent();
  m_store = wk()->xaodStore();

  // as a check, let's see the number of events in our xAOD
  Info("initialize()", "Number of events = %lli", m_event->getEntries() ); // print long long int

  return EL::StatusCode::SUCCESS;
}

const double PI = 3.14159265358979323844;
double mod2pi(double phi){
  while(phi > 2.0*PI) phi -= 2.0*PI;
  while(phi < 0) phi += 2.0*PI;
  return phi;
}
//Have to define custom comparator for PseudoJets in order to have a map from PJs to anything
//Comparison is fuzzy to account for rounding errors
struct VoronoiWeights :: PJcomp {
  bool operator() (const fastjet::PseudoJet& lhs, const fastjet::PseudoJet& rhs) const
  {
    if(lhs.pt()<rhs.pt()-1) return true;
    if(lhs.pt()>rhs.pt()+1) return false;
    if(lhs.eta()<rhs.eta()-.01) return true;
    if(lhs.eta()>rhs.eta()+.01) return false;
    if(mod2pi(lhs.phi())<mod2pi(rhs.phi())-.01) return true; 
    return false;
    //The comparator must be a strict weak ordering. 
    //Might be an issue if there are two clusters right next to each other with the exact same Pt - highly unlikely
  }
};

EL::StatusCode VoronoiWeights :: execute ()
{
  const char* APP_NAME = "VoronoiWeights::execute()";

  const xAOD::EventInfo*                        eventInfo     (nullptr);
  const xAOD::CaloClusterContainer*             in_clusters   (nullptr);
  const xAOD::JetContainer*                     in_jets       (nullptr);
  
  // start grabbing all the containers that we can
  RETURN_CHECK("VoronoiWeights::execute()", HF::retrieve(eventInfo,    m_eventInfo,        m_event, m_store, m_debug), "Could not get the EventInfo container.");
  if(!m_clust.empty())
    RETURN_CHECK("VoronoiWeights::execute()", HF::retrieve(in_clusters,     m_clust,       m_event, m_store, m_debug), "Could not get the clusters container.");
  if(!m_jets.empty())
    RETURN_CHECK("VoronoiWeights::execute()", HF::retrieve(in_jets,     m_jets,       m_event, m_store, m_debug), "Could not get the jets container.");

  clusters.clear();
  for(const auto clust: *in_clusters){
    fastjet::PseudoJet test = fastjet::PseudoJet(clust->p4()); //read in clusters as PseudoJets
    if(test.E() >= 0) clusters.push_back(test);
  }
  std::map<fastjet::PseudoJet,float,PJcomp> ptmap; //Map from PseudoJets to their corrected (Voronoi subtractd) Pts
  if(MakeVoronoiClusters(ptmap) != EL::StatusCode::SUCCESS) Error(APP_NAME,"Error in MakeVoronoiClusters");

  /*std::cout << "Size after: " << ptmap.size() << std::endl;
  for (std::map<fastjet::PseudoJet,float,PJcomp>::iterator it=ptmap.begin(); it!=ptmap.end(); ++it){
      fastjet::PseudoJet cons = it->first;
//Eta: -1.59458; Phi: 2.9978
      if(fabs(cons.eta()+1.59458)<0.1 && fabs(cons.phi()-2.9978)<0.1){
        std::cout << "Pt: " << cons.pt() << "; Eta: " << cons.eta() <<"; Phi: " << cons.phi() << std::endl;
        std::cout << ptmap.count(cons) << "; " << ptmap[cons] << std::endl;
      }
  }*/
  static SG::AuxElement::Decorator< float > correctedPt("correctedPt");
  for(const auto clust: *in_clusters){
    correctedPt(*clust) = 0;
    fastjet::PseudoJet pjclust = fastjet::PseudoJet(clust->p4());
    //std::cout << ptmap.count(pjclust) << std::endl;
    float subpt=0;
    if(ptmap.count(pjclust)==0){ //Should mean E<0. Occasionally means the comparator fucked up (about one cluster every hundred events).
      if(pjclust.E()<0) continue; //If E<0, correctedPt = 0 - don't use that cluster.
      else{ 
        //Problem with the comparator. We have to go back through the map and search manually for the cluster. This means that for the vast majority of clusters the setting of the corrected Pt is O(logn) but for a few rare ones it's O(n).
        //std::cout << "Pt: " << pjclust.pt() << "; Eta: " << pjclust.eta() <<"; Phi: " << pjclust.phi() << std::endl;
        bool found = 0; 
        for (std::map<fastjet::PseudoJet,float,PJcomp>::iterator it=ptmap.begin(); it!=ptmap.end(); ++it){
          fastjet::PseudoJet cons = it->first;
          if(fabs(cons.eta()-pjclust.eta())<0.01 && fabs(mod2pi(cons.phi())-mod2pi(pjclust.phi()))<0.01 && fabs(cons.pt()-pjclust.pt())<1){ //Literally exactly the same comparator. Only the C++ gods know why it didn't work the first time.
            found = 1;
            subpt = it->second;
            break;
          }
        }
        if(!found) Error(APP_NAME,"Cluster with E>0 with no Voronoi area");
      }
    }
    else subpt = ptmap[pjclust];
    if(subpt < 0) continue;
    correctedPt(*clust) = subpt;
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

EL::StatusCode VoronoiWeights::MakeVoronoiClusters(std::map<fastjet::PseudoJet,float,PJcomp>& correctedptmap){
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
      correctedptmap[cons] = correctedPt;
    } // end loop over cons
  } // end loop over jets
  //std::cout << "Size: " << correctedptmap.size() << std::endl;

return EL::StatusCode::SUCCESS;
}

