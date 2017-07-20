#include <vector>
#include <map>
#include <math.h>
#include <ios>

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
//#include "xAODAnaHelpers/tools/ReturnCheck.h"

//ANN:
//#include "ANN/ANN.h"

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
  bool operator() (const std::pair<fastjet::PseudoJet, std::vector<float> >& lhsp, const std::pair<fastjet::PseudoJet, std::vector<float> >& rhsp)
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
  ANA_CHECK(HF::retrieve(eventInfo,    m_eventInfo,        m_event, m_store, msg()));
  if(!m_clust.empty())
    ANA_CHECK(HF::retrieve(in_clusters,     m_clust,       m_event, m_store, msg()));

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

  std::vector< std::pair< fastjet::PseudoJet, std::vector<float> > > ptvec; //vector of pairs of PJs and their corrected pTs
  if(MakeVoronoiClusters(ptvec) != EL::StatusCode::SUCCESS) Error(APP_NAME,"Error in MakeVoronoiClusters");
  std::sort(ptvec.begin(), ptvec.end(), PJcomp());

  int i=0;
  static SG::AuxElement::Decorator< float > voro0Pt("voro0Pt");
  static SG::AuxElement::Decorator< float > voro1Pt("voro1Pt");
  static SG::AuxElement::Decorator< float > spreadPt("spreadPt");
  for(const auto clust: HF::sort_container_pt(in_clusters)){
    voro0Pt(*clust) = 0;
    voro1Pt(*clust) = 0;
    spreadPt(*clust) = 0;
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
      if(m_debug) std::cout << clust->pt() << "\t" << clust->eta() << "\t" << clust->phi() << "\t" << clust->m() << "\t" << ptvec[i].second[0] << "\t" << ptvec[i].second[3] << std::endl;
      voro0Pt(*clust) = ptvec[i].second[1];
      voro1Pt(*clust) = ptvec[i].second[2];
      spreadPt(*clust) = ptvec[i].second[3];
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

EL::StatusCode VoronoiWeights::MakeVoronoiClusters(std::vector< std::pair< fastjet::PseudoJet,std::vector<float> > >& correctedptvec){
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

  int nsigma = 1;
  float rho = bge.rho();
  float sigma = bge.sigma();
  for(unsigned int iJet = 0 ; iJet < inclusiveJets.size() ; iJet++){
    fastjet::PseudoJet jet = inclusiveJets[iJet];
    std::vector<fastjet::PseudoJet> constituents = jet.constituents();
    for(auto cons : constituents){
      float pt = cons.pt();
      float area = cons.area();
      float subPt = pt-rho*area;
      //std::cout << "Area: " << area << "; Rho: " << bge.rho() << "; pt: " << constituents[iCons].pt() << "; corrected: " << correctedPt << std::endl;
      //std::cout << "Pt: " << cons.pt() << "; Eta: " << cons.eta() <<"; Phi: " << cons.phi() << std::endl;
      //fastjet::PseudoJet constituentP;
      /*if(correctedPt<0.) continue;
      constituentP.reset_PtYPhiM(correctedPt, constituents[iCons].rap(), constituents[iCons].phi(), constituents[iCons].m());
      clusters_voronoi.push_back(constituentP);*/
      //correctedptmap[cons] = correctedPt;
      float voro0pt = subPt * (subPt > 0);
      float voro1pt = subPt * (subPt > sqrt(area)*sigma*(float)nsigma);
      std::vector<float> algopts;
      algopts.push_back(subPt);
      algopts.push_back(voro0pt);
      algopts.push_back(voro1pt);
      algopts.push_back(0);
      std::pair <fastjet::PseudoJet,std::vector<float> > pjcptpair (cons,algopts);
      correctedptvec.push_back(pjcptpair);
    } // end loop over cons
  } // end loop over jets
  //std::cout << "Size: " << correctedptmap.size() << std::endl;

  SpreadPt(correctedptvec);

return EL::StatusCode::SUCCESS;
}

void VoronoiWeights::SpreadPt(std::vector< std::pair< fastjet::PseudoJet,std::vector<float> > >& correctedptvec, float spreadr, float alpha){
  const float PI = 3.14159265;
  //default alpha = 2
  //Set up neighbors within spreadr:
  int clusters = correctedptvec.size();
  std::vector<float> spreadPT(clusters);
  std::vector<bool> isPositive(clusters);
  for(int iCl = 0; iCl < clusters; iCl++){
    spreadPT[iCl] = correctedptvec[iCl].second[0];
    isPositive[iCl] = spreadPT[iCl]>0;
  }

  std::vector<std::vector<std::pair<int,float> > > cluster_drs; //for each cluster, list of nearby positive pT clusters and their distances
  for(int iCl = 0; iCl < clusters; iCl++){
    fastjet::PseudoJet icluster = correctedptvec[iCl].first;
    float ieta = icluster.eta();
    float iphi = icluster.phi();
    std::vector<std::pair<int,float> > this_cluster_drs;
    for(int jCl = 0; jCl < clusters; jCl++){
      if(iCl == jCl) continue;
      if(!isPositive[jCl]) continue;
      fastjet::PseudoJet jcluster = correctedptvec[jCl].first;
      float jeta = jcluster.eta();
      float jphi = jcluster.phi();
      float dphi = icluster.delta_phi_to(jcluster);
      float deta = icluster.eta() - jcluster.eta(); //fastjet::pseudojet::delta_R(const PseudoJet& other) gives rap-phi distance
      float dr2 = pow(dphi,2) + pow(deta,2);
      if(dr2 > pow(spreadr,2)) continue;
      std::pair<int,float> jdrpair (jCl,dr2);
      this_cluster_drs.push_back(jdrpair);
    }
    cluster_drs.push_back(this_cluster_drs);
  }

  for(int i = 0; i < clusters; i++){
    if(!(spreadPT[i]<0)) continue; //only spread from negative pT clusters
    //find closest positive PT cluster:
    float sumdR2 = 0;
    //iterate over nearby positive pT clusters
    for(int j=0; j<cluster_drs[i].size(); j++){
      //cout << "j: " << j << " realid: " << realid << " Eta: " << points[realid][0]<< " Phi: " << points[realid][1] << " Pt:" << spreadPT[realid] << " Dist: " << dists[j] << endl;  // dists[j] = dR^2
      float dr = cluster_drs[i][j].second;
      if(dr>0) sumdR2 += 1./(pow(dr,alpha/2));
    }
    //if at least one neighbor
    if(sumdR2 > 0){
      float spreadPT_orig = spreadPT[i];
      //std::cout << "orig: " << spreadPT_orig << std::endl;
      for(int j=0; j<cluster_drs[i].size(); j++){
        float dr = cluster_drs[i][j].second;
        float realid = cluster_drs[i][j].first;
        if(dr>0){
          float weight = (1./pow(dr,alpha/2))/sumdR2;
          //std::cout << dr << "; " << weight << std::endl;
          //std::cout << "Before spreading: " << weight << ";" << weight*spreadPT_orig << ";" << spreadPT[realid] << std::endl;
          if(fabs(weight*spreadPT_orig)>spreadPT[realid]){
            spreadPT[i]+=spreadPT[realid];
            spreadPT[realid]=0;
          }
          else{
            spreadPT[realid]+=weight*spreadPT_orig;
            spreadPT[i]-=weight*spreadPT_orig;
          }
          //std::cout << "After spreading: " << weight << ";" << weight*spreadPT_orig << ";" << spreadPT[realid] << std::endl;
        }
      }
      //std::cout << "final: "  << spreadPT[i] << std::endl;
    }
    //cout << i << ";" << cluster(i,key).Float("correctedPT") << ";" << spreadPT[i]<< endl;
  }

  /*float totalcorrpt=0, totalspreadpt=0;
    for(int i=0; i<clusters; i++){ totalcorrpt+=cluster(i,key).Float("correctedPT"); totalspreadpt+=spreadPT[i];}
    cout << totalcorrpt << ";" << totalspreadpt << endl; //should be the same*/

  for(int iCl = 0; iCl < clusters; iCl++){
    correctedptvec[iCl].second[3] = spreadPT[iCl] * (spreadPT[iCl] > 0);
  }
}
