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



MyxAODAnalysis :: MyxAODAnalysis ():
   m_event(0),
   m_ntupleSvc(0),
   m_eventInfo(0),
   m_vertices(0)
{}

EL::StatusCode MyxAODAnalysis :: setupJob (EL::Job& job)
{
  // let's initialize the algorithm to use the xAODRootAccess package
  job.useXAOD ();
  xAOD::Init(); // call before opening first file

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode MyxAODAnalysis :: histInitialize () {

  // get the output file, create a new TTree and connect it to that output
  // define what braches will go in that tree
  TFile *outputFile = wk()->getOutputFile (outputName);
  tree = new TTree ("tree", "tree");
  tree->SetDirectory (outputFile);
  tree->Branch("EventNumber", &EventNumber); 
  tree->Branch("EventCounter", &m_eventCounter); 
  tree->Branch("NPV", &eventNPV); 
  tree->Branch("Mu", &eventMu); 

  tree->Branch("j0pt","std::vector<float>",&j0pt);
  tree->Branch("j0eta","std::vector<float>",&j0eta);
  tree->Branch("j0phi","std::vector<float>",&j0phi);
  tree->Branch("j0m","std::vector<float>",&j0m);

  tree->Branch("jvoro0pt","std::vector<float>",&jvoro0pt);
  tree->Branch("jvoro0eta","std::vector<float>",&jvoro0eta);
  tree->Branch("jvoro0phi","std::vector<float>",&jvoro0phi);
  tree->Branch("jvoro0m","std::vector<float>",&jvoro0m);

  tree->Branch("truejetpt","std::vector<float>",&truejetpt);
  tree->Branch("truejeteta","std::vector<float>",&truejeteta);
  tree->Branch("truejetphi","std::vector<float>",&truejetphi);
  tree->Branch("truejetm","std::vector<float>",&truejetm);

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode MyxAODAnalysis :: fileExecute () {return EL::StatusCode::SUCCESS;}

EL::StatusCode MyxAODAnalysis :: changeInput (bool /*firstFile*/) {return EL::StatusCode::SUCCESS;}

EL::StatusCode MyxAODAnalysis :: initialize ()
{
  m_event = wk()->xaodEvent();
  m_eventCounter = 0;
  // as a check, let's see the number of events in our xAOD
  Info("initialize()", "Number of events = %lli", m_event->getEntries() ); // print long long int

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode MyxAODAnalysis :: execute ()
{
  const char* APP_NAME = "MyxAODAnalysis::execute()";
  m_eventCounter++;

  if( ! m_event->retrieve( m_eventInfo, "EventInfo").isSuccess() ){
    Error(APP_NAME, "Failed to retrieve event info collection. Exiting." );
    return EL::StatusCode::FAILURE;
   }
  // fill the branches of our trees
  EventNumber = m_eventInfo->eventNumber();

  if (getObjects() != EL::StatusCode::SUCCESS) {
    Error(APP_NAME, "Error in getObjects");
    return EL::StatusCode::FAILURE;
  }
  FillJetVars(jets_AntiKt4LCTopo,j0pt,j0eta,j0phi,j0m);
  FillJetVars(jets_voronoi,jvoro0pt,jvoro0eta,jvoro0phi,jvoro0m);
  FillTruthJetVars(jets_truth,truejetpt,truejeteta,truejetphi,truejetm);

  tree->Fill();
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

EL::StatusCode MyxAODAnalysis::getObjects()
{
  clusters.clear();
  clusters_voronoi.clear();

  jets_AntiKt4LCTopo.clear();
  jets_voronoi.clear();
  jets_truth.clear();

  const char *APP_NAME = "MyxAODAnalysis::getObjects";


  //////////////////////////
  //  Clusters
  //////////////////////////


  if (!m_event->retrieve(m_clust, "CaloCalTopoClusters").isSuccess()) {    // retrieve arguments: container type, container key
     Error(APP_NAME, "Failed to retrieve CaloCalTopoClusters container");
     return EL::StatusCode::FAILURE;
  }

   for(auto *clust: *m_clust){
    fastjet::PseudoJet test = fastjet::PseudoJet( clust->p4());
    if(test.E() >= 0) clusters.push_back(test);
   } 
   MakeJetsWArea(fastjet::antikt_algorithm, 0.4, clusters,jets_AntiKt4LCTopo,true);

   MakeVoronoiClusters();
   MakeJetsWArea(fastjet::antikt_algorithm, 0.4, clusters_voronoi,jets_voronoi,false);

   ///////////
   // Jets
   ///////////
   /*m_lc_jets = 0;

   if (!m_event->retrieve(m_lc_jets, "AntiKt4LCTopoJets").isSuccess()) {    // retrieve arguments: container type, container key
      Error(APP_NAME, "Failed to retrieve AntiKt4LCTopo Jet container");
      return EL::StatusCode::FAILURE;
   }*/
   /*for(auto *xjet: *m_lc_jets){
    fastjet::PseudoJet jet = fastjet::PseudoJet(xjet->p4());
     if(jet.pt() < 5000) continue;
     print(jet);
   } */

   if (!m_event->retrieve(m_truth_jets, "AntiKt4TruthJets").isSuccess()) {    // retrieve arguments: container type, container key
      return EL::StatusCode::FAILURE;
   }
   for(auto *xjet: *m_truth_jets){
    fastjet::PseudoJet jet = fastjet::PseudoJet(xjet->p4());
    jets_truth.push_back(jet);
   } 

   m_vertices = 0;
   if (!m_event->retrieve(m_vertices, "PrimaryVertices").isSuccess()) {    // retrieve arguments: container type, container key
      return EL::StatusCode::FAILURE;
   }

   eventMu = m_eventInfo->averageInteractionsPerCrossing();
   eventNPV = 0;
   for ( auto *ivert : *m_vertices ){
     if ( (ivert)->nTrackParticles() >= 2 ) ++eventNPV;
   }

   //Info("GetObjects()","Event number: %i; NPV: %i",m_eventCounter,eventNPV);   

   m_my_truth_jets = xAOD::shallowCopyContainer(*m_truth_jets);

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode MyxAODAnalysis::MakeVoronoiClusters(){
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

  for(unsigned int iJet = 0 ; iJet < inclusiveJets.size() ; iJet++){
        fastjet::PseudoJet jet = inclusiveJets[iJet];
        std::vector<fastjet::PseudoJet> constituents = jet.constituents();
        for(unsigned int iCons = 0; iCons < constituents.size(); iCons++){
                fastjet::PseudoJet constituentP;
                float area = constituents[iCons].area();
                float rho = bge.rho();
                float correctedPt = constituents[iCons].pt()-rho*area;
       //         if(m_eventCounter==1){std::cout << "Area: " << area << "; Rho: " << bge.rho() << "; pt: " << constituents[iCons].pt() << "; corrected: " << correctedPt << std::endl;}
                if(correctedPt<0.) continue;
                constituentP.reset_PtYPhiM(correctedPt, constituents[iCons].rap(), constituents[iCons].phi(), constituents[iCons].m());
                clusters_voronoi.push_back(constituentP);
        } // end loop over cons
  }// end loop over jets

return EL::StatusCode::SUCCESS;
}

EL::StatusCode MyxAODAnalysis::MakeJetsWArea(const fastjet::JetAlgorithm algo, const double jetR, std::vector<fastjet::PseudoJet> inputConst, std::vector<fastjet::PseudoJet>& out_jets,bool doareasub){

  fastjet::JetDefinition jetDef(algo, jetR,fastjet::E_scheme, fastjet::Best);
  fastjet::AreaDefinition active_area = fastjet::AreaDefinition(fastjet::active_area);
  fastjet::ClusterSequenceArea clustSeq(inputConst, jetDef, active_area);
  std::vector<fastjet::PseudoJet> inclusiveJets = sorted_by_pt(clustSeq.inclusive_jets(5.));
  std::vector<fastjet::PseudoJet> subtractedJets = inclusiveJets;

  //incomplete
  if(doareasub){
    fastjet::Selector jselector = fastjet::SelectorAbsRapRange(0.0,2.1);
    fastjet::JetDefinition jetDefkt(fastjet::kt_algorithm,0.4);
    fastjet::AreaDefinition voronoi_area = fastjet::AreaDefinition(fastjet::voronoi_area,fastjet::VoronoiAreaSpec(0.9));
    fastjet::JetMedianBackgroundEstimator bge(jselector,jetDefkt,voronoi_area);
    fastjet::Subtractor subtractor(&bge);
    bge.set_particles(inputConst);
    //double rho = bge.rho();
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

EL::StatusCode MyxAODAnalysis::FillJetVars(std::vector<fastjet::PseudoJet> jets,
                                            std::vector<float>& jetpt,
                                            std::vector<float>& jeteta,
                                            std::vector<float>& jetphi,
                                            std::vector<float>& jetm)
{
  jetpt.clear();
  jeteta.clear();
  jetphi.clear();
  jetm.clear();

  for(auto jet: jets){
    jetpt.push_back(jet.pt());
    jeteta.push_back(jet.eta());
    jetphi.push_back(jet.phi());
    jetm.push_back(jet.m());
  }

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode MyxAODAnalysis::FillTruthJetVars(std::vector<fastjet::PseudoJet> jets,
                                            std::vector<float>& jetpt,
                                            std::vector<float>& jeteta,
                                            std::vector<float>& jetphi,
                                            std::vector<float>& jetm)
{
  jetpt.clear();
  jeteta.clear();
  jetphi.clear();
  jetm.clear();

  for(auto jet: jets){
    jetpt.push_back(jet.pt());
    jeteta.push_back(jet.eta());
    jetphi.push_back(jet.phi());
    jetm.push_back(jet.m());
  }

  return EL::StatusCode::SUCCESS;
}
void MyxAODAnalysis::print(fastjet::PseudoJet p){
  std::cout << "Eta: " << p.eta() << "; Phi: " << p.phi() << "; Pt: " << p.pt() << std::endl;
}

