#ifndef MyAnalysis_MyxAODAnalysis_H
#define MyAnalysis_MyxAODAnalysis_H

#include <EventLoop/Algorithm.h>
#include "EventLoopAlgs/NTupleSvc.h"
#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TEvent.h"

#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TEvent.h"
#include "xAODRootAccess/TStore.h"
#include "TRegexp.h"
#include <fstream>
#include <iostream>

#ifndef __CINT__
#include "xAODCore/ShallowCopy.h"
#include "xAODEventInfo/EventInfo.h"
#include "xAODTracking/VertexContainer.h"
#include "xAODJet/JetContainer.h"
#include "xAODJet/JetAuxContainer.h"
#include "xAODCaloEvent/CaloClusterContainer.h"
#include "xAODJet/JetConstituentVector.h"
#endif

#include "xAODTracking/VertexContainer.h"
#include "JetCalibTools/JetCalibrationTool.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/JetDefinition.hh"

#include "PATInterfaces/SystematicRegistry.h"

#include <TTree.h>

//#include "MyAnalysis/Candidate.h"
//#include "MyAnalysis/constituentTests.h"

class MyxAODAnalysis : public EL::Algorithm
{
  // put your configuration variables here as public variables.
  // that way they can be set directly from CINT and python.
public:
  // float cutValue;



  // variables that don't get filled at submission time should be
  // protected from being send from the submission node to the worker
  // node (done by the //!)
public:
  // Tree *myTree; //!
  // TH1 *myHist; //!

  std::string outputName;
  TTree *tree; //!
  int EventNumber; //!


  // this is a standard constructor
  MyxAODAnalysis ();

  xAOD::TEvent *m_event; //!
  xAOD::TStore m_store;  //!
  int m_eventCounter; //!

  EL::NTupleSvc *m_ntupleSvc; //!

   // things read from the xAOD and made accessible to functions
#ifndef __CINT__
   const xAOD::EventInfo *m_eventInfo; //!
   const xAOD::VertexContainer *m_vertices; //!
   const xAOD::JetContainer *m_lc_jets; //!
   const xAOD::JetContainer *m_truth_jets;//!
   std::pair< xAOD::JetContainer*, xAOD::ShallowAuxContainer* > m_my_truth_jets;//!
   const xAOD::CaloClusterContainer *m_clust;//!
#endif

  //Jet collections:
  std::vector<fastjet::PseudoJet> jets_AntiKt4LCTopo;//!
  std::vector<fastjet::PseudoJet> jets_truth;//!
  std::vector<fastjet::PseudoJet> jets_voronoi;//!

  //Cluster collections:
  std::vector<fastjet::PseudoJet> clusters;//!
  std::vector<fastjet::PseudoJet> clusters_voronoi;//!
  
  //Branches:
  std::vector<float> j0pt; //!
  std::vector<float> j0eta; //!
  std::vector<float> j0phi; //!
  std::vector<float> j0m; //!

  std::vector<float> jvoro0pt; //!
  std::vector<float> jvoro0eta; //!
  std::vector<float> jvoro0phi; //!
  std::vector<float> jvoro0m; //!

  std::vector<float> truejetpt; //!
  std::vector<float> truejeteta; //!
  std::vector<float> truejetphi; //!
  std::vector<float> truejetm; //!

  //Event info:
  int eventMu;//!
  int eventNPV;//!

   // methods used in the analysis
   virtual EL::StatusCode getObjects();
   virtual EL::StatusCode MakeVoronoiClusters();
   virtual EL::StatusCode MakeJetsWArea(const fastjet::JetAlgorithm, const double,std::vector<fastjet::PseudoJet>, std::vector<fastjet::PseudoJet>&,bool doareasub=true);
   virtual EL::StatusCode FillJetVars(std::vector<fastjet::PseudoJet> jets,
                                             std::vector<float>& jetpt,
                                             std::vector<float>& jeteta,
                                             std::vector<float>& jetphi,
                                             std::vector<float>& jetm);
   virtual EL::StatusCode FillTruthJetVars(std::vector<fastjet::PseudoJet> jets,
                                             std::vector<float>& jetpt,
                                             std::vector<float>& jeteta,
                                             std::vector<float>& jetphi,
                                             std::vector<float>& jetm);
   virtual void print(fastjet::PseudoJet);
   //virtual EL::StatusCode fillOutputTree();
   //virtual EL::StatusCode fillJets();
   //virtual EL::StatusCode NormalClusters();

  // these are the functions inherited from Algorithm
  virtual EL::StatusCode setupJob (EL::Job& job);
  virtual EL::StatusCode fileExecute ();
  virtual EL::StatusCode histInitialize ();
  virtual EL::StatusCode changeInput (bool firstFile);
  virtual EL::StatusCode initialize ();
  virtual EL::StatusCode execute ();
  virtual EL::StatusCode postExecute ();
  virtual EL::StatusCode finalize ();
  virtual EL::StatusCode histFinalize ();

  // this is needed to distribute the algorithm to the workers
  ClassDef(MyxAODAnalysis, 1);
};

#endif
