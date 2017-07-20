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

#include <VoronoiWeightTool/VoronoiWeightTool.h>

//Reclustering
#include <JetReclustering/JetReclusteringTool.h>

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

  int EventNumber; //!
  bool m_debug = false;
  bool m_doLC = false;
  std::array<VoronoiWeightTool*, 3> m_voronoiWeightTools; //!
  std::array<JetReclusteringTool*, 3> m_jetReclusteringTools; //!

  std::string m_inputContainer = "CaloCalTopoClusters",
              m_outputContainer = "VoronoiClusters";

  EL::StatusCode FindRho(const xAOD::CaloClusterContainer* in_clusters,float& rho, float& sigma);
  EL::StatusCode SetRho(DataVector<xAOD::Jet_v1> jets,float event_rho);


  // this is a standard constructor
  MyxAODAnalysis ();

  xAOD::TEvent *m_event; //!
  xAOD::TStore *m_store;  //!
  int m_eventCounter; //!

   // things read from the xAOD and made accessible to functions
#ifndef __CINT__
   const xAOD::EventInfo *m_eventInfo; //!
#endif
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
