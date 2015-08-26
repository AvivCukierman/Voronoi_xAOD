#ifndef MyAnalysis_VoronoiJets_H
#define MyAnalysis_VoronoiJets_H

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

class VoronoiJets : public EL::Algorithm
{
  // put your configuration variables here as public variables.
  // that way they can be set directly from CINT and python.
public:
  bool m_debug = false;
  std::string m_eventInfo      = "EventInfo",
              m_clust          = "CaloCalTopoClusters",
              m_jets          = "AntiKt4LCTopoJets";

   // methods used in the analysis
  virtual EL::StatusCode MakeJetsWArea(const fastjet::JetAlgorithm, const double,const std::vector<fastjet::PseudoJet>&, std::vector<fastjet::PseudoJet>&,bool doareasub=true);
  virtual void print(fastjet::PseudoJet);

  // variables that don't get filled at submission time should be
  // protected from being send from the submission node to the worker
  // node (done by the //!)
private:
  xAOD::TEvent *m_event; //!
  xAOD::TStore *m_store;  //!
  //Cluster collections:
  std::vector<fastjet::PseudoJet> clusters;//!
  std::vector<fastjet::PseudoJet> clusters_voronoi;//!
  std::vector<fastjet::PseudoJet> jets_voronoi;//!
public:
  // this is a standard constructor
  VoronoiJets ();

  
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
  ClassDef(VoronoiJets, 1);
};

#endif
