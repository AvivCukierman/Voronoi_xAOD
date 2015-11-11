#ifndef MyAnalysis_VoronoiWeights_H
#define MyAnalysis_VoronoiWeights_H

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

class VoronoiWeights : public EL::Algorithm
{
  // put your configuration variables here as public variables.
  // that way they can be set directly from CINT and python.
public:
    bool m_debug = false;
    std::string m_eventInfo      = "EventInfo",
                m_clust          = "CaloCalTopoClusters",
                m_jets          = "AntiKt4LCTopoJets";
    bool m_doLC = false;

   // methods used in the analysis
   struct PJcomp;
   virtual EL::StatusCode MakeVoronoiClusters(std::vector< std::pair<fastjet::PseudoJet,std::vector<float> > >&);
   void SpreadPt(std::vector< std::pair< fastjet::PseudoJet,std::vector<float> > >& correctedptvec, float spreadr=0.4, float alpha=2);

private:
  xAOD::TEvent *m_event; //!
  xAOD::TStore *m_store;  //!
  std::vector<fastjet::PseudoJet> clusters; //!

public:
  // this is a standard constructor
  VoronoiWeights ();

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
  ClassDef(VoronoiWeights, 1);
};

#endif
