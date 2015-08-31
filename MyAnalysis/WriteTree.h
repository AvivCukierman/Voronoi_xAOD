#ifndef MyAnalysis_WriteTree_H
#define MyAnalysis_WriteTree_H

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

class WriteTree : public EL::Algorithm
{
  // put your configuration variables here as public variables.
  // that way they can be set directly from CINT and python.
public:
  bool m_debug = false;
  std::string m_eventInfo      = "EventInfo",
              m_jets           = "AntiKt4LCTopoJets",
              m_truth_jets     = "AntiKt4TruthJets",
              m_voronoi_jets   = "AntiKt4VoronoiJets",
              m_vertices       = "PrimaryVertices";

   // methods used in the analysis
  EL::StatusCode FillJetVars(const DataVector<xAOD::Jet_v1>* jets,
                              const DataVector<xAOD::Jet_v1>* tjets,
                                        std::vector<float>& jetpt,
                                        std::vector<float>& jeteta,
                                        std::vector<float>& jetphi,
                                        std::vector<float>& jetmass,
                                        std::vector<float>& jetwidth,
                                        std::vector<float>& jetmindr,
                                        std::vector<float>& tjetpt,
                                        std::vector<float>& tjeteta,
                                        std::vector<float>& tjetphi,
                                        std::vector<float>& tjetmass,
                                        std::vector<float>& tjetwidth,
                                        std::vector<float>& tjetmindr);


  // variables that don't get filled at submission time should be
  // protected from being send from the submission node to the worker
  // node (done by the //!)
private:
  xAOD::TEvent *m_event; //!
  xAOD::TStore *m_store;  //!
  TTree *m_tree; //!

  //everything below here is filled in the tree
  int m_eventNumber; //!
  int m_NPV; //!
  int m_mu; //!

  std::vector<float> m_jvoro0pt; //!
  std::vector<float> m_jvoro0eta; //!
  std::vector<float> m_jvoro0phi; //!
  std::vector<float> m_jvoro0mass; //!
  std::vector<float> m_jvoro0width; //!
  std::vector<float> m_jvoro0mindr; //!
  std::vector<float> m_tjvoro0pt; //!
  std::vector<float> m_tjvoro0eta; //!
  std::vector<float> m_tjvoro0phi; //!
  std::vector<float> m_tjvoro0mass; //!
  std::vector<float> m_tjvoro0width; //!
  std::vector<float> m_tjvoro0mindr; //!
public:
  // this is a standard constructor
  WriteTree ();

  
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
  ClassDef(WriteTree, 1);
};

#endif
