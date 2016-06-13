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
              m_jets           = "AntiKt4NoAreaJets",
              m_truth_jets     = "AntiKt4TruthJets",
              m_voronoi0_jets   = "AntiKt4Voronoi0Jets",
              m_voronoi1_jets   = "AntiKt4Voronoi1Jets",
              m_voronois_jets   = "AntiKt4VoronoiSpreadJets",
              m_vertices       = "PrimaryVertices";
  bool m_doLC = false;

   // methods used in the analysis
  EL::StatusCode FillTJetVars(const DataVector<xAOD::Jet_v1>* tjets,
                                        std::vector<float>& tjetpt,
                                        std::vector<float>& tjeteta,
                                        std::vector<float>& tjetphi,
                                        std::vector<float>& tjetmass,
                                        std::vector<float>& tjetwidth,
                                        std::vector<float>& tjetmindr
                                        );

  EL::StatusCode FillJetVars(const DataVector<xAOD::Jet_v1>* jets,
                              const DataVector<xAOD::Jet_v1>* tjets,
                                        std::vector<float>& jetpt,
                                        std::vector<float>& jeteta,
                                        std::vector<float>& jetphi,
                                        std::vector<float>& jetmass,
                                        std::vector<float>& jetwidth,
                                        std::vector<float>& jetmindr,
                                        std::vector<bool>& jetisPU,
                                        std::vector<float>& cl0pt,
                                        std::vector<float>& tjetpt,
                                        std::vector<float>& tjeteta,
                                        std::vector<float>& tjetphi,
                                        std::vector<float>& tjetmass,
                                        std::vector<float>& tjetwidth,
                                        std::vector<float>& tjetmindr,
                                        bool doareasub=false);


  // variables that don't get filled at submission time should be
  // protected from being send from the submission node to the worker
  // node (done by the //!)
private:
  xAOD::TEvent *m_event; //!
  xAOD::TStore *m_store;  //!
  TTree *m_tree; //!

  //everything below here is filled in the tree
  int m_eventNumber; //!
  float m_eventWeight; //!
  int m_NPV; //!
  int m_mu; //!
  float m_rho; //!

  std::vector<float> m_tjpt; //!
  std::vector<float> m_tjeta; //!
  std::vector<float> m_tjphi; //!
  std::vector<float> m_tjmass; //!
  std::vector<float> m_tjwidth; //!
  std::vector<float> m_tjmindr; //!

  std::vector<float> m_jvoro0pt; //!
  std::vector<float> m_jvoro0eta; //!
  std::vector<float> m_jvoro0phi; //!
  std::vector<float> m_jvoro0mass; //!
  std::vector<float> m_jvoro0width; //!
  std::vector<float> m_jvoro0mindr; //!
  std::vector<bool> m_jvoro0isPU; //!
  std::vector<float> m_jvoro0cl0pt; //!
  std::vector<float> m_tjvoro0pt; //!
  std::vector<float> m_tjvoro0eta; //!
  std::vector<float> m_tjvoro0phi; //!
  std::vector<float> m_tjvoro0mass; //!
  std::vector<float> m_tjvoro0width; //!
  std::vector<float> m_tjvoro0mindr; //!

  std::vector<float> m_jvoro1pt; //!
  std::vector<float> m_jvoro1eta; //!
  std::vector<float> m_jvoro1phi; //!
  std::vector<float> m_jvoro1mass; //!
  std::vector<float> m_jvoro1width; //!
  std::vector<float> m_jvoro1mindr; //!
  std::vector<bool> m_jvoro1isPU; //!
  std::vector<float> m_jvoro1cl0pt; //!
  std::vector<float> m_tjvoro1pt; //!
  std::vector<float> m_tjvoro1eta; //!
  std::vector<float> m_tjvoro1phi; //!
  std::vector<float> m_tjvoro1mass; //!
  std::vector<float> m_tjvoro1width; //!
  std::vector<float> m_tjvoro1mindr; //!

  std::vector<float> m_jvorospt; //!
  std::vector<float> m_jvoroseta; //!
  std::vector<float> m_jvorosphi; //!
  std::vector<float> m_jvorosmass; //!
  std::vector<float> m_jvoroswidth; //!
  std::vector<float> m_jvorosmindr; //!
  std::vector<bool> m_jvorosisPU; //!
  std::vector<float> m_jvoroscl0pt; //!
  std::vector<float> m_tjvorospt; //!
  std::vector<float> m_tjvoroseta; //!
  std::vector<float> m_tjvorosphi; //!
  std::vector<float> m_tjvorosmass; //!
  std::vector<float> m_tjvoroswidth; //!
  std::vector<float> m_tjvorosmindr; //!

  std::vector<float> m_jnoarea0pt; //!
  std::vector<float> m_jnoarea0eta; //!
  std::vector<float> m_jnoarea0phi; //!
  std::vector<float> m_jnoarea0mass; //!
  std::vector<float> m_jnoarea0width; //!
  std::vector<float> m_jnoarea0mindr; //!
  std::vector<bool> m_jnoarea0isPU; //!
  std::vector<float> m_jnoarea0cl0pt; //!
  std::vector<float> m_tjnoarea0pt; //!
  std::vector<float> m_tjnoarea0eta; //!
  std::vector<float> m_tjnoarea0phi; //!
  std::vector<float> m_tjnoarea0mass; //!
  std::vector<float> m_tjnoarea0width; //!
  std::vector<float> m_tjnoarea0mindr; //!

  std::vector<float> m_j0pt; //!
  std::vector<float> m_j0eta; //!
  std::vector<float> m_j0phi; //!
  std::vector<float> m_j0mass; //!
  std::vector<float> m_j0width; //!
  std::vector<float> m_j0mindr; //!
  std::vector<bool> m_j0isPU; //!
  std::vector<float> m_j0cl0pt; //!
  std::vector<float> m_tj0pt; //!
  std::vector<float> m_tj0eta; //!
  std::vector<float> m_tj0phi; //!
  std::vector<float> m_tj0mass; //!
  std::vector<float> m_tj0width; //!
  std::vector<float> m_tj0mindr; //!
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
