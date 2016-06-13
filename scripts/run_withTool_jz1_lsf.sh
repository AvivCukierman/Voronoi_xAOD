#!/bin/bash
python MyAnalysis/scripts/Run_withTool.py \
    --doLC \
    --submitDir "New/LSF_LC_JZ1_jvoro" \
    -w \
    --inputDQ2 \
    --inputFiles "mc15_13TeV.361021.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ1W.merge.AOD.e3569_s2576_s2132_r6765_r6282/" \
    --driver "lsf"
exit 0
