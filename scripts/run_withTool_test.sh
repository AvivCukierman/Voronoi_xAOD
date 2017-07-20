#!/bin/bash
python MyAnalysis/scripts/Run_withTool.py \
    --submitDir "test_LC" \
    --doLC \
    -w \
    --inputFiles "/atlas/local/acukierm/mc15_13TeV.361022.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ2W.merge.AOD.e3668_s2576_s2132_r6765_r6282_tid05771189_00/" \
    --nevents 10 \
    direct

exit 0

#--inputFiles "/atlas/local/acukierm/dijetz1and2/" \
