#!/bin/bash
python MyAnalysis/scripts/Run.py \
    --submitDir "inputDQ2_test" \
    --inputDQ2 \
    --inputFiles "mc15_13TeV.361107.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zmumu.merge.AOD.e3601_s2576_s2132_r6765_r6282/" \
    --driver "direct" \
    -w \
    --nevents 10
exit 0
