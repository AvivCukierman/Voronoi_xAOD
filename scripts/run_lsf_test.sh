#!/bin/bash
python MyAnalysis/scripts/Run.py \
    --submitDir "inputDQ2_test" \
    --inputDQ2 \
    --inputFiles "mc15_13TeV.361020.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ0W.merge.AOD.e3569_s2576_s2132_r6765_r6282/" \
    --driver "direct" \
    --nevents 10
exit 0
