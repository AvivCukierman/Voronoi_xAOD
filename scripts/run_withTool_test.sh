#!/bin/bash
python MyAnalysis/scripts/Run_withTool.py \
    --submitDir "test_LC" \
    -w \
    --inputFiles "/atlas/local/acukierm/dijetz1and2/" \
    --doLC \
    --nevents 10 \
    direct

exit 0

#--inputFiles "/atlas/dq2/user/bnachman/" \
