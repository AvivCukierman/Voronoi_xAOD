#!/bin/bash
python MyAnalysis/scripts/Run.py \
    --submitDir "test_LC" \
    --inputFiles "/atlas/local/acukierm/dijetz1and2/" \
    -w \
    --doLC \
    --nevents 10 \
    direct

exit 0
