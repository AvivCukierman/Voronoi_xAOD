#!/bin/bash
python MyAnalysis/scripts/Run.py \
    --submitDir "test" \
    --inputFiles "/atlas/local/acukierm/dijetz1and2/" \
    -w \
    --nevents 10 \
    direct

exit 0
