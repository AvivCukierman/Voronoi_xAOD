#!/bin/bash
python MyAnalysis/scripts/Run.py \
    --submitDir "inputDQ2_test" \
    --inputDQ2 \
    --driver "lsf" \
    --nevents 10
exit 0
