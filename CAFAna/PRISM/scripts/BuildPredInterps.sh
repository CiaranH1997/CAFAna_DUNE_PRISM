#!/bin/bash

# Build a pred interp script

INPUTDIR=/home/hasnipl/DP_DATA

echo "${INPUTDIR}"

source /home/hasnipl/CLionProjects/CAFAna_DUNE_PRISM/CAFAna/build/Linux/CAFAnaEnv.sh

MakePRISMPredInterps -o ${INPUTDIR}/PredInterps/PRISMState_EProxy_NomHC_25Mar.root \
        -N ${INPUTDIR}/CAF_FHC_PRISM_PROD4-14.root \
        -F ${INPUTDIR}/FD_FHC_nonswap.root \
        --bin-descriptor testopt --syst-descriptor nosyst \
        --no-fakedata-dials -A EProxy
        #--FakeSR "${INPUTDIR}/HCRat.root;bla" \
        #--OA-bin-descriptor "OneNegXBin"
