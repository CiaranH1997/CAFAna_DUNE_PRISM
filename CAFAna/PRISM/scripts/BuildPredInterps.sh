#!/bin/bash

# Build a pred interp script

INPUTDIR=/home/hasnipl/DP_DATA

echo "${INPUTDIR}"

source /home/hasnipl/CLionProjects/NewCode_CAFAna/CAFAna_DUNE_PRISM/CAFAna/build/Linux/CAFAnaEnv.sh

MakePRISMPredInterps -o ${INPUTDIR}/NewCode_PredInterps/PRISMState_EProxy_AltHC_NewProd_10Aug.root \
        -N-nu ${INPUTDIR}/ND_CAF_FHC_FVCut_NewProd8Aug.root \
        -F-nu ${INPUTDIR}/FD_FHC_nonswap.root \
        --bin-descriptor testopt --no-fakedata-dials \
        -A EProxy --syst-descriptor "nosyst" \
        --FakeSR "OnAxis280kA" --no-fakedata-dials \
        --OA-bin-descriptor "OneNegXBin"


