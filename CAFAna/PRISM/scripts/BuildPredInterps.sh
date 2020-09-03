#!/bin/bash

# Build a pred interp script

INPUTDIR=/dune/data/users/chasnip/

echo "${INPUTDIR}"

#source /dune/app/users/chasnip/CAFAna_DUNE_PRISM/CAFAna/build/Linux/CAFAnaEnv.sh

MakePRISMPredInterps -o ${INPUTDIR}/NewCode/PredInterps/PRISMState_AltHC_OldProd_ChargeHadBiasTest_28Aug.root \
        -N-nu ${INPUTDIR}/OffAxisCAFs/CAF_FHC_PRISM_PROD4-14.root \
        -F-nu ${INPUTDIR}/OffAxisCAFs/FD_FHC_nonswap.root \
        --bin-descriptor testopt --no-fakedata-dials \
        -A EProxy --syst-descriptor "nosyst" \
        --FakeSR "OnAxis280kA" \
        --OA-bin-descriptor "OneNegXBin"


