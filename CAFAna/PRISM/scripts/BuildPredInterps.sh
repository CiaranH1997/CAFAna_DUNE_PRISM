#!/bin/bash

# Build a pred interp script

INPUTDIR=/dune/data/users/chasnip/OffAxisCAFs

echo "${INPUTDIR}"

source /dune/app/users/chasnip/CH_DUNE_PRISM/lblpwgtools/CAFAna/build/Linux/CAFAnaEnv.sh

MakePRISMPredInterps -o ${INPUTDIR}/../PredInterps/PRISMState_EProxy_AltHC_16Apr_Flux_XSec_Syst.root \
        -N ${INPUTDIR}/CAF_FHC_PRISM_PROD4-14.root \
        -F ${INPUTDIR}/FD_FHC_nonswap.root \
        --bin-descriptor testopt  \
        -A EProxy --syst-descriptor "nov17flux:nodet" \
        --FakeSR "${INPUTDIR}/HCRat.root;bla" \
        --OA-bin-descriptor "OneNegXBin"

#--no-fakedata-dials
#--syst-descriptor fakedata
#--NueSwap ${INPUTDIR}/FD_FHC_nueswap.root \
#--TauSwap ${INPUTDIR}/FD_FHC_tauswap.root \
#--PRISM-fake-data "MissingProtonFakeData"
