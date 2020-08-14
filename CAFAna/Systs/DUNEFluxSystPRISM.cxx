//
// Created by C. Hasnip on 25/06/2020.
//
#include "DUNEFluxSystPRISM.h"
#ifdef USE_TH2JAGGED
#include "CAFAna/PRISM/EffectiveFluxUncertaintyHelper.h"
static PRISM::EffectiveFluxUncertaintyHelper const *fOffAxisFluxSystHelper =
  nullptr;
#endif

#include "CAFAna/Core/Utilities.h"

#include "StandardRecord/StandardRecord.h"

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"

#include <cassert>

namespace ana {

  //----------------------------------------------------------------------
  DUNEFluxSystPRISM::~DUNEFluxSystPRISM() {}

  //-----------------------------------------------------------------------

  void DUNEFluxSystPRISM::Shift(double sigma, Restorer &restore, caf::StandardRecord *sr,
                                double &weight) const {

    if (abs(sr->dune.nuPDGunosc) == 16) {
      return;
    }

    #ifdef USE_TH2JAGGED

    weight = fOffAxisFluxSystHelper->GetFluxWeight(
      fIdx, sigma, sr->dune.OffAxisFluxBin, sr->dune.OffAxisFluxConfig);
    #endif
  }

  //-----------------------------------------------------------------------

  const DUNEFluxSystPRISM *GetDUNEFluxSystPRISM(unsigned int i, bool applyPenalty) {

    static std::vector<const DUNEFluxSystPRISM *> cache_OffAxis;

    auto c = &cache_OffAxis;
    if (i >= c->size()) {
      c->resize(i + 1);
    }
    if (!c->at(i)) {
      c->at(i) = new DUNEFluxSystPRISM(i, applyPenalty);
    }
    return c->at(i);
  }

  //-----------------------------------------------------------------------

  DUNEFluxSystPRISMVector GetDUNEFluxSystsPRISM(unsigned int N, bool applyPenalty) {

    #ifdef USE_TH2JAGGED
    if (!fOffAxisFluxSystHelper) {
      fOffAxisFluxSystHelper = &PRISM::EffectiveFluxUncertaintyHelper::Get();
    }
    #endif

    DUNEFluxSystPRISMVector ret;
    for (unsigned int i = 0; i < N; i++) {
      ret.push_back(GetDUNEFluxSystPRISM(i, applyPenalty));
    }
    return ret;
  }

}
