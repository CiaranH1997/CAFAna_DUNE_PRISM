//
// Created by C. Hasnip on 25/06/2020.
//

#pragma once

#include "CAFAna/Core/ISyst.h"
#ifdef USE_TH2JAGGED
#include "CAFAna/PRISM/EffectiveFluxUncertaintyHelper.h"
#endif
#include "TString.h"

namespace ana {
  class DUNEFluxSystPRISM : public ISyst {
  public:
    virtual ~DUNEFluxSystPRISM();

    virtual void Shift(double sigma, Restorer &restore, caf::StandardRecord *sr,
                       double &weight) const override;

  protected:
    friend const DUNEFluxSystPRISM *GetDUNEFluxSystPRISM(unsigned int, bool);

    DUNEFluxSystPRISM(int i, bool applyPenalty)
      : ISyst(TString::Format("flux_%s_%i", "OffAxis", i).Data(),
              TString::Format("Flux #%i (%s)", i, "OffAxis").Data(),
              applyPenalty),
        fIdx(i) {}

    int fIdx;
    //fOffAxisFluxSystHelper()
    //PRISM::EffectiveFluxUncertaintyHelper const *fOffAxisFluxSystHelper =
    //  &PRISM::EffectiveFluxUncertaintyHelper::Get();
    //mutable TH1 *fScale[2][2][2][2]; // ND/FD, numu/nue, bar, FHC/RHC

  };

  const DUNEFluxSystPRISM *GetDUNEFluxSystPRISM(unsigned int i, bool applyPenalty = true);

  // Because vector<T*> won't automatically convert to vector<U*> even when U
  // inherits from V.
  struct DUNEFluxSystPRISMVector : public std::vector<const DUNEFluxSystPRISM *> {
    operator std::vector<const ISyst *>() {
      return std::vector<const ISyst *>(begin(), end());
    }
  };

  DUNEFluxSystPRISMVector GetDUNEFluxSystsPRISM(unsigned int N, bool applyPenalty = true);
}