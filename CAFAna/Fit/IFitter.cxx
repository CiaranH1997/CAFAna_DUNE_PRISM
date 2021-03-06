#include "CAFAna/Fit/IFitter.h"

#include "CAFAna/Core/IFitVar.h"
#include "CAFAna/Core/ISyst.h"
#include "CAFAna/Core/Utilities.h"
#include "CAFAna/Experiment/IChiSqExperiment.h"

#include "OscLib/func/IOscCalculator.h"

#include "Utilities/func/MathUtil.h"

#include <algorithm>

#include <cassert>

namespace ana
{

  SystShifts junkShifts = SystShifts();

  //----------------------------------------------------------------------
  IFitter::IFitter(const std::vector<const IFitVar*>& vars,
                   const std::vector<const ISyst*>& systs)
    : fVars(vars), fSysts(systs)
  {
  }

  //----------------------------------------------------------------------
  IFitter::IFitter(const IFitter &other)
      : fVars(other.fVars),
        fSysts(other.fSysts),
        fShifts(std::make_unique<SystShifts>(*other.fShifts))
  {}

  //----------------------------------------------------------------------
  IFitter::~IFitter()
  {
  }

  //----------------------------------------------------------------------
  std::vector<IFitter::SeedPt>
  IFitter::ExpandSeeds(const SeedList& seedPts,
                       const std::vector<SystShifts>& systSeedPts) const
  {
    std::vector<SeedPt> ret;
    for(Seed seed: seedPts.GetSeeds()) ret.emplace_back(seed, SystShifts::Nominal());

    // Now duplicate as many times as required for the syst seeds
    if(!systSeedPts.empty()){
      std::vector<SeedPt> newret;
      for(const SystShifts& s: systSeedPts){
        for(SeedPt pt: ret){
          pt.shift = s;
          newret.push_back(pt);
        }
      }
      ret = newret;
    }

    return ret;
  }

  //----------------------------------------------------------------------
  double IFitter::Fit(osc::IOscCalculatorAdjustable* seed,
                      SystShifts& bestSysts,
                      const SeedList& seedPts,
                      const std::vector<SystShifts> &systSeedPts,
                      Verbosity verb) const
  {
    ValidateSeeds(seed, seedPts, systSeedPts);

    if (verb == kVerbose)
    {
      std::cout << "Finding best fit for";
      for (const auto *v: fVars) std::cout << " " << v->ShortName();
      for (const auto *s: fSysts) std::cout << " " << s->ShortName();
//      if (fSupportsDerivatives) std::cout << " using analytic derivatives";
      std::cout << "..." << std::endl;
    }

    // Do all the actual work. This wrapper function is just so we can have
    // better control of printing.
    const double chi = FitHelper(seed, bestSysts, seedPts, systSeedPts, verb);

    if (verb == kVerbose)
    {
      std::cout << "Best fit";
      for (const auto *v: fVars)
      {
        std::cout << ", " << v->ShortName() << " = " << v->GetValue(seed);
      }
      for (const auto *s: fSysts)
      {
        std::cout << ", " << s->ShortName() << " = " << bestSysts.GetShift(s);
      }
      std::cout << ", LL = " << chi << std::endl;

//      std::cout << "  found with " << fNEval << " evaluations of the likelihood";
//      if (fNEvalFiniteDiff > 0)
//        std::cout << " (" << fNEvalFiniteDiff << " to evaluate gradients numerically)";
//      if (fNEvalGrad > 0)
//        std::cout << " and " << fNEvalGrad << " evaluations of the gradient";
//      std::cout << std::endl;
    }

    return chi;
  }

  //----------------------------------------------------------------------
  double IFitter::FitHelper(osc::IOscCalculatorAdjustable* initseed,
                            SystShifts& bestSysts,
                            const SeedList& seedPts,
                            const std::vector<SystShifts>& systSeedPts,
                            Verbosity verb) const
  {
    // if user passed a derived kind of SystShifts, this preserves it
    fShifts = bestSysts.Copy();
    fShifts->ResetToNominal();

    const std::vector<SeedPt> pts = ExpandSeeds(seedPts, systSeedPts);

    double minchi = 1e10;
    std::vector<double> bestFitPars, bestSystPars;

    for (const SeedPt &pt: pts)
    {
      osc::IOscCalculatorAdjustable* seed = nullptr;
      if(initseed) seed = initseed->Copy();

      pt.fitvars.ResetCalc(seed);

      // be sure to keep any derived class stuff around
      auto shift = fShifts->Copy();
      *shift = pt.shift;
      const double chi = FitHelperSeeded(seed, *shift, verb);
      // todo: this test doesn't make much sense for Bayes LL fitting...
      if (chi < minchi)
      {
        minchi = chi;
        // Store the best fit values of all the parameters we know are being
        // varied.
        bestFitPars.clear();
        bestSystPars.clear();
        for (const auto *v: fVars)
          bestFitPars.push_back(v->GetValue(seed));
        for (const auto *s: fSysts)
          bestSystPars.push_back(shift->GetShift(s));
      }

      delete seed;
    } // end for pt

    // Stuff the results of the actual best fit back into the seeds
    for (unsigned int i = 0; i < fVars.size(); ++i)
      fVars[i]->SetValue(initseed, bestFitPars[i]);
    for (unsigned int i = 0; i < fSysts.size(); ++i)
      bestSysts.SetShift(fSysts[i], bestSystPars[i]);

    return minchi;
  }

  //----------------------------------------------------------------------
  void IFitter::ValidateSeeds(osc::IOscCalculatorAdjustable* seed,
                              const SeedList& seedPts,
                              const std::vector<SystShifts>& systSeedPts) const
  {
    if(!seed && !fVars.empty()){
      std::cout << "ERROR: MinuitFitter::Fit() trying to fit oscillation parameters without an oscillation calculator" << std::endl;
      abort();
    }

    for(const IFitVar* v: seedPts.ActiveFitVars()){
      if (std::find(fVars.begin(), fVars.end(), v) == fVars.end()){
        std::cout << "ERROR MinuitFitter::Fit() trying to seed '"
                  << v->ShortName()
                  << "' which is not part of the fit." << std::endl;
        abort();
      }
    }

    for(const SystShifts& pt: systSeedPts){
      for(const ISyst* s: pt.ActiveSysts()){
        if(std::find(fSysts.begin(), fSysts.end(), s) == fSysts.end()){
          std::cout << "ERROR MinuitFitter::Fit() trying to seed '"
                    << s->ShortName()
                    << "' which is not part of the fit." << std::endl;
          abort();
        }
      }
    }
  }


} // namespace ana
