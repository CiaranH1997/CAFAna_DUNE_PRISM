//
// Created by hasnipl on 29/06/2020.
//

#include "CAFAna/Analysis/AnalysisVars.h"
#include "CAFAna/Analysis/Exposures.h"
#include "CAFAna/Analysis/common_fit_definitions.h"

#include "CAFAna/Core/HistCache.h"
#include "CAFAna/Core/SystShifts.h"

#include "CAFAna/Experiment/ReactorExperiment.h"
#include "CAFAna/Experiment/SingleSampleExperiment.h"

#include "CAFAna/PRISM/PRISMExtrapolator.h"
#include "CAFAna/PRISM/PRISMUtils.h"
#include "CAFAna/PRISM/PredictionPRISM.h"
#include "CAFAna/PRISM/SimpleChi2Experiment.h"

#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/make_ParameterSet.h"

#include "TArrow.h"
#include "TCanvas.h"
#include "TGraph2D.h"
#include "TMarker.h"
#include <map>
#include <utility>

using namespace ana;
using namespace PRISM;

std::map<std::string, PRISMStateBlob> States;

void PRISMScan(fhicl::ParameterSet const &scan) {

  std::string const &state_file = scan.get<std::string>("state_file");

  std::vector<std::string> const &output_file =
    scan.get<std::vector<std::string>>("output_file");
  std::string const &output_dir = scan.get<std::string>("output_dir", "");

  std::string const &varname =
  scan.get<std::string>("projection_name", "EProxy");

  // default to 1 year
  double POT = scan.get<double>("POT_years", 1) * POT120;
  double POT_FD = POT * pot_fd_FVMassFactor;
  std::cout << "POT : " << POT << ", " << POT_FD << std::endl;

  bool use_PRISM = scan.get<bool>("use_PRISM", true);

  auto GOFps = scan.get<fhicl::ParameterSet>("GOF", {});
  bool use_PRISM_ND_stats = GOFps.get<bool>("use_PRISM_ND_stats", true);

  bool fit_nuisance = GOFps.get<bool>("fit_nuisance", false);
  bool poisson_throw = GOFps.get<bool>("poisson_throw", false);

  std::vector<const IFitVar *> free_oscpars = GetOscVars(
    GOFps.get<std::vector<std::string>>("free_osc_params", {}));
  std::vector<ISyst const *> freesysts = ana::GetListOfSysts(
    GOFps.get<std::vector<std::string>>("free_syst_params", {}));

  for (auto &s : freesysts) {
    std::cout << "\t" << s->ShortName() << " free." << std::endl;
  }

  //SystShifts shift;
  double sigma(1);
  std::map<const ISyst *, double> shift_map;
  for (auto &sh : freesysts) {
    shift_map[sh] = sigma;
    //shift.SetShift(sh, sigma);
  }
  SystShifts shift(shift_map);
  for (auto &s : freesysts) {
    std::cout << "\t" << shift.GetShift(s) << " shift." << std::endl;
  }

  auto PRISMps = scan.get<fhicl::ParameterSet>("PRISM", {});

  auto RunPlans = PRISMps.get<fhicl::ParameterSet>("RunPlans", {});
  bool PRISM_write_debug = PRISMps.get<bool>("write_debug", false);

  RunPlan run_plan_nu, run_plan_nub;

  if (RunPlans.has_key("numode")) {
    run_plan_nu =
      make_RunPlan(RunPlans.get<fhicl::ParameterSet>("numode"), POT);
  }
  if (RunPlans.has_key("nubmode")) {
    run_plan_nub =
      make_RunPlan(RunPlans.get<fhicl::ParameterSet>("nubmode"), POT);
  }

  auto params = scan.get<std::vector<fhicl::ParameterSet>>("scan_params");
  std::vector<std::string> param_names;
  std::vector<std::vector<double>> steps;
  for (auto &sp : params) {
    param_names.push_back(sp.get<std::string>("name", ""));
    std::cout << "param name: " << sp.get<std::string>("name", "") << std::endl;
    auto vec_element = sp.get<std::vector<double>>("scan_steps");
    for (auto &v : vec_element)
      std::cout << "scan def: " << v << std::endl;
    steps.push_back(vec_element);
  }

  // Get x and y parameter names
  int nparams = param_names.size();
  std::string xparam_name = param_names[0];
  std::string yparam_name;
  if (param_names.size() > 1) yparam_name = param_names[1];
  bool dmsq32_scan(false), ssth23_scan(false);

  // vector<const IFitVar *> of parameters to scan
  std::vector<const IFitVar *> scan_vars = GetOscVars(param_names);

  // Don't let scan params move freely
  ScrubOscVars(free_oscpars, param_names);

  // set up seed points for multiple fits.
  std::map<const IFitVar *, std::vector<double>> oscSeeds;
  if (std::find(free_oscpars.begin(), free_oscpars.end(),
    &kFitDeltaInPiUnits) != free_oscpars.end()) {
    oscSeeds[&kFitDeltaInPiUnits] = {-1, -0.5, 0, 0.5};
  }
  if (std::find(free_oscpars.begin(), free_oscpars.end(), &kFitSinSqTheta23) !=
    free_oscpars.end()) {
    oscSeeds[&kFitSinSqTheta23] = {0.4, 0.6};
  }

  bool dmsq32_scale = false;
  if (std::find(scan_vars.begin(), scan_vars.end(), &kFitDmSq32NHScaled) !=
      scan_vars.end())
    dmsq32_scale = true;

  osc::IOscCalculatorAdjustable *calc =
    ConfigureCalc(scan.get<fhicl::ParameterSet>("true_osc", {}));

  // Lazy load the state file
  if (!States.count(state_file)) {
    TFile fs(state_file.c_str());
    std::cout << "Loading " << varname << " state from " << state_file
      << std::endl;
    States[state_file] = LoadPRISMState(fs, varname);
    std::cout << "Done!" << std::endl;
    fs.Close();
  }

  PRISMStateBlob &state = States[state_file];

  TFile f(output_file[0].c_str(),
    output_file.size() > 1 ? output_file[1].c_str() : "RECREATE");

  TDirectory *dir = &f;
  if (output_dir.size()) {
    dir = f.mkdir(output_dir.c_str());
  }
  dir->cd();

  PRISMExtrapolator fluxmatcher;
  if (use_PRISM) {
    fluxmatcher.Initialize({
      {"ND_nu", state.MatchPredInterps[kND_nu].get()},
      {"FD_nu", state.MatchPredInterps[kFD_nu_nonswap].get()},
      {"ND_nub", state.MatchPredInterps[kND_nub].get()},
      {"FD_nub", state.MatchPredInterps[kFD_nub_nonswap].get()},
    });

    for (auto const &channel_conditioning :
      PRISMps.get<std::vector<fhicl::ParameterSet>>("match_conditioning")) {

      auto ch = GetMatchChan(channel_conditioning.get<fhicl::ParameterSet>("chan"));

      double chan_reg = channel_conditioning.get<double>("reg_factor", 1E-16);
      std::array<double, 2> chan_energy_range =
        channel_conditioning.get<std::array<double, 2>>("energy_range",
                                                       {0, 4});
      bool chan_is_fake_spec_run =
        channel_conditioning.get<bool>("is_fake_spec_run", false);
      if (chan_is_fake_spec_run) fluxmatcher.SetNoRegAltHC(true);

      fluxmatcher.SetTargetConditioning(ch, chan_reg, (chan_is_fake_spec_run ?
        std::vector<double>{{0},}
        : std::vector<double>{}),
        chan_energy_range);
    }
    if (PRISM_write_debug) {
      fluxmatcher.SetStoreDebugMatches();
    }
    state.PRISM->SetFluxMatcher(&fluxmatcher);
  }

  if (run_plan_nu.GetPlanPOT() > 0) {
    state.PRISM->SetNDRunPlan(run_plan_nu, BeamMode::kNuMode);
  }

  if (run_plan_nub.GetPlanPOT() > 0) {
    state.PRISM->SetNDRunPlan(run_plan_nub, BeamMode::kNuBarMode);
  }

  std::map<std::string, MatchChan> Channels;
  if (scan.is_key_to_sequence("samples")) {
    for (auto const &fs : scan.get<std::vector<fhicl::ParameterSet>>("samples")) {
      auto ch = GetMatchChan(fs);
      Channels[GetMatchChanShortName(ch)] = ch;
    }
  } else {
    auto ch = GetMatchChan(scan.get<fhicl::ParameterSet>("samples"));
    Channels[GetMatchChanShortName(ch)] = ch;
  }

  std::vector<Spectrum> DataSpectra;
  DataSpectra.reserve(Channels.size());

  std::vector<const IChiSqExperiment*> Expts;
  Expts.reserve(Channels.size());

  for (auto const ch : Channels) {

    int osc_from = FluxSpeciesPDG(ch.second.from.chan);
    int osc_to = FluxSpeciesPDG(ch.second.to.chan);
    size_t NDConfig_enum = GetConfigFromNuChan(ch.second.from, true);
    size_t FDConfig_enum = GetConfigFromNuChan(ch.second.to, false);
    size_t FDfdConfig_enum = GetFDConfigFromNuChan(ch.second.to);

    if ((NDConfig_enum == kND_nu) && !run_plan_nu.GetPlanPOT()) {
      std::cout << "[ERROR]: Have ND nu channel, but no numode run plan."
                << std::endl;
      abort();
    }
    if ((NDConfig_enum == kND_nub) && !run_plan_nub.GetPlanPOT()) {
      std::cout << "[ERROR]: Have ND nubar channel, but no numode run plan."
                << std::endl;
      abort();
    }

    TDirectory *chan_dir =
      dir->mkdir(DescribeFDConfig(FDfdConfig_enum).c_str());
    chan_dir->cd();

    DataSpectra.push_back(state.FarDetData_nonswap[FDfdConfig_enum]->Oscillated(
      calc, osc_from, osc_to
    ));

    TH1 *Data = DataSpectra.back().ToTH1(POT_FD);
    Data->Scale(1, "width");
    chan_dir->WriteTObject(Data, "Data_Total");
    Data->SetDirectory(nullptr);

    if (use_PRISM) {
      auto PRISMComponents =
        state.PRISM->PredictPRISMComponents(calc, shift, ch.second);
      TH1 *PRISMPred =
        PRISMComponents.at(PredictionPRISM::kPRISMPred).ToTH1(POT_FD);
      PRISMPred->Scale(1, "width");
      chan_dir->WriteTObject(PRISMPred, "PRISMPred");
      PRISMPred->SetDirectory(nullptr);
    }

    //Expts.emplace_back(new SingleSampleExperiment(state.PRISM.get(),
    //                                              DataSpectra.back().FakeData(POT_FD)));
    Expts.emplace_back(new PRISMChi2Experiment(state.PRISM.get(),
      DataSpectra.back().FakeData(POT_FD), false, POT_FD, ch.second, {0, 8}));
    Expts.emplace_back(new ReactorExperiment(0.088, 0.003));

    const MultiExperiment CombExpts = MultiExperiment(Expts);

    std::vector<double> x_scan;
    std::vector<double> y_scan;
    std::unique_ptr <TH1D> scan_hist_1D;
    std::unique_ptr <TH2D> scan_hist_2D;

    if (nparams == 1) {

      const int NXSteps = steps.at(0).at(0);
      const double x_low_bound = dmsq32_scale ? steps.at(0).at(1) * 1000 : steps.at(0).at(1);
      const double x_high_bound = dmsq32_scale ? steps.at(0).at(2) * 1000 : steps.at(0).at(2);
      std::cout << NXSteps << ", " << x_low_bound << ", " << x_high_bound << std::endl;
      scan_hist_1D = std::make_unique<TH1D>("dchi2_1D", "dchi2_1D",
        NXSteps, x_low_bound, x_high_bound);

      // place the scan points in vector
      for (int i = 0; i < scan_hist_1D->GetNbinsX(); i++) {
        x_scan.emplace_back(scan_hist_1D->GetXaxis()->GetBinCenter(i + 1));
      }
      // else its a 2D hist
    } else {

      const int NXSteps = steps.at(0).at(0);
      const double x_low_bound = steps.at(0).at(1);
      const double x_high_bound = steps.at(0).at(2);
      const int NYSteps = steps.at(1).at(0);
      const double y_low_bound = dmsq32_scale ? steps.at(1).at(1) * 1000 : steps.at(1).at(1);
      const double y_high_bound = dmsq32_scale ? steps.at(1).at(2) * 1000 : steps.at(1).at(2);
      std::cout << NXSteps << ", " << x_low_bound << ", " << x_high_bound << std::endl;
      std::cout << NYSteps << ", " << y_low_bound << ", " << y_high_bound << std::endl;
      scan_hist_2D = std::make_unique<TH2D>("dchi2_2D", "dchi2_2D",
        NXSteps, x_low_bound, x_high_bound,
        NYSteps, y_low_bound, y_high_bound);

      // place scan points in vector
      for (int i = 0; i < scan_hist_2D->GetNbinsX(); i++)
        x_scan.emplace_back(scan_hist_2D->GetXaxis()->GetBinCenter(i + 1));
      for (int i = 0; i < scan_hist_2D->GetNbinsY(); i++)
        y_scan.emplace_back(scan_hist_2D->GetYaxis()->GetBinCenter(i + 1));
    }

    const IFitVar *ssTh23 = &kFitSinSqTheta23;
    const IFitVar *dmsq32 = &kFitDmSq32Scaled;

    if (std::find(scan_vars.begin(), scan_vars.end(), &kFitSinSqTheta23)
        != scan_vars.end()) {
      std::cout << "found ssth23" << std::endl;
      ssth23_scan = true;
    }
    if (std::find(scan_vars.begin(), scan_vars.end(), &kFitDmSq32NHScaled)
        != scan_vars.end()) {
      std::cout << "found dmsq32" << std::endl;
      dmsq32_scan = true;
    }

    //-----------------------
    // Do the minimization 1D
    //-----------------------
    if (nparams == 1) {
      for (const auto &x : x_scan) {
        if (ssth23_scan) {
          ssTh23->SetValue(calc, x);
        }
        else if (dmsq32_scan) {
          dmsq32->SetValue(calc, x);
        }

        std::cout << "dMsq32 = " << calc->GetDmsq32() << std::endl;
        std::cout << "ssth23 = " << calc->GetTh23() << std::endl;

        MinuitFitter fitter(&CombExpts, free_oscpars, freesysts, MinuitFitter::kNormal);
        //SystShifts bestSysts;
        //const SeedList &seedPts = SeedList(); //oscSeeds
        double chi = fitter.Fit(calc, junkShifts, oscSeeds, {}, MinuitFitter::kVerbose);
        // fill hist
        scan_hist_1D->Fill(x, chi);
      }
      // Get minimum ChiSq value (LL value)
      double minchi = 1e10;
      int minx, miny;
      if (nparams == 1) {
        //double minchi = 1e10;
        minx = scan_hist_1D->GetNbinsX() / 2;
        for (int x = 1; x <= scan_hist_1D->GetNbinsX(); ++x) {
          const double chi = scan_hist_1D->GetBinContent(x);
          if (chi < minchi) {
            minchi = chi;
            minx = x;
          }
        }
      }
      if (ssth23_scan) ssTh23->SetValue(calc, scan_hist_1D->GetXaxis()->GetBinCenter(minx));
      else if (dmsq32_scan) dmsq32->SetValue(calc, scan_hist_1D->GetXaxis()->GetBinCenter(minx));
      std::cout << "Bestfit parameter values: " <<
        scan_hist_1D->GetXaxis()->GetBinCenter(minx) << std::endl;
      double BestLL = minchi;
      for (int x = 0; x < scan_hist_1D->GetNbinsX(); x++) {
          scan_hist_1D->SetBinContent(x + 1, scan_hist_1D->GetBinContent(x + 1) - BestLL);
      }
      chan_dir->WriteTObject(scan_hist_1D.release(), "dChi2Scan");
    }
    //-----------------------
    // Do the minimization 2D
    //-----------------------
    else if (param_names.size() > 1) {
      for (const auto &x : x_scan) {
        for (const auto &y : y_scan) {
          std::cout << "x = " << x << ", y = " << y << std::endl;
          ssTh23->SetValue(calc, x);
          dmsq32->SetValue(calc, y);

          MinuitFitter fitter(&CombExpts, free_oscpars, freesysts, MinuitFitter::kNormal);
          SystShifts bestSysts;
          //const SeedList &seedPts = SeedList(); //oscSeeds
          double chi = fitter.Fit(calc, bestSysts, oscSeeds, {}, MinuitFitter::kVerbose);
          // fill hist
          scan_hist_2D->Fill(x, y, chi);
        }
      }
      // Get minimum ChiSq value (LL value)
      double minchi = 1e10;
      int minx = scan_hist_2D->GetNbinsX() / 2;
      int miny = scan_hist_2D->GetNbinsY() / 2;
      for (int x = 1; x <= scan_hist_2D->GetNbinsX(); ++x) {
        for (int y = 1; y <= scan_hist_2D->GetNbinsY(); ++y) {
          const double chi = scan_hist_2D->GetBinContent(x, y);
          if (chi < minchi) {
            minchi = chi;
            minx = x;
            miny = y;
          }
        }
      }

      ssTh23->SetValue(calc, scan_hist_2D->GetXaxis()->GetBinCenter(minx));
      dmsq32->SetValue(calc, scan_hist_2D->GetYaxis()->GetBinCenter(miny));
      std::cout << "Bestfit parameter values: " << ssTh23->GetValue(calc) <<
                ", " << dmsq32->GetValue(calc) << std::endl;
      double BestLL = minchi;
      for (int x = 0; x < scan_hist_2D->GetNbinsX(); x++) {
        for (int y = 0; y < scan_hist_2D->GetNbinsY(); y++) {
          scan_hist_2D->SetBinContent(x + 1, y + 1,
            scan_hist_2D->GetBinContent(x + 1, y + 1) - BestLL);
        }
      }
      chan_dir->WriteTObject(scan_hist_2D.release(), "dChi2Scan");
    }
  }
  // Write total output file
  f.Write();
}

//-------------------------------------------

// main function

int main(int argc, char const *argv[]) {
  // Make sure systs are applied to ND distributions which are per 1 POT.
  setenv("CAFANA_PRED_MINMCSTATS", "0", 1);
  gROOT->SetMustClean(false);

  if (argc != 2) {
    std::cout << "[ERROR]: Expected to be passed the location of a single "
                 "configuration FHiCL file."
              << std::endl;
    return 1;
  }

  fhicl::ParameterSet const &ps = fhicl::make_ParameterSet(argv[1]);

  for (fhicl::ParameterSet const &pred :
       ps.get<std::vector<fhicl::ParameterSet>>("scans")) {
    PRISMScan(pred);
  }
}
