#include "TCanvas.h"
#include "TMarker.h"
#include <cmath>

#include "CAFAna/Analysis/AnalysisVars.h"
#include "CAFAna/Analysis/common_fit_definitions.h"

#include "CAFAna/PRISM/PRISMExtrapolator.h"
#include "CAFAna/PRISM/PRISMUtils.h"
#include "CAFAna/PRISM/PredictionPRISM.h"

#include "CAFAna/Fit/FrequentistSurface.h"
#include "CAFAna/Fit/IFitter.h"
#include "CAFAna/Fit/MinuitFitter.h"
#include "CAFAna/Experiment/SingleSampleExperiment.h"

#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/make_ParameterSet.h"

using namespace ana;

std::map<std::string, PRISMStateBlob> States;
std::map<std::string, std::unique_ptr<PredictionInterp>> FarDetSelStates;

void PRISMPrediction(fhicl::ParameterSet const &pred) {

  std::string const &state_file = pred.get<std::string>("state_file");
  std::vector<std::string> const &output_file =
      pred.get<std::vector<std::string>>("output_file");
  std::string const &output_dir = pred.get<std::string>("output_dir");
  std::string const &varname =
      pred.get<std::string>("projection_name");
  bool isfhc = pred.get<bool>("isFHC", true);
  bool is_fake_spec_run = pred.get<bool>("is_fake_spec_run");

  double reg = pred.get<double>("reg_factor");
  std::array<double, 2> fit_range =
      pred.get<std::array<double, 2>>("fit_range");

  (void)GetListOfSysts();

  osc::IOscCalculatorAdjustable *calc =
      ConfigureCalc(pred.get<fhicl::ParameterSet>("Osc", {}));

  if (!States.count(state_file)) {
    TFile fs(state_file.c_str());
    std::cout << "Loading " << varname << " state from " << state_file
              << std::endl;
    States[state_file] = LoadPRISMState(fs, varname, !isfhc);
    FarDetSelStates[state_file] = PredictionInterp::LoadFrom(
        fs.GetDirectory((std::string("FarDetSel_") + varname +
                         std::string(!isfhc ? "_rhc" : "_fhc"))
                            .c_str()));
    std::cout << "Done!" << std::endl;
    fs.Close();
  }

  SystShifts shift = GetSystShifts(pred.get<fhicl::ParameterSet>("syst", {}));

  PRISMStateBlob &state = States[state_file];

  TFile f(output_file[0].c_str(),
          output_file.size() > 1 ? output_file[1].c_str() : "RECREATE");

  TDirectory *dir = &f;
  if (output_dir.size()) {
    dir = f.mkdir(output_dir.c_str());
  }
  dir->cd();

  int id = 0;
  PRISMExtrapolator fluxmatcher;

  fluxmatcher.InitializeEventRateMatcher(state.NDMatchInterp.get(),
                                         state.FDMatchInterp.get());
  fluxmatcher.SetStoreDebugMatches();

  if (!is_fake_spec_run) {
    fluxmatcher.SetTargetConditioning(reg,
                                      {{0},},
                                      fit_range[0], fit_range[1]);
  } else { // we are doing a fake special run
    fluxmatcher.SetNoRegAltHC();
    fluxmatcher.SetTargetConditioning(reg, {}, fit_range[0], fit_range[1]);
  }

  state.PRISM->SetFluxMatcher(&fluxmatcher);

  //state.PRISM->SetNCCorrection();
  //state.PRISM->SetWSBCorrection();
  //state.PRISM->SetNueCorrection();

  Spectrum PRISMPredEvRateMatchSpec = state.PRISM->PredictSyst(calc, shift);

  double pot = pot_fd * (1.0 / 3.5);

  TH1 *PRISMPredEvRateMatch_h = PRISMPredEvRateMatchSpec.ToTHX(pot);

  PRISMPredEvRateMatch_h->Scale(1, "width");
  PRISMPredEvRateMatch_h->SetTitle(";E_{#nu} (GeV);Pred. FD EvRate per 1 GeV");
  PRISMPredEvRateMatch_h->Write("PRISMPredEvRateMatch");

  for (auto &compspec : state.PRISM->PredictPRISMComponents(calc, shift)) {
    TH1 *comp = compspec.second.ToTHX(pot);
    comp->Scale(1, "width");
    comp->SetTitle(";E_{#nu} (GeV);Pred. FD Contribution per 1 GeV");
    comp->Write((std::string("PRISMPredEvRateMatch_") +
                 PredictionPRISM::GetComponentString(compspec.first))
                    .c_str());
  }

  fluxmatcher.Write(dir->mkdir("PRISMEventRateMatches"));

  TH1 *FarDet_h = state.FarDet->Predict(calc).ToTHX(pot);

  for (int bin_it = 0; bin_it < FarDet_h->GetXaxis()->GetNbins(); ++bin_it) {
    FarDet_h->SetBinError(bin_it + 1,
                          sqrt(FarDet_h->GetBinContent(bin_it + 1)));
  }
  FarDet_h->Scale(1, "width");

  FarDet_h->SetTitle(";E_{#nu} (GeV);FD EvRate");
  FarDet_h->Write("FarDet");

  int disp_pid = isfhc ? 14 : -14;
  TH1 *FarDetData_h =
      state.FarDetData->Oscillated(calc, disp_pid, disp_pid).ToTHX(pot);

  for (int bin_it = 0; bin_it < FarDetData_h->GetXaxis()->GetNbins();
       ++bin_it) {
    FarDetData_h->SetBinError(bin_it + 1,
                              sqrt(FarDetData_h->GetBinContent(bin_it + 1)));
  }
  FarDetData_h->Scale(1, "width");

  FarDetData_h->SetTitle(";E_{#nu} (GeV);FD EvRate");
  FarDetData_h->Write("FarDetData");

  if (FarDetSelStates[state_file]) {
    TH1 *FarDetSelNC_h =
        FarDetSelStates[state_file]
            ->PredictComponent(calc, ana::Flavors::kAll, ana::Current::kNC,
                               ana::Sign::kBoth)
            .ToTHX(pot);
    TH1 *FarDetSelWSB_h =
        FarDetSelStates[state_file]
            ->PredictComponent(calc, ana::Flavors::kAllNuMu, ana::Current::kCC,
                               ana::Sign::kAntiNu)
            .ToTHX(pot);
    FarDetSelNC_h->Scale(1, "width");
    FarDetSelNC_h->Write("FarDetSel_NC");
    FarDetSelWSB_h->Scale(1, "width");
    FarDetSelWSB_h->Write("FarDetSel_WSB");
  }

  TH1 *FarDet_unosc_h = state.FarDet->PredictUnoscillated().ToTHX(pot);
  FarDet_unosc_h->Scale(1, "width");

  FarDet_unosc_h->SetTitle(";E_{#nu} (GeV);FD EvRate");
  FarDet_unosc_h->Write("FarDet_unosc");

  TH1 *NearDet_h = state.NDMatchInterp->Predict(calc).ToTHX(pot);
  NearDet_h->SetTitle(";E_{#nu} (GeV);OffAxis;FD EvRate");
  NearDet_h->Write("NearDet");

  //******************************************
  // Attempt a ChiSq fit using CAFAna method *
  //******************************************

  std::cout << std::endl << "True Osc Parameters: " << std::endl;
  std::cout << "dMsq32 = " << calc->GetDmsq32() << std::endl;
  std::cout << "Theta23 = " << calc->GetTh23() << std::endl;
  std::cout << "Theta13 = " << calc->GetTh13() << std::endl;
  std::cout << "dCP = " << calc->GetdCP() << std::endl;

  std::vector<const ISyst*> allsysts = shift.ActiveSysts();
  std::vector<const IFitVar*> profOscParams = {&kFitDeltaInPiUnits}; //&kFitSinSq2Theta13
  std::vector<double> seed_values;
  for (const IFitVar *v : profOscParams) {
    seed_values.push_back(v->GetValue(calc));
  }

  // PRISM prediction comparison with FD observed 'data'
  const SingleSampleExperiment exptData = SingleSampleExperiment(state.PRISM.get(),
                                                           state.FarDetData->Oscillated(calc, 14, 14).FakeData(pot));

  const IFitVar* Dmsq32 = &kFitDmSq32Scaled;
  const IFitVar* ssTh23 = &kFitSinSqTheta23;
  std::vector<double> Dmsq32_scan;
  std::vector<double> ssTh23_scan;

  const int NXSteps(41);
  const int NYSteps(61);

  std::unique_ptr<TH2D> scan_hist = std::make_unique<TH2D>(
    "dchi2_2DScan", "dchi2", NXSteps, 0.42, 0.62, NYSteps, 2.3, 2.6);
  for (int i = 0; i < NXSteps; i++)
    ssTh23_scan.emplace_back(scan_hist->GetXaxis()->GetBinCenter(i + 1));
  for (int i = 0; i < NYSteps; i++)
    Dmsq32_scan.emplace_back(scan_hist->GetYaxis()->GetBinCenter(i + 1));

  for (const auto &x : ssTh23_scan) {
    for (const auto &y : Dmsq32_scan) {
      ssTh23->SetValue(calc, x);
      Dmsq32->SetValue(calc, y);
      // set profiled params back to see value each iteration
      for (unsigned int i = 0; i < seed_values.size(); i++) {
        profOscParams[i]->SetValue(calc, seed_values[i]);
      }
      MinuitFitter fitter(&exptData, profOscParams, allsysts, MinuitFitter::kNormal);
      fitter.SetFitOpts(MinuitFitter::kNormal);
      SystShifts bestSysts;
      const SeedList &seedPts = SeedList();
      double chi = fitter.Fit(calc, bestSysts, seedPts, {}, MinuitFitter::kVerbose);
      // fill 2D scan with chisq value
      scan_hist->Fill(x, y, chi);
    }
  }
  // Get minimum ChiSq value (LL value)
  double minchi = 1e10;
  int minx = scan_hist->GetNbinsX()/2;
  int miny = scan_hist->GetNbinsY()/2;
  for (int x = 1; x <= scan_hist->GetNbinsX(); ++x) {
    for (int y = 1; y <= scan_hist->GetNbinsY(); ++y) {
      const double chi = scan_hist->GetBinContent(x, y);
      if (chi < minchi) {
        minchi = chi;
        minx = x;
        miny = y;
      }
    }
  }
  ssTh23->SetValue(calc, scan_hist->GetXaxis()->GetBinCenter(minx));
  Dmsq32->SetValue(calc, scan_hist->GetYaxis()->GetBinCenter(miny));
  double BestParamFitX = ssTh23->GetValue(calc);
  double BestParamFitY = Dmsq32->GetValue(calc);
  std::cout << "Best fit parameter values: " << BestParamFitX << ", " << BestParamFitY << std::endl;
  double BestLL = minchi;

  for (int x = 0; x < scan_hist->GetNbinsX(); x++) {
    for (int y = 0; y < scan_hist->GetNbinsY(); y++) {
      scan_hist->SetBinContent(x + 1, y + 1, scan_hist->GetBinContent(x + 1, y + 1) - BestLL);
    }
  }

  std::cout << std::endl << "Best-Fit Osc Parameters: " << std::endl;
  std::cout << "dMsq32 = " << calc->GetDmsq32() << std::endl;
  std::cout << "Theta23 = " << calc->GetTh23() << std::endl;
  std::cout << "Theta13 = " << calc->GetTh13() << std::endl;
  std::cout << "dCP = " << calc->GetdCP() << std::endl;

  scan_hist->Write();
  std::unique_ptr<TCanvas> c = std::make_unique<TCanvas>("c", "Contours", 800, 600);
  scan_hist->Draw("COLZ");

  c->Update();
  c->Write();

  //******************************************

  f.Write();
}

int main(int argc, char const *argv[]) {
  gROOT->SetMustClean(false);

  if (argc != 2) {
    std::cout << "[ERROR]: Expected to be passed the location of a single "
                 "configuration FHiCL file."
              << std::endl;
    return 1;
  }

  fhicl::ParameterSet const &ps = fhicl::make_ParameterSet(argv[1]);

  for (fhicl::ParameterSet const &pred :
       ps.get<std::vector<fhicl::ParameterSet>>("predictions")) {
    PRISMPrediction(pred);
  }
}
