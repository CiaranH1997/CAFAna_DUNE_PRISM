#include "TCanvas.h"
#include "TMarker.h"
#include <cmath>

#include "CAFAna/Analysis/AnalysisVars.h"
#include "CAFAna/Analysis/common_fit_definitions.h"

#include "CAFAna/PRISM/PRISMExtrapolator.h"
#include "CAFAna/PRISM/PRISMUtils.h"
#include "CAFAna/PRISM/PredictionPRISM.h"

#include "CAFAna/Fit/FrequentistSurface.h"
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

  //******************************************
  // Attempt a ChiSq fit using CAFAna method *
  //******************************************

  std::cout << std::endl << "True Osc Parameters: " << std::endl;
  std::cout << "dMsq32 = " << calc->GetDmsq32() << std::endl;
  std::cout << "Theta23 = " << calc->GetTh23() << std::endl;

  TMarker *true_point = new TMarker(pow(sin(calc->GetTh23()), 2), calc->GetDmsq32()*1000, kFullCircle);

  // FD MC prediction comparison with FD observed 'data'
  SingleSampleExperiment exptMC = SingleSampleExperiment(state.FarDet.get(),
                                                         state.FarDetData->Oscillated(calc, 14, 14).FakeData(pot));

  // PRISM prediction comparison with FD observed 'data'
  SingleSampleExperiment exptData = SingleSampleExperiment(state.PRISM.get(),
                                                           state.FarDetData->Oscillated(calc, 14, 14).FakeData(pot));

  FrequentistSurface surfMC = FrequentistSurface(&exptMC, calc,
                                                 &kFitSinSqTheta23, 30, 0.4, 0.6,
                                                 &kFitDmSq32Scaled, 30, 2.2, 2.6);

  FrequentistSurface surfData = FrequentistSurface(&exptData, calc,
                                                 &kFitSinSqTheta23, 30, 0.4, 0.6,
                                                 &kFitDmSq32Scaled, 30, 2.2, 2.6);

  TH2 *crit90sigMC = Gaussian90Percent2D(surfMC);
  TH2 *crit90sigData = Gaussian90Percent2D(surfData);

  TCanvas *c = new TCanvas("c", "Contours", 800, 600);
  surfMC.DrawContour(crit90sigMC, kSolid, kBlue);
  //surfMC.DrawBestFit(kBlue);
  surfData.DrawContour(crit90sigData, kSolid, kRed);
  //surfData.DrawBestFit(kRed);
  true_point->SetMarkerSize(1.5);
  true_point->SetMarkerColor(kBlack);
  true_point->Draw();

  c->Update();
  c->Write();

  //******************************************

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
