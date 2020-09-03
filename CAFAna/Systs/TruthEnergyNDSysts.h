#pragma once

#include "CAFAna/Core/ISyst.h"
#include "StandardRecord/StandardRecord.h"
#include "CAFAna/Core/Utilities.h"

#include "TFile.h"
#include "TH1.h"
#include "TRandom3.h"

#include <cassert>

// Implement ND energy scale uncertainties on truth information
// for use in eRecProxy
// LepE, eP, ePip, ePim, ePi0, eother

namespace ana
{

  // Global true energy scale syst for eRecProxy
  // Don't shift muon energies with this

  class TruthEnergyScaleND : public ISyst {
  public:
    TruthEnergyScaleND() : ISyst("TruthEnergyScaleND", "Global Truth Energy Scale ND Syst") {}
    void Shift(double sigma,
               Restorer& restore, 
               caf::StandardRecord* sr,
               double& weight) const override {
  
      restore.Add(sr->dune.eRecProxy,
                  sr->dune.LepE);

      const double scale = 0.02 * sigma;
      if (!sr->dune.isFD) { // in the ND
        if (sr->dune.isCC && abs(sr->dune.nuPDG) == 14) { // take away muon energy
          sr->dune.eRecProxy += (sr->dune.eRecProxy - sr->dune.LepE) * scale;
        }
        else if (sr->dune.isCC && abs(sr->dune.nuPDG) == 12) { // fine to include electron energy
          sr->dune.eRecProxy += sr->dune.eRecProxy * scale; 
        }
      }
    }
  };

  extern const TruthEnergyScaleND kTruthEnergyScaleND;

  // Total energy scale syst varying with sqrt of the energy
  class TruthEnergySqrtND: public ISyst {
  public:
    TruthEnergySqrtND() : ISyst("TruthEnergySqrtND", "Sqrt Total Energy Scale ND Syst") {}
    void Shift(double sigma,
               Restorer& restore,
               caf::StandardRecord* sr,
               double& weight) const override {
      
      restore.Add(sr->dune.eRecProxy,
                  sr->dune.LepE);

      const double scale = 0.01 * sigma;
      if (!sr->dune.isFD) { // in the ND
        if (sr->dune.isCC && abs(sr->dune.nuPDG) == 14) { // take away muon energy
          sr->dune.eRecProxy += (sr->dune.eRecProxy - sr->dune.LepE) * scale *
            pow((sr->dune.eRecProxy - sr->dune.LepE), 0.5);
        }
        else if (sr->dune.isCC && abs(sr->dune.nuPDG) == 12) { // fine to include electron energy
          sr->dune.eRecProxy += (sr->dune.eRecProxy) * scale *
            pow((sr->dune.eRecProxy), 0.5);
        } 
      }
    }
  };

  extern const TruthEnergySqrtND kTruthEnergySqrtND;

  // Total energy scale syst varying with sqrt of the energy
  class TruthEnergyInvSqrtND: public ISyst {
  public:
    TruthEnergyInvSqrtND() : ISyst("TruthEnergyInvSqrtND", "Inv Sqrt Total Energy Scale ND Syst") {}
    void Shift(double sigma,
               Restorer& restore,
               caf::StandardRecord* sr,
               double& weight) const override {

      restore.Add(sr->dune.eRecProxy,
                  sr->dune.LepE);

      const double scale = 0.02 * sigma;
      if (!sr->dune.isFD) { // in the ND
        if (sr->dune.isCC && abs(sr->dune.nuPDG) == 14) { // take away muon energy
          sr->dune.eRecProxy += (sr->dune.eRecProxy - sr->dune.LepE) * scale *
            pow((sr->dune.eRecProxy - sr->dune.LepE)+0.1, -0.5);
        }
        else if (sr->dune.isCC && abs(sr->dune.nuPDG) == 12) { // fine to include electron energy
          sr->dune.eRecProxy += (sr->dune.eRecProxy) * scale *
            pow((sr->dune.eRecProxy+0.1), -0.5);
        }
      }
    }
  };

  extern const TruthEnergyInvSqrtND kTruthEnergyInvSqrtND;

  //------------------------------------------------------------------------------
  // Electromagnetic

  // Systematic for pi0s and electrons
  // 2.5% on true energy of electrons and pi0s
  class EMTruthUncorrND : public ISyst {
  public:
    EMTruthUncorrND() : ISyst("EMTruthUncorrND", "EM Shower Uncorrelated ND Syst") {}
  
    void Shift(double sigma,
               Restorer& restore,
               caf::StandardRecord* sr, 
               double& weight) const override {

      restore.Add(sr->dune.LepE,
                  sr->dune.ePi0,
                  sr->dune.eRecProxy);

      const double scale = 0.025 * sigma;
      if (!sr->dune.isFD) { // in the ND
        sr->dune.ePi0 += sr->dune.ePi0 * scale;
        sr->dune.eRecProxy += sr->dune.ePi0 * scale;
        if (sr->dune.isCC && abs(sr->dune.nuPDG) == 12) {
          sr->dune.LepE += sr->dune.LepE * scale;
          sr->dune.eRecProxy += sr->dune.LepE * scale;
        }
      }
    }
  };  

  extern const EMTruthUncorrND kEMTruthUncorrND;

  // 2.5% systematic on EM energy for sqrt parameter
  //
  class EMTruthUncorrSqrtND : public ISyst {
  public:
    EMTruthUncorrSqrtND() : ISyst("EMTruthUncorrSqrtND", "EM Shower Uncorrelated ND Syst Sqrt") {}

    void Shift(double sigma,
               Restorer& restore,
               caf::StandardRecord* sr, 
               double& weight) const override {
      
      restore.Add(sr->dune.LepE,
                  sr->dune.ePi0,
                  sr->dune.eRecProxy);

      const double scale = 0.025 * sigma;
      if (!sr->dune.isFD) { // in the ND
        sr->dune.ePi0 += sr->dune.ePi0 * scale * pow(sr->dune.ePi0, 0.5);
        sr->dune.eRecProxy += sr->dune.ePi0 * scale * pow(sr->dune.ePi0, 0.5); 
        if (sr->dune.isCC && abs(sr->dune.nuPDG) == 12) {
          sr->dune.LepE += sr->dune.LepE * scale * pow(sr->dune.LepE, 0.5);
          sr->dune.eRecProxy += sr->dune.LepE * scale * pow(sr->dune.LepE, 0.5);
        }
      }
    }
  };

  extern const EMTruthUncorrSqrtND kEMTruthUncorrSqrtND;

  // 2.5% syst on EM energy for inv sqrt param
  //
  class EMTruthUncorrInvSqrtND : public ISyst {
  public: 
    EMTruthUncorrInvSqrtND() : ISyst("EMTruthUncorrInvSqrtND", "EM Shower Uncorrelated ND Syst Inv Sqrt") {}

    void Shift(double sigma,
               Restorer& restore,
               caf::StandardRecord* sr,
               double& weight) const override {

      restore.Add(sr->dune.LepE,
                  sr->dune.ePi0,
                  sr->dune.eRecProxy);

      const double scale = 0.025 * sigma;
      if (!sr->dune.isFD) { // in the ND
        sr->dune.ePi0 += sr->dune.ePi0 * scale * pow(sr->dune.ePi0+0.1, -0.5);
        sr->dune.eRecProxy += sr->dune.ePi0 * scale * pow(sr->dune.ePi0+0.1, -0.5);
        if (sr->dune.isCC && abs(sr->dune.nuPDG) == 12) {
          sr->dune.LepE += sr->dune.LepE * scale * pow(sr->dune.LepE+0.1, -0.5);
          sr->dune.eRecProxy += sr->dune.LepE * scale * pow(sr->dune.LepE+0.1, -0.5);
        }
      }
    }
  };

  extern const EMTruthUncorrInvSqrtND kEMTruthUncorrInvSqrtND;

  //---------------------------------------------------------------------
  // Charged Hadrons

  // Systematic for charged hadrons: p, pi+, pi-
  // 5% on true energy of charged pions and protons
  class ChargedHadTruthUncorrND : public ISyst {
  public:
    ChargedHadTruthUncorrND() : ISyst("ChargedHadTruthUncorrND", "Charged Hadron Uncorrelated ND Syst") {}

    void Shift(double sigma,
               Restorer& restore,
               caf::StandardRecord* sr,
               double& weight) const override {
    
      restore.Add(sr->dune.ePip,
                  sr->dune.ePim,
                  sr->dune.eP,
                  sr->dune.eRecProxy);

      const double scale = 0.05 * sigma;
      const double sumE = sr->dune.ePip + sr->dune.ePim + sr->dune.eP;
      if (!sr->dune.isFD) { // in the ND
        sr->dune.ePip += sr->dune.ePip * scale;
        sr->dune.ePim += sr->dune.ePim * scale;
        sr->dune.eP += sr->dune.eP * scale;
        sr->dune.eRecProxy += sumE * scale;
      }
    }
  };

  extern const ChargedHadTruthUncorrND kChargedHadTruthUncorrND;

  // 5% syst for charged hadrons sqrt param
  //
  class ChargedHadTruthUncorrSqrtND : public ISyst {
  public: 
    ChargedHadTruthUncorrSqrtND() : ISyst("ChargedHadTruthUncorrSqrtND", "Charged Had Uncorrelated Sqrt ND Syst") {}
    
    void Shift(double sigma,
               Restorer& restore,
               caf::StandardRecord* sr,
               double& weight) const override {

      restore.Add(sr->dune.ePip,
                  sr->dune.ePim,
                  sr->dune.eP,
                  sr->dune.eRecProxy);

      const double scale = 0.05 * sigma;
      const double sumE = sr->dune.ePip + sr->dune.ePim + sr->dune.eP;
      if (!sr->dune.isFD) { // in the ND
        sr->dune.ePip += sr->dune.ePip * scale * pow(sumE, 0.5);
        sr->dune.ePim += sr->dune.ePim * scale * pow(sumE, 0.5);
        sr->dune.eP += sr->dune.eP * scale * pow(sumE, 0.5); 
        sr->dune.eRecProxy += sumE * scale * pow(sumE, 0.5);
      }
    }
  };

  extern const ChargedHadTruthUncorrSqrtND kChargedHadTruthUncorrSqrtND;

  // 5% syst for charged hadrons inv sqrt param
  //
  class ChargedHadTruthUncorrInvSqrtND : public ISyst {
  public:
    ChargedHadTruthUncorrInvSqrtND() : ISyst("ChargedHadTruthUncorrInvSqrtND", "Charged Had Uncorrelated Inv Sqrt ND Syst") {}

    void Shift(double sigma,
               Restorer& restore,
               caf::StandardRecord* sr,
               double& weight) const override {

      restore.Add(sr->dune.ePip,
                  sr->dune.ePim,
                  sr->dune.eP,
                  sr->dune.eRecProxy);

      const double scale = 0.05 * sigma;
      const double sumE = sr->dune.ePip + sr->dune.ePim + sr->dune.eP;
      if (!sr->dune.isFD) { // in the ND
        sr->dune.ePip += sr->dune.ePip * scale * pow(sumE+0.1, -0.5);
        sr->dune.ePim += sr->dune.ePim * scale * pow(sumE+0.1, -0.5);
        sr->dune.eP += sr->dune.eP * scale * pow(sumE+0.1, -0.5);
        sr->dune.eRecProxy += sumE * scale * pow(sumE+0.1, -0.5);
      }
    }
  };

  extern const ChargedHadTruthUncorrInvSqrtND kChargedHadTruthUncorrInvSqrtND;

  //-------------------------------------------------------------------
  // Muons

  // Systematic for Muon 
  // 2% on true energy of muon from CC numu event
  class ETruthScaleMuLArND : public ISyst {
  public:
    ETruthScaleMuLArND() : ISyst("ETruthScaleMuLArND", "Muon Energy Scale ND Syst") {}

    void Shift(double sigma,
               Restorer& restore,
               caf::StandardRecord* sr,
               double& weight) const override {

      restore.Add(sr->dune.LepE,
                  sr->dune.eRecProxy);
 
      const double scale = 0.02 * sigma;
      if (!sr->dune.isFD) { // in the ND
        if (sr->dune.isCC && abs(sr->dune.nuPDG) == 14) { 
          sr->dune.LepE += sr->dune.LepE * scale;
          sr->dune.eRecProxy += sr->dune.LepE * scale;
        }
      }
    } 
  };

  extern const ETruthScaleMuLArND kETruthScaleMuLArND;

  // true energy of muon from CC numu event
  // Sqrt param 2% for LArTPC
  class ETruthScaleMuLArSqrtND : public ISyst {
  public:
    ETruthScaleMuLArSqrtND() : ISyst("ETruthScaleMuLArSqrtND", "Muon E Scale Sqrt ND Syst") {}

    void Shift(double sigma,
               Restorer& restore,
               caf::StandardRecord* sr,
               double& weight) const override {

      restore.Add(sr->dune.LepE,
                  sr->dune.eRecProxy);
 
      // 2% in LArTPC reco by range
      const double scale = 0.02 * sigma; // is 0.005 in FD
      if (!sr->dune.isFD) { // in the ND
        if (sr->dune.isCC && abs(sr->dune.nuPDG) == 14) {
          sr->dune.LepE += sr->dune.LepE * scale * pow(sr->dune.LepE, 0.5);
          sr->dune.eRecProxy += sr->dune.LepE * scale * pow(sr->dune.LepE, 0.5);
        }
      }
    }
  };

  extern const ETruthScaleMuLArSqrtND kETruthScaleMuLArSqrtND;

  // 2% on true energy of muon from CC numu event
  // Inv Sqrt param
  class ETruthScaleMuLArInvSqrtND : public ISyst {
  public:
    ETruthScaleMuLArInvSqrtND() : ISyst("ETruthScaleMuLArInvSqrtND", "Muon E Scale Inv Sqrt ND Syst") {}
  
    void Shift(double sigma,
               Restorer& restore,
               caf::StandardRecord* sr,
               double& weight) const override {

      restore.Add(sr->dune.LepE,
                  sr->dune.eRecProxy);

      const double scale = 0.02 * sigma;
      if (!sr->dune.isFD) { // in the ND
        if (sr->dune.isCC && abs(sr->dune.nuPDG) == 14) {
          sr->dune.LepE += sr->dune.LepE * scale * pow(sr->dune.LepE+0.1, -0.5);
          sr->dune.eRecProxy += sr->dune.LepE * scale * pow(sr->dune.LepE+0.1, -0.5);
        }
      }
    }
  };

  extern const ETruthScaleMuLArInvSqrtND kETruthScaleMuLArInvSqrtND;

  //---------------------------------------------------------------------------------

  // Vector of the truth energy scale systematics
  struct TruthEnergyNDSystVector : public std::vector<const ISyst*> {};

  TruthEnergyNDSystVector GetTrueENDSysts();


}
