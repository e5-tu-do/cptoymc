// STL
#include <cmath>
#include <iostream>
#include <vector>

// ROOT
#include "TBranch.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"

struct Observables {
  double mass_true;
  double time_true;
  int    tag_true;
  double mass;
  double time;
  int    tag_class; // -1: SS&&!OS, 0: !OS&&!SS, 1: OS&&!SS, 10: OS&&SS
  int    tag_OS;
  double eta_OS;
  int    tag_SS;
  double eta_SS;
}; 

struct ParsMass {
  double mean;
  double gamma;
};

struct ParsTimeAndCP {
  double tau;
  double dGamma;
  double dm;
  double Sf;
  double Cf;
  double Df;
  double AP;
};

struct ParsMassResol  {
  double bias;
  double sigma;
};

struct ParsTimeResol {
  double bias;
  double sigma;
};

struct ParsTagging {
  double eff_OS;
  double eff_SS;
  double eff_SSOS;
  double w_OS;
  double dw_OS;
  double w_SS;
  double dw_SS;
};


void resetObs(Observables& obs);

void generateTrueMass(TRandom& rndm, const ParsMass& pars_mass, double& obs_mass_true);
void generateTrueTimeAndTag(TRandom& rndm, const ParsTimeAndCP& pars_time_CP, double& obs_time_true, int& obs_tag_true);

void generateMass(TRandom& rndm, const ParsMassResol& pars_mass_resol, const double obs_mass_true, double& obs_mass);
void generateTime(TRandom& rndm, const ParsTimeResol& pars_time_resol, const double obs_time_true, double& obs_time);
void generateTagAndEta(TRandom& rndm, const ParsTagging& pars_tagging, const int tag_true, int& tag_OS, double& eta_OS, int& tag_SS, double& eta_SS, int& tag_class);

void generateTagAndEtaOS(TRandom& rndm, const ParsTagging& pars_tagging, const int tag_true, int& tag_OS, double& eta_OS);
void generateTagAndEtaSS(TRandom& rndm, const ParsTagging& pars_tagging, const int tag_true, int& tag_SS, double& eta_SS);

double BCPV_PDF(double t, double d, double tau, double dGamma, double dm, double Sf, double Cf, double Df);
double BCPV_PDF_Envelope(double t, double gamma_min, double Sf, double Cf, double Df);

int main() {
  int num_events = 100000;

  //============================================================================
  // Physics Parameters
  //----------------------------------------------------------------------------
  // Mass parameters
  ParsMass pars_mass;
  pars_mass.mean  = 5279.58;
  pars_mass.gamma = 0.;

  //----------------------------------------------------------------------------
  // Lifetime and CPV parameters
  ParsTimeAndCP pars_time_CP;
  pars_time_CP.tau    = 1.519;
  pars_time_CP.dGamma = 0.;
  pars_time_CP.dm     = 0.510;
  pars_time_CP.Sf     = 0.7;
  pars_time_CP.Cf     = 0.0;
  pars_time_CP.Df     = +sqrt(pars_time_CP.Sf*pars_time_CP.Sf+pars_time_CP.Cf*pars_time_CP.Cf);
  pars_time_CP.AP     = 0.0;

  //============================================================================
  // Experimental Parameters
  //----------------------------------------------------------------------------
  // Mass parameters
  ParsMassResol pars_mass_resol;
  pars_mass_resol.bias  = 0.;
  pars_mass_resol.sigma = 8.;

  //----------------------------------------------------------------------------
  // Decay time parameters
  ParsTimeResol pars_time_resol;
  pars_time_resol.bias  = 0.;
  pars_time_resol.sigma = 8.;

  //----------------------------------------------------------------------------
  // Tagging parameters
  ParsTagging pars_tagging;
  pars_tagging.eff_OS   = 0.3292;
  pars_tagging.eff_SS   = 0.0625;
  pars_tagging.eff_SSOS = 0.0436;

  pars_tagging.w_OS    = 0.25;
  pars_tagging.w_SS    = 0.25;

  pars_tagging.dw_OS    = 0.;
  pars_tagging.dw_SS    = 0.;


  //============================================================================
  // Initialise TFile and TTree
  TFile out_file("ToyMC.root","RECREATE");
  TTree out_tree("ToyMCTree","A tree filled with ToyMC observables");

  Observables obs;

  out_tree.Branch("obsMassTrue"       , &(obs.mass_true)  , "obsMassTrue/D"       );
  out_tree.Branch("obsTimeTrue"       , &(obs.time_true)  , "obsTimeTrue/D"       );
  out_tree.Branch("obsTagTrue"        , &(obs.tag_true)   , "obsTagTrue/I"        );
  out_tree.Branch("obsMass"           , &(obs.mass)       , "obsMass/D"           );
  out_tree.Branch("obsTime"           , &(obs.time)       , "obsTime/D"           );
  out_tree.Branch("catTaggedOSSSPion" , &(obs.tag_class)  , "catTaggedOSSSPion/I" );
  out_tree.Branch("obsTagOS"          , &(obs.tag_OS)     , "obsTagOS/I"          );
  out_tree.Branch("obsEtaOS"          , &(obs.eta_OS)     , "obsEtaOS/D"          );
  out_tree.Branch("obsTagSSPion"      , &(obs.tag_SS)     , "obsTagSSPion/I"      );
  out_tree.Branch("obsEtaSSPion"      , &(obs.eta_SS)     , "obsEtaSSPion/D"      );


  //============================================================================
  // Generate
  //----------------------------------------------------------------------------
  TRandom3 rndm(4357);

  for (int i = 0; i < num_events; ++i) {
    resetObs(obs);
    generateTrueMass(rndm, pars_mass, obs.mass_true);
    generateTrueTimeAndTag(rndm, pars_time_CP, obs.time_true, obs.tag_true);

    generateMass(rndm, pars_mass_resol, obs.mass_true, obs.mass);
    generateTime(rndm, pars_time_resol, obs.time_true, obs.time);
    generateTagAndEta(rndm, pars_tagging, obs.tag_true, obs.tag_OS, obs.eta_OS, obs.tag_SS, obs.eta_SS, obs.tag_class);


    out_tree.Fill();
  }


  //============================================================================
  // Finalise TFile
  out_file.Write();
  out_file.Close();

  return 0;
}


void resetObs(Observables& obs) {
  obs.mass_true = -1000.;
  obs.time_true = -1000.;
  obs.tag_true  = 0;
  obs.mass      = -1000.;
  obs.time      = -1000.; 
  obs.tag_class = -10000;
  obs.tag_OS    = 0;
  obs.eta_OS    = 0.5;
  obs.tag_SS    = 0;
  obs.eta_SS    = 0.5;
}

void generateTrueMass(TRandom& rndm, const ParsMass& pars_mass, double& obs_mass_true) {
  obs_mass_true = rndm.BreitWigner(pars_mass.mean, pars_mass.gamma);
}

void generateTrueTimeAndTag(TRandom& rndm, const ParsTimeAndCP& pars_time_CP, double& obs_time_true, int& obs_tag_true) {
  // helper quantities
  double prob_B = (1. - pars_time_CP.AP)/2.;
  double gamma_min = 1./pars_time_CP.tau - std::abs(pars_time_CP.dGamma)/2.;

  // local values of tag and time
  int    val_d = 0;
  double val_t = 0.;

  // now hit and miss, stolen from BDecay!!!
  double val_envelope = 0.;
  double val_pdf      = 0.;
  while(true) {   
    // get initial flavour while taking production asymmetry into account
    val_d = (rndm.Uniform() < prob_B) ? +1 : -1;
  
    // get time
    val_t = -log(rndm.Uniform())/gamma_min;

    // get pdf value and envelope value
    val_pdf      = BCPV_PDF(val_t, val_d, pars_time_CP.tau, pars_time_CP.dGamma, pars_time_CP.dm, pars_time_CP.Sf, pars_time_CP.Cf, pars_time_CP.Df);
    val_envelope = BCPV_PDF_Envelope(val_t, gamma_min, pars_time_CP.Sf, pars_time_CP.Cf, pars_time_CP.Df);

    if (val_envelope < val_pdf) {
      std::cout << "WARNING: Envelope smaller than PDF!" << std::endl;
    }

    if(val_envelope*rndm.Uniform() > val_pdf) continue;
    else break;
  }

  obs_tag_true  = val_d;
  obs_time_true = val_t;
}


void generateMass(TRandom& rndm, const ParsMassResol& pars_mass_resol, const double obs_mass_true, double& obs_mass) {
  obs_mass = obs_mass_true;
  obs_mass += rndm.Gaus(pars_mass_resol.bias, pars_mass_resol.sigma);
}

void generateTime(TRandom& rndm, const ParsTimeResol& pars_time_resol, const double obs_time_true, double& obs_time) {
  obs_time = obs_time_true;
  obs_time += rndm.Gaus(pars_time_resol.bias, pars_time_resol.sigma);
}

void generateTagAndEta(TRandom& rndm, const ParsTagging& pars_tagging, const int tag_true,  
                       int& tag_OS, double& eta_OS, int& tag_SS, double& eta_SS, int& tag_class) {

}

void generateTagAndEtaOS(TRandom& rndm, const ParsTagging& pars_tagging, const int tag_true, int& tag_OS, double& eta_OS) {

}

void generateTagAndEtaSS(TRandom& rndm, const ParsTagging& pars_tagging, const int tag_true, int& tag_SS, double& eta_SS) {

}

double BCPV_PDF(double t, double d, double tau, double dGamma, double dm, double Sf, double Cf, double Df) {
 return exp(-t/tau)*(cosh(dGamma*t/2.)+Df*sinh(dGamma*t/2.)+d*Cf*cos(dm*t)-d*Sf*sin(dm*t));
}

double BCPV_PDF_Envelope(double t, double gamma_min, double Sf, double Cf, double Df) {
  return exp(-t*gamma_min)*(1.+std::abs(Df)+sqrt(Sf*Sf+Cf*Cf));
}
