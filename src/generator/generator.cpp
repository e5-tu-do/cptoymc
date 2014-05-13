#include "generator/generator.h"

// from STL
#include <cmath>
#include <iostream>

// from ROOT
#include "TRandom.h"

// from project
#include "configuration/configuration.h"


namespace cptoymc {
namespace generator {





void GenerateMassBreitWigner(TRandom& rndm, double par_mean, double par_gamma, double& obs_mass_true) {
  obs_mass_true = rndm.BreitWigner(par_mean, par_gamma);
}

void GenerateCPV_P2PV(TRandom& rndm, double par_prod_asym,
                      double par_tau, double par_dGamma, double par_dm,
                      double par_Sf, double par_Cf, double par_Df,
                      double& obs_time_true, int& obs_tag_true) {
  // helper quantities
  double prob_B = (1. - par_prod_asym)/2.;
  double gamma_min = 1./par_tau - std::abs(par_dGamma)/2.;
  
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
    val_pdf      = BCPV_PDF(val_t, val_d, par_tau, par_dGamma, par_dm, par_Sf, par_Cf, par_Df);
    val_envelope = BCPV_PDF_Envelope(val_t, gamma_min, par_Sf, par_Cf, par_Df);
    
    if (val_envelope < val_pdf) {
      std::cout << "WARNING: Envelope smaller than PDF!" << std::endl;
    }
    
    if(val_envelope*rndm.Uniform() > val_pdf) continue;
    else break;
  }
  
  obs_tag_true  = val_d;
  obs_time_true = val_t;
}

double BCPV_PDF(double t, double d, double tau, double dGamma, double dm,
                double Sf, double Cf, double Df) {
  return exp(-t/tau)*(cosh(dGamma*t/2.)+Df*sinh(dGamma*t/2.)+d*Cf*cos(dm*t)-d*Sf*sin(dm*t));
}
  
double BCPV_PDF_Envelope(double t, double gamma_min, double Sf, double Cf, double Df) {
  return exp(-t*gamma_min)*(1.+std::abs(Df)+sqrt(Sf*Sf+Cf*Cf));
}

void GenerateResolSingleGauss(TRandom& rndm, double par_bias, double par_sigma, double obs_true, double& obs_meas) {
  obs_meas = obs_true;
  obs_meas += rndm.Gaus(par_bias, par_sigma);
}

void GenerateEtaFlat(TRandom& rndm, double& obs_eta) {
  obs_eta = rndm.Uniform(0.0,0.5);
}

void GenerateTag(TRandom& rndm, double par_omega, double par_domega, int obs_tag_true, int& obs_tag_meas) {
  if (rndm.Uniform() < (par_omega+(double)obs_tag_true*par_domega/2.) ) {
    obs_tag_meas = -1*obs_tag_true;
  } else {
    obs_tag_meas = obs_tag_true;
  }
}

void GenerateTag(TRandom& rndm,
                 std::function<double(double)>& func_omega,
                 std::function<double(double)>& func_domega,
                 int obs_tag_true, double obs_eta, int& obs_tag_meas) {
  GenerateTag(rndm, func_omega(obs_eta), func_domega(obs_eta), obs_tag_true, obs_tag_meas);
}
  
void GenerateRandomTag(TRandom& rndm, int& obs_tag_meas) {
  obs_tag_meas = (rndm.Uniform() < 0.5) ? +1 : -1;
}




int yieldToGenerate(TRandom& rndm, double yield_exp) {
  return rndm.Poisson(yield_exp);
}


} // namespace cptoymc
} // namespace generator

