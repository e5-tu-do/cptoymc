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


void generate_BSig_CPV_P2VP(TRandom& rndm, const configuration::CompConfig& comp_config, Observables& obs){
  
}



void generateTrueMass(TRandom& rndm, const cptoymc::configuration::ParsMass& pars_mass, double& obs_mass_true) {
  obs_mass_true = rndm.BreitWigner(pars_mass.mean, pars_mass.gamma);
}

void generateTrueTimeAndTag(TRandom& rndm, const cptoymc::configuration::ParsTimeAndCP& pars_time_CP, double& obs_time_true, int& obs_tag_true) {
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

void generateMass(TRandom& rndm, const cptoymc::configuration::ParsMassResol& pars_mass_resol, const double obs_mass_true, double& obs_mass) {
  obs_mass = obs_mass_true;
  obs_mass += rndm.Gaus(pars_mass_resol.bias, pars_mass_resol.sigma);
}

void generateTime(TRandom& rndm, const cptoymc::configuration::ParsTimeResol& pars_time_resol, const double obs_time_true, double& obs_time) {
  obs_time = obs_time_true;
  obs_time += rndm.Gaus(pars_time_resol.bias, pars_time_resol.sigma);
}

void generateTagAndEta(TRandom& rndm, const cptoymc::configuration::ParsTagging& pars_tagging, const int tag_true, int& tag_OS, double& eta_OS, int& tag_SS, double& eta_SS, int& tag_class) {
  
  double random_val = rndm.Uniform();
  
  if (random_val < pars_tagging.eff_OS) {
    generateTagAndEtaOS(rndm, pars_tagging, tag_true, tag_OS, eta_OS);
    tag_SS = 1;
    eta_SS = 0.5;
    tag_class = 1;
  }
  else if (random_val < (pars_tagging.eff_OS + pars_tagging.eff_SS)) {
    generateTagAndEtaSS(rndm, pars_tagging, tag_true, tag_SS, eta_SS);
    tag_OS = 1;
    eta_OS = 0.5;
    tag_class = -1;
  }
  else if (random_val < (pars_tagging.eff_OS + pars_tagging.eff_SS + pars_tagging.eff_SSOS)) {
    generateTagAndEtaOS(rndm, pars_tagging, tag_true, tag_OS, eta_OS);
    generateTagAndEtaSS(rndm, pars_tagging, tag_true, tag_SS, eta_SS);
    tag_class = 10;
  }
  else {
    tag_SS = 1;
    eta_SS = 0.5;
    tag_OS = 1;
    eta_OS = 0.5;
    tag_class = 0;
  }
  
}

void generateTagAndEtaOS(TRandom& rndm, const cptoymc::configuration::ParsTagging& pars_tagging, const int tag_true, int& tag_OS, double& eta_OS) {
  eta_OS = rndm.Uniform(0.0, 0.5);
  if (rndm.Uniform() < (eta_OS+(double)tag_true*pars_tagging.dw_OS/2.)){
    tag_OS = -1*tag_true;
  }
  else {
    tag_OS = tag_true;
  }
}

void generateTagAndEtaSS(TRandom& rndm, const cptoymc::configuration::ParsTagging& pars_tagging, const int tag_true, int& tag_SS, double& eta_SS) {
  eta_SS = rndm.Uniform(0.0, 0.44);
  if (rndm.Uniform() < (eta_SS+(double)tag_true*pars_tagging.dw_SS/2.)) {
    tag_SS = -1*tag_true;
  }
  else {
    tag_SS = tag_true;
  }
}


double BCPV_PDF(double t, double d, double tau, double dGamma, double dm, double Sf, double Cf, double Df) {
  return exp(-t/tau)*(cosh(dGamma*t/2.)+Df*sinh(dGamma*t/2.)+d*Cf*cos(dm*t)-d*Sf*sin(dm*t));
}

double BCPV_PDF_Envelope(double t, double gamma_min, double Sf, double Cf, double Df) {
  return exp(-t*gamma_min)*(1.+std::abs(Df)+sqrt(Sf*Sf+Cf*Cf));
}

int yieldToGenerate(TRandom& rndm, double yield_exp) {
  return rndm.Poisson(yield_exp);
}


} // namespace cptoymc
} // namespace generator

