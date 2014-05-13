#include "generator/CompGenerator.h"

// from STL
#include <iostream>

// from ROOT
#include "TRandom.h"

// from project - generator
#include "generator/generator.h"
#include "generator/Observables.h"

// from project - configuration
#include "configuration/CompConfig.h"



namespace cptoymc {
namespace generator {

CompGenerator::CompGenerator() { }


CompGenerator::~CompGenerator() { }

BSig_CPV_P2VP_Generator::BSig_CPV_P2VP_Generator() :
  CompGenerator(),
  params_mass_{5279.15, 0.},
  params_massresol_{0.,8.},
  params_timeandcp_{1.5,0.,0.5,0.7,0.,0.7,0.},
  params_timeresol_{0.,0.05},
  params_taggingeffs_{0.30,0.06,0.04},
  params_taggingOS_{1.0,0.25,0.25,0.0,0.0},
  params_taggingSS_{1.0,0.25,0.25,0.0,0.0},
  comp_cat_(-1000),
  tag_calib_func_omegaOS_(
    [=](double eta) -> double { return params_taggingOS_.p1*(eta-params_taggingOS_.etabar)+params_taggingOS_.p0; }
  ),
  tag_calib_func_domegaOS_(
    [=](double eta) -> double { return params_taggingOS_.dp1*(eta-params_taggingOS_.etabar)+params_taggingOS_.dp0; } ),
  tag_calib_func_omegaSS_(
    [=](double eta) -> double { return params_taggingSS_.p1*(eta-params_taggingSS_.etabar)+params_taggingSS_.p0; } ),
  tag_calib_func_domegaSS_(
    [=](double eta) -> double { return params_taggingSS_.dp1*(eta-params_taggingSS_.etabar)+params_taggingSS_.dp0; } )
{
  
}

BSig_CPV_P2VP_Generator::~BSig_CPV_P2VP_Generator() {
  
}

void BSig_CPV_P2VP_Generator::Configure(const configuration::CompConfig& comp_config) {

}

void BSig_CPV_P2VP_Generator::GenerateEvent(TRandom& rndm, Observables& observables) {
  observables.comp_cat = comp_cat_;
  
  GenerateMassBreitWigner(rndm, params_mass_.mean, params_mass_.width, observables.mass_true);
  
  GenerateResolSingleGauss(rndm, params_massresol_.bias, params_massresol_.sigma, observables.mass_true, observables.mass_meas);
  
  GenerateCPV_P2PV(rndm, params_timeandcp_.prod_asym, params_timeandcp_.tau,
                   params_timeandcp_.dGamma, params_timeandcp_.dm,
                   params_timeandcp_.Sf, params_timeandcp_.Cf, params_timeandcp_.Df,
                   observables.time_true, observables.tag_true);
  
  GenerateResolSingleGauss(rndm, params_timeresol_.bias, params_timeresol_.sigma,
                           observables.time_true, observables.time_meas);

  
  
  
  double random_val = rndm.Uniform();
  
  if (random_val < params_taggingeffs_.eff_OS) {
    GenerateEtaFlat(rndm, observables.eta_OS);
    GenerateTag(rndm,tag_calib_func_omegaOS_,tag_calib_func_domegaOS_,
                observables.tag_true, observables.eta_OS, observables.tag_OS);
    observables.tag_SS = 1;
    observables.eta_SS = 0.5;
    observables.tag_class = 1;
  }
  else if (random_val < (params_taggingeffs_.eff_OS + params_taggingeffs_.eff_SS)) {
    GenerateEtaFlat(rndm, observables.eta_SS);
    GenerateTag(rndm,tag_calib_func_omegaSS_,tag_calib_func_domegaSS_,
                observables.tag_true, observables.eta_SS, observables.tag_SS);
    observables.tag_OS = 1;
    observables.eta_OS = 0.5;
    observables.tag_class = -1;
  }
  else if (random_val < (   params_taggingeffs_.eff_OS
                          + params_taggingeffs_.eff_SS
                          + params_taggingeffs_.eff_SSOS) ) {
    GenerateEtaFlat(rndm, observables.eta_OS);
    GenerateTag(rndm,tag_calib_func_omegaOS_,tag_calib_func_domegaOS_,
                observables.tag_true, observables.eta_OS, observables.tag_OS);
    GenerateEtaFlat(rndm, observables.eta_SS);
    GenerateTag(rndm,tag_calib_func_omegaSS_,tag_calib_func_domegaSS_,
                observables.tag_true, observables.eta_SS, observables.tag_SS);
    observables.tag_class = 10;
  }
  else {
    observables.tag_SS = 1;
    observables.eta_SS = 0.5;
    observables.tag_OS = 1;
    observables.eta_OS = 0.5;
    observables.tag_class = 0;
  }
  
}



  
} // namespace generator
} // namespace cptoymc
