#include "generator/LLBkg_Generator.h"

// from ROOT
#include "TRandom.h"

// from project - generator
#include "generator/CompGeneratorFactory.h"
#include "generator/generator.h"
#include "generator/Observables.h"

// from project - configuration
#include "configuration/CompConfig.h"

namespace cptoymc {
namespace generator {


LLBkg_Generator::LLBkg_Generator() :
  CompGenerator(),
  params_mass_{0.},
  params_timeandcp_{1.,0.},
  params_timeresol_{0.,0.05},
  params_taggingeffs_{0.30,0.06,0.04},
  params_taggingOS_{0.4},
  params_taggingSS_{0.4},
  comp_cat_(-1001)
{
  
}

LLBkg_Generator::~LLBkg_Generator() {
  
}

void LLBkg_Generator::Configure(const configuration::CompConfig& comp_config) {
  auto config_ptree = comp_config.model_ptree();
  // Mass
  auto sub_config_ptree = config_ptree.get_child("Mass");
  params_mass_.expo  = sub_config_ptree.get("expo", params_mass_.expo);
  
  // TimeAndCP
  sub_config_ptree = config_ptree.get_child("TimeAndCP");
  params_timeandcp_.tau       = sub_config_ptree.get("tau",       params_timeandcp_.tau);
  params_timeandcp_.prod_asym = sub_config_ptree.get("prod_asym", params_timeandcp_.prod_asym);
  
  sub_config_ptree = config_ptree.get_child("TimeResol");
  params_timeresol_.bias  = sub_config_ptree.get("bias" , params_timeresol_.bias  );
  params_timeresol_.sigma = sub_config_ptree.get("sigma", params_timeresol_.sigma );
  
  // Tagging
  sub_config_ptree = config_ptree.get_child("Tagging");
  params_taggingeffs_.eff_OS    = sub_config_ptree.get("eff_OS"   ,params_taggingeffs_.eff_OS   );
  params_taggingeffs_.eff_SS    = sub_config_ptree.get("eff_SS"   ,params_taggingeffs_.eff_SS   );
  params_taggingeffs_.eff_SSOS  = sub_config_ptree.get("eff_SSOS" ,params_taggingeffs_.eff_SSOS );
  
  params_taggingOS_.omega   = sub_config_ptree.get("omega_OS" , params_taggingOS_.omega );
  params_taggingSS_.omega   = sub_config_ptree.get("omega_SS" , params_taggingSS_.omega );
  
}

void LLBkg_Generator::GenerateEvent(TRandom& rndm, Observables& observables) {

}
  
void LLBkg_Generator::GenerateMass(TRandom& rndm, ObservableReal& obs_mass_true, ObservableReal& obs_mass_meas) {
  unsigned int trials = 0;
//  while (trials < max_trials_) {
//    GenerateMassBreitWigner(rndm, params_mass_.mean, params_mass_.width, obs_mass_true.value_);
//    
//    GenerateResolSingleGauss(rndm, params_massresol_.bias, params_massresol_.sigma, obs_mass_true.value(), obs_mass_meas.value_);
//    
//    if (obs_mass_true.HasValidValue() && obs_mass_meas.HasValidValue()) {
//      break;
//    } else {
//      ++trials;
//    }
//    std::cout
//    << "Problem in generation: Maximum trials reached without generating valid values for "
//    << obs_mass_true.dim_name() << " and " << obs_mass_meas.dim_name() << " !!!"
//    << std::endl;
//  }
}

void LLBkg_Generator::GenerateTimeAndTrueTag(TRandom& rndm, ObservableReal& obs_time_true, ObservableInt& obs_tag_true, ObservableReal& obs_time_meas) {
  unsigned int trials = 0;
//  while (trials < max_trials_) {
//    GenerateCPV_P2PV(rndm, params_timeandcp_.prod_asym, params_timeandcp_.tau,
//                     params_timeandcp_.dGamma, params_timeandcp_.dm,
//                     params_timeandcp_.Sf, params_timeandcp_.Cf, params_timeandcp_.Df,
//                     obs_time_true.value_, obs_tag_true.value_);
//    
//    GenerateResolSingleGauss(rndm, params_timeresol_.bias, params_timeresol_.sigma,
//                             obs_time_true.value_, obs_time_meas.value_);
//    
//    
//    if (obs_time_true.HasValidValue() && obs_time_meas.HasValidValue() && obs_tag_true.HasValidValue()) {
//      break;
//    } else {
//      ++trials;
//    }
//    std::cout
//    << "Problem in generation: Maximum trials reached without generating valid values for "
//    << obs_tag_true.dim_name()  << ","
//    << obs_time_true.dim_name() << " and " << obs_time_meas.dim_name() << " !!!"
//    << std::endl;
//  }
}


void LLBkg_Generator::GenerateTagAndEta(TRandom& rndm, const ObservableInt& obs_tag_true,
                                                ObservableInt& obs_tag_OS, ObservableReal& obs_eta_OS,
                                                ObservableInt& obs_tag_SS, ObservableReal& obs_eta_SS,
                                                ObservableInt& obs_tag_class)
{
//  double random_val = rndm.Uniform();
//  
//  if (random_val < params_taggingeffs_.eff_OS) {
//    GenerateEtaFlat(rndm, obs_eta_OS.value_);
//    GenerateTag(rndm,tag_calib_func_omegaOS_,tag_calib_func_domegaOS_,
//                obs_tag_true.value(), obs_eta_OS.value_, obs_tag_OS.value_);
//    obs_tag_SS.value_ = 1;
//    obs_eta_SS.value_ = 0.5;
//    obs_tag_class.value_ = 1;
//  }
//  else if (random_val < (params_taggingeffs_.eff_OS + params_taggingeffs_.eff_SS)) {
//    GenerateEtaFlat(rndm, obs_eta_SS.value_);
//    GenerateTag(rndm,tag_calib_func_omegaSS_,tag_calib_func_domegaSS_,
//                obs_tag_true.value(), obs_eta_SS.value_, obs_tag_SS.value_);
//    obs_tag_OS.value_ = 1;
//    obs_eta_OS.value_ = 0.5;
//    obs_tag_class.value_ = -1;
//  }
//  else if (random_val < (  params_taggingeffs_.eff_OS
//                         + params_taggingeffs_.eff_SS
//                         + params_taggingeffs_.eff_SSOS) ) {
//    GenerateEtaFlat(rndm, obs_eta_OS.value_);
//    GenerateTag(rndm,tag_calib_func_omegaOS_,tag_calib_func_domegaOS_,
//                obs_tag_true.value(), obs_eta_OS.value_, obs_tag_OS.value_);
//    GenerateEtaFlat(rndm, obs_eta_SS.value_);
//    GenerateTag(rndm,tag_calib_func_omegaSS_,tag_calib_func_domegaSS_,
//                obs_tag_true.value(), obs_eta_SS.value_, obs_tag_SS.value_);
//    obs_tag_class.value_ = 10;
//  }
//  else {
//    obs_tag_SS.value_ = 1;
//    obs_eta_SS.value_ = 0.5;
//    obs_tag_OS.value_ = 1;
//    obs_eta_OS.value_ = 0.5;
//    obs_tag_class.value_ = 0;
//  }
}
 

} // namespace generator
} // namespace cptoymc
