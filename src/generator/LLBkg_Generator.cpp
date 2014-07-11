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
  params_taggingOS_{0.4,0.},
  params_taggingSS_{0.4,0.},
  comp_cat_(-1001)
{
  
}

LLBkg_Generator::~LLBkg_Generator() {
  
}

void LLBkg_Generator::Configure(const configuration::CompConfig& comp_config) {
  auto config_ptree = comp_config.model_ptree();

  comp_cat_ = comp_config.comp_cat();
  
  // Mass
  auto sub_config_ptree = config_ptree.get_child("Mass");
  params_mass_.expo  = sub_config_ptree.get("expo", params_mass_.expo);
  
  // TimeAndCP
  sub_config_ptree = config_ptree.get_child("TimeAndCP");
  params_timeandcp_.tau       = sub_config_ptree.get("tau", params_timeandcp_.tau);
  params_timeandcp_.prod_asym = sub_config_ptree.get("AP" , params_timeandcp_.prod_asym);
  
  sub_config_ptree = config_ptree.get_child("TimeResol");
  params_timeresol_.bias  = sub_config_ptree.get("bias" , params_timeresol_.bias  );
  params_timeresol_.sigma = sub_config_ptree.get("sigma", params_timeresol_.sigma );
  
  // Tagging
  sub_config_ptree = config_ptree.get_child("Tagging");
  params_taggingeffs_.eff_OS    = sub_config_ptree.get("eff_OS"   ,params_taggingeffs_.eff_OS   );
  params_taggingeffs_.eff_SS    = sub_config_ptree.get("eff_SS"   ,params_taggingeffs_.eff_SS   );
  params_taggingeffs_.eff_SSOS  = sub_config_ptree.get("eff_SSOS" ,params_taggingeffs_.eff_SSOS );
  
  params_taggingOS_.omega   = sub_config_ptree.get("omega_OS"  , params_taggingOS_.omega  );
  params_taggingOS_.domega  = sub_config_ptree.get("domega_OS" , params_taggingOS_.domega );
  params_taggingSS_.omega   = sub_config_ptree.get("omega_SS"  , params_taggingSS_.omega );
  params_taggingSS_.domega  = sub_config_ptree.get("domega_SS" , params_taggingSS_.domega );
}

void LLBkg_Generator::GenerateEvent(TRandom& rndm, Observables& observables) {
  observables.comp_cat.set_value(comp_cat_);
  
  GenerateMass(rndm, observables.mass_true, observables.mass_meas);
  GenerateTimeAndTrueTag(rndm, observables.time_true, observables.tag_true, observables.time_meas);
  GenerateTagAndEta(rndm, observables.tag_true,
                    observables.tag_OS, observables.eta_OS,
                    observables.tag_SS, observables.eta_SS,
                    observables.tag_class);
}
  
void LLBkg_Generator::GenerateMass(TRandom& rndm, ObservableReal& obs_mass_true, ObservableReal& obs_mass_meas) {
  
  obs_mass_true.value_ = -1000.;
  GenerateExpo(rndm,-1.*params_mass_.expo,obs_mass_meas.value_,obs_mass_meas.min_value(),obs_mass_meas.max_value());
}

void LLBkg_Generator::GenerateTimeAndTrueTag(TRandom& rndm, ObservableReal& obs_time_true, ObservableInt& obs_tag_true, ObservableReal& obs_time_meas) {
  
  // true "tag" according to prodasym
  obs_tag_true.value_ = (rndm.Uniform() < (1. - params_timeandcp_.prod_asym)/2.) ? +1 : -1;
  
  // decay time
  GenerateExpo(rndm,1./params_timeandcp_.tau,obs_time_true.value_,obs_time_true.min_value(),obs_time_true.max_value());
  
  GenerateResolSingleGauss(rndm, params_timeresol_.bias, params_timeresol_.sigma, obs_time_true.value(), obs_time_meas.value_);
}


void LLBkg_Generator::GenerateTagAndEta(TRandom& rndm, const ObservableInt& obs_tag_true,
                                                ObservableInt& obs_tag_OS, ObservableReal& obs_eta_OS,
                                                ObservableInt& obs_tag_SS, ObservableReal& obs_eta_SS,
                                                ObservableInt& obs_tag_class)
{
  double random_val = rndm.Uniform();
  
  if (random_val < params_taggingeffs_.eff_OS) {
    GenerateEtaFlat(rndm, obs_eta_OS.value_);
    GenerateTag(rndm, params_taggingOS_.omega, params_taggingOS_.domega,
                obs_tag_true.value(), obs_tag_OS.value_);
    obs_tag_SS.value_ = 1;
    obs_eta_SS.value_ = 0.5;
    obs_tag_class.value_ = 1;
  }
  else if (random_val < (params_taggingeffs_.eff_OS + params_taggingeffs_.eff_SS)) {
    GenerateEtaFlat(rndm, obs_eta_SS.value_);
    GenerateTag(rndm, params_taggingSS_.omega, params_taggingSS_.domega,
                obs_tag_true.value(), obs_tag_SS.value_);
    obs_tag_OS.value_ = 1;
    obs_eta_OS.value_ = 0.5;
    obs_tag_class.value_ = -1;
  }
  else if (random_val < (  params_taggingeffs_.eff_OS
                         + params_taggingeffs_.eff_SS
                         + params_taggingeffs_.eff_SSOS) ) {
    GenerateEtaFlat(rndm, obs_eta_OS.value_);
    GenerateTag(rndm, params_taggingOS_.omega, params_taggingOS_.domega,
                obs_tag_true.value(), obs_tag_OS.value_);
    GenerateEtaFlat(rndm, obs_eta_SS.value_);
    GenerateTag(rndm, params_taggingSS_.omega, params_taggingSS_.domega,
                obs_tag_true.value(), obs_tag_SS.value_);
    obs_tag_class.value_ = 10;
  }
  else {
    obs_tag_SS.value_ = 1;
    obs_eta_SS.value_ = 0.5;
    obs_tag_OS.value_ = 1;
    obs_eta_OS.value_ = 0.5;
    obs_tag_class.value_ = 0;
  }
}
 

} // namespace generator
} // namespace cptoymc
