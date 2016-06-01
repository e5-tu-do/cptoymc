#include "generator/Bkg_2Final_Generator.h"

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


Bkg_2Final_Generator::Bkg_2Final_Generator() :
  CompGenerator(),
  params_mass_{0.},
  params_timeandcp_{1.,0.,0.},
  params_timeresol_{0.033,0.72,0.,1.0},
  params_taggingeffs_{0.30,0.06,0.04},
  params_taggingOS_{0.5,0.,0.0,-1.0},
  params_taggingSS_{0.5,0.,0.0,-1.0},
  comp_cat_(-1001)
{
}

Bkg_2Final_Generator::~Bkg_2Final_Generator() {

}

void Bkg_2Final_Generator::Configure(const configuration::CompConfig& comp_config) {
  auto config_ptree = comp_config.model_ptree();

  comp_cat_ = comp_config.comp_cat();

  // Mass
  auto sub_config_ptree = config_ptree.get_child("Mass");
  params_mass_.expo  = sub_config_ptree.get("expo", params_mass_.expo);

  // TimeAndCP
  sub_config_ptree = config_ptree.get_child("TimeAndCP");
  params_timeandcp_.tau       = sub_config_ptree.get("tau", params_timeandcp_.tau);
  params_timeandcp_.prod_asym = sub_config_ptree.get("AP" , params_timeandcp_.prod_asym);
  params_timeandcp_.det_asym  = sub_config_ptree.get("AD" , params_timeandcp_.det_asym);

  sub_config_ptree = config_ptree.get_child("TimeResol");
  params_timeresol_.lognormal_m  = sub_config_ptree.get("lognormal_m", params_timeresol_.lognormal_m);
  params_timeresol_.lognormal_k  = sub_config_ptree.get("lognormal_k", params_timeresol_.lognormal_k);
  params_timeresol_.bias         = sub_config_ptree.get("bias"       , params_timeresol_.bias       );
  params_timeresol_.scale        = sub_config_ptree.get("scale"      , params_timeresol_.scale      );

  // Tagging
  sub_config_ptree = config_ptree.get_child("Tagging");
  params_taggingeffs_.eff_OS    = sub_config_ptree.get("eff_OS"   ,params_taggingeffs_.eff_OS   );
  params_taggingeffs_.eff_SS    = sub_config_ptree.get("eff_SS"   ,params_taggingeffs_.eff_SS   );
  params_taggingeffs_.eff_SSOS  = sub_config_ptree.get("eff_SSOS" ,params_taggingeffs_.eff_SSOS );

  params_taggingOS_.omega          = sub_config_ptree.get("omega_OS"         , params_taggingOS_.omega         );
  params_taggingOS_.domega         = sub_config_ptree.get("domega_OS"        , params_taggingOS_.domega        );
  params_taggingOS_.eta_dist_mean  = sub_config_ptree.get("eta_dist_mean_OS" , params_taggingOS_.eta_dist_mean );
  params_taggingOS_.eta_dist_sigma = sub_config_ptree.get("eta_dist_sigma_OS", params_taggingOS_.eta_dist_sigma);

  params_taggingSS_.omega          = sub_config_ptree.get("omega_SS"         , params_taggingSS_.omega         );
  params_taggingSS_.domega         = sub_config_ptree.get("domega_SS"        , params_taggingSS_.domega        );
  params_taggingSS_.eta_dist_mean  = sub_config_ptree.get("eta_dist_mean_SS" , params_taggingSS_.eta_dist_mean );
  params_taggingSS_.eta_dist_sigma = sub_config_ptree.get("eta_dist_sigma_SS", params_taggingSS_.eta_dist_sigma);
}

bool Bkg_2Final_Generator::TryGenerateEvent(TRandom& rndm, Observables& observables) {
  bool gen_success = true;

  observables.comp_cat.set_value(comp_cat_);

  gen_success &= GenerateMass(rndm, observables.mass_true, observables.mass_meas);
  gen_success &= GenerateTimeAndTrueTag(rndm, observables.time_true, observables.timeerror, observables.tag_true,
                                              observables.finalstate, observables.time_meas);
  gen_success &= GenerateTagAndEta(rndm, observables.tag_true,
                                   observables.tag_OS, observables.eta_OS,
                                   observables.tag_SS, observables.eta_SS,
                                   observables.tag_class);
  return gen_success;
}

bool Bkg_2Final_Generator::GenerateMass(TRandom& rndm, ObservableReal& obs_mass_true, ObservableReal& obs_mass_meas) {
  bool gen_success = true;

  obs_mass_true.value_ = -1000.;
  gen_success &= GenerateExpo(rndm,-1.*params_mass_.expo,obs_mass_meas.value_,
                              obs_mass_meas.min_value(),obs_mass_meas.max_value());

  return gen_success;
}

bool Bkg_2Final_Generator::GenerateTimeAndTrueTag(TRandom& rndm, ObservableReal& obs_time_true, ObservableReal& obs_timeerror, ObservableInt&
                                             obs_tag_true, ObservableInt& obs_finalstate, ObservableReal& obs_time_meas) {

  unsigned int trials = 0;
  bool gen_success = true;

  while (trials < max_trials_) {
    gen_success = true;

    // true "tag" according to prodasym
    obs_tag_true.value_ = (rndm.Uniform() < (1. - params_timeandcp_.prod_asym)/2.) ? +1 : -1;

    // finalstate according to detasym
    obs_finalstate.value_ = (rndm.Uniform() < (1. - params_timeandcp_.det_asym)/2.) ? +1 : -1;

    // decay time
    gen_success &= GenerateExpo(rndm,1./params_timeandcp_.tau,obs_time_true.value_,obs_time_true.min_value(),obs_time_true.max_value());

    gen_success &= GenerateLognormal(rndm, params_timeresol_.lognormal_m, params_timeresol_.lognormal_k, obs_timeerror.min_value(), obs_timeerror.max_value(), obs_timeerror.value_);
    gen_success &= GenerateResolSingleGaussPerEvent(rndm, params_timeresol_.bias, params_timeresol_.scale, obs_timeerror.value_, obs_time_true.value(), obs_time_meas.value_);

    if (gen_success && obs_time_true.HasValidValue() && obs_time_meas.HasValidValue() && obs_tag_true.HasValidValue() && obs_finalstate.HasValidValue()) {
      break;
    } else {
      ++trials;
    }
  }
  if (!gen_success && trials >= max_trials_) {
    std::cout
    << "Problem in LLBkg generation for component " << comp_cat_
    << ": Maximum trials reached without generating valid values for "
    << obs_tag_true.dim_name()  << ","
    << obs_time_true.dim_name() << " and " << obs_time_meas.dim_name() << " !!!"
    << std::endl;
  }

  return gen_success;
}


bool Bkg_2Final_Generator::GenerateTagAndEta(TRandom& rndm, const ObservableInt& obs_tag_true,
                                                ObservableInt& obs_tag_OS, ObservableReal& obs_eta_OS,
                                                ObservableInt& obs_tag_SS, ObservableReal& obs_eta_SS,
                                                ObservableInt& obs_tag_class)
{
  bool gen_success = true;
  double random_val = rndm.Uniform();

  if (random_val < params_taggingeffs_.eff_OS) { // generate OS tags and mistags
    // gen_success &= GenerateEtaFlat(rndm, obs_eta_OS.min_value(), obs_eta_OS.max_value(), obs_eta_OS.value_);
    gen_success &= GenerateEtaGauss(rndm, params_taggingOS_.eta_dist_mean, params_taggingOS_.eta_dist_sigma, obs_eta_OS.min_value(), obs_eta_OS.max_value(), obs_eta_OS.value_);
    gen_success &= GenerateTag(rndm, params_taggingOS_.omega, params_taggingOS_.domega,
                               obs_tag_true.GetValueForType("B") , obs_tag_true.GetValueForType("Bb"),
                               obs_tag_OS.GetValueForType("B")   , obs_tag_OS.GetValueForType("Bb"),
                               obs_tag_true.value(), obs_tag_OS.value_);

    obs_tag_SS.value_ = obs_tag_SS.GetValueForType("None");
    obs_eta_SS.value_ = 0.5;
    obs_tag_class.value_ = obs_tag_class.GetValueForType("OSonly");
  }
  else if (random_val < (params_taggingeffs_.eff_OS + params_taggingeffs_.eff_SS)) { // generate SS tags and mistags
    // gen_success &= GenerateEtaFlat(rndm, obs_eta_SS.min_value(), obs_eta_SS.max_value(), obs_eta_SS.value_);
    gen_success &= GenerateEtaGauss(rndm, params_taggingSS_.eta_dist_mean, params_taggingSS_.eta_dist_sigma, obs_eta_SS.min_value(), obs_eta_SS.max_value(), obs_eta_SS.value_);
    gen_success &= GenerateTag(rndm, params_taggingSS_.omega, params_taggingSS_.domega,
                               obs_tag_true.GetValueForType("B") , obs_tag_true.GetValueForType("Bb"),
                               obs_tag_SS.GetValueForType("B")   , obs_tag_SS.GetValueForType("Bb"),
                               obs_tag_true.value(), obs_tag_SS.value_);
    obs_tag_OS.value_ = obs_tag_OS.GetValueForType("None");
    obs_eta_OS.value_ = 0.5;
    obs_tag_class.value_ = obs_tag_class.GetValueForType("SSonly");
  }
  else if (random_val < (  params_taggingeffs_.eff_OS
                         + params_taggingeffs_.eff_SS
                         + params_taggingeffs_.eff_SSOS) ) { // generate overlap tags and mistags
    // gen_success &= GenerateEtaFlat(rndm, obs_eta_OS.min_value(), obs_eta_OS.max_value(), obs_eta_OS.value_);
    gen_success &= GenerateEtaGauss(rndm, params_taggingOS_.eta_dist_mean, params_taggingOS_.eta_dist_sigma, obs_eta_OS.min_value(), obs_eta_OS.max_value(), obs_eta_OS.value_);
    gen_success &= GenerateTag(rndm, params_taggingOS_.omega, params_taggingOS_.domega,
                               obs_tag_true.GetValueForType("B") , obs_tag_true.GetValueForType("Bb"),
                               obs_tag_OS.GetValueForType("B")   , obs_tag_OS.GetValueForType("Bb"),
                               obs_tag_true.value(), obs_tag_OS.value_);


    // gen_success &= GenerateEtaFlat(rndm, obs_eta_SS.min_value(), obs_eta_SS.max_value(), obs_eta_SS.value_);
    gen_success &= GenerateEtaGauss(rndm, params_taggingSS_.eta_dist_mean, params_taggingSS_.eta_dist_sigma, obs_eta_SS.min_value(), obs_eta_SS.max_value(), obs_eta_SS.value_);
    gen_success &= GenerateTag(rndm, params_taggingSS_.omega, params_taggingSS_.domega,
                               obs_tag_true.GetValueForType("B") , obs_tag_true.GetValueForType("Bb"),
                               obs_tag_SS.GetValueForType("B")   , obs_tag_SS.GetValueForType("Bb"),
                               obs_tag_true.value(), obs_tag_SS.value_);
    obs_tag_class.value_ = obs_tag_class.GetValueForType("OSandSS");
  }



  else { // untagged
    obs_tag_SS.value_ = obs_tag_SS.GetValueForType("None");
    obs_eta_SS.value_ = 0.5;
    obs_tag_OS.value_ = obs_tag_OS.GetValueForType("None");
    obs_eta_OS.value_ = 0.5;
    obs_tag_class.value_ = obs_tag_class.GetValueForType("untagged");
  }



  return gen_success;
}


} // namespace generator
} // namespace cptoymc
