#include "generator/BSig_CPV_P2VP_Generator.h"

// from STL
#include <iostream>

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


BSig_CPV_P2VP_Generator::BSig_CPV_P2VP_Generator() :
  CompGenerator(),
  params_mass_{5279.15, 0.},
  params_massresol_{0.,8.},
  params_timeandcp_{1.5,0.,0.5,0.7,0.,0.7,-2.0,-2.0,-2.0,0.,0.},
  params_timeresol_{0.033,0.72,0.,1.0},
  params_taggingeffs_{0.30,0.06,0.04},
  params_taggingOS_{1.0,0.25,0.25,0.0,0.0,0.0,-1.0},
  params_taggingSS_{1.0,0.25,0.25,0.0,0.0,0.0,-1.0},
  comp_cat_(-1000),
  tag_calib_func_omegaOS_(
    [&](double eta) -> double { return params_taggingOS_.p1*(eta-params_taggingOS_.etabar)+params_taggingOS_.p0; } ),
  tag_calib_func_domegaOS_(
    [&](double eta) -> double { return params_taggingOS_.dp1*(eta-params_taggingOS_.etabar)+params_taggingOS_.dp0; } ),
  tag_calib_func_omegaSS_(
    [&](double eta) -> double { return params_taggingSS_.p1*(eta-params_taggingSS_.etabar)+params_taggingSS_.p0; } ),
  tag_calib_func_domegaSS_(
    [&](double eta) -> double { return params_taggingSS_.dp1*(eta-params_taggingSS_.etabar)+params_taggingSS_.dp0; } )
{

}

BSig_CPV_P2VP_Generator::~BSig_CPV_P2VP_Generator() {

}

void BSig_CPV_P2VP_Generator::Configure(const configuration::CompConfig& comp_config) {
  comp_cat_ = comp_config.comp_cat();

  auto config_ptree = comp_config.model_ptree();

  // Mass
  auto sub_config_ptree = config_ptree.get_child("Mass");
  params_mass_.mean  = sub_config_ptree.get("mean", params_mass_.mean);
  params_mass_.width = sub_config_ptree.get("width",params_mass_.width);

  sub_config_ptree = config_ptree.get_child("MassResol");
  params_massresol_.bias  = sub_config_ptree.get("bias",params_massresol_.bias);
  params_massresol_.sigma = sub_config_ptree.get("sigma",params_massresol_.sigma);

  // TimeAndCP
  sub_config_ptree = config_ptree.get_child("TimeAndCP");
  params_timeandcp_.tau       = sub_config_ptree.get("tau",     params_timeandcp_.tau);
  params_timeandcp_.dGamma    = sub_config_ptree.get("dGamma",  params_timeandcp_.dGamma);
  params_timeandcp_.dm        = sub_config_ptree.get("dm",      params_timeandcp_.dm);
  params_timeandcp_.Sf        = sub_config_ptree.get("Sf",      params_timeandcp_.Sf);
  params_timeandcp_.Cf        = sub_config_ptree.get("Cf",      params_timeandcp_.Cf);
  params_timeandcp_.Df        = sub_config_ptree.get("Df",      params_timeandcp_.Df);
  params_timeandcp_.Sfbar     = sub_config_ptree.get("Sfbar",   params_timeandcp_.Sfbar);
  params_timeandcp_.Cfbar     = sub_config_ptree.get("Cfbar",   params_timeandcp_.Cfbar);
  params_timeandcp_.Dfbar     = sub_config_ptree.get("Dfbar",   params_timeandcp_.Dfbar);
  params_timeandcp_.prod_asym = sub_config_ptree.get("AP",      params_timeandcp_.prod_asym);
  params_timeandcp_.det_asym  = sub_config_ptree.get("AD",      params_timeandcp_.det_asym);

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

  params_taggingOS_.p1             = sub_config_ptree.get("p1_OS"            , params_taggingOS_.p1            );
  params_taggingOS_.p0             = sub_config_ptree.get("p0_OS"            , params_taggingOS_.p0            );
  params_taggingOS_.etabar         = sub_config_ptree.get("etabar_OS"        , params_taggingOS_.etabar        );
  params_taggingOS_.dp1            = sub_config_ptree.get("dp1_OS"           , params_taggingOS_.dp1           );
  params_taggingOS_.dp0            = sub_config_ptree.get("dp0_OS"           , params_taggingOS_.dp0           );
  params_taggingOS_.eta_dist_mean  = sub_config_ptree.get("eta_dist_mean_OS" , params_taggingOS_.eta_dist_mean );
  params_taggingOS_.eta_dist_sigma = sub_config_ptree.get("eta_dist_sigma_OS", params_taggingOS_.eta_dist_sigma);

  params_taggingSS_.p1             = sub_config_ptree.get("p1_SS"            , params_taggingSS_.p1            );
  params_taggingSS_.p0             = sub_config_ptree.get("p0_SS"            , params_taggingSS_.p0            );
  params_taggingSS_.etabar         = sub_config_ptree.get("etabar_SS"        , params_taggingSS_.etabar        );
  params_taggingSS_.dp1            = sub_config_ptree.get("dp1_SS"           , params_taggingSS_.dp1           );
  params_taggingSS_.dp0            = sub_config_ptree.get("dp0_SS"           , params_taggingSS_.dp0           );
  params_taggingSS_.eta_dist_mean  = sub_config_ptree.get("eta_dist_mean_SS" , params_taggingSS_.eta_dist_mean );
  params_taggingSS_.eta_dist_sigma = sub_config_ptree.get("eta_dist_sigma_SS", params_taggingSS_.eta_dist_sigma);
}

bool BSig_CPV_P2VP_Generator::TryGenerateEvent(TRandom& rndm, Observables& observables) {
  observables.comp_cat.set_value(comp_cat_);
  bool gen_success = true;
  gen_success &= GenerateMass(rndm, observables.mass_true, observables.mass_meas);
  //gen_success &= GenerateLognormal(rndm, m, k, observables.timeerror);
  gen_success &= GenerateTimeAndTrueTag(rndm, observables.time_true, observables.timeerror, observables.tag_true,
                                              observables.time_meas, observables.finalstate, observables.mix_true);
  gen_success &= GenerateTagAndEta(rndm, observables.tag_true,
                    observables.tag_OS, observables.eta_OS,
                    observables.tag_SS, observables.eta_SS,
                    observables.tag_class);

  return gen_success;
}

bool BSig_CPV_P2VP_Generator::GenerateMass(TRandom& rndm, ObservableReal& obs_mass_true, ObservableReal& obs_mass_meas) {
  bool gen_success = true;
  unsigned int trials = 0;
  while (trials < max_trials_) {
    gen_success = true;
    gen_success &= GenerateMassBreitWigner(rndm, params_mass_.mean, params_mass_.width, obs_mass_true.value_);

    gen_success &= GenerateResolSingleGauss(rndm, params_massresol_.bias, params_massresol_.sigma, obs_mass_true.value(), obs_mass_meas.value_);

    if (gen_success && obs_mass_true.HasValidValue() && obs_mass_meas.HasValidValue()) {
      break;
    } else {
      ++trials;
    }

    gen_success = false;
  }

  if (!gen_success && trials >= max_trials_) {
    std::cout
    << "Problem in BSig_CPV generation for component " << comp_cat_ <<": Maximum trials reached without generating valid values for "
    << obs_mass_true.dim_name() << " and " << obs_mass_meas.dim_name() << " !!!"
    << std::endl;
  }
  return gen_success;
}

bool BSig_CPV_P2VP_Generator::GenerateTimeAndTrueTag(TRandom& rndm, ObservableReal& obs_time_true, ObservableReal& obs_timeerror, ObservableInt& obs_tag_true, ObservableReal& obs_time_meas, ObservableInt& obs_finalstate, ObservableInt& obs_mix_true) {
  unsigned int trials = 0;
  bool gen_success = true;
  while (trials < max_trials_) {
    gen_success = true;
    gen_success &= GenerateCPV_P2PV(rndm, params_timeandcp_.prod_asym, params_timeandcp_.det_asym, params_timeandcp_.tau,
                     params_timeandcp_.dGamma, params_timeandcp_.dm,
                     params_timeandcp_.Sf, params_timeandcp_.Cf, params_timeandcp_.Df,
                     params_timeandcp_.Sfbar, params_timeandcp_.Cfbar, params_timeandcp_.Dfbar,
                     obs_time_true.value_, obs_tag_true.value_, obs_finalstate.value_, obs_mix_true.value_);

    gen_success &= GenerateLognormal(rndm, params_timeresol_.lognormal_m, params_timeresol_.lognormal_k, obs_timeerror.min_value(), obs_timeerror.max_value(), obs_timeerror.value_);
    gen_success &= GenerateResolSingleGaussPerEvent(rndm, params_timeresol_.bias, params_timeresol_.scale, obs_timeerror.value_,
                             obs_time_true.value_, obs_time_meas.value_);

    if (gen_success && obs_time_true.HasValidValue() && obs_time_meas.HasValidValue() && obs_tag_true.HasValidValue() && obs_finalstate.HasValidValue()) {
      break;
    } else {
      ++trials;
    }
  }
  if (!gen_success && trials >= max_trials_) {
    std::cout
    << "Problem in BSig_CPV generation for component " << comp_cat_
    << ": Maximum trials reached without generating valid values for "
    << obs_tag_true.dim_name()  << ","
    << obs_time_true.dim_name() << " and " << obs_time_meas.dim_name() << " !!!"
    << std::endl;
  }
  return gen_success;
}


bool BSig_CPV_P2VP_Generator::GenerateTagAndEta(TRandom& rndm, const ObservableInt& obs_tag_true,
                                            ObservableInt& obs_tag_OS, ObservableReal& obs_eta_OS,
                                            ObservableInt& obs_tag_SS, ObservableReal& obs_eta_SS,
                                            ObservableInt& obs_tag_class)
{
  bool gen_success = true;
  double random_val = rndm.Uniform();

  if (random_val < params_taggingeffs_.eff_OS) { // generate OS tags and mistags
    // gen_success &= GenerateEtaFlat(rndm, obs_eta_OS.min_value(), obs_eta_OS.max_value(), obs_eta_OS.value_);
    gen_success &= GenerateEtaGauss(rndm, params_taggingOS_.eta_dist_mean, params_taggingOS_.eta_dist_sigma, obs_eta_OS.min_value(), obs_eta_OS.max_value(), obs_eta_OS.value_);
    gen_success &= GenerateTag(rndm,tag_calib_func_omegaOS_,tag_calib_func_domegaOS_,
                               obs_tag_true.GetValueForType("B"), obs_tag_true.GetValueForType("Bb"),
                               obs_tag_OS.GetValueForType("B")  , obs_tag_OS.GetValueForType("Bb"),
                               obs_tag_true.value(), obs_eta_OS.value_, obs_tag_OS.value_);
    obs_tag_SS.value_ = obs_tag_SS.GetValueForType("None");
    obs_eta_SS.value_ = 0.5;
    obs_tag_class.value_ = obs_tag_class.GetValueForType("OSonly");
  }
  else if (random_val < (params_taggingeffs_.eff_OS + params_taggingeffs_.eff_SS)) { // generate SS tags and mistags
    // gen_success &= GenerateEtaFlat(rndm, obs_eta_SS.min_value(), obs_eta_SS.max_value(), obs_eta_SS.value_);
    gen_success &= GenerateEtaGauss(rndm, params_taggingSS_.eta_dist_mean, params_taggingSS_.eta_dist_sigma, obs_eta_SS.min_value(), obs_eta_SS.max_value(), obs_eta_SS.value_);
    gen_success &= GenerateTag(rndm,tag_calib_func_omegaSS_,tag_calib_func_domegaSS_,
                               obs_tag_true.GetValueForType("B"), obs_tag_true.GetValueForType("Bb"),
                               obs_tag_SS.GetValueForType("B"), obs_tag_SS.GetValueForType("Bb"),
                               obs_tag_true.value(), obs_eta_SS.value_, obs_tag_SS.value_);
    obs_tag_OS.value_ = obs_tag_OS.GetValueForType("None");
    obs_eta_OS.value_ = 0.5;
    obs_tag_class.value_ = obs_tag_class.GetValueForType("SSonly");
  }
  else if (random_val < (  params_taggingeffs_.eff_OS
                         + params_taggingeffs_.eff_SS
                         + params_taggingeffs_.eff_SSOS) ) { // generate overlap tags and mistags
    // gen_success &= GenerateEtaFlat(rndm, obs_eta_OS.min_value(), obs_eta_OS.max_value(), obs_eta_OS.value_);
    gen_success &= GenerateEtaGauss(rndm, params_taggingOS_.eta_dist_mean, params_taggingOS_.eta_dist_sigma, obs_eta_OS.min_value(), obs_eta_OS.max_value(), obs_eta_OS.value_);
    gen_success &= GenerateTag(rndm,tag_calib_func_omegaOS_,tag_calib_func_domegaOS_,
                               obs_tag_true.GetValueForType("B"), obs_tag_true.GetValueForType("Bb"),
                               obs_tag_OS.GetValueForType("B"), obs_tag_OS.GetValueForType("Bb"),
                               obs_tag_true.value(), obs_eta_OS.value_, obs_tag_OS.value_);

    // gen_success &= GenerateEtaFlat(rndm, obs_eta_SS.min_value(), obs_eta_SS.max_value(), obs_eta_SS.value_);
    gen_success &= GenerateEtaGauss(rndm, params_taggingSS_.eta_dist_mean, params_taggingSS_.eta_dist_sigma, obs_eta_SS.min_value(), obs_eta_SS.max_value(), obs_eta_SS.value_);
    gen_success &= GenerateTag(rndm,tag_calib_func_omegaSS_,tag_calib_func_domegaSS_,
                               obs_tag_true.GetValueForType("B"), obs_tag_true.GetValueForType("Bb"),
                               obs_tag_SS.GetValueForType("B"), obs_tag_SS.GetValueForType("Bb"),
                               obs_tag_true.value(), obs_eta_SS.value_, obs_tag_SS.value_);
    obs_tag_class.value_ = obs_tag_class.GetValueForType("OSandSS");
  }
  else { // untagged
    obs_tag_SS.value_ = obs_tag_SS.GetValueForType("None");
    obs_eta_SS.value_ = 0.5;
    obs_tag_OS.value_ = obs_tag_OS.GetValueForType("None");
    obs_eta_OS.value_ = 0.5;
    obs_tag_class.value_ = obs_tag_class.GetValueForType("untagged");;
  }
  return gen_success;
}

} // namespace generator
} // namespace cptoymc
