#ifndef CPTOYMC_GENERATOR_LLBKG_GENERATOR_H
#define CPTOYMC_GENERATOR_LLBKG_GENERATOR_H

// STL
#include <algorithm>
#include <map>
#include <memory>
#include <string>

// from project
#include "generator/CompGenerator.h"

namespace cptoymc {
namespace generator {

// forward declarations
class ObservableReal;
class ObservableInt;

class LLBkg_Generator : public CompGenerator {
public:
  LLBkg_Generator();
  virtual ~LLBkg_Generator();

  virtual void Configure(const configuration::CompConfig& comp_config);
  virtual bool TryGenerateEvent(TRandom& rndm, Observables& observables);

private:
  bool GenerateMass(TRandom& rndm, ObservableReal& obs_mass_true, ObservableReal& obs_mass_meas);
  bool GenerateTimeAndTrueTag(TRandom& rndm, ObservableReal& obs_time_true, ObservableReal& obs_timeerror,
                              ObservableInt& obs_tag_true, ObservableInt& obs_finalstate, ObservableReal& obs_time_meas);
  bool GenerateTagAndEta(TRandom& rndm, const ObservableInt& obs_tag_true,
                         ObservableInt& obs_tag_OS, ObservableReal& obs_eta_OS,
                         ObservableInt& obs_tag_SS, ObservableReal& obs_eta_SS,
                         ObservableInt& obs_tag_class);

  struct ParamsMass {
    double expo;
  } params_mass_;

  struct ParamsTimeAndCP {
    double tau;
    double prod_asym;
    double det_asym;
  } params_timeandcp_;

  struct ParamsTimeResol {
    double lognormal_m;
    double lognormal_k;
    double bias;
    double scale;
  } params_timeresol_;

  struct ParamsTaggingEffs {
    double eff_OS;
    double eff_SS;
    double eff_SSOS;
  } params_taggingeffs_;

  struct ParamsTagging {
    double omega;
    double domega;
    double eta_dist_mean;
    double eta_dist_sigma;
  } params_taggingOS_, params_taggingSS_;

  int comp_cat_;
};

} // namespace generator
} // namespace cptoymc

#endif // CPTOYMC_GENERATOR_LLBKG_GENERATOR_H
