#ifndef CPTOYMC_GENERATOR_BSIG_CPV_P2VP_GENERATOR_H
#define CPTOYMC_GENERATOR_BSIG_CPV_P2VP_GENERATOR_H

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

  
class BSig_CPV_P2VP_Generator : public CompGenerator {
public:
  BSig_CPV_P2VP_Generator();
  virtual ~BSig_CPV_P2VP_Generator();
  
  virtual void Configure(const configuration::CompConfig& comp_config);
  virtual bool TryGenerateEvent(TRandom& rndm, Observables& observables);
  
private:
  bool GenerateMass(TRandom& rndm, ObservableReal& obs_mass_true, ObservableReal& obs_mass_meas);
  bool GenerateTimeAndTrueTag(TRandom& rndm, ObservableReal& obs_time_true, ObservableReal& obs_timeerror, ObservableInt& obs_tag_true, ObservableReal& obs_time_meas);
  bool GenerateTagAndEta(TRandom& rndm, const ObservableInt& obs_tag_true, const ObservableReal& obs_time_meas,
                         ObservableInt& obs_tag_OS, ObservableReal& obs_eta_OS,
                         ObservableInt& obs_tag_SS, ObservableReal& obs_eta_SS,
                         ObservableInt& obs_tag_class);
  
  struct ParamsMass {
    double mean;
    double width;
  } params_mass_;
  
  struct ParamsMassResol {
    double bias;
    double sigma;
  } params_massresol_;
  
  struct ParamsTimeAndCP {
    double tau;
    double dGamma;
    double dm;
    double Sf;
    double Cf;
    double Df;
    double prod_asym;
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
    double p1;
    double p0;
    double etabar;
    double dp1;
    double dp0;
    double eta_dist_mean;
    double eta_dist_sigma;
  } params_taggingOS_, params_taggingSS_;

  struct ParamsTaggingCorrelation {
    double slope;
    double offset;
  } params_etatimecorrOS_, params_etatimecorrSS_;
  
  int comp_cat_;
  
  std::function<double(double)> tag_calib_func_omegaOS_;
  std::function<double(double)> tag_calib_func_domegaOS_;
  std::function<double(double)> tag_calib_func_omegaSS_;
  std::function<double(double)> tag_calib_func_domegaSS_;
  
};

} // namespace generator
} // namespace cptoymc

#endif // CPTOYMC_GENERATOR_BSIG_CPV_P2VP_GENERATOR_H
