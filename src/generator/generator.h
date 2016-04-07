#ifndef CPTOYMC_GENERATOR_GENERATOR_H
#define CPTOYMC_GENERATOR_GENERATOR_H

// STL
#include <functional>
#include <memory>

class TRandom;
class TTree;

namespace cptoymc {


// forward declarations
namespace configuration {
class CompConfig;
} // namespace configuration

namespace generator {

// forward declarations
class Observables;

// generate functions

bool GenerateExpo(TRandom& rndm, double par_expo, double& obs, double min, double max);

  
// True value generators and helper functions
  
bool GenerateMassBreitWigner(TRandom& rndm, double par_mean, double par_gamma,
                             double& obs_mass_true);

  
bool GenerateCPV_P2PV(TRandom& rndm, double par_prod_asym,
                      double par_tau, double par_dGamma, double par_dm,
                      double par_Sf, double par_Cf, double par_Df,
                      double& obs_time_true, int& obs_tag_true);

double BCPV_PDF(double t, double d, double tau, double dGamma, double dm,
                double Sf, double Cf, double Df);

double BCPV_PDF_Envelope(double t, double gamma_min,
                       double Sf, double Cf, double Df);

bool GenerateLognormal(TRandom& rndm, double m, double k, double min, double max,
                       double& obs_sigma_t);
  
// Resolution generators
bool GenerateResolSingleGauss(TRandom& rndm, double par_bias, double par_sigma,
                              double obs_true, double& obs_meas);
bool GenerateResolSingleGaussPerEvent(TRandom& rndm, double par_bias, double par_scale,
                              double obs_per_event_error, double obs_true, double& obs_meas);
  

// Tagging generators

bool GenerateEtaFlat(TRandom& rndm, double& obs_eta);
bool GenerateEtaFlat(TRandom& rndm, double obs_eta_min, double obs_eta_max, double& obs_eta);
bool GenerateEtaTimeCorrelated(double obs_eta_min, double obs_eta_max, double& obs_eta, const double obs_time_meas, double par_corr_eta_time_slope, double par_corr_eta_time_offset);
bool GenerateEtaGauss(TRandom& rndm, double m, double s, double obs_eta_min, double obs_eta_max, double& obs_eta);
bool GenerateEtaGaussTimeCorrelated(TRandom& rndm, double s, double obs_eta_min, double obs_eta_max, double& obs_eta, const double obs_time_meas, double par_corr_eta_time_slope, double par_corr_eta_time_offset);

bool GenerateTag(TRandom& rndm, double par_omega, double par_domega,
                 const int par_tag_true_B, const int par_tag_true_Bb,
                 const int par_tag_B,      const int par_tag_Bb,
                 int obs_tag_true, int& obs_tag_meas);
  
bool GenerateTag(TRandom& rndm, double par_omega, double par_domega,
                 int obs_tag_true, int& obs_tag_meas);

bool GenerateTag(TRandom& rndm,
                 std::function<double(double)>& func_omega,
                 std::function<double(double)>& func_domega,
                 const int par_tag_true_B, const int par_tag_true_Bb,
                 const int par_tag_B,      const int par_tag_Bb,
                 int obs_tag_true, double obs_eta, int& obs_tag_meas);
  
bool GenerateRandomTag(TRandom& rndm, int& obs_tag_meas);


// void generateTrueTimeAndTag

int yieldToGenerate(TRandom& rndm, double yield_exp);

} // namespace generator
} // namespace cptoymc

#endif // CPTOYMC_GENERATOR_GENERATOR_H
