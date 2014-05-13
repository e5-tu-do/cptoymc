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

// True value generators and helper functions

void GenerateMassBreitWigner(TRandom& rndm, double par_mean, double par_gamma,
                             double& obs_mass_true);

void GenerateCPV_P2PV(TRandom& rndm, double par_prod_asym,
                      double par_tau, double par_dGamma, double par_dm,
                      double par_Sf, double par_Cf, double par_Df,
                      double& obs_time_true, int& obs_tag_true);

double BCPV_PDF(double t, double d, double tau, double dGamma, double dm,
                double Sf, double Cf, double Df);

double BCPV_PDF_Envelope(double t, double gamma_min,
                         double Sf, double Cf, double Df);
  
// Resolution generators
void GenerateResolSingleGauss(TRandom& rndm, double par_bias, double par_sigma,
                           double obs_true, double& obs_meas);
  
// Tagging generators

void GenerateEtaFlat(TRandom& rndm, double& obs_eta);

void GenerateTag(TRandom& rndm, double par_omega, double par_domega,
                 int obs_tag_true, int& obs_tag_meas);

void GenerateTag(TRandom& rndm,
                 std::function<double(double)>& func_omega,
                 std::function<double(double)>& func_domega,
                 int obs_tag_true, double obs_eta, int& obs_tag_meas);
  
void GenerateRandomTag(TRandom& rndm, int& obs_tag_meas);

// Background
//void GenerateMassBkg(TRandom& rndm, const cptoymc::configuration::ParsMass& pars_mass, double& obs_mass_true);
// void generateTrueTimeAndTag

int yieldToGenerate(TRandom& rndm, double yield_exp);

} // namespace generator
} // namespace cptoymc

#endif // CPTOYMC_GENERATOR_GENERATOR_H
