#ifndef CPTOYMC_GENERATOR_GENERATOR_H
#define CPTOYMC_GENERATOR_GENERATOR_H

// STL
#include <memory>

class TRandom;
class TTree;

namespace cptoymc {


// forward declarations
namespace configuration {
struct ParsMass;
struct ParsTimeAndCP;
struct ParsMassResol;
struct ParsTimeResol;
struct ParsTagging;
class CompConfig;
} // namespace configuration

namespace generator {

// forward declarations
class Observables;

// generate functions

// True value generators and helper functions

void GenerateMassBreitWigner(TRandom& rndm, double par_mean, double par_gamma, double& obs_mass_true);

void GenerateCPV_P2PV(TRandom& rndm, const cptoymc::configuration::ParsTimeAndCP& pars_time_CP, double& obs_time_true, int& obs_tag_true);

double BCPV_PDF(double t, double d, double tau, double dGamma, double dm, double Sf, double Cf, double Df);
double BCPV_PDF_Envelope(double t, double gamma_min, double Sf, double Cf, double Df);

  
// Resolution generators
void ResolutionSingleGauss(TRandom& rndm, double par_bias, double par_sigma, double mass_true, double& mass_meas);
  
// Tagging generators

void GenerateEtaFlat(TRandom& rndm, double& obs_eta);

void GenerateTag(TRandom& rndm, double par_omega, double par_domega, int obs_tag_true, int& obs_tag_meas);

// Background
void GenerateMassBkg(TRandom& rndm, const cptoymc::configuration::ParsMass& pars_mass, double& obs_mass_true);
// void generateTrueTimeAndTag

int yieldToGenerate(TRandom& rndm, double yield_exp);

} // namespace generator
} // namespace cptoymc

#endif // CPTOYMC_GENERATOR_GENERATOR_H
