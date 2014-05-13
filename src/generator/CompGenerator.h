#ifndef CPTOYMC_GENERATOR_COMPGENERATOR_H
#define CPTOYMC_GENERATOR_COMPGENERATOR_H


// STL
#include <algorithm>
#include <map>
#include <memory>
#include <string>

// from project
#include "generator/CompGeneratorFactory.h"


class TRandom;
class TTree;



namespace cptoymc {
  
// forward declarations
namespace configuration {
  class CompConfig;
}

namespace generator {
  
// forward declarations
class Observables;

// Generator classes
class CompGenerator {
public:
  CompGenerator();
  virtual ~CompGenerator();

  virtual void Configure(const configuration::CompConfig& comp_config) = 0;
  virtual void GenerateEvent(TRandom& rndm, Observables& observables)  = 0;
  
protected:

};

class BSig_CPV_P2VP_Generator : public CompGenerator {
 public:
  BSig_CPV_P2VP_Generator();
  virtual ~BSig_CPV_P2VP_Generator();
  
  virtual void Configure(const configuration::CompConfig& comp_config);
  virtual void GenerateEvent(TRandom& rndm, Observables& observables);

 private:
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
    double bias;
    double sigma;
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
  } params_taggingOS_, params_taggingSS_;
  
  int comp_cat_;
  
  std::function<double(double)> tag_calib_func_omegaOS_;
  std::function<double(double)> tag_calib_func_domegaOS_;
  std::function<double(double)> tag_calib_func_omegaSS_;
  std::function<double(double)> tag_calib_func_domegaSS_;

};
  
  

class LLBkg_Generator : public CompGenerator {
public:
  LLBkg_Generator();
  virtual ~LLBkg_Generator();

  virtual void Configure(const configuration::CompConfig& comp_config);
  virtual void GenerateEvent(TRandom& rndm, Observables& observables);
};

// register classes to Factory
static CompGeneratorRegistrar<BSig_CPV_P2VP_Generator> registrar("BSig_CPV_P2VP");
  
} // namespace generator
} // namespace cptoymc

    
#endif // CPTOYMC_GENERATOR_COMPGENERATOR_H
