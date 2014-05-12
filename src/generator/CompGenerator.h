#ifndef CPTOYMC_GENERATOR_COMPGENERATOR_H
#define CPTOYMC_GENERATOR_COMPGENERATOR_H


// STL
#include <algorithm>
#include <map>
#include <memory>
#include <string>

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
  CompGenerator(const configuration::CompConfig& comp_config);
  virtual ~CompGenerator();
  
  virtual void generateEvent(TRandom& rndm, Observables& observables) = 0;
  
protected:
  const configuration::CompConfig& comp_config_;
};

class BSig_CPV_P2VP_Generator : public CompGenerator {
 public:
  BSig_CPV_P2VP_Generator(const configuration::CompConfig& comp_config);
  virtual ~BSig_CPV_P2VP_Generator();
  
  virtual void generateEvent(TRandom& rndm, Observables& observables);

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
    double AP;
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
    double eta;
    double dp1;
    double dp0;
  } params_taggingOS_, params_taggingSS_;
  
};
  
  

class LLBkg_Generator : public CompGenerator {
public:
  LLBkg_Generator(const configuration::CompConfig& comp_config);
  virtual ~LLBkg_Generator();
  
  virtual void generateEvent(TRandom& rndm, Observables& observables);
};

} // namespace generator
} // namespace cptoymc

    
#endif // CPTOYMC_GENERATOR_COMPGENERATOR_H
