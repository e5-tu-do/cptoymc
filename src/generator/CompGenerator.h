#ifndef CPTOYMC_GENERATOR_COMPGENERATOR_H
#define CPTOYMC_GENERATOR_COMPGENERATOR_H


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
