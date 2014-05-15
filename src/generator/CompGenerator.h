#ifndef CPTOYMC_GENERATOR_COMPGENERATOR_H
#define CPTOYMC_GENERATOR_COMPGENERATOR_H


class TRandom;

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

} // namespace generator
} // namespace cptoymc

    
#endif // CPTOYMC_GENERATOR_COMPGENERATOR_H
