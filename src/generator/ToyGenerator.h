#ifndef CPTOYMC_GENERATOR_TOYGENERATOR_H
#define CPTOYMC_GENERATOR_TOYGENERATOR_H

// STL
#include <memory>

class TRandom;
class TTree;

namespace cptoymc {
namespace configuration {
// forward declarations
class ToyConfig;

} // namespace configuration
  
namespace generator {
// forward declaration
class Observables;


class ToyGenerator {
public:
  ToyGenerator();
  ToyGenerator(cptoymc::configuration::ToyConfig config);
  ~ToyGenerator();

  std::unique_ptr<TTree> GenerateToy();

private:
  void GenerateBd(TRandom& rndm, const cptoymc::configuration::ToyConfig& config, Observables& observables);
  void GenerateBs(TRandom& rndm, const cptoymc::configuration::ToyConfig& config, Observables& observables);
  void GenerateBkg(TRandom& rndm, const cptoymc::configuration::ToyConfig& config, Observables& observables);
};


} // namespace generator
} // namespace cptoymc

#endif // CPTOYMC_GENERATOR_TOYGENERATOR_H



