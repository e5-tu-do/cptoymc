#ifndef CPTOYMC_GENERATOR_TOYGENERATOR_H
#define CPTOYMC_GENERATOR_TOYGENERATOR_H

// STL
#include <memory>
#include <string>
#include <map>

// from project
#include "Observables.h"

class TRandom;
class TTree;

namespace cptoymc {
namespace configuration {
// forward declarations
class ToyConfig;

} // namespace configuration
  
namespace generator {
// forward declarations
class CompGenerator;
  
  
class ToyGenerator {
public:
  ToyGenerator(const cptoymc::configuration::ToyConfig& config, unsigned int seed = 0);
  ~ToyGenerator();

  void GenerateToy(TTree& out_tree);
  void GenerateToy(TTree& out_tree, unsigned int seed);

private:
  ToyGenerator();
  const configuration::ToyConfig& config_;
  Observables obs_;
  std::map<std::string,std::shared_ptr<CompGenerator>> comp_generators_;
  
  unsigned int seed_;
  TRandom* rndm_;
};


} // namespace generator
} // namespace cptoymc

#endif // CPTOYMC_GENERATOR_TOYGENERATOR_H



