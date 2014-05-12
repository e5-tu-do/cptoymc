#ifndef CPTOYMC_GENERATOR_TOYGENERATOR_H
#define CPTOYMC_GENERATOR_TOYGENERATOR_H

// STL
#include <memory>
#include <string>

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
  ToyGenerator(const cptoymc::configuration::ToyConfig& config);
  ~ToyGenerator();

 void GenerateToy(TTree& out_tree);

private:
  ToyGenerator();
  const configuration::ToyConfig& config_;
  
};


} // namespace generator
} // namespace cptoymc

#endif // CPTOYMC_GENERATOR_TOYGENERATOR_H



