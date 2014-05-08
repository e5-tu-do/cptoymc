#ifndef CPTOYMC_GENERATOR_TOYGENERATOR_H
#define CPTOYMC_GENERATOR_TOYGENERATOR_H

// STL
#include <memory>

class TRandom;
class TTree;

namespace cptoymc {
namespace configuration {
class ToyConfig;
} // namespace configuration
  
namespace generator {

class ToyGenerator {
public:
  ToyGenerator();
  ToyGenerator(cptoymc::configuration::ToyConfig config);
  ~ToyGenerator();

  std::unique_ptr<TTree> GenerateToy();

private:

};


} // namespace generator
} // namespace cptoymc

#endif // CPTOYMC_GENERATOR_GENERATOR_H



#endif // CPTOYMC_GENERATOR_TOYGENERATOR_H
