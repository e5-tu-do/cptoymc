#ifndef CPTOYMC_GENERATOR_COMPGENERATORFACTORY_H
#define CPTOYMC_GENERATOR_COMPGENERATORFACTORY_H

// from STL
#include <algorithm>
#include <map>
#include <string>


namespace cptoymc {
  
// forward declarations
namespace configuration {
  class CompConfig;
}

namespace generator {

class CompGenerator;
  
//// Factory stuff
//template<typename T>
//std::shared_ptr<CompGenerator> createGenerator(const configuration::CompConfig& comp_config ) {
//  return std::make_shared<T>(new T(comp_config));
//}

template<typename T>
CompGenerator* createInstance(const configuration::CompConfig& comp_config ) {
  return new T(comp_config);
}

class GeneratorFactory {
public:
  GeneratorFactory();
  ~GeneratorFactory();
  
  CompGenerator CreateGenerator(const configuration::CompConfig& comp_config) const;
  
private:
  std::map<std::string,std::function<CompGenerator(const configuration::CompConfig&)>> map_generators;
  
};




//map_generators.emplace("BSig_CPV_P2VP",createGenerator<BSig_CPV_P2VP_Generator>);
//
//std::shared_ptr<CompGenerator> GetGenerator(const configuration::CompConfig& comp_config);

} // namespace generator
} // namespace cptoymc


#endif // CPTOYMC_GENERATOR_COMPGENERATORFACTORY_H


