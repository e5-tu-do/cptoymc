#include "generator/CompGeneratorFactory.h"

// from project
#include "generator/CompGenerator.h"

namespace cptoymc {
namespace generator {

//std::shared_ptr<CompGenerator> GetGenerator(const configuration::CompConfig& comp_config) {
//  return map_generators[comp_config.model()](comp_config);
//}

GeneratorFactory::GeneratorFactory() {
  //map_generators.emplace("BSig_CPV_P2VP",&createInstance<BSig_CPV_P2VP_Generator>);
}

GeneratorFactory::~GeneratorFactory() {
  
}

CompGenerator GeneratorFactory::CreateGenerator(const configuration::CompConfig& comp_config) const {
  
}

} // namespace generator
} // namespace cptoymc