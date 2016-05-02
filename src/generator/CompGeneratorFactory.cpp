#include "generator/CompGeneratorFactory.h"

// from STL

// from project
#include "configuration/CompConfig.h"

// from project - component generators
#include "generator/CompGenerator.h"
#include "generator/BSig_CPV_P2VP_Generator.h"
#include "generator/LLBkg_Generator.h"

namespace cptoymc {
namespace generator {


CompGeneratorFactory::CompGeneratorFactory() {
  RegisterGenerator<BSig_CPV_P2VP_Generator>("BSig_CPV_P2VP");
  RegisterGenerator<LLBkg_Generator>("LLBkg");
}

CompGeneratorFactory::~CompGeneratorFactory() { }

CompGeneratorFactory* CompGeneratorFactory::Instance() {
  static CompGeneratorFactory generator_factory;
  return &generator_factory;
}

void CompGeneratorFactory::RegisterGenerator(const std::string& model_name,
                                         std::function<CompGenerator*(void)> generator_factory_function) {
  generator_registry.emplace(model_name,generator_factory_function);
}
  
std::shared_ptr<CompGenerator> CompGeneratorFactory::CreateGenerator(const configuration::CompConfig& comp_config) const {
  CompGenerator* generator_instance = nullptr;

  // find name in the registry and call factory method
  auto gen_it = generator_registry.find(comp_config.model());
  if (gen_it != generator_registry.end()) {
    generator_instance = gen_it->second();
    generator_instance->Configure(comp_config);
  }

  if (generator_instance != nullptr) {
    return std::shared_ptr<CompGenerator>(generator_instance);
  } else {
    std::cout << "Could not find generator model with name " << comp_config.model() << std::endl;
    return nullptr;
  }
}

} // namespace generator
} // namespace cptoymc
