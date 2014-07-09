#ifndef CPTOYMC_GENERATOR_COMPGENERATORFACTORY_H
#define CPTOYMC_GENERATOR_COMPGENERATORFACTORY_H

// from STL
#include <algorithm>
#include <map>
#include <memory>
#include <string>
#include <iostream>

namespace cptoymc {
  
// forward declarations
namespace configuration {
  class CompConfig;
}

namespace generator {

class CompGenerator;
  

class CompGeneratorFactory {
public:
  ~CompGeneratorFactory();

  static CompGeneratorFactory* Instance();
  void RegisterGenerator(const std::string& model_name,
                         std::function<CompGenerator*(void)> generator_factory_function);
  
  std::shared_ptr<CompGenerator> CreateGenerator(const configuration::CompConfig& comp_config) const;


private:
  CompGeneratorFactory();
  std::map<std::string,std::function<CompGenerator*(void)>> generator_registry;
  
  template<class T>
  void RegisterGenerator(const std::string& model_name) {
    RegisterGenerator(model_name, [](void) -> CompGenerator* { return new T();});
  }
  
};

  

} // namespace generator
} // namespace cptoymc


#endif // CPTOYMC_GENERATOR_COMPGENERATORFACTORY_H


