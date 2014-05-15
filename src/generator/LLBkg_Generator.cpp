#include "generator/LLBkg_Generator.h"

// from ROOT
#include "TRandom.h"

// from project - generator
#include "generator/CompGeneratorFactory.h"
#include "generator/generator.h"
#include "generator/Observables.h"

// from project - configuration
#include "configuration/CompConfig.h"

namespace cptoymc {
namespace generator {


LLBkg_Generator::LLBkg_Generator() :
  CompGenerator()
{
  
}

LLBkg_Generator::~LLBkg_Generator() {
  
}

void LLBkg_Generator::Configure(const configuration::CompConfig& comp_config) {
  auto config_ptree = comp_config.model_ptree();
  // Mass
  //auto sub_config_ptree = config_ptree.get_child("Mass");

}

void LLBkg_Generator::GenerateEvent(TRandom& rndm, Observables& observables) {

}
 

} // namespace generator
} // namespace cptoymc
