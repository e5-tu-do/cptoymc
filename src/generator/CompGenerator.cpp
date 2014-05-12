#include "generator/CompGenerator.h"

// from ROOT
#include "TRandom.h"

// from project - configuration
#include "configuration/CompConfig.h"

// from project - generator


namespace cptoymc {
namespace generator {


CompGenerator::CompGenerator(const configuration::CompConfig& comp_config ) :
  comp_config_(comp_config)
{ }

CompGenerator::~CompGenerator() {
}

BSig_CPV_P2VP_Generator::BSig_CPV_P2VP_Generator(const configuration::CompConfig& comp_config) :
CompGenerator(comp_config) {
  
}

BSig_CPV_P2VP_Generator::~BSig_CPV_P2VP_Generator() {
  
}

void BSig_CPV_P2VP_Generator::generateEvent(TRandom& rndm, Observables& observables) {
  
}
  
} // namespace generator
} // namespace cptoymc
