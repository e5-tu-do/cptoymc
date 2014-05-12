#include "generator/CompGenerator.h"

// from ROOT
#include "TRandom.h"

// from project - configuration
#include "configuration/CompConfig.h"

// from project - generator


namespace cptoymc {
namespace generator {


CompGenerator::CompGenerator(const configuration::CompConfig& comp_config) :
  comp_config_(comp_config)
{ }

CompGenerator::~CompGenerator() {
}

BSig_CPV_P2VP_Generator::BSig_CPV_P2VP_Generator(const configuration::CompConfig& comp_config) :
  CompGenerator(comp_config),
  params_mass_{5279.15, 0.},
  params_massresol_{0.,8.},
  params_timeandcp_{1.5,0.,0.5,0.7,0.,0.7,0.},
  params_timeresol_{0.,0.05},
  params_taggingeffs_{0.30,0.06,0.04},
  params_taggingOS_{1.0,0.0,0.0,0.0,0.0},
  params_taggingSS_{1.0,0.0,0.0,0.0,0.0}
{
  
}

BSig_CPV_P2VP_Generator::~BSig_CPV_P2VP_Generator() {
  
}

void BSig_CPV_P2VP_Generator::generateEvent(TRandom& rndm, Observables& observables) {
  
}
  

  
} // namespace generator
} // namespace cptoymc
