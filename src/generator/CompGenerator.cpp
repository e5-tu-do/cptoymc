#include "generator/CompGenerator.h"

// from STL
#include <iostream>

// from ROOT
#include "TRandom.h"

// from project - generator
#include "generator/generator.h"
#include "generator/Observables.h"

// from project - configuration
#include "configuration/CompConfig.h"



namespace cptoymc {
namespace generator {

CompGenerator::CompGenerator() :
  max_trials_(10000)
{ }


CompGenerator::~CompGenerator() { }

void CompGenerator::GenerateEvent(TRandom& rndm, Observables& observables){
  bool success = false;
  for (unsigned int i=0; i<max_trials_; ++i) {
    success = TryGenerateEvent(rndm, observables);
    if (success == true) return;
  }
  std::cout << "Failed to generate event after " << max_trials_ << " trials." << std::endl;
}

  
} // namespace generator
} // namespace cptoymc
