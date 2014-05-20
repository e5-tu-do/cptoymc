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
    max_tries_(10000)
  { }


CompGenerator::~CompGenerator() { }





  
} // namespace generator
} // namespace cptoymc
