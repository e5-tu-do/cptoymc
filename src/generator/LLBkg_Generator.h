#ifndef CPTOYMC_GENERATOR_LLBKG_GENERATOR_H
#define CPTOYMC_GENERATOR_LLBKG_GENERATOR_H

// STL
#include <algorithm>
#include <map>
#include <memory>
#include <string>

// from project
#include "generator/CompGenerator.h"

namespace cptoymc {
namespace generator {

class LLBkg_Generator : public CompGenerator {
public:
  LLBkg_Generator();
  virtual ~LLBkg_Generator();
  
  virtual void Configure(const configuration::CompConfig& comp_config);
  virtual void GenerateEvent(TRandom& rndm, Observables& observables);
};

} // namespace generator
} // namespace cptoymc

#endif // CPTOYMC_GENERATOR_LLBKG_GENERATOR_H
