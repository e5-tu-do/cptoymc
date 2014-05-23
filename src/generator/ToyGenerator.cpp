#include "generator/ToyGenerator.h"

// from STL
#include <iostream>

// from project
#include "configuration/ToyConfig.h"
#include "configuration/CompConfig.h"
#include "generator/generator.h"
#include "generator/Observables.h"

// from project - component generators
#include "generator/CompGenerator.h"
#include "generator/CompGeneratorFactory.h"
#include "generator/BSig_CPV_P2VP_Generator.h"


// from ROOT
#include "TRandom3.h"
#include "TTree.h"

namespace cptoymc {
namespace generator {


ToyGenerator::ToyGenerator(const cptoymc::configuration::ToyConfig& config) :
  config_(config)
{
  
}
  
  
ToyGenerator::~ToyGenerator() {
  
}

void ToyGenerator::GenerateToy(TTree& out_tree) {
  using configuration::CompConfig;
  //CompGeneratorRegistrar<BSig_CPV_P2VP_Generator> registrar("BSig_CPV_P2VP");

  TRandom3 rndm;
  
  std::string tree_name = "ToyMC";
  std::string tree_desc = "ToyMC Tree";
  
  // Prepare Observables
  Observables obs;
  
  
  
  // Prepare Tree
  obs.registerObservableBranches(out_tree);
  
  
  // loop over components, get their yield
  std::string comp_name  = "";
  int comp_cat           = -1000;
  int exp_events_of_comp = 0;
  int num_events_of_comp = 0;
  int num_events_total   = 0;
  
  
  
  for ( auto comp_config : config_.comp_configs() ) {
    comp_name = comp_config.first;
    comp_cat  = comp_config.second.comp_cat();
    exp_events_of_comp = comp_config.second.yield();
    
    num_events_of_comp = generator::yieldToGenerate(rndm, exp_events_of_comp);
    num_events_total += num_events_of_comp;
    std::cout << "Generating " << num_events_of_comp << " events for component "
              << comp_name << " (expected yield: " << exp_events_of_comp << ")"
              << std::endl;
    auto comp_generator = CompGeneratorFactory::Instance()->CreateGenerator(comp_config.second);
    for (int i=0; i < num_events_of_comp; ++i) {
      obs.reset();
      comp_generator->GenerateEvent(rndm, obs);
      out_tree.Fill();
    }
  }
  
}
  
} // namespace cptoymc
} // namespace generator

