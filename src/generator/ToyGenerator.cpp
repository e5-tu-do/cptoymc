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


ToyGenerator::ToyGenerator(const cptoymc::configuration::ToyConfig& config, unsigned int seed) :
  config_(config),
  obs_(),
  comp_generators_(),
  seed_(seed),
  rndm_(new TRandom3(seed_))
{
  // initialize observables
  obs_.Configure(config_.obs_config());
  
  // initialize component generators
  for ( auto comp_config : config_.comp_configs() ) {
    comp_generators_.emplace(comp_config.first,CompGeneratorFactory::Instance()->CreateGenerator(comp_config.second));
  }

}
  
  
ToyGenerator::~ToyGenerator() {
  delete rndm_;
}

void ToyGenerator::GenerateToy(TTree& out_tree) {
  using configuration::CompConfig;
  
  // Prepare Tree
  obs_.Reset();
  obs_.RegisterObservableBranches(out_tree);
  
  
  // loop over components, get their yield
  std::string comp_name  = "";
  int comp_cat           = -1000;
  int exp_events_of_comp = 0;
  int num_events_of_comp = 0;
  int num_events_total   = 0;
  
  for ( auto comp_config : config_.comp_configs() ) {
    comp_name = comp_config.first;
    comp_cat  = comp_config.second.comp_cat();
    exp_events_of_comp = comp_config.second.exp_yield();
    
    //num_events_of_comp = exp_events_of_comp;
    num_events_of_comp = generator::yieldToGenerate(*rndm_, exp_events_of_comp);
    num_events_total += num_events_of_comp;
    std::cout << "Generating " << num_events_of_comp << " events for component "
              << comp_name << " (expected yield: " << exp_events_of_comp << ")"
              << std::endl;
    
    auto comp_generator = comp_generators_.find(comp_name);
    if (comp_generator != comp_generators_.end()) {
      for (int i=0; i < num_events_of_comp; ++i) {
        obs_.Reset();
        comp_generator->second->GenerateEvent(*rndm_, obs_);
        out_tree.Fill();
      }
    } else {
      std:: cout << "Could not find component generator for component " << comp_name << std::endl;
    }
    
  }
  
}
  
void ToyGenerator::GenerateToy(TTree& out_tree, unsigned int seed) {
  rndm_->SetSeed(seed);
  GenerateToy(out_tree);
}
  
} // namespace cptoymc
} // namespace generator

