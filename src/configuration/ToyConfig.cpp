#include "configuration/ToyConfig.h"

// from STDL
#include <iostream>


// from Boost - property_tree
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/info_parser.hpp>

// from project
#include "configuration/CompConfig.h"

namespace cptoymc {
namespace configuration {
  
ToyConfig::ToyConfig() 
{

}

ToyConfig::~ToyConfig() {

}

void ToyConfig::load(const std::string& config_file) {
  using boost::property_tree::ptree;
  ptree pt;

  // parse input file
  boost::property_tree::info_parser::read_info(config_file, pt);

  if (!pt.empty()){ 
    for (auto pt_it: pt) {
      // Prepare Components
      
    }
  } else {
    std::cout << "Tree is empty!" << std::endl;
  }

  // get component_trees
  auto pt_comp = pt.get_child("Components");
  for (auto pt_comp_it: pt_comp) {
    comp_configs_.emplace(
      pt_comp_it.second.get<std::string>("name"),
      CompConfig(pt_comp_it.second.get<std::string>("name"),
                 pt_comp_it.second.get("comp_cat",-1000),
                 pt_comp_it.second.get("yield",0))
    );
  }

  
  
}

} // namespace configuration
} // namespace cptoymc
