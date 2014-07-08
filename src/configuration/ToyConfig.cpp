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

  if (pt.empty()) {
    std::cout << "Tree is empty! Cannot configure the toy generator. Aborting." << std::endl;
    return;
  }

  auto pt_comp = pt.get_child("Components");
  for (auto pt_comp_it : pt_comp) {
    comp_configs_.emplace(
      pt_comp_it.first,
      CompConfig(pt_comp_it.second.get<std::string>("name"),
                 pt_comp_it.second.get("comp_cat",-1000),
                 pt_comp_it.second.get("yield",0),
                 pt_comp_it.second.get<std::string>("model","NoModel"),
                 pt_comp_it.second.get_child("model"))
                 );
  }


}

} // namespace configuration
} // namespace cptoymc
