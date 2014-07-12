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

ToyConfig::ToyConfig() :
  comp_configs_(),
  obs_config_()
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
  
  // Configuration of observables
  auto pt_obs = pt.get_child("Observables");
  obs_config_ = std::unique_ptr<ObsConfig>(new ObsConfig(pt_obs));
  
  // Configuriation of components
  std::cout << "Configuring the generators for the different components." << std::endl;
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
    std::cout << "   Added component " << pt_comp_it.first << std::endl;
  }


}

} // namespace configuration
} // namespace cptoymc
