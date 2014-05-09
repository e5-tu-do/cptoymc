#include "configuration/ToyConfig.h"

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

  // // parse input file
  // boost::property_tree::info_parser::read_info(config_file, pt);

  // if (!pt.empty()){ 
  //   for (auto pt_it: pt) {
  //     std::cout << pt_it.first << std::endl;
  //   }
  // } else {
  //   std::cout << "Tree is empty!" << std::endl;
  // }

  comp_configs_.emplace("Sig_Bd",CompConfig("Sig_Bd",0,10000));
  comp_configs_.emplace("Sig_Bs",CompConfig("Sig_Bs",0,100));
  comp_configs_.emplace("Bkg",   CompConfig("Bkg",0,10000));
  
  
}

} // namespace configuration
} // namespace cptoymc
