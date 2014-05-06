#include "Configuration.h"

// from STL

// from Boost - property_tree
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/info_parser.hpp>


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


  
}

  
} // namespace configuration
} // namespace cptoymc
