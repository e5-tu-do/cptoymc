#include "configuration/CompConfig.h"


namespace cptoymc {
namespace configuration {

  
  
CompConfig::CompConfig(const std::string& name, int comp_cat, int yield,
                       const std::string& model,
                       const boost::property_tree::ptree& model_ptree) :
  name_(name),
  comp_cat_(comp_cat),
  yield_(yield),
  model_(model),
  model_ptree_(model_ptree)
{

}

CompConfig::~CompConfig() {

}

  
} // namespace configuration
} // namespace cptoymc
