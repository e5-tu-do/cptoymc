#include "configuration/CompConfig.h"


namespace cptoymc {
namespace configuration {

  
  
CompConfig::CompConfig(const std::string& name, int comp_cat, double exp_yield,
                       const std::string& model,
                       const boost::property_tree::ptree& model_ptree) :
  name_(name),
  comp_cat_(comp_cat),
  exp_yield_(exp_yield),
  model_(model),
  model_ptree_(model_ptree)
{

}

CompConfig::~CompConfig() {

}

  
} // namespace configuration
} // namespace cptoymc
