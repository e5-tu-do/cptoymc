#include "configuration/configuration.h"

// from STL




namespace cptoymc {
namespace configuration {

CompConfig::CompConfig(const std::string& name, int bkg_cat, int yield) :
  name_(name),
  bkg_cat_(bkg_cat),
  yield_(yield)
{

}

CompConfig::~CompConfig() {

}

  
} // namespace configuration
} // namespace cptoymc
