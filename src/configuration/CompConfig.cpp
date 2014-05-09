#include "configuration/CompConfig.h"

namespace cptoymc {
namespace configuration {

CompConfig::CompConfig(const std::string& name, int comp_cat, int yield) :
  name_(name),
  comp_cat_(comp_cat),
  yield_(yield)
{
  
}

CompConfig::~CompConfig() {

}

  
} // namespace configuration
} // namespace cptoymc
