#ifndef CPTOYMC_CONFIGURATION_TOYCONFIG_H
#define CPTOYMC_CONFIGURATION_TOYCONFIG_H

// STL
#include <map>
#include <string>

// from project
#include "configuration/CompConfig.h"

namespace cptoymc {
namespace configuration {

class ToyConfig {
public:
  ToyConfig();
  ~ToyConfig();

  void load(const std::string& config_file);
  const std::map<std::string,CompConfig>& comp_configs() const {
    return comp_configs_;
  }
  
private:
  std::map<std::string,CompConfig> comp_configs_;
  //boost::property_tree::ptree obs_ptree_;
};

} // namespace configuration
} // namespace cptoymc

#endif // CPTOYMC_CONFIGURATION_TOYCONFIG_H
