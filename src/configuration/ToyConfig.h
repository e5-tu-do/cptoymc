#ifndef CPTOYMC_CONFIGURATION_TOYCONFIG_H
#define CPTOYMC_CONFIGURATION_TOYCONFIG_H

// STL
#include <string>


namespace cptoymc {
namespace configuration {

class ToyConfig {
public:
  ToyConfig();
  ~ToyConfig();

  void load(const std::string& config_file);

private:


};

} // namespace configuration
} // namespace cptoymc

#endif // CPTOYMC_CONFIGURATION_TOYCONFIG_H
