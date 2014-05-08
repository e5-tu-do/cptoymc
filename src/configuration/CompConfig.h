#ifndef CPTOYMC_CONFIGURATION_COMPCONFIG_H
#define CPTOYMC_CONFIGURATION_COMPCONFIG_H 

// STL
#include <string>

namespace cptoymc {
namespace configuration {

class CompConfig {
public:
  CompConfig(const std::string& name, int bkg_cat, int yield);
  ~CompConfig();

  const std::string& name() { return name_; }
  int bkg_cat() { return bkg_cat_; }
  int yield() { return yield_; }

private:
  const std::string name_;
  const int bkg_cat_;
  const int yield_;

};

} // namespace configuration
} // namespace cptoymc

#endif // CPTOYMC_CONFIGURATION_COMPCONFIG_H
