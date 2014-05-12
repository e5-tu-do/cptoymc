#ifndef CPTOYMC_CONFIGURATION_COMPCONFIG_H
#define CPTOYMC_CONFIGURATION_COMPCONFIG_H 

// STL
#include <string>

namespace cptoymc {
namespace configuration {

class CompConfig {
public:
  CompConfig(const std::string& name, int comp_cat, int yield, const std::string& model);
  ~CompConfig();

  const std::string& name() const { return name_; }
  int comp_cat() const { return comp_cat_; }
  int yield() const { return yield_; }
  const std::string& model() const { return model_; }
  
private:
  const std::string name_;
  const int comp_cat_;
  const int yield_;
  std::string model_;

};

} // namespace configuration
} // namespace cptoymc

#endif // CPTOYMC_CONFIGURATION_COMPCONFIG_H
