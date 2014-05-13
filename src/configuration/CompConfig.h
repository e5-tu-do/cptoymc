#ifndef CPTOYMC_CONFIGURATION_COMPCONFIG_H
#define CPTOYMC_CONFIGURATION_COMPCONFIG_H 

// STL
#include <string>

// from Boost - property_tree
#include <boost/property_tree/ptree.hpp>

namespace cptoymc {
namespace configuration {

class CompConfig {
public:
  CompConfig(const std::string& name, int comp_cat, int yield,
             const std::string& model,
             const boost::property_tree::ptree& model_ptree);
  ~CompConfig();

  const std::string& name() const { return name_; }
  int comp_cat() const { return comp_cat_; }
  int yield() const { return yield_; }
  const std::string& model() const { return model_; }
  const boost::property_tree::ptree model_ptree() const { return model_ptree_; }
  
private:
  const std::string name_;
  const int comp_cat_;
  const int yield_;
  std::string model_;
  const boost::property_tree::ptree model_ptree_;

};

} // namespace configuration
} // namespace cptoymc

#endif // CPTOYMC_CONFIGURATION_COMPCONFIG_H
