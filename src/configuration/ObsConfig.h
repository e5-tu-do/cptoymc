#ifndef CPTOYMC_CONFIGURATION_OBSCONFIG_H
#define CPTOYMC_CONFIGURATION_OBSCONFIG_H

// STL
#include <map>
#include <set>
#include <string>
#include <tuple>

// from Boost - property_tree
#include <boost/property_tree/ptree.hpp>

namespace cptoymc {
namespace configuration {

class ObsConfig {
 public:
  ObsConfig(const boost::property_tree::ptree& obs_ptree);
  ~ObsConfig();
  
  const boost::property_tree::ptree model_ptree() const { return obs_ptree_; }
    
 private:
  const boost::property_tree::ptree obs_ptree_;
  
  typedef std::tuple<std::string,std::string,double,double> ObsConfReal;
  typedef std::tuple<std::string,std::string,std::set<int>> ObsConfInt;
  
  std::map<std::string,ObsConfReal> obs_configs_real_;
  std::map<std::string,ObsConfInt>  obs_configs_int_;
};
    
} // namespace configuration
} // namespace cptoymc

#endif // CPTOYMC_CONFIGURATION_OBSCONFIG_H
