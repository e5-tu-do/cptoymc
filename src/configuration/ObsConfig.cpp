#include "configuration/ObsConfig.h"

namespace cptoymc {
namespace configuration {
  
ObsConfig::ObsConfig(const boost::property_tree::ptree& obs_ptree) :
  obs_ptree_(obs_ptree),
  obs_configs_real_(),
  obs_configs_int_()
{
  std::string obs_internal_name;
  std::string obs_name;
  std::string obs_title;
  std::string obs_type;
  
  // for real vars
  double min, max = -10000.;
  
  // for int vars
  std::set<int> allowed_vals;
  
  for (auto pt_obs_it : obs_ptree) {
    obs_internal_name = pt_obs_it.first;
    obs_name  = pt_obs_it.second.get("name",obs_internal_name);
    obs_title = pt_obs_it.second.get("title","No Title Given");
    obs_type  = pt_obs_it.second.get<std::string>("type");

    std::cout << "   Adding Observables " << obs_internal_name
    << " with name " << obs_name;
    if ( obs_type == "Integer") {
      allowed_vals.clear();
      for (auto val : pt_obs_it.second.get_child("range")) {
        allowed_vals.insert(val.second.get<int>(""));
      }
      obs_configs_int_.emplace(obs_internal_name,ObsConfInt(obs_name,obs_title,allowed_vals));
      std::cout << " and allowed values { ";
      for (auto allowed_val : allowed_vals) {
        std::cout << allowed_val << " ";
      }
      std::cout << "}."<< std::endl;;
    }
    else if ( obs_type == "Real" ) {
      min = pt_obs_it.second.get<double>("min");
      max = pt_obs_it.second.get<double>("max");
      obs_configs_real_.emplace(obs_internal_name,ObsConfReal(obs_name,obs_title,min,max));
      std::cout << " and allowed range [ " << min << ", " << max << " ]." << std::endl;
    }
    
  }
}
  
ObsConfig::~ObsConfig() { }

} // namespace configuration
} // namespace cptoymc