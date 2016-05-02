#include "configuration/ObsConfig.h"

#include <iostream>

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
  std::map<std::string,int> allowed_type_vals;
  
  for (auto pt_obs_it : obs_ptree) {
    obs_internal_name = pt_obs_it.first;
    obs_name  = pt_obs_it.second.get("name",obs_internal_name);
    obs_title = pt_obs_it.second.get("title","No Title Given");
    obs_type  = pt_obs_it.second.get<std::string>("type");

    std::cout << "   Preparing observable " << obs_internal_name
    << " with name " << obs_name;
    if ( obs_type == "Integer") {
      allowed_type_vals.clear();
      for (auto val : pt_obs_it.second.get_child("range")) {
        allowed_type_vals.emplace(val.first, val.second.get<int>(""));
      }
      obs_configs_int_.emplace(obs_internal_name,ObsConfInt(obs_name,obs_title,allowed_type_vals));
      std::cout << " and allowed type/value pairs { ";
      for (auto allowed_type_val : allowed_type_vals) {
        std::cout << "{" << allowed_type_val.first << " : " << allowed_type_val.second << "} ";
      }
      std::cout << "}."<< std::endl;;
    }
    else if ( obs_type == "Real" ) {
      min = pt_obs_it.second.get<double>("min");
      max = pt_obs_it.second.get<double>("max");
      obs_configs_real_.emplace(obs_internal_name,ObsConfReal(obs_name,obs_title,min,max));
      std::cout << " and allowed range [ " << min << ", " << max << " ]." << std::endl;
    }
    else {
      std::cout << " with unknown type " << obs_type << "! "
                << "Cannot add this observable!" << std::endl;
    }
    
  }
}
  
ObsConfig::~ObsConfig() { }

  
  
} // namespace configuration
} // namespace cptoymc