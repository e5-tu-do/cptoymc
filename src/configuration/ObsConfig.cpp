#include "configuration/ObsConfig.h"

namespace cptoymc {
namespace configuration {
  
ObsConfig::ObsConfig(const boost::property_tree::ptree& obs_ptree) :
  obs_ptree_(obs_ptree),
  obs_configs_real_(),
  obs_configs_int_()
{
  std::string obs_type = "None";
  double min, max = -10000.;
  std::set<int> allowed_vals;
  
  for (auto pt_obs_it : obs_ptree) {
    obs_type = pt_obs_it.second.get<std::string>("type");
    if ( obs_type == "Integer") {
      
      
//      obs_configs_real_.emplace(
//        pt_obs_it.second.get<std::string>("name"),
//        ObsConfInt{
//          pt_obs_it.second.get<std::string>("name"),
//          pt_obs_it.second.get("title","No title given"),
//          pt_obs_it.second.get<double>(
//        }
//      );
    }
//    comp_configs_.emplace(
//                          pt_comp_it.first,
//                          CompConfig(pt_comp_it.second.get<std::string>("name"),
//                                     pt_comp_it.second.get("comp_cat",-1000),
//                                     pt_comp_it.second.get("yield",0),
//                                     pt_comp_it.second.get<std::string>("model","NoModel"),
//                                     pt_comp_it.second.get_child("model"))
//                          );
    std::cout << "   Added component " << pt_obs_it.first << std::endl;
  }
}
  
ObsConfig::~ObsConfig() { }

} // namespace configuration
} // namespace cptoymc