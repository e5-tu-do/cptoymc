#include "generator/Observables.h"

// from ROOT
#include <TTree.h>

// from Project
#include "configuration/ObsConfig.h"

namespace cptoymc {
namespace generator {

Observable::Observable(const std::string& dim_name, const std::string& var_name, const std::string& var_title) :
  dim_name_(dim_name),
  var_name_(var_name),
  var_title_(var_title)
  { }
  
ObservableReal::ObservableReal(const std::string& dim_name, const std::string& var_name, const std::string& var_title, double value, double min_value, double max_value) :
  Observable(dim_name, var_name, var_title),
  value_(value),
  min_value_(min_value),
  max_value_(max_value),
  var_type_("D")
{ }

ObservableInt::ObservableInt(const std::string& dim_name, const std::string& var_name, const std::string& var_title, int value, const std::set<int>& valid_values) :
  Observable(dim_name, var_name, var_title),
  value_(value),
  valid_values_(valid_values),
  var_type_("I")
  {}
  

Observables::Observables() :
  mass_true("mass_true","obsMassTrue","m'",-1000.,5000.,5500.),
  time_true("time_true","obsTimeTrue","t'",-1000.,0.,18.),
  tag_true("tag_true","obsTagTrue","d'",0,{-1,+1}),
  mass_meas("mass_meas","obsMass","m",-1000.,5000.,5500.),
  time_meas("time_meas","obsTime","t",-1000.,-2.,18.),
  tag_class("tag_class","catTag","catTag",-10000,{-1,0,+1,+10}),
  tag_OS("tag_OS","obsTagOS","d_{\\text{OS}}",0,{-1,0,+1}),
  eta_OS("eta_OS","obsEtaOS","\\eta_{\\text{OS}}",0.5,0.0,0.5),
  tag_SS("tag_SS","obsTagSS","d_{\\text{SS}}",0,{-1,0,+1}),
  eta_SS("eta_SS","obsEtaSS","\\eta_{\\text{SS}}",0.5,0.0,0.5),
  comp_cat("comp_cat","catBkg","catBkg",-10000,{0,1,10,100}),
  observables_real_(),
  observables_int_()
{
  observables_real_.emplace(mass_true.dim_name() , &mass_true );
  observables_real_.emplace(time_true.dim_name() , &time_true );
  observables_real_.emplace(mass_meas.dim_name() , &mass_meas );
  observables_real_.emplace(time_meas.dim_name() , &time_meas );
  observables_real_.emplace(eta_OS.dim_name()    , &eta_OS    );
  observables_real_.emplace(eta_SS.dim_name()    , &eta_SS    );

  observables_int_.emplace(tag_true.dim_name()  , &tag_true  );
  observables_int_.emplace(tag_class.dim_name() , &tag_class );
  observables_int_.emplace(tag_OS.dim_name()    , &tag_OS    );
  observables_int_.emplace(tag_SS.dim_name()    , &tag_SS    );
  observables_int_.emplace(comp_cat.dim_name()  , &comp_cat  );
}

void Observables::Configure(const std::shared_ptr<configuration::ObsConfig> obs_config) {
  // loop over real observables
  for ( auto obs_real_config_entry : obs_config->obs_configs_real()) {
    auto obs_real_entry = observables_real_.find(obs_real_config_entry.first);
    if (obs_real_entry != observables_real_.end()) {
      auto obs_real = obs_real_entry->second;
      auto obs_real_config = obs_real_config_entry.second;
      obs_real->set_var_name(std::get<0>(obs_real_config));
      obs_real->set_var_title(std::get<1>(obs_real_config));
      obs_real->set_range(std::get<2>(obs_real_config), std::get<3>(obs_real_config));
    } else {
      std::cout << "No observable called " << obs_real_config_entry.first << " known to generator." << std::endl;
    }
  }
  // loop over int obervables
  for ( auto obs_int_config_entry : obs_config->obs_configs_int()) {
    auto obs_int_entry = observables_int_.find(obs_int_config_entry.first);
    if (obs_int_entry != observables_int_.end()) {
      auto obs_int = obs_int_entry->second;
      auto obs_int_config = obs_int_config_entry.second;
      obs_int->set_var_name(std::get<0>(obs_int_config));
      obs_int->set_var_title(std::get<1>(obs_int_config));
      obs_int->set_valid_values(std::get<2>(obs_int_config));
    } else {
      std::cout << "No observable called " << obs_int_config_entry.first << " known to generator." << std::endl;
    }
  }
}

void Observables::Reset() {
  mass_true.set_value(-1000.);
  time_true.set_value(-1000.);
  tag_true.set_value(0);
  mass_meas.set_value(-1000.);
  time_meas.set_value(-1000.);
  tag_class.set_value(-10000);
  tag_OS.set_value(0);
  eta_OS.set_value(0.5);
  tag_SS.set_value(0);
  eta_SS.set_value(0.5);
  comp_cat.set_value(-10000);
}

void Observables::RegisterObservableBranches(TTree& out_tree) {
  for (auto& obs_pair : observables_real_) {
    out_tree.Branch(obs_pair.second->var_name().c_str(),
                    &(obs_pair.second->value_),
                    (obs_pair.second->var_name()+"/"+obs_pair.second->var_type()).c_str()
                    );
  }
  for (auto& obs_pair : observables_int_) {
    out_tree.Branch(obs_pair.second->var_name().c_str(),
                    &(obs_pair.second->value_),
                    (obs_pair.second->var_name()+"/"+obs_pair.second->var_type()).c_str()
                    );
  }

}
  
  
} // namespace generator
} // namespace cptoymc

