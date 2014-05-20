#include "generator/Observables.h"

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
  max_value_(max_value)
{ }

ObservableInt::ObservableInt(const std::string& dim_name, const std::string& var_name, const std::string& var_title, int value, const std::set<int>& valid_values) :
  Observable(dim_name, var_name, var_title),
  value_(value),
  valid_values_(valid_values)
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
  comp_cat("comp_cat","catBkg","catBkg",-10000,{0,1,10,100})
{

}

Observables::~Observables() {
  
}


void Observables::reset() {
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

} // namespace generator
} // namespace cptoymc

