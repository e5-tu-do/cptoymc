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
  mass_true(-1000.),
  time_true(-1000.),
  tag_true(0),
  mass_meas(-1000.),
  time_meas(-1000.),
  tag_class(-10000),
  tag_OS(0),
  eta_OS(0.5),
  tag_SS(0),
  eta_SS(0.5),
  comp_cat(-1000)
{

}

Observables::~Observables() {
  
}


void Observables::reset() {
  mass_true = -1000.;
  time_true = -1000.;
  tag_true  = 0;
  mass_meas = -1000.;
  time_meas = -1000.;
  tag_class = -10000;
  tag_OS    = 0;
  eta_OS    = 0.5;
  tag_SS    = 0;
  eta_SS    = 0.5;
  comp_cat  = -1000;
}

} // namespace generator
} // namespace cptoymc

