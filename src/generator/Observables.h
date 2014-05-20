#ifndef CPTOYMC_GENERATOR_OBSERVABLES_H
#define CPTOYMC_GENERATOR_OBSERVABLES_H

// from STL
#include <iostream>
#include <set>
#include <string>

namespace cptoymc {
namespace generator {

class Observable {
public:
  Observable(const std::string& dim_name, const std::string& var_name, const std::string& var_title);
  virtual ~Observable() { };
  
  virtual bool HasValidValue() = 0;

protected:
  std::string dim_name_;
  std::string var_name_;
  std::string var_title_;
};
  
class ObservableReal : public Observable {
public:
  ObservableReal(const std::string& dim_name, const std::string& var_name, const std::string& var_title, double value, double min_value, double max_value);
  virtual ~ObservableReal() { };
  
  void set_range(double min_value, double max_value) {
    min_value_ = min_value;
    max_value_ = max_value;
  }
  
  void set_value(double value) { value_ = value; }
  double value() const { return value_; }
  
  virtual bool HasValidValue() {
    return IsValidValue(value_);
  }
  
  bool IsValidValue(double value) {
    return (value < max_value_) && (value >= min_value_);
  }
  
private:
  double value_;
  double min_value_;
  double max_value_;
};

class ObservableInt : public Observable {
public:
  ObservableInt(const std::string& dim_name, const std::string& var_name, const std::string& var_title, int value, const std::set<int>& valid_values);
  virtual ~ObservableInt() { };
  
  virtual bool HasValidValue() {return IsValidValue(value_);}
  
  bool IsValidValue(int value) {
    return (valid_values_.find(value) != valid_values_.end());
  }
  
private:
  int value_;
  std::set<int> valid_values_;
};
  
class Observables {
public:
  Observables();
  ~Observables();

  void reset();

  double mass_true;
  double time_true;
  int    tag_true;
  double mass_meas;
  double time_meas;
  int    tag_class; // -1: SS&&!OS, 0: !OS&&!SS, 1: OS&&!SS, 10: OS&&SS
  int    tag_OS;
  double eta_OS;
  int    tag_SS;
  double eta_SS;
  int    comp_cat;

};

} // generator
} // cptoymc


#endif // CPTOYMC_GENERATOR_OBSERVABLES_H
