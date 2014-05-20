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

  const std::string& dim_name() const {return dim_name_;}
  const std::string& var_name() const {return var_name_;}
  const std::string& var_title() const {return var_title_;}
  
protected:
  std::string dim_name_;
  std::string var_name_;
  std::string var_title_;
};
  
class ObservableReal : public Observable {
public:
  ObservableReal(const std::string& dim_name, const std::string& var_name, const std::string& var_title, double value, double min_value, double max_value);
  virtual ~ObservableReal() { };
  
  void set_value(double value) { value_ = value; }
  double value() const { return value_; }
  
  void set_range(double min_value, double max_value) {
    min_value_ = min_value;
    max_value_ = max_value;
  }
  
  virtual bool HasValidValue() {
    return IsValidValue(value_);
  }
  
  bool IsValidValue(double value) {
    return (value < max_value_) && (value >= min_value_);
  }
  double value_;
 
private:
  double min_value_;
  double max_value_;
};

class ObservableInt : public Observable {
public:
  ObservableInt(const std::string& dim_name, const std::string& var_name, const std::string& var_title, int value, const std::set<int>& valid_values);
  virtual ~ObservableInt() { };
  
  void set_value(int value) { value_ = value; }
  int value() const { return value_; }
  void set_valid_values(const std::set<int> valid_values) {valid_values_ = valid_values;}
  
  virtual bool HasValidValue() {return IsValidValue(value_);}
  
  bool IsValidValue(int value) {
    return (valid_values_.find(value) != valid_values_.end());
  }

  int value_;
  
private:
  std::set<int> valid_values_;
};
  
class Observables {
public:
  Observables();
  ~Observables() {};

  void reset();

  ObservableReal mass_true;
  ObservableReal time_true;
  ObservableInt  tag_true;
  ObservableReal mass_meas;
  ObservableReal time_meas;
  ObservableInt  tag_class; // -1: SS&&!OS, 0: !OS&&!SS, 1: OS&&!SS, 10: OS&&SS
  ObservableInt  tag_OS;
  ObservableReal eta_OS;
  ObservableInt  tag_SS;
  ObservableReal eta_SS;
  ObservableInt  comp_cat;

};

} // generator
} // cptoymc


#endif // CPTOYMC_GENERATOR_OBSERVABLES_H
