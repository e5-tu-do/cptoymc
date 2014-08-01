#ifndef CPTOYMC_GENERATOR_OBSERVABLES_H
#define CPTOYMC_GENERATOR_OBSERVABLES_H

// from STL
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <string>

// forward declarations
class TTree;

namespace cptoymc {
namespace configuration{
// forward declaration
class ObsConfig;
} // namespace configuration
  
namespace generator {

class Observable {
public:
  Observable(const std::string& dim_name, const std::string& var_name, const std::string& var_title);
  virtual ~Observable() { };
  
  virtual bool HasValidValue() = 0;

  const std::string& dim_name() const {return dim_name_;}
  const std::string& var_name() const {return var_name_;}
  const std::string& var_title() const {return var_title_;}
  
  void set_var_name(const std::string& var_name)   { var_name_  = var_name;}
  void set_var_title(const std::string& var_title) { var_title_ = var_title;}
  
  virtual const std::string& var_type() = 0;
  
protected:
  const std::string dim_name_;
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
 
  virtual const std::string& var_type() {
    return var_type_;
  }
  
  double max_value() const {return max_value_;}
  double min_value() const {return min_value_;}
  
private:
  double min_value_;
  double max_value_;
  const std::string var_type_;
};

class ObservableInt : public Observable {
public:
  ObservableInt(const std::string& dim_name, const std::string& var_name,
                const std::string& var_title, int value, const std::map<std::string, int>& types);
  
  virtual ~ObservableInt() { };
  
  void set_value(int value) { value_ = value; }
  int value() const { return value_; }
  
  //void set_valid_values(const std::set<int> valid_values) {valid_values_ = valid_values;}
  
  void set_valid_types_values(const std::map<std::string,int>& types) { types_ = types; }
  
  virtual bool HasValidValue() {return IsValidValue(value_);}
  
  bool IsValidValue(int value) {
    return (valid_values_.find(value) != valid_values_.end());
  }
  
  bool IsValidKey(const std::string& type_name) {
    return (types_.find(type_name) != types_.end());
  }
  
  int GetValueForType(const std::string& type_name) {
    auto type_value_pair = types_.find(type_name);
    if (type_value_pair != types_.end()) {
      return type_value_pair->second;
    }
    else {
      std::cout << "Cannot find category type " << type_name << " in Category " << dim_name_ << "." << std::endl;
      return -10000;
    }
  }
  
  virtual const std::string& var_type() {
    return var_type_;
  }
  
  int value_;
  
private:
  std::map<std::string,int> types_;
  std::set<int> valid_values_;
  const std::string var_type_;
};
  
class Observables {
public:
  Observables();
  ~Observables() { }

  void Configure(const std::shared_ptr<configuration::ObsConfig> obs_config);
  void Reset();
  void RegisterObservableBranches(TTree& out_tree);

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

private:
  std::map<std::string,ObservableReal*> observables_real_;
  std::map<std::string,ObservableInt*> observables_int_;
  
};

} // generator
} // cptoymc


#endif // CPTOYMC_GENERATOR_OBSERVABLES_H
