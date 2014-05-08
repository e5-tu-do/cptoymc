#ifndef CPTOYMC_CONFIGURATION_H
#define CPTOYMC_CONFIGURATION_H

// from STL
#include <string>

namespace cptoymc {
namespace configuration {



struct ParsMass {
  double mean;
  double gamma;
};

struct ParsTimeAndCP {
  double tau;
  double dGamma;
  double dm;
  double Sf;
  double Cf;
  double Df;
  double AP;
};

struct ParsMassResol  {
  double bias;
  double sigma;
};

struct ParsTimeResol {
  double bias;
  double sigma;
};

struct ParsTagging {
  double eff_OS;
  double eff_SS;
  double eff_SSOS;
  double w_OS;
  double dw_OS;
  double w_SS;
  double dw_SS;
};

class ToyConfig {
public:
  ToyConfig();
  ~ToyConfig();

  void load(const std::string& config_file);

private:


};

class CompConfig {
public:
  CompConfig(const std::string& name, int bkg_cat, int yield);
  ~CompConfig();

  const std::string& name() { return name_; }
  int bkg_cat() { return bkg_cat_; }
  int yield() { return yield_; }

private:
  const std::string name_;
  const int bkg_cat_;
  const int yield_;

};

} // namespace configuration
} // namespace cptoymc

#endif // CPTOYMC_CONFIGURATION_H
