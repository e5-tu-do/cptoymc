#ifndef CPTOYMC_CONFIGURATION_H
#define CPTOYMC_CONFIGURATION_H

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

} // namespace configuration
} // namespace cptoymc

#endif // CPTOYMC_CONFIGURATION_H
