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
  ParsMass        pars_sigBd_mass;
  ParsTimeAndCP   pars_sigBd_timeandcp;
  ParsMassResol   pars_sigBd_mass_resol;  
  ParsTimeResol   pars_sigBd_time_resol;
  ParsTagging     pars_sigBd_tagging;

  ParsMass        pars_sigBs_mass;
  ParsTimeAndCP   pars_sigBs_timeandcp;
  ParsMassResol   pars_sigBs_mass_resol;  
  ParsTimeResol   pars_sigBs_time_resol;
  ParsTagging     pars_sigBs_tagging;

  ParsMass        pars_bkg_mass;
  ParsTimeAndCP   pars_bkg_timeandcp;
  ParsMassResol   pars_bkg_mass_resol;  
  ParsTimeResol   pars_bkg_time_resol;
  ParsTagging     pars_bkg_tagging;

}; // class ToyConfig

// class ComponentConfig {

// }

} // namespace configuration
} // namespace cptoymc

#endif // CPTOYMC_CONFIGURATION_H
