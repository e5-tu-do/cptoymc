#include "generator/generator.h"

// from STL
#include <cmath>
#include <iostream>

// from ROOT
#include "TRandom.h"

// from project
#include "configuration/configuration.h"


namespace cptoymc {
namespace generator {

bool GenerateExpo(TRandom& rndm, double par_expo, double& obs, double min, double max) {
  if (par_expo != 0.) {
    obs = (-1./par_expo)*log(exp(-1.*max*par_expo)+rndm.Uniform()*(exp(-1.*min*par_expo)-exp(-1.*max*par_expo)));
  } else {
    obs = rndm.Uniform(min,max);
  }
  return true;
}


bool GenerateMassBreitWigner(TRandom& rndm, double par_mean, double par_gamma, double& obs_mass_true) {
  obs_mass_true = rndm.BreitWigner(par_mean, par_gamma);
  return true;
}

bool GenerateLognormal(TRandom& rndm, double m, double k, double min, double max,
                       double& obs_sigma_t) {
  unsigned int num_samples(0);

  if (min == max) {
    obs_sigma_t = min;
    return true;
  } else {
    obs_sigma_t = std::exp(std::log(m) + std::log(k)*rndm.Gaus(0,1));
    ++num_samples;

    while (obs_sigma_t > max || obs_sigma_t < min) {
      obs_sigma_t = std::exp(std::log(m) + std::log(k)*rndm.Gaus(0,1));
      ++num_samples;

      if (num_samples % 10000000 == 0) {
        std::cout << "WARNING in GenerateLognormal(rndm, m=" << m << ", k=" << k << ", min=" << min << ", max=" << max << "): Generated " << num_samples << " sample values without one candidate passing. You probably want to check your parameters." << std::endl;
      }
    }

    return true;
  }
}

bool GenerateCPV_P2PV(TRandom& rndm, double par_prod_asym, double par_det_asym,
                      double par_tau, double par_dGamma, double par_dm,
                      double par_Sf, double par_Cf, double par_Df,
                      double par_Sfbar, double par_Cfbar, double par_Dfbar,
                      double& obs_time_true, int& obs_tag_true, int& finalstate,
                      int& mix_true) {
  // helper quantities
  double prob_B = (1. - par_prod_asym)/2.;
  double gamma_min = 1./par_tau - std::abs(par_dGamma)/2.;

  // local values of tag and time
  int    val_d = 0;
  int    val_final = 1;
  double val_t = 0.;

  // now hit and miss, stolen from BDecay!!!
  double val_envelope = 0.;
  double val_pdf      = 0.;
  double val_bar_pdf  = 0.;
  while(true) {
    // get initial flavour while taking production asymmetry into account
    val_d = (rndm.Uniform() < prob_B) ? +1 : -1;

    // get time
    val_t = -log(rndm.Uniform())/gamma_min;

    // get pdf value and envelope value
    // the bar CPV params are on -2.0 on default; this makes noe sense for both finalstates, thus generate only for one finalstate
    if(par_Sfbar == -2 || par_Cfbar == -2 || par_Dfbar == -2){
      val_pdf      = BCPV_PDF(val_t, val_d, par_tau, par_dGamma, par_dm, par_Sf, par_Cf, par_Df);
      val_envelope = BCPV_PDF_Envelope(val_t, gamma_min, par_Sf, par_Cf, par_Df);

      if (val_envelope < val_pdf) {
        std::cout << "WARNING: Envelope smaller than PDF!" << std::endl;
      }

      if(val_envelope*rndm.Uniform() > val_pdf) continue;
      else break;
    }
    else{
      double prob_detect = (1 - par_det_asym)/2;

      val_final = (rndm.Uniform() < prob_detect) ? +1 : -1;

      if(val_final  == 1) {
        val_pdf = BCPV_PDF(val_t, val_d, par_tau, par_dGamma, par_dm, par_Sf, par_Cf, par_Df);
        val_envelope = BCPV_PDF_Envelope(val_t, gamma_min, par_Sf, par_Cf, par_Df);
        if (val_envelope < val_pdf) std::cout << "WARNING: Envelope smaller than PDF!" << std::endl;
        if(val_envelope*rndm.Uniform() > val_pdf) continue;
        else break;
      }
      else {
        val_bar_pdf = BCPV_bar_PDF(val_t, val_d, par_tau, par_dGamma, par_dm, par_Sfbar, par_Cfbar, par_Dfbar);
        val_envelope = BCPV_PDF_Envelope(val_t, gamma_min, par_Sfbar, par_Cfbar, par_Dfbar);
        if (val_envelope < val_bar_pdf) std::cout << "WARNING: Envelope smaller than PDF!" << std::endl;
        if(val_envelope*rndm.Uniform() > val_bar_pdf) continue;
        else break;
      }
    }
  }

  obs_tag_true  = val_d;
  obs_time_true = val_t;
  finalstate    = val_final;
  mix_true      = val_final * val_d;

  return true;
}

double BCPV_PDF(double t, double d, double tau, double dGamma, double dm,
                double Sf, double Cf, double Df) {
  return exp(-t/tau)*(cosh(dGamma*t/2.)+Df*sinh(dGamma*t/2.)+d*Cf*cos(dm*t)-d*Sf*sin(dm*t));
}

double BCPV_bar_PDF(double t, double d, double tau, double dGamma, double dm,
                    double Sfbar, double Cfbar, double Dfbar) {
  return exp(-t/tau)*(cosh(dGamma*t/2.)+Dfbar*sinh(dGamma*t/2.)-d*Cfbar*cos(dm*t)+d*Sfbar*sin(dm*t));
}

double BCPV_PDF_Envelope(double t, double gamma_min, double S, double C, double D) {
  return exp(-t*gamma_min)*(1.+std::abs(D)+sqrt(S*S+C*C));
}

bool GenerateResolSingleGauss(TRandom& rndm, double par_bias, double par_sigma, double obs_true, double& obs_meas) {
  obs_meas = obs_true;
  obs_meas += rndm.Gaus(par_bias, par_sigma);

  return true;
}

bool GenerateResolSingleGaussPerEvent(TRandom& rndm, double par_bias, double par_scale, double obs_per_event_error, double obs_true, double& obs_meas) {
  obs_meas = obs_true;
  obs_meas += rndm.Gaus(par_bias, par_scale*obs_per_event_error);

  return true;
}


bool GenerateEtaFlat(TRandom& rndm, double& obs_eta) {
  return GenerateEtaFlat(rndm,0.0,0.5,obs_eta);
}

bool GenerateEtaFlat(TRandom& rndm, double obs_eta_min, double obs_eta_max, double& obs_eta) {
  obs_eta = rndm.Uniform(obs_eta_min,obs_eta_max);
  return true;
}

bool GenerateEtaGauss(TRandom& rndm, double m, double s, double obs_eta_min, double obs_eta_max, double& obs_eta) {
  unsigned int num_samples(0);

  // s is set to -1.0 on default; a negative Gaussian width does not make sense, thus generate a uniform distribution
  if (s < 0.0) {
    obs_eta = rndm.Uniform(obs_eta_min,obs_eta_max);
    return true;
  } else {
    obs_eta = rndm.Gaus(m,s);
    ++num_samples;

    while (obs_eta > obs_eta_max || obs_eta < obs_eta_min) {
      obs_eta = rndm.Gaus(m,s);
      ++num_samples;

      if (num_samples % 10000000 == 0) {
        std::cout << "WARNING in cptoymc::generator::GenerateEtaGauss(rndm, m=" << m << ", s=" << s << ", obs_eta_min=" << obs_eta_min << ", obs_eta_max=" << obs_eta_max << "): Generated " << num_samples << " sample values without one candidate passing. You probably want to check your parameters." << std::endl;
      }
    }

    return true;
  }
}

bool GenerateTag(TRandom& rndm, double par_omega, double par_domega,
                 const int par_tag_true_B, const int par_tag_true_Bb,
                 const int par_tag_B,      const int par_tag_Bb,
                 int obs_tag_true, int& obs_tag_meas) {
  //if (par_omega > 0.5) par_omega = 0.5;
  if (par_omega < 0.0) par_omega = 0.;

  int correct_tag = 0;

  // always assume that omega_B  = omega + dOmega/2 for true B mesons
  // and that           omega_Bb = omega - dOmega/2 for true Bb mesons
  if (obs_tag_true == par_tag_true_B) {
    par_omega += par_domega/2.; // B meson case
    correct_tag = par_tag_B;
  }
  else if (obs_tag_true == par_tag_true_Bb) {
    par_omega -= par_domega/2.; // Bb meson case
    correct_tag = par_tag_Bb;
  } else {
    std::cout << "Cannot interpret true tag of " << obs_tag_true << ". Failed." << std::endl;
    return false;
  }

  if (rndm.Uniform() < par_omega) {
    obs_tag_meas = -1*correct_tag;
  } else {
    obs_tag_meas = correct_tag;
  }
  return true;
}


bool GenerateTag(TRandom& rndm, double par_omega, double par_domega, int obs_tag_true, int& obs_tag_meas) {
  return GenerateTag(rndm, par_omega, par_domega, +1, -1, +1, -1, obs_tag_true, obs_tag_meas);
}

bool GenerateTag(TRandom& rndm,
                 std::function<double(double)>& func_omega,
                 std::function<double(double)>& func_domega,
                 const int par_tag_true_B, const int par_tag_true_Bb,
                 const int par_tag_B,      const int par_tag_Bb,
                 int obs_tag_true, double obs_eta, int& obs_tag_meas) {
  return GenerateTag(rndm, func_omega(obs_eta), func_domega(obs_eta),
                     par_tag_true_B, par_tag_true_Bb, par_tag_B, par_tag_Bb,
                     obs_tag_true, obs_tag_meas);
}

bool GenerateRandomTag(TRandom& rndm, int& obs_tag_meas) {
  obs_tag_meas = (rndm.Uniform() < 0.5) ? +1 : -1;
  return true;
}




int yieldToGenerate(TRandom& rndm, double yield_exp) {
  return rndm.Poisson(yield_exp);
}


} // namespace cptoymc
} // namespace generator

