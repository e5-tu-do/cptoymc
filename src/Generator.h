#ifndef CPTOYMC_GENERATOR_H
#define CPTOYMC_GENERATOR_H

class TRandom;

namespace cptoymc {

namespace configuration {
struct ParsMass;
struct ParsTimeAndCP;
struct ParsMassResol;
struct ParsTimeResol;
struct ParsTagging;
} // namespace configuration
  
namespace generator {

struct Observables {
  double mass_true;
  double time_true;
  int    tag_true;
  double mass;
  double time;
  int    tag_class; // -1: SS&&!OS, 0: !OS&&!SS, 1: OS&&!SS, 10: OS&&SS
  int    tag_OS;
  double eta_OS;
  int    tag_SS;
  double eta_SS;
};
  
void resetObs(cptoymc::generator::Observables& obs);

void generateTrueMass(TRandom& rndm, const cptoymc::configuration::ParsMass& pars_mass, double& obs_mass_true);
void generateTrueTimeAndTag(TRandom& rndm, const cptoymc::configuration::ParsTimeAndCP& pars_time_CP, double& obs_time_true, int& obs_tag_true);

void generateMass(TRandom& rndm, const cptoymc::configuration::ParsMassResol& pars_mass_resol, const double obs_mass_true, double& obs_mass);
void generateTime(TRandom& rndm, const cptoymc::configuration::ParsTimeResol& pars_time_resol, const double obs_time_true, double& obs_time);
void generateTagAndEta(TRandom& rndm, const cptoymc::configuration::ParsTagging& pars_tagging, const int tag_true, int& tag_OS, double& eta_OS, int& tag_SS, double& eta_SS, int& tag_class);

void generateTagAndEtaOS(TRandom& rndm, const cptoymc::configuration::ParsTagging& pars_tagging, const int tag_true, int& tag_OS, double& eta_OS);
void generateTagAndEtaSS(TRandom& rndm, const cptoymc::configuration::ParsTagging& pars_tagging, const int tag_true, int& tag_SS, double& eta_SS);

double BCPV_PDF(double t, double d, double tau, double dGamma, double dm, double Sf, double Cf, double Df);
double BCPV_PDF_Envelope(double t, double gamma_min, double Sf, double Cf, double Df);

} // namespace generator
} // namespace cptoymc

#endif // CPTOYMC_GENERATOR_H