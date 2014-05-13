#ifndef CPTOYMC_GENERATOR_OBSERVABLES_H
#define CPTOYMC_GENERATOR_OBSERVABLES_H

namespace cptoymc {
namespace generator {

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
