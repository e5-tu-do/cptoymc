#include "generator/Observables.h"

namespace cptoymc {
namespace generator {


Observables::Observables() :
  mass_true(-1000.),
  time_true(-1000.),
  tag_true(0),
  mass(-1000.),
  time(-1000.),
  tag_class(-10000),
  tag_OS(0),
  eta_OS(0.5),
  tag_SS(0),
  eta_SS(0.5),
  bkg_cat(-1000) 
{

}

Observables::~Observables() {
  
}


void Observables::reset() {
  mass_true = -1000.;
  time_true = -1000.;
  tag_true  = 0;
  mass      = -1000.;
  time      = -1000.;
  tag_class = -10000;
  tag_OS    = 0;
  eta_OS    = 0.5;
  tag_SS    = 0;
  eta_SS    = 0.5;
  bkg_cat   = -1000;
}

} // namespace generator
} // namespace cptoymc

