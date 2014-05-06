// STL
#include <cmath>
#include <iostream>
#include <vector>

// ROOT
#include "TBranch.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"

// from Project
#include "Configuration.h"
#include "Generator.h"




int main() {
  using namespace cptoymc::configuration;

  ToyConfig config;
  config.load("../configs/Bd2JpsiKS_dirtyconfig.info");



  // int num_events = 100000;

  // using namespace cptoymc;
  
  // //============================================================================
  // // Physics Parameters
  // //----------------------------------------------------------------------------
  // // Mass parameters
  // configuration::ParsMass pars_mass;
  // pars_mass.mean  = 5279.58;
  // pars_mass.gamma = 0.;

  // //----------------------------------------------------------------------------
  // // Lifetime and CPV parameters
  // configuration::ParsTimeAndCP pars_time_CP;
  // pars_time_CP.tau    = 1.519;
  // pars_time_CP.dGamma = 0.;
  // pars_time_CP.dm     = 0.510;
  // pars_time_CP.Sf     = 0.7;
  // pars_time_CP.Cf     = 0.0;
  // pars_time_CP.Df     = +sqrt(pars_time_CP.Sf*pars_time_CP.Sf+pars_time_CP.Cf*pars_time_CP.Cf);
  // pars_time_CP.AP     = 0.;

  // //============================================================================
  // // Experimental Parameters
  // //----------------------------------------------------------------------------
  // // Mass parameters
  // configuration::ParsMassResol pars_mass_resol;
  // pars_mass_resol.bias  = 0.;
  // pars_mass_resol.sigma = 8.;

  // //----------------------------------------------------------------------------
  // // Decay time parameters
  // configuration::ParsTimeResol pars_time_resol;
  // pars_time_resol.bias  = 0.;
  // pars_time_resol.sigma = 0.05;

  // //----------------------------------------------------------------------------
  // // Tagging parameters
  // configuration::ParsTagging pars_tagging;
  // pars_tagging.eff_OS   = 0.3292;
  // pars_tagging.eff_SS   = 0.0625;
  // pars_tagging.eff_SSOS = 0.0436;

  // pars_tagging.w_OS    = 0.25; // no effect
  // pars_tagging.w_SS    = 0.25; // no effect

  // pars_tagging.dw_OS    = 0.;
  // pars_tagging.dw_SS    = 0.;


  // //============================================================================
  // // Initialise TFile and TTree
  // TFile out_file("ToyMC.root","RECREATE");
  // TTree out_tree("ToyMCTree","A tree filled with ToyMC observables");

  // generator::Observables obs;

  // out_tree.Branch("obsMassTrue"       , &(obs.mass_true)  , "obsMassTrue/D"       );
  // out_tree.Branch("obsTimeTrue"       , &(obs.time_true)  , "obsTimeTrue/D"       );
  // out_tree.Branch("obsTagTrue"        , &(obs.tag_true)   , "obsTagTrue/I"        );
  // out_tree.Branch("obsMass"           , &(obs.mass)       , "obsMass/D"           );
  // out_tree.Branch("obsTime"           , &(obs.time)       , "obsTime/D"           );
  // out_tree.Branch("catTaggedOSSSPion" , &(obs.tag_class)  , "catTaggedOSSSPion/I" );
  // out_tree.Branch("obsTagOS"          , &(obs.tag_OS)     , "obsTagOS/I"          );
  // out_tree.Branch("obsEtaOS"          , &(obs.eta_OS)     , "obsEtaOS/D"          );
  // out_tree.Branch("obsTagSSPion"      , &(obs.tag_SS)     , "obsTagSSPion/I"      );
  // out_tree.Branch("obsEtaSSPion"      , &(obs.eta_SS)     , "obsEtaSSPion/D"      );


  // //============================================================================
  // // Generate
  // //----------------------------------------------------------------------------
  // TRandom3 rndm(4357);

  // for (int i = 0; i < num_events; ++i) {
  //   generator::resetObs(obs);
  //   generator::generateTrueMass(rndm, pars_mass, obs.mass_true);
  //   generator::generateTrueTimeAndTag(rndm, pars_time_CP, obs.time_true, obs.tag_true);

  //   generator::generateMass(rndm, pars_mass_resol, obs.mass_true, obs.mass);
  //   generator::generateTime(rndm, pars_time_resol, obs.time_true, obs.time);
  //   generator::generateTagAndEta(rndm, pars_tagging, obs.tag_true, obs.tag_OS, obs.eta_OS, obs.tag_SS, obs.eta_SS, obs.tag_class);


  //   out_tree.Fill();
  // }


  // //============================================================================
  // // Finalise TFile
  // out_file.Write();
  // out_file.Close();

  return 0;
}


