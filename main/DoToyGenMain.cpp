// STL
#include <cmath>
#include <iostream>
#include <memory>
#include <vector>

// ROOT
#include "TBranch.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"

// from Project
#include "configuration/ToyConfig.h"
#include "generator/ToyGenerator.h"




int main() {
  using cptoymc::configuration::ToyConfig;
  using cptoymc::generator::ToyGenerator;

  ToyConfig config;
  config.load("../configs/Bd2JpsiKS_dirtyconfig.info");

  TFile out_file("ToyMC.root","RECREATE");
  TTree out_tree("ToyMCTree","A tree filled with ToyMC observables");
  
  ToyGenerator generator(config);
  generator.GenerateToy(out_tree,1234);
  
  out_file.Write();
  out_file.Close();

  return 0;
}


