#include <iostream>
#include <sstream>

// from ROOT
#include "TCanvas.h"
#include "TPad.h"
#include "TFile.h"
#include "TROOT.h"
#include "TH1.h"
#include "TLine.h"

//from RooFit
#include "RooCmdArg.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "RooArgList.h"
#include "RooCategory.h"
#include "RooDataHist.h"
#include "RooConstVar.h"
#include "RooLinkedList.h"
#include "RooNumIntConfig.h"
#include "RooBinning.h"

// from RooFit PDFs
#include "RooGaussModel.h"
#include "RooAbsPdf.h"
#include "RooBDecay.h"
#include "RooDecay.h"
#include "RooProdPdf.h"
#include "RooGaussian.h"
#include "RooAddPdf.h"
#include "RooExponential.h"
#include "RooHistPdf.h"
#include "RooExtendPdf.h"
#include "RooSimultaneous.h"
#include "RooBCPGenDecay.h"

// from DooCore
#include "doocore/io/MsgStream.h"
#include "doocore/lutils/lutils.h"
#include "doocore/io/EasyTuple.h"

// from DooFit
#include "doofit/roofit/functions/SinCoeffWithProdAsymm.h"
// #include "doofit/roofit/functions/CosCoeffWithProdAsymm.h"
// #include "doofit/roofit/functions/CosCoeffCombo.h"
#include "doofit/roofit/functions/SinCoeffCombo.h"
#include "doofit/roofit/functions/CoshCoeff.h"
#include "doofit/roofit/functions/CoshCoeffCombo.h"
#include "doofit/roofit/functions/SingleMistagCalibrationWithAsymmetries.h"
#include "doofit/roofit/pdfs/BiasDelta.h"
#include "doofit/config/CommonConfig.h"
#include "doofit/toy/ToyFactoryStd/ToyFactoryStd.h"
#include "doofit/toy/ToyFactoryStd/ToyFactoryStdConfig.h"
#include "doofit/toy/ToyStudyStd/ToyStudyStd.h"
#include "doofit/toy/ToyStudyStd/ToyStudyStdConfig.h"

#include "configuration/ToyConfig.h"
#include "generator/ToyGenerator.h"

using namespace RooFit;
using namespace std;
using namespace doocore::io;
using namespace doocore::lutils;
using namespace doofit::roofit::functions;
using namespace doofit::roofit::pdfs;
using namespace doofit::toy;
using namespace cptoymc::configuration;
using namespace cptoymc::generator;

int main(int argc, char * argv[]){
  
  // Plots
  // gROOT->SetStyle("Plain");
  // setStyle("LHCb");
  // TCanvas canvas("canvas","canvas",800,600);
  // RooPlot* plot_frame;
  // TLatex label(0.7,0.83,"Toys");
  // label.SetNDC();
  
  if (argc < 6)
  {
    cout  <<  "You need to provide:"  <<  endl;
    cout  <<  "the method (g for generation of toys, e for evaluation of toy study)"  <<  endl;
    cout  <<  "the number of toys to generate"  <<  endl;
    cout  <<  "the random seed"  <<  endl;
    cout  <<  "the number of CPUs to fit on"  <<  endl;
    cout  <<  "the name of the config file for the toy production"  <<  endl;
    cout  <<  "the name of the starting values file" <<  endl;
    return 1;
  }

  TString           method = argv[1];
  TString           StringNumToys = argv[2];
  int num_toys = StringNumToys.Atoi();
  TString           StringRandomSeed = argv[3];
  int random_seed = StringRandomSeed.Atoi();
  TString           StringNumCPU = argv[4];
  int num_cpu = StringNumCPU.Atoi();
  
  RooRealVar        obsTime("obsTime","#it{t}",0.,18.,"ps");
  RooCategory       obsTagOS("obsTagOS","Flavour Tag");
  obsTagOS.defineType("B0",1);
  obsTagOS.defineType("B0bar",-1);
  obsTagOS.defineType("untagged",0);
  RooCategory       obsTagSS("obsTagSS","Flavour Tag");
  obsTagSS.defineType("B0",1);
  obsTagSS.defineType("B0bar",-1);
  obsTagSS.defineType("untagged",0);
  
  RooArgSet         observables(obsTime,obsTagOS,obsTagSS,"observables");
  
  // Resolution model
  RooRealVar        parResMean("parResMean","parResMean",0.);
  RooRealVar        parResSigma("parResSigma","parResSigma",0.05);
  RooGaussModel     resGauss("resGauss","Resolution model",obsTime,parResMean,parResSigma);

  // Decay Time PDF
  RooRealVar        parSigTimeTau("parSigTimeTau","#tau",1.5,1.,2.);
  RooRealVar        parSigTimeDeltaM("parSigTimeDeltaM","B0 mixing",0.517);            // PDG2013 (0.510 +- 0.004)/ps
  RooRealVar        parSigTimeDeltaG("parSigTimeDeltaG","#Delta#Gamma",0.);

  // Additional parameters needed for B Decay with CPV
  RooRealVar        parSigTimeSin2b("parSigTimeSin2b","#it{S_{J/#kern[-0.3]{#psi} K_{S}}}",0.7,0.6,0.8);
  RooRealVar        parSigTimeCjpsiKS("parSigTimeCjpsiKS","#it{C_{J/#kern[-0.3]{#psi} K_{S}}}",0,-0.5,0.5);
      
  // Tagging asymmetries
  RooRealVar        parSigEtaDeltaProd("parSigEtaDeltaProd","asymmetry",0.02,-0.5,0.5);
  RooConstVar       parSigEtaMean_OS("parSigEtaMean_OS","parSigEtaMean_OS",0.25);
  RooConstVar       parSigEtaMean_SS("parSigEtaMean_SS","parSigEtaMean_SS",0.25);

  // Decay Time PDF
  // RooBDecay params
  RooConstVar       parSigTimeSinh("parSigTimeSinh","Sh_{f}",0.0);
  CoshCoeffCombo    parSigTimeCosh_Combo("parSigTimeCosh_Combo",obsTagOS,RooConst(0.),RooConst(1.),RooConst(0.),parSigEtaMean_OS,RooConst(0.),RooConst(0.),obsTagSS,RooConst(0.),RooConst(1.),RooConst(0.),parSigEtaMean_SS,RooConst(0.),RooConst(0.),parSigEtaDeltaProd);
  SinCoeffCombo     parSigTimeSin_Combo("parSigTimeSin_Combo",parSigTimeSin2b,obsTagOS,RooConst(0.),RooConst(1.),RooConst(0.),parSigEtaMean_OS,RooConst(0.),RooConst(0.),obsTagSS,RooConst(0.),RooConst(1.),RooConst(0.),parSigEtaMean_SS,RooConst(0.),RooConst(0.),parSigEtaDeltaProd,SinCoeffCombo::kSType);
  SinCoeffCombo     parSigTimeCos_Combo("parSigTimeCos_Combo",parSigTimeCjpsiKS,obsTagOS,RooConst(0.),RooConst(1.),RooConst(0.),parSigEtaMean_OS,RooConst(0.),RooConst(0.),obsTagSS,RooConst(0.),RooConst(1.),RooConst(0.),parSigEtaMean_SS,RooConst(0.),RooConst(0.),parSigEtaDeltaProd,SinCoeffCombo::kCType);

  RooBDecay               pdfSigTime_Combo("pdfSigTime_Combo","P_{S}^{l}(t,d|#eta)",obsTime,parSigTimeTau,parSigTimeDeltaG,parSigTimeCosh_Combo,parSigTimeSinh,parSigTimeCos_Combo,parSigTimeSin_Combo,parSigTimeDeltaM,resGauss,RooBDecay::SingleSided);
  if (method.EqualTo("g") || method.EqualTo("e")) {
            
    doofit::config::CommonConfig cfg_com("common");
    cfg_com.InitializeOptions(argc, argv);
    
    ToyFactoryStdConfig cfg_tfac("toyfac");
    cfg_tfac.InitializeOptions(cfg_com);
    
    ToyStudyStdConfig cfg_tstudy("toystudy");
    cfg_tstudy.InitializeOptions(cfg_tfac);
    
    cfg_com.CheckHelpFlagAndPrintHelp();
        
    cfg_com.PrintAll();
    
    ToyStudyStd tstudy(cfg_com, cfg_tstudy);
    
    RooLinkedList fitting_args;
    fitting_args.Add((TObject*)(new RooCmdArg(NumCPU(num_cpu))));
    fitting_args.Add((TObject*)(new RooCmdArg(Minos(false))));
    fitting_args.Add((TObject*)(new RooCmdArg(Strategy(2))));
    fitting_args.Add((TObject*)(new RooCmdArg(Save(true))));
    fitting_args.Add((TObject*)(new RooCmdArg(Timer(true))));
    fitting_args.Add((TObject*)(new RooCmdArg(Minimizer("Minuit2","migrad"))));
    fitting_args.Add((TObject*)(new RooCmdArg(SumW2Error(false))));
    fitting_args.Add((TObject*)(new RooCmdArg(Extended(false))));
    fitting_args.Add((TObject*)(new RooCmdArg(Optimize(1))));
    // fitting_args.Add((TObject*)(new RooCmdArg(ConditionalObservables(obsEtaOS))));

    ToyConfig     cfg_cptoymc;
    cfg_cptoymc.load(argv[5]);
    ToyGenerator  cptoymc(cfg_cptoymc);
    
    RooDataSet* data = NULL;

    if (method.EqualTo("g")) {
      for (int i = 0; i < num_toys; ++i) {
        cout  <<  i <<  endl;
        try {
          // TFile       outfile("ToyMC.root","recreate");
          TTree       tree("ToyMCTreetree","Tree of generation");
          cptoymc.GenerateToy(tree,random_seed);
          // outfile.Write();
          // outfile.Close();
          data = new RooDataSet("data","Toy MC data",&tree, observables);
          // data->Print("v");
          pdfSigTime_Combo.getParameters(*data)->readFromFile(argv[6]);
          RooFitResult* fit_result = pdfSigTime_Combo.fitTo(*data,fitting_args);
          // fit_result->Print("v");
          tstudy.StoreFitResult(fit_result);
          delete data;
        } catch (...) {
          i--;
        }
      }
    }
    
    if (method.EqualTo("e")) {
      tstudy.ReadFitResults();
      tstudy.EvaluateFitResults();
      tstudy.PlotEvaluatedParameters();
    }
  }
  
  return 0 ;
}