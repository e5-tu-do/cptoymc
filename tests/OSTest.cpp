#include <iostream>
#include <sstream>

// from ROOT
#include "TCanvas.h"
#include "TPad.h"
#include "TFile.h"
#include "TROOT.h"
#include "TH1.h"
#include "TLine.h"
#include "TTree.h"

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

// from DooFit
#include "doofit/roofit/functions/CPCoefficient.h"
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
using namespace doofit::roofit::functions;
using namespace doofit::toy;
using namespace cptoymc::configuration;
using namespace cptoymc::generator;

int main(int argc, char * argv[]){
  
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
  
  RooRealVar        obsTime("obsTime","#it{t}",-2.,18.,"ps");
  RooRealVar        obsMass("obsMass","#it{m_{J/#kern[-0.3]{#psi} K_{S}}}",5230,5330,"MeV/c^{2}");
  RooRealVar        obsEtaOS("obsEtaOS","#eta_{OS}",0.,0.5);
  RooCategory       obsTagOS("obsTagOS","Flavour Tag");
  obsTagOS.defineType("B0",1);
  obsTagOS.defineType("B0bar",-1);
  
  RooArgSet         observables(obsTime,obsMass,obsTagOS,"observables");
  
  // Resolution model
  RooRealVar        parResMean("parResMean","parResMean",0.);
  RooRealVar        parResSigma("parResSigma","parResSigma",0.05);
  RooGaussModel     resGauss("resGauss","Resolution model",obsTime,parResMean,parResSigma);
  
  // Mass model
  RooRealVar        parSigMassMean("parSigMassMean","Bd Mean Mass",5279,5270,5290,"MeV/c^{2}");
  RooRealVar        parSigMassSigma("parSigMassSigma","Sigma of Gaussian Mass",8.0,4.0,12.0,"MeV/c^{2}");
  RooGaussian       pdfSigMass("pdfSigMass","Mass PDF",obsMass,parSigMassMean,parSigMassSigma);

  // Decay Time PDF
  RooRealVar        parSigTimeTau("parSigTimeTau","#tau",1.5,1.,2.);
  RooRealVar        parSigTimeDeltaM("parSigTimeDeltaM","B0 mixing",0.510);            // PDG2013 (0.510 +- 0.004)/ps
  RooRealVar        parSigTimeDeltaG("parSigTimeDeltaG","#Delta#Gamma",0.);

  // Additional parameters needed for B Decay with CPV
  RooRealVar        parSigTimeSin2b("parSigTimeSin2b","#it{S_{J/#kern[-0.3]{#psi} K_{S}}}",0.7,0.6,0.8);
  RooRealVar        parSigTimeCjpsiKS("parSigTimeCjpsiKS","#it{C_{J/#kern[-0.3]{#psi} K_{S}}}",0,-0.5,0.5);
      
  // Tagging asymmetries
  RooRealVar        parSigEtaDeltaProd("parSigEtaDeltaProd","asymmetry",0.02,-0.5,0.5);
  RooConstVar       parSigTimeDelta("parSigTimeDelta","asymmetry",0.);
  RooConstVar       parSigEtaMean_OS("parSigEtaMean_OS","parSigEtaMean_OS",0.25);

  // Decay Time PDF
  // RooBDecay params
  RooConstVar       parSigTimeSinh("parSigTimeSinh","Sh_{f}",0.0);
  CPCoefficient     parSigTimeCosh_OS("parSigTimeCosh_OS",RooConst(1.0),obsTagOS,parSigEtaMean_OS,RooConst(1.),parSigEtaMean_OS,obsEtaOS,RooConst(0.),RooConst(0.),parSigEtaDeltaProd,CPCoefficient::kCosh);
  CPCoefficient     parSigTimeSin_OS("parSigTimeSin_OS",parSigTimeSin2b,obsTagOS,parSigEtaMean_OS,RooConst(1.),parSigEtaMean_OS,obsEtaOS,RooConst(0.),RooConst(0.),parSigEtaDeltaProd,CPCoefficient::kSin);
  CPCoefficient     parSigTimeCos_OS("parSigTimeCos_OS",parSigTimeCjpsiKS,obsTagOS,parSigEtaMean_OS,RooConst(1.),parSigEtaMean_OS,obsEtaOS,RooConst(0.),RooConst(0.),parSigEtaDeltaProd,CPCoefficient::kCos);
  
  RooDecay                pdfSigTimeDecay_OS("pdfSigTimeDecay_OS","P_{S}^{OS}(t)",obsTime,parSigTimeTau,resGauss,RooDecay::SingleSided);
  // RooBCPGenDecay          pdfSigTime_OS("pdfSigTime_OS","P_{S}^{l}(t,d|#eta)",obsTime,obsTagOS,parSigTimeTau,parSigTimeDeltaM,obsEtaOS,parSigTimeCjpsiKS,parSigTimeSin2b,parSigTimeDelta,parSigEtaDeltaProd,resGauss,RooBCPGenDecay::SingleSided);
  RooBDecay               pdfSigTime_OS("pdfSigTime_OS","P_{S}^{l}(t,d|#eta)",obsTime,parSigTimeTau,parSigTimeDeltaG,parSigTimeCosh_OS,parSigTimeSinh,parSigTimeCos_OS,parSigTimeSin_OS,parSigTimeDeltaM,resGauss,RooBDecay::SingleSided);
  
  // Combination of observables
  RooProdPdf              pdfSig_OS("pdfSig_OS","pdfSig_OS",RooArgList(pdfSigTime_OS,pdfSigMass));

  // Background
  // Mass model
  RooRealVar              parBkgMassExponent("parBkgMassExponent","Background Mass Exponent",-0.002,-1.,1.,"1/(MeV/c^{2})");
  RooExponential          pdfBkgMass("pdfBkgMass","Background Mass PDF",obsMass,parBkgMassExponent);

  // Decay time model
  RooRealVar              parBkgTimeTau("parBkgTimeTau","parBkgTimeTau",1,0.5,1.5);
  RooDecay                pdfBkgTime("pdfBkgTime","pdfBkgTime",obsTime,parBkgTimeTau,resGauss,RooDecay::SingleSided);
  
  // Combination of observables
  RooProdPdf              pdfBkg_OS("pdfBkg_OS","pdfBkg_OS",RooArgList(pdfBkgTime,pdfBkgMass));
  
  // Combining signal and background
  RooRealVar              parSigYield_OS("parSigYield_OS","parSigYield_OS",100000,0,200000);
  RooRealVar              parBkgYield_OS("parBkgYield_OS","parBkgYield_OS",100000,0,200000);
  RooAddPdf               pdf_OS("pdf_OS","pdf_OS",RooArgList(pdfSig_OS,pdfBkg_OS),RooArgList(parSigYield_OS,parBkgYield_OS));
  
  if (method.EqualTo("g") || method.EqualTo("e")) {

    doofit::config::CommonConfig cfg_com("common");
    cfg_com.InitializeOptions(argc, argv);
    
    ToyFactoryStdConfig cfg_tfac("toyfac");
    cfg_tfac.InitializeOptions(cfg_com);
    
    ToyStudyStdConfig cfg_tstudy("toystudy");
    cfg_tstudy.InitializeOptions(cfg_tfac);
    
    cfg_com.CheckHelpFlagAndPrintHelp();
        
    // RooWorkspace* ws = new RooWorkspace("ws");
    // ws->import(pdf_OS);
    // ws->defineSet("observables",observables);
    // ws->Print();
    // cfg_tfac.set_workspace(ws);
    // ToyFactoryStd tfac(cfg_com, cfg_tfac);

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
    fitting_args.Add((TObject*)(new RooCmdArg(Extended(true))));
    fitting_args.Add((TObject*)(new RooCmdArg(Optimize(1))));
    fitting_args.Add((TObject*)(new RooCmdArg(ConditionalObservables(obsEtaOS))));
    
    ToyConfig     cfg_cptoymc;
    cfg_cptoymc.load(argv[5]);
    ToyGenerator  cptoymc(cfg_cptoymc);
    
    RooDataSet* data = NULL;

    if (method.EqualTo("g")) {
      for (int i = 0; i < num_toys; ++i) {
        cout  <<  i <<  endl;
        try {
          TTree       tree("ToyMCTreetree","Tree of generation");
          cptoymc.GenerateToy(tree,random_seed);
          data = new RooDataSet("data","Toy MC data",&tree, observables);
          // data = tfac.Generate();
          pdf_OS.getParameters(*data)->readFromFile(argv[6]);
          RooFitResult* fit_result = pdf_OS.fitTo(*data,fitting_args);
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