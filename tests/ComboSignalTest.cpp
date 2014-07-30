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
  RooRealVar        obsEtaSS("obsEtaSS","#eta_{SS#pi}",0.,0.5);
  RooCategory       obsTagSS("obsTagSS","Flavour Tag");
  obsTagSS.defineType("B0",1);
  obsTagSS.defineType("B0bar",-1);
  
  RooArgSet         observables(obsTime,obsMass,obsTagOS,obsTagSS,obsEtaOS,obsEtaSS,"observables");
  
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
  CoshCoeffCombo    parSigTimeCosh_Combo("parSigTimeCosh_Combo",obsTagOS,parSigEtaMean_OS,RooConst(0.),RooConst(0.),obsEtaOS,RooConst(0.),RooConst(0.),obsTagSS,parSigEtaMean_SS,RooConst(0.),RooConst(0.),obsEtaSS,RooConst(0.),RooConst(0.),parSigEtaDeltaProd);
  SinCoeffCombo     parSigTimeSin_Combo("parSigTimeSin_Combo",parSigTimeSin2b,obsTagOS,parSigEtaMean_OS,RooConst(0.),RooConst(0.),obsEtaOS,RooConst(0.),RooConst(0.),obsTagSS,parSigEtaMean_SS,RooConst(0.),RooConst(0.),obsEtaSS,RooConst(0.),RooConst(0.),parSigEtaDeltaProd,SinCoeffCombo::kSType);
  SinCoeffCombo     parSigTimeCos_Combo("parSigTimeCos_Combo",parSigTimeCjpsiKS,obsTagOS,parSigEtaMean_OS,RooConst(0.),RooConst(0.),obsEtaOS,RooConst(0.),RooConst(0.),obsTagSS,parSigEtaMean_SS,RooConst(0.),RooConst(0.),obsEtaSS,RooConst(0.),RooConst(0.),parSigEtaDeltaProd,SinCoeffCombo::kCType);

  RooBDecay         pdfSigTime_Combo("pdfSigTime_Combo","P_{S}^{l}(t,d|#eta)",obsTime,parSigTimeTau,parSigTimeDeltaG,parSigTimeCosh_Combo,parSigTimeSinh,parSigTimeCos_Combo,parSigTimeSin_Combo,parSigTimeDeltaM,resGauss,RooBDecay::SingleSided);

  // Mass PDF
  RooRealVar        parSigMassMean("parSigMassMean","Bd Mean Mass",5279,5270,5290,"MeV/c^{2}");
  RooRealVar        parSigMassSigma("parSigMassSigma","Sigma of Gaussian Mass",8.0,4.0,12.0,"MeV/c^{2}");
  RooGaussian       pdfSigMass("pdfSigMass","Mass PDF",obsMass,parSigMassMean,parSigMassSigma);

  // Full PDF
  RooProdPdf        pdf_Combo("pdf_Combo","pdf_Combo",RooArgList(pdfSigTime_Combo,pdfSigMass));
  RooRealVar        parSigYield_Combo("parSigYield_Combo","parSigYield_Combo",100000,0,200000);
  RooExtendPdf      pdfExtend_Combo("pdfExtend_Combo","pdfExtend_Combo",pdf_Combo,parSigYield_Combo);

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
    fitting_args.Add((TObject*)(new RooCmdArg(Extended(true))));
    fitting_args.Add((TObject*)(new RooCmdArg(Optimize(1))));

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
          pdfExtend_Combo.getParameters(*data)->readFromFile(argv[6]);
          RooFitResult* fit_result = pdfExtend_Combo.fitTo(*data,fitting_args);
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