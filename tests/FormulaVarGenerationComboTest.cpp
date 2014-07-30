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
#include "RooWorkspace.h"
#include "RooFormulaVar.h"

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
#include "RooProduct.h"

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

using namespace RooFit;
using namespace std;
using namespace doocore::io;
using namespace doocore::lutils;
using namespace doofit::roofit::functions;
using namespace doofit::roofit::pdfs;
using namespace doofit::toy;

int main(int argc, char * argv[]){
  
  if (argc < 6)
  {
    cout  <<  "You need to provide:"  <<  endl;
    cout  <<  "the method (g for generation of toys, e for evaluation of toy study)"  <<  endl;
    cout  <<  "the number of toys to generate"  <<  endl;
    cout  <<  "the number of CPUs to fit on"  <<  endl;
    cout  <<  "the name of the starting values file" <<  endl;
    return 1;
  }

  TString           method = argv[1];
  TString           StringNumToys = argv[2];
  int num_toys = StringNumToys.Atoi();
  TString           StringNumCPU = argv[3];
  int num_cpu = StringNumCPU.Atoi();
  
  RooRealVar        obsTime("obsTime","#it{t}",-2.,18.,"ps");
  RooRealVar        obsMass("obsMass","#it{m_{J/#kern[-0.3]{#psi} K_{S}}}",5230,5330,"MeV/c^{2}");
  RooCategory       obsTagOS("obsTagOS","Flavour Tag");
  obsTagOS.defineType("B0",1);
  obsTagOS.defineType("B0bar",-1);
  RooCategory       obsTagSS("obsTagSS","Flavour Tag");
  obsTagSS.defineType("B0",1);
  obsTagSS.defineType("B0bar",-1);
  RooCategory       obsTagOSSS("obsTagOSSS","same or different tag decision");
  obsTagOSSS.defineType("same",1);
  obsTagOSSS.defineType("different",-1);
  
  RooArgSet         observables(obsTime,obsMass,obsTagOS,obsTagSS,obsTagOSSS,"observables");
  
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
  RooConstVar       parSigEtaMean_OS("parSigEtaMean_OS","parSigEtaMean_OS",0.25);
  RooConstVar       parSigEtaMean_SS("parSigEtaMean_SS","parSigEtaMean_SS",0.22);

  // Decay Time PDF
  // RooBDecay params
  RooConstVar       parSigTimeSinh("parSigTimeSinh","Sh_{f}",0.0);
  CoshCoeffCombo    parSigTimeCosh_Combo("parSigTimeCosh_Combo",obsTagOS,parSigEtaMean_OS,RooConst(0.),RooConst(0.),parSigEtaMean_OS,RooConst(0.),RooConst(0.),obsTagSS,parSigEtaMean_SS,RooConst(0.),RooConst(0.),parSigEtaMean_SS,RooConst(0.),RooConst(0.),parSigEtaDeltaProd);
  SinCoeffCombo     parSigTimeSin_Combo("parSigTimeSin_Combo",parSigTimeSin2b,obsTagOS,parSigEtaMean_OS,RooConst(0.),RooConst(0.),parSigEtaMean_OS,RooConst(0.),RooConst(0.),obsTagSS,parSigEtaMean_SS,RooConst(0.),RooConst(0.),parSigEtaMean_SS,RooConst(0.),RooConst(0.),parSigEtaDeltaProd,SinCoeffCombo::kSType);
  SinCoeffCombo     parSigTimeCos_Combo("parSigTimeCos_Combo",parSigTimeCjpsiKS,obsTagOS,parSigEtaMean_OS,RooConst(0.),RooConst(0.),parSigEtaMean_OS,RooConst(0.),RooConst(0.),obsTagSS,parSigEtaMean_SS,RooConst(0.),RooConst(0.),parSigEtaMean_SS,RooConst(0.),RooConst(0.),parSigEtaDeltaProd,SinCoeffCombo::kCType);
  RooBDecay         pdfSigTime_Combo("pdfSigTime_Combo","P_{S}^{l}(t,d|#eta)",obsTime,parSigTimeTau,parSigTimeDeltaG,parSigTimeCosh_Combo,parSigTimeSinh,parSigTimeCos_Combo,parSigTimeSin_Combo,parSigTimeDeltaM,resGauss,RooBDecay::SingleSided);
  
  RooConstVar       parSigEta_same("parSigEta_same","parSigEta_same",0.0859375);
  RooFormulaVar		  parSigTimeCosh_sameGeneration("parSigTimeCosh_sameGeneration","parSigTimeCosh_sameGeneration","1.0 - @0*@1*(1-2.0*@2)",RooArgList(obsTagSS,parSigEtaDeltaProd,parSigEta_same));
  RooFormulaVar		  parSigTimeCos_sameGeneration("parSigTimeCos_sameGeneration","parSigTimeCos_sameGeneration","@0*(@1*(1-2.0*@2)-@3)",RooArgList(parSigTimeCjpsiKS,obsTagSS,parSigEta_same,parSigEtaDeltaProd));
  RooFormulaVar		  parSigTimeSin_sameGeneration("parSigTimeSin_sameGeneration","parSigTimeSin_sameGeneration","-@0*(@1*(1-2.0*@2)-@3)",RooArgList(parSigTimeSin2b,obsTagSS,parSigEta_same,parSigEtaDeltaProd));
  RooBDecay         pdfSigTime_sameGeneration("pdfSigTime_sameGeneration","P_{S}^{l}(t,d|#eta)",obsTime,parSigTimeTau,parSigTimeDeltaG,parSigTimeCosh_sameGeneration,parSigTimeSinh,parSigTimeCos_sameGeneration,parSigTimeSin_sameGeneration,parSigTimeDeltaM,resGauss,RooBDecay::SingleSided);

  RooConstVar       parSigEta_different("parSigEta_different","parSigEta_different",11./24.);
  RooFormulaVar     parSigTimeCosh_differentGeneration("parSigTimeCosh_differentGeneration","parSigTimeCosh_differentGeneration","1.0 - @0*@1*(1-2.0*@2)",RooArgList(obsTagSS,parSigEtaDeltaProd,parSigEta_different));
  RooFormulaVar     parSigTimeCos_differentGeneration("parSigTimeCos_differentGeneration","parSigTimeCos_differentGeneration","@0*(@1*(1-2.0*@2)-@3)",RooArgList(parSigTimeCjpsiKS,obsTagSS,parSigEta_different,parSigEtaDeltaProd));
  RooFormulaVar     parSigTimeSin_differentGeneration("parSigTimeSin_differentGeneration","parSigTimeSin_differentGeneration","-@0*(@1*(1-2.0*@2)-@3)",RooArgList(parSigTimeSin2b,obsTagSS,parSigEta_different,parSigEtaDeltaProd));
  RooBDecay         pdfSigTime_differentGeneration("pdfSigTime_differentGeneration","P_{S}^{l}(t,d|#eta)",obsTime,parSigTimeTau,parSigTimeDeltaG,parSigTimeCosh_differentGeneration,parSigTimeSinh,parSigTimeCos_differentGeneration,parSigTimeSin_differentGeneration,parSigTimeDeltaM,resGauss,RooBDecay::SingleSided);

  // Combination of observables
  RooProdPdf        pdfSig_Combo("pdfSig_Combo","pdfSig_Combo",RooArgList(pdfSigTime_Combo,pdfSigMass));
  RooProdPdf        pdfSig_sameGeneration("pdfSig_sameGeneration","pdfSig_sameGeneration",RooArgList(pdfSigTime_sameGeneration,pdfSigMass));
  RooProdPdf        pdfSig_differentGeneration("pdfSig_differentGeneration","pdfSig_differentGeneration",RooArgList(pdfSigTime_differentGeneration,pdfSigMass));

  // Background
  // Mass model
  RooRealVar              parBkgMassExponent("parBkgMassExponent","Background Mass Exponent",-0.002,-1.,1.,"1/(MeV/c^{2})");
  RooExponential          pdfBkgMass("pdfBkgMass","Background Mass PDF",obsMass,parBkgMassExponent);

  // Decay time model
  RooRealVar              parBkgTimeTau("parBkgTimeTau","parBkgTimeTau",1,0.5,1.5);
  RooDecay                pdfBkgTime("pdfBkgTime","pdfBkgTime",obsTime,parBkgTimeTau,resGauss,RooDecay::SingleSided);

  // Combination of observables
  RooProdPdf              pdfBkg_Combo("pdfBkg_Combo","pdfBkg_Combo",RooArgList(pdfBkgTime,pdfBkgMass));

  // Combining signal and background
  RooRealVar              parSigYield_Combo("parSigYield_Combo","parSigYield_Combo",100000,0,200000);
  RooRealVar              parBkgYield_Combo("parBkgYield_Combo","parBkgYield_Combo",100000,0,200000);
  RooAddPdf               pdf_Combo("pdf_Combo","pdf_Combo",RooArgList(pdfSig_Combo,pdfBkg_Combo),RooArgList(parSigYield_Combo,parBkgYield_Combo));

  RooRealVar              parSigsameTagDecisionYield("parSigsameTagDecisionYield","parSigsameTagDecisionYield",25000,0,50000);
  RooRealVar              parSigdifferentTagDecisionYield("parSigdifferentTagDecisionYield","parSigdifferentTagDecisionYield",25000,0,50000);
  RooConstVar             parBkgsamedifferentTagDecisionRatio("parBkgsamedifferentTagDecisionRatio","parBkgsamedifferentTagDecisionRatio",0.5);
  RooProduct              parBkgsameTagDecisionYield("parBkgsameTagDecisionYield","parBkgsameTagDecisionYield",RooArgList(parBkgYield_Combo,parBkgsamedifferentTagDecisionRatio));
  RooProduct              parBkgdifferentTagDecisionYield("parBkgdifferentTagDecisionYield","parBkgdifferentTagDecisionYield",RooArgList(parBkgYield_Combo,parBkgsamedifferentTagDecisionRatio));

  RooAddPdf               pdf_same("pdf_same","pdf_same",RooArgList(pdfSig_sameGeneration,pdfBkg_Combo),RooArgList(parSigsameTagDecisionYield,parBkgsameTagDecisionYield));
  RooAddPdf               pdf_different("pdf_different","pdf_different",RooArgList(pdfSig_differentGeneration,pdfBkg_Combo),RooArgList(parSigdifferentTagDecisionYield,parBkgdifferentTagDecisionYield));

  RooSimultaneous         pdf_Generation("pdf_Generation","pdf_Generation",obsTagOSSS);
  pdf_Generation.addPdf(pdf_same,"same");
  pdf_Generation.addPdf(pdf_different,"different");

  if (method.EqualTo("g") || method.EqualTo("e")) {
            
    RooWorkspace* ws = new RooWorkspace("ws");
    ws->import(pdf_Generation);
    ws->import(observables);
    ws->defineSet("observables",observables);

    doofit::config::CommonConfig cfg_com("common");
    cfg_com.InitializeOptions(argc, argv);
    
    ToyFactoryStdConfig cfg_tfac("toyfac");
    cfg_tfac.InitializeOptions(cfg_com);
    
    ToyStudyStdConfig cfg_tstudy("toystudy");
    cfg_tstudy.InitializeOptions(cfg_tfac);
    
    cfg_com.CheckHelpFlagAndPrintHelp();
        
    ws->Print();
    
    cfg_tfac.set_workspace(ws);
    
    ToyFactoryStd tfac(cfg_com, cfg_tfac);

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
    
    RooDataSet* data = NULL;

    if (method.EqualTo("g")) {
      for (int i = 0; i < num_toys; ++i) {
        cout  <<  i <<  endl;
        try {
          data = tfac.Generate();
          pdf_Combo.getParameters(*data)->readFromFile(argv[4]);
          RooFitResult* fit_result = pdf_Combo.fitTo(*data,fitting_args);
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