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
  
  if (argc < 5)
  {
    cout  <<  "You need to provide:"  <<  endl;
    cout  <<  "the method (g for generation of toys, e for evaluation of toy study)"  <<  endl;
    cout  <<  "the number of CPUs to fit on"  <<  endl;
    cout  <<  "the name of the config file for the toy production"  <<  endl;
    cout  <<  "the name of the starting values file" <<  endl;
    return 1;
  }

  TString           method;
  if (argc >= 2) {
    method = argv[1];
  }
  int num_cpu;
  if (argc >= 3) {
    TString         StringNumCPU = argv[2];
    num_cpu = StringNumCPU.Atoi();
  }
  
  RooRealVar        obsTime("obsTime","#it{t}",0.,18.,"ps");
  RooRealVar        obsEtaOS("obsEtaOS","#eta_{OS}",0.001,0.501);
  RooRealVar        obsMass("obsMass","#it{m_{J/#kern[-0.3]{#psi} K_{S}}}",5230,5330,"MeV/c^{2}");
  RooCategory       obsTagOS("obsTagOS","Flavour Tag");
  obsTagOS.defineType("B0",1);
  obsTagOS.defineType("B0bar",-1);
  obsTagOS.defineType("untagged",0);
  RooRealVar        obsEtaSSPion("obsEtaSS","#eta_{SS#pi}",0.,0.5);
  RooCategory       obsTagSSPion("obsTagSS","Flavour Tag");
  obsTagSSPion.defineType("B0",1);
  obsTagSSPion.defineType("B0bar",-1);
  obsTagSSPion.defineType("untagged",0);
  
  RooCategory       catTaggedExclOSSSPion("catTag","OS or SSPion tagged");
  catTaggedExclOSSSPion.defineType("OS",1);
  // catTaggedExclOSSSPion.defineType("SSPion",-1);
  // catTaggedExclOSSSPion.defineType("ut",0);
  
  RooArgSet         observables(obsTime,obsMass,obsEtaOS,obsTagOS,obsEtaSSPion,obsTagSSPion,catTaggedExclOSSSPion,"observables");
  
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
  RooRealVar        parSigTimeDeltaM("parSigTimeDeltaM","B0 mixing",0.517);            // PDG2013 (0.510 +- 0.004)/ps
  RooRealVar        parSigTimeDeltaG("parSigTimeDeltaG","#Delta#Gamma",0.);

  // Additional parameters needed for B Decay with CPV
  RooRealVar        parSigTimeSin2b("parSigTimeSin2b","#it{S_{J/#kern[-0.3]{#psi} K_{S}}}",0.7,0.6,0.8);
  RooRealVar        parSigTimeCjpsiKS("parSigTimeCjpsiKS","#it{C_{J/#kern[-0.3]{#psi} K_{S}}}",0,-0.5,0.5);
      
  // Tagging asymmetries
  RooRealVar        parSigEtaDeltaProd("parSigEtaDeltaProd","asymmetry",0.02,-0.5,0.5);
  RooConstVar       parSigTimeDelta("parSigTimeDelta","asymmetry",0.);
  
  // Tagging calibration
  RooConstVar                             parSigEtaP0_OS("parSigEtaP0_OS","parSigEtaP0_OS",0.25);
  RooConstVar                             parSigEtaP1_OS("parSigEtaP1_OS","parSigEtaP1_OS",1.);
  RooConstVar                             parSigEtaDeltaP0_OS("parSigEtaDeltaP0_OS","parSigEtaDeltaP0_OS",0.);
  RooConstVar                             parSigEtaDeltaP1_OS("parSigEtaDeltaP1_OS","parSigEtaDeltaP1_OS",0.);
  RooConstVar                             parSigEtaMean_OS("parSigEtaMean_OS","parSigEtaMean_OS",0.25);
  RooConstVar                             parSigEtaP0_SS("parSigEtaP0_SS","parSigEtaP0_SS",0.);
  RooConstVar                             parSigEtaP1_SS("parSigEtaP1_SS","parSigEtaP1_SS",1.);
  RooConstVar                             parSigEtaMean_SS("parSigEtaMean_SS","parSigEtaMean_SS",0.);
  RooConstVar                             parSigEtaDeltaP0_SS("parSigEtaDeltaP0_SS","parSigEtaDeltaP0_SS",0.);
  RooConstVar                             parSigEtaDeltaP1_SS("parSigEtaDeltaP1_SS","parSigEtaDeltaP1_SS",0.);
  SingleMistagCalibrationWithAsymmetries  parSigTimeOmega_OS_Bd("parSigTimeOmega_OS_Bd","tagging calibration of OS taggers",obsEtaOS,parSigEtaP0_OS,parSigEtaP1_OS,parSigEtaDeltaP0_OS,parSigEtaDeltaP1_OS,parSigEtaMean_OS,SingleMistagCalibrationWithAsymmetries::kBdType);
  SingleMistagCalibrationWithAsymmetries  parSigTimeOmega_OS_Bdb("parSigTimeOmega_OS_Bdb","tagging calibration of OS taggers",obsEtaOS,parSigEtaP0_OS,parSigEtaP1_OS,parSigEtaDeltaP0_OS,parSigEtaDeltaP1_OS,parSigEtaMean_OS,SingleMistagCalibrationWithAsymmetries::kBdbType);

  // Decay Time PDF
  // RooBDecay params
  RooConstVar             parSigTimeSinh("parSigTimeSinh","Sh_{f}",0.0);
  CoshCoeff               parSigTimeCosh_OS("parSigTimeCosh_OS","cosh coefficient OS",parSigTimeOmega_OS_Bd,parSigTimeOmega_OS_Bdb,parSigEtaDeltaProd,obsTagOS);
  SinCoeffWithProdAsymm   parSigTimeSin_OS("parSigTimeSin_OS",parSigTimeSin2b,parSigTimeOmega_OS_Bd,parSigTimeOmega_OS_Bdb,obsTagOS,parSigEtaDeltaProd,SinCoeffWithProdAsymm::kSType);
  SinCoeffWithProdAsymm   parSigTimeCos_OS("parSigTimeCos_OS",parSigTimeCjpsiKS,parSigTimeOmega_OS_Bd,parSigTimeOmega_OS_Bdb,obsTagOS,parSigEtaDeltaProd,SinCoeffWithProdAsymm::kCType);

  // RooFormulaVar           parSigTimeCosh_OS("parSigTimeCosh_OS","cosh coefficient OS","1.0 - @0*@1*(1.0 - 2.0*@2)",RooArgList(parSigEtaDeltaProd,obsTagOS,obsEtaOS));
  // RooFormulaVar           parSigTimeSin_OS("parSigTimeSin_OS","sin coefficient OS","-@0*(@1*(1.0 - @2 - @3)- @4*(1.0 - @1*(@2 - @3)))",RooArgList(parSigTimeSin2b,obsTagOS,parSigTimeOmega_OS_Bd,parSigTimeOmega_OS_Bdb,parSigEtaDeltaProd));
  // RooFormulaVar           parSigTimeCos_OS("parSigTimeCos_OS","cos coefficient OS"," @0*(@1*(1.0 - @2 - @3)- @4*(1.0 - @1*(@2 - @3)))",RooArgList(parSigTimeCjpsiKS,obsTagOS,parSigTimeOmega_OS_Bd,parSigTimeOmega_OS_Bdb,parSigEtaDeltaProd));

  
  CoshCoeffCombo    parSigTimeCosh("parSigTimeCosh",obsTagOS,parSigEtaP0_OS,parSigEtaP1_OS,parSigEtaMean_OS,obsEtaOS,parSigEtaDeltaP0_OS,parSigEtaDeltaP1_OS,obsTagSSPion,parSigEtaP0_SS,parSigEtaP1_SS,parSigEtaMean_SS,obsEtaSSPion,parSigEtaDeltaP0_SS,parSigEtaDeltaP1_SS,parSigEtaDeltaProd);
  SinCoeffCombo     parSigTimeSin("parSigTimeSin",parSigTimeSin2b,obsTagOS,parSigEtaP0_OS,parSigEtaP1_OS,parSigEtaMean_OS,obsEtaOS,parSigEtaDeltaP0_OS,parSigEtaDeltaP1_OS,obsTagSSPion,parSigEtaP0_SS,parSigEtaP1_SS,parSigEtaMean_SS,obsEtaSSPion,parSigEtaDeltaP0_SS,parSigEtaDeltaP1_SS,parSigEtaDeltaProd,SinCoeffCombo::kSType);
  SinCoeffCombo     parSigTimeCos("parSigTimeCos",parSigTimeCjpsiKS,obsTagOS,parSigEtaP0_OS,parSigEtaP1_OS,parSigEtaMean_OS,obsEtaOS,parSigEtaDeltaP0_OS,parSigEtaDeltaP1_OS,obsTagSSPion,parSigEtaP0_SS,parSigEtaP1_SS,parSigEtaMean_SS,obsEtaSSPion,parSigEtaDeltaP0_SS,parSigEtaDeltaP1_SS,parSigEtaDeltaProd,SinCoeffCombo::kCType);

  // Eta PDF
  TH1D*                   histSigEta_OS = new TH1D("histSigEta_OS","histogram of OS tagger",100,obsEtaOS.getMin(),obsEtaOS.getMax());
  for (int i = 0; i <= 100; ++i)
  {
    histSigEta_OS->SetBinContent(i,1);
  }
  RooDataHist             datahistSigEta_OS("datahistSigEta_OS","data histogram of OS tagger",obsEtaOS,histSigEta_OS);
  RooHistPdf              pdfSigEta_OS("pdfSigEta_OS","Signal Eta OS PDF",obsEtaOS,datahistSigEta_OS);

  TH1D*                   histSigEta_SSPion = new TH1D("histSigEta_SSPion","histogram of SSPion tagger",100,obsEtaSSPion.getMin(),obsEtaSSPion.getMax());
  for (int i = 0; i <= 100; ++i)
  {
    histSigEta_SSPion->SetBinContent(i,1);
  }
  RooDataHist             datahistSigEta_SSPion("datahistSigEta_SSPion","data histogram of SSPion tagger",obsEtaSSPion,histSigEta_SSPion);
  RooHistPdf              pdfSigEta_SSPion("pdfSigEta_SSPion","Signal Eta SSPion PDF",obsEtaSSPion,datahistSigEta_SSPion);

  RooRealVar              untaggedvalue("untaggedvalue","untaggedvalue",0.5);
  BiasDelta               pdfSigEta_OS_ut("pdfSigEta_OS_ut","untagged signal eta distribution",obsEtaOS,untaggedvalue);
  
  RooRealVar              parSigEtaOSTaggingFraction("parSigEtaOSTaggingFraction","fraction of OS tagged events in combined sample",0.5,0.,1.);
  RooAddPdf               pdfSigEta_OS_Combo("pdfSigEta_OS_Combo","OS eta PDF including untagged distribution",RooArgList(pdfSigEta_OS,pdfSigEta_OS_ut),parSigEtaOSTaggingFraction);

  RooDecay                pdfSigTimeDecay_OS("pdfSigTimeDecay_OS","P_{S}^{OS}(t)",obsTime,parSigTimeTau,resGauss,RooDecay::SingleSided);
  // RooBCPGenDecay          pdfSigTime_OS("pdfSigTime_OS","P_{S}^{l}(t,d|#eta)",obsTime,obsTagOS,parSigTimeTau,parSigTimeDeltaM,obsEtaOS,parSigTimeCjpsiKS,parSigTimeSin2b,parSigTimeDelta,parSigEtaDeltaProd,resGauss,RooBCPGenDecay::SingleSided);
  RooBDecay               pdfSigTime_OS("pdfSigTime_OS","P_{S}^{l}(t,d|#eta)",obsTime,parSigTimeTau,parSigTimeDeltaG,parSigTimeCosh_OS,parSigTimeSinh,parSigTimeCos_OS,parSigTimeSin_OS,parSigTimeDeltaM,resGauss,RooBDecay::SingleSided);
  RooProdPdf              pdfSigTimeCond_OS("pdfSigTimeCond_OS","pdfSigTimeCond_OS",RooArgList(pdfSigEta_OS),Conditional(pdfSigTime_OS,RooArgSet(obsTime,obsTagOS)));
  RooBDecay               pdfSigTime_Combo("pdfSigTime_Combo","P_{S}^{l}(t,d|#eta)",obsTime,parSigTimeTau,parSigTimeDeltaG,parSigTimeCosh,parSigTimeSinh,parSigTimeCos,parSigTimeSin,parSigTimeDeltaM,resGauss,RooBDecay::SingleSided);
  RooProdPdf              pdfSigTimeCond_Combo("pdfSigTimeCond_Combo","P_{S}^{l}(t,d|#eta)",RooArgList(pdfSigEta_OS_Combo,pdfSigEta_SSPion),Conditional(pdfSigTime_Combo,RooArgSet(obsTime,obsTagOS,obsTagSSPion)));

  RooDecay                pdfSigTime_ut("pdfSigTime_ut","pdfSigTime_ut",obsTime,parSigTimeTau,resGauss,RooDecay::SingleSided);
  
  // Combination of observables
  RooProdPdf              pdfSig_OS("pdfSig_OS","pdfSig_OS",RooArgList(pdfSigTimeCond_OS,pdfSigMass));
  RooProdPdf              pdfSig_SSPion("pdfSig_SSPion","pdfSig_SSPion",RooArgList(pdfSigTimeCond_Combo,pdfSigMass));
  RooProdPdf              pdfSig_ut("pdfSig_ut","pdfSig_ut",RooArgList(pdfSigTime_ut,pdfSigMass));

  // Background
  // Mass model
  RooRealVar              parBkgMassExponent("parBkgMassExponent","Background Mass Exponent",-0.002,-1.,1.,"1/(MeV/c^{2})");
  RooExponential          pdfBkgMass("pdfBkgMass","Background Mass PDF",obsMass,parBkgMassExponent);

  // Decay time model
  RooRealVar              parBkgTimeTau("parBkgTimeTau","parBkgTimeTau",1,0.5,1.5);
  RooDecay                pdfBkgTime("pdfBkgTime","pdfBkgTime",obsTime,parBkgTimeTau,resGauss,RooDecay::SingleSided);

  // Eta model
  TH1D*                   histBkgEta_OS = new TH1D("histBkgEta_OS","histogram of OS tagger",100,obsEtaOS.getMin(),obsEtaOS.getMax());
  for (int i = 0; i <= 100; ++i)
  {
    histBkgEta_OS->SetBinContent(i,1);
  }
  RooDataHist             datahistBkgEta_OS("datahistBkgEta_OS","data histogram of OS tagger",obsEtaOS,histBkgEta_OS);
  RooHistPdf              pdfBkgEta_OS("pdfBkgEta_OS","Background Eta OS PDF",obsEtaOS,datahistBkgEta_OS);

  TH1D*                   histBkgEta_SSPion = new TH1D("histBkgEta_SSPion","histogram of SSPion tagger",100,obsEtaSSPion.getMin(),obsEtaSSPion.getMax());
  for (int i = 0; i <= 100; ++i)
  {
    histBkgEta_SSPion->SetBinContent(i,1);
  }
  RooDataHist             datahistBkgEta_SSPion("datahistBkgEta_SSPion","data histogram of SSPion tagger",obsEtaSSPion,histBkgEta_SSPion);
  RooHistPdf              pdfBkgEta_SSPion("pdfBkgEta_SSPion","Background Eta SSPion PDF",obsEtaSSPion,datahistBkgEta_SSPion);

  BiasDelta               pdfBkgEta_OS_ut("pdfBkgEta_OS_ut","untagged background eta distribution",obsEtaOS,untaggedvalue);
  
  RooRealVar              parBkgEtaOSTaggingFraction("parBkgEtaOSTaggingFraction","fraction of OS tagged events in combined sample",0.5,0.,1.);
  RooAddPdf               pdfBkgEta_OS_Combo("pdfBkgEta_OS_Combo","OS eta PDF including untagged distribution",RooArgList(pdfBkgEta_OS,pdfBkgEta_OS_ut),parBkgEtaOSTaggingFraction);

  // Combination of observables
  RooProdPdf              pdfBkg_OS("pdfBkg_OS","pdfBkg_OS",RooArgList(pdfBkgTime,pdfBkgMass,pdfBkgEta_OS));
  RooProdPdf              pdfBkg_SSPion("pdfBkg_SSPion","pdfBkg_SSPion",RooArgList(pdfBkgTime,pdfBkgMass,pdfBkgEta_OS_Combo,pdfBkgEta_SSPion));
  RooProdPdf              pdfBkg_ut("pdfBkg_ut","pdfBkg_ut",RooArgList(pdfBkgTime,pdfBkgMass));

  // Combining signal and background
  RooRealVar              parSigYield_OS("parSigYield_OS","parSigYield_OS",100000,0,200000);
  RooRealVar              parSigYield_SSPion("parSigYield_SSPion","parSigYield_SSPion",100000,0,200000);
  RooRealVar              parSigYield_ut("parSigYield_ut","parSigYield_ut",100000,0,200000);

  RooRealVar              parBkgYield_OS("parBkgYield_OS","parBkgYield_OS",100000,0,200000);
  RooRealVar              parBkgYield_SSPion("parBkgYield_SSPion","parBkgYield_SSPion",100000,0,200000);
  RooRealVar              parBkgYield_ut("parBkgYield_ut","parBkgYield_ut",100000,0,200000);

  RooAddPdf               pdf_OS("pdf_OS","pdf_OS",RooArgList(pdfSig_OS,pdfBkg_OS),RooArgList(parSigYield_OS,parBkgYield_OS));
  RooAddPdf               pdf_SSPion("pdf_SSPion","pdf_SSPion",RooArgList(pdfSig_SSPion,pdfBkg_SSPion),RooArgList(parSigYield_SSPion,parBkgYield_SSPion));
  RooAddPdf               pdf_ut("pdf_ut","pdf_ut",RooArgList(pdfSig_ut,pdfBkg_ut),RooArgList(parSigYield_ut,parBkgYield_ut));

  RooSimultaneous         simpdf("simpdf","simpdf",catTaggedExclOSSSPion);
  simpdf.addPdf(pdf_OS,"OS");
  simpdf.addPdf(pdf_SSPion,"SSPion");
  simpdf.addPdf(pdf_ut,"ut");

  // RooAbsPdf*              pdf = &simpdf;

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

    ToyConfig     cfg_cptoymc;
    cfg_cptoymc.load(argv[3]);
    ToyGenerator  cptoymc(cfg_cptoymc);
    
    RooDataSet* data = NULL;

    if (method.EqualTo("g")) {
      for (int i = 0; i < 1000; ++i) {
        cout  <<  i <<  endl;
        try {
          TTree       tree("ToyMCTreetree","Tree of generation");
          cptoymc.GenerateToy(tree);
          data = new RooDataSet("data","Toy MC data",&tree, observables);
          data->Print();
          pdfSigTimeCond_OS.getParameters(*data)->readFromFile(argv[4]);
          RooFitResult* fit_result = pdfSigTimeCond_OS.fitTo(*data,fitting_args);
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