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
  RooRealVar        obsEtaOS("obsEtaOS","#eta_{OS}",0.,0.5);
  RooRealVar        obsMass("obsMass","#it{m_{J/#kern[-0.3]{#psi} K_{S}}}",5200,5400,"MeV/c^{2}");
  RooCategory       obsTagOS("obsTagOS","Flavour Tag");
  obsTagOS.defineType("B0",1);
  obsTagOS.defineType("B0bar",-1);
  RooRealVar        obsEtaSS("obsEtaSS","#eta_{SS#pi}",0.,0.5);
  RooCategory       obsTagSS("obsTagSS","Flavour Tag");
  obsTagSS.defineType("B0",1);
  obsTagSS.defineType("B0bar",-1);
  
  RooCategory       catTaggedOSSSPion("catTaggedOSSSPion","OS or SSPion tagged");
  catTaggedOSSSPion.defineType("OS",1);
  catTaggedOSSSPion.defineType("SS",-1);
  catTaggedOSSSPion.defineType("OSSS",10);
  
  RooArgSet         observables(obsTime,obsMass,obsEtaOS,obsTagOS,obsEtaSS,obsTagSS,catTaggedOSSSPion,"observables");
  
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
  RooConstVar                             parSigEtaP0_SS("parSigEtaP0_SS","parSigEtaP0_SS",0.25);
  RooConstVar                             parSigEtaP1_SS("parSigEtaP1_SS","parSigEtaP1_SS",1.);
  RooConstVar                             parSigEtaMean_SS("parSigEtaMean_SS","parSigEtaMean_SS",0.25);
  RooConstVar                             parSigEtaDeltaP0_SS("parSigEtaDeltaP0_SS","parSigEtaDeltaP0_SS",0.);
  RooConstVar                             parSigEtaDeltaP1_SS("parSigEtaDeltaP1_SS","parSigEtaDeltaP1_SS",0.);
  SingleMistagCalibrationWithAsymmetries  parSigTimeOmega_OS_Bd("parSigTimeOmega_OS_Bd","tagging calibration of OS taggers",obsEtaOS,parSigEtaP0_OS,parSigEtaP1_OS,parSigEtaDeltaP0_OS,parSigEtaDeltaP1_OS,parSigEtaMean_OS,SingleMistagCalibrationWithAsymmetries::kBdType);
  SingleMistagCalibrationWithAsymmetries  parSigTimeOmega_OS_Bdb("parSigTimeOmega_OS_Bdb","tagging calibration of OS taggers",obsEtaOS,parSigEtaP0_OS,parSigEtaP1_OS,parSigEtaDeltaP0_OS,parSigEtaDeltaP1_OS,parSigEtaMean_OS,SingleMistagCalibrationWithAsymmetries::kBdbType);
  SingleMistagCalibrationWithAsymmetries  parSigTimeOmega_SS_Bd("parSigTimeOmega_SS_Bd","tagging calibration of SS taggers",obsEtaSS,parSigEtaP0_SS,parSigEtaP1_SS,parSigEtaDeltaP0_SS,parSigEtaDeltaP1_SS,parSigEtaMean_SS,SingleMistagCalibrationWithAsymmetries::kBdType);
  SingleMistagCalibrationWithAsymmetries  parSigTimeOmega_SS_Bdb("parSigTimeOmega_SS_Bdb","tagging calibration of SS taggers",obsEtaSS,parSigEtaP0_SS,parSigEtaP1_SS,parSigEtaDeltaP0_SS,parSigEtaDeltaP1_SS,parSigEtaMean_SS,SingleMistagCalibrationWithAsymmetries::kBdbType);

  // Decay Time PDF
  // RooBDecay params
  RooConstVar             parSigTimeSinh("parSigTimeSinh","Sh_{f}",0.0);
  CoshCoeff               parSigTimeCosh_OS("parSigTimeCosh_OS","cosh coefficient OS",parSigTimeOmega_OS_Bd,parSigTimeOmega_OS_Bdb,parSigEtaDeltaProd,obsTagOS);
  SinCoeffWithProdAsymm   parSigTimeSin_OS("parSigTimeSin_OS",parSigTimeSin2b,parSigTimeOmega_OS_Bd,parSigTimeOmega_OS_Bdb,obsTagOS,parSigEtaDeltaProd,SinCoeffWithProdAsymm::kSType);
  SinCoeffWithProdAsymm   parSigTimeCos_OS("parSigTimeCos_OS",parSigTimeCjpsiKS,parSigTimeOmega_OS_Bd,parSigTimeOmega_OS_Bdb,obsTagOS,parSigEtaDeltaProd,SinCoeffWithProdAsymm::kCType);

  // RooFormulaVar           parSigTimeCosh_OS("parSigTimeCosh_OS","cosh coefficient OS","1.0 - @0*@1*(1.0 - 2.0*@2)",RooArgList(parSigEtaDeltaProd,obsTagOS,obsEtaOS));
  // RooFormulaVar           parSigTimeSin_OS("parSigTimeSin_OS","sin coefficient OS","-@0*(@1*(1.0 - @2 - @3) - @4*(1.0 - @1*(@2 - @3)))",RooArgList(parSigTimeSin2b,obsTagOS,parSigTimeOmega_OS_Bd,parSigTimeOmega_OS_Bdb,parSigEtaDeltaProd));
  // RooFormulaVar           parSigTimeCos_OS("parSigTimeCos_OS","cos coefficient OS"," @0*(@1*(1.0 - @2 - @3) - @4*(1.0 - @1*(@2 - @3)))",RooArgList(parSigTimeCjpsiKS,obsTagOS,parSigTimeOmega_OS_Bd,parSigTimeOmega_OS_Bdb,parSigEtaDeltaProd));

  CoshCoeff               parSigTimeCosh_SS("parSigTimeCosh_SS","cosh coefficient OS",parSigTimeOmega_SS_Bd,parSigTimeOmega_SS_Bdb,parSigEtaDeltaProd,obsTagSS);
  SinCoeffWithProdAsymm   parSigTimeSin_SS("parSigTimeSin_SS",parSigTimeSin2b,parSigTimeOmega_SS_Bd,parSigTimeOmega_SS_Bdb,obsTagSS,parSigEtaDeltaProd,SinCoeffWithProdAsymm::kSType);
  SinCoeffWithProdAsymm   parSigTimeCos_SS("parSigTimeCos_SS",parSigTimeCjpsiKS,parSigTimeOmega_SS_Bd,parSigTimeOmega_SS_Bdb,obsTagSS,parSigEtaDeltaProd,SinCoeffWithProdAsymm::kCType);  
  
  CoshCoeffCombo    parSigTimeCosh_Combo("parSigTimeCosh_Combo",obsTagOS,parSigEtaP0_OS,parSigEtaP1_OS,parSigEtaMean_OS,obsEtaOS,parSigEtaDeltaP0_OS,parSigEtaDeltaP1_OS,obsTagSS,parSigEtaP0_SS,parSigEtaP1_SS,parSigEtaMean_SS,obsEtaSS,parSigEtaDeltaP0_SS,parSigEtaDeltaP1_SS,parSigEtaDeltaProd);
  SinCoeffCombo     parSigTimeSin_Combo("parSigTimeSin_Combo",parSigTimeSin2b,obsTagOS,parSigEtaP0_OS,parSigEtaP1_OS,parSigEtaMean_OS,obsEtaOS,parSigEtaDeltaP0_OS,parSigEtaDeltaP1_OS,obsTagSS,parSigEtaP0_SS,parSigEtaP1_SS,parSigEtaMean_SS,obsEtaSS,parSigEtaDeltaP0_SS,parSigEtaDeltaP1_SS,parSigEtaDeltaProd,SinCoeffCombo::kSType);
  SinCoeffCombo     parSigTimeCos_Combo("parSigTimeCos_Combo",parSigTimeCjpsiKS,obsTagOS,parSigEtaP0_OS,parSigEtaP1_OS,parSigEtaMean_OS,obsEtaOS,parSigEtaDeltaP0_OS,parSigEtaDeltaP1_OS,obsTagSS,parSigEtaP0_SS,parSigEtaP1_SS,parSigEtaMean_SS,obsEtaSS,parSigEtaDeltaP0_SS,parSigEtaDeltaP1_SS,parSigEtaDeltaProd,SinCoeffCombo::kCType);

  // Eta PDF
  TH1D*                   histSigEta_OS = new TH1D("histSigEta_OS","histogram of OS tagger",100,obsEtaOS.getMin(),obsEtaOS.getMax());
  for (int i = 0; i <= 100; ++i)
  {
    histSigEta_OS->SetBinContent(i,1);
  }
  RooDataHist             datahistSigEta_OS("datahistSigEta_OS","data histogram of OS tagger",obsEtaOS,histSigEta_OS);
  RooHistPdf              pdfSigEta_OS("pdfSigEta_OS","Signal Eta OS PDF",obsEtaOS,datahistSigEta_OS);

  TH1D*                   histSigEta_SS = new TH1D("histSigEta_SS","histogram of SSPion tagger",100,obsEtaSS.getMin(),obsEtaSS.getMax());
  for (int i = 0; i <= 100; ++i)
  {
    histSigEta_SS->SetBinContent(i,1);
  }
  RooDataHist             datahistSigEta_SS("datahistSigEta_SS","data histogram of SSPion tagger",obsEtaSS,histSigEta_SS);
  RooHistPdf              pdfSigEta_SS("pdfSigEta_SS","Signal Eta SSPion PDF",obsEtaSS,datahistSigEta_SS);

  RooRealVar              untaggedvalue("untaggedvalue","untaggedvalue",0.5);
  BiasDelta               pdfSigEta_OS_ut("pdfSigEta_OS_ut","untagged signal eta distribution",obsEtaOS,untaggedvalue);
  
  RooDecay                pdfSigTimeDecay_OS("pdfSigTimeDecay_OS","P_{S}^{OS}(t)",obsTime,parSigTimeTau,resGauss,RooDecay::SingleSided);
  // RooBCPGenDecay          pdfSigTime_OS("pdfSigTime_OS","P_{S}^{l}(t,d|#eta)",obsTime,obsTagOS,parSigTimeTau,parSigTimeDeltaM,obsEtaOS,parSigTimeCjpsiKS,parSigTimeSin2b,parSigTimeDelta,parSigEtaDeltaProd,resGauss,RooBCPGenDecay::SingleSided);
  RooBDecay               pdfSigTime_OS("pdfSigTime_OS","P_{S}^{l}(t,d|#eta)",obsTime,parSigTimeTau,parSigTimeDeltaG,parSigTimeCosh_OS,parSigTimeSinh,parSigTimeCos_OS,parSigTimeSin_OS,parSigTimeDeltaM,resGauss,RooBDecay::SingleSided);
  RooProdPdf              pdfSigTimeCond_OS("pdfSigTimeCond_OS","pdfSigTimeCond_OS",RooArgList(pdfSigEta_OS),Conditional(pdfSigTime_OS,RooArgSet(obsTime,obsTagOS)));
  RooBDecay               pdfSigTime_SS("pdfSigTime_SS","P_{S}^{l}(t,d|#eta)",obsTime,parSigTimeTau,parSigTimeDeltaG,parSigTimeCosh_SS,parSigTimeSinh,parSigTimeCos_SS,parSigTimeSin_SS,parSigTimeDeltaM,resGauss,RooBDecay::SingleSided);
  RooProdPdf              pdfSigTimeCond_SS("pdfSigTimeCond_SS","pdfSigTimeCond_OS",RooArgList(pdfSigEta_SS),Conditional(pdfSigTime_SS,RooArgSet(obsTime,obsTagSS)));
  RooBDecay               pdfSigTime_Combo("pdfSigTime_Combo","P_{S}^{l}(t,d|#eta)",obsTime,parSigTimeTau,parSigTimeDeltaG,parSigTimeCosh_Combo,parSigTimeSinh,parSigTimeCos_Combo,parSigTimeSin_Combo,parSigTimeDeltaM,resGauss,RooBDecay::SingleSided);
  RooProdPdf              pdfSigTimeCond_Combo("pdfSigTimeCond_Combo","P_{S}^{l}(t,d|#eta)",RooArgList(pdfSigEta_OS,pdfSigEta_SS),Conditional(pdfSigTime_Combo,RooArgSet(obsTime,obsTagOS,obsTagSS)));

  RooDecay                pdfSigTime_ut("pdfSigTime_ut","pdfSigTime_ut",obsTime,parSigTimeTau,resGauss,RooDecay::SingleSided);
  
  // Combination of observables
  RooProdPdf              pdfSig_OS("pdfSig_OS","pdfSig_OS",RooArgList(pdfSigTimeCond_OS,pdfSigMass));
  RooProdPdf              pdfSig_SS("pdfSig_SS","pdfSig_SS",RooArgList(pdfSigTimeCond_SS,pdfSigMass));
  RooProdPdf              pdfSig_Combo("pdfSig_Combo","pdfSig_Combo",RooArgList(pdfSigTimeCond_Combo,pdfSigMass));
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

  TH1D*                   histBkgEta_SS = new TH1D("histBkgEta_SS","histogram of SSPion tagger",100,obsEtaSS.getMin(),obsEtaSS.getMax());
  for (int i = 0; i <= 100; ++i)
  {
    histBkgEta_SS->SetBinContent(i,1);
  }
  RooDataHist             datahistBkgEta_SS("datahistBkgEta_SS","data histogram of SSPion tagger",obsEtaSS,histBkgEta_SS);
  RooHistPdf              pdfBkgEta_SS("pdfBkgEta_SS","Background Eta SSPion PDF",obsEtaSS,datahistBkgEta_SS);

  BiasDelta               pdfBkgEta_ut("pdfBkgEta_ut","untagged background eta distribution",obsEtaOS,untaggedvalue);
  
  // Combination of observables
  RooProdPdf              pdfBkg_OS("pdfBkg_OS","pdfBkg_OS",RooArgList(pdfBkgTime,pdfBkgMass,pdfBkgEta_OS));
  RooProdPdf              pdfBkg_SS("pdfBkg_SS","pdfBkg_SS",RooArgList(pdfBkgTime,pdfBkgMass,pdfBkgEta_SS));
  RooProdPdf              pdfBkg_Combo("pdfBkg_Combo","pdfBkg_Combo",RooArgList(pdfBkgTime,pdfBkgMass,pdfBkgEta_OS,pdfBkgEta_SS));
  RooProdPdf              pdfBkg_ut("pdfBkg_ut","pdfBkg_ut",RooArgList(pdfBkgTime,pdfBkgMass));

  // Combining signal and background
  RooRealVar              parSigYield_OS("parSigYield_OS","parSigYield_OS",100000,0,200000);
  RooRealVar              parSigYield_SS("parSigYield_SS","parSigYield_SS",100000,0,200000);
  RooRealVar              parSigYield_Combo("parSigYield_Combo","parSigYield_Combo",100000,0,200000);
  RooRealVar              parSigYield_ut("parSigYield_ut","parSigYield_ut",100000,0,200000);

  RooRealVar              parBkgYield_OS("parBkgYield_OS","parBkgYield_OS",100000,0,200000);
  RooRealVar              parBkgYield_SS("parBkgYield_SS","parBkgYield_SS",100000,0,200000);
  RooRealVar              parBkgYield_Combo("parBkgYield_Combo","parBkgYield_Combo",100000,0,200000);
  RooRealVar              parBkgYield_ut("parBkgYield_ut","parBkgYield_ut",100000,0,200000);

  RooAddPdf               pdf_OS("pdf_OS","pdf_OS",RooArgList(pdfSig_OS,pdfBkg_OS),RooArgList(parSigYield_OS,parBkgYield_OS));
  RooAddPdf               pdf_SS("pdf_SS","pdf_SS",RooArgList(pdfSig_SS,pdfBkg_SS),RooArgList(parSigYield_SS,parBkgYield_SS));
  RooAddPdf               pdf_Combo("pdf_Combo","pdf_Combo",RooArgList(pdfSig_Combo,pdfBkg_Combo),RooArgList(parSigYield_Combo,parBkgYield_Combo));
  RooAddPdf               pdf_ut("pdf_ut","pdf_ut",RooArgList(pdfSig_ut,pdfBkg_ut),RooArgList(parSigYield_ut,parBkgYield_ut));

  RooSimultaneous         simpdf("simpdf","simpdf",catTaggedOSSSPion);
  simpdf.addPdf(pdf_OS,"OS");
  simpdf.addPdf(pdf_Combo,"OSSS");
  simpdf.addPdf(pdf_SS,"SS");

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
          simpdf.getParameters(*data)->readFromFile(argv[6]);
          RooFitResult* fit_result = simpdf.fitTo(*data,fitting_args);
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