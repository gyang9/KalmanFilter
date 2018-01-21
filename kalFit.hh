#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <iomanip>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <sstream>
#include <TList.h>
#include <TROOT.h>
#include <TCanvas.h>
#include <TApplication.h>
#include <TH1.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TFile.h>
#include <TRint.h>
#include <TH2.h>
#include <TFormula.h>
#include <TF1.h>
#include <TF2.h>
#include <TMath.h>
#include <Math/DistFunc.h>
#include <TLine.h>
#include <TTree.h>
#include <TGraph.h>
#include <TRandom.h>
#include <TRandom3.h>
#include <TGraphErrors.h>
#include <TVirtualFFT.h>
#include <TFoamIntegrand.h>
#include <TMatrixD.h>
#include <TVectorT.h>
#include <TDecompChol.h>
#include <RooFit.h>
#include "RooGlobalFunc.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooDataSet.h"
#include "RooArgList.h"
#include "RooArgSet.h"
#include "RooGaussian.h"
#include "RooProdPdf.h"
#include "RooWorkspace.h"
#include "RooMinuit.h"
#include "RooNLLVar.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "RooDataSet.h"
#include "RooExtendPdf.h"
#include "RooChi2Var.h"
#include "RooMinuit.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooRandom.h"
#include <RooMsgService.h>
#include <RooHist.h>
#include <RooTrace.h>
#include <RooCategory.h>
#include "RooConstVar.h"
#include "RooBinning.h"
#include "TStopwatch.h"
#include "TFile.h"
#include "TMinuit.h"

#include "RooFit.h"
#include "RooMinuit.h"
#include "RooFitResult.h"
#include "TMinuit.h"
#include <RooRealVar.h>

using namespace RooFit;
//using namespace std;

class KALFIT : public RooAbsReal {

public:

  KALFIT(const char* name);
  KALFIT (const KALFIT & other, const char* name = 0): RooAbsReal(other,name) {};
  virtual TObject* clone(const char* newname) const {return new KALFIT (*this, newname);};
  virtual ~KALFIT();

  KALFIT (const KALFIT & KALFIT );

  //double GetIntegral(TString hyp);

  void STARTFIT(TH2D* h);

    Double_t FillEv(RooListProxy* _pulls, TH1D* Mmean, TH1D* Mvariance, TH1D* MPE, TH1D* MMCS) const;

    TH1D* Mmean;
    TH1D* Mvariance;
    TH1D* MPE;
    TH1D* MMCS;
    
    Double_t ExtraPull(RooListProxy* _pulls) const;

    TMatrixD* prepareMatrix(TH1D* Mmean, TH1D* Mvariance, TH1D* MPE, TH1D* MMCS) const;

    void setSyst(Double_t syst) ;

    RooFormulaVar* Chi2() ;

    void setData(RooListProxy* _pulls) const;

    void setPull(TH1D* pullvecCV) ;
    void setPullUnc(TH1D* pullvecUnc) ;
    Double_t getPullUnc(Int_t pN) ;

    RooRealVar Par1 ;
    RooRealVar Par2 ;
    RooRealVar Par3 ;
    RooRealVar Par4 ;
    RooRealVar Par5 ;
    RooRealVar Par6 ;
    RooRealVar Par7 ;
    RooRealVar Par8 ;
    RooRealVar Par9 ;
    RooRealVar Par10 ;
    RooRealVar Par11 ;
    RooRealVar Par12 ;
    RooRealVar Par13 ;
    RooRealVar Par14 ;
    RooRealVar Par15 ;
    RooRealVar Par16 ;
    RooRealVar Par17 ;
    RooRealVar Par18 ;
    RooRealVar Par19 ;
    RooRealVar Par20 ;
    RooRealVar Par21 ;
    RooRealVar Par22 ;
    RooRealVar Par23 ;
    RooRealVar Par24 ;
    Double_t _par1;
    Double_t _par2;
    Double_t _par3;
    Double_t _par4;
    Double_t _par5;
    Double_t _par6;
    Double_t _par7;
    Double_t _par8;
    Double_t _par9;
    Double_t _par10;
    Double_t _par11;
    Double_t _par12;
    Double_t _par13;
    Double_t _par14;
    Double_t _par15;
    Double_t _par16;
    Double_t _par17;
    Double_t _par18;
    Double_t _par19;
    Double_t _par20;
    Double_t _par21;
    Double_t _par22;
    Double_t _par23;
    Double_t _par24;

    Double_t getPar(int i) ;
    RooRealVar* getParVar(int i) ;
    RooListProxy* getParVar() ;
    RooListProxy* getPullList() const;
    RooArgList _parlist;
    RooListProxy* _pulls;

    //Double_t statL[200];
    //Double_t meanL[200];
    //Double_t varL[200];
    //Double_t mscL[200];
    Double_t currY;
    Double_t currYl;
    Double_t thetaZ;

  virtual  Double_t evaluate() const ;

private:

  void Init();
  void Delete();
  
};

