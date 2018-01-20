
//
// TODO:
// - evaluate PID using p-value for both hypotheses
//   and calculating e.g. the CLs = pvalue_muon / (1-pvalue_prot)
// 


#include <TMath.h>
#include <TString.h>
#include "TMath.h"

#include "RooArgList.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include <TVectorD.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TFrame.h>
#include "kalFit.hh"

#define DEBUG 0 // 0 --> disable debugging
                // 1 --> enable debugging

//using namespace std;


KALFIT::KALFIT(const char* name )
: RooAbsReal(name, name)
{
  this->Init();
  std::cout<<"starting KALFIT.. "<<std::endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

KALFIT::~KALFIT()
{
  this->Delete();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void KALFIT::Init(){

  _pulls     = new RooListProxy("_pulls","_pulls",this);
  RooRealVar* Par1 = new RooRealVar("s12","par1",TMath::ASin(TMath::Sqrt(0.85))/2.,0,100);
  RooRealVar* Par2 = new RooRealVar("s23","par2",TMath::ASin(TMath::Sqrt(0.95))/2.,0,100);
  RooRealVar* Par3 = new RooRealVar("s13","par3",0.1,0,100);
  RooRealVar* Par4 = new RooRealVar("delta","par4",-1.5,-10,10);
  RooRealVar* Par5 = new RooRealVar("dm21","par5",0.000075,-10,10);
  RooRealVar* Par6 = new RooRealVar("dm32","par6",0.00238,-10,10);
  RooRealVar* Par7 = new RooRealVar("dm31","par7",0.00244,-10,10);
  RooRealVar* Par8 = new RooRealVar("numuX","par8",1,0.,100);
  RooRealVar* Par9 = new RooRealVar("nueX","par9",1,0.,100);
  RooRealVar* Par10 = new RooRealVar("numuSel","par10",1,0.,100);
  RooRealVar* Par11 = new RooRealVar("nueSel","par11",1,0.,100);


Par1->setConstant(false);
Par2->setConstant(false);
Par3->setConstant(false);
Par4->setConstant(false);
Par5->setConstant(false);
Par6->setConstant(false);
Par7->setConstant(false);
Par8->setConstant(false);
Par9->setConstant(false);
Par10->setConstant(false);
Par11->setConstant(false);


_parlist.add(*Par1);
_parlist.add(*Par2);
_parlist.add(*Par3);
_parlist.add(*Par4);
_parlist.add(*Par5);
_parlist.add(*Par6);
_parlist.add(*Par7);
_parlist.add(*Par8);
_parlist.add(*Par9);
_parlist.add(*Par10);
_parlist.add(*Par11);
_pulls->add(_parlist);

  this->addServerList(*_pulls);
}

Double_t KALFIT ::FillEv( RooListProxy* _pulls ) const
{

}

Double_t KALFIT ::ExtraPull (RooListProxy* _pulls) const
{
}


Double_t KALFIT ::evaluate() const
{

Double_t matPart = this->FillEv(_pulls);

Double_t extraPull = this -> ExtraPull (_pulls);
Double_t tot = matPart + extraPull; //If needed, add pull terms here.

return tot;

}

Double_t KALFIT ::getPar(int i) {
(((RooAbsReal*)_pulls->at(i))->getVal());
}

RooRealVar* KALFIT ::getParVar(int i) {
return ((RooRealVar*)_pulls->at(i));
}





//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void KALFIT::Delete(){

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......    


void KALFIT::STARTFIT(TH2D* h){

std::cout<<"starting real KALFIT.. "<<std::endl;

}



