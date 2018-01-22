
//
// fitter for B-field in 3DST
// Author: Guang Yang 


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

//#define DEBUG 0 // 0 --> disable debugging
                // 1 --> enable debugging

using namespace std;


KALFIT::KALFIT(const char* name )
: RooAbsReal(name, name)
{

  _pulls     = new RooListProxy("_pulls","_pulls",this);
  RooRealVar* Par1 = new RooRealVar("par1","par1",0,-TMath::Pi(),TMath::Pi());
  RooRealVar* Par2 = new RooRealVar("par2","par2",0,-30,30);
  RooRealVar* Par3 = new RooRealVar("par3","par3",0.1,0,100);
  RooRealVar* Par4 = new RooRealVar("par4","par4",-1.5,-10,10);
  RooRealVar* Par5 = new RooRealVar("par5","par5",0.000075,-10,10);
  RooRealVar* Par6 = new RooRealVar("par6","par6",0.00238,-10,10);
  RooRealVar* Par7 = new RooRealVar("par7","par7",0.00244,-10,10);
  RooRealVar* Par8 = new RooRealVar("par8","par8",1,0.,100);
  RooRealVar* Par9 = new RooRealVar("par9","par9",1,0.,100);
  RooRealVar* Par10 = new RooRealVar("par10","par10",1,0.,100);
  RooRealVar* Par11 = new RooRealVar("par11","par11",1,0.,100);


Par1->setConstant(false);
Par2->setConstant(false);
Par3->setConstant(true);
Par4->setConstant(true);
Par5->setConstant(true);
Par6->setConstant(true);
Par7->setConstant(true);
Par8->setConstant(true);
Par9->setConstant(true);
Par10->setConstant(true);
Par11->setConstant(true);


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
  //this->Init();
  std::cout<<"starting KALFIT.. "<<std::endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

KALFIT::~KALFIT()
{
  //this->Delete();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void KALFIT::Init(){

  _pulls     = new RooListProxy("_pulls","_pulls",this);
  RooRealVar* Par1 = new RooRealVar("par1","par1",0,-TMath::Pi(),TMath::Pi());
  RooRealVar* Par2 = new RooRealVar("par2","par2",0,-30,30);
  RooRealVar* Par3 = new RooRealVar("par3","par3",0.1,0,100);
  RooRealVar* Par4 = new RooRealVar("par4","par4",-1.5,-10,10);
  RooRealVar* Par5 = new RooRealVar("par5","par5",0.000075,-10,10);
  RooRealVar* Par6 = new RooRealVar("par6","par6",0.00238,-10,10);
  RooRealVar* Par7 = new RooRealVar("par7","par7",0.00244,-10,10);
  RooRealVar* Par8 = new RooRealVar("par8","par8",1,0.,100);
  RooRealVar* Par9 = new RooRealVar("par9","par9",1,0.,100);
  RooRealVar* Par10 = new RooRealVar("par10","par10",1,0.,100);
  RooRealVar* Par11 = new RooRealVar("par11","par11",1,0.,100);


Par1->setConstant(false);
Par2->setConstant(false);
Par3->setConstant(true);
Par4->setConstant(true);
Par5->setConstant(true);
Par6->setConstant(true);
Par7->setConstant(true);
Par8->setConstant(true);
Par9->setConstant(true);
Par10->setConstant(true);
Par11->setConstant(true);


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


TMatrixD* KALFIT::prepareMatrix(TH1D* Mmean, TH1D* Mvariance, TH1D* MPE, TH1D* MMCS) const
{
//Int_t nbin = Mmean->GetNbinsX();
Int_t nbin = 10;
/*
Int_t nbin = 0;
for(Int_t loopBB = 0; loopBB<Mmean->GetNbinsX(); loopBB++){
if(Mmean->GetBinContent(loopBB+1)==0 && Mmean->GetBinContent(loopBB+2)==0&& Mmean->GetBinContent(loopBB+3)==0 && Mmean->GetBinContent(loopBB+4)==0 ) nbin = loopBB; break;
}
*/
TMatrixD* convMat = new TMatrixD(nbin,nbin);

Int_t startingPP = 0;
   for(Int_t i=0;i<nbin;i++){
        if(Mmean->GetBinContent(i+1)!=0){
        startingPP = i;
                }
        }

for(Int_t i=0;i<nbin;i++){
(*convMat)(i,i) = MPE->GetBinContent(startingPP+i+1);
}
    for(Int_t i=0;i<nbin;i++){
      for(Int_t j=0;j<nbin; j++){
        (*convMat)(i,j) += Mvariance->GetBinContent(startingPP+i+1)*Mvariance->GetBinContent(startingPP+j+1);
	        (*convMat)(i,j) += MMCS->GetBinContent(startingPP+i+1)*MMCS->GetBinContent(startingPP+j+1);
                                }
                             }
return convMat;
}


Double_t KALFIT ::FillEv( RooListProxy* _pulls ) const 
{
TH1D* Mmean     =this->GetMmean();
TH1D* Mvariance =this->GetMvariance();
TH1D* MPE       =this->GetMPE();
TH1D* MMCS      =this->GetMMCS();

std::cout<<"testing input for fit: "<<Mmean->Integral()<<" "<<Mvariance->Integral()<<" "<<MPE->Integral()<<" "<<MMCS->Integral()<<std::endl;
Int_t nbin = 10;
/*
for(Int_t loopBB = 0; loopBB<Mmean->GetNbinsX(); loopBB++){
if(Mmean->GetBinContent(loopBB+1)==0 && Mmean->GetBinContent(loopBB+2)==0 && Mmean->GetBinContent(loopBB+3)==0 && Mmean->GetBinContent(loopBB+4)==0 ) nbin = loopBB; break;
}
*/
std::cout<<"number of bins for fit: "<<nbin<<std::endl;
TVectorD* fVec = new TVectorD(nbin);
TVectorD* fData = new TVectorD(nbin);

Int_t startingPP = 0;
   for(Int_t i=0;i<nbin;i++){
        if(Mmean->GetBinContent(i+1)!=0){
	startingPP = i;
		}
	}
	

   for(Int_t i=0;i<nbin;i++){ 
	(*fVec)[i] = Mmean->GetBinContent(startingPP+1) + TMath::Tan(((RooAbsReal*)_pulls->at(0))->getVal())* 10*i + ((RooAbsReal*)_pulls->at(1))->getVal()*i;
	(*fData)[i] = Mmean->GetBinContent(startingPP+1);
	}

   fVec->Print();
   fData->Print();

TMatrixD* covMat = this->prepareMatrix(Mmean,Mvariance,MPE,MMCS);

   for(Int_t i=0;i<nbin; i++){
       (*fVec)[i] -= (*fData)[i];
         if( (*covMat)(i,i) ==0 ){(*covMat)(i,i) = 1000000000;}
     }

   covMat->Invert();

   TVectorD mulVec(*fVec);
   mulVec *= (*covMat);

   Double_t fResult = TMath::Abs(mulVec*(*fVec));

   std::cout<<"finished filling in.. "<<fResult<<std::endl;
   return fResult;

}

Double_t KALFIT ::ExtraPull (RooListProxy* _pulls) const
{
}


Double_t KALFIT ::evaluate() const 
{

Double_t matPart = this -> FillEv(_pulls);

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

RooListProxy* KALFIT::getPullList() const
{
return _pulls;
}


void KALFIT::SetHistograms(TH1D* MmeanS, TH1D* MvarianceS, TH1D* MPES, TH1D* MMCSS)
{
Mmean = MmeanS;
Mvariance = MvarianceS;
MPE = MPES;
MMCS = MMCSS;
}

TH1D* KALFIT::GetMmean() const{
return Mmean;
}
TH1D* KALFIT::GetMvariance() const{
return Mvariance;
}
TH1D* KALFIT::GetMPE() const{
return MPE;
}
TH1D* KALFIT::GetMMCS() const{
return MMCS;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void KALFIT::Delete(){

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......    


void KALFIT::STARTFIT(TH2D* h){

std::cout<<"starting real KALFIT.. "<<std::endl;

}

int main(int argc, char* argv[])
{

int event;
double hitLocation[3],hitPE[3],hitT[3],adc[3],loadc[3],Q[3],hitLowQ[3];
int separateS=0,separateF=0;
std::string track1Name;
std::string track2Name;
int n3Dcase = 0;
double trueCos,trueC;
double true3Mom[3];
double trueLen,trueL;
double trueMom,trueM;
double trueCos_muon,trueC_muon;
double trueLen_muon,trueL_muon;
double trueMom_muon,trueM_muon,true3Mom_muon[3],true3M_muon;
double trueCos_pion,trueC_pion;
double trueLen_pion,trueL_pion;
double trueMom_pion,trueM_pion,true3Mom_pion[3],true3M_pion;
double trueCos_proton,trueC_proton;
double trueLen_proton,trueL_proton;
double trueMom_proton,trueM_proton,true3Mom_proton[3],true3M_proton;

int prim, PDG;
double ener;
int muonID = 0;
std::cout<<"now doing sample "<<atoi(argv[1])<< " and distance cut for angular resolution is: "<<atof(argv[2])<<std::endl;

//TFile file(Form("/home/gyang/work/dune-ndx/spack/var/spack/stage/edep-sim-master-rpl36hsf6e7fnmtgw4pd74gyffybpekl/edep-sim/SimChain/canDeleteData/testEvent_particleGun2000MeVMuon_2018_shift_sample24.root"));
TFile file(Form("testEvent_particleGun1000MeVMuon_noMCS_30cone_sample2.root",argv[1]));
TTree* c = (TTree*)file.Get("EDepSimTree");
c->SetBranchAddress("event",&event);
c->SetBranchAddress("hitLocation",&hitLocation);
c->SetBranchAddress("hitPE",&hitPE);
c->SetBranchAddress("hitT",&hitT);
c->SetBranchAddress("hitADC",&adc);
c->SetBranchAddress("hitLowADC",&loadc);
c->SetBranchAddress("hitQ",&Q);
c->SetBranchAddress("hitLowQ",&hitLowQ);
c->SetBranchAddress("hitPrim",&prim);
c->SetBranchAddress("hitPDG",&PDG);
c->SetBranchAddress("hitE",&ener);
c->SetBranchAddress("trueLen",&trueLen);
c->SetBranchAddress("trueMom",&trueMom);
c->SetBranchAddress("true3Mom",&true3Mom);
c->SetBranchAddress("trueCos",&trueCos);

Int_t nevent = c->GetEntries();

Int_t nnevent =200;
TH2F* hist2D_XY_Q[nnevent];
TH2F* hist2D_XZ_Q[nnevent];
TH2F* hist2D_YZ_Q[nnevent];
TH2F* hist2D_XY_PE[nnevent];
TH2F* hist2D_XZ_PE[nnevent];
TH2F* hist2D_YZ_PE[nnevent];
TH2F* hist2D_XY_ADC[nnevent];
TH2F* hist2D_XZ_ADC[nnevent];
TH2F* hist2D_YZ_ADC[nnevent];

for(Int_t i=0;i<nnevent;i++){
hist2D_XY_Q[i] = new TH2F("","",240,0,2400,240,0,2400);
hist2D_XZ_Q[i] = new TH2F("","",240,0,2400,200,0,2000);
hist2D_YZ_Q[i] = new TH2F("","",240,0,2400,200,0,2000);

hist2D_XY_PE[i] = new TH2F("","",240,0,2400,240,0,2400);
hist2D_XZ_PE[i] = new TH2F("","",240,0,2400,200,0,2000);
hist2D_YZ_PE[i] = new TH2F("","",240,0,2400,200,0,2000);

hist2D_XY_ADC[i] = new TH2F("","",240,0,2400,240,0,2400);
hist2D_XZ_ADC[i] = new TH2F("","",240,0,2400,200,0,2000);
hist2D_YZ_ADC[i] = new TH2F("","",240,0,2400,200,0,2000);
}

TH2D* hist2D_XY_E[3][nnevent];
TH2D* hist2D_XZ_E[3][nnevent];
TH2D* hist2D_YZ_E[3][nnevent];
for(Int_t i=0;i<3;i++){
for(Int_t j=0;j<nnevent;j++){
hist2D_XY_E[i][j] = new TH2D("","",240,0,2400,240,0,2400);
hist2D_YZ_E[i][j] = new TH2D("","",240,0,2400,200,0,2000);
hist2D_XZ_E[i][j] = new TH2D("","",240,0,2400,200,0,2000);
}
}

TH2D* muon2DprepCos[nnevent];
for(Int_t j=0;j<nnevent;j++){
muon2DprepCos[j] = new TH2D("","",3,0,3,240,0,2400);
}

int eventS=-1;
int eventE=-2;
int seeProton =0;
bool newEvent=false;
int list[100000]={};
double happened[10000]={};
int initA=1;
double hitInit[3]={};
double lastHenergy = 0;
double hardEnergy = 0;
double hardColl = 0;
bool doLoop = false;
KALFIT* rep;

for(Int_t ii=0;ii<nevent;ii++){

c->GetEntry(ii);

event -= atoi(argv[1])*200;

eventS = event;

if(eventS != eventE){newEvent = true;}
else {newEvent = false;}
if(newEvent) {
std::cout<<"event "<<event<<std::endl;
std::cout<<hitLocation[0]<<" "<<hitLocation[1]<<" "<<hitLocation[2]<<std::endl;
}
if( !newEvent ){
std::cout<<"... "<<hitLocation[1]<<" "<<hitLocation[2]<<" "<<hitPE[1]+hitPE[2]<<std::endl;
hist2D_XY_Q[event]->Fill(hitLocation[0],hitLocation[1],Q[0]+Q[1]);
hist2D_XZ_Q[event]->Fill(hitLocation[0],hitLocation[2],Q[0]+Q[2]);
hist2D_YZ_Q[event]->Fill(hitLocation[1],hitLocation[2],Q[1]+Q[2]);
hist2D_XY_PE[event]->Fill(hitLocation[0],hitLocation[1],hitPE[0]+hitPE[1]);
hist2D_XZ_PE[event]->Fill(hitLocation[0],hitLocation[2],hitPE[0]+hitPE[2]);
hist2D_YZ_PE[event]->Fill(hitLocation[1],hitLocation[2],hitPE[1]+hitPE[2]);
hist2D_XY_ADC[event]->Fill(hitLocation[0],hitLocation[1],adc[0]+adc[1]);
hist2D_XZ_ADC[event]->Fill(hitLocation[0],hitLocation[2],adc[0]+adc[2]);
hist2D_YZ_ADC[event]->Fill(hitLocation[1],hitLocation[2],adc[1]+adc[2]);

if((PDG>100 && PDG<10000) || (PDG<-100 && PDG>-10000)){
if(PDG>0){if(happened[PDG]==1){doLoop==false; } else doLoop == true;}
else{if(happened[-PDG+1]==1){doLoop==false; } else doLoop == true;}
if(lastHenergy != PDG && doLoop) {hardEnergy += trueMom;}
lastHenergy = PDG;
if(PDG>0)happened[PDG]=1;
else happened[-PDG+1]=1;
}
if((PDG>100 && PDG<10000) || (PDG<-100 && PDG>-10000)){
hardColl += ener;
}
if(PDG == 13 || PDG == -13 ){

if(initA==0 && TMath::Sqrt(TMath::Power(hitLocation[0]-hitInit[0],2) + TMath::Power(hitLocation[1]-hitInit[1],2) + TMath::Power(hitLocation[2]-hitInit[2],2))<atof(argv[2]) ){

hist2D_XY_E[0][event]->Fill(hitLocation[0],hitLocation[1],ener);
hist2D_XZ_E[0][event]->Fill(hitLocation[0],hitLocation[2],ener);
hist2D_YZ_E[0][event]->Fill(hitLocation[1],hitLocation[2],ener);
muon2DprepCos[event]->Fill(0.,hitLocation[0],(hitPE[0]+hitPE[1]+hitPE[2])/2);
muon2DprepCos[event]->Fill(1.,hitLocation[1],(hitPE[0]+hitPE[1]+hitPE[2])/2);
muon2DprepCos[event]->Fill(2.,hitLocation[2],(hitPE[0]+hitPE[1]+hitPE[2])/2);

trueCos_muon = trueCos;
trueLen_muon = trueLen;
trueMom_muon = trueMom;
true3Mom_muon[0] = true3Mom[0];
true3Mom_muon[1] = true3Mom[1];
true3Mom_muon[2] = true3Mom[2];

}
if(initA==1) {
hitInit[0]=hitLocation[0];
hitInit[1]=hitLocation[1];
hitInit[2]=hitLocation[2];
initA=0;
}
}
if(PDG == 211 || PDG == -211 ){
hist2D_XY_E[1][event]->Fill(hitLocation[0],hitLocation[1],ener);
hist2D_XZ_E[1][event]->Fill(hitLocation[0],hitLocation[2],ener);
hist2D_YZ_E[1][event]->Fill(hitLocation[1],hitLocation[2],ener);
trueCos_pion = trueCos;
trueLen_pion = trueLen;
trueMom_pion = trueMom;
true3Mom_pion[0] = true3Mom[0];
true3Mom_pion[1] = true3Mom[1];
true3Mom_pion[2] = true3Mom[2];
}
if( PDG == 2212 ){hist2D_XY_E[2][event]->Fill(hitLocation[0],hitLocation[1],ener);
hist2D_XZ_E[2][event]->Fill(hitLocation[0],hitLocation[2],ener);
hist2D_YZ_E[2][event]->Fill(hitLocation[1],hitLocation[2],ener);
trueCos_proton = trueCos;
trueLen_proton = trueLen;
trueMom_proton = trueMom;
true3Mom_proton[0] = true3Mom[0];
true3Mom_proton[1] = true3Mom[1];
true3Mom_proton[2] = true3Mom[2];
}
}
if(newEvent && ii>0){
//KALFIT->STARTFIT(hist2D_YZ_E[0][event]);
/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////


RooFitResult* res;
rep = new KALFIT ("_rep");
char formula[10];

int nbinn = hist2D_YZ_PE[event-1]->GetNbinsY();

double statL[nbinn]={};
double meanL[nbinn]={};
double varL[nbinn]={};
double mscL[nbinn]={};

for(Int_t rbin=0;rbin<nbinn;rbin++){
statL[rbin]=0;meanL[rbin]=0;varL[rbin]=0;mscL[rbin]=0;
}
double currY = 0;
double currYl = 0;
double thetaZ = 0;
double currV = 0;

double inputMom = 1000;

std::cout<<"integrated: "<<hist2D_YZ_PE[event-1]->Integral()<<std::endl;

for(Int_t rbin=0;rbin<nbinn;rbin++)
{ 
for(Int_t rrbin=0;rrbin<hist2D_YZ_PE[event-1]->GetNbinsX();rrbin++){
currY += rrbin*hist2D_YZ_PE[event-1] -> GetBinContent(rrbin+1,rbin+1); 
currYl += hist2D_YZ_PE[event-1] -> GetBinContent(rrbin+1,rbin+1);
}
if(currYl!=0) meanL[rbin] = currY/currYl;

std::cout<<"bin numbers for Y and Z: "<<hist2D_YZ_PE[event-1]->GetNbinsX()<<" "<<nbinn<<std::endl;
std::cout<<"currY and currYl "<<currY<<" "<<currYl<<std::endl;

for(Int_t rrbin=0;rrbin<hist2D_YZ_PE[event-1]->GetNbinsX();rrbin++){
currV += hist2D_YZ_PE[event-1] -> GetBinContent(rrbin+1,rbin)* (rrbin-meanL[rbin])*(rrbin-meanL[rbin]) ;
}
if(currYl!=0) varL[rbin] = currV/currYl;

std::cout<<"currV "<<currV<<std::endl;

statL[rbin] = currYl;

thetaZ = 13.6/ inputMom * 6 * TMath::Sqrt((rbin+1)/40.) * (1 + 0.038 * TMath::Log((rbin+1)/40.));
mscL[rbin] = thetaZ * 10 * rbin;

std::cout<<"statL, thetaZ and mscL "<<currYl<<" "<<thetaZ<<" "<<mscL[rbin]<<std::endl;

currY=0;
currYl=0;
currV=0;

}

std::cout<<"loading status done.. n of bins.. "<<" "<<nbinn<<std::endl;

TH1D* MmeanS = new TH1D("","",nbinn,0,nbinn); 
TH1D* MvarianceS= new TH1D("","",nbinn,0,nbinn);; 
TH1D* MPES= new TH1D("","",nbinn,0,nbinn);; 
TH1D* MMCSS= new TH1D("","",nbinn,0,nbinn);;
for(Int_t rbin=0;rbin<nbinn;rbin++){
MmeanS -> SetBinContent(rbin+1,meanL[rbin]);
MvarianceS -> SetBinContent(rbin+1,varL[rbin]);
MPES -> SetBinContent(rbin+1,statL[rbin]);
MMCSS -> SetBinContent(rbin+1,mscL[rbin]);
}

std::cout<<"setted up the input histograms.."<<std::endl;

rep->SetHistograms(MmeanS,MvarianceS,MPES,MMCSS);
rep->FillEv(rep->getPullList());

RooArgList list("list");
list.add(*rep);
sprintf(formula,"%s","@0");
RooFormulaVar* fcn = new RooFormulaVar("fit","fit",formula,list);

std::cout<<"setting up the fitter step 1.."<<std::endl;

RooMinuit m(*fcn);
m.setStrategy(2);
Double_t callsEDM[2] = {10500., 1.e-6};
Int_t irf = 0;

std::cout<<"input values for pull paras.. "<<rep->getPar(0)<<" "<<rep->getPar(1)<<std::endl;
std::cout<<"setting up the fitter step 2.."<<std::endl;

gMinuit->mnexcm("MIGRAD",callsEDM,2,irf);
std::cout<<"running ?? "<<std::endl;
m.migrad();
//m.hesse();
//m.minos(); 
res = m.save();
double bestFit = res->minNll();

std::cout<<"fitted.. chi2.. "<<bestFit<<" charge value.. "<<rep->getPar(1)<<std::endl;

for(Int_t rbin=0;rbin<nbinn;rbin++){
statL[rbin]=0;meanL[rbin]=0;varL[rbin]=0;mscL[rbin]=0;
}
currY = 0;
currYl = 0;
thetaZ = 0;
currV = 0;
//rep->FillEv(rep->getPullList());

/////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
initA = 1;
}

eventE = event;
}



}

