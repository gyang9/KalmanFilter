#include <TFile.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TTree.h>
#include <TVector3.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TMath.h>
#include <TStopwatch.h>
#include <TPolyMarker3D.h>
#include <string>
#include <iostream>
#include <fstream>

#include "kalFit.hh"

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

TFile file(Form("/home/gyang/work/dune-ndx/spack/var/spack/stage/edep-sim-master-rpl36hsf6e7fnmtgw4pd74gyffybpekl/edep-sim/SimChain/canDeleteData/testEvent_particleGun2000MeVMuon_2018_shift_sample24.root"));
//TFile file(Form("testEvent_%s.root",argv[1]));
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
hist2D_XY_Q[i] = new TH2F("","",200,0,2000,200,0,2000);
hist2D_XZ_Q[i] = new TH2F("","",200,0,2000,240,0,2400);
hist2D_YZ_Q[i] = new TH2F("","",200,0,2000,240,0,2400);

hist2D_XY_PE[i] = new TH2F("","",200,0,2000,200,0,2000);
hist2D_XZ_PE[i] = new TH2F("","",200,0,2000,240,0,2400);
hist2D_YZ_PE[i] = new TH2F("","",200,0,2000,240,0,2400);

hist2D_XY_ADC[i] = new TH2F("","",200,0,2000,200,0,2000);
hist2D_XZ_ADC[i] = new TH2F("","",200,0,2000,240,0,2400);
hist2D_YZ_ADC[i] = new TH2F("","",200,0,2000,240,0,2400);
}

TH2D* hist2D_XY_E[3][nnevent];
TH2D* hist2D_XZ_E[3][nnevent];
TH2D* hist2D_YZ_E[3][nnevent];
for(Int_t i=0;i<3;i++){
for(Int_t j=0;j<nnevent;j++){
hist2D_XY_E[i][j] = new TH2D("","",200,0,2000,200,0,2000);
hist2D_YZ_E[i][j] = new TH2D("","",200,0,2000,240,0,2400);
hist2D_XZ_E[i][j] = new TH2D("","",200,0,2000,240,0,2400);
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
KALFIT* KALFIT;

for(Int_t ii=0;ii<nevent;ii++){

c->GetEntry(ii);

//if(event == 226 || event ==227) continue;
event -= atoi(argv[1])*200;

eventS = event;
//list[eventS] =1;
if(eventS != eventE){newEvent = true;}
else {newEvent = false;}
if(newEvent) {
std::cout<<"event "<<event<<std::endl;
std::cout<<hitLocation[0]<<" "<<hitLocation[1]<<" "<<hitLocation[2]<<std::endl;
}
//if(list[eventS]==0 && newEvent) cout<<"double new "<<endl;
if( !newEvent ){
//if( event ==1){
hist2D_XY_Q[event]->Fill(hitLocation[0],hitLocation[1],Q[0]+Q[1]);
hist2D_XZ_Q[event]->Fill(hitLocation[0],hitLocation[2],Q[0]+Q[2]);
hist2D_YZ_Q[event]->Fill(hitLocation[1],hitLocation[2],Q[1]+Q[2]);
hist2D_XY_PE[event]->Fill(hitLocation[0],hitLocation[1],hitPE[0]+hitPE[1]);
hist2D_XZ_PE[event]->Fill(hitLocation[0],hitLocation[2],hitPE[0]+hitPE[2]);
hist2D_YZ_PE[event]->Fill(hitLocation[1],hitLocation[2],hitPE[1]+hitPE[2]);
hist2D_XY_ADC[event]->Fill(hitLocation[0],hitLocation[1],adc[0]+adc[1]);
hist2D_XZ_ADC[event]->Fill(hitLocation[0],hitLocation[2],adc[0]+adc[2]);
hist2D_YZ_ADC[event]->Fill(hitLocation[1],hitLocation[2],adc[1]+adc[2]);
//cout<<"PDG and energy "<<PDG<<" "<<ener<<endl;

if((PDG>100 && PDG<10000) || (PDG<-100 && PDG>-10000)){
if(PDG>0){if(happened[PDG]==1){doLoop==false; } else doLoop == true;}
else{if(happened[-PDG+1]==1){doLoop==false; } else doLoop == true;}
if(lastHenergy != PDG && doLoop) {hardEnergy += trueMom;}
//cout<<PDG<<" "<<trueMom<<" "<<ener<<endl;
lastHenergy = PDG;
if(PDG>0)happened[PDG]=1;
else happened[-PDG+1]=1;
}

if((PDG>100 && PDG<10000) || (PDG<-100 && PDG>-10000)){
hardColl += ener;
}

//proton2212  pion 211  muon 13  electron 11
if(PDG == 13 || PDG == -13 ){

if(initA==0 && TMath::Sqrt(TMath::Power(hitLocation[0]-hitInit[0],2) + TMath::Power(hitLocation[1]-hitInit[1],2) + TMath::Power(hitLocation[2]-hitInit[2],2))<atof(argv[2]) ){

hist2D_XY_E[0][event]->Fill(hitLocation[0],hitLocation[1],ener);
hist2D_XZ_E[0][event]->Fill(hitLocation[0],hitLocation[2],ener);
hist2D_YZ_E[0][event]->Fill(hitLocation[1],hitLocation[2],ener);

// PE /50 works fine, /200 works bad for getting rid of delta ray.
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
//list[ii]=1;
//eventE = event;
//}

if(newEvent && ii>0){

//KALFIT* KALFIT;
KALFIT->STARTFIT(hist2D_YZ_E[0][event]);
 
initA = 1;

} 

eventE = event;
}
}
