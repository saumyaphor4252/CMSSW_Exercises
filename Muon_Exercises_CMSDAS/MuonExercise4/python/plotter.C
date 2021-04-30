#include <string.h> 
#include "TChain.h"
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "TKey.h"
#include "Riostream.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TF1.h"
#include "TFile.h"
#include "TFitResult.h"
#include "TGraph.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraphErrors.h"
#include "THStack.h"
#include "TLegend.h"
#include "TMinuit.h"
#include "TMath.h"
#include "TProfile.h"
#include "TPaveStats.h"
#include "TPaveText.h"
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"
#include "TGraphAsymmErrors.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include "TString.h"
#include "TFrame.h"
#include "TSystem.h"
#include "TInterpreter.h"

using namespace std;

void plotter(TString file, TString histo1d, TString histo1n, TString XTitle, TString YTitle, double setxmin, double setxmax, int rebin) {

gROOT->ProcessLine(".L tdrstyle.C");
setTDRStyle();

TFile *fileOpened = TFile::Open(file+".root");

TH1F *histoDen = (TH1F*)fileOpened->Get("demo/"+histo1d);
TH1F *histoNum = (TH1F*)fileOpened->Get("demo/"+histo1n);

if(setxmin == 20){
  
  double binn[12]={20,30,40,50,60,70,80,90,100,130,160,180};
  double t = 11;
  
  histoDen->Rebin(t,"NData_den1",binn);
  histoNum->Rebin(t,"NData_num1",binn);
  
  NData_den1->Draw();
  NData_num1->Draw();
  
  histoDen = NData_den1;
  histoNum = NData_num1;
 } else {
 
// rebinning histos
//
histoDen->Rebin(rebin);
histoNum->Rebin(rebin);

}
TString canvName = histo1n;

TCanvas* canv =  new TCanvas(canvName,canvName,0,0,500,500); 

TGraphAsymmErrors *dividend1 = new TGraphAsymmErrors(histoNum,histoDen,"cp");
  dividend1->SetLineColor(kBlack);
  dividend1->SetMarkerColor(kRed);//kPink+10);
  dividend1->SetMarkerStyle(8);
  dividend1->SetMarkerStyle(20);
  dividend1->SetTitle("");
  dividend1->GetYaxis()->SetTitle(YTitle);
  dividend1->GetXaxis()->SetTitle(XTitle);
  dividend1->GetXaxis()->SetRangeUser(setxmin,setxmax);
  dividend1->SetMaximum(1.2); 
  dividend1->SetMinimum(0.);
  dividend1->Draw("AP");
  TLegend *leg_fonchia = new TLegend(0.189516,0.838983,0.689516,0.938559);
  leg_fonchia->SetFillColor(0);
  leg_fonchia->SetLineColor(1);
  leg_fonchia->SetTextSize(0.04);

  leg_fonchia->AddEntry(dividend1, "Z#rightarrow#mu#mu", "p");
  leg_fonchia->Draw("same");
 
  canv->Print(canvName+".pdf",".pdf");
  canv->Print(canvName+".png",".png");
  canv->SaveAs(canvName+".root",".root");






}
