#include <string>
#include <iostream>
#include <algorithm>
#include "TFile.h"
#include "TH1.h"
#include "TList.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"

bool compareStrings (const std::string& l, const std::string& r) {
  return (l==r);
}

void overlayHistos (std::string infile) {

  const char* enumNames[] = {"prompt", "hf", "lf", "other"};    
  std::vector<std::string> histNames;
    
  TFile ifile(infile.c_str());
  ifile.cd("muonexercise2");
  TList* list_of_histos = gDirectory->GetListOfKeys();
  std::cout << "found " << list_of_histos->GetSize() << " objects." << std::endl;
  for (int idx = 0; idx < list_of_histos->GetSize(); idx++) {
    TString histName = list_of_histos->At(idx)->GetName();
    TObjArray* tokens = histName.Tokenize("_");
    //if (tokens->GetEntries() != 3) std::cout << "Invalid histogram name!" << std::endl;
    TString token;
    for(unsigned int jdx=0; jdx<tokens->GetEntries()-1; ++jdx) {
      TObjString* otoken = (TObjString*)tokens->At(jdx);
      if(jdx>0) token += "_";
      token += otoken->GetString();
    }
    histNames.push_back(token.Data());
  }

  std::sort(histNames.begin(), histNames.end());
  histNames.erase(std::unique(histNames.begin(), histNames.end(), compareStrings), histNames.end());

  gStyle->SetOptStat("");    
  for (size_t idx = 0; idx < histNames.size(); idx++) {        
    double xmin = 0.5;
    double ymin = 0.6;
    double xmax = 0.75;
    double ymax = 0.85;
    if (histNames.at(idx) == "npixels") {
      xmin = 0.2;
      xmax = 0.45;
    }
    if (histNames.at(idx) == "d0") {
      xmin = 0.6;
      xmax = 0.85;
    }
    TCanvas c1("c1","c1",600,600);        
    TLegend leg(xmin,ymin,xmax,ymax);
    leg.SetFillColor(0);        
    leg.SetBorderSize(0);        
    leg.SetLineWidth(0);        
    leg.SetLineColor(0);        
    leg.Clear();
    for (int pidx = 0; pidx < 4; pidx++) {
      TH1F *h1 = (TH1F*)gDirectory->Get(Form("%s_%s", histNames.at(idx).c_str(), enumNames[pidx]));
      h1->SetLineColor(pidx+1);
      h1->SetLineWidth(2);
      double max_range = 1.;
      if (histNames.at(idx) == "pt" || histNames.at(idx) == "eta") max_range = 0.3;
      if (histNames.at(idx) == "nchi2" || histNames.at(idx) == "nlayers" || histNames.at(idx) == "puiso" || histNames.at(idx) == "nhits") max_range = 0.6;
      if (histNames.at(idx) == "d0") max_range = 0.8;
      if (histNames.at(idx).find("nvtx")!=std::string::npos) max_range = 0.3;
      if (histNames.at(idx).rfind("_eff")!=std::string::npos)
	h1->Draw("histsame"); //->GetYaxis()->SetRangeUser(0,max_range);
      else 
	h1->DrawNormalized("histsame")->GetYaxis()->SetRangeUser(0,max_range);
      leg.AddEntry(Form("%s_%s", histNames.at(idx).c_str(), enumNames[pidx]), enumNames[pidx], "L");
    }
    leg.Draw();
    c1.Print(Form("%s.png",histNames.at(idx).c_str()));
  }

  return;
}
