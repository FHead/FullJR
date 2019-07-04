//cpp dependencies
#include <iostream>
#include <string>

//ROOT dependencies
#include "TFile.h"
#include "TH2D.h"

int quickHistSum(std::string inFileName, double xMin, double yMin, double xMax, double yMax)
{
  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");
  //  TH2D* hist_p = (TH2D*)inFile_p->Get("akCs4PU3PFFlowJetAnalyzer/response_SmallBins_akCs4PU3PFFlowJetAnalyzer_PbPb_Cent0to10_LightMUAndCHID_ResponseMod0p10_AbsEta0p0to2p0_RecoGenAsymm_h");
  TH2D* hist_p = (TH2D*)inFile_p->Get("akCs4PU3PFFlowJetAnalyzer/response_akCs4PU3PFFlowJetAnalyzer_PbPb_Cent0to10_LightMUAndCHID_ResponseMod0p10_AbsEta0p0to2p0_General_h");

  Double_t sum = 0.0;
  for(Int_t bIX = 0; bIX < hist_p->GetXaxis()->GetNbins(); ++bIX){
    Double_t binCenterX = hist_p->GetXaxis()->GetBinCenter(bIX+1);
    if(binCenterX < xMin) continue;
    if(binCenterX > xMax) continue;

    for(Int_t bIY = 0; bIY < hist_p->GetYaxis()->GetNbins(); ++bIY){
      Double_t binCenterY = hist_p->GetYaxis()->GetBinCenter(bIY+1);
      if(binCenterY < yMin) continue;
      if(binCenterY > yMax) continue;

      sum += hist_p->GetBinContent(bIX+1, bIY+1);
    }
  }

  std::cout << "SUM: " << sum << std::endl;

  inFile_p->Close();
  delete inFile_p;
  
  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 6){
    std::cout << "Usage: ./bin/quickHistSum.exe <inFileName> <xMin> <yMin> <xMax> <yMax>. return 1" << std::endl;
    return 1;
  }
  
  int retVal = 0;
  retVal += quickHistSum(argv[1], std::stod(argv[2]), std::stod(argv[3]), std::stod(argv[4]), std::stod(argv[5]));
  return retVal;
}
