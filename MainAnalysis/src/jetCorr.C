#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "TDatime.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TTree.h"

#include "Utility/include/checkMakeDir.h"
#include "Utility/include/etaPhiFunc.h"
#include "Utility/include/ncollFunctions_5TeV.h"
#include "Utility/include/plotUtilities.h"
#include "Utility/include/stringUtil.h"

int jetCorr(std::string inFileName, std::string pthatFile)
{
  std::vector<double> pthats;
  std::vector<double> pthatWeights;

  std::ifstream inFile(pthatFile.c_str());
  std::string tempStr;
  while(std::getline(inFile, tempStr)){
    std::cout << tempStr << std::endl;
    std::vector<std::string> subStrs;
    tempStr = tempStr + ",";
    while(tempStr.find(",,") != std::string::npos){
      tempStr.replace(tempStr.find(",,"), 2, ",");
    }
    while(tempStr.find("=") != std::string::npos){
      tempStr.replace(tempStr.find("="), 1, ",");
    }
    if(tempStr.size() == 1) continue;
    while(tempStr.find(",") != std::string::npos){
      subStrs.push_back(tempStr.substr(0, tempStr.find(",")));
      tempStr.replace(0, tempStr.find(",")+1, "");
    }
    if(subStrs.size() == 0) continue;
    if(subStrs[0].find("PTHAT") == std::string::npos) continue;

    if(isStrSame("PTHAT", subStrs[0])){
      for(unsigned int sI = 1; sI < subStrs.size(); ++sI){
	pthats.push_back(std::stod(subStrs[sI]));
      }
    }
    else if(isStrSame("PTHATWEIGHTS", subStrs[0])){
      for(unsigned int sI = 1; sI < subStrs.size(); ++sI){
	pthatWeights.push_back(std::stod(subStrs[sI]));
      }
    }
  }
  inFile.close();
  pthats.push_back(100000000);

  TDatime* date = new TDatime();
  const std::string dateStr = std::to_string(date->GetDate());
  delete date;

  const Int_t nCentBins = 4;
  Int_t centBinsLow[nCentBins] = {0, 10, 30, 50};
  Int_t centBinsHi[nCentBins] = {10, 30, 50, 90};

  const Int_t nJets = 6;
  const Int_t jetR[nJets] = {2, 3, 4, 6, 8, 10};
  Int_t r4Pos = -1;
  for(Int_t jI = 0; jI < nJets; ++jI){
    if(jetR[jI] == 4){
      r4Pos = jI;
    }
  }

  checkMakeDir("output");
  checkMakeDir("output/" + dateStr);
  std::string outFileName = inFileName;
  while(outFileName.find("/") != std::string::npos){
    outFileName.replace(0, outFileName.find("/")+1, "");
  }
  outFileName.replace(outFileName.find(".root"), 5, "");

  outFileName = "output/" + dateStr + "/" + outFileName + "_JETCORR_" + dateStr + ".root";

  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");
  TH1D* pthat_p[nCentBins];
  TH1D* pthatWeight_p[nCentBins];
  TH2D* histPtCorr_p[nJets][nCentBins];
  TH2D* histPtCorrWeight_p[nJets][nCentBins];

  const Double_t min = 20;
  const Double_t max = 600;
  Int_t nBins = (max - min)/10;

  for(Int_t cI = 0; cI < nCentBins; ++cI){
    std::string centStr = "Cent" + std::to_string(centBinsLow[cI]) + "to" + std::to_string(centBinsHi[cI]);

    pthat_p[cI] = new TH1D(("pthat_" + centStr + "_h").c_str(), ";p_{T}Hat;Counts", 600, 15, 615);
    pthatWeight_p[cI] = new TH1D(("pthatWeight_" + centStr + "_h").c_str(), ";p_{T}Hat;Counts", 600, 15, 615);

    for(Int_t jI = 0; jI < nJets; ++jI){
      if(jI == r4Pos) continue;
      
      histPtCorr_p[jI][cI] = new TH2D(("histPtCorr_" + centStr + "_R" + std::to_string(jetR[jI]) + "_h").c_str(), (";Max Jet pT R=" + prettyString(jetR[jI]/10., 1, false) + ";Matched Jet pT R=" + prettyString(jetR[r4Pos]/10., 1, false)).c_str(), nBins, min, max, nBins, min, max);
      
      histPtCorrWeight_p[jI][cI] = new TH2D(("histPtCorrWeight_" + centStr + "_R" + std::to_string(jetR[jI]) + "_h").c_str(), (";Max Jet pT R=" + prettyString(jetR[jI]/10., 1, false) + ";Matched Jet pT R=" + prettyString(jetR[r4Pos]/10., 1, false)).c_str(), nBins, min, max, nBins, min, max);
    }
  }

  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");
  TTree* jetTrees_p[nJets];
  TTree* hiTree_p = (TTree*)inFile_p->Get("hiEvtAnalyzer/HiTree");


  for(Int_t jI = 0; jI < nJets; ++jI){
    jetTrees_p[jI] = (TTree*)inFile_p->Get(("akCs" + std::to_string(jetR[jI]) + "PU3PFFlowJetAnalyzer/t").c_str());
  }

  Int_t hiBin;
  Float_t pthat;

  const Int_t nMaxJets = 500;
  Int_t nJt_[nJets];
  Float_t jtpt_[nJets][nMaxJets];
  Float_t jtphi_[nJets][nMaxJets];
  Float_t jteta_[nJets][nMaxJets];

  hiTree_p->SetBranchStatus("*", 0);
  hiTree_p->SetBranchStatus("hiBin", 1);
  hiTree_p->SetBranchStatus("pthat", 1);

  hiTree_p->SetBranchAddress("hiBin", &hiBin);
  hiTree_p->SetBranchAddress("pthat", &pthat);

  for(Int_t jI = 0; jI < nJets; ++jI){
    jetTrees_p[jI]->SetBranchStatus("*", 0);
    jetTrees_p[jI]->SetBranchStatus("nref", 1);
    jetTrees_p[jI]->SetBranchStatus("jtpt", 1);
    jetTrees_p[jI]->SetBranchStatus("jtphi", 1);
    jetTrees_p[jI]->SetBranchStatus("jteta", 1);

    jetTrees_p[jI]->SetBranchAddress("nref", &(nJt_[jI]));
    jetTrees_p[jI]->SetBranchAddress("jtpt", jtpt_[jI]);
    jetTrees_p[jI]->SetBranchAddress("jtphi", jtphi_[jI]);
    jetTrees_p[jI]->SetBranchAddress("jteta", jteta_[jI]);
  }

  Double_t maxPt[nJets];
  Double_t maxPhi[nJets];
  Double_t maxEta[nJets];

  const Int_t nEntries = jetTrees_p[0]->GetEntries();
  for(Int_t entry = 0; entry < nEntries; ++entry){
    hiTree_p->GetEntry(entry);

    //    if(hiBin > 20) continue;

    Int_t centBinPos = -1;
    for(Int_t cI = 0; cI < nCentBins; ++cI){
      if(hiBin/2 >= centBinsLow[cI] && hiBin/2 < centBinsHi[cI]){
	centBinPos = cI;
      }
    }

    if(centBinPos < 0) continue;

    int weightPos = -1;
    for(unsigned int pI = 0; pI < pthats.size(); ++pI){
      if(pthat >= pthats[pI] && pthat < pthats[pI+1]){
	weightPos = pI;
	break;
      }
    }
    double weight = pthatWeights[weightPos];
    weight *= findNcoll_Renorm(hiBin);

    pthat_p[centBinPos]->Fill(pthat);
    pthatWeight_p[centBinPos]->Fill(pthat, weight);

    for(Int_t jI = 0; jI < nJets;++jI){
      maxPt[jI] = -999;
      maxPhi[jI] = -999;
      maxEta[jI] = -999;

      jetTrees_p[jI]->GetEntry(entry);

      for(Int_t pI = 0; pI < nJt_[jI]; ++pI){
	if(TMath::Abs(jteta_[jI][pI]) > 2.) continue;
	if(jtpt_[jI][pI] < min) continue;
	if(jtpt_[jI][pI] < maxPt[jI]) continue;

	maxPt[jI] = jtpt_[jI][pI];
	maxPhi[jI] = jtphi_[jI][pI];
	maxEta[jI] = jteta_[jI][pI];
      }      
    }

    for(Int_t jI = 0; jI < nJets;++jI){
      if(r4Pos == jI) continue;
      if(maxPt[jI] < 0) continue;

      Double_t maxMatchPt = 0;
      //      Double_t maxMatchPhi = -999;
      //      Double_t maxMatchEta = -999;

      for(Int_t pI = 0; pI < nJt_[r4Pos]; ++pI){
        if(TMath::Abs(jteta_[r4Pos][pI]) > 2.) continue;
	if(getDR(jteta_[r4Pos][pI], jtphi_[r4Pos][pI], maxEta[jI], maxPhi[jI]) > 0.4) continue;
	if(jtpt_[r4Pos][pI] < maxMatchPt) continue;
	
	maxMatchPt = jtpt_[r4Pos][pI];
	//	maxMatchPhi = jtphi_[r4Pos][pI];
	//	maxMatchEta = jteta_[r4Pos][pI];
      }
    
      if(maxMatchPt >= min){
	histPtCorr_p[jI][centBinPos]->Fill(maxPt[jI], maxMatchPt);
	histPtCorrWeight_p[jI][centBinPos]->Fill(maxPt[jI], maxMatchPt, weight);
      }
      else{
	histPtCorr_p[jI][centBinPos]->Fill(maxPt[jI], min+0.5);
	histPtCorrWeight_p[jI][centBinPos]->Fill(maxPt[jI], min+0.5, weight);
      }
    }
  }

  inFile_p->Close();
  delete inFile_p;

  outFile_p->cd();

  for(Int_t cI = 0; cI < nCentBins; ++cI){
    pthat_p[cI]->Write("", TObject::kOverwrite);
    delete pthat_p[cI];
    
    pthatWeight_p[cI]->Write("", TObject::kOverwrite);
    delete pthatWeight_p[cI];
    
    for(Int_t jI = 0; jI < nJets; ++jI){
      if(jI == r4Pos) continue;
      
      histPtCorr_p[jI][cI]->Write("", TObject::kOverwrite);
      delete histPtCorr_p[jI][cI];
      
      histPtCorrWeight_p[jI][cI]->Write("", TObject::kOverwrite);
      delete histPtCorrWeight_p[jI][cI];
    }
  }

  outFile_p->Close();
  delete outFile_p;

  return 0;
}


int main(int argc, char* argv[])
{
  if(argc != 3){
    std::cout << "Usage: ./bin/jetCorr.exe <inFileName> <pthatFile>. return 1" << std::endl;
    return 1;
  }

  int retVal = 0;  
  retVal += jetCorr(argv[1], argv[2]);
  return retVal;
}
