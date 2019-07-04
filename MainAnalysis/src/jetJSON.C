#include <algorithm>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "TFile.h"
#include "TMath.h"
#include "TTree.h"

int jetJSON(const std::string inFileName, const std::string jetName, const double pTCut, const double etaCut)
{
  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");
  TTree* jetTree_p = (TTree*)inFile_p->Get(jetName.c_str());
  TTree* hiTree_p = (TTree*)inFile_p->Get("hiEvtAnalyzer/HiTree");

  UInt_t run, lumi;
  ULong64_t evt;

  hiTree_p->SetBranchStatus("*", 0);
  hiTree_p->SetBranchStatus("run", 1);
  hiTree_p->SetBranchStatus("lumi", 1);
  hiTree_p->SetBranchStatus("evt", 1);

  hiTree_p->SetBranchAddress("run", &run);
  hiTree_p->SetBranchAddress("lumi", &lumi);
  hiTree_p->SetBranchAddress("evt", &evt);

  const Int_t nMaxJets = 500;
  Int_t nref_;
  Float_t jtpt_[nMaxJets];
  Float_t jteta_[nMaxJets];

  jetTree_p->SetBranchStatus("*", 0);
  jetTree_p->SetBranchStatus("nref", 1);
  jetTree_p->SetBranchStatus("jtpt", 1);
  jetTree_p->SetBranchStatus("jteta", 1);

  jetTree_p->SetBranchAddress("nref", &nref_);
  jetTree_p->SetBranchAddress("jtpt", jtpt_);
  jetTree_p->SetBranchAddress("jteta", jteta_);

  std::map<UInt_t, std::map<UInt_t, UInt_t> > runLumiMapCount;

  const Int_t nEntries = jetTree_p->GetEntries();
  const Int_t nDiv = nEntries/20;

  std::cout << "Processing " << nEntries << std::endl;
  for(Int_t entry = 0; entry < nEntries; ++entry){
    if(entry%nDiv == 0) std::cout << " Entry " << entry << "/" << nEntries << std::endl;
    jetTree_p->GetEntry(entry);
    hiTree_p->GetEntry(entry);

    bool isEvtGood = false;
    for(Int_t jI = 0; jI < nref_; ++jI){
      if(jtpt_[jI] < pTCut) continue;
      if(TMath::Abs(jteta_[jI]) >= etaCut) continue;

      isEvtGood = true;
      break;
    }

    if(isEvtGood){
      ++((runLumiMapCount[run])[lumi]);
    }
  }

  inFile_p->Close();
  delete inFile_p;

  std::vector<UInt_t> sortedRuns;
  for(auto const & iter : runLumiMapCount){
    sortedRuns.push_back(iter.first);
  }

  std::sort(std::begin(sortedRuns), std::end(sortedRuns));

  std::string totalStr = "{";
  for(auto const & iter : sortedRuns){
    std::map<UInt_t, UInt_t> tempMap = runLumiMapCount[iter];

    std::vector<UInt_t> sortedLumis;
    for(auto const & iter2 : tempMap){
      sortedLumis.push_back(iter2.first);
    }

    std::string lumiStr = "";
    UInt_t startLumi = 1;
    UInt_t prevLumi = 0;
    UInt_t currLumi = 0;
    for(auto const & iter2 : sortedLumis){
      currLumi = iter2;
      if(currLumi - prevLumi > 1){
	if(prevLumi != 0){
	  if(prevLumi - startLumi != 0) lumiStr = lumiStr + std::to_string(startLumi) + "-" + std::to_string(prevLumi) + ", ";
	  else lumiStr = lumiStr + std::to_string(currLumi) + ", ";
	}
        
	startLumi = currLumi;
      }

      prevLumi = currLumi;
    }
    
    if(currLumi - startLumi == 0) lumiStr = lumiStr + std::to_string(currLumi);
    else lumiStr = lumiStr + std::to_string(startLumi) + "-" + std::to_string(currLumi);
    
    //    std::cout << lumiTemp << std::endl;

    //    std::cout << "\"" << iter << "\":";
    while(lumiStr.find(", ") != std::string::npos){
      lumiStr.replace(lumiStr.find(", "), 2, "],[");
    }
    while(lumiStr.find(",") != std::string::npos){
      lumiStr.replace(lumiStr.find(","), 1, "-");
    }
    while(lumiStr.find("-") != std::string::npos){
      lumiStr.replace(lumiStr.find("-"), 1, ", ");
    }
    lumiStr = "[[" + lumiStr + "]]";

    //    std::cout << " " << lumiStr << std::endl;

    totalStr = totalStr + "\"" + std::to_string(iter) + "\":" + lumiStr + ", ";
  }
  totalStr.replace(totalStr.size()-2, 2, "}");
  std::cout << totalStr << std::endl;

  return 0;
}


int main(int argc, char* argv[])
{
  if(argc != 5){
    std::cout << "Usage: ./bin/jetJSON.exe <inFileName> <jetName> <pTCut> <etaCut>" << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += jetJSON(argv[1], argv[2], std::stod(argv[3]), std::stod(argv[4]));
  return retVal;
}
