#include <iostream>
#include <string>
#include <vector>

#include "TFile.h"
#include "TMath.h"
#include "TTree.h"

#include "Utility/include/returnRootFileContentsList.h"

int counts(const std::string inFileName)
{
  
  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");
  const Int_t nMaxJets = 500;
  Int_t nref_;
  Float_t jtpt_[nMaxJets];
  Float_t jteta_[nMaxJets];

  Int_t hiBin_;
  Float_t hiHF_;
  Float_t vz_;

  Int_t HBHENoiseFilterResultRun2Loose_;
  Int_t pprimaryVertexFilter_;
  Int_t phfCoincFilter3_;
  Int_t pclusterCompatibilityFilter_;
  Int_t pBeamScrapingFilter_;

  std::vector<std::string> treeNames = returnRootFileContentsList(inFile_p, "TTree", "JetAna");

  unsigned int pos = 0;
  while(pos < treeNames.size()){
    if(treeNames[pos].find("ak1P") != std::string::npos || treeNames[pos].find("ak5P") != std::string::npos){
      treeNames.erase(treeNames.begin()+pos);
    }
    else ++pos;
  }

  const Int_t nCentBins = 4;
  Int_t centBinsLow[nCentBins] = {0, 10, 30, 50};
  Int_t centBinsHigh[nCentBins] = {10, 30, 50, 90};

  for(unsigned int rI = 0; rI < treeNames.size(); ++rI){
    TTree* tree_p = (TTree*)inFile_p->Get(treeNames[rI].c_str());
    TTree* hiTree_p = (TTree*)inFile_p->Get("hiEvtAnalyzer/HiTree");
    TTree* skimTree_p = (TTree*)inFile_p->Get("skimanalysis/HltTree");

    Double_t secondCut[nCentBins];
    if(treeNames[rI].find("ak2PF") != std::string::npos){
      secondCut[0] = 160;
      secondCut[1] = 160;
      secondCut[2] = 160;
      secondCut[3] = 160;
    }
    else if(treeNames[rI].find("ak3PF") != std::string::npos){
      secondCut[0] = 160;
      secondCut[1] = 160;
      secondCut[2] = 160;
      secondCut[3] = 160;
    }
    else if(treeNames[rI].find("ak4PF") != std::string::npos){
      secondCut[0] = 160;
      secondCut[1] = 160;
      secondCut[2] = 160;
      secondCut[3] = 160;
    }
    else if(treeNames[rI].find("ak6PF") != std::string::npos){
      secondCut[0] = 160;
      secondCut[1] = 160;
      secondCut[2] = 160;
      secondCut[3] = 160;
    }
    else if(treeNames[rI].find("ak8PF") != std::string::npos){
      secondCut[0] = 160;
      secondCut[1] = 160;
      secondCut[2] = 160;
      secondCut[3] = 160;
    }
    else if(treeNames[rI].find("ak10PF") != std::string::npos){
      secondCut[0] = 160;
      secondCut[1] = 160;
      secondCut[2] = 160;
      secondCut[3] = 160;
    }
    else if(treeNames[rI].find("akCs2") != std::string::npos){
      secondCut[0] = 250;
      secondCut[1] = 160;
      secondCut[2] = 160;
      secondCut[3] = 160;
    }
    else if(treeNames[rI].find("akCs3") != std::string::npos){
      secondCut[0] = 300;
      secondCut[1] = 160;
      secondCut[2] = 160;
      secondCut[3] = 160;
    }
    else if(treeNames[rI].find("akCs4") != std::string::npos){
      secondCut[0] = 300;
      secondCut[1] = 160;
      secondCut[2] = 160;
      secondCut[3] = 160;
    }
    else if(treeNames[rI].find("akCs6") != std::string::npos){
      secondCut[0] = 300;
      secondCut[1] = 200;
      secondCut[2] = 200;
      secondCut[3] = 200;
    }
    else if(treeNames[rI].find("akCs8") != std::string::npos){
      secondCut[0] = 350;
      secondCut[1] = 350;
      secondCut[2] = 200;
      secondCut[3] = 200;
    }
    else if(treeNames[rI].find("akCs10") != std::string::npos){
      secondCut[0] = 400;
      secondCut[1] = 350;
      secondCut[2] = 250;
      secondCut[3] = 200;
    }

    tree_p->SetBranchStatus("*", 0);
    tree_p->SetBranchStatus("nref", 1);
    tree_p->SetBranchStatus("jtpt", 1);
    tree_p->SetBranchStatus("jteta", 1);

    tree_p->SetBranchAddress("nref", &nref_);
    tree_p->SetBranchAddress("jtpt", jtpt_);
    tree_p->SetBranchAddress("jteta", jteta_);

    hiTree_p->SetBranchStatus("*", 0);
    hiTree_p->SetBranchStatus("hiBin", 1);
    hiTree_p->SetBranchStatus("vz", 1);
    hiTree_p->SetBranchStatus("hiHF", 1);

    hiTree_p->SetBranchAddress("hiBin", &hiBin_);
    hiTree_p->SetBranchAddress("vz", &vz_);
    hiTree_p->SetBranchAddress("hiHF", &hiHF_);

    skimTree_p->SetBranchStatus("*", 0);
    skimTree_p->SetBranchStatus("HBHENoiseFilterResultRun2Loose", 1);
    if(treeNames[rI].find("Cs") == std::string::npos){
      skimTree_p->SetBranchStatus("pBeamScrapingFilter", 1);
      skimTree_p->SetBranchStatus("pPAprimaryVertexFilter", 1);
    }
    else{
      skimTree_p->SetBranchStatus("pprimaryVertexFilter", 1);
      skimTree_p->SetBranchStatus("phfCoincFilter3", 1);
      skimTree_p->SetBranchStatus("pclusterCompatibilityFilter", 1);
    }

    skimTree_p->SetBranchAddress("HBHENoiseFilterResultRun2Loose", &HBHENoiseFilterResultRun2Loose_);
    if(treeNames[rI].find("Cs") == std::string::npos){
      skimTree_p->SetBranchAddress("pBeamScrapingFilter", &pBeamScrapingFilter_);
      skimTree_p->SetBranchAddress("pPAprimaryVertexFilter", &pprimaryVertexFilter_);     
    }
    else{
      skimTree_p->SetBranchAddress("pprimaryVertexFilter", &pprimaryVertexFilter_);      
      skimTree_p->SetBranchAddress("phfCoincFilter3", &phfCoincFilter3_);
      skimTree_p->SetBranchAddress("pclusterCompatibilityFilter", &pclusterCompatibilityFilter_);
    }


    const Int_t nEntries = tree_p->GetEntries();

    Int_t counts = 0;
    Int_t counts2[nCentBins] = {0, 0, 0, 0};   

    for(Int_t entry = 0; entry < nEntries; ++entry){
      tree_p->GetEntry(entry);
      hiTree_p->GetEntry(entry);
      skimTree_p->GetEntry(entry);

      if(treeNames[rI].find("Cs") != std::string::npos){
	if(!phfCoincFilter3_) continue;
	if(!pclusterCompatibilityFilter_) continue;

	if(hiBin_ >= 180) continue;
	if(hiHF_ >= 5500.) continue;
      }
      else{
	if(!pBeamScrapingFilter_) continue;
      }

      if(!HBHENoiseFilterResultRun2Loose_) continue;
      if(!pprimaryVertexFilter_) continue;
      if(TMath::Abs(vz_) >= 15.) continue;

      Int_t centPos = -1;
      if(treeNames[rI].find("Cs") != std::string::npos){
	for(Int_t cI = 0; cI < nCentBins; ++cI){
	  if(hiBin_ >= centBinsLow[cI]*2 && hiBin_ < centBinsHigh[cI]*2){
	    centPos = cI;
	    break;
	  }
	}
      }
      else centPos = 0;

      bool isGood = false;
      bool isGood2 = false;
      for(Int_t jI = 0; jI < nref_; ++jI){
	if(TMath::Abs(jteta_[jI]) >= 2.) continue;
	if(jtpt_[jI] >= 140.) isGood = true;
	if(jtpt_[jI] >= secondCut[centPos]) isGood2 = true;

	if(isGood && isGood2) break;
      }

      if(isGood) ++counts;
      if(isGood2) ++counts2[centPos];
    }

    //    std::cout << treeNames[rI] << ": " << counts << std::endl;
    if(treeNames[rI].find("Cs") != std::string::npos){
      std::string rParam = treeNames[rI].substr(0, treeNames[rI].find("PU"));
      rParam.replace(0, 4, "");
      if(rParam.size() == 2) rParam = rParam.substr(0, 1) + "." + rParam.substr(1, 1);
      else rParam = "0." + rParam;

      std::cout << "PbPb & High $p_{T}$ Jet100 & anti-$k_{T}$ R=" << rParam << " jet ";

      for(Int_t cI = 0; cI < nCentBins; ++cI){
	//	std::cout << " pT > " << secondCut[cI] << ", " << centBinsLow[cI] << "-" << centBinsHigh[cI] << "%: " << counts2[cI] << std::endl;
	std::cout << "& " << counts2[cI] << " ($p_{T}>" << secondCut[cI] << "$ GeV/c) ";
      }
      std::cout << " \\\\ \\hline" << std::endl;
    }
    else{
      std::string rParam = treeNames[rI].substr(0, treeNames[rI].find("PF"));
      rParam.replace(0, 2, "");
      if(rParam.size() == 2) rParam = rParam.substr(0, 1) + "." + rParam.substr(1, 1);
      else rParam = "0." + rParam;

      std::cout << "pp & High $p_{T}$ Jet80 & anti-$k_{T}$ R=" << rParam << " jet $p_{T}>" << secondCut[0] << "$ GeV/c & " << counts2[0] << "\\\\ \\hline" << std::endl;

      //      std::cout << " pT > " << secondCut[0] << ": " << counts2[0] << std::endl;
    }
  }


  inFile_p->Close();
  delete inFile_p;

  return 0;
}


int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./bin/counts.exe <inFileName>. return 1" << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += counts(argv[1]);
  return retVal;
}
