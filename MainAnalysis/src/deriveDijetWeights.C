//cpp dependencies
#include <fstream>
#include <iostream>
#include <string>

//ROOT dependencies
#include "TFile.h"
#include "TObjArray.h"
#include "TTree.h"

//Non-local dependencies
#include "Utility/include/checkMakeDir.h"
#include "Utility/include/doGlobalDebug.h"
#include "Utility/include/mntToXRootdFileString.h"
#include "Utility/include/returnRootFileContentsList.h"
#include "Utility/include/plotUtilities.h"

//5.02 TeV PYTHIA 6
const Int_t nPtHat = 12;
const Int_t ptHat[nPtHat+1] = {15, 30, 50, 80, 100, 120, 170, 220, 280, 370, 460, 540, 10000};
//const Double_t ptHatCrossSections[nPtHat+1] = {.5269, .03455, .004068, .0004959, .00007096, .00001223, .000003031, .0000007746, .0000001410, .00000003216, .00000001001, 0.00000};
const Double_t ptHatCrossSections[nPtHat+1] = {.5335, .03378, .003778, .0004412, .0001511, .00006147, .00001018, .000002477, .0000006160, .0000001088, .00000002527,  0.000000007865, 0.00000000};

const Int_t nPtHat8 = 13;
const Int_t ptHat8[nPtHat8+1] = {0, 15, 30, 50, 80, 100, 120, 170, 220, 280, 370, 460, 540, 10000};
const Double_t ptHatCrossSections8[nPtHat8+1] = {67.89, .5269, .03455, .004068, .0004959, .000173, .00007096, .00001223, .000003031, .0000007746, .0000001410, .00000003216, .00000001001, 0.000000000};

//2.76TeV
//const Int_t nPtHat = 9;
//const Int_t ptHat[nPtHat+1] = {15, 30, 50, 80, 120, 170, 220, 280, 370, 10000};
//const Double_t ptHatCrossSections[nPtHat+1] = {.5269, .03455, .004068, .0004959, .00007096, .00001223, .000003031, .0000007746, .0000001410, .00000003216, .00000001001, 0.00000};
//const Double_t ptHatCrossSections[nPtHat+1] = {.2034, .01075, .001025, .00009865, .00001129, .000001465, .0000002837, .00000005323, .000000005934, 00000000.000000};


int deriveDijetWeights(const std::string inConfigFileName, const bool isPYTHIA6)
{
  if(!checkFile(inConfigFileName) || inConfigFileName.find(".txt") == std::string::npos){
    std::cout << "Given inConfigFileName \'" << inConfigFileName << "\' is invalid. return 1" << std::endl;
    return 1;
  }

  std::vector<Int_t> pthats;
  std::vector<std::string> fileList;

  std::ifstream file(inConfigFileName.c_str());
  std::string tempStr;
  while(std::getline(file, tempStr)){    
    while(tempStr.find(" ") != std::string::npos){tempStr.replace(tempStr.find(" "), 1, "");}
    if(tempStr.size() == 0) continue;
    if(tempStr.find(".root") != std::string::npos) fileList.push_back(tempStr);
    else{
      if(tempStr.substr(0, std::string("PTHAT=").size()).find("PTHAT=") != std::string::npos){
	tempStr.replace(0, std::string("PTHAT=").size(), "");
	while(tempStr.find(",") != std::string::npos){
	  pthats.push_back(std::stod(tempStr.substr(0, tempStr.find(","))));
	  tempStr.replace(0, tempStr.find(",")+1, "");
	}
	if(tempStr.size() != 0) pthats.push_back(std::stod(tempStr));
      }
      else std::cout << "WARNING: Line in \'" << inConfigFileName << "\', \'" << tempStr << "\' is invalid. check input" << std::endl;
    }
  }

  file.close();

  unsigned int pos = 0;
  while(pthats.size() > pos){
    bool didSwap = false;

    for(unsigned int pI = pos+1; pI < pthats.size(); ++pI){
      if(pthats.at(pI) < pthats.at(pos)){
	int temp = pthats.at(pI);
	pthats.at(pI) = pthats.at(pos);
	pthats.at(pos) = temp;       
	didSwap = true;
      }
    }

    if(!didSwap) ++pos;
  }

  std::vector<Double_t> crossSections;
  std::vector<Int_t> pthatCounts;

  for(unsigned int pI = 0; pI < pthats.size(); pI++){
    Int_t pthatPos = -1;
    Int_t currentPthat = pthats.at(pI);

    if(isPYTHIA6){
      for(Int_t pI2 = 0; pI2 < nPtHat; ++pI2){
	if(pthats.at(pI) == ptHat[pI2]){
	  crossSections.push_back(ptHatCrossSections[pI2]);
	  pthatCounts.push_back(0);
	  pthatPos = pI2;
	  break;
	}
      }
    }
    else{
      for(Int_t pI2 = 0; pI2 < nPtHat8; ++pI2){
        if(pthats.at(pI) == ptHat8[pI2]){
          crossSections.push_back(ptHatCrossSections8[pI2]);
          pthatCounts.push_back(0);
          pthatPos = pI2;
          break;
        }
      }
    }

    if(pthatPos < 0){
      std::cout << "WARNING: Cannot find pthat \'" << currentPthat << "\' in array. please check for valid pthat. return 1" << std::endl;
      return 1;     
    }
  }

  if(isPYTHIA6){
    pthats.push_back(ptHat[nPtHat]);
    pthatCounts.push_back(0);
    crossSections.push_back(ptHatCrossSections[nPtHat]);
  }
  else{
    pthats.push_back(ptHat8[nPtHat8]);
    pthatCounts.push_back(0);
    crossSections.push_back(ptHatCrossSections8[nPtHat8]);
  }

  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl; 

  std::cout << "Processing \'" << fileList.size() << "\' files..." << std::endl;
  for(unsigned int fileIter = 0; fileIter < fileList.size(); fileIter++){
    std::string tempFileName = mntToXRootdFileString(fileList.at(fileIter));
    std::cout << " Processing file " << fileIter << "/" << fileList.size() << ": \'" << fileList.at(fileIter) << "\'" << std::endl;
    
    TFile* inFile_p = NULL;
    inFile_p = TFile::Open(tempFileName.c_str(), "READ");
    if(inFile_p == 0){
      std::cout << "given file \'" << tempFileName << "\' was not able to be read, continue" << std::endl;
      continue;
    }
    else if(inFile_p->IsZombie()){
      std::cout << "given file \'" << tempFileName << "\' was not able to be read, continue" << std::endl;
      continue;
    }
    
    std::vector<std::string> jetTreeList = returnRootFileContentsList(inFile_p, "TTree", "ak");
    int pthatPos = -1;
    for(unsigned int jtI = 0; jtI < jetTreeList.size(); ++jtI){
      TTree* jetTree_p = (TTree*)inFile_p->Get(jetTreeList.at(jtI).c_str());
      TObjArray* branchList = (TObjArray*)jetTree_p->GetListOfBranches();

      for(Int_t bI = 0; bI < branchList->GetEntries(); ++bI){
	std::string tempStr = branchList->At(bI)->GetName();
	if(tempStr.find("pthat") != std::string::npos && tempStr.size() == std::string("pthat").size()){
	  pthatPos = jtI;
	  break;
	}
      }

      if(pthatPos >= 0) break;
    }
    
    
    if(jetTreeList.size() == 0 || pthatPos < 0){
      std::cout << "given file \'" << tempFileName << "\' contains no valid jet trees, continue" << std::endl;
      inFile_p->Close();
      delete inFile_p;
      continue;
    }
    
    TTree* jetTree_p = (TTree*)inFile_p->Get(jetTreeList.at(pthatPos).c_str());
    
    if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl; 
    
    Float_t ptHat_;
    
    jetTree_p->SetBranchStatus("*", 0);
    jetTree_p->SetBranchStatus("pthat", 1);
    
    jetTree_p->SetBranchAddress("pthat", &ptHat_);
    
    const Int_t nEntries = jetTree_p->GetEntries();
    const Int_t printInterval = TMath::Max(1, nEntries/20);

    for(Int_t entry = 0; entry < nEntries; ++entry){
      if(nEntries >= 50000 && entry%printInterval == 0) std::cout << "  Entry: " << entry << "/" << nEntries << std::endl;
      jetTree_p->GetEntry(entry);
            
      for(Int_t ptHatIter = 0; ptHatIter < (Int_t)pthats.size()-1; ptHatIter++){
	if(ptHat_ < pthats.at(ptHatIter+1)){
	  pthatCounts.at(ptHatIter)++;
	  break;
	}
      }
    } 
    
    if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl; 
    
    inFile_p->Close();
    delete inFile_p;
  }

  if(doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  ULong64_t multAdjust = 1;

  std::vector<Double_t> weights;

  bool isFirst = true;
  for(unsigned int iter = 0; iter < pthats.size()-1; iter++){

    std::cout << "Pthat: " << pthats.at(iter) << "-" << pthats.at(iter+1) << std::endl;
    std::cout << " Entries: " << pthatCounts.at(iter) << std::endl;
    Double_t weight = (crossSections.at(iter) - crossSections.at(iter+1))/(Double_t)pthatCounts.at(iter);
    std::cout << " Weight: " << weight << std::endl;

    while(isFirst && weight*multAdjust < 1.){
      multAdjust *= 10;
    }
    isFirst = false;

    std::cout << " Adjusted weight: " << weight*multAdjust << std::endl;

    weights.push_back(weight*multAdjust);
  }

  std::cout << "Pthats size: " << pthats.size() << std::endl;
  std::cout << std::endl;
  std::cout << "Weights for copy-paste: ";

  for(unsigned int iter = 0; iter < pthats.size()-1; iter++){
    std::cout << weights.at(iter) << ", ";
  }
  std::cout << std::endl;

  std::cout << "Renorm weights for copy-paste: ";

  for(unsigned int iter = 0; iter < pthats.size()-1; iter++){
    std::cout << weights.at(iter)/weights.at(0) << ", ";
  }
  std::cout << std::endl;

  std::cout << "Renorm weights2 for copy-paste: ";

  for(unsigned int iter = 0; iter < pthats.size()-1; iter++){
    std::string tempStr = "";
    Double_t val = weights.at(iter)/weights.at(0);

    if(val < 1.){
      int mult = 0;
      while(val < 1.){
	++mult;
	val *= 10;
      }
      
      tempStr = prettyString(val, 6, false);
      tempStr.replace(tempStr.find("."), 1, "");
      for(int i = 0; i < mult-1; ++i){
	tempStr = "0" + tempStr;
      }
      tempStr = "0." + tempStr;
    }
    else if(val < 10) tempStr = prettyString(val, 6, false);
    else if(val < 100) tempStr = prettyString(val, 5, false);
    else if(val < 1000) tempStr = prettyString(val, 4, false);
    else if(val < 10000) tempStr = prettyString(val, 3, false);
    else if(val < 100000) tempStr = prettyString(val, 2, false);
    else if(val < 1000000) tempStr = prettyString(val, 1, false);
    else tempStr = std::to_string((Int_t)val);

    std::cout << tempStr << ", ";
  }
  std::cout << std::endl;

  return 0;
}


int main(int argc, char *argv[])
{
  if(argc != 3){
    std::cout << "Usage: deriveDijetWeights.exe <inConfigFile> <isPYTHIA6>" << std::endl;
    std::cout << "Number of args given: " << argc << std::endl;
    for(int iter = 0; iter < argc; iter++){
      std::cout << "  argv[" << iter << "]: " << argv[iter] << std::endl;
    }
    return -1;
  }

  int retVal = 0;
  retVal = deriveDijetWeights(argv[1], std::stoi(argv[2]));

  std::cout << "Job complete. Return " << retVal << "." << std::endl;
  return retVal;
}
