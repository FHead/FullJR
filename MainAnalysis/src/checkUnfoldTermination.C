//cpp dependencies
#include <iostream>
#include <string>
#include <vector>

//ROOT dependencies
#include "TFile.h"

//Non-local (Utility) dependencies
#include "Utility/include/checkMakeDir.h"

//Local (MainAnalysis) dependencies
#include "MainAnalysis/include/cutPropagator.h"

int checkUnfoldTermination(const std::string inFileName)
{
  if(!checkFile(inFileName) || inFileName.find(".root") == std::string::npos){
    std::cout << "Given inFileName \'" << inFileName << "\' is invalid. return 1" << std::endl;
    return 1;
  }

  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");
  cutPropagator cutProp;
  cutProp.Clean();
  cutProp.GetAllVarFromFile(inFile_p);
 
  std::vector<std::string> histTag = cutProp.GetHistTagBayes();
  std::vector<int> histBestBayes = cutProp.GetHistBestBayes();

  int totalNegPos = 0;
  int totalNegAndHighPos = 0;

  int nPassesSel = 0;
  int nPassesSelNegPos = 0;
  int nPassesSelNegAndHighPos = 0;

  for(unsigned int i = 0; i < histTag.size(); ++i){
    bool passesSel = true;
    if(histTag.at(i).find("AbsEta0p0to2p0") == std::string::npos) passesSel = false;
    else if(histTag.at(i).find("LightMUAndCHID") == std::string::npos) passesSel = false;
    else if(histTag.at(i).find("ResponseMod0p10") == std::string::npos) passesSel = false;

    if(passesSel) nPassesSel++;

    if(histBestBayes.at(i) < 0){
      std::cout << "WARNING: " << histTag.at(i) << " has terminalPos of " << histBestBayes.at(i) << std::endl;
      totalNegPos += 1;
      totalNegAndHighPos += 1;
      if(passesSel){
	nPassesSelNegPos++;
	nPassesSelNegAndHighPos++;
      }
    }
    else if(histBestBayes.at(i) > 20){
      std::cout << "WARNING: " << histTag.at(i) << " has terminalPos of " << histBestBayes.at(i) << std::endl;
      totalNegAndHighPos += 1;
      if(passesSel){
	nPassesSelNegAndHighPos++;
      }
    }
  }
  
  std::cout << "Total hist with zero termination: " << totalNegPos << "/" << histTag.size() << std::endl;
  std::cout << "Total hist with zero or high termination: " << totalNegAndHighPos << "/" << histTag.size() << std::endl;


  std::cout << "Passing hist with zero termination: " << nPassesSelNegPos << "/" << nPassesSel << std::endl;
  std::cout << "Passing hist with zero or high termination: " << nPassesSelNegAndHighPos << "/" << nPassesSel << std::endl;


  inFile_p->Close();
  delete inFile_p;

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./bin/checkUnfoldTermination.exe <inFileName>" << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += checkUnfoldTermination(argv[1]);
  return retVal;
}
