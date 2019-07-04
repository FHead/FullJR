//cpp depdendencies
#include <iostream>
#include <map>
#include <string>
#include <unordered_map>

//ROOT dependencies
#include "TDirectoryFile.h"
#include "TFile.h"
#include "TH1D.h"
#include "TKey.h"
#include "TMath.h"

//Local dependencies

//Non-local dependencies
#include "Utility/include/checkMakeDir.h"
#include "Utility/include/cppWatch.h"

int compareAllHist(const std::string inFileName1, const std::string inFileName2, const std::string stringToRemoveFile2 = "")
{
  cppWatch totalTimer;
  totalTimer.start();
  
  const Int_t nFiles = 2;
  TFile* inFile_p[nFiles];
  std::vector<TDirectoryFile*> dirs_p;
  const std::string fileNames[nFiles] = {inFileName1, inFileName2};
  std::map<std::string, ULong64_t> histNamesToCounts;
  std::map<std::string, TH1D*> histNamesToFile1Contents;
  std::unordered_map<std::string, TH1D*> histNamesToFile2Contents;
  ULong64_t totalKeys = 0;
  
  for(Int_t fI = 0; fI < nFiles; ++fI){
    if(checkFile(fileNames[fI]) && fileNames[fI].find(".root") != std::string::npos) continue;
    std::cout << "File \'" << fileNames[fI] << "\' is not valid. return 1" << std::endl;
    return 1;    
  }  
  
  for(Int_t tI = 0; tI < nFiles; ++tI){
    inFile_p[tI] = new TFile(fileNames[tI].c_str(), "READ");

    TIter next(inFile_p[tI]->GetListOfKeys());
    TKey* key = NULL;

    while((key = (TKey*)next())){
      const std::string name = key->GetName();
      const std::string className = key->GetClassName();
      if(className.find("TH1D") == std::string::npos && className.find("TDirectory") == std::string::npos) continue;


      if(className.find("TH1D") != std::string::npos){
	++(histNamesToCounts[name]);
	if(tI == 0) histNamesToFile1Contents[name] = (TH1D*)inFile_p[tI]->Get(name.c_str());
	else{
	  std::string name2 = name;
	  if(stringToRemoveFile2.size() != 0) name2.replace(name2.find(stringToRemoveFile2), stringToRemoveFile2.size(), "");
	  histNamesToFile2Contents[name2] = (TH1D*)inFile_p[tI]->Get(name.c_str());	
	}
      }
      else{
	dirs_p.push_back((TDirectoryFile*)inFile_p[tI]->Get(name.c_str()));
	TIter next3(dirs_p[dirs_p.size()-1]->GetListOfKeys());
	TKey* key3 = NULL;

	while((key3 = (TKey*)next3())){     
	  const std::string name3 = key3->GetName();
	  const std::string className3 = key3->GetClassName();
	  if(className3.find("TH1D") == std::string::npos) continue;	  

	  ++(histNamesToCounts[name3]);
	  if(tI == 0) histNamesToFile1Contents[name3] = (TH1D*)dirs_p[dirs_p.size()-1]->Get(name3.c_str());
	  else{
	    std::string name4 = name3;
	    if(stringToRemoveFile2.size() != 0) name4.replace(name4.find(stringToRemoveFile2), stringToRemoveFile2.size(), "");
	    histNamesToFile2Contents[name4] = (TH1D*)dirs_p[dirs_p.size()-1]->Get(name3.c_str());
	  }
	}
      }
    }

    if(tI == 0){
      for(auto const & val : histNamesToCounts){
	++totalKeys;
	if(val.second != 1){
	  std::cout << "WARNING: in file 1 \'" << inFileName1 << "\', \'" << val.first << "\' is appearing more than once!" << std::endl;
	}
      }
    }
    else{
      for(auto const & val : histNamesToCounts){
	if(val.second != 2){
	  std::cout << "WARNING: in file 2 \'" << inFileName2 << "\', \'" << val.first << "\' is appearing not twice (" << val.second << ")!" << std::endl;
	}
      }
    }
  }

  std::cout << "Total keys: " << totalKeys << std::endl;
  std::cout << "Checking all hist..." << std::endl;
  
  const ULong64_t nMaxBins = 1000;
  Double_t generalBins1[nMaxBins+1];
  Double_t generalBins2[nMaxBins+1];

  for(auto const & val : histNamesToFile1Contents){
    const Int_t nBins1 = val.second->GetNbinsX();
    TH1D* temp_h = histNamesToFile2Contents[val.first];
    const Int_t nBins2 = temp_h->GetNbinsX();

    if(nBins1 != nBins2){
      std::cout << "Bins have changed (" << nBins1 << "->" << nBins2 << ") between files for \'" << val.first << "\'" << std::endl;
      continue;
    }

    for(Int_t bIX = 0; bIX < nBins1+1; ++bIX){
      generalBins1[bIX] = val.second->GetBinLowEdge(bIX+1);
    }
    for(Int_t bIX = 0; bIX < nBins2+1; ++bIX){
      generalBins2[bIX] = temp_h->GetBinLowEdge(bIX+1);
    }

    bool binsAreGood = true;
    for(Int_t bIX = 0; bIX < nBins1+1; ++bIX){
      if(TMath::Abs(generalBins1[bIX] - generalBins2[bIX]) > 0.01){
	binsAreGood = false;
	break;
      }
    }

    if(!binsAreGood){
      std::cout << "Bins in \'" << val.first << "\' are different between files: " << std::endl;
      std::cout << " ";
      for(Int_t bIX = 0; bIX < nBins1+1; ++bIX){
	std::cout << generalBins1[bIX] << ", ";
      }
      std::cout << generalBins1[nBins1] << "." << std::endl;;
      std::cout << " ";
      for(Int_t bIX = 0; bIX < nBins2+1; ++bIX){
	std::cout << generalBins2[bIX] << ", ";
      }
      std::cout << generalBins2[nBins1] << "." << std::endl;;
      continue;
    }

    for(Int_t bIX = 0; bIX < nBins1; ++bIX){
      generalBins1[bIX] = val.second->GetBinContent(bIX+1);
    }
    for(Int_t bIX = 0; bIX < nBins2; ++bIX){
      generalBins2[bIX] = temp_h->GetBinContent(bIX+1);
    }

    binsAreGood = true;
    for(Int_t bIX = 0; bIX < nBins1; ++bIX){
      if(TMath::Abs(generalBins1[bIX] - generalBins2[bIX]) > 0.01){
	binsAreGood = false;
	break;
      }
    }

    if(!binsAreGood){
      std::cout << "Bin contents in \'" << val.first << "\' are different between files: " << std::endl;
      std::cout << " ";
      for(Int_t bIX = 0; bIX < nBins1; ++bIX){
	std::cout << generalBins1[bIX] << ", ";
      }
      std::cout << generalBins1[nBins1-1] << "." << std::endl;;
      std::cout << " ";
      for(Int_t bIX = 0; bIX < nBins2; ++bIX){
	std::cout << generalBins2[bIX] << ", ";
      }
      std::cout << generalBins2[nBins1-1] << "." << std::endl;;
      continue;
    }

    for(Int_t bIX = 0; bIX < nBins1; ++bIX){
      generalBins1[bIX] = val.second->GetBinError(bIX+1);
    }
    for(Int_t bIX = 0; bIX < nBins2; ++bIX){
      generalBins2[bIX] = temp_h->GetBinError(bIX+1);
    }

    binsAreGood = true;
    for(Int_t bIX = 0; bIX < nBins1; ++bIX){
      if(TMath::Abs(generalBins1[bIX] - generalBins2[bIX]) > 0.01){
	binsAreGood = false;
	break;
      }
    }

    if(!binsAreGood){
      std::cout << "Bin errors in \'" << val.first << "\' are different between files: " << std::endl;
      std::cout << " ";
      for(Int_t bIX = 0; bIX < nBins1; ++bIX){
	std::cout << generalBins1[bIX] << ", ";
      }
      std::cout << generalBins1[nBins1-1] << "." << std::endl;;
      std::cout << " ";
      for(Int_t bIX = 0; bIX < nBins2; ++bIX){
	std::cout << generalBins2[bIX] << ", ";
      }
      std::cout << generalBins2[nBins1-1] << "." << std::endl;;
      continue;
    }    
  }
  
  for(Int_t tI = 0; tI < nFiles; ++tI){
    inFile_p[tI]->Close();
    delete inFile_p[tI];
  }

  totalTimer.stop();
  std::cout << "Final timing report:" << std::endl;
  std::cout << " Total timeCPU: " << totalTimer.totalCPU() << std::endl;
  std::cout << " Total timeWall: " << totalTimer.totalWall() << std::endl;
  std::cout << "compareAllHist Complete!" << std::endl;
  
  return 0;
}

int main(int argc, char* argv[])
{
  if(argc < 3 || argc > 4){
    std::cout << "Usage: ./bin/compareAllHist.exe <inFileName1> <inFileName2> <stringToRemoveFile2>" << std::endl;
    return 1;
  }

  int retVal = 0;
  if(argc == 3) retVal += compareAllHist(argv[1], argv[2]);
  else if(argc == 4) retVal += compareAllHist(argv[1], argv[2], argv[3]);
  return retVal;
}
