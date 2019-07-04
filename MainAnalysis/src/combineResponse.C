//cpp dependencies
#include <iostream>
#include <map>
#include <string>
#include <vector>

//ROOT dependencies
#include "TDatime.h"
#include "TDirectory.h"
#include "TDirectoryFile.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TH1I.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TH2I.h"
#include "TKey.h"
#include "TMath.h"

//Local dependencies
#include "MainAnalysis/include/cutPropagator.h"

//Non-local dependencies
#include "Utility/include/checkMakeDir.h"
#include "Utility/include/returnRootFileContentsList.h"
#include "Utility/include/stringUtil.h"

bool checkAllNames(std::vector<std::string> names1, std::vector<std::string> names2, std::string file1, std::string file2)
{
  for(auto const & name1: names1){
    bool name1Found = false;
    for(auto const & name2: names2){
      if(isStrSame(name1, name2)){
	name1Found = true;
	break;
      }
    }
    if(!name1Found){
      std::cout << "Hist \'" << name1 << "\' from file \'" << file1 << "\' is not found in file \'" << file2 << "\'" << std::endl;
      return false;
    }
  }
  return true;
}

void constructHist1D(TFile* outFile_p, std::vector<TH1D*>* hists_p, std::map<std::string, unsigned int>* histNameToPos, TFile* inFile_p, std::string dirName, bool doNew = false)
{
  TDirectory* dir_p = (TDirectory*)inFile_p->Get(dirName.c_str());
  TIter next(dir_p->GetListOfKeys());
  TKey* key=NULL;

  const Int_t nMaxBins = 400;
  Double_t binsX[nMaxBins+1];  

  while((key = (TKey*)next())){
    const std::string name = key->GetName();
    const std::string name2 = dirName + "/" + name;
    const std::string className = key->GetClassName();

    if(className.find("TH1D") == std::string::npos) continue;
    if(histNameToPos->count(name2) == 0) continue;

    unsigned int pos = (*histNameToPos)[name2];
    TH1D* tempHist_p = (TH1D*)key->ReadObj();

    if(doNew){
      std::string newName = name;
      while(newName.find("/") != std::string::npos){newName.replace(0, newName.find("/")+1, "");}

      const Int_t nBinsX = tempHist_p->GetNbinsX();

      if(nBinsX > nMaxBins){
	std::cout << "WARNING: nBinsX (" << nBinsX << ") less than nMaxBins " << nMaxBins << std::endl;
      }

      for(Int_t bIX = 0; bIX < nBinsX+1; ++bIX){
	binsX[bIX] = tempHist_p->GetXaxis()->GetBinLowEdge(bIX+1);
      }

      std::string title = ";" + std::string(tempHist_p->GetXaxis()->GetTitle()) + ";" + std::string(tempHist_p->GetYaxis()->GetTitle());

      outFile_p->cd();
      outFile_p->cd(dirName.c_str());

      (*hists_p)[pos] = new TH1D(newName.c_str(), title.c_str(), nBinsX, binsX);
    }

    for(Int_t bIX = 0; bIX < tempHist_p->GetNbinsX(); ++bIX){
      Double_t val = (*hists_p)[pos]->GetBinContent(bIX+1);
      Double_t err = (*hists_p)[pos]->GetBinError(bIX+1);
      
      val += tempHist_p->GetBinContent(bIX+1);
      err = TMath::Sqrt(err*err + tempHist_p->GetBinError(bIX+1)*tempHist_p->GetBinError(bIX+1));
      
      (*hists_p)[pos]->SetBinContent(bIX+1, val);
      (*hists_p)[pos]->SetBinError(bIX+1, err);
    } 

    delete tempHist_p;
  }

  delete dir_p;

  return;
}


void constructHist2D(TFile* outFile_p, std::vector<TH2D*>* hists_p, std::map<std::string, unsigned int>* histNameToPos, TFile* inFile_p, std::string dirName, bool doNew = false)
{
  TDirectory* dir_p = (TDirectory*)inFile_p->Get(dirName.c_str());
  TIter next(dir_p->GetListOfKeys());
  TKey* key=NULL;

  const Int_t nMaxBins = 400;
  Double_t binsX[nMaxBins+1];
  Double_t binsY[nMaxBins+1];  

  while((key = (TKey*)next())){
    const std::string name = key->GetName();
    const std::string name2 = dirName + "/" + name;
    const std::string className = key->GetClassName();

    if(className.find("TH2D") == std::string::npos) continue;
    if(histNameToPos->count(name2) == 0) continue;

    unsigned int pos = (*histNameToPos)[name2];
    TH2D* tempHist_p = (TH2D*)key->ReadObj();
      
    if(doNew){
      std::string newName = name;
      while(newName.find("/") != std::string::npos){newName.replace(0, newName.find("/")+1, "");}
      
      const Int_t nBinsX = tempHist_p->GetNbinsX();
      const Int_t nBinsY = tempHist_p->GetNbinsY();

      if(nBinsX > nMaxBins || nBinsY > nMaxBins){
	std::cout << "WARNING: nBinsX or nBinsY (" << nBinsX << " or " << nBinsY << ") less than nMaxBins " << nMaxBins << std::endl;
      }

      for(Int_t bIX = 0; bIX < nBinsX+1; ++bIX){
	binsX[bIX] = tempHist_p->GetXaxis()->GetBinLowEdge(bIX+1);
      }

      for(Int_t bIY = 0; bIY < nBinsY+1; ++bIY){
	binsY[bIY] = tempHist_p->GetYaxis()->GetBinLowEdge(bIY+1);
      }
      std::string title = ";" + std::string(tempHist_p->GetXaxis()->GetTitle()) + ";" + std::string(tempHist_p->GetYaxis()->GetTitle());
      outFile_p->cd();    
      outFile_p->cd(dirName.c_str());    
      
      (*hists_p)[pos] = new TH2D(newName.c_str(), title.c_str(), nBinsX, binsX, nBinsY, binsY);
    }

    for(Int_t bIX = 0; bIX < tempHist_p->GetNbinsX(); ++bIX){
      for(Int_t bIY = 0; bIY < tempHist_p->GetNbinsY(); ++bIY){
	Double_t val = (*hists_p)[pos]->GetBinContent(bIX+1, bIY+1);
	Double_t err = (*hists_p)[pos]->GetBinError(bIX+1, bIY+1);

	val += tempHist_p->GetBinContent(bIX+1, bIY+1);
	err = TMath::Sqrt(err*err + tempHist_p->GetBinError(bIX+1, bIY+1)*tempHist_p->GetBinError(bIX+1, bIY+1));

	(*hists_p)[pos]->SetBinContent(bIX+1, bIY+1, val);
	(*hists_p)[pos]->SetBinError(bIX+1, bIY+1, err);
      }
    }
    delete tempHist_p;
  }

  delete dir_p;
  
  return;
}

int combineResponse(std::string inFileNames, std::string tagStr)
{
  TDatime* date = new TDatime();
  const std::string dateStr = std::to_string(date->GetDate());
  delete date;

  const std::string permaFileNames = inFileNames;

  if(inFileNames.size() != 0){
    if(inFileNames.substr(inFileNames.size()-1, 1).find(",") == std::string::npos) inFileNames = inFileNames + ",";
  }

  std::vector<std::string> fileNames;
  while(inFileNames.find(",") != std::string::npos){
    std::string tempFileName = inFileNames.substr(0, inFileNames.find(","));
    if(!checkFile(tempFileName)) continue;
    if(tempFileName.find(".root") == std::string::npos) continue;
    fileNames.push_back(tempFileName);
    inFileNames.replace(0, inFileNames.find(",")+1, "");
  }

  if(fileNames.size() == 0){
    std::cout << "Given inputs give no valid files. \'" << permaFileNames << "\'. return 1." << std::endl;
    return 1;
  }

  std::string outFileName = "output/" + dateStr + "/combinedResponse_" + tagStr + "_";
  if(checkFile(outFileName + dateStr + ".root")) outFileName = outFileName + "UPDATED_";
  outFileName = outFileName + dateStr + ".root";

  checkMakeDir("output");
  checkMakeDir("output/" + dateStr);  

  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");
  outFile_p->SetBit(TFile::kDevNull);
  TH1::AddDirectory(kFALSE);

  std::vector<TH1D*> hists1D_p;
  std::vector<TH2D*> hists2D_p;
  std::map<std::string, unsigned int> histNameToPos1D;
  std::map<std::string, unsigned int> histNameToPos2D;

  std::map<std::string, std::string> histNameToDirName;

  ULong64_t totalNResponseTH1 = 0;
  ULong64_t totalNResponseTH2 = 0;
  TFile* inFile_p = new TFile(fileNames[0].c_str(), "READ");

  TDirectory* cutDir_p = (TDirectory*)outFile_p->mkdir("cutDir");
  TDirectory* subFileDir_p = (TDirectory*)outFile_p->mkdir("cutDir/fullFiles");
  std::vector<TDirectory*> outDirs_p;
  cutPropagator cutProp;
  cutProp.Clean();
  cutProp.GetAllVarFromFile(inFile_p);  

  std::vector<std::string> responseNamesTH1 = returnRootFileContentsList(inFile_p, "TH1D", "");
  std::vector<std::string> responseNamesTH2 = returnRootFileContentsList(inFile_p, "TH2D", "");

  inFile_p->Close();
  delete inFile_p;

  std::cout << "N responseNamesTH1: " << responseNamesTH1.size() << std::endl;
  for(unsigned int nI = 0; nI < responseNamesTH1.size(); ++nI){
    ++totalNResponseTH1;
    histNameToPos1D[responseNamesTH1[nI]] = nI;

    std::string histName = responseNamesTH1[nI].substr(responseNamesTH1[nI].find("/")+1, responseNamesTH1[nI].size());
    std::string dirName = responseNamesTH1[nI].substr(0, responseNamesTH1[nI].find("/"));
    histNameToDirName[histName] = dirName;
  }
  hists1D_p.reserve(totalNResponseTH1);
  for(unsigned int nI = 0; nI < totalNResponseTH1; ++nI){
    hists1D_p.push_back(NULL);
  }

  std::cout << "N responseNamesTH2: " << responseNamesTH2.size() << std::endl;
  for(unsigned int nI = 0; nI < responseNamesTH2.size(); ++nI){
    ++totalNResponseTH2;
    histNameToPos2D[responseNamesTH2[nI]] = nI;

    std::string histName = responseNamesTH2[nI].substr(responseNamesTH2[nI].find("/")+1, responseNamesTH2[nI].size());
    std::string dirName = responseNamesTH2[nI].substr(0, responseNamesTH2[nI].find("/"));
    histNameToDirName[histName] = dirName;
  }
  hists2D_p.reserve(totalNResponseTH2);
  for(unsigned int nI = 0; nI < totalNResponseTH2; ++nI){
    hists2D_p.push_back(NULL);
  }

  for(unsigned int fI = 0; fI < fileNames.size(); ++fI){
    std::cout << "Processing file " << fI << "/" << fileNames.size() << ", \'" << fileNames[fI] << "\'.." << std::endl;
    inFile_p = new TFile(fileNames[fI].c_str(), "READ");
    
    std::vector<std::string> dirNames = returnRootFileContentsList(inFile_p, "TDirectory", "");
    std::vector<std::string> dirNames2 = returnRootFileContentsList(inFile_p, "TDirectoryFile", "");
    dirNames.insert(std::end(dirNames), std::begin(dirNames2), std::end(dirNames2));

    for(auto const& dir : dirNames){
      if(dir.find("cutDir") != std::string::npos) continue;

      if(fI == 0){
	outDirs_p.push_back(outFile_p->mkdir(dir.c_str()));

	constructHist1D(outFile_p, &hists1D_p, &histNameToPos1D, inFile_p, dir, true);
	constructHist2D(outFile_p, &hists2D_p, &histNameToPos2D, inFile_p, dir, true);
      }
      else{
	constructHist1D(outFile_p, &hists1D_p, &histNameToPos1D, inFile_p, dir, false);
	constructHist2D(outFile_p, &hists2D_p, &histNameToPos2D, inFile_p, dir, false);
      }
    }

    inFile_p->Close();
    delete inFile_p;
  }

  outFile_p->cd();

  for(unsigned int hI = 0; hI < hists1D_p.size(); ++hI){
    outFile_p->cd();
    outFile_p->cd(histNameToDirName[std::string(hists1D_p[hI]->GetName())].c_str());
    hists1D_p[hI]->Write("", TObject::kOverwrite);
    delete hists1D_p[hI];
  }

  for(unsigned int hI = 0; hI < hists2D_p.size(); ++hI){
    outFile_p->cd();
    outFile_p->cd(histNameToDirName[std::string(hists2D_p[hI]->GetName())].c_str());
    hists2D_p[hI]->Write("", TObject::kOverwrite);
    delete hists2D_p[hI];
  }

  if(!cutProp.WriteAllVarToFile(outFile_p, cutDir_p, subFileDir_p)) std::cout << "Warning: Cut writing has failed" << std::endl;

  outFile_p->Close();
  delete outFile_p;

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 3){
    std::cout << "Usage: ./bin/combineResponse.exe <commaSeparatedInFileNames> <tagStr>" << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += combineResponse(argv[1], argv[2]);
  return retVal;
}
