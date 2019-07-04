//cpp depdencies
#include <iostream>
#include <string>
#include <vector>

//ROOT dependencies
#include "TCollection.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TH1I.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TH2I.h"
#include "TKey.h"
#include "TNamed.h"
#include "TTree.h"

//RooUnfold dependencies
#include "src/RooUnfoldResponse.h"

//Local dependencies
#include "MainAnalysis/include/cutPropagator.h"
#include "MainAnalysis/include/doLocalDebug.h"

//Non-local (Utility) dependencies
#include "Utility/include/checkMakeDir.h"
#include "Utility/include/cppWatch.h"
#include "Utility/include/returnRootFileContentsList.h"
#include "Utility/include/stringUtil.h"

int recursiveCombineDir(TFile* inFile_p, const std::string dirName, TFile* outFile_p, const Int_t nDir, TDirectory* dirs_p[])
{
  int retVal = 0;

  inFile_p->cd();
  TDirectoryFile* dir_p = (TDirectoryFile*)inFile_p->Get(dirName.c_str());

  Int_t dirPos = -1;
  Int_t maxVal = -1;
  for(Int_t dI = 0; dI < nDir; ++dI){
    std::string tempName = dirs_p[dI]->GetName();

    if(tempName.find(dirName) != std::string::npos && maxVal < (Int_t )tempName.size()){
      dirPos = dI;
      maxVal = tempName.size();
    }
  }

  inFile_p->cd();
  dir_p->cd();
  TIter next(dir_p->GetListOfKeys());
  TKey* key = NULL;

  while( (key = (TKey*)next()) ){
    const std::string name = key->GetName();
    const std::string className = key->GetClassName();

    if(className.find("TDirectory") != std::string::npos){retVal += recursiveCombineDir(inFile_p, dirName + "/" + name, outFile_p, nDir, dirs_p);}

    if(className.find("TH1I") != std::string::npos){
      TH1I* hist_p = (TH1I*)key->ReadObj();

      outFile_p->cd();
      dirs_p[dirPos]->cd();
      
      hist_p->Write("", TObject::kOverwrite);
      
      inFile_p->cd();
      dir_p->cd();
      //pretty clear memory leak w/o delete here
      delete hist_p;
    }
    else if(className.find("TH1F") != std::string::npos){
      TH1F* hist_p = (TH1F*)key->ReadObj();

      outFile_p->cd();
      dirs_p[dirPos]->cd();
      
      hist_p->Write("", TObject::kOverwrite);
      
      inFile_p->cd();
      dir_p->cd();
      //pretty clear memory leak w/o delete here
      delete hist_p;
    }
    else if(className.find("TH1D") != std::string::npos){
      TH1D* hist_p = (TH1D*)key->ReadObj();

      outFile_p->cd();
      dirs_p[dirPos]->cd();
      
      hist_p->Write("", TObject::kOverwrite);
      
      inFile_p->cd();
      dir_p->cd();
      //pretty clear memory leak w/o delete here
      delete hist_p;
    }
    else if(className.find("TH2I") != std::string::npos){
      TH2I* hist_p = (TH2I*)key->ReadObj();

      outFile_p->cd();
      dirs_p[dirPos]->cd();
      
      hist_p->Write("", TObject::kOverwrite);
      
      inFile_p->cd();
      dir_p->cd();
      //pretty clear memory leak w/o delete here
      delete hist_p;
    }
    else if(className.find("TH2F") != std::string::npos){
      TH2F* hist_p = (TH2F*)key->ReadObj();

      outFile_p->cd();
      dirs_p[dirPos]->cd();
      
      hist_p->Write("", TObject::kOverwrite);
      
      inFile_p->cd();
      dir_p->cd();
      //pretty clear memory leak w/o delete here
      delete hist_p;
    }
    else if(className.find("TH2D") != std::string::npos){
      TH2D* hist_p = (TH2D*)key->ReadObj();

      outFile_p->cd();
      dirs_p[dirPos]->cd();
      
      hist_p->Write("", TObject::kOverwrite);
      
      inFile_p->cd();
      dir_p->cd();
      //pretty clear memory leak w/o delete here
      delete hist_p;
    }
    else if(className.find("RooUnfoldResponse") != std::string::npos){
      RooUnfoldResponse* rooRes_p = (RooUnfoldResponse*)key->ReadObj();

      outFile_p->cd();
      dirs_p[dirPos]->cd();
      
      rooRes_p->Write("", TObject::kOverwrite);
      
      inFile_p->cd();
      dir_p->cd();
      //pretty clear memory leak w/o delete here
      delete rooRes_p;
    }
  }

  return retVal;
}

int combineFiles(const std::string outFileName, std::vector<std::string> fileList)
{
  cppWatch globalWatch;
  globalWatch.start();

  for(unsigned int fI = 0; fI < fileList.size(); ++fI){
    if(!checkFile(fileList.at(fI)) || fileList.at(fI).find(".root") == std::string::npos){
      std::cout << "File \'" << fileList.at(fI) << "\' is invalid. return 1" << std::endl;
      return 1;
    }
  }

  cutPropagator prevProp;
  prevProp.Clean();

  std::vector<std::string> combinedHistTagBayes;
  std::vector<int> combinedBestBayes;

  std::vector<std::string> combinedHistTagSVD;
  std::vector<int> combinedBestSVD;

  std::vector<std::string> tdirNameVect;
  std::vector<std::vector<std::string> > tdirNamePerFile;

  std::cout << "Combining " << fileList.size() << " files..." << std::endl;
  
  std::vector<std::string> jtAlgosStr;
  std::vector<double> minJtPt;
  std::vector<double> multiJtPt;
  std::vector<int> recoTruncPos;

  const Int_t nMaxJtAlgos = 10;
  const Int_t nMaxCentBins = 4;

  const Int_t nMaxPtBins = 500;
  Int_t nGenPtBins[nMaxJtAlgos][nMaxCentBins];
  Int_t nRecoPtBins[nMaxJtAlgos][nMaxCentBins];
  Double_t genPtBins[nMaxJtAlgos][nMaxCentBins][nMaxPtBins];
  Double_t recoPtBins[nMaxJtAlgos][nMaxCentBins][nMaxPtBins];
  
  for(Int_t jI = 0; jI < nMaxJtAlgos; ++jI){
    for(Int_t cI = 0; cI < nMaxCentBins; ++cI){
      nGenPtBins[jI][cI] = -1;
      nRecoPtBins[jI][cI] = -1;
    }
  }

  std::vector<std::string> globalAlgos;
  std::vector<int> centBinsLow;
  std::vector<int> centBinsHi;

  for(unsigned int fI = 0; fI < fileList.size(); ++fI){
    std::cout << " " << fI << "/" << fileList.size() << ": " << fileList.at(fI) << std::endl;
    TFile* inFile_p = new TFile(fileList.at(fI).c_str(), "READ");

    cutPropagator cutPropCurrent;
    cutPropCurrent.Clean();    
    cutPropCurrent.GetAllVarFromFile(inFile_p);
    
    std::vector<std::string> tempHistTagBayes = cutPropCurrent.GetHistTagBayes();
    std::vector<int> tempHistBestBayes = cutPropCurrent.GetHistBestBayes();

    std::vector<std::string> tempHistTagSVD = cutPropCurrent.GetHistTagSVD();
    std::vector<int> tempHistBestSVD = cutPropCurrent.GetHistBestSVD();

    if(tempHistTagBayes.size() != 0){
      combinedHistTagBayes.insert(combinedHistTagBayes.end(), tempHistTagBayes.begin(), tempHistTagBayes.end());
      combinedBestBayes.insert(combinedBestBayes.end(), tempHistBestBayes.begin(), tempHistBestBayes.end());
    }

    if(tempHistTagSVD.size() != 0){
      combinedHistTagSVD.insert(combinedHistTagSVD.end(), tempHistTagSVD.begin(), tempHistTagSVD.end());
      combinedBestSVD.insert(combinedBestSVD.end(), tempHistBestSVD.begin(), tempHistBestSVD.end());
    }
  
    std::vector<std::string> tempAlgos = cutPropCurrent.GetJtAlgos();
    std::vector<double> tempMinJtPt = cutPropCurrent.GetMinJtPtCut();
    std::vector<double> tempMultiJtPt = cutPropCurrent.GetMultiJtPtCut();
    std::vector<int> tempRecoTruncPos = cutPropCurrent.GetRecoTruncPos();
    centBinsLow = cutPropCurrent.GetCentBinsLow();
    centBinsHi = cutPropCurrent.GetCentBinsHi();

    for(unsigned int algoI = 0; algoI < tempAlgos.size(); ++algoI){
      bool isAlgoFound = false;

      globalAlgos.push_back(tempAlgos[algoI]);
      Int_t rVal = getRVal(tempAlgos[algoI]);

      for(unsigned int cI = 0; cI < centBinsLow.size(); ++cI){
	std::string centBinsStr = "Cent" + std::to_string(centBinsLow[cI]) + "to" + std::to_string(centBinsHi[cI]);
      
	if(nGenPtBins[globalAlgos.size()-1][cI] < 0){
	  nGenPtBins[globalAlgos.size()-1][cI] = cutPropCurrent.GetGenNBinsFromRValCent(rVal, centBinsStr);
	  cutPropCurrent.GetGenBinsFromRValCent(rVal, centBinsStr, genPtBins[globalAlgos.size()-1][cI]);
	}
	if(nRecoPtBins[globalAlgos.size()-1][cI] < 0){
	  nRecoPtBins[globalAlgos.size()-1][cI] = cutPropCurrent.GetRecoNBinsFromRValCent(rVal, centBinsStr);
	  cutPropCurrent.GetRecoBinsFromRValCent(rVal, centBinsStr, recoPtBins[globalAlgos.size()-1][cI]);
	}
      }

      for(unsigned int jtI = 0; jtI < jtAlgosStr.size(); ++jtI){
	if(isStrSame(jtAlgosStr[jtI], tempAlgos[algoI])){
	  isAlgoFound = true;
	  break;
	}
      }

      if(!isAlgoFound){
	jtAlgosStr.push_back(tempAlgos[algoI]);
	minJtPt.push_back(tempMinJtPt[algoI]);
	multiJtPt.push_back(tempMultiJtPt[algoI]);
	recoTruncPos.push_back(tempRecoTruncPos[algoI]);
      }      
    }

    for(unsigned int algoI = 0; algoI < globalAlgos.size(); ++algoI){
      Int_t rVal = getRVal(globalAlgos[algoI]);

      for(unsigned int cI = 0; cI < centBinsLow.size(); ++cI){
	std::string centBinsStr = "Cent" + std::to_string(centBinsLow[cI]) + "to" + std::to_string(centBinsHi[cI]);
	cutPropCurrent.SetGenPtBins(rVal, centBinsStr, prevProp.GetGenPtBins(rVal, centBinsStr));
	cutPropCurrent.SetRecoPtBins(rVal, centBinsStr, prevProp.GetRecoPtBins(rVal, centBinsStr));
      }
    }
    
    std::vector<std::string> classList;
    std::vector<std::string> contentsList = returnRootFileContentsList(inFile_p, "", "", -1, &(classList));

    tdirNamePerFile.push_back({});

    for(unsigned int i = 0; i < contentsList.size(); ++i){
      if(classList.at(i).find("TDirectory") == std::string::npos) continue;
      
      bool isSame = false;
      for(unsigned int j = 0; j < tdirNameVect.size(); ++j){
	if(isStrSame(contentsList.at(i), tdirNameVect.at(j))){
	  isSame = true;
	  break;
	}
      }

      if(!isSame) tdirNameVect.push_back(contentsList.at(i));

      tdirNamePerFile.at(tdirNamePerFile.size()-1).push_back(contentsList.at(i));
    }
  
    if(fI != 0){
      bool testProp = cutPropCurrent.CheckPropagatorsMatch(prevProp, true, true);
      if(!testProp){
	std::cout << "Warning: \'" << fileList.at(fI-1) << "\', and \'" << fileList.at(fI) << "\' doesn't match." << std::endl;
      }
    }
    else{
      prevProp.Clean();
      prevProp.GetAllVarFromFile(inFile_p);
    }

    inFile_p->Close();
    delete inFile_p;
  }

  unsigned int pos = 0;
  while(pos < tdirNameVect.size()){
    bool isSwap = false;
    for(unsigned int dI = pos+1; dI < tdirNameVect.size(); ++dI){
      if(tdirNameVect.at(pos).find("/") != std::string::npos && tdirNameVect.at(dI).find("/") == std::string::npos){
	std::string tempStr = tdirNameVect.at(pos);
        tdirNameVect.at(pos) = tdirNameVect.at(dI);
        tdirNameVect.at(dI) = tempStr;
        isSwap = true;
      }
    }

    for(unsigned int dI = pos+1; dI < tdirNameVect.size(); ++dI){
      if(tdirNameVect.at(dI).find("/") != std::string::npos) continue;

      if(getRVal(tdirNameVect.at(dI)) < getRVal(tdirNameVect.at(pos))){
	std::string tempStr = tdirNameVect.at(pos);
        tdirNameVect.at(pos) = tdirNameVect.at(dI);
        tdirNameVect.at(dI) = tempStr;
        isSwap = true;
      }
    }
    
    if(!isSwap) pos++;
  }

  pos = 0;
  while(pos < tdirNameVect.size()){
    bool isFound = false;
    for(unsigned int dI = pos+1; dI < tdirNameVect.size(); ++dI){
      if(tdirNameVect.at(pos).find(tdirNameVect.at(dI)) != std::string::npos){
	std::string tempStr = tdirNameVect.at(pos);
	tdirNameVect.at(pos) = tdirNameVect.at(dI);
	tdirNameVect.at(dI) = tempStr;
	isFound = true;
      }
    }

    if(!isFound) pos++;
  }

  std::cout << "Reproducing " << tdirNameVect.size() << " directories..." << std::endl;
  for(unsigned int i = 0; i < tdirNameVect.size(); ++i){
    std::cout << " " << i << "/" << tdirNameVect.size() << ": " << tdirNameVect.at(i) << std::endl;
  }

  const Int_t nDir = tdirNameVect.size();
  TDirectory* dirs_p[nDir];

  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");
  //Following two lines necessary else you get absolutely clobbered on deletion timing. See threads here:
  //https://root-forum.cern.ch/t/closing-root-files-is-dead-slow-is-this-a-normal-thing/5273/16
  //https://root-forum.cern.ch/t/tfile-speed/17549/25
  //Bizarre
  outFile_p->SetBit(TFile::kDevNull);
  TH1::AddDirectory(kFALSE);

  
  int cutDirPos = -1;
  int subDirPos = -1;
  int unfoldDirBayesPos = -1;
  int unfoldDirSVDPos = -1;
  
  for(Int_t dI = 0; dI < nDir; ++dI){
    outFile_p->cd();
    dirs_p[dI] = NULL;
    if(tdirNameVect.at(dI).find("/") == std::string::npos) dirs_p[dI] = (TDirectory*)outFile_p->mkdir(tdirNameVect.at(dI).c_str());
    else{
      Int_t prePos = -1;
      Int_t tempMaxVal = -1;
      for(Int_t dI2 = 0; dI2 < dI; ++dI2){
	if(tdirNameVect.at(dI).find(tdirNameVect.at(dI2)) != std::string::npos && tempMaxVal < (int)tdirNameVect.at(dI2).size()){
	  prePos = dI2;
	  tempMaxVal = tdirNameVect.at(dI2).size();
	}
      }
    
      std::string tempDirName = tdirNameVect.at(dI);
      tempDirName.replace(0, tempDirName.rfind("/")+1, "");
      std::cout << "Prepos: " << prePos << ", " << dirs_p[prePos]->GetName() << std::endl;
      dirs_p[prePos]->cd();
      dirs_p[dI] = (TDirectory*)dirs_p[prePos]->mkdir(tempDirName.c_str());
    }
    dirs_p[dI]->cd();

    if(isStrSame("cutDir", tdirNameVect.at(dI))) cutDirPos = dI;
    else if(tdirNameVect.at(dI).find("subDir") != std::string::npos) subDirPos = dI;
    else if(tdirNameVect.at(dI).find("unfoldDirBayes") != std::string::npos) unfoldDirBayesPos = dI;
    else if(tdirNameVect.at(dI).find("unfoldDirSVD") != std::string::npos) unfoldDirSVDPos = dI;
  }

  if(unfoldDirBayesPos >= 0){
    std::cout << "UnfoldDirBayesPos: " << unfoldDirBayesPos << std::endl;
    std::cout << " " << dirs_p[unfoldDirBayesPos]->GetName() << std::endl;
  }
  else if(unfoldDirSVDPos >= 0){
    std::cout << "UnfoldDirSVDPos: " << unfoldDirSVDPos << std::endl;
    std::cout << " " << dirs_p[unfoldDirSVDPos]->GetName() << std::endl;
  }
  else{
    std::cout << "No unfolding dir found warning" << std::endl;
  }

  outFile_p->cd();

  cppWatch writingWatch;
  writingWatch.start();

  int retVal = 0;

  std::cout << "Writing..." << std::endl;

  for(unsigned int fI = 0; fI < fileList.size(); ++fI){
    TFile* inFile_p = new TFile(fileList.at(fI).c_str(), "READ");

    TIter next(inFile_p->GetListOfKeys());
    TKey* key = NULL;
    
    while( (key = (TKey*)next() ) ){
      const std::string name = key->GetName();
      const std::string className = key->GetClassName();
      
      if(className.find("TDirectory") != std::string::npos){retVal += recursiveCombineDir(inFile_p, name, outFile_p, nDir, dirs_p);}
      
      if(className.find("TH1I") != std::string::npos){      
	TH1I* hist_p = (TH1I*)key->ReadObj();
      
	outFile_p->cd();
	hist_p->Write("", TObject::kOverwrite);
	//pretty clear memory leak w/o delete here
	delete hist_p;
      }
      else if(className.find("TH1F") != std::string::npos){      
	TH1F* hist_p = (TH1F*)key->ReadObj();
      
	outFile_p->cd();
	hist_p->Write("", TObject::kOverwrite);
	//pretty clear memory leak w/o delete here
	delete hist_p;
      }
      else if(className.find("TH1D") != std::string::npos){      
	TH1D* hist_p = (TH1D*)key->ReadObj();
      
	outFile_p->cd();
	hist_p->Write("", TObject::kOverwrite);
	//pretty clear memory leak w/o delete here
	delete hist_p;
      }
      else if(className.find("TH2I") != std::string::npos){      
	TH2I* hist_p = (TH2I*)key->ReadObj();
      
	outFile_p->cd();
	hist_p->Write("", TObject::kOverwrite);
	//pretty clear memory leak w/o delete here
	delete hist_p;
      }
      else if(className.find("TH2F") != std::string::npos){      
	TH2F* hist_p = (TH2F*)key->ReadObj();
      
	outFile_p->cd();
	hist_p->Write("", TObject::kOverwrite);
	//pretty clear memory leak w/o delete here
	delete hist_p;
      }
      else if(className.find("TH2D") != std::string::npos){      
	TH2D* hist_p = (TH2D*)key->ReadObj();
      
	outFile_p->cd();
	hist_p->Write("", TObject::kOverwrite);
	//pretty clear memory leak w/o delete here
	delete hist_p;
      }
      else if(className.find("RooUnfoldResponse") != std::string::npos){      
	RooUnfoldResponse* rooRes_p = (RooUnfoldResponse*)key->ReadObj();
      
	outFile_p->cd();
	rooRes_p->Write("", TObject::kOverwrite);
	//pretty clear memory leak w/o delete here
	delete rooRes_p;
      }
    }

    inFile_p->Close();
    delete inFile_p;
  }

  if(retVal != 0) std::cout << "Uhoh rec nonzero something bad" << std::endl;

  writingWatch.stop();
  
  std::cout << "Write time: " << writingWatch.totalWall() << std::endl;

  prevProp.SetNJtAlgos(jtAlgosStr.size());
  prevProp.SetJtAlgos(jtAlgosStr);
  prevProp.SetMinJtPtCut(minJtPt);
  prevProp.SetMultiJtPtCut(multiJtPt);
  prevProp.SetRecoTruncPos(recoTruncPos);

  prevProp.SetNHistDim(combinedHistTagBayes.size());
  prevProp.SetHistTagBayes(combinedHistTagBayes);
  prevProp.SetHistBestBayes(combinedBestBayes);

  prevProp.SetHistTagSVD(combinedHistTagSVD);
  prevProp.SetHistBestSVD(combinedBestSVD);


  for(unsigned int algoI = 0; algoI < globalAlgos.size(); ++algoI){
    Int_t rVal = getRVal(globalAlgos[algoI]);

    for(unsigned int cI = 0; cI < centBinsLow.size(); ++cI){
      std::string centBinsStr = "Cent" + std::to_string(centBinsLow[cI]) + "to" + std::to_string(centBinsHi[cI]);
      
      prevProp.SetGenPtBins(rVal, centBinsStr, nGenPtBins[algoI][cI], genPtBins[algoI][cI]);
      prevProp.SetRecoPtBins(rVal, centBinsStr, nRecoPtBins[algoI][cI], recoPtBins[algoI][cI]);
    }
  }

  std::cout << " Best bayes: " << combinedHistTagBayes.size() << ", " << combinedBestBayes.size() << std::endl;
  std::cout << "  unfoldDirBayesPos: " << unfoldDirBayesPos << std::endl;
  
  if(unfoldDirBayesPos < 0 || unfoldDirSVDPos < 0){
    if(!prevProp.WriteAllVarToFile(outFile_p, dirs_p[cutDirPos], dirs_p[subDirPos], NULL)) std::cout << "Warning: Cut writing has failed" << std::endl;
  }
  else{
    if(!prevProp.WriteAllVarToFile(outFile_p, dirs_p[cutDirPos], dirs_p[subDirPos], dirs_p[unfoldDirBayesPos], dirs_p[unfoldDirSVDPos])) std::cout << "Warning: Cut writing has failed" << std::endl;
  }

  outFile_p->Close();
  delete outFile_p;

  globalWatch.stop();

  std::cout << "Job complete." << std::endl;
  std::cout << "Total time: " << globalWatch.totalWall() << std::endl;

  return 0;
}


int main(int argc, char* argv[])
{
  if(argc < 3){
    std::cout << "Usage: ./bin/combineFiles.exe <outFileName> <Long> <List> <Of> <Input> <Files> ..." << std::endl;
    return 1;
  }

  std::vector<std::string> fileList;
  for(int i = 2; i < argc; ++i){
    fileList.push_back(argv[i]);
  }

  int retVal = 0;
  retVal += combineFiles(argv[1], fileList);
  return retVal;
}
