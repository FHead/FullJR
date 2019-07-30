#ifndef SCALEERRORTOOL_H
#define SCALEERRORTOOL_H

//cpp dependencies
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <fstream>

//ROOT dependencies
#include "TMath.h"
#include "TH2F.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TStyle.h"

//Local (Utility) dependencies
#include "Utility/include/checkMakeDir.h"
#include "Utility/include/doGlobalDebug.h"
#include "Utility/include/plotUtilities.h"
#include "Utility/include/enums.h"

class scaleErrorTool{
 public:
  scaleErrorTool();
  scaleErrorTool(std::string inErrorFile);
  ~scaleErrorTool(){};

  bool Init();
  bool Init(std::string inErrorFile);

  std::string getValidString(std::string inFullLine, std::vector<std::string> inCompVect);
  unsigned int getValidPos(std::string compStr, std::vector<std::string> inCompVect);
  std::string getValidStringFromPos(unsigned int pos, std::vector<std::string> inCompVect);
  unsigned int getKey(std::string fullLine);
  unsigned int getKey(std::string centStr, std::string absEtaStr, std::string rStr, std::string flowStr);
  unsigned int getKey(CentType centType, AbsEtaType absEtaType, RType rType, FlowType flowType);

  std::string getStringFromKey(unsigned int key);

  double getMuDataMinusMC(std::string centStr, std::string absEtaStr, std::string rStr, std::string flowStr);
  double getMuDataMinusMCErr(std::string centStr, std::string absEtaStr, std::string rStr, std::string flowStr);
  double getSigDataMinusMC(std::string centStr, std::string absEtaStr, std::string rStr, std::string flowStr);
  double getSigDataMinusMCErr(std::string centStr, std::string absEtaStr, std::string rStr, std::string flowStr);

  double getMuDataMinusMC(int centVal, double absEtaVal, int rVal, std::string flowStr);
  double getMuDataMinusMCErr(int centVal, double absEtaVal, int rVal, std::string flowStr);
  double getSigDataMinusMC(int centVal, double absEtaVal, int rVal, std::string flowStr);
  double getSigDataMinusMCErr(int centVal, double absEtaVal, int rVal, std::string flowStr);

  double getMuDataMinusMC(int centVal, double absEtaVal, int rVal, FlowType flowType);
  double getMuDataMinusMCErr(int centVal, double absEtaVal, int rVal, FlowType flowType);
  double getSigDataMinusMC(int centVal, double absEtaVal, int rVal, FlowType flowType);
  double getSigDataMinusMCErr(int centVal, double absEtaVal, int rVal, FlowType flowType);

  double getMuDataMinusMC(CentType centType, AbsEtaType absEtaType, RType rType, FlowType flowType);
  double getMuDataMinusMCErr(CentType centType, AbsEtaType absEtaType, RType rType, FlowType flowType);
  double getSigDataMinusMC(CentType centType, AbsEtaType absEtaType, RType rType, FlowType flowType);
  double getSigDataMinusMCErr(CentType centType, AbsEtaType absEtaType, RType rType, FlowType flowType);

  double getMuDataMinusMC(std::string fullLine);
  double getMuDataMinusMCErr(std::string fullLine);
  double getSigDataMinusMC(std::string fullLine);
  double getSigDataMinusMCErr(std::string fullLine);

  CentType   getCentType(std::string centStr)     { return static_cast<CentType>(getValidPos(centStr, centStrV)); }
  AbsEtaType getAbsEtaType(std::string absEtaStr) { return static_cast<AbsEtaType>(getValidPos(absEtaStr, absEtaStrV)); }
  RType      getRType(std::string rStr)           { return static_cast<RType>(getValidPos(rStr, rStrV)); }
  FlowType   getFlowType(std::string flowStr)     { return static_cast<FlowType>(getValidPos(flowStr, flowStrV)); }

  void makeErrorMaps();

  std::string errorFile_;

  //  std::vector<std::string> centStrV = {"0-10%", "10-30%", "30-50%", "50-100%"};//0 thru 3
  std::vector<std::string> centStrV = {"0-5%", "5-10%", "10-20%", "20-30%", "30-40%", "40-50%", "50-70%", "70-100%"};//0 thru 7
  std::vector<std::string> absEtaStrV = {"AbsEta0p0to1p0", "AbsEta1p0to2p0"};//0 thru 1
  std::vector<std::string> rStrV = {"R0p3", "R0p4", "R0p6", "R0p8", "R1p0"}; //0 thru 4
  std::vector<std::string> flowStrV = {"NoFlow", "FlowDefaultInRho"};//0 thru 1

  //  std::vector<int> centValLowV = {0, 10, 30, 50};//0 thru 3
  //  std::vector<int> centValHiV = {10, 30, 50, 90};//0 thru 3
  std::vector<int> centValLowV = {0,  5, 10, 20, 30, 40, 50, 70};//0 thru 7
  std::vector<int> centValHiV = { 5, 10, 20, 30, 40, 50, 70, 100};//0 thru 7
  std::vector<double> absEtaValLowV = {0., 1.};//0 thru 1
  std::vector<double> absEtaValHiV = {1., 2.};//0 thru 1
  std::vector<int> rValV = {3, 4, 6, 8, 10}; //0 thru 4

  std::map<unsigned int, double> muDataMinusMC;
  std::map<unsigned int, double> sigDataMinusMC;
  std::map<unsigned int, double> muDataMinusMCErr;
  std::map<unsigned int, double> sigDataMinusMCErr;

  Bool_t initialized = false;
};

scaleErrorTool::scaleErrorTool()
{
  errorFile_ = "";
  return;
}

scaleErrorTool::scaleErrorTool(std::string inErrorFile)
{
  errorFile_ = inErrorFile;
  return;
}

bool scaleErrorTool::Init()
{
  if(errorFile_.size() == 0  || errorFile_.find(".txt") == std::string::npos){
    std::cout << "SCALEERRORTOOL: Error file string \'" << errorFile_ << "\' is invalid. return false" << std::endl;
    initialized = false;
    return initialized;
  }

  std::ifstream file(errorFile_.c_str());
  std::string tempStr;
  while(std::getline(file, tempStr)){
    if(tempStr.find("AbsEta0p0to2p0") != std::string::npos) continue;
    if(tempStr.size() == 0) continue;
    
    unsigned int keyVal = getKey(tempStr);

    bool mcMinusData = true;
    if(tempStr.find("Data,MC") != std::string::npos) mcMinusData = false;

    tempStr.replace(0, tempStr.find("=")+1, "");
    double nominalMu = std::stod(tempStr.substr(0,tempStr.find("#")));
    if(mcMinusData) nominalMu *= -1;
    muDataMinusMC[keyVal] = nominalMu;
    tempStr.replace(0, tempStr.find("#pm")+3, "");
    muDataMinusMCErr[keyVal] = std::stod(tempStr.substr(0,tempStr.find(",")));

    tempStr.replace(0, tempStr.find("=")+1, "");
    double nominalSig = std::stod(tempStr.substr(0,tempStr.find("#")));
    if(mcMinusData) nominalSig *= -1;
    sigDataMinusMC[keyVal] = nominalSig;
    tempStr.replace(0, tempStr.find("#pm")+3, "");
    sigDataMinusMCErr[keyVal] = std::stod(tempStr);
  }

  std::cout << "SCALE TOOL ERROR INIT: " << std::endl;
  for(std::map<unsigned int, double>::iterator mI = muDataMinusMC.begin(); mI != muDataMinusMC.end(); ++mI){
    std::cout << " " << mI->first << ", " << getStringFromKey(mI->first) << ", " << mI->second << std::endl;
  }

  file.close();

  initialized = true;
  return initialized;
}


std::string scaleErrorTool::getValidString(std::string inFullLine, std::vector<std::string> inCompVect)
{
  std::string retStr = "";

  for(unsigned int i = 0; i < inCompVect.size(); ++i){
    if(inFullLine.find(inCompVect.at(i)) != std::string::npos){
      retStr = inCompVect.at(i);
      break;
    }
  }

  return retStr;
}

unsigned int scaleErrorTool::getValidPos(std::string compStr, std::vector<std::string> inCompVect)
{
  unsigned int retVal = 9;

  for(unsigned int i = 0; i < inCompVect.size(); ++i){
    if(compStr.find(inCompVect.at(i)) != std::string::npos && compStr.size() == inCompVect.at(i).size()){
      retVal = i;
      break;
    }
  }

  return retVal;
}


std::string scaleErrorTool::getValidStringFromPos(unsigned int pos, std::vector<std::string> inCompVect)
{
  return inCompVect.at(pos);
}


unsigned int scaleErrorTool::getKey(std::string fullLine)
{
  std::string centStr = getValidString(fullLine, centStrV);
  std::string absEtaStr = getValidString(fullLine, absEtaStrV);
  std::string rStr = getValidString(fullLine, rStrV);
  std::string flowStr = getValidString(fullLine, flowStrV);
  
  unsigned int retKey = 99999;
  
  if(centStr.size() == 0 || absEtaStr.size() == 0 || rStr.size() == 0 || flowStr.size() == 0){
    std::cout << "Missing key in line \'" << fullLine << "\'. Init return 9999 key" << std::endl;
    std::cout << " Cent: " << centStr << std::endl;
    std::cout << " AbsEta: " << absEtaStr << std::endl;
    std::cout << " R: " << rStr << std::endl;
    std::cout << " Flow: " << flowStr << std::endl;
  }
  else retKey = getKey(centStr, absEtaStr, rStr, flowStr);

  return retKey;
}

unsigned int scaleErrorTool::getKey(std::string centStr, std::string absEtaStr, std::string rStr, std::string flowStr)
{
  return getKey(getCentType(centStr), getAbsEtaType(absEtaStr), getRType(rStr), getFlowType(flowStr));
}
  
unsigned int scaleErrorTool::getKey(CentType centType, AbsEtaType absEtaType, RType rType, FlowType flowType)
{
  unsigned int retKey = 0;
  retKey += centType;
  retKey += 10 * absEtaType;
  retKey += 100 * rType;
  retKey += 1000 * flowType;

  return retKey;
}


std::string scaleErrorTool::getStringFromKey(unsigned int key)
{
  std::string keyStr = getValidStringFromPos(key%10, centStrV);
  key /= 10;
  keyStr = keyStr + ", " + getValidStringFromPos(key%10, absEtaStrV);
  key /= 10;
  keyStr = keyStr + ", " + getValidStringFromPos(key%10, rStrV);
  key /= 10;
  keyStr = keyStr + ", " + getValidStringFromPos(key%10, flowStrV);

  return keyStr;
}


double scaleErrorTool::getMuDataMinusMC(std::string centStr, std::string absEtaStr, std::string rStr, std::string flowStr)
{
  unsigned int key = getKey(centStr, absEtaStr, rStr, flowStr);
  return muDataMinusMC[key];
}

double scaleErrorTool::getMuDataMinusMCErr(std::string centStr, std::string absEtaStr, std::string rStr, std::string flowStr)
{
  unsigned int key = getKey(centStr, absEtaStr, rStr, flowStr);
  return muDataMinusMCErr[key];
}

double scaleErrorTool::getSigDataMinusMC(std::string centStr, std::string absEtaStr, std::string rStr, std::string flowStr)
{
  unsigned int key = getKey(centStr, absEtaStr, rStr, flowStr);
  return sigDataMinusMC[key];
}

double scaleErrorTool::getSigDataMinusMCErr(std::string centStr, std::string absEtaStr, std::string rStr, std::string flowStr)
{
  unsigned int key = getKey(centStr, absEtaStr, rStr, flowStr);
  return sigDataMinusMCErr[key];
}

double scaleErrorTool::getMuDataMinusMC(int centVal, double absEtaVal, int rVal, std::string flowStr)
{
   FlowType flowType = getFlowType(flowStr);
   return getMuDataMinusMC(centVal, absEtaVal, rVal, flowType);
}

double scaleErrorTool::getMuDataMinusMCErr(int centVal, double absEtaVal, int rVal, std::string flowStr)
{
   FlowType flowType = getFlowType(flowStr);
   return getMuDataMinusMCErr(centVal, absEtaVal, rVal, flowType);
}

double scaleErrorTool::getSigDataMinusMC(int centVal, double absEtaVal, int rVal, std::string flowStr)
{
   FlowType flowType = getFlowType(flowStr);
   return getSigDataMinusMC(centVal, absEtaVal, rVal, flowType);
}

double scaleErrorTool::getSigDataMinusMCErr(int centVal, double absEtaVal, int rVal, std::string flowStr)
{
   FlowType flowType = getFlowType(flowStr);
   return getSigDataMinusMCErr(centVal, absEtaVal, rVal, flowType);
}
  
double scaleErrorTool::getMuDataMinusMC(int centVal, double absEtaVal, int rVal, FlowType flowType)
{
   CentType centType;
   AbsEtaType absEtaType;
   RType rType;

   absEtaVal = TMath::Abs(absEtaVal);

   for(unsigned int cI = 0; cI < centValLowV.size(); ++cI){
      if(centValLowV.at(cI) <= centVal && centVal < centValHiV.at(cI)){
         centType = static_cast<CentType>(cI);
         break;
      }
   }

   for(unsigned int cI = 0; cI < absEtaValLowV.size(); ++cI){
      if(absEtaValLowV.at(cI) <= absEtaVal && absEtaVal < absEtaValHiV.at(cI)){
         absEtaType = static_cast<AbsEtaType>(cI);
         break;
      }
   }

   for(unsigned int cI = 0; cI < rValV.size(); ++cI){
      if(rValV.at(cI) == rVal){
         rType = static_cast<RType>(cI);
         break;
      }
   }

   return getMuDataMinusMC(centType, absEtaType, rType, flowType);
}

double scaleErrorTool::getMuDataMinusMCErr(int centVal, double absEtaVal, int rVal, FlowType flowType)
{
   CentType centType;
   AbsEtaType absEtaType;
   RType rType;

   absEtaVal = TMath::Abs(absEtaVal);

   for(unsigned int cI = 0; cI < centValLowV.size(); ++cI){
      if(centValLowV.at(cI) <= centVal && centVal < centValHiV.at(cI)){
         centType = static_cast<CentType>(cI);
         break;
      }
   }

   for(unsigned int cI = 0; cI < absEtaValLowV.size(); ++cI){
      if(absEtaValLowV.at(cI) <= absEtaVal && absEtaVal < absEtaValHiV.at(cI)){
         absEtaType = static_cast<AbsEtaType>(cI);
         break;
      }
   }

   for(unsigned int cI = 0; cI < rValV.size(); ++cI){
      if(rValV.at(cI) == rVal){
         rType = static_cast<RType>(cI);
         break;
      }
   }

   return getMuDataMinusMCErr(centType, absEtaType, rType, flowType);
}

double scaleErrorTool::getSigDataMinusMC(int centVal, double absEtaVal, int rVal, FlowType flowType)
{
   CentType centType;
   AbsEtaType absEtaType;
   RType rType;

   absEtaVal = TMath::Abs(absEtaVal);

   for(unsigned int cI = 0; cI < centValLowV.size(); ++cI){
      if(centValLowV.at(cI) <= centVal && centVal < centValHiV.at(cI)){
         centType = static_cast<CentType>(cI);
         break;
      }
   }

   for(unsigned int cI = 0; cI < absEtaValLowV.size(); ++cI){
      if(absEtaValLowV.at(cI) <= absEtaVal && absEtaVal < absEtaValHiV.at(cI)){
         absEtaType = static_cast<AbsEtaType>(cI);
         break;
      }
   }

   for(unsigned int cI = 0; cI < rValV.size(); ++cI){
      if(rValV.at(cI) == rVal){
         rType = static_cast<RType>(cI);
         break;
      }
   }

   return getSigDataMinusMC(centType, absEtaType, rType, flowType);
}

double scaleErrorTool::getSigDataMinusMCErr(int centVal, double absEtaVal, int rVal, FlowType flowType)
{
   CentType centType;
   AbsEtaType absEtaType;
   RType rType;

   absEtaVal = TMath::Abs(absEtaVal);

   for(unsigned int cI = 0; cI < centValLowV.size(); ++cI){
      if(centValLowV.at(cI) <= centVal && centVal < centValHiV.at(cI)){
         centType = static_cast<CentType>(cI);
         break;
      }
   }

   for(unsigned int cI = 0; cI < absEtaValLowV.size(); ++cI){
      if(absEtaValLowV.at(cI) <= absEtaVal && absEtaVal < absEtaValHiV.at(cI)){
         absEtaType = static_cast<AbsEtaType>(cI);
         break;
      }
   }

   for(unsigned int cI = 0; cI < rValV.size(); ++cI){
      if(rValV.at(cI) == rVal){
         rType = static_cast<RType>(cI);
         break;
      }
   }

   return getSigDataMinusMCErr(centType, absEtaType, rType, flowType);
}

double scaleErrorTool::getMuDataMinusMC(CentType centType, AbsEtaType absEtaType, RType rType, FlowType flowType)
{
  unsigned int key = getKey(centType, absEtaType, rType, flowType);
  return muDataMinusMC[key];
}

double scaleErrorTool::getMuDataMinusMCErr(CentType centType, AbsEtaType absEtaType, RType rType, FlowType flowType)
{
  unsigned int key = getKey(centType, absEtaType, rType, flowType);
  return muDataMinusMCErr[key];
}

double scaleErrorTool::getSigDataMinusMC(CentType centType, AbsEtaType absEtaType, RType rType, FlowType flowType)
{
  unsigned int key = getKey(centType, absEtaType, rType, flowType);
  return sigDataMinusMC[key];
}

double scaleErrorTool::getSigDataMinusMCErr(CentType centType, AbsEtaType absEtaType, RType rType, FlowType flowType)
{
  unsigned int key = getKey(centType, absEtaType, rType, flowType);
  return sigDataMinusMCErr[key];
}

double scaleErrorTool::getMuDataMinusMC(std::string fullLine)
{
  std::string centStr = getValidString(fullLine, centStrV);
  std::string absEtaStr = getValidString(fullLine, absEtaStrV);
  std::string rStr = getValidString(fullLine, rStrV);
  std::string flowStr = getValidString(fullLine, flowStrV);

  unsigned int key = getKey(centStr, absEtaStr, rStr, flowStr);
  return muDataMinusMC[key];
}

double scaleErrorTool::getMuDataMinusMCErr(std::string fullLine)
{
  std::string centStr = getValidString(fullLine, centStrV);
  std::string absEtaStr = getValidString(fullLine, absEtaStrV);
  std::string rStr = getValidString(fullLine, rStrV);
  std::string flowStr = getValidString(fullLine, flowStrV);

  unsigned int key = getKey(centStr, absEtaStr, rStr, flowStr);
  return muDataMinusMCErr[key];
}

double scaleErrorTool::getSigDataMinusMC(std::string fullLine)
{
  std::string centStr = getValidString(fullLine, centStrV);
  std::string absEtaStr = getValidString(fullLine, absEtaStrV);
  std::string rStr = getValidString(fullLine, rStrV);
  std::string flowStr = getValidString(fullLine, flowStrV);

  unsigned int key = getKey(centStr, absEtaStr, rStr, flowStr);
  return sigDataMinusMC[key];
}

double scaleErrorTool::getSigDataMinusMCErr(std::string fullLine)
{
  std::string centStr = getValidString(fullLine, centStrV);
  std::string absEtaStr = getValidString(fullLine, absEtaStrV);
  std::string rStr = getValidString(fullLine, rStrV);
  std::string flowStr = getValidString(fullLine, flowStrV);

  unsigned int key = getKey(centStr, absEtaStr, rStr, flowStr);
  return sigDataMinusMCErr[key];
}


void scaleErrorTool::makeErrorMaps()
{
  if(!initialized){
    std::cout << "WARNING: scaleErrorTool is not initialized, cannot create error maps. return 1" << std::endl;
    return;
  }

  std::cout << "scaleErrorTool is creating error maps..." << std::endl;
  
  //  const Int_t nCentBins = 4;
  //  const Double_t centBins[nCentBins+1] = {0, 10, 30, 50, 100};
  const Int_t nCentBins = 7;
  const Double_t centBins[nCentBins+1] = {0, 5, 10, 20, 30, 40, 50, 70};
  const Int_t nRBins = 5;
  const Double_t rBins[nRBins+1] = {0.25, 0.35, 0.45, 0.7, 0.9, 1.1};
  const Double_t rBins2[nRBins] = {0.3, 0.4, 0.6, 0.8, 1.0};

  std::vector<std::string> typeVal = {"mu", "mu", "mu", "mu", "sigma", "sigma", "sigma", "sigma"};
  std::vector<std::string> etaVal = {"AbsEta0p0to1p0", "AbsEta1p0to2p0", "AbsEta0p0to1p0", "AbsEta1p0to2p0", "AbsEta0p0to1p0", "AbsEta1p0to2p0", "AbsEta0p0to1p0", "AbsEta1p0to2p0"};
  std::vector<std::string> etaLabel;
  for(unsigned int etaI = 0; etaI < etaVal.size(); ++etaI){
    std::string etaStr = etaVal.at(etaI);
    etaStr.replace(0,6,"");
    while(etaStr.find("p") != std::string::npos){
      etaStr.replace(etaStr.find("p"), 1, ".");     
    }
    etaStr.replace(etaStr.find("to"), 2, " < |#eta| < ");
    etaLabel.push_back(etaStr);
  }

  std::vector<std::string> flowVal = {"NoFlow", "NoFlow", "FlowDefaultInRho", "FlowDefaultInRho", "NoFlow", "NoFlow", "FlowDefaultInRho", "FlowDefaultInRho"};
  std::vector<std::string> flowLabel;
  for(unsigned int fI = 0; fI < flowVal.size(); ++fI){
    if(flowVal.at(fI).find("No") != std::string::npos) flowLabel.push_back("No Flow Corr.");
    else flowLabel.push_back("Flow Corr.");
  }

  std::vector<TH2F*> histVect;
  histVect.reserve(typeVal.size());
  for(unsigned int hI = 0; hI < typeVal.size(); ++hI){
    histVect.push_back(NULL);
    histVect.at(hI) = new TH2F((typeVal.at(hI) + "ErrorMap_" + etaVal.at(hI) + "_" + flowVal.at(hI) + "_h").c_str(), ";Centrality%;Jet or Random Cone R", nCentBins, centBins, nRBins, rBins);
    histVect.at(hI)->GetXaxis()->CenterTitle();
    histVect.at(hI)->GetYaxis()->CenterTitle();
  }

  Double_t maxValMu = -1;
  Double_t maxValSig = -1;

  for(unsigned int hI = 0; hI < histVect.size(); ++hI){
    for(Int_t bIX = 0; bIX < histVect.at(hI)->GetNbinsX(); ++bIX){
      const std::string centStr = std::to_string((int)centBins[bIX]) + "-" + std::to_string((int)centBins[bIX+1]) + "%";
      
      for(Int_t bIY = 0; bIY < histVect.at(hI)->GetNbinsY(); ++bIY){
	const std::string rStr = rStrV.at(bIY);
	
	Double_t val = -1;
	if(typeVal.at(hI).find("mu") != std::string::npos){
	  val = TMath::Abs(getMuDataMinusMC(centStr, etaVal.at(hI), rStr, flowVal.at(hI)));
	  if(val > maxValMu) maxValMu = val;
	}
	else{
	  val = TMath::Abs(getSigDataMinusMC(centStr, etaVal.at(hI), rStr, flowVal.at(hI)));
	  if(val > maxValSig) maxValSig = val;
	}
	
	histVect.at(hI)->SetBinContent(bIX+1, bIY+1, val);
      }
    }
  }

  TDatime* date = new TDatime();
  std::string dirName1 = "pdfDir/";
  std::string dirName2 = "pdfDir/" + std::to_string(date->GetDate()) + "/";
  checkMakeDir(dirName1);
  checkMakeDir(dirName2);

  TLatex* label_p = new TLatex();
  label_p->SetNDC();
  label_p->SetTextFont(43);
  label_p->SetTextSize(16);

  for(unsigned int hI = 0; hI < histVect.size(); ++hI){
    TCanvas* canv_p = new TCanvas("temp_p", "", 500, 400);
    canv_p->SetLeftMargin(canv_p->GetLeftMargin()*1.2);
    canv_p->SetBottomMargin(canv_p->GetLeftMargin());
    canv_p->SetTopMargin(canv_p->GetLeftMargin()*3./4.);
    canv_p->SetRightMargin(canv_p->GetLeftMargin());

    canv_p->cd();

    if(typeVal.at(hI).find("mu") != std::string::npos) histVect.at(hI)->SetMaximum(maxValMu);
    else histVect.at(hI)->SetMaximum(maxValSig);

    histVect.at(hI)->SetMinimum(0.0);

    histVect.at(hI)->GetYaxis()->SetLabelColor(0);
    histVect.at(hI)->SetMarkerSize(histVect.at(hI)->GetMarkerSize()*4./5.);
    histVect.at(hI)->SetTickLength(0.0, "xyz");
    histVect.at(hI)->DrawCopy("COLZ TEXT");
    gStyle->SetOptStat(0);
    gStyle->SetPaintTextFormat("1.3f");

    label_p->DrawLatex(.15, .95, ("|#Delta#" + typeVal.at(hI) + "|_{Data,MC}, " + etaLabel.at(hI) + ", " + flowLabel.at(hI)).c_str());

    label_p->SetNDC(0);
    label_p->SetTextSize(14);
    
    for(Int_t rI = 0; rI < nRBins; ++rI){
      label_p->DrawLatex(-4, rBins2[rI]-.025, prettyString(rBins2[rI], 1, false).c_str());
    }

    label_p->SetNDC();
    label_p->SetTextSize(16);


    std::string saveName = histVect.at(hI)->GetName();
    saveName = dirName2 + saveName + "_" + std::to_string(date->GetDate()) + ".pdf";
    quietSaveAs(canv_p, saveName);
    delete canv_p;
  }

  for(unsigned int hI = 0; hI < histVect.size(); ++hI){
    delete histVect.at(hI);
  }
  histVect.clear();
  delete date;
  delete label_p;

  std::cout << "Error maps complete." << std::endl;

  return;
}

#endif
