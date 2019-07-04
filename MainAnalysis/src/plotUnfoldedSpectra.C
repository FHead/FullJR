//cpp dependencies
#include <iostream>
#include <string>
#include <vector>

//ROOT dependencies
#include "TBox.h"
#include "TCanvas.h"
#include "TCollection.h"
#include "TDatime.h"
#include "TDirectory.h"
#include "TError.h"
#include "TFile.h"
#include "TH1D.h"
#include "TKey.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TMath.h"
#include "TPad.h"
#include "TStyle.h"

//Local dependencies
#include "MainAnalysis/include/cutPropagator.h"
#include "MainAnalysis/include/doLocalDebug.h"
#include "MainAnalysis/include/smallOrLargeR.h"
#include "MainAnalysis/include/systFunctions.h"
#include "MainAnalysis/include/texSlideCreator.h"

//Non-local dependencies
#include "Utility/include/checkMakeDir.h"
#include "Utility/include/cppWatch.h"
#include "Utility/include/doGlobalDebug.h"
#include "Utility/include/getLinBins.h"
#include "Utility/include/getLogBins.h"
#include "Utility/include/histDefUtility.h"
#include "Utility/include/kirchnerPalette.h"
#include "Utility/include/lumiAndTAAUtil.h"
#include "Utility/include/returnRootFileContentsList.h"
#include "Utility/include/vanGoghPalette.h"

void doRatioPlotting(std::vector<TH1D*> inHist_p, int divPos, const std::string xStr, const std::string yStr, const int xMin, const int xMax, const std::string labelStr, std::string saveName)
{
  TDatime* date = new TDatime();
  const std::string dateStr = std::to_string(date->GetDate());
  delete date;

  TCanvas* canv_p = new TCanvas("convCanv_p", "", 450, 450);
  canv_p->SetTopMargin(0.01);
  canv_p->SetRightMargin(0.01);
  canv_p->SetLeftMargin(0.14);
  canv_p->SetBottomMargin(0.14);
  canv_p->cd();
 
  TLatex* label_p = new TLatex();
  label_p->SetTextFont(43);
  label_p->SetTextSize(16);
	  
  TH1D* tempHist_p = new TH1D("tempHist_p", (";" + xStr + ";" + yStr).c_str(), 10, xMin, xMax);

  centerTitles(tempHist_p);
  tempHist_p->Draw("HIST");

  Double_t maxVal = -1;
  Double_t minVal = 1000000;
	  
  const Int_t nTempStyle = 3;
  const Int_t tempStyle[nTempStyle] = {20, 21, 34};

  kirchnerPalette kPalette;  
  const Int_t nTempColor = 4;
  const Int_t tempColor[nTempColor] = {kPalette.getColor(0), kPalette.getColor(1), kPalette.getColor(2), kPalette.getColor(3)};

  const Int_t nClones = inHist_p.size();
  std::vector<TH1D*> clones_p;
  clones_p.reserve(nClones);
  for(Int_t bI = 0; bI < nClones; ++bI){
    clones_p.push_back((TH1D*)inHist_p.at(bI)->Clone(("clone_" + std::to_string(bI)).c_str()));
  }
	
  for(unsigned int bI = 0; bI < clones_p.size(); ++bI){
    clones_p.at(bI)->SetMarkerStyle(tempStyle[bI%nTempStyle]);
    clones_p.at(bI)->SetMarkerSize(0.8);
    clones_p.at(bI)->SetMarkerColor(tempColor[bI%nTempColor]);
    clones_p.at(bI)->SetLineColor(tempColor[bI%nTempColor]);
    
    if(bI == (unsigned int)divPos) continue;
    
    clones_p.at(bI)->Divide(clones_p.at(divPos));
    
    for(Int_t bIX = 0; bIX < clones_p.at(bI)->GetNbinsX(); ++bIX){
      Double_t binCenter = clones_p.at(bI)->GetBinCenter(bIX+1);
      
      clones_p.at(bI)->SetBinError(bIX+1, 0.0);
      if(binCenter >= xMin) continue;
      if(binCenter <= xMax) continue;
      
      clones_p.at(bI)->SetBinContent(bIX+1, 0.0);
    }
    
    Double_t tempMaxVal = clones_p.at(bI)->GetMaximum();
    Double_t tempMinVal = getMinGTZero(clones_p.at(bI));
    if(tempMaxVal > maxVal) maxVal = tempMaxVal;
    if(tempMinVal < minVal) minVal = tempMinVal;
  }
  
  clones_p.at(divPos)->Divide(clones_p.at(divPos));
  
  for(Int_t bIX = 0; bIX < clones_p.at(divPos)->GetNbinsX(); ++bIX){
    Double_t binCenter = clones_p.at(divPos)->GetBinCenter(bIX+1);
    
    clones_p.at(divPos)->SetBinError(bIX+1, 0.0);
    if(binCenter >= 200.) continue;
    if(binCenter <= 1000.) continue;
    
    clones_p.at(divPos)->SetBinContent(bIX+1, 0.0);
  }

  const std::string maxValStr = "Max: " + prettyString(maxVal, 3, false); 
  const std::string minValStr = "Min: " + prettyString(minVal, 3, false);
 
  Double_t interval = maxVal - minVal;
  maxVal += interval/10.;
  minVal = TMath::Max(0., minVal - interval/10.);
  
  if(maxVal < 1.0001 && minVal > .9999){
    maxVal = 1.0001;
    minVal = 0.9999;
  } 

  tempHist_p->SetMaximum(maxVal);
  tempHist_p->SetMinimum(minVal);
  
  for(Int_t bI = 0; bI < nClones; ++bI){
    clones_p.at(bI)->DrawCopy("HIST E1 P SAME");
  }
  
  const Int_t nTempBins = 20;
  Double_t tempBins[nTempBins+1];

  getLinBins(minVal, maxVal, nTempBins, tempBins);
  label_p->DrawLatex(250, tempBins[19], labelStr.c_str());
  label_p->DrawLatex(250, tempBins[18], maxValStr.c_str());
  label_p->DrawLatex(250, tempBins[17], minValStr.c_str());

  gStyle->SetOptStat(0);
  gPad->SetLogx();
	  
  saveName = "pdfDir/" + dateStr + "/PlotUnfold/" + saveName;
  quietSaveAs(canv_p, saveName);
  
  delete tempHist_p;
  delete canv_p;
  delete label_p;

  return;
}


int plotUnfoldedSpectra(const std::string inFileNamePP, const std::string inFileNamePbPb)
{
  std::cout << __LINE__ << std::endl;
  std::vector<std::string> slideTitles;
  std::vector<std::vector<std::string > > pdfPerSlide;
  std::vector<std::string> slideTitlesMain;
  std::vector<std::vector<std::string > > pdfPerSlideMain;

  double total = 0;
  cppWatch toPPFile;
  cppWatch ppFile;
  cppWatch pbpbFile;
  cppWatch comparisonOfFiles;

  toPPFile.start();

  std::cout << "Initializing files..." << std::endl;

  if(!checkFile(inFileNamePP) || !checkFile(inFileNamePbPb)){
    if(!checkFile(inFileNamePP)) std::cout << "inFileNamePP, \'" << inFileNamePP << "\', is invalid. return 1" << std::endl;
    if(!checkFile(inFileNamePbPb)) std::cout << "inFileNamePbPb, \'" << inFileNamePbPb << "\', is invalid. return 1" << std::endl;
    return 1;
  }

  kirchnerPalette kPalette;

  TDatime* date = new TDatime();
  const std::string dateStr = std::to_string(date->GetDate());
  delete date;

  const std::string dirStr = "pdfDir/" + dateStr + "/PlotUnfold";
  checkMakeDir("pdfDir");
  checkMakeDir("pdfDir/" + dateStr);
  checkMakeDir(dirStr);

  toPPFile.stop();
  total += toPPFile.totalWall();
  std::cout << "To PP File time: " << toPPFile.totalWall() << "/" << total << std::endl;
  ppFile.start();

  TFile* inFilePP_p = new TFile(inFileNamePP.c_str(), "READ");
  cutPropagator cutPropPP;
  cutPropPP.Clean();
  cutPropPP.GetAllVarFromFile(inFilePP_p);

  std::vector<std::string> histTagPP = cutPropPP.GetHistTag();
  std::vector<int> histBestBayesPP = cutPropPP.GetHistBestBayes();

  std::map<std::string, int> histTagMapPP;
  for(unsigned int i = 0; i < histTagPP.size(); ++i){
    histTagMapPP[histTagPP.at(i)] = histBestBayesPP.at(i);
  }

  std::vector<std::string> jetPPList = returnRootFileContentsList(inFilePP_p, "TDirectoryFile", "JetAnalyzer", 1);

  std::cout << "JetPPList: " << std::endl;
  for(unsigned int jI = 0; jI < jetPPList.size(); ++jI){
    std::cout << " " << jI << "/" << jetPPList.size() << ": " << jetPPList.at(jI) << std::endl;
  }

  ppFile.stop();
  total += ppFile.totalWall();
  std::cout << "PP File time: " << ppFile.totalWall() << "/" << total << std::endl;
  pbpbFile.start();

  TFile* inFilePbPb_p = new TFile(inFileNamePbPb.c_str(), "READ");
  cutPropagator cutPropPbPb;
  cutPropPbPb.Clean();
  cutPropPbPb.GetAllVarFromFile(inFilePbPb_p);
  std::vector<std::string> jetPbPbList = returnRootFileContentsList(inFilePbPb_p, "TDirectoryFile", "JetAnalyzer", 1);

  std::vector<std::string> histTagPbPb = cutPropPbPb.GetHistTag();
  std::vector<int> histBestBayesPbPb = cutPropPbPb.GetHistBestBayes();
  std::map<std::string, int> histTagMapPbPb;
  for(unsigned int i = 0; i < histTagPbPb.size(); ++i){
    histTagMapPbPb[histTagPbPb.at(i)] = histBestBayesPbPb.at(i);
  }

  std::cout << "JetPbPbList: " << std::endl;
  for(unsigned int jI = 0; jI < jetPbPbList.size(); ++jI){
    std::cout << " " << jI << "/" << jetPbPbList.size() << ": " << jetPbPbList.at(jI) << std::endl;
  }

  pbpbFile.stop();
  total += pbpbFile.totalWall();
  std::cout << "PBPB File time: " << pbpbFile.totalWall() << "/" << total << std::endl;
  comparisonOfFiles.start();

  if(!cutPropPP.GetIsPP() || cutPropPbPb.GetIsPP() || !cutPropPP.CheckPropagatorsMatch(cutPropPbPb, true, false)){
    if(!cutPropPP.GetIsPP()) std::cout << "inFileNamePP \'" << inFileNamePP << "\' is not pp, return 1" << std::endl;
    if(cutPropPbPb.GetIsPP()) std::cout << "inFileNamePbPb \'" << inFileNamePbPb << "\' is not PbPb, return 1" << std::endl;
    if(!cutPropPP.CheckPropagatorsMatch(cutPropPbPb, true, false)) std::cout << "Cut propagators of pp and pbpb do not match. inspect. return 1" << std::endl;
  
    inFilePP_p->Close();
    delete inFilePP_p;
    
    inFilePbPb_p->Close();
    delete inFilePbPb_p;
    
    return 1;
  }

  const Int_t nMaxJtPtBins = 50;
  const Int_t nGenJtPtSmallBinsSmallRCent0to10 = cutPropPbPb.GetNGenJtPtSmallBinsSmallRCent0to10();
  const Int_t nGenJtPtLargeBinsSmallRCent0to10 = cutPropPbPb.GetNGenJtPtLargeBinsSmallRCent0to10();
  const Int_t nGenJtPtSmallBinsLargeRCent0to10 = cutPropPbPb.GetNGenJtPtSmallBinsLargeRCent0to10();
  const Int_t nGenJtPtLargeBinsLargeRCent0to10 = cutPropPbPb.GetNGenJtPtLargeBinsLargeRCent0to10();
  const Int_t nRecoJtPtBinsSmallRCent0to10 = cutPropPbPb.GetNRecoJtPtBinsSmallRCent0to10();
  const Int_t nRecoJtPtBinsLargeRCent0to10 = cutPropPbPb.GetNRecoJtPtBinsLargeRCent0to10();

  const Int_t nGenJtPtSmallBinsSmallRCent10to30 = cutPropPbPb.GetNGenJtPtSmallBinsSmallRCent10to30();
  const Int_t nGenJtPtLargeBinsSmallRCent10to30 = cutPropPbPb.GetNGenJtPtLargeBinsSmallRCent10to30();
  const Int_t nGenJtPtSmallBinsLargeRCent10to30 = cutPropPbPb.GetNGenJtPtSmallBinsLargeRCent10to30();
  const Int_t nGenJtPtLargeBinsLargeRCent10to30 = cutPropPbPb.GetNGenJtPtLargeBinsLargeRCent10to30();
  const Int_t nRecoJtPtBinsSmallRCent10to30 = cutPropPbPb.GetNRecoJtPtBinsSmallRCent10to30();
  const Int_t nRecoJtPtBinsLargeRCent10to30 = cutPropPbPb.GetNRecoJtPtBinsLargeRCent0to10();

  const Int_t nGenJtPtSmallBinsSmallRCent30to50 = cutPropPbPb.GetNGenJtPtSmallBinsSmallRCent30to50();
  const Int_t nGenJtPtLargeBinsSmallRCent30to50 = cutPropPbPb.GetNGenJtPtLargeBinsSmallRCent30to50();
  const Int_t nGenJtPtSmallBinsLargeRCent30to50 = cutPropPbPb.GetNGenJtPtSmallBinsLargeRCent30to50();
  const Int_t nGenJtPtLargeBinsLargeRCent30to50 = cutPropPbPb.GetNGenJtPtLargeBinsLargeRCent30to50();
  const Int_t nRecoJtPtBinsSmallRCent30to50 = cutPropPbPb.GetNRecoJtPtBinsSmallRCent30to50();
  const Int_t nRecoJtPtBinsLargeRCent30to50 = cutPropPbPb.GetNRecoJtPtBinsLargeRCent30to50();

  const Int_t nGenJtPtSmallBinsSmallRCent50to90 = cutPropPbPb.GetNGenJtPtSmallBinsSmallRCent50to90();
  const Int_t nGenJtPtLargeBinsSmallRCent50to90 = cutPropPbPb.GetNGenJtPtLargeBinsSmallRCent50to90();
  const Int_t nGenJtPtSmallBinsLargeRCent50to90 = cutPropPbPb.GetNGenJtPtSmallBinsLargeRCent50to90();
  const Int_t nGenJtPtLargeBinsLargeRCent50to90 = cutPropPbPb.GetNGenJtPtLargeBinsLargeRCent50to90();
  const Int_t nRecoJtPtBinsSmallRCent50to90 = cutPropPbPb.GetNRecoJtPtBinsSmallRCent50to90();
  const Int_t nRecoJtPtBinsLargeRCent50to90 = cutPropPbPb.GetNRecoJtPtBinsLargeRCent50to90();

  std::vector<Double_t> genJtPtSmallBinsSmallRTemp = cutPropPbPb.GetGenJtPtSmallBinsSmallR();
  std::vector<Double_t> genJtPtLargeBinsSmallRTemp = cutPropPbPb.GetGenJtPtLargeBinsSmallR();
  std::vector<Double_t> genJtPtSmallBinsLargeRTemp = cutPropPbPb.GetGenJtPtSmallBinsLargeR();
  std::vector<Double_t> genJtPtLargeBinsLargeRTemp = cutPropPbPb.GetGenJtPtLargeBinsLargeR();
  std::vector<Double_t> recoJtPtBinsSmallRTemp = cutPropPbPb.GetRecoJtPtBinsSmallR();
  std::vector<Double_t> recoJtPtBinsLargeRTemp = cutPropPbPb.GetRecoJtPtBinsLargeR();

  smallOrLargeR rReader;
  if(!rReader.CheckNGenJtPtSmallBinsSmallRCent0to10(nGenJtPtSmallBinsSmallRCent0to10)) return 1;
  if(!rReader.CheckNGenJtPtSmallBinsLargeRCent0to10(nGenJtPtSmallBinsLargeRCent0to10)) return 1;
  if(!rReader.CheckNGenJtPtLargeBinsSmallRCent0to10(nGenJtPtLargeBinsSmallRCent0to10)) return 1;
  if(!rReader.CheckNGenJtPtLargeBinsLargeRCent0to10(nGenJtPtLargeBinsLargeRCent0to10)) return 1;
  if(!rReader.CheckNRecoJtPtBinsSmallRCent0to10(nRecoJtPtBinsSmallRCent0to10)) return 1;
  if(!rReader.CheckNRecoJtPtBinsLargeRCent0to10(nRecoJtPtBinsLargeRCent0to10)) return 1;

  if(!rReader.CheckNGenJtPtSmallBinsSmallRCent10to30(nGenJtPtSmallBinsSmallRCent10to30)) return 1;
  if(!rReader.CheckNGenJtPtSmallBinsLargeRCent10to30(nGenJtPtSmallBinsLargeRCent10to30)) return 1;
  if(!rReader.CheckNGenJtPtLargeBinsSmallRCent10to30(nGenJtPtLargeBinsSmallRCent10to30)) return 1;
  if(!rReader.CheckNGenJtPtLargeBinsLargeRCent10to30(nGenJtPtLargeBinsLargeRCent10to30)) return 1;
  if(!rReader.CheckNRecoJtPtBinsSmallRCent10to30(nRecoJtPtBinsSmallRCent10to30)) return 1;
  if(!rReader.CheckNRecoJtPtBinsLargeRCent10to30(nRecoJtPtBinsLargeRCent10to30)) return 1;

  if(!rReader.CheckNGenJtPtSmallBinsSmallRCent30to50(nGenJtPtSmallBinsSmallRCent30to50)) return 1;
  if(!rReader.CheckNGenJtPtSmallBinsLargeRCent30to50(nGenJtPtSmallBinsLargeRCent30to50)) return 1;
  if(!rReader.CheckNGenJtPtLargeBinsSmallRCent30to50(nGenJtPtLargeBinsSmallRCent30to50)) return 1;
  if(!rReader.CheckNGenJtPtLargeBinsLargeRCent30to50(nGenJtPtLargeBinsLargeRCent30to50)) return 1;
  if(!rReader.CheckNRecoJtPtBinsSmallRCent30to50(nRecoJtPtBinsSmallRCent30to50)) return 1;
  if(!rReader.CheckNRecoJtPtBinsLargeRCent30to50(nRecoJtPtBinsLargeRCent30to50)) return 1;

  if(!rReader.CheckNGenJtPtSmallBinsSmallRCent50to90(nGenJtPtSmallBinsSmallRCent50to90)) return 1;
  if(!rReader.CheckNGenJtPtSmallBinsLargeRCent50to90(nGenJtPtSmallBinsLargeRCent50to90)) return 1;
  if(!rReader.CheckNGenJtPtLargeBinsSmallRCent50to90(nGenJtPtLargeBinsSmallRCent50to90)) return 1;
  if(!rReader.CheckNGenJtPtLargeBinsLargeRCent50to90(nGenJtPtLargeBinsLargeRCent50to90)) return 1;
  if(!rReader.CheckNRecoJtPtBinsSmallRCent50to90(nRecoJtPtBinsSmallRCent50to90)) return 1;
  if(!rReader.CheckNRecoJtPtBinsLargeRCent50to90(nRecoJtPtBinsLargeRCent50to90)) return 1;

  if(!rReader.CheckGenJtPtSmallBinsSmallR(genJtPtSmallBinsSmallRTemp)) return 1;
  if(!rReader.CheckGenJtPtSmallBinsLargeR(genJtPtSmallBinsLargeRTemp)) return 1;
  if(!rReader.CheckGenJtPtLargeBinsSmallR(genJtPtLargeBinsSmallRTemp)) return 1;
  if(!rReader.CheckGenJtPtLargeBinsLargeR(genJtPtLargeBinsLargeRTemp)) return 1;
  if(!rReader.CheckRecoJtPtBinsSmallR(recoJtPtBinsSmallRTemp)) return 1;
  if(!rReader.CheckRecoJtPtBinsLargeR(recoJtPtBinsLargeRTemp)) return 1;

  comparisonOfFiles.stop();
  total += comparisonOfFiles.totalWall();
  std::cout << "Comp File time: " << comparisonOfFiles.totalWall() << "/" << total << std::endl;

  Int_t minForForLoop = 10000000;
  if(doLocalDebug || doGlobalDebug) minForForLoop = 1;

  std::cout << "Extracting cuts..." << std::endl;

  const Int_t nAdditionalSyst = 4;
  const std::string additionalSyst[nAdditionalSyst] = {"LumiUp", "LumiDown", "TAAUp", "TAADown"};

  const Int_t nJtMax = 10;
  const Int_t nJtPP = jetPPList.size();
  const Int_t nJtPbPb = jetPbPbList.size();

  if(nJtPP > nJtMax){
    std::cout << "nJtPP \'" << nJtPP << "\' greater than nJtMax \'" << nJtMax << "\'. return 1" << std::endl;
    return 1;
  }

  if(nJtPbPb > nJtMax){
    std::cout << "nJtPbPb \'" << nJtPbPb << "\' greater than nJtMax \'" << nJtMax << "\'. return 1" << std::endl;
    return 1;
  }

  const Int_t nMaxSyst = 20;
  const Int_t nSystOrig = 1;//cutPropPbPb.GetNSyst() - 1; // Minus 1 because removing priorflat
  //  const Int_t nSyst = TMath::Min(minForForLoop, nSystOrig + nAdditionalSyst);
  const Int_t nSyst = nSystOrig + nAdditionalSyst; 
  std::vector<std::string> systStr = {""};//cutPropPbPb.GetSystStr();
  Int_t jerDataPos = -1;

  if(nSyst > nMaxSyst){
    std::cout << "nSyst \'" << nSyst << "\' greater than nMaxSyst \'" << nMaxSyst << "\'. return 1" << std::endl;
    return 1;
  }

  for(unsigned int sI = 0; sI < systStr.size(); ++sI){
    if(isStrSame(systStr.at(sI), "PriorFlat")){
      systStr.erase(systStr.begin()+sI);
      break;
    }
  }

  for(Int_t sI = 0; sI < nAdditionalSyst; ++sI){
    systStr.push_back(additionalSyst[sI]);
  }

  std::cout << "Processing full systematics: " << std::endl;
  
  for(Int_t sI = 0; sI < nSyst; ++sI){
    std::cout << " " << sI << "/" << nSyst << ": " << systStr.at(sI) << std::endl;
    if(isStrSame(systStr.at(sI), "JERData")) jerDataPos = sI - 1;
  }

  if(jerDataPos >= 0) std::cout << "JERDATAPOS: " << jerDataPos << std::endl;
  else std::cout << "WARNING: JERDATAPOS NOT FOUND" << std::endl;

  const Int_t nBayesCap = 50;
  const Int_t nBayes = TMath::Min(nBayesCap, cutPropPbPb.GetNBayes());
  const Int_t nBigBayesSymm = cutPropPbPb.GetNBigBayesSymm();
  std::vector<int> bayesVal = cutPropPbPb.GetBayesVal();
  const Int_t nSuperBayes = cutPropPbPb.GetNSuperBayes();

  if(nBayes > nBayesCap){
    std::cout << "nBayes \'" << nBayes << "\' greater than nBayesCap \'" << nBayesCap << "\'. return 1" << std::endl;
    return 1;
  }
  
  const Int_t nMaxResponseMod = 3;
  const Int_t nResponseMod = TMath::Min(minForForLoop, cutPropPbPb.GetNResponseMod());
  //  const Int_t nResponseMod = cutPropPbPb.GetNResponseMod();
  std::vector<double> responseMod = cutPropPbPb.GetResponseMod();

  if(nResponseMod > nMaxResponseMod){
    std::cout << "nResponseMod \'" << nResponseMod << "\' greater than nMaxResponseMod \'" << nMaxResponseMod << "\'. return 1" << std::endl;
    return 1;
  }
  
  /*
  const Int_t nResponseMod = 1;
  std::vector<double> responseMod = {0.10};
  */

  const Int_t nMaxCentBins = 4;
  const Int_t nCentBins = cutPropPbPb.GetNCentBins();

  if(nCentBins > nMaxCentBins){
    std::cout << "nCentBins \'" << nCentBins << "\' greater than nMaxCentBins \'" << nMaxCentBins << "\'. return 1" << std::endl;
    return 1;
  }

  std::vector<Int_t> centBinsLow = cutPropPbPb.GetCentBinsLow();
  std::vector<Int_t> centBinsHi = cutPropPbPb.GetCentBinsHi();
  std::vector<Double_t> centBinsScalingFact;
  for(Int_t cI = 0; cI < nCentBins; ++cI){
    centBinsScalingFact.push_back(TMath::Power(10., cI+1.));
  }

  //  const Int_t nJtPtBins = cutPropPbPb.GetNJtPtBins();
  //  std::vector<Double_t> jtPtBinsTemp = cutPropPbPb.GetJtPtBins();

  //  const Int_t nJtAbsEtaBins = TMath::Min(minForForLoop, cutPropPbPb.GetNJtAbsEtaBins());

  
  const Int_t nMaxJtAbsEtaBins = 6;
  const Int_t nJtAbsEtaBins = cutPropPbPb.GetNJtAbsEtaBins();

  if(nJtAbsEtaBins > nMaxJtAbsEtaBins){
    std::cout << "nJtAbsEtaBins \'" << nJtAbsEtaBins << "\' greater than nMaxJtAbsEtaBins \'" << nMaxJtAbsEtaBins << "\'. return 1" << std::endl;
    return 1;
  }

  std::vector<Double_t> jtAbsEtaBinsLow = cutPropPbPb.GetJtAbsEtaBinsLow();
  std::vector<Double_t> jtAbsEtaBinsHi = cutPropPbPb.GetJtAbsEtaBinsHi();
  
  /*
  const Int_t nJtAbsEtaBins = 1;
  std::vector<Double_t> jtAbsEtaBinsLow = {0.0};
  std::vector<Double_t> jtAbsEtaBinsHi = {2.0};
  */

  //  const Int_t nID = cutPropPbPb.GetNID();

  const Int_t nMaxID = 6;
  const Int_t nID = TMath::Min(minForForLoop, cutPropPbPb.GetNID());
  std::vector<std::string> idStr = cutPropPbPb.GetIdStr();

  if(nID > nMaxID){
    std::cout << "nID \'" << nID << "\' greater than nMaxID \'" << nMaxID << "\'. return 1" << std::endl;
    return 1;
  }
  
  /*
  const Int_t nID = 1;
  std::vector<std::string> idStr = {"LightMUAndCHID"};
  */

  std::vector<int> jetPPMatchedToPbPb;
  for(Int_t jI = 0; jI < nJtPbPb; ++jI){
    const Int_t rValPbPb = getRVal(jetPbPbList.at(jI));

    bool isMatched = false;
    for(Int_t jI2 = 0; jI2 < nJtPP; ++jI2){
      const Int_t rValPP = getRVal(jetPPList.at(jI2));
      
      if(rValPbPb == rValPP){
	jetPPMatchedToPbPb.push_back(jI2);
	isMatched = true;
	break;
      }
    }

    if(!isMatched){
      std::cout << "Warning: " << jetPbPbList.at(jI) << " has no match." << std::endl;
      std::cout << " Options: ";
      for(Int_t jI2 = 0; jI2 < nJtPP; ++jI2){
	std::cout << jetPPList.at(jI2) << ", ";
      }
      std::cout << std::endl;
    }

  }

  if(jetPPMatchedToPbPb.size() != jetPbPbList.size()){
    std::cout << "Size of pp algos, " << jetPPMatchedToPbPb.size() << ",  matched to jetPbPbList of different size, " << jetPbPbList.size() << ". return 1" << std::endl;
    
    inFilePP_p->Close();
    delete inFilePP_p;
    
    inFilePbPb_p->Close();
    delete inFilePbPb_p;    

    return 1;
  }

  std::cout << "Processing the following pairs of jets: " << std::endl;
  for(Int_t jI = 0; jI < nJtPbPb; ++jI){
    std::cout << " " << jI << "/" << nJtPbPb << ": " << jetPbPbList.at(jI) << ", " << jetPPList.at(jetPPMatchedToPbPb.at(jI)) << std::endl;
  }

  if(doLocalDebug || doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  if(doLocalDebug || doGlobalDebug){
    std::cout << "nJtPbPb*nCentBins*nID*nResponseMod*nJtAbsEtaBins*nSyst*nBayes=Total" << std::endl;
    Int_t Total = nJtPbPb*nCentBins*nID*nResponseMod*nJtAbsEtaBins*nSyst*nBayes;
    std::cout << nJtPbPb << "*" << nCentBins << "*" << nID << "*" << nResponseMod << "*" << nJtAbsEtaBins << "*" << nSyst << "*" << nBayes << "=" << Total << std::endl;
  }


  std::string outFileName = inFileNamePP;
  while(outFileName.find("/") != std::string::npos){outFileName.replace(0, outFileName.find("/")+1, "");}
  if(outFileName.find(".txt") != std::string::npos) outFileName.replace(outFileName.find(".txt"), std::string(".txt").size(), "");
  else if(outFileName.find(".root") != std::string::npos) outFileName.replace(outFileName.find(".root"), std::string(".root").size(), "");

  std::string outFileName2 = inFileNamePbPb;
  while(outFileName2.find("/") != std::string::npos){outFileName2.replace(0, outFileName2.find("/")+1, "");}
  if(outFileName2.find(".txt") != std::string::npos) outFileName2.replace(outFileName2.find(".txt"), std::string(".txt").size(), "");
  else if(outFileName2.find(".root") != std::string::npos) outFileName2.replace(outFileName2.find(".root"), std::string(".root").size(), "");

  const Int_t sizeToTruncName = 40;
  while(outFileName.size() > sizeToTruncName){outFileName = outFileName.substr(0,outFileName.size()-1);}
  while(outFileName2.size() > sizeToTruncName){outFileName2 = outFileName2.substr(0,outFileName2.size()-1);}

  std::string tagStr = "";
  if(nJtPbPb == 1){
    tagStr = jetPbPbList.at(0) + "_" + jetPPList.at(jetPPMatchedToPbPb.at(0)) + "_";
  }

  outFileName = "output/plotUnfoldedSpectra_NSuperBayes" + std::to_string(nSuperBayes) + "_" + tagStr + dateStr + ".root";

  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");
  //Following two lines necessary else you get absolutely clobbered on deletion timing. See threads here:
  //https://root-forum.cern.ch/t/closing-root-files-is-dead-slow-is-this-a-normal-thing/5273/16
  //https://root-forum.cern.ch/t/tfile-speed/17549/25
  //Bizarre
  outFile_p->SetBit(TFile::kDevNull);
  TH1::AddDirectory(kFALSE);
 
  TDirectory* dirPP_p[nJtMax];
  TDirectory* dirPbPb_p[nJtMax];
  
  for(Int_t tI = 0; tI < nJtPbPb; ++tI){
    std::string tempStr = jetPbPbList.at(tI);
    if(tempStr.find("/") != std::string::npos) tempStr.replace(tempStr.find("/"), tempStr.size() - tempStr.find("/"), "");
    dirPbPb_p[tI] = (TDirectory*)outFile_p->mkdir(tempStr.c_str());
    dirPbPb_p[tI]->cd();
  }

  TH1D* jtPtUnfolded_RecoGenAsymm_PbPb_h[nJtMax][nMaxCentBins][nMaxID][nMaxResponseMod][nMaxJtAbsEtaBins][nMaxSyst];
  Int_t bayesPosPbPb[nJtMax][nMaxCentBins][nMaxID][nMaxResponseMod][nMaxJtAbsEtaBins][nMaxSyst];
  
  for(Int_t tI = 0; tI < nJtPP; ++tI){
    std::string tempStr = jetPPList.at(tI);
    if(tempStr.find("/") != std::string::npos) tempStr.replace(tempStr.find("/"), tempStr.size() - tempStr.find("/"), "");
    dirPP_p[tI] = (TDirectory*)outFile_p->mkdir(tempStr.c_str());
    dirPP_p[tI]->cd();
  }

  TH1D* jtPtUnfolded_RecoGenAsymm_PP_h[nJtMax][nMaxID][nMaxResponseMod][nMaxJtAbsEtaBins][nMaxSyst];
  Int_t bayesPosPP[nJtMax][nMaxID][nMaxResponseMod][nMaxJtAbsEtaBins][nMaxSyst];

  if(doLocalDebug || doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  std::cout << "Grabbing histograms..." << std::endl;

  if(doLocalDebug || doGlobalDebug) std::cout << "idI,mI,aI,sI,bI,tI,cI" << std::endl;

  cppWatch histGrabTot;
  cppWatch histGrab;
  cppWatch histClone;

  histGrabTot.start();

  inFilePbPb_p->cd();
  for(unsigned int tI = 0; tI < jetPbPbList.size(); ++tI){
    cppWatch histGrabLocal;
    histGrabLocal.start();
    histGrab.start();
    std::cout << " Grabbing PbPb " << tI << "/" << jetPbPbList.size() << ": " << jetPbPbList.at(tI) << std::endl;

    TDirectory* dir_p = (TDirectory*)inFilePbPb_p->Get(jetPbPbList.at(tI).c_str());
    TIter nextPbPb(dir_p->GetListOfKeys());
    TKey* key = NULL;

    while( (key = (TKey*)nextPbPb()) ){
      const std::string name = key->GetName();
      const std::string className = key->GetClassName();

      if(className.find("TH1") == std::string::npos) continue;

      int idPos = -1;
      int modPos = -1;
      int absEtaPos = -1;
      int systPos = -1;
      int bayesPos = -1;
      int centPos = -1;

      for(unsigned int idI = 0; idI < idStr.size(); ++idI){
	if(name.find("_" + idStr.at(idI) + "_") != std::string::npos){
	  idPos = idI;
	  break;
	}
      }

      for(int mI = 0; mI < nResponseMod; ++mI){
	if(name.find("_ResponseMod" + prettyString(responseMod.at(mI), 2, true) + "_") != std::string::npos){
	  modPos = mI;
	  break;
	}
      }

      for(int aI = 0; aI < nJtAbsEtaBins; ++aI){
	const std::string jtAbsEtaStr = "AbsEta" + prettyString(jtAbsEtaBinsLow.at(aI), 1, true) + "to" + prettyString(jtAbsEtaBinsHi.at(aI), 1, true);	
	if(name.find(jtAbsEtaStr) != std::string::npos){
	  absEtaPos = aI;
	  break;
	}
      }

      for(Int_t sI = 1; sI < nSystOrig; ++sI){
	if(name.find(systStr.at(sI)) != std::string::npos){
	  systPos = sI;
	  break;
	}
      }
      if(systPos < 0 && name.find("PriorFlat") == std::string::npos) systPos = 0;

      for(Int_t bI = 0; bI < nBayes; ++bI){
	const std::string bayesStr = "_Bayes" + std::to_string(bayesVal.at(bI)) + "_";
	if(name.find(bayesStr) != std::string::npos){
	  bayesPos = bI;
	  break;
	}
      }

      for(Int_t cI = 0; cI < nCentBins; ++cI){		
	const std::string centStr = "Cent" + std::to_string(centBinsLow.at(cI)) + "to" + std::to_string(centBinsHi.at(cI));
	if(name.find(centStr) != std::string::npos){
	  centPos = cI;
	  break;
	}
      }


      if(idPos < 0) std::cout << "Warning: Cannot find idPos for \'" << name << "\'. Code will likely break" << std::endl;
      if(modPos < 0) std::cout << "Warning: Cannot find modPos for \'" << name << "\'. Code will likely break" << std::endl;
      if(absEtaPos < 0) std::cout << "Warning: Cannot find absEtaPos for \'" << name << "\'. Code will likely break" << std::endl;
      if(systPos < 0){
	if(name.find("PriorFlat") == std::string::npos) std::cout << "Warning: Cannot find systPos for \'" << name << "\'. Code will likely break" << std::endl;
	else continue;
      }
      if(bayesPos < 0) std::cout << "Warning: Cannot find bayesPos for \'" << name << "\'. Code will likely break" << std::endl;
      if(centPos < 0) std::cout << "Warning: Cannot find centPos for \'" << name << "\'. Code will likely break" << std::endl;

      if(bayesPos != 3) continue;

      jtPtUnfolded_RecoGenAsymm_PbPb_h[tI][centPos][idPos][modPos][absEtaPos][systPos] = (TH1D*)key->ReadObj();   
    }

    //Building best bayes
    for(Int_t cI = 0; cI < nCentBins; ++cI){
      const std::string centStr = "Cent" + std::to_string(centBinsLow.at(cI)) + "to" + std::to_string(centBinsHi.at(cI));
      
      for(unsigned int idI = 0; idI < idStr.size(); ++idI){
	for(Int_t mI = 0; mI < nResponseMod; ++mI){
	  const std::string resStr = "ResponseMod" + prettyString(responseMod.at(mI), 2, true);
	  
	  for(Int_t aI = 0; aI < nJtAbsEtaBins; ++aI){
	    const std::string jtAbsEtaStr = "AbsEta" + prettyString(jtAbsEtaBinsLow.at(aI), 1, true) + "to" + prettyString(jtAbsEtaBinsHi.at(aI), 1, true);	
	    
	    for(Int_t sI = 0; sI < nSystOrig; ++sI){
	      const std::string tempHistTag = jetPbPbList.at(tI) + "_" + centStr + "_" + idStr.at(idI) + "_" + resStr + "_" + jtAbsEtaStr + "_" + systStr.at(sI);
	      if(histTagMapPbPb.count(tempHistTag) == 0){
		std::cout << "WARNING: tag \'" << tempHistTag << "\' not found. Setting to -1" << std::endl;
		bayesPosPbPb[tI][cI][idI][mI][aI][sI] = -1;
	      }
	      else{
		int tempSet = histTagMapPbPb[tempHistTag];
		std::cout << "tag \'" << tempHistTag << "\' found. Setting to " << tempSet << std::endl;
		bayesPosPbPb[tI][cI][idI][mI][aI][sI] = tempSet;
	      }	      
	    }
	  }
	}
      }
    }
  
    histGrabLocal.stop();
    histGrab.stop();

    cppWatch histCloneLocal;
    histCloneLocal.start();
    histClone.start();
    
    std::cout << "  Cloning..." << std::endl;

    for(Int_t cI = 0; cI < nCentBins; ++cI){
      const std::string centStr = "Cent" + std::to_string(centBinsLow.at(cI)) + "to" + std::to_string(centBinsHi.at(cI));

      for(unsigned int idI = 0; idI < idStr.size(); ++idI){
	for(Int_t mI = 0; mI < nResponseMod; ++mI){
	  const std::string resStr = "ResponseMod" + prettyString(responseMod.at(mI), 2, true);
	  
	  for(Int_t aI = 0; aI < nJtAbsEtaBins; ++aI){
	    const std::string jtAbsEtaStr = "AbsEta" + prettyString(jtAbsEtaBinsLow.at(aI), 1, true) + "to" + prettyString(jtAbsEtaBinsHi.at(aI), 1, true);	
	    
	    for(Int_t sI = nSystOrig; sI < nSyst; ++sI){
	      std::string tempSystStr = "_" + systStr.at(sI) + "_";
	      while(tempSystStr.find("__") != std::string::npos){tempSystStr.replace(tempSystStr.find("__"), 2, "_");}

	      bayesPosPbPb[tI][cI][idI][mI][aI][sI] = bayesPosPbPb[tI][cI][idI][mI][aI][0];

	      for(Int_t bI = 0; bI < nBayes; ++bI){
		const std::string bayesStr = "Bayes" + std::to_string(bayesVal.at(bI));
		
		const std::string newName = "jtPtUnfolded_RecoGenAsymm_" + jetPbPbList.at(tI) + "_" + centStr + "_" + idStr.at(idI) + "_" + resStr + "_" + jtAbsEtaStr + tempSystStr + bayesStr + "_h";
		
		std::cout << " Name to grab: " << newName << std::endl;
		jtPtUnfolded_RecoGenAsymm_PbPb_h[tI][cI][idI][mI][aI][sI] = (TH1D*)jtPtUnfolded_RecoGenAsymm_PbPb_h[tI][cI][idI][mI][aI][0]->Clone(newName.c_str());
		std::cout << "Grab successful"  << std::endl;
	      }
	    }
	  }
	}
      }
    }

    histCloneLocal.stop();
    histClone.stop();

    std::cout << "  Local histgrab, clone times: "<< histGrabLocal.totalWall() << ", " << histCloneLocal.totalWall() << std::endl;
  }

  

  inFilePP_p->cd();
  for(unsigned int tI = 0; tI < jetPPList.size(); ++tI){
    cppWatch histGrabLocal;
    histGrabLocal.start();
    histGrab.start();
    std::cout << " Grabbing PP " << tI << "/" << jetPPList.size() << ": " << jetPPList.at(tI) << std::endl;

    TDirectory* dir_p = (TDirectory*)inFilePP_p->Get(jetPPList.at(tI).c_str());
    TIter nextPP(dir_p->GetListOfKeys());
    TKey* key = NULL;

    while( (key = (TKey*)nextPP()) ){
      const std::string name = key->GetName();
      const std::string className = key->GetClassName();

      if(className.find("TH1") == std::string::npos) continue;

      int idPos = -1;
      int modPos = -1;
      int absEtaPos = -1;
      int systPos = -1;
      int bayesPos = -1;

      for(unsigned int idI = 0; idI < idStr.size(); ++idI){
	if(name.find("_" + idStr.at(idI) + "_") != std::string::npos){
	  idPos = idI;
	  break;
	}
      }

      for(int mI = 0; mI < nResponseMod; ++mI){
	if(name.find("_ResponseMod" + prettyString(responseMod.at(mI), 2, true) + "_") != std::string::npos){
	  modPos = mI;
	  break;
	}
      }

      for(int aI = 0; aI < nJtAbsEtaBins; ++aI){
	const std::string jtAbsEtaStr = "AbsEta" + prettyString(jtAbsEtaBinsLow.at(aI), 1, true) + "to" + prettyString(jtAbsEtaBinsHi.at(aI), 1, true);	
	if(name.find(jtAbsEtaStr) != std::string::npos){
	  absEtaPos = aI;
	  break;
	}
      }

      for(Int_t sI = 1; sI < nSystOrig; ++sI){
	if(name.find(systStr.at(sI)) != std::string::npos){
	  systPos = sI;
	  break;
	}
      }
      if(systPos < 0 && name.find("PriorFlat") == std::string::npos) systPos = 0;

      for(Int_t bI = 0; bI < nBayes; ++bI){
	const std::string bayesStr = "_Bayes" + std::to_string(bayesVal.at(bI)) + "_";
	if(name.find(bayesStr) != std::string::npos){
	  bayesPos = bI;
	  break;
	}
      }

      if(idPos < 0) std::cout << "Warning: Cannot find idPos for \'" << name << "\'. Code will likely break" << std::endl;
      if(modPos < 0) std::cout << "Warning: Cannot find modPos for \'" << name << "\'. Code will likely break" << std::endl;
      if(absEtaPos < 0) std::cout << "Warning: Cannot find absEtaPos for \'" << name << "\'. Code will likely break" << std::endl;
      if(systPos < 0){
	if(name.find("PriorFlat") == std::string::npos) std::cout << "Warning: Cannot find systPos for \'" << name << "\'. Code will likely break" << std::endl;
	else continue;
      }
      if(bayesPos < 0) std::cout << "Warning: Cannot find bayesPos for \'" << name << "\'. Code will likely break" << std::endl;

      if(bayesPos != 3) continue;

      jtPtUnfolded_RecoGenAsymm_PP_h[tI][idPos][modPos][absEtaPos][systPos] = (TH1D*)key->ReadObj();   
    }

    //Building best bayes
    for(unsigned int idI = 0; idI < idStr.size(); ++idI){
      for(Int_t mI = 0; mI < nResponseMod; ++mI){
	const std::string resStr = "ResponseMod" + prettyString(responseMod.at(mI), 2, true);
	
	for(Int_t aI = 0; aI < nJtAbsEtaBins; ++aI){
	  const std::string jtAbsEtaStr = "AbsEta" + prettyString(jtAbsEtaBinsLow.at(aI), 1, true) + "to" + prettyString(jtAbsEtaBinsHi.at(aI), 1, true);	
	  
	  for(Int_t sI = 0; sI < nSystOrig; ++sI){
	    const std::string tempHistTag = jetPPList.at(tI) + "_PP_" + idStr.at(idI) + "_" + resStr + "_" + jtAbsEtaStr + "_" + systStr.at(sI);
	    if(histTagMapPP.count(tempHistTag) == 0){
	      std::cout << "WARNING: tag \'" << tempHistTag << "\' not found. Setting to -1" << std::endl;
	      bayesPosPP[tI][idI][mI][aI][sI] = -1;
	    }
	    else{
	      int tempSet = histTagMapPP[tempHistTag];
	      std::cout << "tag \'" << tempHistTag << "\' found. Setting to " << tempSet << std::endl;
	      bayesPosPP[tI][idI][mI][aI][sI] = tempSet;
	    }
	  }
	}
      }
    }

    histGrabLocal.stop();
    histGrab.stop();
    
    cppWatch histCloneLocal;
    histCloneLocal.start();
    histClone.start();

    std::cout << "  Cloning..." << std::endl;

    for(unsigned int idI = 0; idI < idStr.size(); ++idI){
      for(Int_t mI = 0; mI < nResponseMod; ++mI){
	const std::string resStr = "ResponseMod" + prettyString(responseMod.at(mI), 2, true);
	
	for(Int_t aI = 0; aI < nJtAbsEtaBins; ++aI){
	  const std::string jtAbsEtaStr = "AbsEta" + prettyString(jtAbsEtaBinsLow.at(aI), 1, true) + "to" + prettyString(jtAbsEtaBinsHi.at(aI), 1, true);	
	  
	  for(Int_t sI = nSystOrig; sI < nSyst; ++sI){
	    const std::string tempSystStr = systStr.at(sI) + "_";

	    bayesPosPP[tI][idI][mI][aI][sI] = bayesPosPP[tI][idI][mI][aI][0];
	    
	    for(Int_t bI = 0; bI < nBayes; ++bI){
	      const std::string bayesStr = "Bayes" + std::to_string(bayesVal.at(bI));
	      
	      const std::string newName = "jtPtUnfolded_RecoGenAsymm_" + jetPPList.at(tI) + "_PP_" + idStr.at(idI) + "_" + resStr + "_" + jtAbsEtaStr + "_" + tempSystStr + bayesStr + "_h";
	      
	      jtPtUnfolded_RecoGenAsymm_PP_h[tI][idI][mI][aI][sI] = (TH1D*)jtPtUnfolded_RecoGenAsymm_PP_h[tI][idI][mI][aI][0]->Clone(newName.c_str());
	    }
	  }
	}
      }
    }

    histCloneLocal.stop();
    histClone.stop();

    std::cout << "  Local histgrab, clone times: " << histGrabLocal.totalWall() << ", " << histCloneLocal.totalWall() << std::endl;
  }

  
  histGrabTot.stop();
  std::cout << "Total histGrabTot: " << histGrabTot.totalWall() << std::endl;
  std::cout << "  just grabbing: " << histGrab.totalWall() << std::endl;
  std::cout << "  just cloning: " << histClone.totalWall() << std::endl;


  if(doLocalDebug || doGlobalDebug){
    std::cout << "Started return" << std::endl;
    return 1;
  }


  const Double_t lumiFactor = getLumiFactor();
  const Double_t nMBEvents = getNMBEvents();

  std::cout << "Scaling histograms..." << std::endl;

  for(Int_t idI = 0; idI < nID; ++idI){
    for(Int_t mI = 0; mI < nResponseMod; ++mI){
      for(Int_t aI = 0; aI < nJtAbsEtaBins; ++aI){
	const Double_t etaBinWidth = TMath::Abs(jtAbsEtaBinsHi[aI] - jtAbsEtaBinsLow[aI]);

	for(Int_t sI = 0; sI < nSyst; ++sI){
	  for(Int_t bI = 0; bI < nBayes; ++bI){
	    for(Int_t tI = 0; tI < nJtPbPb; ++tI){
	      for(Int_t cI = 0; cI < nCentBins; ++cI){		
		const std::string centStr = "Cent" + std::to_string(centBinsLow.at(cI)) + "to" + std::to_string(centBinsHi.at(cI));
		const Double_t centBinWidth = TMath::Abs(centBinsHi[cI] - centBinsLow[cI])/100.;

		Double_t tempTAAFactor = getTAAScaleFactorMB(centStr);
		if(isStrSame(systStr[sI], "TAAUp")) tempTAAFactor += tempTAAFactor*getTAAScaleFactorUp(centStr);
		else if(isStrSame(systStr[sI], "TAADown")) tempTAAFactor -= tempTAAFactor*getTAAScaleFactorDown(centStr);

		Double_t totalPbPbFactor = tempTAAFactor*nMBEvents*2.*etaBinWidth*centBinWidth;

		jtPtUnfolded_RecoGenAsymm_PbPb_h[tI][cI][idI][mI][aI][sI]->Scale(1./totalPbPbFactor);
		divBinWidth(jtPtUnfolded_RecoGenAsymm_PbPb_h[tI][cI][idI][mI][aI][sI]);
	      }
	    }

	    for(Int_t tI = 0; tI < nJtPP; ++tI){
	      Double_t tempLumiFactor = lumiFactor;
	      if(isStrSame(systStr[sI], "LumiUp")) tempLumiFactor += tempLumiFactor*getLumiPercentError();
	      else if(isStrSame(systStr[sI], "LumiDown")) tempLumiFactor -= tempLumiFactor*getLumiPercentError();	      

	      Double_t totalPPFactor = tempLumiFactor*2.*etaBinWidth;
	      jtPtUnfolded_RecoGenAsymm_PP_h[tI][idI][mI][aI][sI]->Scale(1./totalPPFactor);
	      divBinWidth(jtPtUnfolded_RecoGenAsymm_PP_h[tI][idI][mI][aI][sI]);
	    }
	  }
	}
      }
    }
  }

  std::cout << "Writing histograms..." << std::endl;

  outFile_p->cd();
  for(Int_t tI = 0; tI < nJtPbPb; ++tI){
    dirPbPb_p[tI]->cd();

    for(Int_t cI = 0; cI < nCentBins; ++cI){		
      for(Int_t idI = 0; idI < nID; ++idI){
	for(Int_t mI = 0; mI < nResponseMod; ++mI){
	  for(Int_t aI = 0; aI < nJtAbsEtaBins; ++aI){
	    for(Int_t sI = 0; sI < nSyst; ++sI){
	      for(Int_t bI = 0; bI < nBayes; ++bI){
		std::string newName = jtPtUnfolded_RecoGenAsymm_PbPb_h[tI][cI][idI][mI][aI][sI]->GetName();
		newName.replace(newName.find("_h"), 2, "_Rescaled_h");
		jtPtUnfolded_RecoGenAsymm_PbPb_h[tI][cI][idI][mI][aI][sI]->Write(newName.c_str(), TObject::kOverwrite);
	      }
	    }
	  }
	}
      }
    }
  }

  for(Int_t tI = 0; tI < nJtPP; ++tI){
    dirPP_p[tI]->cd();
    
    for(Int_t idI = 0; idI < nID; ++idI){
      for(Int_t mI = 0; mI < nResponseMod; ++mI){
	for(Int_t aI = 0; aI < nJtAbsEtaBins; ++aI){
	  for(Int_t sI = 0; sI < nSyst; ++sI){
	    for(Int_t bI = 0; bI < nBayes; ++bI){
	      std::string newName = jtPtUnfolded_RecoGenAsymm_PP_h[tI][idI][mI][aI][sI]->GetName();
	      newName.replace(newName.find("_h"), 2, "_Rescaled_h");
	      jtPtUnfolded_RecoGenAsymm_PP_h[tI][idI][mI][aI][sI]->Write(newName.c_str(), TObject::kOverwrite);
	    }
	  }
	}
      }
    }
  }

  std::cout << "Plotting... " << std::endl;
  const std::string plotID = "LightMUAndCHID";
  const std::string plotAbsEtaStr = "AbsEta0p0to2p0";
  const Int_t plotBayesVal = 100;

  Int_t idPos = -1;
  Int_t absEtaPos = -1;
  Int_t bayesPos = -1;

  for(Int_t idI = 0; idI < nID; ++idI){
    if(isStrSame(idStr.at(idI), plotID)){
      idPos = idI;
      break;
    }
  }

  for(Int_t aI = 0; aI < nJtAbsEtaBins; ++aI){
    const std::string jtAbsEtaStr = "AbsEta" + prettyString(jtAbsEtaBinsLow.at(aI), 1, true) + "to" + prettyString(jtAbsEtaBinsHi.at(aI), 1, true);	
    
    if(isStrSame(jtAbsEtaStr, plotAbsEtaStr)){
      absEtaPos = aI;
      break;
    }
  }

  for(Int_t bI = 0; bI < nBayes; ++bI){
    if(bayesVal.at(bI) == plotBayesVal){
      bayesPos = bI;
      break;
    }
  }
  
  if(idPos == -1 || absEtaPos == -1 || bayesPos == -1){
    std::cout << "Plotting pos not found, skip plotting, return 1" << std::endl;

    std::cout << "Closing files..." << std::endl;
    
    outFile_p->Close();
    delete outFile_p;
  
    inFilePP_p->Close();
    delete inFilePP_p;

    inFilePbPb_p->Close();
    delete inFilePbPb_p;

    return 1;
  }

  std::vector<std::vector<std::string> > systToCombo = {{"JECUpMC", "JECDownMC", "JECMC"}, {"JECUpData", "JECDownData", "JECData"}, {"JECUpUE", "JECDownUE", "JECUE"}, {"PriorUp1PowerPthat", "PriorDown1PowerPthat", "Prior1PowerPthat"}, {"LumiUp", "LumiDown", "Lumi"}, {"TAAUp", "TAADown", "TAA"}};

  const std::string idNameStr = idStr.at(idPos);
  const std::string plotBayesStr = "Bayes" + std::to_string(plotBayesVal);

  const Int_t nSystSmooth = 2;
  const std::string systSmooth[nSystSmooth] = {"Unsmoothed", "Smoothed"};
  const bool systSmoothBool[nSystSmooth] = {false, true};

  //Check convergence bands...
  std::cout << "Do convergence bands" << std::endl;
  for(Int_t tI = 0; tI < nJtPbPb; ++tI){
    const Int_t rVal = getRVal(jetPbPbList.at(tI));
    const Int_t ppPos = jetPPMatchedToPbPb.at(tI);
    const std::string rValStr = std::to_string(rVal);

    const bool isSmallR = rReader.GetIsSmallR(rVal);

    std::string rStr = "R=";
    if(getRVal(jetPbPbList.at(tI)) < 10) rStr = rStr + "0." + std::to_string(getRVal(jetPbPbList.at(tI)));
    else rStr = rStr + "1.0";

    const std::string jetStrPbPb = jetPbPbList.at(tI).substr(0, jetPbPbList.at(tI).find("JetAna"));
    const std::string jetStrPP = jetPPList.at(ppPos).substr(0, jetPPList.at(ppPos).find("JetAna"));

    for(Int_t cI = 0; cI < nCentBins; ++ cI){
      const std::string centStr = "Cent" + std::to_string(centBinsLow.at(cI)) + "to" + std::to_string(centBinsHi.at(cI));
      const std::string centStr2 = std::to_string(centBinsLow.at(cI)) + "-" + std::to_string(centBinsHi.at(cI)) + "%";
      
      const Int_t nRecoJtPtBins = rReader.GetSmallOrLargeRNBins(isSmallR, false, false, centStr);
      Double_t recoJtPtBins[nMaxJtPtBins+1];
      rReader.GetSmallOrLargeRBins(isSmallR, false, nRecoJtPtBins+1, recoJtPtBins, false);
      
      for(Int_t mI = 0; mI < nResponseMod; ++mI){
	//start pp closures
	slideTitles.push_back("Conv. " + jetStrPP + " (ModRes " + prettyString(responseMod[mI],2,false) + ")");
	pdfPerSlide.push_back({});
	
	for(Int_t sI = 0; sI < nSyst; ++sI){	  
	  const Int_t totBigBayes = 1 + 2*nBigBayesSymm;
	  std::vector<TH1D*> bigBayesClones_p;
	  bigBayesClones_p.reserve(totBigBayes);
	  for(Int_t bI = 0; bI < totBigBayes; ++bI){
	    //	  Int_t binPos = nBayes - 1 -2*nBigBayesSymm + bI;
	    bigBayesClones_p.push_back((TH1D*)jtPtUnfolded_RecoGenAsymm_PP_h[ppPos][idPos][mI][absEtaPos][sI]->Clone(("bigBayesClone_" + std::to_string(bI)).c_str()));
	  }      
	  
	  std::string saveName = "convBayes_NBigBayesSymm" + std::to_string(nBigBayesSymm) + "_" + jetPPList.at(ppPos) + "_R" + rValStr + "_PP_" + centStr + "_" + idNameStr + "_ResponseMod" + prettyString(responseMod[mI], 2, true) + "_" + plotAbsEtaStr + "_" + systStr.at(sI) + "_" + dateStr + ".pdf";
	  pdfPerSlide.at(pdfPerSlide.size()-1).push_back(saveName);
	  
	  doRatioPlotting(bigBayesClones_p, nBigBayesSymm, "Jet p_{T} (GeV/c)", "Ratio (Bayes Iter 100 #pm" + std::to_string(nBigBayesSymm) + ")", recoJtPtBins[0], 1000., systStr.at(sI), saveName);
	  for(unsigned int ucI = 0; ucI < bigBayesClones_p.size(); ++ucI){
	    delete bigBayesClones_p.at(ucI);
	    bigBayesClones_p.at(ucI) = NULL;
	  }
	}
	
	slideTitles.push_back("BestComp. " + jetStrPP + " (ModRes " + prettyString(responseMod[mI],2,false) + ")");
	pdfPerSlide.push_back({});
	
	for(Int_t sI = 0; sI < nSyst; ++sI){	  
	  std::vector<TH1D*> clones_p;
	  
	  Int_t binPos = nBayes - 1 -nBigBayesSymm;	  
	  clones_p.push_back((TH1D*)jtPtUnfolded_RecoGenAsymm_PP_h[ppPos][idPos][mI][absEtaPos][sI]->Clone(("clone_" + std::to_string(binPos)).c_str()));
	  
	  if(bayesPosPP[ppPos][idPos][mI][absEtaPos][sI] < 0) bayesPosPP[ppPos][idPos][mI][absEtaPos][sI] = 25;
      
	  std::cout << jtPtUnfolded_RecoGenAsymm_PP_h[ppPos][idPos][mI][absEtaPos][sI]->GetName() << std::endl;
	  Int_t tempBayesPos = bayesPosPP[ppPos][idPos][mI][absEtaPos][sI];
	  tempBayesPos = 3;
	  
	  clones_p.push_back((TH1D*)jtPtUnfolded_RecoGenAsymm_PP_h[ppPos][idPos][mI][absEtaPos][sI]->Clone(("clone_" + std::to_string(tempBayesPos)).c_str()));
	
	  std::string saveName = "convBestBayes_NBigBayesSymm" + std::to_string(nBigBayesSymm) + "_" + jetPPList.at(ppPos) + "_R" + rValStr + "_PP_" + idNameStr + "_ResponseMod" + prettyString(responseMod[mI], 2, true) + "_" + plotAbsEtaStr + "_" + systStr.at(sI) + "_" + dateStr + ".pdf";
	  pdfPerSlide.at(pdfPerSlide.size()-1).push_back(saveName);
	  
	  std::string labelStr = systStr.at(sI) + ", Bayes" + std::to_string(tempBayesPos);
	  doRatioPlotting(clones_p, 0, "Jet p_{T} (GeV/c)", "Ratio (Bayes Iter 100 w/ " + std::to_string(tempBayesPos) + ")", recoJtPtBins[0], 1000., labelStr, saveName);
	
	  for(unsigned int ucI = 0; ucI < clones_p.size(); ++ucI){
	    delete clones_p.at(ucI);
	    clones_p.at(ucI) = NULL;
	  }      
	}
	//end pp closure ratios
	
	slideTitles.push_back("Conv. " + jetStrPbPb + " (ModRes " + prettyString(responseMod[mI],2,false) + ", " + centStr2 + ")");
	pdfPerSlide.push_back({});

	for(Int_t sI = 0; sI < nSyst; ++sI){	  
	  const Int_t totBigBayes = 1 + 2*nBigBayesSymm;
	  std::vector<TH1D*> bigBayesClones_p;
	  bigBayesClones_p.reserve(totBigBayes);
	  for(Int_t bI = 0; bI < totBigBayes; ++bI){
	    //	    Int_t binPos = nBayes - 1 -2*nBigBayesSymm + bI;
	    bigBayesClones_p.push_back((TH1D*)jtPtUnfolded_RecoGenAsymm_PbPb_h[tI][cI][idPos][mI][absEtaPos][sI]->Clone(("bigBayesClone_" + std::to_string(bI)).c_str()));
	  }      

	  std::string saveName = "convBayes_NBigBayesSymm" + std::to_string(nBigBayesSymm) + "_" + jetPbPbList.at(tI) + "_R" + rValStr + "_PbPb_" + centStr + "_" + idNameStr + "_ResponseMod" + prettyString(responseMod[mI], 2, true) + "_" + plotAbsEtaStr + "_" + systStr.at(sI) + "_" + dateStr + ".pdf";
	  pdfPerSlide.at(pdfPerSlide.size()-1).push_back(saveName);

	  doRatioPlotting(bigBayesClones_p, nBigBayesSymm, "Jet p_{T} (GeV/c)", "Ratio (Bayes Iter 100 #pm" + std::to_string(nBigBayesSymm) + ")", recoJtPtBins[0], 1000., systStr.at(sI), saveName);
	  for(unsigned int ucI = 0; ucI < bigBayesClones_p.size(); ++ucI){
	    delete bigBayesClones_p.at(ucI);
	    bigBayesClones_p.at(ucI) = NULL;
	  }
	}
      
	slideTitles.push_back("BestComp. " + jetStrPbPb + " (ModRes " + prettyString(responseMod[mI],2,false) + ", " + centStr2 + ")");
	pdfPerSlide.push_back({});

	for(Int_t sI = 0; sI < nSyst; ++sI){	  
	  std::vector<TH1D*> clones_p;

	  Int_t binPos = nBayes - 1 -nBigBayesSymm;	  
	  clones_p.push_back((TH1D*)jtPtUnfolded_RecoGenAsymm_PbPb_h[tI][cI][idPos][mI][absEtaPos][sI]->Clone(("clone_" + std::to_string(binPos)).c_str()));
	  
	  if(bayesPosPbPb[tI][cI][idPos][mI][absEtaPos][sI] < 0) bayesPosPbPb[tI][cI][idPos][mI][absEtaPos][sI] = 25;

	  Int_t tempBayesPos = bayesPosPbPb[tI][cI][idPos][mI][absEtaPos][sI];
	  tempBayesPos = 3;

	  clones_p.push_back((TH1D*)jtPtUnfolded_RecoGenAsymm_PbPb_h[tI][cI][idPos][mI][absEtaPos][sI]->Clone(("clone_" + std::to_string(tempBayesPos)).c_str()));

	  std::string saveName = "convBestBayes_NBigBayesSymm" + std::to_string(nBigBayesSymm) + "_" + jetPbPbList.at(tI) + "_R" + rValStr + "_PbPb_" + centStr + "_" + idNameStr + "_ResponseMod" + prettyString(responseMod[mI], 2, true) + "_" + plotAbsEtaStr + "_" + systStr.at(sI) + "_" + dateStr + ".pdf";
	  pdfPerSlide.at(pdfPerSlide.size()-1).push_back(saveName);

	  std::string labelStr = systStr.at(sI) + ", Bayes" + std::to_string(tempBayesPos);
	  doRatioPlotting(clones_p, 0, "Jet p_{T} (GeV/c)", "Ratio (Bayes Iter 100 w/ " + std::to_string(tempBayesPos) + ")", recoJtPtBins[0], 1000., labelStr, saveName);

	  for(unsigned int ucI = 0; ucI < clones_p.size(); ++ucI){
	    delete clones_p.at(ucI);
	    clones_p.at(ucI) = NULL;
	  }

	}
      }
    }
  }

  //Actual spectra plots...
  for(Int_t tI = 0; tI < nJtPbPb; ++tI){
    const Int_t rVal = getRVal(jetPbPbList.at(tI));
    const Int_t ppPos = jetPPMatchedToPbPb.at(tI);
    const std::string rValStr = std::to_string(rVal);
    
    std::string rStr = "R=";
    if(getRVal(jetPbPbList.at(tI)) < 10) rStr = rStr + "0." + std::to_string(getRVal(jetPbPbList.at(tI)));
    else rStr = rStr + "1.0";
    
    Double_t xMinVal = 200;
    Double_t xPointMinVal = xMinVal;
    if(rVal >= 8) xPointMinVal = 300;
    Double_t xMaxVal = 1200;
    Double_t xPointMaxVal = 1000;
    Double_t xPointMaxVal5090 = 600;
    
    TLegend* legRes_p = new TLegend(0.7, 0.65, 0.95, 0.88);
    legRes_p->SetBorderSize(0);
    legRes_p->SetFillColor(0);
    legRes_p->SetFillStyle(0);
    legRes_p->SetTextFont(43);
    legRes_p->SetTextSize(16);
    
    const Int_t nPad = 2;
    TCanvas* resModCompCanv_p[nSystSmooth];
    TPad* resModCompPad_p[nSystSmooth][nPad];
    TH1D* resModHist_p[nSystSmooth];
    TH1D* resModHistRat_p[nSystSmooth];
    
    const Double_t yPadFrac = 0.35;
    const Double_t yPadLow[nPad] = {yPadFrac, 0.0};
    const Double_t yPadHi[nPad] = {1.0, yPadFrac};
    const Double_t padLeftMarg = 0.14;
    
    for(Int_t usI = 0; usI < nSystSmooth; ++usI){
      resModCompCanv_p[usI] = new TCanvas(("resModCompCanv_" + systSmooth[usI]).c_str(), "", 450, 450);
      resModCompCanv_p[usI]->SetTopMargin(0.001);
      resModCompCanv_p[usI]->SetRightMargin(0.001);
      resModCompCanv_p[usI]->SetLeftMargin(0.001);
      resModCompCanv_p[usI]->SetBottomMargin(0.001);    

      for(Int_t pI = 0; pI < nPad; ++pI){
	resModCompCanv_p[usI]->cd();
	resModCompPad_p[usI][pI] = new TPad(("resModCompPad_" + systSmooth[usI] + "_" + std::to_string(pI)).c_str(), "", 0.0, yPadLow[pI], 1.0, yPadHi[pI]);
	resModCompPad_p[usI][pI]->SetRightMargin(0.001);
	resModCompPad_p[usI][pI]->SetLeftMargin(padLeftMarg);
	
	if(pI == 0) resModCompPad_p[usI][pI]->SetTopMargin(0.01);
	else resModCompPad_p[usI][pI]->SetTopMargin(0.001);
	
	if(pI == 0) resModCompPad_p[usI][pI]->SetBottomMargin(0.001);
	else resModCompPad_p[usI][pI]->SetBottomMargin(padLeftMarg/yPadFrac);
	
	resModCompPad_p[usI][pI]->Draw("SAME");
      }
      
      resModHist_p[usI] = new TH1D(("resModHist_" + systSmooth[usI]).c_str(), ";Jet p_{T} (GeV/c);Fractional Systematic", 10, xPointMinVal, xPointMaxVal);
      centerTitles(resModHist_p[usI]);
      
      resModHist_p[usI]->SetMaximum(0.0);
      resModHist_p[usI]->SetMinimum(-1);
      resModCompPad_p[usI][0]->cd();
      resModHist_p[usI]->GetYaxis()->SetTitleFont(43);
      resModHist_p[usI]->GetYaxis()->SetTitleSize(16);
      resModHist_p[usI]->GetYaxis()->SetLabelFont(43);
      resModHist_p[usI]->GetYaxis()->SetLabelSize(10);
      resModHist_p[usI]->Draw("HIST");
      
      
      resModHistRat_p[usI] = new TH1D(("resModHistRat_" + systSmooth[usI]).c_str(), ";Jet p_{T} (GeV/c);Old/New", 10, xPointMinVal, xPointMaxVal);
      centerTitles(resModHistRat_p[usI]);
      
      resModHistRat_p[usI]->SetMaximum(1.1);
      resModHistRat_p[usI]->SetMinimum(0.4);
      resModCompPad_p[usI][1]->cd();
      
      resModHistRat_p[usI]->GetYaxis()->SetTitleFont(43);
      resModHistRat_p[usI]->GetYaxis()->SetTitleSize(16);
      resModHistRat_p[usI]->GetYaxis()->SetLabelFont(43);
      resModHistRat_p[usI]->GetYaxis()->SetLabelSize(10);
      
      resModHistRat_p[usI]->GetXaxis()->SetTitleFont(43);
      resModHistRat_p[usI]->GetXaxis()->SetTitleSize(16);
      resModHistRat_p[usI]->GetXaxis()->SetLabelFont(43);
      resModHistRat_p[usI]->GetXaxis()->SetLabelSize(10);
      
      resModHistRat_p[usI]->GetXaxis()->SetTitleOffset(2.0);
      
      resModHistRat_p[usI]->Draw("HIST");
    }
    
    std::vector<TH1D*> prevJERData;
    
    for(Int_t mI = 0; mI < nResponseMod; ++mI){
      const std::string responseStr = prettyString(responseMod.at(mI), 2, true);
      
      for(Int_t usI = 0; usI < nSystSmooth; ++usI){
	TLatex* label_p = new TLatex();
	label_p->SetTextFont(43);
	label_p->SetTextSize(16);
	
	Int_t nYBins = 20;
	Double_t yBins[nYBins];
	
	TCanvas* spectCanv_p = new TCanvas("spectCanv_p", "", 450, 450);
	spectCanv_p->SetTopMargin(0.08);
	spectCanv_p->SetRightMargin(0.01);
	spectCanv_p->SetLeftMargin(0.14);
	spectCanv_p->SetBottomMargin(0.14);
	
	const std::string yAxisTitle = "#frac{d^{2}#sigma^{pp}}{dp_{T}d#eta} or #frac{1}{N_{evt}}#frac{1}{#LTTAA#GT}#frac{d^{2}N^{PbPb}}{dp_{T}d#eta} [#frac{nb}{GeV/c}]";
	TH1D* tempHist_p = new TH1D("tempHist_p", (";Jet p_{T} (GeV/c);" + yAxisTitle).c_str(), 10, xMinVal, xMaxVal);
	centerTitles(tempHist_p);
	
	tempHist_p->GetXaxis()->SetTitleOffset(1.8); 
	tempHist_p->GetYaxis()->SetTitleOffset(1.6);
	
	for(Int_t bIX = 0; bIX < jtPtUnfolded_RecoGenAsymm_PP_h[ppPos][idPos][mI][absEtaPos][0]->GetNbinsX(); ++bIX){
	  bool binIsBad = jtPtUnfolded_RecoGenAsymm_PP_h[ppPos][idPos][mI][absEtaPos][0]->GetBinCenter(bIX+1) < xPointMinVal;
	  binIsBad = binIsBad || jtPtUnfolded_RecoGenAsymm_PP_h[ppPos][idPos][mI][absEtaPos][0]->GetBinCenter(bIX+1) > xPointMaxVal;
	  if(binIsBad){
	    jtPtUnfolded_RecoGenAsymm_PP_h[ppPos][idPos][mI][absEtaPos][0]->SetBinContent(bIX+1, 0.0);
	    jtPtUnfolded_RecoGenAsymm_PP_h[ppPos][idPos][mI][absEtaPos][0]->SetBinError(bIX+1, 0.0);
	  }
	}
	
	Double_t minVal = getMinGTZero(jtPtUnfolded_RecoGenAsymm_PP_h[ppPos][idPos][mI][absEtaPos][0]);
	Double_t maxVal = getMax(jtPtUnfolded_RecoGenAsymm_PP_h[ppPos][idPos][mI][absEtaPos][0]);

	for(Int_t cI = 0; cI < nCentBins; ++cI){	  
	  Double_t xPointMaxValUsed = xPointMaxVal;
	  if(centBinsLow[cI] >= 50) xPointMaxValUsed = xPointMaxVal5090;
	  
	  for(Int_t sI = 0; sI < nSyst; ++sI){
	    scaleCentralAndErrorValues(jtPtUnfolded_RecoGenAsymm_PbPb_h[tI][cI][idPos][mI][absEtaPos][sI], centBinsScalingFact.at(cI));
	  }
      	  
	  for(Int_t bIX = 0; bIX < jtPtUnfolded_RecoGenAsymm_PbPb_h[tI][cI][idPos][mI][absEtaPos][0]->GetNbinsX(); ++bIX){
	    bool binIsBad = jtPtUnfolded_RecoGenAsymm_PbPb_h[tI][cI][idPos][mI][absEtaPos][0]->GetBinCenter(bIX+1) < xPointMinVal;
	    binIsBad = binIsBad || jtPtUnfolded_RecoGenAsymm_PbPb_h[tI][cI][idPos][mI][absEtaPos][0]->GetBinCenter(bIX+1) > xPointMaxValUsed;
	    if(binIsBad){
	      jtPtUnfolded_RecoGenAsymm_PbPb_h[tI][cI][idPos][mI][absEtaPos][0]->SetBinContent(bIX+1, 0.0);
	      jtPtUnfolded_RecoGenAsymm_PbPb_h[tI][cI][idPos][mI][absEtaPos][0]->SetBinError(bIX+1, 0.0);
	    }
	  }
	  
	  Double_t tempMinVal = getMinGTZero(jtPtUnfolded_RecoGenAsymm_PbPb_h[tI][cI][idPos][mI][absEtaPos][0]);
	  Double_t tempMaxVal = getMax(jtPtUnfolded_RecoGenAsymm_PbPb_h[tI][cI][idPos][mI][absEtaPos][0]);
	  if(minVal > tempMinVal) minVal = tempMinVal;
	  if(maxVal < tempMaxVal) maxVal = tempMaxVal;
	}
	
	maxVal *= 100.;
	minVal /= 100.;
	
	tempHist_p->SetMaximum(maxVal);
	tempHist_p->SetMinimum(minVal);
	
	tempHist_p->GetYaxis()->SetLabelFont(43);
	tempHist_p->GetYaxis()->SetLabelSize(11);
	
	std::cout << "Title offset: " << tempHist_p->GetYaxis()->GetTitleOffset() << std::endl;
	
	getLogBins(minVal, maxVal, nYBins, yBins);
	
	tempHist_p->DrawCopy("HIST E1 P");
	
	TLegend* leg_p = new TLegend(0.7, 0.65, 0.95, 0.88);
	leg_p->SetBorderSize(0);
	leg_p->SetFillColor(0);
	leg_p->SetFillStyle(0);
	leg_p->SetTextFont(43);
	leg_p->SetTextSize(16);
	
	jtPtUnfolded_RecoGenAsymm_PP_h[ppPos][idPos][mI][absEtaPos][0]->SetMarkerColor(kPalette.getColor(getColorPosFromCent("", true)));
	jtPtUnfolded_RecoGenAsymm_PP_h[ppPos][idPos][mI][absEtaPos][0]->SetLineColor(kPalette.getColor(getColorPosFromCent("", true)));
	jtPtUnfolded_RecoGenAsymm_PP_h[ppPos][idPos][mI][absEtaPos][0]->SetMarkerStyle(getStyleFromCent("", true));
	jtPtUnfolded_RecoGenAsymm_PP_h[ppPos][idPos][mI][absEtaPos][0]->SetMarkerSize(1);
	
	std::vector<TH1D*> systHistVectPP;
	systHistVectPP.reserve(nSyst-1);
	for(Int_t sI = 1; sI < nSyst; ++sI){
	  systHistVectPP.push_back((TH1D*)jtPtUnfolded_RecoGenAsymm_PP_h[ppPos][idPos][mI][absEtaPos][sI]->Clone(("systPP_" + systStr.at(sI)).c_str()));
	}      
	
	Int_t slideSubVal = 0;
	if(usI == 0) pdfPerSlide.push_back({});
	else slideSubVal = nCentBins+1;
	
	std::vector<std::string> tempSystStr(systStr.begin()+1, systStr.end());
	std::vector<double> systValVectPP = getSyst(dirStr, jtPtUnfolded_RecoGenAsymm_PP_h[ppPos][idPos][mI][absEtaPos][0], systHistVectPP, tempSystStr, xPointMinVal, xPointMaxVal, &(pdfPerSlide.at(pdfPerSlide.size()- 1 - slideSubVal)), systSmoothBool[usI], systToCombo);
    
	resModCompCanv_p[usI]->cd();
	resModCompPad_p[usI][0]->cd();
	Int_t resStyle = -1;
	Int_t resCol = kPalette.getColor(getColorPosFromCent("", true));
	std::string drawStr = "HIST P SAME";
	if(mI == 0) resStyle = getStyleFromCent("", true);
	else resStyle = getAltStyleFromCent("", true);
	
	if(jerDataPos >= 0){
	  systHistVectPP.at(jerDataPos)->SetMarkerSize(0.8);
	  systHistVectPP.at(jerDataPos)->SetMarkerStyle(resStyle);
	  systHistVectPP.at(jerDataPos)->SetMarkerColor(resCol);
	  systHistVectPP.at(jerDataPos)->SetLineColor(resCol);
	  
	  systHistVectPP.at(jerDataPos)->DrawCopy(drawStr.c_str());
	  
	  if(systHistVectPP.at(jerDataPos)->GetMaximum() > resModHist_p[usI]->GetMaximum()) resModHist_p[usI]->SetMaximum(systHistVectPP.at(jerDataPos)->GetMaximum());

	  if(mI == 0) prevJERData.push_back((TH1D*)systHistVectPP.at(jerDataPos)->Clone((std::string(systHistVectPP.at(jerDataPos)->GetName()) + "_CLONE").c_str()));
	  else{
	    if(usI == 0) systHistVectPP.at(jerDataPos)->Divide(prevJERData.at(0));
	    else if(usI == 1) systHistVectPP.at(jerDataPos)->Divide(prevJERData.at(1 + nCentBins));
	    
	    resModCompCanv_p[usI]->cd();
	    resModCompPad_p[usI][1]->cd();
	    
	    systHistVectPP.at(jerDataPos)->DrawCopy("HIST P SAME");	  
	  }
	}
   
	if(usI == 0) slideTitles.push_back("PP Systematic (" + jetPPList.at(ppPos) + ", " + responseStr +  ")");
	
	for(unsigned int sI = 0; sI < systHistVectPP.size(); ++sI){
	  delete systHistVectPP.at(sI);
	}
	     
	drawSyst(spectCanv_p, jtPtUnfolded_RecoGenAsymm_PP_h[ppPos][idPos][mI][absEtaPos][0], &systValVectPP, xPointMinVal, xPointMaxVal);
	
	for(Int_t cI = 0; cI < nCentBins; ++cI){
	  const std::string centStr = "Cent" + std::to_string(centBinsLow.at(cI)) + "to" + std::to_string(centBinsHi.at(cI));
	  const std::string centStr2 = std::to_string(centBinsLow.at(cI)) + "-" + std::to_string(centBinsHi.at(cI)) + "%";
	  
	  Double_t xPointMaxValUsed = xPointMaxVal;
	  if(centBinsLow[cI] >= 50) xPointMaxValUsed = xPointMaxVal5090;

	  jtPtUnfolded_RecoGenAsymm_PbPb_h[tI][cI][idPos][mI][absEtaPos][0]->SetMarkerColor(kPalette.getColor(getColorPosFromCent(centStr, false)));
	  jtPtUnfolded_RecoGenAsymm_PbPb_h[tI][cI][idPos][mI][absEtaPos][0]->SetLineColor(kPalette.getColor(getColorPosFromCent(centStr, false)));
	  jtPtUnfolded_RecoGenAsymm_PbPb_h[tI][cI][idPos][mI][absEtaPos][0]->SetMarkerStyle(getStyleFromCent(centStr, false));
	  jtPtUnfolded_RecoGenAsymm_PbPb_h[tI][cI][idPos][mI][absEtaPos][0]->SetMarkerSize(1);
	  
	  std::vector<TH1D*> systHistVectPbPb;
	  systHistVectPbPb.reserve(nSyst-1);
	  for(Int_t sI = 1; sI < nSyst; ++sI){
	    systHistVectPbPb.push_back((TH1D*)jtPtUnfolded_RecoGenAsymm_PbPb_h[tI][cI][idPos][mI][absEtaPos][sI]->Clone(("systPbPb_" + systStr.at(sI)).c_str()));
	  }      
	  
	  Int_t centSlideSubVal = 0;
	  if(usI == 0) pdfPerSlide.push_back({});
	  else centSlideSubVal = slideSubVal - 1 - cI;
	  
	  std::vector<double> systValVectPbPb = getSyst(dirStr, jtPtUnfolded_RecoGenAsymm_PbPb_h[tI][cI][idPos][mI][absEtaPos][0], systHistVectPbPb, tempSystStr, xPointMinVal, xPointMaxValUsed, &(pdfPerSlide.at(pdfPerSlide.size() - 1 - centSlideSubVal)), systSmoothBool[usI], systToCombo);
	  if(usI == 0) slideTitles.push_back(centStr2 + " Systematic (" + jetPbPbList.at(tI) + ", " + responseStr +  ")");
	  
	  resModCompCanv_p[usI]->cd();
	  resModCompPad_p[usI][0]->cd();
	  resStyle = -1;
	  resCol = kPalette.getColor(getColorPosFromCent(centStr, false));
	  drawStr = "HIST P SAME";
	  if(mI == 0) resStyle = getStyleFromCent(centStr, false);
	  else resStyle = getAltStyleFromCent(centStr, false);
	  
	  if(jerDataPos >= 0){
	    systHistVectPbPb.at(jerDataPos)->SetMarkerSize(0.8);
	    systHistVectPbPb.at(jerDataPos)->SetMarkerStyle(resStyle);
	    systHistVectPbPb.at(jerDataPos)->SetMarkerColor(resCol);
	    systHistVectPbPb.at(jerDataPos)->SetLineColor(resCol);
	    
	    systHistVectPbPb.at(jerDataPos)->DrawCopy(drawStr.c_str());
	    
	    if(systHistVectPbPb.at(jerDataPos)->GetMaximum() > resModHist_p[usI]->GetMaximum()) resModHist_p[usI]->SetMaximum(systHistVectPbPb.at(jerDataPos)->GetMaximum());
	    
	    if(mI == 0) prevJERData.push_back((TH1D*)systHistVectPbPb.at(jerDataPos)->Clone((std::string(systHistVectPbPb.at(jerDataPos)->GetName()) + "_CLONE").c_str()));	  
	    else{
	      if(usI == 0) systHistVectPbPb.at(jerDataPos)->Divide(prevJERData.at(0 + cI + 1));
	      else if(usI == 1) systHistVectPbPb.at(jerDataPos)->Divide(prevJERData.at(1 + nCentBins + cI + 1));
	      
	      resModCompCanv_p[usI]->cd();
	      resModCompPad_p[usI][1]->cd();
	      
	      systHistVectPbPb.at(jerDataPos)->DrawCopy("HIST P SAME");
	    }
	  }
	  
	  for(unsigned int sI = 0; sI < systHistVectPbPb.size(); ++sI){
	    delete systHistVectPbPb.at(sI);
	  }
	  drawSyst(spectCanv_p, jtPtUnfolded_RecoGenAsymm_PbPb_h[tI][cI][idPos][mI][absEtaPos][0], &systValVectPbPb, xPointMinVal, xPointMaxValUsed);	
	}
	
	jtPtUnfolded_RecoGenAsymm_PP_h[ppPos][idPos][mI][absEtaPos][0]->DrawCopy("HIST E1 P SAME");

	for(Int_t cI = 0; cI < nCentBins; ++cI){
	  jtPtUnfolded_RecoGenAsymm_PbPb_h[tI][cI][idPos][mI][absEtaPos][0]->DrawCopy("HIST E1 P SAME");
	}
	
	for(Int_t cI = nCentBins-1; cI >= 0; --cI){
	  const std::string centLegStr = std::to_string(centBinsLow.at(cI)) + "-" + std::to_string(centBinsHi.at(cI)) + "% x 10^{" + std::to_string(cI+1) + "}";
	  const std::string centLegStr2 = std::to_string(centBinsLow.at(cI)) + "-" + std::to_string(centBinsHi.at(cI)) + "%";
	  
	  jtPtUnfolded_RecoGenAsymm_PbPb_h[tI][cI][idPos][mI][absEtaPos][0]->SetFillColorAlpha(jtPtUnfolded_RecoGenAsymm_PbPb_h[tI][cI][idPos][mI][absEtaPos][0]->GetMarkerColor(), .25);      
	  leg_p->AddEntry(jtPtUnfolded_RecoGenAsymm_PbPb_h[tI][cI][idPos][mI][absEtaPos][0], centLegStr.c_str(), "P L F");	
	  if(usI == 0 && mI == 0) legRes_p->AddEntry(jtPtUnfolded_RecoGenAsymm_PbPb_h[tI][cI][idPos][mI][absEtaPos][0], centLegStr2.c_str(), "P L");	
	}
	
	jtPtUnfolded_RecoGenAsymm_PP_h[ppPos][idPos][mI][absEtaPos][0]->SetFillColorAlpha(jtPtUnfolded_RecoGenAsymm_PP_h[ppPos][idPos][mI][absEtaPos][0]->GetMarkerColor(), .25);      
	leg_p->AddEntry(jtPtUnfolded_RecoGenAsymm_PP_h[ppPos][idPos][mI][absEtaPos][0], "pp", "P L F");
	if(usI == 0 && mI == 0) legRes_p->AddEntry(jtPtUnfolded_RecoGenAsymm_PP_h[ppPos][idPos][mI][absEtaPos][0], "pp", "P L");
	
	gPad->SetLogy();
	gPad->SetLogx();
	gStyle->SetOptStat(0);
	
	drawWhiteBox(900, 1100, minVal/10, minVal*.9);

	label_p->DrawLatex(200, yBins[0]/4., "200");
	label_p->DrawLatex(400, yBins[0]/4., "400");
	label_p->DrawLatex(600, yBins[0]/4., "600");
	label_p->DrawLatex(800, yBins[0]/4., "800");
	label_p->DrawLatex(1000, yBins[0]/4., "1000");
        
	label_p->DrawLatex(250, yBins[19], ("anti-k_{T} " + rStr).c_str());
	label_p->DrawLatex(250, yBins[18], "|#eta_{jets}| < 2");
	
	label_p->SetNDC();
	
	label_p->DrawLatex(0.1, 0.95, "#bf{CMS Preliminary}");
	label_p->DrawLatex(0.4, 0.95, "27.4 pb^{-1} pp + 404 #mub^{-1} PbPb (5.02 TeV)");
	
	leg_p->Draw("SAME");
	
	gPad->RedrawAxis();
	gPad->SetTicks(1,2);
	
	std::string saveName = "spectra_" + jetPbPbList.at(tI) + "_R" + rValStr + "_" + idNameStr + "_" + plotAbsEtaStr + "_" + plotBayesStr + "_" + responseStr + "_" + systSmooth[usI] + "_" + dateStr + ".pdf";

	if(usI == 0){
	  slideTitles.push_back("Spectra (" + responseStr + ")");
	  pdfPerSlide.push_back({});
	  
	  slideTitlesMain.push_back("Spectra (" + responseStr + ")");
	  pdfPerSlideMain.push_back({});
	}
	
	pdfPerSlide.at(pdfPerSlide.size() - 1).push_back(saveName);
	pdfPerSlideMain.at(pdfPerSlideMain.size() - 1).push_back(saveName);

	saveName = "pdfDir/" + dateStr + "/PlotUnfold/" + saveName;
	
	quietSaveAs(spectCanv_p, saveName);
	
	delete spectCanv_p;
	delete label_p;
	delete leg_p;
	
	for(Int_t cI = 0; cI < nCentBins; ++cI){
	  for(Int_t sI = 0; sI < nSyst; ++sI){
	    scaleCentralAndErrorValues(jtPtUnfolded_RecoGenAsymm_PbPb_h[tI][cI][idPos][mI][absEtaPos][sI], 1./centBinsScalingFact.at(cI));
	  }
	}
      }
    }
    
    for(Int_t usI = 0; usI < nSystSmooth; ++usI){
      TLatex* label_p = new TLatex();
      label_p->SetTextFont(43);
      label_p->SetTextSize(16);
      
      resModCompCanv_p[usI]->cd();
      resModCompPad_p[usI][0]->cd();
      Double_t maxVal = resModHist_p[usI]->GetMaximum()*1.2;
      resModHist_p[usI]->SetMaximum(maxVal);
      resModHist_p[usI]->SetMinimum(0.0);

      label_p->DrawLatex(500, .12, rStr.c_str());
      
      legRes_p->Draw("SAME");
      
      gPad->SetLogx();
      gPad->Modified();
      
      resModCompCanv_p[usI]->cd();
      resModCompPad_p[usI][1]->cd();
      gPad->SetLogx();
      
      drawWhiteBox(900, 1100, 0.2, .39);
      
      label_p->DrawLatex(200, .3, "200");
      label_p->DrawLatex(400, .3, "400");
      label_p->DrawLatex(600, .3, "600");
      label_p->DrawLatex(800, .3, "800");
      
      std::string saveName = "resErrComp_" + jetPbPbList.at(tI) + "_R" + rValStr + "_" + idNameStr + "_" + plotAbsEtaStr + "_" + plotBayesStr + "_" + systSmooth[usI] + "_"+ dateStr + ".pdf";

      if(usI == 0){
	slideTitles.push_back("Resolution Error Comp.");
	pdfPerSlide.push_back({});
      }
      
      pdfPerSlide.at(pdfPerSlide.size() - 1).push_back(saveName);
      
      saveName = "pdfDir/" + dateStr + "/PlotUnfold/" + saveName;
      
      quietSaveAs(resModCompCanv_p[usI], saveName);
      
      for(Int_t pI = 0; pI < nPad; ++pI){delete resModCompPad_p[usI][pI];}
      
      delete resModCompCanv_p[usI];
      delete resModHist_p[usI];
    }
    
    delete legRes_p;
    
    for(unsigned int i = 0; i < prevJERData.size(); ++i){
      delete prevJERData.at(i);
      prevJERData.at(i) = NULL;
    }
  }

  //Traditional RAA

  //R scan RAA

  std::cout << "Closing files..." << std::endl;
  
  const Int_t nHistDim = cutPropPbPb.GetNHistDim() + cutPropPP.GetNHistDim();
  histTagPbPb.insert(histTagPbPb.end(), histTagPP.begin(), histTagPP.end());
  histBestBayesPbPb.insert(histBestBayesPbPb.end(), histBestBayesPP.begin(), histBestBayesPP.end());
  
  std::cout << "nHistDim: " << nHistDim << ", " << cutPropPbPb.GetNHistDim() << ", " << cutPropPP.GetNHistDim() << std::endl;
  std::cout << " vect: " << histTagPbPb.size() << ", " << histTagPP.size() << std::endl;
  std::cout << " vect: " << histBestBayesPbPb.size() << ", " << histBestBayesPP.size() << std::endl;
  
  if(nHistDim != (int)(histTagPbPb.size())){
    std::cout << "Warning: Error in combination of histTags from pp and pbpb" << std::endl;
  }
  if(nHistDim != (int)(histBestBayesPbPb.size())){
    std::cout << "Warning: Error in combination of histTags from pp and pbpb" << std::endl;
  }
  
  cutPropPbPb.SetNHistDim(nHistDim);
  cutPropPbPb.SetHistTag(histTagPbPb);
  cutPropPbPb.SetHistBestBayes(histBestBayesPbPb);
  
  TDirectory* cutDir_p = (TDirectory*)outFile_p->mkdir("cutDir");
  TDirectory* subDir_p = (TDirectory*)cutDir_p->mkdir("subDir");
  TDirectory* unfoldDir_p = (TDirectory*)cutDir_p->mkdir("unfoldDir");
  
  if(!cutPropPbPb.WriteAllVarToFile(outFile_p, cutDir_p, subDir_p, unfoldDir_p)) std::cout << "Warning: Cut writing has failed" << std::endl;  
  
  outFile_p->Close();
  delete outFile_p;
  
  inFilePP_p->Close();
  delete inFilePP_p;

  inFilePbPb_p->Close();
  delete inFilePbPb_p;

  texSlideCreator tex;
  tex.Clean();
  tex.Init(outFileName);
  tex.InitDir("pdfDir/" + dateStr + "/PlotUnfold/");

  if(nJtPbPb == 1) tex.InitTag("AllPlots_" + jetPbPbList.at(0));
  else tex.InitTag("AllPlots");

  tex.SetAuthor("Christopher McGinn");
  tex.SetSlideTitles(slideTitles);
  tex.SetSlidePdfs(pdfPerSlide);
  if(!(tex.CreateTexSlides())){
    std::cout << "Warning: .tex slide creation failed" << std::endl;
  }


  texSlideCreator texMain;
  texMain.Clean();
  texMain.Init(outFileName);
  if(nJtPbPb == 1) texMain.InitTag("MainPlots_" + jetPbPbList.at(0));
  else texMain.InitTag("MainPlots");
  texMain.SetAuthor("Christopher McGinn");
  texMain.SetSlideTitles(slideTitlesMain);
  texMain.SetSlidePdfs(pdfPerSlideMain);
  if(!(texMain.CreateTexSlides())){
    std::cout << "Warning: .tex slide creation failed" << std::endl;
  }

  std::cout << "Job complete" << std::endl;

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 3){
    std::cout << "Usage: ./bin/plotUnfoldedSpectra.exe <inFileNamePP> <inFileNamePbPb>" << std::endl;
    return 1;
  }

  int retVal = 0;
  std::cout << __LINE__ << std::endl;
  retVal += plotUnfoldedSpectra(argv[1], argv[2]);
  return retVal;
}
