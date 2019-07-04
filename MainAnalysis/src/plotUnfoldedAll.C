//cpp dependencies
#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

//ROOT dependencies
#include "TCanvas.h"
#include "TDatime.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TH1D.h"
#include "TKey.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMath.h"
#include "TPad.h"
#include "TStyle.h"
#include "TTree.h"

//Local dependencies
#include "MainAnalysis/include/canvNDCToXY.h"
#include "MainAnalysis/include/cutPropagator.h"
#include "MainAnalysis/include/macroHistToSubsetHist.h"
#include "MainAnalysis/include/smallOrLargeR.h"
#include "MainAnalysis/include/systFunctions.h"

//Non-Local dependencies
#include "Utility/include/checkMakeDir.h"
#include "Utility/include/fileUtilities.h"
#include "Utility/include/histDefUtility.h"
#include "Utility/include/kirchnerPalette.h"
#include "Utility/include/lumiAndTAAUtil.h"
#include "Utility/include/plotUtilities.h"
#include "Utility/include/returnRootFileContentsList.h"
#include "Utility/include/vanGoghPalette.h"

template <class T>
bool compWithWarning(const std::string typeStr, T a, T b)
{
  bool val = a == b;
  if(!val) std::cout << typeStr << " \'" << a << "\' from one file does not match \'" << b << "\' from other. fail";
  return val;
}

ULong64_t getKey(ULong64_t binsI, ULong64_t jI, ULong64_t cI, ULong64_t idI, ULong64_t rI, ULong64_t aI, ULong64_t sI){return binsI + jI*10 + 1000*cI + 100000*idI + 10000000*rI + 1000000000*aI + 100000000000*sI;}

ULong64_t getKey(ULong64_t binsI, ULong64_t jI, ULong64_t cI, ULong64_t idI, ULong64_t rI, ULong64_t aI, ULong64_t sI, ULong64_t bI){return binsI + 10*jI + 1000*cI + 100000*idI + 10000000*rI + 1000000000*aI + 100000000000*sI + 100000000000000*bI;}

int posInStrVect(std::string strToCheck, std::string front, std::vector<std::string> vect, std::string back)
{
  int pos = -1;
  for(unsigned i = 0; i < vect.size(); ++i){
    if(vect[i].size() == 0) continue;
    if(strToCheck.find(front + vect[i] + back) != std::string::npos){
      pos = i;
      break;
    }
  }

  return pos;
}

int posInStrVectExact(std::string inStr, std::vector<std::string> inVect)
{
  Int_t pos = -1;
  for(unsigned int i = 0; i < inVect.size(); ++i){
    if(isStrSame(inStr, inVect[i])){pos = i; break;}
  }
  return pos;
}

void defineCanv(TCanvas* canv_p)
{
  canv_p->SetTopMargin(0.08);
  canv_p->SetRightMargin(0.01);
  canv_p->SetLeftMargin(0.18);
  canv_p->SetBottomMargin(0.1);

  return;
}

template <class T>
bool setGeneralPtBins(const Int_t nGeneralPtBins, Double_t generalPtBins[], std::vector<T> binVect)
{
  if(nGeneralPtBins+1 < (Int_t)binVect.size()){
    std::cout << "testGeneralPtBins input binVect size \'" << binVect.size() << "\' is greater than the generalBins max+1, \'" << nGeneralPtBins << "+1\'. return false"  << std::endl;
    return false;
  }
  for(unsigned int vI = 0; vI < binVect.size(); ++vI){generalPtBins[vI] = binVect[vI];}
  return true;
}

void divHistByWidth(TH1D* hist_p)
{
  for(Int_t bIX = 0; bIX < hist_p->GetNbinsX(); ++bIX){
    Double_t binWidth = hist_p->GetBinWidth(bIX+1);
    hist_p->SetBinContent(bIX+1, hist_p->GetBinContent(bIX+1)/binWidth);
    hist_p->SetBinError(bIX+1, hist_p->GetBinError(bIX+1)/binWidth);
  }
  return;
}

void divHistByWidth(std::vector<TH1D*> hist_p)
{
  for(unsigned int sI = 0; sI < hist_p.size(); ++sI){
    if(hist_p[sI] == NULL) continue;

    for(Int_t bIX = 0; bIX < hist_p[sI]->GetNbinsX(); ++bIX){
      Double_t binWidth = hist_p[sI]->GetBinWidth(bIX+1);
      hist_p[sI]->SetBinContent(bIX+1, hist_p[sI]->GetBinContent(bIX+1)/binWidth);
      hist_p[sI]->SetBinError(bIX+1, hist_p[sI]->GetBinError(bIX+1)/binWidth);
    }
  }
  return;
}

void scaleHist(TH1D* hist_p, Double_t scaleFactor)
{
  for(Int_t bIX = 0; bIX < hist_p->GetNbinsX(); ++bIX){
    hist_p->SetBinContent(bIX+1, hist_p->GetBinContent(bIX+1)*scaleFactor);
    hist_p->SetBinError(bIX+1, hist_p->GetBinError(bIX+1)*scaleFactor);
  }
  return;
}

void scaleVect(std::vector<Double_t>* vect, Double_t scaleFactor)
{
  for(ULong64_t vI = 0; vI < (ULong64_t)vect->size(); ++vI){   
    (*vect)[vI] *= scaleFactor;
  }
  return;
}

Double_t getHistMax(TH1D* hist_p)
{
  Double_t max = hist_p->GetMinimum();
  for(Int_t bIX = 0; bIX < hist_p->GetNbinsX(); ++bIX){
    if(max < hist_p->GetBinContent(bIX+1)) max = hist_p->GetBinContent(bIX+1);
  }
  return max;
}

Double_t getHistMin(TH1D* hist_p)
{
  Double_t min = hist_p->GetMaximum();
  for(Int_t bIX = 0; bIX < hist_p->GetNbinsX(); ++bIX){
    if(min > hist_p->GetBinContent(bIX+1)) min = hist_p->GetBinContent(bIX+1);
  }
  return min;
}

Double_t getHistMinGTZero(TH1D* hist_p)
{
  Double_t min = hist_p->GetMaximum();
  for(Int_t bIX = 0; bIX < hist_p->GetNbinsX(); ++bIX){
    if(hist_p->GetBinContent(bIX+1) <= 0) continue;
    if(min > hist_p->GetBinContent(bIX+1)) min = hist_p->GetBinContent(bIX+1);
  }
  return min;
}

void createRAA(TH1D* pbpbHist_p, TH1D* ppHist_p)
{
  for(Int_t bIX = 0; bIX < pbpbHist_p->GetNbinsX(); ++bIX){
    Double_t val = pbpbHist_p->GetBinContent(bIX+1)/ppHist_p->GetBinContent(bIX+1);
    Double_t relErr1 = pbpbHist_p->GetBinError(bIX+1)/pbpbHist_p->GetBinContent(bIX+1);
    Double_t relErr2 = ppHist_p->GetBinError(bIX+1)/ppHist_p->GetBinContent(bIX+1);
    Double_t relErr = TMath::Sqrt(relErr1*relErr1 + relErr2*relErr2);
    Double_t err = val*relErr;

    pbpbHist_p->SetBinContent(bIX+1, val);
    pbpbHist_p->SetBinError(bIX+1, err);
  }
 
  return;
}

template <class T>
bool binRCheck(bool maxOrMin, std::string maxOrMinStr, Double_t maxOrMinVal, Int_t nBins, std::vector<T> bins)
{
  if(!maxOrMin){
    std::cout << maxOrMinStr << " largeR val \'" << maxOrMinVal << "\' is not found in bins: ";
    for(Int_t gI = 0; gI < nBins; ++gI){
      std::cout << bins[gI] << ", ";
    }
    std::cout << bins[nBins] << ". return 1" << std::endl;
    return false;
  }
  return true;
}

void labelStandard(TLatex* label_p)
{
  label_p->SetNDC();
  label_p->SetTextFont(43);
  label_p->SetTextSize(16);
 
  return;
}

void legStandard(TLegend* leg_p)
{
  leg_p->SetBorderSize(0);
  leg_p->SetFillColor(0);
  leg_p->SetFillStyle(0);
  leg_p->SetTextFont(43);
  leg_p->SetTextSize(16);
  leg_p->SetTextColor(1);
 
  return;
}

void setupSyst(const ULong64_t nReducedSyst, std::vector<std::vector<Int_t> > nBins, std::vector<std::vector<Int_t> > nBinsPP, const ULong64_t nJtAlgos, const ULong64_t nCentBins, std::vector<std::vector<std::vector<std::vector<Double_t> > > >* systPbPb, std::vector<std::vector<std::vector<Double_t> > >* systSumPbPb, std::vector<std::vector<std::vector<std::vector<Double_t> > > >* systPP=NULL, std::vector<std::vector<std::vector<Double_t> > >* systSumPP=NULL)
{
  
  for(ULong64_t jI = 0; jI < (ULong64_t)nJtAlgos; ++jI){
    systSumPbPb->push_back({});
    systPbPb->push_back({});

    for(ULong64_t cI = 0; cI < (ULong64_t)nCentBins; ++cI){
      (*systSumPbPb)[jI].push_back({});
      (*systPbPb)[jI].push_back({});
      
      for(ULong64_t sI = 0; sI < nReducedSyst; ++sI){
	(*systPbPb)[jI][cI].push_back({});
      }
    }
  }

  for(ULong64_t jI = 0; jI < (ULong64_t)nJtAlgos; ++jI){
    for(ULong64_t cI = 0; cI < (ULong64_t)nCentBins; ++cI){
      for(Int_t bI = 0; bI < nBins[jI][cI]; ++bI){
	(*systSumPbPb)[jI][cI].push_back(0.0);
      }

      for(ULong64_t sI = 0; sI < nReducedSyst; ++sI){
	for(Int_t bI = 0; bI < nBins[jI][cI]; ++bI){
	  (*systPbPb)[jI][cI][sI].push_back(1000000000000.);
	}
      }
    }
  }

  if(systPP != NULL){
    for(ULong64_t jI = 0; jI < (ULong64_t)nJtAlgos; ++jI){
      systSumPP->push_back({});
      systPP->push_back({});
      for(ULong64_t cI = 0; cI < (ULong64_t)nCentBins; ++cI){
	(*systSumPP)[jI].push_back({});
	(*systPP)[jI].push_back({});
		
	for(ULong64_t sI = 0; sI < nReducedSyst; ++sI){
	  (*systPP)[jI][cI].push_back({});
	}
      }
    }

    for(ULong64_t jI = 0; jI < (ULong64_t)nJtAlgos; ++jI){
      for(ULong64_t cI = 0; cI < (ULong64_t)nCentBins; ++cI){
	for(Int_t bI = 0; bI < nBinsPP[jI][cI]; ++bI){
	  (*systSumPP)[jI][cI].push_back(0.0);
	}
	
	for(ULong64_t sI = 0; sI < nReducedSyst; ++sI){
	  for(Int_t bI = 0; bI < nBinsPP[jI][cI]; ++bI){
	    (*systPP)[jI][cI][sI].push_back(1000000000000.);
	  }
	}
      }
    }
  }

  return;
}

void defaultPlotSet(TH1D* hist_p, std::string centStr = "")
{
  bool isPP = true;
  if(centStr.size() != 0) isPP = false;
 
  kirchnerPalette kPalette;

  hist_p->SetMarkerColor(kPalette.getColor(getColorPosFromCent(centStr, isPP)));
  hist_p->SetLineColor(kPalette.getColor(getColorPosFromCent(centStr, isPP)));
  hist_p->SetMarkerStyle(getStyleFromCent(centStr, isPP));
  hist_p->SetMarkerSize(1);
  centerTitles(hist_p);

  hist_p->GetXaxis()->SetTitleOffset(hist_p->GetXaxis()->GetTitleOffset()*1.3);
  hist_p->GetYaxis()->SetTitleOffset(2.2);
  hist_p->GetXaxis()->SetLabelColor(0);

  return;
}

void defaultPlotSet(TGraphAsymmErrors* hist_p, std::string centStr = "")
{
  bool isPP = true;
  if(centStr.size() != 0) isPP = false;
 
  kirchnerPalette kPalette;

  hist_p->SetMarkerColor(kPalette.getColor(getColorPosFromCent(centStr, isPP)));
  hist_p->SetLineColor(kPalette.getColor(getColorPosFromCent(centStr, isPP)));
  hist_p->SetMarkerStyle(getStyleFromCent(centStr, isPP));
  hist_p->SetMarkerSize(1);
  centerTitles(hist_p);

  hist_p->GetXaxis()->SetTitleOffset(hist_p->GetXaxis()->GetTitleOffset()*1.3);
  hist_p->GetYaxis()->SetTitleOffset(2.2);
  hist_p->GetXaxis()->SetLabelColor(0);

  hist_p->GetXaxis()->SetTitleFont(43);
  hist_p->GetYaxis()->SetTitleFont(43);
  hist_p->GetXaxis()->SetTitleSize(16);
  hist_p->GetYaxis()->SetTitleSize(16);

  hist_p->GetXaxis()->SetLabelFont(43);
  hist_p->GetYaxis()->SetLabelFont(43);
  hist_p->GetXaxis()->SetLabelSize(14);
  hist_p->GetYaxis()->SetLabelSize(14);

  return;
}

void defaultPlotSetR(TH1D* hist_p, std::string rStr = "")
{ 
  kirchnerPalette kPalette;

  hist_p->SetMarkerColor(kPalette.getColor(getColorPosFromAlgo(rStr)));
  hist_p->SetLineColor(kPalette.getColor(getColorPosFromAlgo(rStr)));
  hist_p->SetMarkerStyle(getStyleFromAlgo(rStr));
  hist_p->SetMarkerSize(1);
  centerTitles(hist_p);

  hist_p->GetXaxis()->SetTitleOffset(hist_p->GetXaxis()->GetTitleOffset()*1.3);
  hist_p->GetYaxis()->SetTitleOffset(2.2);
  hist_p->GetXaxis()->SetLabelColor(0);

  hist_p->GetXaxis()->SetTitleFont(43);
  hist_p->GetYaxis()->SetTitleFont(43);
  hist_p->GetXaxis()->SetTitleSize(16);
  hist_p->GetYaxis()->SetTitleSize(16);

  hist_p->GetXaxis()->SetLabelFont(43);
  hist_p->GetYaxis()->SetLabelFont(43);
  hist_p->GetXaxis()->SetLabelSize(14);
  hist_p->GetYaxis()->SetLabelSize(14);

  return;
}

void writeAll(TFile* outFile_p, std::vector<TH1D*> hist_p)
{
  outFile_p->cd();
  for(auto const & hist : hist_p){
    if(hist == NULL) continue;
    hist->Write("", TObject::kOverwrite);
  }
  return;
}

void deleteAll(TFile* outFile_p, std::vector<TH1D*>* hist_p)
{
  outFile_p->cd();
  for(unsigned int hI = 0; hI < hist_p->size(); ++hI){
    if((*hist_p)[hI] == NULL) continue;
    delete (*hist_p)[hI];
  }
  return;
}

std::string labelAbsEta(std::string absEtaString)
{
  absEtaString.replace(absEtaString.find("AbsEta"), 6, "");
  while(absEtaString.find("p") != std::string::npos){absEtaString.replace(absEtaString.find("p"), 1, ".");}
  absEtaString.replace(absEtaString.find("to"), 2, " < |#eta_{Jet}| < ");
  while(absEtaString.substr(0,1).find(" ") != std::string::npos){absEtaString.replace(0,1,"");}
  if(absEtaString.substr(0,3).find("0.0") != std::string::npos){
    unsigned int endPos = absEtaString.find("<");
    bool goodToRemove = true;
    for(unsigned int iter = 3; iter < endPos; ++iter){
      std::string tempStr = absEtaString.substr(iter,1);
      if(tempStr.find(" ") != std::string::npos) continue;
      if(tempStr.find("0") != std::string::npos) continue;
      if(tempStr.find("<") != std::string::npos) continue;

      goodToRemove = false;
      break;
    }

    if(goodToRemove) absEtaString.replace(0, absEtaString.find("<")+1, "");
  }
  while(absEtaString.substr(0,1).find(" ") != std::string::npos){absEtaString.replace(0,1,"");}

  unsigned int startPos = absEtaString.find("<");
  unsigned int endPos = absEtaString.rfind("<");

  if(startPos != endPos){
    std::string tempStr = absEtaString.substr(0, startPos);
    absEtaString.replace(0, startPos, "");
    while(isStrSame(tempStr.substr(tempStr.size()-1, 1), " ") || isStrSame(tempStr.substr(tempStr.size()-1, 1), "0")){
      tempStr = tempStr.substr(0,tempStr.size()-1);
    }
    tempStr = tempStr.substr(0,tempStr.size()-1);
    absEtaString = tempStr + " " + absEtaString;
  }

  std::string tempStr = absEtaString.substr(endPos, absEtaString.size());
  absEtaString.replace(endPos, absEtaString.size(), "");
  while(isStrSame(tempStr.substr(tempStr.size()-1, 1), " ") || isStrSame(tempStr.substr(tempStr.size()-1, 1), "0")){
    tempStr = tempStr.substr(0,tempStr.size()-1);
  }
  tempStr = tempStr.substr(0,tempStr.size()-1);
  absEtaString = absEtaString + tempStr;
 
  return absEtaString;
}

void drawXLabels(TPad* canv_p, TH1D* hist_p, TLatex* label_p, const Int_t nXVals, const Double_t xVals[])
{
  canvNDCToXY labelAid(canv_p, hist_p);

  canv_p->cd();	 

  for(Int_t xI = 0; xI < nXVals; ++xI){
    if(xVals[xI] < hist_p->GetBinLowEdge(1)) continue;
    if(xVals[xI] >= hist_p->GetBinLowEdge(hist_p->GetNbinsX()+1)) continue;

    if(xVals[0] > 1){
      label_p->DrawLatex(labelAid.getXRelFromAbs(xVals[xI], canv_p->GetLogx()), canv_p->GetBottomMargin()-0.12, std::to_string((Int_t)xVals[xI]).c_str());
    }
    else{
      label_p->DrawLatex(labelAid.getXRelFromAbs(xVals[xI], canv_p->GetLogx()), canv_p->GetBottomMargin()-0.12, prettyString(xVals[xI], 1, false).c_str());
    }
  }    

  return;
}
void drawAllLabels(TCanvas* canv_p, TPad* pad_p, TH1D* hist_p, TLatex* label_p, const Int_t nXVals, const Double_t xVals[], const std::string rValStr, const std::string absEtaStr, const std::string pTStr ="", const bool drawPPOnly = false)
{
  canvNDCToXY labelAid(canv_p, hist_p);

  canv_p->cd();	 
  if(pad_p != NULL) pad_p->cd();

  for(Int_t xI = 0; xI < nXVals; ++xI){
    if(xVals[xI] < hist_p->GetBinLowEdge(1)) continue;
    if(xVals[xI] >= hist_p->GetBinLowEdge(hist_p->GetNbinsX()+1)) continue;

    if(pad_p == NULL){
      if(xVals[0] > 1){
	label_p->DrawLatex(labelAid.getXRelFromAbs(xVals[xI], canv_p->GetLogx()), canv_p->GetBottomMargin()-0.035, std::to_string((Int_t)xVals[xI]).c_str());
      }
      else{
	label_p->DrawLatex(labelAid.getXRelFromAbs(xVals[xI], canv_p->GetLogx()), canv_p->GetBottomMargin()-0.035, prettyString(xVals[xI], 1, false).c_str());
      }
    }  
  }

  if(pad_p == NULL){
    label_p->DrawLatex(gPad->GetLeftMargin()+0.04, 1.0-gPad->GetTopMargin()-0.05, "#bf{CMS Work In-Progress}");
    if(!drawPPOnly) label_p->DrawLatex(gPad->GetLeftMargin(), 1.0-gPad->GetTopMargin()+0.02, "#sqrt{s_{NN}} = 5.02 TeV, PbPb 404 #mub^{-1}, pp 27.4 pb^{-1}");
    else label_p->DrawLatex(gPad->GetLeftMargin(), 1.0-gPad->GetTopMargin()+0.02, "#sqrt{s_{NN}} = 5.02 TeV, pp 27.4 pb^{-1}");
    if(rValStr.size() != 0) label_p->DrawLatex(0.7, 1.0-gPad->GetTopMargin()-0.05, ("anti-k_{t} R=" + rValStr + " jets").c_str());
    else if(pTStr.size() != 0) label_p->DrawLatex(0.7, 1.0-gPad->GetTopMargin()-0.05, pTStr.c_str());

    label_p->DrawLatex(0.7, 1.0-gPad->GetTopMargin()-0.10, labelAbsEta(absEtaStr).c_str());
  }
  else{
    label_p->DrawLatex(pad_p->GetLeftMargin()+0.04, 1.0-pad_p->GetTopMargin()-0.05, "#bf{CMS Work In-Progress}");
    if(!drawPPOnly) label_p->DrawLatex(pad_p->GetLeftMargin(), 1.0-pad_p->GetTopMargin()+0.02, "#sqrt{s_{NN}} = 5.02 TeV, PbPb 404 #mub^{-1}, pp 27.4 pb^{-1}");
    else label_p->DrawLatex(pad_p->GetLeftMargin(), 1.0-pad_p->GetTopMargin()+0.02, "#sqrt{s_{NN}} = 5.02 TeV, pp 27.4 pb^{-1}");
    if(rValStr.size() != 0) label_p->DrawLatex(0.7, 1.0-pad_p->GetTopMargin()-0.05, ("anti-k_{t} R=" + rValStr + " jets").c_str());
    else if(pTStr.size() != 0) label_p->DrawLatex(0.7, 1.0-pad_p->GetTopMargin()-0.05, pTStr.c_str());

    label_p->DrawLatex(0.7, 1.0-pad_p->GetTopMargin()-0.10, labelAbsEta(absEtaStr).c_str());
  }

  return;
}


void plotSyst(TH1D* hist_p, std::vector<std::vector<std::vector<Double_t> > > systVect, std::vector<std::string> systStr, const Int_t centPos, std::string saveStr)
{
  TDatime* date = new TDatime();
  const std::string dateStr = std::to_string(date->GetDate());
  delete date;

  const Int_t nSyst = systStr.size();
  kirchnerPalette kPalette;
  const Int_t nColors = 6;
  const Int_t nLineStyles = 2;
  if(nSyst > nColors*nLineStyles) std::cout << "WARNING: nSyst \'" << nSyst << "\' is greater than nColors*nLineStyles=" << nColors << "*" << nLineStyles << "=" << nColors*nLineStyles << std::endl;

  TCanvas* canv_p = new TCanvas("systCanv_p", "", 450, 450);
  defineCanv(canv_p);
  canv_p->SetTopMargin(0.01);

  const Int_t nArbitraryBins = 100;
  const Int_t nBins = hist_p->GetNbinsX();
  Double_t bins[nArbitraryBins+1];
  for(Int_t bIX = 0; bIX < nBins+1; ++bIX){
    bins[bIX] = hist_p->GetXaxis()->GetBinLowEdge(bIX+1);
  }
  TH1D* systHist_p = new TH1D("systHist_h", (";" + std::string(hist_p->GetXaxis()->GetTitle()) + ";Fraction Error").c_str(), nBins, bins);

  std::vector<TH1D*> systLeg_p;
  systLeg_p.reserve(nSyst);
  
  TLegend* leg_p = new TLegend(0.25, 0.70, 0.4, 0.98);
  legStandard(leg_p);
  leg_p->SetTextSize(14);

  bool isDrawn = false;

  Double_t total[nArbitraryBins];
  Double_t totalMinusTAALumi[nArbitraryBins];

  for(Int_t aI = 0; aI < nArbitraryBins; ++aI){
    total[aI] = 0.;
    totalMinusTAALumi[aI] = 0.;
  }

  bool goodNominal = hist_p->Integral() >= TMath::Power(10,-20);
  goodNominal = goodNominal && saveStr.find("_Smooth") != std::string::npos;

  std::string histName = hist_p->GetName();
  std::string yStr = "% Error";
  if(saveStr.find("spectraRat") != std::string::npos) yStr = "Spectra Double Ratio" + yStr;
  else if(saveStr.find("SpectraRat") != std::string::npos) yStr = "Spectra Double Ratio" + yStr;
  else if(saveStr.find("SPECTRARAT") != std::string::npos) yStr = "Spectra Double Ratio" + yStr;
  else if(saveStr.find("spectra") != std::string::npos) yStr = "Spectra" + yStr;
  else if(saveStr.find("Spectra") != std::string::npos) yStr = "Spectra" + yStr;
  else if(saveStr.find("SPECTRA") != std::string::npos) yStr = "Spectra" + yStr;
  else if(saveStr.find("rraa") != std::string::npos) yStr = "RAA Ratio" + yStr;
  else if(saveStr.find("RRaa") != std::string::npos) yStr = "RAA Ratio" + yStr;
  else if(saveStr.find("RRAA") != std::string::npos) yStr = "RAA Ratio" + yStr;
  else if(saveStr.find("raa") != std::string::npos) yStr = "RAA" + yStr;
  else if(saveStr.find("Raa") != std::string::npos) yStr = "RAA" + yStr;
  else if(saveStr.find("RAA") != std::string::npos) yStr = "RAA" + yStr;    

  if(saveStr.find("_Smooth") != std::string::npos) yStr = yStr + " (w/ smoothing)";
  else yStr = yStr + " (No smoothing)";

  std::string centStr = "PP";
  std::string algoStr = histName;

  std::cout << "LINE: " << __LINE__ << std::endl;

  if(histName.find("PbPb") != std::string::npos){
    std::cout << "RUNNING " << __LINE__ << std::endl;
    centStr = histName.substr(histName.find("Cent")+4, histName.size());
    centStr = centStr.substr(0, centStr.find("_"));
    centStr.replace(centStr.find("to"), 2, "-");
    centStr = centStr + "%";

    algoStr.replace(0, algoStr.find("akCs")+4, "");
    algoStr.replace(algoStr.find("PU"), algoStr.size(), "");
    if(algoStr.size() == 2) algoStr = algoStr.substr(0,1) + "." + algoStr.substr(1,1);
    else algoStr = "0." + algoStr;

    algoStr = "R=" + algoStr;
  }
  else if(histName.find("rraa") != std::string::npos){
    std::cout << "LINE: " << __LINE__ << ", " << algoStr << std::endl;

    centStr = histName.substr(histName.find("Cent")+4, histName.size());
    centStr = centStr.substr(0, centStr.find("_"));
    algoStr = algoStr.replace(0, algoStr.find(centStr), "");
    centStr.replace(centStr.find("to"), 2, "-");
    centStr = centStr + "%";

    std::cout << "LINE: " << __LINE__ << std::endl;

    algoStr = algoStr.replace(0, algoStr.find("_")+1, "");
    algoStr = algoStr.substr(0, algoStr.find("_"));

    std::cout << "LINE: " << __LINE__ << ", " << algoStr << std::endl;

    while(algoStr.find("p") != std::string::npos){
      algoStr.replace(algoStr.find("p"), 1, ".");
    }
    std::cout << "LINE: " << __LINE__ << ", " << algoStr << std::endl;
    algoStr.replace(algoStr.find("to"), 2, " < p_{T} < ");
    std::cout << "LINE: " << __LINE__ << std::endl;

  }
  else{
    std::cout << "LINE: " << __LINE__ << std::endl;

    algoStr.replace(0, algoStr.find("ak")+2, "");
    algoStr.replace(algoStr.find("PF"), algoStr.size(), "");

    if(algoStr.size() == 2) algoStr = algoStr.substr(0,1) + "." + algoStr.substr(1,1);
    else algoStr = "0." + algoStr;
    algoStr = "R=" + algoStr;
  }

  std::string typeOfPlot = "spectraRatio";
  if(yStr.find("RAA Ratio") != std::string::npos) typeOfPlot = "raaRatio";
  else if(yStr.find("RAA") != std::string::npos) typeOfPlot = "raa";
  else if(yStr.find("Spectra%") != std::string::npos) typeOfPlot = "spectra";

  if(goodNominal) std::cout << "SYSTTOP," << typeOfPlot << "," << algoStr << "," << centStr << std::endl;
  if(goodNominal){
    std::cout << "SYST,Type,";
    for(Int_t bIX = 0; bIX < nBins; ++bIX){
      std::cout << bins[bIX] << "-" << bins[bIX+1] << ",";
    }
    std::cout << std::endl;
  }

  for(Int_t sI = 0; sI < nSyst; ++sI){
    if(systStr[sI].size() == 0) continue;
    
    if(goodNominal) std::cout << "SYST," << systStr[sI] << ",";
    for(Int_t bIX = 0; bIX < nBins; ++bIX){
      Double_t val = (systVect[centPos][sI])[bIX];
      total[bIX] = TMath::Sqrt(total[bIX]*total[bIX] + val*val);
      
      if(systStr[sI].find("Lumi") == std::string::npos){
	if(systStr[sI].find("TAA") == std::string::npos){
	  totalMinusTAALumi[bIX] = TMath::Sqrt(totalMinusTAALumi[bIX]*totalMinusTAALumi[bIX] + val*val);
	}
      }      

      if(hist_p->GetBinContent(bIX+1) > 0){
	systHist_p->SetBinContent(bIX+1, val/hist_p->GetBinContent(bIX+1));
	systHist_p->SetBinError(bIX+1, 0.);
      }
      else{
	systHist_p->SetBinContent(bIX+1, 0.0);
	systHist_p->SetBinError(bIX+1, 0.0);
      }
      
      if(goodNominal) std::cout << systHist_p->GetBinContent(bIX+1) << ",";
    }
    if(goodNominal) std::cout << std::endl;

    
    if(saveStr.find("rraa") != std::string::npos) systHist_p->SetMaximum(0.15);
    else if(saveStr.find("raa") != std::string::npos) systHist_p->SetMaximum(0.15);
    else if(saveStr.find("spectraRat") != std::string::npos) systHist_p->SetMaximum(0.07);
    else if(saveStr.find("spectra") != std::string::npos) systHist_p->SetMaximum(0.3);
    else systHist_p->SetMaximum(0.5);

    systHist_p->SetMinimum(0.0);

    systHist_p->SetMarkerSize(0.01);
    systHist_p->SetMarkerColor(kPalette.getColor(sI%nColors));
    systHist_p->SetLineColor(kPalette.getColor(sI%nColors));
    systHist_p->SetLineWidth(3);
    systHist_p->SetLineStyle(1 + sI/nColors);

    centerTitles(systHist_p);

    systHist_p->GetYaxis()->SetTitle(yStr.c_str());

    if(!isDrawn){
      systHist_p->DrawCopy("HIST");      
      isDrawn = true;
    }
    else systHist_p->DrawCopy("HIST  SAME");
  
    systLeg_p.push_back(new TH1D(("clone_" + systStr[sI]).c_str(), "", nBins, bins));
    systLeg_p[systLeg_p.size()-1]->SetMarkerSize(0.01);
    systLeg_p[systLeg_p.size()-1]->SetMarkerColor(kPalette.getColor(sI%nColors));
    systLeg_p[systLeg_p.size()-1]->SetLineColor(kPalette.getColor(sI%nColors));
    systLeg_p[systLeg_p.size()-1]->SetLineWidth(3);
    systLeg_p[systLeg_p.size()-1]->SetLineStyle(1 + sI/nColors);

    leg_p->AddEntry(systLeg_p[systLeg_p.size()-1], systStr[sI].c_str(), "L");
  }

  systLeg_p.push_back(new TH1D("systLegTotal", "", nBins, bins));
  for(Int_t bIX = 0; bIX < nBins; ++bIX){
    if(hist_p->GetBinContent(bIX+1) > 0){
      systLeg_p[systLeg_p.size()-1]->SetBinContent(bIX+1, total[bIX]/hist_p->GetBinContent(bIX+1));
      systLeg_p[systLeg_p.size()-1]->SetBinError(bIX+1, 0.);
    }
    else{
      systLeg_p[systLeg_p.size()-1]->SetBinContent(bIX+1, 0.0);
      systLeg_p[systLeg_p.size()-1]->SetBinError(bIX+1, 0.0);
    }

    systLeg_p[systLeg_p.size()-1]->SetLineStyle(1);
    systLeg_p[systLeg_p.size()-1]->SetLineColor(1);
    systLeg_p[systLeg_p.size()-1]->SetMarkerColor(1);
  }

  if(goodNominal) std::cout << "SYST,Total,";
    for(Int_t bIX = 0; bIX < nBins; ++bIX){      
      if(goodNominal) std::cout << systLeg_p[systLeg_p.size()-1]->GetBinContent(bIX+1) << ",";
    }
    if(goodNominal) std::cout << std::endl;


  systLeg_p[systLeg_p.size()-1]->DrawCopy("HIST SAME");
  leg_p->AddEntry(systLeg_p[systLeg_p.size()-1], "Total", "L");

  systLeg_p.push_back(new TH1D("systLegTotal2", "", nBins, bins));
  for(Int_t bIX = 0; bIX < nBins; ++bIX){
    if(hist_p->GetBinContent(bIX+1) > 0){
      systLeg_p[systLeg_p.size()-1]->SetBinContent(bIX+1, totalMinusTAALumi[bIX]/hist_p->GetBinContent(bIX+1));
      systLeg_p[systLeg_p.size()-1]->SetBinError(bIX+1, 0.);
    }
    else{
      systLeg_p[systLeg_p.size()-1]->SetBinContent(bIX+1, 0.0);
      systLeg_p[systLeg_p.size()-1]->SetBinError(bIX+1, 0.);
    }

    systLeg_p[systLeg_p.size()-1]->SetLineStyle(2);
    systLeg_p[systLeg_p.size()-1]->SetLineColor(1);
    systLeg_p[systLeg_p.size()-1]->SetMarkerColor(1);
  }

  canv_p->cd();

  systLeg_p[systLeg_p.size()-1]->DrawCopy("HIST SAME");
  leg_p->AddEntry(systLeg_p[systLeg_p.size()-1], "Total - Lumi/TAA", "L");
  leg_p->Draw("SAME");

  TLatex* label_p = new TLatex();
  label_p->SetNDC();
  label_p->SetTextFont(43);
  label_p->SetTextSize(22);


  std::cout << "LINE: " << __LINE__ << std::endl;

  label_p->DrawLatex(0.6, 0.9, algoStr.c_str());
  label_p->DrawLatex(0.6, 0.8, centStr.c_str());

  std::cout << "LINE: " << __LINE__ << std::endl;

  delete label_p;

  std::cout << "LINE: " << __LINE__ <<  std::endl;

  
  checkMakeDir("pdfDir");
  checkMakeDir("pdfDir/" + dateStr);
  quietSaveAs(canv_p, "pdfDir/" + dateStr + "/syst_" + saveStr + "_" + dateStr + ".pdf");
  delete canv_p;

  std::cout << "LINE: " << __LINE__ <<  std::endl;

  delete systHist_p;

  std::cout << "LINE: " << __LINE__ <<  std::endl;

  for(unsigned int sI = 0; sI < systLeg_p.size(); ++sI){
    delete systLeg_p[sI];
  }  
  systLeg_p.clear();

  std::cout << "LINE: " << __LINE__ <<  std::endl;

  delete leg_p;

  std::cout << "LINE: " << __LINE__ <<  std::endl;

  return;
}


void getNewBBPos(std::string systStr, ULong64_t &idKey, ULong64_t &idKeyBBPP, ULong64_t &idKeyBBPbPb)
{
  std::string newPosPP = std::to_string(idKeyBBPP - idKey);
  while(newPosPP.substr(newPosPP.size()-1, 1).find("0") != std::string::npos){
    newPosPP.replace(newPosPP.size()-1, 1, "");
  }
  std::string newBBPP = std::to_string(idKeyBBPP);
  while(newBBPP.substr(0, 1).find("0") == std::string::npos){
    newBBPP.replace(0, 1, "");
  }
  int newPosPP2 = std::stoi(newPosPP);
  if(systStr.find("PriorAdjUp") != std::string::npos) ++newPosPP2;
  else --newPosPP2;
  newPosPP = std::to_string(newPosPP2);
  newBBPP = newPosPP + newBBPP;
  idKeyBBPP = std::stol(newBBPP);
  
  std::string newPosPbPb = std::to_string(idKeyBBPbPb - idKey);
  while(newPosPbPb.substr(newPosPbPb.size()-1, 1).find("0") != std::string::npos){
    newPosPbPb.replace(newPosPbPb.size()-1, 1, "");
  }
  std::string newBBPbPb = std::to_string(idKeyBBPbPb);
  while(newBBPbPb.substr(0, 1).find("0") == std::string::npos){
    newBBPbPb.replace(0, 1, "");
  }
  int newPosPbPb2 = std::stoi(newPosPbPb);
  if(systStr.find("PriorAdjUp") != std::string::npos) ++newPosPbPb2;
  else --newPosPbPb2;
  newPosPbPb = std::to_string(newPosPbPb2);
  newBBPbPb = newPosPbPb + newBBPbPb;
  idKeyBBPbPb = std::stol(newBBPbPb);

  return;
}

std::string getRStrFromInt(int r)
{
  std::string rStr = std::to_string(r);
  if(rStr.size() == 2) rStr = rStr.substr(0, 1) + "." + rStr.substr(1, 1);
  else rStr = "0." + rStr;
  return rStr;
}

int getRFromJetString(std::string inStr)
{
  if(inStr.find("Cs") != std::string::npos){
    inStr.replace(0, inStr.find("Cs")+2, "");
    if(inStr.find("PU") != std::string::npos) inStr.replace(inStr.find("PU"), inStr.size(), "");
  }
  else{
    inStr.replace(0, inStr.find("ak")+2, "");
    if(inStr.find("PF") != std::string::npos) inStr.replace(inStr.find("PF"), inStr.size(), "");
  }

  return std::stoi(inStr);
}



int plotUnfoldedAll(const std::string inFileNamePP, const std::string inFileNamePbPb, const std::string inATLASFileName="", const std::string tagStr = "", const std::string overrideFile = "")
{
  if(!fileIsGood(inFileNamePP, ".root")) return 1;
  if(!fileIsGood(inFileNamePbPb, ".root")) return 1;
  if(inATLASFileName.size() != 0 && !fileIsGood(inATLASFileName, ".root")) return 1;
   
  TDatime* date = new TDatime();
  const std::string dateStr = std::to_string(date->GetDate());
  delete date;

  checkMakeDir("pdfDir");
  checkMakeDir("pdfDir/" + dateStr);

  gStyle->SetOptStat(0);

  kirchnerPalette kPalette;

  const Double_t lumiFactor = getLumiFactor();
  const Double_t nMBEvents = getNMBEvents();

  checkMakeDir("output");
  checkMakeDir("output/" + dateStr);
  std::string outFileName = "output/" + dateStr + "/plotUnfoldedAll_";
  if(tagStr.size() != 0) outFileName = outFileName + tagStr + "_";
  if(checkFile(outFileName + dateStr + ".root")) outFileName = outFileName + "UPDATED_" + dateStr + ".root";
  else outFileName = outFileName + dateStr + ".root";

  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");
  //Following two lines necessary else you get absolutely clobbered on deletion timing. See threads here:
  //https://root-forum.cern.ch/t/closing-root-files-is-dead-slow-is-this-a-normal-thing/5273/16
  //https://root-forum.cern.ch/t/tfile-speed/17549/25
  //Bizarre
  outFile_p->SetBit(TFile::kDevNull);
  TH1::AddDirectory(kFALSE);
 
  const Int_t nFiles = 2;
  TFile* inFile_p[nFiles];
  std::vector<bool> isPbPb = {true, false};
  std::vector<std::string> fileNames = {inFileNamePbPb, inFileNamePP};
  std::vector<cutPropagator*> cutProps;

  std::vector<std::string> jtAlgosPbPb;
  std::vector<std::string> jtAlgosPP;
  std::vector<std::string> dirListPbPb;
  std::vector<std::string> dirListPP;

  std::vector<std::string> histTagsPbPb;
  std::vector<std::string> histTagsPP;
  std::vector<int> histBestBayesPbPb;
  std::vector<int> histBestBayesPP;
  std::map<std::string, int> histTagsToBestBayesPbPb;
  std::map<std::string, int> histTagsToBestBayesPP;
  std::map<ULong64_t, ULong64_t> histKeyToBestBayesKeyPbPb;
  std::map<ULong64_t, ULong64_t> histKeyToBestBayesKeyPP;

  std::vector<int> histBestSvdPbPb;
  std::vector<int> histBestSvdPP;
  std::map<std::string, int> histTagsToBestSvdPbPb;
  std::map<std::string, int> histTagsToBestSvdPP;
  std::map<ULong64_t, ULong64_t> histKeyToBestSvdKeyPbPb;
  std::map<ULong64_t, ULong64_t> histKeyToBestSvdKeyPP;

  const Int_t nMaxJtAlgos = 10;
  const Int_t nMaxCentBins = 4;
  Int_t nCentBins = -1;
  std::vector<Int_t> centBinsLow;
  std::vector<Int_t> centBinsHi;
  std::vector<std::string> centBinsStr;
  std::vector<Double_t> centBinsWidth;

  Int_t nID = -1;
  std::vector<std::string> idStr;
  std::vector<bool> goodID;

  const Int_t nSmooth = 2;
  std::vector<std::string> smoothStr = {"Smooth", "NotSmooth"};
  std::vector<bool> smoothBool = {true, false};

  const Int_t nRedStat = 2;
  std::vector<std::string> redStatStr = {"RedStat", "NotRedStat"};
  std::vector<bool> redStatBool = {true, false};

  const Int_t nRatio = 2;
  std::vector<std::string> ratioStr = {"NoRatio", "Ratio"};
  std::vector<bool> ratioBool = {false, true};

  Int_t nResponseMod = -1;
  std::vector<double> responseMod;
  std::vector<std::string> responseModStr;
  std::vector<bool> goodResponseMod;

  Int_t nJtAbsEtaBins = -1; 
  std::vector<double> jtAbsEtaBinsLow;
  std::vector<double> jtAbsEtaBinsHi;
  std::vector<std::string> jtAbsEtaBinsStr;
  std::vector<bool> goodJtAbsEtaBins;
  std::vector<double> jtAbsEtaBinsWidth;

  Int_t nSystInFile = -1;
  std::vector<std::string> systStrInFile;
  Int_t nSystOutFile = 7;
  //  std::vector<std::string> systStrOutFile = {"LumiUp", "LumiDown", "TAAUp", "TAADown"};
  //  std::map<std::string, std::string> systToCombo = {{"JECUpMC", "JECMC"}, {"JECDownMC", "JECMC"}, {"JECUpData", "JECData"}, {"JECDownData", "JECData"}, {"JECUpUE", "JECUE"}, {"JECDownUE", "JECUE"}, {"PriorUp1PowerPthat", "Prior1PowerPthat"}, {"PriorDown1PowerPthat", "Prior1PowerPthat"}, {"LumiUp", "Lumi"}, {"LumiDown", "Lumi"}, {"TAAUp", "TAA"}, {"TAADown", "TAA"}};
 
  std::vector<std::string> systStrOutFile = {"SVD", "PriorAdjUp", "PriorAdjDown", "LumiUp", "LumiDown", "TAAUp", "TAADown"};
  std::map<std::string, std::string> systToCombo = {{"JECUpMC", "JECMC"}, {"JECDownMC", "JECMC"}, {"JECUpData", "JECData"}, {"JECDownData", "JECData"}, {"JECUpUE", "JECUE"}, {"JECDownUE", "JECUE"}, {"PriorUp1PowerPthat", "Prior1PowerPthat"}, {"PriorDown1PowerPthat", "Prior1PowerPthat"}, {"PriorAdjUp", "PriorAdj"}, {"PriorAdjDown", "PriorAdj"}, {"LumiUp", "Lumi"}, {"LumiDown", "Lumi"}, {"TAAUp", "TAA"}, {"TAADown", "TAA"}, {"MatrixStatRandom", "MatrixStat"}, {"MatrixStatUp", "MatrixStat"}, {"MatrixStatDown", "MatrixStat"}};

  std::vector<std::string> reducedSystStr;

  std::vector<std::string> corrSystStr = {"JECUpMC", "JECDownMC", "JECUpData", "JECDownData", "JECUpUE", "JECDownUE", "JERData"};

  std::vector<std::string> doubleCorrSystStr = {"JECUpMC", "JECDownMC", "JECUpData", "JECDownData", "JECUpUE", "JECDownUE", "JERData", "LumiUp", "LumiDown", "TAAUp", "TAADown"};
 
  Int_t nBayes = -1;

  Int_t nSVD = -1;

  std::vector<int> bayesVal;
  std::vector<std::string> bayesValStr;

  std::vector<std::string> svdValStr;

  //Use this to do any or all pt binnings
  const Int_t nMaxGeneralPtBins = 300;
  //  Int_t nGeneralPtBins;

  Int_t nRRAABins = -1;
  Double_t rraaBins[nMaxGeneralPtBins+1];
  //  Double_t generalPtBins[nMaxGeneralPtBins+1];

  std::vector<double> jtPtBins;
      
  Int_t rVals[nMaxJtAlgos];
  std::string centStrs[nMaxCentBins];
  
  Int_t nGenJtPtBinsPP[nMaxJtAlgos][nMaxCentBins];
  Double_t genJtPtBinsPP[nMaxJtAlgos][nMaxCentBins][nMaxGeneralPtBins+1];

  Int_t nGenJtPtBins[nMaxJtAlgos][nMaxCentBins];
  Double_t genJtPtBins[nMaxJtAlgos][nMaxCentBins][nMaxGeneralPtBins+1];

  Int_t nGenJtPtBinsRRAA[nMaxJtAlgos][nMaxCentBins];
  Double_t genJtPtBinsRRAA[nMaxJtAlgos][nMaxCentBins][nMaxGeneralPtBins+1];
  
  std::cout << "Processing inputs..." << std::endl;
  unsigned int pos = 0;
  for(auto const & fName : fileNames){
    inFile_p[pos] = new TFile(fName.c_str(), "READ");

    cutProps.push_back(new cutPropagator());
    cutProps[pos]->Clean();
    cutProps[pos]->GetAllVarFromFile(inFile_p[pos]);       
   
    if(isPbPb[pos]){
      dirListPbPb = returnRootFileContentsList(inFile_p[pos], "TDirectoryFile", "JetAnalyzer", 1);
      histTagsPbPb = cutProps[pos]->GetHistTagBayes();
      histBestBayesPbPb = cutProps[pos]->GetHistBestBayes();
      histBestSvdPbPb = cutProps[pos]->GetHistBestSVD();

      unsigned int sI = 0;
      while(sI < histTagsPbPb.size()){
	if(histTagsPbPb[sI].find("Flat") != std::string::npos){
	  histTagsPbPb.erase(histTagsPbPb.begin()+sI);
	  histBestBayesPbPb.erase(histBestBayesPbPb.begin()+sI);
	  histBestSvdPbPb.erase(histBestSvdPbPb.begin()+sI);
	}
	else ++sI;
      }

      nCentBins = cutProps[pos]->GetNCentBins();  
      centBinsLow = cutProps[pos]->GetCentBinsLow();
      centBinsHi = cutProps[pos]->GetCentBinsHi();
      jtAlgosPbPb = cutProps[pos]->GetJtAlgos();

      std::vector<std::string> tempDirList;
      for(unsigned int jtI = 0; jtI < jtAlgosPbPb.size(); ++jtI){
	Int_t jtPos = -1;

	for(unsigned int jtI2 = 0; jtI2 < dirListPbPb.size(); ++jtI2){
	  if(jtAlgosPbPb[jtI].find(dirListPbPb[jtI2]) != std::string::npos){
	    jtPos = jtI2;
	    break;
	  }
	}

	if(jtPos < 0){
	  std::cout << "NO JTPOS FOUND. RETURN 1" << std::endl;
	  return 1;
	}

	tempDirList.push_back(dirListPbPb[jtPos]);
      }

      dirListPbPb = tempDirList;
    }     
    else{
      dirListPP = returnRootFileContentsList(inFile_p[pos], "TDirectoryFile", "JetAnalyzer", 1);
      histTagsPP = cutProps[pos]->GetHistTagBayes();
      histBestBayesPP = cutProps[pos]->GetHistBestBayes();
      histBestSvdPP = cutProps[pos]->GetHistBestSVD();

      unsigned int sI = 0;
      while(sI < histTagsPP.size()){
	if(histTagsPP[sI].find("Flat") != std::string::npos){
	  histTagsPP.erase(histTagsPP.begin()+sI);
	  histBestBayesPP.erase(histBestBayesPP.begin()+sI);
	  histBestSvdPP.erase(histBestSvdPP.begin()+sI);
	}
	else ++sI;
      }


      jtAlgosPP = cutProps[pos]->GetJtAlgos();     

      std::vector<std::string> tempDirList;
      for(unsigned int jtI = 0; jtI < jtAlgosPP.size(); ++jtI){
	Int_t jtPos = -1;

	for(unsigned int jtI2 = 0; jtI2 < dirListPP.size(); ++jtI2){
	  if(jtAlgosPP[jtI].find(dirListPP[jtI2]) != std::string::npos){
	    jtPos = jtI2;
	    break;
	  }
	}

	if(jtPos < 0){
	  std::cout << "NO JTPOS FOUND. RETURN 1" << std::endl;
	  return 1;
	}

	tempDirList.push_back(dirListPP[jtPos]);
      }

      dirListPP = tempDirList;
    }
  
    if(pos == 0){
      nID = cutProps[pos]->GetNID();
      nResponseMod = cutProps[pos]->GetNResponseMod();
      nSystInFile = cutProps[pos]->GetNSyst();
      nJtAbsEtaBins = cutProps[pos]->GetNJtAbsEtaBins();
      nBayes = cutProps[pos]->GetNBayes();
      nSVD = cutProps[pos]->GetNSVD();

      idStr = cutProps[pos]->GetIdStr();
      responseMod = cutProps[pos]->GetResponseMod();
      jtAbsEtaBinsLow = cutProps[pos]->GetJtAbsEtaBinsLow();
      jtAbsEtaBinsHi = cutProps[pos]->GetJtAbsEtaBinsHi();      
      systStrInFile = cutProps[pos]->GetSystStr();
      
      bayesVal = cutProps[pos]->GetBayesVal();

      //      nGeneralPtBins = cutProps[pos]->GetNGeneralBins();
      std::vector<double> tempGeneralPtBins = cutProps[pos]->GetGeneralBins();
    }
    else{
      if(!compWithWarning("nID", nID, cutProps[pos]->GetNID())) return 1;
      if(!compWithWarning("nResponseMod", nResponseMod, cutProps[pos]->GetNResponseMod())) return 1;
      if(!compWithWarning("nSystInFile", nSystInFile, cutProps[pos]->GetNSyst())) return 1;
      if(!compWithWarning("nJtAbsEtaBins", nJtAbsEtaBins, cutProps[pos]->GetNJtAbsEtaBins())) return 1;

      bool propMatch = cutProps[0]->CheckPropagatorsMatch(*(cutProps[pos]), true, false, false);

      if(!propMatch) std::cout << " Mismatch in propagator 0 to " << pos << std::endl;
      else std::cout << " Good propagator " << pos << std::endl;
    }
  

    inFile_p[pos]->Close();
    delete inFile_p[pos];
    ++pos;
  }

  if(nCentBins > nMaxCentBins){
    std::cout << "nCentBins \'" << nCentBins << "\' is less than nMaxCentBins \'" << nMaxCentBins << "\'. return 1" << std::endl;
    return 1;
  }

  std::cout << __LINE__ << std::endl;
  std::vector<int> jtAlgosR;
  std::vector<int> jtAlgosROrdered;
  for(unsigned int dI = 0; dI < dirListPbPb.size(); ++dI){
    std::cout << dI << ", " << __LINE__ << std::endl;
    std::cout << "DIR: " << dirListPbPb[dI] << ", " << jtAlgosPbPb[dI] << std::endl;
    std::cout << "DIR: " << dirListPP[dI] << ", " << jtAlgosPP[dI] << std::endl;
    
    jtAlgosR.push_back(getRFromJetString(dirListPbPb[dI]));
    jtAlgosROrdered.push_back(dI);
  }


  for(unsigned int jI = 0; jI < jtAlgosROrdered.size(); ++jI){
    for(unsigned int jI2 = jI+1; jI2 < jtAlgosROrdered.size(); ++jI2){
      if(jtAlgosR[jI] > jtAlgosR[jI2]){
	int oldR = jtAlgosR[jI];
	int oldPos = jtAlgosROrdered[jI];
	
	jtAlgosR[jI] = jtAlgosR[jI2];
	jtAlgosROrdered[jI] = jtAlgosROrdered[jI2];

	jtAlgosR[jI2] = oldR;
	jtAlgosROrdered[jI2] = oldPos;
      }
    }
  }
  
  

  smallOrLargeR rReader;

  for(unsigned int jI = 0; jI < dirListPbPb.size(); ++jI){
    std::string tempStr = dirListPbPb[jI];
    if(tempStr.find("/") != std::string::npos) tempStr.replace(tempStr.find("/"), tempStr.size() - tempStr.find("/"), "");
    const Int_t rVal = getRVal(tempStr);

    rVals[jI] = rVal;
    
    for(Int_t cI = 0; cI < nCentBins; ++cI){
      std::string centStr = "Cent" + std::to_string(centBinsLow[cI]) + "to" + std::to_string(centBinsHi[cI]);	  

      centStrs[cI] = centStr;

      for(unsigned int fI = 0; fI < fileNames.size(); ++fI){
	if(!isPbPb[fI]){
	  nGenJtPtBinsPP[jI][cI] = cutProps[fI]->GetGenNBinsFromRValCent(rVal, centStr);
	  cutProps[fI]->GetGenBinsFromRValCent(rVal, centStr, genJtPtBinsPP[jI][cI]);
	}
      }
    }
  }

  for(unsigned int jI = 0; jI < dirListPbPb.size(); ++jI){
    for(Int_t cI = 0; cI < nCentBins; ++cI){

      nGenJtPtBins[jI][cI] = nGenJtPtBinsPP[jI][cI];
      for(Int_t gI = 0; gI < nGenJtPtBins[jI][cI]+1; ++gI){
	genJtPtBins[jI][cI][gI] = genJtPtBinsPP[jI][cI][gI];
      }
    }
  }

  if(overrideFile.size() != 0){
    std::cout << "RUNNING OVERRIDE" << std::endl;
    std::ifstream inFile(overrideFile.c_str());
    std::string tempStr;
    while(std::getline(inFile, tempStr)){
      if(tempStr.size() == 0) continue;
      
      tempStr = tempStr + ",";
      while(tempStr.find(",,") != std::string::npos){
	tempStr.replace(tempStr.find(",,"), 2, ",");
      }
      
      std::cout << tempStr << std::endl;

      std::vector<std::string> subStrs;
      while(tempStr.find(",") != std::string::npos){
	subStrs.push_back(tempStr.substr(0, tempStr.find(",")));
	tempStr.replace(0, tempStr.find(",")+1, "");
      }

      if(subStrs.size() == 0) continue;
      if(subStrs[0].find("#") != std::string::npos) continue;

      if(subStrs[0].find("rraa") != std::string::npos){
	for(unsigned int bIX = 1; bIX < subStrs.size(); ++bIX){
	  rraaBins[bIX-1] = std::stod(subStrs[bIX]);
	}
	nRRAABins = subStrs.size()-1;
	continue;
      }

      Int_t rPos = -1;
      Int_t centPos = -1;

      for(unsigned int rI = 0; rI < jtAlgosPbPb.size(); ++rI){
	if(std::stoi(subStrs[0]) == rVals[rI]){
	  rPos = rI;
	  break;
	}
      }

      for(int cI = 0; cI < nCentBins; ++cI){
	if(subStrs[1].find(centStrs[cI]) != std::string::npos){
	  centPos = cI;
	  break;
	}
      }

      std::vector<double> newBins;
      for(unsigned int gI = 2; gI < subStrs.size(); ++gI){
	newBins.push_back(std::stod(subStrs[gI]));
      }
    
      nGenJtPtBins[rPos][centPos] = newBins.size()-1;

      for(unsigned int gI = 0; gI < newBins.size(); ++gI){
	genJtPtBins[rPos][centPos][gI] = newBins[gI];
      }

      std::cout << rPos << ", " << centPos << std::endl;
      for(unsigned int gI = 0; gI < newBins.size(); ++gI){
	std::cout << genJtPtBins[rPos][centPos][gI] << std::endl;
      }
    }    

    inFile.close();
  }


  for(unsigned int rI = 0; rI < jtAlgosPbPb.size(); ++rI){
    std::vector<Double_t> compositeBins;

    for(Int_t cI = 0; cI < nCentBins; ++cI){
      for(Int_t gI = 0; gI < nGenJtPtBins[rI][cI]+1; ++gI){
	bool isFound = false;

	for(unsigned int dI = 0; dI < compositeBins.size(); ++dI){
	  if(TMath::Abs(genJtPtBins[rI][cI][gI] - compositeBins[dI]) > 1.) continue;

	  isFound = true;
	  break;
	}

	if(!isFound) compositeBins.push_back(genJtPtBins[rI][cI][gI]);
      }
    }

    std::sort(std::begin(compositeBins), std::end(compositeBins));

    
    nGenJtPtBinsPP[rI][0] = compositeBins.size()-1;
    for(unsigned int cI = 0; cI < compositeBins.size(); ++cI){
      genJtPtBinsPP[rI][0][cI] = compositeBins[cI];
    }
    
  }

  if(nRRAABins < 0){
    for(unsigned int jI = 0; jI < jtAlgosPbPb.size(); ++jI){
      for(Int_t cI = 0; cI < nCentBins; ++cI){
	nGenJtPtBinsRRAA[jI][cI] = nGenJtPtBins[jI][cI];

	for(Int_t bIX = 0; bIX < nGenJtPtBins[jI][cI]; ++bIX){
	  genJtPtBinsRRAA[jI][cI][bIX] = genJtPtBins[jI][cI][bIX];
	}
      }
    }
    nRRAABins = -1;
  }
  else{
    for(unsigned int jI = 0; jI < jtAlgosPbPb.size(); ++jI){
      for(Int_t cI = 0; cI < nCentBins; ++cI){
	
	nGenJtPtBinsRRAA[jI][cI] = 0;
	for(Int_t bIX = 0; bIX < nGenJtPtBins[jI][cI]+1; ++bIX){

	  bool found = false;
	  for(Int_t bIX2 = 0; bIX2 < nRRAABins; ++bIX2){
	    if(TMath::Abs(rraaBins[bIX2] - genJtPtBins[jI][cI][bIX]) > 1.) continue;
	    found = true;
	    break;
	  }
	
	  if(!found) continue;
       
	  genJtPtBinsRRAA[jI][cI][nGenJtPtBinsRRAA[jI][cI]] = genJtPtBins[jI][cI][bIX];
	  ++nGenJtPtBinsRRAA[jI][cI];
	  found = true;
	}
	
	--nGenJtPtBinsRRAA[jI][cI];
      }
    }
  }


  std::cout << "RRAABInS: " << std::endl;
  for(unsigned int jI = 0; jI < jtAlgosPbPb.size(); ++jI){
    for(Int_t cI = 0; cI < nCentBins; ++cI){
      std::cout << jI << ", " << cI << ": " << nGenJtPtBinsRRAA[jI][cI] << std::endl;
    }
  }

  std::cout << "Checking binning against rReader.." << std::endl;
  std::cout << " All bins good!" << std::endl;
  
  const Int_t nXVals = 4;
  const Double_t xVals[nXVals] = {200, 400, 600, 800};  


  /*  
  const Int_t nXValsReduced = 3;
  const Int_t xValsReduced[nXValsReduced] = {300, 400, 500};
  */

  for(unsigned int sI = 0; sI < systStrInFile.size(); ++sI){
    std::cout << " " << systStrInFile[sI] << std::endl;
    if(systStrInFile[sI].find("Flat") != std::string::npos){
      systStrInFile.erase(systStrInFile.begin()+sI);
      --nSystInFile;      
      break;
    }
  }

  for(unsigned int sI = 0; sI < systStrInFile.size(); ++sI){
    if(systToCombo.count(systStrInFile[sI]) > 0){
      std::string systTemp = systToCombo[systStrInFile[sI]];
      bool stringIsFound = false;
      for(auto const & sI2 : reducedSystStr){
	if(isStrSame(sI2, systTemp)) stringIsFound = true;
      }
      if(!stringIsFound) reducedSystStr.push_back(systToCombo[systStrInFile[sI]]);
    }
    else reducedSystStr.push_back(systStrInFile[sI]);
  }
 
  for(unsigned int sI = 0; sI < systStrOutFile.size(); ++sI){
    if(systToCombo.count(systStrOutFile[sI]) > 0){
      std::string systTemp = systToCombo[systStrOutFile[sI]];
      bool stringIsFound = false;
      for(auto const & sI2 : reducedSystStr){
	if(isStrSame(sI2, systTemp)) stringIsFound = true;
      }
      if(!stringIsFound) reducedSystStr.push_back(systToCombo[systStrOutFile[sI]]);
    }
    else reducedSystStr.push_back(systStrOutFile[sI]);
  }

  std::vector<std::string> systStrTotal = systStrInFile;
  systStrTotal.insert(std::end(systStrTotal), std::begin(systStrOutFile), std::end(systStrOutFile));

  Int_t nSystTotal = nSystInFile + nSystOutFile;
  Int_t nReducedSyst = reducedSystStr.size();

  Int_t lumiPos = -1;
  Int_t taaPos = -1;

  for(Int_t rI = 0; rI < nReducedSyst; ++rI){
    if(reducedSystStr[rI].find("Lumi") != std::string::npos) lumiPos = rI;
    else if(reducedSystStr[rI].find("TAA") != std::string::npos) taaPos = rI;
  }

  for(unsigned int hI = 0; hI < histTagsPbPb.size(); ++hI){
    if(histTagsPbPb[hI].find("nHistDim") != std::string::npos) continue;
    histTagsToBestBayesPbPb[histTagsPbPb[hI]] = histBestBayesPbPb[hI];
  }
  for(unsigned int hI = 0; hI < histTagsPbPb.size(); ++hI){
    if(histTagsPbPb[hI].find("nHistDim") != std::string::npos) continue;
    histTagsToBestSvdPbPb[histTagsPbPb[hI]] = histBestSvdPbPb[hI];
  }

  for(unsigned int hI = 0; hI < histTagsPP.size(); ++hI){
    if(histTagsPP[hI].find("nHistDim") != std::string::npos) continue;
    histTagsToBestBayesPP[histTagsPP[hI]] = histBestBayesPP[hI];
  }
  for(unsigned int hI = 0; hI < histTagsPP.size(); ++hI){
    if(histTagsPP[hI].find("nHistDim") != std::string::npos) continue;
    histTagsToBestSvdPP[histTagsPP[hI]] = histBestSvdPP[hI];

    std::cout << "PP TAG: " << histTagsPP[hI] << std::endl;
  }

  for(auto const & res : responseMod){responseModStr.push_back("ResponseMod" + prettyString(res, 2, true));}

  for(unsigned int aI = 0; aI < jtAbsEtaBinsLow.size(); ++aI){
    jtAbsEtaBinsStr.push_back("AbsEta" + prettyString(jtAbsEtaBinsLow[aI], 1, true) + "to" + prettyString(jtAbsEtaBinsHi[aI], 1, true));
  }

  for(auto const & b : bayesVal){bayesValStr.push_back("Bayes" + std::to_string(b));}

  for(Int_t sI = 0; sI < nSVD; ++sI){
    svdValStr.push_back("Svd" + std::to_string(sI+1));
  }

  for(unsigned int cI = 0; cI < centBinsLow.size(); ++cI){
    centBinsStr.push_back("Cent" + std::to_string(centBinsLow[cI]) + "to" + std::to_string(centBinsHi[cI]));
    centBinsWidth.push_back(centBinsHi[cI] - centBinsLow[cI]);
  }

  std::cout << "GENBINS: " << std::endl;
  for(unsigned int jI = 0; jI < jtAlgosPbPb.size(); ++jI){
    for(Int_t cI = 0; cI < nCentBins; ++cI){
      //      std::cout << jI << ", " << cI << std::endl;
      std::cout << jtAlgosPbPb[jI] << ", " << centBinsStr[cI] << ", " << nGenJtPtBins[jI][cI] << ", " << jI << ", " << cI << std::endl;
      for(Int_t bIX = 0; bIX < nGenJtPtBins[jI][cI]+1; ++bIX){
	std::cout << " " << genJtPtBins[jI][cI][bIX] << ",";
      }
      std::cout << std::endl;
    }
  }

  for(unsigned int aI = 0; aI < jtAbsEtaBinsLow.size(); ++aI){
    jtAbsEtaBinsStr.push_back("JtAbsEta" + std::to_string(jtAbsEtaBinsLow[aI]) + "to" + std::to_string(jtAbsEtaBinsHi[aI]));
    jtAbsEtaBinsWidth.push_back(jtAbsEtaBinsHi[aI] - jtAbsEtaBinsLow[aI]);
  }


  const UInt_t nJtAlgos = jtAlgosPbPb.size();
  /*
  Int_t r2Pos = -1;
  for(unsigned int rI = 0; rI < jtAlgosPbPb.size(); ++rI){
    if(jtAlgosPbPb[rI].find("akCs3") != std::string::npos){
      r2Pos = rI;
      break;
    }
  }
  */
  std::vector<Int_t> rValI;
  std::vector<Double_t> rValD;
  Double_t rValMin = 100;
  Double_t rValMax = -1;

  std::vector<std::string> rValStr;

  std::cout << "NBAYES: " << nBayes << std::endl;

  std::map<ULong64_t, ULong64_t> keyToVectPos;
  std::map<ULong64_t, ULong64_t> keyToVectPosSvd;
  UInt_t nKey = 0;
  UInt_t nKeySvd = 0;
  for(ULong64_t jI = 0; jI < (ULong64_t)nJtAlgos; ++jI){
    rValI.push_back(getRVal(dirListPbPb.at(jI)));
    rValD.push_back(((Double_t)rValI[jI])/10.);
    if(rValI[jI] < rValMin) rValMin = rValI[jI]; 
    if(rValI[jI] > rValMax) rValMax = rValI[jI]; 

    rValStr.push_back(getRValStr(dirListPbPb.at(jI)));

    for(ULong64_t idI = 0; idI < (ULong64_t)nID; ++idI){
      goodID.push_back(false);
      for(ULong64_t rI = 0; rI < (ULong64_t)nResponseMod; ++rI){
	goodResponseMod.push_back(false);
	for(ULong64_t aI = 0; aI < (ULong64_t)nJtAbsEtaBins; ++aI){
	  goodJtAbsEtaBins.push_back(false);
	  for(ULong64_t sI = 0; sI < (ULong64_t)nSystInFile; ++sI){

	    for(ULong64_t bI = 0; bI < (ULong64_t)nBayes; ++bI){
	      for(ULong64_t cI = 0; cI < (ULong64_t)nCentBins; ++cI){
		ULong64_t key = getKey(0, jI, cI, idI, rI, aI, sI, bayesVal[bI]);

		if(sI == 13){
		  std::cout << "KEY CHECK: " << 0 << ", " << jI << ", " << cI << ", " << idI << ", " << rI << ", " << aI << ", " << sI << ", " << bayesVal[bI] << std::endl;
		  std::cout << " KEY: " << key << ", " << nKey << std::endl;
		}
		keyToVectPos[key] = nKey;
		
		++nKey;
	      }	     
	    }

	    for(ULong64_t bI = 0; bI < (ULong64_t)nSVD; ++bI){
	      for(ULong64_t cI = 0; cI < (ULong64_t)nCentBins; ++cI){
		ULong64_t key = getKey(0, jI, cI, idI, rI, aI, sI, bI+1);
		keyToVectPosSvd[key] = nKeySvd;
		
		++nKeySvd;
	      }	     
	    }

	  }
	}
      }
    }
  }

  std::cout << "CHECK: " << keyToVectPos[704000000000] << std::endl;
  
  Int_t r2Pos = -1;
  Int_t r4Pos = -1;
  Int_t r10Pos = -1;
  std::string r2Str="";
  for(unsigned int rI = 0; rI < rValI.size(); ++rI){
    if(rValI[rI] == 2){
      r2Pos = rI;
      r2Str = std::to_string(rValI[rI]);
    }
    if(rValI[rI] == 4) r4Pos = rI;
    if(rValI[rI] == 10) r10Pos = rI;
  }
  
  rValMin /= 10.;
  rValMax /= 10.;
  rValMin -= 0.05;
  rValMax += 0.05;

  const Int_t nRBinsMax = 30;
  Double_t rBins[nRBinsMax+1];
  Int_t nRBins = (rValMax - rValMin)/0.1;
  getLinBins(rValMin, rValMax, nRBins, rBins);

  const Int_t nRValLabels = 5;
  Double_t rValLabels[nRValLabels] = {0.2, 0.4, 0.6, 0.8, 1.0};  

  std::cout << "RBins: " << std::endl;
  for(Int_t rI = 0; rI < nRBins+1; ++rI){
    std::cout << " " << rI << ": " << rBins[rI] << std::endl;
  }

  std::cout << "nKeys: " << nKey << std::endl;
  std::vector<TH1D*> histVectPbPb;
  std::vector<TH1D*> histVectPP;

  std::vector<TH1D*> histVectPbPbSvd;
  std::vector<TH1D*> histVectPPSvd;

  histVectPbPb.reserve(nKey);
  for(ULong64_t kI = 0; kI < nKey; ++kI){histVectPbPb.push_back(NULL);}

  histVectPP.reserve(nKey);
  for(ULong64_t kI = 0; kI < nKey; ++kI){histVectPP.push_back(NULL);}

  histVectPbPbSvd.reserve(nKeySvd);
  for(ULong64_t kI = 0; kI < nKeySvd; ++kI){histVectPbPbSvd.push_back(NULL);}

  histVectPPSvd.reserve(nKeySvd);
  for(ULong64_t kI = 0; kI < nKeySvd; ++kI){histVectPPSvd.push_back(NULL);}

  std::cout << "Getting all hists...." << std::endl;
  for(Int_t fI = 0; fI < nFiles; ++fI){
    std::vector<std::string> dirs;
    if(isPbPb[fI]) dirs = dirListPbPb;
    else dirs = dirListPP;

    inFile_p[fI] = new TFile(fileNames[fI].c_str(), "READ");

    for(auto const & dir : dirs){
      TDirectory* dir_p = (TDirectory*)inFile_p[fI]->Get(dir.c_str());
      TIter next(dir_p->GetListOfKeys());
      TKey* key = NULL;
   
      int algoPos = -1;
      std::vector<std::string> algos;
      if(isPbPb[fI]) algos = jtAlgosPbPb;
      else algos = jtAlgosPP;

      for(UInt_t algoI = 0; algoI < algos.size(); ++algoI){
	if(algos[algoI].find(dir) != std::string::npos){
	  algoPos = algoI;
	  break;
	}
      }
   
      while((key = (TKey*)next())){
	const std::string name = key->GetName();
	const std::string className = key->GetClassName();
	if(className.find("TH1") == std::string::npos) continue;
	
	if(name.find("Flat") != std::string::npos) continue;

	int binsPos = 0;
	int idPos = posInStrVect(name, "_", idStr, "_");
	int modPos = posInStrVect(name, "_", responseModStr, "_");
	int absEtaPos = posInStrVect(name, "_", jtAbsEtaBinsStr, "_");	int systPos = posInStrVect(name, "_", systStrInFile, "_");
	if(systPos < 0 && name.find("PriorFlat") == std::string::npos) systPos = 0;
	int bayesPos = posInStrVect(name, "_", bayesValStr, "_");
	int svdPos = posInStrVect(name, "_", svdValStr, "_");
	int centPos = posInStrVect(name, "_", centBinsStr, "_");

	//	if(name.find("Svd") != std::string::npos) continue;

	if(algoPos < 0) std::cout << __LINE__ << ": Missing algoPos in name \'" << name << "\'" << std::endl;
	if(idPos < 0) std::cout << __LINE__ << ": Missing idPos in name \'" << name << "\'" << std::endl;
	if(binsPos < 0) std::cout << __LINE__ << ": Missing binsPos in name \'" << name << "\'" << std::endl;
	if(modPos < 0) std::cout << __LINE__ << ": Missing modPos in name \'" << name << "\'" << std::endl;
	if(absEtaPos < 0) std::cout << __LINE__ << ": Missing absEtaPos in name \'" << name << "\'" << std::endl;
	if(systPos < 0) std::cout << __LINE__ << ": Missing systPos in name \'" << name << "\'" << std::endl;
	if(bayesPos < 0 && svdPos < 0) std::cout << __LINE__ << ": Missing bayesPos/svdPos in name \'" << name << "\'" << std::endl;
	if(centPos < 0) std::cout << __LINE__ << ": Missing centPos in name \'" << name << "\'" << std::endl;	

	goodID[idPos] = true;
	goodResponseMod[modPos] = true;
	goodJtAbsEtaBins[absEtaPos] = true;

	ULong64_t idKey = 0;
	ULong64_t vectPos = 0;

	if(bayesPos >= 0){
	  idKey = getKey(binsPos, algoPos, centPos, idPos, modPos, absEtaPos, systPos, bayesVal[bayesPos]);
	  vectPos = keyToVectPos[idKey];
	  if(isPbPb[fI]) histVectPbPb[vectPos] = (TH1D*)key->ReadObj();
	  else histVectPP[vectPos] = (TH1D*)key->ReadObj();	
	}
	else if(svdPos >= 0){
	  idKey = getKey(binsPos, algoPos, centPos, idPos, modPos, absEtaPos, systPos, svdPos+1);
	  vectPos = keyToVectPosSvd[idKey];
	  if(isPbPb[fI]) histVectPbPbSvd[vectPos] = (TH1D*)key->ReadObj();
	  else histVectPPSvd[vectPos] = (TH1D*)key->ReadObj();	
	}

      }
    }
  }

  std::cout << "Processing best bayes PbPb..." << std::endl;
  for(auto const & tagToBB : histTagsToBestBayesPbPb){
    int binsPos = 0;
    int algoPos = posInStrVect(tagToBB.first, "", dirListPbPb, "_");
    int idPos = posInStrVect(tagToBB.first, "_", idStr, "_");
    int modPos = posInStrVect(tagToBB.first, "_", responseModStr, "_");
    int absEtaPos = posInStrVect(tagToBB.first, "_", jtAbsEtaBinsStr, "_");
    int systPos = posInStrVect(tagToBB.first, "_", systStrInFile, "");
    if(systPos < 0 && tagToBB.first.find("PriorFlat") == std::string::npos) systPos = 0;
    int centPos = posInStrVect(tagToBB.first, "_", centBinsStr, "_");
    int bayesPos = tagToBB.second;   

    //    if(name.find("Svd") != std::string::npos) continue;
    
    if(binsPos < 0) std::cout << __LINE__ << ": Missing binsPos in tagToBB.first \'" << tagToBB.first << "\'" << std::endl;
    if(algoPos < 0) std::cout << __LINE__ << ": Missing algoPos in tagToBB.first \'" << tagToBB.first << "\'" << std::endl;
    if(idPos < 0) std::cout << __LINE__ << ": Missing idPos in tagToBB.first \'" << tagToBB.first << "\'" << std::endl;
    if(modPos < 0) std::cout << __LINE__ << ": Missing modPos in tagToBB.first \'" << tagToBB.first << "\'" << std::endl;
    if(absEtaPos < 0) std::cout << __LINE__ << ": Missing absEtaPos in tagToBB.first \'" << tagToBB.first << "\'" << std::endl;
    if(systPos < 0) std::cout << __LINE__ << ": Missing systPos in tagToBB.first \'" << tagToBB.first << "\'" << std::endl;
    if(centPos < 0) std::cout << __LINE__ << ": Missing centPos in tagToBB.first \'" << tagToBB.first << "\'" << std::endl;	
    if(bayesPos < 0){
      std::cout << __LINE__ << ": Missing bayesPos in \'" << tagToBB.first << "\', " << tagToBB.second << ". set to default 10" << std::endl;
      bayesPos = 10;
    }
  
    ULong64_t idKey = getKey(binsPos, algoPos, centPos, idPos, modPos, absEtaPos, systPos);
    ULong64_t idKeyBB = getKey(binsPos, algoPos, centPos, idPos, modPos, absEtaPos, systPos, bayesPos);  

    histKeyToBestBayesKeyPbPb[idKey] = idKeyBB;
    if(centPos == 0 && systPos == 7 && absEtaPos == 4 && idPos == 0 && algoPos == 0){
      std::cout << "FOUND THE KEY IN PBPB IDKEY, IDKEYBB " << idKey << ", " << idKeyBB << std::endl;
    }    
  }

  for(auto const & tagToBB : histTagsToBestSvdPbPb){
    int binsPos = 0;
    int algoPos = posInStrVect(tagToBB.first, "", dirListPbPb, "_");
    int idPos = posInStrVect(tagToBB.first, "_", idStr, "_");
    int modPos = posInStrVect(tagToBB.first, "_", responseModStr, "_");
    int absEtaPos = posInStrVect(tagToBB.first, "_", jtAbsEtaBinsStr, "_");
    int systPos = posInStrVect(tagToBB.first, "_", systStrInFile, "");
    if(systPos < 0 && tagToBB.first.find("PriorFlat") == std::string::npos) systPos = 0;
    int centPos = posInStrVect(tagToBB.first, "_", centBinsStr, "_");
    int svdPos = tagToBB.second;   

    //    if(name.find("Svd") != std::string::npos) continue;
    
    if(binsPos < 0) std::cout << __LINE__ << ": Missing binsPos in tagToBB.first \'" << tagToBB.first << "\'" << std::endl;
    if(algoPos < 0) std::cout << __LINE__ << ": Missing algoPos in tagToBB.first \'" << tagToBB.first << "\'" << std::endl;
    if(idPos < 0) std::cout << __LINE__ << ": Missing idPos in tagToBB.first \'" << tagToBB.first << "\'" << std::endl;
    if(modPos < 0) std::cout << __LINE__ << ": Missing modPos in tagToBB.first \'" << tagToBB.first << "\'" << std::endl;
    if(absEtaPos < 0) std::cout << __LINE__ << ": Missing absEtaPos in tagToBB.first \'" << tagToBB.first << "\'" << std::endl;
    if(systPos < 0) std::cout << __LINE__ << ": Missing systPos in tagToBB.first \'" << tagToBB.first << "\'" << std::endl;
    if(centPos < 0) std::cout << __LINE__ << ": Missing centPos in tagToBB.first \'" << tagToBB.first << "\'" << std::endl;	
    if(svdPos < 0){
      std::cout << __LINE__ << ": Missing svdPos in \'" << tagToBB.first << "\', " << tagToBB.second << ". set to default 10" << std::endl;
      svdPos = 1;
    }
  
    ULong64_t idKey = getKey(binsPos, algoPos, centPos, idPos, modPos, absEtaPos, systPos);
    ULong64_t idKeyBB = getKey(binsPos, algoPos, centPos, idPos, modPos, absEtaPos, systPos, svdPos);  

    histKeyToBestSvdKeyPbPb[idKey] = idKeyBB;
    if(centPos == 0 && systPos == 7 && absEtaPos == 4 && idPos == 0 && algoPos == 0){
      std::cout << "FOUND THE KEY IN PBPB IDKEY, IDKEYBB " << idKey << ", " << idKeyBB << std::endl;
    }    
  }

  std::cout << "Processing best bayes PP..." << std::endl;
  for(auto const & tagToBB : histTagsToBestBayesPP){
    int binsPos = 0;
    int algoPos = posInStrVect(tagToBB.first, "", dirListPP, "_");
    int idPos = posInStrVect(tagToBB.first, "_", idStr, "_");
    int modPos = posInStrVect(tagToBB.first, "_", responseModStr, "_");
    int absEtaPos = posInStrVect(tagToBB.first, "_", jtAbsEtaBinsStr, "_");
    int systPos = posInStrVect(tagToBB.first, "_", systStrInFile, "");
    if(systPos < 0 && tagToBB.first.find("PriorFlat") == std::string::npos) systPos = 0;
    int centPos = posInStrVect(tagToBB.first, "_", centBinsStr, "_");

    if(tagToBB.first.find("MatrixStat") != std::string::npos){
      std::cout << "TAG CHECK: " << tagToBB.first << ", " << systPos << std::endl;
    }

    int bayesPos = tagToBB.second;
   
    if(binsPos < 0) std::cout << __LINE__ << ": Missing binsPos in tagToBB.first \'" << tagToBB.first << "\'" << std::endl;
    if(algoPos < 0) std::cout << __LINE__ << ": Missing algoPos in tagToBB.first \'" << tagToBB.first << "\'" << std::endl;
    if(idPos < 0) std::cout << __LINE__ << ": Missing idPos in tagToBB.first \'" << tagToBB.first << "\'" << std::endl;
    if(modPos < 0) std::cout << __LINE__ << ": Missing modPos in tagToBB.first \'" << tagToBB.first << "\'" << std::endl;
    if(absEtaPos < 0) std::cout << __LINE__ << ": Missing absEtaPos in tagToBB.first \'" << tagToBB.first << "\'" << std::endl;
    if(systPos < 0) std::cout << __LINE__ << ": Missing systPos in tagToBB.first \'" << tagToBB.first << "\'" << std::endl;
    if(centPos < 0) std::cout << __LINE__ << ": Missing centPos in tagToBB.first \'" << tagToBB.first << "\'" << std::endl;	
    if(bayesPos < 0){
      std::cout << __LINE__ << ": Missing bayesPos in \'" << tagToBB.first << "\'. set to default 10" << std::endl;
      bayesPos = 10;
    }

    ULong64_t idKey = getKey(binsPos, algoPos, centPos, idPos, modPos, absEtaPos, systPos);

    ULong64_t idKeyBB = getKey(binsPos, algoPos, centPos, idPos, modPos, absEtaPos, systPos, bayesPos);
    if(300604000003000 == idKeyBB){
      std::cout << "KEY FOUND BB PP, " << idKeyBB << ", " << idKey  << std::endl;
      std::cout << "KEY FOUND PP: " << binsPos << ", " << algoPos << ", " << centPos << ", " << idPos << ", " << modPos << ", " << absEtaPos << ", " << systPos << ", " << (ULong64_t)bayesPos << std::endl;
    }
    histKeyToBestBayesKeyPP[idKey] = idKeyBB;
  }

  std::cout << "Processing best svd PP..." << std::endl;
  for(auto const & tagToBB : histTagsToBestSvdPP){
    int binsPos = 0;
    int algoPos = posInStrVect(tagToBB.first, "", dirListPP, "_");
    int idPos = posInStrVect(tagToBB.first, "_", idStr, "_");
    int modPos = posInStrVect(tagToBB.first, "_", responseModStr, "_");
    int absEtaPos = posInStrVect(tagToBB.first, "_", jtAbsEtaBinsStr, "_");
    int systPos = posInStrVect(tagToBB.first, "_", systStrInFile, "");
    if(systPos < 0 && tagToBB.first.find("PriorFlat") == std::string::npos) systPos = 0;
    int centPos = posInStrVect(tagToBB.first, "_", centBinsStr, "_");

    int svdPos = tagToBB.second;
   
    if(binsPos < 0) std::cout << __LINE__ << ": Missing binsPos in tagToBB.first \'" << tagToBB.first << "\'" << std::endl;
    if(algoPos < 0) std::cout << __LINE__ << ": Missing algoPos in tagToBB.first \'" << tagToBB.first << "\'" << std::endl;
    if(idPos < 0) std::cout << __LINE__ << ": Missing idPos in tagToBB.first \'" << tagToBB.first << "\'" << std::endl;
    if(modPos < 0) std::cout << __LINE__ << ": Missing modPos in tagToBB.first \'" << tagToBB.first << "\'" << std::endl;
    if(absEtaPos < 0) std::cout << __LINE__ << ": Missing absEtaPos in tagToBB.first \'" << tagToBB.first << "\'" << std::endl;
    if(systPos < 0) std::cout << __LINE__ << ": Missing systPos in tagToBB.first \'" << tagToBB.first << "\'" << std::endl;
    if(centPos < 0) std::cout << __LINE__ << ": Missing centPos in tagToBB.first \'" << tagToBB.first << "\'" << std::endl;	
    if(svdPos < 0){
      std::cout << __LINE__ << ": Missing svdPos in \'" << tagToBB.first << "\'. set to default 10" << std::endl;
      svdPos = 10;
    }

    ULong64_t idKey = getKey(binsPos, algoPos, centPos, idPos, modPos, absEtaPos, systPos);
    ULong64_t idKeyBB = getKey(binsPos, algoPos, centPos, idPos, modPos, absEtaPos, systPos, svdPos);

    histKeyToBestSvdKeyPP[idKey] = idKeyBB;
  }

  std::vector<int> spectraPos;
  
  std::cout << "Traditional Spectra" << std::endl;
  for(ULong64_t jI = 0; jI < (ULong64_t)nJtAlgos; ++jI){
    std::cout << "Jet Algo: " << jtAlgosPbPb[jI] << "," << jtAlgosPP[jI] << std::endl;

    for(Int_t smoothI = 0; smoothI < nSmooth; ++smoothI){
      for(ULong64_t idI = 0; idI < (ULong64_t)nID; ++idI){
	if(!goodID[idI]) continue;
	
	for(ULong64_t rI = 0; rI < (ULong64_t)nResponseMod; ++rI){
	  if(!goodResponseMod[rI]) continue;
	  
	  for(ULong64_t aI = 0; aI < (ULong64_t)nJtAbsEtaBins; ++aI){
	    if(!goodJtAbsEtaBins[aI]) continue;	
	    
	    TLatex* label_p = new TLatex();
	    labelStandard(label_p);
	    
	    TLegend* leg_p = new TLegend(0.20, 0.12, 0.4, 0.35);
	    legStandard(leg_p);
	    
	    TCanvas* canv_p = new TCanvas("canv_p", "", 450, 450);	 
	    defineCanv(canv_p);
	    
	    std::vector<TH1D*> histVectPbPbClones;
	    histVectPbPbClones.reserve(nSystTotal*nCentBins);
	    std::vector<TH1D*> histVectPPClones;	 
	    histVectPPClones.reserve(nSystTotal*nCentBins);
	    
	    //Grab the subset histograms
	    for(ULong64_t sI = 0; sI < (ULong64_t)nSystTotal; ++sI){
	      ULong64_t sIPos = sI;
	      if(sI >= (ULong64_t)nSystInFile) sIPos = 0;
	      
	      for(ULong64_t cI = 0; cI < (ULong64_t)nCentBins; ++cI){
		ULong64_t idKey = getKey(0, jI, cI, idI, rI, aI, sIPos);
		ULong64_t idKeyBBPP = histKeyToBestBayesKeyPP[idKey];
		ULong64_t idKeyBBPbPb = histKeyToBestBayesKeyPbPb[idKey];

		//                ULong64_t key = getKey(0, jI, cI, idI, rI, aI, sI, bayesVal[bI]);
				
		std::cout << systStrTotal[sI] << ", " << __LINE__ << std::endl;
		
		if(systStrTotal[sI].find("PriorAdj") != std::string::npos){
		  getNewBBPos(systStrTotal[sI], idKey, idKeyBBPP, idKeyBBPbPb);
		}
		else if(systStrTotal[sI].find("SVD") != std::string::npos){
		  idKeyBBPP = histKeyToBestSvdKeyPP[idKey];
		  idKeyBBPbPb = histKeyToBestSvdKeyPbPb[idKey];		  
		}
		/*
		1304000000000, 13/22
		  0, 0, 0, 0, 0, 4, 13
0
		  MatrixStatUp, 1794, 0
		*/
		std::cout << idKey << ", " << sI << "/" << nSystTotal << std::endl;
		std::cout << 0 << ", " << jI << ", " << cI << ", " << idI << ", " << rI << ", " << aI << ", " << sIPos << std::endl;
		std::cout << idKeyBBPP << std::endl;
		ULong64_t vectPos = keyToVectPos[idKeyBBPP];
		std::cout << systStrTotal[sI] << ", " << __LINE__ << ", " << vectPos << std::endl;
		std::string tempName = "";
		if(systStrTotal[sI].find("SVD") != std::string::npos){
		  vectPos = keyToVectPosSvd[idKeyBBPP];
		  tempName = std::string(histVectPPSvd[vectPos]->GetName()) + "_Clone_" + systStrTotal[sI];		
		}
		else tempName = std::string(histVectPP[vectPos]->GetName()) + "_Clone_" + systStrTotal[sI];		

		std::cout << systStrTotal[sI] << ", " << __LINE__ << std::endl;
		
		histVectPPClones.push_back(new TH1D(tempName.c_str(), ";Jet p_{T} [GeV];#frac{1}{#LTTAA#GT} #frac{1}{N_{evt}} #frac{d^{2}N_{Jet}}{dp_{T}d#eta}  and  #frac{d^{2}#sigma_{Jet}}{dp_{T}d#eta} [nb/GeV]", nGenJtPtBinsPP[jI][0], genJtPtBinsPP[jI][0]));	

		std::cout << systStrTotal[sI] << ", " << __LINE__ << std::endl;
		
		  
		if(systStrTotal[sI].find("SVD") != std::string::npos){
		  if(!macroHistToSubsetHist(histVectPPSvd[vectPos], histVectPPClones[histVectPPClones.size()-1])) return 1;
		}
		else{
		  if(!macroHistToSubsetHist(histVectPP[vectPos], histVectPPClones[histVectPPClones.size()-1])) return 1;
		}
		std::cout << systStrTotal[sI] << ", " << __LINE__ << std::endl;
		
		if(systStrTotal[sI].size() == 0 && cI == 0){
		  std::cout << "INIT PRINT" << std::endl;
		  histVectPP[vectPos]->Print("ALL");
		  spectraPos.push_back(vectPos);
		}
		
		vectPos = keyToVectPos[idKeyBBPbPb];
		if(systStrTotal[sI].find("SVD") != std::string::npos){
		  vectPos = keyToVectPosSvd[idKeyBBPbPb];
		  tempName = std::string(histVectPbPbSvd[vectPos]->GetName()) + "_Clone_" + systStrTotal[sI];
		}
		else tempName = std::string(histVectPbPb[vectPos]->GetName()) + "_Clone_" + systStrTotal[sI];

		std::cout << systStrTotal[sI] << ", " << __LINE__ << std::endl;
		
		histVectPbPbClones.push_back(new TH1D(tempName.c_str(), ";Jet p_{T} [GeV];#frac{1}{#LTTAA#GT} #frac{1}{N_{evt}} #frac{d^{2}N_{Jet}}{dp_{T}d#eta}  and  #frac{d^{2}#sigma_{Jet}}{dp_{T}d#eta} [nb/GeV]", nGenJtPtBins[jI][cI], genJtPtBins[jI][cI]));

		  
		if(systStrTotal[sI].find("SVD") != std::string::npos){
		  if(!macroHistToSubsetHist(histVectPbPbSvd[vectPos], histVectPbPbClones[histVectPbPbClones.size()-1])) return 1;
		}
		else{
		  if(!macroHistToSubsetHist(histVectPbPb[vectPos], histVectPbPbClones[histVectPbPbClones.size()-1])) return 1;

		}
		std::cout << systStrTotal[sI] << ", " << __LINE__ << std::endl;

	      }
	    }
	    
	    //Div by bin widths;
	    divHistByWidth(histVectPPClones);
	    divHistByWidth(histVectPbPbClones);
	    
	    std::cout << "LINE: " << __LINE__ << std::endl;
	    
	    Double_t scaleFactor = 1.;
	    for(ULong64_t cI = 0; cI < (ULong64_t)nCentBins; ++cI){
	      
	      for(ULong64_t sI = 0; sI < (ULong64_t)nSystTotal; ++sI){	 
		Double_t lumiFactorToApply = lumiFactor;		
		if(isStrSame(systStrTotal[sI], "LumiUp")) lumiFactorToApply += getLumiAbsError();
		else if(isStrSame(systStrTotal[sI], "LumiDown")) lumiFactorToApply -= getLumiAbsError();		
		scaleHist(histVectPPClones[cI + sI*nCentBins], 1./(2.*jtAbsEtaBinsWidth[aI]*lumiFactorToApply));
	      }
	      
	      scaleFactor *= 10.;
	      
	      Double_t totalFactor = nMBEvents*2.*jtAbsEtaBinsWidth[aI]*centBinsWidth[cI]/(100.*scaleFactor);
	      
	      for(ULong64_t sI = 0; sI < (ULong64_t)nSystTotal; ++sI){
		Double_t taaFactor = getTAAScaleFactorNB(centBinsStr[cI]);
		if(isStrSame(systStrTotal[sI], "TAAUp")) taaFactor += taaFactor*getTAAScaleFactorUp(centBinsStr[cI]);
		else if(isStrSame(systStrTotal[sI], "TAADown")) taaFactor -= taaFactor*getTAAScaleFactorUp(centBinsStr[cI]);		
		scaleHist(histVectPbPbClones[cI + sI*nCentBins], 1./(totalFactor*taaFactor));
	      }
	    }
	    
	    //Setup systematics
	    std::vector<std::vector<std::vector<std::vector<Double_t> > > > reducedSystPP, reducedSystPbPb;
	    std::vector<std::vector<std::vector<Double_t> > > reducedSystSumPP, reducedSystSumPbPb;
	    
	    std::vector<std::vector<Int_t> > nGenJtPtBinsVect, nGenJtPtBinsPPVect;
	    nGenJtPtBinsVect.push_back({});
	    nGenJtPtBinsPPVect.push_back({});
	    for(Int_t cI = 0; cI < nCentBins; ++cI){
	      nGenJtPtBinsVect[0].push_back(nGenJtPtBins[jI][cI]);
	      nGenJtPtBinsPPVect[0].push_back(nGenJtPtBinsPP[jI][cI]);
	    }
	    
	    //	    std::vector<std::vector<Double_t> > reducedSystPP, reducedSystPbPb, reducedSystSumPP, reducedSystSumPbPb;
	    std::cout << "NGENJTPTBINSPP: " << nGenJtPtBinsPPVect[0][0] << std::endl;
	    setupSyst(nReducedSyst, nGenJtPtBinsVect, nGenJtPtBinsPPVect, 1, nCentBins, &reducedSystPbPb, &reducedSystSumPbPb, &reducedSystPP, &reducedSystSumPP);
	    
	    if(smoothBool[smoothI]){
	      for(ULong64_t sI = 0; sI < (ULong64_t)nSystTotal; ++sI){
		if(systStrTotal[sI].size() == 0) continue;
		
		for(ULong64_t cI = 0; cI < (ULong64_t)nCentBins; ++cI){
		  if(nGenJtPtBins[jI][cI] < 3) continue;
		  
		  smoothErrors(histVectPPClones[cI + sI*nCentBins], histVectPPClones[cI]);
		  smoothErrors(histVectPbPbClones[cI + sI*nCentBins], histVectPbPbClones[cI]);
		}
	      }
	    }
	    
	    //Extract systematics
	    for(ULong64_t sI = 0; sI < (ULong64_t)nSystTotal; ++sI){
	      ULong64_t pos = 0;
	      
	      if(systToCombo.count(systStrTotal[sI]) > 0) pos = posInStrVectExact(systToCombo[systStrTotal[sI]], reducedSystStr);
	      else pos = posInStrVectExact(systStrTotal[sI], reducedSystStr);		
	      
	      for(ULong64_t cI = 0; cI < (ULong64_t)nCentBins; ++cI){
		for(Int_t bI = 0; bI < nGenJtPtBinsPP[jI][cI]; ++bI){
		  Double_t delta = TMath::Abs(histVectPPClones[cI + sI*nCentBins]->GetBinContent(bI+1) - histVectPPClones[cI]->GetBinContent(bI+1));
		  (reducedSystPP[0][cI][pos])[bI] = TMath::Min((reducedSystPP[0][cI][pos])[bI], delta);
		}
		
		for(Int_t bI = 0; bI < nGenJtPtBins[jI][cI]; ++bI){
		  Double_t delta = TMath::Abs(histVectPbPbClones[cI + sI*nCentBins]->GetBinContent(bI+1) - histVectPbPbClones[cI]->GetBinContent(bI+1));
		  
		  (reducedSystPbPb[0][cI][pos])[bI] = TMath::Min((reducedSystPbPb[0][cI][pos])[bI], delta);
		}
	      }
	    }
	    
	    for(ULong64_t cI = 0; cI < (ULong64_t)nCentBins; ++cI){
	      for(Int_t bI = 0; bI < nGenJtPtBinsPP[jI][cI]; ++bI){
		for(ULong64_t sI = 0; sI < (ULong64_t)reducedSystStr.size(); ++sI){
		  if(reducedSystStr[sI].find("PriorFlat") != std::string::npos) continue;
		  
		  (reducedSystSumPP[0][cI])[bI] = TMath::Sqrt((reducedSystSumPP[0][cI])[bI]*(reducedSystSumPP[0][cI])[bI] + (reducedSystPP[0][cI][sI])[bI]*(reducedSystPP[0][cI][sI])[bI]);
		}
	      }
	      
	      for(Int_t bI = 0; bI < nGenJtPtBins[jI][cI]; ++bI){
		for(ULong64_t sI = 0; sI < (ULong64_t)reducedSystStr.size(); ++sI){
		  if(reducedSystStr[sI].find("PriorFlat") != std::string::npos) continue;	  
		  (reducedSystSumPbPb[0][cI])[bI] = TMath::Sqrt((reducedSystSumPbPb[0][cI])[bI]*(reducedSystSumPbPb[0][cI])[bI] + (reducedSystPbPb[0][cI][sI])[bI]*(reducedSystPbPb[0][cI][sI])[bI]);
		}
	      }
	    }
	    
	    //Grab Max and Min vals...	 
	    Double_t max = 10000;
	    Double_t min = 0.000001;
	    
	    histVectPPClones[0]->SetMaximum(max);
	    histVectPPClones[0]->SetMinimum(min);
	    
	    canv_p->cd();
	    defaultPlotSet(histVectPPClones[0], "");	   
	    
	    const Int_t nBinsScale = 400;
	    Double_t binsScale[nBinsScale+1];
	    getLogBins(100, 2000, nBinsScale, binsScale);
	    Int_t lowPos = -1;
	    Int_t hiPos = -1;
	    Double_t lowestVal = 10000;
	    Double_t highestVal = -1;
	    
	    for(UInt_t jI2 = 0; jI2 < nJtAlgos; ++jI2){
	      if(lowestVal > genJtPtBinsPP[jI2][0][0]) lowestVal = genJtPtBinsPP[jI2][0][0];
	      if(highestVal < genJtPtBinsPP[jI2][0][nGenJtPtBinsPP[jI2][0]]) highestVal = genJtPtBinsPP[jI2][0][nGenJtPtBinsPP[jI2][0]];
	    }
	    
	    for(Int_t bI = 0; bI < nBinsScale; ++bI){
	      if(binsScale[bI] <= lowestVal && lowestVal < binsScale[bI+1]) lowPos = bI;
	      if(binsScale[bI] <= highestVal && highestVal < binsScale[bI+1]) hiPos = bI;
	    } 
	    
	    lowPos -= 20;
	    hiPos += 20;	    
	    
	    std::cout << "BINSRAA " << binsScale[lowPos] << ", " << binsScale[hiPos] << ", " << genJtPtBinsPP[jI][0][0] << ", " << genJtPtBinsPP[jI][0][nGenJtPtBinsPP[jI][0]]  << std::endl;
	    
	    TH1D* plotHist_p = new TH1D("plotHist_h", ";;", 10, binsScale[lowPos], binsScale[hiPos]);
	    plotHist_p->SetMaximum(max);
	    plotHist_p->SetMinimum(min);
	    defaultPlotSet(plotHist_p, "");
	    plotHist_p->GetXaxis()->SetTitle(histVectPPClones[0]->GetXaxis()->GetTitle());
	    plotHist_p->GetYaxis()->SetTitle(histVectPPClones[0]->GetYaxis()->GetTitle());
	      
	    plotHist_p->DrawCopy("HIST E1 P");
	    
	    histVectPPClones[0]->DrawCopy("HIST E1 P SAME");	 
	    canvNDCToXY labelAid(canv_p, histVectPPClones[0]);
	    gPad->SetLogx();	   
	      
	    std::string saveStr = "spectra_" + dirListPP[jI] + "_PP_" + centBinsStr[0] + "_" + idStr[idI] + "_" + responseModStr[rI] + "_" + jtAbsEtaBinsStr[aI] + "_" + smoothStr[smoothI];	     
	    
	    drawSyst(canv_p, NULL, histVectPPClones[0], &(reducedSystSumPP[0][0]), histVectPPClones[0]->GetBinLowEdge(1), histVectPPClones[0]->GetBinLowEdge(histVectPPClones[0]->GetNbinsX()+1));
	    
	    plotSyst(histVectPPClones[0], reducedSystPP[0], reducedSystStr, 0, saveStr);
	      
	    histVectPPClones[0]->DrawCopy("HIST E1 P SAME");	 
	    
	    for(ULong64_t cI = 0; cI < (ULong64_t)nCentBins; ++cI){
	      defaultPlotSet(histVectPbPbClones[cI], centBinsStr[cI]);
	      histVectPbPbClones[cI]->DrawCopy("HIST E1 P SAME");
	      drawSyst(canv_p, NULL, histVectPbPbClones[cI], &(reducedSystSumPbPb[0][cI]), histVectPbPbClones[cI]->GetBinLowEdge(1), histVectPbPbClones[cI]->GetBinLowEdge(histVectPbPbClones[cI]->GetNbinsX()+1));
	      
	      saveStr = "spectra_" + dirListPbPb[jI] + "_PbPb_" + centBinsStr[cI] + "_" + idStr[idI] + "_" + responseModStr[rI] + "_" + jtAbsEtaBinsStr[aI] + "_" + smoothStr[smoothI];
	      plotSyst(histVectPbPbClones[cI], reducedSystPbPb[0], reducedSystStr, cI, saveStr);
	      
	      histVectPbPbClones[cI]->DrawCopy("HIST E1 P SAME");
	    }
	    
	    for(ULong64_t cI = 0; cI < (ULong64_t)nCentBins; ++cI){
	      ULong64_t pos = nCentBins - 1 - cI;
	      std::string centStr = std::to_string(centBinsLow[pos]) + "-" + std::to_string(centBinsHi[pos]) + "%";
	      leg_p->AddEntry(histVectPbPbClones[pos], (centStr + " #times 10^{" + std::to_string((Int_t)pos+1) + "}").c_str(), "P L");
	    }
	    
	    leg_p->AddEntry(histVectPPClones[0], "pp #times 10^{0}", "P L");
	    
	    canv_p->cd();
	    gPad->SetLogy();
	    
	    drawAllLabels(canv_p, NULL, plotHist_p, label_p, nXVals, xVals, rValStr[jI], jtAbsEtaBinsStr[aI]);
	  	    
	    delete plotHist_p;

	    leg_p->Draw("SAME");
	    TFile* inATLASFile_p = NULL;
	    TFile* ianFile_p = NULL;
	    TH1D* ianHist_p = NULL;
	    TH1D* histPP_p = NULL;

	    if(inATLASFileName.size() != 0 && rValI[jI] == 4 && false){
	      inATLASFile_p = new TFile(inATLASFileName.c_str(), "READ");
	     
	      TH1D* tempHistPP_p = (TH1D*)inATLASFile_p->Get("Table 5/Hist1D_y1");
	      std::vector<int> atlasInts = {5, 6, 7, 8, 9};
	      std::vector<double> atlasEtas = {0.3, 0.5, 0.4, 0.4, 0.5};

	      Int_t nAtlasBins = -1;
	      Double_t atlasBins[nMaxGeneralPtBins+1];
	      Int_t startPosATLAS = 0;	      
	      
	      for(Int_t bIX = 0; bIX < tempHistPP_p->GetNbinsX()+1; ++bIX){
		if(tempHistPP_p->GetBinLowEdge(bIX+1) < genJtPtBins[jI][0][0] && tempHistPP_p->GetBinLowEdge(bIX+1) < genJtPtBins[jI][nCentBins-1][0]){
		  ++startPosATLAS;
		  continue;
		}

		atlasBins[nAtlasBins+1] = tempHistPP_p->GetBinLowEdge(bIX+1);
		++nAtlasBins;
	      }

	      for(Int_t bIX = nAtlasBins+1; bIX >= 1; --bIX){
		atlasBins[bIX] = atlasBins[bIX-1];
	      }
	      atlasBins[0] = tempHistPP_p->GetBinLowEdge(startPosATLAS);
	      --startPosATLAS;
	      ++nAtlasBins;

	      histPP_p = new TH1D("histPP_p", "", nAtlasBins, atlasBins);
	      for(unsigned int aeI = 0; aeI < atlasEtas.size(); ++aeI){
		tempHistPP_p = (TH1D*)inATLASFile_p->Get(("Table " + std::to_string(atlasInts[aeI]) + "/Hist1D_y1").c_str());

		for(Int_t bIX = startPosATLAS; bIX < tempHistPP_p->GetNbinsX(); ++bIX){
		  Double_t val1 = tempHistPP_p->GetBinContent(bIX+1)*atlasEtas[aeI];
		  Double_t val2 = histPP_p->GetBinContent(bIX+1-startPosATLAS);
		  Double_t err1 = tempHistPP_p->GetBinError(bIX+1)*atlasEtas[aeI];
		  Double_t err2 = histPP_p->GetBinError(bIX+1-startPosATLAS);
		  Double_t val = val1 + val2;
		  Double_t err = TMath::Sqrt(err1*err1 + err2*err2);

		  histPP_p->SetBinContent(bIX+1-startPosATLAS, val);
		  histPP_p->SetBinError(bIX+1-startPosATLAS, err);
		}
	      }

	      TH1D* tempHistPbPb5090_p = (TH1D*)inATLASFile_p->Get("Table 11/Hist1D_y1");
	      TH1D* histPbPb5090_p = new TH1D("histPbPb5090_h", "", nAtlasBins, atlasBins);
	      macroHistToSubsetHist(tempHistPbPb5090_p, histPbPb5090_p, true);
	      
	      histPP_p->SetMarkerSize(1);
	      histPP_p->SetMarkerStyle(24);
	      histPP_p->SetMarkerColor(1);
	      histPP_p->SetFillColor(0);
	      histPP_p->SetLineColor(1);
	      histPP_p->SetLineWidth(1);
	      for(Int_t bIX = 0; bIX < histPP_p->GetNbinsX(); ++bIX){
		histPP_p->SetBinError(bIX+1, 0);
	      }

	      for(Int_t bIX = 0; bIX < nAtlasBins; ++bIX){
		Double_t val = histPP_p->GetBinContent(bIX+1);
		Double_t err = histPP_p->GetBinError(bIX+1);

		val /= 2.1;
		err /= 2.1;

		histPP_p->SetBinContent(bIX+1, val);
		histPP_p->SetBinError(bIX+1, err);
	      }

	      scaleHist(histPP_p, 2.1/2.0);
	      std::cout << "ATLAS PRINT: " << std::endl;
	      histPP_p->Print("ALL");
	      histPP_p->DrawCopy("HIST E1 SAME");
	      histPP_p->DrawCopy("HIST E1 P SAME");
	      
	      histPbPb5090_p->SetMarkerSize(1);
	      histPbPb5090_p->SetMarkerStyle(25);
	      histPbPb5090_p->SetMarkerColor(1);
	      histPbPb5090_p->SetFillColor(0);
	      histPbPb5090_p->SetLineColor(1);
	      histPbPb5090_p->SetLineWidth(1);
	      for(Int_t bIX = 0; bIX < histPbPb5090_p->GetNbinsX(); ++bIX){
		histPbPb5090_p->SetBinError(bIX+1, 0);
	      }
	      
	      scaleFactor = 1;
	      for(Int_t scaleI = 0; scaleI < nCentBins; ++scaleI){
		scaleFactor *= 10;
	      }
	      scaleHist(histPbPb5090_p, scaleFactor*2.8/2.0);
	      histPbPb5090_p->DrawCopy("HIST E1  SAME");	    
	      histPbPb5090_p->DrawCopy("HIST E1 P SAME");	    
	      delete histPbPb5090_p;

	      leg_p->AddEntry(histPP_p, "Matched ATLAS PP, 50-90%", "L P");

	      const std::string cmsPPIanFile = "/afs/cern.ch/work/c/cmcginn/private/Projects/FullJR/Apr1to7_2019/BayesianUnfoldedData_wPY8_FULLRECO/ak4PFJets_wjtID_anabins_Bayes_PY8_FullRECO_LOMC_";
	      std::vector<std::string> cmsPPIanEtas = {"00eta05", "05eta10", "10eta15", "15eta20"};

	      Int_t nIanBins = -1;
	      Double_t ianBins[nMaxGeneralPtBins+1];
	      Int_t startPosIan = 0;

	      ianFile_p = new TFile((cmsPPIanFile + cmsPPIanEtas[0] + ".root").c_str(), "READ");
	      std::cout << "IAN: " << genJtPtBins[jI][0][0] << ", " << genJtPtBins[jI][nCentBins-1][0] << std::endl;
	      TH1D* tempHist_p = (TH1D*)ianFile_p->Get("Data_unf");
	      for(Int_t bIX = 0; bIX < tempHist_p->GetNbinsX()+1; ++bIX){
		if(tempHist_p->GetBinLowEdge(bIX+1) < genJtPtBins[jI][0][0] && tempHist_p->GetBinLowEdge(bIX+1) < genJtPtBins[jI][nCentBins-1][0]){
		  ++startPosIan;
		  continue;
		}

		ianBins[nIanBins+1] = tempHist_p->GetBinLowEdge(bIX+1);
		++nIanBins;
	      }	      
	      
	      for(Int_t bIX = nIanBins+1; bIX >= 1; --bIX){
		ianBins[bIX] = ianBins[bIX-1];
	      }
	      ianBins[0] = tempHist_p->GetBinLowEdge(startPosIan);
	      --startPosIan;
	      ++nIanBins;

	      std::cout << "IAN BINS: ";
	      for(Int_t bIX = 0; bIX < nIanBins+1; ++bIX){
		std::cout << ianBins[bIX] << ", ";
	      }
	      std::cout << std::endl;

	      ianHist_p = new TH1D("ianHist_h", ";;", nIanBins, ianBins);
	      for(Int_t bIX = 0; bIX < nIanBins; ++bIX){
		ianHist_p->SetBinContent(bIX+1, 0.0);
	      }

	      for(unsigned int eI = 0; eI < cmsPPIanEtas.size(); ++eI){
		ianFile_p->Close();
		delete ianFile_p;

		ianFile_p = new TFile((cmsPPIanFile + cmsPPIanEtas[eI] + ".root").c_str(), "READ");
		tempHist_p = (TH1D*)ianFile_p->Get("Data_unf");
		for(Int_t bIX = startPosIan; bIX < nIanBins+startPosIan; ++bIX){
		  Double_t val1 = tempHist_p->GetBinContent(bIX+1)*0.5;
		  Double_t val2 = ianHist_p->GetBinContent(bIX+1-startPosIan);
		  Double_t err1 = tempHist_p->GetBinError(bIX+1)*0.5;
		  Double_t err2 = ianHist_p->GetBinError(bIX+1-startPosIan);
		  Double_t val = val1 + val2;
		  Double_t err = TMath::Sqrt(err1*err1 + err2*err2);

		  ianHist_p->SetBinContent(bIX+1-startPosIan, val);
		  ianHist_p->SetBinError(bIX+1-startPosIan, err);
		}	       
	      }

	      for(Int_t bIX = 0; bIX < nIanBins; ++bIX){
		Double_t val = ianHist_p->GetBinContent(bIX+1);
		Double_t err = ianHist_p->GetBinError(bIX+1);

		val /= 2.0;
		err /= 2.0;

		ianHist_p->SetBinContent(bIX+1, val);
		ianHist_p->SetBinError(bIX+1, err);
	      }
	      
	      ianHist_p->SetMarkerColor(2);
	      ianHist_p->SetLineColor(2);
	      ianHist_p->SetMarkerSize(1);
	      ianHist_p->SetMarkerStyle(28);
	      ianHist_p->DrawCopy("HIST E1 P SAME");
	      leg_p->AddEntry(ianHist_p, "CMS-SMP PP", "L P");
	    }
	    
	    gPad->RedrawAxis();

	    const std::string saveName = "pdfDir/" + dateStr + "/spectra_" + dirListPbPb[jI] + "_" + idStr[idI] + "_" + responseModStr[rI] + "_" + jtAbsEtaBinsStr[aI] + "_" + smoothStr[smoothI] + "_" + dateStr+ ".pdf";
	    std::cout << "SAVE: " << saveName << std::endl;	    
	    quietSaveAs(canv_p, saveName);
	    delete canv_p;
	    delete leg_p;
	    
	    if(inATLASFile_p != NULL){
	      inATLASFile_p->Close();
	      delete inATLASFile_p;

	      delete ianHist_p;

	      ianFile_p->Close();
	      delete ianFile_p;
	    }
	    
	    delete label_p;
	    
	    outFile_p->cd();
	    
	    scaleFactor = 1.;
	    for(ULong64_t cI = 0; cI < (ULong64_t)nCentBins; ++cI){
	      scaleFactor *= 10.;
	      for(ULong64_t sI = 0; sI < (ULong64_t)nSystTotal; ++sI){
		scaleHist(histVectPbPbClones[cI + sI*nCentBins], 1./scaleFactor);
	      }
	    }
	    
	    writeAll(outFile_p, histVectPPClones);
	    writeAll(outFile_p, histVectPbPbClones);
	    deleteAll(outFile_p, &histVectPPClones);
	    deleteAll(outFile_p, &histVectPbPbClones);	   
	    histVectPPClones.clear();
	    histVectPbPbClones.clear();
	  }
	}
      }
    }
  } 
  std::cout << "End Traditional Spectra" << std::endl;
 
  unsigned int counter = 0;
  while(spectraPos.size() > counter){
    unsigned int counter2 = counter+1;
    while(spectraPos.size() > counter2){
      if(spectraPos[counter] == spectraPos[counter2]){
	spectraPos.erase(spectraPos.begin()+counter2);
      }
      else ++counter2;
    }

    std::cout << spectraPos[counter] << ": " << histVectPP[spectraPos[counter]]->GetName() << std::endl;
    ++counter;
  }
  std::cout << "COUNTER: " << counter << std::endl;


  for(UInt_t jI = 0; jI < nJtAlgos; ++jI){
    TFile* outFile2_p = new TFile(("output/" + dateStr + "/unfoldedPP_R" + std::to_string(rValI[jI]) + "_" + dateStr + ".root").c_str(), "RECREATE");
  
    for(unsigned int sI = 0; sI < spectraPos.size(); ++sI){
      std::string tempStr = histVectPP[spectraPos[sI]]->GetName();
      int rValTemp = getRFromJetString(tempStr);
     
      if(rValTemp == rValI[jI]){
	histVectPP[spectraPos[sI]]->Write("", TObject::kOverwrite);
	std::cout << "FINAL PRINT: " << std::endl;
	histVectPP[spectraPos[sI]]->Print("ALL");
      }
    }
    
    outFile2_p->Close();
    delete outFile2_p;
  }


  std::cout << "Tradtiional RAA" << std::endl;
  for(ULong64_t jI = 0; jI < (ULong64_t)nJtAlgos; ++jI){
    std::cout << "Jet Algo: " << jtAlgosPbPb[jI] << ", " << jtAlgosPP[jI] << std::endl;

    for(Int_t smoothI = 0; smoothI < nSmooth; ++smoothI){	
      for(ULong64_t idI = 0; idI < (ULong64_t)nID; ++idI){
	if(!goodID[idI]) continue;
	
	for(ULong64_t rI = 0; rI < (ULong64_t)nResponseMod; ++rI){
	  if(!goodResponseMod[rI]) continue;
	  
	  for(ULong64_t aI = 0; aI < (ULong64_t)nJtAbsEtaBins; ++aI){
	    if(!goodJtAbsEtaBins[aI]) continue;	
	    
	    TLatex* label_p = new TLatex();
	    labelStandard(label_p);
	    
	    TLegend* leg_p = new TLegend(0.25, 0.12, 0.4, 0.35);
	    legStandard(leg_p);

	    TLegend* leg2_p = new TLegend(0.45, 0.12, 0.8, 0.35);
	    legStandard(leg2_p);

	    TLegend* leg3_p = new TLegend(0.65, 0.12, 1.0, 0.35);
	    legStandard(leg3_p);
	    
	    TCanvas* canv_p = new TCanvas("canv_p", "", 450, 450);	 
	    defineCanv(canv_p);
	    
	    std::vector<TH1D*> histVectPbPbClones;
	    histVectPbPbClones.reserve(nSystTotal*nCentBins);
	    std::vector<TH1D*> histVectPPClones;	 
	    histVectPPClones.reserve(nSystTotal*nCentBins);

	    std::vector<TH1D*> histFillClones;
	    for(Int_t bI = 0; bI < nCentBins+1; ++bI){
	      histFillClones.push_back(new TH1D(("histFillClones" + std::to_string(bI)).c_str(), "", 1, 0, 10));
	      histFillClones[bI]->SetLineColor(0);
	      if(bI == 0) histFillClones[bI]->SetFillColorAlpha(kPalette.getColor(getColorPosFromCent("", true)), .25);
	      else histFillClones[bI]->SetFillColorAlpha(kPalette.getColor(getColorPosFromCent(centBinsStr[bI-1], false)), .25);
	    }
	    
	    //Grab the subset histograms
	    for(ULong64_t sI = 0; sI < (ULong64_t)nSystTotal; ++sI){
	      ULong64_t sIPos = sI;
	      if(sI >= (ULong64_t)nSystInFile) sIPos = 0;
	      
	      for(ULong64_t cI = 0; cI < (ULong64_t)nCentBins; ++cI){
		ULong64_t idKey = getKey(0, jI, cI, idI, rI, aI, sIPos);
		ULong64_t idKeyBBPP = histKeyToBestBayesKeyPP[idKey];
		ULong64_t idKeyBBPbPb = histKeyToBestBayesKeyPbPb[idKey];

		if(systStrTotal[sI].find("PriorAdj") != std::string::npos){
		  getNewBBPos(systStrTotal[sI], idKey, idKeyBBPP, idKeyBBPbPb);
		}
		else if(systStrTotal[sI].find("SVD") != std::string::npos){
		  idKeyBBPP = histKeyToBestSvdKeyPP[idKey];
		  idKeyBBPbPb = histKeyToBestSvdKeyPbPb[idKey];		  
		}


		ULong64_t vectPos = keyToVectPos[idKeyBBPP];
		std::string tempName = "";
		if(systStrTotal[sI].find("SVD") != std::string::npos){
		  vectPos = keyToVectPosSvd[idKeyBBPP];
		  tempName = "raa_" + std::string(histVectPPSvd[vectPos]->GetName()) + "_Clone_" + systStrTotal[sI];
		}
		else tempName = "raa_" + std::string(histVectPP[vectPos]->GetName()) + "_Clone_" + systStrTotal[sI];

		histVectPPClones.push_back(new TH1D(tempName.c_str(), ";Jet p_{T} [GeV];R_{AA}", nGenJtPtBins[jI][cI], genJtPtBins[jI][cI]));

		if(systStrTotal[sI].find("SVD") != std::string::npos){
		  if(!macroHistToSubsetHist(histVectPPSvd[vectPos], histVectPPClones[histVectPPClones.size()-1])) return 1;	     
		}
		else{
		  if(!macroHistToSubsetHist(histVectPP[vectPos], histVectPPClones[histVectPPClones.size()-1])) return 1;	     
		}

		vectPos = keyToVectPos[idKeyBBPbPb];

		if(systStrTotal[sI].find("SVD") != std::string::npos){
		  vectPos = keyToVectPosSvd[idKeyBBPbPb];
		  tempName = "raa_" + std::string(histVectPbPbSvd[vectPos]->GetName()) + "_Clone_" + systStrTotal[sI];
		}
		else{
		  tempName = "raa_" + std::string(histVectPbPb[vectPos]->GetName()) + "_Clone_" + systStrTotal[sI];
		}

		histVectPbPbClones.push_back(new TH1D(tempName.c_str(), ";Jet p_{T} [GeV];R_{AA}", nGenJtPtBins[jI][cI], genJtPtBins[jI][cI]));
		
		if(systStrTotal[sI].find("SVD") != std::string::npos){
		  if(!macroHistToSubsetHist(histVectPbPbSvd[vectPos], histVectPbPbClones[histVectPbPbClones.size()-1])) return 1;
		}
		else{
		  if(!macroHistToSubsetHist(histVectPbPb[vectPos], histVectPbPbClones[histVectPbPbClones.size()-1])) return 1;
		}

		
	      }
	    }
	    
	    //Div by bin widths;
	    divHistByWidth(histVectPPClones);
	    divHistByWidth(histVectPbPbClones);
	    
	    for(ULong64_t cI = 0; cI < (ULong64_t)nCentBins; ++cI){
	      for(ULong64_t sI = 0; sI < (ULong64_t)nSystTotal; ++sI){
		Double_t tempTAAFactor = getTAAScaleFactorNB(centBinsStr[cI]);	     
		if(isStrSame(systStrTotal[sI], "TAAUp")) tempTAAFactor += tempTAAFactor*getTAAScaleFactorUp(centBinsStr[cI]);
		else if(isStrSame(systStrTotal[sI], "TAADown")) tempTAAFactor -= tempTAAFactor*getTAAScaleFactorUp(centBinsStr[cI]);		
		
		Double_t totalFactor = tempTAAFactor*nMBEvents*2.*jtAbsEtaBinsWidth[aI]*centBinsWidth[cI]/100.;
		scaleHist(histVectPbPbClones[cI + sI*nCentBins], 1./totalFactor);
		
		Double_t lumiFactorToApply = lumiFactor;		
		if(isStrSame(systStrTotal[sI], "LumiUp")) lumiFactorToApply += getLumiAbsError();
		else if(isStrSame(systStrTotal[sI], "LumiDown")) lumiFactorToApply -= getLumiAbsError();		
		
		scaleHist(histVectPPClones[cI + sI*nCentBins], 1./(2.*jtAbsEtaBinsWidth[aI]*lumiFactorToApply));
		
		bool isSystCorr = false;
		for(ULong64_t sI2 = 0; sI2 < (ULong64_t)corrSystStr.size(); ++sI2){
		  if(isStrSame(corrSystStr[sI2], systStrTotal[sI])){
		    isSystCorr = true;
		    break;
		  }
		}
		
		if(isSystCorr) createRAA(histVectPbPbClones[cI + sI*nCentBins], histVectPPClones[cI + sI*nCentBins]);
		else if(systStrTotal[sI].find("Lumi") != std::string::npos) createRAA(histVectPbPbClones[cI + sI*nCentBins], histVectPPClones[cI + sI*nCentBins]);
		else createRAA(histVectPbPbClones[cI + sI*nCentBins], histVectPPClones[cI]);
		
	      }
	    }
	    
	    //Setup systematics
	    std::vector<std::vector<std::vector<std::vector<Double_t> > > > reducedSystPbPb;
	    std::vector<std::vector<std::vector<Double_t> > > reducedSystSumPbPb;

	    std::vector<std::vector<Int_t> > nGenJtPtBinsVect;
	    nGenJtPtBinsVect.push_back({});
	    for(Int_t cI = 0; cI < nCentBins; ++cI){
	      nGenJtPtBinsVect[0].push_back(nGenJtPtBins[jI][cI]);
	    }

	    setupSyst(nReducedSyst, nGenJtPtBinsVect, nGenJtPtBinsVect, 1, nCentBins, &reducedSystPbPb, &reducedSystSumPbPb);
	    
	    if(smoothBool[smoothI]){
	      for(ULong64_t sI = 0; sI < (ULong64_t)nSystTotal; ++sI){
		if(systStrTotal[sI].size() == 0) continue;
		
		for(ULong64_t cI = 0; cI < (ULong64_t)nCentBins; ++cI){
		  if(histVectPbPbClones[cI]->GetNbinsX() < 3) continue;
		  
		  smoothErrors(histVectPbPbClones[cI + sI*nCentBins], histVectPbPbClones[cI]);
		}
	      }
	    }
	    
	      
	    //Extract systematics
	    for(ULong64_t sI = 0; sI < (ULong64_t)nSystTotal; ++sI){
	      ULong64_t pos = 0;
	      if(systToCombo.count(systStrTotal[sI]) > 0) pos = posInStrVectExact(systToCombo[systStrTotal[sI]], reducedSystStr);
	      else pos = posInStrVectExact(systStrTotal[sI], reducedSystStr);
	  
	      for(ULong64_t cI = 0; cI < (ULong64_t)nCentBins; ++cI){
		for(Int_t bI = 0; bI < nGenJtPtBins[jI][cI]; ++bI){
		  Double_t delta = TMath::Abs(histVectPbPbClones[cI + sI*nCentBins]->GetBinContent(bI+1) - histVectPbPbClones[cI]->GetBinContent(bI+1));
		  
		  (reducedSystPbPb[0][cI][pos])[bI] = TMath::Min((reducedSystPbPb[0][cI][pos])[bI], delta);
		}
	      }
	    }
	    
	    for(ULong64_t cI = 0; cI < (ULong64_t)nCentBins; ++cI){
	      for(Int_t bI = 0; bI < nGenJtPtBins[jI][cI]; ++bI){	   
		for(ULong64_t sI = 0; sI < (ULong64_t)reducedSystStr.size(); ++sI){
		  if(reducedSystStr[sI].find("PriorFlat") != std::string::npos) continue;
		  else if(reducedSystStr[sI].find("TAA") != std::string::npos) continue;
		  else if(reducedSystStr[sI].find("Lumi") != std::string::npos) continue;
		  
		  (reducedSystSumPbPb[0][cI])[bI] = TMath::Sqrt((reducedSystSumPbPb[0][cI])[bI]*(reducedSystSumPbPb[0][cI])[bI] + (reducedSystPbPb[0][cI][sI])[bI]*(reducedSystPbPb[0][cI][sI])[bI]);
		}
	      }
	    }
	    
	    //Grab Max and Min vals...	 
	    Double_t min = 0.0;
	    Double_t max = 1.2;
	    histVectPbPbClones[0]->SetMaximum(max);
	    histVectPbPbClones[0]->SetMinimum(min);
	    
	    canv_p->cd();
	    defaultPlotSet(histVectPbPbClones[0], centBinsStr[0]);


	    const Int_t nBinsScale = 400;
	    Double_t binsScale[nBinsScale+1];
	    getLogBins(100, 2000, nBinsScale, binsScale);
	    Int_t lowPos = -1;
	    Int_t hiPos = -1;
	    Double_t lowestVal = 10000;
	    Double_t highestVal = -1;

	    /*
	    for(Int_t bI = 0; bI < nBinsScale; ++bI){
	      if(binsScale[bI] <= genJtPtBinsPP[jI][0][0] && genJtPtBinsPP[jI][0][0] < binsScale[bI+1]) lowPos = bI;
	      if(binsScale[bI] <= genJtPtBinsPP[jI][0][nGenJtPtBinsPP[jI][0]] && genJtPtBinsPP[jI][0][nGenJtPtBinsPP[jI][0]] < binsScale[bI+1]) hiPos = bI;
	      }*/	   

	    for(UInt_t jI2 = 0; jI2 < nJtAlgos; ++jI2){
	      if(lowestVal > genJtPtBinsPP[jI2][0][0]) lowestVal = genJtPtBinsPP[jI2][0][0];
	      if(highestVal < genJtPtBinsPP[jI2][0][nGenJtPtBinsPP[jI2][0]]) highestVal = genJtPtBinsPP[jI2][0][nGenJtPtBinsPP[jI2][0]];
	    }

	    for(Int_t bI = 0; bI < nBinsScale; ++bI){
	      if(binsScale[bI] <= lowestVal && lowestVal < binsScale[bI+1]) lowPos = bI;
	      if(binsScale[bI] <= highestVal && highestVal < binsScale[bI+1]) hiPos = bI;
	    } 

	    lowPos -= 20;
	    hiPos += 20;

	    std::cout << "BINSRAA " << binsScale[lowPos] << ", " << binsScale[hiPos] << ", " << genJtPtBinsPP[jI][0][0] << ", " << genJtPtBinsPP[jI][0][nGenJtPtBinsPP[jI][0]]  << std::endl;

	    TH1D* plotHist_p = new TH1D("plotHist_h", ";;", 10, binsScale[lowPos], binsScale[hiPos]);
	    plotHist_p->SetMaximum(max);
	    plotHist_p->SetMinimum(min);
	    defaultPlotSet(plotHist_p, "");
	    plotHist_p->GetXaxis()->SetTitle(histVectPbPbClones[0]->GetXaxis()->GetTitle());
	    plotHist_p->GetYaxis()->SetTitle(histVectPbPbClones[0]->GetYaxis()->GetTitle());
	    plotHist_p->DrawCopy("HIST E1 P");
   
	    histVectPbPbClones[0]->DrawCopy("HIST E1 P SAME");	 
	    TLine* line_p = new TLine();
	    line_p->SetLineStyle(2);
	    line_p->DrawLine(plotHist_p->GetXaxis()->GetBinLowEdge(1), 1., plotHist_p->GetXaxis()->GetBinLowEdge(plotHist_p->GetNbinsX()+1), 1.);
	    delete line_p;
	    histVectPbPbClones[0]->DrawCopy("HIST E1 P SAME");	 
	    gPad->SetLogx();
	    canvNDCToXY labelAid(canv_p, histVectPbPbClones[0]);
	    
	    drawSyst(canv_p, NULL, histVectPbPbClones[0], &(reducedSystSumPbPb[0][0]), histVectPbPbClones[0]->GetBinLowEdge(1), histVectPbPbClones[0]->GetBinLowEdge(histVectPbPbClones[0]->GetNbinsX()+1));
	    histVectPbPbClones[0]->DrawCopy("HIST E1 P SAME");
    
	    for(ULong64_t cI = 1; cI < (ULong64_t)nCentBins; ++cI){
	      defaultPlotSet(histVectPbPbClones[cI], centBinsStr[cI]);
	      if(nGenJtPtBins[jI][cI] == 0) continue;
	      
	      histVectPbPbClones[cI]->DrawCopy("HIST E1 P SAME");
	      drawSyst(canv_p, NULL, histVectPbPbClones[cI], &(reducedSystSumPbPb[0][cI]), histVectPbPbClones[cI]->GetBinLowEdge(1), histVectPbPbClones[cI]->GetBinLowEdge(histVectPbPbClones[cI]->GetNbinsX()+1));
	      histVectPbPbClones[cI]->DrawCopy("HIST E1 P SAME");
	    }
	    
	    for(ULong64_t cI = 0; cI < (ULong64_t)nCentBins; ++cI){
	      std::string saveStr = "raa_" + dirListPbPb[jI] + "_PbPb_" + centBinsStr[cI] + "_" + idStr[idI] + "_" + responseModStr[rI] + "_" + jtAbsEtaBinsStr[aI] + "_" + smoothStr[smoothI];
	      plotSyst(histVectPbPbClones[cI], reducedSystPbPb[0], reducedSystStr, cI, saveStr);
	    }
	  
	    reducedSystPbPb[0][0][lumiPos][0] /= histVectPbPbClones[0]->GetBinContent(1);
	  	    
	    drawLumiOrTAA(canv_p, histVectPbPbClones[0], plotHist_p->GetXaxis()->GetBinLowEdge(1), plotHist_p->GetXaxis()->GetBinLowEdge(plotHist_p->GetNbinsX()+1), reducedSystPbPb[0][0][lumiPos], "Lumi");

	    unsigned int tempPos = 0;
	    for(Int_t cI = nCentBins-1; cI >= 0; --cI){	    
	      if(nGenJtPtBins[jI][cI] == 0) continue;
	      
	      reducedSystPbPb[0][cI][taaPos][0] /= histVectPbPbClones[cI]->GetBinContent(1);
	      
	      drawLumiOrTAA(canv_p, histVectPbPbClones[cI], plotHist_p->GetXaxis()->GetBinLowEdge(1), plotHist_p->GetXaxis()->GetBinLowEdge(plotHist_p->GetNbinsX()+1), reducedSystPbPb[0][cI][taaPos], "TAA", tempPos);
	      ++tempPos;
	    }
	    
	    leg3_p->AddEntry(histFillClones[0], "Lumi", "F");
	    leg3_p->AddEntry(histFillClones[0], "", "");
	    leg3_p->AddEntry(histFillClones[0], "", "");
	    leg3_p->AddEntry(histFillClones[0], "", "");

	    for(ULong64_t cI = 0; cI < (ULong64_t)nCentBins; ++cI){
	      if(cI == 0) leg2_p->AddEntry(histFillClones[nCentBins-cI], "TAA", "F");
	      else leg2_p->AddEntry(histFillClones[nCentBins-cI], "", "F");
	    }

	    for(ULong64_t cI = 0; cI < (ULong64_t)nCentBins; ++cI){
	      if(nGenJtPtBins[jI][cI] == 0) continue;
	      ULong64_t pos = nCentBins - 1 - cI;
	      std::string centStr = std::to_string(centBinsLow[pos]) + "-" + std::to_string(centBinsHi[pos]) + "%";
	      leg_p->AddEntry(histVectPbPbClones[pos], centStr.c_str(), "P L");
	    }
	    
	    drawAllLabels(canv_p, NULL, plotHist_p, label_p, nXVals, xVals, rValStr[jI], jtAbsEtaBinsStr[aI]);
	    
	    delete plotHist_p;

	    leg_p->Draw("SAME");
	    leg2_p->Draw("SAME");
	    leg3_p->Draw("SAME");
	    gPad->RedrawAxis();
	    const std::string saveName = "pdfDir/" + dateStr + "/raa_" + dirListPbPb[jI] + "_" + idStr[idI] + "_" + responseModStr[rI] + "_" + jtAbsEtaBinsStr[aI] + "_" + smoothStr[smoothI] + "_" + dateStr+ ".pdf";
	    quietSaveAs(canv_p, saveName);
	    delete canv_p;
	    delete leg_p;
	    delete leg2_p;
	    delete leg3_p;
	    
	    delete label_p;
	    
	    outFile_p->cd();
	    
	    
	    writeAll(outFile_p, histVectPbPbClones);
	    deleteAll(outFile_p, &histVectPPClones);
	    deleteAll(outFile_p, &histVectPbPbClones);	   
	    histVectPPClones.clear();
	    histVectPbPbClones.clear();	   
	    deleteAll(outFile_p, &histFillClones);	   
	    histFillClones.clear();
	  }
	}
      }
    }
  }
   
  std::cout << "End Traditional RAA" << std::endl;
  
  if(r2Pos >= 0 && nRRAABins > 0 && r10Pos >= 0){    
    //Spectra RAT


    std::cout << __LINE__ << std::endl;
    for(Int_t redI = 0; redI < nRedStat; ++redI){
      for(Int_t ratI = 0; ratI < nRatio; ++ratI){
	for(Int_t smoothI = 0; smoothI < nSmooth; ++smoothI){
	  for(ULong64_t idI = 0; idI < (ULong64_t)nID; ++idI){
	    if(!goodID[idI]) continue;
	    
	    
	    for(ULong64_t rI = 0; rI < (ULong64_t)nResponseMod; ++rI){
	      if(!goodResponseMod[rI]) continue;
	    
	      std::cout << __LINE__ << std::endl;
	      
	      for(ULong64_t aI = 0; aI < (ULong64_t)nJtAbsEtaBins; ++aI){
		if(!goodJtAbsEtaBins[aI]) continue;	
		
		Int_t nCurrPtBins = 0;
		Double_t currPtBins[nMaxGeneralPtBins+1];
		std::vector<double> tempPtBins;
		
		for(Int_t jI = 0; jI < (Int_t)nJtAlgos; ++jI){
		  for(Int_t cI = 0; cI < nCentBins; ++cI){
		    
		    for(Int_t gI = 0; gI < nGenJtPtBinsPP[jI][cI]+1; ++gI){
		      bool isFound = false;
		      
		      for(unsigned int pI = 0; pI < tempPtBins.size(); ++pI){
			if(TMath::Abs(genJtPtBinsPP[jI][cI][gI] - tempPtBins[pI]) < 1){
			  isFound = true;
			  break;
			}
		      }
		      
		      if(!isFound) tempPtBins.push_back(genJtPtBinsPP[jI][cI][gI]);
		    }
		    
		  }
		}
		
		std::sort(std::begin(tempPtBins), std::end(tempPtBins));
		
		std::cout << "TEMPPTBINS: " << std::endl;
		for(unsigned int tI = 0; tI < tempPtBins.size(); ++tI){
		  std::cout << tempPtBins[tI] << std::endl;
		}
		nCurrPtBins = tempPtBins.size()-4;
		for(unsigned int pI = 1; pI < tempPtBins.size()-2; ++pI){
		  currPtBins[pI-1] = tempPtBins[pI];
		}
		currPtBins[nCurrPtBins] = tempPtBins[tempPtBins.size()-2];
		
		std::cout << "CURRPTBINS: " << std::endl;
		for(int pI = 0; pI < nCurrPtBins+1; ++pI){
		  std::cout << currPtBins[pI] << std::endl;
		}
		
		std::vector<TH1D*> histVectPPClones;	 
		histVectPPClones.reserve(nJtAlgos*nSystTotal*nCentBins);	    	  	 
		
		std::cout << __LINE__ << std::endl;
		
		for(ULong64_t sI = 0; sI < (ULong64_t)nSystTotal; ++sI){
		  ULong64_t sIPos = sI;
		  if(sI >= (ULong64_t)nSystInFile) sIPos = 0;
		  
		  std::cout << __LINE__ << std::endl;
		  for(ULong64_t cI = 0; cI < (ULong64_t)nCentBins; ++cI){
		    for(ULong64_t jI = 0; jI < (ULong64_t)nJtAlgos; ++jI){
		      
		      ULong64_t idKey = getKey(0, jI, cI, idI, rI, aI, sIPos);
		      ULong64_t idKeyBBPP = histKeyToBestBayesKeyPP[idKey];
		      std::cout << __LINE__ << std::endl;
		      ULong64_t idKeyDummy = idKeyBBPP;
		      if(systStrTotal[sI].find("PriorAdj") != std::string::npos){
			getNewBBPos(systStrTotal[sI], idKey, idKeyBBPP, idKeyDummy);
		      }
		      else if(systStrTotal[sI].find("SVD") != std::string::npos){
			idKeyBBPP = histKeyToBestSvdKeyPP[idKey];
		      }
		      std::cout << __LINE__ << std::endl;
		      
		      ULong64_t vectPos = keyToVectPos[idKeyBBPP];
		      std::string tempName = "";
		      
		      std::cout << __LINE__ << std::endl;
		      
		      if(systStrTotal[sI].find("SVD") != std::string::npos){
			std::cout << __LINE__ << std::endl;
			vectPos = keyToVectPosSvd[idKeyBBPP];		    
			std::cout << __LINE__ << std::endl;
			tempName = "raa_" + std::string(histVectPPSvd[vectPos]->GetName()) + "_Clone_" + systStrTotal[sI];
			std::cout << __LINE__ << std::endl;
		      }
		      else{
			std::cout << cI << ", " << jI << std::endl;
			std::cout << __LINE__ << ", " << vectPos << ", " << idKeyBBPP << ", " << systStrTotal[sI] << std::endl;
			
			tempName = "raa_" + std::string(histVectPP[vectPos]->GetName()) + "_Clone_" + systStrTotal[sI];
			std::cout << __LINE__ << std::endl;
		      }
		      
		      std::cout << __LINE__ << std::endl;
		      
		      histVectPPClones.push_back(new TH1D(tempName.c_str(), ";Jet p_{T} [GeV];#frac{d^{2}N^{R}}{dp_{T}d#eta}/#frac{d^{2}N^{R=1.0}}{dp_{T}d#eta}", nCurrPtBins, currPtBins));
		      
		      if(systStrTotal[sI].find("SVD") != std::string::npos){
			if(!macroHistToSubsetHist(histVectPPSvd[vectPos], histVectPPClones[histVectPPClones.size()-1])) return 1;	     
		      }
		      else{
			if(!macroHistToSubsetHist(histVectPP[vectPos], histVectPPClones[histVectPPClones.size()-1])) return 1;	     
			
			if(sI == 0){
			  std::cout << "INITIAL CHECK: " << std::endl;
			  histVectPPClones[histVectPPClones.size()-1]->Print("ALL");
			}
		      }
		    }
		  }
		}
		
		//Div by bin widths;
		divHistByWidth(histVectPPClones);
		
		for(ULong64_t jI = 0; jI < (ULong64_t)nJtAlgos; ++jI){
		  for(ULong64_t cI = 0; cI < (ULong64_t)nCentBins; ++cI){
		    for(ULong64_t sI = 0; sI < (ULong64_t)nSystTotal; ++sI){
		      
		      Double_t lumiFactorToApply = lumiFactor;		
		      if(isStrSame(systStrTotal[sI], "LumiUp")) lumiFactorToApply += getLumiAbsError();
		      else if(isStrSame(systStrTotal[sI], "LumiDown")) lumiFactorToApply -= getLumiAbsError();		
		      
		      scaleHist(histVectPPClones[jI + cI*nJtAlgos + sI*nCentBins*nJtAlgos], 1./(2.*jtAbsEtaBinsWidth[aI]*lumiFactorToApply));
		    }
		  }
		}
		
		for(ULong64_t jI = 0; jI < (ULong64_t)nJtAlgos; ++jI){
		  for(ULong64_t cI = 0; cI < (ULong64_t)nCentBins; ++cI){
		    for(ULong64_t sI = 0; sI < (ULong64_t)nSystTotal; ++sI){
		      bool isSystDoubleCorr = false;
		      for(ULong64_t sI2 = 0; sI2 < (ULong64_t)doubleCorrSystStr.size(); ++sI2){
			if(isStrSame(doubleCorrSystStr[sI2], systStrTotal[sI])){
			  isSystDoubleCorr = true;
			  break;
			}
		      }
		      
		      if(jI != (ULong64_t)r10Pos){
			if(sI == 0){
			  std::cout << "BEFORE RAA: " << std::endl;
			  histVectPPClones[jI + cI*nJtAlgos + sI*nCentBins*nJtAlgos]->Print("ALL");
			}
			
			std::cout << "DOUBLE CORR: " << systStrTotal[sI] << ", " << isSystDoubleCorr << std::endl;
			
			if(isSystDoubleCorr) createRAA(histVectPPClones[jI + cI*nJtAlgos + sI*nCentBins*nJtAlgos], histVectPPClones[r10Pos + cI*nJtAlgos + sI*nCentBins*nJtAlgos]);
			else createRAA(histVectPPClones[jI + cI*nJtAlgos + sI*nCentBins*nJtAlgos], histVectPPClones[r10Pos + cI*nJtAlgos]);
			
			if(isStrSame(redStatStr[redI], "RedStat")){
			  for(Int_t bIX = 0; bIX < histVectPPClones[jI + cI*nJtAlgos + sI*nCentBins*nJtAlgos]->GetXaxis()->GetNbins(); ++bIX){
			    Double_t tempErr = histVectPPClones[jI + cI*nJtAlgos + sI*nCentBins*nJtAlgos]->GetBinError(bIX+1);
			    
			    if(rValI[jI] == 8) tempErr *= 0.29;
			    else if(rValI[jI] == 6) tempErr *= 0.4425;
			    else if(rValI[jI] == 4) tempErr *= 0.5625;
			    else if(rValI[jI] == 3) tempErr *= 0.605;
			    else if(rValI[jI] == 2) tempErr *= 0.703;
			    histVectPPClones[jI + cI*nJtAlgos + sI*nCentBins*nJtAlgos]->SetBinError(bIX+1, tempErr);		      
			  }
			}
			if(sI == 0){
			  std::cout << "AFTER RAA: " << std::endl;
			  histVectPPClones[jI + cI*nJtAlgos + sI*nCentBins*nJtAlgos]->Print("ALL");
			}
			
		      }
		    }
		  }
		}	          	 
		
		std::cout << __LINE__ << std::endl;
		std::vector<std::vector<std::vector<std::vector<Double_t> > > > reducedSystPP;
		std::vector<std::vector<std::vector<Double_t> > > reducedSystSumPP;
		
		std::vector<std::vector<Int_t> > nGenJtPtBinsVect;
		for(ULong64_t jI = 0; jI < (ULong64_t)nJtAlgos; ++jI){
		  nGenJtPtBinsVect.push_back({});	    
		  for(Int_t cI = 0; cI < nCentBins; ++cI){
		    nGenJtPtBinsVect[jI].push_back(nCurrPtBins);
		  }
		}
		
		std::cout << __LINE__ << std::endl;
		setupSyst(nReducedSyst, nGenJtPtBinsVect, nGenJtPtBinsVect, nJtAlgos, nCentBins, &reducedSystPP, &reducedSystSumPP);
		
		std::cout << __LINE__ << std::endl;
		
		if(smoothBool[smoothI]){
		  for(ULong64_t sI = 0; sI < (ULong64_t)nSystTotal; ++sI){
		    if(systStrTotal[sI].size() == 0) continue;
		    
		    for(ULong64_t cI = 0; cI < (ULong64_t)nCentBins; ++cI){
		      for(ULong64_t jI = 0; jI < (ULong64_t)nJtAlgos; ++jI){
			
			if(histVectPPClones[jI + cI*nJtAlgos]->GetNbinsX() < 3) continue;
			smoothErrors(histVectPPClones[jI + cI*nJtAlgos + sI*nCentBins*nJtAlgos], histVectPPClones[jI + cI*nJtAlgos]);
		      }
		    }
		  }
		}
		
		std::cout << __LINE__ << std::endl;
		
		//Extract systematics
		for(ULong64_t sI = 0; sI < (ULong64_t)nSystTotal; ++sI){
		  
		  ULong64_t pos = 0;
		  if(systToCombo.count(systStrTotal[sI]) > 0) pos = posInStrVectExact(systToCombo[systStrTotal[sI]], reducedSystStr);
		  else pos = posInStrVectExact(systStrTotal[sI], reducedSystStr);
		  
		  for(ULong64_t jI = 0; jI < (ULong64_t)nJtAlgos; ++jI){
		    for(ULong64_t cI = 0; cI < (ULong64_t)nCentBins; ++cI){
		      
		      
		      for(Int_t bI = 0; bI < nCurrPtBins; ++bI){
			Double_t delta = TMath::Abs(histVectPPClones[jI + cI*nJtAlgos + sI*nCentBins*nJtAlgos]->GetBinContent(bI+1) - histVectPPClones[jI + cI*nJtAlgos]->GetBinContent(bI+1));
		      
			(reducedSystPP[jI][cI][pos])[bI] = TMath::Min((reducedSystPP[jI][cI][pos])[bI], delta);
		      }		  
		    }
		  }
		}
		
		std::cout << __LINE__ << std::endl;
		
		for(ULong64_t jI = 0; jI < (ULong64_t)nJtAlgos; ++jI){
		  for(ULong64_t cI = 0; cI < (ULong64_t)nCentBins; ++cI){
		    for(Int_t bI = 0; bI < nCurrPtBins; ++bI){
		      for(ULong64_t sI = 0; sI < (ULong64_t)reducedSystStr.size(); ++sI){
			(reducedSystSumPP[jI][cI])[bI] = TMath::Sqrt((reducedSystSumPP[jI][cI])[bI]*(reducedSystSumPP[jI][cI])[bI] + (reducedSystPP[jI][cI][sI])[bI]*(reducedSystPP[jI][cI][sI])[bI]);
		      }
		    }		
		  }
		}
		
		TLatex* label_p = new TLatex();
		labelStandard(label_p);
		
		TLegend* leg_p = NULL;
		if(ratioBool[ratI]) leg_p = new TLegend(0.25, 0.02, 0.4, 0.25);
		else leg_p = new TLegend(0.25, 0.12, 0.4, 0.35);
		legStandard(leg_p);
		
		TLegend* leg2_p = NULL;
		if(ratioBool[ratI]) leg2_p = new TLegend(0.45, 0.02, 0.6, 0.25);
		else leg2_p = new TLegend(0.45, 0.12, 0.6, 0.35);
		legStandard(leg2_p);
		
		TLegend* leg3_p = NULL;
		if(ratioBool[ratI]) leg3_p = new TLegend(0.65, 0.02, 0.8, 0.25);
		else leg3_p = new TLegend(0.65, 0.12, 0.8, 0.35);
		legStandard(leg3_p);
		
		TCanvas* canv_p = new TCanvas("canv_p", "", 450, 450);	 
		defineCanv(canv_p);		

		const Int_t nPads = 2;
		TPad* pads_p[nPads] = {NULL, NULL};
		Double_t ySplit = 0.35;

		std::cout << "LINE: " << __LINE__ << std::endl;
		if(ratioBool[ratI]){
		  Double_t origTopMargin = canv_p->GetTopMargin();
		  Double_t origLeftMargin = canv_p->GetLeftMargin();
		  Double_t origRightMargin = canv_p->GetRightMargin();
		  Double_t origBottomMargin = canv_p->GetBottomMargin();

		  canv_p->SetTopMargin(0.01);
		  canv_p->SetRightMargin(0.01);
		  canv_p->SetLeftMargin(0.01);
		  canv_p->SetBottomMargin(0.01);

		  std::cout << "ORIG MARG: " << origTopMargin << ", " << origBottomMargin << ", " << origLeftMargin << ", " << origRightMargin << std::endl;

		  canv_p->cd();
		  pads_p[0] = new TPad("pad0", "", 0.0, ySplit, 1.0, 1.0);
		  pads_p[0]->SetRightMargin(origRightMargin);
		  pads_p[0]->SetTopMargin(origTopMargin/(1.-ySplit));
		  pads_p[0]->SetLeftMargin(origLeftMargin);
		  pads_p[0]->SetBottomMargin(0.001);
		  pads_p[0]->Draw("SAME");

		  canv_p->cd();
		  pads_p[1] = new TPad("pad1", "", 0.0, 0.0, 1.0, ySplit);
		  pads_p[1]->SetRightMargin(origRightMargin);
		  pads_p[1]->SetTopMargin(.001);
		  pads_p[1]->SetLeftMargin(origLeftMargin);
		  pads_p[1]->SetBottomMargin(origBottomMargin/(ySplit));
		  pads_p[1]->Draw("SAME");

		  pads_p[0]->cd();
		}

		std::cout << "LINE: " << __LINE__ << std::endl;
		
		bool isDrawn = false;
		for(ULong64_t cI = 0; cI < (ULong64_t)1; ++cI){		
		  for(ULong64_t sI = 0; sI < (ULong64_t)nSystTotal; ++sI){
		    for(ULong64_t jI = 0; jI < (ULong64_t)nJtAlgos; ++jI){
		      if(r10Pos == (Int_t)jI) continue;
		      //canv_p->cd();
		    		    
		      defaultPlotSetR(histVectPPClones[jI + cI*nJtAlgos + sI*nJtAlgos*nCentBins], jtAlgosPP[jI]);		
		      if(sI == 0){
			if(!isDrawn){
			  if(ratioBool[ratI]){
			    histVectPPClones[jI + cI*nJtAlgos + sI*nJtAlgos*nCentBins]->GetXaxis()->SetTitleColor(0);
			    histVectPPClones[jI + cI*nJtAlgos + sI*nJtAlgos*nCentBins]->GetXaxis()->SetLabelColor(0);
			  }

			  histVectPPClones[jI + cI*nJtAlgos + sI*nJtAlgos*nCentBins]->SetMinimum(0.0);
			  histVectPPClones[jI + cI*nJtAlgos + sI*nJtAlgos*nCentBins]->SetMaximum(1.2);

			  histVectPPClones[jI + cI*nJtAlgos + sI*nJtAlgos*nCentBins]->GetXaxis()->SetTickSize(0.03/(1.-ySplit));


//			  for(Int_t bIX = 0; bIX < histVectPPClones[jI + cI*nJtAlgos + sI*nJtAlgos*nCentBins]->GetNbinsX(); ++bIX){
//			    histVectPPClones[jI + cI*nJtAlgos + sI*nJtAlgos*nCentBins]->SetBinContent(bIX+1, -100);
//			    histVectPPClones[jI + cI*nJtAlgos + sI*nJtAlgos*nCentBins]->SetBinError(bIX+1, 0.0);
//
//			  }

			  histVectPPClones[jI + cI*nJtAlgos + sI*nJtAlgos*nCentBins]->DrawCopy("HIST E1 P");
			  isDrawn = true;
			}
			else histVectPPClones[jI + cI*nJtAlgos + sI*nJtAlgos*nCentBins]->DrawCopy("HIST E1 P SAME");		      		      
		      }		    
		    }
		  }

		std::cout << "LINE: " << __LINE__ << std::endl;
		
		  for(ULong64_t jI = 0; jI < (ULong64_t)nJtAlgos; ++jI){
		    int orderPos = jtAlgosROrdered[jI];
		    if(orderPos == r10Pos) continue;
		    
		    int rVal = getRFromJetString(jtAlgosPbPb[orderPos]);
		    std::string rStr = getRStrFromInt(rVal);
		    
		    leg_p->AddEntry(histVectPPClones[orderPos], ("R="+rStr).c_str(), "P L");		
		  }
		
		  for(ULong64_t jI = 0; jI < (ULong64_t)nJtAlgos; ++jI){
		    if(r10Pos == (Int_t)jI) continue;
		  
		    std::cout << jtAlgosPP[jI] << std::endl;
		    
		    drawSyst(canv_p, pads_p[0], histVectPPClones[jI], &(reducedSystSumPP[jI][0]), histVectPPClones[jI]->GetBinLowEdge(1), histVectPPClones[jI]->GetBinLowEdge(histVectPPClones[jI]->GetNbinsX()+1));
		    
		    std::string saveStr = "spectraRat_" + dirListPP[jI] + "_PP_" + centBinsStr[0] + "_" + idStr[idI] + "_" + responseModStr[rI] + "_" + jtAbsEtaBinsStr[aI] + "_" + smoothStr[smoothI];	     
		    
		    plotSyst(histVectPPClones[jI], reducedSystPP[jI], reducedSystStr, 0, saveStr);
		  }
		
		  std::cout << "LINE: " << __LINE__ << std::endl;

		  for(ULong64_t jI = 0; jI < (ULong64_t)nJtAlgos; ++jI){
		    if(r10Pos == (Int_t)jI) continue;
		    //		    canv_p->cd();

		    
		    histVectPPClones[jI]->Print("ALL");
		    histVectPPClones[jI]->DrawCopy("HIST E1 P SAME");
		  }
		  
		  
		  std::vector<TH1D*> theoryCompPYT6;
		  std::vector<TH1D*> theoryCompPYT8;
		  for(ULong64_t jI = 0; jI < (ULong64_t)nJtAlgos; ++jI){
		    theoryCompPYT6.push_back(new TH1D(("theoryCompPYT6_" + std::to_string(jI)).c_str(), "", nCurrPtBins, currPtBins));
		  
		    TFile* theoryFile_p = new TFile(("/data/cmcginn/FullJR/Response/20190502/20190502/combinedResponse_ak" + std::to_string(rValI[jI]) + "_20190502.root").c_str(), "READ");
		    TH1D* tempTheory_p = (TH1D*)theoryFile_p->Get(("ak" + std::to_string(rValI[jI])+ "PFJetAnalyzer/genJtPt_ak" + std::to_string(rValI[jI]) + "PFJetAnalyzer_PP_Cent0to10_LightMUAndCHID_ResponseMod0p00_AbsEta0p0to2p0_General_h").c_str());
		    
		    macroHistToSubsetHist(tempTheory_p, theoryCompPYT6[theoryCompPYT6.size()-1], true);
		    theoryFile_p->Close();
		    delete theoryFile_p;
		  }
		  
		  for(ULong64_t jI = 0; jI < (ULong64_t)nJtAlgos; ++jI){
		    theoryCompPYT8.push_back(new TH1D(("theoryCompPYT8_" + std::to_string(jI)).c_str(), "", nCurrPtBins, currPtBins));
		    
		    TFile* theoryFile_p = new TFile(("/data/cmcginn/FullJR/Response/20190603/combinedResponse_ak" + std::to_string(rValI[jI]) + "_PYTHIA8_20190603.root").c_str(), "READ");
		    TH1D* tempTheory_p = (TH1D*)theoryFile_p->Get(("ak" + std::to_string(rValI[jI])+ "PFJetAnalyzer/genJtPt_ak" + std::to_string(rValI[jI]) + "PFJetAnalyzer_PP_Cent0to10_LightMUAndCHID_ResponseMod0p00_AbsEta0p0to2p0_General_h").c_str());
		    
		    macroHistToSubsetHist(tempTheory_p, theoryCompPYT8[theoryCompPYT8.size()-1], true);
		    theoryFile_p->Close();
		    delete theoryFile_p;
		  }

		  if(ratioBool[ratI]) pads_p[0]->cd();
		  for(ULong64_t jI = 0; jI < (ULong64_t)nJtAlgos; ++jI){
		    if(jI == (ULong64_t)r10Pos) continue;
		    
		    createRAA(theoryCompPYT6[jI], theoryCompPYT6[r10Pos]);
		    defaultPlotSetR(theoryCompPYT6[jI], jtAlgosPP[jI]);		
		    theoryCompPYT6[jI]->SetLineStyle(2);
		    theoryCompPYT6[jI]->SetLineWidth(3);
		    theoryCompPYT6[jI]->DrawCopy("HIST C SAME");
		  }
	       
			  
		  if(ratioBool[ratI]){
		    pads_p[1]->cd();

		    bool isDrawn = false;
		    for(ULong64_t jI = 0; jI < (ULong64_t)nJtAlgos; ++jI){
		      if(jI == (ULong64_t)r10Pos) continue;
		      if(jI != (ULong64_t)r4Pos && jI != (ULong64_t)r2Pos) continue;

		      for(Int_t bIX = 0; bIX < theoryCompPYT6[jI]->GetNbinsX(); ++bIX){
			Double_t val = theoryCompPYT6[jI]->GetBinContent(bIX+1)/histVectPPClones[jI]->GetBinContent(bIX+1);
			Double_t relErr = val*histVectPPClones[jI]->GetBinError(bIX+1)/histVectPPClones[jI]->GetBinContent(bIX+1);

			theoryCompPYT6[jI]->SetBinContent(bIX+1, val);
			theoryCompPYT6[jI]->SetBinError(bIX+1, relErr);
		      }		      


		      theoryCompPYT6[jI]->GetYaxis()->SetNdivisions(404);
		      theoryCompPYT6[jI]->GetYaxis()->SetTitle("#frac{Theory}{Data} R=0.2,0.4");
		      theoryCompPYT6[jI]->GetXaxis()->SetTitle("Jet p_{T} [GeV]");
		      theoryCompPYT6[jI]->GetXaxis()->SetTitleOffset(3.25);
		      centerTitles(theoryCompPYT6[jI]);
		      theoryCompPYT6[jI]->SetMaximum(1.25);
		      theoryCompPYT6[jI]->SetMinimum(0.75);

		      if(!isDrawn){
			isDrawn = true;
			pads_p[1]->SetLogx();

			theoryCompPYT6[jI]->GetXaxis()->SetTickSize(0.03/ySplit);
			theoryCompPYT6[jI]->GetYaxis()->SetTitleSize(theoryCompPYT6[jI]->GetYaxis()->GetTitleSize()-2);
			theoryCompPYT6[jI]->DrawCopy("HIST C");

			drawXLabels(pads_p[1], theoryCompPYT6[jI], label_p, nXVals, xVals);

		      }
		      else theoryCompPYT6[jI]->DrawCopy("HIST C SAME");
		    }

		    pads_p[0]->cd();
		  }
		  
		  for(ULong64_t jI = 0; jI < (ULong64_t)nJtAlgos; ++jI){
		    if(jI == (ULong64_t)r10Pos) continue;
		    //		    if(jI != (ULong64_t)r4Pos && jI != (ULong64_t)r8Pos) continue;
		    
		    createRAA(theoryCompPYT8[jI], theoryCompPYT8[r10Pos]);
		    defaultPlotSetR(theoryCompPYT8[jI], jtAlgosPP[jI]);		
		    theoryCompPYT8[jI]->SetLineStyle(1);
		    theoryCompPYT8[jI]->SetLineWidth(3);
		    theoryCompPYT8[jI]->DrawCopy("HIST C SAME");
		  }

		  if(ratioBool[ratI]){
		    pads_p[1]->cd();

		    for(ULong64_t jI = 0; jI < (ULong64_t)nJtAlgos; ++jI){
		      if(jI == (ULong64_t)r10Pos) continue;
		      if(jI != (ULong64_t)r4Pos && jI != (ULong64_t)r2Pos) continue;

		      for(Int_t bIX = 0; bIX < theoryCompPYT8[jI]->GetNbinsX(); ++bIX){
			Double_t val = theoryCompPYT8[jI]->GetBinContent(bIX+1)/histVectPPClones[jI]->GetBinContent(bIX+1);
			Double_t relErr = val*histVectPPClones[jI]->GetBinError(bIX+1)/histVectPPClones[jI]->GetBinContent(bIX+1);

			theoryCompPYT8[jI]->SetBinContent(bIX+1, val);
			theoryCompPYT8[jI]->SetBinError(bIX+1, relErr);
		      }		      


		      theoryCompPYT8[jI]->GetYaxis()->SetNdivisions(505);
		      theoryCompPYT8[jI]->GetYaxis()->SetTitle("Theory/Data");
		      theoryCompPYT8[jI]->SetMaximum(1.25);
		      theoryCompPYT8[jI]->SetMinimum(0.75);

		      theoryCompPYT8[jI]->DrawCopy("HIST C SAME");
		    }
		    

		    pads_p[0]->cd();
		  }
		  
		  bool isDonePYT6 = false;
		  bool isDonePYT8 = false;
		  for(ULong64_t jI = 0; jI < (ULong64_t)nJtAlgos; ++jI){
		    int orderPos = jtAlgosROrdered[jI];
		    if(orderPos == r10Pos) continue;
		    
		    int rVal = getRFromJetString(jtAlgosPbPb[orderPos]);
		    std::string rStr = getRStrFromInt(rVal);
		    
		    if(!isDonePYT6){
		      leg2_p->AddEntry(theoryCompPYT6[orderPos], "PYTHIA6", "L");
		      isDonePYT6 = true;
		    }
		    else leg2_p->AddEntry(theoryCompPYT6[orderPos], "", "L");
		    
		    if(!isDonePYT8){
		      leg3_p->AddEntry(theoryCompPYT8[orderPos], "PYTHIA8", "L");
		      isDonePYT8 = true;
		    }
		    else leg3_p->AddEntry(theoryCompPYT8[orderPos], "", "L");
		  }	      

		std::cout << "LINE: " << __LINE__ << std::endl;
		  
		  TLine* line_p = new TLine();
		  line_p->SetLineStyle(2);
		  line_p->DrawLine(histVectPPClones[0]->GetXaxis()->GetBinLowEdge(1), 1., histVectPPClones[0]->GetXaxis()->GetBinLowEdge(histVectPPClones[0]->GetNbinsX()+1), 1.);
		  
		  leg_p->Draw("SAME");	      	    	      
		  leg2_p->Draw("SAME");	      	    	      
		  leg3_p->Draw("SAME");	      	    	      
		  
		  //		  canv_p->cd();	      
		  gPad->SetLogx();
		  gPad->RedrawAxis();
		  if(ratioBool[ratI]){
		    pads_p[0]->cd();
		    pads_p[0]->RedrawAxis();

		    drawWhiteBox(160, 199.5, -0.01, 0.05);
		      

		    pads_p[1]->cd();
		    pads_p[1]->RedrawAxis();

		    line_p->DrawLine(histVectPPClones[0]->GetXaxis()->GetBinLowEdge(1), 1., histVectPPClones[0]->GetXaxis()->GetBinLowEdge(histVectPPClones[0]->GetNbinsX()+1), 1.);
		  }
		  delete line_p;

		  drawAllLabels(canv_p, pads_p[0], histVectPPClones[0], label_p, nXVals, xVals, "", jtAbsEtaBinsStr[aI], "", true);
		  
		  checkMakeDir("pdfDir");
		  checkMakeDir("pdfDir/" + dateStr);
		  
		  const std::string saveName = "pdfDir/" + dateStr + "/spectraRat_PP_" + idStr[idI] + "_" + responseModStr[rI] + "_" + jtAbsEtaBinsStr[aI] + "_" + smoothStr[smoothI] + "_" + redStatStr[redI] + "_" + ratioStr[ratI] + "_" + dateStr+ ".pdf";
		  
		  quietSaveAs(canv_p, saveName);
		  std::cout << "LINE: " << __LINE__ << ", " << ratI << std::endl;
		  
		  for(unsigned int tI = 0; tI < theoryCompPYT6.size(); ++tI){
		    delete theoryCompPYT6[tI];
		  }
		  
		  for(unsigned int tI = 0; tI < theoryCompPYT8.size(); ++tI){
		    delete theoryCompPYT8[tI];
		  }

		  delete canv_p;
		  delete leg_p;
		  delete leg2_p;
		  delete leg3_p;
		  delete label_p;	      	
		}
		
		deleteAll(outFile_p, &histVectPPClones);
		histVectPPClones.clear();
		
	      }
	    }
	  }
	}
      }
    }
 
    //RRAA
    for(Int_t redI = 0; redI < nRedStat; ++redI){
      for(Int_t smoothI = 0; smoothI < nSmooth; ++smoothI){
	for(ULong64_t idI = 0; idI < (ULong64_t)nID; ++idI){
	  if(!goodID[idI]) continue;
	  
	  for(ULong64_t rI = 0; rI < (ULong64_t)nResponseMod; ++rI){
	    if(!goodResponseMod[rI]) continue;
	    
	    for(ULong64_t aI = 0; aI < (ULong64_t)nJtAbsEtaBins; ++aI){
	      if(!goodJtAbsEtaBins[aI]) continue;	
	      
	      std::vector<TH1D*> histVectPbPbClones;
	      std::vector<TH1D*> histVectPbPbDenomClones;
	      histVectPbPbClones.reserve(nJtAlgos*nSystTotal*nCentBins);
	      histVectPbPbDenomClones.reserve(nSystTotal*nCentBins);
	      std::vector<TH1D*> histVectPPClones;	 
	      std::vector<TH1D*> histVectPPDenomClones;	 
	      histVectPPClones.reserve(nJtAlgos*nSystTotal*nCentBins);
	      histVectPPDenomClones.reserve(nSystTotal*nCentBins);
	      
	      std::vector<double> globalBins;
	      bool goodBins[nMaxGeneralPtBins][nMaxJtAlgos][nMaxCentBins];
	      
	      for(ULong64_t jI = 0; jI < nJtAlgos; ++jI){
		for(Int_t cI = 0; cI < nCentBins; ++cI){
		  
		  for(Int_t gI = 0; gI < nGenJtPtBinsRRAA[jI][cI]+1; ++gI){
		    bool binEdgeFound = false;
		    
		    for(unsigned int bI = 0; bI < globalBins.size(); ++bI){
		      if(TMath::Abs(globalBins[bI] - genJtPtBinsRRAA[jI][cI][gI]) < 1){
			binEdgeFound = true;
			break;
		      }
		    }
		    
		    if(!binEdgeFound) globalBins.push_back(genJtPtBinsRRAA[jI][cI][gI]);
		  }
		}
	      }
	      
	      std::sort(std::begin(globalBins), std::end(globalBins));
	      
	      for(unsigned int bI = 0; bI < globalBins.size(); ++bI){
		for(ULong64_t jI = 0; jI < nJtAlgos; ++jI){
		  for(Int_t cI = 0; cI < nCentBins; ++cI){
		    goodBins[bI][jI][cI] = true;
		    
		    bool isFound = false;
		    for(Int_t gI = 0; gI < nGenJtPtBinsRRAA[jI][cI]+1; ++gI){
		      if(TMath::Abs(globalBins[bI] - genJtPtBinsRRAA[jI][cI][gI]) < 1.){
			isFound = true;
			break;
		      }
		    }
		    
		    if(!isFound) goodBins[bI][jI][cI] = false;
		  }
		}
	      }
	      
	      std::cout << "GLOBAL BINS: ";
	      for(unsigned int bI = 0; bI < globalBins.size(); ++bI){
		std::cout << bI << ": " << globalBins[bI] << std::endl;
	      }
	      
	      for(ULong64_t sI = 0; sI < (ULong64_t)nSystTotal; ++sI){
		ULong64_t sIPos = sI;
		if(sI >= (ULong64_t)nSystInFile) sIPos = 0;
		
		for(ULong64_t cI = 0; cI < (ULong64_t)nCentBins; ++cI){
		  for(ULong64_t jI = 0; jI < (ULong64_t)nJtAlgos; ++jI){
		    
		    ULong64_t idKey = getKey(0, jI, cI, idI, rI, aI, sIPos);
		    ULong64_t idKeyBBPP = histKeyToBestBayesKeyPP[idKey];
		    ULong64_t idKeyBBPbPb = histKeyToBestBayesKeyPbPb[idKey];
		    
		    if(systStrTotal[sI].find("PriorAdj") != std::string::npos){
		      getNewBBPos(systStrTotal[sI], idKey, idKeyBBPP, idKeyBBPbPb);
		    }
		    else if(systStrTotal[sI].find("SVD") != std::string::npos){
		      idKeyBBPP = histKeyToBestSvdKeyPP[idKey];
		      idKeyBBPbPb = histKeyToBestSvdKeyPbPb[idKey];		  
		    }
		    
		    
		    ULong64_t vectPos = keyToVectPos[idKeyBBPP];
		    std::string tempName = "";
		    
		    if(systStrTotal[sI].find("SVD") != std::string::npos){
		      vectPos = keyToVectPosSvd[idKeyBBPP];		    
		      tempName = "raa_" + std::string(histVectPPSvd[vectPos]->GetName()) + "_Clone_" + systStrTotal[sI];
		    }
		    else{
		      tempName = "raa_" + std::string(histVectPP[vectPos]->GetName()) + "_Clone_" + systStrTotal[sI];
		    }
		    
		    histVectPPClones.push_back(new TH1D(tempName.c_str(), ";Jet p_{T} [GeV];R_{AA}", nGenJtPtBinsRRAA[jI][cI], genJtPtBinsRRAA[jI][cI]));
		    
		    if(jI == (ULong64_t)r2Pos) histVectPPDenomClones.push_back(new TH1D((tempName + "_DENOM").c_str(), ";Jet p_{T} [GeV];R_{AA}", nGenJtPtBinsRRAA[jI][cI], genJtPtBinsRRAA[jI][cI]));
		    
		    if(systStrTotal[sI].find("SVD") != std::string::npos){
		      if(!macroHistToSubsetHist(histVectPPSvd[vectPos], histVectPPClones[histVectPPClones.size()-1])) return 1;	     
		      
		      if(jI == (ULong64_t)r2Pos){
			if(!macroHistToSubsetHist(histVectPPSvd[vectPos], histVectPPDenomClones[histVectPPDenomClones.size()-1])) return 1;
		      }	     
		    }
		    else{
		      if(!macroHistToSubsetHist(histVectPP[vectPos], histVectPPClones[histVectPPClones.size()-1])) return 1;	     
		      
		      if(jI == (ULong64_t)r2Pos){
			if(!macroHistToSubsetHist(histVectPP[vectPos], histVectPPDenomClones[histVectPPDenomClones.size()-1])) return 1;	     
		      }
		    }
		    
		    vectPos = keyToVectPos[idKeyBBPbPb];
		    if(systStrTotal[sI].find("SVD") != std::string::npos){
		      vectPos = keyToVectPosSvd[idKeyBBPbPb];
		      tempName = "raa_" + std::string(histVectPbPbSvd[vectPos]->GetName()) + "_Clone_" + systStrTotal[sI];
		    }
		    else{
		      tempName = "raa_" + std::string(histVectPbPb[vectPos]->GetName()) + "_Clone_" + systStrTotal[sI];
		    }
		    
		    histVectPbPbClones.push_back(new TH1D(tempName.c_str(), ";Jet p_{T} [GeV];R_{AA}", nGenJtPtBinsRRAA[jI][cI], genJtPtBinsRRAA[jI][cI]));
		    
		    if(jI == (ULong64_t)r2Pos) histVectPbPbDenomClones.push_back(new TH1D((tempName + "_DENOM").c_str(), ";Jet p_{T} [GeV];R_{AA}", nGenJtPtBinsRRAA[jI][cI], genJtPtBinsRRAA[jI][cI]));
		    
		    
		    if(systStrTotal[sI].find("SVD") != std::string::npos){
		      if(!macroHistToSubsetHist(histVectPbPbSvd[vectPos], histVectPbPbClones[histVectPbPbClones.size()-1])) return 1;
		      
		      if(jI == (ULong64_t)r2Pos){
			if(!macroHistToSubsetHist(histVectPbPbSvd[vectPos], histVectPbPbDenomClones[histVectPbPbDenomClones.size()-1])) return 1;
		      }
		    }
		    else{
		      if(!macroHistToSubsetHist(histVectPbPb[vectPos], histVectPbPbClones[histVectPbPbClones.size()-1])) return 1;
		      
		      if(jI == (ULong64_t)r2Pos){
			if(!macroHistToSubsetHist(histVectPbPb[vectPos], histVectPbPbDenomClones[histVectPbPbDenomClones.size()-1])) return 1;
		      }
		    }
		  }
		}
	      }
	      
	      //Div by bin widths;
	      divHistByWidth(histVectPPClones);
	      divHistByWidth(histVectPbPbClones);
	      divHistByWidth(histVectPPDenomClones);
	      divHistByWidth(histVectPbPbDenomClones);
	    
	      for(ULong64_t cI = 0; cI < (ULong64_t)nCentBins; ++cI){
		for(ULong64_t sI = 0; sI < (ULong64_t)nSystTotal; ++sI){
		  Double_t tempTAAFactor = getTAAScaleFactorNB(centBinsStr[cI]);	     
		  if(isStrSame(systStrTotal[sI], "TAAUp")) tempTAAFactor += tempTAAFactor*getTAAScaleFactorUp(centBinsStr[cI]);
		  else if(isStrSame(systStrTotal[sI], "TAADown")) tempTAAFactor -= tempTAAFactor*getTAAScaleFactorUp(centBinsStr[cI]);		
		  
		  Double_t totalFactor = tempTAAFactor*nMBEvents*2.*jtAbsEtaBinsWidth[aI]*centBinsWidth[cI]/100.;
		  scaleHist(histVectPbPbDenomClones[cI + sI*nCentBins], 1./totalFactor);
		  
		  
		  Double_t lumiFactorToApply = lumiFactor;		
		  if(isStrSame(systStrTotal[sI], "LumiUp")) lumiFactorToApply += getLumiAbsError();
		  else if(isStrSame(systStrTotal[sI], "LumiDown")) lumiFactorToApply -= getLumiAbsError();		
		  
		  scaleHist(histVectPPDenomClones[cI + sI*nCentBins], 1./(2.*jtAbsEtaBinsWidth[aI]*lumiFactorToApply));
		}
	      }
	      
	      
	      for(ULong64_t jI = 0; jI < (ULong64_t)nJtAlgos; ++jI){
		for(ULong64_t cI = 0; cI < (ULong64_t)nCentBins; ++cI){
		  for(ULong64_t sI = 0; sI < (ULong64_t)nSystTotal; ++sI){
		    Double_t tempTAAFactor = getTAAScaleFactorNB(centBinsStr[cI]);	     
		    if(isStrSame(systStrTotal[sI], "TAAUp")) tempTAAFactor += tempTAAFactor*getTAAScaleFactorUp(centBinsStr[cI]);
		    else if(isStrSame(systStrTotal[sI], "TAADown")) tempTAAFactor -= tempTAAFactor*getTAAScaleFactorUp(centBinsStr[cI]);		
		    
		    Double_t totalFactor = tempTAAFactor*nMBEvents*2.*jtAbsEtaBinsWidth[aI]*centBinsWidth[cI]/100.;
		    scaleHist(histVectPbPbClones[jI + cI*nJtAlgos + sI*nCentBins*nJtAlgos], 1./totalFactor);
		    
		    
		    Double_t lumiFactorToApply = lumiFactor;		
		    if(isStrSame(systStrTotal[sI], "LumiUp")) lumiFactorToApply += getLumiAbsError();
		    else if(isStrSame(systStrTotal[sI], "LumiDown")) lumiFactorToApply -= getLumiAbsError();		
		    
		    scaleHist(histVectPPClones[jI + cI*nJtAlgos + sI*nCentBins*nJtAlgos], 1./(2.*jtAbsEtaBinsWidth[aI]*lumiFactorToApply));
		    
		    bool isSystDoubleCorr = false;
		    for(ULong64_t sI2 = 0; sI2 < (ULong64_t)doubleCorrSystStr.size(); ++sI2){
		      if(isStrSame(doubleCorrSystStr[sI2], systStrTotal[sI])){
			isSystDoubleCorr = true;
			break;
		      }
		    }
		    
		    if(isSystDoubleCorr){
		      createRAA(histVectPbPbClones[jI + cI*nJtAlgos + sI*nCentBins*nJtAlgos], histVectPbPbDenomClones[cI + sI*nCentBins]);
		      createRAA(histVectPPClones[jI + cI*nJtAlgos + sI*nCentBins*nJtAlgos], histVectPPDenomClones[cI + sI*nCentBins]);
		    }
		    else{
		      createRAA(histVectPbPbClones[jI + cI*nJtAlgos + sI*nCentBins*nJtAlgos], histVectPbPbDenomClones[cI]);		  
		      createRAA(histVectPPClones[jI + cI*nJtAlgos + sI*nCentBins*nJtAlgos], histVectPPDenomClones[cI]);		  
		    }


		    if(isStrSame(redStatStr[redI], "RedStat")){
		      for(Int_t bIX = 0; bIX < histVectPbPbClones[jI + cI*nJtAlgos + sI*nCentBins*nJtAlgos]->GetXaxis()->GetNbins(); ++bIX){
			Double_t tempErr = histVectPbPbClones[jI + cI*nJtAlgos + sI*nCentBins*nJtAlgos]->GetBinError(bIX+1);
			
			if(rValI[jI] == 3) tempErr *= 0.38;
			else if(rValI[jI] == 4) tempErr *= 0.472;
			else if(rValI[jI] == 6) tempErr *= 0.5917;
			else if(rValI[jI] == 8) tempErr *= 0.648;
			else if(rValI[jI] == 10) tempErr *= 0.701;
			histVectPbPbClones[jI + cI*nJtAlgos + sI*nCentBins*nJtAlgos]->SetBinError(bIX+1, tempErr);		      
		      }
		      
		      for(Int_t bIX = 0; bIX < histVectPPClones[jI + cI*nJtAlgos + sI*nCentBins*nJtAlgos]->GetXaxis()->GetNbins(); ++bIX){
			Double_t tempErr = histVectPPClones[jI + cI*nJtAlgos + sI*nCentBins*nJtAlgos]->GetBinError(bIX+1);
			
			if(rValI[jI] == 3) tempErr *= 0.38;
			else if(rValI[jI] == 4) tempErr *= 0.472;
			else if(rValI[jI] == 6) tempErr *= 0.5917;
			else if(rValI[jI] == 8) tempErr *= 0.648;
			else if(rValI[jI] == 10) tempErr *= 0.701;
			histVectPPClones[jI + cI*nJtAlgos + sI*nCentBins*nJtAlgos]->SetBinError(bIX+1, tempErr);		      
		      }
		    }
		  }
		}
	      }
	      
	      for(ULong64_t jI = 0; jI < (ULong64_t)nJtAlgos; ++jI){
		if(jI == (ULong64_t)r2Pos) continue;
		
		for(ULong64_t cI = 0; cI < (ULong64_t)nCentBins; ++cI){
		  for(ULong64_t sI = 0; sI < (ULong64_t)nSystTotal; ++sI){
		    //		    if(systStrTotal[sI].size() == 0) continue;
		    
		    bool isSystCorr = false;
		    for(ULong64_t sI2 = 0; sI2 < (ULong64_t)doubleCorrSystStr.size(); ++sI2){
		      if(isStrSame(doubleCorrSystStr[sI2], systStrTotal[sI])){
			isSystCorr = true;
			break;
		      }
		    }
		    
		    if(isSystCorr) createRAA(histVectPbPbClones[jI + cI*nJtAlgos + sI*nCentBins*nJtAlgos], histVectPPClones[jI + cI*nJtAlgos + sI*nCentBins*nJtAlgos]);
		    else createRAA(histVectPbPbClones[jI + cI*nJtAlgos + sI*nCentBins*nJtAlgos], histVectPPClones[jI + cI*nJtAlgos]);		  
		  }
		}
	      }	   
	      
	      for(ULong64_t cI = 0; cI < (ULong64_t)nCentBins; ++cI){
		TH1D* temp_p = new TH1D("temp_h", "", nGenJtPtBinsRRAA[r2Pos][cI], genJtPtBinsRRAA[r2Pos][cI]);
		macroHistToSubsetHist(histVectPbPbClones[r2Pos + cI*nJtAlgos], temp_p, true);
		
		for(ULong64_t sI = 0; sI < (ULong64_t)nSystTotal; ++sI){
		  //		  if(systStrTotal[sI].size() == 0) continue;
		  
		  bool isSystCorr = false;
		  for(ULong64_t sI2 = 0; sI2 < (ULong64_t)doubleCorrSystStr.size(); ++sI2){
		    if(isStrSame(doubleCorrSystStr[sI2], systStrTotal[sI])){
		      isSystCorr = true;
		      break;
		    }
		  }
		  
		  if(isSystCorr) createRAA(histVectPbPbClones[r2Pos + cI*nJtAlgos + sI*nCentBins*nJtAlgos], histVectPbPbClones[r2Pos + cI*nJtAlgos + sI*nCentBins*nJtAlgos]);
		  else createRAA(histVectPbPbClones[r2Pos + cI*nJtAlgos + sI*nCentBins*nJtAlgos], temp_p);
		}
		delete temp_p;
	      }
	      
	      
	      std::vector<std::vector<std::vector<std::vector<Double_t> > > > reducedSystPbPb;
	      std::vector<std::vector<std::vector<Double_t> > > reducedSystSumPbPb;
	      
	      std::vector<std::vector<Int_t> > nGenJtPtBinsVect;
	      for(ULong64_t jI = 0; jI < (ULong64_t)nJtAlgos; ++jI){
		nGenJtPtBinsVect.push_back({});	    
		for(Int_t cI = 0; cI < nCentBins; ++cI){
		  nGenJtPtBinsVect[jI].push_back(nGenJtPtBinsRRAA[jI][cI]);
		}
	      }
	      
	      setupSyst(nReducedSyst, nGenJtPtBinsVect, nGenJtPtBinsVect, nJtAlgos, nCentBins, &reducedSystPbPb, &reducedSystSumPbPb);
	      
	      if(smoothBool[smoothI]){
		for(ULong64_t sI = 0; sI < (ULong64_t)nSystTotal; ++sI){
		  if(systStrTotal[sI].size() == 0) continue;
		  
		  for(ULong64_t cI = 0; cI < (ULong64_t)nCentBins; ++cI){
		    for(ULong64_t jI = 0; jI < (ULong64_t)nJtAlgos; ++jI){
		      
		      if(histVectPbPbClones[jI + cI*nJtAlgos]->GetNbinsX() < 3) continue;
		      smoothErrors(histVectPbPbClones[jI + cI*nJtAlgos + sI*nCentBins*nJtAlgos], histVectPbPbClones[jI + cI*nJtAlgos]);
		    }
		  }
		}
	      }
	      
	      //Extract systematics
	      for(ULong64_t sI = 0; sI < (ULong64_t)nSystTotal; ++sI){
		
		ULong64_t pos = 0;
		if(systToCombo.count(systStrTotal[sI]) > 0) pos = posInStrVectExact(systToCombo[systStrTotal[sI]], reducedSystStr);
		else pos = posInStrVectExact(systStrTotal[sI], reducedSystStr);
	      
		for(ULong64_t jI = 0; jI < (ULong64_t)nJtAlgos; ++jI){
		  for(ULong64_t cI = 0; cI < (ULong64_t)nCentBins; ++cI){
		    
		    
		    for(Int_t bI = 0; bI < nGenJtPtBinsRRAA[jI][cI]; ++bI){
		      Double_t delta = TMath::Abs(histVectPbPbClones[jI + cI*nJtAlgos + sI*nCentBins*nJtAlgos]->GetBinContent(bI+1) - histVectPbPbClones[jI + cI*nJtAlgos]->GetBinContent(bI+1));
		      
		      (reducedSystPbPb[jI][cI][pos])[bI] = TMath::Min((reducedSystPbPb[jI][cI][pos])[bI], delta);
		    }		  
		  }
		}
	      }
	      
	      for(ULong64_t jI = 0; jI < (ULong64_t)nJtAlgos; ++jI){
		for(ULong64_t cI = 0; cI < (ULong64_t)nCentBins; ++cI){
		  for(Int_t bI = 0; bI < nGenJtPtBinsRRAA[jI][cI]; ++bI){
		    for(ULong64_t sI = 0; sI < (ULong64_t)reducedSystStr.size(); ++sI){
		      (reducedSystSumPbPb[jI][cI])[bI] = TMath::Sqrt((reducedSystSumPbPb[jI][cI])[bI]*(reducedSystSumPbPb[jI][cI])[bI] + (reducedSystPbPb[jI][cI][sI])[bI]*(reducedSystPbPb[jI][cI][sI])[bI]);
		    }
		  }		
		}
	      }
	      
	      for(unsigned int bI = 0; bI < globalBins.size()-1; ++bI){
		if(!goodBins[bI]) continue;
		if(!goodBins[bI+1]) continue;
		
		Double_t globalBinVal = (globalBins[bI] + globalBins[bI+1])/2.;
		std::string globalBinsStr = prettyString(globalBins[bI], 1, true) + "to" + prettyString(globalBins[bI+1], 1, true);
		std::string globalBinsStr2 = std::to_string((int)globalBins[bI]) + " < p_{T} < " + std::to_string((int)globalBins[bI+1]);
		
		TLatex* label_p = new TLatex();
		labelStandard(label_p);
		
		TLegend* leg_p = new TLegend(0.25, 0.12, 0.4, 0.35);
		legStandard(leg_p);
		
		TCanvas* canv_p = new TCanvas("canv_p", "", 450, 450);	 
		defineCanv(canv_p);		
		
		std::vector<TH1D*> rraa_p;
		rraa_p.reserve(nCentBins*nSystTotal);
		std::vector<TGraphAsymmErrors*> rraaTGA_p;
		rraaTGA_p.reserve(nCentBins*nSystTotal);
		
		bool isDrawn = false;
		for(ULong64_t cI = 0; cI < (ULong64_t)nCentBins; ++cI){
		  
		  for(ULong64_t sI = 0; sI < (ULong64_t)nSystTotal; ++sI){
		    rraa_p.push_back(new TH1D(("rraa_" + centBinsStr[cI] + "_" + globalBinsStr + "_" + systStrTotal[sI]).c_str(), (";Jet R;R_{AA}^{R}/R_{AA}^{R=0." + r2Str + "}").c_str(), nRBins, rBins));
		    rraaTGA_p.push_back(NULL);
		    centerTitles(rraa_p[rraa_p.size()-1]);
		    
		    for(ULong64_t jI = 0; jI < (ULong64_t)nJtAlgos; ++jI){
		      Int_t rBinPos = -1;
		      Int_t globalBinPos = -1;
		      
		      if(!goodBins[bI][jI][cI]) continue;
		      
		      for(Int_t bIX = 0; bIX < nRBins; ++bIX){
			if(rValD[jI] >= rBins[bIX] && rValD[jI] < rBins[bIX+1]){
			  rBinPos = bIX;
			  break;
			}
		      }
		      
		      for(Int_t bIX = 0; bIX < nGenJtPtBinsRRAA[jI][cI]; ++bIX){
			if(genJtPtBinsRRAA[jI][cI][bIX] <= globalBinVal && globalBinVal < genJtPtBinsRRAA[jI][cI][bIX+1]){
			  globalBinPos = bIX;
			  break;
			}		      
		      }
		      
		      if(jI != (ULong64_t)r2Pos){
			std::cout << "RBINPOS: " << rBinPos << ", " << r2Pos << ", " << rValD[jI] << ", " << dirListPbPb.size() << std::endl;
			std::cout << " " << dirListPbPb[jI] << std::endl;
			rraa_p[sI + cI*nSystTotal]->SetBinContent(rBinPos+1, histVectPbPbClones[jI + cI*nJtAlgos + sI*nJtAlgos*nCentBins]->GetBinContent(globalBinPos+1));
			rraa_p[sI + cI*nSystTotal]->SetBinError(rBinPos+1, histVectPbPbClones[jI + cI*nJtAlgos + sI*nJtAlgos*nCentBins]->GetBinError(globalBinPos+1));
		      }
		    }
		    
		    canv_p->cd();
		    
		    defaultPlotSet(rraa_p[sI + cI*nSystTotal], centBinsStr[cI]);
		    rraaTGA_p[sI + cI*nSystTotal] = new TGraphAsymmErrors(rraa_p[sI + cI*nSystTotal]);
		    defaultPlotSet(rraaTGA_p[sI + cI*nSystTotal], centBinsStr[cI]);
		    
		    for(Int_t bIX = 0; bIX < rraa_p[sI + cI*nSystTotal]->GetXaxis()->GetNbins(); ++bIX){
		      rraaTGA_p[sI + cI*nSystTotal]->SetPointEXlow(bIX, 0);
		      rraaTGA_p[sI + cI*nSystTotal]->SetPointEXhigh(bIX, 0);
		    }
		    
		    if(sI == 0){
		      if(!isDrawn){
			rraa_p[sI + cI*nSystTotal]->SetMinimum(0.45);
			rraa_p[sI + cI*nSystTotal]->SetMaximum(1.55);
			//		      rraa_p[sI + cI*nSystTotal]->DrawCopy("HIST E1 P");
			TH1D* tempHist_p = (TH1D*)rraa_p[sI + cI*nSystTotal]->Clone("tempHistClone");
			tempHist_p->Scale(1./1000000000.);
			
			tempHist_p->SetMaximum(1.55);
			tempHist_p->SetMinimum(0.45);
			tempHist_p->DrawCopy();
			delete tempHist_p;
			
			
			rraaTGA_p[sI + cI*nSystTotal]->Draw("P");
			isDrawn = true;
		      }
		      else{
			//		      rraa_p[sI + cI*nSystTotal]->DrawCopy("HIST E1 P SAME");
			
			rraaTGA_p[sI + cI*nSystTotal]->Draw("P");
		      }
		      
		      std::string centStr = std::to_string(centBinsLow[cI]) + "-" + std::to_string(centBinsHi[cI]) + "%";
		      if(rraa_p[cI*nSystTotal]->Integral() > TMath::Power(10,-10)){
			leg_p->AddEntry(rraa_p[cI*nSystTotal], centStr.c_str(), "P L");
		      }
		      
		    }
		  }
		}
		
		TLine* line_p = new TLine();
		line_p->SetLineStyle(2);
		line_p->DrawLine(rraa_p[0]->GetXaxis()->GetBinLowEdge(1), 1., rraa_p[0]->GetXaxis()->GetBinLowEdge(rraa_p[0]->GetNbinsX()+1), 1.);
		delete line_p;
		
		for(Int_t cI = 0; cI < nCentBins; ++cI){
		  bool isGood = false;
		  for(Int_t jI = 0; jI < (Int_t)nJtAlgos; ++jI){
		    if(!goodBins[bI][jI][cI]) continue;
		    
		    isGood = true;
		  }
		  if(!isGood) continue;
		  
	     
		  std::string centStr = std::to_string(centBinsLow[cI]) + "-" + std::to_string(centBinsHi[cI]) + "%";
		  //	leg_p->AddEntry(rraa_p[cI*nSystTotal], centStr.c_str(), "P L");
		}
		
		leg_p->Draw("SAME");
	      
		std::vector<std::vector< std::vector<Double_t> > > reducedSystPbPb2;
		for(ULong64_t cI = 0; cI < (ULong64_t)nCentBins; ++cI){
		  reducedSystPbPb2.push_back({});
		  for(ULong64_t sI = 0; sI < (ULong64_t)reducedSystStr.size(); ++sI){
		    reducedSystPbPb2[cI].push_back({});
		  }
		}
		
		for(ULong64_t cI = 0; cI < (ULong64_t)nCentBins; ++cI){
		  std::vector<Double_t> tempSyst;
		  
		  for(Int_t rI = 0; rI < nRBins; ++rI){
		    Int_t jetBinPos = -1;
		    Int_t globalBinPos = -1;
		    
		    for(Int_t jI = 0; jI < (Int_t)nJtAlgos; ++jI){
		      if(rBins[rI] <= rValD[jI] && rValD[jI] < rBins[rI+1]){
			jetBinPos = jI;
			break;
		      }
		    }
		    
		    for(Int_t bIX = 0; bIX < nGenJtPtBinsRRAA[jetBinPos][cI]; ++bIX){
		      if(genJtPtBinsRRAA[jetBinPos][cI][bIX] <= globalBinVal && globalBinVal < genJtPtBinsRRAA[jetBinPos][cI][bIX+1]){
			globalBinPos = bIX;
			break;
		      }		      
		    }
		    
		    if(jetBinPos >= 0){
		      tempSyst.push_back(reducedSystSumPbPb[jetBinPos][cI][globalBinPos]);
		      
		      for(unsigned int sI = 0; sI < reducedSystStr.size(); ++sI){
			reducedSystPbPb2[cI][sI].push_back(reducedSystPbPb[jetBinPos][cI][sI][globalBinPos]);
		      }		      
		    }
		    else{
		      tempSyst.push_back(0.0);		  
		      for(unsigned int sI = 0; sI < reducedSystStr.size(); ++sI){
			reducedSystPbPb2[cI][sI].push_back(0.0);
		      }
		    }
		  }
		  
		  drawSyst(canv_p, NULL, rraa_p[cI*nSystTotal], &tempSyst, rraa_p[cI*nSystTotal]->GetBinLowEdge(1), rraa_p[cI*nSystTotal]->GetBinLowEdge(rraa_p[cI]->GetNbinsX()+1), true);	      
		}
		
		for(ULong64_t cI = 0; cI < (ULong64_t)nCentBins; ++cI){
		  rraaTGA_p[cI*nSystTotal]->Draw("P");
		}
		
		for(ULong64_t cI = 0; cI < (ULong64_t)nCentBins; ++cI){
		  std::string saveStr = "rraa_" + globalBinsStr + "_PbPb_" + centBinsStr[cI] + "_" + idStr[idI] + "_" + responseModStr[rI] + "_" + jtAbsEtaBinsStr[aI] + "_" + smoothStr[smoothI];
		  plotSyst(rraa_p[cI*nSystTotal], reducedSystPbPb2, reducedSystStr, cI, saveStr);
		}
	      	
		canv_p->cd();	      
		gPad->RedrawAxis();
		drawAllLabels(canv_p, NULL, rraa_p[0], label_p, nRValLabels, rValLabels, "", jtAbsEtaBinsStr[aI], globalBinsStr2);
		
		checkMakeDir("pdfDir");
		checkMakeDir("pdfDir/" + dateStr);
		
		const std::string saveName = "pdfDir/" + dateStr + "/rraa_" + globalBinsStr + "_" + idStr[idI] + "_" + responseModStr[rI] + "_" + jtAbsEtaBinsStr[aI] + "_" + smoothStr[smoothI] + "_" + redStatStr[redI] + "_" + dateStr+ ".pdf";
		
		quietSaveAs(canv_p, saveName);
		
		for(ULong64_t cI = 0; cI < (ULong64_t)nCentBins; ++cI){
		  for(ULong64_t sI = 0; sI < (ULong64_t)nSystTotal; ++sI){
		    delete rraa_p[sI + cI*nSystTotal];
		    delete rraaTGA_p[sI + cI*nSystTotal];
		  }
		}
		
		
		delete canv_p;
		delete leg_p;
		delete label_p;
	      }	   	    
	      
	      deleteAll(outFile_p, &histVectPPClones);
	      deleteAll(outFile_p, &histVectPbPbClones);
	      deleteAll(outFile_p, &histVectPPDenomClones);
	      deleteAll(outFile_p, &histVectPbPbDenomClones);
	      histVectPPClones.clear();
	      histVectPbPbClones.clear();	   	    
	      histVectPPDenomClones.clear();
	      histVectPbPbDenomClones.clear();	   	    
	    }
	  }
	}
      }
    }
  }
 
  for(Int_t fI = 0; fI < nFiles; ++fI){
    inFile_p[fI]->Close();
    delete inFile_p[fI];
  }
  
  std::cout << "Infile closed Files" << std::endl;  

  for(unsigned int i = 0; i < cutProps.size(); ++i){
    delete (cutProps[i]);
  }

  std::cout << "Props cut" << std::endl;  

  outFile_p->Close();
  delete outFile_p;
 
  std::cout << "Outfile closed Files" << std::endl;  

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc < 3 || argc > 6){
    std::cout << "Usage: ./bin/plotUnfoldedAll.exe <inFileNamePP> <inFileNamePbPb> <inATLASFileName-opt> <tagStr-opt> <fileBinOverride>. return 1" << std::endl;
    return 1;
  }

  int retVal = 0;
  if(argc == 3) retVal += plotUnfoldedAll(argv[1], argv[2]);
  else if(argc == 4) retVal += plotUnfoldedAll(argv[1], argv[2], argv[3]);
  else if(argc == 5) retVal += plotUnfoldedAll(argv[1], argv[2], argv[3], argv[4]);
  else if(argc == 6) retVal += plotUnfoldedAll(argv[1], argv[2], argv[3], argv[4], argv[5]);
  return retVal;
}
