#ifndef SYSTFUNCTIONS_H
#define SYSTFUNCTIONS_H

//cpp dependencies
#include <string>
#include <vector>

//ROOT dependencies
#include "TBox.h"
#include "TCanvas.h"
#include "TDatime.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TPad.h"
#include "TStyle.h"

//Non-local (Utility, etc.) dependencies
#include "Utility/include/checkMakeDir.h"
#include "Utility/include/getLinBins.h"
#include "Utility/include/getLogBins.h"
#include "Utility/include/histDefUtility.h"
#include "Utility/include/kirchnerPalette.h"
#include "Utility/include/lumiAndTAAUtil.h"
#include "Utility/include/plotUtilities.h"
#include "Utility/include/vanGoghPalette.h"

std::vector<double> getSyst(const std::string dirStr, TH1D* nominal_p, std::vector<TH1D*> syst_p, std::vector<std::string> systStr, Double_t minXVal, Double_t maxXVal, std::vector<std::string>* plotNames, bool doSmoothing, std::vector<std::vector<std::string > > systToCombo)
{
  std::cout << "Starting getSyst: " << std::endl;
  std::cout << "SystHist input: " << std::endl;
  for(unsigned int i = 0; i < syst_p.size(); ++i){
    std::cout << syst_p.at(i)->GetName() << std::endl;
  }
  for(unsigned int sI = 0; sI < systStr.size(); ++sI){
    std::cout << systStr.at(sI) << std::endl;
  }

  //Setting some style and color for TH1

  vanGoghPalette vg;
  const Int_t nColor = 5;
  const Int_t colors[nColor] = {vg.getColor(0), vg.getColor(1), vg.getColor(2), vg.getColor(3), vg.getColor(4)};

  const Int_t nStyle = 4;
  Int_t styles[nStyle] = {20, 21, 47, 34};

  TDatime* date = new TDatime();
  const std::string dateStr = std::to_string(date->GetDate());
  delete date;

  double smallDelta = 0.000001;
  std::vector<double> systVect;
  std::vector<double> nomValVect;

  //Building parallel systematic histograms

  Int_t tempNBins = 0;
  for(Int_t bIX = 0; bIX < nominal_p->GetNbinsX(); ++bIX){
    Double_t binCenter = nominal_p->GetBinCenter(bIX+1);
    if(binCenter > maxXVal || binCenter < minXVal) continue;

    tempNBins++;
  }

  const Int_t nBins = tempNBins;
  Double_t bins[nBins+1];

  Int_t binFillPos = 0;
  for(Int_t bIX = 0; bIX < nominal_p->GetNbinsX(); ++bIX){
    Double_t binCenter = nominal_p->GetBinCenter(bIX+1);
    if(binCenter < minXVal) continue;
    if(binCenter > maxXVal){
      bins[binFillPos] = nominal_p->GetBinLowEdge(bIX+1);
      break;
    }

    bins[binFillPos] = nominal_p->GetBinLowEdge(bIX+1);
    binFillPos++;
  }

  std::cout << "CHECK BINS:" << std::endl;
  for(Int_t bIX = 0; bIX < nBins+1; ++bIX){
    std::cout << bins[bIX] << ", ";
  }
  std::cout << std::endl;

  const Int_t nSystHist = syst_p.size();
  bool systToBeCombined[nSystHist];
  std::vector<std::vector<int> > postToCombine;
  std::string systLegStr[nSystHist];
  bool systToBeSkipped[nSystHist];

  for(Int_t sI = 0; sI < nSystHist; ++sI){
    systToBeCombined[sI] = false;
    systLegStr[sI] = systStr.at(sI);
    systToBeSkipped[sI] = false;
    postToCombine.push_back({});

    for(unsigned int cI = 0; cI < systToCombo.size(); ++cI){
      if(isStrSame(systStr.at(sI), systToCombo.at(cI).at(0))){
        systToBeCombined[sI] = true;
        systLegStr[sI] = systToCombo.at(cI).at(systToCombo.at(cI).size()-1);

	for(Int_t sI2 = 0; sI2 < nSystHist; ++sI2){
          bool doCombine = false;
          for(unsigned int skI = 1; skI < systToCombo.at(cI).size()-1; ++skI){
            if(isStrSame(systStr.at(sI2), systToCombo.at(cI).at(skI))){
              doCombine = true;
              break;
            }
          }

          if(doCombine) postToCombine.at(postToCombine.size()-1).push_back(sI2);
	}

        break;
      }
      else{
        bool doBreak = false;
	for(unsigned int skI = 1; skI < systToCombo.at(cI).size()-1; ++skI){
          if(isStrSame(systStr.at(sI), systToCombo.at(cI).at(skI))){
            systToBeSkipped[sI] = true;
            doBreak = true;
            break;
          }
        }
        if(doBreak) break;
      }
    }
  }

  TH1D* systHist_p[nSystHist];

  for(Int_t sI = 0; sI < nSystHist; ++sI){
    systHist_p[sI] = new TH1D(("systHist_" + std::to_string(sI)).c_str(), "", nBins, bins);

    Int_t binFillPos = 0;
    for(Int_t bIX = 0; bIX < syst_p.at(sI)->GetNbinsX(); ++bIX){
      Double_t binCenter = nominal_p->GetBinCenter(bIX+1);
      if(binCenter > maxXVal || binCenter < minXVal) continue;

      Double_t binVal = TMath::Abs(syst_p.at(sI)->GetBinContent(bIX+1) - nominal_p->GetBinContent(bIX+1))/(nominal_p->GetBinContent(bIX+1));
      Double_t binErr = nominal_p->GetBinContent(bIX+1)*nominal_p->GetBinContent(bIX+1)/(nominal_p->GetBinError(bIX+1)*nominal_p->GetBinError(bIX+1));
      binErr += syst_p.at(sI)->GetBinContent(bIX+1)*syst_p.at(sI)->GetBinContent(bIX+1)/(syst_p.at(sI)->GetBinError(bIX+1)*syst_p.at(sI)->GetBinError(bIX+1));
      binErr = binVal*TMath::Sqrt(binErr);

      systHist_p[sI]->SetBinContent(binFillPos+1, binVal);
      systHist_p[sI]->SetBinError(binFillPos+1, binErr);
      binFillPos++;
    }

    if(doSmoothing) systHist_p[sI]->Smooth();

    binFillPos = 0;

    for(Int_t bIX = 0; bIX < syst_p.at(sI)->GetNbinsX(); ++bIX){
      Double_t binCenter = nominal_p->GetBinCenter(bIX+1);
      if(binCenter > maxXVal || binCenter < minXVal) continue;

      Double_t binVal = systHist_p[sI]->GetBinContent(binFillPos+1)*nominal_p->GetBinContent(bIX+1);
      Double_t binErr = systHist_p[sI]->GetBinError(binFillPos+1)*nominal_p->GetBinContent(bIX+1);

      syst_p.at(sI)->SetBinContent(bIX+1, binVal);
      syst_p.at(sI)->SetBinError(bIX+1, binErr);
      binFillPos++;
    }

    delete systHist_p[sI];
  }


  for(unsigned int sI = 0; sI < syst_p.size(); ++sI){
    if(!systToBeCombined[sI]) continue;

    for(unsigned int cI = 0; cI < postToCombine.at(sI).size(); ++cI){
      std::cout << "Combining " << systStr.at(sI) << " and " << systStr.at(postToCombine.at(sI).at(cI)) << std::endl;

      for(Int_t bIX = 0; bIX < syst_p.at(sI)->GetNbinsX(); ++bIX){
	Double_t binCenter = syst_p.at(sI)->GetBinCenter(bIX+1);
        if(binCenter > maxXVal || binCenter < minXVal) continue;

	Double_t binC1 = syst_p.at(sI)->GetBinContent(bIX+1);
        Double_t binC2 = syst_p.at(postToCombine.at(sI).at(cI))->GetBinContent(bIX+1);

	if(binC2 > binC1){
          syst_p.at(sI)->SetBinContent(bIX+1, binC2);
          syst_p.at(sI)->SetBinError(bIX+1, syst_p.at(sI)->GetBinError(bIX+1)*binC2/binC1);
        }
      }
    }
  }

  for(Int_t bIX = 0; bIX < nominal_p->GetNbinsX(); ++bIX){
    Double_t binCenter = nominal_p->GetBinCenter(bIX+1);
    if(binCenter > maxXVal || binCenter < minXVal){
      for(unsigned int sI = 0; sI < syst_p.size(); ++sI){
	syst_p.at(sI)->SetBinContent(bIX+1, 0.0);
        syst_p.at(sI)->SetBinError(bIX+1, 0.0);
      }
      continue;
    }

    Double_t tempVal = 0.0;
    for(unsigned int sI = 0; sI < syst_p.size(); ++sI){
      Double_t binVal = syst_p.at(sI)->GetBinContent(bIX+1)/nominal_p->GetBinContent(bIX+1);
      Double_t tempRelErr = syst_p.at(sI)->GetBinError(bIX+1)/nominal_p->GetBinContent(bIX+1);

      if(!systToBeSkipped[sI]) tempVal = TMath::Sqrt(tempVal*tempVal + syst_p.at(sI)->GetBinContent(bIX+1)*syst_p.at(sI)->GetBinContent(bIX+1));

      syst_p.at(sI)->SetBinContent(bIX+1, binVal);
      syst_p.at(sI)->SetBinError(bIX+1, tempRelErr);
      if(syst_p.at(sI)->GetBinContent(bIX+1) <= 0){
        syst_p.at(sI)->SetBinContent(bIX+1, smallDelta);
        syst_p.at(sI)->SetBinError(bIX+1, smallDelta);
      }
    }

    std::cout << "Pushing back: " << tempVal << ", " << nominal_p->GetBinLowEdge(bIX+1) << "-" << nominal_p->GetBinLowEdge(bIX+2) << std::endl;
    systVect.push_back(tempVal);
    nomValVect.push_back(nominal_p->GetBinContent(bIX+1));
  }

  double maxVal = -1;

  for(unsigned int sI = 0; sI < syst_p.size(); ++sI){
    if(systToBeSkipped[sI]) continue;
    if(syst_p.at(sI)->GetMaximum() > maxVal) maxVal = syst_p.at(sI)->GetMaximum();
  }

  const Int_t nMax = 2;
  const Double_t max[nMax] = {maxVal, 0.25};

  for(Int_t nI = 0; nI < nMax; ++nI){
    TLegend* leg_p = new TLegend(0.2, 0.52, 0.4, 0.99);
    leg_p->SetBorderSize(0);
    leg_p->SetFillColor(0);
    leg_p->SetFillStyle(0);
    leg_p->SetTextFont(43);
    leg_p->SetTextSize(16);

    TLatex* label_p = new TLatex();
    label_p->SetTextFont(43);
    label_p->SetTextSize(16);


    TCanvas* canv_p = new TCanvas("canv_p", "", 450, 450);
    canv_p->SetTopMargin(0.01);
    canv_p->SetRightMargin(0.001);
    canv_p->SetLeftMargin(0.12);
    canv_p->SetBottomMargin(0.12);

    TH1D* tempHist_p = new TH1D("tempHist_p", ";Jet p_{T} (GeV);Fractional Systetmatics", nBins, bins);
    tempHist_p->GetXaxis()->SetTitleOffset(1.6);
    tempHist_p->SetMarkerStyle(1);
    tempHist_p->SetMarkerSize(0.001);
    tempHist_p->SetMarkerColor(1);
    tempHist_p->SetLineColor(1);

    canv_p->cd();
    centerTitles(tempHist_p);
    tempHist_p->SetMaximum(1.5*max[nI]);
    tempHist_p->DrawCopy("HIST  P");

    gStyle->SetOptStat(0);
    gPad->SetLogx();

    double systMaxVal = -1;
    unsigned int systMaxPos = -1;
    double systSubMaxVal = -1;
    unsigned int systSubMaxPos = -1;

    for(unsigned int sI = 0; sI < syst_p.size(); ++sI){
      if(systToBeSkipped[sI]) continue;

      if(systMaxVal < syst_p.at(sI)->GetMaximum()){
        systSubMaxVal = systMaxVal;
        systSubMaxPos = systMaxPos;

        systMaxVal = syst_p.at(sI)->GetMaximum();
	systMaxPos = sI;
      }
      else if(systSubMaxVal < syst_p.at(sI)->GetMaximum()){
	systSubMaxVal = syst_p.at(sI)->GetMaximum();
        systSubMaxPos = sI;
      }
    }

    for(unsigned int sI = 0; sI < syst_p.size(); ++sI){
      if(systToBeSkipped[sI]) continue;

      Int_t binFillPos = 0;
      for(Int_t bIX = 0; bIX < syst_p.at(sI)->GetNbinsX(); ++bIX){
        Double_t binCenter = syst_p.at(sI)->GetBinCenter(bIX+1);
        if(binCenter > maxXVal || binCenter < minXVal) continue;

        Double_t binVal = tempHist_p->GetBinContent(binFillPos+1)*tempHist_p->GetBinContent(binFillPos+1);
        binVal += syst_p.at(sI)->GetBinContent(bIX+1)*syst_p.at(sI)->GetBinContent(bIX+1);
        binVal = TMath::Sqrt(binVal);

        tempHist_p->SetBinContent(binFillPos+1, binVal);
        tempHist_p->SetBinError(binFillPos+1, smallDelta);

        binFillPos++;
      }

      syst_p.at(sI)->SetMarkerColor(colors[sI%nColor]);
      syst_p.at(sI)->SetMarkerStyle(styles[sI%nStyle]);
      syst_p.at(sI)->SetLineColor(colors[sI%nColor]);
      syst_p.at(sI)->SetMarkerSize(0.8);

      syst_p.at(sI)->DrawCopy("HIST  P SAME");
      syst_p.at(sI)->DrawCopy("HIST  SAME");

      std::string legStr = systLegStr[sI];
      if(systMaxPos == sI) legStr = legStr + " (Dominant)";
      else if(systSubMaxPos == sI) legStr = legStr + " (Sub-dominant)";
      leg_p->AddEntry(syst_p.at(sI), legStr.c_str(), "P L");
    }

    tempHist_p->DrawCopy("HIST  SAME");

    if(nI == 0){
      std::cout << "Quick closure check: " << std::endl;
      for(Int_t bIX = 0; bIX < nBins; ++bIX){
	std::string binLowStr = prettyString(tempHist_p->GetBinLowEdge(bIX+1), 1, false);
	std::string binHiStr = prettyString(tempHist_p->GetBinLowEdge(bIX+2), 1, false);
      }
    }

    drawWhiteBox(900, 1100, -max[nI]*9./10., -max[nI]/100.);

    double yPos = 0 - max[nI]/20.;

    label_p->DrawLatex(200, yPos, "200");
    label_p->DrawLatex(400, yPos, "400");
    label_p->DrawLatex(600, yPos, "600");
    label_p->DrawLatex(800, yPos, "800");

    leg_p->Draw("SAME");

    std::string smoothStr = "Unsmoothed";
    if(doSmoothing) smoothStr = "Smoothed";

    std::string maxValStr = "NormalMax";
    if(nI != 0) maxValStr = "ZOOM";

    Int_t nLabelBins = 20;
    Double_t labelBins[nLabelBins+1];
    getLinBins(0, max[nI], nLabelBins, labelBins);

    label_p->DrawLatex(500, labelBins[19], maxValStr.c_str());
    label_p->DrawLatex(500, labelBins[17], smoothStr.c_str());

    gPad->RedrawAxis();
    gPad->SetTicks(1,2);


    std::string nominalName = nominal_p->GetName();
    std::string ppCentTag = "PP";
    if(nominalName.find("Cent") != std::string::npos){
      ppCentTag = nominalName.substr(nominalName.find("Cent"), nominalName.size());
      ppCentTag.replace(0, 4, "");
      ppCentTag.replace(ppCentTag.find("to"), 2, "-");
      ppCentTag.replace(ppCentTag.find("_"), ppCentTag.size(), "");
      ppCentTag = "PbPb " + ppCentTag + "%";
    }
    label_p->DrawLatex(500, labelBins[15], ppCentTag.c_str());

    std::string dirName = "pdfDir/" + dateStr;
    if(checkDir(dirStr)) dirName = dirStr;
    checkMakeDir(dirName);
    const std::string saveName = "systErr_" + nominalName + "_" + maxValStr + "_" + smoothStr + "_" + dateStr +  ".pdf";

    quietSaveAs(canv_p, (dirName + "/" + saveName));
    plotNames->push_back(saveName);

    delete tempHist_p;
    delete canv_p;
    delete leg_p;
    delete label_p;
  }

  return systVect;
}


void smoothErrors(TH1D* syst_p, TH1D* nom_p)
{
  const Int_t nBinsMax = 100;
  Double_t bins[nBinsMax];
  Int_t nBins = syst_p->GetNbinsX();

  for(Int_t bI = 0; bI < nBins+1; ++bI){
    bins[bI] = syst_p->GetBinLowEdge(bI+1);
  }

  TH1D* temp_p = new TH1D("temp_h", "", nBins, bins);

  for(Int_t bI = 0; bI < nBins; ++bI){
    temp_p->SetBinContent(bI+1, TMath::Abs(syst_p->GetBinContent(bI+1) - nom_p->GetBinContent(bI+1))/nom_p->GetBinContent(bI+1));
    temp_p->SetBinError(bI+1, 0);
  }

  temp_p->Smooth();

  for(Int_t bI = 0; bI < nBins;++bI){
    Double_t newVal = temp_p->GetBinContent(bI+1)*nom_p->GetBinContent(bI+1);
    if(nom_p->GetBinContent(bI+1) < syst_p->GetBinContent(bI+1)) newVal = nom_p->GetBinContent(bI+1) + newVal;
    else newVal = nom_p->GetBinContent(bI+1) - newVal;

    syst_p->SetBinContent(bI+1, newVal);
  }

  delete temp_p;
  
  return;
}


void drawSyst(TCanvas* canv_p, TPad* pad_p, TH1D* nominal_p, std::vector<double>* syst_, Double_t minXVal, Double_t maxXVal, bool skipZeros = false)
{
  canv_p->cd();
  if(pad_p != NULL) pad_p->cd();

  TBox* tempBox_p = new TBox();
  tempBox_p->SetFillColorAlpha(nominal_p->GetMarkerColor(), .25);

  std::string nominalName = nominal_p->GetName();

  Int_t binFillPos = 0;
  for(Int_t bIX = 0; bIX < nominal_p->GetNbinsX(); ++bIX){
    Double_t binCenter = nominal_p->GetBinCenter(bIX+1);
    Double_t binLowEdge = nominal_p->GetBinLowEdge(bIX+1);
    Double_t binHiEdge = nominal_p->GetBinLowEdge(bIX+2);
    Double_t binContent = nominal_p->GetBinContent(bIX+1);

    if(nominalName.find("rraa") != std::string::npos){
      binLowEdge += 0.035;
      binHiEdge -= 0.035;
    }

    if(binCenter < minXVal) continue;
    if(binCenter > maxXVal) continue;

    if(skipZeros && nominal_p->GetBinContent(bIX+1) <= 0.000001){
      ++binFillPos;
      continue;
    }

    tempBox_p->DrawBox(binLowEdge, binContent - syst_->at(binFillPos), binHiEdge, binContent + syst_->at(binFillPos));

    ++binFillPos;
  }

  delete tempBox_p;

  return;
}


void drawLumiOrTAA(TCanvas* canv_p, TH1D* hist_p, Double_t min, Double_t max, std::vector<double> syst, std::string lumiOrTAA, int pos = -1)
{

  kirchnerPalette kPalette;

  const Int_t nBins = 200;
  Double_t bins[nBins+1];
  if(canv_p->GetLogx()) getLogBins(min, max, nBins, bins);
  else getLinBins(min, max, nBins, bins);
  
  TBox* tempBox_p = new TBox();
  if(lumiOrTAA.find("Lumi") != std::string::npos){
    tempBox_p->SetFillColorAlpha(kPalette.getColor(getColorPosFromCent("", true)), .25);

    std::cout << "LUMIORTAA: " << bins[6] << "-" << bins[24] << ", " << min << ", " << max << std::endl;
    tempBox_p->DrawBox(bins[6], 1 - syst[0], bins[24], 1 + syst[0]);
  }
  else if(lumiOrTAA.find("TAA") != std::string::npos){
    tempBox_p->SetFillColorAlpha(hist_p->GetMarkerColor(), .25);
    tempBox_p->DrawBox(bins[26 + pos*20], 1 - syst[0], bins[44 + pos*20], 1 + syst[0]);
    std::cout << "LUMIORTAA: " << bins[26 + pos*20] << "-" << bins[44 + pos*20] << ", " << min << ", " << max << std::endl;
    std::cout << " LUMIORTAA: " << 26 + pos*20 << "-" << 44 + pos*20  << std::endl;
  }

  delete tempBox_p;

  return;
}


#endif
