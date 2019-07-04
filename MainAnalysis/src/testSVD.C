#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TDatime.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLatex.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TStyle.h"

#include "MainAnalysis/include/macroHistToSubsetHist.h"
#include "Utility/include/checkMakeDir.h"
#include "Utility/include/histDefUtility.h"
#include "Utility/include/kirchnerPalette.h"

#include "src/RooUnfoldBayes.h"
#include "src/RooUnfoldResponse.h"
#include "src/RooUnfoldSvd.h"

std::vector<double> stringToVectDouble(std::string inStr)
{
  std::vector<double> vect;

  while(inStr.find(" ") != std::string::npos){
    inStr.replace(inStr.find(" "), 1, "");
  }
 
  if(inStr.size() != 0){
    inStr = inStr + ",";
    while(inStr.find(",,") != std::string::npos){
      inStr.replace(inStr.find(",,"), 2, ",");
    }

    while(inStr.find(",") != std::string::npos){
      vect.push_back(std::stod(inStr.substr(0, inStr.find(","))));
      inStr.replace(0, inStr.find(",")+1, "");
    }
  } 
  
  return vect;
}

std::vector<std::string> stringToVectString(std::string inStr)
{
  while(inStr.find(" ") != std::string::npos){
    inStr.replace(inStr.find(" "), 1, "");
  }

  std::vector<std::string> vect;

  if(inStr.size() != 0){
    inStr = inStr + ",";
    while(inStr.find(",,") != std::string::npos){
      inStr.replace(inStr.find(",,"), 2, ",");
    }

    while(inStr.find(",") != std::string::npos){
      vect.push_back(inStr.substr(0, inStr.find(",")));
      inStr.replace(0, inStr.find(",")+1, "");
    }
  } 
  
  return vect;
}


void vectToArr(std::vector<double> vect, Double_t arr[])
{
  for(unsigned int vI = 0; vI < vect.size(); ++vI){
    arr[vI] = vect[vI];
  }
  
  return;
}

int testSVD(std::string inResponseName, std::string inDataName, std::string genBinsStr, std::string recoBinsStr, std::string inAlgoStr, std::string inCentStr)
{
  while(genBinsStr.find(" ") != std::string::npos){
    genBinsStr.replace(genBinsStr.find(" "), 1, "");
  }

  while(recoBinsStr.find(" ") != std::string::npos){
    recoBinsStr.replace(recoBinsStr.find(" "), 1, "");
  }

  TDatime* date = new TDatime();
  const std::string dateStr = std::to_string(date->GetDate());
  delete date;

  kirchnerPalette kPalette;
  const Int_t nColors = 6;

  const Int_t nStyles = 4;
  const Int_t styles[nStyles] = {21, 34, 47, 33};
  const Int_t styles2[nStyles] = {25, 28, 46, 27};
  
  std::vector<double> genBinsVect =stringToVectDouble(genBinsStr);
  std::vector<double> recoBinsVect =stringToVectDouble(recoBinsStr);

  const Int_t nMaxBins = 1000;
  Int_t nGenBins = genBinsVect.size()-1;  
  Int_t nRecoBins = recoBinsVect.size()-1;
  Double_t genBins[nMaxBins+1];
  Double_t genBins2[nMaxBins+1];
  Double_t recoBins[nMaxBins+1];
  vectToArr(genBinsVect, genBins);
  vectToArr(recoBinsVect, recoBins);

  for(Int_t gI = 1; gI < nGenBins; ++gI){
    genBins2[gI-1] = genBins[gI];
  }
  
  std::vector<std::string> algoStr;
  std::vector<std::string> centStr;
  std::vector<double> cParams;
  std::vector<double> sParams;
  std::vector<double> nParams;
  
  std::ifstream inFile("tables/csnParams_UPDATE_20190510.txt");
  std::string tempStr;

  int curPos = 0;
  int algoCentPos = -1;
  while(std::getline(inFile, tempStr)){
    std::vector<std::string> tempVect = stringToVectString(tempStr);
    if(tempVect[0].find("ALGO") != std::string::npos) continue;

    algoStr.push_back(tempVect[0]);
    centStr.push_back(tempVect[1]);
    cParams.push_back(std::stod(tempVect[2]));
    sParams.push_back(std::stod(tempVect[3]));
    nParams.push_back(std::stod(tempVect[4]));
    if(inAlgoStr.find(tempVect[0]) != std::string::npos){
      if(inCentStr.find(tempVect[1]) != std::string::npos){
	algoCentPos = curPos;
      }
    }
    
    ++curPos;
  } 
  inFile.close();

  Double_t cParam = cParams[algoCentPos];
  Double_t sParam = sParams[algoCentPos];
  Double_t nParam = nParams[algoCentPos];

  std::cout << inAlgoStr << ", " << inCentStr << ", " << cParam << ", " << sParam << ", " << nParam << std::endl;
  
  //  return 1;

  TFile* responseFile_p = new TFile(inResponseName.c_str(), "READ");
  TH2D* response_p = (TH2D*)responseFile_p->Get((inAlgoStr + "JetAnalyzer/response_" + inAlgoStr + "JetAnalyzer_PbPb_" + inCentStr + "_LightMUAndCHID_ResponseMod0p10_AbsEta0p0to2p0_General_h").c_str());

  Int_t nResBinsY = response_p->GetYaxis()->GetNbins();
  Int_t nResBinsX = response_p->GetXaxis()->GetNbins();
  Double_t resBinsY[nMaxBins+1];
  Double_t resBinsX[nMaxBins+1];
  for(Int_t rI = 0; rI < nResBinsX+1; ++rI){
    resBinsX[rI] = response_p->GetXaxis()->GetBinLowEdge(rI+1);
  }
  
  for(Int_t rI = 0; rI < nResBinsY+1; ++rI){
    resBinsY[rI] = response_p->GetYaxis()->GetBinLowEdge(rI+1);
  }

  std::vector<double> createBinComposite;
  for(Int_t rI = 0; rI < nResBinsY+1; ++rI){
    
    if(resBinsY[rI] >= 200 && ((Int_t)resBinsY[rI])%10 != 0)  continue;
    if(resBinsY[rI] >= 700 && resBinsY[rI] < 1000 &&  ((Int_t)resBinsY[rI])%20 != 0)  continue;
    if(resBinsY[rI] >= 1000 && ((Int_t)resBinsY[rI])%100 != 0)  continue;
    
    createBinComposite.push_back(resBinsY[rI]);
  }

  for(unsigned int cI = 0; cI < createBinComposite.size(); ++cI){
    std::cout << createBinComposite[cI] << ", ";
  }
  std::cout << std::endl;

  nResBinsY = createBinComposite.size()-1;
  nResBinsX = createBinComposite.size()-1;

  for(unsigned int cI = 0; cI < createBinComposite.size(); ++cI){
    resBinsX[cI] = createBinComposite[cI];
    resBinsY[cI] = createBinComposite[cI];
  }
  TH2D* responseCreate_p = new TH2D("responseCreate_h", "", nResBinsX, resBinsX, nResBinsY, resBinsY);
  macroHistToSubsetHist(response_p, responseCreate_p, true);

  TH1D* res1D_p[nMaxBins];
  for(Int_t rI = 0; rI < nResBinsY; ++rI){
    res1D_p[rI] = new TH1D(("res1D_" + std::to_string(rI) + "_h").c_str(), "", nResBinsX, resBinsX);

    for(Int_t bIX = 0; bIX < nResBinsX; ++bIX){
      res1D_p[rI]->SetBinContent(bIX+1, responseCreate_p->GetBinContent(bIX+1, rI+1));
      res1D_p[rI]->SetBinError(bIX+1, responseCreate_p->GetBinError(bIX+1, rI+1));
    }

    if(res1D_p[rI]->Integral() < TMath::Power(10, -20)){
      std::cout <<  "EMPTY: " << resBinsY[rI] << "-" << resBinsY[rI+1] << std::endl;
    }
  }

  
  TFile* dataFile_p = new TFile(inDataName.c_str(), "READ");
  TH1D* data_p = (TH1D*)dataFile_p->Get((inAlgoStr + "JetAnalyzer/jtPtRaw_General_" + inAlgoStr + "JetAnalyzer_PbPb_" + inCentStr + "_LightMUAndCHID_AbsEta0p0to2p0_h").c_str());
  TH1D* dataRed_p = new TH1D("dataRed_p", "", nRecoBins, recoBins);

  macroHistToSubsetHist(data_p, dataRed_p);
  
  const Int_t nSVD = 20;
  const Int_t nToy = 50;

  checkMakeDir("output");
  checkMakeDir("output/" + dateStr);
  
  TFile* outFile_p = new TFile(("output/" + dateStr + "/svd_" + inAlgoStr + "_" + inCentStr + ".root").c_str(), "RECREATE");
  TH1D* genHist_p[nToy];
  //  TH1D* genHistDiv_p = new TH1D("genHistDiv_h", "", nGenBins, genBins);
  TH1D* recoHist_p[nToy];
  TH1D* unfoldHist_p[nToy][nSVD];

  TH1D* unfoldChi2Prob_h[nSVD];
  
  //  genHistDiv_p->Sumw2();
  
  Int_t nSVDTrue = TMath::Min(nSVD, nGenBins);

  TH1D* unfoldChi2Prob_Mean_h = new TH1D("unfoldChi2Prob_Mean_h", ";kReg;#LTProb#GT", nSVDTrue, 0.5, nSVDTrue+0.5);
  TH1D* unfoldChi2Prob_Sigma_h = new TH1D("unfoldChi2Prob_Sigma_h", ";kReg;#sigma(Prob)", nSVDTrue, 0.5, nSVDTrue+0.5);

  TH1D* dVect_Mean_p = new TH1D("dVect_Mean_p", ";kReg;#LTDVal#GT", nSVDTrue, 0.5, nSVDTrue+0.5);
  TH1D* dVect_Sigma_p = new TH1D("dVect_Sigma_p", ";kReg;#sigma(DVal)", nSVDTrue, 0.5, nSVDTrue+0.5);

  std::vector<std::vector<double> > dVectVals;
  std::vector<double> dVectMeans;
  std::vector<double> dVectSigmas;

  for(Int_t dI = 0; dI < nSVDTrue; ++dI){
    dVectVals.push_back({});
    dVectMeans.push_back(0.0);
    dVectSigmas.push_back(0.0);
  }
  
  centerTitles({unfoldChi2Prob_Mean_h, unfoldChi2Prob_Sigma_h, dVect_Mean_p, dVect_Sigma_p});
  
   
  for(Int_t tI = 0; tI < nToy; ++tI){
    genHist_p[tI] = new TH1D(("genHist_Toy" + std::to_string(tI) + "_h").c_str(), "", nGenBins, genBins);
    recoHist_p[tI] = new TH1D(("recoHist_Toy" + std::to_string(tI) + "_h").c_str(), "", nRecoBins, recoBins);

    for(Int_t sI = 0; sI < nSVDTrue; ++sI){  
      unfoldHist_p[tI][sI] = new TH1D(("unfoldHist_Toy" + std::to_string(tI) + "_SVD" + std::to_string(sI+1) + "_h").c_str(), "", nGenBins, genBins);

      if(tI == 0){
	unfoldChi2Prob_h[sI] = new TH1D(("unfoldChi2Prob_Svd" + std::to_string(sI+1) + "_h").c_str(), "", 25, 0.0, 1.0);
      }
    }
  }
  
  TRandom3* randGen_p = new TRandom3(0);


  TH2D* randHist_p = new TH2D("randHist_h", "", nRecoBins, recoBins, nGenBins, genBins);
  TH1D* randGenHist_p = new TH1D("randGenHist_h", "", nGenBins, genBins);
  TH1D* randGenHistDiv_p = new TH1D("randGenHistDiv_h", "", nGenBins, genBins);
  randHist_p->Sumw2();
  randGenHist_p->Sumw2();
  randGenHistDiv_p->Sumw2();
  Int_t nFill = 0;


  Double_t targetPower = 6.5;
  if(inAlgoStr.find("akCs2") != std::string::npos){
    if(inCentStr.find("Cent50to90") != std::string::npos) targetPower = 6.3;
    else if(inCentStr.find("Cent30to50") != std::string::npos) targetPower = 6.4;
    else if(inCentStr.find("Cent10to30") != std::string::npos) targetPower = 6.5;
    else if(inCentStr.find("Cent0to10") != std::string::npos) targetPower = 6.9;
  }
  else if(inAlgoStr.find("akCs3") != std::string::npos){
    if(inCentStr.find("Cent50to90") != std::string::npos) targetPower = 6.3;
    else if(inCentStr.find("Cent30to50") != std::string::npos) targetPower = 6.4;
    else if(inCentStr.find("Cent10to30") != std::string::npos) targetPower = 6.5;
    else if(inCentStr.find("Cent0to10") != std::string::npos) targetPower = 6.9;
  }
  else if(inAlgoStr.find("akCs4") != std::string::npos){
    if(inCentStr.find("Cent50to90") != std::string::npos) targetPower = 6.4;
    else if(inCentStr.find("Cent30to50") != std::string::npos) targetPower = 6.5;
    else if(inCentStr.find("Cent10to30") != std::string::npos) targetPower = 6.7;
    else if(inCentStr.find("Cent0to10") != std::string::npos) targetPower = 6.9;
  }
  else if(inAlgoStr.find("akCs6") != std::string::npos){
    if(inCentStr.find("Cent50to90") != std::string::npos) targetPower = 6.4;
    else if(inCentStr.find("Cent30to50") != std::string::npos) targetPower = 6.5;
    else if(inCentStr.find("Cent10to30") != std::string::npos) targetPower = 6.7;
    else if(inCentStr.find("Cent0to10") != std::string::npos) targetPower = 6.9;
  }
  else if(inAlgoStr.find("akCs8") != std::string::npos){
    if(inCentStr.find("Cent50to90") != std::string::npos) targetPower = 6.6;
    else if(inCentStr.find("Cent30to50") != std::string::npos) targetPower = 6.6;
    else if(inCentStr.find("Cent10to30") != std::string::npos) targetPower = 6.9;
    else if(inCentStr.find("Cent0to10") != std::string::npos) targetPower = 7.1;
  }
  else if(inAlgoStr.find("akCs10") != std::string::npos){
    if(inCentStr.find("Cent50to90") != std::string::npos) targetPower = 6.6;
    else if(inCentStr.find("Cent30to50") != std::string::npos) targetPower = 6.6;
    else if(inCentStr.find("Cent10to30") != std::string::npos) targetPower = 6.9;
    else if(inCentStr.find("Cent0to10") != std::string::npos) targetPower = 7.1;
  }
  std::cout << "PARAMS: " << cParam << ", " << sParam << ", " << nParam << std::endl;
  std::cout << "POWERS: " << targetPower << std::endl;

  Double_t weightPower = targetPower - 1 - 2;
  
  const Int_t nTarget = 1000000;
  while(nFill < nTarget){    
    Double_t val = randGen_p->Uniform(0,1);
    Double_t genPt = genBins[0]/TMath::Power(val, 1./2);   
    Double_t weight = TMath::Power(genBins[0]/genPt, weightPower);
    
    if(genPt < genBins[0]) continue;
    if(genPt >= genBins[nGenBins]) continue;


    Double_t res = TMath::Sqrt(cParam*cParam + sParam*sParam/genPt + nParam*nParam/(genPt*genPt));

    Double_t recoPt = genPt*randGen_p->Gaus(1., res);

    
    Int_t pos = -1;
    for(Int_t rI = 0; rI < nResBinsY; ++rI){
      if(genPt >= resBinsY[rI] && genPt < resBinsY[rI+1]){
	pos = rI;
	break;
      }
    }
    //    std::cout << "POS: " << pos << ", " << genPt << ", " << resBinsY[pos] << "-" << resBinsY[pos+1] << std::endl;
    recoPt = res1D_p[pos]->GetRandom();


    if(recoPt < recoBins[0]) continue;
    if(recoPt >= recoBins[nRecoBins]) continue;


    //    std::cout << recoPt << ", " << genPt << std::endl;
    randHist_p->Fill(recoPt, genPt, weight);
    randGenHist_p->Fill(genPt, weight);
    randGenHistDiv_p->Fill(genPt, weight);

    if(nFill%10000 == 0){
      std::cout << "nFill: " << nFill << "/" << nTarget << std::endl;
    }
    
    ++nFill;
  }
  
  for(Int_t tI = 0; tI < nToy; ++tI){
    while(genHist_p[tI]->Integral() < dataRed_p->Integral()){    
      Double_t genPt, recoPt;
      randHist_p->GetRandom2(recoPt, genPt);    
      
      genHist_p[tI]->Fill(genPt);
      recoHist_p[tI]->Fill(recoPt);      
    }
  }
 
  for(Int_t tI = 0; tI < nToy; ++tI){
    std::cout << "UNFOLDING TOY: " << tI << "/" << nToy << std::endl;
    for(Int_t sI = 0; sI < nSVDTrue; ++sI){
      TH2D* tempRes_p = new TH2D("tempRes_p", "", nRecoBins, recoBins, nGenBins, genBins);      
      TH1D* tempGen_p = new TH1D("tempGen_p", "", nGenBins, genBins);
      TH1D* tempReco_p = new TH1D("tempReco_p", "", nRecoBins, recoBins);
      macroHistToSubsetHist(response_p, tempRes_p, true);
      
      Double_t integral = tempRes_p->Integral();
      for(Int_t bIY = 0; bIY < tempRes_p->GetYaxis()->GetNbins(); ++bIY){
	Double_t total = 0.0;

	for(Int_t bIX = 0; bIX < tempRes_p->GetXaxis()->GetNbins(); ++bIX){
	  total += tempRes_p->GetBinContent(bIX+1, bIY+1);
	}

	for(Int_t bIX = 0; bIX < tempRes_p->GetXaxis()->GetNbins(); ++bIX){
	  Double_t val = tempRes_p->GetBinContent(bIX+1, bIY+1)/total;
	  Double_t err = tempRes_p->GetBinError(bIX+1, bIY+1)/total;
	  tempRes_p->SetBinContent(bIX+1, bIY+1, val);
	  tempRes_p->SetBinError(bIX+1, bIY+1, err);
	}	
      }
      tempRes_p->Scale(integral/tempRes_p->Integral());

      macroHistToSubsetHistX(tempRes_p, tempReco_p, true);
      macroHistToSubsetHistY(tempRes_p, tempGen_p, true);
      
      RooUnfoldResponse* rooRes_p = new RooUnfoldResponse(tempReco_p, tempGen_p, tempRes_p, "rooRes_p");
      RooUnfoldSvd* rooSvd_p = new RooUnfoldSvd(rooRes_p, recoHist_p[tI], TMath::Min(nSVDTrue-1, TMath::Max(3, nSVDTrue-2)));      
      TH1D* tempUnfold_p = (TH1D*)rooSvd_p->Hreco(RooUnfold::kCovToy);
      
      delete rooRes_p;
      delete rooSvd_p;
      
      for(Int_t bIY = 0; bIY < tempRes_p->GetYaxis()->GetNbins(); ++bIY){
	Double_t total = 0.0;

        for(Int_t bIX = 0; bIX < tempRes_p->GetXaxis()->GetNbins(); ++bIX){
          total += tempRes_p->GetBinContent(bIX+1, bIY+1);
	}

	Double_t factor = tempUnfold_p->GetBinContent(bIY+1)/total;

	for(Int_t bIX = 0; bIX < tempRes_p->GetXaxis()->GetNbins(); ++bIX){
          Double_t val = tempRes_p->GetBinContent(bIX+1, bIY+1)*factor;
          Double_t err = tempRes_p->GetBinError(bIX+1, bIY+1)*factor;

          tempRes_p->SetBinContent(bIX+1, bIY+1, val);
          tempRes_p->SetBinError(bIX+1, bIY+1, err);
        }
      }

      rooRes_p = new RooUnfoldResponse(tempReco_p, tempGen_p, tempRes_p, "rooRes_p");
      rooSvd_p = new RooUnfoldSvd(rooRes_p, recoHist_p[tI], sI+1);
      tempUnfold_p = (TH1D*)rooSvd_p->Hreco(RooUnfold::kCovToy);
      

      for(Int_t gI = 0; gI < tempUnfold_p->GetNbinsX(); ++gI){
	unfoldHist_p[tI][sI]->SetBinContent(gI+1, tempUnfold_p->GetBinContent(gI+1));
	unfoldHist_p[tI][sI]->SetBinError(gI+1, tempUnfold_p->GetBinError(gI+1));
      }

      if(sI == 0){
	TH1D* tempDVect_p = (TH1D*)rooSvd_p->GetD();
	for(Int_t bIX = 0; bIX < tempDVect_p->GetNbinsX(); ++bIX){
	  Double_t val = dVect_Mean_p->GetBinContent(bIX+1) + tempDVect_p->GetBinContent(bIX+1);
	  dVect_Mean_p->SetBinContent(bIX+1, val);
	  dVect_Mean_p->SetBinError(bIX+1, 0.0);

	  dVectVals[bIX].push_back(tempDVect_p->GetBinContent(bIX+1));
	  dVectMeans[bIX] += tempDVect_p->GetBinContent(bIX+1)/((Double_t)nToy);
	}	
      }
      
      TH1D* temp2_p = new TH1D("temp2_h", "", nGenBins-2, genBins2);
      TH1D* temp3_p = new TH1D("temp3_h", "", nGenBins-2, genBins2);
      macroHistToSubsetHist(genHist_p[tI], temp2_p);
      macroHistToSubsetHist(tempUnfold_p, temp3_p);
      
      Double_t prob = temp3_p->Chi2Test(temp2_p, "WW");
      
      unfoldChi2Prob_h[sI]->Fill(prob);

      delete temp2_p;
      delete temp3_p;
      
      delete rooSvd_p;
      delete rooRes_p;
      delete tempReco_p;
      delete tempGen_p;
      delete tempRes_p;
    }
  }

  
  delete randGen_p;
  

  for(Int_t sI = 0; sI < nSVDTrue; ++sI){  
    unfoldChi2Prob_Mean_h->SetBinContent(sI+1, unfoldChi2Prob_h[sI]->GetMean());
    unfoldChi2Prob_Mean_h->SetBinError(sI+1, unfoldChi2Prob_h[sI]->GetMeanError());

    unfoldChi2Prob_Sigma_h->SetBinContent(sI+1, unfoldChi2Prob_h[sI]->GetStdDev());
    unfoldChi2Prob_Sigma_h->SetBinError(sI+1, unfoldChi2Prob_h[sI]->GetStdDevError());
  }

  outFile_p->cd();


  TCanvas* canv_p = new TCanvas("canv_p", "", 450*2, 450);
  canv_p->SetTopMargin(0.01);
  canv_p->SetRightMargin(0.01);
  canv_p->SetLeftMargin(0.01);
  canv_p->SetBottomMargin(0.01);

  canv_p->Divide(2, 1);

  canv_p->cd(1);
  gPad->SetTopMargin(0.01);
  gPad->SetRightMargin(0.01);
  gPad->SetLeftMargin(0.14);
  gPad->SetBottomMargin(0.14);
  
  checkMakeDir("pdfDir");
  checkMakeDir("pdfDir/" + dateStr);

  unfoldChi2Prob_Mean_h->SetMarkerStyle(20);
  unfoldChi2Prob_Mean_h->SetMarkerSize(1);
  unfoldChi2Prob_Mean_h->SetMarkerColor(1);
  unfoldChi2Prob_Mean_h->SetLineColor(1);

  unfoldChi2Prob_Sigma_h->SetMarkerStyle(21);
  unfoldChi2Prob_Sigma_h->SetMarkerSize(1);
  unfoldChi2Prob_Sigma_h->SetMarkerColor(4);
  unfoldChi2Prob_Sigma_h->SetLineColor(4);

  unfoldChi2Prob_Mean_h->GetYaxis()->SetTitle("#LTProb#GT or #sigma(Prob)");

  unfoldChi2Prob_Mean_h->SetMaximum(1.0);
  unfoldChi2Prob_Mean_h->SetMinimum(0.0);
  unfoldChi2Prob_Mean_h->DrawCopy("HIST E1 P");
  unfoldChi2Prob_Sigma_h->DrawCopy("HIST E1 P SAME");
  gStyle->SetOptStat(0);

  TLatex* label_p = new TLatex();
  label_p->SetNDC();
  label_p->SetTextSize(12);
  label_p->SetTextFont(43);

  label_p->DrawLatex(0.6, 0.9, inAlgoStr.c_str());
  label_p->DrawLatex(0.6, 0.8, inCentStr.c_str());

  canv_p->cd(2);
  gPad->SetTopMargin(0.01);
  gPad->SetRightMargin(0.01);
  gPad->SetLeftMargin(0.14);
  gPad->SetBottomMargin(0.14);  
  
  for(Int_t bIX = 0; bIX < dVect_Mean_p->GetNbinsX(); ++bIX){
    dVect_Mean_p->SetBinContent(bIX+1, dVect_Mean_p->GetBinContent(bIX+1)/((Double_t)nToy));
    std::cout << dVectMeans[bIX] << "/" << dVect_Mean_p->GetBinContent(bIX+1);
    
    for(unsigned int sI = 0; sI < dVectVals[bIX].size(); ++sI){
      Double_t val = dVectVals[bIX][sI] - dVectMeans[bIX];
      val *= val;

      dVectSigmas[bIX] += val;      
    }
    dVectSigmas[bIX] /= ((Double_t)nToy-1.);

    dVectSigmas[bIX] = TMath::Sqrt(dVectSigmas[bIX]);

    dVect_Mean_p->SetBinError(bIX+1, dVectSigmas[bIX]/TMath::Sqrt((Double_t)nToy));    

    dVect_Sigma_p->SetBinContent(bIX+1, dVectSigmas[bIX]);    
  }

  dVect_Mean_p->SetMarkerStyle(20);
  dVect_Mean_p->SetMarkerSize(1);
  dVect_Mean_p->SetMarkerColor(1);
  dVect_Mean_p->SetLineColor(1);

  dVect_Sigma_p->SetMarkerStyle(21);
  dVect_Sigma_p->SetMarkerSize(1);
  dVect_Sigma_p->SetMarkerColor(4);
  dVect_Sigma_p->SetLineColor(4);

  dVect_Mean_p->GetYaxis()->SetTitle("#LTDVal#GT or #sigma(DVal)");
  dVect_Mean_p->SetMinimum(0.5);
  dVect_Mean_p->SetMaximum(200.);
  dVect_Mean_p->DrawCopy("HIST E1 P");
  dVect_Sigma_p->DrawCopy("HIST E1 P SAME");

  label_p->DrawLatex(0.3, 0.9, genBinsStr.c_str());
  label_p->DrawLatex(0.3, 0.8, recoBinsStr.c_str());
   

  dVect_Mean_p->GetYaxis()->SetTitle("#LTDVal#GT");

  gPad->SetLogy();
  gStyle->SetOptStat(0);
  
  canv_p->SaveAs(("pdfDir/" + dateStr + "/bestSVD_" + inAlgoStr + "_" + inCentStr + "_" + dateStr + ".pdf").c_str());
  delete canv_p;


  std::cout << __LINE__ << std::endl;
  for(Int_t sI = 0; sI < nSVDTrue; ++sI){ 
    std::cout << __LINE__ << std::endl;
    canv_p = new TCanvas("canv_p", "", 450*3, 450);
    canv_p->SetTopMargin(0.01);
    canv_p->SetRightMargin(0.01);
    canv_p->SetLeftMargin(0.01);
    canv_p->SetBottomMargin(0.01);

    std::cout << __LINE__ << std::endl;

    canv_p->Divide(3, 1);
    
    canv_p->cd();
    canv_p->cd(1);
    gPad->SetTopMargin(0.01);
    gPad->SetRightMargin(0.01);
    gPad->SetLeftMargin(0.14);
    gPad->SetBottomMargin(0.14);

    std::cout << __LINE__ << std::endl;

    for(Int_t tI = 0; tI < TMath::Min(1, nToy); ++tI){
      TH1D* temp_p = new TH1D("temp_h", "", nGenBins, genBins);
      macroHistToSubsetHist(genHist_p[tI], temp_p);
      temp_p->Divide(unfoldHist_p[tI][sI]);

      temp_p->SetMarkerStyle(styles[tI%nStyles]);
      temp_p->SetMarkerSize(1);
      temp_p->SetMarkerColor(kPalette.getColor(tI%nColors));
      temp_p->SetLineColor(kPalette.getColor(tI%nColors));

      temp_p->SetMaximum(2.0);
      temp_p->SetMinimum(0.0);

      if(tI == 0) temp_p->DrawCopy("HIST E1 P");
      else temp_p->DrawCopy("HIST E1 P SAME");

      delete temp_p;
    }
  
    canv_p->cd(1);
    label_p->DrawLatex(0.6, 0.9, inAlgoStr.c_str());
    label_p->DrawLatex(0.6, 0.82, inCentStr.c_str());
    label_p->DrawLatex(0.6, 0.74, ("SVD=" + std::to_string(sI+1)).c_str());

    canv_p->cd(2);
    gPad->SetTopMargin(0.01);
    gPad->SetRightMargin(0.01);
    gPad->SetLeftMargin(0.14);
    gPad->SetBottomMargin(0.14);

    dataRed_p->SetMarkerStyle(20);
    dataRed_p->SetLineColor(1);
    dataRed_p->SetMarkerColor(1);
    dataRed_p->SetMarkerSize(1);
  
    if(sI == 0){
      dataRed_p->SetMaximum(dataRed_p->GetMaximum()*2.);
      dataRed_p->SetMinimum(0.5);
    }
    dataRed_p->DrawCopy("HIST E1 P");
  
    for(Int_t tI = 0; tI < TMath::Min(1, nToy); ++tI){
      recoHist_p[tI]->SetMarkerStyle(styles[tI%nStyles]);
      recoHist_p[tI]->SetMarkerSize(1);
      recoHist_p[tI]->SetMarkerColor(kPalette.getColor(tI%nColors));
      recoHist_p[tI]->SetLineColor(kPalette.getColor(tI%nColors));

      recoHist_p[tI]->DrawCopy("HIST E1 P SAME");

      genHist_p[tI]->SetMarkerStyle(styles2[tI%nStyles]);
      genHist_p[tI]->SetMarkerSize(1);
      genHist_p[tI]->SetMarkerColor(kPalette.getColor(tI%nColors));
      genHist_p[tI]->SetLineColor(kPalette.getColor(tI%nColors));

      genHist_p[tI]->DrawCopy("HIST E1 P SAME");
    }
    gPad->SetLogy();
    dataRed_p->DrawCopy("HIST E1 P SAME");

    canv_p->cd(3);
    gPad->SetTopMargin(0.01);
    gPad->SetRightMargin(0.01);
    gPad->SetLeftMargin(0.14);
    gPad->SetBottomMargin(0.14);

    for(Int_t tI = 0; tI < TMath::Min(3, nToy); ++tI){
      TH1D* temp_p = new TH1D("temp_h", "", nRecoBins, recoBins);
      macroHistToSubsetHist(recoHist_p[tI], temp_p);
      

      temp_p->SetMarkerStyle(styles[tI%nStyles]);
      temp_p->SetMarkerSize(1);
      temp_p->SetMarkerColor(kPalette.getColor(tI%nColors));
      temp_p->SetLineColor(kPalette.getColor(tI%nColors));
    
      temp_p->Divide(dataRed_p);
      temp_p->SetMaximum(2.0);
      temp_p->SetMinimum(0.5);

      if(tI == 0) temp_p->DrawCopy("HIST E1 P");
      else temp_p->DrawCopy("HIST E1 P SAME");

      delete temp_p;
    }

    std::cout << __LINE__ << std::endl;
    canv_p->SaveAs(("pdfDir/" + dateStr + "/toyRats_" + inAlgoStr + "_" + inCentStr + "_SVD" + std::to_string(sI+1) + "_" + dateStr + ".pdf").c_str());
    std::cout << __LINE__ << std::endl; 
    delete canv_p;
    std::cout << __LINE__ << std::endl;
  }

  delete label_p;

  std::cout << __LINE__ << std::endl;

  
  unfoldChi2Prob_Mean_h->GetYaxis()->SetTitle("#LTProb#GT");

  unfoldChi2Prob_Mean_h->Write("", TObject::kOverwrite);
  delete unfoldChi2Prob_Mean_h;

  unfoldChi2Prob_Sigma_h->Write("", TObject::kOverwrite);
  delete unfoldChi2Prob_Sigma_h;

  dVect_Mean_p->Write("", TObject::kOverwrite);
  delete dVect_Mean_p;

  dVect_Sigma_p->Write("", TObject::kOverwrite);
  delete dVect_Sigma_p;

  delete dataRed_p;
  
  randHist_p->Write("", TObject::kOverwrite);
  randGenHist_p->Write("", TObject::kOverwrite);

  for(Int_t gI = 0; gI < nGenBins; ++gI){
    Double_t width = randGenHistDiv_p->GetBinWidth(gI+1);
    Double_t val = randGenHistDiv_p->GetBinContent(gI+1)/width;
    Double_t err = randGenHistDiv_p->GetBinError(gI+1)/width;

    randGenHistDiv_p->SetBinContent(gI+1, val);
    randGenHistDiv_p->SetBinError(gI+1, err);
  }
    std::cout << __LINE__ << std::endl;
  randGenHistDiv_p->Write("", TObject::kOverwrite);
  delete randGenHistDiv_p;
    std::cout << __LINE__ << std::endl;

  
  for(Int_t tI = 0; tI < nToy; ++tI){
    genHist_p[tI]->Write("", TObject::kOverwrite);
    recoHist_p[tI]->Write("", TObject::kOverwrite);
  }
    std::cout << __LINE__ << std::endl;
  
  for(Int_t tI = 0; tI < nToy; ++tI){
    for(Int_t sI = 0; sI < nSVDTrue; ++sI){
      unfoldHist_p[tI][sI]->Write("", TObject::kOverwrite);
      delete unfoldHist_p[tI][sI];
    }
  }
    std::cout << __LINE__ << std::endl;
  
  for(Int_t tI = 0; tI < nToy; ++tI){   
    delete genHist_p[tI];
    delete recoHist_p[tI];
  }
    std::cout << __LINE__ << std::endl;

  for(Int_t sI = 0; sI < nSVDTrue; ++sI){  
    unfoldChi2Prob_h[sI]->Write("", TObject::kOverwrite);
    delete unfoldChi2Prob_h[sI];
  }
    std::cout << __LINE__ << std::endl;

  outFile_p->Close();
  delete outFile_p;
  
  dataFile_p->Close();
  delete dataFile_p;

  responseFile_p->Close();
  delete responseFile_p;

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 7){
    std::cout << "Usage: ./bin/testSVD.exe <inResponseName> <inDataName> <genBins> <recoBins> <inAlgoStr> <inCentStr>" << std::endl;
    return 1;
  }

  int retVal = 0;;
  retVal += testSVD(argv[1], argv[2], argv[3], argv[4], argv[5], argv[6]);
  return retVal;
}
