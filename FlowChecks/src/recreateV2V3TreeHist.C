#include <iostream>
#include <string>
#include <vector>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TMath.h"
#include "TDatime.h"
#include "TF1.h"
#include "Math/Minimizer.h"

#include "Utility/include/cppWatch.h"
#include "Utility/include/doGlobalDebug.h"

int recreateV2V3TreeHist(const std::string inFileName)
{
  TDatime* date = new TDatime();
  const std::string dateStr = std::to_string(date->GetDate());
  delete date;

  TFile* outFile_p = new TFile(("output/v2CrossCheck_TreeHist_" + dateStr + ".root").c_str(), "RECREATE");

  const Int_t nCentBins = 11;
  const Int_t centBinsLow[nCentBins] = {5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55};
  const Int_t centBinsHi[nCentBins] = {10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60};
  const Double_t centBins[nCentBins+1] = {5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60};

  TH1F* v2Raw_PF_h[nCentBins];
  TH1F* v2RawCorr_PF_h[nCentBins];

  TH1F* v2Obs_PF_h[nCentBins];
  TH1F* v2ObsCorr_PF_h[nCentBins];

  TH1F* v2Raw_Trk_h[nCentBins];
  TH1F* v2RawCorr_Trk_h[nCentBins];
  
  TH1F* v2Obs_Trk_h[nCentBins];
  TH1F* v2ObsCorr_Trk_h[nCentBins];

  TH1F* v2Raw_EByE_h[nCentBins];
  TH1F* v2RawCorr_EByE_h[nCentBins];  
  TH1F* v2Obs_EByE_h[nCentBins];
  TH1F* v2ObsCorr_EByE_h[nCentBins];
  TH1F* v2Fit_EByE_h[nCentBins];
  TH1F* v2FitCorr_EByE_h[nCentBins];  
  TH1F* v2FitV4_EByE_h[nCentBins];
  TH1F* v2FitV4Corr_EByE_h[nCentBins];  


  TH1F* v2Raw_Mean_PF_h = new TH1F("v2Raw_Mean_PF_h", ";Centrality (%);#LTv_{2}^{raw}#GT", nCentBins, centBins);
  TH1F* v2Raw_Sigma_PF_h = new TH1F("v2Raw_Sigma_PF_h", ";Centrality (%);#sigma(v_{2}^{raw})", nCentBins, centBins);
  TH1F* v2RawCorr_Mean_PF_h = new TH1F("v2RawCorr_Mean_PF_h", ";Centrality (%);#LTv_{2}^{raw}#GT", nCentBins, centBins);
  TH1F* v2RawCorr_Sigma_PF_h = new TH1F("v2RawCorr_Sigma_PF_h", ";Centrality (%);#sigma(v_{2}^{raw})", nCentBins, centBins);

  TH1F* v2Obs_Mean_PF_h = new TH1F("v2Obs_Mean_PF_h", ";Centrality (%);#LTv_{2}^{obs}#GT", nCentBins, centBins);
  TH1F* v2Obs_Sigma_PF_h = new TH1F("v2Obs_Sigma_PF_h", ";Centrality (%);#sigma(v_{2}^{obs})", nCentBins, centBins);
  TH1F* v2ObsCorr_Mean_PF_h = new TH1F("v2ObsCorr_Mean_PF_h", ";Centrality (%);#LTv_{2}^{obs}#GT", nCentBins, centBins);
  TH1F* v2ObsCorr_Sigma_PF_h = new TH1F("v2ObsCorr_Sigma_PF_h", ";Centrality (%);#sigma(v_{2}^{obs})", nCentBins, centBins);

  TH1F* v2Raw_Mean_Trk_h = new TH1F("v2Raw_Mean_Trk_h", ";Centrality (%);#LTv_{2}^{raw}#GT", nCentBins, centBins);
  TH1F* v2Raw_Sigma_Trk_h = new TH1F("v2Raw_Sigma_Trk_h", ";Centrality (%);#sigma(v_{2}^{raw})", nCentBins, centBins);
  TH1F* v2RawCorr_Mean_Trk_h = new TH1F("v2RawCorr_Mean_Trk_h", ";Centrality (%);#LTv_{2}^{raw}#GT", nCentBins, centBins);
  TH1F* v2RawCorr_Sigma_Trk_h = new TH1F("v2RawCorr_Sigma_Trk_h", ";Centrality (%);#sigma(v_{2}^{raw})", nCentBins, centBins);

  TH1F* v2Obs_Mean_Trk_h = new TH1F("v2Obs_Mean_Trk_h", ";Centrality (%);#LTv_{2}^{obs}#GT", nCentBins, centBins);
  TH1F* v2Obs_Sigma_Trk_h = new TH1F("v2Obs_Sigma_Trk_h", ";Centrality (%);#sigma(v_{2}^{obs})", nCentBins, centBins);
  TH1F* v2ObsCorr_Mean_Trk_h = new TH1F("v2ObsCorr_Mean_Trk_h", ";Centrality (%);#LTv_{2}^{obs}#GT", nCentBins, centBins);
  TH1F* v2ObsCorr_Sigma_Trk_h = new TH1F("v2ObsCorr_Sigma_Trk_h", ";Centrality (%);#sigma(v_{2}^{obs})", nCentBins, centBins);


  TH1F* v2Raw_Mean_EByE_h = new TH1F("v2Raw_Mean_EByE_h", ";Centrality (%);#LTv_{2}^{raw}#GT", nCentBins, centBins);
  TH1F* v2Raw_Sigma_EByE_h = new TH1F("v2Raw_Sigma_EByE_h", ";Centrality (%);#sigma(v_{2}^{raw})", nCentBins, centBins);
  TH1F* v2RawCorr_Mean_EByE_h = new TH1F("v2RawCorr_Mean_EByE_h", ";Centrality (%);#LTv_{2}^{raw}#GT", nCentBins, centBins);
  TH1F* v2RawCorr_Sigma_EByE_h = new TH1F("v2RawCorr_Sigma_EByE_h", ";Centrality (%);#sigma(v_{2}^{raw})", nCentBins, centBins);

  TH1F* v2Obs_Mean_EByE_h = new TH1F("v2Obs_Mean_EByE_h", ";Centrality (%);#LTv_{2}^{obs}#GT", nCentBins, centBins);
  TH1F* v2Obs_Sigma_EByE_h = new TH1F("v2Obs_Sigma_EByE_h", ";Centrality (%);#sigma(v_{2}^{obs})", nCentBins, centBins);
  TH1F* v2ObsCorr_Mean_EByE_h = new TH1F("v2ObsCorr_Mean_EByE_h", ";Centrality (%);#LTv_{2}^{obs}#GT", nCentBins, centBins);
  TH1F* v2ObsCorr_Sigma_EByE_h = new TH1F("v2ObsCorr_Sigma_EByE_h", ";Centrality (%);#sigma(v_{2}^{obs})", nCentBins, centBins);

  TH1F* v2Fit_Mean_EByE_h = new TH1F("v2Fit_Mean_EByE_h", ";Centrality (%);#LTv_{2}^{raw}#GT", nCentBins, centBins);
  TH1F* v2Fit_Sigma_EByE_h = new TH1F("v2Fit_Sigma_EByE_h", ";Centrality (%);#sigma(v_{2}^{raw})", nCentBins, centBins);
  TH1F* v2FitCorr_Mean_EByE_h = new TH1F("v2FitCorr_Mean_EByE_h", ";Centrality (%);#LTv_{2}^{raw}#GT", nCentBins, centBins);
  TH1F* v2FitCorr_Sigma_EByE_h = new TH1F("v2FitCorr_Sigma_EByE_h", ";Centrality (%);#sigma(v_{2}^{raw})", nCentBins, centBins);

  
  for(Int_t cI = 0; cI < nCentBins; ++cI){
    const std::string centStr = "Cent" + std::to_string(centBinsLow[cI]) + "to" + std::to_string(centBinsHi[cI]);
    v2Raw_PF_h[cI] = new TH1F(("v2Raw_" + centStr + "_PF_h").c_str(), ";v_{2};Counts", 150, 0.0, 0.6);
    v2RawCorr_PF_h[cI] = new TH1F(("v2RawCorr_" + centStr + "_PF_h").c_str(), ";v_{2};Counts", 150, 0.0, 0.6);
    v2Obs_PF_h[cI] = new TH1F(("v2Obs_" + centStr + "_PF_h").c_str(), ";v_{2};Counts", 150, 0.0, 0.6);
    v2ObsCorr_PF_h[cI] = new TH1F(("v2ObsCorr_" + centStr + "_PF_h").c_str(), ";v_{2};Counts", 150, 0.0, 0.6);

    v2Raw_Trk_h[cI] = new TH1F(("v2Raw_" + centStr + "_Trk_h").c_str(), ";v_{2};Counts", 150, 0.0, 0.6);
    v2RawCorr_Trk_h[cI] = new TH1F(("v2RawCorr_" + centStr + "_Trk_h").c_str(), ";v_{2};Counts", 150, 0.0, 0.6);
    v2Obs_Trk_h[cI] = new TH1F(("v2Obs_" + centStr + "_Trk_h").c_str(), ";v_{2};Counts", 150, 0.0, 0.6);
    v2ObsCorr_Trk_h[cI] = new TH1F(("v2ObsCorr_" + centStr + "_Trk_h").c_str(), ";v_{2};Counts", 150, 0.0, 0.6);

    v2Raw_EByE_h[cI] = new TH1F(("v2Raw_" + centStr + "_EByE_h").c_str(), ";v_{2};Counts", 150, 0.0, 0.6);
    v2RawCorr_EByE_h[cI] = new TH1F(("v2RawCorr_" + centStr + "_EByE_h").c_str(), ";v_{2};Counts", 150, 0.0, 0.6);
    v2Obs_EByE_h[cI] = new TH1F(("v2Obs_" + centStr + "_EByE_h").c_str(), ";v_{2};Counts", 150, 0.0, 0.6);
    v2ObsCorr_EByE_h[cI] = new TH1F(("v2ObsCorr_" + centStr + "_EByE_h").c_str(), ";v_{2};Counts", 150, 0.0, 0.6);
    v2Fit_EByE_h[cI] = new TH1F(("v2Fit_" + centStr + "_EByE_h").c_str(), ";v_{2};Counts", 150, 0.0, 0.6);
    v2FitCorr_EByE_h[cI] = new TH1F(("v2FitCorr_" + centStr + "_EByE_h").c_str(), ";v_{2};Counts", 150, 0.0, 0.6);
    v2FitV4_EByE_h[cI] = new TH1F(("v2FitV4_" + centStr + "_EByE_h").c_str(), ";v_{2};Counts", 150, 0.0, 0.6);
    v2FitV4Corr_EByE_h[cI] = new TH1F(("v2FitV4Corr_" + centStr + "_EByE_h").c_str(), ";v_{2};Counts", 150, 0.0, 0.6);
  }

  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");
  TTree* inTree_p = (TTree*)inFile_p->Get("v2V3Tree");
  TObjArray* listOfBranches = (TObjArray*)inTree_p->GetListOfBranches();
  std::vector<std::string> vectorOfBranches;
  for(Int_t oI = 0; oI < listOfBranches->GetEntries(); ++oI){
    std::string branch = listOfBranches->At(oI)->GetName();
    vectorOfBranches.push_back(branch);
  }

  bool hasTrk = false;
  bool hasPF = false;
  bool hasEByE = false;

  for(unsigned int tI = 0; tI < vectorOfBranches.size(); ++tI){
    if(vectorOfBranches.at(tI).find("trk") != std::string::npos){
      hasTrk = true;
      break;
    }

    if(vectorOfBranches.at(tI).find("pf") != std::string::npos){
      hasPF = true;
      break;
    }

    if(vectorOfBranches.at(tI).find("eByE") != std::string::npos){
      hasEByE = true;
      break;
    }
  }

  Int_t hiBin_;
  Float_t hiEvt2Plane_;
  Float_t hiEvt3Plane_;
  Float_t v2FromTree_;
  std::vector<float>* pfPhi_p=NULL;
  std::vector<float>* pfWeight_p=NULL;
  std::vector<float>* trkPt_p=NULL;
  std::vector<float>* trkPhi_p=NULL;
  std::vector<float>* trkWeight_p=NULL;  
  std::vector<float>* eByEPt_p=NULL;
  std::vector<float>* eByEPhi_p=NULL;
  std::vector<float>* eByEWeight_p=NULL;  
  
  inTree_p->SetBranchAddress("hiBin", &hiBin_);
  inTree_p->SetBranchAddress("hiEvt2Plane", &hiEvt2Plane_);
  inTree_p->SetBranchAddress("hiEvt3Plane", &hiEvt3Plane_);
  inTree_p->SetBranchAddress("v2FromTree", &v2FromTree_);
  if(hasPF){
    inTree_p->SetBranchAddress("pfPhi", &pfPhi_p);
    inTree_p->SetBranchAddress("pfWeight", &pfWeight_p);
  }
  if(hasTrk){
    inTree_p->SetBranchAddress("trkPt", &trkPt_p);
    inTree_p->SetBranchAddress("trkPhi", &trkPhi_p);
    inTree_p->SetBranchAddress("trkWeight", &trkWeight_p);
  }
  if(hasEByE){
    inTree_p->SetBranchAddress("eByEPt", &eByEPt_p);
    inTree_p->SetBranchAddress("eByEPhi", &eByEPhi_p);
    inTree_p->SetBranchAddress("eByEWeight", &eByEWeight_p);
  }

  const Int_t nEntries = TMath::Min((Int_t)inTree_p->GetEntries(), (Int_t)200000);

  double globalN[nCentBins];
  double globalV2XRawPF[nCentBins];
  double globalV2YRawPF[nCentBins];
  double globalV2XRawPFCorr[nCentBins];
  double globalV2YRawPFCorr[nCentBins];

  double globalV2XRawTrk[nCentBins];
  double globalV2YRawTrk[nCentBins];
  double globalV2XRawTrkCorr[nCentBins];
  double globalV2YRawTrkCorr[nCentBins];

  double globalV2XRawEByE[nCentBins];
  double globalV2YRawEByE[nCentBins];
  double globalV2XRawEByECorr[nCentBins];
  double globalV2YRawEByECorr[nCentBins];

  for(Int_t cI = 0; cI < nCentBins; ++cI){
    globalN[cI] = 0.;
    globalV2XRawPF[cI] = 0.;
    globalV2YRawPF[cI] = 0.;
    globalV2XRawPFCorr[cI] = 0.;
    globalV2YRawPFCorr[cI] = 0.;

    globalV2XRawTrk[cI] = 0.;
    globalV2YRawTrk[cI] = 0.;
    globalV2XRawTrkCorr[cI] = 0.;
    globalV2YRawTrkCorr[cI] = 0.;

    globalV2XRawEByE[cI] = 0.;
    globalV2YRawEByE[cI] = 0.;
    globalV2XRawEByECorr[cI] = 0.;
    globalV2YRawEByECorr[cI] = 0.;
  }


  ROOT::Math::MinimizerOptions::SetDefaultStrategy(0);
  ROOT::Math::MinimizerOptions::SetDefaultTolerance(10);

  std::string flowFitForm = "[0]*(1. + 2.*([1]*TMath::Cos(2.*(x - [2])) + [3]*TMath::Cos(3.*(x - [4]))))";
  std::string flowFitFormV4 = "[0]*(1. + 2.*([1]*TMath::Cos(2.*(x - [2])) + [3]*TMath::Cos(3.*(x - [4])) + [5]*TMath::Cos(4.*(x - [6])) ))";

  TF1* flowFit_p = new TF1("flowFit_p", flowFitForm.c_str(), -TMath::Pi(), TMath::Pi());
  TF1* flowFitCorr_p = new TF1("flowFitCorr_p", flowFitForm.c_str(), -TMath::Pi(), TMath::Pi());
  TF1* flowFitV4_p = new TF1("flowFitV4_p", flowFitFormV4.c_str(), -TMath::Pi(), TMath::Pi());
  TF1* flowFitV4Corr_p = new TF1("flowFitV4Corr_p", flowFitFormV4.c_str(), -TMath::Pi(), TMath::Pi());

  
  cppWatch allLoop;
  cppWatch fitLoop;
  cppWatch fitLoopA;
  cppWatch fitLoopB;
  cppWatch fitLoopC;

  allLoop.start();
  for(Int_t entry = 0; entry < nEntries; ++entry){
    if(entry%10000 == 0) std::cout << " Entry " << entry << "/" << nEntries << std::endl;
    inTree_p->GetEntry(entry);

    int centPos = -1;
    for(Int_t cI = 0; cI < nCentBins; ++cI){
      if(centBinsLow[cI] <= hiBin_/2 && centBinsHi[cI] > hiBin_/2){
	centPos = cI;
	break;
      }
    }

    if(centPos < 0) continue;

    double v2xRawPF = 0.;
    double v2yRawPF = 0.;

    double v2xRawPFCorr = 0.;
    double v2yRawPFCorr = 0.;

    double v2xRawTrk = 0.;
    double v2yRawTrk = 0.;

    double v2xRawTrkCorr = 0.;
    double v2yRawTrkCorr = 0.;

    double v2xRawEByE = 0.;
    double v2yRawEByE = 0.;

    double v2xRawEByECorr = 0.;
    double v2yRawEByECorr = 0.;

    double weightPF = 0.;
    double weightTrk = 0.;
    double weightEByE = 0.;
  
    if(hasPF){
      for(unsigned pfI = 0; pfI < pfPhi_p->size(); pfI++){
	double tempWeightPF = pfWeight_p->at(pfI);
	double deltaEventPhi = pfPhi_p->at(pfI) - hiEvt2Plane_;
	
	v2xRawPF += TMath::Cos(2*(deltaEventPhi));
	v2yRawPF += TMath::Sin(2*(deltaEventPhi));
	
	v2xRawPFCorr += tempWeightPF*TMath::Cos(2*(deltaEventPhi));
	v2yRawPFCorr += tempWeightPF*TMath::Sin(2*(deltaEventPhi));
	
	weightPF += tempWeightPF;
      }
    }

    Int_t trkCounter = 0;
    if(hasTrk){    
      for(unsigned tI = 0; tI < trkPhi_p->size(); tI++){
	//	if(trkPt_p->at(tI) >= 2.4) continue;
	trkCounter++;
	double tempWeightTrk = trkWeight_p->at(tI);
	double deltaEventPhi = trkPhi_p->at(tI) - hiEvt2Plane_;
	
	v2xRawTrk += TMath::Cos(2*(deltaEventPhi));
	v2yRawTrk += TMath::Sin(2*(deltaEventPhi));

	v2xRawTrkCorr += tempWeightTrk*TMath::Cos(2*(deltaEventPhi));
	v2yRawTrkCorr += tempWeightTrk*TMath::Sin(2*(deltaEventPhi));
	
	weightTrk += tempWeightTrk;
      }
    }

    Int_t eByECounter = 0;
    if(hasEByE){    
      
      const Int_t nEByE = eByEPhi_p->size();
      const int nPhiBins = std::fmax(10, nEByE/30);

      double eventPlane2Cos = 0;
      double eventPlane2Sin = 0;

      double eventPlane3Cos = 0;
      double eventPlane3Sin = 0;

      double eventPlane4Cos = 0;
      double eventPlane4Sin = 0;

      double eventPlane2CosCorr = 0;
      double eventPlane2SinCorr = 0;

      double eventPlane3CosCorr = 0;
      double eventPlane3SinCorr = 0;

      double eventPlane4CosCorr = 0;
      double eventPlane4SinCorr = 0;

      TH1F* phi_h = new TH1F("phi_h", ";#phi;Counts (.3 < p_{T} < 3.)", nPhiBins, -TMath::Pi(), TMath::Pi());
      TH1F* phiCorr_h = new TH1F("phiCorr_h", ";#phi;Counts Corr. (.3 < p_{T} < 3.)", nPhiBins, -TMath::Pi(), TMath::Pi());
      

      for(unsigned tI = 0; tI < eByEPhi_p->size(); tI++){
	//	if(eByEPt_p->at(tI) >= 2.4) continue;
	eByECounter++;
	double tempWeightEByE = eByEWeight_p->at(tI);
	//	double deltaEventPhi = eByEPhi_p->at(tI) - hiEvt2Plane_;
	double deltaEventPhi = eByEPhi_p->at(tI);
	
	v2xRawEByE += TMath::Cos(2*(deltaEventPhi));
	v2yRawEByE += TMath::Sin(2*(deltaEventPhi));

	v2xRawEByECorr += tempWeightEByE*TMath::Cos(2*(deltaEventPhi));
	v2yRawEByECorr += tempWeightEByE*TMath::Sin(2*(deltaEventPhi));
	
	weightEByE += tempWeightEByE;
	
	phi_h->Fill(deltaEventPhi);
	phiCorr_h->Fill(deltaEventPhi, tempWeightEByE);

	eventPlane2Cos += std::cos(2*deltaEventPhi);
	eventPlane2Sin += std::sin(2*deltaEventPhi);

	eventPlane3Cos += std::cos(3*deltaEventPhi);
	eventPlane3Sin += std::sin(3*deltaEventPhi);

	eventPlane4Cos += std::cos(4*deltaEventPhi);
	eventPlane4Sin += std::sin(4*deltaEventPhi);

	eventPlane2CosCorr += tempWeightEByE*std::cos(2*deltaEventPhi);
	eventPlane2SinCorr += tempWeightEByE*std::sin(2*deltaEventPhi);

	eventPlane3CosCorr += tempWeightEByE*std::cos(3*deltaEventPhi);
	eventPlane3SinCorr += tempWeightEByE*std::sin(3*deltaEventPhi);       

	eventPlane4CosCorr += tempWeightEByE*std::cos(4*deltaEventPhi);
	eventPlane4SinCorr += tempWeightEByE*std::sin(4*deltaEventPhi);       
      }

      
      double eventPlane2 = std::atan2(eventPlane2Sin, eventPlane2Cos)/2.;
      double eventPlane3 = std::atan2(eventPlane3Sin, eventPlane3Cos)/3.;
      double eventPlane4 = std::atan2(eventPlane4Sin, eventPlane4Cos)/4.;

      double eventPlane2Corr = std::atan2(eventPlane2SinCorr, eventPlane2CosCorr)/2.;
      double eventPlane3Corr = std::atan2(eventPlane3SinCorr, eventPlane3CosCorr)/3.;
      double eventPlane4Corr = std::atan2(eventPlane4SinCorr, eventPlane4CosCorr)/4.;


      flowFit_p->SetParameter(0, phi_h->GetMaximum()*3./4.);
      flowFit_p->SetParameter(1, 0);
      flowFit_p->SetParameter(2, eventPlane2);
      flowFit_p->SetParLimits(2, eventPlane2-0.0001, eventPlane2 + 0.0001);

      flowFit_p->SetParameter(3, 0);
      flowFit_p->SetParameter(4, eventPlane3);
      flowFit_p->SetParLimits(4, eventPlane3 - 0.0001, eventPlane3 + 0.0001);

      fitLoopA.stop();
      fitLoopB.start();

      phi_h->Fit("flowFit_p", "Q N 0", "", -TMath::Pi(), TMath::Pi());

      fitLoopB.stop();
      fitLoopC.start();

      v2Fit_EByE_h[centPos]->Fill(flowFit_p->GetParameter(1));


      flowFitCorr_p->SetParameter(0, phiCorr_h->GetMaximum()*3./4.);
      flowFitCorr_p->SetParameter(1, 0);
      flowFit_p->SetParameter(2, eventPlane2Corr);
      flowFit_p->SetParLimits(2, eventPlane2Corr-0.0001, eventPlane2Corr + 0.0001);

      flowFit_p->SetParameter(3, 0);
      flowFit_p->SetParameter(4, eventPlane3Corr);
      flowFit_p->SetParLimits(4, eventPlane3Corr - 0.0001, eventPlane3Corr + 0.0001);

      fitLoopA.stop();
      fitLoopB.start();

      phiCorr_h->Fit("flowFitCorr_p", "Q N 0", "", -TMath::Pi(), TMath::Pi());

      v2FitCorr_EByE_h[centPos]->Fill(flowFitCorr_p->GetParameter(1));


      flowFitV4_p->SetParameter(0, phi_h->GetMaximum()*3./4.);
      flowFitV4_p->SetParameter(1, 0);
      flowFitV4_p->SetParameter(2, eventPlane2);
      flowFitV4_p->SetParLimits(2, eventPlane2-0.0001, eventPlane2 + 0.0001);

      flowFitV4_p->SetParameter(3, 0);
      flowFitV4_p->SetParameter(4, eventPlane3);
      flowFitV4_p->SetParLimits(4, eventPlane3 - 0.0001, eventPlane3 + 0.0001);

      flowFitV4_p->SetParameter(5, 0);
      flowFitV4_p->SetParameter(6, eventPlane4);
      flowFitV4_p->SetParLimits(6, eventPlane4 - 0.0001, eventPlane4 + 0.0001);

      fitLoopA.stop();
      fitLoopB.start();

      phi_h->Fit("flowFitV4_p", "Q N 0", "", -TMath::Pi(), TMath::Pi());

      fitLoopB.stop();
      fitLoopC.start();

      v2FitV4_EByE_h[centPos]->Fill(flowFitV4_p->GetParameter(1));


      flowFitV4Corr_p->SetParameter(0, phiCorr_h->GetMaximum()*3./4.);
      flowFitV4Corr_p->SetParameter(1, 0);
      flowFitV4Corr_p->SetParameter(2, eventPlane2Corr);
      flowFitV4Corr_p->SetParLimits(2, eventPlane2Corr-0.0001, eventPlane2Corr + 0.0001);

      flowFitV4Corr_p->SetParameter(3, 0);
      flowFitV4Corr_p->SetParameter(4, eventPlane3Corr);
      flowFitV4Corr_p->SetParLimits(4, eventPlane3Corr - 0.0001, eventPlane3Corr + 0.0001);

      flowFitV4Corr_p->SetParameter(5, 0);
      flowFitV4Corr_p->SetParameter(6, eventPlane4Corr);
      flowFitV4Corr_p->SetParLimits(6, eventPlane4Corr - 0.0001, eventPlane4Corr + 0.0001);

      fitLoopA.stop();
      fitLoopB.start();

      phiCorr_h->Fit("flowFitV4Corr_p", "Q N 0", "", -TMath::Pi(), TMath::Pi());

      v2FitV4Corr_EByE_h[centPos]->Fill(flowFitV4Corr_p->GetParameter(1));
 
      fitLoopC.stop();
      fitLoop.stop();
      delete phi_h;      
      delete phiCorr_h;      

    }
  
    if(hasPF){
      v2xRawPF /= (double)pfPhi_p->size();
      v2yRawPF /= (double)pfPhi_p->size();
      
      v2xRawPFCorr /= weightPF;
      v2yRawPFCorr /= weightPF;
      
      double v2RawPF = TMath::Sqrt(v2xRawPF*v2xRawPF + v2yRawPF*v2yRawPF);
      double v2RawPFCorr = TMath::Sqrt(v2xRawPFCorr*v2xRawPFCorr + v2yRawPFCorr*v2yRawPFCorr);
      
      v2Raw_PF_h[centPos]->Fill(v2RawPF);
      v2RawCorr_PF_h[centPos]->Fill(v2RawPFCorr);
      
      globalN[centPos] += 1;
      globalV2XRawPF[centPos] += v2xRawPF;
      globalV2YRawPF[centPos] += v2yRawPF;
      globalV2XRawPFCorr[centPos] += v2xRawPFCorr;
      globalV2YRawPFCorr[centPos] += v2yRawPFCorr;
    }

    if(hasTrk){
      v2xRawTrk /= (double)trkCounter;
      v2yRawTrk /= (double)trkCounter;
      
      v2xRawTrkCorr /= weightTrk;
      v2yRawTrkCorr /= weightTrk;
      
      double v2RawTrk = TMath::Sqrt(v2xRawTrk*v2xRawTrk + v2yRawTrk*v2yRawTrk);
      double v2RawTrkCorr = TMath::Sqrt(v2xRawTrkCorr*v2xRawTrkCorr + v2yRawTrkCorr*v2yRawTrkCorr);
      
      v2Raw_Trk_h[centPos]->Fill(v2RawTrk);
      v2RawCorr_Trk_h[centPos]->Fill(v2RawTrkCorr);
      
      globalN[centPos] += 1;
      globalV2XRawTrk[centPos] += v2xRawTrk;
      globalV2YRawTrk[centPos] += v2yRawTrk;
      globalV2XRawTrkCorr[centPos] += v2xRawTrkCorr;
      globalV2YRawTrkCorr[centPos] += v2yRawTrkCorr;
    }

    if(hasEByE){
      v2xRawEByE /= (double)eByECounter;
      v2yRawEByE /= (double)eByECounter;
      
      v2xRawEByECorr /= weightEByE;
      v2yRawEByECorr /= weightEByE;
      

      //      std::cout << "hiBin,X,Y,XCorr,YCorr: " << hiBin_ << ", " << v2xRawEByE << ", " << v2yRawEByE << ", " << v2xRawEByECorr << ", " << v2yRawEByECorr << std::endl;

      double v2RawEByE = TMath::Sqrt(v2xRawEByE*v2xRawEByE + v2yRawEByE*v2yRawEByE);
      double v2RawEByECorr = TMath::Sqrt(v2xRawEByECorr*v2xRawEByECorr + v2yRawEByECorr*v2yRawEByECorr);
      
      v2Raw_EByE_h[centPos]->Fill(v2RawEByE);
      v2RawCorr_EByE_h[centPos]->Fill(v2RawEByECorr);
      
      globalN[centPos] += 1;
      globalV2XRawEByE[centPos] += v2xRawEByE;
      globalV2YRawEByE[centPos] += v2yRawEByE;
      globalV2XRawEByECorr[centPos] += v2xRawEByECorr;
      globalV2YRawEByECorr[centPos] += v2yRawEByECorr;
    }
  }
  allLoop.stop();

  delete flowFit_p;
  delete flowFitCorr_p;
  delete flowFitV4_p;
  delete flowFitV4Corr_p;
  

  std::cout << "All loop: " << allLoop.total() << std::endl;
  std::cout << " Fit loop: " << fitLoop.total() << ", " << fitLoop.total()/allLoop.total() << std::endl;
  std::cout << "  Fit loop A: " << fitLoopA.total() << ", " << fitLoopA.total()/allLoop.total() << std::endl;
  std::cout << "  Fit loop B: " << fitLoopB.total() << ", " << fitLoopB.total()/allLoop.total() << std::endl;
  std::cout << "  Fit loop C: " << fitLoopC.total() << ", " << fitLoopC.total()/allLoop.total() << std::endl;
  

  for(Int_t cI = 0; cI < nCentBins; ++cI){
    globalV2XRawPF[cI] /= globalN[cI];
    globalV2YRawPF[cI] /= globalN[cI];
    globalV2XRawPFCorr[cI] /= globalN[cI];
    globalV2YRawPFCorr[cI] /= globalN[cI];

    if(hasTrk){
      globalV2XRawTrk[cI] /= globalN[cI];
      globalV2YRawTrk[cI] /= globalN[cI];
      globalV2XRawTrkCorr[cI] /= globalN[cI];
      globalV2YRawTrkCorr[cI] /= globalN[cI];
    }

    if(hasEByE){
      globalV2XRawEByE[cI] /= globalN[cI];
      globalV2YRawEByE[cI] /= globalN[cI];
      globalV2XRawEByECorr[cI] /= globalN[cI];
      globalV2YRawEByECorr[cI] /= globalN[cI];
    }

    std::cout << "Cent, N: " << centBinsLow[cI] << "-" << centBinsHi[cI] << "%, " << globalN[cI] << std::endl;

    if(hasPF){
      std::cout << "  Global v2xRawPF, v2yRawPF, N: " << globalV2XRawPF[cI] << ", " << globalV2YRawPF[cI] << std::endl;
      std::cout << "  Global v2xRawPFCorr, v2yRawPFCorr, N: " << globalV2XRawPFCorr[cI] << ", " << globalV2YRawPFCorr[cI] << std::endl;
    }

    if(hasTrk){  
      std::cout << "  Global v2xRawTrk, v2yRawTrk, N: " << globalV2XRawTrk[cI] << ", " << globalV2YRawTrk[cI] << std::endl;
      std::cout << "  Global v2xRawTrkCorr, v2yRawTrkCorr, N: " << globalV2XRawTrkCorr[cI] << ", " << globalV2YRawTrkCorr[cI] << std::endl;
    }

    if(hasEByE){  
      std::cout << "  Global v2xRawEByE, v2yRawEByE, N: " << globalV2XRawEByE[cI] << ", " << globalV2YRawEByE[cI] << std::endl;
      std::cout << "  Global v2xRawEByECorr, v2yRawEByECorr, N: " << globalV2XRawEByECorr[cI] << ", " << globalV2YRawEByECorr[cI] << std::endl;
    }
  }



  for(Int_t entry = 0; entry < nEntries; ++entry){
    if(entry%10000 == 0) std::cout << " Entry " << entry << "/" << nEntries << std::endl;
    inTree_p->GetEntry(entry);

    int centPos = -1;
    for(Int_t cI = 0; cI < nCentBins; ++cI){
      if(centBinsLow[cI] <= hiBin_/2 && centBinsHi[cI] > hiBin_/2){
	centPos = cI;
	break;
      }
    }

    if(centPos < 0) continue;

    double v2xObsPF = 0.;
    double v2yObsPF = 0.;

    double v2xObsPFCorr = 0.;
    double v2yObsPFCorr = 0.;

    double weightPF = 0.;

    if(hasPF){
      for(unsigned pfI = 0; pfI < pfPhi_p->size(); pfI++){
	double tempWeightPF = pfWeight_p->at(pfI);
	double deltaEventPhi = pfPhi_p->at(pfI) - hiEvt2Plane_;
	
	v2xObsPF += TMath::Cos(2*(deltaEventPhi));
	v2yObsPF += TMath::Sin(2*(deltaEventPhi));
	
	v2xObsPFCorr += tempWeightPF*TMath::Cos(2*(deltaEventPhi));
	v2yObsPFCorr += tempWeightPF*TMath::Sin(2*(deltaEventPhi));
	
	weightPF += tempWeightPF;
      }
    }


    double v2xObsTrk = 0.;
    double v2yObsTrk = 0.;

    double v2xObsTrkCorr = 0.;
    double v2yObsTrkCorr = 0.;

    double weightTrk = 0.;

    Int_t trkCounter = 0;
    if(hasTrk){
      for(unsigned tI = 0; tI < trkPhi_p->size(); tI++){
	//	if(trkPt_p->at(tI) >= 2.4) continue;
	trkCounter++;
	double tempWeightTrk = trkWeight_p->at(tI);
	double deltaEventPhi = trkPhi_p->at(tI) - hiEvt2Plane_;
	
	v2xObsTrk += TMath::Cos(2*(deltaEventPhi));
	v2yObsTrk += TMath::Sin(2*(deltaEventPhi));
	
	v2xObsTrkCorr += tempWeightTrk*TMath::Cos(2*(deltaEventPhi));
	v2yObsTrkCorr += tempWeightTrk*TMath::Sin(2*(deltaEventPhi));
	
	weightTrk += tempWeightTrk;
      }
    }

    double v2xObsEByE = 0.;
    double v2yObsEByE = 0.;

    double v2xObsEByECorr = 0.;
    double v2yObsEByECorr = 0.;

    double weightEByE = 0.;

    Int_t eByECounter = 0;
    if(hasEByE){
      for(unsigned tI = 0; tI < eByEPhi_p->size(); tI++){
	//	if(eByEPt_p->at(tI) >= 2.4) continue;
	eByECounter++;
	double tempWeightEByE = eByEWeight_p->at(tI);
	//	double deltaEventPhi = eByEPhi_p->at(tI) - hiEvt2Plane_;
	double deltaEventPhi = eByEPhi_p->at(tI);
	
	v2xObsEByE += TMath::Cos(2*(deltaEventPhi));
	v2yObsEByE += TMath::Sin(2*(deltaEventPhi));
	
	v2xObsEByECorr += tempWeightEByE*TMath::Cos(2*(deltaEventPhi));
	v2yObsEByECorr += tempWeightEByE*TMath::Sin(2*(deltaEventPhi));
	
	weightEByE += tempWeightEByE;
      }
    }

    if(hasPF){
      v2xObsPF /= (double)pfPhi_p->size();
      v2yObsPF /= (double)pfPhi_p->size();
      
      v2xObsPFCorr /= weightPF;
      v2yObsPFCorr /= weightPF;
      
      v2xObsPF -= globalV2XRawPF[centPos];
      v2yObsPF -= globalV2YRawPF[centPos];
      
      v2xObsPFCorr -= globalV2XRawPFCorr[centPos];
      v2yObsPFCorr -= globalV2YRawPFCorr[centPos];
      
      double v2ObsPF = TMath::Sqrt(v2xObsPF*v2xObsPF + v2yObsPF*v2yObsPF);
      double v2ObsPFCorr = TMath::Sqrt(v2xObsPFCorr*v2xObsPFCorr + v2yObsPFCorr*v2yObsPFCorr);
      
      v2Obs_PF_h[centPos]->Fill(v2ObsPF);
      v2ObsCorr_PF_h[centPos]->Fill(v2ObsPFCorr);
    }

    if(hasTrk){
      v2xObsTrk /= (double)trkCounter;
      v2yObsTrk /= (double)trkCounter;
      
      v2xObsTrkCorr /= weightTrk;
      v2yObsTrkCorr /= weightTrk;
      
      v2xObsTrk -= globalV2XRawTrk[centPos];
      v2yObsTrk -= globalV2YRawTrk[centPos];
      
      v2xObsTrkCorr -= globalV2XRawTrkCorr[centPos];
      v2yObsTrkCorr -= globalV2YRawTrkCorr[centPos];
      
      double v2ObsTrk = TMath::Sqrt(v2xObsTrk*v2xObsTrk + v2yObsTrk*v2yObsTrk);
      double v2ObsTrkCorr = TMath::Sqrt(v2xObsTrkCorr*v2xObsTrkCorr + v2yObsTrkCorr*v2yObsTrkCorr);
      
      v2Obs_Trk_h[centPos]->Fill(v2ObsTrk);
      v2ObsCorr_Trk_h[centPos]->Fill(v2ObsTrkCorr);
    }

    if(hasEByE){
      v2xObsEByE /= (double)eByECounter;
      v2yObsEByE /= (double)eByECounter;
      
      v2xObsEByECorr /= weightEByE;
      v2yObsEByECorr /= weightEByE;
      
      v2xObsEByE -= globalV2XRawEByE[centPos];
      v2yObsEByE -= globalV2YRawEByE[centPos];
      
      v2xObsEByECorr -= globalV2XRawEByECorr[centPos];
      v2yObsEByECorr -= globalV2YRawEByECorr[centPos];
      
      double v2ObsEByE = TMath::Sqrt(v2xObsEByE*v2xObsEByE + v2yObsEByE*v2yObsEByE);
      double v2ObsEByECorr = TMath::Sqrt(v2xObsEByECorr*v2xObsEByECorr + v2yObsEByECorr*v2yObsEByECorr);
      
      v2Obs_EByE_h[centPos]->Fill(v2ObsEByE);
      v2ObsCorr_EByE_h[centPos]->Fill(v2ObsEByECorr);
    }
  }


  inFile_p->Close();
  delete inFile_p;

  outFile_p->cd();

  if(hasPF){
    for(Int_t cI = 0; cI < nCentBins; ++cI){  
      v2Raw_PF_h[cI]->Write("", TObject::kOverwrite);
      v2RawCorr_PF_h[cI]->Write("", TObject::kOverwrite);
      
      v2Obs_PF_h[cI]->Write("", TObject::kOverwrite);
      v2ObsCorr_PF_h[cI]->Write("", TObject::kOverwrite);
      
      v2Raw_Mean_PF_h->SetBinContent(cI+1, v2Raw_PF_h[cI]->GetMean());
      v2Raw_Mean_PF_h->SetBinError(cI+1, v2Raw_PF_h[cI]->GetMeanError());
      v2Raw_Sigma_PF_h->SetBinContent(cI+1, v2Raw_PF_h[cI]->GetStdDev());
      v2Raw_Sigma_PF_h->SetBinError(cI+1, v2Raw_PF_h[cI]->GetStdDevError());
      
      v2RawCorr_Mean_PF_h->SetBinContent(cI+1, v2RawCorr_PF_h[cI]->GetMean());
      v2RawCorr_Mean_PF_h->SetBinError(cI+1, v2RawCorr_PF_h[cI]->GetMeanError());
      v2RawCorr_Sigma_PF_h->SetBinContent(cI+1, v2RawCorr_PF_h[cI]->GetStdDev());
      v2RawCorr_Sigma_PF_h->SetBinError(cI+1, v2RawCorr_PF_h[cI]->GetStdDevError());
      
      v2Obs_Mean_PF_h->SetBinContent(cI+1, v2Obs_PF_h[cI]->GetMean());
      v2Obs_Mean_PF_h->SetBinError(cI+1, v2Obs_PF_h[cI]->GetMeanError());
      v2Obs_Sigma_PF_h->SetBinContent(cI+1, v2Obs_PF_h[cI]->GetStdDev());
      v2Obs_Sigma_PF_h->SetBinError(cI+1, v2Obs_PF_h[cI]->GetStdDevError());
      
      v2ObsCorr_Mean_PF_h->SetBinContent(cI+1, v2ObsCorr_PF_h[cI]->GetMean());
      v2ObsCorr_Mean_PF_h->SetBinError(cI+1, v2ObsCorr_PF_h[cI]->GetMeanError());
      v2ObsCorr_Sigma_PF_h->SetBinContent(cI+1, v2ObsCorr_PF_h[cI]->GetStdDev());
      v2ObsCorr_Sigma_PF_h->SetBinError(cI+1, v2ObsCorr_PF_h[cI]->GetStdDevError());
      
      delete v2Raw_PF_h[cI];
      delete v2RawCorr_PF_h[cI];
      
      delete v2Obs_PF_h[cI];
      delete v2ObsCorr_PF_h[cI];
    }

    v2Raw_Mean_PF_h->Write("", TObject::kOverwrite);
    v2Raw_Sigma_PF_h->Write("", TObject::kOverwrite);
    
    v2RawCorr_Mean_PF_h->Write("", TObject::kOverwrite);
    v2RawCorr_Sigma_PF_h->Write("", TObject::kOverwrite);
    
    v2Obs_Mean_PF_h->Write("", TObject::kOverwrite);
    v2Obs_Sigma_PF_h->Write("", TObject::kOverwrite);
    
    v2ObsCorr_Mean_PF_h->Write("", TObject::kOverwrite);
    v2ObsCorr_Sigma_PF_h->Write("", TObject::kOverwrite);
    
    delete v2Raw_Mean_PF_h;
    delete v2Raw_Sigma_PF_h;
    delete v2RawCorr_Mean_PF_h;
    delete v2RawCorr_Sigma_PF_h;
    
    delete v2Obs_Mean_PF_h;
    delete v2Obs_Sigma_PF_h;
    delete v2ObsCorr_Mean_PF_h;
    delete v2ObsCorr_Sigma_PF_h;
  }

  for(Int_t cI = 0; cI < nCentBins; ++cI){  
    if(hasTrk){
      v2Raw_Trk_h[cI]->Write("", TObject::kOverwrite);
      v2RawCorr_Trk_h[cI]->Write("", TObject::kOverwrite);
      
      v2Obs_Trk_h[cI]->Write("", TObject::kOverwrite);
      v2ObsCorr_Trk_h[cI]->Write("", TObject::kOverwrite);

      v2Raw_Mean_Trk_h->SetBinContent(cI+1, v2Raw_Trk_h[cI]->GetMean());
      v2Raw_Mean_Trk_h->SetBinError(cI+1, v2Raw_Trk_h[cI]->GetMeanError());
      v2Raw_Sigma_Trk_h->SetBinContent(cI+1, v2Raw_Trk_h[cI]->GetStdDev());
      v2Raw_Sigma_Trk_h->SetBinError(cI+1, v2Raw_Trk_h[cI]->GetStdDevError());
      
      v2RawCorr_Mean_Trk_h->SetBinContent(cI+1, v2RawCorr_Trk_h[cI]->GetMean());
      v2RawCorr_Mean_Trk_h->SetBinError(cI+1, v2RawCorr_Trk_h[cI]->GetMeanError());
      v2RawCorr_Sigma_Trk_h->SetBinContent(cI+1, v2RawCorr_Trk_h[cI]->GetStdDev());
      v2RawCorr_Sigma_Trk_h->SetBinError(cI+1, v2RawCorr_Trk_h[cI]->GetStdDevError());
      
      v2Obs_Mean_Trk_h->SetBinContent(cI+1, v2Obs_Trk_h[cI]->GetMean());
      v2Obs_Mean_Trk_h->SetBinError(cI+1, v2Obs_Trk_h[cI]->GetMeanError());
      v2Obs_Sigma_Trk_h->SetBinContent(cI+1, v2Obs_Trk_h[cI]->GetStdDev());
      v2Obs_Sigma_Trk_h->SetBinError(cI+1, v2Obs_Trk_h[cI]->GetStdDevError());
      
      v2ObsCorr_Mean_Trk_h->SetBinContent(cI+1, v2ObsCorr_Trk_h[cI]->GetMean());
      v2ObsCorr_Mean_Trk_h->SetBinError(cI+1, v2ObsCorr_Trk_h[cI]->GetMeanError());
      v2ObsCorr_Sigma_Trk_h->SetBinContent(cI+1, v2ObsCorr_Trk_h[cI]->GetStdDev());
      v2ObsCorr_Sigma_Trk_h->SetBinError(cI+1, v2ObsCorr_Trk_h[cI]->GetStdDevError());
    }

    delete v2Raw_Trk_h[cI];
    delete v2RawCorr_Trk_h[cI];

    delete v2Obs_Trk_h[cI];
    delete v2ObsCorr_Trk_h[cI];
  }

  v2Raw_Mean_Trk_h->Write("", TObject::kOverwrite);
  v2Raw_Sigma_Trk_h->Write("", TObject::kOverwrite);

  v2RawCorr_Mean_Trk_h->Write("", TObject::kOverwrite);
  v2RawCorr_Sigma_Trk_h->Write("", TObject::kOverwrite);

  v2Obs_Mean_Trk_h->Write("", TObject::kOverwrite);
  v2Obs_Sigma_Trk_h->Write("", TObject::kOverwrite);

  v2ObsCorr_Mean_Trk_h->Write("", TObject::kOverwrite);
  v2ObsCorr_Sigma_Trk_h->Write("", TObject::kOverwrite);

  if(hasTrk){
    delete v2Raw_Mean_Trk_h;
    delete v2Raw_Sigma_Trk_h;
    delete v2RawCorr_Mean_Trk_h;
    delete v2RawCorr_Sigma_Trk_h;
    
    delete v2Obs_Mean_Trk_h;
    delete v2Obs_Sigma_Trk_h;
    delete v2ObsCorr_Mean_Trk_h;
    delete v2ObsCorr_Sigma_Trk_h;
  }


  for(Int_t cI = 0; cI < nCentBins; ++cI){  
    if(hasEByE){
      v2Raw_EByE_h[cI]->Write("", TObject::kOverwrite);
      v2RawCorr_EByE_h[cI]->Write("", TObject::kOverwrite);

      v2Fit_EByE_h[cI]->Write("", TObject::kOverwrite);
      v2FitCorr_EByE_h[cI]->Write("", TObject::kOverwrite);
      v2FitV4_EByE_h[cI]->Write("", TObject::kOverwrite);
      v2FitV4Corr_EByE_h[cI]->Write("", TObject::kOverwrite);
      
      v2Obs_EByE_h[cI]->Write("", TObject::kOverwrite);
      v2ObsCorr_EByE_h[cI]->Write("", TObject::kOverwrite);

      v2Raw_Mean_EByE_h->SetBinContent(cI+1, v2Raw_EByE_h[cI]->GetMean());
      v2Raw_Mean_EByE_h->SetBinError(cI+1, v2Raw_EByE_h[cI]->GetMeanError());
      v2Raw_Sigma_EByE_h->SetBinContent(cI+1, v2Raw_EByE_h[cI]->GetStdDev());
      v2Raw_Sigma_EByE_h->SetBinError(cI+1, v2Raw_EByE_h[cI]->GetStdDevError());
      
      v2RawCorr_Mean_EByE_h->SetBinContent(cI+1, v2RawCorr_EByE_h[cI]->GetMean());
      v2RawCorr_Mean_EByE_h->SetBinError(cI+1, v2RawCorr_EByE_h[cI]->GetMeanError());
      v2RawCorr_Sigma_EByE_h->SetBinContent(cI+1, v2RawCorr_EByE_h[cI]->GetStdDev());
      v2RawCorr_Sigma_EByE_h->SetBinError(cI+1, v2RawCorr_EByE_h[cI]->GetStdDevError());

      v2Fit_Mean_EByE_h->SetBinContent(cI+1, v2Fit_EByE_h[cI]->GetMean());
      v2Fit_Mean_EByE_h->SetBinError(cI+1, v2Fit_EByE_h[cI]->GetMeanError());
      v2Fit_Sigma_EByE_h->SetBinContent(cI+1, v2Fit_EByE_h[cI]->GetStdDev());
      v2Fit_Sigma_EByE_h->SetBinError(cI+1, v2Fit_EByE_h[cI]->GetStdDevError());
      
      v2FitCorr_Mean_EByE_h->SetBinContent(cI+1, v2FitCorr_EByE_h[cI]->GetMean());
      v2FitCorr_Mean_EByE_h->SetBinError(cI+1, v2FitCorr_EByE_h[cI]->GetMeanError());
      v2FitCorr_Sigma_EByE_h->SetBinContent(cI+1, v2FitCorr_EByE_h[cI]->GetStdDev());
      v2FitCorr_Sigma_EByE_h->SetBinError(cI+1, v2FitCorr_EByE_h[cI]->GetStdDevError());
      
      v2Obs_Mean_EByE_h->SetBinContent(cI+1, v2Obs_EByE_h[cI]->GetMean());
      v2Obs_Mean_EByE_h->SetBinError(cI+1, v2Obs_EByE_h[cI]->GetMeanError());
      v2Obs_Sigma_EByE_h->SetBinContent(cI+1, v2Obs_EByE_h[cI]->GetStdDev());
      v2Obs_Sigma_EByE_h->SetBinError(cI+1, v2Obs_EByE_h[cI]->GetStdDevError());
      
      v2ObsCorr_Mean_EByE_h->SetBinContent(cI+1, v2ObsCorr_EByE_h[cI]->GetMean());
      v2ObsCorr_Mean_EByE_h->SetBinError(cI+1, v2ObsCorr_EByE_h[cI]->GetMeanError());
      v2ObsCorr_Sigma_EByE_h->SetBinContent(cI+1, v2ObsCorr_EByE_h[cI]->GetStdDev());
      v2ObsCorr_Sigma_EByE_h->SetBinError(cI+1, v2ObsCorr_EByE_h[cI]->GetStdDevError());
    }

    delete v2Raw_EByE_h[cI];
    delete v2RawCorr_EByE_h[cI];

    delete v2Fit_EByE_h[cI];
    delete v2FitCorr_EByE_h[cI];

    delete v2FitV4_EByE_h[cI];
    delete v2FitV4Corr_EByE_h[cI];

    delete v2Obs_EByE_h[cI];
    delete v2ObsCorr_EByE_h[cI];
  }

  v2Raw_Mean_EByE_h->Write("", TObject::kOverwrite);
  v2Raw_Sigma_EByE_h->Write("", TObject::kOverwrite);

  v2RawCorr_Mean_EByE_h->Write("", TObject::kOverwrite);
  v2RawCorr_Sigma_EByE_h->Write("", TObject::kOverwrite);

  v2Fit_Mean_EByE_h->Write("", TObject::kOverwrite);
  v2Fit_Sigma_EByE_h->Write("", TObject::kOverwrite);

  v2FitCorr_Mean_EByE_h->Write("", TObject::kOverwrite);
  v2FitCorr_Sigma_EByE_h->Write("", TObject::kOverwrite);

  v2Obs_Mean_EByE_h->Write("", TObject::kOverwrite);
  v2Obs_Sigma_EByE_h->Write("", TObject::kOverwrite);

  v2ObsCorr_Mean_EByE_h->Write("", TObject::kOverwrite);
  v2ObsCorr_Sigma_EByE_h->Write("", TObject::kOverwrite);

  if(hasEByE){
    delete v2Raw_Mean_EByE_h;
    delete v2Raw_Sigma_EByE_h;
    delete v2Fit_Mean_EByE_h;
    delete v2Fit_Sigma_EByE_h;
    delete v2RawCorr_Mean_EByE_h;
    delete v2RawCorr_Sigma_EByE_h;
    delete v2FitCorr_Mean_EByE_h;
    delete v2FitCorr_Sigma_EByE_h;
    
    delete v2Obs_Mean_EByE_h;
    delete v2Obs_Sigma_EByE_h;
    delete v2ObsCorr_Mean_EByE_h;
    delete v2ObsCorr_Sigma_EByE_h;
  }

  outFile_p->Close();
  delete outFile_p;

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "USAGE: ./recreateV2V3TreeHist.exe <inFileName>" << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += recreateV2V3TreeHist(argv[1]);
  return retVal;
}
