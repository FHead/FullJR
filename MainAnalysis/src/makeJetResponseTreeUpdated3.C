//cpp dependencies
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

//ROOT dependencies
#include "TDirectory.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLatex.h"
#include "TMath.h"
#include "TNamed.h"
#include "TRandom3.h"
#include "TTree.h"

//RooUnfold dependencies
#include "src/RooUnfoldResponse.h"

//Local FullJR (MainAnalysis) dependencies
#include "MainAnalysis/include/cutPropagator.h"
#include "MainAnalysis/include/doLocalDebug.h"
#include "MainAnalysis/include/flatWeightReader.h"
#include "MainAnalysis/include/smallOrLargeR.h"

//Non-local FullJR (Utility, etc.) dependencies
#include "Utility/include/checkMakeDir.h"
#include "Utility/include/cppWatch.h"
#include "Utility/include/doGlobalDebug.h"
#include "Utility/include/etaPhiFunc.h"
#include "Utility/include/getLinBins.h"
#include "Utility/include/goodGlobalSelection.h"
#include "Utility/include/histDefUtility.h"
#include "Utility/include/mntToXRootdFileString.h"
#include "Utility/include/ncollFunctions_5TeV.h"
#include "Utility/include/plotUtilities.h"
#include "Utility/include/scaleErrorTool.h"
#include "Utility/include/specialHYDJETEventExclude.h"
#include "Utility/include/stringUtil.h"
#include "Utility/include/CustomAssert.h"
#include "Utility/include/CommandLine.h"

template <typename T> std::string to_string_with_precision(const T a_value, const int n);
int FindBin(int N, double Min, double Max, double V);
int GetBin(int C = 0, int I = 0, int R = 0, int E = 0, int S = 0);
int makeJetResponseTree(const std::string inName, const std::string treeName,
   bool isPP = false, double inEntryFrac = 1.,
   const bool doRooResponse = false, const bool doSystReduced = false);
int main(int argc, char *argv[]);

template <typename T> std::string to_string_with_precision(const T a_value, const int n)
{
   std::ostringstream out;
   out << std::setprecision(n) << a_value;
   return out.str();
}

int FindBin(int N, double Min, double Max, double V)
{
   int B = (V - Min) / (Max - Min) * N;
   if(B < 0)
      B = 0;
   if(B >= N)
      B = N - 1;
   return B;
}

int GetBin(int C, int I, int R, int E, int S)
{
   const int NC = 4;
   const int NI = 1;
   const int NR = 2;
   const int NE = 5;
   // const int NS = 14;

   return (((S * NE + E) * NR + R) * NI + I) * NC + C;
}

int makeJetResponseTree(const std::string inName, const std::string treeName,
   bool isPP, double inEntryFrac,
   const bool doRooResponse, const bool doSystReduced)
{
   if(doLocalDebug || doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

   cppWatch totalRunWatch;
   cppWatch fileLoopWatch;
   cppWatch writeLoopWatch;
   cppWatch deleteLoopWatch;
   totalRunWatch.start();

   std::vector<std::string> fileList;
   std::vector<double> pthats;
   std::vector<double> pthatWeights;
   std::string workDir = ".";

   if(inName.find(".root") != std::string::npos)
   {
      fileList.push_back(inName);
      pthats.push_back(1.);
      pthatWeights.push_back(1.);
   }
   else if(inName.find(".txt") != std::string::npos)
   {
      std::ifstream file(inName.c_str());
      std::string tempStr;

      while(std::getline(file, tempStr))
      {
         while(tempStr.find(" ") != std::string::npos){tempStr.replace(tempStr.find(" "), 1, "");}
         if(tempStr.size() == 0) continue;     
         if(tempStr.find(".root") != std::string::npos) fileList.push_back(tempStr);
         else{
            if(tempStr.substr(0, std::string("PTHAT=").size()).find("PTHAT=") != std::string::npos){
               tempStr.replace(0, std::string("PTHAT=").size(), "");
               while(tempStr.find(",") != std::string::npos){
                  pthats.push_back(std::stod(tempStr.substr(0, tempStr.find(","))));
                  tempStr.replace(0, tempStr.find(",")+1, "");
               }
               if(tempStr.size() != 0) pthats.push_back(std::stod(tempStr));
            }
            else if(tempStr.substr(0, std::string("PTHATWEIGHTS=").size()).find("PTHATWEIGHTS=") != std::string::npos){
               tempStr.replace(0, std::string("PTHATWEIGHTS=").size(), "");
               while(tempStr.find(",") != std::string::npos){
                  pthatWeights.push_back(std::stod(tempStr.substr(0, tempStr.find(","))));
                  tempStr.replace(0, tempStr.find(",")+1, "");
               }
               if(tempStr.size() != 0) pthatWeights.push_back(std::stod(tempStr));	
            }
            else if(tempStr.substr(0, std::string("ISPP=").size()).find("ISPP=") != std::string::npos){
               tempStr.replace(0,tempStr.find("=")+1, "");
               isPP = std::stoi(tempStr);
            }
            else if(tempStr.substr(0, std::string("WORKDIR=").size()).find("WORKDIR=") != std::string::npos)
            {
               tempStr.replace(0,tempStr.find("=")+1, "");
               workDir = tempStr;
            }
            else std::cout << "WARNING: Line in \'" << inName << "\', \'" << tempStr << "\' is invalid. check input" << std::endl;
         }     
      }

      file.close();
   }
   else
      Assert(false, "Input filename invalid");

   Assert(fileList.size() > 0,     "is gives no valid root files.");
   Assert(pthats.size() > 0,       "contains no pthat list.");
   Assert(pthatWeights.size() > 0, "contains no pthatWeights list.");

   if(doLocalDebug || doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

   //Post possible returns, start your random number generator
   TRandom3 *randGen_p = new TRandom3(0);

   unsigned int pos = 0;
   while(pos < pthats.size())
   {
      bool isMoved = false;
      for(unsigned int pI = pos+1; pI < pthats.size(); ++pI)
      {
         if(pthats[pI] < pthats[pos])
         {
            double pthatTemp = pthats[pos];
            double pthatWeightTemp = pthatWeights[pos];

            pthats[pos] = pthats[pI];
            pthatWeights[pos] = pthatWeights[pI];

            pthats[pI] = pthatTemp;
            pthatWeights[pI] = pthatWeightTemp;

            isMoved = true;
         }
      }

      if(!isMoved) pos++;
   }

   std::cout << "Pthats and weights: " << std::endl;
   for(unsigned int pI = 0; pI < pthats.size(); ++pI){
      std::cout << " " << pI << "/" << pthats.size() << ": " << pthats[pI] << ", " << pthatWeights[pI] << std::endl;
   }

   std::cout << "isPP: " << isPP << std::endl;

   if(doLocalDebug || doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

   TFile *inFile_p = nullptr;

   Int_t posR4Temp = -1;

   Int_t rValI = getRVal(treeName);
   Double_t rValD = getRVal(treeName) / 10;

   smallOrLargeR rReader;

   if(treeName.find("akCs4") != std::string::npos && posR4Temp < 0)
      posR4Temp = 0;
   else if(treeName.find("ak4") != std::string::npos && posR4Temp < 0)
      posR4Temp = 0;
   
   const Int_t posR4 = posR4Temp;
   const Int_t posGeneral = 0;

   std::cout << "PosR4: " << posR4 << std::endl;
   std::cout << "PosGeneral: " << posGeneral << std::endl;

   const Int_t nCentBinsPerma = 4;
   const Int_t centBinsLowPerma[nCentBinsPerma] = {0, 10, 30, 50};
   const Int_t centBinsHiPerma[nCentBinsPerma] = {10, 30, 50, 90};

   Int_t nCentBinsTemp = nCentBinsPerma;

   const Int_t nMaxCentBins = 4;
   const Int_t nCentBins = nCentBinsTemp;

   Assert(nCentBins <= nMaxCentBins, "nCentBins is larger than nMaxCentBins.");

   std::vector<Int_t> centBinsLow, centBinsHi;
   if(isPP && false)
   { 
      centBinsLow.push_back(0);
      centBinsHi.push_back(100);
   }
   else
   {
      for(Int_t cI = 0; cI < nCentBinsPerma; ++cI)
      {
         centBinsLow.push_back(centBinsLowPerma[cI]);
         centBinsHi.push_back(centBinsHiPerma[cI]);
      }
   }

   const Int_t nPthatBins = 100;
   const Float_t pthatLow = pthats[0];
   const Float_t pthatHi = pthats[pthats.size()-1]*2.;
   Double_t pthatBins[nPthatBins+1];
   getLinBins(pthatLow, pthatHi, nPthatBins, pthatBins);

   if(doLocalDebug || doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;  

   const Int_t nCentBins2 = 100;
   const Float_t centBinsLow2 = 0;
   const Float_t centBinsHi2 = 100;
   Double_t centBins2[nCentBins2+1];
   getLinBins(centBinsLow2, centBinsHi2, nCentBins2, centBins2);

   const double fracParaFills = 0.1;

   TDatime *date = new TDatime();
   //  const std::string dateStr = std::to_string(date->GetDate()) + "_" + std::to_string(date->GetHour());
   const std::string dateStr = std::to_string(date->GetDate());
   delete date;
   const std::string fullPath = std::getenv("FULLJRDIR");

   std::string inDir = workDir;
   inDir = inDir.substr(0, inDir.rfind("/"));

   std::string outFileName = inName;
   while(outFileName.find("/") != std::string::npos){outFileName.replace(0, outFileName.find("/")+1, "");}
   if(outFileName.find(".root") != std::string::npos) outFileName.replace(outFileName.find(".root"), std::string(".root").size(), "");
   else if(outFileName.find(".txt") != std::string::npos) outFileName.replace(outFileName.find(".txt"), std::string(".txt").size(), "");
   outFileName = "output/" + dateStr + "/" + outFileName + "_FracNEntries" + prettyString(inEntryFrac, 2, true) + "_JetResponse_";

   if(inDir.find("/mnt") == std::string::npos){
      outFileName = inDir + "/" + outFileName;
      if(checkFile(outFileName + dateStr + ".root")) outFileName = outFileName + "UPDATED_";
      outFileName = outFileName + dateStr + ".root";
      checkMakeDir(inDir + "/output");
      checkMakeDir(inDir + "/output/" + dateStr);
   }
   else{
      if(checkFile(outFileName + dateStr + ".root")) outFileName = outFileName + "UPDATED_";
      outFileName = outFileName + dateStr + ".root";
      checkMakeDir("output");
      checkMakeDir("output/" + dateStr);
   }

   const Double_t jtAbsEtaMax = 2.;

   const Int_t nJtAbsEtaBins = 5;
   const Double_t jtAbsEtaBinsLow[nJtAbsEtaBins] = {0.0, 0.5, 1.0, 1.5, 0.0};
   const Double_t jtAbsEtaBinsHi[nJtAbsEtaBins] = {0.5, 1.0, 1.5, 2.0, 2.0};

   const Double_t jecVarMC = 0.02;
   const Double_t jerVarMC = 0.07;

   const Double_t jecVarData = 0.02;
   const Int_t nResponseMod = 2;
   const Double_t responseMod[nResponseMod] = {0.00, 0.10};
   const Double_t jerVarData[nResponseMod] = {0.15, 0.10};

   const int nR = rReader.GetNR();
   std::vector<int> rVals = rReader.GetRVals();

   const int nPtBinsMax = 200;
   int nGenJtPtBins[nMaxCentBins];
   double genJtPtBins[nMaxCentBins][nPtBinsMax+1];
   int nRecoJtPtBins[nMaxCentBins];
   double recoJtPtBins[nMaxCentBins][nPtBinsMax+1];

   for(Int_t cI = 0; cI < nCentBins; ++cI)
   {
      std::string centStr = "Cent" + std::to_string(centBinsLow[cI]) + "to" + std::to_string(centBinsHi[cI]);

      int B = GetBin(cI);
      nGenJtPtBins[B] = rReader.GetGenNBinsFromRValCent(rValI, centStr);
      nRecoJtPtBins[B] = rReader.GetRecoNBinsFromRValCent(rValI, centStr);

      rReader.GetGenBinsFromRValCent(rValI, centStr, genJtPtBins[B]);
      rReader.GetRecoBinsFromRValCent(rValI, centStr, recoJtPtBins[B]);      
   }

   Int_t bigJetRecoTrunc = -1;
   for(Int_t jI = 0; jI < nRecoJtPtBins[0]; ++jI)
   {
      Double_t val = (recoJtPtBins[0][jI] + recoJtPtBins[0][jI+1]) / 2.;
      if(val > 200.)
      {
         bigJetRecoTrunc = jI + 1;
         break;
      }
   }

   Double_t minJtPtCut = 20;
   Double_t multiJtPtCut = 50;
   Int_t recoTruncPos = 1;

   bool isBigJt = (treeName.find("ak8") != std::string::npos)
      || (treeName.find("ak10") != std::string::npos)
      || (treeName.find("akCs8") != std::string::npos)
      || (treeName.find("akCs10") != std::string::npos);
   if(isBigJt)
   {
      multiJtPtCut = 100;
      minJtPtCut = 20.;
      recoTruncPos = bigJetRecoTrunc;
   }

   //FULL ID taken from here https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2016, FullLight and FullTight correspond to previous levels of cuts + the Loose and TightLepVeto versions, respectively
   /*
      const Int_t nID = 5;
      const std::string idStr[nID] = {"NoID", "LightMUID", "LightMUAndCHID", "FullLight", "FullTight"};
      const Double_t jtPfCHMFCutLow[nID] = {0.0, 0.0, 0.00, 0.00, 0.00};
      const Double_t jtPfCHMFCutHi[nID] = {1.0, 1.0, 0.80, 0.80, 0.80}; // THIS CUT IS DROPPED FROM 0.9 to 0.8 to account for a few jets still screwing up response. This is a temporary solution. Please instead check the event for a compatible nonsense muon that this cut is meant to address
      const Double_t jtPfMUMFCutLow[nID] = {0.0, 0.0, 0.00, 0.00, 0.00};
      const Double_t jtPfMUMFCutHi[nID] = {1.0, 0.60, 0.60, 0.60, 0.60};
      const Double_t jtPfNHFCutLow[nID] = {0.0, 0.0, 0.0, 0.0, 0.0};
      const Double_t jtPfNHFCutHi[nID] = {1.0, 1.0, 1.0, 0.99, 0.90};
      const Double_t jtPfNEFCutLow[nID] = {0.0, 0.0, 0.0, 0.0, 0.0};
      const Double_t jtPfNEFCutHi[nID] = {1.0, 1.0, 1.0, 0.99, 0.90};
      const Int_t jtPfMinMult[nID] = {1, 1, 1, 2, 2};
      const Double_t jtPfMUFCutLow[nID] = {0.0, 0.0, 0.0, 0.0, 0.0};
      const Double_t jtPfMUFCutHi[nID] = {1.0, 1.0, 1.0, 1.0, 0.8};
      const Double_t jtPfCHFCutLow[nID] = {0.0, 0.0, 0.0, 0.0000001, 0.0000001};
      const Double_t jtPfCHFCutHi[nID] = {1.0, 1.0, 1.0, 1.0, 1.0};
      const Int_t jtPfMinChgMult[nID] = {0, 0, 0, 1, 1};
      const Double_t jtPfCEFCutLow[nID] = {0.0, 0.0, 0.0, 0.0, 0.0};
      const Double_t jtPfCEFCutHi[nID] = {1.0, 1.0, 1.0, 0.99, 0.90};
      */
   const Int_t nID = 1;
   const std::string idStr[nID] = {"LightMUAndCHID"};
   const Double_t jtPfCHMFCutLow[nID] = {0.00,};
   const Double_t jtPfCHMFCutHi[nID] = {0.80}; // THIS CUT IS DROPPED FROM 0.9 to 0.8 to account for a few jets still screwing up response. This is a temporary solution. Please instead check the event for a compatible nonsense muon that this cut is meant to address
   const Double_t jtPfMUMFCutLow[nID] = {0.00};
   const Double_t jtPfMUMFCutHi[nID] = {0.60};
   const Double_t jtPfNHFCutLow[nID] = {0.0};
   const Double_t jtPfNHFCutHi[nID] = {1.0};
   const Double_t jtPfNEFCutLow[nID] = {0.0};
   const Double_t jtPfNEFCutHi[nID] = {1.0};
   const Int_t jtPfMinMult[nID] = {1};
   const Double_t jtPfMUFCutLow[nID] = {0.0};
   const Double_t jtPfMUFCutHi[nID] = {1.0};
   const Double_t jtPfCHFCutLow[nID] = {0.0};
   const Double_t jtPfCHFCutHi[nID] = {1.0};
   const Int_t jtPfMinChgMult[nID] = {0};
   const Double_t jtPfCEFCutLow[nID] = {0.0};
   const Double_t jtPfCEFCutHi[nID] = {1.0};


   enum SystType {None, JECUpMC, JECDownMC, JECUpData, JECDownData, JECUpUE, JECDownUE, JERMC, JERData, Fake, PriorUp1PowerPthat, PriorDown1PowerPthat, PriorFlat, MatrixStat};

   const Int_t nSyst = 14;
   const std::string systStr[nSyst] = {"", "JECUpMC", "JECDownMC", "JECUpData", "JECDownData", "JECUpUE", "JECDownUE", "JERMC", "JERData", "Fake", "PriorUp1PowerPthat", "PriorDown1PowerPthat", "PriorFlat", "MatrixStat"};
   const SystType systType[nSyst] = {None, JECUpMC, JECDownMC, JECUpData, JECDownData, JECUpUE, JECDownUE, JERMC, JERData, Fake, PriorUp1PowerPthat, PriorDown1PowerPthat, PriorFlat, MatrixStat};
   Int_t priorFlatPos = -1;
   for(Int_t sI = 0; sI < nSyst; ++sI)
   {
      if(systType[sI] == PriorFlat)
      {
         priorFlatPos = sI;
         break;
      }
   }

   Int_t nSystReduced = nSyst;
   if(doSystReduced) nSystReduced = 3;

   const std::string rcDiffFileName = "MainAnalysis/tables/rcDifferences_20180418.txt";
   scaleErrorTool scaleErr((fullPath + "/" + rcDiffFileName).c_str());
   scaleErr.Init();

   cutPropagator cutProp;
   cutProp.Clean();
   cutProp.SetInFileNames({inName});
   cutProp.SetInFullFileNames(fileList);
   cutProp.SetIsPP(isPP);
   cutProp.SetRCDiffFileName(rcDiffFileName);
   cutProp.SetJtAbsEtaMax(jtAbsEtaMax);
   cutProp.SetJECVarMC(jecVarMC);
   cutProp.SetJERVarMC(jerVarMC);
   cutProp.SetJECVarData(jecVarData);
   cutProp.SetNResponseMod(nResponseMod);
   cutProp.SetResponseMod(nResponseMod, responseMod);
   cutProp.SetJERVarData(nResponseMod, jerVarData);
   cutProp.SetNR(nR);
   cutProp.SetRVals(rVals);

   for(Int_t cI = 0; cI < nCentBins; ++cI)
   {
      std::string centStr = "Cent" + std::to_string(centBinsLow[cI]) + "to" + std::to_string(centBinsHi[cI]);
      cutProp.SetGenPtBins(rValI, centStr, nGenJtPtBins[cI], genJtPtBins[cI]);
      cutProp.SetRecoPtBins(rValI, centStr, nRecoJtPtBins[cI], recoJtPtBins[cI]);
   }

   cutProp.SetNGeneralBins(rReader.GetNGeneralBins());
   cutProp.SetGeneralBins(rReader.GetGeneralBins());

   cutProp.SetNJtAlgos(1);
   cutProp.SetJtAlgos({treeName});
   cutProp.SetMinJtPtCut({minJtPtCut});
   cutProp.SetMultiJtPtCut({multiJtPtCut});
   cutProp.SetRecoTruncPos({recoTruncPos});
   cutProp.SetNJtAbsEtaBins(nJtAbsEtaBins);
   cutProp.SetJtAbsEtaBinsLow(nJtAbsEtaBins, jtAbsEtaBinsLow);
   cutProp.SetJtAbsEtaBinsHi(nJtAbsEtaBins, jtAbsEtaBinsHi);
   cutProp.SetNPthats(pthats.size());
   cutProp.SetPthats(pthats);
   cutProp.SetPthatWeights(pthatWeights);
   cutProp.SetNCentBins(nCentBins);
   cutProp.SetCentBinsLow(centBinsLow);
   cutProp.SetCentBinsHi(centBinsHi);
   cutProp.SetNID(nID);
   cutProp.SetIdStr(nID, idStr);
   cutProp.SetJtPfCHMFCutLow(nID, jtPfCHMFCutLow);
   cutProp.SetJtPfCHMFCutHi(nID, jtPfCHMFCutHi);
   cutProp.SetJtPfMUMFCutLow(nID, jtPfMUMFCutLow);
   cutProp.SetJtPfMUMFCutHi(nID, jtPfMUMFCutHi);
   cutProp.SetJtPfNHFCutLow(nID, jtPfNHFCutLow);
   cutProp.SetJtPfNHFCutHi(nID, jtPfNHFCutHi);
   cutProp.SetJtPfNEFCutLow(nID, jtPfNEFCutLow);
   cutProp.SetJtPfNEFCutHi(nID, jtPfNEFCutHi);
   cutProp.SetJtPfMUFCutLow(nID, jtPfMUFCutLow);
   cutProp.SetJtPfMUFCutHi(nID, jtPfMUFCutHi);
   cutProp.SetJtPfCHFCutLow(nID, jtPfCHFCutLow);
   cutProp.SetJtPfCHFCutHi(nID, jtPfCHFCutHi);
   cutProp.SetJtPfCEFCutLow(nID, jtPfCEFCutLow);
   cutProp.SetJtPfCEFCutHi(nID, jtPfCEFCutHi);
   cutProp.SetJtPfMinMult(nID, jtPfMinMult);
   cutProp.SetJtPfMinChgMult(nID, jtPfMinChgMult);
   cutProp.SetNSyst(nSystReduced);
   cutProp.SetSystStr(nSystReduced, systStr);

   const std::string flatWeightNamePbPb = "MainAnalysis/tables/Pythia6_Dijet_pp502_Hydjet_Cymbal_MB_PbPb_MCDijet_20180521_ExcludeTop4_ExcludeToFrac_Frac0p7_Full_5Sigma_20180608_SVM_FlatGenJetResponse_20180828_9.root";

   const std::string flatWeightNamePP = "MainAnalysis/tables/Pythia6_Dijet_pp502_MCDijet_20180712_ExcludeTop4_ExcludeToFrac_Frac0p7_Full_5Sigma_20180712_SVM_FlatGenJetResponse_20180828_9.root";

   std::string flatWeightName = fullPath + "/";
   if(isPP) flatWeightName = flatWeightName + flatWeightNamePP;
   else flatWeightName = flatWeightName + flatWeightNamePbPb;
   std::cout << "CutProp: " << cutProp.GetNJtAlgos() << std::endl;
   flatWeightReader flatWeight(flatWeightName, cutProp);

   const Int_t nResponseBins = 300;
   Double_t responseBins[nResponseBins+1];
   getLinBins(0., 9., nResponseBins, responseBins);

   TFile *outFile_p = new TFile(outFileName.c_str(), "RECREATE");
   //Following two lines necessary else you get absolutely clobbered on deletion timing. See threads here:
   //https://root-forum.cern.ch/t/closing-root-files-is-dead-slow-is-this-a-normal-thing/5273/16
   //https://root-forum.cern.ch/t/tfile-speed/17549/25
   //Bizarre
   outFile_p->SetBit(TFile::kDevNull);
   TH1::AddDirectory(kFALSE);

   TDirectory *generalDir_p = (TDirectory *)outFile_p->mkdir("generalHistDir");
   generalDir_p->cd();
   TH1D *pthat_h = new TH1D("pthat_h", ";p_{T} Hat;Counts (Unweighted)", nPthatBins, pthatLow, pthatHi);
   TH1D *pthatWeighted_h = new TH1D("pthatWeighted_h", ";p_{T} Hat;Counts (Weighted)", nPthatBins, pthatLow, pthatHi);
   TH1D *pthatFullWeighted_h = new TH1D("pthatFullWeighted_h", ";p_{T} Hat;Counts (Full Weighted)", nPthatBins, pthatLow, pthatHi);
   TH1D *pthatFullRatio_h = new TH1D("pthatFullRatio_h", ";p_{T} Hat;Ratio", nPthatBins, pthatLow, pthatHi);

   centerTitles({pthat_h, pthatWeighted_h, pthatFullWeighted_h, pthatFullRatio_h});
   setSumW2({pthat_h, pthatWeighted_h, pthatFullWeighted_h, pthatFullRatio_h});

   TH1D *centrality_h = NULL;
   TH1D *centralityWeighted_h = NULL;
   TH1D *centralityFullWeighted_h = NULL;
   TH1D *centralityFullRatio_h = NULL;

   if(isPP == false)
   {
      centrality_h = new TH1D("centrality_h", ";Centrality (%);Counts (Unweighted)", nCentBins2, centBins2);
      centralityWeighted_h = new TH1D("centralityWeighted_h", ";Centrality (%);Counts (Weighted)", nCentBins2, centBins2);
      centralityFullWeighted_h = new TH1D("centralityFullWeighted_h", ";Centrality (%);Counts (Full Weighted)", nCentBins2, centBins2);
      centralityFullRatio_h = new TH1D("centralityFullRatio_h", ";Centrality (%);Ratio", nCentBins2, centBins2);

      centerTitles({centrality_h, centralityWeighted_h, centralityFullWeighted_h, centralityFullRatio_h});
      setSumW2({centrality_h, centralityWeighted_h, centralityFullWeighted_h, centralityFullRatio_h});
   }

   Int_t nNormBins = 200;
   Double_t normBins[nNormBins+1];
   getLinBins(0,3,nNormBins, normBins);

   Int_t nGeneralBins = rReader.GetNGeneralBins();

   Double_t generalBins[nPtBinsMax+1];
   rReader.GetGeneralBins(nGeneralBins, generalBins);

   Double_t minGenJtPt = 999999999;
   Double_t minRecoJtPt = 999999999;

   std::string dirName = treeName;
   dirName = dirName.substr(0, dirName.find("/"));

   for(Int_t bIX = 0; bIX < nGeneralBins+1; ++bIX)
   {
      if(minRecoJtPt > generalBins[bIX]) minRecoJtPt = generalBins[bIX];
      if(minGenJtPt > generalBins[bIX]) minGenJtPt = generalBins[bIX];
   }

   for(Int_t cI = 0; cI < nCentBins; ++cI)
   {
      for(Int_t bIX = 0; bIX < nRecoJtPtBins[cI] + 1; ++bIX)
         if(minRecoJtPt > recoJtPtBins[cI][bIX]) minRecoJtPt = recoJtPtBins[cI][bIX];
      for(Int_t bIX = 0; bIX < nGenJtPtBins[cI] + 1; ++bIX)
         if(minGenJtPt > genJtPtBins[cI][bIX]) minGenJtPt = genJtPtBins[cI][bIX]; 
   }

   std::cout << "MINGENJTPT: " << minGenJtPt << ", " << minJtPtCut << std::endl;

   TDirectory *dir_p = nullptr;

   TH1D *deltaPtOrig_h[nMaxCentBins];
   TH1D *deltaPhiOrig_h[nMaxCentBins];
   TH1D *deltaEtaOrig_h[nMaxCentBins];
   TH1D *deltaROrig_h[nMaxCentBins];

   TH1D *deltaPtReplace_h[nMaxCentBins];
   TH1D *deltaPhiReplace_h[nMaxCentBins];
   TH1D *deltaEtaReplace_h[nMaxCentBins];
   TH1D *deltaRReplace_h[nMaxCentBins];

   TH1D *deltaPtNotReplace_h[nMaxCentBins];
   TH1D *deltaPhiNotReplace_h[nMaxCentBins];
   TH1D *deltaEtaNotReplace_h[nMaxCentBins];
   TH1D *deltaRNotReplace_h[nMaxCentBins];

   TH1D *recoJtPt_h[nMaxCentBins*nID*nResponseMod*nJtAbsEtaBins];
   TH1D *recoJtPt_RecoTrunc_h[nMaxCentBins*nID*nResponseMod*nJtAbsEtaBins];
   TH1D *recoJtPt_NoTruth_h[nMaxCentBins*nID*nResponseMod*nJtAbsEtaBins];

   TH1D *recoJtPt_ParaFills_h[nMaxCentBins*nID*nResponseMod*nJtAbsEtaBins];
   TH1D *recoJtPt_RecoTrunc_ParaFills_h[nMaxCentBins*nID*nResponseMod*nJtAbsEtaBins];
   TH1D *recoJtPt_NoTruth_ParaFills_h[nMaxCentBins*nID*nResponseMod*nJtAbsEtaBins];

   TH1D *recoJtPt_GoodGen_h[nMaxCentBins*nID*nResponseMod*nJtAbsEtaBins*nSyst];
   TH1D *recoJtPt_GoodGen_ParaFills_h[nMaxCentBins*nID*nResponseMod*nJtAbsEtaBins*nSyst];

   TH1D *recoJtPt_General_h[nMaxCentBins*nID*nResponseMod*nJtAbsEtaBins*nSyst];
   TH1D *recoJtPt_General_AllReco_h[nMaxCentBins*nID*nResponseMod*nJtAbsEtaBins*nSyst];
   TH1D *recoJtPt_General_ParaFills_h[nMaxCentBins*nID*nResponseMod*nJtAbsEtaBins*nSyst];

   TH1D *genJtPt_h[nMaxCentBins*nID*nResponseMod*nJtAbsEtaBins];
   TH1D *genJtPt_ParaFills_h[nMaxCentBins*nID*nResponseMod*nJtAbsEtaBins];
   TH1D *genJtPt_All_h[nMaxCentBins*nID*nResponseMod*nJtAbsEtaBins];
   TH1D *genJtPt_GoodReco_h[nMaxCentBins*nID*nResponseMod*nJtAbsEtaBins*nSyst];
   TH1D *genJtPt_GoodReco_ParaFills_h[nMaxCentBins*nID*nResponseMod*nJtAbsEtaBins*nSyst];
   TH1D *genJtPt_General_h[nMaxCentBins*nID*nResponseMod*nJtAbsEtaBins*nSyst];
   TH1D *genJtPt_General_AllGen_h[nMaxCentBins*nID*nResponseMod*nJtAbsEtaBins*nSyst];
   TH1D *genJtPt_General_ParaFills_h[nMaxCentBins*nID*nResponseMod*nJtAbsEtaBins*nSyst];
   TH2D *response_RecoGenSymm_h[nMaxCentBins*nID*nResponseMod*nJtAbsEtaBins*nSyst];
   TH2D *response_RecoGenAsymm_h[nMaxCentBins*nID*nResponseMod*nJtAbsEtaBins*nSyst];

   TH2D *response_General_h[nMaxCentBins*nID*nResponseMod*nJtAbsEtaBins*nSyst];
   TH2D *response_General_Half1_h[nMaxCentBins*nID*nResponseMod*nJtAbsEtaBins*nSyst];
   TH2D *response_General_Half2_h[nMaxCentBins*nID*nResponseMod*nJtAbsEtaBins*nSyst];
   TH2D *response_General_ParaFills_h[nMaxCentBins*nID*nResponseMod*nJtAbsEtaBins*nSyst];
   TH2D *responseNorm_General_h[nMaxCentBins*nID*nResponseMod*nJtAbsEtaBins*nSyst];

   RooUnfoldResponse *rooResponse_RecoGenAsymm_h[nMaxCentBins*nID*nResponseMod*nJtAbsEtaBins*nSyst];
   TH1D *genJtPt_CheckPriorFlat_h[nMaxCentBins*nID*nResponseMod*nJtAbsEtaBins];

   outFile_p->cd();
   dirName = treeName;
   dirName = dirName.substr(0, dirName.find("/"));

   dir_p = (TDirectory *)outFile_p->mkdir(dirName.c_str());

   for(Int_t cI = 0; cI < nCentBins; ++cI)
   {
      std::string centStr = "Cent" + std::to_string(centBinsLow[cI]) + "to" + std::to_string(centBinsHi[cI]);      
      if(isPP) centStr = "PP_" + centStr;
      else centStr = "PbPb_" + centStr;

      int B = GetBin(cI);

      deltaPtOrig_h[cI] = new TH1D(("deltaPtOrig_" + dirName + "_" + centStr + "_h").c_str(), ";#Delta pT/Reco Jet pT;Counts", 100, 0, 2);
      deltaPhiOrig_h[cI] = new TH1D(("deltaPhiOrig_" + dirName + "_" + centStr + "_h").c_str(), ";#Delta #phi;Counts", 100, -0.5, 0.5);
      deltaEtaOrig_h[cI] = new TH1D(("deltaEtaOrig_" + dirName + "_" + centStr + "_h").c_str(), ";#Delta #eta;Counts", 100, -0.5, 0.5);
      deltaROrig_h[cI] = new TH1D(("deltaROrig_" + dirName + "_" + centStr + "_h").c_str(), ";#Delta R;Counts", 100, -0.5, 0.5);

      deltaPtReplace_h[cI] = new TH1D(("deltaPtReplace_" + dirName + "_" + centStr + "_h").c_str(), ";#Delta pT/Reco Jet pT;Counts", 100, 0, 2);
      deltaPhiReplace_h[cI] = new TH1D(("deltaPhiReplace_" + dirName + "_" + centStr + "_h").c_str(), ";#Delta #phi;Counts", 100, -0.5, 0.5);
      deltaEtaReplace_h[cI] = new TH1D(("deltaEtaReplace_" + dirName + "_" + centStr + "_h").c_str(), ";#Delta #eta;Counts", 100, -0.5, 0.5);
      deltaRReplace_h[cI] = new TH1D(("deltaRReplace_" + dirName + "_" + centStr + "_h").c_str(), ";#Delta R;Counts", 100, -0.5, 0.5);

      deltaPtNotReplace_h[cI] = new TH1D(("deltaPtNotReplace_" + dirName + "_" + centStr + "_h").c_str(), ";#Delta pT/Reco Jet pT;Counts", 100, 0, 2);
      deltaPhiNotReplace_h[cI] = new TH1D(("deltaPhiNotReplace_" + dirName + "_" + centStr + "_h").c_str(), ";#Delta #phi;Counts", 100, -0.5, 0.5);
      deltaEtaNotReplace_h[cI] = new TH1D(("deltaEtaNotReplace_" + dirName + "_" + centStr + "_h").c_str(), ";#Delta #eta;Counts", 100, -0.5, 0.5);
      deltaRNotReplace_h[cI] = new TH1D(("deltaRNotReplace_" + dirName + "_" + centStr + "_h").c_str(), ";#Delta R;Counts", 100, -0.5, 0.5);

      for(Int_t iI = 0; iI < nID; ++iI)
      {
         for(Int_t mI = 0; mI < nResponseMod; ++mI)
         {
            std::string resStr = "ResponseMod" + prettyString(responseMod[mI], 2, true);

            for(Int_t aI = 0; aI < nJtAbsEtaBins; ++aI)
            {
               const std::string jtAbsEtaStr = "AbsEta" + prettyString(jtAbsEtaBinsLow[aI], 1, true)
                  + "to" + prettyString(jtAbsEtaBinsHi[aI], 1, true);

               std::string nameStr = dirName + "_" + centStr + "_" + idStr[iI] + "_" + resStr
                  + "_" + jtAbsEtaStr;

               int BB = GetBin(cI, iI, mI, aI);
               recoJtPt_h[BB] = new TH1D(("recoJtPt_" + nameStr + "_h").c_str(),
                     ";Reco. Jet p_{T};Counts (Weighted)", nRecoJtPtBins[B], recoJtPtBins[B]);
               recoJtPt_RecoTrunc_h[BB] = new TH1D(("recoJtPt_" + nameStr + "_RecoTrunc_h").c_str(),
                     ";Reco. Jet p_{T};Counts (Weighted)", nRecoJtPtBins[B], recoJtPtBins[B]);
               recoJtPt_NoTruth_h[BB] = new TH1D(("recoJtPt_" + nameStr + "_NoTruth_h").c_str(),
                     ";Reco. Jet p_{T};Counts (Weighted)", nRecoJtPtBins[B], recoJtPtBins[B]);

               recoJtPt_ParaFills_h[BB] = new TH1D(("recoJtPt_" + nameStr + "_ParaFills_h").c_str(),
                     ";Reco. Jet p_{T};Counts (Weighted)", nRecoJtPtBins[B], recoJtPtBins[B]);
               recoJtPt_RecoTrunc_ParaFills_h[BB] = new TH1D(("recoJtPt_" + nameStr + "_RecoTrunc_ParaFills_h").c_str(),
                     ";Reco. Jet p_{T};Counts (Weighted)", nRecoJtPtBins[B], recoJtPtBins[B]);
               recoJtPt_NoTruth_ParaFills_h[BB] = new TH1D(("recoJtPt_" + nameStr + "_NoTruth_ParaFills_h").c_str(),
                     ";Reco. Jet p_{T};Counts (Weighted)", nRecoJtPtBins[B], recoJtPtBins[B]);

               genJtPt_h[BB] = new TH1D(("genJtPt_" + nameStr + "_h").c_str(),
                     ";Gen. Jet p_{T};Counts (Weighted)", nGenJtPtBins[B], genJtPtBins[B]);

               genJtPt_ParaFills_h[BB] = new TH1D(("genJtPt_" + nameStr + "_ParaFills_h").c_str(),
                     ";Gen. Jet p_{T};Counts (Weighted)", nGenJtPtBins[B], genJtPtBins[B]);;	    	  

               genJtPt_All_h[BB] = new TH1D(("genJtPt_" + nameStr + "_All_h").c_str(),
                     ";Gen. Jet p_{T};Counts (Weighted)", nGenJtPtBins[B], genJtPtBins[B]);
               genJtPt_CheckPriorFlat_h[BB] = new TH1D(("rooResponse_" + nameStr + "_CheckPriorFlat_h").c_str(),
                     ";Gen. p_{T} (Weighted Flat);Counts (Weighted)", nGenJtPtBins[B], genJtPtBins[B]);

               std::vector<TH1*> tempVect = {recoJtPt_RecoTrunc_h[BB], recoJtPt_NoTruth_h[BB], recoJtPt_h[BB], recoJtPt_RecoTrunc_ParaFills_h[BB], recoJtPt_NoTruth_ParaFills_h[BB], recoJtPt_ParaFills_h[BB], genJtPt_h[BB], genJtPt_ParaFills_h[BB], genJtPt_All_h[BB], genJtPt_CheckPriorFlat_h[BB]};
               centerTitles(tempVect);
               setSumW2(tempVect);	    	    	    	  

               for(Int_t sI = 0; sI < nSystReduced; ++sI)
               {
                  std::string tempSysStr =  "_" + systStr[sI] + "_";	 
                  while(tempSysStr.find("__") != std::string::npos)
                     tempSysStr.replace(tempSysStr.find("__"), 2, "_");

                  int BB = GetBin(cI, iI, mI, aI, sI);

                  recoJtPt_GoodGen_h[BB] = new TH1D(("recoJtPt_" + nameStr + tempSysStr + "GoodGen_h").c_str(),
                        ";Reco. Jet p_{T};Counts (Weighted)", nRecoJtPtBins[B], recoJtPtBins[B]);

                  recoJtPt_GoodGen_ParaFills_h[BB] = new TH1D(("recoJtPt_" + nameStr + tempSysStr + "GoodGen_ParaFills_h").c_str(),
                        ";Reco. Jet p_{T};Counts (Weighted)", nRecoJtPtBins[B], recoJtPtBins[B]);

                  recoJtPt_General_h[BB] = new TH1D(("recoJtPt_" + nameStr + tempSysStr + "General_h").c_str(),
                        ";Reco. Jet p_{T};Counts (Weighted)", 148, 20, 1500);
                  recoJtPt_General_AllReco_h[BB] = new TH1D(("recoJtPt_" + nameStr + tempSysStr + "General_AllReco_h").c_str(),
                        ";Reco. Jet p_{T};Counts (Weighted)", 148, 20, 1500);
                  recoJtPt_General_ParaFills_h[BB] = new TH1D(("recoJtPt_" + nameStr + tempSysStr + "General_ParaFills_h").c_str(),
                        ";Reco. Jet p_{T};Counts (Weighted)", 148, 20, 1500);

                  genJtPt_General_h[BB] = new TH1D(("genJtPt_" + nameStr + tempSysStr + "General_h").c_str(),
                        ";Gen. Jet p_{T};Counts (Weighted)", 148, 20, 1500);

                  genJtPt_General_AllGen_h[BB] = new TH1D(("genJtPt_" + nameStr + tempSysStr + "General_AllGen_h").c_str(),
                        ";Gen. Jet p_{T};Counts (Weighted)", 148, 20, 1500);

                  genJtPt_General_ParaFills_h[BB] = new TH1D(("genJtPt_" + nameStr + tempSysStr + "General_ParaFills_h").c_str(),
                        ";Gen. Jet p_{T};Counts (Weighted)", 148, 20, 1500);

                  // std::cout << "response_General_h " << " " << nGeneralBins << std::endl;
                  response_General_h[BB] = new TH2D(("response_" + nameStr + tempSysStr + "General_h").c_str(),
                        ";Reco. Jet p_{T};Gen. Jet p_{T}", 148, 20, 1500, 148, 20, 1500);

                  response_General_Half1_h[BB] = new TH2D(("response_" + nameStr + tempSysStr + "General_Half1_h").c_str(),
                        ";Reco. Jet p_{T};Gen. Jet p_{T}", 148, 20, 1500, 148, 20, 1500);
                  response_General_Half2_h[BB] = new TH2D(("response_" + nameStr + tempSysStr + "General_Half2_h").c_str(),
                        ";Reco. Jet p_{T};Gen. Jet p_{T}", 148, 20, 1500, 148, 20, 1500);

                  response_General_ParaFills_h[BB] = new TH2D(("response_" + nameStr + tempSysStr + "General_ParaFills_h").c_str(),
                        ";Reco. Jet p_{T};Gen. Jet p_{T}", 148, 20, 1500, 148, 20, 1500);

                  responseNorm_General_h[BB] = new TH2D(("responseNorm_" + nameStr + tempSysStr + "General_h").c_str(),
                        ";Reco. Jet p_{T};Gen. Jet p_{T}", nNormBins, 0, 3, 148, 20, 1500);

                  std::vector<TH1*> tempVect = {recoJtPt_GoodGen_h[BB], recoJtPt_GoodGen_ParaFills_h[BB], recoJtPt_General_h[BB], recoJtPt_General_AllReco_h[BB], recoJtPt_General_ParaFills_h[BB], genJtPt_General_h[BB], genJtPt_General_AllGen_h[BB], genJtPt_General_ParaFills_h[BB], response_General_h[BB], response_General_Half1_h[BB], response_General_Half2_h[BB], response_General_ParaFills_h[BB], responseNorm_General_h[BB]};
                  centerTitles(tempVect);
                  setSumW2(tempVect);

                  genJtPt_GoodReco_h[BB] = new TH1D(("genJtPt_" + nameStr + tempSysStr + "GoodReco_h").c_str(),
                        ";Gen. Jet p_{T};Counts (Weighted)", nGenJtPtBins[B], genJtPtBins[B]);

                  genJtPt_GoodReco_ParaFills_h[BB] = new TH1D(("genJtPt_" + nameStr + tempSysStr + "GoodReco_ParaFills_h").c_str(),
                        ";Gen. Jet p_{T};Counts (Weighted)", nGenJtPtBins[B], genJtPtBins[B]);

                  // std::cout << "response_RecoGenSymm_h " << " " << nRecoJtPtBins[B] << " " << nGenJtPtBins[B] << std::endl;
                  response_RecoGenSymm_h[BB] = new TH2D(("response_" + nameStr + tempSysStr + "RecoGenSymm_h").c_str(),
                        ";Reco. Jet p_{T};Gen. Jet p_{T}", nRecoJtPtBins[B], recoJtPtBins[B], nGenJtPtBins[B], genJtPtBins[B]);

                  response_RecoGenAsymm_h[BB] = new TH2D(("response_" + nameStr + tempSysStr + "RecoGenAsymm_h").c_str(),
                        ";Reco. Jet p_{T};Gen. Jet p_{T}", nRecoJtPtBins[B], recoJtPtBins[B], nGenJtPtBins[B], genJtPtBins[B]);

                  if(doRooResponse){
                     rooResponse_RecoGenAsymm_h[BB] = new RooUnfoldResponse(("rooResponse_" + nameStr + tempSysStr + "RecoGenAsymm_h").c_str(), "");

                     rooResponse_RecoGenAsymm_h[BB]->Setup(recoJtPt_GoodGen_h[BB], genJtPt_GoodReco_h[BB]);
                  }

                  tempVect = {genJtPt_GoodReco_h[BB], genJtPt_GoodReco_ParaFills_h[BB], response_RecoGenSymm_h[BB], response_RecoGenAsymm_h[BB]};
                  centerTitles(tempVect);
                  setSumW2(tempVect);
               }
            }
         }
      }
   }

   const Int_t nMaxJet = 500;
   goodGlobalSelection globalSel;
   globalSel.setIsPbPb(!isPP);

   specialHYDJETEventExclude specialSel;

   Int_t nFileLoopEvt = 0;
   fileLoopWatch.start();

   std::vector<int> centPos;
   if(isPP){for(Int_t cI = 0; cI < nCentBins; ++cI){centPos.push_back(cI);}}
   else centPos.push_back(-1);

   Int_t nref_;
   Float_t jtpt_[nMaxJet];
   Float_t rawpt_[nMaxJet];
   Float_t jteta_[nMaxJet];
   Float_t jtphi_[nMaxJet];
   Float_t refpt_[nMaxJet];
   Float_t refeta_[nMaxJet];
   Float_t refphi_[nMaxJet];
   Float_t jtPfCHF_[nMaxJet];
   Float_t jtPfCEF_[nMaxJet];
   Float_t jtPfNHF_[nMaxJet];
   Float_t jtPfNEF_[nMaxJet];
   Float_t jtPfMUF_[nMaxJet];
   Float_t jtPfCHMF_[nMaxJet];
   Float_t jtPfCEMF_[nMaxJet];
   Float_t jtPfNHMF_[nMaxJet];
   Float_t jtPfNEMF_[nMaxJet];
   Float_t jtPfMUMF_[nMaxJet];
   Int_t jtPfCHM_[nMaxJet];
   Int_t jtPfCEM_[nMaxJet];
   Int_t jtPfNHM_[nMaxJet];
   Int_t jtPfNEM_[nMaxJet];
   Int_t jtPfMUM_[nMaxJet];

   Int_t ngen_;
   Float_t genpt_[nMaxJet];
   Float_t genphi_[nMaxJet];
   Float_t geneta_[nMaxJet];
   Int_t gensubid_[nMaxJet];

   for(unsigned int fI = 0; fI < fileList.size(); ++fI)
   {
      std::cout << "Processing file " << fI << "/" << fileList.size() << ": \'" << fileList[fI] << "\'" << std::endl;

      inFile_p = TFile::Open(mntToXRootdFileString(fileList[fI]).c_str(), "READ");
      TTree *jetTrees_p = nullptr;

      jetTrees_p = (TTree *)inFile_p->Get(treeName.c_str());

      jetTrees_p->SetBranchStatus("*", 0);
      jetTrees_p->SetBranchStatus("nref", 1);
      jetTrees_p->SetBranchStatus("jtpt", 1);
      jetTrees_p->SetBranchStatus("rawpt", 1);
      jetTrees_p->SetBranchStatus("jteta", 1);
      jetTrees_p->SetBranchStatus("jtphi", 1);
      jetTrees_p->SetBranchStatus("refpt", 1);
      jetTrees_p->SetBranchStatus("refphi", 1);
      jetTrees_p->SetBranchStatus("refeta", 1);
      jetTrees_p->SetBranchStatus("jtPfCHF", 1);
      jetTrees_p->SetBranchStatus("jtPfCEF", 1);
      jetTrees_p->SetBranchStatus("jtPfNHF", 1);
      jetTrees_p->SetBranchStatus("jtPfNEF", 1);
      jetTrees_p->SetBranchStatus("jtPfMUF", 1);
      jetTrees_p->SetBranchStatus("jtPfCHMF", 1);
      jetTrees_p->SetBranchStatus("jtPfCEMF", 1);
      jetTrees_p->SetBranchStatus("jtPfNHMF", 1);
      jetTrees_p->SetBranchStatus("jtPfNEMF", 1);
      jetTrees_p->SetBranchStatus("jtPfMUMF", 1);
      jetTrees_p->SetBranchStatus("jtPfCHM", 1);
      jetTrees_p->SetBranchStatus("jtPfCEM", 1);
      jetTrees_p->SetBranchStatus("jtPfNHM", 1);
      jetTrees_p->SetBranchStatus("jtPfNEM", 1);
      jetTrees_p->SetBranchStatus("jtPfMUM", 1);

      jetTrees_p->SetBranchStatus("ngen", 1);
      jetTrees_p->SetBranchStatus("genpt", 1);
      jetTrees_p->SetBranchStatus("geneta", 1);
      jetTrees_p->SetBranchStatus("genphi", 1);
      jetTrees_p->SetBranchStatus("gensubid", 1);

      jetTrees_p->SetBranchAddress("nref", &(nref_));
      jetTrees_p->SetBranchAddress("jtpt", jtpt_);
      jetTrees_p->SetBranchAddress("rawpt", rawpt_);
      jetTrees_p->SetBranchAddress("jteta", jteta_);
      jetTrees_p->SetBranchAddress("jtphi", jtphi_);
      jetTrees_p->SetBranchAddress("refpt", refpt_);
      jetTrees_p->SetBranchAddress("refphi", refphi_);
      jetTrees_p->SetBranchAddress("refeta", refeta_);
      jetTrees_p->SetBranchAddress("jtPfCHF", jtPfCHF_);
      jetTrees_p->SetBranchAddress("jtPfCEF", jtPfCEF_);
      jetTrees_p->SetBranchAddress("jtPfNHF", jtPfNHF_);
      jetTrees_p->SetBranchAddress("jtPfNEF", jtPfNEF_);
      jetTrees_p->SetBranchAddress("jtPfMUF", jtPfMUF_);
      jetTrees_p->SetBranchAddress("jtPfCHMF", jtPfCHMF_);
      jetTrees_p->SetBranchAddress("jtPfCEMF", jtPfCEMF_);
      jetTrees_p->SetBranchAddress("jtPfNHMF", jtPfNHMF_);
      jetTrees_p->SetBranchAddress("jtPfNEMF", jtPfNEMF_);
      jetTrees_p->SetBranchAddress("jtPfMUMF", jtPfMUMF_);
      jetTrees_p->SetBranchAddress("jtPfCHM", jtPfCHM_);
      jetTrees_p->SetBranchAddress("jtPfCEM", jtPfCEM_);
      jetTrees_p->SetBranchAddress("jtPfNHM", jtPfNHM_);
      jetTrees_p->SetBranchAddress("jtPfNEM", jtPfNEM_);
      jetTrees_p->SetBranchAddress("jtPfMUM", jtPfMUM_);

      jetTrees_p->SetBranchAddress("ngen", &(ngen_));
      jetTrees_p->SetBranchAddress("genpt", genpt_);
      jetTrees_p->SetBranchAddress("geneta", geneta_);
      jetTrees_p->SetBranchAddress("genphi", genphi_);
      jetTrees_p->SetBranchAddress("gensubid", gensubid_);

      Float_t pthat_;

      jetTrees_p->SetBranchStatus("pthat", 1);
      jetTrees_p->SetBranchAddress("pthat", &pthat_);

      Float_t vz_;
      Float_t hiHF_;
      Int_t hiBin_;
      unsigned int run_, lumi_;
      unsigned long long evt_;

      Int_t HBHENoiseFilterResultRun2Loose_ = -1;
      Int_t pprimaryVertexFilter_ = -1;
      Int_t pBeamScrapingFilter_ = -1;
      Int_t phfCoincFilter3_ = -1;
      Int_t pclusterCompatibilityFilter_ = -1;

      TTree *hiTree_p = (TTree*)inFile_p->Get("hiEvtAnalyzer/HiTree");

      hiTree_p->SetBranchStatus("*", 0);
      hiTree_p->SetBranchStatus("hiBin", 1);
      hiTree_p->SetBranchStatus("vz", 1);
      hiTree_p->SetBranchStatus("hiHF", 1);
      hiTree_p->SetBranchStatus("run", 1);
      hiTree_p->SetBranchStatus("lumi", 1);
      hiTree_p->SetBranchStatus("evt", 1);

      hiTree_p->SetBranchAddress("hiBin", &hiBin_);
      hiTree_p->SetBranchAddress("vz", &vz_);
      hiTree_p->SetBranchAddress("hiHF", &hiHF_);
      hiTree_p->SetBranchAddress("run", &run_);
      hiTree_p->SetBranchAddress("lumi", &lumi_);
      hiTree_p->SetBranchAddress("evt", &evt_);

      TTree *skimTree_p = (TTree*)inFile_p->Get("skimanalysis/HltTree");

      skimTree_p->SetBranchStatus("*", 0);
      skimTree_p->SetBranchStatus("HBHENoiseFilterResultRun2Loose", 1);
      skimTree_p->SetBranchAddress("HBHENoiseFilterResultRun2Loose", &HBHENoiseFilterResultRun2Loose_);

      if(!isPP){
         skimTree_p->SetBranchStatus("pprimaryVertexFilter", 1);
         skimTree_p->SetBranchStatus("phfCoincFilter3", 1);
         skimTree_p->SetBranchStatus("pclusterCompatibilityFilter", 1);

         skimTree_p->SetBranchAddress("pprimaryVertexFilter", &pprimaryVertexFilter_);
         skimTree_p->SetBranchAddress("phfCoincFilter3", &phfCoincFilter3_);
         skimTree_p->SetBranchAddress("pclusterCompatibilityFilter", &pclusterCompatibilityFilter_);
      }
      else{
         skimTree_p->SetBranchStatus("pBeamScrapingFilter", 1);
         skimTree_p->SetBranchStatus("pPAprimaryVertexFilter", 1);

         skimTree_p->SetBranchAddress("pBeamScrapingFilter", &pBeamScrapingFilter_);
         skimTree_p->SetBranchAddress("pPAprimaryVertexFilter", &pprimaryVertexFilter_);
      }

      const Int_t nEntries = TMath::Min((Int_t)1000000000, (Int_t)jetTrees_p->GetEntries());
      const Double_t entryFrac = inEntryFrac;
      const Int_t nEntriesToProcess = entryFrac*nEntries;
      const Int_t printInterval = TMath::Max(1, nEntriesToProcess/20);

      nFileLoopEvt += nEntriesToProcess;

      cppWatch subEntryWatch, nonJet, jetLoop1, jetLoop2, jetLoop3, jetLoop3Sub1, jetLoop3Sub2, jetLoop3Sub3, jetLoop3Sub3Suba, jetLoop3Sub3Subb, jetLoop3Sub3Subc, jetLoop3Sub3Sub0;

      subEntryWatch.clear();
      nonJet.clear();
      jetLoop1.clear();
      jetLoop2.clear();
      jetLoop3.clear();

      subEntryWatch.start();
      std::cout << "Processing " << nEntriesToProcess << "... (" << prettyString(entryFrac, 2, false) << " fraction of " << nEntries << ")" << std::endl;
      for(Int_t entry = 0; entry < nEntriesToProcess; ++entry)
      {
         if(nEntriesToProcess >= 50000 && entry % printInterval == 0)
         {
            std::cout << " Entry: " << entry << "/" << nEntriesToProcess << std::endl;
            subEntryWatch.stop();
            std::cout << "  Timing: " << subEntryWatch.totalWall() << std::endl;
            subEntryWatch.clear();
            subEntryWatch.start();

            nonJet.stop();

            std::cout << "   nonjet: " << nonJet.totalWall() << std::endl;	
            std::cout << "   jetloop1: " << jetLoop1.totalWall() << std::endl;
            std::cout << "   jetloop2: " << jetLoop2.totalWall() << std::endl;
            std::cout << "   jetloop3: " << jetLoop3.totalWall() << std::endl;
            std::cout << "    jetloop3Sub1: " << jetLoop3Sub1.totalWall() << std::endl;
            std::cout << "    jetloop3Sub2: " << jetLoop3Sub2.totalWall() << std::endl;
            std::cout << "    jetloop3Sub3: " << jetLoop3Sub3.totalWall() << std::endl;
            std::cout << "     jetloop3Sub3Sub0: " << jetLoop3Sub3Sub0.totalWall() << std::endl;
            std::cout << "     jetloop3Sub3Suba: " << jetLoop3Sub3Suba.totalWall() << std::endl;
            std::cout << "     jetloop3Sub3Subb: " << jetLoop3Sub3Subb.totalWall() << std::endl;
            std::cout << "     jetloop3Sub3Subc: " << jetLoop3Sub3Subc.totalWall() << std::endl;

            nonJet.clear();
            nonJet.start();

            jetLoop1.clear();
            jetLoop2.clear();
            jetLoop3.clear();
            jetLoop3Sub1.clear();
            jetLoop3Sub2.clear();
            jetLoop3Sub3.clear();
            jetLoop3Sub3Sub0.clear();
            jetLoop3Sub3Suba.clear();
            jetLoop3Sub3Subb.clear();
            jetLoop3Sub3Subc.clear();
         }

         hiTree_p->GetEntry(entry);
         skimTree_p->GetEntry(entry);

         globalSel.setVz(vz_);
         globalSel.setHiHF(hiHF_);
         globalSel.setPprimaryVertexFilter(pprimaryVertexFilter_);
         globalSel.setPBeamScrapingFilter(pBeamScrapingFilter_);
         globalSel.setPhfCoincFilter3(phfCoincFilter3_);
         globalSel.setHBHENoiseFilterResultRun2Loose(HBHENoiseFilterResultRun2Loose_);
         globalSel.setPclusterCompatibilityFilter(pclusterCompatibilityFilter_);

         if(!globalSel.isGood()) continue;

         jetTrees_p->GetEntry(entry);

         Double_t ncollWeight_ = 1.;      
         if(!isPP)
         {
            centPos[0] = -1;
            for(Int_t cI = 0; cI < nCentBins; ++cI)
            {
               if(centBinsLow[cI] * 2 <= hiBin_ && hiBin_ < centBinsHi[cI]*2)
               {
                  centPos[0] = cI;
                  break;
               }
            }
            if(centPos[0] < 0)
               continue;

            bool badJetSpecialSel = specialSel.CheckEventBadJet(ngen_, genpt_, genphi_, geneta_, gensubid_, entry);
            if(badJetSpecialSel)
               continue;

            ncollWeight_ = findNcoll_Renorm(hiBin_);
         }

         Double_t pthatWeight_ = -1;
         for(unsigned int pI = 0; pI < pthats.size() - 1; ++pI){
            if(pthats[pI] <= pthat_ && pthat_ < pthats[pI+1]){
               pthatWeight_ = pthatWeights[pI];
               break;
            }
         }
         if(pthat_ > pthats[pthats.size()-1]) pthatWeight_ = pthatWeights[pthatWeights.size()-1];

         if(pthatWeight_ < 0){
            std::cout << "WARNING - NO WEIGHT FOR pthat \'" << pthat_ << "\'. Set to 1" << std::endl;
            pthatWeight_ = 1.;
         }

         Double_t fullWeight_ = ncollWeight_*pthatWeight_;
         pthat_h->Fill(pthat_);
         pthatWeighted_h->Fill(pthat_, pthatWeight_);
         pthatFullWeighted_h->Fill(pthat_, fullWeight_);

         if(!isPP){
            centrality_h->Fill(hiBin_/2.);
            centralityWeighted_h->Fill(hiBin_/2., ncollWeight_);
            centralityFullWeighted_h->Fill(hiBin_/2., fullWeight_);
         }

         Int_t scaleHiBin = 160;
         if(!isPP) scaleHiBin = hiBin_/2;

         Bool_t isPara = randGen_p->Uniform(0., 1.) < fracParaFills;

         std::string algoName = treeName;
         //	if(algoName.find("/") != std::string::npos) algoName = algoName.substr(0, algoName.find("/"));

         const Int_t nUsed = ngen_;
         bool isUsed[nUsed];
         for(Int_t gI = 0; gI < nUsed; ++gI)
            isUsed[gI] = false;

         nonJet.stop();

         //WE HAVE TO CLEAN BS REFPT MATCHES
         for(Int_t jI = 0; jI < nref_; ++jI)
         {
            if(jtpt_[jI] < 20) continue;
            if(TMath::Abs(jteta_[jI]) > 3.) continue;
            if(refpt_[jI] < 0) continue;
            if(refpt_[jI] / jtpt_[jI] > 0.25) continue;

            Double_t deltaPTInit = refpt_[jI] / jtpt_[jI];
            Double_t deltaPhiInit = (jtphi_[jI] - refphi_[jI]);
            Double_t deltaEtaInit = (jteta_[jI] - refeta_[jI]);
            Double_t deltaRInit = getDR(jteta_[jI], jtphi_[jI], refeta_[jI], refphi_[jI]);
            bool isReplace = false;

            for(Int_t gI = 0; gI < ngen_; ++gI)
            {
               //	    if(gensubid_[gI] == 0) continue;
               Double_t deltaR = getDR(jteta_[jI], jtphi_[jI], geneta_[gI], genphi_[gI]);
               Double_t deltaPT = genpt_[gI] / jtpt_[jI];
               Double_t deltaPhi = (jtphi_[jI] - genphi_[gI]);   // WTF
               Double_t deltaEta = (jteta_[jI] - geneta_[gI]);

               if(deltaR < 0.15 && genpt_[gI] / jtpt_[jI] > 0.5)
               {
                  if(centPos.size() != 0)
                  {
                     int B = GetBin(centPos[0]);
                     deltaPtOrig_h[B]->Fill(deltaPTInit);
                     deltaPhiOrig_h[B]->Fill(deltaPhiInit);
                     deltaEtaOrig_h[B]->Fill(deltaEtaInit);
                     deltaROrig_h[B]->Fill(deltaRInit);

                     deltaPtReplace_h[B]->Fill(deltaPT);
                     deltaPhiReplace_h[B]->Fill(deltaPhi);
                     deltaEtaReplace_h[B]->Fill(deltaEta);
                     deltaRReplace_h[B]->Fill(deltaR);
                  }

                  refpt_[jI] = -999;
                  refphi_[jI] = -999;
                  refeta_[jI] = -999;

                  isReplace = true;
                  break;
               }
            }

            if(!isReplace && jtpt_[jI] > 100 && TMath::Abs(jteta_[jI]) < 2.)
            {
               if(centPos.size() != 0)
               {
                  int B = GetBin(centPos[0]);
                  deltaPtNotReplace_h[B]->Fill(deltaPTInit);
                  deltaPhiNotReplace_h[B]->Fill(deltaPhiInit);
                  deltaEtaNotReplace_h[B]->Fill(deltaEtaInit);
                  deltaRNotReplace_h[B]->Fill(deltaRInit);
               }
            }
         }

         for(Int_t jI = 0; jI < nref_; ++jI)
         {
            jetLoop1.start();

            if(TMath::Abs(jteta_[jI]) > jtAbsEtaMax)
               continue;

            Double_t refptTemp = refpt_[jI];
            if(refptTemp < 0)
            {
               for(Int_t gI = 0; gI < ngen_; ++gI)
               {
                  if(gensubid_[gI] == 0) continue;
                  else if(isUsed[gI]) continue;
                  else if(getDR(jteta_[jI], jtphi_[jI], geneta_[gI], genphi_[gI]) < 0.2 + TMath::Min(0.0, 0.3 - rValD)){
                     isUsed[gI] = true;
                     refptTemp = genpt_[gI];
                     break;
                  }
               }
            }

            jetLoop1.stop();

            if(refptTemp < minGenJtPt) continue;

            std::vector<bool> passesID;
            for(Int_t iI = 0; iI < nID; ++iI)
            {
               bool pass = jtPfCHMFCutLow[iI] <= jtPfCHMF_[jI] && jtPfCHMF_[jI] <= jtPfCHMFCutHi[iI];
               pass = pass && jtPfMUMFCutLow[iI] <= jtPfMUMF_[jI] && jtPfMUMF_[jI] <= jtPfMUMFCutHi[iI];
               pass = pass && jtPfNHFCutLow[iI] <= jtPfNHF_[jI] && jtPfNHF_[jI] <= jtPfNHFCutHi[iI];
               pass = pass && jtPfNEFCutLow[iI] <= jtPfNEF_[jI] && jtPfNEF_[jI] <= jtPfNEFCutHi[iI];
               pass = pass && jtPfMUFCutLow[iI] <= jtPfMUF_[jI] && jtPfMUF_[jI] <= jtPfMUFCutHi[iI];
               pass = pass && jtPfCHFCutLow[iI] <= jtPfCHF_[jI] && jtPfCHF_[jI] <= jtPfCHFCutHi[iI];
               pass = pass && jtPfCEFCutLow[iI] <= jtPfCEF_[jI] && jtPfCEF_[jI] <= jtPfCEFCutHi[iI];
               pass = pass && jtPfCEM_[jI] + jtPfNEM_[jI] + jtPfCHM_[jI] + jtPfNHM_[jI] + jtPfMUM_[jI] >= jtPfMinMult[iI];
               pass = pass && jtPfCHM_[jI] >= jtPfMinChgMult[iI];

               passesID.push_back(pass);
            }

            std::vector<int> jtAbsEtaPoses;
            for(Int_t aI = 0; aI < nJtAbsEtaBins; ++aI){
               if(TMath::Abs(jteta_[jI]) >= jtAbsEtaBinsLow[aI] && TMath::Abs(jteta_[jI]) < jtAbsEtaBinsHi[aI]){
                  jtAbsEtaPoses.push_back(aI);
               }
            }

            for(Int_t mI = 0; mI < nResponseMod; ++mI){
               jetLoop2.start();

               Double_t jtPtFillVal[nSyst];
               double fullWeight2[nSyst];

               bool oneRecoIsGood = false;

               for(Int_t sI = 0; sI < nSystReduced; ++sI){	  
                  jtPtFillVal[sI] = jtpt_[jI] + (jtpt_[jI] - refpt_[jI]) * responseMod[mI];
                  fullWeight2[sI] = fullWeight_;

                  if((systType[sI] == JECUpMC)) jtPtFillVal[sI] += jtPtFillVal[sI]*jecVarMC;
                  else if((systType[sI] == JECDownMC)) jtPtFillVal[sI] -= jtPtFillVal[sI]*jecVarMC;
                  else if((systType[sI] == JECUpData)) jtPtFillVal[sI] += jtPtFillVal[sI]*jecVarData;
                  else if((systType[sI] == JECDownData)) jtPtFillVal[sI] -= jtPtFillVal[sI]*jecVarData;
                  else if((systType[sI] == JECUpUE) || (systType[sI] == JECDownUE) )
                  {
                     Double_t tempScale =  jtPtFillVal[sI]/rawpt_[jI];
                     Double_t tempRawPt = rawpt_[jI];
                     if((systType[sI] == JECUpUE)) tempRawPt += TMath::Abs(scaleErr.getMuDataMinusMC(scaleHiBin, jteta_[jI], rValI, FlowDefaultInRho));
                     else tempRawPt -= TMath::Abs(scaleErr.getMuDataMinusMC(scaleHiBin, jteta_[jI], rValI, FlowDefaultInRho));

                     jtPtFillVal[sI] = tempRawPt*tempScale;
                  }
                  else if((systType[sI] == JERMC)) jtPtFillVal[sI] += (jtPtFillVal[sI] - refpt_[jI])*jerVarMC;
                  else if((systType[sI] == JERData)) jtPtFillVal[sI] += (jtPtFillVal[sI] - refpt_[jI])*jerVarData[mI];
                  else if((systType[sI] == PriorUp1PowerPthat)) fullWeight2[sI] *= pthat_/pthatLow;
                  else if((systType[sI] == PriorDown1PowerPthat)) fullWeight2[sI] *= pthatLow/pthat_;
                  else if((systType[sI] == PriorFlat)) fullWeight2[sI] *= 1./pthatWeight_;
                  /*
                     else if((systType[sI] == PriorFlat) && refpt_[jI] > 0) fullWeight2[sI] *= flatWeight.getJtWeight(algoName == hiBin_/2 == refpt_[jI] == jteta_[jI]);
                     else if((systType[sI] == PriorNoWeight) && refpt_[jI] > 0) fullWeight2[sI] /= pthatWeight_;
                     */
                  if(jtPtFillVal[sI] > minRecoJtPt) oneRecoIsGood = true;
               }

               jetLoop2.stop();

               if(!oneRecoIsGood) continue;

               jetLoop3.start();

               for(auto const & cent : centPos){
                  jetLoop3Sub1.start();

                  int B = GetBin(cent);
                  bool goodTruth = (refpt_[jI] >= genJtPtBins[B][0] && refpt_[jI] < genJtPtBins[B][nGenJtPtBins[B]] && refpt_[jI] > minJtPtCut);

                  bool fillLowTruth = false;
                  bool fillHighTruth = false;
                  if(!goodTruth && refpt_[jI] > 0.0){
                     if(refpt_[jI] < genJtPtBins[B][0]) fillLowTruth = true;
                     else if(refpt_[jI] >= genJtPtBins[B][nGenJtPtBins[B]]) fillHighTruth = true;
                  }

                  bool goodTruthGeneral = (refpt_[jI] >= generalBins[0] && refpt_[jI] < generalBins[nGeneralBins] && refpt_[jI] > minJtPtCut);  

                  std::vector<Int_t> goodRecoPos;
                  std::vector<Int_t> goodRecoTruncPos;
                  std::vector<Int_t> goodRecoGeneralPos;

                  for(Int_t sI = 0; sI < nSystReduced; ++sI)
                  {
                     if(jtPtFillVal[sI] >= genJtPtBins[B][0] && jtPtFillVal[sI] < genJtPtBins[B][nGenJtPtBins[B]])
                        goodRecoPos.push_back(sI);
                     if(jtPtFillVal[sI] >= recoJtPtBins[B][0] && jtPtFillVal[sI] < recoJtPtBins[B][nRecoJtPtBins[B]])
                        goodRecoTruncPos.push_back(sI);
                     if(jtPtFillVal[sI] >= generalBins[0] && jtPtFillVal[sI] < generalBins[nGeneralBins])
                        goodRecoGeneralPos.push_back(sI);
                  }

                  jetLoop3Sub1.stop();
                  jetLoop3Sub2.start();

                  if(refpt_[jI] < 0)
                  {
                     int B = GetBin(cent);
                     bool goodTruth2 = true;
                     if((refptTemp >= genJtPtBins[B][0]) == false)
                        goodTruth2 = false;
                     if((refptTemp < genJtPtBins[B][nGenJtPtBins[B]]) == false)
                        goodTruth2 = false;
                     if((refptTemp > minJtPtCut) == false)
                        goodTruth2 = false;

                     if(!goodTruth2)
                     {
                        for(unsigned int aI = 0; aI < jtAbsEtaPoses.size(); ++aI)
                        {
                           for(unsigned int iI = 0; iI < passesID.size(); ++iI)
                           {
                              if(!passesID[iI]) continue;

                              B = GetBin(cent, iI, mI, jtAbsEtaPoses[aI]);
                              if(isPara)
                                 recoJtPt_NoTruth_ParaFills_h[B]->Fill(jtPtFillVal[0], fullWeight_);
                              else
                                 recoJtPt_NoTruth_h[B]->Fill(jtPtFillVal[0], fullWeight_);
                           }
                        }
                     }	    
                  }

                  jetLoop3Sub2.stop();
                  jetLoop3Sub3.start();	      

                  jetLoop3Sub3Sub0.start();	      
                  if(!isPara)
                  {
                     for(unsigned int iI = 0; iI < passesID.size(); ++iI)
                     {
                        if(!passesID[iI]) continue;
                        for(unsigned int aI = 0; aI < jtAbsEtaPoses.size(); ++aI)
                        {
                           for(unsigned int sI = 0; sI < goodRecoGeneralPos.size(); ++sI)
                           {
                              int B = GetBin(cent, iI, mI, jtAbsEtaPoses[aI], goodRecoGeneralPos[sI]);
                              double W = fullWeight2[goodRecoGeneralPos[sI]];
                              recoJtPt_General_AllReco_h[B]->Fill(refpt_[jI], W);
                           }
                        }
                     }
                  }

                  if(goodTruthGeneral)
                  {
                     for(unsigned int iI = 0; iI < passesID.size(); ++iI)
                     {
                        if(!passesID[iI]) continue;
                        for(unsigned int aI = 0; aI < jtAbsEtaPoses.size(); ++aI)
                        {
                           if(!isPara)
                           {
                              for(Int_t sI = 0; sI < nSystReduced; ++sI)
                              {
                                 int B = GetBin(cent, iI, mI, jtAbsEtaPoses[aI], sI);
                                 genJtPt_General_AllGen_h[B]->Fill(refpt_[jI], fullWeight2[sI]);
                              }
                           }

                           for(unsigned int sI = 0; sI < goodRecoGeneralPos.size(); ++sI)
                           {
                              int B = GetBin(cent, iI, mI, jtAbsEtaPoses[aI], goodRecoGeneralPos[sI]);
                              double W = fullWeight2[goodRecoGeneralPos[sI]];
                              // int BRefPT = FindBin(148, 20, 1500, refpt_[jI]);
                              // int BJTPT = FindBin(148, 20, 1500, jtPtFillVal[goodRecoGeneralPos[sI]]);
                              if(/*isPara*/ false)
                              {
                                 genJtPt_General_ParaFills_h[B]->Fill(refpt_[jI], W);
                                 response_General_ParaFills_h[B]->Fill(jtPtFillVal[goodRecoGeneralPos[sI]], refpt_[jI], W);
                                 recoJtPt_General_ParaFills_h[B]->Fill(jtPtFillVal[goodRecoGeneralPos[sI]], W);
                              }
                              else{	       
                                 genJtPt_General_h[B]->Fill(refpt_[jI], W);
                                 response_General_h[B]->Fill(jtPtFillVal[goodRecoGeneralPos[sI]], refpt_[jI], W);

                                 if(randGen_p->Uniform(0.0, 1.0) < 0.5)
                                    response_General_Half1_h[B]->Fill(jtPtFillVal[goodRecoGeneralPos[sI]], refpt_[jI], W);
                                 else
                                    response_General_Half2_h[B]->Fill(jtPtFillVal[goodRecoGeneralPos[sI]], refpt_[jI], W);

                                 double Norm = jtPtFillVal[goodRecoGeneralPos[sI]] / refpt_[jI];
                                 responseNorm_General_h[B]->Fill(Norm, refpt_[jI], W);
                                 recoJtPt_General_h[B]->Fill(jtPtFillVal[goodRecoGeneralPos[sI]], W);			
                              }
                           }
                        }
                     }

                  }

                  jetLoop3Sub3Sub0.stop();	      

                  if(goodTruth || fillLowTruth || fillHighTruth)
                  {
                     Double_t tempRefPt = refpt_[jI];

                     int B = GetBin(cent);
                     if(fillLowTruth)       tempRefPt = genJtPtBins[B][0] + 1;
                     else if(fillHighTruth) tempRefPt = genJtPtBins[B][nGenJtPtBins[B]] - 1;

                     for(unsigned int iI = 0; iI < passesID.size(); ++iI)
                     {
                        if(!passesID[iI]) continue;
                        for(unsigned int aI = 0; aI < jtAbsEtaPoses.size(); ++aI)
                        {
                           if(!isPara)
                           {
                              jetLoop3Sub3Suba.start();	      

                              int B = GetBin(cent, iI, mI, jtAbsEtaPoses[aI]);
                              genJtPt_All_h[B]->Fill(tempRefPt, fullWeight_);

                              for(unsigned int sI = 0; sI < goodRecoPos.size(); ++sI)
                              {
                                 B = GetBin(cent, iI, mI, jtAbsEtaPoses[aI], goodRecoPos[sI]);
                                 response_RecoGenSymm_h[B]->Fill(jtPtFillVal[goodRecoPos[sI]], tempRefPt, fullWeight2[goodRecoPos[sI]]);
                              }  

                              jetLoop3Sub3Suba.stop();	      		      
                              jetLoop3Sub3Subb.start();	      

                              for(unsigned int sI = 0; sI < goodRecoTruncPos.size(); ++sI)
                              {
                                 B = GetBin(cent, iI, mI, jtAbsEtaPoses[aI], goodRecoTruncPos[sI]);
                                 genJtPt_GoodReco_h[B]->Fill(tempRefPt, fullWeight2[goodRecoTruncPos[sI]]);
                                 response_RecoGenAsymm_h[B]->Fill(jtPtFillVal[goodRecoTruncPos[sI]],
                                       tempRefPt, fullWeight2[goodRecoTruncPos[sI]]);
                                 if(doRooResponse)
                                    rooResponse_RecoGenAsymm_h[B]->Fill(jtPtFillVal[goodRecoTruncPos[sI]],
                                          tempRefPt, fullWeight2[goodRecoTruncPos[sI]]);

                                 recoJtPt_GoodGen_h[B]->Fill(jtPtFillVal[goodRecoTruncPos[sI]],
                                       fullWeight2[goodRecoTruncPos[sI]]);
                              }

                              jetLoop3Sub3Subb.stop();	      
                              jetLoop3Sub3Subc.start();	      

                              if(priorFlatPos >= 0)
                              {
                                 B = GetBin(cent, iI, mI, aI);
                                 genJtPt_CheckPriorFlat_h[B]->Fill(tempRefPt, fullWeight2[priorFlatPos]);
                              }

                              if(goodRecoPos.size() != 0)
                              {
                                 B = GetBin(cent, iI, mI, jtAbsEtaPoses[aI]);
                                 if(goodRecoPos[0] == 0)
                                    recoJtPt_h[B]->Fill(jtPtFillVal[0], fullWeight_);
                              }

                              if(goodRecoTruncPos.size() != 0)
                              {
                                 B = GetBin(cent, iI, mI, jtAbsEtaPoses[aI]);
                                 if(goodRecoTruncPos[0] == 0)
                                    recoJtPt_RecoTrunc_h[B]->Fill(jtPtFillVal[0], fullWeight_);			
                              }

                              B = GetBin(cent, iI, mI, jtAbsEtaPoses[aI]);
                              genJtPt_h[B]->Fill(tempRefPt, fullWeight_);

                              jetLoop3Sub3Subc.stop();	      
                           }		  
                           else
                           {
                              for(unsigned int sI = 0; sI < goodRecoTruncPos.size(); ++sI)
                              {
                                 int B = GetBin(cent, iI, mI, jtAbsEtaPoses[aI], goodRecoTruncPos[sI]);
                                 double W = fullWeight2[goodRecoTruncPos[sI]];
                                 genJtPt_GoodReco_ParaFills_h[B]->Fill(tempRefPt, W);
                                 recoJtPt_GoodGen_ParaFills_h[B]->Fill(jtPtFillVal[goodRecoTruncPos[sI]], W);
                              }		  

                              if(goodRecoPos.size() != 0)
                              {
                                 int B = GetBin(cent, iI, mI, jtAbsEtaPoses[aI]);
                                 if(goodRecoPos[0] == 0)
                                    recoJtPt_ParaFills_h[B]->Fill(jtPtFillVal[0], fullWeight_);
                              }

                              if(goodRecoTruncPos.size() != 0)
                              {
                                 int B = GetBin(cent, iI, mI, jtAbsEtaPoses[aI]);
                                 if(goodRecoTruncPos[0] == 0)
                                    recoJtPt_RecoTrunc_ParaFills_h[B]->Fill(jtPtFillVal[0], fullWeight_);
                              }

                              int B = GetBin(cent, iI, mI, jtAbsEtaPoses[aI]);
                              genJtPt_ParaFills_h[B]->Fill(tempRefPt, fullWeight_);
                           }
                        }		
                     }	     	      
                  }
                  jetLoop3Sub3.stop();
               }	  
               jetLoop3.stop();
            }	
         }

         nonJet.start();
      }

      inFile_p->Close();
      delete inFile_p;
      inFile_p = NULL;
   }

   fileLoopWatch.stop();
   writeLoopWatch.start();

   outFile_p->cd();
   generalDir_p->cd();
   pthat_h->Write("", TObject::kOverwrite);
   pthatWeighted_h->Write("", TObject::kOverwrite);
   pthatFullWeighted_h->Write("", TObject::kOverwrite);
   pthatFullRatio_h->Divide(pthatFullWeighted_h, pthatWeighted_h);
   pthatFullRatio_h->Write("", TObject::kOverwrite);

   if(!isPP)
   {
      centrality_h->Write("", TObject::kOverwrite);
      centralityWeighted_h->Write("", TObject::kOverwrite);
      centralityFullWeighted_h->Write("", TObject::kOverwrite);
      centralityFullRatio_h->Divide(centralityFullWeighted_h, centralityWeighted_h);
      centralityFullRatio_h->Write("", TObject::kOverwrite);
   }

   if(doLocalDebug || doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

   outFile_p->cd();
   dir_p->cd();

   for(Int_t cI = 0; cI < nCentBins; ++cI)
   {
      int B = GetBin(cI);

      deltaPtOrig_h[B]->Write("", TObject::kOverwrite);
      deltaPhiOrig_h[B]->Write("", TObject::kOverwrite);
      deltaEtaOrig_h[B]->Write("", TObject::kOverwrite);
      deltaROrig_h[B]->Write("", TObject::kOverwrite);

      deltaPtReplace_h[B]->Write("", TObject::kOverwrite);
      deltaPhiReplace_h[B]->Write("", TObject::kOverwrite);
      deltaEtaReplace_h[B]->Write("", TObject::kOverwrite);
      deltaRReplace_h[B]->Write("", TObject::kOverwrite);

      deltaPtNotReplace_h[B]->Write("", TObject::kOverwrite);
      deltaPhiNotReplace_h[B]->Write("", TObject::kOverwrite);
      deltaEtaNotReplace_h[B]->Write("", TObject::kOverwrite);
      deltaRNotReplace_h[B]->Write("", TObject::kOverwrite);

      for(Int_t iI = 0; iI < nID; ++iI)
      {
         for(Int_t mI = 0; mI < nResponseMod; ++mI)
         {
            for(Int_t aI = 0; aI < nJtAbsEtaBins; ++aI)
            {
               B = GetBin(cI, iI, mI, aI);

               recoJtPt_h[B]->Write("", TObject::kOverwrite);
               recoJtPt_RecoTrunc_h[B]->Write("", TObject::kOverwrite);
               recoJtPt_NoTruth_h[B]->Write("", TObject::kOverwrite);
               recoJtPt_ParaFills_h[B]->Write("", TObject::kOverwrite);
               recoJtPt_RecoTrunc_ParaFills_h[B]->Write("", TObject::kOverwrite);
               recoJtPt_NoTruth_ParaFills_h[B]->Write("", TObject::kOverwrite);

               genJtPt_h[B]->Write("", TObject::kOverwrite);
               genJtPt_ParaFills_h[B]->Write("", TObject::kOverwrite);
               genJtPt_All_h[B]->Write("", TObject::kOverwrite);

               genJtPt_CheckPriorFlat_h[B]->Write("", TObject::kOverwrite);

               for(Int_t sI = 0; sI < nSystReduced; ++sI)
               {
                  B = GetBin(cI, iI, mI, aI, sI);

                  recoJtPt_GoodGen_h[B]->Write("", TObject::kOverwrite);
                  recoJtPt_GoodGen_ParaFills_h[B]->Write("", TObject::kOverwrite);

                  recoJtPt_General_h[B]->Write("", TObject::kOverwrite);
                  recoJtPt_General_AllReco_h[B]->Write("", TObject::kOverwrite);

                  genJtPt_General_h[B]->Write("", TObject::kOverwrite);
                  genJtPt_General_AllGen_h[B]->Write("", TObject::kOverwrite);
                  response_General_h[B]->Write("", TObject::kOverwrite);
                  response_General_Half1_h[B]->Write("", TObject::kOverwrite);
                  response_General_Half2_h[B]->Write("", TObject::kOverwrite);
                  responseNorm_General_h[B]->Write("", TObject::kOverwrite);

                  recoJtPt_General_ParaFills_h[B]->Write("", TObject::kOverwrite);
                  genJtPt_General_ParaFills_h[B]->Write("", TObject::kOverwrite);
                  response_General_ParaFills_h[B]->Write("", TObject::kOverwrite);

                  genJtPt_GoodReco_h[B]->Write("", TObject::kOverwrite);
                  genJtPt_GoodReco_ParaFills_h[B]->Write("", TObject::kOverwrite);
                  response_RecoGenSymm_h[B]->Write("", TObject::kOverwrite);
                  response_RecoGenAsymm_h[B]->Write("", TObject::kOverwrite);
                  if(doRooResponse)
                     rooResponse_RecoGenAsymm_h[B]->Write("", TObject::kOverwrite);
               }
            }
         }
      }
   }

   writeLoopWatch.stop();

   deleteLoopWatch.start();

   if(doLocalDebug || doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

   outFile_p->cd();
   generalDir_p->cd();

   delete pthat_h;
   delete pthatWeighted_h;
   delete pthatFullWeighted_h;
   delete pthatFullRatio_h;

   if(!isPP){
      delete centrality_h;
      delete centralityWeighted_h;
      delete centralityFullWeighted_h;
      delete centralityFullRatio_h;
   }

   if(doLocalDebug || doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;


   deleteLoopWatch.stop();

   outFile_p->cd();

   TDirectory *cutDir_p = (TDirectory*)outFile_p->mkdir("cutDir");
   TDirectory *subDir_p = (TDirectory*)cutDir_p->mkdir("subDir");

   if(doLocalDebug || doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

   cutProp.SetNJtAlgos(1);
   cutProp.SetJtAlgos({treeName});
   cutProp.SetMinJtPtCut({minJtPtCut});
   cutProp.SetMultiJtPtCut({multiJtPtCut});
   cutProp.SetRecoTruncPos({recoTruncPos});

   if(doLocalDebug || doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;
   if(!cutProp.WriteAllVarToFile(outFile_p, cutDir_p, subDir_p)) std::cout << "Warning: Cut writing has failed" << std::endl;

   if(doLocalDebug || doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

   outFile_p->Close();
   delete outFile_p;

   if(doLocalDebug || doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

   delete randGen_p;

   totalRunWatch.stop();
   const double fileLoopWatchTotal = fileLoopWatch.totalWall();
   const double totalRunWatchTotal = totalRunWatch.totalWall();
   const double writeLoopWatchTotal = writeLoopWatch.totalWall();
   const double deleteLoopWatchTotal = deleteLoopWatch.totalWall();
   const double nonFileLoopTotal = totalRunWatch.totalWall() - fileLoopWatch.totalWall(); 
   const double nJetTrees = 1;

   std::cout << "File loop watch: " << fileLoopWatch.totalWall() << std::endl;
   std::cout << " Per event: " << fileLoopWatchTotal/nFileLoopEvt << std::endl;
   std::cout << " Per jetTree: " << fileLoopWatchTotal/nJetTrees << std::endl;
   std::cout << " Per event x jetTree: " << fileLoopWatchTotal/(nJetTrees*nFileLoopEvt) << std::endl;

   std::cout << "Total run watch: " << totalRunWatch.totalWall() << std::endl;
   std::cout << " File loop fraction: " << fileLoopWatchTotal/totalRunWatchTotal << std::endl;
   std::cout << " Non-file loop num: " << totalRunWatch.totalWall() - fileLoopWatch.totalWall() << std::endl;
   std::cout << " Write loop num: " << writeLoopWatch.totalWall() << std::endl;
   std::cout << "  Fraction of non-file loop: " << writeLoopWatchTotal/nonFileLoopTotal << std::endl;
   std::cout << " Delete loop num: " << deleteLoopWatch.totalWall() << std::endl;
   std::cout << "  Fraction of non-file loop: " << deleteLoopWatchTotal/nonFileLoopTotal << std::endl;

   return 0;
}

int main(int argc, char *argv[])
{
   CommandLine CL(argc, argv);

   if(CL.Get("h", "X") != "X")
   {
      std::cout << "Usage ./bin/makeJetResponseTreeUpdated2.exe" << std::endl;
      std::cout << "   --h                   display help" << std::endl;
      std::cout << "   --input FileName      input filename" << std::endl;
      std::cout << "   --tree TreeName       what tree to run on" << std::endl;
      std::cout << "   [--isPP false]        is PP?" << std::endl;
      std::cout << "   [--fraction 1.0]      fraction to process" << std::endl;
      std::cout << "   [--roo false]         do RooResponse?" << std::endl;
      std::cout << "   [--reducesys false]   reduced systematics?" << std::endl;
      return 0;
   }

   std::string Input = CL.Get("input");
   std::string TreeName = CL.Get("tree");
   bool IsPP = CL.GetBool("isPP", false);
   double Fraction = CL.GetDouble("fraction", 1.0);
   bool DoRoo = CL.GetBool("roo", false);
   bool ReduceSys = CL.GetBool("reducesys", false);

   int retVal = makeJetResponseTree(Input, TreeName, IsPP, Fraction, DoRoo, ReduceSys);
   std::cout << "Job complete. Return " << retVal << "." << std::endl;

   return retVal;
}
