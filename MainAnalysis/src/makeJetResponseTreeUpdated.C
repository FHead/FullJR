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
#include "Utility/include/returnRootFileContentsList.h"
#include "Utility/include/scaleErrorTool.h"
#include "Utility/include/specialHYDJETEventExclude.h"
#include "Utility/include/stringUtil.h"

   template <typename T>
std::string to_string_with_precision(const T a_value, const int n)
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

int makeJetResponseTree(const std::string inName, bool isPP = false, double inEntryFrac = 1., const bool doRooResponse = false, const bool doSystReduced = false)
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

   if(inName.find(".root") != std::string::npos){
      fileList.push_back(inName);
      pthats.push_back(1.);
      pthatWeights.push_back(1.);
   }
   else if(inName.find(".txt") != std::string::npos){
      std::ifstream file(inName.c_str());
      std::string tempStr;

      while(std::getline(file, tempStr)){
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
   else{
      std::cout << "Given inName \'" << inName << "\' is invalid. return 1" << std::endl;
      return 1;
   }

   if(fileList.size() == 0){
      std::cout << "Given inName \'" << inName << "\' is gives no valid root files. return 1" << std::endl;
      return 1;
   }
   else if(pthats.size() == 0){
      std::cout << "Given inName \'" << inName << "\' contains no pthat list. return 1" << std::endl;
      return 1;
   }
   else if(pthatWeights.size() == 0){
      std::cout << "Given inName \'" << inName << "\' contains no pthatWeights list. return 1" << std::endl;
      return 1;
   }

   if(doLocalDebug || doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

   //Post possible returns, start your random number generator
   TRandom3* randGen_p = new TRandom3(0);

   unsigned int pos = 0;
   while(pos < pthats.size()){
      bool isMoved = false;
      for(unsigned int pI = pos+1; pI < pthats.size(); ++pI){

         if(pthats[pI] < pthats[pos]){
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

   TFile* inFile_p = TFile::Open(mntToXRootdFileString(fileList[0]).c_str(), "READ");
   std::vector<std::string> responseTrees = returnRootFileContentsList(inFile_p, "TTree", "JetAna");  

   pos = 0;
   while(responseTrees.size() > pos){
      //For testing, uncomment to exclude all but R=0.4 trees
      //    if(responseTrees[pos].find("akCs4") == std::string::npos) responseTrees.erase(responseTrees.begin()+pos);

      if(isPP){
         if(responseTrees[pos].find("akCs") != std::string::npos) responseTrees.erase(responseTrees.begin()+pos);
         else if(responseTrees[pos].find("akPu") != std::string::npos) responseTrees.erase(responseTrees.begin()+pos);
         else ++pos;
      }
      else{
         if(responseTrees[pos].find("akCs3PF") != std::string::npos) responseTrees.erase(responseTrees.begin()+pos);
         else if(responseTrees[pos].find("akCs4PF") != std::string::npos) responseTrees.erase(responseTrees.begin()+pos); 
         else if(responseTrees[pos].find("akPu3PF") != std::string::npos) responseTrees.erase(responseTrees.begin()+pos);
         else if(responseTrees[pos].find("akPu4PF") != std::string::npos) responseTrees.erase(responseTrees.begin()+pos);
         else if(responseTrees[pos].find("akCs3PU3PFJet") != std::string::npos) responseTrees.erase(responseTrees.begin()+pos);
         else if(responseTrees[pos].find("akCs4PU3PFJet") != std::string::npos) responseTrees.erase(responseTrees.begin()+pos);
         else if(responseTrees[pos].find("akCs6PU3PFJet") != std::string::npos) responseTrees.erase(responseTrees.begin()+pos);
         else if(responseTrees[pos].find("akCs8PU3PFJet") != std::string::npos) responseTrees.erase(responseTrees.begin()+pos);
         else if(responseTrees[pos].find("akCs10PU3PFJet") != std::string::npos) responseTrees.erase(responseTrees.begin()+pos);
         else ++pos;
      }
   }

   if(doLocalDebug || doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

   inFile_p->Close();
   delete inFile_p;
   inFile_p = NULL;

   Int_t posR4Temp = -1;
   Int_t posGeneralTemp = -1;
   const Int_t nMaxTrees = 10;
   const Int_t nTrees = responseTrees.size(); 

   Int_t rValI[nMaxTrees];
   Double_t rValD[nMaxTrees];

   smallOrLargeR rReader;

   for(unsigned int tI = 0; tI < responseTrees.size(); ++tI){
      std::string dirName = responseTrees[tI];
      rValI[tI] = getRVal(dirName);
      rValD[tI] = (getRVal(dirName))/10.;
   }

   if(nTrees > nMaxTrees){
      std::cout << "nTrees \'" << nTrees << "\' is greater than nMaxTrees \'" << nMaxTrees << "\'. return 1" << std::endl;
      return 1;
   }

   std::cout << "Making response matrices for the following " << nTrees << " jet trees: " << std::endl;
   for(int jI = 0; jI < nTrees; ++jI){
      std::cout << " " << jI << "/" << nTrees << ": " << responseTrees[jI] << std::endl;

      if(responseTrees[jI].find("akCs4") != std::string::npos && posR4Temp < 0) posR4Temp = jI;
      else if(responseTrees[jI].find("ak4") != std::string::npos && posR4Temp < 0) posR4Temp = jI;
      else posGeneralTemp = jI;
   }
   if(posGeneralTemp < 0 && nTrees == 1) posGeneralTemp = 0;

   const Int_t posR4 = posR4Temp;
   const Int_t posGeneral = posGeneralTemp;

   std::cout << "PosR4: " << posR4 << std::endl;
   std::cout << "PosGeneral: " << posGeneral << std::endl;

   const Int_t nCentBinsPerma = 4;
   const Int_t centBinsLowPerma[nCentBinsPerma] = {0, 10, 30, 50};
   const Int_t centBinsHiPerma[nCentBinsPerma] = {10, 30, 50, 90};

   Int_t nCentBinsTemp = nCentBinsPerma;//1;
   //  if(!isPP) nCentBinsTemp = nCentBinsPerma;

   const Int_t nMaxCentBins = 4;
   const Int_t nCentBins = nCentBinsTemp;

   if(nCentBins > nMaxCentBins){
      std::cout << "nCentBins \'" << nCentBins << "\' is greater than nMaxCentBins \'" << nMaxCentBins << "\'. return 1" << std::endl;
      return 1;
   }

   std::vector<Int_t> centBinsLow, centBinsHi;
   if(isPP && false){ // WE WANT PP DEPENDENT CENTRALITY
      centBinsLow.push_back(0);
      centBinsHi.push_back(100);
   }
   else{
      for(Int_t cI = 0; cI < nCentBinsPerma; ++cI){
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

   TDatime* date = new TDatime();
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

   /*
      const Int_t nRecoJtPtBins = 10;
      const Double_t jtPtBins[nRecoJtPtBins+1] = {100., 200., 300., 400., 500., 600., 700., 800., 900., 1000., 1100.};
      */

   const int nR = rReader.GetNR();
   std::vector<int> rVals = rReader.GetRVals();

   const int nPtBinsMax = 200;
   int nGenJtPtBins[nMaxTrees][nMaxCentBins];
   double genJtPtBins[nMaxTrees][nMaxCentBins][nPtBinsMax+1];
   int nRecoJtPtBins[nMaxTrees][nMaxCentBins];
   double recoJtPtBins[nMaxTrees][nMaxCentBins][nPtBinsMax+1];

   for(Int_t jI = 0; jI < nTrees; ++jI){
      for(Int_t cI = 0; cI < nCentBins; ++cI){
         std::string centStr = "Cent" + std::to_string(centBinsLow[cI]) + "to" + std::to_string(centBinsHi[cI]);

         nGenJtPtBins[jI][cI] = rReader.GetGenNBinsFromRValCent(rValI[jI], centStr);
         nRecoJtPtBins[jI][cI] = rReader.GetRecoNBinsFromRValCent(rValI[jI], centStr);

         rReader.GetGenBinsFromRValCent(rValI[jI], centStr, genJtPtBins[jI][cI]);
         rReader.GetRecoBinsFromRValCent(rValI[jI], centStr, recoJtPtBins[jI][cI]);      
      }
   }

   for(Int_t jI = 0; jI < nTrees; ++jI){   
      for(Int_t cI = 0; cI < nCentBins; ++cI){
         std::string centStr = "Cent" + std::to_string(centBinsLow[cI]) + "to" + std::to_string(centBinsHi[cI]);

         std::cout << "jI, tree: " << jI << ", " << centStr << std::endl;

         for(Int_t bI = 0; bI < nGenJtPtBins[jI][cI]+1; ++bI){
            std::cout << " " << genJtPtBins[jI][cI][bI] << ",";
         }
         std::cout << std::endl;

         for(Int_t bI = 0; bI < nRecoJtPtBins[jI][cI]+1; ++bI){
            std::cout << " " << recoJtPtBins[jI][cI][bI] << ",";
         }
         std::cout << std::endl;
      }
   }

   Int_t bigJetRecoTrunc = -1;
   for(Int_t jI = 0; jI < nRecoJtPtBins[0][0]; ++jI){
      Double_t val = (recoJtPtBins[0][0][jI] + recoJtPtBins[0][0][jI+1])/2.;
      if(val > 200.){
         bigJetRecoTrunc = jI+1;
         break;
      }
   }

   std::vector<Double_t> minJtPtCut;
   std::vector<Double_t> multiJtPtCut;
   std::vector<Int_t> recoTruncPos;

   for(Int_t jI = 0; jI < nTrees; ++jI){
      multiJtPtCut.push_back(50.);
      minJtPtCut.push_back(20.);
      recoTruncPos.push_back(1);

      bool isBigJt = responseTrees[jI].find("ak8") != std::string::npos || responseTrees[jI].find("ak10") != std::string::npos || responseTrees[jI].find("akCs8") != std::string::npos || responseTrees[jI].find("akCs10") != std::string::npos;
      if(isBigJt){
         multiJtPtCut[jI] = 100;
         minJtPtCut[jI] = 20.;
         recoTruncPos[jI] = bigJetRecoTrunc;
      }
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
   for(Int_t sI = 0; sI < nSyst; ++sI){
      if(systType[sI] == PriorFlat){
         priorFlatPos = sI;
         break;
      }
   }

   Int_t nSystReduced = nSyst;
   if(doSystReduced) nSystReduced = 3;
   /*
      const Int_t nSyst = 1;
      const std::string systStr[nSyst] = {""};
      */

   //LightMUAndCHID == jtPfCHMF < 0.9 && jtPfMUMF < 0.6

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

   for(Int_t jI = 0; jI < nTrees; ++jI){
      for(Int_t cI = 0; cI < nCentBins; ++cI){
         std::string centStr = "Cent" + std::to_string(centBinsLow[cI]) + "to" + std::to_string(centBinsHi[cI]);

         cutProp.SetGenPtBins(rValI[jI], centStr, nGenJtPtBins[jI][cI], genJtPtBins[jI][cI]);
         cutProp.SetRecoPtBins(rValI[jI], centStr, nRecoJtPtBins[jI][cI], recoJtPtBins[jI][cI]);
      }
   }

   cutProp.SetNGeneralBins(rReader.GetNGeneralBins());
   cutProp.SetGeneralBins(rReader.GetGeneralBins());

   cutProp.SetNJtAlgos(nTrees);
   cutProp.SetJtAlgos(responseTrees);
   cutProp.SetMinJtPtCut(minJtPtCut);
   cutProp.SetMultiJtPtCut(multiJtPtCut);
   cutProp.SetRecoTruncPos(recoTruncPos);
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

   TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");
   //Following two lines necessary else you get absolutely clobbered on deletion timing. See threads here:
   //https://root-forum.cern.ch/t/closing-root-files-is-dead-slow-is-this-a-normal-thing/5273/16
   //https://root-forum.cern.ch/t/tfile-speed/17549/25
   //Bizarre
   outFile_p->SetBit(TFile::kDevNull);
   TH1::AddDirectory(kFALSE);

   TDirectory* generalDir_p = (TDirectory*)outFile_p->mkdir("generalHistDir");
   generalDir_p->cd();
   TH1D* pthat_h = new TH1D("pthat_h", ";p_{T} Hat;Counts (Unweighted)", nPthatBins, pthatLow, pthatHi);
   TH1D* pthatWeighted_h = new TH1D("pthatWeighted_h", ";p_{T} Hat;Counts (Weighted)", nPthatBins, pthatLow, pthatHi);
   TH1D* pthatFullWeighted_h = new TH1D("pthatFullWeighted_h", ";p_{T} Hat;Counts (Full Weighted)", nPthatBins, pthatLow, pthatHi);
   TH1D* pthatFullRatio_h = new TH1D("pthatFullRatio_h", ";p_{T} Hat;Ratio", nPthatBins, pthatLow, pthatHi);

   centerTitles({pthat_h, pthatWeighted_h, pthatFullWeighted_h, pthatFullRatio_h});
   setSumW2({pthat_h, pthatWeighted_h, pthatFullWeighted_h, pthatFullRatio_h});

   TH1D* centrality_h = NULL;
   TH1D* centralityWeighted_h = NULL;
   TH1D* centralityFullWeighted_h = NULL;
   TH1D* centralityFullRatio_h = NULL;

   if(!isPP){
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

   for(Int_t tI = 0; tI < nTrees; ++tI){
      std::string dirName = responseTrees[tI];
      dirName = dirName.substr(0, dirName.find("/"));

      for(Int_t bIX = 0; bIX < nGeneralBins+1; ++bIX){
         if(minRecoJtPt > generalBins[bIX]) minRecoJtPt = generalBins[bIX];
         if(minGenJtPt > generalBins[bIX]) minGenJtPt = generalBins[bIX];
      }

      for(Int_t cI = 0; cI < nCentBins; ++cI){
         for(Int_t bIX = 0; bIX < nRecoJtPtBins[tI][cI]+1; ++bIX){
            if(minRecoJtPt > recoJtPtBins[tI][cI][bIX]) minRecoJtPt = recoJtPtBins[tI][cI][bIX];
         }

         for(Int_t bIX = 0; bIX < nGenJtPtBins[tI][cI]+1; ++bIX){
            if(minGenJtPt > genJtPtBins[tI][cI][bIX]) minGenJtPt = genJtPtBins[tI][cI][bIX]; 
         }
      }
   }  

   std::cout << "MINGENJTPT: " << minGenJtPt << ", " << minJtPtCut[0] << std::endl;

   TDirectory* dir_p[nMaxTrees] = {NULL};

   TH1D* deltaPtOrig_h[nMaxTrees][nMaxCentBins];
   TH1D* deltaPhiOrig_h[nMaxTrees][nMaxCentBins];
   TH1D* deltaEtaOrig_h[nMaxTrees][nMaxCentBins];
   TH1D* deltaROrig_h[nMaxTrees][nMaxCentBins];

   TH1D* deltaPtReplace_h[nMaxTrees][nMaxCentBins];
   TH1D* deltaPhiReplace_h[nMaxTrees][nMaxCentBins];
   TH1D* deltaEtaReplace_h[nMaxTrees][nMaxCentBins];
   TH1D* deltaRReplace_h[nMaxTrees][nMaxCentBins];

   TH1D* deltaPtNotReplace_h[nMaxTrees][nMaxCentBins];
   TH1D* deltaPhiNotReplace_h[nMaxTrees][nMaxCentBins];
   TH1D* deltaEtaNotReplace_h[nMaxTrees][nMaxCentBins];
   TH1D* deltaRNotReplace_h[nMaxTrees][nMaxCentBins];

   TH1D* recoJtPt_h[nMaxTrees][nMaxCentBins][nID][nResponseMod][nJtAbsEtaBins];
   TH1D* recoJtPt_RecoTrunc_h[nMaxTrees][nMaxCentBins][nID][nResponseMod][nJtAbsEtaBins];
   TH1D* recoJtPt_NoTruth_h[nMaxTrees][nMaxCentBins][nID][nResponseMod][nJtAbsEtaBins];

   TH1D* recoJtPt_ParaFills_h[nMaxTrees][nMaxCentBins][nID][nResponseMod][nJtAbsEtaBins];
   TH1D* recoJtPt_RecoTrunc_ParaFills_h[nMaxTrees][nMaxCentBins][nID][nResponseMod][nJtAbsEtaBins];
   TH1D* recoJtPt_NoTruth_ParaFills_h[nMaxTrees][nMaxCentBins][nID][nResponseMod][nJtAbsEtaBins];

   TH1D* recoJtPt_GoodGen_h[nMaxTrees][nMaxCentBins][nID][nResponseMod][nJtAbsEtaBins][nSyst];
   TH1D* recoJtPt_GoodGen_ParaFills_h[nMaxTrees][nMaxCentBins][nID][nResponseMod][nJtAbsEtaBins][nSyst];

   TH1D* recoJtPt_General_h[nMaxTrees][nMaxCentBins][nID][nResponseMod][nJtAbsEtaBins][nSyst];
   TH1D* recoJtPt_General_AllReco_h[nMaxTrees][nMaxCentBins][nID][nResponseMod][nJtAbsEtaBins][nSyst];
   TH1D* recoJtPt_General_ParaFills_h[nMaxTrees][nMaxCentBins][nID][nResponseMod][nJtAbsEtaBins][nSyst];

   TH1D* genJtPt_h[nMaxTrees][nMaxCentBins][nID][nResponseMod][nJtAbsEtaBins];
   TH1D* genJtPt_ParaFills_h[nMaxTrees][nMaxCentBins][nID][nResponseMod][nJtAbsEtaBins];
   TH1D* genJtPt_All_h[nMaxTrees][nMaxCentBins][nID][nResponseMod][nJtAbsEtaBins];
   TH1D* genJtPt_GoodReco_h[nMaxTrees][nMaxCentBins][nID][nResponseMod][nJtAbsEtaBins][nSyst];
   TH1D* genJtPt_GoodReco_ParaFills_h[nMaxTrees][nMaxCentBins][nID][nResponseMod][nJtAbsEtaBins][nSyst];
   TH1D* genJtPt_General_h[nMaxTrees][nMaxCentBins][nID][nResponseMod][nJtAbsEtaBins][nSyst];
   TH1D* genJtPt_General_AllGen_h[nMaxTrees][nMaxCentBins][nID][nResponseMod][nJtAbsEtaBins][nSyst];
   TH1D* genJtPt_General_ParaFills_h[nMaxTrees][nMaxCentBins][nID][nResponseMod][nJtAbsEtaBins][nSyst];
   TH2D* response_RecoGenSymm_h[nMaxTrees][nMaxCentBins][nID][nResponseMod][nJtAbsEtaBins][nSyst];
   TH2D* response_RecoGenAsymm_h[nMaxTrees][nMaxCentBins][nID][nResponseMod][nJtAbsEtaBins][nSyst];

   TH2D* response_General_h[nMaxTrees][nMaxCentBins][nID][nResponseMod][nJtAbsEtaBins][nSyst];
   TH2D* response_General_Half1_h[nMaxTrees][nMaxCentBins][nID][nResponseMod][nJtAbsEtaBins][nSyst];
   TH2D* response_General_Half2_h[nMaxTrees][nMaxCentBins][nID][nResponseMod][nJtAbsEtaBins][nSyst];
   TH2D* response_General_ParaFills_h[nMaxTrees][nMaxCentBins][nID][nResponseMod][nJtAbsEtaBins][nSyst];
   TH2D* responseNorm_General_h[nMaxTrees][nMaxCentBins][nID][nResponseMod][nJtAbsEtaBins][nSyst];

   RooUnfoldResponse* rooResponse_RecoGenAsymm_h[nMaxTrees][nMaxCentBins][nID][nResponseMod][nJtAbsEtaBins][nSyst];
   TH1D* genJtPt_CheckPriorFlat_h[nMaxTrees][nMaxCentBins][nID][nResponseMod][nJtAbsEtaBins];

   for(Int_t dI = 0; dI < nTrees; ++dI){
      if(nTrees == 2 && dI == posR4 && !isPP) continue;

      outFile_p->cd();
      std::string dirName = responseTrees[dI];
      dirName = dirName.substr(0, dirName.find("/"));

      dir_p[dI] = (TDirectory*)outFile_p->mkdir(dirName.c_str());

      for(Int_t cI = 0; cI < nCentBins; ++cI){
         std::string centStr = "Cent" + std::to_string(centBinsLow[cI]) + "to" + std::to_string(centBinsHi[cI]);      
         if(isPP) centStr = "PP_" + centStr;
         else centStr = "PbPb_" + centStr;

         deltaPtOrig_h[dI][cI] = new TH1D(("deltaPtOrig_" + dirName + "_" + centStr + "_h").c_str(), ";#Delta pT/Reco Jet pT;Counts", 100, 0, 2);
         deltaPhiOrig_h[dI][cI] = new TH1D(("deltaPhiOrig_" + dirName + "_" + centStr + "_h").c_str(), ";#Delta #phi;Counts", 100, -0.5, 0.5);
         deltaEtaOrig_h[dI][cI] = new TH1D(("deltaEtaOrig_" + dirName + "_" + centStr + "_h").c_str(), ";#Delta #eta;Counts", 100, -0.5, 0.5);
         deltaROrig_h[dI][cI] = new TH1D(("deltaROrig_" + dirName + "_" + centStr + "_h").c_str(), ";#Delta R;Counts", 100, -0.5, 0.5);

         deltaPtReplace_h[dI][cI] = new TH1D(("deltaPtReplace_" + dirName + "_" + centStr + "_h").c_str(), ";#Delta pT/Reco Jet pT;Counts", 100, 0, 2);
         deltaPhiReplace_h[dI][cI] = new TH1D(("deltaPhiReplace_" + dirName + "_" + centStr + "_h").c_str(), ";#Delta #phi;Counts", 100, -0.5, 0.5);
         deltaEtaReplace_h[dI][cI] = new TH1D(("deltaEtaReplace_" + dirName + "_" + centStr + "_h").c_str(), ";#Delta #eta;Counts", 100, -0.5, 0.5);
         deltaRReplace_h[dI][cI] = new TH1D(("deltaRReplace_" + dirName + "_" + centStr + "_h").c_str(), ";#Delta R;Counts", 100, -0.5, 0.5);

         deltaPtNotReplace_h[dI][cI] = new TH1D(("deltaPtNotReplace_" + dirName + "_" + centStr + "_h").c_str(), ";#Delta pT/Reco Jet pT;Counts", 100, 0, 2);
         deltaPhiNotReplace_h[dI][cI] = new TH1D(("deltaPhiNotReplace_" + dirName + "_" + centStr + "_h").c_str(), ";#Delta #phi;Counts", 100, -0.5, 0.5);
         deltaEtaNotReplace_h[dI][cI] = new TH1D(("deltaEtaNotReplace_" + dirName + "_" + centStr + "_h").c_str(), ";#Delta #eta;Counts", 100, -0.5, 0.5);
         deltaRNotReplace_h[dI][cI] = new TH1D(("deltaRNotReplace_" + dirName + "_" + centStr + "_h").c_str(), ";#Delta R;Counts", 100, -0.5, 0.5);

         for(Int_t iI = 0; iI < nID; ++iI){
            for(Int_t mI = 0; mI < nResponseMod; ++mI){
               std::string resStr = "ResponseMod" + prettyString(responseMod[mI], 2, true);

               for(Int_t aI = 0; aI < nJtAbsEtaBins; ++aI){
                  const std::string jtAbsEtaStr = "AbsEta" + prettyString(jtAbsEtaBinsLow[aI], 1, true) + "to" + prettyString(jtAbsEtaBinsHi[aI], 1, true);

                  recoJtPt_h[dI][cI][iI][mI][aI] = new TH1D(("recoJtPt_" + dirName + "_" + centStr + "_" + idStr[iI] + "_" + resStr + "_" + jtAbsEtaStr + "_h").c_str(), ";Reco. Jet p_{T};Counts (Weighted)", nRecoJtPtBins[dI][cI], recoJtPtBins[dI][cI]);
                  recoJtPt_RecoTrunc_h[dI][cI][iI][mI][aI] = new TH1D(("recoJtPt_" + dirName + "_" + centStr + "_" + idStr[iI] + "_" + resStr + "_" + jtAbsEtaStr + "_RecoTrunc_h").c_str(), ";Reco. Jet p_{T};Counts (Weighted)", nRecoJtPtBins[dI][cI], recoJtPtBins[dI][cI]);
                  recoJtPt_NoTruth_h[dI][cI][iI][mI][aI] = new TH1D(("recoJtPt_" + dirName + "_" + centStr + "_" + idStr[iI] + "_" + resStr + "_" + jtAbsEtaStr + "_NoTruth_h").c_str(), ";Reco. Jet p_{T};Counts (Weighted)", nRecoJtPtBins[dI][cI], recoJtPtBins[dI][cI]);

                  recoJtPt_ParaFills_h[dI][cI][iI][mI][aI] = new TH1D(("recoJtPt_" + dirName + "_" + centStr + "_" + idStr[iI] + "_" + resStr + "_" + jtAbsEtaStr + "_ParaFills_h").c_str(), ";Reco. Jet p_{T};Counts (Weighted)", nRecoJtPtBins[dI][cI], recoJtPtBins[dI][cI]);
                  recoJtPt_RecoTrunc_ParaFills_h[dI][cI][iI][mI][aI] = new TH1D(("recoJtPt_" + dirName + "_" + centStr + "_" + idStr[iI] + "_" + resStr + "_" + jtAbsEtaStr + "_RecoTrunc_ParaFills_h").c_str(), ";Reco. Jet p_{T};Counts (Weighted)", nRecoJtPtBins[dI][cI], recoJtPtBins[dI][cI]);
                  recoJtPt_NoTruth_ParaFills_h[dI][cI][iI][mI][aI] = new TH1D(("recoJtPt_" + dirName + "_" + centStr + "_" + idStr[iI] + "_" + resStr + "_" + jtAbsEtaStr + "_NoTruth_ParaFills_h").c_str(), ";Reco. Jet p_{T};Counts (Weighted)", nRecoJtPtBins[dI][cI], recoJtPtBins[dI][cI]);

                  genJtPt_h[dI][cI][iI][mI][aI] = new TH1D(("genJtPt_" + dirName + "_" + centStr + "_" + idStr[iI] + "_" + resStr + "_" + jtAbsEtaStr + "_h").c_str(), ";Gen. Jet p_{T};Counts (Weighted)", nGenJtPtBins[dI][cI], genJtPtBins[dI][cI]);

                  genJtPt_ParaFills_h[dI][cI][iI][mI][aI] = new TH1D(("genJtPt_" + dirName + "_" + centStr + "_" + idStr[iI] + "_" + resStr + "_" + jtAbsEtaStr + "_ParaFills_h").c_str(), ";Gen. Jet p_{T};Counts (Weighted)", nGenJtPtBins[dI][cI], genJtPtBins[dI][cI]);;	    	  

                  genJtPt_All_h[dI][cI][iI][mI][aI] = new TH1D(("genJtPt_" + dirName + "_" + centStr + "_" + idStr[iI] + "_" + resStr + "_" + jtAbsEtaStr + "_All_h").c_str(), ";Gen. Jet p_{T};Counts (Weighted)", nGenJtPtBins[dI][cI], genJtPtBins[dI][cI]);
                  genJtPt_CheckPriorFlat_h[dI][cI][iI][mI][aI] = new TH1D(("rooResponse_" + dirName + "_" + centStr + "_" + idStr[iI] + "_" + resStr + "_" + jtAbsEtaStr + "_CheckPriorFlat_h").c_str(), ";Gen. p_{T} (Weighted Flat);Counts (Weighted)", nGenJtPtBins[dI][cI], genJtPtBins[dI][cI]);


                  for(Int_t sI = 0; sI < nSystReduced; ++sI){
                     std::string tempSysStr =  "_" + systStr[sI] + "_";	 
                     while(tempSysStr.find("__") != std::string::npos){tempSysStr.replace(tempSysStr.find("__"), 2, "_");}

                     recoJtPt_GoodGen_h[dI][cI][iI][mI][aI][sI] = new TH1D(("recoJtPt_" + dirName + "_" + centStr + "_" + idStr[iI] + "_" + resStr + "_" + jtAbsEtaStr + tempSysStr + "GoodGen_h").c_str(), ";Reco. Jet p_{T};Counts (Weighted)", nRecoJtPtBins[dI][cI], recoJtPtBins[dI][cI]);

                     recoJtPt_GoodGen_ParaFills_h[dI][cI][iI][mI][aI][sI] = new TH1D(("recoJtPt_" + dirName + "_" + centStr + "_" + idStr[iI] + "_" + resStr + "_" + jtAbsEtaStr + tempSysStr + "GoodGen_ParaFills_h").c_str(), ";Reco. Jet p_{T};Counts (Weighted)", nRecoJtPtBins[dI][cI], recoJtPtBins[dI][cI]);

                     recoJtPt_General_h[dI][cI][iI][mI][aI][sI] = new TH1D(("recoJtPt_" + dirName + "_" + centStr + "_" + idStr[iI] + "_" + resStr + "_" + jtAbsEtaStr + tempSysStr + "General_h").c_str(), ";Reco. Jet p_{T};Counts (Weighted)", 148, 20, 1500);
                     recoJtPt_General_AllReco_h[dI][cI][iI][mI][aI][sI] = new TH1D(("recoJtPt_" + dirName + "_" + centStr + "_" + idStr[iI] + "_" + resStr + "_" + jtAbsEtaStr + tempSysStr + "General_AllReco_h").c_str(), ";Reco. Jet p_{T};Counts (Weighted)", 148, 20, 1500);
                     recoJtPt_General_ParaFills_h[dI][cI][iI][mI][aI][sI] = new TH1D(("recoJtPt_" + dirName + "_" + centStr + "_" + idStr[iI] + "_" + resStr + "_" + jtAbsEtaStr + tempSysStr + "General_ParaFills_h").c_str(), ";Reco. Jet p_{T};Counts (Weighted)", 148, 20, 1500);
                     
                     genJtPt_General_h[dI][cI][iI][mI][aI][sI] = new TH1D(("genJtPt_" + dirName + "_" + centStr + "_" + idStr[iI] + "_" + resStr + "_" + jtAbsEtaStr + tempSysStr + "General_h").c_str(), ";Gen. Jet p_{T};Counts (Weighted)", 148, 20, 1500);

                     genJtPt_General_AllGen_h[dI][cI][iI][mI][aI][sI] = new TH1D(("genJtPt_" + dirName + "_" + centStr + "_" + idStr[iI] + "_" + resStr + "_" + jtAbsEtaStr + tempSysStr + "General_AllGen_h").c_str(), ";Gen. Jet p_{T};Counts (Weighted)", 148, 20, 1500);

                     genJtPt_General_ParaFills_h[dI][cI][iI][mI][aI][sI] = new TH1D(("genJtPt_" + dirName + "_" + centStr + "_" + idStr[iI] + "_" + resStr + "_" + jtAbsEtaStr + tempSysStr + "General_ParaFills_h").c_str(), ";Gen. Jet p_{T};Counts (Weighted)", 148, 20, 1500);

                     // std::cout << "response_General_h " << " " << nGeneralBins << std::endl;
                     response_General_h[dI][cI][iI][mI][aI][sI] = new TH2D(("response_" + dirName + "_" + centStr + "_" + idStr[iI] + "_" + resStr + "_" + jtAbsEtaStr + tempSysStr + "General_h").c_str(), ";Reco. Jet p_{T};Gen. Jet p_{T}", 148, 20, 1500, 148, 20, 1500);

                     response_General_Half1_h[dI][cI][iI][mI][aI][sI] = new TH2D(("response_" + dirName + "_" + centStr + "_" + idStr[iI] + "_" + resStr + "_" + jtAbsEtaStr + tempSysStr + "General_Half1_h").c_str(), ";Reco. Jet p_{T};Gen. Jet p_{T}", 148, 20, 1500, 148, 20, 1500);
                     response_General_Half2_h[dI][cI][iI][mI][aI][sI] = new TH2D(("response_" + dirName + "_" + centStr + "_" + idStr[iI] + "_" + resStr + "_" + jtAbsEtaStr + tempSysStr + "General_Half2_h").c_str(), ";Reco. Jet p_{T};Gen. Jet p_{T}", 148, 20, 1500, 148, 20, 1500);

                     response_General_ParaFills_h[dI][cI][iI][mI][aI][sI] = new TH2D(("response_" + dirName + "_" + centStr + "_" + idStr[iI] + "_" + resStr + "_" + jtAbsEtaStr + tempSysStr + "General_ParaFills_h").c_str(), ";Reco. Jet p_{T};Gen. Jet p_{T}", 148, 20, 1500, 148, 20, 1500);

                     responseNorm_General_h[dI][cI][iI][mI][aI][sI] = new TH2D(("responseNorm_" + dirName + "_" + centStr + "_" + idStr[iI] + "_" + resStr + "_" + jtAbsEtaStr + tempSysStr + "General_h").c_str(), ";Reco. Jet p_{T};Gen. Jet p_{T}", nNormBins, 0, 3, 148, 20, 1500);

                     std::vector<TH1*> tempVect = {recoJtPt_GoodGen_h[dI][cI][iI][mI][aI][sI], recoJtPt_GoodGen_ParaFills_h[dI][cI][iI][mI][aI][sI], recoJtPt_General_h[dI][cI][iI][mI][aI][sI], recoJtPt_General_AllReco_h[dI][cI][iI][mI][aI][sI], recoJtPt_General_ParaFills_h[dI][cI][iI][mI][aI][sI], genJtPt_General_h[dI][cI][iI][mI][aI][sI], genJtPt_General_AllGen_h[dI][cI][iI][mI][aI][sI], genJtPt_General_ParaFills_h[dI][cI][iI][mI][aI][sI], response_General_h[dI][cI][iI][mI][aI][sI], response_General_Half1_h[dI][cI][iI][mI][aI][sI], response_General_Half2_h[dI][cI][iI][mI][aI][sI], response_General_ParaFills_h[dI][cI][iI][mI][aI][sI], responseNorm_General_h[dI][cI][iI][mI][aI][sI]};
                     centerTitles(tempVect);
                     setSumW2(tempVect);

                     genJtPt_GoodReco_h[dI][cI][iI][mI][aI][sI] = new TH1D(("genJtPt_" + dirName + "_" + centStr + "_" + idStr[iI] + "_" + resStr + "_" + jtAbsEtaStr + tempSysStr + "GoodReco_h").c_str(), ";Gen. Jet p_{T};Counts (Weighted)", nGenJtPtBins[dI][cI], genJtPtBins[dI][cI]);

                     genJtPt_GoodReco_ParaFills_h[dI][cI][iI][mI][aI][sI] = new TH1D(("genJtPt_" + dirName + "_" + centStr + "_" + idStr[iI] + "_" + resStr + "_" + jtAbsEtaStr + tempSysStr + "GoodReco_ParaFills_h").c_str(), ";Gen. Jet p_{T};Counts (Weighted)", nGenJtPtBins[dI][cI], genJtPtBins[dI][cI]);

                     // std::cout << "response_RecoGenSymm_h " << " " << nRecoJtPtBins[dI][cI] << " " << nGenJtPtBins[dI][cI] << std::endl;
                     response_RecoGenSymm_h[dI][cI][iI][mI][aI][sI] = new TH2D(("response_" + dirName + "_" + centStr + "_" + idStr[iI] + "_" + resStr + "_" + jtAbsEtaStr + tempSysStr + "RecoGenSymm_h").c_str(), ";Reco. Jet p_{T};Gen. Jet p_{T}", nRecoJtPtBins[dI][cI], recoJtPtBins[dI][cI], nGenJtPtBins[dI][cI], genJtPtBins[dI][cI]);

                     response_RecoGenAsymm_h[dI][cI][iI][mI][aI][sI] = new TH2D(("response_" + dirName + "_" + centStr + "_" + idStr[iI] + "_" + resStr + "_" + jtAbsEtaStr + tempSysStr + "RecoGenAsymm_h").c_str(), ";Reco. Jet p_{T};Gen. Jet p_{T}", nRecoJtPtBins[dI][cI], recoJtPtBins[dI][cI], nGenJtPtBins[dI][cI], genJtPtBins[dI][cI]);

                     if(doRooResponse){
                        rooResponse_RecoGenAsymm_h[dI][cI][iI][mI][aI][sI] = new RooUnfoldResponse(("rooResponse_" + dirName + "_" + centStr + "_" + idStr[iI] + "_" + resStr + "_" + jtAbsEtaStr + tempSysStr + "RecoGenAsymm_h").c_str(), "");

                        rooResponse_RecoGenAsymm_h[dI][cI][iI][mI][aI][sI]->Setup(recoJtPt_GoodGen_h[dI][cI][iI][mI][aI][sI], genJtPt_GoodReco_h[dI][cI][iI][mI][aI][sI]);
                     }

                     tempVect = {genJtPt_GoodReco_h[dI][cI][iI][mI][aI][sI], genJtPt_GoodReco_ParaFills_h[dI][cI][iI][mI][aI][sI], response_RecoGenSymm_h[dI][cI][iI][mI][aI][sI], response_RecoGenAsymm_h[dI][cI][iI][mI][aI][sI]};
                     centerTitles(tempVect);
                     setSumW2(tempVect);
                  }


                  std::vector<TH1*> tempVect = {recoJtPt_RecoTrunc_h[dI][cI][iI][mI][aI], recoJtPt_NoTruth_h[dI][cI][iI][mI][aI], recoJtPt_h[dI][cI][iI][mI][aI], recoJtPt_RecoTrunc_ParaFills_h[dI][cI][iI][mI][aI], recoJtPt_NoTruth_ParaFills_h[dI][cI][iI][mI][aI], recoJtPt_ParaFills_h[dI][cI][iI][mI][aI], genJtPt_h[dI][cI][iI][mI][aI], genJtPt_ParaFills_h[dI][cI][iI][mI][aI], genJtPt_All_h[dI][cI][iI][mI][aI], genJtPt_CheckPriorFlat_h[dI][cI][iI][mI][aI]};
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

   Int_t nref_[nMaxTrees];
   Float_t jtpt_[nMaxTrees][nMaxJet];
   Float_t rawpt_[nMaxTrees][nMaxJet];
   Float_t jteta_[nMaxTrees][nMaxJet];
   Float_t jtphi_[nMaxTrees][nMaxJet];
   Float_t refpt_[nMaxTrees][nMaxJet];
   Float_t refeta_[nMaxTrees][nMaxJet];
   Float_t refphi_[nMaxTrees][nMaxJet];
   Float_t jtPfCHF_[nMaxTrees][nMaxJet];
   Float_t jtPfCEF_[nMaxTrees][nMaxJet];
   Float_t jtPfNHF_[nMaxTrees][nMaxJet];
   Float_t jtPfNEF_[nMaxTrees][nMaxJet];
   Float_t jtPfMUF_[nMaxTrees][nMaxJet];
   Float_t jtPfCHMF_[nMaxTrees][nMaxJet];
   Float_t jtPfCEMF_[nMaxTrees][nMaxJet];
   Float_t jtPfNHMF_[nMaxTrees][nMaxJet];
   Float_t jtPfNEMF_[nMaxTrees][nMaxJet];
   Float_t jtPfMUMF_[nMaxTrees][nMaxJet];
   Int_t jtPfCHM_[nMaxTrees][nMaxJet];
   Int_t jtPfCEM_[nMaxTrees][nMaxJet];
   Int_t jtPfNHM_[nMaxTrees][nMaxJet];
   Int_t jtPfNEM_[nMaxTrees][nMaxJet];
   Int_t jtPfMUM_[nMaxTrees][nMaxJet];

   Int_t ngen_[nMaxTrees];
   Float_t genpt_[nMaxTrees][nMaxJet];
   Float_t genphi_[nMaxTrees][nMaxJet];
   Float_t geneta_[nMaxTrees][nMaxJet];
   Int_t gensubid_[nMaxTrees][nMaxJet];

   for(unsigned int fI = 0; fI < fileList.size(); ++fI){
      std::cout << "Processing file " << fI << "/" << fileList.size() << ": \'" << fileList[fI] << "\'" << std::endl;

      inFile_p = TFile::Open(mntToXRootdFileString(fileList[fI]).c_str(), "READ");
      TTree* jetTrees_p[nMaxTrees] = {NULL};

      for(Int_t tI = 0; tI < nTrees; ++tI){
         jetTrees_p[tI] = (TTree*)inFile_p->Get(responseTrees[tI].c_str());

         if(nTrees != 2 || tI != posR4){
            jetTrees_p[tI]->SetBranchStatus("*", 0);
            jetTrees_p[tI]->SetBranchStatus("nref", 1);
            jetTrees_p[tI]->SetBranchStatus("jtpt", 1);
            jetTrees_p[tI]->SetBranchStatus("rawpt", 1);
            jetTrees_p[tI]->SetBranchStatus("jteta", 1);
            jetTrees_p[tI]->SetBranchStatus("jtphi", 1);
            jetTrees_p[tI]->SetBranchStatus("refpt", 1);
            jetTrees_p[tI]->SetBranchStatus("refphi", 1);
            jetTrees_p[tI]->SetBranchStatus("refeta", 1);
            jetTrees_p[tI]->SetBranchStatus("jtPfCHF", 1);
            jetTrees_p[tI]->SetBranchStatus("jtPfCEF", 1);
            jetTrees_p[tI]->SetBranchStatus("jtPfNHF", 1);
            jetTrees_p[tI]->SetBranchStatus("jtPfNEF", 1);
            jetTrees_p[tI]->SetBranchStatus("jtPfMUF", 1);
            jetTrees_p[tI]->SetBranchStatus("jtPfCHMF", 1);
            jetTrees_p[tI]->SetBranchStatus("jtPfCEMF", 1);
            jetTrees_p[tI]->SetBranchStatus("jtPfNHMF", 1);
            jetTrees_p[tI]->SetBranchStatus("jtPfNEMF", 1);
            jetTrees_p[tI]->SetBranchStatus("jtPfMUMF", 1);
            jetTrees_p[tI]->SetBranchStatus("jtPfCHM", 1);
            jetTrees_p[tI]->SetBranchStatus("jtPfCEM", 1);
            jetTrees_p[tI]->SetBranchStatus("jtPfNHM", 1);
            jetTrees_p[tI]->SetBranchStatus("jtPfNEM", 1);
            jetTrees_p[tI]->SetBranchStatus("jtPfMUM", 1);
         }

         jetTrees_p[tI]->SetBranchStatus("ngen", 1);
         jetTrees_p[tI]->SetBranchStatus("genpt", 1);
         jetTrees_p[tI]->SetBranchStatus("geneta", 1);
         jetTrees_p[tI]->SetBranchStatus("genphi", 1);
         jetTrees_p[tI]->SetBranchStatus("gensubid", 1);

         if(nTrees != 2 || tI != posR4){
            jetTrees_p[tI]->SetBranchAddress("nref", &(nref_[tI]));
            jetTrees_p[tI]->SetBranchAddress("jtpt", jtpt_[tI]);
            jetTrees_p[tI]->SetBranchAddress("rawpt", rawpt_[tI]);
            jetTrees_p[tI]->SetBranchAddress("jteta", jteta_[tI]);
            jetTrees_p[tI]->SetBranchAddress("jtphi", jtphi_[tI]);
            jetTrees_p[tI]->SetBranchAddress("refpt", refpt_[tI]);
            jetTrees_p[tI]->SetBranchAddress("refphi", refphi_[tI]);
            jetTrees_p[tI]->SetBranchAddress("refeta", refeta_[tI]);
            jetTrees_p[tI]->SetBranchAddress("jtPfCHF", jtPfCHF_[tI]);
            jetTrees_p[tI]->SetBranchAddress("jtPfCEF", jtPfCEF_[tI]);
            jetTrees_p[tI]->SetBranchAddress("jtPfNHF", jtPfNHF_[tI]);
            jetTrees_p[tI]->SetBranchAddress("jtPfNEF", jtPfNEF_[tI]);
            jetTrees_p[tI]->SetBranchAddress("jtPfMUF", jtPfMUF_[tI]);
            jetTrees_p[tI]->SetBranchAddress("jtPfCHMF", jtPfCHMF_[tI]);
            jetTrees_p[tI]->SetBranchAddress("jtPfCEMF", jtPfCEMF_[tI]);
            jetTrees_p[tI]->SetBranchAddress("jtPfNHMF", jtPfNHMF_[tI]);
            jetTrees_p[tI]->SetBranchAddress("jtPfNEMF", jtPfNEMF_[tI]);
            jetTrees_p[tI]->SetBranchAddress("jtPfMUMF", jtPfMUMF_[tI]);
            jetTrees_p[tI]->SetBranchAddress("jtPfCHM", jtPfCHM_[tI]);
            jetTrees_p[tI]->SetBranchAddress("jtPfCEM", jtPfCEM_[tI]);
            jetTrees_p[tI]->SetBranchAddress("jtPfNHM", jtPfNHM_[tI]);
            jetTrees_p[tI]->SetBranchAddress("jtPfNEM", jtPfNEM_[tI]);
            jetTrees_p[tI]->SetBranchAddress("jtPfMUM", jtPfMUM_[tI]);
         }

         jetTrees_p[tI]->SetBranchAddress("ngen", &(ngen_[tI]));
         jetTrees_p[tI]->SetBranchAddress("genpt", genpt_[tI]);
         jetTrees_p[tI]->SetBranchAddress("geneta", geneta_[tI]);
         jetTrees_p[tI]->SetBranchAddress("genphi", genphi_[tI]);
         jetTrees_p[tI]->SetBranchAddress("gensubid", gensubid_[tI]);
      }

      Float_t pthat_;

      jetTrees_p[posGeneral]->SetBranchStatus("pthat", 1);
      jetTrees_p[posGeneral]->SetBranchAddress("pthat", &pthat_);

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

      TTree* hiTree_p = (TTree*)inFile_p->Get("hiEvtAnalyzer/HiTree");

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

      TTree* skimTree_p = (TTree*)inFile_p->Get("skimanalysis/HltTree");

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

      const Int_t nEntries = TMath::Min((Int_t)1000000000, (Int_t)jetTrees_p[posGeneral]->GetEntries());
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
      for(Int_t entry = 0; entry < nEntriesToProcess; ++entry){
         if(nEntriesToProcess >= 50000 && entry%printInterval == 0){
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

         for(Int_t tI = 0; tI < nTrees; ++tI){jetTrees_p[tI]->GetEntry(entry);}

         Double_t ncollWeight_ = 1.;      
         if(!isPP){
            centPos[0] = -1;
            for(Int_t cI = 0; cI < nCentBins; ++cI){
               if(centBinsLow[cI]*2 <= hiBin_ && hiBin_ < centBinsHi[cI]*2){
                  centPos[0] = cI;
                  break;
               }
            }
            if(centPos[0] < 0) continue;

            bool badJetSpecialSel = specialSel.CheckEventBadJet(ngen_[posR4], genpt_[posR4], genphi_[posR4], geneta_[posR4], gensubid_[posR4], entry);
            if(badJetSpecialSel) continue;

            ncollWeight_ = findNcoll_Renorm(hiBin_);
         }

         Double_t pthatWeight_ = -1;
         for(unsigned int pI = 0; pI < pthats.size()-1; ++pI){
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

         for(Int_t tI = 0; tI < nTrees; ++tI){
            if(nTrees == 2 && posR4 == tI) continue;

            std::string algoName = responseTrees[tI];
            //	if(algoName.find("/") != std::string::npos) algoName = algoName.substr(0, algoName.find("/"));

            const Int_t nUsed = ngen_[tI];
            bool isUsed[nUsed];
            for(Int_t gI = 0; gI < nUsed; ++gI){isUsed[gI] = false;}

            nonJet.stop();

            //WE HAVE TO CLEAN BS REFPT MATCHES
            for(Int_t jI = 0; jI < nref_[tI]; ++jI){
               if(jtpt_[tI][jI] < 20) continue;
               if(TMath::Abs(jteta_[tI][jI]) > 3.) continue;
               if(refpt_[tI][jI] < 0) continue;
               if(refpt_[tI][jI]/jtpt_[tI][jI] > 0.25) continue;

               Double_t deltaPTInit = refpt_[tI][jI]/jtpt_[tI][jI];
               Double_t deltaPhiInit = (jtphi_[tI][jI] - refphi_[tI][jI]);
               Double_t deltaEtaInit = (jteta_[tI][jI] - refeta_[tI][jI]);
               Double_t deltaRInit = getDR(jteta_[tI][jI], jtphi_[tI][jI], refeta_[tI][jI], refphi_[tI][jI]);
               bool isReplace = false;

               for(Int_t gI = 0; gI < ngen_[tI]; ++gI){
                  //	    if(gensubid_[tI][gI] == 0) continue;
                  Double_t deltaR = getDR(jteta_[tI][jI], jtphi_[tI][jI], geneta_[tI][gI], genphi_[tI][gI]);
                  Double_t deltaPT = genpt_[tI][gI]/jtpt_[tI][jI];
                  Double_t deltaPhi = (jtphi_[tI][jI] - genphi_[tI][gI]);
                  Double_t deltaEta = (jteta_[tI][jI] - geneta_[tI][gI]);

                  if(deltaR < 0.15 && genpt_[tI][gI]/jtpt_[tI][jI] > 0.5){
                     //	      std::cout << "WARNING CLEANING JET (pt,phi,eta): " << jtpt_[tI][jI] << ", " << jtphi_[tI][jI] << ", " << jteta_[tI][jI] << " matched to " << refpt_[tI][jI] << ", " << refphi_[tI][jI] << ", " << refeta_[tI][jI] << "(deltaRInit " << deltaRInit << "), better match found is " << genpt_[tI][gI] << ", " << genphi_[tI][gI] << ", " << geneta_[tI][gI] << "(deltaR " << deltaR << ")" << std::endl;

                     if(centPos.size() != 0){
                        deltaPtOrig_h[tI][centPos[0]]->Fill(deltaPTInit);
                        deltaPhiOrig_h[tI][centPos[0]]->Fill(deltaPhiInit);
                        deltaEtaOrig_h[tI][centPos[0]]->Fill(deltaEtaInit);
                        deltaROrig_h[tI][centPos[0]]->Fill(deltaRInit);

                        deltaPtReplace_h[tI][centPos[0]]->Fill(deltaPT);
                        deltaPhiReplace_h[tI][centPos[0]]->Fill(deltaPhi);
                        deltaEtaReplace_h[tI][centPos[0]]->Fill(deltaEta);
                        deltaRReplace_h[tI][centPos[0]]->Fill(deltaR);
                     }

                     refpt_[tI][jI] = -999;
                     refphi_[tI][jI] = -999;
                     refeta_[tI][jI] = -999;

                     isReplace = true;
                     break;
                  }
               }

               if(!isReplace && jtpt_[tI][jI] > 100 && TMath::Abs(jteta_[tI][jI]) < 2.){
                  //	    std::cout << "WARNING JET " << jtpt_[tI][jI] << ", " << jtphi_[tI][jI] << ", " << jteta_[tI][jI] << " matched to " << refpt_[tI][jI] << ", " << refphi_[tI][jI] << ", " << refeta_[tI][jI] << "(deltaRInit " << deltaRInit << ") NOT REPLACED IN EVENT " << entry << std::endl;	    

                  if(centPos.size() != 0){
                     deltaPtNotReplace_h[tI][centPos[0]]->Fill(deltaPTInit);
                     deltaPhiNotReplace_h[tI][centPos[0]]->Fill(deltaPhiInit);
                     deltaEtaNotReplace_h[tI][centPos[0]]->Fill(deltaEtaInit);
                     deltaRNotReplace_h[tI][centPos[0]]->Fill(deltaRInit);
                  }

                  for(Int_t gI = 0; gI < ngen_[tI]; ++gI){
                     //	      std::cout << " WARNING JET POSSIBLE MATCH: " << gI << ": " << genpt_[tI][gI] << ", " << genphi_[tI][gI] << ", " << geneta_[tI][gI] << ", " << gensubid_[tI][gI] << std::endl;
                  }
               }
            }

            for(Int_t jI = 0; jI < nref_[tI]; ++jI){
               jetLoop1.start();

               if(TMath::Abs(jteta_[tI][jI]) > jtAbsEtaMax) continue;
               Double_t refptTemp = refpt_[tI][jI];
               if(refptTemp < 0){
                  for(Int_t gI = 0; gI < ngen_[tI]; ++gI){
                     if(gensubid_[tI][gI] == 0) continue;
                     else if(isUsed[gI]) continue;
                     else if(getDR(jteta_[tI][jI], jtphi_[tI][jI], geneta_[tI][gI], genphi_[tI][gI]) < 0.2 + TMath::Min(0.0, 0.3 - rValD[tI])){
                        isUsed[gI] = true;
                        refptTemp = genpt_[tI][gI];
                        break;
                     }
                  }
               }

               jetLoop1.stop();

               if(refptTemp < minGenJtPt) continue;

               std::vector<bool> passesID;
               for(Int_t iI = 0; iI < nID; ++iI){
                  bool pass = jtPfCHMFCutLow[iI] <= jtPfCHMF_[tI][jI] && jtPfCHMF_[tI][jI] <= jtPfCHMFCutHi[iI];
                  pass = pass && jtPfMUMFCutLow[iI] <= jtPfMUMF_[tI][jI] && jtPfMUMF_[tI][jI] <= jtPfMUMFCutHi[iI];
                  pass = pass && jtPfNHFCutLow[iI] <= jtPfNHF_[tI][jI] && jtPfNHF_[tI][jI] <= jtPfNHFCutHi[iI];
                  pass = pass && jtPfNEFCutLow[iI] <= jtPfNEF_[tI][jI] && jtPfNEF_[tI][jI] <= jtPfNEFCutHi[iI];
                  pass = pass && jtPfMUFCutLow[iI] <= jtPfMUF_[tI][jI] && jtPfMUF_[tI][jI] <= jtPfMUFCutHi[iI];
                  pass = pass && jtPfCHFCutLow[iI] <= jtPfCHF_[tI][jI] && jtPfCHF_[tI][jI] <= jtPfCHFCutHi[iI];
                  pass = pass && jtPfCEFCutLow[iI] <= jtPfCEF_[tI][jI] && jtPfCEF_[tI][jI] <= jtPfCEFCutHi[iI];
                  pass = pass && jtPfCEM_[tI][jI] + jtPfNEM_[tI][jI] + jtPfCHM_[tI][jI] + jtPfNHM_[tI][jI] + jtPfMUM_[tI][jI] >= jtPfMinMult[iI];
                  pass = pass && jtPfCHM_[tI][jI] >= jtPfMinChgMult[iI];

                  passesID.push_back(pass);
               }

               std::vector<int> jtAbsEtaPoses;
               for(Int_t aI = 0; aI < nJtAbsEtaBins; ++aI){
                  if(TMath::Abs(jteta_[tI][jI]) >= jtAbsEtaBinsLow[aI] && TMath::Abs(jteta_[tI][jI]) < jtAbsEtaBinsHi[aI]){
                     jtAbsEtaPoses.push_back(aI);
                  }
               }

               for(Int_t mI = 0; mI < nResponseMod; ++mI){
                  jetLoop2.start();

                  Double_t jtPtFillVal[nSyst];
                  double fullWeight2[nSyst];

                  bool oneRecoIsGood = false;

                  for(Int_t sI = 0; sI < nSystReduced; ++sI){	  
                     jtPtFillVal[sI] = jtpt_[tI][jI] + (jtpt_[tI][jI] - refpt_[tI][jI]) * responseMod[mI];
                     fullWeight2[sI] = fullWeight_;

                     if((systType[sI] == JECUpMC)) jtPtFillVal[sI] += jtPtFillVal[sI]*jecVarMC;
                     else if((systType[sI] == JECDownMC)) jtPtFillVal[sI] -= jtPtFillVal[sI]*jecVarMC;
                     else if((systType[sI] == JECUpData)) jtPtFillVal[sI] += jtPtFillVal[sI]*jecVarData;
                     else if((systType[sI] == JECDownData)) jtPtFillVal[sI] -= jtPtFillVal[sI]*jecVarData;
                     else if((systType[sI] == JECUpUE) || (systType[sI] == JECDownUE) )
                     {
                        Double_t tempScale =  jtPtFillVal[sI]/rawpt_[tI][jI];
                        Double_t tempRawPt = rawpt_[tI][jI];
                        if((systType[sI] == JECUpUE)) tempRawPt += TMath::Abs(scaleErr.getMuDataMinusMC(scaleHiBin, jteta_[tI][jI], rValI[tI], FlowDefaultInRho));
                        else tempRawPt -= TMath::Abs(scaleErr.getMuDataMinusMC(scaleHiBin, jteta_[tI][jI], rValI[tI], FlowDefaultInRho));

                        jtPtFillVal[sI] = tempRawPt*tempScale;
                     }
                     else if((systType[sI] == JERMC)) jtPtFillVal[sI] += (jtPtFillVal[sI] - refpt_[tI][jI])*jerVarMC;
                     else if((systType[sI] == JERData)) jtPtFillVal[sI] += (jtPtFillVal[sI] - refpt_[tI][jI])*jerVarData[mI];
                     else if((systType[sI] == PriorUp1PowerPthat)) fullWeight2[sI] *= pthat_/pthatLow;
                     else if((systType[sI] == PriorDown1PowerPthat)) fullWeight2[sI] *= pthatLow/pthat_;
                     else if((systType[sI] == PriorFlat)) fullWeight2[sI] *= 1./pthatWeight_;
                     /*
                        else if((systType[sI] == PriorFlat) && refpt_[tI][jI] > 0) fullWeight2[sI] *= flatWeight.getJtWeight(algoName == hiBin_/2 == refpt_[tI][jI] == jteta_[tI][jI]);
                        else if((systType[sI] == PriorNoWeight) && refpt_[tI][jI] > 0) fullWeight2[sI] /= pthatWeight_;
                        */
                     if(jtPtFillVal[sI] > minRecoJtPt) oneRecoIsGood = true;
                  }

                  jetLoop2.stop();

                  if(!oneRecoIsGood) continue;

                  jetLoop3.start();

                  for(auto const & cent : centPos){
                     jetLoop3Sub1.start();

                     bool goodTruth = (refpt_[tI][jI] >= genJtPtBins[tI][cent][0] && refpt_[tI][jI] < genJtPtBins[tI][cent][nGenJtPtBins[tI][cent]] && refpt_[tI][jI] > minJtPtCut[tI]);

                     bool fillLowTruth = false;
                     bool fillHighTruth = false;
                     if(!goodTruth && refpt_[tI][jI] > 0.0){
                        if(refpt_[tI][jI] < genJtPtBins[tI][cent][0]) fillLowTruth = true;
                        else if(refpt_[tI][jI] >= genJtPtBins[tI][cent][nGenJtPtBins[tI][cent]]) fillHighTruth = true;
                     }

                     bool goodTruthGeneral = (refpt_[tI][jI] >= generalBins[0] && refpt_[tI][jI] < generalBins[nGeneralBins] && refpt_[tI][jI] > minJtPtCut[tI]);  

                     std::vector<Int_t> goodRecoPos;
                     std::vector<Int_t> goodRecoTruncPos;
                     std::vector<Int_t> goodRecoGeneralPos;

                     for(Int_t sI = 0; sI < nSystReduced; ++sI){
                        if(jtPtFillVal[sI] >= genJtPtBins[tI][cent][0] && jtPtFillVal[sI] < genJtPtBins[tI][cent][nGenJtPtBins[tI][cent]]) goodRecoPos.push_back(sI);
                        if(jtPtFillVal[sI] >= recoJtPtBins[tI][cent][0] && jtPtFillVal[sI] < recoJtPtBins[tI][cent][nRecoJtPtBins[tI][cent]]) goodRecoTruncPos.push_back(sI);
                        if(jtPtFillVal[sI] >= generalBins[0] && jtPtFillVal[sI] < generalBins[nGeneralBins]) goodRecoGeneralPos.push_back(sI);
                     }

                     jetLoop3Sub1.stop();
                     jetLoop3Sub2.start();

                     if(refpt_[tI][jI] < 0){
                        bool goodTruth2 = (refptTemp >= genJtPtBins[tI][cent][0] && refptTemp < genJtPtBins[tI][cent][nGenJtPtBins[tI][cent]] && refptTemp > minJtPtCut[tI]);
                        if(!goodTruth2){
                           for(unsigned int aI = 0; aI < jtAbsEtaPoses.size(); ++aI){
                              for(unsigned int iI = 0; iI < passesID.size(); ++iI){
                                 if(!passesID[iI]) continue;

                                 if(isPara) recoJtPt_NoTruth_ParaFills_h[tI][cent][iI][mI][jtAbsEtaPoses[aI]]->Fill(jtPtFillVal[0], fullWeight_);
                                 else recoJtPt_NoTruth_h[tI][cent][iI][mI][jtAbsEtaPoses[aI]]->Fill(jtPtFillVal[0], fullWeight_);
                              }
                           }
                        }	    
                     }

                     jetLoop3Sub2.stop();
                     jetLoop3Sub3.start();	      

                     jetLoop3Sub3Sub0.start();	      
                     if(!isPara){
                        for(unsigned int iI = 0; iI < passesID.size(); ++iI){
                           if(!passesID[iI]) continue;
                           for(unsigned int aI = 0; aI < jtAbsEtaPoses.size(); ++aI){
                              for(unsigned int sI = 0; sI < goodRecoGeneralPos.size(); ++sI){
                                 recoJtPt_General_AllReco_h[tI][cent][iI][mI][jtAbsEtaPoses[aI]][goodRecoGeneralPos[sI]]->Fill(refpt_[tI][jI], fullWeight2[goodRecoGeneralPos[sI]]);
                              }
                           }
                        }
                     }
                     
                     if(goodTruthGeneral){
                        for(unsigned int iI = 0; iI < passesID.size(); ++iI){
                           if(!passesID[iI]) continue;
                           for(unsigned int aI = 0; aI < jtAbsEtaPoses.size(); ++aI){

                              if(!isPara){
                                 for(Int_t sI = 0; sI < nSystReduced; ++sI){
                                    genJtPt_General_AllGen_h[tI][cent][iI][mI][jtAbsEtaPoses[aI]][sI]->Fill(refpt_[tI][jI], fullWeight2[sI]);
                                 }
                              }

                              for(unsigned int sI = 0; sI < goodRecoGeneralPos.size(); ++sI){
                                 int BRefPT = FindBin(148, 20, 1500, refpt_[tI][jI]);
                                 int BJTPT = FindBin(148, 20, 1500, jtPtFillVal[goodRecoGeneralPos[sI]]);
                                 if(/*isPara*/ false){
                                    genJtPt_General_ParaFills_h[tI][cent][iI][mI][jtAbsEtaPoses[aI]][goodRecoGeneralPos[sI]]->AddBinContent(BRefPT, fullWeight2[goodRecoGeneralPos[sI]]);
                                    int B = response_General_ParaFills_h[tI][cent][iI][mI][jtAbsEtaPoses[aI]][goodRecoGeneralPos[sI]]->GetBin(BJTPT, BRefPT);
                                    response_General_ParaFills_h[tI][cent][iI][mI][jtAbsEtaPoses[aI]][goodRecoGeneralPos[sI]]->AddBinContent(B, fullWeight2[goodRecoGeneralPos[sI]]);
                                    recoJtPt_General_ParaFills_h[tI][cent][iI][mI][jtAbsEtaPoses[aI]][goodRecoGeneralPos[sI]]->AddBinContent(BJTPT, fullWeight2[goodRecoGeneralPos[sI]]);
                                 }
                                 else{	       
                                    genJtPt_General_h[tI][cent][iI][mI][jtAbsEtaPoses[aI]][goodRecoGeneralPos[sI]]->AddBinContent(BRefPT, fullWeight2[goodRecoGeneralPos[sI]]);
                                    int B =  response_General_h[tI][cent][iI][mI][jtAbsEtaPoses[aI]][goodRecoGeneralPos[sI]]->GetBin(BJTPT, BRefPT);
                                       response_General_h[tI][cent][iI][mI][jtAbsEtaPoses[aI]][goodRecoGeneralPos[sI]]->AddBinContent(B, fullWeight2[goodRecoGeneralPos[sI]]);

                                    if(randGen_p->Uniform(0.0, 1.0) < 0.5) response_General_Half1_h[tI][cent][iI][mI][jtAbsEtaPoses[aI]][goodRecoGeneralPos[sI]]->AddBinContent(B, fullWeight2[goodRecoGeneralPos[sI]]);
                                    else response_General_Half2_h[tI][cent][iI][mI][jtAbsEtaPoses[aI]][goodRecoGeneralPos[sI]]->AddBinContent(B, fullWeight2[goodRecoGeneralPos[sI]]);

                                    responseNorm_General_h[tI][cent][iI][mI][jtAbsEtaPoses[aI]][goodRecoGeneralPos[sI]]->Fill(jtPtFillVal[goodRecoGeneralPos[sI]]/refpt_[tI][jI], refpt_[tI][jI], fullWeight2[goodRecoGeneralPos[sI]]);
                                    recoJtPt_General_h[tI][cent][iI][mI][jtAbsEtaPoses[aI]][goodRecoGeneralPos[sI]]->AddBinContent(BJTPT, fullWeight2[goodRecoGeneralPos[sI]]);			
                                 }
                              }
                           }
                        }

                     }

                     jetLoop3Sub3Sub0.stop();	      

                     if(goodTruth || fillLowTruth || fillHighTruth){
                        Double_t tempRefPt = refpt_[tI][jI];

                        if(fillLowTruth) tempRefPt = genJtPtBins[tI][cent][0]+1;
                        else if(fillHighTruth) tempRefPt = genJtPtBins[tI][cent][nGenJtPtBins[tI][cent]]-1;

                        for(unsigned int iI = 0; iI < passesID.size(); ++iI){
                           if(!passesID[iI]) continue;
                           for(unsigned int aI = 0; aI < jtAbsEtaPoses.size(); ++aI){

                              if(!isPara){
                                 jetLoop3Sub3Suba.start();	      

                                 genJtPt_All_h[tI][cent][iI][mI][jtAbsEtaPoses[aI]]->Fill(tempRefPt, fullWeight_);

                                 for(unsigned int sI = 0; sI < goodRecoPos.size(); ++sI){
                                    response_RecoGenSymm_h[tI][cent][iI][mI][jtAbsEtaPoses[aI]][goodRecoPos[sI]]->Fill(jtPtFillVal[goodRecoPos[sI]], tempRefPt, fullWeight2[goodRecoPos[sI]]);
                                 }  

                                 jetLoop3Sub3Suba.stop();	      		      
                                 jetLoop3Sub3Subb.start();	      

                                 for(unsigned int sI = 0; sI < goodRecoTruncPos.size(); ++sI){
                                    genJtPt_GoodReco_h[tI][cent][iI][mI][jtAbsEtaPoses[aI]][goodRecoTruncPos[sI]]->Fill(tempRefPt, fullWeight2[goodRecoTruncPos[sI]]);
                                    response_RecoGenAsymm_h[tI][cent][iI][mI][jtAbsEtaPoses[aI]][goodRecoTruncPos[sI]]->Fill(jtPtFillVal[goodRecoTruncPos[sI]], tempRefPt, fullWeight2[goodRecoTruncPos[sI]]);
                                    if(doRooResponse) rooResponse_RecoGenAsymm_h[tI][cent][iI][mI][jtAbsEtaPoses[aI]][goodRecoTruncPos[sI]]->Fill(jtPtFillVal[goodRecoTruncPos[sI]], tempRefPt, fullWeight2[goodRecoTruncPos[sI]]);

                                    recoJtPt_GoodGen_h[tI][cent][iI][mI][jtAbsEtaPoses[aI]][goodRecoTruncPos[sI]]->Fill(jtPtFillVal[goodRecoTruncPos[sI]], fullWeight2[goodRecoTruncPos[sI]]);
                                 }

                                 jetLoop3Sub3Subb.stop();	      
                                 jetLoop3Sub3Subc.start();	      

                                 if(priorFlatPos >= 0){
                                    genJtPt_CheckPriorFlat_h[tI][cent][iI][mI][aI]->Fill(tempRefPt, fullWeight2[priorFlatPos]);
                                 }

                                 if(goodRecoPos.size() != 0){
                                    if(goodRecoPos[0] == 0) recoJtPt_h[tI][cent][iI][mI][jtAbsEtaPoses[aI]]->Fill(jtPtFillVal[0], fullWeight_);
                                 }

                                 if(goodRecoTruncPos.size() != 0){
                                    if(goodRecoTruncPos[0] == 0) recoJtPt_RecoTrunc_h[tI][cent][iI][mI][jtAbsEtaPoses[aI]]->Fill(jtPtFillVal[0], fullWeight_);			
                                 }

                                 genJtPt_h[tI][cent][iI][mI][jtAbsEtaPoses[aI]]->Fill(tempRefPt, fullWeight_);

                                 jetLoop3Sub3Subc.stop();	      
                              }		  
                              else{
                                 for(unsigned int sI = 0; sI < goodRecoTruncPos.size(); ++sI){
                                    genJtPt_GoodReco_ParaFills_h[tI][cent][iI][mI][jtAbsEtaPoses[aI]][goodRecoTruncPos[sI]]->Fill(tempRefPt, fullWeight2[goodRecoTruncPos[sI]]);
                                    recoJtPt_GoodGen_ParaFills_h[tI][cent][iI][mI][jtAbsEtaPoses[aI]][goodRecoTruncPos[sI]]->Fill(jtPtFillVal[goodRecoTruncPos[sI]], fullWeight2[goodRecoTruncPos[sI]]);
                                 }		  

                                 if(goodRecoPos.size() != 0){
                                    if(goodRecoPos[0] == 0) recoJtPt_ParaFills_h[tI][cent][iI][mI][jtAbsEtaPoses[aI]]->Fill(jtPtFillVal[0], fullWeight_);
                                 }

                                 if(goodRecoTruncPos.size() != 0){		     
                                    if(goodRecoTruncPos[0] == 0) recoJtPt_RecoTrunc_ParaFills_h[tI][cent][iI][mI][jtAbsEtaPoses[aI]]->Fill(jtPtFillVal[0], fullWeight_);
                                 }

                                 genJtPt_ParaFills_h[tI][cent][iI][mI][jtAbsEtaPoses[aI]]->Fill(tempRefPt, fullWeight_);
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

   if(!isPP){
      centrality_h->Write("", TObject::kOverwrite);
      centralityWeighted_h->Write("", TObject::kOverwrite);
      centralityFullWeighted_h->Write("", TObject::kOverwrite);
      centralityFullRatio_h->Divide(centralityFullWeighted_h, centralityWeighted_h);
      centralityFullRatio_h->Write("", TObject::kOverwrite);
   }

   if(doLocalDebug || doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

   for(Int_t dI = 0; dI < nTrees; ++dI){
      if(nTrees == 2 && dI == posR4 && !isPP) continue;

      outFile_p->cd();
      dir_p[dI]->cd();

      for(Int_t cI = 0; cI < nCentBins; ++cI){
         deltaPtOrig_h[dI][cI]->Write("", TObject::kOverwrite);
         deltaPhiOrig_h[dI][cI]->Write("", TObject::kOverwrite);
         deltaEtaOrig_h[dI][cI]->Write("", TObject::kOverwrite);
         deltaROrig_h[dI][cI]->Write("", TObject::kOverwrite);

         deltaPtReplace_h[dI][cI]->Write("", TObject::kOverwrite);
         deltaPhiReplace_h[dI][cI]->Write("", TObject::kOverwrite);
         deltaEtaReplace_h[dI][cI]->Write("", TObject::kOverwrite);
         deltaRReplace_h[dI][cI]->Write("", TObject::kOverwrite);

         deltaPtNotReplace_h[dI][cI]->Write("", TObject::kOverwrite);
         deltaPhiNotReplace_h[dI][cI]->Write("", TObject::kOverwrite);
         deltaEtaNotReplace_h[dI][cI]->Write("", TObject::kOverwrite);
         deltaRNotReplace_h[dI][cI]->Write("", TObject::kOverwrite);

         for(Int_t iI = 0; iI < nID; ++iI){
            for(Int_t mI = 0; mI < nResponseMod; ++mI){
               for(Int_t aI = 0; aI < nJtAbsEtaBins; ++aI){
                  recoJtPt_h[dI][cI][iI][mI][aI]->Write("", TObject::kOverwrite);
                  recoJtPt_RecoTrunc_h[dI][cI][iI][mI][aI]->Write("", TObject::kOverwrite);
                  recoJtPt_NoTruth_h[dI][cI][iI][mI][aI]->Write("", TObject::kOverwrite);
                  recoJtPt_ParaFills_h[dI][cI][iI][mI][aI]->Write("", TObject::kOverwrite);
                  recoJtPt_RecoTrunc_ParaFills_h[dI][cI][iI][mI][aI]->Write("", TObject::kOverwrite);
                  recoJtPt_NoTruth_ParaFills_h[dI][cI][iI][mI][aI]->Write("", TObject::kOverwrite);

                  genJtPt_h[dI][cI][iI][mI][aI]->Write("", TObject::kOverwrite);
                  genJtPt_ParaFills_h[dI][cI][iI][mI][aI]->Write("", TObject::kOverwrite);
                  genJtPt_All_h[dI][cI][iI][mI][aI]->Write("", TObject::kOverwrite);

                  for(Int_t sI = 0; sI < nSystReduced; ++sI){
                     recoJtPt_GoodGen_h[dI][cI][iI][mI][aI][sI]->Write("", TObject::kOverwrite);
                     recoJtPt_GoodGen_ParaFills_h[dI][cI][iI][mI][aI][sI]->Write("", TObject::kOverwrite);

                     recoJtPt_General_h[dI][cI][iI][mI][aI][sI]->Write("", TObject::kOverwrite);
                     recoJtPt_General_AllReco_h[dI][cI][iI][mI][aI][sI]->Write("", TObject::kOverwrite);

                     genJtPt_General_h[dI][cI][iI][mI][aI][sI]->Write("", TObject::kOverwrite);
                     genJtPt_General_AllGen_h[dI][cI][iI][mI][aI][sI]->Write("", TObject::kOverwrite);
                     response_General_h[dI][cI][iI][mI][aI][sI]->Write("", TObject::kOverwrite);
                     response_General_Half1_h[dI][cI][iI][mI][aI][sI]->Write("", TObject::kOverwrite);
                     response_General_Half2_h[dI][cI][iI][mI][aI][sI]->Write("", TObject::kOverwrite);
                     responseNorm_General_h[dI][cI][iI][mI][aI][sI]->Write("", TObject::kOverwrite);

                     recoJtPt_General_ParaFills_h[dI][cI][iI][mI][aI][sI]->Write("", TObject::kOverwrite);
                     genJtPt_General_ParaFills_h[dI][cI][iI][mI][aI][sI]->Write("", TObject::kOverwrite);
                     response_General_ParaFills_h[dI][cI][iI][mI][aI][sI]->Write("", TObject::kOverwrite);

                     genJtPt_GoodReco_h[dI][cI][iI][mI][aI][sI]->Write("", TObject::kOverwrite);
                     genJtPt_GoodReco_ParaFills_h[dI][cI][iI][mI][aI][sI]->Write("", TObject::kOverwrite);
                     response_RecoGenSymm_h[dI][cI][iI][mI][aI][sI]->Write("", TObject::kOverwrite);
                     response_RecoGenAsymm_h[dI][cI][iI][mI][aI][sI]->Write("", TObject::kOverwrite);
                     if(doRooResponse) rooResponse_RecoGenAsymm_h[dI][cI][iI][mI][aI][sI]->Write("", TObject::kOverwrite);
                  }

                  genJtPt_CheckPriorFlat_h[dI][cI][iI][mI][aI]->Write("", TObject::kOverwrite);
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

   for(Int_t dI = 0; dI < nTrees; ++dI){
      if(nTrees == 2 && dI == posR4 && !isPP) continue;

      outFile_p->cd();
      dir_p[dI]->cd();

      for(Int_t cI = 0; cI < nCentBins; ++cI){
         delete deltaPtOrig_h[dI][cI];
         delete deltaPhiOrig_h[dI][cI];
         delete deltaEtaOrig_h[dI][cI];
         delete deltaROrig_h[dI][cI];

         delete deltaPtReplace_h[dI][cI];
         delete deltaPhiReplace_h[dI][cI];
         delete deltaEtaReplace_h[dI][cI];
         delete deltaRReplace_h[dI][cI];

         delete deltaPtNotReplace_h[dI][cI];
         delete deltaPhiNotReplace_h[dI][cI];
         delete deltaEtaNotReplace_h[dI][cI];
         delete deltaRNotReplace_h[dI][cI];

         for(Int_t iI = 0; iI < nID; ++iI){
            for(Int_t mI = 0; mI < nResponseMod; ++mI){
               for(Int_t aI = 0; aI < nJtAbsEtaBins; ++aI){

                  delete recoJtPt_RecoTrunc_h[dI][cI][iI][mI][aI];
                  delete recoJtPt_NoTruth_h[dI][cI][iI][mI][aI];
                  delete recoJtPt_h[dI][cI][iI][mI][aI];

                  delete recoJtPt_RecoTrunc_ParaFills_h[dI][cI][iI][mI][aI];
                  delete recoJtPt_NoTruth_ParaFills_h[dI][cI][iI][mI][aI];
                  delete recoJtPt_ParaFills_h[dI][cI][iI][mI][aI];

                  delete genJtPt_h[dI][cI][iI][mI][aI];
                  delete genJtPt_ParaFills_h[dI][cI][iI][mI][aI];	   	    
                  delete genJtPt_All_h[dI][cI][iI][mI][aI];

                  for(Int_t sI = 0; sI < nSystReduced; ++sI){
                     delete recoJtPt_GoodGen_h[dI][cI][iI][mI][aI][sI];
                     delete recoJtPt_GoodGen_ParaFills_h[dI][cI][iI][mI][aI][sI];

                     delete recoJtPt_General_h[dI][cI][iI][mI][aI][sI];
                     delete recoJtPt_General_AllReco_h[dI][cI][iI][mI][aI][sI];

                     delete genJtPt_General_h[dI][cI][iI][mI][aI][sI];
                     delete genJtPt_General_AllGen_h[dI][cI][iI][mI][aI][sI];

                     delete response_General_h[dI][cI][iI][mI][aI][sI];
                     delete response_General_Half1_h[dI][cI][iI][mI][aI][sI];
                     delete response_General_Half2_h[dI][cI][iI][mI][aI][sI];
                     delete responseNorm_General_h[dI][cI][iI][mI][aI][sI];

                     delete recoJtPt_General_ParaFills_h[dI][cI][iI][mI][aI][sI];
                     delete genJtPt_General_ParaFills_h[dI][cI][iI][mI][aI][sI];
                     delete response_General_ParaFills_h[dI][cI][iI][mI][aI][sI];

                     delete genJtPt_GoodReco_h[dI][cI][iI][mI][aI][sI];
                     delete genJtPt_GoodReco_ParaFills_h[dI][cI][iI][mI][aI][sI];
                     delete response_RecoGenSymm_h[dI][cI][iI][mI][aI][sI];
                     delete response_RecoGenAsymm_h[dI][cI][iI][mI][aI][sI];
                     if(doRooResponse) delete rooResponse_RecoGenAsymm_h[dI][cI][iI][mI][aI][sI];
                  }

                  delete genJtPt_CheckPriorFlat_h[dI][cI][iI][mI][aI];
               }
            }
         }
      }
   }    

   deleteLoopWatch.stop();

   outFile_p->cd();

   TDirectory* cutDir_p = (TDirectory*)outFile_p->mkdir("cutDir");
   TDirectory* subDir_p = (TDirectory*)cutDir_p->mkdir("subDir");

   if(doLocalDebug || doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

   if(doLocalDebug || doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

   if(nTrees == 2 && !isPP){
      responseTrees.erase(responseTrees.begin() + posR4);
      minJtPtCut.erase(minJtPtCut.begin() + posR4);
      multiJtPtCut.erase(multiJtPtCut.begin() + posR4);
      recoTruncPos.erase(recoTruncPos.begin() + posR4);

      cutProp.SetNJtAlgos(nTrees-1);
      cutProp.SetJtAlgos(responseTrees);
      cutProp.SetMinJtPtCut(minJtPtCut);
      cutProp.SetMultiJtPtCut(multiJtPtCut);
      cutProp.SetRecoTruncPos(recoTruncPos);
   }
   else{
      cutProp.SetNJtAlgos(nTrees);
      cutProp.SetJtAlgos(responseTrees);
      cutProp.SetMinJtPtCut(minJtPtCut);
      cutProp.SetMultiJtPtCut(multiJtPtCut);
      cutProp.SetRecoTruncPos(recoTruncPos);
   }

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
   const double nJetTrees = responseTrees.size();

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

int main(int argc, char* argv[])
{
   if(argc < 2 && argc > 6){
      std::cout << "Usage ./bin/makeJetResponseTree.exe <inName> <isPP-Opt> <entryFrac-Opt> <doRooResponse-opt> <doSystReduced-opt>" << std::endl;
      return 1;
   }

   int retVal = 0;
   if(argc == 2) retVal += makeJetResponseTree(argv[1]);
   else if(argc == 3) retVal += makeJetResponseTree(argv[1], std::stoi(argv[2]));
   else if(argc == 4) retVal += makeJetResponseTree(argv[1], std::stoi(argv[2]), std::stod(argv[3]));
   else if(argc == 5) retVal += makeJetResponseTree(argv[1], std::stoi(argv[2]), std::stod(argv[3]), std::stoi(argv[4]));
   else if(argc == 6) retVal += makeJetResponseTree(argv[1], std::stoi(argv[2]), std::stod(argv[3]), std::stoi(argv[4]), std::stoi(argv[5]));

   std::cout << "Job complete. Return " << retVal << "." << std::endl;
   return retVal;
}
