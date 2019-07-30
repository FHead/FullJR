#ifndef CUTPROPAGATOR_H
#define CUTPROPAGATOR_H

//cpp dependencies
#include <string>
#include <iostream>
#include <vector>
#include <sstream>
#include <iomanip>

//ROOT dependencies
#include "TFile.h"
#include "TDirectory.h"
#include "TNamed.h"
#include "TCollection.h"
#include "TList.h"
#include "TKey.h"
#include "TString.h"

//Local FullJR (MainAnalysis) dependencies
#include "MainAnalysis/include/doLocalDebug.h"

//Non-local FullJR (Utility, etc.) dependencies
#include "Utility/include/doGlobalDebug.h"
#include "Utility/include/plotUtilities.h"
#include "Utility/include/returnRootFileContentsList.h"
#include "Utility/include/stringUtil.h"

class cutPropagator
{
 public:
  std::vector<std::string> inFileNames;
  std::vector<std::string> inFullFileNames;

  bool isPP;

  double jtAbsEtaMax;

  std::string rcDiffFileName;
  std::string flatPriorFileName;
  double jecVarMC;
  double jerVarMC;

  double jecVarData;

  double delta = 0.0001;

  int nResponseMod;
  std::vector<double> responseMod;
  std::vector<double> jerVarData;

  int nJtAlgos;
  std::vector<std::string> jtAlgos;
  std::vector<double> minJtPtCut;
  std::vector<double> multiJtPtCut;
  std::vector<int> recoTruncPos;

  int nR;
  std::vector<int> rVals;

  std::vector<double> recoJtPtBinsR1Cent0to10;
  std::vector<double> recoJtPtBinsR1Cent10to30;
  std::vector<double> recoJtPtBinsR1Cent30to50;
  std::vector<double> recoJtPtBinsR1Cent50to90;

  std::vector<double> genJtPtBinsR1Cent0to10;
  std::vector<double> genJtPtBinsR1Cent10to30;
  std::vector<double> genJtPtBinsR1Cent30to50;
  std::vector<double> genJtPtBinsR1Cent50to90;

  std::vector<double> recoJtPtBinsR2Cent0to10;
  std::vector<double> recoJtPtBinsR2Cent10to30;
  std::vector<double> recoJtPtBinsR2Cent30to50;
  std::vector<double> recoJtPtBinsR2Cent50to90;

  std::vector<double> genJtPtBinsR2Cent0to10;
  std::vector<double> genJtPtBinsR2Cent10to30;
  std::vector<double> genJtPtBinsR2Cent30to50;
  std::vector<double> genJtPtBinsR2Cent50to90;

  std::vector<double> recoJtPtBinsR3Cent0to10;
  std::vector<double> recoJtPtBinsR3Cent10to30;
  std::vector<double> recoJtPtBinsR3Cent30to50;
  std::vector<double> recoJtPtBinsR3Cent50to90;

  std::vector<double> genJtPtBinsR3Cent0to10;
  std::vector<double> genJtPtBinsR3Cent10to30;
  std::vector<double> genJtPtBinsR3Cent30to50;
  std::vector<double> genJtPtBinsR3Cent50to90;

  std::vector<double> recoJtPtBinsR4Cent0to10;
  std::vector<double> recoJtPtBinsR4Cent10to30;
  std::vector<double> recoJtPtBinsR4Cent30to50;
  std::vector<double> recoJtPtBinsR4Cent50to90;

  std::vector<double> genJtPtBinsR4Cent0to10;
  std::vector<double> genJtPtBinsR4Cent10to30;
  std::vector<double> genJtPtBinsR4Cent30to50;
  std::vector<double> genJtPtBinsR4Cent50to90;

  std::vector<double> recoJtPtBinsR5Cent0to10;
  std::vector<double> recoJtPtBinsR5Cent10to30;
  std::vector<double> recoJtPtBinsR5Cent30to50;
  std::vector<double> recoJtPtBinsR5Cent50to90;

  std::vector<double> genJtPtBinsR5Cent0to10;
  std::vector<double> genJtPtBinsR5Cent10to30;
  std::vector<double> genJtPtBinsR5Cent30to50;
  std::vector<double> genJtPtBinsR5Cent50to90;

  std::vector<double> recoJtPtBinsR6Cent0to10;
  std::vector<double> recoJtPtBinsR6Cent10to30;
  std::vector<double> recoJtPtBinsR6Cent30to50;
  std::vector<double> recoJtPtBinsR6Cent50to90;

  std::vector<double> genJtPtBinsR6Cent0to10;
  std::vector<double> genJtPtBinsR6Cent10to30;
  std::vector<double> genJtPtBinsR6Cent30to50;
  std::vector<double> genJtPtBinsR6Cent50to90;

  std::vector<double> recoJtPtBinsR8Cent0to10;
  std::vector<double> recoJtPtBinsR8Cent10to30;
  std::vector<double> recoJtPtBinsR8Cent30to50;
  std::vector<double> recoJtPtBinsR8Cent50to90;

  std::vector<double> genJtPtBinsR8Cent0to10;
  std::vector<double> genJtPtBinsR8Cent10to30;
  std::vector<double> genJtPtBinsR8Cent30to50;
  std::vector<double> genJtPtBinsR8Cent50to90;

  std::vector<double> recoJtPtBinsR10Cent0to10;
  std::vector<double> recoJtPtBinsR10Cent10to30;
  std::vector<double> recoJtPtBinsR10Cent30to50;
  std::vector<double> recoJtPtBinsR10Cent50to90;

  std::vector<double> genJtPtBinsR10Cent0to10;
  std::vector<double> genJtPtBinsR10Cent10to30;
  std::vector<double> genJtPtBinsR10Cent30to50;
  std::vector<double> genJtPtBinsR10Cent50to90;

  int nGeneralBins;
  std::vector<double> generalBins;

  int nJtAbsEtaBins;
  std::vector<double> jtAbsEtaBinsLow;
  std::vector<double> jtAbsEtaBinsHi;
  
  int nPthats;
  std::vector<double> pthats;
  std::vector<double> pthatWeights;

  int nCentBins;
  std::vector<int> centBinsLow;
  std::vector<int> centBinsHi;

  int nID;
  std::vector<std::string> idStr;
  std::vector<double> jtPfCHMFCutLow;
  std::vector<double> jtPfCHMFCutHi;
  std::vector<double> jtPfMUMFCutLow;
  std::vector<double> jtPfMUMFCutHi;
  std::vector<double> jtPfNHFCutLow;
  std::vector<double> jtPfNHFCutHi;
  std::vector<double> jtPfNEFCutLow;
  std::vector<double> jtPfNEFCutHi;
  std::vector<double> jtPfMUFCutLow;
  std::vector<double> jtPfMUFCutHi;
  std::vector<double> jtPfCHFCutLow;
  std::vector<double> jtPfCHFCutHi;
  std::vector<double> jtPfCEFCutLow;
  std::vector<double> jtPfCEFCutHi;
  std::vector<int> jtPfMinMult;
  std::vector<int> jtPfMinChgMult;

  int nSyst;
  std::vector<std::string> systStr;

  int nBayes;
  int nBigBayesSymm;
  std::vector<int> bayesVal;

  int nSuperBayes;

  int nSVD;

  int nHistDim;
  std::vector<std::string> histTagBayes;
  std::vector<std::string> histTagSVD;
  std::vector<int> histBestBayes;
  std::vector<int> histBestSVD;

  void Clean();
  bool GetAllVarFromFile(TFile* inFile_p);
  std::string VectToString(std::vector<double> inVect, int prec);
  bool WriteAllVarToFile(TFile* inFile_p, TDirectory* inDir_p, TDirectory* inSubDir_p, TDirectory* unfoldDirBayes_p, TDirectory* unfoldDirSVD_p);
  bool CheckPropagatorsMatch(cutPropagator inCutProp, bool doBothMCOrBothData, bool doBothPPOrBothPbPb, bool skipJtAlgo);

  bool CheckDouble(double in1, double in2);
  bool CheckInt(int in1, int in2);
  bool CheckBool(bool in1, bool in2);
  bool CheckString(std::string in1, std::string in2);
  bool CheckVectDouble(std::vector<double> in1, std::vector<double> in2);
  bool CheckVectDouble(std::vector<double> in1, std::vector<double> in2, std::string valStr);
  bool CheckVectInt(std::vector<int> in1, std::vector<int> in2);
  bool CheckVectString(std::vector<std::string> in1, std::vector<std::string> in2);

  bool CheckJtAbsEtaMax(double inJtAbsEtaMax);
  bool CheckJtAbsEtaMax(cutPropagator inCutProp);

  bool CheckNR(int inNR);
  bool CheckNR(cutPropagator cutProp);
  bool CheckRVals(std::vector<int> inRVals);
  bool CheckRVals(cutPropagator cutProp);

  bool CheckNBins(int inBins, int compBins, std::string binStr);
  bool CheckStartEndVals(double inVals, double compVals, std::string valStr);

  bool CheckGenPtBinsR1Cent0to10(std::vector<double> inGenPtBinsR1Cent0to10);
  bool CheckGenPtBinsR1Cent0to10(cutPropagator inCutProp);
  bool CheckRecoPtBinsR1Cent0to10(std::vector<double> inRecoPtBinsR1Cent0to10);
  bool CheckRecoPtBinsR1Cent0to10(cutPropagator inCutProp);

  bool CheckGenPtBinsR1Cent10to30(std::vector<double> inGenPtBinsR1Cent10to30);
  bool CheckGenPtBinsR1Cent10to30(cutPropagator inCutProp);
  bool CheckRecoPtBinsR1Cent10to30(std::vector<double> inRecoPtBinsR1Cent10to30);
  bool CheckRecoPtBinsR1Cent10to30(cutPropagator inCutProp);

  bool CheckGenPtBinsR1Cent30to50(std::vector<double> inGenPtBinsR1Cent30to50);
  bool CheckGenPtBinsR1Cent30to50(cutPropagator inCutProp);
  bool CheckRecoPtBinsR1Cent30to50(std::vector<double> inRecoPtBinsR1Cent30to50);
  bool CheckRecoPtBinsR1Cent30to50(cutPropagator inCutProp);

  bool CheckGenPtBinsR1Cent50to90(std::vector<double> inGenPtBinsR1Cent50to90);
  bool CheckGenPtBinsR1Cent50to90(cutPropagator inCutProp);
  bool CheckRecoPtBinsR1Cent50to90(std::vector<double> inRecoPtBinsR1Cent50to90);
  bool CheckRecoPtBinsR1Cent50to90(cutPropagator inCutProp);

  bool CheckGenPtBinsR2Cent0to10(std::vector<double> inGenPtBinsR2Cent0to10);
  bool CheckGenPtBinsR2Cent0to10(cutPropagator inCutProp);
  bool CheckRecoPtBinsR2Cent0to10(std::vector<double> inRecoPtBinsR2Cent0to10);
  bool CheckRecoPtBinsR2Cent0to10(cutPropagator inCutProp);

  bool CheckGenPtBinsR2Cent10to30(std::vector<double> inGenPtBinsR2Cent10to30);
  bool CheckGenPtBinsR2Cent10to30(cutPropagator inCutProp);
  bool CheckRecoPtBinsR2Cent10to30(std::vector<double> inRecoPtBinsR2Cent10to30);
  bool CheckRecoPtBinsR2Cent10to30(cutPropagator inCutProp);

  bool CheckGenPtBinsR2Cent30to50(std::vector<double> inGenPtBinsR2Cent30to50);
  bool CheckGenPtBinsR2Cent30to50(cutPropagator inCutProp);
  bool CheckRecoPtBinsR2Cent30to50(std::vector<double> inRecoPtBinsR2Cent30to50);
  bool CheckRecoPtBinsR2Cent30to50(cutPropagator inCutProp);

  bool CheckGenPtBinsR2Cent50to90(std::vector<double> inGenPtBinsR2Cent50to90);
  bool CheckGenPtBinsR2Cent50to90(cutPropagator inCutProp);
  bool CheckRecoPtBinsR2Cent50to90(std::vector<double> inRecoPtBinsR2Cent50to90);
  bool CheckRecoPtBinsR2Cent50to90(cutPropagator inCutProp);

  bool CheckGenPtBinsR3Cent0to10(std::vector<double> inGenPtBinsR3Cent0to10);
  bool CheckGenPtBinsR3Cent0to10(cutPropagator inCutProp);
  bool CheckRecoPtBinsR3Cent0to10(std::vector<double> inRecoPtBinsR3Cent0to10);
  bool CheckRecoPtBinsR3Cent0to10(cutPropagator inCutProp);

  bool CheckGenPtBinsR3Cent10to30(std::vector<double> inGenPtBinsR3Cent10to30);
  bool CheckGenPtBinsR3Cent10to30(cutPropagator inCutProp);
  bool CheckRecoPtBinsR3Cent10to30(std::vector<double> inRecoPtBinsR3Cent10to30);
  bool CheckRecoPtBinsR3Cent10to30(cutPropagator inCutProp);

  bool CheckGenPtBinsR3Cent30to50(std::vector<double> inGenPtBinsR3Cent30to50);
  bool CheckGenPtBinsR3Cent30to50(cutPropagator inCutProp);
  bool CheckRecoPtBinsR3Cent30to50(std::vector<double> inRecoPtBinsR3Cent30to50);
  bool CheckRecoPtBinsR3Cent30to50(cutPropagator inCutProp);

  bool CheckGenPtBinsR3Cent50to90(std::vector<double> inGenPtBinsR3Cent50to90);
  bool CheckGenPtBinsR3Cent50to90(cutPropagator inCutProp);
  bool CheckRecoPtBinsR3Cent50to90(std::vector<double> inRecoPtBinsR3Cent50to90);
  bool CheckRecoPtBinsR3Cent50to90(cutPropagator inCutProp);

  bool CheckGenPtBinsR4Cent0to10(std::vector<double> inGenPtBinsR4Cent0to10);
  bool CheckGenPtBinsR4Cent0to10(cutPropagator inCutProp);
  bool CheckRecoPtBinsR4Cent0to10(std::vector<double> inRecoPtBinsR4Cent0to10);
  bool CheckRecoPtBinsR4Cent0to10(cutPropagator inCutProp);

  bool CheckGenPtBinsR4Cent10to30(std::vector<double> inGenPtBinsR4Cent10to30);
  bool CheckGenPtBinsR4Cent10to30(cutPropagator inCutProp);
  bool CheckRecoPtBinsR4Cent10to30(std::vector<double> inRecoPtBinsR4Cent10to30);
  bool CheckRecoPtBinsR4Cent10to30(cutPropagator inCutProp);

  bool CheckGenPtBinsR4Cent30to50(std::vector<double> inGenPtBinsR4Cent30to50);
  bool CheckGenPtBinsR4Cent30to50(cutPropagator inCutProp);
  bool CheckRecoPtBinsR4Cent30to50(std::vector<double> inRecoPtBinsR4Cent30to50);
  bool CheckRecoPtBinsR4Cent30to50(cutPropagator inCutProp);

  bool CheckGenPtBinsR4Cent50to90(std::vector<double> inGenPtBinsR4Cent50to90);
  bool CheckGenPtBinsR4Cent50to90(cutPropagator inCutProp);
  bool CheckRecoPtBinsR4Cent50to90(std::vector<double> inRecoPtBinsR4Cent50to90);
  bool CheckRecoPtBinsR4Cent50to90(cutPropagator inCutProp);


  bool CheckGenPtBinsR5Cent0to10(std::vector<double> inGenPtBinsR5Cent0to10);
  bool CheckGenPtBinsR5Cent0to10(cutPropagator inCutProp);
  bool CheckRecoPtBinsR5Cent0to10(std::vector<double> inRecoPtBinsR5Cent0to10);
  bool CheckRecoPtBinsR5Cent0to10(cutPropagator inCutProp);

  bool CheckGenPtBinsR5Cent10to30(std::vector<double> inGenPtBinsR5Cent10to30);
  bool CheckGenPtBinsR5Cent10to30(cutPropagator inCutProp);
  bool CheckRecoPtBinsR5Cent10to30(std::vector<double> inRecoPtBinsR5Cent10to30);
  bool CheckRecoPtBinsR5Cent10to30(cutPropagator inCutProp);

  bool CheckGenPtBinsR5Cent30to50(std::vector<double> inGenPtBinsR5Cent30to50);
  bool CheckGenPtBinsR5Cent30to50(cutPropagator inCutProp);
  bool CheckRecoPtBinsR5Cent30to50(std::vector<double> inRecoPtBinsR5Cent30to50);
  bool CheckRecoPtBinsR5Cent30to50(cutPropagator inCutProp);

  bool CheckGenPtBinsR5Cent50to90(std::vector<double> inGenPtBinsR5Cent50to90);
  bool CheckGenPtBinsR5Cent50to90(cutPropagator inCutProp);
  bool CheckRecoPtBinsR5Cent50to90(std::vector<double> inRecoPtBinsR5Cent50to90);
  bool CheckRecoPtBinsR5Cent50to90(cutPropagator inCutProp);

  bool CheckGenPtBinsR6Cent0to10(std::vector<double> inGenPtBinsR6Cent0to10);
  bool CheckGenPtBinsR6Cent0to10(cutPropagator inCutProp);
  bool CheckRecoPtBinsR6Cent0to10(std::vector<double> inRecoPtBinsR6Cent0to10);
  bool CheckRecoPtBinsR6Cent0to10(cutPropagator inCutProp);

  bool CheckGenPtBinsR6Cent10to30(std::vector<double> inGenPtBinsR6Cent10to30);
  bool CheckGenPtBinsR6Cent10to30(cutPropagator inCutProp);
  bool CheckRecoPtBinsR6Cent10to30(std::vector<double> inRecoPtBinsR6Cent10to30);
  bool CheckRecoPtBinsR6Cent10to30(cutPropagator inCutProp);

  bool CheckGenPtBinsR6Cent30to50(std::vector<double> inGenPtBinsR6Cent30to50);
  bool CheckGenPtBinsR6Cent30to50(cutPropagator inCutProp);
  bool CheckRecoPtBinsR6Cent30to50(std::vector<double> inRecoPtBinsR6Cent30to50);
  bool CheckRecoPtBinsR6Cent30to50(cutPropagator inCutProp);

  bool CheckGenPtBinsR6Cent50to90(std::vector<double> inGenPtBinsR6Cent50to90);
  bool CheckGenPtBinsR6Cent50to90(cutPropagator inCutProp);
  bool CheckRecoPtBinsR6Cent50to90(std::vector<double> inRecoPtBinsR6Cent50to90);
  bool CheckRecoPtBinsR6Cent50to90(cutPropagator inCutProp);

  bool CheckGenPtBinsR8Cent0to10(std::vector<double> inGenPtBinsR8Cent0to10);
  bool CheckGenPtBinsR8Cent0to10(cutPropagator inCutProp);
  bool CheckRecoPtBinsR8Cent0to10(std::vector<double> inRecoPtBinsR8Cent0to10);
  bool CheckRecoPtBinsR8Cent0to10(cutPropagator inCutProp);

  bool CheckGenPtBinsR8Cent10to30(std::vector<double> inGenPtBinsR8Cent10to30);
  bool CheckGenPtBinsR8Cent10to30(cutPropagator inCutProp);
  bool CheckRecoPtBinsR8Cent10to30(std::vector<double> inRecoPtBinsR8Cent10to30);
  bool CheckRecoPtBinsR8Cent10to30(cutPropagator inCutProp);

  bool CheckGenPtBinsR8Cent30to50(std::vector<double> inGenPtBinsR8Cent30to50);
  bool CheckGenPtBinsR8Cent30to50(cutPropagator inCutProp);
  bool CheckRecoPtBinsR8Cent30to50(std::vector<double> inRecoPtBinsR8Cent30to50);
  bool CheckRecoPtBinsR8Cent30to50(cutPropagator inCutProp);

  bool CheckGenPtBinsR8Cent50to90(std::vector<double> inGenPtBinsR8Cent50to90);
  bool CheckGenPtBinsR8Cent50to90(cutPropagator inCutProp);
  bool CheckRecoPtBinsR8Cent50to90(std::vector<double> inRecoPtBinsR8Cent50to90);
  bool CheckRecoPtBinsR8Cent50to90(cutPropagator inCutProp);

  bool CheckGenPtBinsR10Cent0to10(std::vector<double> inGenPtBinsR10Cent0to10);
  bool CheckGenPtBinsR10Cent0to10(cutPropagator inCutProp);
  bool CheckRecoPtBinsR10Cent0to10(std::vector<double> inRecoPtBinsR10Cent0to10);
  bool CheckRecoPtBinsR10Cent0to10(cutPropagator inCutProp);

  bool CheckGenPtBinsR10Cent10to30(std::vector<double> inGenPtBinsR10Cent10to30);
  bool CheckGenPtBinsR10Cent10to30(cutPropagator inCutProp);
  bool CheckRecoPtBinsR10Cent10to30(std::vector<double> inRecoPtBinsR10Cent10to30);
  bool CheckRecoPtBinsR10Cent10to30(cutPropagator inCutProp);

  bool CheckGenPtBinsR10Cent30to50(std::vector<double> inGenPtBinsR10Cent30to50);
  bool CheckGenPtBinsR10Cent30to50(cutPropagator inCutProp);
  bool CheckRecoPtBinsR10Cent30to50(std::vector<double> inRecoPtBinsR10Cent30to50);
  bool CheckRecoPtBinsR10Cent30to50(cutPropagator inCutProp);

  bool CheckGenPtBinsR10Cent50to90(std::vector<double> inGenPtBinsR10Cent50to90);
  bool CheckGenPtBinsR10Cent50to90(cutPropagator inCutProp);
  bool CheckRecoPtBinsR10Cent50to90(std::vector<double> inRecoPtBinsR10Cent50to90);
  bool CheckRecoPtBinsR10Cent50to90(cutPropagator inCutProp);


  bool CheckNGeneralBins(int nGeneralBins);
  bool CheckNGeneralBins(cutPropagator inCutProp);

  bool CheckNJtAbsEtaBins(int inNJtAbsEtaBins);
  bool CheckNJtAbsEtaBins(cutPropagator inCutProp);
  bool CheckNID(int inNID);
  bool CheckNID(cutPropagator inCutProp);

  bool CheckGeneralBins(std::vector<double> generalBins);
  bool CheckGeneralBins(cutPropagator inCutProp);

  bool CheckJtAbsEtaBinsLow(std::vector<double> inJtAbsEtaBinsLow);
  bool CheckJtAbsEtaBinsLow(cutPropagator inCutProp);
  bool CheckJtAbsEtaBinsHi(std::vector<double> inJtAbsEtaBinsHi);
  bool CheckJtAbsEtaBinsHi(cutPropagator inCutProp);
  bool CheckIdStr(std::vector<std::string> inIDStr);
  bool CheckIdStr(cutPropagator inCutProp);
  bool CheckJtPfCHMFCutLow(std::vector<double> inJtPfCHMFCutLow);
  bool CheckJtPfCHMFCutLow(cutPropagator inCutProp);
  bool CheckJtPfCHMFCutHi(std::vector<double> inJtPfCHMFCutHi);
  bool CheckJtPfCHMFCutHi(cutPropagator inCutProp);
  bool CheckJtPfMUMFCutLow(std::vector<double> inJtPfMUMFCutLow);
  bool CheckJtPfMUMFCutLow(cutPropagator inCutProp);
  bool CheckJtPfMUMFCutHi(std::vector<double> inJtPfMUMFCutHi);
  bool CheckJtPfMUMFCutHi(cutPropagator inCutProp);
  bool CheckJtPfNHFCutLow(std::vector<double> inJtPfNHFCutLow);
  bool CheckJtPfNHFCutLow(cutPropagator inCutProp);
  bool CheckJtPfNHFCutHi(std::vector<double> inJtPfNHFCutHi);
  bool CheckJtPfNHFCutHi(cutPropagator inCutProp);
  bool CheckJtPfNEFCutLow(std::vector<double> inJtPfNEFCutLow);
  bool CheckJtPfNEFCutLow(cutPropagator inCutProp);
  bool CheckJtPfNEFCutHi(std::vector<double> inJtPfNEFCutHi);
  bool CheckJtPfNEFCutHi(cutPropagator inCutProp);
  bool CheckJtPfMUFCutLow(std::vector<double> inJtPfMUFCutLow);
  bool CheckJtPfMUFCutLow(cutPropagator inCutProp);
  bool CheckJtPfMUFCutHi(std::vector<double> inJtPfMUFCutHi);
  bool CheckJtPfMUFCutHi(cutPropagator inCutProp);
  bool CheckJtPfCHFCutLow(std::vector<double> inJtPfCHFCutLow);
  bool CheckJtPfCHFCutLow(cutPropagator inCutProp);
  bool CheckJtPfCHFCutHi(std::vector<double> inJtPfCHFCutHi);
  bool CheckJtPfCHFCutHi(cutPropagator inCutProp);
  bool CheckJtPfCEFCutLow(std::vector<double> inJtPfCEFCutLow);
  bool CheckJtPfCEFCutLow(cutPropagator inCutProp);
  bool CheckJtPfCEFCutHi(std::vector<double> inJtPfCEFCutHi);
  bool CheckJtPfCEFCutHi(cutPropagator inCutProp);
  bool CheckJtPfMinMult(std::vector<int> inJtPfMinMult);
  bool CheckJtPfMinMult(cutPropagator inCutProp);
  bool CheckJtPfMinChgMult(std::vector<int> inJtPfMinChgMult);
  bool CheckJtPfMinChgMult(cutPropagator inCutProp);
  bool CheckNSyst(int inNSyst);
  bool CheckNSyst(cutPropagator inCutProp);
  bool CheckSystStr(std::vector<std::string> inSystStr);
  bool CheckSystStr(cutPropagator inCutProp);
  bool CheckNBayes(int inNBayes);
  bool CheckNBayes(cutPropagator inCutProp);
  bool CheckNSVD(int inNSVD);
  bool CheckNSVD(cutPropagator inCutProp);
  bool CheckNBigBayesSymm(int inNBigBayesSymm);
  bool CheckNBigBayesSymm(cutPropagator inCutProp);
  bool CheckBayesVal(std::vector<int> inBayesVal);
  bool CheckBayesVal(cutPropagator inCutProp);
  bool CheckNSuperBayes(int inNSuperBayes);
  bool CheckNSuperBayes(cutPropagator inCutProp);
  bool CheckNHistDim(int inNHistDim);
  bool CheckNHistDim(cutPropagator inCutProp);
  bool CheckHistBestBayes(std::vector<int> inHistBestBayes);
  bool CheckHistBestBayes(cutPropagator inCutProp);
  bool CheckHistBestSVD(std::vector<int> inHistBestSVD);
  bool CheckHistBestSVD(cutPropagator inCutProp);
  bool CheckRCDiffFileName(std::string inRCDiffFileName);
  bool CheckRCDiffFileName(cutPropagator inCutProp);
  bool CheckPriorFlatFileName(std::string inPriorFlatFileName);
  bool CheckPriorFlatFileName(cutPropagator inCutProp);
  bool CheckJECVarMC(double inJECVarMC);
  bool CheckJECVarMC(cutPropagator inCutProp);
  bool CheckJERVarMC(double inJERVarMC);
  bool CheckJERVarMC(cutPropagator inCutProp);
  bool CheckJECVarData(double inJECVarData);
  bool CheckJECVarData(cutPropagator inCutProp);
  bool CheckNResponseMod(int inNResponseMod);
  bool CheckNResponseMod(cutPropagator inCutProp);
  bool CheckResponseMod(std::vector<double> inResponseMod);
  bool CheckResponseMod(cutPropagator inCutProp);
  bool CheckJERVarData(std::vector<double> inJERVarData);
  bool CheckJERVarData(cutPropagator inCutProp);
  bool CheckNPthats(int inNPthats);
  bool CheckNPthats(cutPropagator inCutProp);
  bool CheckPthats(std::vector<double> inPthats);
  bool CheckPthats(cutPropagator inCutProp);
  bool CheckPthatWeights(std::vector<double> inPthatWeights);
  bool CheckPthatWeights(cutPropagator inCutProp);
  bool CheckNJtAlgos(int inNJtAlgos);
  bool CheckNJtAlgos(cutPropagator inCutProp);
  bool CheckJtAlgos(std::vector<std::string> inJtAlgos);
  bool CheckJtAlgos(cutPropagator inCutProp);
  bool CheckMinJtPtCut(std::vector<double> inMinJtPtCut);
  bool CheckMinJtPtCut(cutPropagator inCutProp);
  bool CheckMultiJtPtCut(std::vector<double> inMultiJtPtCut);
  bool CheckMultiJtPtCut(cutPropagator inCutProp);
  bool CheckRecoTruncPos(std::vector<int> inRecoTruncPos);
  bool CheckRecoTruncPos(cutPropagator inCutProp);
  bool CheckIsPP(bool inIsPP);
  bool CheckIsPP(cutPropagator inCutProp);
  bool CheckNCentBins(int inNCentBins);
  bool CheckNCentBins(cutPropagator inCutProp);
  bool CheckCentBinsLow(std::vector<int> inCentBinsLow);
  bool CheckCentBinsLow(cutPropagator inCutProp);
  bool CheckCentBinsHi(std::vector<int> inCentBinsHi);
  bool CheckCentBinsHi(cutPropagator inCutProp);

  std::vector<std::string> GetInFileNames(){return inFileNames;}
  std::vector<std::string> GetInFullFileNames(){return inFullFileNames;}

  bool GetIsPP(){return isPP;}

  double GetJtAbsEtaMax(){return jtAbsEtaMax;}

  std::string GetRCDiffFileName(){return rcDiffFileName;}
  std::string GetPriorFlatFileName(){return flatPriorFileName;}

  double GetJECVarMC(){return jecVarMC;}
  double GetJERVarMC(){return jerVarMC;}

  double GetJECVarData(){return jecVarData;}

  int GetNResponseMod(){return nResponseMod;}
  std::vector<double> GetResponseMod(){return responseMod;}
  std::vector<double> GetJERVarData(){return jerVarData;}

  int GetNJtAlgos(){return nJtAlgos;}
  std::vector<std::string> GetJtAlgos(){return jtAlgos;}
  std::vector<double> GetMinJtPtCut(){return minJtPtCut;}
  std::vector<double> GetMultiJtPtCut(){return multiJtPtCut;}
  std::vector<int> GetRecoTruncPos(){return recoTruncPos;}

  int GetNR(){return nR;}
  std::vector<int> GetRVals(){return rVals;}

  std::vector<double> GetGenPtBinsR1Cent0to10(){return genJtPtBinsR1Cent0to10;}
  std::vector<double> GetRecoPtBinsR1Cent0to10(){return recoJtPtBinsR1Cent0to10;}
  std::vector<double> GetGenPtBinsR1Cent10to30(){return genJtPtBinsR1Cent10to30;}
  std::vector<double> GetRecoPtBinsR1Cent10to30(){return recoJtPtBinsR1Cent10to30;}
  std::vector<double> GetGenPtBinsR1Cent30to50(){return genJtPtBinsR1Cent30to50;}
  std::vector<double> GetRecoPtBinsR1Cent30to50(){return recoJtPtBinsR1Cent30to50;}
  std::vector<double> GetGenPtBinsR1Cent50to90(){return genJtPtBinsR1Cent50to90;}
  std::vector<double> GetRecoPtBinsR1Cent50to90(){return recoJtPtBinsR1Cent50to90;}

  std::vector<double> GetGenPtBinsR2Cent0to10(){return genJtPtBinsR2Cent0to10;}
  std::vector<double> GetRecoPtBinsR2Cent0to10(){return recoJtPtBinsR2Cent0to10;}
  std::vector<double> GetGenPtBinsR2Cent10to30(){return genJtPtBinsR2Cent10to30;}
  std::vector<double> GetRecoPtBinsR2Cent10to30(){return recoJtPtBinsR2Cent10to30;}
  std::vector<double> GetGenPtBinsR2Cent30to50(){return genJtPtBinsR2Cent30to50;}
  std::vector<double> GetRecoPtBinsR2Cent30to50(){return recoJtPtBinsR2Cent30to50;}
  std::vector<double> GetGenPtBinsR2Cent50to90(){return genJtPtBinsR2Cent50to90;}
  std::vector<double> GetRecoPtBinsR2Cent50to90(){return recoJtPtBinsR2Cent50to90;}

  std::vector<double> GetGenPtBinsR3Cent0to10(){return genJtPtBinsR3Cent0to10;}
  std::vector<double> GetRecoPtBinsR3Cent0to10(){return recoJtPtBinsR3Cent0to10;}
  std::vector<double> GetGenPtBinsR3Cent10to30(){return genJtPtBinsR3Cent10to30;}
  std::vector<double> GetRecoPtBinsR3Cent10to30(){return recoJtPtBinsR3Cent10to30;}
  std::vector<double> GetGenPtBinsR3Cent30to50(){return genJtPtBinsR3Cent30to50;}
  std::vector<double> GetRecoPtBinsR3Cent30to50(){return recoJtPtBinsR3Cent30to50;}
  std::vector<double> GetGenPtBinsR3Cent50to90(){return genJtPtBinsR3Cent50to90;}
  std::vector<double> GetRecoPtBinsR3Cent50to90(){return recoJtPtBinsR3Cent50to90;}

  std::vector<double> GetGenPtBinsR4Cent0to10(){return genJtPtBinsR4Cent0to10;}
  std::vector<double> GetRecoPtBinsR4Cent0to10(){return recoJtPtBinsR4Cent0to10;}
  std::vector<double> GetGenPtBinsR4Cent10to30(){return genJtPtBinsR4Cent10to30;}
  std::vector<double> GetRecoPtBinsR4Cent10to30(){return recoJtPtBinsR4Cent10to30;}
  std::vector<double> GetGenPtBinsR4Cent30to50(){return genJtPtBinsR4Cent30to50;}
  std::vector<double> GetRecoPtBinsR4Cent30to50(){return recoJtPtBinsR4Cent30to50;}
  std::vector<double> GetGenPtBinsR4Cent50to90(){return genJtPtBinsR4Cent50to90;}
  std::vector<double> GetRecoPtBinsR4Cent50to90(){return recoJtPtBinsR4Cent50to90;}

  std::vector<double> GetGenPtBinsR5Cent0to10(){return genJtPtBinsR5Cent0to10;}
  std::vector<double> GetRecoPtBinsR5Cent0to10(){return recoJtPtBinsR5Cent0to10;}
  std::vector<double> GetGenPtBinsR5Cent10to30(){return genJtPtBinsR5Cent10to30;}
  std::vector<double> GetRecoPtBinsR5Cent10to30(){return recoJtPtBinsR5Cent10to30;}
  std::vector<double> GetGenPtBinsR5Cent30to50(){return genJtPtBinsR5Cent30to50;}
  std::vector<double> GetRecoPtBinsR5Cent30to50(){return recoJtPtBinsR5Cent30to50;}
  std::vector<double> GetGenPtBinsR5Cent50to90(){return genJtPtBinsR5Cent50to90;}
  std::vector<double> GetRecoPtBinsR5Cent50to90(){return recoJtPtBinsR5Cent50to90;}

  std::vector<double> GetGenPtBinsR6Cent0to10(){return genJtPtBinsR6Cent0to10;}
  std::vector<double> GetRecoPtBinsR6Cent0to10(){return recoJtPtBinsR6Cent0to10;}
  std::vector<double> GetGenPtBinsR6Cent10to30(){return genJtPtBinsR6Cent10to30;}
  std::vector<double> GetRecoPtBinsR6Cent10to30(){return recoJtPtBinsR6Cent10to30;}
  std::vector<double> GetGenPtBinsR6Cent30to50(){return genJtPtBinsR6Cent30to50;}
  std::vector<double> GetRecoPtBinsR6Cent30to50(){return recoJtPtBinsR6Cent30to50;}
  std::vector<double> GetGenPtBinsR6Cent50to90(){return genJtPtBinsR6Cent50to90;}
  std::vector<double> GetRecoPtBinsR6Cent50to90(){return recoJtPtBinsR6Cent50to90;}

  std::vector<double> GetGenPtBinsR8Cent0to10(){return genJtPtBinsR8Cent0to10;}
  std::vector<double> GetRecoPtBinsR8Cent0to10(){return recoJtPtBinsR8Cent0to10;}
  std::vector<double> GetGenPtBinsR8Cent10to30(){return genJtPtBinsR8Cent10to30;}
  std::vector<double> GetRecoPtBinsR8Cent10to30(){return recoJtPtBinsR8Cent10to30;}
  std::vector<double> GetGenPtBinsR8Cent30to50(){return genJtPtBinsR8Cent30to50;}
  std::vector<double> GetRecoPtBinsR8Cent30to50(){return recoJtPtBinsR8Cent30to50;}
  std::vector<double> GetGenPtBinsR8Cent50to90(){return genJtPtBinsR8Cent50to90;}
  std::vector<double> GetRecoPtBinsR8Cent50to90(){return recoJtPtBinsR8Cent50to90;}

  std::vector<double> GetGenPtBinsR10Cent0to10(){return genJtPtBinsR10Cent0to10;}
  std::vector<double> GetRecoPtBinsR10Cent0to10(){return recoJtPtBinsR10Cent0to10;}
  std::vector<double> GetGenPtBinsR10Cent10to30(){return genJtPtBinsR10Cent10to30;}
  std::vector<double> GetRecoPtBinsR10Cent10to30(){return recoJtPtBinsR10Cent10to30;}
  std::vector<double> GetGenPtBinsR10Cent30to50(){return genJtPtBinsR10Cent30to50;}
  std::vector<double> GetRecoPtBinsR10Cent30to50(){return recoJtPtBinsR10Cent30to50;}
  std::vector<double> GetGenPtBinsR10Cent50to90(){return genJtPtBinsR10Cent50to90;}
  std::vector<double> GetRecoPtBinsR10Cent50to90(){return recoJtPtBinsR10Cent50to90;}


  std::vector<double> GetGenPtBins(int rVal, std::string centStr);
  std::vector<double> GetRecoPtBins(int rVal, std::string centStr);

  int GetGenNBinsFromRValCent(int rVal, std::string centStr);
  void GetGenBinsFromRValCent(int rVal, std::string centStr, double outBins[]);
  int GetRecoNBinsFromRValCent(int rVal, std::string centStr);
  void GetRecoBinsFromRValCent(int rVal, std::string centStr, double outBins[]);

  int GetNGeneralBins(){return nGeneralBins;}
  std::vector<double> GetGeneralBins(){return generalBins;}

  int GetNJtAbsEtaBins(){return nJtAbsEtaBins;}
  std::vector<double> GetJtAbsEtaBinsLow(){return jtAbsEtaBinsLow;}
  std::vector<double> GetJtAbsEtaBinsHi(){return jtAbsEtaBinsHi;}

  int GetNPthats(){return nPthats;}
  std::vector<double> GetPthats(){return pthats;}
  std::vector<double> GetPthatWeights(){return pthatWeights;}

  int GetNCentBins(){return nCentBins;}
  std::vector<int> GetCentBinsLow(){return centBinsLow;}
  std::vector<int> GetCentBinsHi(){return centBinsHi;}

  int GetNID(){return nID;}
  std::vector<std::string> GetIdStr(){return idStr;}
  std::vector<double> GetJtPfCHMFCutLow(){return jtPfCHMFCutLow;}
  std::vector<double> GetJtPfCHMFCutHi(){return jtPfCHMFCutHi;}
  std::vector<double> GetJtPfMUMFCutLow(){return jtPfMUMFCutLow;}
  std::vector<double> GetJtPfMUMFCutHi(){return jtPfMUMFCutHi;}
  std::vector<double> GetJtPfNHFCutLow(){return jtPfNHFCutLow;}
  std::vector<double> GetJtPfNHFCutHi(){return jtPfNHFCutHi;}
  std::vector<double> GetJtPfNEFCutLow(){return jtPfNEFCutLow;}
  std::vector<double> GetJtPfNEFCutHi(){return jtPfNEFCutHi;}
  std::vector<double> GetJtPfMUFCutLow(){return jtPfMUFCutLow;}
  std::vector<double> GetJtPfMUFCutHi(){return jtPfMUFCutHi;}
  std::vector<double> GetJtPfCHFCutLow(){return jtPfCHFCutLow;}
  std::vector<double> GetJtPfCHFCutHi(){return jtPfCHFCutHi;}
  std::vector<double> GetJtPfCEFCutLow(){return jtPfCEFCutLow;}
  std::vector<double> GetJtPfCEFCutHi(){return jtPfCEFCutHi;}
  std::vector<int> GetJtPfMinMult(){return jtPfMinMult;}
  std::vector<int> GetJtPfMinChgMult(){return jtPfMinChgMult;}

  int GetNSyst(){return nSyst;}
  std::vector<std::string> GetSystStr(){return systStr;}

  int GetNBayes(){return nBayes;}
  int GetNBigBayesSymm(){return nBigBayesSymm;}
  std::vector<int> GetBayesVal(){return bayesVal;}

  int GetNSuperBayes(){return nSuperBayes;}

  int GetNSVD(){return nSVD;}

  int GetNHistDim(){return nHistDim;}
  std::vector<std::string> GetHistTagBayes(){return histTagBayes;}
  std::vector<std::string> GetHistTagSVD(){return histTagSVD;}
  std::vector<int> GetHistBestBayes(){return histBestBayes;}
  std::vector<int> GetHistBestSVD(){return histBestSVD;}

  void SetInFileNames(std::vector<std::string> inInFileNames){inFileNames = inInFileNames; return;}
  void SetInFullFileNames(std::vector<std::string> inInFullFileNames){inFullFileNames = inInFullFileNames; return;}

  void SetIsPP(bool inIsPP){isPP = inIsPP; return;}

  void SetJtAbsEtaMax(double inJtAbsEtaMax){jtAbsEtaMax = inJtAbsEtaMax; return;}

  void SetRCDiffFileName(std::string inRCDiffFileName){rcDiffFileName = inRCDiffFileName; return;}
  void SetPriorFlatFileName(std::string inPriorFlatFileName){flatPriorFileName = inPriorFlatFileName; return;}
  
  void SetJECVarMC(double inJECVarMC){jecVarMC = inJECVarMC; return;}
  void SetJERVarMC(double inJERVarMC){jerVarMC = inJERVarMC; return;}

  void SetJECVarData(double inJECVarData){jecVarData = inJECVarData; return;}

  void SetNResponseMod(int inNResponseMod){nResponseMod = inNResponseMod; return;}
  void SetResponseMod(std::vector<double> inResponseMod){responseMod = inResponseMod; return;}
  void SetResponseMod(int inN, const Double_t inResponseMod[]);
  void SetJERVarData(std::vector<double> inJERVarData){jerVarData = inJERVarData; return;}
  void SetJERVarData(int inN, const Double_t inJERVarData[]);

  void SetNJtAlgos(int inNJtAlgos){nJtAlgos = inNJtAlgos; return;}
  void SetJtAlgos(std::vector<std::string> inJtAlgos){jtAlgos = inJtAlgos; return;}
  void SetJtAlgos(int inN, std::string inJtAlgos[]);
  void SetMinJtPtCut(std::vector<double> inMinJtPtCut){minJtPtCut = inMinJtPtCut; return;}
  void SetMinJtPtCut(int inN, double inMinJtPtCut[]);
  void SetMultiJtPtCut(std::vector<double> inMultiJtPtCut){multiJtPtCut = inMultiJtPtCut; return;}
  void SetMultiJtPtCut(int inN, double inMultiJtPtCut[]);
  void SetRecoTruncPos(std::vector<int> inRecoTruncPos){recoTruncPos = inRecoTruncPos; return;}
  void SetRecoTruncPos(int inN, int inRecoTruncPos[]);


  void SetNR(int inNR){nR = inNR; return;}
  void SetRVals(std::vector<int> inRVals){rVals = inRVals; return;}
  void SetRVals(int inN, const Int_t inRVals[]);


  void SetGenPtBinsR1Cent0to10(std::vector<double> inGenPtBinsR1Cent0to10){genJtPtBinsR1Cent0to10 = inGenPtBinsR1Cent0to10; return;}
  void SetRecoPtBinsR1Cent0to10(std::vector<double> inRecoPtBinsR1Cent0to10){recoJtPtBinsR1Cent0to10 = inRecoPtBinsR1Cent0to10; return;}
  void SetGenPtBinsR1Cent10to30(std::vector<double> inGenPtBinsR1Cent10to30){genJtPtBinsR1Cent10to30 = inGenPtBinsR1Cent10to30; return;}
  void SetRecoPtBinsR1Cent10to30(std::vector<double> inRecoPtBinsR1Cent10to30){recoJtPtBinsR1Cent10to30 = inRecoPtBinsR1Cent10to30; return;}
  void SetGenPtBinsR1Cent30to50(std::vector<double> inGenPtBinsR1Cent30to50){genJtPtBinsR1Cent30to50 = inGenPtBinsR1Cent30to50; return;}
  void SetRecoPtBinsR1Cent30to50(std::vector<double> inRecoPtBinsR1Cent30to50){recoJtPtBinsR1Cent30to50 = inRecoPtBinsR1Cent30to50; return;}
  void SetGenPtBinsR1Cent50to90(std::vector<double> inGenPtBinsR1Cent50to90){genJtPtBinsR1Cent50to90 = inGenPtBinsR1Cent50to90; return;}
  void SetRecoPtBinsR1Cent50to90(std::vector<double> inRecoPtBinsR1Cent50to90){recoJtPtBinsR1Cent50to90 = inRecoPtBinsR1Cent50to90; return;}

  void SetGenPtBinsR2Cent0to10(std::vector<double> inGenPtBinsR2Cent0to10){genJtPtBinsR2Cent0to10 = inGenPtBinsR2Cent0to10; return;}
  void SetRecoPtBinsR2Cent0to10(std::vector<double> inRecoPtBinsR2Cent0to10){recoJtPtBinsR2Cent0to10 = inRecoPtBinsR2Cent0to10; return;}
  void SetGenPtBinsR2Cent10to30(std::vector<double> inGenPtBinsR2Cent10to30){genJtPtBinsR2Cent10to30 = inGenPtBinsR2Cent10to30; return;}
  void SetRecoPtBinsR2Cent10to30(std::vector<double> inRecoPtBinsR2Cent10to30){recoJtPtBinsR2Cent10to30 = inRecoPtBinsR2Cent10to30; return;}
  void SetGenPtBinsR2Cent30to50(std::vector<double> inGenPtBinsR2Cent30to50){genJtPtBinsR2Cent30to50 = inGenPtBinsR2Cent30to50; return;}
  void SetRecoPtBinsR2Cent30to50(std::vector<double> inRecoPtBinsR2Cent30to50){recoJtPtBinsR2Cent30to50 = inRecoPtBinsR2Cent30to50; return;}
  void SetGenPtBinsR2Cent50to90(std::vector<double> inGenPtBinsR2Cent50to90){genJtPtBinsR2Cent50to90 = inGenPtBinsR2Cent50to90; return;}
  void SetRecoPtBinsR2Cent50to90(std::vector<double> inRecoPtBinsR2Cent50to90){recoJtPtBinsR2Cent50to90 = inRecoPtBinsR2Cent50to90; return;}

  void SetGenPtBinsR3Cent0to10(std::vector<double> inGenPtBinsR3Cent0to10){genJtPtBinsR3Cent0to10 = inGenPtBinsR3Cent0to10; return;}
  void SetRecoPtBinsR3Cent0to10(std::vector<double> inRecoPtBinsR3Cent0to10){recoJtPtBinsR3Cent0to10 = inRecoPtBinsR3Cent0to10; return;}
  void SetGenPtBinsR3Cent10to30(std::vector<double> inGenPtBinsR3Cent10to30){genJtPtBinsR3Cent10to30 = inGenPtBinsR3Cent10to30; return;}
  void SetRecoPtBinsR3Cent10to30(std::vector<double> inRecoPtBinsR3Cent10to30){recoJtPtBinsR3Cent10to30 = inRecoPtBinsR3Cent10to30; return;}
  void SetGenPtBinsR3Cent30to50(std::vector<double> inGenPtBinsR3Cent30to50){genJtPtBinsR3Cent30to50 = inGenPtBinsR3Cent30to50; return;}
  void SetRecoPtBinsR3Cent30to50(std::vector<double> inRecoPtBinsR3Cent30to50){recoJtPtBinsR3Cent30to50 = inRecoPtBinsR3Cent30to50; return;}
  void SetGenPtBinsR3Cent50to90(std::vector<double> inGenPtBinsR3Cent50to90){genJtPtBinsR3Cent50to90 = inGenPtBinsR3Cent50to90; return;}
  void SetRecoPtBinsR3Cent50to90(std::vector<double> inRecoPtBinsR3Cent50to90){recoJtPtBinsR3Cent50to90 = inRecoPtBinsR3Cent50to90; return;}

  void SetGenPtBinsR4Cent0to10(std::vector<double> inGenPtBinsR4Cent0to10){genJtPtBinsR4Cent0to10 = inGenPtBinsR4Cent0to10; return;}
  void SetRecoPtBinsR4Cent0to10(std::vector<double> inRecoPtBinsR4Cent0to10){recoJtPtBinsR4Cent0to10 = inRecoPtBinsR4Cent0to10; return;}
  void SetGenPtBinsR4Cent10to30(std::vector<double> inGenPtBinsR4Cent10to30){genJtPtBinsR4Cent10to30 = inGenPtBinsR4Cent10to30; return;}
  void SetRecoPtBinsR4Cent10to30(std::vector<double> inRecoPtBinsR4Cent10to30){recoJtPtBinsR4Cent10to30 = inRecoPtBinsR4Cent10to30; return;}
  void SetGenPtBinsR4Cent30to50(std::vector<double> inGenPtBinsR4Cent30to50){genJtPtBinsR4Cent30to50 = inGenPtBinsR4Cent30to50; return;}
  void SetRecoPtBinsR4Cent30to50(std::vector<double> inRecoPtBinsR4Cent30to50){recoJtPtBinsR4Cent30to50 = inRecoPtBinsR4Cent30to50; return;}
  void SetGenPtBinsR4Cent50to90(std::vector<double> inGenPtBinsR4Cent50to90){genJtPtBinsR4Cent50to90 = inGenPtBinsR4Cent50to90; return;}
  void SetRecoPtBinsR4Cent50to90(std::vector<double> inRecoPtBinsR4Cent50to90){recoJtPtBinsR4Cent50to90 = inRecoPtBinsR4Cent50to90; return;}

  void SetGenPtBinsR5Cent0to10(std::vector<double> inGenPtBinsR5Cent0to10){genJtPtBinsR5Cent0to10 = inGenPtBinsR5Cent0to10; return;}
  void SetRecoPtBinsR5Cent0to10(std::vector<double> inRecoPtBinsR5Cent0to10){recoJtPtBinsR5Cent0to10 = inRecoPtBinsR5Cent0to10; return;}
  void SetGenPtBinsR5Cent10to30(std::vector<double> inGenPtBinsR5Cent10to30){genJtPtBinsR5Cent10to30 = inGenPtBinsR5Cent10to30; return;}
  void SetRecoPtBinsR5Cent10to30(std::vector<double> inRecoPtBinsR5Cent10to30){recoJtPtBinsR5Cent10to30 = inRecoPtBinsR5Cent10to30; return;}
  void SetGenPtBinsR5Cent30to50(std::vector<double> inGenPtBinsR5Cent30to50){genJtPtBinsR5Cent30to50 = inGenPtBinsR5Cent30to50; return;}
  void SetRecoPtBinsR5Cent30to50(std::vector<double> inRecoPtBinsR5Cent30to50){recoJtPtBinsR5Cent30to50 = inRecoPtBinsR5Cent30to50; return;}
  void SetGenPtBinsR5Cent50to90(std::vector<double> inGenPtBinsR5Cent50to90){genJtPtBinsR5Cent50to90 = inGenPtBinsR5Cent50to90; return;}
  void SetRecoPtBinsR5Cent50to90(std::vector<double> inRecoPtBinsR5Cent50to90){recoJtPtBinsR5Cent50to90 = inRecoPtBinsR5Cent50to90; return;}

  void SetGenPtBinsR6Cent0to10(std::vector<double> inGenPtBinsR6Cent0to10){genJtPtBinsR6Cent0to10 = inGenPtBinsR6Cent0to10; return;}
  void SetRecoPtBinsR6Cent0to10(std::vector<double> inRecoPtBinsR6Cent0to10){recoJtPtBinsR6Cent0to10 = inRecoPtBinsR6Cent0to10; return;}
  void SetGenPtBinsR6Cent10to30(std::vector<double> inGenPtBinsR6Cent10to30){genJtPtBinsR6Cent10to30 = inGenPtBinsR6Cent10to30; return;}
  void SetRecoPtBinsR6Cent10to30(std::vector<double> inRecoPtBinsR6Cent10to30){recoJtPtBinsR6Cent10to30 = inRecoPtBinsR6Cent10to30; return;}
  void SetGenPtBinsR6Cent30to50(std::vector<double> inGenPtBinsR6Cent30to50){genJtPtBinsR6Cent30to50 = inGenPtBinsR6Cent30to50; return;}
  void SetRecoPtBinsR6Cent30to50(std::vector<double> inRecoPtBinsR6Cent30to50){recoJtPtBinsR6Cent30to50 = inRecoPtBinsR6Cent30to50; return;}
  void SetGenPtBinsR6Cent50to90(std::vector<double> inGenPtBinsR6Cent50to90){genJtPtBinsR6Cent50to90 = inGenPtBinsR6Cent50to90; return;}
  void SetRecoPtBinsR6Cent50to90(std::vector<double> inRecoPtBinsR6Cent50to90){recoJtPtBinsR6Cent50to90 = inRecoPtBinsR6Cent50to90; return;}

  void SetGenPtBinsR8Cent0to10(std::vector<double> inGenPtBinsR8Cent0to10){genJtPtBinsR8Cent0to10 = inGenPtBinsR8Cent0to10; return;}
  void SetRecoPtBinsR8Cent0to10(std::vector<double> inRecoPtBinsR8Cent0to10){recoJtPtBinsR8Cent0to10 = inRecoPtBinsR8Cent0to10; return;}
  void SetGenPtBinsR8Cent10to30(std::vector<double> inGenPtBinsR8Cent10to30){genJtPtBinsR8Cent10to30 = inGenPtBinsR8Cent10to30; return;}
  void SetRecoPtBinsR8Cent10to30(std::vector<double> inRecoPtBinsR8Cent10to30){recoJtPtBinsR8Cent10to30 = inRecoPtBinsR8Cent10to30; return;}
  void SetGenPtBinsR8Cent30to50(std::vector<double> inGenPtBinsR8Cent30to50){genJtPtBinsR8Cent30to50 = inGenPtBinsR8Cent30to50; return;}
  void SetRecoPtBinsR8Cent30to50(std::vector<double> inRecoPtBinsR8Cent30to50){recoJtPtBinsR8Cent30to50 = inRecoPtBinsR8Cent30to50; return;}
  void SetGenPtBinsR8Cent50to90(std::vector<double> inGenPtBinsR8Cent50to90){genJtPtBinsR8Cent50to90 = inGenPtBinsR8Cent50to90; return;}
  void SetRecoPtBinsR8Cent50to90(std::vector<double> inRecoPtBinsR8Cent50to90){recoJtPtBinsR8Cent50to90 = inRecoPtBinsR8Cent50to90; return;}

  void SetGenPtBinsR10Cent0to10(std::vector<double> inGenPtBinsR10Cent0to10){genJtPtBinsR10Cent0to10 = inGenPtBinsR10Cent0to10; return;}
  void SetRecoPtBinsR10Cent0to10(std::vector<double> inRecoPtBinsR10Cent0to10){recoJtPtBinsR10Cent0to10 = inRecoPtBinsR10Cent0to10; return;}
  void SetGenPtBinsR10Cent10to30(std::vector<double> inGenPtBinsR10Cent10to30){genJtPtBinsR10Cent10to30 = inGenPtBinsR10Cent10to30; return;}
  void SetRecoPtBinsR10Cent10to30(std::vector<double> inRecoPtBinsR10Cent10to30){recoJtPtBinsR10Cent10to30 = inRecoPtBinsR10Cent10to30; return;}
  void SetGenPtBinsR10Cent30to50(std::vector<double> inGenPtBinsR10Cent30to50){genJtPtBinsR10Cent30to50 = inGenPtBinsR10Cent30to50; return;}
  void SetRecoPtBinsR10Cent30to50(std::vector<double> inRecoPtBinsR10Cent30to50){recoJtPtBinsR10Cent30to50 = inRecoPtBinsR10Cent30to50; return;}
  void SetGenPtBinsR10Cent50to90(std::vector<double> inGenPtBinsR10Cent50to90){genJtPtBinsR10Cent50to90 = inGenPtBinsR10Cent50to90; return;}
  void SetRecoPtBinsR10Cent50to90(std::vector<double> inRecoPtBinsR10Cent50to90){recoJtPtBinsR10Cent50to90 = inRecoPtBinsR10Cent50to90; return;}


  void SetGenPtBins(int rVal, std::string centStr, std::vector<double> inGenPtBins);
  void SetRecoPtBins(int rVal, std::string centStr, std::vector<double> inRecoPtBins);
  void SetGenPtBins(int rVal, std::string centStr, int inNGenPtBins, double inGenPtBins[]);
  void SetRecoPtBins(int rVal, std::string centStr, int inNRecoPtBins, double inGenPtBins[]);

  void SetNGeneralBins(int inNGeneralBins){nGeneralBins = inNGeneralBins; return;}

  void SetGeneralBins(std::vector<double> inGeneralBins){generalBins = inGeneralBins; return;}
  void SetGeneralBins(int inN, const Double_t inGeneralBins[]);


  void SetNJtAbsEtaBins(int inNJtAbsEtaBins){nJtAbsEtaBins = inNJtAbsEtaBins; return;}
  void SetJtAbsEtaBinsLow(std::vector<double> inJtAbsEtaBinsLow){jtAbsEtaBinsLow = inJtAbsEtaBinsLow; return;}
  void SetJtAbsEtaBinsLow(int inN, const Double_t inJtAbsEtaBinsLow[]);  
  void SetJtAbsEtaBinsHi(std::vector<double> inJtAbsEtaBinsHi){jtAbsEtaBinsHi = inJtAbsEtaBinsHi; return;}
  void SetJtAbsEtaBinsHi(int inN, const Double_t inJtAbsEtaBinsHi[]);  

  void SetNPthats(int inNPthats){nPthats = inNPthats; return;}
  void SetPthats(std::vector<double> inPthats){pthats = inPthats; return;}
  void SetPthats(int inN, const Double_t inPthats[]);
  void SetPthatWeights(std::vector<double> inPthatWeights){pthatWeights = inPthatWeights; return;}
  void SetPthatWeights(int inN, const Double_t inPthatWeights[]);

  void SetNCentBins(int inNCentBins){nCentBins = inNCentBins; return;}
  void SetCentBinsLow(std::vector<int> inCentBinsLow){centBinsLow = inCentBinsLow; return;}
  void SetCentBinsLow(int inN, const Int_t inCentBinsLow[]);
  void SetCentBinsHi(std::vector<int> inCentBinsHi){centBinsHi = inCentBinsHi; return;}
  void SetCentBinsHi(int inN, const Int_t inCentBinsHi[]);

  void SetNID(int inNID){nID = inNID; return;}
  void SetIdStr(std::vector<std::string> inIdStr){idStr = inIdStr; return;};
  void SetIdStr(int inN, const std::string inIdStr[]);
  void SetJtPfCHMFCutLow(std::vector<double> inJtPfCHMFCutLow){jtPfCHMFCutLow = inJtPfCHMFCutLow; return;}
  void SetJtPfCHMFCutLow(int inN, const Double_t inJtPfCHMFCutLow[]);
  void SetJtPfCHMFCutHi(std::vector<double> inJtPfCHMFCutHi){jtPfCHMFCutHi = inJtPfCHMFCutHi; return;}
  void SetJtPfCHMFCutHi(int inN, const Double_t inJtPfCHMFCutHi[]);
  void SetJtPfMUMFCutLow(std::vector<double> inJtPfMUMFCutLow){jtPfMUMFCutLow = inJtPfMUMFCutLow; return;}
  void SetJtPfMUMFCutLow(int inN, const Double_t inJtPfMUMFCutLow[]);
  void SetJtPfMUMFCutHi(std::vector<double> inJtPfMUMFCutHi){jtPfMUMFCutHi = inJtPfMUMFCutHi; return;}
  void SetJtPfMUMFCutHi(int inN, const Double_t inJtPfMUMFCutHi[]);
  void SetJtPfNHFCutLow(std::vector<double> inJtPfNHFCutLow){jtPfNHFCutLow = inJtPfNHFCutLow; return;}
  void SetJtPfNHFCutLow(int inN, const Double_t inJtPfNHFCutLow[]);
  void SetJtPfNHFCutHi(std::vector<double> inJtPfNHFCutHi){jtPfNHFCutHi = inJtPfNHFCutHi; return;}
  void SetJtPfNHFCutHi(int inN, const Double_t inJtPfNHFCutHi[]);
  void SetJtPfNEFCutLow(std::vector<double> inJtPfNEFCutLow){jtPfNEFCutLow = inJtPfNEFCutLow; return;}
  void SetJtPfNEFCutLow(int inN, const Double_t inJtPfNEFCutLow[]);
  void SetJtPfNEFCutHi(std::vector<double> inJtPfNEFCutHi){jtPfNEFCutHi = inJtPfNEFCutHi; return;}
  void SetJtPfNEFCutHi(int inN, const Double_t inJtPfNEFCutHi[]);
  void SetJtPfMUFCutLow(std::vector<double> inJtPfMUFCutLow){jtPfMUFCutLow = inJtPfMUFCutLow; return;}
  void SetJtPfMUFCutLow(int inN, const Double_t inJtPfMUFCutLow[]);
  void SetJtPfMUFCutHi(std::vector<double> inJtPfMUFCutHi){jtPfMUFCutHi = inJtPfMUFCutHi; return;}
  void SetJtPfMUFCutHi(int inN, const Double_t inJtPfMUFCutHi[]);
  void SetJtPfCHFCutLow(std::vector<double> inJtPfCHFCutLow){jtPfCHFCutLow = inJtPfCHFCutLow; return;}
  void SetJtPfCHFCutLow(int inN, const Double_t inJtPfCHFCutLow[]);
  void SetJtPfCHFCutHi(std::vector<double> inJtPfCHFCutHi){jtPfCHFCutHi = inJtPfCHFCutHi; return;}
  void SetJtPfCHFCutHi(int inN, const Double_t inJtPfCHFCutHi[]);
  void SetJtPfCEFCutLow(std::vector<double> inJtPfCEFCutLow){jtPfCEFCutLow = inJtPfCEFCutLow; return;}
  void SetJtPfCEFCutLow(int inN, const Double_t inJtPfCEFCutLow[]);
  void SetJtPfCEFCutHi(std::vector<double> inJtPfCEFCutHi){jtPfCEFCutHi = inJtPfCEFCutHi; return;}
  void SetJtPfCEFCutHi(int inN, const Double_t inJtPfCEFCutHi[]);
  void SetJtPfMinMult(std::vector<int> inJtPfMinMult){jtPfMinMult = inJtPfMinMult; return;}
  void SetJtPfMinMult(int inN, const Int_t inJtPfMinMult[]);
  void SetJtPfMinChgMult(std::vector<int> inJtPfMinChgMult){jtPfMinChgMult = inJtPfMinChgMult; return;}
  void SetJtPfMinChgMult(int inN, const Int_t inJtPfMinChgMult[]);

  void SetNSyst(int inNSyst){nSyst = inNSyst; return;}
  void SetSystStr(std::vector<std::string> inSystStr){systStr = inSystStr; return;};
  void SetSystStr(int inN, const std::string inSystStr[]);

  void SetNBayes(int inNBayes){nBayes = inNBayes; return;}
  void SetNBigBayesSymm(int inNBigBayesSymm){nBigBayesSymm = inNBigBayesSymm; return;}
  void SetBayesVal(std::vector<int> inBayesVal){bayesVal = inBayesVal; return;};
  void SetBayesVal(int inN, const int inBayesVal[]);
  void SetNSuperBayes(int inNSuperBayes){nSuperBayes = inNSuperBayes; return;}

  void SetNSVD(int inNSVD){nSVD = inNSVD; return;}

  void SetNHistDim(int inNHistDim){nHistDim = inNHistDim; return;}
  void SetHistTagBayes(std::vector<std::string> inHistTagBayes){histTagBayes = inHistTagBayes; return;};
  void SetHistTagBayes(int inN, const std::string inHistTagBayes[]);

  void SetHistTagSVD(std::vector<std::string> inHistTagSVD){histTagSVD = inHistTagSVD; return;};
  void SetHistTagSVD(int inN, const std::string inHistTagSVD[]);

  void SetHistBestBayes(std::vector<int> inHistBestBayes){histBestBayes = inHistBestBayes; return;};
  void SetHistBestSVD(std::vector<int> inHistBestSVD){histBestSVD = inHistBestSVD; return;};
  void SetHistBestBayes(int inN, const int inHistBestBayes[]);
  void SetHistBestSVD(int inN, const int inHistBestSVD[]);

  std::vector<int> StringToIntVect(std::string inStr);
  std::vector<double> StringToDoubleVect(std::string inStr);
  std::vector<std::string> StringToStringVect(std::string inStr);

  std::string to_string_with_precision(double a_value, const int n);
};


void cutPropagator::Clean()
{
  inFileNames.clear();
  inFullFileNames.clear();

  isPP = false;

  jtAbsEtaMax = -99;

  rcDiffFileName = "";
  flatPriorFileName = "";

  jecVarMC = -99;
  jerVarMC = -99;

  jecVarData = -99;

  nResponseMod = -1;
  responseMod.clear();
  jerVarData.clear();

  nJtAlgos = -1;
  jtAlgos.clear();
  minJtPtCut.clear();
  multiJtPtCut.clear();
  recoTruncPos.clear();

  nR = -1;
  rVals.clear();

  genJtPtBinsR1Cent0to10.clear();
  recoJtPtBinsR1Cent0to10.clear();
  genJtPtBinsR1Cent10to30.clear();
  recoJtPtBinsR1Cent10to30.clear();
  genJtPtBinsR1Cent30to50.clear();
  recoJtPtBinsR1Cent30to50.clear();
  genJtPtBinsR1Cent50to90.clear();
  recoJtPtBinsR1Cent50to90.clear();

  genJtPtBinsR2Cent0to10.clear();
  recoJtPtBinsR2Cent0to10.clear();
  genJtPtBinsR2Cent10to30.clear();
  recoJtPtBinsR2Cent10to30.clear();
  genJtPtBinsR2Cent30to50.clear();
  recoJtPtBinsR2Cent30to50.clear();
  genJtPtBinsR2Cent50to90.clear();
  recoJtPtBinsR2Cent50to90.clear();

  genJtPtBinsR3Cent0to10.clear();
  recoJtPtBinsR3Cent0to10.clear();
  genJtPtBinsR3Cent10to30.clear();
  recoJtPtBinsR3Cent10to30.clear();
  genJtPtBinsR3Cent30to50.clear();
  recoJtPtBinsR3Cent30to50.clear();
  genJtPtBinsR3Cent50to90.clear();
  recoJtPtBinsR3Cent50to90.clear();

  genJtPtBinsR4Cent0to10.clear();
  recoJtPtBinsR4Cent0to10.clear();
  genJtPtBinsR4Cent10to30.clear();
  recoJtPtBinsR4Cent10to30.clear();
  genJtPtBinsR4Cent30to50.clear();
  recoJtPtBinsR4Cent30to50.clear();
  genJtPtBinsR4Cent50to90.clear();
  recoJtPtBinsR4Cent50to90.clear();

  genJtPtBinsR5Cent0to10.clear();
  recoJtPtBinsR5Cent0to10.clear();
  genJtPtBinsR5Cent10to30.clear();
  recoJtPtBinsR5Cent10to30.clear();
  genJtPtBinsR5Cent30to50.clear();
  recoJtPtBinsR5Cent30to50.clear();
  genJtPtBinsR5Cent50to90.clear();
  recoJtPtBinsR5Cent50to90.clear();

  genJtPtBinsR6Cent0to10.clear();
  recoJtPtBinsR6Cent0to10.clear();
  genJtPtBinsR6Cent10to30.clear();
  recoJtPtBinsR6Cent10to30.clear();
  genJtPtBinsR6Cent30to50.clear();
  recoJtPtBinsR6Cent30to50.clear();
  genJtPtBinsR6Cent50to90.clear();
  recoJtPtBinsR6Cent50to90.clear();

  genJtPtBinsR8Cent0to10.clear();
  recoJtPtBinsR8Cent0to10.clear();
  genJtPtBinsR8Cent10to30.clear();
  recoJtPtBinsR8Cent10to30.clear();
  genJtPtBinsR8Cent30to50.clear();
  recoJtPtBinsR8Cent30to50.clear();
  genJtPtBinsR8Cent50to90.clear();
  recoJtPtBinsR8Cent50to90.clear();

  genJtPtBinsR10Cent0to10.clear();
  recoJtPtBinsR10Cent0to10.clear();
  genJtPtBinsR10Cent10to30.clear();
  recoJtPtBinsR10Cent10to30.clear();
  genJtPtBinsR10Cent30to50.clear();
  recoJtPtBinsR10Cent30to50.clear();
  genJtPtBinsR10Cent50to90.clear();
  recoJtPtBinsR10Cent50to90.clear();

  nGeneralBins = -1;

  generalBins.clear();

  nJtAbsEtaBins = -1;
  jtAbsEtaBinsLow.clear();
  jtAbsEtaBinsHi.clear();

  nPthats = -1;
  pthats.clear();
  pthatWeights.clear();

  nCentBins = -1;
  centBinsLow.clear();
  centBinsHi.clear();

  nID = -1;
  idStr.clear();
  jtPfCHMFCutLow.clear();
  jtPfCHMFCutHi.clear();
  jtPfMUMFCutLow.clear();
  jtPfMUMFCutHi.clear();
  jtPfNHFCutLow.clear();
  jtPfNHFCutHi.clear();
  jtPfNEFCutLow.clear();
  jtPfNEFCutHi.clear();
  jtPfMUFCutLow.clear();
  jtPfMUFCutHi.clear();
  jtPfCHFCutLow.clear();
  jtPfCHFCutHi.clear();
  jtPfCEFCutLow.clear();
  jtPfCEFCutHi.clear();
  jtPfMinMult.clear();
  jtPfMinChgMult.clear();

  nSyst = -1;
  systStr.clear();

  nBayes = -1;
  nBigBayesSymm = -1;
  bayesVal.clear();

  nSuperBayes = -1;

  nSVD = -1;

  nHistDim = -1;
  histTagBayes.clear();
  histTagSVD.clear();
  histBestBayes.clear();
  histBestSVD.clear();

  return;
}

bool cutPropagator::GetAllVarFromFile(TFile* inFile_p)
{
  inFile_p->cd();

  std::vector<std::string> cutDirCheck = returnRootFileContentsList(inFile_p, "TDirectoryFile", "cutDir", 1);
  if(cutDirCheck.size() == 0){
    std::cout << "cutPropagator::GetAllVarFromFile - given inFile_p \'" << inFile_p->GetName() << "\' contains no cutDir. return 1" << std::endl;
    return false;
  }

  std::vector<std::string> cutDirTNameds = returnTDirContentsList(inFile_p, "cutDir", "TNamed", "cutDir", 0, -1);
  std::vector<std::string> cutDirTDir = returnTDirContentsList(inFile_p, "cutDir", "TDirectoryFile", "", 0, -1);
  //returnRootFileContentsList(inFile_p, "TNamed", "cutDir");
  unsigned int pos = 0;
  while(pos < cutDirTNameds.size()){
    if(cutDirTNameds.at(pos).find("cutDir/") != std::string::npos) ++pos;
    else cutDirTNameds.erase(cutDirTNameds.begin()+pos);
  }

  if(cutDirTNameds.size() == 0){
    std::cout << "cutPropagator::GetAllVarFromFile - given inFile_p \'" << inFile_p->GetName() << "\' contains empty cutDir. return 1" << std::endl;
    return false;
  }

  //  std::cout << "Printing " << cutDirTNameds.size() << " TNameds..." << std::endl;
  for(unsigned int cI = 0; cI < cutDirTNameds.size(); ++cI){
    std::string tempStr = cutDirTNameds.at(cI);
    while(tempStr.find("/") != std::string::npos){tempStr.replace(0, tempStr.find("/")+1, "");}

    //    std::cout << " " << cI << "/" << cutDirTNameds.size() << ": " << cutDirTNameds.at(cI) << std::endl;
    bool isHistBayesBest = cutDirTNameds.at(cI).find("/unfoldDirBayes/") != std::string::npos;
    bool isHistSVDBest = cutDirTNameds.at(cI).find("/unfoldDirSVD/") != std::string::npos;

    if(isHistBayesBest){      
      histTagBayes.push_back(tempStr);
      histBestBayes.push_back(std::stoi(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle()));
      
      continue;
    }

    if(isHistSVDBest){    
      histTagSVD.push_back(tempStr);
      histBestSVD.push_back(std::stoi(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle()));
      
      continue;
    }


    if(isStrSame("nCentBins", tempStr)) nCentBins = std::stoi(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("centBinsLow", tempStr)) centBinsLow = StringToIntVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("centBinsHi", tempStr)) centBinsHi = StringToIntVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("jtAbsEtaMax", tempStr)) jtAbsEtaMax = std::stof(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("rcDiffFileName", tempStr)) rcDiffFileName = ((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle();
    else if(isStrSame("flatPriorFileName", tempStr)) flatPriorFileName = ((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle();
    else if(isStrSame("jecVarMC", tempStr)) jecVarMC = std::stof(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("jerVarMC", tempStr)) jerVarMC = std::stof(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("jecVarData", tempStr)) jecVarData = std::stof(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("nResponseMod", tempStr)) nResponseMod = std::stoi(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("nJtAlgos", tempStr)) nJtAlgos = std::stoi(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("jtAlgos", tempStr)) jtAlgos = StringToStringVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("minJtPtCut", tempStr)) minJtPtCut = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("multiJtPtCut", tempStr)) multiJtPtCut = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("recoTruncPos", tempStr)) recoTruncPos = StringToIntVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("nR", tempStr)) nR = std::stoi(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("rVals", tempStr)) rVals = StringToIntVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("genJtPtBinsR1Cent0to10", tempStr)) genJtPtBinsR1Cent0to10 = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("recoJtPtBinsR1Cent0to10", tempStr)) recoJtPtBinsR1Cent0to10 = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("genJtPtBinsR1Cent10to30", tempStr)) genJtPtBinsR1Cent10to30 = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("recoJtPtBinsR1Cent10to30", tempStr)) recoJtPtBinsR1Cent10to30 = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("genJtPtBinsR1Cent30to50", tempStr)) genJtPtBinsR1Cent30to50 = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("recoJtPtBinsR1Cent30to50", tempStr)) recoJtPtBinsR1Cent30to50 = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("genJtPtBinsR1Cent50to90", tempStr)) genJtPtBinsR1Cent50to90 = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("recoJtPtBinsR1Cent50to90", tempStr)) recoJtPtBinsR1Cent50to90 = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("genJtPtBinsR2Cent0to10", tempStr)) genJtPtBinsR2Cent0to10 = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("recoJtPtBinsR2Cent0to10", tempStr)) recoJtPtBinsR2Cent0to10 = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("genJtPtBinsR2Cent10to30", tempStr)) genJtPtBinsR2Cent10to30 = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("recoJtPtBinsR2Cent10to30", tempStr)) recoJtPtBinsR2Cent10to30 = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("genJtPtBinsR2Cent30to50", tempStr)) genJtPtBinsR2Cent30to50 = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("recoJtPtBinsR2Cent30to50", tempStr)) recoJtPtBinsR2Cent30to50 = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("genJtPtBinsR2Cent50to90", tempStr)) genJtPtBinsR2Cent50to90 = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("recoJtPtBinsR2Cent50to90", tempStr)) recoJtPtBinsR2Cent50to90 = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("genJtPtBinsR3Cent0to10", tempStr)) genJtPtBinsR3Cent0to10 = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("recoJtPtBinsR3Cent0to10", tempStr)) recoJtPtBinsR3Cent0to10 = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("genJtPtBinsR3Cent10to30", tempStr)) genJtPtBinsR3Cent10to30 = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("recoJtPtBinsR3Cent10to30", tempStr)) recoJtPtBinsR3Cent10to30 = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("genJtPtBinsR3Cent30to50", tempStr)) genJtPtBinsR3Cent30to50 = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("recoJtPtBinsR3Cent30to50", tempStr)) recoJtPtBinsR3Cent30to50 = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("genJtPtBinsR3Cent50to90", tempStr)) genJtPtBinsR3Cent50to90 = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("recoJtPtBinsR3Cent50to90", tempStr)) recoJtPtBinsR3Cent50to90 = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("genJtPtBinsR4Cent0to10", tempStr)) genJtPtBinsR4Cent0to10 = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("recoJtPtBinsR4Cent0to10", tempStr)) recoJtPtBinsR4Cent0to10 = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("genJtPtBinsR4Cent10to30", tempStr)) genJtPtBinsR4Cent10to30 = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("recoJtPtBinsR4Cent10to30", tempStr)) recoJtPtBinsR4Cent10to30 = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("genJtPtBinsR4Cent30to50", tempStr)) genJtPtBinsR4Cent30to50 = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("recoJtPtBinsR4Cent30to50", tempStr)) recoJtPtBinsR4Cent30to50 = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("genJtPtBinsR4Cent50to90", tempStr)) genJtPtBinsR4Cent50to90 = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("recoJtPtBinsR4Cent50to90", tempStr)) recoJtPtBinsR4Cent50to90 = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("genJtPtBinsR5Cent0to10", tempStr)) genJtPtBinsR5Cent0to10 = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("recoJtPtBinsR5Cent0to10", tempStr)) recoJtPtBinsR5Cent0to10 = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("genJtPtBinsR5Cent10to30", tempStr)) genJtPtBinsR5Cent10to30 = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("recoJtPtBinsR5Cent10to30", tempStr)) recoJtPtBinsR5Cent10to30 = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("genJtPtBinsR5Cent30to50", tempStr)) genJtPtBinsR5Cent30to50 = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("recoJtPtBinsR5Cent30to50", tempStr)) recoJtPtBinsR5Cent30to50 = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("genJtPtBinsR5Cent50to90", tempStr)) genJtPtBinsR5Cent50to90 = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("recoJtPtBinsR5Cent50to90", tempStr)) recoJtPtBinsR5Cent50to90 = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("genJtPtBinsR6Cent0to10", tempStr)) genJtPtBinsR6Cent0to10 = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("recoJtPtBinsR6Cent0to10", tempStr)) recoJtPtBinsR6Cent0to10 = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("genJtPtBinsR6Cent10to30", tempStr)) genJtPtBinsR6Cent10to30 = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("recoJtPtBinsR6Cent10to30", tempStr)) recoJtPtBinsR6Cent10to30 = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("genJtPtBinsR6Cent30to50", tempStr)) genJtPtBinsR6Cent30to50 = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("recoJtPtBinsR6Cent30to50", tempStr)) recoJtPtBinsR6Cent30to50 = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("genJtPtBinsR6Cent50to90", tempStr)) genJtPtBinsR6Cent50to90 = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("recoJtPtBinsR6Cent50to90", tempStr)) recoJtPtBinsR6Cent50to90 = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("genJtPtBinsR8Cent0to10", tempStr)) genJtPtBinsR8Cent0to10 = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("recoJtPtBinsR8Cent0to10", tempStr)) recoJtPtBinsR8Cent0to10 = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("genJtPtBinsR8Cent10to30", tempStr)) genJtPtBinsR8Cent10to30 = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("recoJtPtBinsR8Cent10to30", tempStr)) recoJtPtBinsR8Cent10to30 = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("genJtPtBinsR8Cent30to50", tempStr)) genJtPtBinsR8Cent30to50 = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("recoJtPtBinsR8Cent30to50", tempStr)) recoJtPtBinsR8Cent30to50 = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("genJtPtBinsR8Cent50to90", tempStr)) genJtPtBinsR8Cent50to90 = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("recoJtPtBinsR8Cent50to90", tempStr)) recoJtPtBinsR8Cent50to90 = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("genJtPtBinsR10Cent0to10", tempStr)) genJtPtBinsR10Cent0to10 = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("recoJtPtBinsR10Cent0to10", tempStr)) recoJtPtBinsR10Cent0to10 = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("genJtPtBinsR10Cent10to30", tempStr)) genJtPtBinsR10Cent10to30 = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("recoJtPtBinsR10Cent10to30", tempStr)) recoJtPtBinsR10Cent10to30 = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("genJtPtBinsR10Cent30to50", tempStr)) genJtPtBinsR10Cent30to50 = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("recoJtPtBinsR10Cent30to50", tempStr)) recoJtPtBinsR10Cent30to50 = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("genJtPtBinsR10Cent50to90", tempStr)) genJtPtBinsR10Cent50to90 = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("recoJtPtBinsR10Cent50to90", tempStr)) recoJtPtBinsR10Cent50to90 = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("nGeneralBins", tempStr)) nGeneralBins = std::stoi(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());    
    else if(isStrSame("nJtAbsEtaBins", tempStr)) nJtAbsEtaBins = std::stoi(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("inFileNames", tempStr)) inFileNames = StringToStringVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(tempStr.find("inFullFileNames_") != std::string::npos){
      std::string tempStr2 = ((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle();
      inFullFileNames.push_back(tempStr2);
    }
    else if(isStrSame("isPP", tempStr)) isPP = std::stoi(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("nPthats", tempStr)) nPthats = std::stoi(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("pthats", tempStr)) pthats = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("pthatWeights", tempStr)) pthatWeights = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("responseMod", tempStr)) responseMod = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("jerVarData", tempStr)) jerVarData = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("generalBins", tempStr)) generalBins = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("jtAbsEtaBinsLow", tempStr)) jtAbsEtaBinsLow = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("jtAbsEtaBinsHi", tempStr)) jtAbsEtaBinsHi = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("nID", tempStr)) nID = std::stoi(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("idStr", tempStr)) idStr = StringToStringVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("jtPfCHMFCutLow", tempStr)) jtPfCHMFCutLow = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("jtPfCHMFCutHi", tempStr)) jtPfCHMFCutHi = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("jtPfMUMFCutLow", tempStr)) jtPfMUMFCutLow = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("jtPfMUMFCutHi", tempStr)) jtPfMUMFCutHi = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("jtPfNHFCutLow", tempStr)) jtPfNHFCutLow = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("jtPfNHFCutHi", tempStr)) jtPfNHFCutHi = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("jtPfNEFCutLow", tempStr)) jtPfNEFCutLow = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("jtPfNEFCutHi", tempStr)) jtPfNEFCutHi = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("jtPfMUFCutLow", tempStr)) jtPfMUFCutLow = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("jtPfMUFCutHi", tempStr)) jtPfMUFCutHi = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("jtPfCHFCutLow", tempStr)) jtPfCHFCutLow = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("jtPfCHFCutHi", tempStr)) jtPfCHFCutHi = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("jtPfCEFCutLow", tempStr)) jtPfCEFCutLow = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("jtPfCEFCutHi", tempStr)) jtPfCEFCutHi = StringToDoubleVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("jtPfMinMult", tempStr)) jtPfMinMult = StringToIntVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("jtPfMinChgMult", tempStr)) jtPfMinChgMult = StringToIntVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("nSyst", tempStr)) nSyst = std::stoi(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("systStr", tempStr)) systStr = StringToStringVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("nBayes", tempStr)) nBayes = std::stoi(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("nBigBayesSymm", tempStr)) nBigBayesSymm = std::stoi(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("bayesVal", tempStr)) bayesVal = StringToIntVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("nSuperBayes", tempStr)) nSuperBayes = std::stoi(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("nSVD", tempStr)) nSVD = std::stoi(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("nHistDim", tempStr)) nHistDim = std::stoi(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("histTagBayes", tempStr)) histTagBayes = StringToStringVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("histTagSVD", tempStr)) histTagSVD = StringToStringVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("histBestBayes", tempStr)) histBestBayes = StringToIntVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else if(isStrSame("histBestSVD", tempStr)) histBestSVD = StringToIntVect(((TNamed*)inFile_p->Get(cutDirTNameds.at(cI).c_str()))->GetTitle());
    else{
      std::cout << "WARNING: TNAMED \'" << tempStr << "\' is unaccounted for in cutPropagator. Consider fixing" << std::endl;
    }
  }
  
  //  std::cout << "Handling cutPropDirs: " << std::endl;
  for(unsigned int dI = 0; dI < cutDirTDir.size(); ++dI){
    //    std::cout << " " << dI << "/" << cutDirTDir.size() << ": " << cutDirTDir.at(dI) << std::endl;
    TDirectory* tempDir_p = (TDirectory*)inFile_p->Get(cutDirTDir.at(dI).c_str());
    TIter next(tempDir_p->GetListOfKeys());
    TKey* tempKey = NULL;

    while( (tempKey = (TKey*)next()) ){
      const std::string name = tempKey->GetName();
      const std::string className = tempKey->GetClassName();

      if(className.find("TNamed") == std::string::npos) continue;

      TNamed* tempName = (TNamed*)tempKey->ReadObj();
      if(cutDirTDir.at(dI).find("unfoldDirBayes") != std::string::npos){
	histTagBayes.push_back(std::string(tempName->GetName()));
	histBestBayes.push_back(std::stoi(std::string(tempName->GetTitle())));
      }
      else if(cutDirTDir.at(dI).find("unfoldDirSVD") != std::string::npos){
	histTagSVD.push_back(std::string(tempName->GetName()));
	histBestSVD.push_back(std::stoi(std::string(tempName->GetTitle())));
      }
      else if(cutDirTDir.at(dI).find("subDir") != std::string::npos){
	inFullFileNames.push_back(std::string(tempName->GetTitle()));
      }
    }

  }

  return true;
}


std::string cutPropagator::VectToString(std::vector<double> inVect, int prec)
{
  std::string retStr = "";
  for(unsigned int vI = 0; vI < inVect.size(); ++vI){
    retStr = retStr + prettyString(inVect[vI], prec, false) + ",";
  }
  return retStr;
}


bool cutPropagator::WriteAllVarToFile(TFile* inFile_p, TDirectory* inDir_p, TDirectory* inSubDir_p, TDirectory* unfoldDirBayes_p = NULL, TDirectory* unfoldDirSVD_p = NULL)
{
  if(inDir_p == NULL){
    std::cout << "cutPropafator::WriteAllVarToFile - Given indir is null. return false" << std::endl;
    return false;
  }
  else if(inSubDir_p == NULL){
    std::cout << "cutPropafator::WriteAllVarToFile - Given insubdir is null. return false" << std::endl;
    return false;
  }
  else if(std::string(inDir_p->GetName()).size() != std::string("cutDir").size()){
    std::cout << "cutPropafator::WriteAllVarToFile - Given indir is not \'cutDir\'. return false" << std::endl;
    return false;
  }
  else if(std::string(inDir_p->GetName()).find("cutDir") == std::string::npos){
    std::cout << "cutPropafator::WriteAllVarToFile - Given indir is not \'cutDir\'. return false" << std::endl;
    return false;
  }

  if(doLocalDebug || doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  inFile_p->cd();
  inDir_p->cd();

  std::string inFileNames2 = "";
  for(unsigned int i = 0; i < inFileNames.size(); ++i){
    inFileNames2 = inFileNames2 + inFileNames.at(i) + ",";
  }

  if(doLocalDebug || doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  std::string responseModStr = "";
  std::string jerVarDataStr = "";
  for(int jI = 0; jI < nResponseMod; ++jI){
    responseModStr = responseModStr + std::to_string(responseMod.at(jI)) + ",";
    jerVarDataStr = jerVarDataStr + std::to_string(jerVarData.at(jI)) + ",";
  }

  if(doLocalDebug || doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  std::string rValsStr = "";
  for(Int_t rI = 0; rI < nR; ++rI){
    rValsStr = rValsStr + std::to_string(rVals.at(rI)) + ",";
  }

  std::string generalBinsStr = "";
  for(unsigned int jI = 0; jI < generalBins.size(); ++jI){
    generalBinsStr = generalBinsStr + std::to_string(generalBins.at(jI)) + ",";
  }


  std::string jtAlgosStr = "";
  std::string minJtPtCutStr = "";
  std::string multiJtPtCutStr = "";
  std::string recoTruncPosStr = "";

  for(int jI = 0; jI < nJtAlgos; ++jI){
    if(jtAlgos.size() != 0) jtAlgosStr = jtAlgosStr + jtAlgos.at(jI) + ",";
    if(minJtPtCut.size() != 0) minJtPtCutStr = minJtPtCutStr + prettyString(minJtPtCut.at(jI), 1, false) + ",";
    if(multiJtPtCut.size() != 0) multiJtPtCutStr = multiJtPtCutStr + prettyString(multiJtPtCut.at(jI), 1, false) + ",";
    if(recoTruncPos.size() != 0) recoTruncPosStr = recoTruncPosStr + std::to_string((int)(recoTruncPos.at(jI))) + ",";
  }

  std::string jtAbsEtaBinsLowStr = "";
  std::string jtAbsEtaBinsHiStr = "";
  for(int jI = 0; jI < nJtAbsEtaBins; ++jI){
    jtAbsEtaBinsLowStr = jtAbsEtaBinsLowStr + std::to_string(jtAbsEtaBinsLow.at(jI)) + ",";
    jtAbsEtaBinsHiStr = jtAbsEtaBinsHiStr + std::to_string(jtAbsEtaBinsHi.at(jI)) + ",";
  }

  std::string nPthatsStr = std::to_string(nPthats);
  std::string pthatsStr = "";
  for(unsigned int jI = 0; jI < pthats.size(); ++jI){
    pthatsStr = pthatsStr + std::to_string(pthats.at(jI)) + ",";
  }

  std::string pthatWeightsStr = "";
  for(unsigned int jI = 0; jI < pthatWeights.size(); ++jI){
    pthatWeightsStr = pthatWeightsStr + to_string_with_precision(pthatWeights.at(jI), 12) + ",";
    std::cout << "Weight orig, string: " << pthatWeights.at(jI) << ", " << to_string_with_precision(pthatWeights.at(jI), 12) << std::endl;
  }

  std::string centBinsLowStr = "";
  for(int jI = 0; jI < nCentBins; ++jI){
    centBinsLowStr = centBinsLowStr + std::to_string(centBinsLow.at(jI)) + ",";
  }

  std::string centBinsHiStr = "";
  for(int jI = 0; jI < nCentBins; ++jI){
    centBinsHiStr = centBinsHiStr + std::to_string(centBinsHi.at(jI)) + ",";
  }

  if(doLocalDebug || doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  std::string nIDStr = std::to_string(nID);
  std::string idStr2 = "";
  std::string jtPfCHMFCutLowStr = "";
  std::string jtPfCHMFCutHiStr = "";
  std::string jtPfMUMFCutLowStr = "";
  std::string jtPfMUMFCutHiStr = "";
  std::string jtPfNHFCutLowStr = "";
  std::string jtPfNHFCutHiStr = "";
  std::string jtPfNEFCutLowStr = "";
  std::string jtPfNEFCutHiStr = "";
  std::string jtPfMUFCutLowStr = "";
  std::string jtPfMUFCutHiStr = "";
  std::string jtPfCHFCutLowStr = "";
  std::string jtPfCHFCutHiStr = "";
  std::string jtPfCEFCutLowStr = "";
  std::string jtPfCEFCutHiStr = "";
  std::string jtPfMinMultStr = "";
  std::string jtPfMinChgMultStr = "";

  for(int iI = 0; iI < nID; ++iI){
    idStr2 = idStr2 + idStr.at(iI) + ",";
    jtPfCHMFCutLowStr = jtPfCHMFCutLowStr + std::to_string(jtPfCHMFCutLow.at(iI)) + ",";
    jtPfCHMFCutHiStr = jtPfCHMFCutHiStr + std::to_string(jtPfCHMFCutHi.at(iI)) + ",";
    jtPfMUMFCutLowStr = jtPfMUMFCutLowStr + std::to_string(jtPfMUMFCutLow.at(iI)) + ",";
    jtPfMUMFCutHiStr = jtPfMUMFCutHiStr + std::to_string(jtPfMUMFCutHi.at(iI)) + ",";
    jtPfNHFCutLowStr = jtPfNHFCutLowStr + std::to_string(jtPfNHFCutLow.at(iI)) + ",";
    jtPfNHFCutHiStr = jtPfNHFCutHiStr + std::to_string(jtPfNHFCutHi.at(iI)) + ",";
    jtPfNEFCutLowStr = jtPfNEFCutLowStr + std::to_string(jtPfNEFCutLow.at(iI)) + ",";
    jtPfNEFCutHiStr = jtPfNEFCutHiStr + std::to_string(jtPfNEFCutHi.at(iI)) + ",";
    jtPfMUFCutLowStr = jtPfMUFCutLowStr + std::to_string(jtPfMUFCutLow.at(iI)) + ",";
    jtPfMUFCutHiStr = jtPfMUFCutHiStr + std::to_string(jtPfMUFCutHi.at(iI)) + ",";
    jtPfCHFCutLowStr = jtPfCHFCutLowStr + std::to_string(jtPfCHFCutLow.at(iI)) + ",";
    jtPfCHFCutHiStr = jtPfCHFCutHiStr + std::to_string(jtPfCHFCutHi.at(iI)) + ",";
    jtPfCEFCutLowStr = jtPfCEFCutLowStr + std::to_string(jtPfCEFCutLow.at(iI)) + ",";
    jtPfCEFCutHiStr = jtPfCEFCutHiStr + std::to_string(jtPfCEFCutHi.at(iI)) + ",";
    jtPfMinMultStr = jtPfMinMultStr + std::to_string(jtPfMinMult.at(iI)) + ",";
    jtPfMinChgMultStr = jtPfMinChgMultStr + std::to_string(jtPfMinChgMult.at(iI)) + ",";
  }

  std::string nSystStr = std::to_string(nSyst);
  std::string systStr2 = "";

  for(int sI = 0; sI < nSyst; ++sI){
    systStr2 = systStr2 + systStr.at(sI) + ",";
  }

  if(doLocalDebug || doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  std::string nBayesStr = std::to_string(nBayes);
  std::string nBigBayesSymmStr = std::to_string(nBigBayesSymm);
  std::string nSuperBayesStr = std::to_string(nSuperBayes);
  std::string bayesVal2 = "";

  std::string nSVDStr = std::to_string(nSVD);

  for(int sI = 0; sI < nBayes; ++sI){
    bayesVal2 = bayesVal2 + std::to_string(bayesVal.at(sI)) + ",";
  }

  std::string nHistDimStr = std::to_string(nHistDim);
  

  if(doLocalDebug || doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;


  std::string genJtPtBinsR1Cent0to10Str = VectToString(genJtPtBinsR1Cent0to10, 1);
  std::string recoJtPtBinsR1Cent0to10Str = VectToString(recoJtPtBinsR1Cent0to10, 1);
  std::string genJtPtBinsR1Cent10to30Str = VectToString(genJtPtBinsR1Cent10to30, 1);
  std::string recoJtPtBinsR1Cent10to30Str = VectToString(recoJtPtBinsR1Cent10to30, 1);
  std::string genJtPtBinsR1Cent30to50Str = VectToString(genJtPtBinsR1Cent30to50, 1);
  std::string recoJtPtBinsR1Cent30to50Str = VectToString(recoJtPtBinsR1Cent30to50, 1);
  std::string genJtPtBinsR1Cent50to90Str = VectToString(genJtPtBinsR1Cent50to90, 1);
  std::string recoJtPtBinsR1Cent50to90Str = VectToString(recoJtPtBinsR1Cent50to90, 1);

  std::string genJtPtBinsR2Cent0to10Str = VectToString(genJtPtBinsR2Cent0to10, 1);
  std::string recoJtPtBinsR2Cent0to10Str = VectToString(recoJtPtBinsR2Cent0to10, 1);
  std::string genJtPtBinsR2Cent10to30Str = VectToString(genJtPtBinsR2Cent10to30, 1);
  std::string recoJtPtBinsR2Cent10to30Str = VectToString(recoJtPtBinsR2Cent10to30, 1);
  std::string genJtPtBinsR2Cent30to50Str = VectToString(genJtPtBinsR2Cent30to50, 1);
  std::string recoJtPtBinsR2Cent30to50Str = VectToString(recoJtPtBinsR2Cent30to50, 1);
  std::string genJtPtBinsR2Cent50to90Str = VectToString(genJtPtBinsR2Cent50to90, 1);
  std::string recoJtPtBinsR2Cent50to90Str = VectToString(recoJtPtBinsR2Cent50to90, 1);

  std::string genJtPtBinsR3Cent0to10Str = VectToString(genJtPtBinsR3Cent0to10, 1);
  std::string recoJtPtBinsR3Cent0to10Str = VectToString(recoJtPtBinsR3Cent0to10, 1);
  std::string genJtPtBinsR3Cent10to30Str = VectToString(genJtPtBinsR3Cent10to30, 1);
  std::string recoJtPtBinsR3Cent10to30Str = VectToString(recoJtPtBinsR3Cent10to30, 1);
  std::string genJtPtBinsR3Cent30to50Str = VectToString(genJtPtBinsR3Cent30to50, 1);
  std::string recoJtPtBinsR3Cent30to50Str = VectToString(recoJtPtBinsR3Cent30to50, 1);
  std::string genJtPtBinsR3Cent50to90Str = VectToString(genJtPtBinsR3Cent50to90, 1);
  std::string recoJtPtBinsR3Cent50to90Str = VectToString(recoJtPtBinsR3Cent50to90, 1);

  std::string genJtPtBinsR4Cent0to10Str = VectToString(genJtPtBinsR4Cent0to10, 1);
  std::string recoJtPtBinsR4Cent0to10Str = VectToString(recoJtPtBinsR4Cent0to10, 1);
  std::string genJtPtBinsR4Cent10to30Str = VectToString(genJtPtBinsR4Cent10to30, 1);
  std::string recoJtPtBinsR4Cent10to30Str = VectToString(recoJtPtBinsR4Cent10to30, 1);
  std::string genJtPtBinsR4Cent30to50Str = VectToString(genJtPtBinsR4Cent30to50, 1);
  std::string recoJtPtBinsR4Cent30to50Str = VectToString(recoJtPtBinsR4Cent30to50, 1);
  std::string genJtPtBinsR4Cent50to90Str = VectToString(genJtPtBinsR4Cent50to90, 1);
  std::string recoJtPtBinsR4Cent50to90Str = VectToString(recoJtPtBinsR4Cent50to90, 1);

  std::string genJtPtBinsR5Cent0to10Str = VectToString(genJtPtBinsR5Cent0to10, 1);
  std::string recoJtPtBinsR5Cent0to10Str = VectToString(recoJtPtBinsR5Cent0to10, 1);
  std::string genJtPtBinsR5Cent10to30Str = VectToString(genJtPtBinsR5Cent10to30, 1);
  std::string recoJtPtBinsR5Cent10to30Str = VectToString(recoJtPtBinsR5Cent10to30, 1);
  std::string genJtPtBinsR5Cent30to50Str = VectToString(genJtPtBinsR5Cent30to50, 1);
  std::string recoJtPtBinsR5Cent30to50Str = VectToString(recoJtPtBinsR5Cent30to50, 1);
  std::string genJtPtBinsR5Cent50to90Str = VectToString(genJtPtBinsR5Cent50to90, 1);
  std::string recoJtPtBinsR5Cent50to90Str = VectToString(recoJtPtBinsR5Cent50to90, 1);

  std::string genJtPtBinsR6Cent0to10Str = VectToString(genJtPtBinsR6Cent0to10, 1);
  std::string recoJtPtBinsR6Cent0to10Str = VectToString(recoJtPtBinsR6Cent0to10, 1);
  std::string genJtPtBinsR6Cent10to30Str = VectToString(genJtPtBinsR6Cent10to30, 1);
  std::string recoJtPtBinsR6Cent10to30Str = VectToString(recoJtPtBinsR6Cent10to30, 1);
  std::string genJtPtBinsR6Cent30to50Str = VectToString(genJtPtBinsR6Cent30to50, 1);
  std::string recoJtPtBinsR6Cent30to50Str = VectToString(recoJtPtBinsR6Cent30to50, 1);
  std::string genJtPtBinsR6Cent50to90Str = VectToString(genJtPtBinsR6Cent50to90, 1);
  std::string recoJtPtBinsR6Cent50to90Str = VectToString(recoJtPtBinsR6Cent50to90, 1);

  std::string genJtPtBinsR8Cent0to10Str = VectToString(genJtPtBinsR8Cent0to10, 1);
  std::string recoJtPtBinsR8Cent0to10Str = VectToString(recoJtPtBinsR8Cent0to10, 1);
  std::string genJtPtBinsR8Cent10to30Str = VectToString(genJtPtBinsR8Cent10to30, 1);
  std::string recoJtPtBinsR8Cent10to30Str = VectToString(recoJtPtBinsR8Cent10to30, 1);
  std::string genJtPtBinsR8Cent30to50Str = VectToString(genJtPtBinsR8Cent30to50, 1);
  std::string recoJtPtBinsR8Cent30to50Str = VectToString(recoJtPtBinsR8Cent30to50, 1);
  std::string genJtPtBinsR8Cent50to90Str = VectToString(genJtPtBinsR8Cent50to90, 1);
  std::string recoJtPtBinsR8Cent50to90Str = VectToString(recoJtPtBinsR8Cent50to90, 1);

  std::string genJtPtBinsR10Cent0to10Str = VectToString(genJtPtBinsR10Cent0to10, 1);
  std::string recoJtPtBinsR10Cent0to10Str = VectToString(recoJtPtBinsR10Cent0to10, 1);
  std::string genJtPtBinsR10Cent10to30Str = VectToString(genJtPtBinsR10Cent10to30, 1);
  std::string recoJtPtBinsR10Cent10to30Str = VectToString(recoJtPtBinsR10Cent10to30, 1);
  std::string genJtPtBinsR10Cent30to50Str = VectToString(genJtPtBinsR10Cent30to50, 1);
  std::string recoJtPtBinsR10Cent30to50Str = VectToString(recoJtPtBinsR10Cent30to50, 1);
  std::string genJtPtBinsR10Cent50to90Str = VectToString(genJtPtBinsR10Cent50to90, 1);
  std::string recoJtPtBinsR10Cent50to90Str = VectToString(recoJtPtBinsR10Cent50to90, 1);


  std::vector<TNamed*> nameList;  
  nameList.push_back(new TNamed("inFileNames", inFileNames2.c_str()));
  nameList.push_back(new TNamed("isPP", std::to_string(isPP).c_str()));
  nameList.push_back(new TNamed("jtAbsEtaMax", std::to_string(jtAbsEtaMax).c_str()));
  nameList.push_back(new TNamed("rcDiffFileName", rcDiffFileName.c_str()));
  nameList.push_back(new TNamed("flatPriorFileName", flatPriorFileName.c_str()));
  nameList.push_back(new TNamed("jecVarMC", std::to_string(jecVarMC).c_str()));
  nameList.push_back(new TNamed("jerVarMC", std::to_string(jerVarMC).c_str()));
  nameList.push_back(new TNamed("jecVarData", std::to_string(jecVarData).c_str()));
  nameList.push_back(new TNamed("nResponseMod", std::to_string(nResponseMod).c_str()));
  nameList.push_back(new TNamed("responseMod", responseModStr.c_str()));
  nameList.push_back(new TNamed("jerVarData", jerVarDataStr.c_str()));
  nameList.push_back(new TNamed("nJtAlgos", std::to_string(nJtAlgos).c_str()));
  nameList.push_back(new TNamed("jtAlgos", jtAlgosStr.c_str()));
  nameList.push_back(new TNamed("minJtPtCut", minJtPtCutStr.c_str()));
  nameList.push_back(new TNamed("multiJtPtCut", multiJtPtCutStr.c_str()));
  nameList.push_back(new TNamed("recoTruncPos", recoTruncPosStr.c_str()));
  nameList.push_back(new TNamed("nR", std::to_string(nR).c_str()));
  nameList.push_back(new TNamed("rVals", rValsStr.c_str()));
  nameList.push_back(new TNamed("genJtPtBinsR1Cent0to10", genJtPtBinsR1Cent0to10Str.c_str()));
  nameList.push_back(new TNamed("recoJtPtBinsR1Cent0to10", recoJtPtBinsR1Cent0to10Str.c_str()));
  nameList.push_back(new TNamed("genJtPtBinsR1Cent10to30", genJtPtBinsR1Cent10to30Str.c_str()));
  nameList.push_back(new TNamed("recoJtPtBinsR1Cent10to30", recoJtPtBinsR1Cent10to30Str.c_str()));
  nameList.push_back(new TNamed("genJtPtBinsR1Cent30to50", genJtPtBinsR1Cent30to50Str.c_str()));
  nameList.push_back(new TNamed("recoJtPtBinsR1Cent30to50", recoJtPtBinsR1Cent30to50Str.c_str()));
  nameList.push_back(new TNamed("genJtPtBinsR1Cent50to90", genJtPtBinsR1Cent50to90Str.c_str()));
  nameList.push_back(new TNamed("recoJtPtBinsR1Cent50to90", recoJtPtBinsR1Cent50to90Str.c_str()));
  nameList.push_back(new TNamed("genJtPtBinsR2Cent0to10", genJtPtBinsR2Cent0to10Str.c_str()));
  nameList.push_back(new TNamed("recoJtPtBinsR2Cent0to10", recoJtPtBinsR2Cent0to10Str.c_str()));
  nameList.push_back(new TNamed("genJtPtBinsR2Cent10to30", genJtPtBinsR2Cent10to30Str.c_str()));
  nameList.push_back(new TNamed("recoJtPtBinsR2Cent10to30", recoJtPtBinsR2Cent10to30Str.c_str()));
  nameList.push_back(new TNamed("genJtPtBinsR2Cent30to50", genJtPtBinsR2Cent30to50Str.c_str()));
  nameList.push_back(new TNamed("recoJtPtBinsR2Cent30to50", recoJtPtBinsR2Cent30to50Str.c_str()));
  nameList.push_back(new TNamed("genJtPtBinsR2Cent50to90", genJtPtBinsR2Cent50to90Str.c_str()));
  nameList.push_back(new TNamed("recoJtPtBinsR2Cent50to90", recoJtPtBinsR2Cent50to90Str.c_str()));

  nameList.push_back(new TNamed("genJtPtBinsR3Cent0to10", genJtPtBinsR3Cent0to10Str.c_str()));
  nameList.push_back(new TNamed("recoJtPtBinsR3Cent0to10", recoJtPtBinsR3Cent0to10Str.c_str()));
  nameList.push_back(new TNamed("genJtPtBinsR3Cent10to30", genJtPtBinsR3Cent10to30Str.c_str()));
  nameList.push_back(new TNamed("recoJtPtBinsR3Cent10to30", recoJtPtBinsR3Cent10to30Str.c_str()));
  nameList.push_back(new TNamed("genJtPtBinsR3Cent30to50", genJtPtBinsR3Cent30to50Str.c_str()));
  nameList.push_back(new TNamed("recoJtPtBinsR3Cent30to50", recoJtPtBinsR3Cent30to50Str.c_str()));
  nameList.push_back(new TNamed("genJtPtBinsR3Cent50to90", genJtPtBinsR3Cent50to90Str.c_str()));
  nameList.push_back(new TNamed("recoJtPtBinsR3Cent50to90", recoJtPtBinsR3Cent50to90Str.c_str()));

  nameList.push_back(new TNamed("genJtPtBinsR4Cent0to10", genJtPtBinsR4Cent0to10Str.c_str()));
  nameList.push_back(new TNamed("recoJtPtBinsR4Cent0to10", recoJtPtBinsR4Cent0to10Str.c_str()));
  nameList.push_back(new TNamed("genJtPtBinsR4Cent10to30", genJtPtBinsR4Cent10to30Str.c_str()));
  nameList.push_back(new TNamed("recoJtPtBinsR4Cent10to30", recoJtPtBinsR4Cent10to30Str.c_str()));
  nameList.push_back(new TNamed("genJtPtBinsR4Cent30to50", genJtPtBinsR4Cent30to50Str.c_str()));
  nameList.push_back(new TNamed("recoJtPtBinsR4Cent30to50", recoJtPtBinsR4Cent30to50Str.c_str()));
  nameList.push_back(new TNamed("genJtPtBinsR4Cent50to90", genJtPtBinsR4Cent50to90Str.c_str()));
  nameList.push_back(new TNamed("recoJtPtBinsR4Cent50to90", recoJtPtBinsR4Cent50to90Str.c_str()));

  nameList.push_back(new TNamed("genJtPtBinsR5Cent0to10", genJtPtBinsR5Cent0to10Str.c_str()));
  nameList.push_back(new TNamed("recoJtPtBinsR5Cent0to10", recoJtPtBinsR5Cent0to10Str.c_str()));
  nameList.push_back(new TNamed("genJtPtBinsR5Cent10to30", genJtPtBinsR5Cent10to30Str.c_str()));
  nameList.push_back(new TNamed("recoJtPtBinsR5Cent10to30", recoJtPtBinsR5Cent10to30Str.c_str()));
  nameList.push_back(new TNamed("genJtPtBinsR5Cent30to50", genJtPtBinsR5Cent30to50Str.c_str()));
  nameList.push_back(new TNamed("recoJtPtBinsR5Cent30to50", recoJtPtBinsR5Cent30to50Str.c_str()));
  nameList.push_back(new TNamed("genJtPtBinsR5Cent50to90", genJtPtBinsR5Cent50to90Str.c_str()));
  nameList.push_back(new TNamed("recoJtPtBinsR5Cent50to90", recoJtPtBinsR5Cent50to90Str.c_str()));

  nameList.push_back(new TNamed("genJtPtBinsR6Cent0to10", genJtPtBinsR6Cent0to10Str.c_str()));
  nameList.push_back(new TNamed("recoJtPtBinsR6Cent0to10", recoJtPtBinsR6Cent0to10Str.c_str()));
  nameList.push_back(new TNamed("genJtPtBinsR6Cent10to30", genJtPtBinsR6Cent10to30Str.c_str()));
  nameList.push_back(new TNamed("recoJtPtBinsR6Cent10to30", recoJtPtBinsR6Cent10to30Str.c_str()));
  nameList.push_back(new TNamed("genJtPtBinsR6Cent30to50", genJtPtBinsR6Cent30to50Str.c_str()));
  nameList.push_back(new TNamed("recoJtPtBinsR6Cent30to50", recoJtPtBinsR6Cent30to50Str.c_str()));
  nameList.push_back(new TNamed("genJtPtBinsR6Cent50to90", genJtPtBinsR6Cent50to90Str.c_str()));
  nameList.push_back(new TNamed("recoJtPtBinsR6Cent50to90", recoJtPtBinsR6Cent50to90Str.c_str()));

  nameList.push_back(new TNamed("genJtPtBinsR8Cent0to10", genJtPtBinsR8Cent0to10Str.c_str()));
  nameList.push_back(new TNamed("recoJtPtBinsR8Cent0to10", recoJtPtBinsR8Cent0to10Str.c_str()));
  nameList.push_back(new TNamed("genJtPtBinsR8Cent10to30", genJtPtBinsR8Cent10to30Str.c_str()));
  nameList.push_back(new TNamed("recoJtPtBinsR8Cent10to30", recoJtPtBinsR8Cent10to30Str.c_str()));
  nameList.push_back(new TNamed("genJtPtBinsR8Cent30to50", genJtPtBinsR8Cent30to50Str.c_str()));
  nameList.push_back(new TNamed("recoJtPtBinsR8Cent30to50", recoJtPtBinsR8Cent30to50Str.c_str()));
  nameList.push_back(new TNamed("genJtPtBinsR8Cent50to90", genJtPtBinsR8Cent50to90Str.c_str()));
  nameList.push_back(new TNamed("recoJtPtBinsR8Cent50to90", recoJtPtBinsR8Cent50to90Str.c_str()));

  nameList.push_back(new TNamed("genJtPtBinsR10Cent0to10", genJtPtBinsR10Cent0to10Str.c_str()));
  nameList.push_back(new TNamed("recoJtPtBinsR10Cent0to10", recoJtPtBinsR10Cent0to10Str.c_str()));
  nameList.push_back(new TNamed("genJtPtBinsR10Cent10to30", genJtPtBinsR10Cent10to30Str.c_str()));
  nameList.push_back(new TNamed("recoJtPtBinsR10Cent10to30", recoJtPtBinsR10Cent10to30Str.c_str()));
  nameList.push_back(new TNamed("genJtPtBinsR10Cent30to50", genJtPtBinsR10Cent30to50Str.c_str()));
  nameList.push_back(new TNamed("recoJtPtBinsR10Cent30to50", recoJtPtBinsR10Cent30to50Str.c_str()));
  nameList.push_back(new TNamed("genJtPtBinsR10Cent50to90", genJtPtBinsR10Cent50to90Str.c_str()));
  nameList.push_back(new TNamed("recoJtPtBinsR10Cent50to90", recoJtPtBinsR10Cent50to90Str.c_str()));

  nameList.push_back(new TNamed("nGeneralBins", std::to_string(nGeneralBins).c_str()));

  nameList.push_back(new TNamed("generalBins", generalBinsStr.c_str()));

  nameList.push_back(new TNamed("nJtAbsEtaBins", std::to_string(nJtAbsEtaBins).c_str()));
  nameList.push_back(new TNamed("jtAbsEtaBinsLow", jtAbsEtaBinsLowStr.c_str()));
  nameList.push_back(new TNamed("jtAbsEtaBinsHi", jtAbsEtaBinsHiStr.c_str()));
  nameList.push_back(new TNamed("nPthats", nPthatsStr.c_str()));
  nameList.push_back(new TNamed("pthats", pthatsStr.c_str()));
  nameList.push_back(new TNamed("pthatWeights", pthatWeightsStr.c_str()));
  nameList.push_back(new TNamed("nCentBins", std::to_string(nCentBins).c_str()));
  nameList.push_back(new TNamed("centBinsLow", centBinsLowStr.c_str()));
  nameList.push_back(new TNamed("centBinsHi", centBinsHiStr.c_str()));
  nameList.push_back(new TNamed("nID", nIDStr.c_str()));
  nameList.push_back(new TNamed("idStr", idStr2.c_str()));
  nameList.push_back(new TNamed("jtPfCHMFCutLow", jtPfCHMFCutLowStr.c_str()));
  nameList.push_back(new TNamed("jtPfCHMFCutHi", jtPfCHMFCutHiStr.c_str()));
  nameList.push_back(new TNamed("jtPfMUMFCutLow", jtPfMUMFCutLowStr.c_str()));
  nameList.push_back(new TNamed("jtPfMUMFCutHi", jtPfMUMFCutHiStr.c_str()));
  nameList.push_back(new TNamed("jtPfNHFCutLow", jtPfNHFCutLowStr.c_str()));
  nameList.push_back(new TNamed("jtPfNHFCutHi", jtPfNHFCutHiStr.c_str()));
  nameList.push_back(new TNamed("jtPfNEFCutLow", jtPfNEFCutLowStr.c_str()));
  nameList.push_back(new TNamed("jtPfNEFCutHi", jtPfNEFCutHiStr.c_str()));
  nameList.push_back(new TNamed("jtPfMUFCutLow", jtPfMUFCutLowStr.c_str()));
  nameList.push_back(new TNamed("jtPfMUFCutHi", jtPfMUFCutHiStr.c_str()));
  nameList.push_back(new TNamed("jtPfCHFCutLow", jtPfCHFCutLowStr.c_str()));
  nameList.push_back(new TNamed("jtPfCHFCutHi", jtPfCHFCutHiStr.c_str()));
  nameList.push_back(new TNamed("jtPfCEFCutLow", jtPfCEFCutLowStr.c_str()));
  nameList.push_back(new TNamed("jtPfCEFCutHi", jtPfCEFCutHiStr.c_str()));
  nameList.push_back(new TNamed("jtPfMinMult", jtPfMinMultStr.c_str()));
  nameList.push_back(new TNamed("jtPfMinChgMult", jtPfMinChgMultStr.c_str()));
  nameList.push_back(new TNamed("nSyst", nSystStr.c_str()));
  nameList.push_back(new TNamed("systStr", systStr2.c_str()));
  nameList.push_back(new TNamed("nBayes", nBayesStr.c_str()));
  nameList.push_back(new TNamed("nSVD", nSVDStr.c_str()));
  nameList.push_back(new TNamed("nBigBayesSymm", nBigBayesSymmStr.c_str()));
  nameList.push_back(new TNamed("nSuperBayes", nSuperBayesStr.c_str()));
  nameList.push_back(new TNamed("bayesVal", bayesVal2.c_str()));
  nameList.push_back(new TNamed("nHistDim", nHistDimStr.c_str()));
  
  for(auto const & name : nameList){
    name->Write("", TObject::kOverwrite);
  }

  if(doLocalDebug || doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  inFile_p->cd();
  inDir_p->cd();
  inSubDir_p->cd();  


  if(doLocalDebug || doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  std::vector<TNamed> inFullFileNames2;
  for(unsigned int i = 0; i < inFullFileNames.size(); ++i){
    inFullFileNames2.push_back(TNamed(("inFullFileNames_" + std::to_string(i)).c_str(), inFullFileNames.at(i).c_str()));
  }

  if(doLocalDebug || doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  for(unsigned int i = 0; i < inFullFileNames2.size(); ++i){
    inFullFileNames2.at(i).Write("", TObject::kOverwrite);
  }

  if(doLocalDebug || doGlobalDebug) std::cout << __FILE__ << ", " << __LINE__ << std::endl;

  if(unfoldDirBayes_p != NULL){
    std::cout << "Attempting a write of all unfold obj" << std::endl;

    inFile_p->cd();
    inDir_p->cd();

    std::cout << " " << unfoldDirBayes_p->GetName() << std::endl;

    unfoldDirBayes_p->cd();
    //    nHistDimName.Write("", TObject::kOverwrite);
    
    for(Int_t hI = 0; hI < nHistDim; ++hI){
      unfoldDirBayes_p->cd();
      TNamed tempTName(histTagBayes.at(hI).c_str(), std::to_string(histBestBayes.at(hI)).c_str());
      tempTName.Write("", TObject::kOverwrite);
    }
  }

  if(unfoldDirSVD_p != NULL){
    std::cout << "Attempting a write of all unfold obj" << std::endl;

    inFile_p->cd();
    inDir_p->cd();

    std::cout << " " << unfoldDirSVD_p->GetName() << std::endl;

    unfoldDirSVD_p->cd();
    //    nHistDimName.Write("", TObject::kOverwrite);
    
    for(Int_t hI = 0; hI < nHistDim; ++hI){
      unfoldDirSVD_p->cd();
      TNamed tempTName(histTagSVD.at(hI).c_str(), std::to_string(histBestSVD.at(hI)).c_str());
      tempTName.Write("", TObject::kOverwrite);
    }
  }


  return true;
}


bool cutPropagator::CheckPropagatorsMatch(cutPropagator inCutProp, bool doBothMCOrBothData, bool doBothPPOrBothPbPb, bool skipJtAlgo = false)
{
  if(!CheckJtAbsEtaMax(inCutProp)) return false;
  if(!CheckNR(inCutProp)) return false;
  if(!CheckRVals(inCutProp)) return false;


  if(!CheckGenPtBinsR1Cent0to10(inCutProp)) return false;
  if(!CheckRecoPtBinsR1Cent0to10(inCutProp)) return false;
  if(!CheckGenPtBinsR1Cent10to30(inCutProp)) return false;
  if(!CheckRecoPtBinsR1Cent10to30(inCutProp)) return false;
  if(!CheckGenPtBinsR1Cent30to50(inCutProp)) return false;
  if(!CheckRecoPtBinsR1Cent30to50(inCutProp)) return false;
  if(!CheckGenPtBinsR1Cent50to90(inCutProp)) return false;
  if(!CheckRecoPtBinsR1Cent50to90(inCutProp)) return false;

  if(!CheckGenPtBinsR2Cent0to10(inCutProp)) return false;
  if(!CheckRecoPtBinsR2Cent0to10(inCutProp)) return false;
  if(!CheckGenPtBinsR2Cent10to30(inCutProp)) return false;
  if(!CheckRecoPtBinsR2Cent10to30(inCutProp)) return false;
  if(!CheckGenPtBinsR2Cent30to50(inCutProp)) return false;
  if(!CheckRecoPtBinsR2Cent30to50(inCutProp)) return false;
  if(!CheckGenPtBinsR2Cent50to90(inCutProp)) return false;
  if(!CheckRecoPtBinsR2Cent50to90(inCutProp)) return false;

  if(!CheckGenPtBinsR3Cent0to10(inCutProp)) return false;
  if(!CheckRecoPtBinsR3Cent0to10(inCutProp)) return false;
  if(!CheckGenPtBinsR3Cent10to30(inCutProp)) return false;
  if(!CheckRecoPtBinsR3Cent10to30(inCutProp)) return false;
  if(!CheckGenPtBinsR3Cent30to50(inCutProp)) return false;
  if(!CheckRecoPtBinsR3Cent30to50(inCutProp)) return false;
  if(!CheckGenPtBinsR3Cent50to90(inCutProp)) return false;
  if(!CheckRecoPtBinsR3Cent50to90(inCutProp)) return false;

  if(!CheckGenPtBinsR4Cent0to10(inCutProp)) return false;
  if(!CheckRecoPtBinsR4Cent0to10(inCutProp)) return false;
  if(!CheckGenPtBinsR4Cent10to30(inCutProp)) return false;
  if(!CheckRecoPtBinsR4Cent10to30(inCutProp)) return false;
  if(!CheckGenPtBinsR4Cent30to50(inCutProp)) return false;
  if(!CheckRecoPtBinsR4Cent30to50(inCutProp)) return false;
  if(!CheckGenPtBinsR4Cent50to90(inCutProp)) return false;
  if(!CheckRecoPtBinsR4Cent50to90(inCutProp)) return false;

  if(!CheckGenPtBinsR5Cent0to10(inCutProp)) return false;
  if(!CheckRecoPtBinsR5Cent0to10(inCutProp)) return false;
  if(!CheckGenPtBinsR5Cent10to30(inCutProp)) return false;
  if(!CheckRecoPtBinsR5Cent10to30(inCutProp)) return false;
  if(!CheckGenPtBinsR5Cent30to50(inCutProp)) return false;
  if(!CheckRecoPtBinsR5Cent30to50(inCutProp)) return false;
  if(!CheckGenPtBinsR5Cent50to90(inCutProp)) return false;
  if(!CheckRecoPtBinsR5Cent50to90(inCutProp)) return false;

  if(!CheckGenPtBinsR6Cent0to10(inCutProp)) return false;
  if(!CheckRecoPtBinsR6Cent0to10(inCutProp)) return false;
  if(!CheckGenPtBinsR6Cent10to30(inCutProp)) return false;
  if(!CheckRecoPtBinsR6Cent10to30(inCutProp)) return false;
  if(!CheckGenPtBinsR6Cent30to50(inCutProp)) return false;
  if(!CheckRecoPtBinsR6Cent30to50(inCutProp)) return false;
  if(!CheckGenPtBinsR6Cent50to90(inCutProp)) return false;
  if(!CheckRecoPtBinsR6Cent50to90(inCutProp)) return false;

  if(!CheckGenPtBinsR8Cent0to10(inCutProp)) return false;
  if(!CheckRecoPtBinsR8Cent0to10(inCutProp)) return false;
  if(!CheckGenPtBinsR8Cent10to30(inCutProp)) return false;
  if(!CheckRecoPtBinsR8Cent10to30(inCutProp)) return false;
  if(!CheckGenPtBinsR8Cent30to50(inCutProp)) return false;
  if(!CheckRecoPtBinsR8Cent30to50(inCutProp)) return false;
  if(!CheckGenPtBinsR8Cent50to90(inCutProp)) return false;
  if(!CheckRecoPtBinsR8Cent50to90(inCutProp)) return false;

  if(!CheckGenPtBinsR10Cent0to10(inCutProp)) return false;
  if(!CheckRecoPtBinsR10Cent0to10(inCutProp)) return false;
  if(!CheckGenPtBinsR10Cent10to30(inCutProp)) return false;
  if(!CheckRecoPtBinsR10Cent10to30(inCutProp)) return false;
  if(!CheckGenPtBinsR10Cent30to50(inCutProp)) return false;
  if(!CheckRecoPtBinsR10Cent30to50(inCutProp)) return false;
  if(!CheckGenPtBinsR10Cent50to90(inCutProp)) return false;
  if(!CheckRecoPtBinsR10Cent50to90(inCutProp)) return false;

  if(!CheckNGeneralBins(inCutProp)) return false;

  if(!CheckNJtAbsEtaBins(inCutProp)) return false;
  if(!CheckNID(inCutProp)) return false;

  if(!CheckGeneralBins(inCutProp)) return false;

  if(!CheckJtAbsEtaBinsLow(inCutProp)) return false;
  if(!CheckJtAbsEtaBinsHi(inCutProp)) return false;
  if(!CheckIdStr(inCutProp)) return false;
  if(!CheckJtPfCHMFCutLow(inCutProp)) return false;
  if(!CheckJtPfCHMFCutHi(inCutProp)) return false;
  if(!CheckJtPfMUMFCutLow(inCutProp)) return false;
  if(!CheckJtPfMUMFCutHi(inCutProp)) return false;
  if(!CheckJtPfNHFCutLow(inCutProp)) return false;
  if(!CheckJtPfNHFCutHi(inCutProp)) return false;
  if(!CheckJtPfNEFCutLow(inCutProp)) return false;
  if(!CheckJtPfNEFCutHi(inCutProp)) return false;
  if(!CheckJtPfMUFCutLow(inCutProp)) return false;
  if(!CheckJtPfMUFCutHi(inCutProp)) return false;
  if(!CheckJtPfCHFCutLow(inCutProp)) return false;
  if(!CheckJtPfCHFCutHi(inCutProp)) return false;
  if(!CheckJtPfCEFCutLow(inCutProp)) return false;
  if(!CheckJtPfCEFCutHi(inCutProp)) return false;
  if(!CheckJtPfMinMult(inCutProp)) return false;
  if(!CheckJtPfMinChgMult(inCutProp)) return false;

  if(!CheckNBayes(inCutProp)) return false;
  if(!CheckNBigBayesSymm(inCutProp)) return false;
  if(!CheckBayesVal(inCutProp)) return false;
  if(!CheckNSuperBayes(inCutProp)) return false;

  if(doBothMCOrBothData){
    if(!CheckNSyst(inCutProp)) return false;
    if(!CheckSystStr(inCutProp)) return false;
    if(!CheckRCDiffFileName(inCutProp)) return false;
    if(!CheckPriorFlatFileName(inCutProp)) return false;
    if(!CheckJECVarMC(inCutProp)) return false;
    if(!CheckJERVarMC(inCutProp)) return false;
    if(!CheckJECVarData(inCutProp)) return false;
    if(!CheckNResponseMod(inCutProp)) return false;
    if(!CheckResponseMod(inCutProp)) return false;
    if(!CheckJERVarData(inCutProp)) return false;
    if(!CheckNPthats(inCutProp)) return false;
    if(!CheckPthats(inCutProp)) return false;
    if(!CheckPthatWeights(inCutProp)) return false;
  }

  if(doBothPPOrBothPbPb){
    if(skipJtAlgo){
      if(!CheckNJtAlgos(inCutProp)) return false;
      if(!CheckJtAlgos(inCutProp)) return false;
      if(!CheckMinJtPtCut(inCutProp)) return false;
      if(!CheckMultiJtPtCut(inCutProp)) return false;
      if(!CheckRecoTruncPos(inCutProp)) return false;
    }
    if(!CheckIsPP(inCutProp)) return false;
    if(!CheckNCentBins(inCutProp)) return false;
    if(!CheckCentBinsLow(inCutProp)) return false;
    if(!CheckCentBinsHi(inCutProp)) return false;
  }

  return true;
}

bool cutPropagator::CheckDouble(double in1, double in2)
{
  if(in1 + delta < in2) return false;
  else if(in1 - delta > in2) return false;
  return true;
}

bool cutPropagator::CheckInt(int in1, int in2){return in1 == in2;}

bool cutPropagator::CheckBool(bool in1, bool in2){return in1 == in2;}

bool cutPropagator::CheckString(std::string in1, std::string in2){return isStrSame(in1, in2);}

bool cutPropagator::CheckVectDouble(std::vector<double> in1, std::vector<double> in2)
{
  if(in1.size() != in2.size()) return false;
  for(unsigned int i = 0; i < in1.size(); ++i){
    if(!CheckDouble(in1.at(i), in2.at(i))) return false;
  }
  return true;
}

bool cutPropagator::CheckVectDouble(std::vector<double> in1, std::vector<double> in2, std::string valStr)
{
  bool checkVal = CheckVectDouble(in1, in2);
  if(!checkVal) std::cout << "cutPropagator check failed on " << valStr << std::endl;   
  return checkVal;
}

bool cutPropagator::CheckVectInt(std::vector<int> in1, std::vector<int> in2)
{
  if(in1.size() != in2.size()) return false;
  for(unsigned int i = 0; i < in1.size(); ++i){
    if(!CheckInt(in1.at(i), in2.at(i))) return false;
  }
  return true;
}

bool cutPropagator::CheckVectString(std::vector<std::string> in1, std::vector<std::string> in2)
{
  if(in1.size() != in2.size()) return false;
  for(unsigned int i = 0; i < in1.size(); ++i){
    if(!CheckString(in1.at(i), in2.at(i))) return false;
  }
  return true;
}

//Series of checks
bool cutPropagator::CheckJtAbsEtaMax(double inJtAbsEtaMax)
{
  bool checkVal = CheckDouble(jtAbsEtaMax, inJtAbsEtaMax);
  if(!checkVal) std::cout << "cutPropagator check failed on jtAbsEtaMax" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckJtAbsEtaMax(cutPropagator inCutProp){return CheckJtAbsEtaMax(inCutProp.GetJtAbsEtaMax());}

bool cutPropagator::CheckNR(int inNR)
{
  bool checkVal = CheckInt(nR, inNR);
  if(!checkVal) std::cout << "cutPropagator check failed on nR" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckNR(cutPropagator inCutProp){return CheckNR(inCutProp.GetNR());}

bool cutPropagator::CheckRVals(std::vector<int> inRVals)
{
  bool checkVal = CheckVectInt(rVals, inRVals);
  if(!checkVal) std::cout << "cutPropagator check failed on rVals" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckRVals(cutPropagator inCutProp){return CheckRVals(inCutProp.GetRVals());}

bool cutPropagator::CheckNBins(int inBins, int compBins, std::string binStr)
{
  bool checkVal = CheckInt(compBins, inBins);
  if(!checkVal) std::cout << "cutPropagator check failed on " << binStr << std::endl;
  return checkVal;
}

bool cutPropagator::CheckStartEndVals(double inVal, double compVal, std::string valStr)
{
  bool checkVal = CheckDouble(compVal, inVal);
  if(!checkVal) std::cout << "cutPropagator check failed on " << valStr << std::endl;
  return checkVal;
}


bool cutPropagator::CheckGenPtBinsR1Cent0to10(std::vector<double> inGenPtBinsR1Cent0to10){return CheckVectDouble(inGenPtBinsR1Cent0to10, genJtPtBinsR1Cent0to10, "genJtPtBinsR1Cent0to10");}
bool cutPropagator::CheckGenPtBinsR1Cent0to10(cutPropagator inCutProp){return CheckGenPtBinsR1Cent0to10(inCutProp.GetGenPtBinsR1Cent0to10());}
bool cutPropagator::CheckRecoPtBinsR1Cent0to10(std::vector<double> inRecoPtBinsR1Cent0to10){return CheckVectDouble(inRecoPtBinsR1Cent0to10, recoJtPtBinsR1Cent0to10, "recoJtPtBinsR1Cent0to10");}
bool cutPropagator::CheckRecoPtBinsR1Cent0to10(cutPropagator inCutProp){return CheckRecoPtBinsR1Cent0to10(inCutProp.GetRecoPtBinsR1Cent0to10());}
bool cutPropagator::CheckGenPtBinsR1Cent10to30(std::vector<double> inGenPtBinsR1Cent10to30){return CheckVectDouble(inGenPtBinsR1Cent10to30, genJtPtBinsR1Cent10to30, "genJtPtBinsR1Cent10to30");}
bool cutPropagator::CheckGenPtBinsR1Cent10to30(cutPropagator inCutProp){return CheckGenPtBinsR1Cent10to30(inCutProp.GetGenPtBinsR1Cent10to30());}
bool cutPropagator::CheckRecoPtBinsR1Cent10to30(std::vector<double> inRecoPtBinsR1Cent10to30){return CheckVectDouble(inRecoPtBinsR1Cent10to30, recoJtPtBinsR1Cent10to30, "recoJtPtBinsR1Cent10to30");}
bool cutPropagator::CheckRecoPtBinsR1Cent10to30(cutPropagator inCutProp){return CheckRecoPtBinsR1Cent10to30(inCutProp.GetRecoPtBinsR1Cent10to30());}
bool cutPropagator::CheckGenPtBinsR1Cent30to50(std::vector<double> inGenPtBinsR1Cent30to50){return CheckVectDouble(inGenPtBinsR1Cent30to50, genJtPtBinsR1Cent30to50, "genJtPtBinsR1Cent30to50");}
bool cutPropagator::CheckGenPtBinsR1Cent30to50(cutPropagator inCutProp){return CheckGenPtBinsR1Cent30to50(inCutProp.GetGenPtBinsR1Cent30to50());}
bool cutPropagator::CheckRecoPtBinsR1Cent30to50(std::vector<double> inRecoPtBinsR1Cent30to50){return CheckVectDouble(inRecoPtBinsR1Cent30to50, recoJtPtBinsR1Cent30to50, "recoJtPtBinsR1Cent30to50");}
bool cutPropagator::CheckRecoPtBinsR1Cent30to50(cutPropagator inCutProp){return CheckRecoPtBinsR1Cent30to50(inCutProp.GetRecoPtBinsR1Cent30to50());}
bool cutPropagator::CheckGenPtBinsR1Cent50to90(std::vector<double> inGenPtBinsR1Cent50to90){return CheckVectDouble(inGenPtBinsR1Cent50to90, genJtPtBinsR1Cent50to90, "genJtPtBinsR1Cent50to90");}
bool cutPropagator::CheckGenPtBinsR1Cent50to90(cutPropagator inCutProp){return CheckGenPtBinsR1Cent50to90(inCutProp.GetGenPtBinsR1Cent50to90());}
bool cutPropagator::CheckRecoPtBinsR1Cent50to90(std::vector<double> inRecoPtBinsR1Cent50to90){return CheckVectDouble(inRecoPtBinsR1Cent50to90, recoJtPtBinsR1Cent50to90, "recoJtPtBinsR1Cent50to90");}
bool cutPropagator::CheckRecoPtBinsR1Cent50to90(cutPropagator inCutProp){return CheckRecoPtBinsR1Cent50to90(inCutProp.GetRecoPtBinsR1Cent50to90());}

bool cutPropagator::CheckGenPtBinsR2Cent0to10(std::vector<double> inGenPtBinsR2Cent0to10){return CheckVectDouble(inGenPtBinsR2Cent0to10, genJtPtBinsR2Cent0to10, "genJtPtBinsR2Cent0to10");}
bool cutPropagator::CheckGenPtBinsR2Cent0to10(cutPropagator inCutProp){return CheckGenPtBinsR2Cent0to10(inCutProp.GetGenPtBinsR2Cent0to10());}
bool cutPropagator::CheckRecoPtBinsR2Cent0to10(std::vector<double> inRecoPtBinsR2Cent0to10){return CheckVectDouble(inRecoPtBinsR2Cent0to10, recoJtPtBinsR2Cent0to10, "recoJtPtBinsR2Cent0to10");}
bool cutPropagator::CheckRecoPtBinsR2Cent0to10(cutPropagator inCutProp){return CheckRecoPtBinsR2Cent0to10(inCutProp.GetRecoPtBinsR2Cent0to10());}
bool cutPropagator::CheckGenPtBinsR2Cent10to30(std::vector<double> inGenPtBinsR2Cent10to30){return CheckVectDouble(inGenPtBinsR2Cent10to30, genJtPtBinsR2Cent10to30, "genJtPtBinsR2Cent10to30");}
bool cutPropagator::CheckGenPtBinsR2Cent10to30(cutPropagator inCutProp){return CheckGenPtBinsR2Cent10to30(inCutProp.GetGenPtBinsR2Cent10to30());}
bool cutPropagator::CheckRecoPtBinsR2Cent10to30(std::vector<double> inRecoPtBinsR2Cent10to30){return CheckVectDouble(inRecoPtBinsR2Cent10to30, recoJtPtBinsR2Cent10to30, "recoJtPtBinsR2Cent10to30");}
bool cutPropagator::CheckRecoPtBinsR2Cent10to30(cutPropagator inCutProp){return CheckRecoPtBinsR2Cent10to30(inCutProp.GetRecoPtBinsR2Cent10to30());}
bool cutPropagator::CheckGenPtBinsR2Cent30to50(std::vector<double> inGenPtBinsR2Cent30to50){return CheckVectDouble(inGenPtBinsR2Cent30to50, genJtPtBinsR2Cent30to50, "genJtPtBinsR2Cent30to50");}
bool cutPropagator::CheckGenPtBinsR2Cent30to50(cutPropagator inCutProp){return CheckGenPtBinsR2Cent30to50(inCutProp.GetGenPtBinsR2Cent30to50());}
bool cutPropagator::CheckRecoPtBinsR2Cent30to50(std::vector<double> inRecoPtBinsR2Cent30to50){return CheckVectDouble(inRecoPtBinsR2Cent30to50, recoJtPtBinsR2Cent30to50, "recoJtPtBinsR2Cent30to50");}
bool cutPropagator::CheckRecoPtBinsR2Cent30to50(cutPropagator inCutProp){return CheckRecoPtBinsR2Cent30to50(inCutProp.GetRecoPtBinsR2Cent30to50());}
bool cutPropagator::CheckGenPtBinsR2Cent50to90(std::vector<double> inGenPtBinsR2Cent50to90){return CheckVectDouble(inGenPtBinsR2Cent50to90, genJtPtBinsR2Cent50to90, "genJtPtBinsR2Cent50to90");}
bool cutPropagator::CheckGenPtBinsR2Cent50to90(cutPropagator inCutProp){return CheckGenPtBinsR2Cent50to90(inCutProp.GetGenPtBinsR2Cent50to90());}
bool cutPropagator::CheckRecoPtBinsR2Cent50to90(std::vector<double> inRecoPtBinsR2Cent50to90){return CheckVectDouble(inRecoPtBinsR2Cent50to90, recoJtPtBinsR2Cent50to90, "recoJtPtBinsR2Cent50to90");}
bool cutPropagator::CheckRecoPtBinsR2Cent50to90(cutPropagator inCutProp){return CheckRecoPtBinsR2Cent50to90(inCutProp.GetRecoPtBinsR2Cent50to90());}

bool cutPropagator::CheckGenPtBinsR3Cent0to10(std::vector<double> inGenPtBinsR3Cent0to10){return CheckVectDouble(inGenPtBinsR3Cent0to10, genJtPtBinsR3Cent0to10, "genJtPtBinsR3Cent0to10");}
bool cutPropagator::CheckGenPtBinsR3Cent0to10(cutPropagator inCutProp){return CheckGenPtBinsR3Cent0to10(inCutProp.GetGenPtBinsR3Cent0to10());}
bool cutPropagator::CheckRecoPtBinsR3Cent0to10(std::vector<double> inRecoPtBinsR3Cent0to10){return CheckVectDouble(inRecoPtBinsR3Cent0to10, recoJtPtBinsR3Cent0to10, "recoJtPtBinsR3Cent0to10");}
bool cutPropagator::CheckRecoPtBinsR3Cent0to10(cutPropagator inCutProp){return CheckRecoPtBinsR3Cent0to10(inCutProp.GetRecoPtBinsR3Cent0to10());}
bool cutPropagator::CheckGenPtBinsR3Cent10to30(std::vector<double> inGenPtBinsR3Cent10to30){return CheckVectDouble(inGenPtBinsR3Cent10to30, genJtPtBinsR3Cent10to30, "genJtPtBinsR3Cent10to30");}
bool cutPropagator::CheckGenPtBinsR3Cent10to30(cutPropagator inCutProp){return CheckGenPtBinsR3Cent10to30(inCutProp.GetGenPtBinsR3Cent10to30());}
bool cutPropagator::CheckRecoPtBinsR3Cent10to30(std::vector<double> inRecoPtBinsR3Cent10to30){return CheckVectDouble(inRecoPtBinsR3Cent10to30, recoJtPtBinsR3Cent10to30, "recoJtPtBinsR3Cent10to30");}
bool cutPropagator::CheckRecoPtBinsR3Cent10to30(cutPropagator inCutProp){return CheckRecoPtBinsR3Cent10to30(inCutProp.GetRecoPtBinsR3Cent10to30());}
bool cutPropagator::CheckGenPtBinsR3Cent30to50(std::vector<double> inGenPtBinsR3Cent30to50){return CheckVectDouble(inGenPtBinsR3Cent30to50, genJtPtBinsR3Cent30to50, "genJtPtBinsR3Cent30to50");}
bool cutPropagator::CheckGenPtBinsR3Cent30to50(cutPropagator inCutProp){return CheckGenPtBinsR3Cent30to50(inCutProp.GetGenPtBinsR3Cent30to50());}
bool cutPropagator::CheckRecoPtBinsR3Cent30to50(std::vector<double> inRecoPtBinsR3Cent30to50){return CheckVectDouble(inRecoPtBinsR3Cent30to50, recoJtPtBinsR3Cent30to50, "recoJtPtBinsR3Cent30to50");}
bool cutPropagator::CheckRecoPtBinsR3Cent30to50(cutPropagator inCutProp){return CheckRecoPtBinsR3Cent30to50(inCutProp.GetRecoPtBinsR3Cent30to50());}
bool cutPropagator::CheckGenPtBinsR3Cent50to90(std::vector<double> inGenPtBinsR3Cent50to90){return CheckVectDouble(inGenPtBinsR3Cent50to90, genJtPtBinsR3Cent50to90, "genJtPtBinsR3Cent50to90");}
bool cutPropagator::CheckGenPtBinsR3Cent50to90(cutPropagator inCutProp){return CheckGenPtBinsR3Cent50to90(inCutProp.GetGenPtBinsR3Cent50to90());}
bool cutPropagator::CheckRecoPtBinsR3Cent50to90(std::vector<double> inRecoPtBinsR3Cent50to90){return CheckVectDouble(inRecoPtBinsR3Cent50to90, recoJtPtBinsR3Cent50to90, "recoJtPtBinsR3Cent50to90");}
bool cutPropagator::CheckRecoPtBinsR3Cent50to90(cutPropagator inCutProp){return CheckRecoPtBinsR3Cent50to90(inCutProp.GetRecoPtBinsR3Cent50to90());}

bool cutPropagator::CheckGenPtBinsR4Cent0to10(std::vector<double> inGenPtBinsR4Cent0to10){return CheckVectDouble(inGenPtBinsR4Cent0to10, genJtPtBinsR4Cent0to10, "genJtPtBinsR4Cent0to10");}
bool cutPropagator::CheckGenPtBinsR4Cent0to10(cutPropagator inCutProp){return CheckGenPtBinsR4Cent0to10(inCutProp.GetGenPtBinsR4Cent0to10());}
bool cutPropagator::CheckRecoPtBinsR4Cent0to10(std::vector<double> inRecoPtBinsR4Cent0to10){return CheckVectDouble(inRecoPtBinsR4Cent0to10, recoJtPtBinsR4Cent0to10, "recoJtPtBinsR4Cent0to10");}
bool cutPropagator::CheckRecoPtBinsR4Cent0to10(cutPropagator inCutProp){return CheckRecoPtBinsR4Cent0to10(inCutProp.GetRecoPtBinsR4Cent0to10());}
bool cutPropagator::CheckGenPtBinsR4Cent10to30(std::vector<double> inGenPtBinsR4Cent10to30){return CheckVectDouble(inGenPtBinsR4Cent10to30, genJtPtBinsR4Cent10to30, "genJtPtBinsR4Cent10to30");}
bool cutPropagator::CheckGenPtBinsR4Cent10to30(cutPropagator inCutProp){return CheckGenPtBinsR4Cent10to30(inCutProp.GetGenPtBinsR4Cent10to30());}
bool cutPropagator::CheckRecoPtBinsR4Cent10to30(std::vector<double> inRecoPtBinsR4Cent10to30){return CheckVectDouble(inRecoPtBinsR4Cent10to30, recoJtPtBinsR4Cent10to30, "recoJtPtBinsR4Cent10to30");}
bool cutPropagator::CheckRecoPtBinsR4Cent10to30(cutPropagator inCutProp){return CheckRecoPtBinsR4Cent10to30(inCutProp.GetRecoPtBinsR4Cent10to30());}
bool cutPropagator::CheckGenPtBinsR4Cent30to50(std::vector<double> inGenPtBinsR4Cent30to50){return CheckVectDouble(inGenPtBinsR4Cent30to50, genJtPtBinsR4Cent30to50, "genJtPtBinsR4Cent30to50");}
bool cutPropagator::CheckGenPtBinsR4Cent30to50(cutPropagator inCutProp){return CheckGenPtBinsR4Cent30to50(inCutProp.GetGenPtBinsR4Cent30to50());}
bool cutPropagator::CheckRecoPtBinsR4Cent30to50(std::vector<double> inRecoPtBinsR4Cent30to50){return CheckVectDouble(inRecoPtBinsR4Cent30to50, recoJtPtBinsR4Cent30to50, "recoJtPtBinsR4Cent30to50");}
bool cutPropagator::CheckRecoPtBinsR4Cent30to50(cutPropagator inCutProp){return CheckRecoPtBinsR4Cent30to50(inCutProp.GetRecoPtBinsR4Cent30to50());}
bool cutPropagator::CheckGenPtBinsR4Cent50to90(std::vector<double> inGenPtBinsR4Cent50to90){return CheckVectDouble(inGenPtBinsR4Cent50to90, genJtPtBinsR4Cent50to90, "genJtPtBinsR4Cent50to90");}
bool cutPropagator::CheckGenPtBinsR4Cent50to90(cutPropagator inCutProp){return CheckGenPtBinsR4Cent50to90(inCutProp.GetGenPtBinsR4Cent50to90());}
bool cutPropagator::CheckRecoPtBinsR4Cent50to90(std::vector<double> inRecoPtBinsR4Cent50to90){return CheckVectDouble(inRecoPtBinsR4Cent50to90, recoJtPtBinsR4Cent50to90, "recoJtPtBinsR4Cent50to90");}
bool cutPropagator::CheckRecoPtBinsR4Cent50to90(cutPropagator inCutProp){return CheckRecoPtBinsR4Cent50to90(inCutProp.GetRecoPtBinsR4Cent50to90());}

bool cutPropagator::CheckGenPtBinsR5Cent0to10(std::vector<double> inGenPtBinsR5Cent0to10){return CheckVectDouble(inGenPtBinsR5Cent0to10, genJtPtBinsR5Cent0to10, "genJtPtBinsR5Cent0to10");}
bool cutPropagator::CheckGenPtBinsR5Cent0to10(cutPropagator inCutProp){return CheckGenPtBinsR5Cent0to10(inCutProp.GetGenPtBinsR5Cent0to10());}
bool cutPropagator::CheckRecoPtBinsR5Cent0to10(std::vector<double> inRecoPtBinsR5Cent0to10){return CheckVectDouble(inRecoPtBinsR5Cent0to10, recoJtPtBinsR5Cent0to10, "recoJtPtBinsR5Cent0to10");}
bool cutPropagator::CheckRecoPtBinsR5Cent0to10(cutPropagator inCutProp){return CheckRecoPtBinsR5Cent0to10(inCutProp.GetRecoPtBinsR5Cent0to10());}
bool cutPropagator::CheckGenPtBinsR5Cent10to30(std::vector<double> inGenPtBinsR5Cent10to30){return CheckVectDouble(inGenPtBinsR5Cent10to30, genJtPtBinsR5Cent10to30, "genJtPtBinsR5Cent10to30");}
bool cutPropagator::CheckGenPtBinsR5Cent10to30(cutPropagator inCutProp){return CheckGenPtBinsR5Cent10to30(inCutProp.GetGenPtBinsR5Cent10to30());}
bool cutPropagator::CheckRecoPtBinsR5Cent10to30(std::vector<double> inRecoPtBinsR5Cent10to30){return CheckVectDouble(inRecoPtBinsR5Cent10to30, recoJtPtBinsR5Cent10to30, "recoJtPtBinsR5Cent10to30");}
bool cutPropagator::CheckRecoPtBinsR5Cent10to30(cutPropagator inCutProp){return CheckRecoPtBinsR5Cent10to30(inCutProp.GetRecoPtBinsR5Cent10to30());}
bool cutPropagator::CheckGenPtBinsR5Cent30to50(std::vector<double> inGenPtBinsR5Cent30to50){return CheckVectDouble(inGenPtBinsR5Cent30to50, genJtPtBinsR5Cent30to50, "genJtPtBinsR5Cent30to50");}
bool cutPropagator::CheckGenPtBinsR5Cent30to50(cutPropagator inCutProp){return CheckGenPtBinsR5Cent30to50(inCutProp.GetGenPtBinsR5Cent30to50());}
bool cutPropagator::CheckRecoPtBinsR5Cent30to50(std::vector<double> inRecoPtBinsR5Cent30to50){return CheckVectDouble(inRecoPtBinsR5Cent30to50, recoJtPtBinsR5Cent30to50, "recoJtPtBinsR5Cent30to50");}
bool cutPropagator::CheckRecoPtBinsR5Cent30to50(cutPropagator inCutProp){return CheckRecoPtBinsR5Cent30to50(inCutProp.GetRecoPtBinsR5Cent30to50());}
bool cutPropagator::CheckGenPtBinsR5Cent50to90(std::vector<double> inGenPtBinsR5Cent50to90){return CheckVectDouble(inGenPtBinsR5Cent50to90, genJtPtBinsR5Cent50to90, "genJtPtBinsR5Cent50to90");}
bool cutPropagator::CheckGenPtBinsR5Cent50to90(cutPropagator inCutProp){return CheckGenPtBinsR5Cent50to90(inCutProp.GetGenPtBinsR5Cent50to90());}
bool cutPropagator::CheckRecoPtBinsR5Cent50to90(std::vector<double> inRecoPtBinsR5Cent50to90){return CheckVectDouble(inRecoPtBinsR5Cent50to90, recoJtPtBinsR5Cent50to90, "recoJtPtBinsR5Cent50to90");}
bool cutPropagator::CheckRecoPtBinsR5Cent50to90(cutPropagator inCutProp){return CheckRecoPtBinsR5Cent50to90(inCutProp.GetRecoPtBinsR5Cent50to90());}

bool cutPropagator::CheckGenPtBinsR6Cent0to10(std::vector<double> inGenPtBinsR6Cent0to10){return CheckVectDouble(inGenPtBinsR6Cent0to10, genJtPtBinsR6Cent0to10, "genJtPtBinsR6Cent0to10");}
bool cutPropagator::CheckGenPtBinsR6Cent0to10(cutPropagator inCutProp){return CheckGenPtBinsR6Cent0to10(inCutProp.GetGenPtBinsR6Cent0to10());}
bool cutPropagator::CheckRecoPtBinsR6Cent0to10(std::vector<double> inRecoPtBinsR6Cent0to10){return CheckVectDouble(inRecoPtBinsR6Cent0to10, recoJtPtBinsR6Cent0to10, "recoJtPtBinsR6Cent0to10");}
bool cutPropagator::CheckRecoPtBinsR6Cent0to10(cutPropagator inCutProp){return CheckRecoPtBinsR6Cent0to10(inCutProp.GetRecoPtBinsR6Cent0to10());}
bool cutPropagator::CheckGenPtBinsR6Cent10to30(std::vector<double> inGenPtBinsR6Cent10to30){return CheckVectDouble(inGenPtBinsR6Cent10to30, genJtPtBinsR6Cent10to30, "genJtPtBinsR6Cent10to30");}
bool cutPropagator::CheckGenPtBinsR6Cent10to30(cutPropagator inCutProp){return CheckGenPtBinsR6Cent10to30(inCutProp.GetGenPtBinsR6Cent10to30());}
bool cutPropagator::CheckRecoPtBinsR6Cent10to30(std::vector<double> inRecoPtBinsR6Cent10to30){return CheckVectDouble(inRecoPtBinsR6Cent10to30, recoJtPtBinsR6Cent10to30, "recoJtPtBinsR6Cent10to30");}
bool cutPropagator::CheckRecoPtBinsR6Cent10to30(cutPropagator inCutProp){return CheckRecoPtBinsR6Cent10to30(inCutProp.GetRecoPtBinsR6Cent10to30());}
bool cutPropagator::CheckGenPtBinsR6Cent30to50(std::vector<double> inGenPtBinsR6Cent30to50){return CheckVectDouble(inGenPtBinsR6Cent30to50, genJtPtBinsR6Cent30to50, "genJtPtBinsR6Cent30to50");}
bool cutPropagator::CheckGenPtBinsR6Cent30to50(cutPropagator inCutProp){return CheckGenPtBinsR6Cent30to50(inCutProp.GetGenPtBinsR6Cent30to50());}
bool cutPropagator::CheckRecoPtBinsR6Cent30to50(std::vector<double> inRecoPtBinsR6Cent30to50){return CheckVectDouble(inRecoPtBinsR6Cent30to50, recoJtPtBinsR6Cent30to50, "recoJtPtBinsR6Cent30to50");}
bool cutPropagator::CheckRecoPtBinsR6Cent30to50(cutPropagator inCutProp){return CheckRecoPtBinsR6Cent30to50(inCutProp.GetRecoPtBinsR6Cent30to50());}
bool cutPropagator::CheckGenPtBinsR6Cent50to90(std::vector<double> inGenPtBinsR6Cent50to90){return CheckVectDouble(inGenPtBinsR6Cent50to90, genJtPtBinsR6Cent50to90, "genJtPtBinsR6Cent50to90");}
bool cutPropagator::CheckGenPtBinsR6Cent50to90(cutPropagator inCutProp){return CheckGenPtBinsR6Cent50to90(inCutProp.GetGenPtBinsR6Cent50to90());}
bool cutPropagator::CheckRecoPtBinsR6Cent50to90(std::vector<double> inRecoPtBinsR6Cent50to90){return CheckVectDouble(inRecoPtBinsR6Cent50to90, recoJtPtBinsR6Cent50to90, "recoJtPtBinsR6Cent50to90");}
bool cutPropagator::CheckRecoPtBinsR6Cent50to90(cutPropagator inCutProp){return CheckRecoPtBinsR6Cent50to90(inCutProp.GetRecoPtBinsR6Cent50to90());}

bool cutPropagator::CheckGenPtBinsR8Cent0to10(std::vector<double> inGenPtBinsR8Cent0to10){return CheckVectDouble(inGenPtBinsR8Cent0to10, genJtPtBinsR8Cent0to10, "genJtPtBinsR8Cent0to10");}
bool cutPropagator::CheckGenPtBinsR8Cent0to10(cutPropagator inCutProp){return CheckGenPtBinsR8Cent0to10(inCutProp.GetGenPtBinsR8Cent0to10());}
bool cutPropagator::CheckRecoPtBinsR8Cent0to10(std::vector<double> inRecoPtBinsR8Cent0to10){return CheckVectDouble(inRecoPtBinsR8Cent0to10, recoJtPtBinsR8Cent0to10, "recoJtPtBinsR8Cent0to10");}
bool cutPropagator::CheckRecoPtBinsR8Cent0to10(cutPropagator inCutProp){return CheckRecoPtBinsR8Cent0to10(inCutProp.GetRecoPtBinsR8Cent0to10());}
bool cutPropagator::CheckGenPtBinsR8Cent10to30(std::vector<double> inGenPtBinsR8Cent10to30){return CheckVectDouble(inGenPtBinsR8Cent10to30, genJtPtBinsR8Cent10to30, "genJtPtBinsR8Cent10to30");}
bool cutPropagator::CheckGenPtBinsR8Cent10to30(cutPropagator inCutProp){return CheckGenPtBinsR8Cent10to30(inCutProp.GetGenPtBinsR8Cent10to30());}
bool cutPropagator::CheckRecoPtBinsR8Cent10to30(std::vector<double> inRecoPtBinsR8Cent10to30){return CheckVectDouble(inRecoPtBinsR8Cent10to30, recoJtPtBinsR8Cent10to30, "recoJtPtBinsR8Cent10to30");}
bool cutPropagator::CheckRecoPtBinsR8Cent10to30(cutPropagator inCutProp){return CheckRecoPtBinsR8Cent10to30(inCutProp.GetRecoPtBinsR8Cent10to30());}
bool cutPropagator::CheckGenPtBinsR8Cent30to50(std::vector<double> inGenPtBinsR8Cent30to50){return CheckVectDouble(inGenPtBinsR8Cent30to50, genJtPtBinsR8Cent30to50, "genJtPtBinsR8Cent30to50");}
bool cutPropagator::CheckGenPtBinsR8Cent30to50(cutPropagator inCutProp){return CheckGenPtBinsR8Cent30to50(inCutProp.GetGenPtBinsR8Cent30to50());}
bool cutPropagator::CheckRecoPtBinsR8Cent30to50(std::vector<double> inRecoPtBinsR8Cent30to50){return CheckVectDouble(inRecoPtBinsR8Cent30to50, recoJtPtBinsR8Cent30to50, "recoJtPtBinsR8Cent30to50");}
bool cutPropagator::CheckRecoPtBinsR8Cent30to50(cutPropagator inCutProp){return CheckRecoPtBinsR8Cent30to50(inCutProp.GetRecoPtBinsR8Cent30to50());}
bool cutPropagator::CheckGenPtBinsR8Cent50to90(std::vector<double> inGenPtBinsR8Cent50to90){return CheckVectDouble(inGenPtBinsR8Cent50to90, genJtPtBinsR8Cent50to90, "genJtPtBinsR8Cent50to90");}
bool cutPropagator::CheckGenPtBinsR8Cent50to90(cutPropagator inCutProp){return CheckGenPtBinsR8Cent50to90(inCutProp.GetGenPtBinsR8Cent50to90());}
bool cutPropagator::CheckRecoPtBinsR8Cent50to90(std::vector<double> inRecoPtBinsR8Cent50to90){return CheckVectDouble(inRecoPtBinsR8Cent50to90, recoJtPtBinsR8Cent50to90, "recoJtPtBinsR8Cent50to90");}
bool cutPropagator::CheckRecoPtBinsR8Cent50to90(cutPropagator inCutProp){return CheckRecoPtBinsR8Cent50to90(inCutProp.GetRecoPtBinsR8Cent50to90());}

bool cutPropagator::CheckGenPtBinsR10Cent0to10(std::vector<double> inGenPtBinsR10Cent0to10){return CheckVectDouble(inGenPtBinsR10Cent0to10, genJtPtBinsR10Cent0to10, "genJtPtBinsR10Cent0to10");}
bool cutPropagator::CheckGenPtBinsR10Cent0to10(cutPropagator inCutProp){return CheckGenPtBinsR10Cent0to10(inCutProp.GetGenPtBinsR10Cent0to10());}
bool cutPropagator::CheckRecoPtBinsR10Cent0to10(std::vector<double> inRecoPtBinsR10Cent0to10){return CheckVectDouble(inRecoPtBinsR10Cent0to10, recoJtPtBinsR10Cent0to10, "recoJtPtBinsR10Cent0to10");}
bool cutPropagator::CheckRecoPtBinsR10Cent0to10(cutPropagator inCutProp){return CheckRecoPtBinsR10Cent0to10(inCutProp.GetRecoPtBinsR10Cent0to10());}
bool cutPropagator::CheckGenPtBinsR10Cent10to30(std::vector<double> inGenPtBinsR10Cent10to30){return CheckVectDouble(inGenPtBinsR10Cent10to30, genJtPtBinsR10Cent10to30, "genJtPtBinsR10Cent10to30");}
bool cutPropagator::CheckGenPtBinsR10Cent10to30(cutPropagator inCutProp){return CheckGenPtBinsR10Cent10to30(inCutProp.GetGenPtBinsR10Cent10to30());}
bool cutPropagator::CheckRecoPtBinsR10Cent10to30(std::vector<double> inRecoPtBinsR10Cent10to30){return CheckVectDouble(inRecoPtBinsR10Cent10to30, recoJtPtBinsR10Cent10to30, "recoJtPtBinsR10Cent10to30");}
bool cutPropagator::CheckRecoPtBinsR10Cent10to30(cutPropagator inCutProp){return CheckRecoPtBinsR10Cent10to30(inCutProp.GetRecoPtBinsR10Cent10to30());}
bool cutPropagator::CheckGenPtBinsR10Cent30to50(std::vector<double> inGenPtBinsR10Cent30to50){return CheckVectDouble(inGenPtBinsR10Cent30to50, genJtPtBinsR10Cent30to50, "genJtPtBinsR10Cent30to50");}
bool cutPropagator::CheckGenPtBinsR10Cent30to50(cutPropagator inCutProp){return CheckGenPtBinsR10Cent30to50(inCutProp.GetGenPtBinsR10Cent30to50());}
bool cutPropagator::CheckRecoPtBinsR10Cent30to50(std::vector<double> inRecoPtBinsR10Cent30to50){return CheckVectDouble(inRecoPtBinsR10Cent30to50, recoJtPtBinsR10Cent30to50, "recoJtPtBinsR10Cent30to50");}
bool cutPropagator::CheckRecoPtBinsR10Cent30to50(cutPropagator inCutProp){return CheckRecoPtBinsR10Cent30to50(inCutProp.GetRecoPtBinsR10Cent30to50());}
bool cutPropagator::CheckGenPtBinsR10Cent50to90(std::vector<double> inGenPtBinsR10Cent50to90){return CheckVectDouble(inGenPtBinsR10Cent50to90, genJtPtBinsR10Cent50to90, "genJtPtBinsR10Cent50to90");}
bool cutPropagator::CheckGenPtBinsR10Cent50to90(cutPropagator inCutProp){return CheckGenPtBinsR10Cent50to90(inCutProp.GetGenPtBinsR10Cent50to90());}
bool cutPropagator::CheckRecoPtBinsR10Cent50to90(std::vector<double> inRecoPtBinsR10Cent50to90){return CheckVectDouble(inRecoPtBinsR10Cent50to90, recoJtPtBinsR10Cent50to90, "recoJtPtBinsR10Cent50to90");}
bool cutPropagator::CheckRecoPtBinsR10Cent50to90(cutPropagator inCutProp){return CheckRecoPtBinsR10Cent50to90(inCutProp.GetRecoPtBinsR10Cent50to90());}



bool cutPropagator::CheckNGeneralBins(int inNGeneralBins){return CheckNBins(inNGeneralBins, nGeneralBins, "nGeneralBins");}
bool cutPropagator::CheckNGeneralBins(cutPropagator inCutProp){return CheckNGeneralBins(inCutProp.GetNGeneralBins());}


bool cutPropagator::CheckNJtAbsEtaBins(int inNJtAbsEtaBins)
{
  bool checkVal = CheckInt(nJtAbsEtaBins, inNJtAbsEtaBins);
  if(!checkVal) std::cout << "cutPropagator check failed on nJtAbsEtaBins" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckNJtAbsEtaBins(cutPropagator inCutProp){return CheckNJtAbsEtaBins(inCutProp.GetNJtAbsEtaBins());}


bool cutPropagator::CheckNID(int inNID)
{
  bool checkVal = CheckInt(nID, inNID);
  if(!checkVal) std::cout << "cutPropagator check failed on nID" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckNID(cutPropagator inCutProp){return CheckNID(inCutProp.GetNID());}


bool cutPropagator::CheckGeneralBins(std::vector<double> inGeneralBins)
{
  bool checkVal = CheckVectDouble(generalBins, inGeneralBins);
  if(!checkVal) std::cout << "cutPropagator check failed on generalBins" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckGeneralBins(cutPropagator inCutProp){return CheckGeneralBins(inCutProp.GetGeneralBins());}


bool cutPropagator::CheckJtAbsEtaBinsLow(std::vector<double> inJtAbsEtaBinsLow)
{
  bool checkVal = CheckVectDouble(jtAbsEtaBinsLow, inJtAbsEtaBinsLow);
  if(!checkVal) std::cout << "cutPropagator check failed on jtAbsEtaBinsLow" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckJtAbsEtaBinsLow(cutPropagator inCutProp){return CheckJtAbsEtaBinsLow(inCutProp.GetJtAbsEtaBinsLow());}

bool cutPropagator::CheckJtAbsEtaBinsHi(std::vector<double> inJtAbsEtaBinsHi)
{
  bool checkVal = CheckVectDouble(jtAbsEtaBinsHi, inJtAbsEtaBinsHi);
  if(!checkVal) std::cout << "cutPropagator check failed on jtAbsEtaBinsHi" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckJtAbsEtaBinsHi(cutPropagator inCutProp){return CheckJtAbsEtaBinsHi(inCutProp.GetJtAbsEtaBinsHi());}

bool cutPropagator::CheckIdStr(std::vector<std::string> inIDStr)
{
  bool checkVal = CheckVectString(idStr, inIDStr);
  if(!checkVal) std::cout << "cutPropagator check failed on idStr" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckIdStr(cutPropagator inCutProp){return CheckIdStr(inCutProp.GetIdStr());}

bool cutPropagator::CheckJtPfCHMFCutLow(std::vector<double> inJtPfCHMFCutLow)
{
  bool checkVal = CheckVectDouble(jtPfCHMFCutLow, inJtPfCHMFCutLow);
  if(!checkVal) std::cout << "cutPropagator check failed on jtPfCHMFCutLow" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckJtPfCHMFCutLow(cutPropagator inCutProp){return CheckJtPfCHMFCutLow(inCutProp.GetJtPfCHMFCutLow());}

bool cutPropagator::CheckJtPfCHMFCutHi(std::vector<double> inJtPfCHMFCutHi)
{
  bool checkVal = CheckVectDouble(jtPfCHMFCutHi, inJtPfCHMFCutHi);
  if(!checkVal) std::cout << "cutPropagator check failed on jtPfCHMFCutHi" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckJtPfCHMFCutHi(cutPropagator inCutProp){return CheckJtPfCHMFCutHi(inCutProp.GetJtPfCHMFCutHi());}

bool cutPropagator::CheckJtPfMUMFCutLow(std::vector<double> inJtPfMUMFCutLow)
{
  bool checkVal = CheckVectDouble(jtPfMUMFCutLow, inJtPfMUMFCutLow);
  if(!checkVal) std::cout << "cutPropagator check failed on jtPfMUMFCutLow" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckJtPfMUMFCutLow(cutPropagator inCutProp){return CheckJtPfMUMFCutLow(inCutProp.GetJtPfMUMFCutLow());}

bool cutPropagator::CheckJtPfMUMFCutHi(std::vector<double> inJtPfMUMFCutHi)
{
  bool checkVal = CheckVectDouble(jtPfMUMFCutHi, inJtPfMUMFCutHi);
  if(!checkVal) std::cout << "cutPropagator check failed on jtPfMUMFCutHi" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckJtPfMUMFCutHi(cutPropagator inCutProp){return CheckJtPfMUMFCutHi(inCutProp.GetJtPfMUMFCutHi());}

bool cutPropagator::CheckJtPfNHFCutLow(std::vector<double> inJtPfNHFCutLow)
{
  bool checkVal = CheckVectDouble(jtPfNHFCutLow, inJtPfNHFCutLow);
  if(!checkVal) std::cout << "cutPropagator check failed on jtPfNHFCutLow" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckJtPfNHFCutLow(cutPropagator inCutProp){return CheckJtPfNHFCutLow(inCutProp.GetJtPfNHFCutLow());}

bool cutPropagator::CheckJtPfNHFCutHi(std::vector<double> inJtPfNHFCutHi)
{
  bool checkVal = CheckVectDouble(jtPfNHFCutHi, inJtPfNHFCutHi);
  if(!checkVal) std::cout << "cutPropagator check failed on jtPfNHFCutHi" << std::endl;
 return checkVal;
}
bool cutPropagator::CheckJtPfNHFCutHi(cutPropagator inCutProp){return CheckJtPfNHFCutHi(inCutProp.GetJtPfNHFCutHi());}

bool cutPropagator::CheckJtPfNEFCutLow(std::vector<double> inJtPfNEFCutLow)
{
  bool checkVal = CheckVectDouble(jtPfNEFCutLow, inJtPfNEFCutLow);
  if(!checkVal) std::cout << "cutPropagator check failed on jtPfNEFCutLow" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckJtPfNEFCutLow(cutPropagator inCutProp){return CheckJtPfNEFCutLow(inCutProp.GetJtPfNEFCutLow());}

bool cutPropagator::CheckJtPfNEFCutHi(std::vector<double> inJtPfNEFCutHi)
{
  bool checkVal = CheckVectDouble(jtPfNEFCutHi, inJtPfNEFCutHi);
  if(!checkVal) std::cout << "cutPropagator check failed on jtPfNEFCutHi" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckJtPfNEFCutHi(cutPropagator inCutProp){return CheckJtPfNEFCutHi(inCutProp.GetJtPfNEFCutHi());}

bool cutPropagator::CheckJtPfMUFCutLow(std::vector<double> inJtPfMUFCutLow)
{
  bool checkVal = CheckVectDouble(jtPfMUFCutLow, inJtPfMUFCutLow);
  if(!checkVal) std::cout << "cutPropagator check failed on jtPfMUFCutLow" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckJtPfMUFCutLow(cutPropagator inCutProp){return CheckJtPfMUFCutLow(inCutProp.GetJtPfMUFCutLow());}

bool cutPropagator::CheckJtPfMUFCutHi(std::vector<double> inJtPfMUFCutHi)
{
  bool checkVal = CheckVectDouble(jtPfMUFCutHi, inJtPfMUFCutHi);
  if(!checkVal) std::cout << "cutPropagator check failed on jtPfMUFCutHi" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckJtPfMUFCutHi(cutPropagator inCutProp){return CheckJtPfMUFCutHi(inCutProp.GetJtPfMUFCutHi());}

bool cutPropagator::CheckJtPfCHFCutLow(std::vector<double> inJtPfCHFCutLow)
{
  bool checkVal = CheckVectDouble(jtPfCHFCutLow, inJtPfCHFCutLow);
  if(!checkVal) std::cout << "cutPropagator check failed on jtPfCHFCutLow" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckJtPfCHFCutLow(cutPropagator inCutProp){return CheckJtPfCHFCutLow(inCutProp.GetJtPfCHFCutLow());}

bool cutPropagator::CheckJtPfCHFCutHi(std::vector<double> inJtPfCHFCutHi)
{
  bool checkVal = CheckVectDouble(jtPfCHFCutHi, inJtPfCHFCutHi);
  if(!checkVal) std::cout << "cutPropagator check failed on jtPfCHFCutHi" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckJtPfCHFCutHi(cutPropagator inCutProp){return CheckJtPfCHFCutHi(inCutProp.GetJtPfCHFCutHi());}

bool cutPropagator::CheckJtPfCEFCutLow(std::vector<double> inJtPfCEFCutLow)
{
  bool checkVal = CheckVectDouble(jtPfCEFCutLow, inJtPfCEFCutLow);
  if(!checkVal) std::cout << "cutPropagator check failed on jtPfCEFCutLow" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckJtPfCEFCutLow(cutPropagator inCutProp){return CheckJtPfCEFCutLow(inCutProp.GetJtPfCEFCutLow());}

bool cutPropagator::CheckJtPfCEFCutHi(std::vector<double> inJtPfCEFCutHi)
{
  bool checkVal = CheckVectDouble(jtPfCEFCutHi, inJtPfCEFCutHi);
  if(!checkVal) std::cout << "cutPropagator check failed on jtPfCEFCutHi" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckJtPfCEFCutHi(cutPropagator inCutProp){return CheckJtPfCEFCutHi(inCutProp.GetJtPfCEFCutHi());}

bool cutPropagator::CheckJtPfMinMult(std::vector<int> inJtPfMinMult)
{
  bool checkVal = CheckVectInt(jtPfMinMult, inJtPfMinMult);
  if(!checkVal) std::cout << "cutPropagator check failed on jtPfMinMult" << std::endl;
 return checkVal;
}
bool cutPropagator::CheckJtPfMinMult(cutPropagator inCutProp){return CheckJtPfMinMult(inCutProp.GetJtPfMinMult());}

bool cutPropagator::CheckJtPfMinChgMult(std::vector<int> inJtPfMinChgMult)
{
  bool checkVal = CheckVectInt(jtPfMinChgMult, inJtPfMinChgMult);
  if(!checkVal) std::cout << "cutPropagator check failed on jtPfMinChgMult" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckJtPfMinChgMult(cutPropagator inCutProp){return CheckJtPfMinChgMult(inCutProp.GetJtPfMinChgMult());}

bool cutPropagator::CheckNSyst(int inNSyst)
{
  bool checkVal = CheckInt(inNSyst, nSyst);
  if(!checkVal) std::cout << "cutPropagator check failed on nSyst" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckNSyst(cutPropagator inCutProp){return CheckNSyst(inCutProp.GetNSyst());}

bool cutPropagator::CheckSystStr(std::vector<std::string> inSystStr)
{
  bool checkVal = CheckVectString(inSystStr, systStr);
  if(!checkVal) std::cout << "cutPropagator check failed on systStr" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckSystStr(cutPropagator inCutProp){return CheckSystStr(inCutProp.GetSystStr());}

bool cutPropagator::CheckNBayes(int inNBayes)
{
  bool checkVal = CheckInt(inNBayes, nBayes);
  if(!checkVal) std::cout << "cutPropagator check failed on nBayes" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckNBayes(cutPropagator inCutProp){return CheckNBayes(inCutProp.GetNBayes());}

bool cutPropagator::CheckNSVD(int inNSVD)
{
  bool checkVal = CheckInt(inNSVD, nSVD);
  if(!checkVal) std::cout << "cutPropagator check failed on nSVD" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckNSVD(cutPropagator inCutProp){return CheckNSVD(inCutProp.GetNSVD());}

bool cutPropagator::CheckNBigBayesSymm(int inNBigBayesSymm)
{
  bool checkVal = CheckInt(inNBigBayesSymm, nBigBayesSymm);
  if(!checkVal) std::cout << "cutPropagator check failed on nBigBayesSymm" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckNBigBayesSymm(cutPropagator inCutProp){return CheckNBigBayesSymm(inCutProp.GetNBigBayesSymm());}

bool cutPropagator::CheckBayesVal(std::vector<int> inBayesVal)
{
  bool checkVal = CheckVectInt(inBayesVal, bayesVal);
  if(!checkVal) std::cout << "cutPropagator check failed on bayesVal" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckBayesVal(cutPropagator inCutProp){return CheckBayesVal(inCutProp.GetBayesVal());}

bool cutPropagator::CheckNSuperBayes(int inNSuperBayes)
{
  bool checkVal = CheckInt(inNSuperBayes, nSuperBayes);
  if(!checkVal) std::cout << "cutPropagator check failed on nSuperBayes" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckNSuperBayes(cutPropagator inCutProp){return CheckNSuperBayes(inCutProp.GetNSuperBayes());}

bool cutPropagator::CheckRCDiffFileName(std::string inRCDiffFileName)
{
  bool checkVal = CheckString(inRCDiffFileName, rcDiffFileName);
  if(!checkVal) std::cout << "cutPropagator check failed on rcDiffFileName" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckRCDiffFileName(cutPropagator inCutProp){return CheckRCDiffFileName(inCutProp.GetRCDiffFileName());}

bool cutPropagator::CheckPriorFlatFileName(std::string inPriorFlatFileName)
{
  bool checkVal = CheckString(inPriorFlatFileName, flatPriorFileName);
  if(!checkVal) std::cout << "cutPropagator check failed on flatPriorFileName" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckPriorFlatFileName(cutPropagator inCutProp){return CheckPriorFlatFileName(inCutProp.GetPriorFlatFileName());}

bool cutPropagator::CheckJECVarMC(double inJECVarMC)
{
  bool checkVal = CheckDouble(inJECVarMC, jecVarMC);
  if(!checkVal) std::cout << "cutPropagator check failed on jecVarMC" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckJECVarMC(cutPropagator inCutProp){return CheckJECVarMC(inCutProp.GetJECVarMC());}

bool cutPropagator::CheckJERVarMC(double inJERVarMC)
{
  bool checkVal = CheckDouble(inJERVarMC, jerVarMC);
  if(!checkVal) std::cout << "cutPropagator check failed on jerVarMC" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckJERVarMC(cutPropagator inCutProp){return CheckJERVarMC(inCutProp.GetJERVarMC());}

bool cutPropagator::CheckJECVarData(double inJECVarData)
{
  bool checkVal = CheckDouble(inJECVarData, jecVarData);
  if(!checkVal) std::cout << "cutPropagator check failed on jecVarData" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckJECVarData(cutPropagator inCutProp){return CheckJECVarData(inCutProp.GetJECVarData());}

bool cutPropagator::CheckNResponseMod(int inNResponseMod)
{
  bool checkVal = CheckInt(inNResponseMod, nResponseMod);
  if(!checkVal) std::cout << "cutPropagator check failed on nResponseMod" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckNResponseMod(cutPropagator inCutProp){return CheckNResponseMod(inCutProp.GetNResponseMod());}

bool cutPropagator::CheckResponseMod(std::vector<double> inResponseMod)
{
  bool checkVal = CheckVectDouble(inResponseMod, responseMod);
  if(!checkVal) std::cout << "cutPropagator check failed on responseMod" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckResponseMod(cutPropagator inCutProp){return CheckResponseMod(inCutProp.GetResponseMod());}

bool cutPropagator::CheckJERVarData(std::vector<double> inJERVarData)
{
  bool checkVal = CheckVectDouble(inJERVarData, jerVarData);
  if(!checkVal) std::cout << "cutPropagator check failed on jerVarData" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckJERVarData(cutPropagator inCutProp){return CheckJERVarData(inCutProp.GetJERVarData());}

bool cutPropagator::CheckNPthats(int inNPthats)
{
  bool checkVal = CheckInt(inNPthats, nPthats);
  if(!checkVal) std::cout << "cutPropagator check failed on nPthats" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckNPthats(cutPropagator inCutProp){return CheckNPthats(inCutProp.GetNPthats());}

bool cutPropagator::CheckPthats(std::vector<double> inPthats)
{
  bool checkVal = CheckVectDouble(inPthats, pthats);
  if(!checkVal) std::cout << "cutPropagator check failed on pthats" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckPthats(cutPropagator inCutProp){return CheckPthats(inCutProp.GetPthats());}

bool cutPropagator::CheckPthatWeights(std::vector<double> inPthatWeights)
{
  bool checkVal = CheckVectDouble(inPthatWeights, pthatWeights);
  if(!checkVal) std::cout << "cutPropagator check failed on pthatWeights" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckPthatWeights(cutPropagator inCutProp){return CheckPthatWeights(inCutProp.GetPthatWeights());}

bool cutPropagator::CheckNJtAlgos(int inNJtAlgos)
{
  bool checkVal = CheckInt(inNJtAlgos, nJtAlgos);
  if(!checkVal) std::cout << "cutPropagator check failed on nJtAlgos" << std::endl;
 return checkVal;
}
bool cutPropagator::CheckNJtAlgos(cutPropagator inCutProp){return CheckNJtAlgos(inCutProp.GetNJtAlgos());}

bool cutPropagator::CheckJtAlgos(std::vector<std::string> inJtAlgos)
{
  bool checkVal = CheckVectString(inJtAlgos, jtAlgos);
  if(!checkVal) std::cout << "cutPropagator check failed on jtAlgos" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckJtAlgos(cutPropagator inCutProp){return CheckJtAlgos(inCutProp.GetJtAlgos());}

bool cutPropagator::CheckMinJtPtCut(std::vector<double> inMinJtPtCut)
{
  bool checkVal = CheckVectDouble(inMinJtPtCut, minJtPtCut);
  if(!checkVal) std::cout << "cutPropagator check failed on minJtPtCut" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckMinJtPtCut(cutPropagator inCutProp){return CheckMinJtPtCut(inCutProp.GetMinJtPtCut());}

bool cutPropagator::CheckMultiJtPtCut(std::vector<double> inMultiJtPtCut)
{
  bool checkVal = CheckVectDouble(inMultiJtPtCut, multiJtPtCut);
  if(!checkVal) std::cout << "cutPropagator check failed on multiJtPtCut" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckMultiJtPtCut(cutPropagator inCutProp){return CheckMultiJtPtCut(inCutProp.GetMultiJtPtCut());}

bool cutPropagator::CheckRecoTruncPos(std::vector<int> inRecoTruncPos)
{
  bool checkVal = CheckVectInt(inRecoTruncPos, recoTruncPos);
  if(!checkVal) std::cout << "cutPropagator check failed on recoTruncPos" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckRecoTruncPos(cutPropagator inCutProp){return CheckRecoTruncPos(inCutProp.GetRecoTruncPos());}

bool cutPropagator::CheckIsPP(bool inIsPP)
{
  bool checkVal = CheckBool(inIsPP, isPP);
  if(!checkVal) std::cout << "cutPropagator check failed on isPP" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckIsPP(cutPropagator inCutProp){return CheckIsPP(inCutProp.GetIsPP());}

bool cutPropagator::CheckNCentBins(int inNCentBins)
{
  bool checkVal = CheckInt(inNCentBins, nCentBins);
  if(!checkVal) std::cout << "cutPropagator check failed on nCentBins" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckNCentBins(cutPropagator inCutProp){return CheckNCentBins(inCutProp.GetNCentBins());}

bool cutPropagator::CheckCentBinsLow(std::vector<int> inCentBinsLow)
{
  bool checkVal = CheckVectInt(inCentBinsLow, centBinsLow);
  if(!checkVal) std::cout << "cutPropagator check failed on centBinsLow" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckCentBinsLow(cutPropagator inCutProp){return CheckCentBinsLow(inCutProp.GetCentBinsLow());}

bool cutPropagator::CheckCentBinsHi(std::vector<int> inCentBinsHi)
{
  bool checkVal = CheckVectInt(inCentBinsHi, centBinsHi);
  if(!checkVal) std::cout << "cutPropagator check failed on centBinsHi" << std::endl;
  return checkVal;
}
bool cutPropagator::CheckCentBinsHi(cutPropagator inCutProp){return CheckCentBinsHi(inCutProp.GetCentBinsHi());}
//end checks

int cutPropagator::GetGenNBinsFromRValCent(int rVal, std::string centStr)
{
  return GetGenPtBins(rVal, centStr).size()-1;
}

int cutPropagator::GetRecoNBinsFromRValCent(int rVal, std::string centStr)
{
  return GetRecoPtBins(rVal, centStr).size()-1;
}

void cutPropagator::GetGenBinsFromRValCent(int rVal, std::string centStr, double outBins[])
{
  std::vector<double> bins = GetGenPtBins(rVal,centStr);

  for(unsigned int bI = 0; bI < bins.size(); ++bI){
    outBins[bI] = bins[bI];
  }
  return;
}

void cutPropagator::GetRecoBinsFromRValCent(int rVal, std::string centStr, double outBins[])
{
  std::vector<double> bins = GetRecoPtBins(rVal,centStr);

  for(unsigned int bI = 0; bI < bins.size(); ++bI){
    outBins[bI] = bins[bI];
  }
  return;
}


std::vector<double> cutPropagator::GetGenPtBins(int rVal, std::string centStr)
{
  if(rVal == 1){
    if(centStr.find("Cent0to10") != std::string::npos) return GetGenPtBinsR1Cent0to10();
    else if(centStr.find("Cent10to30") != std::string::npos) return GetGenPtBinsR1Cent10to30();
    else if(centStr.find("Cent30to50") != std::string::npos) return GetGenPtBinsR1Cent30to50();
    else if(centStr.find("Cent50to90") != std::string::npos) return GetGenPtBinsR1Cent50to90();
  }
  else if(rVal == 2){
    if(centStr.find("Cent0to10") != std::string::npos) return GetGenPtBinsR2Cent0to10();
    else if(centStr.find("Cent10to30") != std::string::npos) return GetGenPtBinsR2Cent10to30();
    else if(centStr.find("Cent30to50") != std::string::npos) return GetGenPtBinsR2Cent30to50();
    else if(centStr.find("Cent50to90") != std::string::npos) return GetGenPtBinsR2Cent50to90();
  }
  else if(rVal == 3){
    if(centStr.find("Cent0to10") != std::string::npos) return GetGenPtBinsR3Cent0to10();
    else if(centStr.find("Cent10to30") != std::string::npos) return GetGenPtBinsR3Cent10to30();
    else if(centStr.find("Cent30to50") != std::string::npos) return GetGenPtBinsR3Cent30to50();
    else if(centStr.find("Cent50to90") != std::string::npos) return GetGenPtBinsR3Cent50to90();
  }
  else if(rVal == 4){
    if(centStr.find("Cent0to10") != std::string::npos) return GetGenPtBinsR4Cent0to10();
    else if(centStr.find("Cent10to30") != std::string::npos) return GetGenPtBinsR4Cent10to30();
    else if(centStr.find("Cent30to50") != std::string::npos) return GetGenPtBinsR4Cent30to50();
    else if(centStr.find("Cent50to90") != std::string::npos) return GetGenPtBinsR4Cent50to90();
  }
  else if(rVal == 5){
    if(centStr.find("Cent0to10") != std::string::npos) return GetGenPtBinsR5Cent0to10();
    else if(centStr.find("Cent10to30") != std::string::npos) return GetGenPtBinsR5Cent10to30();
    else if(centStr.find("Cent30to50") != std::string::npos) return GetGenPtBinsR5Cent30to50();
    else if(centStr.find("Cent50to90") != std::string::npos) return GetGenPtBinsR5Cent50to90();
  }
  else if(rVal == 6){
    if(centStr.find("Cent0to10") != std::string::npos) return GetGenPtBinsR6Cent0to10();
    else if(centStr.find("Cent10to30") != std::string::npos) return GetGenPtBinsR6Cent10to30();
    else if(centStr.find("Cent30to50") != std::string::npos) return GetGenPtBinsR6Cent30to50();
    else if(centStr.find("Cent50to90") != std::string::npos) return GetGenPtBinsR6Cent50to90();
  }
  else if(rVal == 8){
    if(centStr.find("Cent0to10") != std::string::npos) return GetGenPtBinsR8Cent0to10();
    else if(centStr.find("Cent10to30") != std::string::npos) return GetGenPtBinsR8Cent10to30();
    else if(centStr.find("Cent30to50") != std::string::npos) return GetGenPtBinsR8Cent30to50();
    else if(centStr.find("Cent50to90") != std::string::npos) return GetGenPtBinsR8Cent50to90();
  }
  else if(rVal == 10){
    if(centStr.find("Cent0to10") != std::string::npos) return GetGenPtBinsR10Cent0to10();
    else if(centStr.find("Cent10to30") != std::string::npos) return GetGenPtBinsR10Cent10to30();
    else if(centStr.find("Cent30to50") != std::string::npos) return GetGenPtBinsR10Cent30to50();
    else if(centStr.find("Cent50to90") != std::string::npos) return GetGenPtBinsR10Cent50to90();
  }

  return {};
}

std::vector<double> cutPropagator::GetRecoPtBins(int rVal, std::string centStr)
{
  if(rVal == 1){
    if(centStr.find("Cent0to10") != std::string::npos) return GetRecoPtBinsR1Cent0to10();
    else if(centStr.find("Cent10to30") != std::string::npos) return GetRecoPtBinsR1Cent10to30();
    else if(centStr.find("Cent30to50") != std::string::npos) return GetRecoPtBinsR1Cent30to50();
    else if(centStr.find("Cent50to90") != std::string::npos) return GetRecoPtBinsR1Cent50to90();
  }
  else if(rVal == 2){
    if(centStr.find("Cent0to10") != std::string::npos) return GetRecoPtBinsR2Cent0to10();
    else if(centStr.find("Cent10to30") != std::string::npos) return GetRecoPtBinsR2Cent10to30();
    else if(centStr.find("Cent30to50") != std::string::npos) return GetRecoPtBinsR2Cent30to50();
    else if(centStr.find("Cent50to90") != std::string::npos) return GetRecoPtBinsR2Cent50to90();
  }
  else if(rVal == 3){
    if(centStr.find("Cent0to10") != std::string::npos) return GetRecoPtBinsR3Cent0to10();
    else if(centStr.find("Cent10to30") != std::string::npos) return GetRecoPtBinsR3Cent10to30();
    else if(centStr.find("Cent30to50") != std::string::npos) return GetRecoPtBinsR3Cent30to50();
    else if(centStr.find("Cent50to90") != std::string::npos) return GetRecoPtBinsR3Cent50to90();
  }
  else if(rVal == 4){
    if(centStr.find("Cent0to10") != std::string::npos) return GetRecoPtBinsR4Cent0to10();
    else if(centStr.find("Cent10to30") != std::string::npos) return GetRecoPtBinsR4Cent10to30();
    else if(centStr.find("Cent30to50") != std::string::npos) return GetRecoPtBinsR4Cent30to50();
    else if(centStr.find("Cent50to90") != std::string::npos) return GetRecoPtBinsR4Cent50to90();
  }
  else if(rVal == 5){
    if(centStr.find("Cent0to10") != std::string::npos) return GetRecoPtBinsR5Cent0to10();
    else if(centStr.find("Cent10to30") != std::string::npos) return GetRecoPtBinsR5Cent10to30();
    else if(centStr.find("Cent30to50") != std::string::npos) return GetRecoPtBinsR5Cent30to50();
    else if(centStr.find("Cent50to90") != std::string::npos) return GetRecoPtBinsR5Cent50to90();
  }
  else if(rVal == 6){
    if(centStr.find("Cent0to10") != std::string::npos) return GetRecoPtBinsR6Cent0to10();
    else if(centStr.find("Cent10to30") != std::string::npos) return GetRecoPtBinsR6Cent10to30();
    else if(centStr.find("Cent30to50") != std::string::npos) return GetRecoPtBinsR6Cent30to50();
    else if(centStr.find("Cent50to90") != std::string::npos) return GetRecoPtBinsR6Cent50to90();
  }
  else if(rVal == 8){
    if(centStr.find("Cent0to10") != std::string::npos) return GetRecoPtBinsR8Cent0to10();
    else if(centStr.find("Cent10to30") != std::string::npos) return GetRecoPtBinsR8Cent10to30();
    else if(centStr.find("Cent30to50") != std::string::npos) return GetRecoPtBinsR8Cent30to50();
    else if(centStr.find("Cent50to90") != std::string::npos) return GetRecoPtBinsR8Cent50to90();
  }
  else if(rVal == 10){
    if(centStr.find("Cent0to10") != std::string::npos) return GetRecoPtBinsR10Cent0to10();
    else if(centStr.find("Cent10to30") != std::string::npos) return GetRecoPtBinsR10Cent10to30();
    else if(centStr.find("Cent30to50") != std::string::npos) return GetRecoPtBinsR10Cent30to50();
    else if(centStr.find("Cent50to90") != std::string::npos) return GetRecoPtBinsR10Cent50to90();
  }

  return {};
}


void cutPropagator::SetResponseMod(int inN, const Double_t inResponseMod[])
{
  for(int i = 0; i < inN; ++i){
    responseMod.push_back(inResponseMod[i]);
  }

  return;
}

void cutPropagator::SetJERVarData(int inN, const Double_t inJERVarData[])
{
  for(int i = 0; i < inN; ++i){
    jerVarData.push_back(inJERVarData[i]);
  }

  return;
}

void cutPropagator::SetJtAlgos(int inN, std::string inJtAlgos[])
{
  for(int i = 0; i < inN; ++i){
    jtAlgos.push_back(inJtAlgos[i]);
  }

  return;
}

void cutPropagator::SetMinJtPtCut(int inN, Double_t inMinJtPtCut[])
{
  for(int i = 0; i < inN; ++i){
    minJtPtCut.push_back(inMinJtPtCut[i]);
  }

  return;
}

void cutPropagator::SetMultiJtPtCut(int inN, Double_t inMultiJtPtCut[])
{
  for(int i = 0; i < inN; ++i){
    multiJtPtCut.push_back(inMultiJtPtCut[i]);
  }

  return;
}

void cutPropagator::SetRecoTruncPos(int inN, Int_t inRecoTruncPos[])
{
  for(int i = 0; i < inN; ++i){
    recoTruncPos.push_back(inRecoTruncPos[i]);
  }

  return;
}

void cutPropagator::SetRVals(int inN, const Int_t inRVals[])
{
  for(int i = 0; i < inN; ++i){
    rVals.push_back(inRVals[i]);
  }

  return;
}

void cutPropagator::SetGenPtBins(int rVal, std::string centStr, int inNGenPtBins, double inGenPtBins[])
{
  std::vector<double> bins;
  for(int i = 0; i < inNGenPtBins+1; ++i){
    bins.push_back(inGenPtBins[i]);
  }
  SetGenPtBins(rVal, centStr, bins);
  return;
}

void cutPropagator::SetRecoPtBins(int rVal, std::string centStr, int inNRecoPtBins, double inRecoPtBins[])
{
  std::vector<double> bins;
  for(int i = 0; i < inNRecoPtBins+1; ++i){
    bins.push_back(inRecoPtBins[i]);
  }
  SetRecoPtBins(rVal, centStr, bins);
  return;
}

void cutPropagator::SetGenPtBins(int rVal, std::string centStr, std::vector<double> inGenPtBins)
{
  if(rVal == 1){
    if(centStr.find("Cent0to10") != std::string::npos) SetGenPtBinsR1Cent0to10(inGenPtBins);
    else if(centStr.find("Cent10to30") != std::string::npos) SetGenPtBinsR1Cent10to30(inGenPtBins);
    else if(centStr.find("Cent30to50") != std::string::npos) SetGenPtBinsR1Cent30to50(inGenPtBins);
    else if(centStr.find("Cent50to90") != std::string::npos) SetGenPtBinsR1Cent50to90(inGenPtBins);
  }
  else if(rVal == 2){
    if(centStr.find("Cent0to10") != std::string::npos) SetGenPtBinsR2Cent0to10(inGenPtBins);
    else if(centStr.find("Cent10to30") != std::string::npos) SetGenPtBinsR2Cent10to30(inGenPtBins);
    else if(centStr.find("Cent30to50") != std::string::npos) SetGenPtBinsR2Cent30to50(inGenPtBins);
    else if(centStr.find("Cent50to90") != std::string::npos) SetGenPtBinsR2Cent50to90(inGenPtBins);
  }
  else if(rVal == 3){
    if(centStr.find("Cent0to10") != std::string::npos) SetGenPtBinsR3Cent0to10(inGenPtBins);
    else if(centStr.find("Cent10to30") != std::string::npos) SetGenPtBinsR3Cent10to30(inGenPtBins);
    else if(centStr.find("Cent30to50") != std::string::npos) SetGenPtBinsR3Cent30to50(inGenPtBins);
    else if(centStr.find("Cent50to90") != std::string::npos) SetGenPtBinsR3Cent50to90(inGenPtBins);
  }
  else if(rVal == 4){
    if(centStr.find("Cent0to10") != std::string::npos) SetGenPtBinsR4Cent0to10(inGenPtBins);
    else if(centStr.find("Cent10to30") != std::string::npos) SetGenPtBinsR4Cent10to30(inGenPtBins);
    else if(centStr.find("Cent30to50") != std::string::npos) SetGenPtBinsR4Cent30to50(inGenPtBins);
    else if(centStr.find("Cent50to90") != std::string::npos) SetGenPtBinsR4Cent50to90(inGenPtBins);
  }
  else if(rVal == 5){
    if(centStr.find("Cent0to10") != std::string::npos) SetGenPtBinsR5Cent0to10(inGenPtBins);
    else if(centStr.find("Cent10to30") != std::string::npos) SetGenPtBinsR5Cent10to30(inGenPtBins);
    else if(centStr.find("Cent30to50") != std::string::npos) SetGenPtBinsR5Cent30to50(inGenPtBins);
    else if(centStr.find("Cent50to90") != std::string::npos) SetGenPtBinsR5Cent50to90(inGenPtBins);
  }
  else if(rVal == 6){
    if(centStr.find("Cent0to10") != std::string::npos) SetGenPtBinsR6Cent0to10(inGenPtBins);
    else if(centStr.find("Cent10to30") != std::string::npos) SetGenPtBinsR6Cent10to30(inGenPtBins);
    else if(centStr.find("Cent30to50") != std::string::npos) SetGenPtBinsR6Cent30to50(inGenPtBins);
    else if(centStr.find("Cent50to90") != std::string::npos) SetGenPtBinsR6Cent50to90(inGenPtBins);
  }
  else if(rVal == 8){
    if(centStr.find("Cent0to10") != std::string::npos) SetGenPtBinsR8Cent0to10(inGenPtBins);
    else if(centStr.find("Cent10to30") != std::string::npos) SetGenPtBinsR8Cent10to30(inGenPtBins);
    else if(centStr.find("Cent30to50") != std::string::npos) SetGenPtBinsR8Cent30to50(inGenPtBins);
    else if(centStr.find("Cent50to90") != std::string::npos) SetGenPtBinsR8Cent50to90(inGenPtBins);
  }
  else if(rVal == 10){
    if(centStr.find("Cent0to10") != std::string::npos) SetGenPtBinsR10Cent0to10(inGenPtBins);
    else if(centStr.find("Cent10to30") != std::string::npos) SetGenPtBinsR10Cent10to30(inGenPtBins);
    else if(centStr.find("Cent30to50") != std::string::npos) SetGenPtBinsR10Cent30to50(inGenPtBins);
    else if(centStr.find("Cent50to90") != std::string::npos) SetGenPtBinsR10Cent50to90(inGenPtBins);
  }

  return;
}


void cutPropagator::SetRecoPtBins(int rVal, std::string centStr, std::vector<double> inRecoPtBins)
{
  if(rVal == 1){
    if(centStr.find("Cent0to10") != std::string::npos) SetRecoPtBinsR1Cent0to10(inRecoPtBins);
    else if(centStr.find("Cent10to30") != std::string::npos) SetRecoPtBinsR1Cent10to30(inRecoPtBins);
    else if(centStr.find("Cent30to50") != std::string::npos) SetRecoPtBinsR1Cent30to50(inRecoPtBins);
    else if(centStr.find("Cent50to90") != std::string::npos) SetRecoPtBinsR1Cent50to90(inRecoPtBins);
  }
  else if(rVal == 2){
    if(centStr.find("Cent0to10") != std::string::npos) SetRecoPtBinsR2Cent0to10(inRecoPtBins);
    else if(centStr.find("Cent10to30") != std::string::npos) SetRecoPtBinsR2Cent10to30(inRecoPtBins);
    else if(centStr.find("Cent30to50") != std::string::npos) SetRecoPtBinsR2Cent30to50(inRecoPtBins);
    else if(centStr.find("Cent50to90") != std::string::npos) SetRecoPtBinsR2Cent50to90(inRecoPtBins);
  }
  else if(rVal == 3){
    if(centStr.find("Cent0to10") != std::string::npos) SetRecoPtBinsR3Cent0to10(inRecoPtBins);
    else if(centStr.find("Cent10to30") != std::string::npos) SetRecoPtBinsR3Cent10to30(inRecoPtBins);
    else if(centStr.find("Cent30to50") != std::string::npos) SetRecoPtBinsR3Cent30to50(inRecoPtBins);
    else if(centStr.find("Cent50to90") != std::string::npos) SetRecoPtBinsR3Cent50to90(inRecoPtBins);
  }
  else if(rVal == 4){
    if(centStr.find("Cent0to10") != std::string::npos) SetRecoPtBinsR4Cent0to10(inRecoPtBins);
    else if(centStr.find("Cent10to30") != std::string::npos) SetRecoPtBinsR4Cent10to30(inRecoPtBins);
    else if(centStr.find("Cent30to50") != std::string::npos) SetRecoPtBinsR4Cent30to50(inRecoPtBins);
    else if(centStr.find("Cent50to90") != std::string::npos) SetRecoPtBinsR4Cent50to90(inRecoPtBins);
  }
  else if(rVal == 5){
    if(centStr.find("Cent0to10") != std::string::npos) SetRecoPtBinsR5Cent0to10(inRecoPtBins);
    else if(centStr.find("Cent10to30") != std::string::npos) SetRecoPtBinsR5Cent10to30(inRecoPtBins);
    else if(centStr.find("Cent30to50") != std::string::npos) SetRecoPtBinsR5Cent30to50(inRecoPtBins);
    else if(centStr.find("Cent50to90") != std::string::npos) SetRecoPtBinsR5Cent50to90(inRecoPtBins);
  }
  else if(rVal == 6){
    if(centStr.find("Cent0to10") != std::string::npos) SetRecoPtBinsR6Cent0to10(inRecoPtBins);
    else if(centStr.find("Cent10to30") != std::string::npos) SetRecoPtBinsR6Cent10to30(inRecoPtBins);
    else if(centStr.find("Cent30to50") != std::string::npos) SetRecoPtBinsR6Cent30to50(inRecoPtBins);
    else if(centStr.find("Cent50to90") != std::string::npos) SetRecoPtBinsR6Cent50to90(inRecoPtBins);
  }
  else if(rVal == 8){
    if(centStr.find("Cent0to10") != std::string::npos) SetRecoPtBinsR8Cent0to10(inRecoPtBins);
    else if(centStr.find("Cent10to30") != std::string::npos) SetRecoPtBinsR8Cent10to30(inRecoPtBins);
    else if(centStr.find("Cent30to50") != std::string::npos) SetRecoPtBinsR8Cent30to50(inRecoPtBins);
    else if(centStr.find("Cent50to90") != std::string::npos) SetRecoPtBinsR8Cent50to90(inRecoPtBins);
  }
  else if(rVal == 10){
    if(centStr.find("Cent0to10") != std::string::npos) SetRecoPtBinsR10Cent0to10(inRecoPtBins);
    else if(centStr.find("Cent10to30") != std::string::npos) SetRecoPtBinsR10Cent10to30(inRecoPtBins);
    else if(centStr.find("Cent30to50") != std::string::npos) SetRecoPtBinsR10Cent30to50(inRecoPtBins);
    else if(centStr.find("Cent50to90") != std::string::npos) SetRecoPtBinsR10Cent50to90(inRecoPtBins);
  }

  return;
}


void cutPropagator::SetGeneralBins(int inN, const Double_t inGeneralBins[])
{
  for(int i = 0; i < inN; ++i){
    generalBins.push_back(inGeneralBins[i]);
  }

  return;
}

void cutPropagator::SetJtAbsEtaBinsLow(int inN, const Double_t inJtAbsEtaBinsLow[])
{
  for(int i = 0; i < inN; ++i){
    jtAbsEtaBinsLow.push_back(inJtAbsEtaBinsLow[i]);
  }

  return;
}


void cutPropagator::SetJtAbsEtaBinsHi(int inN, const Double_t inJtAbsEtaBinsHi[])
{
  for(int i = 0; i < inN; ++i){
    jtAbsEtaBinsHi.push_back(inJtAbsEtaBinsHi[i]);
  }

  return;
}


void cutPropagator::SetPthats(int inN, const Double_t inPthats[])
{
  for(int i = 0; i < inN; ++i){
    pthats.push_back(inPthats[i]);
  }

  return;
}


void cutPropagator::SetPthatWeights(int inN, const Double_t inPthatWeights[])
{
  for(int i = 0; i < inN; ++i){
    pthatWeights.push_back(inPthatWeights[i]);
  }

  return;
}


void cutPropagator::SetCentBinsLow(int inN, const Int_t inCentBinsLow[])
{
  for(int i = 0; i < inN; ++i){
    centBinsLow.push_back(inCentBinsLow[i]);
  }

  return;
}


void cutPropagator::SetCentBinsHi(int inN, const Int_t inCentBinsHi[])
{
  for(int i = 0; i < inN; ++i){
    centBinsHi.push_back(inCentBinsHi[i]);
  }

  return;
}


void cutPropagator::SetIdStr(int inN, const std::string inIdStr[])
{
  for(int i = 0; i < inN; ++i){
    idStr.push_back(inIdStr[i]);
  }

  return;
}

void cutPropagator::SetJtPfCHMFCutLow(int inN, const Double_t inJtPfCHMFCutLow[])
{
  for(int i = 0; i < inN; ++i){
    jtPfCHMFCutLow.push_back(inJtPfCHMFCutLow[i]);
  }

  return;
}

void cutPropagator::SetJtPfCHMFCutHi(int inN, const Double_t inJtPfCHMFCutHi[])
{
  for(int i = 0; i < inN; ++i){
    jtPfCHMFCutHi.push_back(inJtPfCHMFCutHi[i]);
  }

  return;
}


void cutPropagator::SetJtPfMUMFCutLow(int inN, const Double_t inJtPfMUMFCutLow[])
{
  for(int i = 0; i < inN; ++i){
    jtPfMUMFCutLow.push_back(inJtPfMUMFCutLow[i]);
  }

  return;
}

void cutPropagator::SetJtPfMUMFCutHi(int inN, const Double_t inJtPfMUMFCutHi[])
{
  for(int i = 0; i < inN; ++i){
    jtPfMUMFCutHi.push_back(inJtPfMUMFCutHi[i]);
  }

  return;
}


void cutPropagator::SetJtPfNHFCutLow(int inN, const Double_t inJtPfNHFCutLow[])
{
  for(int i = 0; i < inN; ++i){
    jtPfNHFCutLow.push_back(inJtPfNHFCutLow[i]);
  }

  return;
}

void cutPropagator::SetJtPfNHFCutHi(int inN, const Double_t inJtPfNHFCutHi[])
{
  for(int i = 0; i < inN; ++i){
    jtPfNHFCutHi.push_back(inJtPfNHFCutHi[i]);
  }

  return;
}

void cutPropagator::SetJtPfNEFCutLow(int inN, const Double_t inJtPfNEFCutLow[])
{
  for(int i = 0; i < inN; ++i){
    jtPfNEFCutLow.push_back(inJtPfNEFCutLow[i]);
  }

  return;
}

void cutPropagator::SetJtPfNEFCutHi(int inN, const Double_t inJtPfNEFCutHi[])
{
  for(int i = 0; i < inN; ++i){
    jtPfNEFCutHi.push_back(inJtPfNEFCutHi[i]);
  }

  return;
}

void cutPropagator::SetJtPfMUFCutLow(int inN, const Double_t inJtPfMUFCutLow[])
{
  for(int i = 0; i < inN; ++i){
    jtPfMUFCutLow.push_back(inJtPfMUFCutLow[i]);
  }

  return;
}

void cutPropagator::SetJtPfMUFCutHi(int inN, const Double_t inJtPfMUFCutHi[])
{
  for(int i = 0; i < inN; ++i){
    jtPfMUFCutHi.push_back(inJtPfMUFCutHi[i]);
  }

  return;
}

void cutPropagator::SetJtPfCHFCutLow(int inN, const Double_t inJtPfCHFCutLow[])
{
  for(int i = 0; i < inN; ++i){
    jtPfCHFCutLow.push_back(inJtPfCHFCutLow[i]);
  }

  return;
}

void cutPropagator::SetJtPfCHFCutHi(int inN, const Double_t inJtPfCHFCutHi[])
{
  for(int i = 0; i < inN; ++i){
    jtPfCHFCutHi.push_back(inJtPfCHFCutHi[i]);
  }

  return;
}

void cutPropagator::SetJtPfCEFCutLow(int inN, const Double_t inJtPfCEFCutLow[])
{
  for(int i = 0; i < inN; ++i){
    jtPfCEFCutLow.push_back(inJtPfCEFCutLow[i]);
  }

  return;
}

void cutPropagator::SetJtPfCEFCutHi(int inN, const Double_t inJtPfCEFCutHi[])
{
  for(int i = 0; i < inN; ++i){
    jtPfCEFCutHi.push_back(inJtPfCEFCutHi[i]);
  }

  return;
}

void cutPropagator::SetJtPfMinMult(int inN, const Int_t inJtPfMinMult[])
{
  for(int i = 0; i < inN; ++i){
    jtPfMinMult.push_back(inJtPfMinMult[i]);
  }

  return;
}

void cutPropagator::SetJtPfMinChgMult(int inN, const Int_t inJtPfMinChgMult[])
{
  for(int i = 0; i < inN; ++i){
    jtPfMinChgMult.push_back(inJtPfMinChgMult[i]);
  }

  return;
}

void cutPropagator::SetSystStr(int inN, const std::string inSystStr[])
{
  for(int i = 0; i < inN; ++i){
    systStr.push_back(inSystStr[i]);
  }

  return;
}

void cutPropagator::SetBayesVal(int inN, const int inBayesVal[])
{
  for(int i = 0; i < inN; ++i){
    bayesVal.push_back(inBayesVal[i]);
  }

  return;
}


void cutPropagator::SetHistTagBayes(int inN, const std::string inHistTagBayes[])
{
  for(int i = 0; i < inN; ++i){
    histTagBayes.push_back(inHistTagBayes[i]);
  }

  return;
}

void cutPropagator::SetHistTagSVD(int inN, const std::string inHistTagSVD[])
{
  for(int i = 0; i < inN; ++i){
    histTagSVD.push_back(inHistTagSVD[i]);
  }

  return;
}


void cutPropagator::SetHistBestBayes(int inN, const int inHistBestBayes[])
{
  for(int i = 0; i < inN; ++i){
    histBestBayes.push_back(inHistBestBayes[i]);
  }

  return;
}

void cutPropagator::SetHistBestSVD(int inN, const int inHistBestSVD[])
{
  for(int i = 0; i < inN; ++i){
    histBestSVD.push_back(inHistBestSVD[i]);
  }

  return;
}

std::vector<int> cutPropagator::StringToIntVect(std::string inStr)
{
  std::vector<int> intVect;
  while(inStr.find(",") != std::string::npos){
    intVect.push_back(std::stoi(inStr.substr(0, inStr.find(","))));
    inStr.replace(0, inStr.find(",")+1, "");
  }
  if(inStr.size() != 0) intVect.push_back(std::stoi(inStr));
  return intVect;
}

std::vector<double> cutPropagator::StringToDoubleVect(std::string inStr)
{
  std::vector<double> doubleVect;
  while(inStr.find(",") != std::string::npos){
    doubleVect.push_back(std::stod(inStr.substr(0, inStr.find(","))));
    inStr.replace(0, inStr.find(",")+1, "");
  }
  if(inStr.size() != 0) doubleVect.push_back(std::stod(inStr));
  return doubleVect;
}

std::vector<std::string> cutPropagator::StringToStringVect(std::string inStr)
{
  std::vector<std::string> stringVect;
  while(inStr.find(",") != std::string::npos){
    stringVect.push_back(inStr.substr(0, inStr.find(",")));
    inStr.replace(0, inStr.find(",")+1, "");
  }
  if(inStr.size() != 0) stringVect.push_back(inStr);
  return stringVect;
}


std::string cutPropagator::to_string_with_precision(double a_value, const int n)
{
  std::ostringstream out;
  out << std::setprecision(n) << a_value;
  return out.str();
}

#endif
