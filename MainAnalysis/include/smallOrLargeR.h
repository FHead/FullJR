#ifndef SMALLORLARGER_H
#define SMALLORLARGER_H

#include <iostream>
#include <vector>

class smallOrLargeR
{
 public:
  const double deltaBin = 0.1;
  //Chosen loose since bins will be effective integer but double because of TH1 initialization
  static const int nR = 6;
  std::vector<int> rVals = {2, 3, 4, 6, 8, 10};
   
  std::vector<double> genPtBinsR2Cent0to10 = {20, 100, 150, 200, 300, 400, 500, 630, 1000, 1500};
  std::vector<double> recoPtBinsR2Cent0to10 = {200, 300, 400, 500, 630, 1000};
  std::vector<double> genPtBinsR2Cent10to30 = {20, 100, 150, 200, 300, 400, 500, 630, 1000, 1500};
  std::vector<double> recoPtBinsR2Cent10to30 = {200, 300, 400, 500, 630, 1000};
  std::vector<double> genPtBinsR2Cent30to50 = {20, 100, 150, 200, 300, 400, 500, 630, 1000};
  std::vector<double> recoPtBinsR2Cent30to50 = {200, 300, 400, 500, 630};
  std::vector<double> genPtBinsR2Cent50to90 = {20, 100, 150, 200, 300, 400, 500, 630};
  std::vector<double> recoPtBinsR2Cent50to90 = {200, 300, 400, 500};

  std::vector<double> genPtBinsR3Cent0to10 = {20, 100, 150, 200, 300, 400, 500, 630, 1000, 1500};
  std::vector<double> recoPtBinsR3Cent0to10 = {200, 300, 400, 500, 630, 1000};
  std::vector<double> genPtBinsR3Cent10to30 = {20, 100, 150, 200, 300, 400, 500, 630, 1000, 1500};
  std::vector<double> recoPtBinsR3Cent10to30 = {200, 300, 400, 500, 630, 1000};
  std::vector<double> genPtBinsR3Cent30to50 = {20, 100, 150, 200, 300, 400, 500, 630, 1000};
  std::vector<double> recoPtBinsR3Cent30to50 = {200, 300, 400, 500, 630};
  std::vector<double> genPtBinsR3Cent50to90 = {20, 100, 150, 200, 300, 400, 500, 630};
  std::vector<double> recoPtBinsR3Cent50to90 = {200, 300, 400, 500};

  std::vector<double> genPtBinsR4Cent0to10 = {20, 100, 150, 200, 300, 400, 500, 630, 1000, 1500};
  std::vector<double> recoPtBinsR4Cent0to10 = {200, 300, 400, 500, 630, 1000};
  std::vector<double> genPtBinsR4Cent10to30 = {20, 100, 150, 200, 300, 400, 500, 630, 1000, 1500};
  std::vector<double> recoPtBinsR4Cent10to30 = {200, 300, 400, 500, 630, 1000};
  std::vector<double> genPtBinsR4Cent30to50 = {20, 100, 150, 200, 300, 400, 500, 630, 1000};
  std::vector<double> recoPtBinsR4Cent30to50 = {200, 300, 400, 500, 630};
  std::vector<double> genPtBinsR4Cent50to90 = {20, 100, 150, 200, 300, 400, 500, 630};
  std::vector<double> recoPtBinsR4Cent50to90 = {200, 300, 400, 500};

  std::vector<double> genPtBinsR6Cent0to10 = {20, 100, 150, 200, 300, 400, 500, 630, 1000, 1500};
  std::vector<double> recoPtBinsR6Cent0to10 = {200, 300, 400, 500, 630, 1000};
  std::vector<double> genPtBinsR6Cent10to30 = {20, 100, 150, 200, 300, 400, 500, 630, 1000, 1500};
  std::vector<double> recoPtBinsR6Cent10to30 = {200, 300, 400, 500, 630, 1000};
  std::vector<double> genPtBinsR6Cent30to50 = {20, 100, 150, 200, 300, 400, 500, 630, 1000};
  std::vector<double> recoPtBinsR6Cent30to50 = {200, 300, 400, 500, 630};
  std::vector<double> genPtBinsR6Cent50to90 = {20, 100, 150, 200, 300, 400, 500, 630};
  std::vector<double> recoPtBinsR6Cent50to90 = {200, 300, 400, 500};

  std::vector<double> genPtBinsR8Cent0to10 = {20, 100, 150, 200, 300, 400, 500, 630, 1000, 1500};
  std::vector<double> recoPtBinsR8Cent0to10 = {300, 400, 500, 630, 1000};
  std::vector<double> genPtBinsR8Cent10to30 = {20, 100, 150, 200, 300, 400, 500, 630, 1000, 1500};
  std::vector<double> recoPtBinsR8Cent10to30 = {200, 300, 400, 500, 630, 1000};
  std::vector<double> genPtBinsR8Cent30to50 = {20, 100, 150, 200, 300, 400, 500, 630, 1000};
  std::vector<double> recoPtBinsR8Cent30to50 = {200, 300, 400, 500, 630};
  std::vector<double> genPtBinsR8Cent50to90 = {20, 100, 150, 200, 300, 400, 500, 630};
  std::vector<double> recoPtBinsR8Cent50to90 = {200, 300, 400, 500};

  std::vector<double> genPtBinsR10Cent0to10 = {20, 100, 150, 200, 300, 400, 500, 630, 1000, 1500};
  std::vector<double> recoPtBinsR10Cent0to10 = {300, 400, 500, 630, 1000};
  std::vector<double> genPtBinsR10Cent10to30 = {20, 100, 150, 200, 300, 400, 500, 630, 1000, 1500};
  std::vector<double> recoPtBinsR10Cent10to30 = {300, 400, 500, 630, 1000};
  std::vector<double> genPtBinsR10Cent30to50 = {20, 100, 150, 200, 300, 400, 500, 630, 1000};
  std::vector<double> recoPtBinsR10Cent30to50 = {200, 300, 400, 500, 630};
  std::vector<double> genPtBinsR10Cent50to90 = {20, 100, 150, 200, 300, 400, 500, 630};
  std::vector<double> recoPtBinsR10Cent50to90 = {200, 300, 400, 500};


  static const int generalBinsLow = 20;
  static const int generalBinsHigh = 1500;
  static const int generalBinsInterval = 10;

  static const int nGeneralBins = (generalBinsHigh - generalBinsLow)/generalBinsInterval;
  std::vector<double> generalBins;

  smallOrLargeR();

  int GetNR(){return nR;}
  std::vector<int> GetRVals(){return rVals;}

  std::vector<double> GetGenPtBinsR2Cent0to10(){return genPtBinsR2Cent0to10;}
  std::vector<double> GetRecoPtBinsR2Cent0to10(){return recoPtBinsR2Cent0to10;}
  std::vector<double> GetGenPtBinsR2Cent10to30(){return genPtBinsR2Cent10to30;}
  std::vector<double> GetRecoPtBinsR2Cent10to30(){return recoPtBinsR2Cent10to30;}
  std::vector<double> GetGenPtBinsR2Cent30to50(){return genPtBinsR2Cent30to50;}
  std::vector<double> GetRecoPtBinsR2Cent30to50(){return recoPtBinsR2Cent30to50;}
  std::vector<double> GetGenPtBinsR2Cent50to90(){return genPtBinsR2Cent50to90;}
  std::vector<double> GetRecoPtBinsR2Cent50to90(){return recoPtBinsR2Cent50to90;}

  std::vector<double> GetGenPtBinsR3Cent0to10(){return genPtBinsR3Cent0to10;}
  std::vector<double> GetRecoPtBinsR3Cent0to10(){return recoPtBinsR3Cent0to10;}
  std::vector<double> GetGenPtBinsR3Cent10to30(){return genPtBinsR3Cent10to30;}
  std::vector<double> GetRecoPtBinsR3Cent10to30(){return recoPtBinsR3Cent10to30;}
  std::vector<double> GetGenPtBinsR3Cent30to50(){return genPtBinsR3Cent30to50;}
  std::vector<double> GetRecoPtBinsR3Cent30to50(){return recoPtBinsR3Cent30to50;}
  std::vector<double> GetGenPtBinsR3Cent50to90(){return genPtBinsR3Cent50to90;}
  std::vector<double> GetRecoPtBinsR3Cent50to90(){return recoPtBinsR3Cent50to90;}

  std::vector<double> GetGenPtBinsR4Cent0to10(){return genPtBinsR4Cent0to10;}
  std::vector<double> GetRecoPtBinsR4Cent0to10(){return recoPtBinsR4Cent0to10;}
  std::vector<double> GetGenPtBinsR4Cent10to30(){return genPtBinsR4Cent10to30;}
  std::vector<double> GetRecoPtBinsR4Cent10to30(){return recoPtBinsR4Cent10to30;}
  std::vector<double> GetGenPtBinsR4Cent30to50(){return genPtBinsR4Cent30to50;}
  std::vector<double> GetRecoPtBinsR4Cent30to50(){return recoPtBinsR4Cent30to50;}
  std::vector<double> GetGenPtBinsR4Cent50to90(){return genPtBinsR4Cent50to90;}
  std::vector<double> GetRecoPtBinsR4Cent50to90(){return recoPtBinsR4Cent50to90;}

  std::vector<double> GetGenPtBinsR6Cent0to10(){return genPtBinsR6Cent0to10;}
  std::vector<double> GetRecoPtBinsR6Cent0to10(){return recoPtBinsR6Cent0to10;}
  std::vector<double> GetGenPtBinsR6Cent10to30(){return genPtBinsR6Cent10to30;}
  std::vector<double> GetRecoPtBinsR6Cent10to30(){return recoPtBinsR6Cent10to30;}
  std::vector<double> GetGenPtBinsR6Cent30to50(){return genPtBinsR6Cent30to50;}
  std::vector<double> GetRecoPtBinsR6Cent30to50(){return recoPtBinsR6Cent30to50;}
  std::vector<double> GetGenPtBinsR6Cent50to90(){return genPtBinsR6Cent50to90;}
  std::vector<double> GetRecoPtBinsR6Cent50to90(){return recoPtBinsR6Cent50to90;}

  std::vector<double> GetGenPtBinsR8Cent0to10(){return genPtBinsR8Cent0to10;}
  std::vector<double> GetRecoPtBinsR8Cent0to10(){return recoPtBinsR8Cent0to10;}
  std::vector<double> GetGenPtBinsR8Cent10to30(){return genPtBinsR8Cent10to30;}
  std::vector<double> GetRecoPtBinsR8Cent10to30(){return recoPtBinsR8Cent10to30;}
  std::vector<double> GetGenPtBinsR8Cent30to50(){return genPtBinsR8Cent30to50;}
  std::vector<double> GetRecoPtBinsR8Cent30to50(){return recoPtBinsR8Cent30to50;}
  std::vector<double> GetGenPtBinsR8Cent50to90(){return genPtBinsR8Cent50to90;}
  std::vector<double> GetRecoPtBinsR8Cent50to90(){return recoPtBinsR8Cent50to90;}

  std::vector<double> GetGenPtBinsR10Cent0to10(){return genPtBinsR10Cent0to10;}
  std::vector<double> GetRecoPtBinsR10Cent0to10(){return recoPtBinsR10Cent0to10;}
  std::vector<double> GetGenPtBinsR10Cent10to30(){return genPtBinsR10Cent10to30;}
  std::vector<double> GetRecoPtBinsR10Cent10to30(){return recoPtBinsR10Cent10to30;}
  std::vector<double> GetGenPtBinsR10Cent30to50(){return genPtBinsR10Cent30to50;}
  std::vector<double> GetRecoPtBinsR10Cent30to50(){return recoPtBinsR10Cent30to50;}
  std::vector<double> GetGenPtBinsR10Cent50to90(){return genPtBinsR10Cent50to90;}
  std::vector<double> GetRecoPtBinsR10Cent50to90(){return recoPtBinsR10Cent50to90;}

  std::vector<double> GetGenPtBins(int rVal, std::string centStr);
  std::vector<double> GetRecoPtBins(int rVal, std::string centStr);

  int GetGenNBinsFromRValCent(int rVal, std::string centStr);
  void GetGenBinsFromRValCent(int rVal, std::string centStr, double outBins[]);
  int GetRecoNBinsFromRValCent(int rVal, std::string centStr);
  void GetRecoBinsFromRValCent(int rVal, std::string centStr, double outBins[]);

  int GetNGeneralBins(){return nGeneralBins;}
  std::vector<double> GetGeneralBins(){return generalBins;}
  void GetGeneralBins(Int_t nBins, Double_t bins[]);

  bool CheckIntVect(std::vector<int> vect1, std::vector<int> vect2);
  bool CheckNBins(int inNBins, int compNBins, std::string nBinsStr);
  void SizeIsLTOrEQ(int inNBins, std::vector<double> inBins, std::string nBinsStr, std::string binsStr);  

  bool CheckDoubleVect(std::vector<double> vect1, std::vector<double> vect2);
};

smallOrLargeR::smallOrLargeR()
{ 
  generalBins.reserve(nGeneralBins+1);
  
  for(Int_t bI = 0; bI < nGeneralBins+1; ++bI){
    generalBins.push_back(generalBinsLow + bI*generalBinsInterval);
  }
  
  return;
}

int smallOrLargeR::GetGenNBinsFromRValCent(int rVal, std::string centStr)
{
  return GetGenPtBins(rVal, centStr).size()-1;
}

int smallOrLargeR::GetRecoNBinsFromRValCent(int rVal, std::string centStr)
{
  return GetRecoPtBins(rVal, centStr).size()-1;
}

void smallOrLargeR::GetGenBinsFromRValCent(int rVal, std::string centStr, double outBins[])
{
  std::vector<double> bins = GetGenPtBins(rVal,centStr);

  for(unsigned int bI = 0; bI < bins.size(); ++bI){
    outBins[bI] = bins[bI];
  }
  return;
}

void smallOrLargeR::GetRecoBinsFromRValCent(int rVal, std::string centStr, double outBins[])
{
  std::vector<double> bins = GetRecoPtBins(rVal,centStr);

  for(unsigned int bI = 0; bI < bins.size(); ++bI){
    outBins[bI] = bins[bI];
  }
  return;
}
 
std::vector<double> smallOrLargeR::GetGenPtBins(int rVal, std::string centStr)
{
  if(rVal == 2){
    if(centStr.find("Cent0to10") != std::string::npos) return genPtBinsR2Cent0to10;
    else if(centStr.find("Cent10to30") != std::string::npos) return genPtBinsR2Cent10to30;
    else if(centStr.find("Cent30to50") != std::string::npos) return genPtBinsR2Cent30to50;
    else if(centStr.find("Cent50to90") != std::string::npos) return genPtBinsR2Cent50to90;
  }
  else if(rVal == 3){
    if(centStr.find("Cent0to10") != std::string::npos) return genPtBinsR3Cent0to10;
    else if(centStr.find("Cent10to30") != std::string::npos) return genPtBinsR3Cent10to30;
    else if(centStr.find("Cent30to50") != std::string::npos) return genPtBinsR3Cent30to50;
    else if(centStr.find("Cent50to90") != std::string::npos) return genPtBinsR3Cent50to90;
  }
  else if(rVal == 4){
    if(centStr.find("Cent0to10") != std::string::npos) return genPtBinsR4Cent0to10;
    else if(centStr.find("Cent10to30") != std::string::npos) return genPtBinsR4Cent10to30;
    else if(centStr.find("Cent30to50") != std::string::npos) return genPtBinsR4Cent30to50;
    else if(centStr.find("Cent50to90") != std::string::npos) return genPtBinsR4Cent50to90;
  }
  else if(rVal == 6){
    if(centStr.find("Cent0to10") != std::string::npos) return genPtBinsR6Cent0to10;
    else if(centStr.find("Cent10to30") != std::string::npos) return genPtBinsR6Cent10to30;
    else if(centStr.find("Cent30to50") != std::string::npos) return genPtBinsR6Cent30to50;
    else if(centStr.find("Cent50to90") != std::string::npos) return genPtBinsR6Cent50to90;
  }
  else if(rVal == 8){
    if(centStr.find("Cent0to10") != std::string::npos) return genPtBinsR8Cent0to10;
    else if(centStr.find("Cent10to30") != std::string::npos) return genPtBinsR8Cent10to30;
    else if(centStr.find("Cent30to50") != std::string::npos) return genPtBinsR8Cent30to50;
    else if(centStr.find("Cent50to90") != std::string::npos) return genPtBinsR8Cent50to90;
  }
  else if(rVal == 10){
    if(centStr.find("Cent0to10") != std::string::npos) return genPtBinsR10Cent0to10;
    else if(centStr.find("Cent10to30") != std::string::npos) return genPtBinsR10Cent10to30;
    else if(centStr.find("Cent30to50") != std::string::npos) return genPtBinsR10Cent30to50;
    else if(centStr.find("Cent50to90") != std::string::npos) return genPtBinsR10Cent50to90;
  }

  return {};
}

std::vector<double> smallOrLargeR::GetRecoPtBins(int rVal, std::string centStr)
{
  if(rVal == 2){
    if(centStr.find("Cent0to10") != std::string::npos) return recoPtBinsR2Cent0to10;
    else if(centStr.find("Cent10to30") != std::string::npos) return recoPtBinsR2Cent10to30;
    else if(centStr.find("Cent30to50") != std::string::npos) return recoPtBinsR2Cent30to50;
    else if(centStr.find("Cent50to90") != std::string::npos) return recoPtBinsR2Cent50to90;
  }
  else if(rVal == 3){
    if(centStr.find("Cent0to10") != std::string::npos) return recoPtBinsR3Cent0to10;
    else if(centStr.find("Cent10to30") != std::string::npos) return recoPtBinsR3Cent10to30;
    else if(centStr.find("Cent30to50") != std::string::npos) return recoPtBinsR3Cent30to50;
    else if(centStr.find("Cent50to90") != std::string::npos) return recoPtBinsR3Cent50to90;
  }
  else if(rVal == 4){
    if(centStr.find("Cent0to10") != std::string::npos) return recoPtBinsR4Cent0to10;
    else if(centStr.find("Cent10to30") != std::string::npos) return recoPtBinsR4Cent10to30;
    else if(centStr.find("Cent30to50") != std::string::npos) return recoPtBinsR4Cent30to50;
    else if(centStr.find("Cent50to90") != std::string::npos) return recoPtBinsR4Cent50to90;
  }
  else if(rVal == 6){
    if(centStr.find("Cent0to10") != std::string::npos) return recoPtBinsR6Cent0to10;
    else if(centStr.find("Cent10to30") != std::string::npos) return recoPtBinsR6Cent10to30;
    else if(centStr.find("Cent30to50") != std::string::npos) return recoPtBinsR6Cent30to50;
    else if(centStr.find("Cent50to90") != std::string::npos) return recoPtBinsR6Cent50to90;
  }
  else if(rVal == 8){
    if(centStr.find("Cent0to10") != std::string::npos) return recoPtBinsR8Cent0to10;
    else if(centStr.find("Cent10to30") != std::string::npos) return recoPtBinsR8Cent10to30;
    else if(centStr.find("Cent30to50") != std::string::npos) return recoPtBinsR8Cent30to50;
    else if(centStr.find("Cent50to90") != std::string::npos) return recoPtBinsR8Cent50to90;
  }
  else if(rVal == 10){
    if(centStr.find("Cent0to10") != std::string::npos) return recoPtBinsR10Cent0to10;
    else if(centStr.find("Cent10to30") != std::string::npos) return recoPtBinsR10Cent10to30;
    else if(centStr.find("Cent30to50") != std::string::npos) return recoPtBinsR10Cent30to50;
    else if(centStr.find("Cent50to90") != std::string::npos) return recoPtBinsR10Cent50to90;
  }

  return {};
}


void smallOrLargeR::GetGeneralBins(Int_t nBins, Double_t bins[])
{
  for(Int_t bI = 0; bI < nBins+1; ++bI){
    bins[bI] = generalBins[bI];
  }
  return;
}

bool smallOrLargeR::CheckIntVect(std::vector<int> vect1, std::vector<int> vect2)
{
  bool isGood = vect1.size() == vect2.size();
  if(isGood){
    for(unsigned int i = 0; i < vect2.size(); ++i){
      if(vect2.at(i) != vect1.at(i)){
	isGood = false;
	break;
      }
    }
  }
  return isGood;
}

bool smallOrLargeR::CheckNBins(int inNBins, int compNBins, std::string nBinsStr)
{
  bool isGood = inNBins == compNBins;
  if(!isGood) std::cout << "WARNING: " << nBinsStr << " propagated \'" << inNBins << "\' doesn't match rReader \'" << compNBins << "\'. return false" << std::endl;
  return isGood;
}

void smallOrLargeR::SizeIsLTOrEQ(int inNBins, std::vector<double> inBins, std::string nBinsStr, std::string binsStr)
{
  if(((unsigned int)inNBins) > inBins.size()){
    std::cout << "Warning in smallOrLargeR: " << nBinsStr << ", " << inNBins << ", is incompatible w/ " << binsStr << ", " << inBins.size() << "." << std::endl;
  }
  return;
}


bool smallOrLargeR::CheckDoubleVect(std::vector<double> vect1, std::vector<double> vect2)
{
  bool isGood = vect1.size() == vect2.size();
  if(isGood){
    for(unsigned int i = 0; i < vect2.size(); ++i){
      if(TMath::Abs(vect2.at(i) - vect1.at(i)) >= deltaBin){
	isGood = false;
	break;
      }
    }
  }
  return isGood;
}


#endif
