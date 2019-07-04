#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "TH2D.h"

void csvCheck(std::string csvFileName, TH2D** hist_p, std::string name)
{
  std::ifstream file(csvFileName.c_str());
  std::string tempStr;
  
  const Int_t nMaxBins = 100;
  Double_t nBins = -1;
  Double_t bins[nMaxBins+1];

  std::vector<std::vector<std::string> > lines;
  while(std::getline(file, tempStr)){
    if(tempStr.size() == 0) continue;
    tempStr = tempStr + ",";
    while(tempStr.find(",,") != std::string::npos){
      tempStr.replace(tempStr.find(",,"), 2, ",");
    }
    if(tempStr.size() == 1) continue;
    
    lines.push_back({});
    while(tempStr.find(",") != std::string::npos){
      lines[lines.size()-1].push_back(tempStr.substr(0, tempStr.find(",")));
      tempStr.replace(0, tempStr.find(",")+1, "");
    }
  }

  file.close();

  std::vector<double> tempBins;
  for(unsigned int bI = 0; bI < lines.size(); ++bI){
    //    std::cout << "'" << lines[bI][0] << "'" << std::endl;
    double binVal = std::stod(lines[bI][1]);
    bool found = false;
    for(unsigned int tI = 0; tI < tempBins.size(); ++tI){
      if(TMath::Abs(tempBins[tI] - binVal) < 1.){
	found = true;
	break;
      }
    }    

    if(!found) tempBins.push_back(binVal);
  }

  std::sort(std::begin(tempBins), std::end(tempBins));
  
  for(unsigned int tI = 0; tI < tempBins.size(); ++tI){
    bins[tI] = tempBins[tI];
    ++nBins;
  }

  for(Int_t gI = 0; gI < nBins+1; ++gI){
    std::cout << bins[gI] << std::endl;
  }
  
  (*hist_p) = new TH2D(name.c_str(), "", nBins, bins, nBins, bins);
  
  for(unsigned int lI = 0; lI < lines.size(); ++lI){
    Int_t binX = -1;
    Int_t binY = -1;
    
    Double_t binValX = std::stod(lines[lI][1]);
    Double_t binValY = std::stod(lines[lI][3]);

    Double_t val = std::stod(lines[lI][5]);

    for(Int_t bI = 0; bI < nBins; ++bI){
      if(TMath::Abs(bins[bI] - binX) < 1.) binX = bI;
      if(TMath::Abs(bins[bI] - binY) < 1.) binY = bI;
    }

    (*hist_p)->SetBinContent(binX+1, binY+1, val);
  }

  //  std::cout << hist_p << std::endl;

  return;
}
