#include <iostream>
#include <string>

#include "TFile.h"
#include "TTree.h"

int getWeirdJets(std::string inFileName)
{
  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");

  inFile_p->Close();
  delete inFile_p;

  return 0;
}

int main()
{

  int retVal = 0;
  retVal += getWeightJets(argv[1]);
  return retVal;
}
