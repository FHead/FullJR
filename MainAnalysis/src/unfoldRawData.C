//cpp dependencies
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

//ROOT dependencies
#include "TCanvas.h"
#include "TDatime.h"
#include "TDirectory.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMatrix.h"
#include "TPad.h"
#include "TRandom3.h"
#include "TStyle.h"
#include "TTree.h"

//RooUnfold dependencies
#include "src/RooUnfoldBayes.h"
#include "src/RooUnfoldResponse.h"
#include "src/RooUnfoldSvd.h"

//Local dependencies
#include "MainAnalysis/include/canvNDCToXY.h"
#include "MainAnalysis/include/cutPropagator.h"
#include "MainAnalysis/include/macroHistToSubsetHist.h"
#include "MainAnalysis/include/reweightPrior.h"
#include "MainAnalysis/include/smallOrLargeR.h"
#include "MainAnalysis/include/texSlideCreator.h"

//Non-local FullJR dependencies (Utility, etc.)
#include "Utility/include/checkMakeDir.h"
#include "Utility/include/cppWatch.h"
#include "Utility/include/getLinBins.h"
#include "Utility/include/histDefUtility.h"
#include "Utility/include/lumiAndTAAUtil.h"
#include "Utility/include/mntToXRootdFileString.h"
#include "Utility/include/plotUtilities.h"
#include "Utility/include/returnRootFileContentsList.h"
#include "Utility/include/vanGoghPalette.h"

void correctForFakes(TH1D* rawHist_p, std::vector<double> binEdges, std::vector<double> fakeFactor)
{
   double deltaBinEdge = 0.1;

   std::vector<int> corrPos;
   for(unsigned int bI = 0; bI < binEdges.size()-1; ++bI){
      int tempPos = -1;

      for(Int_t bIX = 0; bIX < rawHist_p->GetNbinsX(); ++bIX){
         if(TMath::Abs(binEdges[bI] - rawHist_p->GetBinLowEdge(bIX+1)) > deltaBinEdge) continue;
         if(TMath::Abs(binEdges[bI+1] - rawHist_p->GetBinLowEdge(bIX+2)) > deltaBinEdge) continue;

         tempPos = bIX+1;
         break;
      }

      if(tempPos == -1){
         std::cout << "Warning! - corrPos = -1, mismatch in binnings." << std::endl;
         std::cout << " Input fake bin: " << binEdges[bI] << "-" << binEdges[bI+1] << std::endl;
         std::cout << " Histogram options: ";
         for(Int_t bIX = 0; bIX < rawHist_p->GetNbinsX(); ++bIX){
            std::cout << rawHist_p->GetBinLowEdge(bIX+1) << ", ";
         }
         std::cout << std::endl;
      }
      corrPos.push_back(tempPos);
   }

   for(unsigned int cI = 0; cI < corrPos.size(); ++cI){
      if(corrPos[cI] == -1) continue;

      double binVal = rawHist_p->GetBinContent(corrPos[cI]);
      binVal -= binVal*fakeFactor[cI];
      double binErr = TMath::Sqrt(binVal);

      rawHist_p->SetBinContent(corrPos[cI], binVal);
      rawHist_p->SetBinError(corrPos[cI], binErr);
   }

   return;
}

void getPearsTMatrix(RooUnfoldBayes* bayes_p, TH2D** covarianceToPlot, TH2D** pearsonToPlot)
{
   std::cout << "LINE: " << __LINE__ << ", " << bayes_p << std::endl;
   TMatrixD tempCovBayes = (TMatrixD)bayes_p->Ereco(RooUnfold::kCovToy);
   //  TMatrixD tempCovBayes = (TMatrixD)bayes_p->Ereco(RooUnfold::kErrors);
   std::cout << "LINE: " << __LINE__ << std::endl;
   TMatrixD* pearsonCoefsBayes = (TMatrixD*)tempCovBayes.Clone("pearsonCoefsBayes");

   std::cout << "LINE: " << __LINE__ << std::endl;

   for(Int_t rI = 0; rI < pearsonCoefsBayes->GetNrows(); ++rI){
      for(Int_t cI = 0; cI < pearsonCoefsBayes->GetNcols(); ++cI){
         Double_t valBayes = tempCovBayes(rI, cI);
         bool isGoodDiag = tempCovBayes(rI, rI) > 0.0 && tempCovBayes(cI, cI) > 0.0;
         bool isGoodVal = valBayes > 0.0;

         std::cout << "LINE: " << __LINE__ << std::endl;
         if(isGoodDiag) valBayes /= TMath::Sqrt(tempCovBayes(rI, rI)*tempCovBayes(cI, cI));
         else if(!isGoodVal) valBayes = 0;
         else{
            std::cout << "Warning diag is zero but off diag val is non-zero" << std::endl;
         }

         (*(pearsonCoefsBayes))(rI, cI) = valBayes;
      }
   }

   std::cout << "LINE: " << __LINE__ << std::endl;

   if(*covarianceToPlot == NULL) (*covarianceToPlot) = new TH2D(tempCovBayes);
   else{
      TH2D* tempCov = new TH2D(tempCovBayes);
      for(Int_t bIX = 0; bIX < tempCov->GetXaxis()->GetNbins(); ++bIX){
         for(Int_t bIY = 0; bIY < tempCov->GetYaxis()->GetNbins(); ++bIY){
            Double_t val = tempCov->GetBinContent(bIX+1, bIY+1);

            std::cout << "LINE: " << __LINE__ << std::endl;

            (*covarianceToPlot)->SetBinContent(bIX+1, bIY+1, val);
            (*covarianceToPlot)->SetBinError(bIX+1, bIY+1, 0.0);
         }
      }
   }
   (*pearsonToPlot) = new TH2D(*pearsonCoefsBayes);

   std::cout << "LINE: " << __LINE__ << std::endl;

   return;
}


void solve2By2(Double_t row1ans, Double_t row1col1, Double_t row1col2, Double_t row2ans, Double_t row2col1, Double_t row2col2, Double_t* solcol1, Double_t* solcol2)
{
   (*solcol1) = (row1ans*row2col2 - row2ans*row2col1)/(row1col1*row2col2 - row1col2*row2col1);
   (*solcol2) = (row1ans - (*solcol1)*row1col1)/row2col1;

   return;
}

void do2By2Reweight(TH1D* data_p, TH2D* responseMatrix_p, int nToys)
{
   std::cout << nToys << std::endl;
   /*
      std::vector<Double_t> factors;

      for(Int_t bIY = 0; bIY < responseMatrix_p->GetYaxis()->GetNbins(); ++bIY){
      Double_t total = 0.0;
      for(Int_t bIX = 0; bIX < responseMatrix_p->GetXaxis()->GetNbins(); ++bIX){
      total += responseMatrix_p->GetBinContent(bIX+1, bIY+1);
      }

      for(Int_t bIX = 0; bIX < responseMatrix_p->GetXaxis()->GetNbins(); ++bIX){
      Double_t val = responseMatrix_p->GetBinContent(bIX+1, bIY+1)/total;
      Double_t err = responseMatrix_p->GetBinError(bIX+1, bIY+1)/total;

      responseMatrix_p->SetBinContent(bIX+1, bIY+1, val);
      responseMatrix_p->SetBinError(bIX+1, bIY+1, err);
      }
      }

      for(Int_t bIY = 0; bIY < responseMatrix_p->GetYaxis()->GetNbins(); ++bIY){
      Double_t sol1 = data_p->GetBinContent(bIY+1);
      Double_t sol2 = data_p->GetBinContent(bIY+2);

      for(unsigned int fI = 0; fI < factors.size(); ++fI){
      sol1 -= factors[fI]*responseMatrix_p->GetBinContent(bIY+1, fI+1);
      sol2 -= factors[fI]*responseMatrix_p->GetBinContent(bIY+2, fI+1);
      }

      Double_t tempFactorA = -1;
      Double_t tempFactorB = -1;
      Double_t row1col1 = responseMatrix_p->GetBinContent(factors.size()+1, bIY+1);
      Double_t row1col2 = responseMatrix_p->GetBinContent(factors.size()+2, bIY+1);
      Double_t row2col1 = responseMatrix_p->GetBinContent(factors.size()+1, bIY+2);
      Double_t row2col2 = responseMatrix_p->GetBinContent(factors.size()+2, bIY+2);

      solve2By2(sol1, row1col1, row1col2, sol2, row2col1, row2col2, &tempFactorA, &tempFactorB);

      factors.push_back(tempFactorA);
      if(bIY == responseMatrix_p->GetYaxis()->GetNbins()-2){
      factors.push_back(tempFactorB);
      break;
      }
      }

      for(Int_t bIY = 0; bIY < responseMatrix_p->GetYaxis()->GetNbins(); ++bIY){
      for(Int_t bIX = 0; bIX < responseMatrix_p->GetXaxis()->GetNbins(); ++bIX){
      Double_t val = responseMatrix_p->GetBinContent(bIX+1, bIY+1)*factors[bIY];
      Double_t err = responseMatrix_p->GetBinError(bIX+1, bIY+1)*factors[bIY];

      responseMatrix_p->SetBinContent(bIX+1, bIY+1, val);
      responseMatrix_p->SetBinError(bIX+1, bIY+1, err);
      }
      }
      */


   const Int_t nBinsMax = 100;
   Int_t nGenBins = -1;
   Int_t nRecoBins = -1;;
   Double_t genBins[nBinsMax];
   Double_t recoBins[nBinsMax];

   for(Int_t bI = 0; bI < responseMatrix_p->GetYaxis()->GetNbins()+1; ++bI){
      genBins[bI] = responseMatrix_p->GetYaxis()->GetBinLowEdge(bI+1);
      ++nGenBins;
   }

   for(Int_t bI = 0; bI < data_p->GetNbinsX()+1; ++bI){
      recoBins[bI] = data_p->GetBinLowEdge(bI+1);
      ++nRecoBins;
   }

   TH1D* reco_p = new TH1D("reco_p", "", nRecoBins, recoBins);
   TH1D* gen_p = new TH1D("gen_p", "", nGenBins, genBins);

   macroHistToSubsetHistX(responseMatrix_p, reco_p, true);
   macroHistToSubsetHistY(responseMatrix_p, gen_p, true);

   RooUnfoldResponse* rooRes_p = NULL;
   RooUnfoldBayes* bayes_p = NULL;

   for(Int_t bI = 0; bI < 2; ++bI){
      rooRes_p = new RooUnfoldResponse(reco_p, gen_p, responseMatrix_p, "rooRes_p");
      bayes_p = new RooUnfoldBayes(rooRes_p, data_p, 10, false, "temp");
      bayes_p->SetVerbose(-1);
      bayes_p->SetNToys(nToys);
      TH1D* unfold_h = (TH1D*)bayes_p->Hreco(RooUnfold::kCovToy);
      //TH1D* unfold_h = (TH1D*)bayes_p->Hreco(RooUnfold::kErrors);

      for(Int_t bIY = 0; bIY < responseMatrix_p->GetYaxis()->GetNbins(); ++bIY){
         Double_t total = 0.0;
         for(Int_t bIX = 0; bIX < responseMatrix_p->GetXaxis()->GetNbins(); ++bIX){
            total += responseMatrix_p->GetBinContent(bIX+1, bIY+1);
         }

         if(total <= TMath::Power(10,-20)) continue;
         total = unfold_h->GetBinContent(bIY+1)/total;

         for(Int_t bIX = 0; bIX < responseMatrix_p->GetXaxis()->GetNbins(); ++bIX){
            responseMatrix_p->SetBinContent(bIX+1, bIY+1, responseMatrix_p->GetBinContent(bIX+1, bIY+1)*total);
            responseMatrix_p->SetBinError(bIX+1, bIY+1, responseMatrix_p->GetBinError(bIX+1, bIY+1)*total);
         }
      }

      delete unfold_h;
      delete rooRes_p;
      delete bayes_p;

      responseMatrix_p->Scale(data_p->Integral()/responseMatrix_p->Integral());

      macroHistToSubsetHistX(responseMatrix_p, reco_p, true);
      macroHistToSubsetHistY(responseMatrix_p, gen_p, true);
   }

   rooRes_p = new RooUnfoldResponse(reco_p, gen_p, responseMatrix_p, "rooRes_p");
   bayes_p = new RooUnfoldBayes(rooRes_p, data_p, 10, false, "temp");
   bayes_p->SetVerbose(-1);
   bayes_p->SetNToys(nToys);

   TH1D* unfold_h = (TH1D*)bayes_p->Hreco(RooUnfold::kCovToy);
   //TH1D* unfold_h = (TH1D*)bayes_p->Hreco(RooUnfold::kErrors);

   for(Int_t bIY = 0; bIY < responseMatrix_p->GetYaxis()->GetNbins(); ++bIY){
      Double_t total = 0.0;
      for(Int_t bIX = 0; bIX < responseMatrix_p->GetXaxis()->GetNbins(); ++bIX){
         total += responseMatrix_p->GetBinContent(bIX+1, bIY+1);
      }
      if(total <= TMath::Power(10,-20)) continue;
      total = unfold_h->GetBinContent(bIY+1)/total;

      for(Int_t bIX = 0; bIX < responseMatrix_p->GetXaxis()->GetNbins(); ++bIX){
         responseMatrix_p->SetBinContent(bIX+1, bIY+1, responseMatrix_p->GetBinContent(bIX+1, bIY+1)*total);
         responseMatrix_p->SetBinError(bIX+1, bIY+1, responseMatrix_p->GetBinError(bIX+1, bIY+1)*total);
      }
   }

   responseMatrix_p->Scale(data_p->Integral()/responseMatrix_p->Integral());
   macroHistToSubsetHistY(responseMatrix_p, gen_p, true);


   delete reco_p;
   delete gen_p;
   delete rooRes_p;
   delete bayes_p;
   delete unfold_h;

   return;
}

void doUnfold(TH2D* response_p, TH1D* prior_p, TH1D* data_p, Int_t bIter, RooUnfoldResponse** rooRes_p, RooUnfoldBayes** bayes_p, TH1D* outUnfold_p, int nToys, bool doSuper = false)
{
   std::cout << nToys << std::endl;

   const Int_t nBinsMax = 100;
   Int_t nGenBins = -1;
   Int_t nRecoBins = -1;
   Double_t genBins[nBinsMax];
   Double_t recoBins[nBinsMax];

   for(Int_t bI = 0; bI < prior_p->GetNbinsX()+1; ++bI){
      genBins[bI] = prior_p->GetBinLowEdge(bI+1);
      ++nGenBins;
   }

   for(Int_t bI = 0; bI < data_p->GetNbinsX()+1; ++bI){
      recoBins[bI] = data_p->GetBinLowEdge(bI+1);
      ++nRecoBins;
   }

   TH2D* responseReduced_p = new TH2D("responseReduced_h", "", nRecoBins, recoBins, nGenBins, genBins);
   macroHistToSubsetHist(response_p, responseReduced_p, true, true, true);

   for(Int_t bIY = 0; bIY < responseReduced_p->GetYaxis()->GetNbins(); ++bIY){
      Double_t total = 0.0;
      for(Int_t bIX = 0; bIX < responseReduced_p->GetXaxis()->GetNbins(); ++bIX){
         total += responseReduced_p->GetBinContent(bIX+1, bIY+1);
      }
      if(total <= TMath::Power(10,-20)) continue;

      Double_t scale = prior_p->GetBinContent(bIY+1)/total;

      for(Int_t bIX = 0; bIX < responseReduced_p->GetXaxis()->GetNbins(); ++bIX){
         Double_t val = responseReduced_p->GetBinContent(bIX+1, bIY+1)*scale;
         Double_t err = responseReduced_p->GetBinError(bIX+1, bIY+1)*scale;

         responseReduced_p->SetBinContent(bIX+1, bIY+1, val);
         responseReduced_p->SetBinError(bIX+1, bIY+1, err);
      }
   }

   responseReduced_p->Scale(data_p->Integral()/responseReduced_p->Integral());

   TH1D* reco_p = new TH1D("reco_p", "", nRecoBins, recoBins);
   TH1D* gen_p = new TH1D("gen_p", "", nGenBins, genBins);

   macroHistToSubsetHistX(responseReduced_p, reco_p, true);
   macroHistToSubsetHistY(responseReduced_p, gen_p, true);
   TH1D* unfold_h = NULL;

   if(doSuper){
      for(Int_t bI = 0; bI < 1; ++bI){
         (*rooRes_p) = new RooUnfoldResponse(reco_p, gen_p, responseReduced_p, "rooRes_p");

         (*bayes_p) = new RooUnfoldBayes(*rooRes_p, data_p, 3, false, "temp");
         (*bayes_p)->SetVerbose(-1);
         (*bayes_p)->SetNToys(nToys);
         unfold_h = (TH1D*)(*bayes_p)->Hreco(RooUnfold::kCovToy);
         //      unfold_h = (TH1D*)(*bayes_p)->Hreco(RooUnfold::kErrors);

         delete (*rooRes_p);
         delete (*bayes_p);

         for(Int_t bIY = 0; bIY < responseReduced_p->GetYaxis()->GetNbins(); ++bIY){
            Double_t total = 0.0;
            for(Int_t bIX = 0; bIX < responseReduced_p->GetXaxis()->GetNbins(); ++bIX){
               total += responseReduced_p->GetBinContent(bIX+1, bIY+1);
            }

            if(total <= TMath::Power(10,-20)) continue;
            Double_t scale = unfold_h->GetBinContent(bIY+1)/total;

            for(Int_t bIX = 0; bIX < responseReduced_p->GetXaxis()->GetNbins(); ++bIX){
               Double_t val = responseReduced_p->GetBinContent(bIX+1, bIY+1)*scale;
               Double_t err = responseReduced_p->GetBinError(bIX+1, bIY+1)*scale;

               responseReduced_p->SetBinContent(bIX+1, bIY+1, val);
               responseReduced_p->SetBinError(bIX+1, bIY+1, err);
            }
         }

         delete unfold_h;
         responseReduced_p->Scale(data_p->Integral()/responseReduced_p->Integral());

         macroHistToSubsetHistX(responseReduced_p, reco_p, true);
         macroHistToSubsetHistY(responseReduced_p, gen_p, true);
      }
   }

   (*rooRes_p) = new RooUnfoldResponse(reco_p, gen_p, responseReduced_p, "rooRes_p");

   (*bayes_p) = new RooUnfoldBayes(*rooRes_p, data_p, bIter, false, "temp");
   (*bayes_p)->SetVerbose(-1);
   (*bayes_p)->SetNToys(nToys);
   unfold_h = (TH1D*)(*bayes_p)->Hreco(RooUnfold::kCovToy);
   //unfold_h = (TH1D*)(*bayes_p)->Hreco(RooUnfold::kErrors);

   for(Int_t bIX = 0; bIX < unfold_h->GetNbinsX(); ++bIX){
      outUnfold_p->SetBinContent(bIX+1, unfold_h->GetBinContent(bIX+1));
      outUnfold_p->SetBinError(bIX+1, unfold_h->GetBinError(bIX+1));
   }

   delete unfold_h;

   delete reco_p;
   delete gen_p;
   delete responseReduced_p;

   return;
}

void doUnfoldSvd(TH2D* response_p, TH1D* prior_p, TH1D* data_p, Int_t kReg, RooUnfoldResponse** rooRes_p, RooUnfoldSvd** svd_p, TH1D* outUnfold_p, int nToys)
{
   std::cout << nToys << std::endl;

   const Int_t nBinsMax = 100;
   Int_t nGenBins = -1;
   Int_t nRecoBins = -1;
   Double_t genBins[nBinsMax];
   Double_t recoBins[nBinsMax];

   for(Int_t bI = 0; bI < prior_p->GetNbinsX()+1; ++bI){
      genBins[bI] = prior_p->GetBinLowEdge(bI+1);
      ++nGenBins;
   }

   for(Int_t bI = 0; bI < data_p->GetNbinsX()+1; ++bI){
      recoBins[bI] = data_p->GetBinLowEdge(bI+1);
      ++nRecoBins;
   }

   TH2D* responseReduced_p = new TH2D("responseReduced_h", "", nRecoBins, recoBins, nGenBins, genBins);
   macroHistToSubsetHist(response_p, responseReduced_p, true, true, true);

   for(Int_t bIY = 0; bIY < responseReduced_p->GetYaxis()->GetNbins(); ++bIY){
      Double_t total = 0.0;
      for(Int_t bIX = 0; bIX < responseReduced_p->GetXaxis()->GetNbins(); ++bIX){
         total += responseReduced_p->GetBinContent(bIX+1, bIY+1);
      }
      if(total <= TMath::Power(10,-20)) continue;

      Double_t scale = prior_p->GetBinContent(bIY+1)/total;

      for(Int_t bIX = 0; bIX < responseReduced_p->GetXaxis()->GetNbins(); ++bIX){
         Double_t val = responseReduced_p->GetBinContent(bIX+1, bIY+1)*scale;
         Double_t err = responseReduced_p->GetBinError(bIX+1, bIY+1)*scale;

         responseReduced_p->SetBinContent(bIX+1, bIY+1, val);
         responseReduced_p->SetBinError(bIX+1, bIY+1, err);
      }
   }

   responseReduced_p->Scale(data_p->Integral()/responseReduced_p->Integral());

   TH1D* reco_p = new TH1D("reco_p", "", nRecoBins, recoBins);
   TH1D* gen_p = new TH1D("gen_p", "", nGenBins, genBins);

   macroHistToSubsetHistX(responseReduced_p, reco_p, true);
   macroHistToSubsetHistY(responseReduced_p, gen_p, true);
   TH1D* unfold_h = NULL;
   (*rooRes_p) = new RooUnfoldResponse(reco_p, gen_p, responseReduced_p, "rooRes_p");

   (*svd_p) = new RooUnfoldSvd(*rooRes_p, data_p, kReg);
   (*svd_p)->SetVerbose(-1);
   (*svd_p)->SetNToys(nToys);
   unfold_h = (TH1D*)(*svd_p)->Hreco(RooUnfold::kCovToy);
   //unfold_h = (TH1D*)(*svd_p)->Hreco(RooUnfold::kErrors);

   for(Int_t bIX = 0; bIX < unfold_h->GetNbinsX(); ++bIX){
      outUnfold_p->SetBinContent(bIX+1, unfold_h->GetBinContent(bIX+1));
      outUnfold_p->SetBinError(bIX+1, unfold_h->GetBinError(bIX+1));
   }

   delete unfold_h;

   delete reco_p;
   delete gen_p;
   delete responseReduced_p;

   return;
}


void setPrior(TH2D* res_p, TH1D* prior_p)
{
   for(Int_t bIY = 0; bIY < res_p->GetYaxis()->GetNbins(); ++bIY){
      Double_t total = 0.0;
      for(Int_t bIX = 0; bIX < res_p->GetXaxis()->GetNbins(); ++bIX){
         total += res_p->GetBinContent(bIX+1, bIY+1);
      }

      if(total <= TMath::Power(10, -20)) continue;

      total = prior_p->GetBinContent(bIY+1)/total;
      for(Int_t bIX = 0; bIX < res_p->GetXaxis()->GetNbins(); ++bIX){
         Double_t val = res_p->GetBinContent(bIX+1, bIY+1)*total;
         Double_t err = res_p->GetBinError(bIX+1, bIY+1)*total;

         res_p->SetBinContent(bIX+1, bIY+1, val);
         res_p->SetBinError(bIX+1, bIY+1, err);
      }
   }

   return;
}


void setPrior(TH2D* res_p, TH2D* prior_p)
{
   const Int_t nMaxBins = 100;
   Int_t nBins = res_p->GetYaxis()->GetNbins();
   Double_t bins[nMaxBins+1];
   for(Int_t bIY = 0; bIY < nBins+1; ++bIY){
      bins[bIY] = prior_p->GetYaxis()->GetBinLowEdge(bIY+1);
   }

   TH1D* priorTemp_p = new TH1D("priorTemp_h", "", nBins, bins);
   for(Int_t bIY = 0; bIY < prior_p->GetYaxis()->GetNbins(); ++bIY){
      Double_t total = 0.0;
      for(Int_t bIX = 0; bIX < prior_p->GetXaxis()->GetNbins(); ++bIX){
         total += prior_p->GetBinContent(bIX+1, bIY+1);
      }

      priorTemp_p->SetBinContent(bIY+1, total);
      priorTemp_p->SetBinError(bIY+1, 0.0);
   }
   setPrior(res_p, priorTemp_p);
   delete priorTemp_p;

   return;
}

void cleanMatrix(TH2D* inMat_p)
{
   /*
      for(Int_t bIX = 0; bIX < inMat_p->GetXaxis()->GetNbins(); ++bIX){
      Double_t total = 0.0;
      for(Int_t bIY = 0; bIY < inMat_p->GetYaxis()->GetNbins(); ++bIY){
      total += inMat_p->GetBinContent(bIX+1, bIY+1);
      }

      if(total <= TMath::power(10, -20)) continue;

      Int_t breakPos = -1;
      Double_t total2 = 0.0;
      for(Int_t bIY = 0; bIY < inMat_p->GetYaxis()->GetNbins(); ++bIY){
      total2 += inMat_p->GetBinContent(bIX+1, bIY+1);

      if(total2/total >= 0.01){
      breakPos = bIY;
      break;
      }
      }
      for(Int_t bIY = 0; bIY < breakPos; ++bIY){
      inMat_p->SetBinContent(bIX+1, bIY+1, 0.0);
      inMat_p->SetBinErro(bIX+1, bIY+1, 0.0);
      }
      }
      */

   for(Int_t bIX = 0; bIX < inMat_p->GetXaxis()->GetNbins(); ++bIX){ 
      for(Int_t bIY = 0; bIY < inMat_p->GetYaxis()->GetNbins(); ++bIY){
         Double_t val = inMat_p->GetBinContent(bIX+1, bIY+1);

         if(val <= TMath::Power(10, -20)) continue;

         Double_t err = inMat_p->GetBinError(bIX+1, bIY+1);    

         if(err/val >= 1./TMath::Sqrt(10.)){
            inMat_p->SetBinContent(bIX+1, bIY+1, 0.0);
            inMat_p->SetBinError(bIX+1, bIY+1, 0.0);
         }
      }        
   }

   return;
}

double getVectMean(std::vector<double> vect)
{
   double mean = 0.0;
   for(unsigned int mI = 0; mI < vect.size(); ++mI){
      mean += vect[mI];
   }
   mean /= (double)vect.size();
   return mean;
}

double getVectSigma(std::vector<double> vect)
{
   double mean = getVectMean(vect);
   double sigma = 0.0;
   for(unsigned int mI = 0; mI < vect.size(); ++mI){
      sigma += (vect[mI] - mean)*(vect[mI] - mean);
   }
   return TMath::Sqrt(sigma/(double)(vect.size()-1));
}

void ratioPrint(TH1D* hist1_p, TH1D* hist2_p)
{
   for(Int_t bIX = 0; bIX < hist1_p->GetNbinsX(); ++bIX){
      Double_t val = hist1_p->GetBinContent(bIX+1)/hist2_p->GetBinContent(bIX+1);
      std::cout << prettyString(val, 4, false) << ", ";
   }
   std::cout << std::endl;


   return;
}

void drawResponse(TH2D* res_p, std::string saveStr)
{
   const Int_t nMaxBins = 500;
   Int_t nXBins = res_p->GetXaxis()->GetNbins();
   Int_t nYBins = res_p->GetYaxis()->GetNbins();
   Double_t xBins[nMaxBins+1];
   Double_t yBins[nMaxBins+1];

   for(Int_t xI = 0; xI < nXBins+1; ++xI){
      xBins[xI] = res_p->GetXaxis()->GetBinLowEdge(xI+1);
   }
   for(Int_t yI = 0; yI < nYBins+1; ++yI){
      yBins[yI] = res_p->GetYaxis()->GetBinLowEdge(yI+1);
   }

   TH2D* resNorm_p = new TH2D("resNorm_p", ";Reconstructed p_{T} (GeV/c);Truth p_{T} (GeV/C)", nXBins, xBins, nYBins, yBins);
   TH2D* resXNorm_p = new TH2D("resXNorm_p", ";Reconstructed p_{T} (GeV/c) (Norm to Unity in Y);Truth p_{T} (GeV/C)", nXBins, xBins, nYBins, yBins);
   TH2D* resYNorm_p = new TH2D("resYNorm_p", ";Reconstructed p_{T} (GeV/c);Truth p_{T} (GeV/C) (Norm to Unity in X)", nXBins, xBins, nYBins, yBins);

   centerTitles({resNorm_p, resXNorm_p, resYNorm_p});

   macroHistToSubsetHist(res_p, resNorm_p, true);
   macroHistToSubsetHist(res_p, resXNorm_p, true);
   macroHistToSubsetHist(res_p, resYNorm_p, true);

   for(Int_t xI = 0; xI < nXBins; ++xI){
      Double_t total = 0.0;

      for(Int_t yI = 0; yI < nYBins; ++yI){
         total += resXNorm_p->GetBinContent(xI+1, yI+1);
      }

      for(Int_t yI = 0; yI < nYBins; ++yI){
         Double_t val = resXNorm_p->GetBinContent(xI+1, yI+1)/total;
         Double_t err = resXNorm_p->GetBinError(xI+1, yI+1)/total;

         resXNorm_p->SetBinContent(xI+1, yI+1, val);
         resXNorm_p->SetBinError(xI+1, yI+1, err);
      }
   }

   for(Int_t yI = 0; yI < nYBins; ++yI){
      Double_t total = 0.0;

      for(Int_t xI = 0; xI < nXBins; ++xI){
         total += resYNorm_p->GetBinContent(xI+1, yI+1);
      }

      for(Int_t xI = 0; xI < nXBins; ++xI){
         Double_t val = resYNorm_p->GetBinContent(xI+1, yI+1)/total;
         Double_t err = resYNorm_p->GetBinError(xI+1, yI+1)/total;

         resYNorm_p->SetBinContent(xI+1, yI+1, val);
         resYNorm_p->SetBinError(xI+1, yI+1, err);
      }
   }

   TCanvas* canv_p = new TCanvas("canv_p", "", 450*3, 450);
   canv_p->Divide(3,1);

   canv_p->cd();
   canv_p->SetTopMargin(0.01);
   canv_p->SetLeftMargin(0.01);
   canv_p->SetBottomMargin(0.01);
   canv_p->SetRightMargin(0.01);

   canv_p->cd(1);
   gPad->SetTopMargin(0.10);
   gPad->SetLeftMargin(0.10);
   gPad->SetBottomMargin(0.10);
   gPad->SetRightMargin(0.10);

   TLatex* label_p = new TLatex();
   label_p->SetNDC(0);
   label_p->SetTextFont(43);
   label_p->SetTextSize(18);

   resNorm_p->Scale(1000./resNorm_p->GetBinContent(1, 1));
   resNorm_p->SetMinimum(0.0005);
   resNorm_p->DrawCopy("COLZ TEXT");
   gStyle->SetPaintTextFormat("1.3f");

   /*  
       for(Int_t bIX = 0; bIX < nXBins; ++bIX){
       for(Int_t bIY = 0; bIY < nYBins; ++bIY){
       Double_t centY = (yBins[bIY] + yBins[bIY+1])/2.;
       Double_t centX = (xBins[bIX]*3./4. + xBins[bIX+1]/4.);

       Double_t val = resNorm_p->GetBinContent(bIX+1, bIY+1);
       std::string tempStr = "";
       if(val >= 1000) tempStr = prettyString(val, 1, false);
       else if(val >= 100) tempStr = prettyString(val, 2, false);
       else if(val >= 10) tempStr = prettyString(val, 3, false);
       else if(val >= 1) tempStr = prettyString(val, 4, false);
       else tempStr = prettyString(val, 5, false);

       label_p->DrawLatex(centX, centY, tempStr.c_str());
       }
       }
       */

   gStyle->SetOptStat(0);
   gPad->SetLogx();
   gPad->SetLogy();
   gPad->SetLogz();
   //  gStyle->SetPaintTextFormat("1.3f");

   canv_p->cd(2);
   gPad->SetTopMargin(0.10);
   gPad->SetLeftMargin(0.10);
   gPad->SetBottomMargin(0.10);
   gPad->SetRightMargin(0.10);

   resXNorm_p->DrawCopy("COLZ TEXT");
   gStyle->SetOptStat(0);
   gPad->SetLogx();
   gPad->SetLogy();
   gStyle->SetPaintTextFormat("1.3f");

   canv_p->cd(3);
   gPad->SetTopMargin(0.10);
   gPad->SetLeftMargin(0.10);
   gPad->SetBottomMargin(0.10);
   gPad->SetRightMargin(0.10);

   resYNorm_p->DrawCopy("COLZ TEXT");
   gStyle->SetOptStat(0);
   gPad->SetLogx();
   gPad->SetLogy();
   gStyle->SetPaintTextFormat("1.3f");

   canv_p->cd();

   std::string rStr = res_p->GetName();
   std::string centStr = rStr;
   std::string etaStr = rStr;
   std::string pytHyd = "";
   //  std::string centStr = rStr;

   if(rStr.find("Cs") != std::string::npos){
      rStr.replace(0, rStr.find("akCs")+4, "");
      rStr.replace(rStr.find("PU"), rStr.size(), "");

      centStr.replace(0, centStr.find("Cent")+4, "");
      centStr.replace(centStr.find("_"), centStr.size(), "");
      centStr.replace(centStr.find("to"), 2, "-");
      centStr = centStr + "%";
      pytHyd = "PYTHIA6+HYDJET";
   }
   else{
      rStr.replace(0, rStr.find("ak")+2, "");
      rStr.replace(rStr.find("PF"), rStr.size(), "");

      centStr = "PP";
      pytHyd = "PYTHIA6";
   }

   if(rStr.size() == 2) rStr = "R=" + rStr.substr(0, 1) + "." + rStr.substr(1, 1);
   else rStr = "R=0." + rStr.substr(0, 1);

   etaStr.replace(0, etaStr.find("AbsEta")+6, "");
   etaStr.replace(etaStr.find("_"), etaStr.size(), "");

   while(etaStr.find("p") != std::string::npos){
      etaStr.replace(etaStr.find("p"), 1, ".");
   }

   etaStr.replace(etaStr.find("to"), 2, " < |#eta| < ");

   label_p->SetNDC();
   label_p->DrawLatex(0.1, 0.95, rStr.c_str());
   label_p->DrawLatex(0.15, 0.95, centStr.c_str());
   label_p->DrawLatex(0.2, 0.95, etaStr.c_str());
   label_p->DrawLatex(0.3, 0.95, pytHyd.c_str());
   //  label_p->DrawLatex(0.5, 0.95, "RESPONSE MATRIX");

   delete label_p;

   quietSaveAs(canv_p, saveStr);
   delete canv_p;
   delete resNorm_p;
   delete resXNorm_p;
   delete resYNorm_p;

   return;
}


int svdTerm(std::string jetStr, std::string centStr)
{
   int bestSvd = - 1;

   if(jetStr.find("akCs10") != std::string::npos){
      if(centStr.find("Cent0to10") != std::string::npos) bestSvd = 3;
      else if(centStr.find("Cent10to30") != std::string::npos) bestSvd = 4;
      else if(centStr.find("Cent30to50") != std::string::npos) bestSvd = 4;
      else if(centStr.find("Cent50to90") != std::string::npos) bestSvd = 4;
   }
   else if(jetStr.find("akCs8") != std::string::npos){
      if(centStr.find("Cent0to10") != std::string::npos) bestSvd = 4;
      else if(centStr.find("Cent10to30") != std::string::npos) bestSvd = 4;
      else if(centStr.find("Cent30to50") != std::string::npos) bestSvd = 4;
      else if(centStr.find("Cent50to90") != std::string::npos) bestSvd = 4;
   }
   else if(jetStr.find("akCs6") != std::string::npos){
      if(centStr.find("Cent0to10") != std::string::npos) bestSvd = 4;
      else if(centStr.find("Cent10to30") != std::string::npos) bestSvd = 5;
      else if(centStr.find("Cent30to50") != std::string::npos) bestSvd = 4;
      else if(centStr.find("Cent50to90") != std::string::npos) bestSvd = 4;
   }
   else if(jetStr.find("akCs4") != std::string::npos){
      if(centStr.find("Cent0to10") != std::string::npos) bestSvd = 4;
      else if(centStr.find("Cent10to30") != std::string::npos) bestSvd = 5;
      else if(centStr.find("Cent30to50") != std::string::npos) bestSvd = 4;
      else if(centStr.find("Cent50to90") != std::string::npos) bestSvd = 4;
   }
   else if(jetStr.find("akCs3") != std::string::npos){
      if(centStr.find("Cent0to10") != std::string::npos) bestSvd = 4;
      else if(centStr.find("Cent10to30") != std::string::npos) bestSvd = 6;
      else if(centStr.find("Cent30to50") != std::string::npos) bestSvd = 4;
      else if(centStr.find("Cent50to90") != std::string::npos) bestSvd = 4;
   }
   else if(jetStr.find("akCs2") != std::string::npos){
      if(centStr.find("Cent0to10") != std::string::npos) bestSvd = 4;
      else if(centStr.find("Cent10to30") != std::string::npos) bestSvd = 6;
      else if(centStr.find("Cent30to50") != std::string::npos) bestSvd = 4;
      else if(centStr.find("Cent50to90") != std::string::npos) bestSvd = 4;
   }

   return bestSvd;
}


int unfoldRawData(const std::string inDataFileName, const std::string inResponseName, const std::string selectJtAlgo = "", std::string overrideFile = "", const bool doClean = false, const int nToys = 1000)
{
   std::string cleanStr = "Clean";
   if(!doClean) cleanStr = "NoClean";

   std::string toyStr = "NToy" + std::to_string(nToys);

   std::string addedTagStr = "";

   cppWatch totalRunWatch;
   cppWatch doUnfoldWatch;
   Int_t unfoldCounter = 0;

   totalRunWatch.start();

   TRandom3* randGen_p = new TRandom3(0);

   std::vector<int> overrideInts;
   std::vector<std::string> overrideCents;  
   std::vector<std::vector<double> > overrideGenBins;
   std::vector<std::vector<double> > overrideRecoBins;

   bool doOverride = overrideFile.size() != 0 && overrideFile.find(".txt") != std::string::npos && checkFile(overrideFile);

   if(doOverride){
      std::ifstream inFile(overrideFile.c_str());
      std::string tempStr;

      while(std::getline(inFile, tempStr)){
         if(tempStr.size() == 0) continue;
         if(tempStr.substr(0, 1).find("#") != std::string::npos) continue;
         std::vector<std::string> subStrings;

         tempStr = tempStr + ",";
         while(tempStr.find(",,") != std::string::npos){
            tempStr.replace(tempStr.find(",,"), 2, ",");
         }
         while(tempStr.find(" ") != std::string::npos){
            tempStr.replace(tempStr.find(" "), 1, "");
         }
         if(tempStr.size() == 1) continue;

         while(tempStr.find(",") != std::string::npos){
            subStrings.push_back(tempStr.substr(0, tempStr.find(",")));
            tempStr.replace(0, tempStr.find(",")+1, "");
         }

         for(unsigned int sI = 0; sI < subStrings.size(); ++sI){
            std::cout << subStrings[sI] << ",";
         }
         std::cout << std::endl;

         if(subStrings[0].find("TAG") != std::string::npos){
            addedTagStr = subStrings[1] + "_";
            continue;
         }

         int rVal = std::stoi(subStrings[1]);
         std::string centStr = subStrings[2];

         int pos = -1;
         for(unsigned int rI = 0; rI < overrideInts.size(); ++rI){
            if(rVal != overrideInts[rI]) continue;
            if(!isStrSame(overrideCents[rI], centStr)) continue;

            pos = rI;
            break;
         }

         if(pos < 0){
            overrideInts.push_back(rVal);
            overrideCents.push_back(centStr);
         }

         std::vector<double> bins;
         for(unsigned int rI = 3; rI < subStrings.size(); ++rI){
            bins.push_back(std::stod(subStrings[rI]));
         }

         if(subStrings[0].find("Gen") != std::string::npos){
            if(pos < 0){
               overrideGenBins.push_back(bins); 
               overrideRecoBins.push_back({}); 
            }
            else overrideGenBins[pos] = bins;
         }
         else{
            if(pos < 0){
               overrideRecoBins.push_back(bins); 
               overrideGenBins.push_back({}); 
            }
            else overrideRecoBins[pos] = bins;
         }
      } 
      inFile.close();
   }

   vanGoghPalette vg;
   const Int_t nXVals = 5;
   const Int_t xVals[nXVals] = {200, 400, 600, 800, 1000};

   TDatime* date = new TDatime();
   const std::string dateStr = std::to_string(date->GetDate());
   const std::string dateStr2 = std::to_string(date->GetYear()) + "." + std::to_string(date->GetMonth()) + "." + std::to_string(date->GetDay());
   delete date;

   TFile* responseFile_p = new TFile(inResponseName.c_str(), "READ");
   std::vector<std::string> responseJetDirList = returnRootFileContentsList(responseFile_p, "TDirectoryFile", "JetAnalyzer");
   std::cout << "Printing " << responseJetDirList.size() << " response jets..." << std::endl;

   checkMakeDir("pdfDir");
   checkMakeDir("pdfDir/" + dateStr);
   checkMakeDir("pdfDir/" + dateStr + "/Unfold");

   //editing
   std::vector<std::vector<std::string> > slideTitlesPerAlgo;
   std::vector<std::vector<std::vector<std::string> > > pdfPerSlidePerAlgo;

   for(unsigned int jI = 0; jI < responseJetDirList.size(); ++jI){
      std::cout << " " << jI << "/" << responseJetDirList.size() << ": " << responseJetDirList[jI] << std::endl;

      std::string tempStr = responseJetDirList[jI];
      if(tempStr.find("/") != std::string::npos) tempStr.replace(tempStr.find("/"), tempStr.size() - tempStr.find("/"), "");

      checkMakeDir("pdfDir/" + dateStr + "/Unfold_" + tempStr);
      slideTitlesPerAlgo.push_back({});
      pdfPerSlidePerAlgo.push_back({});
   }

   cutPropagator cutPropResponse;
   cutPropResponse.Clean();
   cutPropResponse.GetAllVarFromFile(responseFile_p);

   responseFile_p->Close();
   delete responseFile_p;

   TFile* dataFile_p = new TFile(inDataFileName.c_str(), "READ");
   std::vector<std::string> dataJetDirList = returnRootFileContentsList(dataFile_p, "TDirectoryFile", "JetAnalyzer");

   cutPropagator cutPropData;
   cutPropData.Clean();
   cutPropData.GetAllVarFromFile(dataFile_p);

   dataFile_p->Close();
   delete dataFile_p;

   if(doOverride){
      for(unsigned int oI = 0; oI < overrideInts.size(); ++oI){
         cutPropData.SetGenPtBins(overrideInts[oI], overrideCents[oI], overrideGenBins[oI]);
         cutPropData.SetRecoPtBins(overrideInts[oI], overrideCents[oI], overrideRecoBins[oI]);

         cutPropResponse.SetGenPtBins(overrideInts[oI], overrideCents[oI], overrideGenBins[oI]);
         cutPropResponse.SetRecoPtBins(overrideInts[oI], overrideCents[oI], overrideRecoBins[oI]);
      }
   }

   if(!cutPropResponse.CheckPropagatorsMatch(cutPropData, false, true)){
      std::cout << "unfoldRawData - Cuts listed in data file \'" << inDataFileName << "\' and response file \'" << inResponseName << "\' do not match. return 1" << std::endl;
      //    return 1;
   }

   Int_t valForForLoops = 100000000;
   if(doLocalDebug || doGlobalDebug){
      std::cout << "DOLOCALDEBUG or DOGLOBALDEBUG in unfoldRawData: Setting all histogram array sizes to 1 for faster processing" << std::endl;
      valForForLoops = 1000000;
   }

   const Int_t isDataPP = cutPropData.GetIsPP();
   const Int_t isResponsePP = cutPropResponse.GetIsPP();

   const Int_t nMaxSyst = 20;
   Int_t nSyst = TMath::Min(valForForLoops, cutPropResponse.GetNSyst());
   std::vector<std::string> systStr = cutPropResponse.GetSystStr();


   systStr.push_back("MatrixStatUp");
   systStr.push_back("MatrixStatDown");
   ++nSyst;
   ++nSyst;

   //  ++nSyst;
   //  systStr.push_back("MatrixStat");

   Int_t mainPos = -1;
   std::vector<Int_t> priorPos;
   Int_t matrixStatPos = -1;
   Int_t matrixStatUpPos = -1;
   Int_t matrixStatDownPos = -1;

   for(unsigned int sI = 0; sI < systStr.size(); ++sI){
      if(systStr[sI].size() == 0) mainPos = sI;
      else if(systStr[sI].find("Prior") != std::string::npos) priorPos.push_back(sI);
      else if(isStrSame(systStr[sI], "MatrixStatUp")) matrixStatUpPos = sI;
      else if(isStrSame(systStr[sI], "MatrixStatDown")) matrixStatDownPos = sI;
      else if(isStrSame(systStr[sI], "MatrixStat")){
         systStr[sI] = "MatrixStatRandom";
         matrixStatPos = sI;
      }
   }

   if(nSyst > nMaxSyst){
      std::cout << "nSyst \'" << nSyst << "\' is greater than nMaxSyst \'" << nMaxSyst << "\'. return 1" << std::endl;
      return 1;
   }

   const Int_t nResponseMod = 1;
   std::vector<double> responseMod = {0.10};
   std::vector<double> jerVarData = {0.10};

   /*
      const Int_t nResponseMod = TMath::Min(valForForLoops, cutPropResponse.GetNResponseMod());
      std::vector<double> responseMod = cutPropResponse.GetResponseMod();
      std::vector<double> jerVarData = cutPropResponse.GetJERVarData();
      */

   const Int_t nMaxCentBins = 4;
   const Int_t nCentBins = TMath::Min(valForForLoops, cutPropData.GetNCentBins());
   if(nCentBins > nMaxCentBins){
      std::cout << "nCentBins \'" << nCentBins << "\' is greater than nMaxCentBins \'" << nMaxCentBins << "\'. return 1" << std::endl;
      return 1;
   }

   std::vector<Int_t> centBinsLow = cutPropData.GetCentBinsLow();
   std::vector<Int_t> centBinsHi = cutPropData.GetCentBinsHi();

   const Int_t nMaxJtPtBins = 50;

   const Int_t nJtAbsEtaBins = 1;
   std::vector<Double_t> jtAbsEtaBinsLowTemp = {0.0};
   std::vector<Double_t> jtAbsEtaBinsHiTemp = {2.0};

   /*
      const Int_t nJtAbsEtaBins = TMath::Min(valForForLoops, cutPropData.GetNJtAbsEtaBins());
      std::vector<Double_t> jtAbsEtaBinsLowTemp = cutPropData.GetJtAbsEtaBinsLow();
      std::vector<Double_t> jtAbsEtaBinsHiTemp = cutPropData.GetJtAbsEtaBinsHi();
      */

   /*
      const Int_t nID = TMath::Min(valForForLoops, cutPropData.GetNID());
      std::vector<std::string> idStr = cutPropData.GetIdStr();
      */

   const Int_t nID = 1;
   std::vector<std::string> idStr = {"LightMUAndCHID"};

   Double_t jtAbsEtaBinsLow[nJtAbsEtaBins];
   Double_t jtAbsEtaBinsHi[nJtAbsEtaBins];
   std::cout << "nJtAbsEtaBins: ";
   for(Int_t jI = 0; jI < nJtAbsEtaBins; ++jI){
      jtAbsEtaBinsLow[jI] = jtAbsEtaBinsLowTemp[jI];
      jtAbsEtaBinsHi[jI] = jtAbsEtaBinsHiTemp[jI];
      std::cout << " " << jtAbsEtaBinsLow[jI] << "-" << jtAbsEtaBinsHi[jI] << ",";
   }
   std::cout << std::endl;

   std::cout << "nCentBins: " << nCentBins << std::endl;
   for(Int_t cI = 0; cI < nCentBins; ++cI){
      std::cout << " " << cI << "/" << nCentBins << ": " << centBinsLow[cI] << "-" << centBinsHi[cI] << std::endl;
   }

   std::cout << "Raw data from file: \'" << inDataFileName << "\'" << std::endl;
   std::cout << "Response data from file: \'" << inResponseName << "\'" << std::endl;

   //Reduce to match jet dirs
   unsigned int pos = 0;
   while(pos < dataJetDirList.size()){
      bool isFound = false;
      for(unsigned int i = 0; i < responseJetDirList.size(); ++i){
         if(dataJetDirList[pos].size() != responseJetDirList[i].size()) continue;
         if(dataJetDirList[pos].find(responseJetDirList[i]) == std::string::npos) continue;

         isFound = true;
         break;
      }

      if(isFound) ++pos;
      else dataJetDirList.erase(dataJetDirList.begin()+pos);
   }

   pos = 0;
   while(pos < responseJetDirList.size()){
      bool isFound = false;
      for(unsigned int i = 0; i < dataJetDirList.size(); ++i){
         if(responseJetDirList[pos].size() != dataJetDirList[i].size()) continue;
         if(responseJetDirList[pos].find(dataJetDirList[i]) == std::string::npos) continue;

         isFound = true;
         break;
      }

      if(isFound) ++pos;
      else responseJetDirList.erase(responseJetDirList.begin()+pos);
   }

   std::cout << "Shared jets to process: " << std::endl;
   for(unsigned int i = 0; i < responseJetDirList.size(); ++i){
      std::cout << " " << i << "/" << responseJetDirList.size() << ": " << responseJetDirList[i] << std::endl;
   }

   if(selectJtAlgo.size() != 0){
      std::cout << "Restricting to \'" << selectJtAlgo << "\'." << std::endl;
      unsigned int pos = 0;
      while(pos < responseJetDirList.size()){
         if(responseJetDirList[pos].find(selectJtAlgo) == std::string::npos) responseJetDirList.erase(responseJetDirList.begin()+pos);
         else ++pos;
      }

      if(responseJetDirList.size() == 0){
         std::cout << "No jet dirs after selection for " << selectJtAlgo << std::endl;
         return 1;
      }
   }

   smallOrLargeR rReader;

   const Int_t nMaxDataJet = 10;
   const Int_t nDataJet = responseJetDirList.size();

   if(nDataJet > nMaxDataJet){
      std::cout << "nDataJet \'" << nDataJet << "\' is greater than nMaxDataJet \'" << nMaxDataJet << "\'. return 1" << std::endl;
      return 1;
   }

   std::string outFileName = inDataFileName;
   while(outFileName.find("/") != std::string::npos){outFileName.replace(0, outFileName.find("/")+1, "");}
   if(outFileName.find(".txt") != std::string::npos) outFileName.replace(outFileName.find(".txt"), std::string(".txt").size(), "");
   else if(outFileName.find(".root") != std::string::npos) outFileName.replace(outFileName.find(".root"), std::string(".root").size(), "");

   std::string debugStr = "";
   if(doLocalDebug || doGlobalDebug) debugStr = "DEBUG_";

   std::string selectJtAlgoStr = selectJtAlgo;
   if(selectJtAlgoStr.size() != 0) selectJtAlgoStr = selectJtAlgoStr + "_";

   checkMakeDir("output");
   checkMakeDir("output/" + dateStr);

   outFileName = outFileName.substr(0, outFileName.find("ProcessRawData"));
   outFileName = "output/" + dateStr + "/" + outFileName + "_UnfoldRawData_" + selectJtAlgoStr + cleanStr + "_" + toyStr + "_" + debugStr + addedTagStr + dateStr + ".root";

   while(outFileName.find("__") != std::string::npos){outFileName.replace(outFileName.find("__"), 2, "_");}

   TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");
   //Following two lines necessary else you get absolutely clobbered on deletion timing. See threads here:
   //https://root-forum.cern.ch/t/closing-root-files-is-dead-slow-is-this-a-normal-thing/5273/16
   //https://root-forum.cern.ch/t/tfile-speed/17549/25
   //Bizarre
   outFile_p->SetBit(TFile::kDevNull);
   TH1::AddDirectory(kFALSE);

   const Int_t nSvd = 10;

   const Int_t nBayes = 14;
   const Int_t nBigBayesSymm = 3;
   const Int_t nBayesBig = nBayes - (2*nBigBayesSymm + 1);
   Int_t temp100Pos = -1;
   Int_t bayesVal[nBayes];
   for(Int_t bI = 0; bI < nBayes; ++bI){
      //    if(bI >= nBayesBig) bayesVal[bI] = 100 + 1 + nBigBayesSymm - (nBayes - bI);
      //    else bayesVal[bI] = bI+1;
      bayesVal[bI] = bI+1;

      //    if(bayesVal[bI] == 100) temp100Pos = bI;
      if(bayesVal[bI] == nBayes-nBigBayesSymm) temp100Pos = bI;
   }


   const Int_t bayes100Pos = temp100Pos;

   std::cout << "Bayes: " << std::endl;
   for(Int_t bI = 0; bI < nBayes; ++bI){
      std::cout << " " << bI << "/" << nBayes << ": " << bayesVal[bI] << std::endl;
   }

   std::cout << "LINE: " << __LINE__ << std::endl;

   const Int_t nHistDim = nDataJet*nCentBins*nID*nResponseMod*nJtAbsEtaBins*nSyst;
   std::vector<std::string> histTag;
   std::vector<int> histBestSvd;
   std::vector<double> histBestSvdProb;
   std::vector<int> histBestBayes;
   std::vector<int> histBestBayesPos;

   std::cout << "LINE: " << __LINE__ << std::endl;

   cutPropData.SetNBayes(nBayes);
   cutPropData.SetNSVD(nSvd);
   cutPropData.SetNBigBayesSymm(nBigBayesSymm);
   cutPropData.SetBayesVal(nBayes, bayesVal);

   cutPropData.SetNResponseMod(nResponseMod);
   cutPropData.SetResponseMod(responseMod);
   cutPropData.SetJERVarData(jerVarData);

   cutPropData.SetNSyst(nSyst);
   cutPropData.SetSystStr(systStr);

   std::cout << "LINE: " << __LINE__ << std::endl;

   int nGenJtPtBins[nMaxDataJet][nMaxCentBins];
   double genJtPtBins[nMaxDataJet][nMaxCentBins][nMaxJtPtBins+1];
   int nRecoJtPtBins[nMaxDataJet][nMaxCentBins];
   double recoJtPtBins[nMaxDataJet][nMaxCentBins][nMaxJtPtBins+1];

   std::cout << "LINE: " << __LINE__ << std::endl;

   //  int nRecoJtPtBins2[nMaxDataJet][nMaxCentBins];
   //  double recoJtPtBins2[nMaxDataJet][nMaxCentBins][nMaxJtPtBins+1];

   double ptBinLow[nMaxDataJet][nMaxCentBins];
   double ptBinHigh[nMaxDataJet][nMaxCentBins];  

   std::cout << "LINE: " << __LINE__ << std::endl;

   for(Int_t jI = 0; jI < nDataJet; ++jI){
      std::string tempStr = responseJetDirList[jI];
      if(tempStr.find("/") != std::string::npos) tempStr.replace(tempStr.find("/"), tempStr.size() - tempStr.find("/"), "");
      const Int_t rVal = getRVal(tempStr);

      for(Int_t cI = 0; cI < nCentBins; ++cI){
         std::string centStr = "Cent" + std::to_string(centBinsLow[cI]) + "to" + std::to_string(centBinsHi[cI]);
         nGenJtPtBins[jI][cI] = cutPropData.GetGenNBinsFromRValCent(rVal, centStr);
         cutPropData.GetGenBinsFromRValCent(rVal, centStr, genJtPtBins[jI][cI]);
         nRecoJtPtBins[jI][cI] = cutPropData.GetRecoNBinsFromRValCent(rVal, centStr);
         cutPropData.GetRecoBinsFromRValCent(rVal, centStr, recoJtPtBins[jI][cI]);
      }
   }

   std::cout << "LINE: " << __LINE__ << std::endl;

   /*
      for(Int_t jI = 0; jI < nDataJet; ++jI){
      for(Int_t cI = 0; cI < nCentBins; ++cI){

      for(Int_t rI = 1; rI < nRecoJtPtBins[jI][cI]; ++rI){
      recoJtPtBins2[jI][cI][rI-1] = recoJtPtBins[jI][cI][rI];
      }
      nRecoJtPtBins2[jI][cI] = nRecoJtPtBins[jI][cI] - 2;


      }
      */

   std::cout << "LINE: " << __LINE__ << std::endl;

   const Int_t nBayesDraw = TMath::Min(valForForLoops, 4);
   TDirectory* dir_p[nMaxDataJet];
   TH1D* jtPtUnfoldedBayes_RecoGenAsymm_h[nMaxDataJet][nMaxCentBins][nID][nResponseMod][nJtAbsEtaBins][nMaxSyst][nBayes];
   TH1D* jtPtUnfoldedBayesMC_RecoGenAsymm_h[nMaxDataJet][nMaxCentBins][nID][nResponseMod][nJtAbsEtaBins][nMaxSyst][nBayes];
   TH1D* jtPtUnfoldedMCTruth_RecoGenAsymm_h[nMaxDataJet][nMaxCentBins][nID][nResponseMod][nJtAbsEtaBins][nMaxSyst];
   RooUnfoldBayes* bayes_p[nMaxDataJet][nMaxCentBins][nID][nResponseMod][nJtAbsEtaBins][nMaxSyst][nBayes];
   RooUnfoldBayes* bayesMC_p[nMaxDataJet][nMaxCentBins][nID][nResponseMod][nJtAbsEtaBins][nMaxSyst][nBayes];

   TH1D* jtPtUnfoldedSvd_RecoGenAsymm_h[nMaxDataJet][nMaxCentBins][nID][nResponseMod][nJtAbsEtaBins][nMaxSyst][nSvd];
   TH1D* jtPtUnfoldedSvdMC_RecoGenAsymm_h[nMaxDataJet][nMaxCentBins][nID][nResponseMod][nJtAbsEtaBins][nMaxSyst][nSvd];
   RooUnfoldSvd* svd_p[nMaxDataJet][nMaxCentBins][nID][nResponseMod][nJtAbsEtaBins][nMaxSyst][nSvd];
   RooUnfoldSvd* svdMC_p[nMaxDataJet][nMaxCentBins][nID][nResponseMod][nJtAbsEtaBins][nMaxSyst][nSvd];
   Double_t lineChi2Svd[nMaxDataJet][nMaxCentBins][nID][nResponseMod][nJtAbsEtaBins][nMaxSyst][nSvd];

   std::cout << "LINE: " << __LINE__ << std::endl;

   Int_t histTermPos[nMaxDataJet][nMaxCentBins][nID][nResponseMod][nJtAbsEtaBins][nMaxSyst];
   const Int_t nCov = 3;
   TH2D* covariance_h[nMaxDataJet][nMaxCentBins][nID][nResponseMod][nJtAbsEtaBins][nCov];
   TH2D* covariance_Mean_h[nMaxDataJet][nMaxCentBins][nID][nResponseMod][nJtAbsEtaBins];
   TH2D* covariance_Sigma_h[nMaxDataJet][nMaxCentBins][nID][nResponseMod][nJtAbsEtaBins];
   TH2D* covariance_SigmaOverMean_h[nMaxDataJet][nMaxCentBins][nID][nResponseMod][nJtAbsEtaBins];

   TH2D* response_Unfold_h[nMaxDataJet][nMaxCentBins][nID][nResponseMod][nJtAbsEtaBins][nMaxSyst];

   std::cout << "LINE: " << __LINE__ << std::endl;

   for(Int_t jI = 0; jI < nDataJet; ++jI){
      std::string tempStr = responseJetDirList[jI];
      if(tempStr.find("/") != std::string::npos) tempStr.replace(tempStr.find("/"), tempStr.size() - tempStr.find("/"), "");
      dir_p[jI] = (TDirectory*)outFile_p->mkdir(tempStr.c_str());

      std::cout << "LINE: " << __LINE__ << std::endl;

      for(Int_t cI = 0; cI < nCentBins; ++cI){
         std::string centStr = "Cent" + std::to_string(centBinsLow[cI]) + "to" + std::to_string(centBinsHi[cI]);

         if(isDataPP) centStr = "PP_" + centStr;
         else centStr = "PbPb_" + centStr;

         for(Int_t idI = 0; idI < nID; ++idI){
            for(Int_t mI = 0; mI < nResponseMod; ++mI){
               const std::string resStr = "ResponseMod" + prettyString(responseMod[mI], 2, true);

               for(Int_t aI = 0; aI < nJtAbsEtaBins; ++aI){
                  const std::string jtAbsEtaStr = "AbsEta" + prettyString(jtAbsEtaBinsLow[aI], 1, true) + "to" + prettyString(jtAbsEtaBinsHi[aI], 1, true);


                  for(Int_t covI = 0; covI < nCov; ++covI){
                     std::string covName = "cov_" + tempStr + "_" + centStr + "_" + idStr[idI] + "_" + resStr + "_" + jtAbsEtaStr + "_Cov" + std::to_string(covI) + "_h";
                     covariance_h[jI][cI][idI][mI][aI][covI] = new TH2D(covName.c_str(), "", nGenJtPtBins[jI][cI], genJtPtBins[jI][cI], nGenJtPtBins[jI][cI], genJtPtBins[jI][cI]);
                  }
                  std::string covName = "cov_" + tempStr + "_" + centStr + "_" + idStr[idI] + "_" + resStr + "_" + jtAbsEtaStr + "_Mean_h";
                  covariance_Mean_h[jI][cI][idI][mI][aI] = new TH2D(covName.c_str(), "", nGenJtPtBins[jI][cI], genJtPtBins[jI][cI], nGenJtPtBins[jI][cI], genJtPtBins[jI][cI]);

                  covName = "cov_" + tempStr + "_" + centStr + "_" + idStr[idI] + "_" + resStr + "_" + jtAbsEtaStr + "_Sigma_h";
                  covariance_Sigma_h[jI][cI][idI][mI][aI] = new TH2D(covName.c_str(), "", nGenJtPtBins[jI][cI], genJtPtBins[jI][cI], nGenJtPtBins[jI][cI], genJtPtBins[jI][cI]);

                  covName = "cov_" + tempStr + "_" + centStr + "_" + idStr[idI] + "_" + resStr + "_" + jtAbsEtaStr + "_SigmaOverMean_h";
                  covariance_SigmaOverMean_h[jI][cI][idI][mI][aI] = new TH2D(covName.c_str(), "", nGenJtPtBins[jI][cI], genJtPtBins[jI][cI], nGenJtPtBins[jI][cI], genJtPtBins[jI][cI]);


                  std::cout << "LINE: " << __LINE__ << std::endl;

                  for(Int_t sI = 0; sI < nSyst; ++sI){
                     std::string tempSystStr = "_" + systStr[sI] + "_";
                     while(tempSystStr.find("__") != std::string::npos){tempSystStr.replace(tempSystStr.find("__"), 2, "_");}

                     histTermPos[jI][cI][idI][mI][aI][sI] = -1;

                     jtPtUnfoldedMCTruth_RecoGenAsymm_h[jI][cI][idI][mI][aI][sI] = new TH1D(("jtPtUnfoldedMCTruth_RecoGenAsymm_" + tempStr + "_" + centStr + "_" + idStr[idI] + "_" + resStr + "_" + jtAbsEtaStr + tempSystStr + "h").c_str(), ";Unfolded MC Jet p_{T};Counts", nGenJtPtBins[jI][cI], genJtPtBins[jI][cI]);

                     response_Unfold_h[jI][cI][idI][mI][aI][sI] = new TH2D(("response_" + tempStr + "_" + centStr + "_" + idStr[idI] + "_" + resStr + "_" + jtAbsEtaStr + tempSystStr + "Unfold_h").c_str(), "", nRecoJtPtBins[jI][cI], recoJtPtBins[jI][cI], nGenJtPtBins[jI][cI], genJtPtBins[jI][cI]);

                     for(Int_t bI = 0; bI < nBayes; ++bI){
                        std::string bayesStr = "Bayes" + std::to_string(bayesVal[bI]);

                        std::cout << "LINE: " << __LINE__ << std::endl;

                        jtPtUnfoldedBayes_RecoGenAsymm_h[jI][cI][idI][mI][aI][sI][bI] = new TH1D(("jtPtUnfoldedBayes_RecoGenAsymm_" + tempStr + "_" + centStr + "_" + idStr[idI] + "_" + resStr + "_" + jtAbsEtaStr + tempSystStr + bayesStr + "_h").c_str(), ";Unfolded Jet p_{T};Counts", nGenJtPtBins[jI][cI], genJtPtBins[jI][cI]);

                        jtPtUnfoldedBayesMC_RecoGenAsymm_h[jI][cI][idI][mI][aI][sI][bI] = new TH1D(("jtPtUnfoldedBayesMC_RecoGenAsymm_" + tempStr + "_" + centStr + "_" + idStr[idI] + "_" + resStr + "_" + jtAbsEtaStr + tempSystStr + bayesStr + "_h").c_str(), ";Unfolded MC Jet p_{T};Counts", nGenJtPtBins[jI][cI], genJtPtBins[jI][cI]);
                        std::vector<TH1*> tempVect = {jtPtUnfoldedBayes_RecoGenAsymm_h[jI][cI][idI][mI][aI][sI][bI], jtPtUnfoldedBayesMC_RecoGenAsymm_h[jI][cI][idI][mI][aI][sI][bI]};

                        centerTitles(tempVect);
                        setSumW2(tempVect);
                     }

                     for(Int_t bI = 0; bI < TMath::Min(nSvd, nGenJtPtBins[jI][cI]); ++bI){
                        std::string svdStr = "Svd" + std::to_string(bI+1);

                        jtPtUnfoldedSvd_RecoGenAsymm_h[jI][cI][idI][mI][aI][sI][bI] = new TH1D(("jtPtUnfoldedSvd_RecoGenAsymm_" + tempStr + "_" + centStr + "_" + idStr[idI] + "_" + resStr + "_" + jtAbsEtaStr + tempSystStr + svdStr + "_h").c_str(), ";Unfolded Jet p_{T};Counts", nGenJtPtBins[jI][cI], genJtPtBins[jI][cI]);

                        jtPtUnfoldedSvdMC_RecoGenAsymm_h[jI][cI][idI][mI][aI][sI][bI] = new TH1D(("jtPtUnfoldedSvdMC_RecoGenAsymm_" + tempStr + "_" + centStr + "_" + idStr[idI] + "_" + resStr + "_" + jtAbsEtaStr + tempSystStr + svdStr + "_h").c_str(), ";Unfolded MC Jet p_{T};Counts", nGenJtPtBins[jI][cI], genJtPtBins[jI][cI]);
                        std::vector<TH1*> tempVect = {jtPtUnfoldedSvd_RecoGenAsymm_h[jI][cI][idI][mI][aI][sI][bI], jtPtUnfoldedSvdMC_RecoGenAsymm_h[jI][cI][idI][mI][aI][sI][bI]};

                        centerTitles(tempVect);
                        setSumW2(tempVect);
                     }

                  }
               }
            }
         }
      }
   }

   std::cout << "LINE: " << __LINE__ << std::endl;

   responseFile_p = new TFile(inResponseName.c_str(), "READ");
   //  RooUnfoldResponse* rooResponse_RecoGenAsymm_h[nMaxDataJet][nMaxCentBins][nID][nResponseMod][nJtAbsEtaBins][nMaxSyst][nSmallLargeBins];
   RooUnfoldResponse* rooResponse_Reweighted_h[nMaxDataJet][nMaxCentBins][nID][nResponseMod][nJtAbsEtaBins][nMaxSyst];

   TH2D* response_RecoGenAsymm_h[nMaxDataJet][nMaxCentBins][nID][nResponseMod][nJtAbsEtaBins][nMaxSyst];  
   TH2D* response_General_h[nMaxDataJet][nMaxCentBins][nID][nResponseMod][nJtAbsEtaBins][nMaxSyst];
   TH1D* recoJtPt_GoodGen_h[nMaxDataJet][nMaxCentBins][nID][nResponseMod][nJtAbsEtaBins][nMaxSyst];
   //  TH1D* genJtPt_GoodReco_h[nMaxDataJet][nMaxCentBins][nID][nResponseMod][nJtAbsEtaBins][nMaxSyst];

   std::cout << "LINE: " << __LINE__ << std::endl;

   for(Int_t jI = 0; jI < nDataJet; ++jI){
      std::string tempStr = responseJetDirList[jI];
      if(tempStr.find("/") != std::string::npos) tempStr.replace(tempStr.find("/"), tempStr.size() - tempStr.find("/"), "");

      std::cout << "LINE: " << __LINE__ << std::endl;

      for(Int_t cI = 0; cI < nCentBins; ++cI){
         std::string centStr = "Cent" + std::to_string(centBinsLow[cI]) + "to" + std::to_string(centBinsHi[cI]);
         if(isResponsePP) centStr = "PP_" + centStr;
         else centStr = "PbPb_" + centStr;

         std::cout << "LINE: " << __LINE__ << std::endl;

         for(Int_t iI = 0; iI < nID; ++iI){
            for(Int_t mI = 0; mI < nResponseMod; ++mI){
               const std::string resStr = "ResponseMod" + prettyString(responseMod[mI], 2, true);

               for(Int_t aI = 0; aI < nJtAbsEtaBins; ++aI){
                  const std::string jtAbsEtaStr = "AbsEta" + prettyString(jtAbsEtaBinsLow[aI], 1, true) + "to" + prettyString(jtAbsEtaBinsHi[aI], 1, true);


                  std::cout << "LINE: " << __LINE__ << std::endl;
                  for(Int_t sI = 0; sI < nSyst; ++sI){
                     std::string tempSystStr = "_" + systStr[sI] + "_";
                     while(tempSystStr.find("__") != std::string::npos){tempSystStr.replace(tempSystStr.find("__"), 2, "_");}


                     std::cout << "LINE: " << __LINE__ << ", " << tempSystStr << std::endl;

                     if(matrixStatPos == sI || matrixStatUpPos == sI || matrixStatDownPos == sI){
                        recoJtPt_GoodGen_h[jI][cI][iI][mI][aI][sI] = (TH1D*)recoJtPt_GoodGen_h[jI][cI][iI][mI][aI][mainPos]->Clone(("recoJtPt_" + tempStr + "_" + centStr + "_" + idStr[iI] + "_" + resStr + "_" + jtAbsEtaStr + tempSystStr + "GoodGen_h").c_str());

                        response_General_h[jI][cI][iI][mI][aI][sI] = (TH2D*)response_General_h[jI][cI][iI][mI][aI][mainPos]->Clone(("response_" + tempStr + "_" + centStr + "_" + idStr[iI] + "_" + resStr + "_" + jtAbsEtaStr + tempSystStr + "General_h").c_str());

                        response_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI] = (TH2D*)response_RecoGenAsymm_h[jI][cI][iI][mI][aI][mainPos]->Clone(("response_" + tempStr + "_" + centStr + "_" + idStr[iI] + "_" + resStr + "_" + jtAbsEtaStr + tempSystStr + "RecoGenAsymm_h").c_str());
                     }
                     else{
                        recoJtPt_GoodGen_h[jI][cI][iI][mI][aI][sI] = (TH1D*)responseFile_p->Get((tempStr + "/recoJtPt_" + tempStr + "_" + centStr + "_" + idStr[iI] + "_" + resStr + "_" + jtAbsEtaStr + tempSystStr + "GoodGen_h").c_str());

                        response_General_h[jI][cI][iI][mI][aI][sI] = (TH2D*)responseFile_p->Get((tempStr + "/response_" + tempStr + "_" + centStr + "_" + idStr[iI] + "_" + resStr + "_" + jtAbsEtaStr + tempSystStr + "General_h").c_str());	     	      

                        response_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI] = (TH2D*)responseFile_p->Get((tempStr + "/response_" + tempStr + "_" + centStr + "_" + idStr[iI] + "_" + resStr + "_" + jtAbsEtaStr + tempSystStr + "RecoGenAsymm_h").c_str());	     	      

                     }	    
                  }
               }
            }
         }
      }
   }

   std::cout << "LINE: " << __LINE__ << std::endl;

   dataFile_p = new TFile(inDataFileName.c_str(), "READ");
   TH1D* jtPtRaw_RecoGenAsymm_h[nMaxDataJet][nMaxCentBins][nID][nJtAbsEtaBins];
   TH1D* jtPtRaw_General_h[nMaxDataJet][nMaxCentBins][nID][nJtAbsEtaBins];

   for(Int_t jI = 0; jI < nDataJet; ++jI){
      std::string tempStr = responseJetDirList[jI];
      if(tempStr.find("/") != std::string::npos) tempStr.replace(tempStr.find("/"), tempStr.size() - tempStr.find("/"), "");

      for(Int_t cI = 0; cI < nCentBins; ++cI){
         std::string centStr = "Cent" + std::to_string(centBinsLow[cI]) + "to" + std::to_string(centBinsHi[cI]);
         if(isDataPP) centStr = "PP_" + centStr;
         else centStr = "PbPb_" + centStr;

         for(Int_t iI = 0; iI < nID; ++iI){
            for(Int_t aI = 0; aI < nJtAbsEtaBins; ++aI){
               const std::string jtAbsEtaStr = "AbsEta" + prettyString(jtAbsEtaBinsLow[aI], 1, true) + "to" + prettyString(jtAbsEtaBinsHi[aI], 1, true);

               //	  const std::string name = tempStr + "/jtPtRaw_RecoGenAsymm_" + tempStr + "_" + centStr + "_" + idStr[iI] + "_" + jtAbsEtaStr + "_h";
               const std::string name = "jtPtRaw_RecoGenAsymm_" + tempStr + "_" + centStr + "_" + idStr[iI] + "_" + jtAbsEtaStr + "_h";
               const std::string nameGeneral = tempStr + "/jtPtRaw_General_" + tempStr + "_" + centStr + "_" + idStr[iI] + "_" + jtAbsEtaStr + "_h";

               //	  jtPtRaw_RecoGenAsymm_h[jI][cI][iI][aI] = (TH1D*)dataFile_p->Get(name.c_str());
               jtPtRaw_RecoGenAsymm_h[jI][cI][iI][aI] = new TH1D(name.c_str(), "", nRecoJtPtBins[jI][cI], recoJtPtBins[jI][cI]);
               jtPtRaw_General_h[jI][cI][iI][aI] = (TH1D*)dataFile_p->Get(nameGeneral.c_str());
            }
         }
      }
   }



   for(Int_t jI = 0; jI < nDataJet; ++jI){
      for(Int_t cI = 0; cI < nCentBins; ++cI){
         ptBinLow[jI][cI] = genJtPtBins[jI][cI][0];
         if(ptBinLow[jI][cI] < recoJtPtBins[jI][cI][0]) ptBinLow[jI][cI] = recoJtPtBins[jI][cI][0];

         ptBinHigh[jI][cI] = genJtPtBins[jI][cI][nGenJtPtBins[jI][cI]];
         if(ptBinHigh[jI][cI] > recoJtPtBins[jI][cI][nRecoJtPtBins[jI][cI]]) ptBinHigh[jI][cI] = recoJtPtBins[jI][cI][nRecoJtPtBins[jI][cI]];
      }
   }


   std::cout << "Start Unfolding..." << std::endl;
   for(Int_t jI = 0; jI < nDataJet; ++jI){
      std::string tempStr = responseJetDirList[jI];
      if(tempStr.find("/") != std::string::npos) tempStr.replace(tempStr.find("/"), tempStr.size() - tempStr.find("/"), "");

      std::cout << " Unfolding " << jI << "/" << nDataJet << ": " << tempStr << std::endl;

      const Int_t rVal = getRVal(tempStr);

      for(Int_t cI = 0; cI < nCentBins; ++cI){
         if(!isDataPP) std::cout << "  " << centBinsLow[cI] << "-" << centBinsHi[cI] << "%..." << std::endl;

         std::string centStr = "Cent" + std::to_string(centBinsLow[cI]) + "to" + std::to_string(centBinsHi[cI]);


         for(Int_t iI = 0; iI < nID; ++iI){
            for(Int_t mI = 0; mI < nResponseMod; ++mI){	 
               const std::string resStr = "ResponseMod" + prettyString(responseMod[mI], 2, true);
               for(Int_t aI = 0; aI < nJtAbsEtaBins; ++aI){
                  const std::string jtAbsEtaStr = "AbsEta" + prettyString(jtAbsEtaBinsLow[aI], 1, true) + "to" + prettyString(jtAbsEtaBinsHi[aI], 1, true);

                  for(Int_t sI = 0; sI < nSyst; ++sI){
                     std::string tempSystStr = "_" + systStr[sI] + "_";
                     while(tempSystStr.find("__") != std::string::npos){tempSystStr.replace(tempSystStr.find("__"), 2, "_");}

                     std::cout << "LINE: " << __LINE__ << ", " << systStr[sI] << std::endl;	      

                     Int_t sIPos = sI;
                     for(unsigned int pI = 0; pI < priorPos.size(); ++pI){
                        if(priorPos[pI] == sI){
                           sIPos = mainPos;
                           break;
                        }
                     }

                     std::string histName = jtPtUnfoldedBayes_RecoGenAsymm_h[jI][cI][iI][mI][aI][sIPos][0]->GetName();

                     //	      TH2D* initRes_p = (TH2D*)response_RecoGenAsymm_h[jI][cI][iI][mI][aI][sIPos]->Clone("initRes_p");
                     //	      TH2D* initRes_p = new TH2D("initRes_p", "", nRecoJtPtBins[jI][cI], recoJtPtBins[jI][cI], nGenJtPtBins[jI][cI], genJtPtBins[jI][cI]);
                     //	      response_Unfold_h[jI][cI][iI][mI][aI][sI] = new TH2D((tempStr + "/response_" + tempStr + "_" + centStr + "_" + idStr[iI] + "_" + resStr + "_" + jtAbsEtaStr + tempSystStr + "General_h").c_str(), "", nRecoJtPtBins[jI][cI], recoJtPtBins[jI][cI], nGenJtPtBins[jI][cI], genJtPtBins[jI][cI]);

                     std::cout << "LINE: " << __LINE__ << ", " << systStr[sI] << std::endl;	 	      
                     TH2D* initResMC_p = new TH2D("initResMC_p", "", nRecoJtPtBins[jI][cI], recoJtPtBins[jI][cI], nGenJtPtBins[jI][cI], genJtPtBins[jI][cI]);
                     TH2D* initResMC2_p = new TH2D("initResMC2_p", "", nRecoJtPtBins[jI][cI], recoJtPtBins[jI][cI], nGenJtPtBins[jI][cI], genJtPtBins[jI][cI]);

                     std::cout << "LINE: " << __LINE__ << ", " << systStr[sI] << std::endl;	 

                     if(doClean){
                        std::cout << "PRIORFIT" << std::endl;

                        const Int_t nPriorBins = 14;
                        Double_t priorBins[nPriorBins+1] = {300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000};

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

                        std::cout << "LINE: " << __LINE__ << ", " << systStr[sI] << std::endl;	 

                        TH1D* tempPrior_p = new TH1D("tempPrior_h", "", nPriorBins, priorBins);
                        macroHistToSubsetHistY(response_General_h[jI][cI][iI][mI][aI][sIPos], tempPrior_p, true);
                        std::cout << "LINE: " << __LINE__ << ", " << systStr[sI] << std::endl;	 
                        tempPrior_p->Print("ALL");
                        for(Int_t bIX = 0; bIX < nPriorBins; ++bIX){
                           Double_t width = tempPrior_p->GetBinWidth(bIX+1);
                           Double_t val = tempPrior_p->GetBinContent(bIX+1)/width;
                           Double_t err = tempPrior_p->GetBinError(bIX+1)/width;
                           tempPrior_p->SetBinContent(bIX+1, val);
                           tempPrior_p->SetBinError(bIX+1, err);
                        }

                        std::cout << "LINE: " << __LINE__ << ", " << systStr[sI] << std::endl;	 
                        TF1* priorFit_p = new TF1("priorFit_p", "[0]*TMath::Power(325/x, [1])", 300, 1000);
                        priorFit_p->SetParameter(0, tempPrior_p->GetBinContent(1));
                        priorFit_p->SetParameter(1, 8);

                        tempPrior_p->Fit("priorFit_p", "M", "", 300, 1000);
                        std::cout << "PARAMS: " << priorFit_p->GetParameter(0) << ", " << priorFit_p->GetParameter(1) << std::endl;

                        double newPar1 = getPower(rVal, centStr, isDataPP);
                        std::cout << "NEWPAR: " << newPar1 << std::endl;

                        std::cout << "LINE: " << __LINE__ << ", " << systStr[sI] << std::endl;	 
                        tempPrior_p->DrawCopy("HIST E1 P");
                        priorFit_p->Draw("SAME");
                        gStyle->SetOptStat(0);
                        gPad->SetLogy();

                        canv_p->cd(2);
                        gPad->SetTopMargin(0.01);
                        gPad->SetRightMargin(0.01);
                        gPad->SetLeftMargin(0.14);
                        gPad->SetBottomMargin(0.14);

                        Double_t newPar0 = tempPrior_p->GetBinContent(3)/TMath::Power(1./tempPrior_p->GetBinCenter(3), newPar1);

                        for(Int_t bIX = 0; bIX < nPriorBins; ++bIX){
                           Double_t val = tempPrior_p->GetBinContent(bIX+1);
                           Double_t err = tempPrior_p->GetBinError(bIX+1);
                           Double_t val2 = priorFit_p->GetParameter(0)*TMath::Power(325/tempPrior_p->GetBinCenter(bIX+1), priorFit_p->GetParameter(1));

                           tempPrior_p->SetBinContent(bIX+1, val/val2);
                           tempPrior_p->SetBinError(bIX+1, err/val2);
                        }

                        std::cout << "LINE: " << __LINE__ << ", " << systStr[sI] << std::endl;	 

                        tempPrior_p->DrawCopy("HIST E1 P");
                        gStyle->SetOptStat(0);

                        TF1* priorFitNew_p = new TF1("priorFit_p", "[0]*TMath::Power(1./x, [1])", 300, 1000);
                        priorFitNew_p->SetParameter(1, newPar1);

                        Double_t par1 = priorFit_p->GetParameter(1);

                        priorFitNew_p->SetParameter(0, newPar0);
                        priorFitNew_p->SetMarkerStyle(20);
                        priorFitNew_p->SetMarkerSize(1);
                        priorFitNew_p->SetMarkerColor(1);
                        priorFitNew_p->SetLineColor(1);
                        priorFitNew_p->SetLineStyle(2);
                        canv_p->cd(1);
                        priorFitNew_p->Draw("SAME");

                        const std::string saveName = "priorFit_" + tempStr + "_" + centStr + "_" + idStr[iI] + "_" + resStr + "_" + jtAbsEtaStr + "_" + tempSystStr +  "AllBayes_RecoGenAsymm_" + cleanStr + "_" + toyStr + "_" + debugStr + addedTagStr + dateStr + ".pdf";
                        canv_p->SaveAs(("pdfDir/" + dateStr + "/Unfold_" + tempStr + "/" + saveName).c_str());

                        delete canv_p;

                        //		Double_t pos300 = -1;
                        Double_t center300 = -1;
                        Double_t total300 = 0.0;
                        for(Int_t bIY = 0; bIY < response_General_h[jI][cI][iI][mI][aI][sIPos]->GetYaxis()->GetNbins(); ++bIY){
                           Double_t lowEdge = response_General_h[jI][cI][iI][mI][aI][sIPos]->GetYaxis()->GetBinLowEdge(bIY+1);
                           Double_t highEdge = response_General_h[jI][cI][iI][mI][aI][sIPos]->GetYaxis()->GetBinLowEdge(bIY+2);

                           if(lowEdge <= 301 && 301 < highEdge){
                              //		    pos300 = bIY;
                              center300 = response_General_h[jI][cI][iI][mI][aI][sIPos]->GetYaxis()->GetBinCenter(bIY+1);

                              for(Int_t bIX = 0; bIX < response_General_h[jI][cI][iI][mI][aI][sIPos]->GetXaxis()->GetNbins(); ++bIX){
                                 total300 += response_General_h[jI][cI][iI][mI][aI][sIPos]->GetBinContent(bIX+1, bIY+1);
                              }		    

                              break;
                           }
                        }

                        Double_t par0 = total300/TMath::Power(1./center300, par1);
                        newPar0 = total300/TMath::Power(1./center300, newPar1);

                        for(Int_t bIY = 0; bIY < response_General_h[jI][cI][iI][mI][aI][sIPos]->GetYaxis()->GetNbins(); ++bIY){
                           Double_t center = response_General_h[jI][cI][iI][mI][aI][sIPos]->GetYaxis()->GetBinCenter(bIY+1);
                           Double_t factor = newPar0*TMath::Power(1./center, newPar1)/(par0*TMath::Power(1./center, par1));

                           for(Int_t bIX = 0; bIX < response_General_h[jI][cI][iI][mI][aI][sIPos]->GetXaxis()->GetNbins(); ++bIX){
                              Double_t val = response_General_h[jI][cI][iI][mI][aI][sIPos]->GetBinContent(bIX+1, bIY+1)*factor;
                              Double_t err = response_General_h[jI][cI][iI][mI][aI][sIPos]->GetBinError(bIX+1, bIY+1)*factor;

                              response_General_h[jI][cI][iI][mI][aI][sIPos]->SetBinContent(bIX+1, bIY+1, val);
                              response_General_h[jI][cI][iI][mI][aI][sIPos]->SetBinError(bIX+1, bIY+1, err);
                           }
                        }		

                        delete tempPrior_p;
                        delete priorFit_p;
                        delete priorFitNew_p;
                     }

                     macroHistToSubsetHist(response_General_h[jI][cI][iI][mI][aI][sIPos], response_Unfold_h[jI][cI][iI][mI][aI][sI], true, false, true);
                     macroHistToSubsetHist(response_General_h[jI][cI][iI][mI][aI][sIPos], initResMC_p, true, false, true);
                     macroHistToSubsetHist(response_General_h[jI][cI][iI][mI][aI][sI], initResMC2_p, true, false, true);

                     if(sI == matrixStatPos || sI == matrixStatUpPos || sI == matrixStatDownPos){
                        Double_t minVal = 10000000000000;

                        for(Int_t bIX = 0; bIX < response_Unfold_h[jI][cI][iI][mI][aI][sI]->GetXaxis()->GetNbins(); ++bIX){
                           for(Int_t bIY = 0; bIY < response_Unfold_h[jI][cI][iI][mI][aI][sI]->GetYaxis()->GetNbins(); ++bIY){
                              if(response_Unfold_h[jI][cI][iI][mI][aI][sI]->GetBinContent(bIX+1, bIY+1) < minVal && response_Unfold_h[jI][cI][iI][mI][aI][sI]->GetBinContent(bIX+1, bIY+1) > TMath::Power(10, -20)){
                                 minVal = response_Unfold_h[jI][cI][iI][mI][aI][sI]->GetBinContent(bIX+1, bIY+1);
                              }
                           }
                        }


                        for(Int_t bIX = 0; bIX < response_Unfold_h[jI][cI][iI][mI][aI][sI]->GetXaxis()->GetNbins(); ++bIX){
                           for(Int_t bIY = 0; bIY < response_Unfold_h[jI][cI][iI][mI][aI][sI]->GetYaxis()->GetNbins(); ++bIY){
                              Double_t cent = response_Unfold_h[jI][cI][iI][mI][aI][sI]->GetBinContent(bIX+1, bIY+1);
                              Double_t err = response_Unfold_h[jI][cI][iI][mI][aI][sI]->GetBinError(bIX+1, bIY+1);		    

                              Double_t newVal = 100;
                              Double_t newErr = 100;

                              if(cent < TMath::Power(10,-20)){
                                 newVal = minVal;
                                 newErr = minVal;

                                 if(sI == matrixStatDownPos){		      
                                    newVal = 0;
                                    newErr = 0;
                                 }
                                 else if(sI == matrixStatPos){
                                    newVal = TMath::Max(0., randGen_p->Gaus(0, minVal));
                                    newErr = minVal*newVal;
                                 }
                              }
                              else{		    
                                 newVal = cent + err;
                                 newErr = err*newVal/cent;
                                 if(sI == matrixStatDownPos){		      
                                    newVal = TMath::Max(0., cent - err);		     
                                    newErr = err*newVal/cent;
                                 }
                                 else if(sI == matrixStatPos){
                                    newVal = TMath::Max(0., randGen_p->Gaus(cent, err));
                                    newErr = err*newVal/cent;
                                 }
                              }		      

                              response_Unfold_h[jI][cI][iI][mI][aI][sI]->SetBinContent(bIX+1, bIY+1, newVal);
                              response_Unfold_h[jI][cI][iI][mI][aI][sI]->SetBinError(bIX+1, bIY+1, newErr);
                           }
                        }

                        for(Int_t bIX = 0; bIX < initResMC_p->GetXaxis()->GetNbins(); ++bIX){
                           for(Int_t bIY = 0; bIY < initResMC_p->GetYaxis()->GetNbins(); ++bIY){
                              Double_t cent = initResMC_p->GetBinContent(bIX+1, bIY+1);
                              Double_t err = initResMC_p->GetBinError(bIX+1, bIY+1);		    
                              Double_t newVal = 100;
                              Double_t newErr = 100;

                              if(cent < TMath::Power(10,-20)){
                                 newVal = minVal;
                                 newErr = minVal;

                                 if(sI == matrixStatDownPos){		      
                                    newVal = 0;
                                    newErr = 0;
                                 }
                                 else if(sI == matrixStatPos){
                                    newVal = TMath::Max(0., randGen_p->Gaus(0, minVal));
                                    newErr = minVal*newVal;
                                 }
                              }
                              else{
                                 newVal = cent + err;
                                 newErr = err*newVal/cent;

                                 if(sI == matrixStatDownPos){		      
                                    newVal = TMath::Max(0., cent - err);		     
                                    newErr = err*newVal/cent;
                                 }
                                 else if(sI == matrixStatPos){
                                    newVal = TMath::Max(0., randGen_p->Gaus(newVal, newErr));
                                    newErr = err*newVal/cent;
                                 }
                              }

                              initResMC_p->SetBinContent(bIX+1, bIY+1, newVal);
                              initResMC_p->SetBinError(bIX+1, bIY+1, newErr);
                           }
                        }

                        for(Int_t bIX = 0; bIX < initResMC2_p->GetXaxis()->GetNbins(); ++bIX){
                           for(Int_t bIY = 0; bIY < initResMC2_p->GetYaxis()->GetNbins(); ++bIY){
                              Double_t cent = initResMC2_p->GetBinContent(bIX+1, bIY+1);
                              Double_t err = initResMC2_p->GetBinError(bIX+1, bIY+1);		    
                              Double_t newVal = 100;
                              Double_t newErr = 100;   

                              if(cent < TMath::Power(10,-20)){
                                 newVal = minVal;
                                 newErr = minVal;

                                 if(sI == matrixStatDownPos){		      
                                    newVal = 0;
                                    newErr = 0;
                                 }
                                 else if(sI == matrixStatPos){
                                    newVal = TMath::Max(0., randGen_p->Gaus(0, minVal));
                                    newErr = minVal*newVal;
                                 }
                              }
                              else{
                                 newVal = cent + err;
                                 newErr = err*newVal/cent;
                                 if(sI == matrixStatDownPos){		      
                                    newVal = TMath::Max(0., cent - err);		     
                                    newErr = err*newVal/cent;
                                 }
                                 else if(sI == matrixStatPos){
                                    newVal = TMath::Max(0., randGen_p->Gaus(newVal, newErr));
                                    newErr = err*newVal/cent;
                                 }
                              }

                              initResMC2_p->SetBinContent(bIX+1, bIY+1, newVal);
                              initResMC2_p->SetBinError(bIX+1, bIY+1, newErr);
                           }
                        }
                     }

                     if(/*doClean*/ false){
                        cleanMatrix(response_Unfold_h[jI][cI][iI][mI][aI][sI]);
                        cleanMatrix(initResMC_p);
                        cleanMatrix(initResMC2_p);
                     }	      	      

                     std::cout << "LINE: " << __LINE__ << ", " << systStr[sI] << std::endl;	 

                     /*	      
                              TH1D* initMeas_p = (TH1D*)recoJtPt_GoodGen_h[jI][cI][iI][mI][aI][sIPos]->Clone("initMeas_p");
                              TH1D* initTrue_p = (TH1D*)genJtPt_GoodReco_h[jI][cI][iI][mI][aI][sIPos]->Clone("initTrue_p");
                              TH1D* rawClone_p = (TH1D*)jtPtRaw_RecoGenAsymm_h[jI][cI][iI][aI]->Clone("rawClone_p"); */
                     std::cout << "LINE: " << __LINE__ << ", " << systStr[sI] << std::endl;	 

                     TH1D* initMeas_p = new TH1D("initMeas_p", "", nRecoJtPtBins[jI][cI], recoJtPtBins[jI][cI]);
                     TH1D* initTrue_p = new TH1D("initTrue_p", "", nGenJtPtBins[jI][cI], genJtPtBins[jI][cI]);

                     TH1D* initMeasMC_p = new TH1D("initMeasMC_p", "", nRecoJtPtBins[jI][cI], recoJtPtBins[jI][cI]);
                     TH1D* initTrueMC_p = new TH1D("initTrueMC_p", "", nGenJtPtBins[jI][cI], genJtPtBins[jI][cI]);

                     TH1D* rawClone_p = new TH1D("rawClone_p", "", nRecoJtPtBins[jI][cI], recoJtPtBins[jI][cI]);
                     TH1D* rawCloneMC_p = new TH1D("rawCloneMC_p", "", nRecoJtPtBins[jI][cI], recoJtPtBins[jI][cI]);

                     std::cout << "LINE: " << __LINE__ << ", " << systStr[sI] << std::endl;	 

                     macroHistToSubsetHist(jtPtRaw_General_h[jI][cI][iI][aI], rawClone_p, false);	      
                     macroHistToSubsetHistX(initResMC_p, rawCloneMC_p, true);	      
                     macroHistToSubsetHistY(initResMC_p, jtPtUnfoldedMCTruth_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI], true);	      

                     if(sI == 0){
                        macroHistToSubsetHist(jtPtRaw_General_h[jI][cI][iI][aI], jtPtRaw_RecoGenAsymm_h[jI][cI][iI][aI], false);
                     }

                     // Using a clone so we can modify for Fake err;

                     setPrior(response_Unfold_h[jI][cI][iI][mI][aI][sI], initResMC2_p);
                     //	      do2By2Reweight(rawClone_p, response_Unfold_h[jI][cI][iI][mI][aI][sI]);
                     macroHistToSubsetHistY(response_Unfold_h[jI][cI][iI][mI][aI][sI], initTrue_p, true);
                     macroHistToSubsetHistX(response_Unfold_h[jI][cI][iI][mI][aI][sI], initMeas_p, true);			    

                     setPrior(initResMC_p, initResMC2_p);
                     //	      do2By2Reweight(rawCloneMC_p, initResMC_p);
                     macroHistToSubsetHistY(initResMC_p, initTrueMC_p, true);
                     macroHistToSubsetHistX(initResMC_p, initMeasMC_p, true);		

                     //NOTE: ERRORS HERE ARE HARD CODED RELATIVE VALUES BASED ON PLOTS IN AN

                     if(isStrSame(systStr[sI], "Fake") && !isDataPP){
                        if(tempStr.find("akCs10PU3PFFlow") != std::string::npos){
                           if(centBinsLow[cI] == 0 && centBinsHi[cI] == 10) correctForFakes(rawClone_p, {100, 150, 200, 250}, {0.425, 0.3, 0.21});
                           else if(centBinsLow[cI] == 10 && centBinsHi[cI] == 30) correctForFakes(rawClone_p, {100, 150, 200}, {0.08, 0.075});
                        }
                        else if(tempStr.find("akCs8PU3PFFlow") != std::string::npos){
                           if(centBinsLow[cI] == 0 && centBinsHi[cI] == 10) correctForFakes(rawClone_p, {100, 150, 200, 250}, {0.31, 0.16, 0.04});
                           else if(centBinsLow[cI] == 10 && centBinsHi[cI] == 30) correctForFakes(rawClone_p, {100, 150}, {0.04}); 
                        }
                        else if(tempStr.find("akCs6PU3PFFlow") != std::string::npos){
                           if(centBinsLow[cI] == 0 && centBinsHi[cI] == 10) correctForFakes(rawClone_p, {100, 150, 200}, {0.07, 0.07}); //taking max since the two bins are not descending and statistically compatible
                        }
                        else if(tempStr.find("akCs4PU3PFFlow") != std::string::npos){
                           if(centBinsLow[cI] == 0 && centBinsHi[cI] == 10) correctForFakes(rawClone_p, {100, 150}, {0.07});
                        }
                     }	      

                     const Int_t nBinsMax = 100;
                     const Int_t nBinsTemp = initMeas_p->GetNbinsX();
                     Double_t binsTemp[nBinsMax];

                     for(Int_t bI = 0; bI < nBinsTemp+1; ++bI){
                        binsTemp[bI] = initMeas_p->GetBinLowEdge(bI+1);
                     }

                     TH1D* tempHist_p = new TH1D("tempHist_p", "", nBinsTemp, binsTemp);
                     tempHist_p->Divide(rawClone_p, initMeas_p);

                     for(Int_t bI = 0; bI < nBayes; ++bI){	    				
                        RooUnfoldResponse* rooResClone_p = NULL;
                        RooUnfoldResponse* rooResCloneMC_p = NULL;

                        TH1D* unfold_h = new TH1D("unfold_h", "", nGenJtPtBins[jI][cI], genJtPtBins[jI][cI]);
                        TH1D* unfoldMC_h = new TH1D("unfoldMC_h", "", nGenJtPtBins[jI][cI], genJtPtBins[jI][cI]);
                        bayes_p[jI][cI][iI][mI][aI][sI][bI] = NULL;
                        bayesMC_p[jI][cI][iI][mI][aI][sI][bI] = NULL;		
                        doUnfoldWatch.start();		  		

                        if(sI == 0){
                           const std::string saveName = "pdfDir/" + dateStr + "/Unfold_" + tempStr + "/responseMat_" + tempStr + "_" + centStr + "_" + idStr[iI] + "_" + resStr + "_" + jtAbsEtaStr + "_" + tempSystStr +  "AllBayes_RecoGenAsymm_" + cleanStr + "_" + toyStr + "_" + debugStr + addedTagStr + dateStr + ".pdf";
                           drawResponse(response_Unfold_h[jI][cI][iI][mI][aI][sI], saveName);
                        }		

                        std::cout << "LINE: " << __LINE__ << std::endl;

                        doUnfold(response_Unfold_h[jI][cI][iI][mI][aI][sI], initTrue_p, rawClone_p, bayesVal[bI], &rooResClone_p, &(bayes_p[jI][cI][iI][mI][aI][sI][bI]), unfold_h, nToys, false);		
                        std::cout << "Final CHI2: " << bayes_p[jI][cI][iI][mI][aI][sI][bI]->getChi2() << std::endl;

                        std::cout << "LINE: " << __LINE__ << std::endl;

                        doUnfold(initResMC_p, initTrueMC_p, rawCloneMC_p, bayesVal[bI], &rooResCloneMC_p, &(bayesMC_p[jI][cI][iI][mI][aI][sI][bI]), unfoldMC_h, nToys, false);

                        ++unfoldCounter;


                        std::cout << "LINE: " << __LINE__ << std::endl;

                        if(bayesVal[bI] == 4 && sI == 0){
                           std::cout << "COVARIANCE" << std::endl;
                           TH1D* unfoldTemp_h = new TH1D("unfoldTemp_h", "", nGenJtPtBins[jI][cI], genJtPtBins[jI][cI]);
                           std::cout << "LINE: " << __LINE__ << std::endl;

                           for(Int_t covI = 0; covI < nCov; ++covI){
                              RooUnfoldResponse* rooResCloneTemp_p = NULL;
                              RooUnfoldBayes* bayesTemp_p = NULL;
                              doUnfold(response_Unfold_h[jI][cI][iI][mI][aI][sI], initTrue_p, rawClone_p, bayesVal[bI], &rooResCloneTemp_p, &bayesTemp_p, unfoldTemp_h, nToys, false);

                              TH2D* pearsonToPlot = NULL;
                              getPearsTMatrix(bayesTemp_p, &covariance_h[jI][cI][iI][mI][aI][covI], &pearsonToPlot);

                              delete pearsonToPlot;
                              delete bayesTemp_p;
                              delete rooResCloneTemp_p;
                           }
                           std::cout << "LINE: " << __LINE__ << std::endl;

                           delete unfoldTemp_h;
                           std::cout << "END COVARIANCE" << std::endl;
                        }

                        std::cout << "LINE: " << __LINE__ << std::endl;

                        doUnfoldWatch.stop();

                        if(bI == 0){
                           std::string reweightName = response_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI]->GetName();
                           reweightName.replace(reweightName.find("RecoGenAsymm"), std::string("RecoGenAsymm").size(), "Reweighted");
                           rooResponse_Reweighted_h[jI][cI][iI][mI][aI][sI] = new RooUnfoldResponse(initMeas_p, initTrue_p, response_Unfold_h[jI][cI][iI][mI][aI][sI], reweightName.c_str());
                        }

                        std::cout << "LINE: " << __LINE__ << std::endl;

                        for(Int_t bIX = 0; bIX < unfold_h->GetNbinsX(); ++bIX){
                           jtPtUnfoldedBayes_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][bI]->SetBinContent(bIX+1, unfold_h->GetBinContent(bIX+1));
                           jtPtUnfoldedBayes_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][bI]->SetBinError(bIX+1, unfold_h->GetBinError(bIX+1));
                        }

                        std::cout << "LINE: " << __LINE__ << std::endl;

                        for(Int_t bIX = 0; bIX < unfold_h->GetNbinsX(); ++bIX){
                           jtPtUnfoldedBayesMC_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][bI]->SetBinContent(bIX+1, unfoldMC_h->GetBinContent(bIX+1));
                           jtPtUnfoldedBayesMC_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][bI]->SetBinError(bIX+1, unfoldMC_h->GetBinError(bIX+1));
                        }

                        delete rooResClone_p;
                        delete rooResCloneMC_p;
                        delete unfold_h;
                        delete unfoldMC_h;
                     }

                     std::cout << "LINE: " << __LINE__ << std::endl;

                     for(Int_t bI = 0; bI < TMath::Min(nGenJtPtBins[jI][cI], nSvd); ++bI){	    				
                        RooUnfoldResponse* rooResClone_p = NULL;
                        RooUnfoldResponse* rooResCloneMC_p = NULL;

                        TH1D* unfold_h = new TH1D("unfold_h", "", nGenJtPtBins[jI][cI], genJtPtBins[jI][cI]);
                        TH1D* unfoldMC_h = new TH1D("unfoldMC_h", "", nGenJtPtBins[jI][cI], genJtPtBins[jI][cI]);
                        svd_p[jI][cI][iI][mI][aI][sI][bI] = NULL;
                        svdMC_p[jI][cI][iI][mI][aI][sI][bI] = NULL;		
                        doUnfoldWatch.start();		  		

                        doUnfoldSvd(response_Unfold_h[jI][cI][iI][mI][aI][sI], initTrue_p, rawClone_p, bI+1, &rooResClone_p, &(svd_p[jI][cI][iI][mI][aI][sI][bI]), unfold_h, nToys);
                        doUnfoldSvd(initResMC_p, initTrueMC_p, rawCloneMC_p, bI+1, &rooResCloneMC_p, &(svdMC_p[jI][cI][iI][mI][aI][sI][bI]), unfoldMC_h, nToys);


                        for(Int_t bIX = 0; bIX < unfold_h->GetNbinsX(); ++bIX){
                           jtPtUnfoldedSvd_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][bI]->SetBinContent(bIX+1, unfold_h->GetBinContent(bIX+1));
                           jtPtUnfoldedSvd_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][bI]->SetBinError(bIX+1, unfold_h->GetBinError(bIX+1));
                        }

                        for(Int_t bIX = 0; bIX < unfold_h->GetNbinsX(); ++bIX){
                           jtPtUnfoldedSvdMC_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][bI]->SetBinContent(bIX+1, unfoldMC_h->GetBinContent(bIX+1));
                           jtPtUnfoldedSvdMC_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][bI]->SetBinError(bIX+1, unfoldMC_h->GetBinError(bIX+1));
                        }

                        delete rooResClone_p;
                        delete rooResCloneMC_p;
                        delete unfold_h;
                        delete unfoldMC_h;
                     }

                     std::cout << "LINE: " << __LINE__ << std::endl;

                     delete initMeas_p;
                     delete initTrue_p;

                     std::cout << "LINE: " << __LINE__ << std::endl;

                     delete initResMC_p;
                     delete initResMC2_p;
                     delete initMeasMC_p;
                     delete initTrueMC_p;

                     std::cout << "LINE: " << __LINE__ << std::endl;

                     delete rawClone_p;	      
                     delete rawCloneMC_p;	      
                     std::cout << "LINE: " << __LINE__ << std::endl;
                  }
                  std::cout << "LINE: " << __LINE__ << std::endl;
               }
               std::cout << "LINE: " << __LINE__ << std::endl;
            }
            std::cout << "LINE: " << __LINE__ << std::endl;
         }
         std::cout << "LINE: " << __LINE__ << std::endl;
      }
      std::cout << "LINE: " << __LINE__ << std::endl;
   }

   std::cout << "LINE: " << __LINE__ << std::endl;

   std::cout << "Unfolding complete." << std::endl;

   dataFile_p->Close();
   delete dataFile_p;

   responseFile_p->Close();
   delete responseFile_p;

   outFile_p->cd();

   const Double_t yPadFrac = 0.45;
   const Double_t marg = 0.12;

   const Int_t nStyles = 6;
   const Int_t styles[nStyles] = {24, 25, 27, 28, 46, 44};

   const Int_t nColors = 5;
   const Int_t colors[nColors] = {1, vg.getColor(0), vg.getColor(1), vg.getColor(2), vg.getColor(4)};

   std::vector<std::vector<std::string> > pdfNames;

   Int_t nPDFTotal = 0;

   for(Int_t jI = 0; jI < nDataJet; ++jI){
      outFile_p->cd();
      dir_p[jI]->cd();

      std::string tempStr = responseJetDirList[jI];
      if(tempStr.find("/") != std::string::npos) tempStr.replace(tempStr.find("/"), tempStr.size() - tempStr.find("/"), "");


      Double_t lowPtTruncVal = 200.;
      Bool_t isBigJet = tempStr.find("ak8") != std::string::npos || tempStr.find("ak10") != std::string::npos || tempStr.find("akCs8") != std::string::npos || tempStr.find("akCs10") != std::string::npos;
      if(isBigJet) lowPtTruncVal = 300.;

      for(Int_t cI = 0; cI < nCentBins; ++cI){
         std::string centStr = "Cent" + std::to_string(centBinsLow[cI]) + "to" + std::to_string(centBinsHi[cI]);
         std::string centStr2 = "PP";
         if(!isDataPP){
            centStr = "PbPb_" + centStr;
            centStr2 = std::to_string(centBinsLow[cI]) + "-" + std::to_string(centBinsHi[cI]) + "%";
         }
         else centStr = "PP_" + centStr;


         for(Int_t iI = 0; iI < nID; ++iI){
            for(Int_t mI = 0; mI < nResponseMod; ++mI){	
               const std::string resStr = "ResponseMod" + prettyString(responseMod[mI], 2, true);

               for(Int_t aI = 0; aI < nJtAbsEtaBins; ++aI){
                  const std::string jtAbsEtaStr = "AbsEta" + prettyString(jtAbsEtaBinsLow[aI], 1, true) + "to" + prettyString(jtAbsEtaBinsHi[aI], 1, true);
                  //	    if(jtAbsEtaStr.find("AbsEta0p0to2p0") == std::string::npos) continue;
                  const std::string jtAbsEtaStr2 = prettyString(jtAbsEtaBinsLow[aI], 1, false) + "<|#eta|<" +  prettyString(jtAbsEtaBinsHi[aI], 1, false);


                  pdfNames.push_back({});

                  for(Int_t sI = 0; sI < nSyst; ++sI){	    	   
                     //	      if(systStr[sI].find("MatrixStat") != std::string::npos) continue;

                     const std::string tempSystStr = systStr[sI] + "_";

                     TLegend* leg_p = new TLegend(0.15, 0.05, 0.5, 0.6);
                     leg_p->SetBorderSize(0.0);
                     leg_p->SetFillStyle(0);
                     leg_p->SetFillColor(0);
                     leg_p->SetTextFont(43);
                     leg_p->SetTextSize(14);

                     TLatex* label_p = new TLatex();
                     label_p->SetTextFont(43);
                     label_p->SetTextSize(14);
                     label_p->SetNDC();                  

                     const Int_t nPads = 12;
                     const Double_t nPadsX = 4;
                     const Double_t nPadsY = 2;

                     Double_t padWidth = 1./nPadsX;
                     Double_t padHeight = 1./nPadsY;

                     const Double_t padsXLow[nPads] = {0*padWidth, 
                        0*padWidth, 
                        0*padWidth, 
                        1*padWidth, 
                        1*padWidth, 
                        2*padWidth, 
                        3*padWidth, 
                        0*padWidth, 
                        1*padWidth, 
                        2*padWidth,
                        3*padWidth,
                        2*padWidth};
                     std::cout << "LINE: " << __LINE__ << std::endl;

                     const Double_t padsXHi[nPads] = {1*padWidth, 
                        1*padWidth, 
                        1*padWidth, 
                        2*padWidth, 
                        2*padWidth, 
                        3*padWidth, 
                        4*padWidth, 
                        1*padWidth, 
                        2*padWidth, 
                        3*padWidth, 
                        4*padWidth,
                        3*padWidth};

                     std::cout << "LINE: " << __LINE__ << std::endl;
                     const Double_t padsYLow[nPads] = {0.5 + 0.5*yPadFrac,
                        0.5 + 0.5*yPadFrac - (0.5*yPadFrac - marg*nPadsY/nPadsX)/2, 
                        1*padHeight, 
                        0.5 + 0.5*yPadFrac - (0.5*yPadFrac - marg*nPadsY/nPadsX)/2, 
                        1*padHeight,
                        1*padHeight, 
                        1*padHeight, 
                        0*padHeight, 
                        0*padHeight, 
                        0*padHeight + 0.5*yPadFrac, 
                        0*padHeight,
                        0*padHeight};

                     std::cout << "LINE: " << __LINE__ << std::endl;
                     const Double_t padsYHi[nPads] = {2*padHeight,
                        0.5 + 0.5*yPadFrac, 
                        0.5 + 0.5*yPadFrac - (0.5*yPadFrac - marg*nPadsY/nPadsX)/2, 
                        2*padHeight, 
                        0.5 + 0.5*yPadFrac - (0.5*yPadFrac - marg*nPadsY/nPadsX)/2,
                        2*padHeight,
                        2*padHeight, 
                        1*padHeight, 
                        1*padHeight,
                        1*padHeight, 
                        1*padHeight,
                        0*padHeight + 0.5*yPadFrac};

                     const Double_t bottomMarg[nPads] = {0.001, 0.001, marg*(nPadsY/nPadsX)/(padsYHi[2] - padsYLow[2]), 0.001, marg*(nPadsY/nPadsX)/(padsYHi[4] - padsYLow[4]), marg*(nPadsY/nPadsX)/(padsYHi[5] - padsYLow[5]), marg*(nPadsY/nPadsX)/(padsYHi[6] - padsYLow[6]), marg*(nPadsY/nPadsX)/(padsYHi[7] - padsYLow[7]), marg*(nPadsY/nPadsX)/(padsYHi[8] - padsYLow[8]), 0.001, marg*(nPadsY/nPadsX)/(padsYHi[10] - padsYLow[10]), marg*(nPadsY/nPadsX)/(padsYHi[11] - padsYLow[11])};

                     TCanvas* canv_p = new TCanvas("canv_p", "", 450*4, 450*2);
                     TPad* pads[nPads];
                     for(Int_t padI = 0; padI < nPads; ++padI){
                        canv_p->cd();
                        pads[padI] = new TPad(("pad" + std::to_string(padI)).c_str(), "", padsXLow[padI], padsYLow[padI], padsXHi[padI], padsYHi[padI]);
                        pads[padI]->SetLeftMargin(marg);

                        Double_t topMargin = 0.0001;
                        if(padI != 1 && padI != 2 && padI != 4 && padI != 11){
                           topMargin = 0.03/(padsYHi[padI] - padsYLow[padI]);
                        }
                        //		topMargin = 0.0001;
                        pads[padI]->SetTopMargin(topMargin);
                        pads[padI]->SetBottomMargin(bottomMarg[padI]);
                        pads[padI]->SetRightMargin(0.002);
                        pads[padI]->Draw();
                     }

                     //start here for terminalpos purposes
                     canv_p->cd();
                     pads[2]->cd();

                     int terminalPos = -1;
                     int terminalPos5 = -1;
                     int terminalPos5AndPears = -1;
                     int terminalChi2 = -1;
                     //	      int terminalChi2AndPears = -1;
                     double terminalSumMin = 999999;
                     double terminalSumMinChi2 = 999999;
                     std::cout << "LINE: " << __LINE__ << std::endl;

                     TH1D* bandValLow_p = (TH1D*)jtPtUnfoldedBayes_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][bayes100Pos]->Clone("bandValLow");
                     TH1D* bandValHi_p = (TH1D*)jtPtUnfoldedBayes_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][bayes100Pos]->Clone("bandValHi");
                     std::cout << "LINE: " << __LINE__ << std::endl;

                     for(Int_t bI = bayes100Pos-nBigBayesSymm; bI <= bayes100Pos+nBigBayesSymm; ++bI){
                        for(Int_t bIX = 0; bIX < bandValLow_p->GetNbinsX(); ++bIX){
                           if(jtPtUnfoldedBayes_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][bI]->GetBinContent(bIX+1) < bandValLow_p->GetBinContent(bIX+1)){
                              bandValLow_p->SetBinContent(bIX+1, jtPtUnfoldedBayes_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][bI]->GetBinContent(bIX+1));
                           }
                           std::cout << "LINE: " << __LINE__ << std::endl;
                           if(jtPtUnfoldedBayes_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][bI]->GetBinContent(bIX+1) > bandValHi_p->GetBinContent(bIX+1)){
                              bandValHi_p->SetBinContent(bIX+1, jtPtUnfoldedBayes_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][bI]->GetBinContent(bIX+1));
                           }
                        }
                     }

                     std::cout << "LINE: " << __LINE__ << std::endl;
                     TH1D* clones_p[nBayes];

                     Double_t max = -1;
                     Double_t min = 100;
                     for(Int_t bI = 0; bI < nBayes; ++bI){
                        std::cout << "LINE: " << __LINE__ << std::endl;
                        clones_p[bI] = (TH1D*)jtPtUnfoldedBayes_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][bI]->Clone(("clone_" + std::to_string(bI)).c_str());

                        bool sub1Perc = true;
                        bool sub5Perc = true;

                        for(Int_t bIX = 0; bIX < clones_p[bI]->GetNbinsX(); ++bIX){
                           double binCenter = (clones_p[bI]->GetBinLowEdge(bIX+1) + clones_p[bI]->GetBinLowEdge(bIX+2))/2.;
                           double binContent = clones_p[bI]->GetBinContent(bIX+1);
                           double bandLowContent = bandValLow_p->GetBinContent(bIX+1);
                           double bandHiContent = bandValHi_p->GetBinContent(bIX+1);

                           if(binCenter > lowPtTruncVal && binCenter < 1000.){
                              if(sub1Perc){
                                 if(binContent/bandHiContent > 1.01) sub1Perc = false;
                                 else if(binContent/bandLowContent < .99) sub1Perc = false;
                              }
                              if(sub5Perc){
                                 if(binContent/bandHiContent > 1.05) sub5Perc = false;
                                 else if(binContent/bandLowContent < .95) sub5Perc = false;
                              }
                           }
                        }
                        std::cout << "LINE: " << __LINE__ << std::endl;

                        if(sub1Perc && terminalPos < 0 && bI < nBayesBig){
                           terminalPos = bI;
                           //		    histTermPos[jI][cI][iI][mI][aI][sI] = terminalPos;
                        }
                        std::cout << "LINE: " << __LINE__ << std::endl;
                        if(sub5Perc && terminalPos5 < 0 && bI < nBayesBig) terminalPos5 = bI;
                        std::cout << "LINE: " << __LINE__ << std::endl;
                        if(sub5Perc && bI < nBayesBig){
                           std::cout << "LINE: " << __LINE__ << std::endl;
                           TH2D* covarianceToPlot = NULL;
                           TH2D* pearsonToPlot = NULL;
                           std::cout << "LINE: " << __LINE__ << ", " << sI << ", " << bI << std::endl;
                           getPearsTMatrix(bayes_p[jI][cI][iI][mI][aI][sI][bI], &covarianceToPlot, &pearsonToPlot);

                           std::cout << "LINE: " << __LINE__ << std::endl;
                           bool keepMinOppDiag = true;
                           Double_t sum = 0;
                           for(Int_t bIX = 0; bIX < pearsonToPlot->GetNbinsX(); ++bIX){
                              for(Int_t bIY = 0; bIY < pearsonToPlot->GetNbinsY(); ++bIY){
                                 if(bIX == bIY) continue;//continue on diagonal

                                 Double_t val = pearsonToPlot->GetBinContent(bIX+1, bIY+1);
                                 if(bIX == 0 && bIY == pearsonToPlot->GetNbinsY()-1 && TMath::Abs(val) > 0.15) keepMinOppDiag = false;

                                 sum += TMath::Abs(val);
                              }
                           }

                           if(sum < terminalSumMin && keepMinOppDiag){
                              terminalSumMin = sum;
                              terminalPos5AndPears = bI;
                              //		    histTermPos[jI][cI][iI][mI][aI][sI] = bI;//terminalPos5AndPears;
                           }

                           delete pearsonToPlot;
                        }
                        std::cout << "LINE: " << __LINE__ << std::endl;

                        double pVal = -1;
                        if(bI >= 1){
                           double chi2 = bayes_p[jI][cI][iI][mI][aI][sI][bI]->getChi2();		  
                           int nPDF = nGenJtPtBins[jI][cI]-1;		  

                           //		  pVal = TMath::Prob(chi2, jtPtUnfoldedBayes_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][bI]->GetXaxis()->GetNbins()-1);
                           //		  TH1D* histNom_p = (TH1D*)jtPtUnfoldedBayes_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][bI]->Clone("histNom_h");
                           //		  TH1D* histPrev_p = (TH1D*)jtPtUnfoldedBayes_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][bI-1]->Clone("histPrev_h");

                           TH1D* histNom_p = new TH1D("histNom_p", "", nGenJtPtBins[jI][cI], genJtPtBins[jI][cI]);
                           TH1D* histPrev_p = new TH1D("histPrev_p", "", nGenJtPtBins[jI][cI], genJtPtBins[jI][cI]);
                           macroHistToSubsetHist(jtPtUnfoldedBayes_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][bI], histNom_p, true);
                           macroHistToSubsetHist(jtPtUnfoldedBayes_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][bI-1], histPrev_p, true);

                           bool ratioSub1 = true;

                           std::cout << "LINE: " << __LINE__ << std::endl;

                           for(Int_t bIX = 1; bIX < histNom_p->GetNbinsX()-1; ++bIX){
                              Double_t val = ((Double_t)histNom_p->GetBinContent(bIX+1))/((Double_t)histPrev_p->GetBinContent(bIX+1));

                              if(histNom_p->GetBinLowEdge(bIX+1) < recoJtPtBins[jI][cI][0]) continue;
                              if(histNom_p->GetBinLowEdge(bIX+2) > recoJtPtBins[jI][cI][nRecoJtPtBins[jI][cI]]) continue;
                              //		    std::cout << "CHECK RATIO: " << val << ", " << bI << std::endl;

                              if(val < 1.01 && val > 0.99) continue;

                              ratioSub1 = false;
                              break;
                           }


                           std::cout << "LINE: " << __LINE__ << std::endl;

                           if(ratioSub1){
                              std::cout << "RATIO IS SUB 1 (BAYES=" << bI << "): " << std::endl;
                              ratioPrint(histNom_p, histPrev_p);
                              histNom_p->Print("ALL");
                              histPrev_p->Print("ALL");
                           }
                           else{
                              std::cout << "RATIO IS NOT SUB 1 (BAYES=" << bI << "): " << std::endl;
                              ratioPrint(histNom_p, histPrev_p);
                              histNom_p->Print("ALL");
                              histPrev_p->Print("ALL");
                           }


                           for(Int_t bIX = 0; bIX < histNom_p->GetXaxis()->GetNbins(); ++bIX){
                              Double_t count1 = histNom_p->GetBinContent(bIX+1);
                              Double_t count2 = histPrev_p->GetBinContent(bIX+1);

                              Double_t diff = TMath::Sqrt(TMath::Abs(count1 - count2));

                              histNom_p->SetBinError(bIX+1, diff);
                              histPrev_p->SetBinError(bIX+1, diff);
                              /*
                                 if(err1 > err2){
                                 Double_t newErr = err1 - err2;

                              //		      Double_t newErr = TMath::Sqrt(err1*err1 - err2*err2);
                              histNom_p->SetBinError(bIX+1, newErr);
                              histPrev_p->SetBinError(bIX+1, newErr*err2/err1);
                              }
                              else{
                              Double_t newErr = err2 - err1;
                              //		      Double_t newErr = TMath::Sqrt(err2*err2 - err1*err1);
                              histNom_p->SetBinError(bIX+1, newErr);
                              histPrev_p->SetBinError(bIX+1, newErr*err1/err2);
                              }
                              */
                           }

                           pVal = histNom_p->Chi2Test(histPrev_p, "WW");
                           pVal = TMath::Prob(chi2, nPDF);

                           /*		  
                                   if(ratioSub1){
                                   pVal = .9999;		  
                           //		    std::cout << "RATIO IS SUB 1" << std::endl;
                           }
                           else pVal = 0.0;
                           */

                           delete histNom_p;
                           delete histPrev_p;
                        }

                        if(bI >= 1 && terminalChi2 < 0 && pVal >= 0.95){
                           terminalChi2 = bI;
                           histTermPos[jI][cI][iI][mI][aI][sI] = bI;
                           std::cout << "TERMINATION AT: " << bI << ", " << bayesVal[bI] << std::endl;
                        }

                        if(bI >= 1 && pVal >= 0.95){
                           TH2D* covarianceToPlot = NULL;
                           TH2D* pearsonToPlot = NULL;
                           getPearsTMatrix(bayes_p[jI][cI][iI][mI][aI][sI][bI], &covarianceToPlot, &pearsonToPlot);

                           std::cout << "LINE: " << __LINE__ << std::endl;
                           bool keepMinOppDiag = true;
                           Double_t sum = 0;
                           for(Int_t bIX = 0; bIX < pearsonToPlot->GetNbinsX(); ++bIX){
                              for(Int_t bIY = 0; bIY < pearsonToPlot->GetNbinsY(); ++bIY){
                                 if(bIX == bIY) continue;//continue on diagonal

                                 Double_t val = pearsonToPlot->GetBinContent(bIX+1, bIY+1);
                                 if(bIX == 0 && bIY == pearsonToPlot->GetNbinsY()-1 && TMath::Abs(val) > 0.15) keepMinOppDiag = false;

                                 sum += TMath::Abs(val);
                              }
                           }

                           if(sum < terminalSumMinChi2 && keepMinOppDiag){
                              terminalSumMinChi2 = sum;
                              //		    terminalChi2AndPears = bI;

                              //		    histTermPos[jI][cI][iI][mI][aI][sI] = bI;
                              //		    histTermPos[jI][cI][iI][mI][aI][sI] = bI;//terminalPos5AndPears;
                           }

                           delete pearsonToPlot;
                        }

                        std::cout << "LINE: " << __LINE__ << std::endl;

                        clones_p[bI]->Divide(jtPtUnfoldedBayes_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][bayes100Pos]);

                        std::cout << "LINE: " << __LINE__ << std::endl;
                        std::cout << "PRINT CHECK 1: " << std::endl;
                        jtPtUnfoldedBayes_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][bayes100Pos]->Print("ALL");
                        std::cout << "PRINT CHECK 2: " << std::endl;
                        jtPtUnfoldedBayes_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][bI]->Print("ALL");

                        std::cout << "LINE: " << __LINE__ << std::endl;

                        Double_t tempMax = getMax(clones_p[bI]);
                        Double_t tempMin = getMinGTZero(clones_p[bI]);
                        std::cout << "LINE: " << __LINE__ << std::endl;
                        if(tempMax > max) max = tempMax;
                        if(tempMin < min) min = tempMin;
                        std::cout << "LINE: " << __LINE__ << std::endl;
                     }
                     std::cout << "LINE: " << __LINE__ << std::endl;

                     //	      histTermPos[jI][cI][iI][mI][aI][sI] = 3;

                     delete bandValLow_p;
                     delete bandValHi_p;

                     Double_t interval = (max - min)/10.;
                     max += interval;
                     min -= interval;

                     const Int_t nTempBins = 10;
                     Double_t tempBins[nTempBins+1];
                     getLinBins(min, max, nTempBins, tempBins);

                     if(max > 1.15){
                        max = 1.15;
                        interval = (max - min)/10.;
                     }
                     if(min < .85){
                        min = 0.85;
                        interval = (max - min)/10.;
                     }

                     max = 1.15;
                     min = 0.85;
                     std::cout << "LINE: " << __LINE__ << std::endl;

                     for(Int_t bI = 1; bI < nBayes; ++bI){		
                        clones_p[bI]->SetMaximum(max);
                        clones_p[bI]->SetMinimum(min);

                        clones_p[bI]->GetYaxis()->SetTitle("Ratio w/ Bayes 100");
                        clones_p[bI]->GetYaxis()->SetTitleSize(9);
                        clones_p[bI]->GetYaxis()->SetLabelSize(9);
                        clones_p[bI]->GetYaxis()->SetTitleOffset(clones_p[bI]->GetYaxis()->GetTitleOffset()*1.5);
                        clones_p[bI]->GetYaxis()->SetNdivisions(505);

                        if(bI < nBayesDraw){
                           Int_t bayesPos = bI;
                           if(bI == nBayesDraw-1 && histTermPos[jI][cI][iI][mI][aI][sI] >= nBayesDraw) bayesPos = histTermPos[jI][cI][iI][mI][aI][sI];

                           clones_p[bayesPos]->SetMarkerStyle(styles[bayesPos%nStyles]);
                           clones_p[bayesPos]->SetMarkerColor(colors[bayesPos%nColors]);
                           clones_p[bayesPos]->SetLineColor(colors[bayesPos%nColors]);
                           clones_p[bayesPos]->SetMarkerSize(1.0);

                           clones_p[bayesPos]->GetXaxis()->SetTitleFont(43);
                           clones_p[bayesPos]->GetYaxis()->SetTitleFont(43);
                           clones_p[bayesPos]->GetXaxis()->SetLabelFont(43);
                           clones_p[bayesPos]->GetYaxis()->SetLabelFont(43);

                           clones_p[bayesPos]->GetXaxis()->SetTitleSize(14);
                           clones_p[bayesPos]->GetYaxis()->SetTitleSize(14);
                           clones_p[bayesPos]->GetXaxis()->SetLabelSize(14);
                           clones_p[bayesPos]->GetYaxis()->SetLabelSize(14);

                           clones_p[bayesPos]->GetXaxis()->SetTitleOffset(5.);
                           clones_p[bayesPos]->GetYaxis()->SetTitleOffset(1.5);



                           if(bayesPos == 1) clones_p[bayesPos]->DrawCopy("HIST E1 P");
                           else clones_p[bayesPos]->DrawCopy("HIST E1 P SAME");
                        }
                     }

                     bool doLogX = false;
                     if(jtPtUnfoldedBayes_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][0]->GetBinWidth(1)*3 < jtPtUnfoldedBayes_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][0]->GetBinWidth(jtPtUnfoldedBayes_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][0]->GetNbinsX()-1)) doLogX = true;

                     gStyle->SetOptStat(0);
                     if(doLogX) gPad->SetLogx();

                     canvNDCToXY labelAid(pads[2], clones_p[0]);
                     for(unsigned int xI = 0; xI < nXVals; ++xI){
                        if(xVals[xI] < clones_p[0]->GetBinLowEdge(1)) continue;
                        if(xVals[xI] >= clones_p[0]->GetBinLowEdge(clones_p[0]->GetNbinsX()+1)) continue;

                        canv_p->cd();
                        pads[2]->cd();
                        label_p->DrawLatex(labelAid.getXRelFromAbs(xVals[xI], pads[2]->GetLogx()), pads[2]->GetBottomMargin(), std::to_string(xVals[xI]).c_str());
                     }		  

                     for(Int_t bI = 0; bI < nBayes; ++bI){
                        delete clones_p[bI];
                        clones_p[bI] = NULL;
                     }

                     drawWhiteBox(900, 1100, .00, min*.999);

                     canv_p->cd();	       
                     pads[0]->cd();

                     min = 1000000000;
                     max = -1;

                     for(Int_t bI = 0; bI < nBayesDraw; ++bI){
                        Int_t bayesPos = bI;
                        if(bI == nBayesDraw-1 && histTermPos[jI][cI][iI][mI][aI][sI] >= nBayesDraw) bayesPos = histTermPos[jI][cI][iI][mI][aI][sI];

                        Double_t tempMax = getMax(jtPtUnfoldedBayes_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][bayesPos]);
                        Double_t tempMin = getMinGTZero(jtPtUnfoldedBayes_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][bayesPos]);

                        if(tempMax > max) max = tempMax;
                        if(tempMin < min) min = tempMin;
                     }

                     if(min < TMath::Power(10,-8)) min = TMath::Power(10,-8);

                     Double_t globalMin = getNearestFactor10Down(min, 3);
                     Double_t globalMax = getNearestFactor10Up(max, 2);

                     globalMin = TMath::Power(10, -7);
                     globalMax = TMath::Power(10, 6);

                     for(Int_t bI = 0; bI < nBayes; ++bI){
                        jtPtUnfoldedBayes_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][bI]->SetMaximum(globalMax);
                        jtPtUnfoldedBayes_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][bI]->SetMinimum(globalMin);
                        std::cout << "LINE: " << __LINE__ << std::endl;


                        jtPtUnfoldedBayes_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][bI]->GetXaxis()->SetTitleFont(43);
                        jtPtUnfoldedBayes_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][bI]->GetYaxis()->SetTitleFont(43);
                        jtPtUnfoldedBayes_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][bI]->GetXaxis()->SetLabelFont(43);
                        jtPtUnfoldedBayes_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][bI]->GetYaxis()->SetLabelFont(43);

                        jtPtUnfoldedBayes_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][bI]->GetXaxis()->SetTitleSize(14);
                        jtPtUnfoldedBayes_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][bI]->GetYaxis()->SetTitleSize(14);
                        jtPtUnfoldedBayes_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][bI]->GetXaxis()->SetLabelSize(14);
                        jtPtUnfoldedBayes_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][bI]->GetYaxis()->SetLabelSize(14);

                        jtPtUnfoldedBayes_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][bI]->GetXaxis()->SetTitleOffset(5.);
                        jtPtUnfoldedBayes_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][bI]->GetYaxis()->SetTitleOffset(1.5);

                        if(bI < nBayesDraw){
                           Int_t bayesPos = bI;
                           if(bI == nBayesDraw-1 && histTermPos[jI][cI][iI][mI][aI][sI] >= nBayesDraw) bayesPos = histTermPos[jI][cI][iI][mI][aI][sI];

                           //		  jtPtUnfoldedBayes_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][bayesPos]->SetMaximum();
                           //		  jtPtUnfoldedBayes_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][bayesPos]->SetMinimum();

                           jtPtUnfoldedBayes_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][bayesPos]->SetMarkerStyle(styles[bayesPos%nStyles]);
                           jtPtUnfoldedBayes_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][bayesPos]->SetMarkerColor(colors[bayesPos%nColors]);
                           jtPtUnfoldedBayes_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][bayesPos]->SetLineColor(colors[bayesPos%nColors]);
                           jtPtUnfoldedBayes_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][bayesPos]->SetMarkerSize(1.0);

                           jtPtUnfoldedBayesMC_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][bI]->SetMarkerStyle(styles[bI%nStyles]);
                           jtPtUnfoldedBayesMC_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][bI]->SetMarkerColor(colors[bI%nColors]);
                           jtPtUnfoldedBayesMC_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][bI]->SetLineColor(colors[bI%nColors]);
                           jtPtUnfoldedBayesMC_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][bI]->SetMarkerSize(1.0);
                           std::cout << "LINE: " << __LINE__ << std::endl;

                           if(bayesPos == 0){
                              std::cout << "UNFOLD CHECK: " << jtPtUnfoldedBayes_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][bayesPos]->GetMinimum() << ", " << jtPtUnfoldedBayes_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][bayesPos]->GetMaximum() << std::endl;

                              jtPtUnfoldedBayes_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][bayesPos]->DrawCopy("HIST E1 P");

                              jtPtRaw_RecoGenAsymm_h[jI][cI][iI][aI]->SetMarkerStyle(1);
                              jtPtRaw_RecoGenAsymm_h[jI][cI][iI][aI]->SetMarkerSize(0.001);
                              jtPtRaw_RecoGenAsymm_h[jI][cI][iI][aI]->SetMarkerColor(0);
                              jtPtRaw_RecoGenAsymm_h[jI][cI][iI][aI]->SetLineColor(1);
                              jtPtRaw_RecoGenAsymm_h[jI][cI][iI][aI]->SetLineWidth(2);
                              jtPtRaw_RecoGenAsymm_h[jI][cI][iI][aI]->DrawCopy("HIST E1 SAME");

                              std::cout << "LINE: " << __LINE__ << std::endl;

                              Double_t rescale = jtPtRaw_RecoGenAsymm_h[jI][cI][iI][aI]->Integral()/recoJtPt_GoodGen_h[jI][cI][iI][mI][aI][sI]->Integral();

                              recoJtPt_GoodGen_h[jI][cI][iI][mI][aI][sI]->Scale(rescale);
                              recoJtPt_GoodGen_h[jI][cI][iI][mI][aI][sI]->SetMarkerStyle(1);
                              recoJtPt_GoodGen_h[jI][cI][iI][mI][aI][sI]->SetMarkerSize(0.001);
                              recoJtPt_GoodGen_h[jI][cI][iI][mI][aI][sI]->SetMarkerColor(0);
                              recoJtPt_GoodGen_h[jI][cI][iI][mI][aI][sI]->SetLineColor(4);
                              recoJtPt_GoodGen_h[jI][cI][iI][mI][aI][sI]->SetLineStyle(2);
                              recoJtPt_GoodGen_h[jI][cI][iI][mI][aI][sI]->SetLineWidth(2);
                              //		    recoJtPt_GoodGen_h[jI][cI][iI][mI][aI][sI]->DrawCopy("HIST E1 SAME");

                              std::cout << "LINE: " << __LINE__ << std::endl;

                              TH2D* initRes_p = new TH2D("initRes_p", "", nRecoJtPtBins[jI][cI], recoJtPtBins[jI][cI], nGenJtPtBins[jI][cI], genJtPtBins[jI][cI]);
                              macroHistToSubsetHist(response_General_h[jI][cI][iI][mI][aI][sI], initRes_p, true, false, true);

                              TH2D* initResInt_p = new TH2D("initResInt_p", "", 1, recoJtPtBins[jI][cI], nGenJtPtBins[jI][cI], genJtPtBins[jI][cI]);
                              macroHistToSubsetHist(response_General_h[jI][cI][iI][mI][aI][sI], initRes_p, true, false, true);
                              macroHistToSubsetHist(response_General_h[jI][cI][iI][mI][aI][sI], initResInt_p, true, false, true);

                              std::cout << "LINE: " << __LINE__ << std::endl;
                              TH1D* initResProj_p = new TH1D("initResProj_p", "", nGenJtPtBins[jI][cI], genJtPtBins[jI][cI]);
                              macroHistToSubsetHistY(initRes_p, initResProj_p, true);
                              TH1D* initResIntProj_p = new TH1D("initResIntProj_p", "", nGenJtPtBins[jI][cI], genJtPtBins[jI][cI]);
                              macroHistToSubsetHistY(initResInt_p, initResIntProj_p, true);

                              rescale = jtPtRaw_RecoGenAsymm_h[jI][cI][iI][aI]->Integral()/initResProj_p->Integral();

                              initResProj_p->Scale(rescale);
                              initResProj_p->SetMarkerStyle(1);
                              initResProj_p->SetMarkerSize(0.001);
                              initResProj_p->SetMarkerColor(0);
                              initResProj_p->SetLineColor(4);
                              initResProj_p->SetLineStyle(2);
                              initResProj_p->SetLineWidth(2);
                              initResProj_p->DrawCopy("HIST E1 SAME");


                              initResIntProj_p->Scale(1./initResIntProj_p->Integral());
                              Double_t total = 0.0;
                              Double_t mark = -1;
                              for(Int_t bI = initResIntProj_p->GetNbinsX(); bI >= 0; --bI){
                                 total += initResIntProj_p->GetBinContent(bI+1);
                                 if(total >= .99){
                                    mark = initResIntProj_p->GetBinCenter(bI+1);
                                    break;
                                 }
                              }

                              TLine* quickLine_p = new TLine();
                              quickLine_p->SetLineStyle(2);
                              quickLine_p->SetLineColor(4);
                              quickLine_p->DrawLine(mark, min, mark, max);
                              delete quickLine_p;

                              leg_p->AddEntry(jtPtRaw_RecoGenAsymm_h[jI][cI][iI][aI], "Folded", "L");
                              //		    leg_p->AddEntry(initResProj_p, "Prior", "L");
                              delete initRes_p;
                              delete initResInt_p;
                              delete initResProj_p;
                              delete initResIntProj_p;
                           }
                           else jtPtUnfoldedBayes_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][bayesPos]->DrawCopy("HIST E1 P SAME"); 		
                           leg_p->AddEntry(jtPtUnfoldedBayes_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][bayesPos], ("Bayes=" + std::to_string(bayesVal[bayesPos])).c_str(), "P L");

                        }
                     }

                     std::cout << "LINE: " << __LINE__ << std::endl;

                     canv_p->cd();
                     pads[0]->cd();
                     gStyle->SetOptStat(0);
                     gPad->SetLogy();
                     gPad->SetTicks(1,2);

                     if(doLogX) gPad->SetLogx();

                     leg_p->Draw("SAME");

                     label_p->DrawLatex(0.60, 0.82, tempStr.c_str());
                     label_p->DrawLatex(0.60, 0.76, centStr2.c_str());
                     label_p->DrawLatex(0.60, 0.70, jtAbsEtaStr2.c_str());
                     label_p->DrawLatex(0.60, 0.64, resStr.c_str());
                     label_p->DrawLatex(0.60, 0.58, idStr[iI].c_str());
                     label_p->DrawLatex(0.60, 0.52, systStr[sI].c_str());

                     canv_p->cd();
                     pads[1]->cd();

                     max = -1;
                     min = 100;

                     std::cout << "LINE: " << __LINE__ << std::endl;

                     for(Int_t bI = 0; bI < nBayes; ++bI){
                        clones_p[bI] = (TH1D*)jtPtUnfoldedBayes_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][bI]->Clone(("clone_" + std::to_string(bI)).c_str());
                        clones_p[bI]->Divide(jtPtUnfoldedBayes_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][0]);
                        Double_t tempMax = getMax(clones_p[bI]);
                        Double_t tempMin = getMinGTZero(clones_p[bI]);
                        if(max < tempMax) max = tempMax;
                        if(min > tempMin) min = tempMin;
                     }

                     interval = (max - min)/10.;
                     max += interval;
                     min -= interval;

                     if(max > 1.15) max = 1.15;
                     if(min < .85) min = 0.85;

                     for(Int_t bI = 0; bI < nBayes; ++bI){
                        clones_p[bI]->SetMaximum(max);
                        clones_p[bI]->SetMinimum(min);

                        clones_p[bI]->GetYaxis()->SetTitle("Ratio w/ Bayes0");
                        clones_p[bI]->GetYaxis()->SetTitleSize(9);
                        clones_p[bI]->GetYaxis()->SetLabelSize(9);
                        clones_p[bI]->GetYaxis()->SetTitleOffset(clones_p[bI]->GetYaxis()->GetTitleOffset()*1.5);
                        clones_p[bI]->GetYaxis()->SetNdivisions(505);	       		

                        if(bI < nBayesDraw){
                           Int_t bayesPos = bI;
                           if(bI == nBayesDraw-1 && histTermPos[jI][cI][iI][mI][aI][sI] >= nBayesDraw) bayesPos = histTermPos[jI][cI][iI][mI][aI][sI];

                           if(bayesPos == 0){
                              Double_t min = clones_p[bayesPos]->GetMinimum();
                              Double_t max = clones_p[bayesPos]->GetMaximum();

                              clones_p[bayesPos]->SetMinimum(0.9);
                              clones_p[bayesPos]->SetMaximum(1.1);

                              clones_p[bayesPos]->GetYaxis()->SetTitleFont(43);
                              clones_p[bayesPos]->GetXaxis()->SetTitleFont(43);
                              clones_p[bayesPos]->GetYaxis()->SetTitleSize(14);
                              clones_p[bayesPos]->GetXaxis()->SetTitleSize(14);

                              clones_p[bayesPos]->DrawCopy("HIST E1 P");

                              clones_p[bayesPos]->SetMinimum(min);
                              clones_p[bayesPos]->SetMaximum(max);
                           }
                           else clones_p[bayesPos]->DrawCopy("HIST E1 P SAME");
                        }
                     }

                     gStyle->SetOptStat(0);
                     if(doLogX) gPad->SetLogx();

                     for(Int_t bI = 0; bI < nBayes; ++bI){
                        delete clones_p[bI];
                        clones_p[bI] = NULL;
                     }


                     std::cout << "LINE: " << __LINE__ << std::endl;

                     canv_p->cd();
                     pads[3]->cd();

                     //Refold EDITING HERE		
                     jtPtRaw_RecoGenAsymm_h[jI][cI][iI][aI]->SetMaximum(globalMax*20.);
                     jtPtRaw_RecoGenAsymm_h[jI][cI][iI][aI]->SetMinimum(globalMin/20.);
                     jtPtRaw_RecoGenAsymm_h[jI][cI][iI][aI]->DrawCopy("HIST E1");

                     //	      label_p->DrawLatex(0.35, .94, "Refolding + Comp. to data");

                     std::cout << "LINE: " << __LINE__ << ", " << systStr[sI] << std::endl;

                     for(Int_t bI = 0; bI < nBayesDraw; ++bI){
                        Int_t bayesPos = bI;
                        if(bI == nBayesDraw-1 && histTermPos[jI][cI][iI][mI][aI][sI] >= nBayesDraw) bayesPos = histTermPos[jI][cI][iI][mI][aI][sI];

                        std::cout << "LINE: " << __LINE__ << std::endl;

                        canv_p->cd();
                        pads[3]->cd();

                        TH1D* tempClone = (TH1D*)rooResponse_Reweighted_h[jI][cI][iI][mI][aI][sI]->ApplyToTruth(jtPtUnfoldedBayes_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][bayesPos], "tempClone");
                        std::cout << "LINE: " << __LINE__ << std::endl;

                        centerTitles(tempClone);
                        tempClone->SetMarkerStyle(styles[bayesPos%nStyles]);
                        tempClone->SetMarkerColor(colors[bayesPos%nColors]);
                        tempClone->SetLineColor(colors[bayesPos%nColors]);
                        tempClone->SetMarkerSize(1.0);

                        tempClone->GetXaxis()->SetTitleFont(43);
                        tempClone->GetYaxis()->SetTitleFont(43);
                        tempClone->GetXaxis()->SetLabelFont(43);
                        tempClone->GetYaxis()->SetLabelFont(43);

                        tempClone->GetXaxis()->SetTitleSize(14);
                        tempClone->GetYaxis()->SetTitleSize(14);
                        tempClone->GetXaxis()->SetLabelSize(14);
                        tempClone->GetYaxis()->SetLabelSize(14);

                        tempClone->GetXaxis()->SetTitleOffset(5.);
                        tempClone->GetYaxis()->SetTitleOffset(1.5);

                        tempClone->SetTitle("");

                        tempClone->DrawCopy("HIST E1 P SAME");

                        canv_p->cd();
                        pads[4]->cd();

                        tempClone->Divide(jtPtRaw_RecoGenAsymm_h[jI][cI][iI][aI]);

                        if(bayesPos == 0){
                           tempClone->SetMaximum(1.25);
                           tempClone->SetMinimum(0.75);
                           tempClone->GetYaxis()->SetNdivisions(505);
                           tempClone->GetYaxis()->SetTitle("Ratio");
                           tempClone->DrawCopy("HIST E1 P");
                        }
                        else tempClone->DrawCopy("HIST E1 P SAME");

                        delete tempClone;
                     }

                     std::cout << "LINE: " << __LINE__ << std::endl;

                     drawWhiteBox(900, 1100, .00, min*.999);

                     label_p->SetNDC(0);
                     //	      label_p->DrawLatex(100, min - interval*3, "100");
                     if(!isBigJet) label_p->DrawLatex(200, .77, "200");
                     label_p->DrawLatex(400, .77, "400");
                     label_p->DrawLatex(600, .77, "600");
                     label_p->DrawLatex(1000, .77, "1000");
                     label_p->SetNDC(1);

                     std::cout << "LINE: " << __LINE__ << ", " << histTermPos[jI][cI][iI][mI][aI][sI] << std::endl;

                     const Int_t nPearsToDraw = 3;
                     for(Int_t bI = 0; bI < nPearsToDraw; ++bI){
                        Int_t bayesPos = bI+1;
                        if(bI == nPearsToDraw-1){
                           if(histTermPos[jI][cI][iI][mI][aI][sI] >= nPearsToDraw) bayesPos = histTermPos[jI][cI][iI][mI][aI][sI];
                           else continue;
                        }

                        std::cout << "LINE: " << __LINE__ << std::endl;

                        canv_p->cd();
                        pads[5+bI]->cd();

                        std::cout << "LINE: " << __LINE__ << std::endl;

                        TH2D* pearsonToPlot = NULL;
                        TH2D* covarianceToPlot = NULL;
                        getPearsTMatrix(bayes_p[jI][cI][iI][mI][aI][sI][bayesPos], &covarianceToPlot, &pearsonToPlot);

                        std::cout << "LINE: " << __LINE__ << std::endl;

                        pearsonToPlot->SetMarkerSize(1.75);		  
                        gStyle->SetPaintTextFormat("1.3f");

                        std::cout << "LINE: " << __LINE__ << std::endl;

                        std::vector<Float_t> tempVals;
                        for(Int_t bIX = 0; bIX < pearsonToPlot->GetNbinsX(); ++bIX){
                           for(Int_t bIY = 0; bIY < pearsonToPlot->GetNbinsY(); ++bIY){
                              Double_t val = pearsonToPlot->GetBinContent(bIX+1, bIY+1);
                              tempVals.push_back(val);
                              pearsonToPlot->SetBinContent(bIX+1, bIY+1, TMath::Abs(val));
                           }
                        }

                        std::cout << "LINE: " << __LINE__ << std::endl;

                        //		pearsonToPlot->SetMaximum(1.05);
                        //		pearsonToPlot->SetMinimum(-0.05);
                        //		pearsonToPlot->DrawCopy("COL");
                        pearsonToPlot->DrawCopy("COLZ");
                        //		std::cout << "COVARIANCE PRINCE" << std::endl;
                        //		pearsonToPlot->Print("ALL");

                        std::cout << "LINE: " << __LINE__ << std::endl;

                        Int_t size = label_p->GetTextSize();
                        label_p->SetTextSize(3*size/4);

                        label_p->SetNDC(0);
                        unsigned int tempPos = 0;


                        std::cout << "LINE: " << __LINE__ << std::endl;

                        const Int_t nBinsMax = 100;
                        Int_t nGenBins = -1;
                        Double_t genBins[nBinsMax];

                        std::cout << "LINE: " << __LINE__ << std::endl;

                        for(Int_t bIX = 0; bIX < jtPtUnfoldedBayes_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][bI]->GetNbinsX()+1; ++bIX){
                           genBins[bIX] = jtPtUnfoldedBayes_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][bI]->GetBinLowEdge(bIX+1);
                           ++nGenBins;
                        }


                        std::cout << "LINE: " << __LINE__ << std::endl;

                        for(Int_t bIX = 0; bIX < pearsonToPlot->GetNbinsX(); ++bIX){
                           Double_t xLow = pearsonToPlot->GetXaxis()->GetBinLowEdge(bIX+1);

                           Double_t xLow2 = genBins[bIX];
                           Double_t xHigh2 = genBins[bIX+1];
                           Double_t xCent = pearsonToPlot->GetXaxis()->GetBinCenter(bIX+1);
                           for(Int_t bIY = 0; bIY < pearsonToPlot->GetNbinsY(); ++bIY){
                              Double_t yLow = genBins[bIY];
                              Double_t yHigh = genBins[bIY+1];
                              Double_t yCent = pearsonToPlot->GetYaxis()->GetBinCenter(bIY+1);
                              if(TMath::Abs(tempVals[tempPos]) > 0.15) label_p->DrawLatex((xCent + xLow)/2., yCent, ("#color[2]{" + prettyString(tempVals[tempPos], 3, false) + "}").c_str());
                              else label_p->DrawLatex((xCent + xLow)/2., yCent, prettyString(tempVals[tempPos], 3, false).c_str());

                              if(bI == 0 && sI == 0) std::cout << "TEMP VALS," << xLow2 << "," << xHigh2 << "," << yLow << "," << yHigh << ", " << tempVals[tempPos] << std::endl;
                              ++tempPos;
                           }
                        }
                        std::cout << "LINE: " << __LINE__ << std::endl;

                        if(bI == 0 && sI == 0) std::cout << "INTEGRAL: " << pearsonToPlot->Integral() << std::endl;

                        label_p->SetTextSize(size);
                        label_p->SetNDC();
                        //		label_p->DrawLatex(0.35, .94, ("#bf{Pearson, Iteration " + std::to_string(bayesVal[bayesPos])+ "}").c_str());

                        std::cout << "LINE: " << __LINE__ << std::endl;

                        //	      gPad->SetLogx();
                        //	      gPad->SetLogy();

                        drawWhiteBox(-0.6, pearsonToPlot->GetXaxis()->GetBinLowEdge(pearsonToPlot->GetNbinsX()+1)+.25, -0.6, -.01);
                        drawWhiteBox(-0.6, -.01, -0.6, pearsonToPlot->GetXaxis()->GetBinLowEdge(pearsonToPlot->GetNbinsX()+1)+.25);

                        label_p->SetNDC(0);

                        std::cout << "LINE: " << __LINE__ << std::endl;

                        Int_t binLow = -1;
                        Int_t binHigh = -1;

                        for(Int_t bIY = 0; bIY < jtPtUnfoldedBayes_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][bI]->GetNbinsX()+1; ++bIY){
                           Int_t binLowEdge = jtPtUnfoldedBayes_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][bI]->GetXaxis()->GetBinLowEdge(bIY+1);

                           if(binLow == -1){
                              if(binLowEdge >= ptBinLow[jI][cI]-1) binLow = bIY;
                           }
                           if(binHigh == -1){
                              if(binLowEdge >= ptBinHigh[jI][cI]-1) binHigh = bIY;
                           }

                           label_p->DrawLatex(-0.7, 0.05 + bIY, std::to_string(int(binLowEdge)).c_str());
                           label_p->DrawLatex(0.05 + bIY, -0.3, std::to_string(int(binLowEdge)).c_str());
                        }

                        std::cout << "LINE: " << __LINE__ << std::endl;
                        label_p->SetNDC(1);

                        TLine* line_p = new TLine();
                        line_p->SetLineStyle(2);
                        line_p->SetLineColor(2);
                        line_p->SetLineWidth(line_p->GetLineWidth()*2);

                        line_p->DrawLine(binLow, binLow, binLow, binHigh);
                        line_p->DrawLine(binLow, binHigh, binHigh, binHigh);
                        line_p->DrawLine(binLow, binLow, binHigh, binLow);
                        line_p->DrawLine(binHigh, binLow, binHigh, binHigh);

                        delete line_p;
                        delete pearsonToPlot;
                     }

                     std::cout << "LINE: " << __LINE__ << std::endl;

                     canv_p->cd();
                     std::cout << "LINE: " << __LINE__ << std::endl;
                     pads[8]->cd();
                     std::cout << "LINE: " << __LINE__ << std::endl;

                     std::string centBinsStr = "Cent" + std::to_string(centBinsLow[cI]) + "to" + std::to_string(centBinsHi[cI]);
                     std::cout << responseJetDirList[jI] << ", " << centBinsStr << std::endl;
                     std::cout << "LINE: " << __LINE__ << std::endl;

                     Int_t tempBestSVD = TMath::Min(nGenJtPtBins[jI][cI], svdTerm(responseJetDirList[jI], centBinsStr));
                     std::cout << "LINE: " << __LINE__ << ", " << tempBestSVD << std::endl;

                     std::cout << "LINE: " << __LINE__ << std::endl;

                     if(tempBestSVD > 0){
                        TH1D* tempSVD_p = svd_p[jI][cI][iI][mI][aI][sI][tempBestSVD-1]->GetD();
                        tempSVD_p->SetMarkerStyle(20);
                        tempSVD_p->SetMarkerSize(1);
                        tempSVD_p->SetMarkerColor(1);
                        tempSVD_p->SetLineColor(1);

                        tempSVD_p->SetMinimum(0.2);

                        tempSVD_p->DrawCopy("HIST E1 P");
                        gPad->SetLogy();
                     }
                     std::cout << "LINE: " << __LINE__ << std::endl;


                     canv_p->cd();
                     pads[9]->cd();

                     std::cout << "LINE: " << __LINE__ << std::endl;

                     for(Int_t bI = 0; bI < nPearsToDraw; ++bI){

                        if(bI == 0){
                           jtPtUnfoldedBayesMC_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][bI]->DrawCopy("HIST E1 P");
                           gPad->SetLogx(); 
                           gPad->SetLogy();
                        }
                        else jtPtUnfoldedBayesMC_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][bI]->DrawCopy("HIST E1 P SAME");
                     }

                     std::cout << "LINE: " << __LINE__ << std::endl;

                     jtPtUnfoldedMCTruth_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI]->SetLineColor(4);
                     jtPtUnfoldedMCTruth_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI]->SetLineStyle(2);
                     jtPtUnfoldedMCTruth_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI]->SetLineStyle(4);
                     jtPtUnfoldedMCTruth_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI]->GetYaxis()->SetTitleFont(43);
                     jtPtUnfoldedMCTruth_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI]->GetXaxis()->SetTitleFont(43);
                     jtPtUnfoldedMCTruth_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI]->GetYaxis()->SetTitleSize(14);
                     jtPtUnfoldedMCTruth_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI]->GetXaxis()->SetTitleSize(14);
                     jtPtUnfoldedMCTruth_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI]->DrawCopy("HIST E1 P SAME");

                     std::cout << "LINE: " << __LINE__ << std::endl;

                     label_p->DrawLatex(0.2, 0.4, "Monte Carlo closure check");

                     std::cout << "LINE: " << __LINE__ << std::endl;

                     canv_p->cd();
                     pads[11]->cd();
                     jtPtUnfoldedMCTruth_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI]->Divide(jtPtUnfoldedBayesMC_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][3]);
                     jtPtUnfoldedMCTruth_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI]->SetMaximum(1.2);
                     jtPtUnfoldedMCTruth_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI]->SetMinimum(0.8);
                     jtPtUnfoldedMCTruth_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI]->SetLineColor(1);
                     jtPtUnfoldedMCTruth_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI]->SetLineStyle(1);
                     jtPtUnfoldedMCTruth_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI]->GetYaxis()->SetTitle("Divide");

                     jtPtUnfoldedMCTruth_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI]->GetYaxis()->SetTitleFont(43);
                     jtPtUnfoldedMCTruth_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI]->GetXaxis()->SetTitleFont(43);
                     jtPtUnfoldedMCTruth_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI]->GetYaxis()->SetTitleSize(14);
                     jtPtUnfoldedMCTruth_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI]->GetXaxis()->SetTitleSize(14);
                     centerTitles(jtPtUnfoldedMCTruth_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI]);
                     jtPtUnfoldedMCTruth_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI]->DrawCopy("HIST E1 P");
                     gPad->SetLogx();

                     std::cout << "LINE: " << __LINE__ << std::endl;


                     canv_p->cd();
                     pads[10]->cd();

                     TH2D* initRes_p = new TH2D("initRes_p", "", nRecoJtPtBins[jI][cI], recoJtPtBins[jI][cI], nGenJtPtBins[jI][cI], genJtPtBins[jI][cI]);
                     macroHistToSubsetHist(response_General_h[jI][cI][iI][mI][aI][sI], initRes_p, true, false, true);
                     if(/*doClean*/false){
                        cleanMatrix(initRes_p);
                     }


                     for(Int_t bIY = 0; bIY < initRes_p->GetYaxis()->GetNbins(); ++bIY){
                        Double_t total = 0.0;

                        for(Int_t bIX = 0; bIX < initRes_p->GetXaxis()->GetNbins(); ++bIX){
                           total += initRes_p->GetBinContent(bIX+1, bIY+1);
                        }

                        if(total <= 0.0) continue;

                        for(Int_t bIX = 0; bIX < initRes_p->GetXaxis()->GetNbins(); ++bIX){
                           Double_t val = initRes_p->GetBinContent(bIX+1, bIY+1)/total;
                           Double_t err = initRes_p->GetBinError(bIX+1, bIY+1)/total;

                           initRes_p->SetBinContent(bIX+1, bIY+1, val);
                           initRes_p->SetBinError(bIX+1, bIY+1, err);
                        }
                     }

                     initRes_p->GetXaxis()->SetTitle("Reconstructed Jet p_{T} (GeV/c)");
                     initRes_p->GetYaxis()->SetTitle("Truth Jet p_{T} (GeV/c)");
                     centerTitles(initRes_p);
                     initRes_p->DrawCopy("COL TEXT");
                     gPad->SetLogx();
                     gPad->SetLogy();
                     gStyle->SetPaintTextFormat("1.3f");

                     delete initRes_p;

                     canv_p->cd();
                     pads[3]->cd();
                     gPad->SetLogy();
                     if(doLogX) gPad->SetLogx();

                     canv_p->cd();
                     pads[4]->cd();
                     if(doLogX) gPad->SetLogx();

                     std::cout << "LINE: " << __LINE__ << std::endl;

                     //	      label_p->SetNDC(0);
                     canv_p->cd();
                     pads[0]->cd();

                     const std::string tempHistTag = tempStr + "_" + centStr + "_" + idStr[iI] + "_ResponseMod" + prettyString(responseMod[mI], 2, true) + "_" + jtAbsEtaStr + "_" + systStr[sI];
                     histTag.push_back(tempHistTag);

                     if(terminalChi2 >= 0){
                        histBestBayes.push_back(bayesVal[terminalChi2]);
                        histBestBayesPos.push_back(terminalChi2);
                     }
                     else{
                        histBestBayes.push_back(-1);
                        histBestBayesPos.push_back(-1);
                     }	    

                     if(false && terminalPos5 >= 0) label_p->DrawLatex(0.25, 0.33, ("Term. 5% at Bayes=" + std::to_string(bayesVal[terminalPos5])).c_str());
                     else if(false) label_p->DrawLatex(0.25, 0.33, "Doesn't term. at 5% level");

                     if(false && terminalPos5AndPears >= 0) label_p->DrawLatex(0.25, 0.26, ("Term. 5%+PearsMin at Bayes=" + std::to_string(bayesVal[terminalPos5AndPears])).c_str());
                     else if(false) label_p->DrawLatex(0.25, 0.26, "Doesn't term. at 5% level+Pears");

                     if(terminalChi2 >= 0) label_p->DrawLatex(0.35, 0.12, ("Term. Chi2 at Bayes=" + std::to_string(bayesVal[terminalChi2])).c_str());
                     else label_p->DrawLatex(0.35, 0.12, "Doesn't term. at Chi2");

                     label_p->DrawLatex(0.35, 0.2, ("Term. SVD at kReg=" + std::to_string(svdTerm(responseJetDirList[jI], centStr))).c_str());


                     if(false && terminalPos >= 0){
                        label_p->DrawLatex(0.25, 0.19, ("Term. 1% at Bayes=" + std::to_string(bayesVal[terminalPos])).c_str());
                        Double_t maxDelta = 0;
                        Double_t lastDelta = 0;

                        for(Int_t bIX = 0; bIX < jtPtUnfoldedBayes_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][terminalPos]->GetNbinsX(); ++bIX){
                           double center = (jtPtUnfoldedBayes_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][terminalPos]->GetBinLowEdge(bIX+1) + jtPtUnfoldedBayes_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][terminalPos]->GetBinLowEdge(bIX+2))/2.;

                           if(center < lowPtTruncVal) continue;
                           if(center > 1000.) continue;

                           double content1 = jtPtUnfoldedBayes_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][terminalPos]->GetBinContent(bIX+1);
                           double content1Max = jtPtUnfoldedBayes_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][nBayes-1]->GetBinContent(bIX+1);
                           double content1MaxMin1 = jtPtUnfoldedBayes_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][nBayes-2]->GetBinContent(bIX+1);

                           if(content1 == 0 && content1Max != 0){
                              maxDelta = 100;
                           }
                           else if(TMath::Abs(content1 - content1Max)/content1 > maxDelta){
                              maxDelta = TMath::Abs(content1 - content1Max)/content1;
                           }

                           if(content1MaxMin1 == 0 && content1Max != 0){
                              lastDelta = 100;
                           }
                           else if(TMath::Abs(content1MaxMin1 - content1Max)/content1MaxMin1 > lastDelta){
                              lastDelta = TMath::Abs(content1MaxMin1 - content1Max)/content1MaxMin1;
                           }
                        }

                        label_p->DrawLatex(0.25, 0.12, ("MaxDelta: " + prettyString(maxDelta*100, 2, false) + "%").c_str());
                        label_p->DrawLatex(0.25, 0.05, ("LastDelta: " + prettyString(lastDelta*100, 2, false) + "%").c_str());
                     }
                     else if(false) label_p->DrawLatex(0.25, 0.19, ("Doesn't term. for nBayes=" + std::to_string(nBayes+1)).c_str());

                     canv_p->cd();
                     pads[2]->cd();

                     label_p->SetNDC(0);

                     //	      label_p->DrawLatex(100, min - interval*2, "100");
                     //	      label_p->DrawLatex(200, min - interval*2, "200");
                     //	      label_p->DrawLatex(400, min - interval*2, "400");
                     //	      label_p->DrawLatex(600, min - interval*2, "600");
                     //	      label_p->DrawLatex(1000, min - interval*2, "1000");

                     for(Int_t pI = 0; pI < nPads; ++pI){
                        canv_p->cd();
                        pads[pI]->cd();
                        gPad->RedrawAxis();
                        gPad->SetTicks(1, 2);
                     }

                     const std::string saveName = "jtPtUnfoldedBayes_" + tempStr + "_" + centStr + "_" + idStr[iI] + "_" + resStr + "_" + jtAbsEtaStr + "_" + tempSystStr +  "AllBayes_RecoGenAsymm_" + cleanStr + "_" + toyStr + "_" + debugStr + addedTagStr + dateStr + ".pdf";

                     std::string titleStr = tempStr + ", " + centStr2 + ", " + idStr[iI] + ", $" + jtAbsEtaStr2 + "$, " + tempSystStr;
                     while(titleStr.find("#") != std::string::npos){titleStr.replace(titleStr.find("#"), 1, "\\");}
                     slideTitlesPerAlgo[jI].push_back(titleStr);
                     pdfPerSlidePerAlgo[jI].push_back({});
                     pdfPerSlidePerAlgo[jI][pdfPerSlidePerAlgo[jI].size()-1].push_back(saveName);
                     label_p->SetTextSize(18);
                     label_p->SetNDC();
                     for(Int_t padI = 0; padI < nPads; ++padI){
                        pads[padI]->cd();
                        if(padI != 1 && padI != 2 && padI != 4 && padI != 11){
                           if(padI == 0) label_p->DrawLatex(0.15, 0.95, "Unfolded Distributions");
                           else if(padI == 3) label_p->DrawLatex(0.15, 0.95, "Refolded Distributions");
                           else if(padI == 5) label_p->DrawLatex(0.15, 0.95, "PEARSON COEFFICIENTS ITER 2");
                           else if(padI == 6) label_p->DrawLatex(0.15, 0.95, "PEARSON COEFFICIENTS ITER 3");
                           else if(padI == 7) label_p->DrawLatex(0.15, 0.95, "PEARSON COEFFICIENTS TERMINATING");
                           else if(padI == 8) label_p->DrawLatex(0.15, 0.95, "SVD Unfolding d-vector");
                           else if(padI == 9) label_p->DrawLatex(0.15, 0.95, "MC Unfolding Distributions");
                           else if(padI == 10) label_p->DrawLatex(0.15, 0.95, "Response Matrix, Unity in Y");
                        }
                     }

                     const std::string finalSaveName = "pdfDir/" + dateStr + "/Unfold_" + tempStr + "/" + saveName;
                     quietSaveAs(canv_p, finalSaveName);
                     ++nPDFTotal;

                     for(Int_t pI = 0; pI < nPads; ++pI){
                        delete pads[pI];
                     }

                     delete canv_p;
                     delete leg_p;
                     delete label_p;
                  }
               }

               for(Int_t aI = 0; aI < nJtAbsEtaBins; ++aI){
                  const std::string jtAbsEtaStr = "AbsEta" + prettyString(jtAbsEtaBinsLow[aI], 1, true) + "to" + prettyString(jtAbsEtaBinsHi[aI], 1, true);

                  std::vector<std::vector<std::vector<double> > > vals;

                  for(Int_t bIX = 0; bIX < covariance_h[jI][cI][iI][mI][aI][0]->GetXaxis()->GetNbins(); ++bIX){
                     vals.push_back({});
                     for(Int_t bIY = 0; bIY < covariance_h[jI][cI][iI][mI][aI][0]->GetYaxis()->GetNbins(); ++bIY){
                        vals[bIX].push_back({});
                     }
                  }

                  for(Int_t covI = 0; covI < nCov; ++covI){
                     covariance_h[jI][cI][iI][mI][aI][covI]->Write("", TObject::kOverwrite);

                     for(Int_t bIX = 0; bIX < covariance_h[jI][cI][iI][mI][aI][covI]->GetXaxis()->GetNbins(); ++bIX){
                        for(Int_t bIY = 0; bIY < covariance_h[jI][cI][iI][mI][aI][covI]->GetYaxis()->GetNbins(); ++bIY){
                           double tempVal = covariance_h[jI][cI][iI][mI][aI][covI]->GetBinContent(bIX+1, bIY+1);
                           vals[bIX][bIY].push_back(tempVal);
                        }
                     }

                     delete covariance_h[jI][cI][iI][mI][aI][covI];
                  }

                  for(unsigned int bIX = 0; bIX < vals.size(); ++bIX){
                     for(unsigned int bIY = 0; bIY < vals[bIX].size(); ++bIY){
                        std::vector<double> tempVect = vals[bIX][bIY];
                        double mean = getVectMean(tempVect);
                        double sigma = getVectSigma(tempVect);

                        covariance_Mean_h[jI][cI][iI][mI][aI]->SetBinContent(bIX+1, bIY+1, mean);
                        covariance_Mean_h[jI][cI][iI][mI][aI]->SetBinError(bIX+1, bIY+1, 0.0);

                        covariance_Sigma_h[jI][cI][iI][mI][aI]->SetBinContent(bIX+1, bIY+1, sigma);
                        covariance_Sigma_h[jI][cI][iI][mI][aI]->SetBinError(bIX+1, bIY+1, 0.0);

                        covariance_SigmaOverMean_h[jI][cI][iI][mI][aI]->SetBinContent(bIX+1, bIY+1, sigma/mean);
                        covariance_SigmaOverMean_h[jI][cI][iI][mI][aI]->SetBinError(bIX+1, bIY+1, 0.0);
                     }
                  }

                  covariance_Mean_h[jI][cI][iI][mI][aI]->Write("", TObject::kOverwrite);
                  delete covariance_Mean_h[jI][cI][iI][mI][aI];

                  covariance_Sigma_h[jI][cI][iI][mI][aI]->Write("", TObject::kOverwrite);
                  delete covariance_Sigma_h[jI][cI][iI][mI][aI];

                  covariance_SigmaOverMean_h[jI][cI][iI][mI][aI]->Write("", TObject::kOverwrite);
                  delete covariance_SigmaOverMean_h[jI][cI][iI][mI][aI];

                  for(Int_t sI = 0; sI < nSyst; ++sI){
                     const std::string tempSystStr = systStr[sI] + "_";
                     std::cout << response_Unfold_h[jI][cI][iI][mI][aI][sI]->GetName() << std::endl;

                     response_Unfold_h[jI][cI][iI][mI][aI][sI]->Write("", TObject::kOverwrite);
                     delete response_Unfold_h[jI][cI][iI][mI][aI][sI];

                     delete jtPtUnfoldedMCTruth_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI];

                     TCanvas* canv_p = new TCanvas("canv_p", "", 450, 450);
                     TLegend* leg_p = new TLegend(0.2, 0.65, 0.5, 0.95);
                     leg_p->SetBorderSize(0);
                     leg_p->SetFillColor(0);
                     leg_p->SetFillStyle(0);
                     leg_p->SetTextFont(43);
                     leg_p->SetTextSize(12);

                     for(Int_t bI = 1; bI < nBayes; ++bI){
                        TH1D* tempHist_p = (TH1D*)jtPtUnfoldedBayes_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][bI]->Clone("tempHist");


                        tempHist_p->Divide(jtPtUnfoldedBayes_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][bI-1]);
                        tempHist_p->SetMaximum(1.2);
                        tempHist_p->SetMinimum(0.8);

                        tempHist_p->SetMarkerStyle(styles[bI%nStyles]);
                        tempHist_p->SetMarkerSize(1);
                        tempHist_p->SetMarkerColor(vg.getColor(bI%nColors));
                        tempHist_p->SetLineColor(vg.getColor(bI%nColors));

                        jtPtUnfoldedBayes_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][bI]->SetMarkerStyle(styles[bI%nStyles]);
                        jtPtUnfoldedBayes_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][bI]->SetMarkerSize(1);
                        jtPtUnfoldedBayes_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][bI]->SetMarkerColor(vg.getColor(bI%nColors));
                        jtPtUnfoldedBayes_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][bI]->SetLineColor(vg.getColor(bI%nColors));

                        if(bI == 1) tempHist_p->DrawCopy("HIST E1 P");
                        else tempHist_p->DrawCopy("SAME HIST E1 P");

                        std::string legStr = "Bayes=" + std::to_string(bayesVal[bI]);
                        if(histTermPos[jI][cI][iI][mI][aI][sI] == bI) legStr = legStr + ",*";

                        leg_p->AddEntry(jtPtUnfoldedBayes_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][bI], legStr.c_str(), "P L");

                        delete tempHist_p;
                     }

                     leg_p->Draw("SAME");
                     std::string saveName = "iterRatioBayes_" + tempStr + "_" + centStr + "_" + idStr[iI] + "_" + resStr + "_" + jtAbsEtaStr + "_" + tempSystStr +  "AllBayes_RecoGenAsymm_" + cleanStr + "_" + toyStr + "_" + debugStr + addedTagStr + dateStr + ".pdf";;
                     saveName = "pdfDir/" + dateStr + "/Unfold_" + tempStr + "/" + saveName;
                     quietSaveAs(canv_p, saveName);

                     delete leg_p;
                     delete canv_p;


                     canv_p = new TCanvas("canv_p", "", 450, 450);
                     leg_p = new TLegend(0.2, 0.65, 0.5, 0.95);
                     leg_p->SetBorderSize(0);
                     leg_p->SetFillColor(0);
                     leg_p->SetFillStyle(0);
                     leg_p->SetTextFont(43);
                     leg_p->SetTextSize(12);

                     lineChi2Svd[jI][cI][iI][mI][aI][sI][0] = 999;

                     Double_t maxProb = 0.00000000;
                     Int_t maxProbPos = -1;

                     Int_t sigma2Pos = -1;
                     Double_t sigma2Prob = 0.000;

                     for(Int_t svdI = 1; svdI < TMath::Min(nGenJtPtBins[jI][cI], nSvd); ++svdI){
                        Double_t newGenBins[nMaxJtPtBins+1];
                        for(Int_t gI = 0; gI < nGenJtPtBins[jI][cI]-1; ++gI){
                           newGenBins[gI] = genJtPtBins[jI][cI][gI+1];
                        }

                        TH1D* tempHist_p = new TH1D("tempHist", "", nGenJtPtBins[jI][cI]-2, newGenBins);
                        TH1D* tempHist2_p = new TH1D("tempHist2", "", nGenJtPtBins[jI][cI]-2, newGenBins);
                        macroHistToSubsetHist(jtPtUnfoldedSvd_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][svdI], tempHist_p);
                        macroHistToSubsetHist(jtPtUnfoldedSvd_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][svdI-1], tempHist2_p);


                        std::cout << "1 FIT" << std::endl;
                        tempHist_p->Divide(tempHist2_p);
                        TF1* lineFit_p = new TF1("lineFit_p", "[0]", genJtPtBins[jI][cI][0], genJtPtBins[jI][cI][nGenJtPtBins[jI][cI]]);
                        tempHist_p->Print("ALL");
                        lineFit_p->SetParameter(0, 1.);
                        lineFit_p->SetParLimits(0, .99, 1.01);
                        tempHist_p->Fit("lineFit_p", "M N", "",  genJtPtBins[jI][cI][0], genJtPtBins[jI][cI][nGenJtPtBins[jI][cI]]);
                        lineChi2Svd[jI][cI][iI][mI][aI][sI][svdI] = lineFit_p->GetProb();	       
                        if(lineChi2Svd[jI][cI][iI][mI][aI][sI][svdI] > 0.95 && sigma2Pos < 0){
                           sigma2Pos = svdI;
                           sigma2Prob = lineChi2Svd[jI][cI][iI][mI][aI][sI][svdI];
                        }

                        if(lineChi2Svd[jI][cI][iI][mI][aI][sI][svdI] > maxProb){
                           maxProb = lineChi2Svd[jI][cI][iI][mI][aI][sI][svdI];
                           maxProbPos = svdI;
                        }

                        tempHist_p->SetMaximum(1.5);
                        tempHist_p->SetMinimum(0.5);

                        tempHist_p->SetMarkerStyle(styles[svdI%nStyles]);
                        tempHist_p->SetMarkerSize(1);
                        tempHist_p->SetMarkerColor(vg.getColor(svdI%nColors));
                        tempHist_p->SetLineColor(vg.getColor(svdI%nColors));

                        jtPtUnfoldedSvd_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][svdI]->SetMarkerStyle(styles[svdI%nStyles]);
                        jtPtUnfoldedSvd_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][svdI]->SetMarkerSize(1);
                        jtPtUnfoldedSvd_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][svdI]->SetMarkerColor(vg.getColor(svdI%nColors));
                        jtPtUnfoldedSvd_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][svdI]->SetLineColor(vg.getColor(svdI%nColors));

                        if(svdI == 1) tempHist_p->DrawCopy("HIST E1 P");
                        else tempHist_p->DrawCopy("SAME HIST E1 P");

                        std::string legStr = "Svd=" + std::to_string(svdI+1) + ", " + prettyString(lineChi2Svd[jI][cI][iI][mI][aI][sI][svdI], 4, false);
                        //		if(histTermPos[jI][cI][iI][mI][aI][sI] == svdI) legStr = legStr + ",*";

                        leg_p->AddEntry(jtPtUnfoldedSvd_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][svdI], legStr.c_str(), "P L");

                        delete lineFit_p;
                        delete tempHist_p;
                     }

                     if(sigma2Pos >= 0){
                        histBestSvd.push_back(sigma2Pos+1);
                        histBestSvdProb.push_back(sigma2Prob);
                     }
                     else if(maxProbPos >= 0 && maxProb > 0.05){
                        histBestSvd.push_back(maxProbPos+1);
                        histBestSvdProb.push_back(maxProb);
                     }
                     else{
                        histBestSvd.push_back(-1);
                        histBestSvdProb.push_back(0.0);
                     }

                     std::string centBinsStr = "Cent" + std::to_string(centBinsLow[cI]) + "to" + std::to_string(centBinsHi[cI]);

                     /*
                        if(responseJetDirList[jI].find("akCs10") != std::string::npos){
                        if(centBinsStr.find("Cent0to10") != std::string::npos) histBestSvd[histBestSvd.size()-1] = 3;
                        else if(centBinsStr.find("Cent10to30") != std::string::npos) histBestSvd[histBestSvd.size()-1] = 4;
                        else if(centBinsStr.find("Cent30to50") != std::string::npos) histBestSvd[histBestSvd.size()-1] = 4;
                        else if(centBinsStr.find("Cent50to90") != std::string::npos) histBestSvd[histBestSvd.size()-1] = 4;
                        }
                        else if(responseJetDirList[jI].find("akCs8") != std::string::npos){
                        if(centBinsStr.find("Cent0to10") != std::string::npos) histBestSvd[histBestSvd.size()-1] = 4;
                        else if(centBinsStr.find("Cent10to30") != std::string::npos) histBestSvd[histBestSvd.size()-1] = 4;
                        else if(centBinsStr.find("Cent30to50") != std::string::npos) histBestSvd[histBestSvd.size()-1] = 4;
                        else if(centBinsStr.find("Cent50to90") != std::string::npos) histBestSvd[histBestSvd.size()-1] = 4;
                        }
                        else if(responseJetDirList[jI].find("akCs6") != std::string::npos){
                        if(centBinsStr.find("Cent0to10") != std::string::npos) histBestSvd[histBestSvd.size()-1] = 4;
                        else if(centBinsStr.find("Cent10to30") != std::string::npos) histBestSvd[histBestSvd.size()-1] = 5;
                        else if(centBinsStr.find("Cent30to50") != std::string::npos) histBestSvd[histBestSvd.size()-1] = 4;
                        else if(centBinsStr.find("Cent50to90") != std::string::npos) histBestSvd[histBestSvd.size()-1] = 4;
                        }
                        else if(responseJetDirList[jI].find("akCs4") != std::string::npos){
                        if(centBinsStr.find("Cent0to10") != std::string::npos) histBestSvd[histBestSvd.size()-1] = 4;
                        else if(centBinsStr.find("Cent10to30") != std::string::npos) histBestSvd[histBestSvd.size()-1] = 5;
                        else if(centBinsStr.find("Cent30to50") != std::string::npos) histBestSvd[histBestSvd.size()-1] = 4;
                        else if(centBinsStr.find("Cent50to90") != std::string::npos) histBestSvd[histBestSvd.size()-1] = 4;
                        }
                        else if(responseJetDirList[jI].find("akCs3") != std::string::npos){
                        if(centBinsStr.find("Cent0to10") != std::string::npos) histBestSvd[histBestSvd.size()-1] = 4;
                        else if(centBinsStr.find("Cent10to30") != std::string::npos) histBestSvd[histBestSvd.size()-1] = 6;
                        else if(centBinsStr.find("Cent30to50") != std::string::npos) histBestSvd[histBestSvd.size()-1] = 4;
                        else if(centBinsStr.find("Cent50to90") != std::string::npos) histBestSvd[histBestSvd.size()-1] = 4;
                        }
                        else if(responseJetDirList[jI].find("akCs2") != std::string::npos){
                        if(centBinsStr.find("Cent0to10") != std::string::npos) histBestSvd[histBestSvd.size()-1] = 4;
                        else if(centBinsStr.find("Cent10to30") != std::string::npos) histBestSvd[histBestSvd.size()-1] = 6;
                        else if(centBinsStr.find("Cent30to50") != std::string::npos) histBestSvd[histBestSvd.size()-1] = 4;
                        else if(centBinsStr.find("Cent50to90") != std::string::npos) histBestSvd[histBestSvd.size()-1] = 4;
                        }
                        */
                     histBestSvd[histBestSvd.size()-1] = svdTerm(responseJetDirList[jI], centBinsStr);

                     leg_p->Draw("SAME");
                     saveName = "iterRatioSvd_" + tempStr + "_" + centStr + "_" + idStr[iI] + "_" + resStr + "_" + jtAbsEtaStr + "_" + tempSystStr +  "AllSvd_RecoGenAsymm_" + cleanStr + "_" + toyStr + "_" + debugStr + addedTagStr + dateStr + ".pdf";;
                     saveName = "pdfDir/" + dateStr + "/Unfold_" + tempStr + "/" + saveName;
                     quietSaveAs(canv_p, saveName);

                     delete leg_p;
                     delete canv_p;

                     canv_p = new TCanvas("canv_p", "", 450, 450);
                     leg_p = new TLegend(0.2, 0.65, 0.5, 0.95);
                     leg_p->SetBorderSize(0);
                     leg_p->SetFillColor(0);
                     leg_p->SetFillStyle(0);
                     leg_p->SetTextFont(43);
                     leg_p->SetTextSize(12);

                     if(sI != 0){
                        for(Int_t bI = 1; bI < nBayes; ++bI){
                           TH1D* tempHist_p = (TH1D*)jtPtUnfoldedBayes_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][bI]->Clone("tempHist");

                           Int_t tempPos = histTermPos[jI][cI][iI][mI][aI][0];
                           std::cout << "RATIO iterRatioSyst " << histTermPos[jI][cI][iI][mI][aI][0] << std::endl;
                           if(tempPos == -1){
                              std::cout << "POS MISSING PICK 5" << std::endl;
                              tempPos = 5;
                           }


                           std::cout << __LINE__ << std::endl;
                           ratioPrint(tempHist_p, jtPtUnfoldedBayes_RecoGenAsymm_h[jI][cI][iI][mI][aI][0][tempPos]);
                           tempHist_p->Divide(jtPtUnfoldedBayes_RecoGenAsymm_h[jI][cI][iI][mI][aI][0][tempPos]);
                           tempHist_p->SetMaximum(1.5);
                           tempHist_p->SetMinimum(0.5);
                           std::cout << __LINE__ << std::endl;

                           tempHist_p->SetMarkerStyle(styles[bI%nStyles]);
                           tempHist_p->SetMarkerSize(1);
                           tempHist_p->SetMarkerColor(vg.getColor(bI%nColors));
                           tempHist_p->SetLineColor(vg.getColor(bI%nColors));

                           jtPtUnfoldedBayes_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][bI]->SetMarkerStyle(styles[bI%nStyles]);
                           jtPtUnfoldedBayes_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][bI]->SetMarkerSize(1);
                           jtPtUnfoldedBayes_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][bI]->SetMarkerColor(vg.getColor(bI%nColors));
                           jtPtUnfoldedBayes_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][bI]->SetLineColor(vg.getColor(bI%nColors));

                           std::cout << __LINE__ << std::endl;

                           if(bI == 1) tempHist_p->DrawCopy("HIST E1 P");
                           else tempHist_p->DrawCopy("SAME HIST E1 P");

                           std::string legStr = "Bayes=" + std::to_string(bayesVal[bI]);
                           if(histTermPos[jI][cI][iI][mI][aI][sI] == bI) legStr = legStr + ",*";
                           if(histTermPos[jI][cI][iI][mI][aI][0] == bI) legStr = legStr + ",**";

                           leg_p->AddEntry(jtPtUnfoldedBayes_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][bI], legStr.c_str(), "P L");

                           std::cout << __LINE__ << std::endl;
                           delete tempHist_p;
                        }

                        leg_p->Draw("SAME");
                        std::string saveName = "iterRatioBayesSyst_" + tempStr + "_" + centStr + "_" + idStr[iI] + "_" + resStr + "_" + jtAbsEtaStr + "_" + tempSystStr +  "AllBayes_RecoGenAsymm_" + cleanStr + "_" + toyStr + "_" + debugStr + addedTagStr + dateStr + ".pdf";;
                        saveName = "pdfDir/" + dateStr + "/Unfold_" + tempStr + "/" + saveName;
                        quietSaveAs(canv_p, saveName);
                     }

                     delete leg_p;
                     delete canv_p;

                     leg_p = new TLegend(0.2, 0.65, 0.5, 0.95);
                     leg_p->SetBorderSize(0);
                     leg_p->SetFillColor(0);
                     leg_p->SetFillStyle(0);
                     leg_p->SetTextFont(43);
                     leg_p->SetTextSize(12);

                     canv_p = new TCanvas("canv_p", "", 450, 450);
                     canv_p->SetTopMargin(0.01);
                     canv_p->SetRightMargin(0.01);
                     canv_p->SetLeftMargin(0.14);
                     canv_p->SetBottomMargin(0.14);

                     Double_t newGenBins[nMaxJtPtBins+1];
                     for(Int_t gI = 0; gI < nGenJtPtBins[jI][cI]-1; ++gI){
                        newGenBins[gI] = genJtPtBins[jI][cI][gI+1];
                     }

                     TH1D* tempHist_p = new TH1D("tempHist_h", ";pT;Ratio SVD/Bayes", nGenJtPtBins[jI][cI]-2, newGenBins);
                     TH1D* tempHist2_p[nSvd];

                     Int_t nSvdTrue = TMath::Min(nSvd, nGenJtPtBins[jI][cI]);

                     Int_t divSvd = histBestSvd[histBestSvd.size()-1]-1;
                     if(divSvd < 0){
                        divSvd = TMath::Min(nSvd, nGenJtPtBins[jI][cI]) - 1;
                        histBestSvd[histBestSvd.size()-1] = TMath::Min(nSvd, nGenJtPtBins[jI][cI]) - 1;
                     }



                     const std::string tempHistTag = tempStr + "_" + centStr + "_" + idStr[iI] + "_ResponseMod" + prettyString(responseMod[mI], 2, true) + "_" + jtAbsEtaStr + "_" + systStr[sI];


                     Int_t pos = -1;

                     for(unsigned int tIX = 0; tIX < histTag.size(); ++tIX){
                        if(isStrSame(histTag[tIX], tempHistTag)){
                           pos = tIX;
                           break;
                        }
                     }

                     Int_t divBayes = histBestBayesPos[pos];
                     if(divBayes < 0) divBayes = 10;



                     macroHistToSubsetHist(jtPtUnfoldedBayes_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][divBayes], tempHist_p);

                     for(Int_t dI = 0; dI < nSvdTrue; ++dI){
                        tempHist2_p[dI] = new TH1D("tempHist2_h", ";pT;Ratio SVD/Bayes", nGenJtPtBins[jI][cI]-2, newGenBins);

                        macroHistToSubsetHist(jtPtUnfoldedSvd_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][dI], tempHist2_p[dI]);
                        tempHist2_p[dI]->Divide(tempHist_p);
                        tempHist2_p[dI]->SetMaximum(1.35);
                        tempHist2_p[dI]->SetMinimum(0.65);

                        tempHist2_p[dI]->SetMarkerStyle(styles[dI%nStyles]);
                        tempHist2_p[dI]->SetMarkerSize(1);
                        tempHist2_p[dI]->SetMarkerColor(vg.getColor(dI%nColors));
                        tempHist2_p[dI]->SetLineColor(vg.getColor(dI%nColors));

                        TF1* lineFit_p = new TF1("lineFit_p", "[0]", genJtPtBins[jI][cI][0], genJtPtBins[jI][cI][nGenJtPtBins[jI][cI]]);
                        tempHist_p->Print("ALL");
                        lineFit_p->SetParameter(0, 1.);
                        lineFit_p->SetParLimits(0, .99, 1.01);
                        tempHist2_p[dI]->Fit("lineFit_p", "M N", "", tempHist_p->GetBinLowEdge(1), tempHist_p->GetBinLowEdge(tempHist_p->GetNbinsX()+1));


                        leg_p->AddEntry(tempHist2_p[dI], ("SVD=" + std::to_string(dI+1) + "," + prettyString(lineFit_p->GetProb(), 3, false)).c_str(), "P L");

                        delete lineFit_p;

                        if(dI == 0) tempHist2_p[dI]->DrawCopy("HIST E1 P");
                        else tempHist2_p[dI]->DrawCopy("HIST E1 P SAME");
                     }

                     TLatex* label_p = new TLatex();
                     label_p->SetTextFont(43);
                     label_p->SetTextSize(14);
                     label_p->SetNDC();                  

                     label_p->DrawLatex(0.6, 0.9, ("SVD=" + std::to_string(divSvd+1)).c_str());
                     label_p->DrawLatex(0.6, 0.8, ("Bayes=" + std::to_string(bayesVal[divBayes])).c_str());

                     leg_p->Draw("SAME");

                     saveName = "ratioBayesSVD" + tempStr + "_" + centStr + "_" + idStr[iI] + "_" + resStr + "_" + jtAbsEtaStr + "_" + tempSystStr +  "AllBayes_RecoGenAsymm_" + cleanStr + "_" + toyStr + "_" + debugStr + addedTagStr + dateStr + ".pdf";;
                     saveName = "pdfDir/" + dateStr + "/Unfold_" + tempStr + "/" + saveName;
                     if(sI == 0) quietSaveAs(canv_p, saveName);

                     delete tempHist_p;
                     for(Int_t dI = 0; dI < nSvdTrue; ++dI){
                        delete tempHist2_p[dI];
                     }
                     delete canv_p;
                     delete label_p;
                     delete leg_p;

                     for(Int_t bI = 0; bI < nBayes; ++bI){		
                        jtPtUnfoldedBayes_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][bI]->Write("", TObject::kOverwrite);

                        delete jtPtUnfoldedBayesMC_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][bI];
                        delete bayes_p[jI][cI][iI][mI][aI][sI][bI];
                        delete bayesMC_p[jI][cI][iI][mI][aI][sI][bI];
                     }

                     for(Int_t bI = 0; bI < TMath::Min(nGenJtPtBins[jI][cI], nSvd); ++bI){		
                        jtPtUnfoldedSvd_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][bI]->Write("", TObject::kOverwrite);

                        delete jtPtUnfoldedSvd_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][bI];
                        delete jtPtUnfoldedSvdMC_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][bI];
                        delete svd_p[jI][cI][iI][mI][aI][sI][bI];
                        delete svdMC_p[jI][cI][iI][mI][aI][sI][bI];
                     }
                  }

                  for(Int_t sI = 0; sI < nSyst; ++sI){
                     for(Int_t bI = 0; bI < nBayes; ++bI){
                        delete jtPtUnfoldedBayes_RecoGenAsymm_h[jI][cI][iI][mI][aI][sI][bI];
                     }
                  }
               }
            }
         }
      }
   }

   /*
   //TEX FILE
   x  const std::string textWidth = "0.24";
   std::string texFileName = outFileName;
   texFileName.replace(texFileName.find(".root"), std::string(".root").size(), ".tex");
   texFileName.replace(0, std::string("output").size(), "pdfDir/" + dateStr + "/Unfold");

   std::ofstream texFile(texFileName.c_str());

   texFile << "\\RequirePackage{xspace}" << std::endl;
   texFile << "\\RequirePackage{amsmath}" << std::endl;
   texFile << std::endl;

   texFile << "\\documentclass[xcolor=dvipsnames]{beamer}" << std::endl;
   texFile << "\\usetheme{Warsaw}" << std::endl;
   texFile << "\\setbeamercolor{structure}{fg=NavyBlue!90!NavyBlue}" << std::endl;
   texFile << "\\setbeamercolor{footlinecolor}{fg=white,bg=lightgray}" << std::endl;
   texFile << std::endl;

   texFile << "\\newcommand{\\pt}{\\ensuremath{p_{\\mathrm{T}}}\\xspace}" << std::endl;
   texFile << std::endl;

   texFile << "\\setbeamersize{text margin left=3pt,text margin right=3pt}" << std::endl;
   texFile << std::endl;

   texFile << "\\setbeamerfont{frametitle}{size=\\tiny}" << std::endl;
   texFile << "\\setbeamertemplate{frametitle}" << std::endl;
   texFile << "{" << std::endl;
   texFile << "    \\nointerlineskip" << std::endl;
   texFile << "    \\begin{beamercolorbox}[sep=0.05cm, ht=1.0em, wd=\\paperwidth]{frametitle}" << std::endl;
   texFile << "        \\vbox{}\\vskip-2ex%" << std::endl;
   texFile << "        \\strut\\insertframetitle\\strut" << std::endl;
   texFile << "        \\vskip-0.8ex%" << std::endl;
   texFile << "    \\end{beamercolorbox}" << std::endl;
   texFile << "}" << std::endl;
   texFile << std::endl;

   texFile << "\\setbeamertemplate{footline}{%" << std::endl;
   texFile << "  \\begin{beamercolorbox}[sep=.4em,wd=\\paperwidth,leftskip=0.5cm,rightskip=0.5cm]{footlinecolor}" << std::endl;
   texFile << "    \\hspace{0.075cm}%" << std::endl;
   texFile << "    \\hfill\\insertauthor \\hfill\\insertpagenumber" << std::endl;
   texFile << "  \\end{beamercolorbox}%" << std::endl;
   texFile << "}" << std::endl;
   texFile << "\\setbeamertemplate{navigation symbols}{}" << std::endl;
   texFile << std::endl;

   texFile << "\\setbeamertemplate{itemize item}[circle]" << std::endl;
   texFile << "\\setbeamertemplate{itemize subitem}[circle]" << std::endl;
   texFile << "\\setbeamertemplate{itemize subsubitem}[circle]" << std::endl;
   texFile << "\\setbeamercolor{itemize item}{fg=black}" << std::endl;
   texFile << "\\setbeamercolor{itemize subitem}{fg=black}" << std::endl;
   texFile << "\\setbeamercolor{itemize subsubitem}{fg=black}" << std::endl;
   texFile << std::endl;

   texFile << "\\definecolor{links}{HTML}{00BFFF}" << std::endl;
   texFile << "\\hypersetup{colorlinks,linkcolor=,urlcolor=links}" << std::endl;
   texFile << std::endl;

   texFile << "\\author[CM]{Chris McGinn}" << std::endl;
   texFile << std::endl;

   texFile << "\\begin{document}" << std::endl;
   texFile << "\\begin{frame}" << std::endl;
   texFile << "\\frametitle{\\centerline{Unfolding termination (" << dateStr2 << ")}}" << std::endl;
   texFile << " \\begin{itemize}" << std::endl;
   texFile << "  \\fontsize{10}{10}\\selectfont" << std::endl;
   texFile << "  \\item{Placeholder}" << std::endl;
   texFile << "  \\begin{itemize}" << std::endl;
   texFile << "   \\fontsize{10}{10}\\selectfont" << std::endl;
   texFile << "   \\item{Placeholder}" << std::endl;
   texFile << "  \\end{itemize}" << std::endl;
   texFile << " \\end{itemize}" << std::endl;
   texFile << "\\end{frame}" << std::endl;


   if(doLocalDebug) std::cout << __FILE__ << ", " << "LINE: " << __LINE__ << std::endl;

   for(unsigned int spI = 0; spI < pdfNames.size(); ++spI){
      if(pdfNames.at(spI).at(0).find("AbsEta0p0to2p0") == std::string::npos) continue;

      std::string centStr = "PP";

      if(!isDataPP){
         centStr = pdfNames.at(spI).at(0).substr(pdfNames.at(spI).at(0).find("Cent"), pdfNames.at(spI).at(0).size());
         centStr.replace(centStr.find("_"), centStr.size(), "");
      }

      std::string jtStr = pdfNames.at(spI).at(0).substr(pdfNames.at(spI).at(0).find("_ak")+1, pdfNames.at(spI).at(0).size());
      jtStr.replace(jtStr.find("_"), jtStr.size(), "");

      std::string resStr = pdfNames.at(spI).at(0).substr(pdfNames.at(spI).at(0).find("_Response")+1, pdfNames.at(spI).at(0).size());
      resStr.replace(resStr.find("_"), resStr.size(), "");

      std::string idStr = "";
      if(pdfNames.at(spI).at(0).find("_NoID") != std::string::npos) idStr = "NoID";
      else if(pdfNames.at(spI).at(0).find("_LightMUID") != std::string::npos) idStr = "LightMUID";
      else if(pdfNames.at(spI).at(0).find("_LightMUAndCHID") != std::string::npos) idStr = "LightMUAndCHID";
      else if(pdfNames.at(spI).at(0).find("_FullLight") != std::string::npos) idStr = "FullLight";
      else if(pdfNames.at(spI).at(0).find("_FullTight") != std::string::npos) idStr = "FullTight";

      if(idStr.find("LightMUAndCHID") == std::string::npos) continue;
      if(resStr.find("0p10") == std::string::npos) continue;

      texFile << "\\begin{frame}" << std::endl;
      texFile << "\\frametitle{\\centerline{" << jtStr << ", " << centStr << ", " << resStr << ", " << idStr << "}}" << std::endl;

      for(unsigned int mI = 0; mI < pdfNames.at(spI).size(); ++mI){
         texFile << "\\includegraphics[width=" << textWidth << "\\textwidth]{" << pdfNames.at(spI).at(mI) << "}";
         if(mI == 3 || mI == 7) texFile << "\\\\";
         texFile << std::endl;
      }

      texFile << "\\begin{itemize}" << std::endl;
      texFile << "\\fontsize{8}{8}\\selectfont" << std::endl;
      texFile << "\\item{test}" << std::endl;
      texFile << "\\end{itemize}" << std::endl;
      texFile << "\\end{frame}" << std::endl;
   }


   texFile << "\\end{document}" << std::endl;
   texFile << std::endl;

   texFile.close();
   */

      std::cout << "Total nPDF: " << nPDFTotal << std::endl;

   for(Int_t jI = 0; jI < nDataJet; ++jI){
      texSlideCreator tex;

      std::string tempStr = responseJetDirList[jI];
      if(tempStr.find("/") != std::string::npos) tempStr.replace(tempStr.find("/"), tempStr.size() - tempStr.find("/"), "");

      tex.Clean();
      tex.Init(outFileName);
      tex.InitDir("pdfDir/" + dateStr + "/Unfold_" + tempStr + "/");
      tex.SetAuthor("Chris McGinn");

      tex.SetSlideTitles(slideTitlesPerAlgo[jI]);
      tex.SetSlidePdfs(pdfPerSlidePerAlgo[jI]);
      if(!(tex.CreateTexSlides())){
         std::cout << "Warning: .tex slide creation failed" << std::endl;
      }
   }

   outFile_p->cd();
   TDirectory* cutDir_p = (TDirectory*)outFile_p->mkdir("cutDir");
   TDirectory* subDir_p = (TDirectory*)cutDir_p->mkdir("subDir");
   TDirectory* unfoldDirBayes_p = (TDirectory*)cutDir_p->mkdir("unfoldDirBayes");
   TDirectory* unfoldDirSVD_p = (TDirectory*)cutDir_p->mkdir("unfoldDirSVD");

   if(nHistDim != (int)histTag.size()) std::cout << "WARNING: nHistDim (" << nHistDim << ") != histTag.size() (" << histTag.size() << ")" << std::endl;
   if(nHistDim != (int)histBestBayes.size()) std::cout << "WARNING: nHistDim (" << nHistDim << ") != histBestBayes.size() (" << histBestBayes.size() << ")" << std::endl;
   if(nHistDim != (int)histBestSvd.size()) std::cout << "WARNING: nHistDim (" << nHistDim << ") != histBestSvd.size() (" << histBestSvd.size() << ")" << std::endl;

   cutPropData.SetNHistDim(nHistDim);
   cutPropData.SetHistTagBayes(histTag);
   cutPropData.SetHistBestBayes(histBestBayes);

   cutPropData.SetHistTagSVD(histTag);
   cutPropData.SetHistBestSVD(histBestSvd);

   if(!cutPropData.WriteAllVarToFile(outFile_p, cutDir_p, subDir_p, unfoldDirBayes_p, unfoldDirSVD_p)) std::cout << "Warning: Cut writing has failed" << std::endl;

   outFile_p->Close();
   delete outFile_p;

   totalRunWatch.stop();

   std::cout << "TOTAL RUN: " << totalRunWatch.totalWall() << std::endl;
   std::cout << " TOTAL UNFOLD: " << doUnfoldWatch.totalWall() << std::endl;
   std::cout << " TOTAL UNFOLDS: " << unfoldCounter << std::endl;

   delete randGen_p;

   return 0;
}

int main(int argc, char* argv[])
{
   if(argc < 3 || argc > 7){
      std::cout << "Usage: ./bin/unfoldRawData.exe <inDataFileName> <inResponseName> <selectJtAlgo-Opt> <overrideFile-Opt> <doClean-Opt> <nToys-Opt>" << std::endl;
      return 1;
   }

   int retVal = 0;
   if(argc == 3) retVal += unfoldRawData(argv[1], argv[2]);
   else if(argc == 4) retVal += unfoldRawData(argv[1], argv[2], argv[3]);
   else if(argc == 5) retVal += unfoldRawData(argv[1], argv[2], argv[3], argv[4]);
   else if(argc == 6) retVal += unfoldRawData(argv[1], argv[2], argv[3], argv[4], std::stoi(argv[5]));
   else if(argc == 7) retVal += unfoldRawData(argv[1], argv[2], argv[3], argv[4], std::stoi(argv[5]), std::stoi(argv[6]));

   std::cout << "Job complete. Return " << retVal << "." << std::endl;
   return retVal;
}
