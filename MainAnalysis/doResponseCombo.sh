#!/bin/bash

files=(output/Pythia6_Dijet_pp502_MCDijet_20180712_ExcludeTop4_ExcludeToFrac_Frac0p7_Full_5Sigma_20180712_SVM_ak10PFJetAnalyzer_FracNEntries1p00_JetResponse_20180823.root output/Pythia6_Dijet_pp502_MCDijet_20180712_ExcludeTop4_ExcludeToFrac_Frac0p7_Full_5Sigma_20180712_SVM_ak3PFJetAnalyzer_FracNEntries1p00_JetResponse_20180823.root output/Pythia6_Dijet_pp502_MCDijet_20180712_ExcludeTop4_ExcludeToFrac_Frac0p7_Full_5Sigma_20180712_SVM_ak4PFJetAnalyzer_FracNEntries1p00_JetResponse_20180823.root output/Pythia6_Dijet_pp502_MCDijet_20180712_ExcludeTop4_ExcludeToFrac_Frac0p7_Full_5Sigma_20180712_SVM_ak6PFJetAnalyzer_FracNEntries1p00_JetResponse_20180823.root output/Pythia6_Dijet_pp502_MCDijet_20180712_ExcludeTop4_ExcludeToFrac_Frac0p7_Full_5Sigma_20180712_SVM_ak8PFJetAnalyzer_FracNEntries1p00_JetResponse_20180823.root)

comboStr=""

for i in ${files[@]}
do
    comboStr="$comboStr $i"
done

#echo $comboStr

./bin/combineFiles.exe output/Pythia6_Dijet_pp502_MCDijet_20180712_ExcludeTop4_ExcludeToFrac_Frac0p7_Full_5Sigma_20180712_SVM_AllAlgos_FracNEntries1p00_JetResponse_20180823.root $comboStr 