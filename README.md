To run:

source setFullJREnv.sh

This will setup some necessary environmental variables for each directory to build. Note that directories, while interdependent in code, are built independently. So if one wants only to run or modify a specific subset of the code, this is done by working within the desired directory.

If source setFullJREnv.sh gives empty paths (likely if not working on CERN computing resources), fill in the appropriate path and source.

To build a directory, cd into directory and run 'make'
To clean, cd into directory and run 'make clean'

NOTE: Check carefully that your ROOT matches your RooUnfold

To run mainline analysis (Until otherwise stated, assuming you are working under MainAnalysis):

Preamble (Necessary files):

Need four sets of forest files:
 * PbPb Data forest
 * pp Data Forest
 * PbPb MC files (some selection of pthats)
 * pp MC files (some selection of pthats)

Step 1 (making input txt files):

Since the MC files will most likely be a set of pthats, dump them in a text file w/ the following format:

ISPP=<ANSWER 1 or 0>
PTHATS=<Comma separated list of pthat for each sample, preferably ordered>
PTHATWEIGHTS=<Comma separated list of pthatWeights corresponding to each value - if not available just set all to one, we will derive lates>
/LINE/WITH/FULL/FILE/PATH/AND/NAME.root
/LINE/WITH/FULL/FILE/PATH/AND/NAME.root
etc..

Step 2 (pthatweights):

If we dont have pthat weights for our MC, we start with derivation of pthatWeights

./bin/deriveDijetWeights.exe <inTxtFileName> <isPYTHIA6 boolean>

where the first argument is the txt file you derived in Step 1. The second argument is for whether you use PYTHIA6 or PYTHIA8 (most likely you are using PYTHIA6 but double check in any case)

This job will spit out pthatWeights, renormalized such that the lowest pthat has a weight of 1 (i.e. identical to counts, convenient for debug), for smooth pthat distribution. Insert these weights in your txt file from step one. They will be used throughout.

Step 3 (Flat weights):

If we want to derive weights for a flat spectra for use during the response matrix construction phase, we have to create a table.

With our text files w/ proper pthatWeights, do

./bin/makeDeriveFlatResponse.exe <inTxtName> <isPP-opt>

We dont have to add the last argument, as the isPP bool is already in the text file.

This will produce a file in the output directory called '<inTxtName>_FlatGenJetResponse_DATEPRODUCED_HOURPRODUCED.root'

mv this file into tables, and in the src/makeJetResponseTree.C find the flatTableReader initiliazation and insert this new file name w/ path starting from 'MainAnalysis', i.e. 'MainAnalysis/tables/<name>'

Step 4 (Building response matrices):

We can now derive response matrices from the MC if we choose. It will probably be more convenient to split the files into single jet algos and process those in parallel. This can be done w/ Repository here:

https://github.com/cfmcginn/GeneralTree

When splitting into single algos, make sure you also keep the R=0.4 jet tree in PbPb as this is used to exclude certain events where hydjet produces jets on ~300 GeV or above (no need to confused the response matrices w/ these rare events)

Either with your split set of trees w/ redefined text files or with the full set of jet algos in a single file, now do:

./bin/makeJetResponseTree.exe <inTxtName> <isPP-opt> <fractionalNumberOfEvents>

where the second argument is optional and handled in the text file input (but probably specify anyways), and the third is an option to do a reduced fraction processing for debuggging that does not require a change to the pthat weights.

This job will produce a response matrix file in outputs with name 'output/<inTxtName>_FracNEntries<fractionalNumberOfEvents>_JetResponse_TODAYSDATE.root'

If you choose to split processing by algo, keep the matrices separate for now. We will do recombination after parallelized unfolding.

Step 5 (Processing Data):

Now we need to process the raw data w/ the following command

./bin/processRawData.exe <inDataFileName> <inResponseName> <isPP-opt> <tagStr-opt>

If you opted to split the response files to run in parallel, I recommend you stick to that splitting now, and give the algo of each job as <tagStr-opt>. Here is an example:

./bin/processRawData.exe /data/cmcginn/Forests/pp2015Data/HIHardProbes/HiForestAOD_HighPtJet80_HLTJet80_LargeROR_PtCut110_AbsEta5_20180627_11Lumi_180627_122059_355_OutOf355_MERGED.root output/Pythia6_Dijet_pp502_MCDijet_20180712_ExcludeTop4_ExcludeToFrac_Frac0p7_Full_5Sigma_20180712_SVM_ak10PFJetAnalyzer_FracNEntries1p00_JetResponse_20180827.root 1 ak10PFJetAnalyzer >& logs/processPP_ak10PFJetAnalyzer.log &

This produces output file 'output/$DATE/HiForestAOD_HighPtJet80_HLTJet80_LargeRO_Pythia6_Dijet_pp502_MCDijet_20180712_Exc_ProcessRawData_ak10PFJetAnalyzer_20180827.root'

Step 6 (Unfolding):
Now that the raw data is processed, we run unfolding of pp and pbpb jtalgos in parallel
./bin/unfoldRawData.exe <inDataFileName> <inResponseName> <selectJtAlgo-Opt>
where inDataFileName is the product of step5, inResponseName is the same file w/ response matrices, and selectJtAlgo useful if you are doing all jet algos in separate files

here is an example command
./bin/unfoldRawData.exe output/HiForestAOD_HighPtJet80_HLTJet80_LargeRO_Pythia6_Dijet_pp502_MCDijet_20180712_Exc_ProcessRawData_ak10PFJetAnalyzer_20180827.root output/Pythia6_Dijet_pp502_MCDijet_20180712_ExcludeTop4_ExcludeToFrac_Frac0p7_Full_5Sigma_20180712_SVM_ak10PFJetAnalyzer_FracNEntries1p00_JetResponse_20180827.root ak10PF >& logs/unfoldPP_ak10PF.log &

which will produce the output file 'output/$DATE/HiForestAOD_HighPtJet80_HLTJet80_LargeRO_Pythia6_Dijet_pp502_MCDijet_20180712_Exc_UnfoldRawData_NSuperBayes1_ak10PF_20180827.root'

NSuperBayes is currently experimental and hardcoded

Step 7 (combining):

At this point it is highly convenient to combine files into single for PbPb and pp before plotting. Since TNamed is not well handled by hadd, a custom script is in place for this. Do:
./bin/combineFiles.exe <outFileName> <Long> <List> <Of> <Input> <Files> 

Step 8 (start plotting):
To do spectra plotting do:
./bin/plotUnfoldedAll.exe <ppFile> <pbpbFile> <HEP data for ATLAS comp> <tagStr Use CombinedAlgos> <tables/overrideBinsPlot.txt - see tables dir, allows for arbitrary rebins> 

Using bash scripts:
bash bash/runResponse.sh # Step 4, long
bash bash/runProcess.sh # Step 5, quick
bash bash/runUnfold.sh # Step 6, quick-to-long depending on how much unfolding is done
bash bash/runCombine.sh # Step 7, quick
bash bash/runPlotUnfold.sh # Step 8, quick