CXX = g++
#O3 for max optimization (go to 0 for debug)
CXXFLAGS = -Wall -Werror -O0 -Wextra -Wno-unused-local-typedefs -Wno-deprecated-declarations -std=c++11 -g
ifeq "$(GCCVERSION)" "1"
  CXXFLAGS += -Wno-error=misleading-indentation
endif

ifndef FULLJRDIR
$(error FULLJRDIR is not set at all. Do 'source setFullJREnv.sh')
endif

INCLUDE=-I $(FULLJRDIR)
ROOT=`root-config --cflags --glibs`

MKDIR_BIN=mkdir -p $(PWD)/bin
MKDIR_OUTPUT=mkdir -p $(PWD)/output
MKDIR_PDF=mkdir -p $(PWD)/pdfDir

CHECK_FULLJRDIR=printf "\nFULLJRDIR is set to: \'$(FULLJRDIR)\'\n   !!!IF FULLJRDIR IS EMPTY THE BUILD WILL FAIL!!!\n   Run 'source \
setStudyMultEnv.sh' in top level StudyMult directory to fix.\n\n"

ROOUNF=-I $(ROOUNFPATH) -L $(ROOUNFPATH) -lRooUnfold

#all: mkdirBin mkdirPdf mkdirOutput checkFullJRDir bin/makeRawRAATree.exe bin/makeJetResponseTree.exe bin/makeJetResponseTreeUpdated.exe bin/makeDeriveFlatResponse.exe bin/validateJetResponse.exe bin/plotJetResponse.exe bin/makeMuonFakeCheck.exe bin/checkGoodJetBadJetPF.exe bin/processRawData.exe bin/unfoldRawData.exe bin/unfoldRawDataTest.exe bin/deriveDijetWeights.exe bin/checkGenHydSpectra.exe bin/plotUnfoldedSpectra.exe bin/plotUnfoldedAll.exe bin/checkUnfoldTermination.exe bin/combineFiles.exe bin/testRooUnfResUnderflow.exe bin/checkFileNClassContents.exe bin/compareAllHist.exe bin/lumiComp.exe bin/extractNMB.exe bin/checkOverlap.exe bin/comboTest.exe bin/combineResponse.exe bin/splitFiles.exe bin/unfoldTest.exe bin/validateSplit.exe bin/extractBins.exe bin/extractNSigma.exe bin/extractSigmaFit.exe bin/quickHistSum.exe #bin/plotUnfoldedRAA.exe

all: mkdirBin mkdirPdf mkdirOutput checkFullJRDir  bin/splitFiles.exe bin/makeRawRAATree.exe bin/makeJetResponseTree.exe bin/makeJetResponseTreeUpdated.exe  bin/makeJetResponseTreeUpdated2.exe bin/makeJetResponseTreeUpdated3.exe bin/makeDeriveFlatResponse.exe bin/validateJetResponse.exe bin/processRawData.exe bin/unfoldRawData.exe bin/unfoldRawDataTest.exe bin/deriveDijetWeights.exe bin/checkGenHydSpectra.exe bin/plotUnfoldedAll.exe bin/checkUnfoldTermination.exe bin/combineFiles.exe bin/testRooUnfResUnderflow.exe bin/checkFileNClassContents.exe bin/compareAllHist.exe bin/lumiComp.exe bin/extractNMB.exe bin/checkOverlap.exe bin/comboTest.exe bin/combineResponse.exe bin/unfoldTest.exe bin/validateSplit.exe bin/extractBins.exe bin/extractNSigma.exe bin/extractSigmaFit.exe bin/quickHistSum.exe bin/toyATLASCMS.exe bin/toyATLASCMSAlt.exe bin/jetCorr.exe bin/testSVD.exe bin/counts.exe bin/jetJSON.exe #bin/plotUnfoldedRAA.exe

mkdirBin:
	$(MKDIR_BIN)

mkdirOutput:
	$(MKDIR_OUTPUT)

mkdirPdf:
	$(MKDIR_PDF)

checkFullJRDir:
	@$(CHECK_FULLJRDIR)

bin/makeRawRAATree.exe: src/makeRawRAATree.C
	$(CXX) $(CXXFLAGS) src/makeRawRAATree.C $(ROOT) $(INCLUDE) -o bin/makeRawRAATree.exe

bin/makeJetResponseTree.exe: src/makeJetResponseTree.C
	$(CXX) $(CXXFLAGS) src/makeJetResponseTree.C $(ROOT) $(ROOUNF) $(INCLUDE) -o bin/makeJetResponseTree.exe

bin/makeJetResponseTreeUpdated.exe: src/makeJetResponseTreeUpdated.C
	$(CXX) $(CXXFLAGS) src/makeJetResponseTreeUpdated.C $(ROOT) $(ROOUNF) $(INCLUDE) -o bin/makeJetResponseTreeUpdated.exe

bin/makeJetResponseTreeUpdated2.exe: src/makeJetResponseTreeUpdated2.C
	$(CXX) $(CXXFLAGS) src/makeJetResponseTreeUpdated2.C $(ROOT) $(ROOUNF) $(INCLUDE) -o bin/makeJetResponseTreeUpdated2.exe

bin/makeJetResponseTreeUpdated3.exe: src/makeJetResponseTreeUpdated3.C
	$(CXX) $(CXXFLAGS) src/makeJetResponseTreeUpdated3.C $(ROOT) $(ROOUNF) $(INCLUDE) -o bin/makeJetResponseTreeUpdated3.exe

bin/makeDeriveFlatResponse.exe: src/makeDeriveFlatResponse.C
	$(CXX) $(CXXFLAGS) src/makeDeriveFlatResponse.C $(ROOT) $(INCLUDE) -o bin/makeDeriveFlatResponse.exe

bin/validateJetResponse.exe: src/validateJetResponse.C
	$(CXX) $(CXXFLAGS) src/validateJetResponse.C $(ROOT) $(ROOUNF) $(INCLUDE) -o bin/validateJetResponse.exe

bin/plotJetResponse.exe: src/plotJetResponse.C
	$(CXX) $(CXXFLAGS) src/plotJetResponse.C $(ROOT) $(INCLUDE) -o bin/plotJetResponse.exe

bin/plotUnfoldedSpectra.exe: src/plotUnfoldedSpectra.C
	$(CXX) $(CXXFLAGS) src/plotUnfoldedSpectra.C $(ROOT) $(INCLUDE) -o bin/plotUnfoldedSpectra.exe

bin/plotUnfoldedAll.exe: src/plotUnfoldedAll.C
	$(CXX) $(CXXFLAGS) src/plotUnfoldedAll.C $(ROOT) $(INCLUDE) -o bin/plotUnfoldedAll.exe

bin/plotUnfoldedRAA.exe: src/plotUnfoldedRAA.C
	$(CXX) $(CXXFLAGS) src/plotUnfoldedRAA.C $(ROOT) $(INCLUDE) -o bin/plotUnfoldedRAA.exe

bin/makeMuonFakeCheck.exe: src/makeMuonFakeCheck.C
	$(CXX) $(CXXFLAGS) src/makeMuonFakeCheck.C $(ROOT) $(INCLUDE) -o bin/makeMuonFakeCheck.exe

bin/checkGoodJetBadJetPF.exe: src/checkGoodJetBadJetPF.C
	$(CXX) $(CXXFLAGS) src/checkGoodJetBadJetPF.C $(ROOT) $(INCLUDE) -o bin/checkGoodJetBadJetPF.exe

bin/checkGenHydSpectra.exe: src/checkGenHydSpectra.C
	$(CXX) $(CXXFLAGS) src/checkGenHydSpectra.C $(ROOT) $(INCLUDE) -o bin/checkGenHydSpectra.exe

bin/processRawData.exe: src/processRawData.C
	$(CXX) $(CXXFLAGS) src/processRawData.C $(ROOT) $(INCLUDE) -o bin/processRawData.exe

bin/unfoldRawData.exe: src/unfoldRawData.C
	$(CXX) $(CXXFLAGS) src/unfoldRawData.C $(ROOT) $(ROOUNF) $(INCLUDE) -o bin/unfoldRawData.exe

bin/unfoldRawDataTest.exe: src/unfoldRawDataTest.C
	$(CXX) $(CXXFLAGS) src/unfoldRawDataTest.C $(ROOT) $(ROOUNF) $(INCLUDE) -o bin/unfoldRawDataTest.exe

bin/deriveDijetWeights.exe: src/deriveDijetWeights.C
	$(CXX) $(CXXFLAGS) src/deriveDijetWeights.C $(ROOT) $(INCLUDE) -o bin/deriveDijetWeights.exe

bin/checkUnfoldTermination.exe: src/checkUnfoldTermination.C
	$(CXX) $(CXXFLAGS) src/checkUnfoldTermination.C $(ROOT) $(INCLUDE) -o bin/checkUnfoldTermination.exe

bin/combineFiles.exe: src/combineFiles.C
	$(CXX) $(CXXFLAGS) src/combineFiles.C $(ROOT) $(ROOUNF) $(INCLUDE) -o bin/combineFiles.exe

bin/testRooUnfResUnderflow.exe: src/testRooUnfResUnderflow.C
	$(CXX) $(CXXFLAGS) src/testRooUnfResUnderflow.C $(ROOT) $(ROOUNF) $(INCLUDE) -o bin/testRooUnfResUnderflow.exe

bin/checkFileNClassContents.exe: src/checkFileNClassContents.C
	$(CXX) $(CXXFLAGS) src/checkFileNClassContents.C $(ROOT) $(INCLUDE) -o bin/checkFileNClassContents.exe

bin/compareAllHist.exe: src/compareAllHist.C
	$(CXX) $(CXXFLAGS) src/compareAllHist.C $(ROOT) $(INCLUDE) -o bin/compareAllHist.exe

bin/lumiComp.exe: src/lumiComp.C
	$(CXX) $(CXXFLAGS) src/lumiComp.C $(ROOT) $(INCLUDE) -o bin/lumiComp.exe

bin/extractNMB.exe: src/extractNMB.C
	$(CXX) $(CXXFLAGS) src/extractNMB.C $(ROOT) $(INCLUDE) -o bin/extractNMB.exe

bin/checkOverlap.exe: src/checkOverlap.C
	$(CXX) $(CXXFLAGS) src/checkOverlap.C $(ROOT) $(INCLUDE) -o bin/checkOverlap.exe

bin/comboTest.exe: src/comboTest.C
	$(CXX) $(CXXFLAGS) src/comboTest.C $(ROOT) $(INCLUDE) -o bin/comboTest.exe

bin/combineResponse.exe: src/combineResponse.C
	$(CXX) $(CXXFLAGS) src/combineResponse.C $(ROOT) $(INCLUDE) -o bin/combineResponse.exe

bin/splitFiles.exe: src/splitFiles.C
	$(CXX) $(CXXFLAGS) src/splitFiles.C $(ROOT) $(INCLUDE) -o bin/splitFiles.exe

bin/unfoldTest.exe: src/unfoldTest.C
	$(CXX) $(CXXFLAGS) src/unfoldTest.C $(ROOT) $(ROOUNF) $(INCLUDE) -o bin/unfoldTest.exe

bin/validateSplit.exe: src/validateSplit.C
	$(CXX) $(CXXFLAGS) src/validateSplit.C $(ROOT) $(INCLUDE) -o bin/validateSplit.exe

bin/extractBins.exe: src/extractBins.C
	$(CXX) $(CXXFLAGS) src/extractBins.C $(ROOT) $(INCLUDE) -o bin/extractBins.exe

bin/extractNSigma.exe: src/extractNSigma.C
	$(CXX) $(CXXFLAGS) src/extractNSigma.C $(ROOT) $(INCLUDE) -o bin/extractNSigma.exe

bin/extractSigmaFit.exe: src/extractSigmaFit.C
	$(CXX) $(CXXFLAGS) src/extractSigmaFit.C $(ROOT) $(INCLUDE) -o bin/extractSigmaFit.exe

bin/quickHistSum.exe: src/quickHistSum.C
	$(CXX) $(CXXFLAGS) src/quickHistSum.C $(ROOT) $(INCLUDE) -o bin/quickHistSum.exe

bin/toyATLASCMS.exe: src/toyATLASCMS.C
	$(CXX) $(CXXFLAGS) src/toyATLASCMS.C $(ROOT) $(INCLUDE) $(ROOUNF) -o bin/toyATLASCMS.exe

bin/jetCorr.exe: src/jetCorr.C
	$(CXX) $(CXXFLAGS) src/jetCorr.C $(ROOT) $(INCLUDE) $(ROOUNF) -o bin/jetCorr.exe

bin/toyATLASCMSAlt.exe: src/toyATLASCMSAlt.C
	$(CXX) $(CXXFLAGS) src/toyATLASCMSAlt.C $(ROOT) $(INCLUDE) $(ROOUNF) -o bin/toyATLASCMSAlt.exe

bin/testSVD.exe: src/testSVD.C
	$(CXX) $(CXXFLAGS) src/testSVD.C $(ROOT) $(INCLUDE) $(ROOUNF) -o bin/testSVD.exe

bin/counts.exe: src/counts.C
	$(CXX) $(CXXFLAGS) src/counts.C $(ROOT) $(INCLUDE) $(ROOUNF) -o bin/counts.exe

bin/jetJSON.exe: src/jetJSON.C
	$(CXX) $(CXXFLAGS) src/jetJSON.C $(ROOT) $(INCLUDE) -o bin/jetJSON.exe


clean:
	rm -f ./*~
	rm -f ./#*#
	rm -f bash/*~
	rm -f bash/#*#
	rm -f bin/*.exe
	rm -rf bin
	rm -f configs/*~
	rm -f configs/#*#
	rm -f include/*~
	rm -f include/#*#
	rm -f src/*~
	rm -f src/#*#
