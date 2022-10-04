#include <TSystem>
//#include "RooFit.h"

class StMaker;
class StChain;
class StPicoDstMaker;
class StPicoEvent;
class StRefMultCorr;
class StEpdGeom;
class StEpdEpInfo;
class StEpdEpFinder;


//StChain *chain;

void readPicoDst(const Char_t *inputFile="test.list", int jobindex, int run=11, float energy=200.0, Char_t *ListDir="./datalist/")
{
	int iJob = jobindex;
	cout << "current job id: " << iJob << endl;

	TString StrOutName = Form("output_%.6d.root", iJob);
	const Char_t *outputFile = StrOutName.Data();

	Long64_t nEvents = 1000000000;
	//nEvents = 50000;

	//Load all the System libraries
	gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
	loadSharedLibraries();

 	gROOT->LoadMacro("lMuDst.C");
	TString outputKFParticleQA = Form("KFParticleQA_%.6d.root", iJob);
	lMuDst(-1,inputFile,"ry2016,RpicoDst,kfpInter,mysql,quiet,nodefault",outputKFParticleQA.Data());

	//gSystem->AddIncludePath("-lRooFitCore");
	//gSystem->Load("libRooFit");
	gSystem->Load("StUtilities");
	gSystem->Load("StRefMultCorr");
	gSystem->Load("StPicoEvent");
	gSystem->Load("StPicoDstMaker");
	//gSystem->Load("StEpdUtil");
	gSystem->Load("StEpdUtil.so");
	//gSystem->Load("StdEdxPull_cxx.so");
	//gSystem->Load("StLambdaDecayPair");
	//gSystem->Load("StKFParticleAnalysisMaker");
	//gSystem->Load("KFParticle");

	StKFParticleAnalysisMaker *anaMaker = new StKFParticleAnalysisMaker("ana", outputFile);
	anaMaker->setRunEnergyAndListDir(run,energy,ListDir);            

	// Loop over the links in the chain
	Int_t iInit = chain -> Init() ;
	if (iInit) chain->Fatal(iInit,"on init");
	cout<<"chain->Init();"<<endl;

	// KFParticle
	StKFParticleInterface::instance()->CleanLowPVTrackEvents();
	//StKFParticleInterface::instance()->UseHFTTracksOnly();
	StKFParticleInterface::instance()->SetSoftKaonPIDMode();
	StKFParticleInterface::instance()->SetSoftTofPidMode();
	StKFParticleInterface::instance()->SetChiPrimaryCut(10);
	StKFParticleInterface::instance()->SetChiPrimaryCut2D(3); // default >3
	StKFParticleInterface::instance()->SetChi2Cut2D(10);      // default <10
	StKFParticleInterface::instance()->SetLCut(1.0); // default >5.0
	// StKFParticleInterface::instance()->SetLdLCut2D(3); // default >5.0, not sure where it is used
	// StKFParticleInterface::instance()->SetMaxDistanceBetweenParticlesCut(1.5); // default <1
	// StKFParticleInterface::instance()->SetLdLCutXiOmega(3); // default >10
	// StKFParticleInterface::instance()->SetChi2TopoCutXiOmega(10); // default <5
	// StKFParticleInterface::instance()->SetChi2CutXiOmega(10); // default <6

	//Add decays to the reconstruction list
	StKFParticleInterface::instance()->AddDecayToReconstructionList( 3122);
	StKFParticleInterface::instance()->AddDecayToReconstructionList(-3122);
	StKFParticleInterface::instance()->AddDecayToReconstructionList( 3334);
	StKFParticleInterface::instance()->AddDecayToReconstructionList(-3334);
	//StKFParticleInterface::instance()->AddDecayToReconstructionList( 333); // test for Ding

	// StPicoDstMaker & chain
	StPicoDstMaker* maker = (StPicoDstMaker *) StMaker::GetTopChain()->Maker("PicoDst");
	if (! maker) return;
	maker->SetStatus("*",1);
	TChain *tree = maker->chain();
	Long64_t nentries = tree->GetEntries();
	if (nentries <= 0) return;
	Long64_t EventsToRun = TMath::Min(nEvents,nentries);
	cout << nentries << " events in chain " << EventsToRun << " will be read." << endl;

	time_t time_start;
	time_t time_now;
	time(&time_start);

	Int_t istat = 0, i = 1;
	while(i<=EventsToRun && istat!=2) {
		//if(i%5000==0){cout << endl; cout << "== Event " << i << " start ==" << endl;}
		chain->Clear();
		istat = chain->Make(i);
		if(i%100==0) {
			cout << endl; 
			cout << "== Event " << i << " finish == " << flush;
			time(&time_now);
			int time_diff = (int)difftime(time_now, time_start);
			cout << time_diff/60 << "min " << time_diff%60 << "s: " << 1.0*time_diff/i << "s/event" << endl;
		}
		if (istat == 2)
			cout << "Last  event processed. Status = " << istat << endl;
		if (istat == 3)
			cout << "Error event processed. Status = " << istat << endl;
		i++;
	}


	cout << "****************************************** " << endl;
	cout << "Work done... now its time to close up shop!"<< endl;
	cout << "****************************************** " << endl;
	chain->Finish();
	cout << "****************************************** " << endl;
	cout << "total number of events  " << EventsToRun << endl;
	cout << "****************************************** " << endl;

	delete chain;


}
