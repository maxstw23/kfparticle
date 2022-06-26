#include "StKFParticleAnalysisMaker.h"
#include "PhysicalConstants.h"
#include "phys_constants.h"
#include "StBTofUtil/tofPathLength.hh"
#include "StThreeVectorF.hh"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TFile.h"
#include "TProfile.h"
#include "TLorentzVector.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
#include "StMessMgr.h"
#include <algorithm>
#include <TMath.h>

#include "KFVertex.h"
#include "KFParticle.h"
#include "KFParticleSIMD.h"
#include "KFPTrack.h"
#include "KFParticleTopoReconstructor.h"
#include "StKFParticleInterface.h"
#include "StKFParticlePerformanceInterface.h"

#include "StTrackHelix.h"
#include "StLambdaDecayPair.h"
#include "MyToolkit.h"

#define pi                 TMath::Pi()
#define LambdaPdgMass      1.11568
#define ProtonPdgMass      0.938272
#define PionPdgMass        0.139570
#define KaonPdgMass		   0.493677
#define LambdaPdg          3122
#define OmegaPdg           3334
#define KaonPdg			   321
#define ProtonPdg          2212
#define PionPdg           -211


//-----------------------------------------------------------------------------
ClassImp(StKFParticleAnalysisMaker)

//-----------------------------------------------------------------------------
StKFParticleAnalysisMaker::StKFParticleAnalysisMaker(const char* name, const char* outName)
: StMaker(name)
{
	mOutName = outName;

	mRun    =0;
	mEnergy =0;
	mListDir="./";
}

//----------------------------------------------------------------------------- 
StKFParticleAnalysisMaker::~StKFParticleAnalysisMaker() {
	SafeDelete(KFParticleInterface);
	SafeDelete(KFParticlePerformanceInterface);
}

//----------------------------------------------------------------------------- 
Int_t StKFParticleAnalysisMaker::openFile()
{    
	// ======= StRefMultCorr ======= //
	// set up StRefMultCorr for centrality
	mRefMultCorr = CentralityMaker::instance()->getRefMultCorr() ;
	// ======= StRefMultCorr end ======= //

	cout << "-----------------------------------" << endl;
	cout << "------- user's files loaded -------" << endl;
	cout << "-----------------------------------" << endl;

	return kStOK;
}

//----------------------------------------------------------------------------- 
Int_t StKFParticleAnalysisMaker::Init() {

        cout << mOutName << endl;
        const char *sDir = mOutName.Data();
        int lDir = mOutName.Length();
        int iDir = lDir-6;
        while(!std::isdigit(sDir[iDir])) iDir --;
        while(std::isdigit(sDir[iDir]))  iDir --;
        mJob = std::atoi(&sDir[iDir+1]);
        cout << "current job id: " << mJob << endl;

	PI = M_PI;
	twoPI = 2*M_PI;

	badList.clear();
	runList.clear();

	if(!readRunList())return kStFatal;
	if(!readBadList())return kStFatal;

	DeclareHistograms();
	Int_t openFileStatus = openFile();
	if(openFileStatus == kStFatal) return kStFatal;

	TFile *f = GetTFile(); // These two lines need to be HERE (though I don't know /why/)- don't throw in another function
	if(f){f->cd(); BookVertexPlots();}

	return kStOK;
}

//----------------------------------------------------------------------------- 
Int_t StKFParticleAnalysisMaker::Finish() {
	if(mOutName!="") {
		TFile *fout = new TFile(mOutName.Data(),"RECREATE");
		fout->cd();
		WriteHistograms();
		fout->Close();
	}
	return kStOK;
}

//-----------------------------------------------------------------------------
void StKFParticleAnalysisMaker::DeclareHistograms() {

	const int nRuns=runList.size();
	int mRunL=runList.at(0);
	int mRunH=runList.at(nRuns-1);
	int mDayL=(int) ((mRunL)/1000%1000);
	int mDayH=(int) ((mRunH)/1000%1000)+1;
	const int nDays=mDayH-mDayL;
	int nTmps=nDays;//>40?40:nDays;
	mStps=ceil((nDays*1.0)/(nTmps*1.0));
	const int nTims=ceil((nDays*1.0)/(mStps*1.0));

	cout<<nDays<<" "<<mDayL<<" "<<mDayH<<" "<<nTims<<" "<<mStps<<endl;

	hNRefMult = new TH1F("RefMult" , "Reference Multiplicity" , 1000, 0.0, 1000.0 ) ;
	hNRefMultA= new TH1F("RefMultA", "Reference MultiplicityA", 1000, 0.0, 1000.0 ) ;
	hNRefMultB= new TH1F("RefMultB", "Reference MultiplicityB", 1000, 0.0, 1000.0 ) ;
	hVertexXY = new TH2F("VertexXY", "Vertex XY Position", 200, -10.0, 10.0 , 200, -10., 10 ) ;
	hVertexZ  = new TH1F("VertexZ" , "Event Vertex Z Position", 200, -100.0, 100.0 ) ;
	hVertex2D = new TH2F("Vertex2D", "VertexZ vs VPD Vz", 200, -100.0, 100.0 , 200, -100., 100 ) ;
	hDiffVz   = new TH1F("VertexZdiff" , "VertexZ-VPDVz diff", 100, -10.0, 10.0 ) ;
	hcent     = new TH1F("centrality","centrality"  ,nCent,0.,nCent);
	hcentw    = new TH1F("centralityw","centralityw",nCent,0.,nCent);

	hcentRefM = new TProfile("hcentRefM","hcentRefM",nCent,0.,nCent,0,1000);
	hcentRefW = new TProfile("hcentRefW","hcentRefW",nCent,0.,nCent,0,1000);

	// xiatong's global kaon QA
	hgKpluspdEdx       = new TH2D("hgKpluspdEdx", "Global K+ dE/dx vs p", 5000, 0., 50., 1000, 0., 10.);
	hgKplusdEdxErr     = new TH1D("hgKplusdEdxErr", "Global K+ dE/dx error", 100, 0., 1.); // may need to adjust
	hgKplusnsigma      = new TH1D("hgKplusnsigma", "Global K+ n_{#sigma}", 600, -6., 6.);
	hgKpluspinvbeta    = new TH2D("hgKpluspinvbeta", "Global K+ inverse beta vs p", 5000, 0., 50., 1000, 0., 10.);
	hgKplusm2          = new TH1D("hgKplusm2", "Global K+ m^2", 500, -1., 4.);
	hgKplusbtofYlocal  = new TH1D("hgKplusbtofYlocal", "Global K+ BTOF Ylocal", 1000, -5., 5.);
	hgKplusp           = new TH1D("hgKplusp", "Global K+ momentum", 500, 0., 10.);
	hgKpluspT          = new TH1D("hgKpluspT", "Global K+ transver momentum", 500, 0., 10.);
	hgKplusDCAtoPV     = new TH1D("hgKplusDCAtoPV", "Global K+ DCA to PV", 500, 0., 10.);
	hgKplusDCAtoO      = new TH1D("hgKplusDCAtoO", "Global K+ DCA to Omega", 500, 0., 10.);
	hgKminuspdEdx      = new TH2D("hgKminuspdEdx", "Global K- dE/dx vs p", 5000, 0., 50., 1000, 0., 10.);
	hgKminusdEdxErr    = new TH1D("hgKminusdEdxErr", "Global K- dE/dx error", 100, 0., 1.); // may need to adjust
	hgKminusnsigma     = new TH1D("hgKminusdEdx", "Global K- n_{#sigma}", 500, 0., 5.);
	hgKminuspinvbeta   = new TH2D("hgKminuspinvbeta", "Global K- inverse beta vs p", 5000, 0., 50., 1000, 0., 10.);
	hgKminusm2         = new TH1D("hgKminusdEdx", "Global K- m^2", 500, 0., 5.);
	hgKminusbtofYlocal = new TH1D("hgKminusbtofYlocal", "Global K+ BTOF Ylocal", 1000, -5., 5.);
	hgKminusp          = new TH1D("hgKminusp", "Global K- momentum", 500, 0., 10.);
	hgKminuspT         = new TH1D("hgKminuspT", "Global K- transver momentum", 500, 0., 10.);
	hgKminusDCAtoPV    = new TH1D("hgKminusDCAtoPV", "Global K- DCA to PV", 500, 0., 10.);
	hgKminusDCAtoO     = new TH1D("hgKminusDCAtoO", "Global K- DCA to Omega", 500, 0., 10.);

	// xiatong's analysis
	hCorrKplusO     = new TH1D("hCorrKplusO"    , "K^{+}-#Omega^{-} Correlation"      , 5000, 0.0, 50.0);
    hCorrKplusObar  = new TH1D("hCorrKplusObar" , "K^{+}-#bar{#Omega^{+}} Correlation", 5000, 0.0, 50.0);
    hCorrKminusO    = new TH1D("hCorrKminusO"   , "K^{-}-#Omega^{-} Correlation"      , 5000, 0.0, 50.0);
    hCorrKminusObar = new TH1D("hCorrKminusObar", "K^{-}-#bar{#Omega^{+}} Correlation", 5000, 0.0, 50.0);

	cout << "----------------------------------" << endl;
	cout << "------- histograms claimed -------" << endl;
	cout << "----------------------------------" << endl;

	return;
}

//-----------------------------------------------------------------------------
void StKFParticleAnalysisMaker::WriteHistograms() {

	///////////////////
	hNRefMult ->Write();  
	hNRefMultA->Write();  
	hNRefMultB->Write();  
	hVertexXY ->Write();  
	hVertexZ  ->Write();  
	hVertex2D ->Write();
	hDiffVz   ->Write();
	hcent     ->Write();  
	hcentw    ->Write();  

	hcentRefM ->Write();
	hcentRefW ->Write();

	hCorrKplusO    ->Write();
    hCorrKplusObar ->Write();
    hCorrKminusO   ->Write();
    hCorrKminusObar->Write();

	hgKpluspdEdx      ->Write();
	hgKplusdEdxErr    ->Write();
	hgKplusnsigma     ->Write();
	hgKpluspinvbeta   ->Write();
	hgKplusm2         ->Write();
	hgKplusbtofYlocal ->Write();
	hgKplusp          ->Write();
	hgKpluspT         ->Write();
	hgKplusDCAtoPV    ->Write();
	hgKplusDCAtoO     ->Write();
	hgKminuspdEdx     ->Write();
	hgKminusdEdxErr   ->Write();
	hgKminusnsigma    ->Write();
	hgKminuspinvbeta  ->Write();
	hgKminusm2        ->Write();
	hgKminusbtofYlocal->Write();
	hgKminusp         ->Write();
	hgKminuspT        ->Write();
	hgKminusDCAtoPV   ->Write();
	hgKminusDCAtoO    ->Write();

	return;
}

//-----------------------------------------------------------------------------
void StKFParticleAnalysisMaker::setRunEnergyAndListDir(int run, double energy,char ListDir[256])
{
	mRun    =run;
	mEnergy =energy;
	mListDir=ListDir;
}

//-----------------------------------------------------------------------------
bool StKFParticleAnalysisMaker::readBadList()
{ 
	if(0) return kTRUE;
	TString inf=mListDir + "/badList/";
	inf += Form("badrun%dList%.1f.list",mRun,mEnergy);
	ifstream inrun; 
	inrun.open(inf);
	if ( inrun.fail() ) {
		cout<< "cannot open " << inf.Data() << endl;
		return kFALSE;
	}
	Int_t runid;
	while ( inrun >> runid ) { badList.push_back(runid); }
	inrun.close();
	sort(badList.begin(),badList.end());

	vector<int>::iterator it;
	it = std::unique (badList.begin(), badList.end());
	badList.resize( std::distance(badList.begin(),it) );

	cout <<"badrun list :" <<inf.Data() << " loaded." << endl;
	cout <<"Total       :" <<badList.size()<< " bad runs. "<< endl;
	return kTRUE;
}

//------------------------------------------------------------------------------
bool StKFParticleAnalysisMaker::readRunList()
{ 
	if(0) return kTRUE;
	TString inf=mListDir + "/runList/";
	inf += Form("run%dList%.1f.list",mRun,mEnergy);
	ifstream inrun; 
	inrun.open(inf);
	if ( inrun.fail() ) {
		cout<< "cannot open " << inf.Data() << endl;
		return kFALSE;
	}
	Int_t runid;
	while ( inrun >> runid ) { runList.push_back(runid); }
	inrun.close();
	sort(runList.begin(),runList.end());

	vector<int>::iterator it;
	it = std::unique (runList.begin(), runList.end());
	runList.resize( std::distance(runList.begin(),it) );

	cout <<"Run list :" <<inf.Data() << " loaded. "<< endl;
	cout <<"Total    :" <<runList.size()<< " runs. "<< endl;

	if(runList.size()<1){cout<<"no run number found!!!"<<endl; return kFALSE;}

	return kTRUE;
}

//-----------------------------------------------------------------------------
bool StKFParticleAnalysisMaker::removeBadID(int runnumber)
{
	for (std::vector<int>::iterator it=badList.begin(); it!=badList.end(); ++it) { 
		if(runnumber==*it) {
			return kTRUE;
		}
	}
	return kFALSE;
}  

//-----------------------------------------------------------------------------
Int_t StKFParticleAnalysisMaker::CheckrunNumber(int runnumber)
{    
	int pointer=-999; 
	int id=0; 
	for (std::vector<int>::iterator it=runList.begin(); it!=runList.end(); ++it) { 
		if(runnumber==*it) pointer=id;
		id++;
	}

	if(pointer==-999) cout << "Run number are not found! "<< runnumber << endl;
	return pointer;
} 

//----------------------------------------------------------------------------- 
Int_t StKFParticleAnalysisMaker::findCentrality(int mult)
{
	// NOT in use
	int centrality=0;
	float centFull[9] = {30,46,65,89,118,154,200,262,305};
	if      (mult>centFull[8]) centrality=9;
	else if (mult>centFull[7]) centrality=8;
	else if (mult>centFull[6]) centrality=7;
	else if (mult>centFull[5]) centrality=6;
	else if (mult>centFull[4]) centrality=5;
	else if (mult>centFull[3]) centrality=4;
	else if (mult>centFull[2]) centrality=3;
	else if (mult>centFull[1]) centrality=2;
	else if (mult>centFull[0]) centrality=1;

	return centrality;
}

//----------------------------------------------------------------------------- 
void StKFParticleAnalysisMaker::Clear(Option_t *opt) {
}

//----------------------------------------------------------------------------- 
Int_t StKFParticleAnalysisMaker::Make() 
{
    	PicoDst = StPicoDst::instance(); 		
	StPicoDst* mPicoDst = PicoDst;
	if(!mPicoDst) {
		LOG_WARN << " No PicoDst! Skip! " << endm;
		return kStOK;
	}

	//     pass event  
	/////////////////////////////////////////////////////////
	StPicoEvent* mEvent= (StPicoEvent*) mPicoDst->event(); 
	if(!mEvent)return kStOK;

	const int  runID    = mEvent->runId();
	const int  evtID    = mEvent->eventId();
	const int  refMult  = mEvent->refMult();
	const int grefMult  = mEvent->grefMult();
	const int  ranking  = mEvent->ranking();
	const int tofMatch  = mEvent->nBTOFMatch();

	const double magnet = mEvent->bField();

	if(               (!mEvent->isTrigger(610001))
			&&(!mEvent->isTrigger(610011))
			&&(!mEvent->isTrigger(610021))
			&&(!mEvent->isTrigger(610031))
			&&(!mEvent->isTrigger(610041))
			&&(!mEvent->isTrigger(610051))
	  )return kStOK; 

	const TVector3 Vertex3D=mEvent->primaryVertex();
	const double VertexX = Vertex3D.x(); 
	const double VertexY = Vertex3D.y(); 
	const double VertexZ = Vertex3D.z(); 
	const double VertexR = sqrt(VertexX*VertexX + VertexY*VertexY);
	const double vpdVz   = mEvent->vzVpd();

	//event cut
	//if(refMult <=2 || refMult > 1000) return kStOK;
	if(removeBadID(runID)) return kStOK;            
	if(mRefMultCorr->isBadRun(runID)) return kStOK; // reject bad run of StRefMultCorr
	if(!mRefMultCorr->passnTofMatchRefmultCut(1.*refMult, 1.*tofMatch)) return kStOK; // reject pileup of StRefMultCorr

	hVertex2D ->Fill(VertexZ,vpdVz);
	hDiffVz   ->Fill(VertexZ-vpdVz); 

	const double DVz = VertexZ-vpdVz;

	if(fabs(VertexZ) > 80) return kStOK; 
	if(sqrt(pow(VertexX,2.)+pow(VertexY,2.))>2.0) return kStOK; 
	if(fabs(VertexZ-vpdVz)>3.) return kStOK;       // no vpd cut in low energy?

	//check run number
	int runnumberPointer = -999;
	runnumberPointer = CheckrunNumber(runID);
	if(runnumberPointer==-999) return kStOK;

	int dayPointer = (int)((runID)/1000%1000);
	int mRunL=runList.at(0);
	int mDayL=(int) ((mRunL)/1000%1000);
	dayPointer -= mDayL;
	int timePointer = dayPointer/mStps;

	//StRefMultCorr
	mRefMultCorr->init(runID);
	mRefMultCorr->initEvent(refMult,VertexZ,mEvent->ZDCx());
	double refmultWght = mRefMultCorr->getWeight();
	double refmultCorr = mRefMultCorr->getRefMultCorr() ;
	int centrality     = mRefMultCorr->getCentralityBin9();  // 0 - 8  be careful !!!!!!!! 
	//double mWght    = 1;
	//double mult_corr= refMult;
	//int centrality  = findCentrality(mult_corr)-1;  // 0-8 
	if( centrality<0||centrality>=(nCent-1)) return kStOK;
	int cent = centrality+1;  

	double mWght = refmultWght;
	double mult_corr = refmultCorr;

	///////////////////////////
	hNRefMult ->Fill(grefMult);
	hNRefMultA->Fill(mult_corr);
	hNRefMultB->Fill(mult_corr,mWght);
	hVertexXY ->Fill(VertexX,VertexY);
	hVertexZ  ->Fill(VertexZ); 
	hcent     ->Fill(cent);         // 1-9
	hcentw    ->Fill(cent,mWght);   // 1-9

	///////////////
	hcentRefM    ->Fill(cent,mult_corr);     
	hcentRefW    ->Fill(cent,mult_corr,mWght);  
	hcentRefM    ->Fill(0.,mult_corr);     
	hcentRefW    ->Fill(0.,mult_corr,mWght);  
	///////////////

// ======= KFParticle ======= //
	vector<StLambdaDecayPair> KFParticleLambdaDecayPair;

	SetupKFParticle();
	if (InterfaceCantProcessEvent) return kStOK;
	for (int iKFParticle=0; iKFParticle < KFParticlePerformanceInterface->GetNReconstructedParticles(); iKFParticle++){ 
		const KFParticle particle = KFParticleInterface->GetParticles()[iKFParticle]; 
		int upQ; if (particle.GetPDG() == LambdaPdg) upQ = 1; else if (particle.GetPDG() == -1*LambdaPdg) upQ = -1; else continue;
		int eLambda = -(upQ-1)/2; // 0 if Lambda, 1 if AntiLambda

		SetDaughterTrackPointers(iKFParticle);
		if (ProtonTrackIndex == -99999 || PionTrackIndex == -99999) continue; if(!ProtonTrack) continue; if(!PionTrack) continue;

		double dmass = -999; // just a placeholder
		TLorentzVector p4Pair, p4Proton; // just a placeholder
		StLambdaDecayPair TmpLambdaDecayPair(p4Pair, p4Proton, ProtonTrackIndex, PionTrackIndex, (eLambda==0), dmass);
		KFParticleLambdaDecayPair.push_back(TmpLambdaDecayPair);
	} // End loop over KFParticles

	// correlation function loop  
  	Int_t nTracks = mPicoDst->numberOfTracks( );
	for (Int_t iTrack = 0; iTrack < nTracks; iTrack++) 
	{
    	StPicoTrack *track = mPicoDst->track(iTrack);
    	if (! track)            continue;
    	if (! track->charge())  continue;
    	if (  track->nHitsFit() < 15) continue;

		// TOF Info
		bool hasTOF = false;
		int tofindex = track->bTofPidTraitsIndex();
		if (tofindex >= 0) 
		{
			int tofflag = (mPicoDst->btofPidTraits(tofindex))->btofMatchFlag();
			float tof = (mPicoDst->btofPidTraits(tofindex))->btof();
			float BtofYLocal = (mPicoDst->btofPidTraits(tofindex))->btofYLocal();
			if (track->charge() > 0) hgKplusbtofYlocal->Fill(BtofYLocal);
			if (track->charge() < 0) hgKminusbtofYlocal->Fill(BtofYLocal);
			if((tofflag >= 1) && (tof > 0) && (BtofYLocal > -1.8) && (BtofYLocal < 1.8)) hasTOF = true;
		}

		// fill QA
		StPicoPhysicalHelix helix = track->helix(magnet);
		TVector3 pkaon = helix.momentum(magnet);
		std::cout << pkaon.Mag() << std::endl;
		if (track->charge() > 0) //K+
		{
			hgKpluspdEdx    ->Fill(pkaon.Mag(), track->dEdx());
			hgKplusdEdxErr ->Fill(track->dEdxError());
			hgKplusnsigma  ->Fill(track->nSigmaKaon());
			hgKplusp       ->Fill(track->gMom().Mag());
			hgKpluspT      ->Fill(track->gMom().Perp());
			hgKplusDCAtoPV ->Fill(track->gDCA(Vertex3D).Mag());
			if (hasTOF)
			{
				float beta = (mPicoDst->btofPidTraits(tofindex))->btofBeta();
				hgKpluspinvbeta->Fill(pkaon.Mag(), 1./beta);
				hgKplusm2  ->Fill(pkaon.Mag2()*(1.0 / beta / beta - 1.0));
			}
		}
		if (track->charge() < 0) //K-
		{
			hgKminuspdEdx   ->Fill(pkaon.Mag(), track->dEdx());
			hgKminusdEdxErr ->Fill(track->dEdxError());
			hgKminusnsigma  ->Fill(track->nSigmaKaon());
			hgKminusp       ->Fill(track->gMom().Mag());
			hgKminuspT      ->Fill(track->gMom().Perp());
			hgKminusDCAtoPV ->Fill(track->gDCA(Vertex3D).Mag());
			if (hasTOF)
			{
				float beta = (mPicoDst->btofPidTraits(tofindex))->btofBeta();
				hgKminuspinvbeta->Fill(pkaon.Mag(), 1./beta);
				hgKminusm2  ->Fill(pkaon.Mag2()*(1.0 / beta / beta - 1.0));
			}
		}

		// start with the most basic PID cut
		if (track->nSigmaKaon() > 2 || track->nSigmaProton() < 2 || track->nSigmaPion() < 2) continue;

		// Omega loop
		const int kaonindex = track->id();
		for (int iKFParticle=0; iKFParticle < KFParticlePerformanceInterface->GetNReconstructedParticles(); iKFParticle++)
		{ 
			const KFParticle particle = KFParticleInterface->GetParticles()[iKFParticle]; 
			int upQ; if (particle.GetPDG() == OmegaPdg) upQ = 1; else if (particle.GetPDG() == -1*OmegaPdg) upQ = -1; else continue;	
			if (IsKaonOmegaDaughter(iKFParticle, kaonindex)) continue;

			// pair-wise cut to be considered
			/* */

			// Omega momentum at DCA to PV
			TVector3 pOmega(particle.GetPx(), particle.GetPy(), particle.GetPz());
			TVector3 xOmega(particle.GetX(), particle.GetY(), particle.GetZ());
			StPicoPhysicalHelix helixOmega(pOmega, xOmega, magnet, particle.GetQ());
            double pathlength = helixOmega.pathLength(Vertex3D, false);
            TVector3 pOmega_tb = helixOmega.momentumAt(pathlength, magnet); 

			// k*
			TLorentzVector lv1(pOmega_tb, particle.GetMass());
			TLorentzVector lv2(track->gMom(), KaonPdgMass);
			TLorentzVector P = lv1 + lv2;
			TVector3 pair_beta = P.BoostVector();
			lv1.Boost((-1)*pair_beta); 	
			lv2.Boost((-1)*pair_beta); 		
			if (track->charge() > 0 && particle.GetQ() < 0) hCorrKplusO    ->Fill(0.5*(lv1-lv2).Vect().Mag());
			if (track->charge() > 0 && particle.GetQ() > 0) hCorrKplusObar ->Fill(0.5*(lv1-lv2).Vect().Mag());
			if (track->charge() < 0 && particle.GetQ() < 0) hCorrKminusO   ->Fill(0.5*(lv1-lv2).Vect().Mag());
			if (track->charge() < 0 && particle.GetQ() > 0) hCorrKminusObar->Fill(0.5*(lv1-lv2).Vect().Mag());

		} // End loop over KFParticles
		
	}

// ======= KFParticle end ======= //

// ======= Lambda loop ======= //
	for(int j=0; j<KFParticleLambdaDecayPair.size(); j++) {
		int i = KFParticleLambdaDecayPair[j].get_idxProton();
		int k = KFParticleLambdaDecayPair[j].get_idxPion();
		if(k == i) continue;

		StPicoTrack* mTrackI = (StPicoTrack*)mPicoDst->track(i);

		int    mchgI = mTrackI->charge();
		int    mhitI = mTrackI->nHitsFit();
		double mMomI = mTrackI->gMom().Mag();
		double mp0xI = mTrackI->gMom().X();
		double mp0yI = mTrackI->gMom().Y();
		double mp0zI = mTrackI->gMom().Z();
		double mpt0I = mTrackI->gMom().Perp();
		double mphiI = mTrackI->gMom().Phi();
		double metaI = mTrackI->gMom().PseudoRapidity();
		double mdcaI = mTrackI->gDCA(Vertex3D).Mag();

		if(mphiI<0) mphiI += 2*M_PI;
		if(mphiI>=2*M_PI) mphiI -= 2*M_PI;

		StPicoTrack* mTrackK = (StPicoTrack*)mPicoDst->track(k);

		int    mchgK = mTrackK->charge();
		int    mhitK = mTrackK->nHitsFit();
		double mMomK = mTrackK->gMom().Mag();
		double mp0xK = mTrackK->gMom().X();
		double mp0yK = mTrackK->gMom().Y();
		double mp0zK = mTrackK->gMom().Z();
		double mpt0K = mTrackK->gMom().Perp();
		double mphiK = mTrackK->gMom().Phi();
		double metaK = mTrackK->gMom().PseudoRapidity();
		double mdcaK = mTrackK->gDCA(Vertex3D).Mag();

		if(mphiK<0) mphiK += 2*M_PI;
		if(mphiK>=2*M_PI) mphiK -= 2*M_PI;

		bool isPP = mchgI>0 && mchgK>0;
		bool isNN = mchgI<0 && mchgK<0;
		bool isPN = mchgI>0 && mchgK<0;
		bool isNP = mchgI<0 && mchgK>0;
		bool isSelfLambda = isPN;
		bool isAntiLambda = isNP;

		// remove the SS cases
		if(isPP || isNN) continue;

		//reconstruction of V0, the parent particle
		TVector3 xv0, op1, op2;
		double dca1to2 = closestDistance(mTrackI, mTrackK, magnet, Vertex3D, xv0, op1, op2);
		TVector3 pv0 = op1 + op2;
		TVector3 xv0toPV = xv0 - Vertex3D;
		double rdotp = xv0toPV.Dot(pv0);
		double dcav0toPV = rdotp*rdotp/pv0.Mag2();
		dcav0toPV = sqrt(xv0toPV.Mag2() - dcav0toPV);
		double v0decaylength = xv0toPV.Mag();
		double v0cosrdotp = rdotp/v0decaylength/pv0.Mag();

		TLorentzVector p4ProtonI, p4PionK;
		p4ProtonI.SetPxPyPzE(op1.X(), op1.Y(), op1.Z(), sqrt(op1.Mag2() + pmass*pmass));
		p4PionK.SetPxPyPzE(  op2.X(), op2.Y(), op2.Z(), sqrt(op2.Mag2() + pimass*pimass));

		// ProtonI & PionK
		TLorentzVector p4Pair = p4ProtonI + p4PionK;
		double massPair = p4Pair.M();
		double ptPair   = p4Pair.Pt();
		double etaPair  = p4Pair.Eta();
		double phiPair  = p4Pair.Phi();
	}
// ======= Lambda loop ends ======= //

	/////////////////////////////////////////////////////////
	return kStOK;

}

//------------------------------------------------------------
bool StKFParticleAnalysisMaker::isGoodObs(double Obs) {

	bool mGoodObs = true;
	mGoodObs = mGoodObs && !std::isnan(Obs);
	mGoodObs = mGoodObs && !std::isinf(Obs);
	//mGoodObs = mGoodObs && fabs(Obs)<2.1;

	return mGoodObs;
}

//------------------------------------------------------------
void StKFParticleAnalysisMaker::BookVertexPlots()
{
  KFParticleInterface = new StKFParticleInterface;
  bool storeMCHistograms = false;
  KFParticlePerformanceInterface = new StKFParticlePerformanceInterface(KFParticleInterface->GetTopoReconstructor(), storeMCHistograms);
}

//------------------------------------------------------------
void StKFParticleAnalysisMaker::SetupKFParticle(){
	int maxGBTrackIndex = -1; //find max global track index
	for (unsigned int iTrack = 0; iTrack < PicoDst->numberOfTracks(); iTrack++){
		StPicoTrack *track = PicoDst->track(iTrack);
		if ( !track ) continue;
		if ( track->id() > maxGBTrackIndex ) maxGBTrackIndex = track->id();
	}
	vector<KFMCTrack> mcTracks(0);
	vector<int> triggeredTracks;
	vector<int> mcIndices(maxGBTrackIndex+1);
	for (unsigned int iIndex = 0; iIndex < mcIndices.size(); iIndex++) mcIndices[iIndex] = -1;
	if (maxGBTrackIndex > 0)  KFParticleInterface->ResizeTrackPidVectors(maxGBTrackIndex+1);

	if ( !KFParticleInterface->ProcessEvent(PicoDst, triggeredTracks) ) InterfaceCantProcessEvent = true; else InterfaceCantProcessEvent = false;

	trackMap.resize(maxGBTrackIndex+1, -1); //make a map from trackID to track index in global track array
	for(unsigned int iTrack = 0; iTrack < PicoDst->numberOfTracks(); iTrack++)
	{
		StPicoTrack *track = PicoDst->track(iTrack); if (!track) continue;
		int index = track->id();
		trackMap[index] = iTrack;
	}

	KFParticlePerformanceInterface->SetMCTracks(mcTracks);
	KFParticlePerformanceInterface->SetMCIndexes(mcIndices);    
	KFParticlePerformanceInterface->SetCentralityBin(-1);
	KFParticlePerformanceInterface->SetCentralityWeight(1.);
	KFParticlePerformanceInterface->SetPrintEffFrequency(100000);
	KFParticlePerformanceInterface->PerformanceAnalysis();
} // void SetupKFParticle

//------------------------------------------------------------
void StKFParticleAnalysisMaker::SetDaughterTrackPointers(int iKFParticle){ // Get Daughter Tracks to calculate decay variables manually
	const KFParticle particle = KFParticleInterface->GetParticles()[iKFParticle];
	int upQ; if (particle.GetPDG() == LambdaPdg) upQ = 1; else if (particle.GetPDG() == -1*LambdaPdg) upQ = -1; else upQ = 0;
	ProtonTrackIndex = -99999, PionTrackIndex = -99999;
	for(int iDaughter=0; iDaughter < particle.NDaughters(); iDaughter++){ 
		const int daughterId = particle.DaughterIds()[iDaughter]; 
		const KFParticle daughter = KFParticleInterface->GetParticles()[daughterId]; 
		const int globalTrackId = daughter.DaughterIds()[0];
		int trackIndex = trackMap[globalTrackId];    
		if (daughter.GetPDG() == upQ*ProtonPdg) ProtonTrackIndex = trackIndex;
		else if (daughter.GetPDG() == upQ*PionPdg) PionTrackIndex = trackIndex;
	}  // iDaughter
	ProtonTrack = PicoDst->track(ProtonTrackIndex); PionTrack = PicoDst->track(PionTrackIndex);	
} // void SetDaughterTrackPointers

bool StKFParticleAnalysisMaker::IsKaonOmegaDaughter(int iKFParticle, int kaonTrackId)
{
	const KFParticle particle = KFParticleInterface->GetParticles()[iKFParticle];
	if (fabs(particle.GetPDG()) != OmegaPdg) return false;
	for(int iDaughter=0; iDaughter < particle.NDaughters(); iDaughter++)
	{ 
		const int daughterId = particle.DaughterIds()[iDaughter]; 
		const KFParticle daughter = KFParticleInterface->GetParticles()[daughterId]; 
		if (fabs(daughter.GetPDG()) != KaonPdg) continue;
		const int globalTrackId = daughter.DaughterIds()[0];
		
		if (globalTrackId == kaonTrackId) return true;
	}  // iDaughter
	return false;
}
