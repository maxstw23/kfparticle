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
#include "StBichsel/StdEdxPull.h"
#include "MyToolkit.h"

#define pi                 TMath::Pi()
#define OmegaPdgMass	   1.67245
#define OmegaMassSigma     0.0021
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
	nOmegaEvtProcessed = 0;

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

	// xiatong's global track QA
	hgpdEdx      = new TH2D("hgpdEdx", "Global dE/dx vs p", 1000, 0., 10., 1000, 0., 10.);
	hgdEdxErr    = new TH1D("hgdEdxErr", "Global dE/dx error", 100, 0., 1.);
	hgpinvbeta   = new TH2D("hgpinvbeta", "Global inverse beta vs p", 1000, 0., 10., 1000, 0., 10.);
	hgm2         = new TH1D("hgm2", "Global m^2", 500, -1., 4.);
	hgpm2        = new TH2D("hgpm2", "Global m^2 vs p", 1000, 0., 10., 500, -1., 4.);
	hgm2nSigmaKaon   = new TH2D("hgm2nSigmaKaon", "Global m2 vs nSigmaKaon", 2000, -10, 10, 1200, -1., 5.);
	hgm2nSigmaPion   = new TH2D("hgm2nSigmaPion", "Global m2 vs nSigmaPion", 2000, -10, 10, 1200, -1., 5.);
	hgm2nSigmaProton = new TH2D("hgm2nSigmaProton", "Global m2 vs nSigmaProton", 2000, -10, 10, 1200, -1., 5.);
	hgp          = new TH1D("hgp", "Global momentum", 1000, 0., 10.); 
	hgpT         = new TH1D("hgpT", "Global transverse momentum", 1000, 0., 10.);
	hgDCAtoPV    = new TH1D("hgDCAtoPV", "Global DCA to PV", 500, 0., 10.);
	hgbtofYlocal = new TH1D("hgbtofYlocal", "Global K+ BTOF Ylocal", 1000, -5., 5.);
	// 2D pid
	char temp[200];
	for (int i = 0; i < 15; i++)
	{	
		sprintf(temp, "hgPID2D_pt_%d", i);
		hgPID2D_pt[i] = new TH2D(temp, temp, 2000, -10, 10, 2000, -10, 10);
		sprintf(temp, "hgzTPC_pt_%d", i);
		hgzTPC_pt[i] = new TH1D(temp, temp, 2000, -10, 10);
	}

	// kaon QA
	hgKpdEdx       = new TH2D("hgKpdEdx", "Global K dE/dx vs p", 1000, 0., 10., 1000, 0., 10.);
	hgKpinvbeta    = new TH2D("hgKpinvbeta", "Global K inverse beta vs p", 1000, 0., 10., 1000, 0., 10.);
	hgKm2          = new TH1D("hgKm2", "Global K m^2", 500, -1., 4.);
	hgKpm2         = new TH2D("hgKpm2", "Global K m^2 vs p", 1000, 0., 10., 500, -1., 4.);
	hgKp           = new TH1D("hgKp", "Global K momentum", 1000, 0., 10.);
	hgKpT          = new TH1D("hgKpT", "Global K transver momentum", 1000, 0., 10.);
	hgKDCAtoPV     = new TH1D("hgKDCAtoPV", "Global K DCA to PV", 500, 0., 10.);
	hgKDCAtoO      = new TH1D("hgKDCAtoO", "Global K DCA to Omega", 500, 0., 10.);
	hgKpionpdEdx   = new TH2D("hgKpionpdEdx", "Misidentified kaon dEdx", 1000, 0., 10., 1000, 0., 10.);
	hgKptnSigma    = new TH2D("hgKptnSigma", "Kaon Pt vs nSigmaKaon", 2000, -10, 10, 1000, 0., 10.);

	// Omega QA
	hOmegaM   = new TH1D("hOmegaM", "Omega Invariant Mass", 1400, 1., 2.4);
	hOmegap   = new TH1D("hOmegap", "Omega Momentum", 1000, 0., 10.);
	hOmegapt  = new TH1D("hOmegapt", "Omega Transverse Momentum", 1000, 0., 10.);
	hOmegay   = new TH1D("hOmegay", "Omega Rapidity", 1000, -5., 5.);
	hOmegaypt = new TH2D("hOmegaypt", "Omega Rapidity vs pT", 1000, 0., 10., 1000, -5., 5.);
	hOmegaphi = new TH1D("hOmegaphi", "Omega Phi", 1000, -pi, pi);
	hOmegaDL  = new TH1D("hOmegaDL", "Omega Decay Length", 1000, 0., 10.);
	hOmegabarM   = new TH1D("hOmegabarM", "Omegabar Invariant Mass", 1400, 1., 2.4);
	hOmegabarp   = new TH1D("hOmegabarp", "Omegabar Momentum", 1000, 0., 10.);
	hOmegabarpt  = new TH1D("hOmegabarpt", "Omegabar Transverse Momentum", 1000, 0., 10.);
	hOmegabary   = new TH1D("hOmegabary", "Omegabar Rapidity", 1000, -5., 5.);
	hOmegabarypt = new TH2D("hOmegabarypt", "Omegabar Rapidity vs pT", 1000, 0., 10., 1000, -5., 5.);
	hOmegabarphi = new TH1D("hOmegabarphi", "Omegabar Phi", 1000, -pi, pi);
	hOmegabarDL  = new TH1D("hOmegabarDL", "Omegabar Decay Length", 1000, 0., 10.);

	// Lambda QA
	hLambdaM   = new TH1D("hLambdaM", "Lambda Invariant Mass", 1200, 0.6, 1.8);
	hLambdap   = new TH1D("hLambdap", "Lambda Momentum", 1000, 0., 10.);
	hLambdapt  = new TH1D("hLambdapt", "Lambda Transverse Momentum", 1000, 0., 10.);
	hLambday   = new TH1D("hLambday", "Lambda Rapidity", 1000, -5., 5.);
	hLambdaphi = new TH1D("hLambdaphi", "Lambda Phi", 1000, -pi, pi);
	hLambdaDL  = new TH1D("hLambdaDL", "Lambda Decay Length", 1000, 0., 10.);
	hLambdabarM   = new TH1D("hLambdabarM", "Lambdabar Invariant Mass", 1200, 0.6, 1.8);
	hLambdabarp   = new TH1D("hLambdabarp", "Lambdabar Momentum", 1000, 0., 10.);
	hLambdabarpt  = new TH1D("hLambdabarpt", "Lambdabar Transverse Momentum", 1000, 0., 10.);
	hLambdabary   = new TH1D("hLambdabary", "Lambdabar Rapidity", 1000, -5., 5.);
	hLambdabarphi = new TH1D("hLambdabarphi", "Lambdabar Phi", 1000, -pi, pi);
	hLambdabarDL  = new TH1D("hLambdabarDL", "Lambdabar Decay Length", 1000, 0., 10.);

	// Daughter kaon QA
	hDauKplusp    = new TH1D("hDauKplusp", "Daughter K+ Momentum", 1000, 0., 10.);   
	hDauKpluspt   = new TH1D("hDauKpluspt", "Daughter K+ Transverse Momentum", 1000, 0., 10.);   
	hDauKplusy    = new TH1D("hDauKplusy", "Daughter K+ Rapidity", 1000, -5., 5.);    
	hDauKplusphi  = new TH1D("hDauKplusphi", "Daughter K+ Phi", 1000, -pi, pi);     
	hDauKplusnSigma = new TH1D("hDauKplusnSigma", "Daughter K+ nSigma", 1000, -10, 10);
	hDauKminusp   = new TH1D("hDauKminusp", "Daughter K- Momentum", 1000, 0., 10.);     
	hDauKminuspt  = new TH1D("hDauKminuspt", "Daughter K- Transverse Momentum", 1000, 0., 10.);     
	hDauKminusy   = new TH1D("hDauKminusy", "Daughter K- Rapidity", 1000, -5., 5.);     
	hDauKminusphi = new TH1D("hDauKminusphi", "Daughter K- Phi", 1000, -pi, pi);
	hDauKminusnSigma = new TH1D("hDauKminusnSigma", "Daughter K+ nSigma", 1000, -10, 10);

	// Daughter Lambda QA
	hDauLambdaM   = new TH1D("hDauLambdaM", "Daughter Lambda Invariant Mass", 1200, 0.6, 1.8);
	hDauLambdap   = new TH1D("hDauLambdap", "Daughter Lambda Momentum", 1000, 0., 10.);
	hDauLambdapt  = new TH1D("hDauLambdapt", "Daughter Lambda Transverse Momentum", 1000, 0., 10.);
	hDauLambday   = new TH1D("hDauLambday", "Daughter Lambda Rapidity", 1000, -5., 5.);
	hDauLambdaphi = new TH1D("hDauLambdaphi", "Daughter Lambda Phi", 1000, -pi, pi);
	hDauLambdabarM   = new TH1D("hDauLambdabarM", "Daughter Lambdabar Invariant Mass", 1200, 0.6, 1.8);
	hDauLambdabarp   = new TH1D("hDauLambdabarp", "Daughter Lambdabar Momentum", 1000, 0., 10.);
	hDauLambdabarpt  = new TH1D("hDauLambdabarpt", "Daughter Lambdabar Transverse Momentum", 1000, 0., 10.);
	hDauLambdabary   = new TH1D("hDauLambdabary", "Daughter Lambdabar Rapidity", 1000, -5., 5.);
	hDauLambdabarphi = new TH1D("hDauLambdabarphi", "Daughter Lambdabar Phi", 1000, -pi, pi);

	// Daughter proton and pion
	hDauProtonp   = new TH1D("hDauProtonp", "Daughter Proton Momentum", 1000, 0., 10.);
	hDauProtonpt  = new TH1D("hDauProtonpt", "Daughter Proton Transverse Momentum", 1000, 0., 10.);
	hDauPionp     = new TH1D("hDauPionp", "Daughter Pion Momentum", 1000, 0., 10.);
	hDauPionpt    = new TH1D("hDauPionpt", "Daughter Pion Transverse Momentum", 1000, 0., 10.);
	hDauProtonbarp   = new TH1D("hDauProtonbarp", "Daughter Protonbar Momentum", 1000, 0., 10.);
	hDauProtonbarpt  = new TH1D("hDauProtonbarpt", "Daughter Protonbar Transverse Momentum", 1000, 0., 10.);
	hDauPionbarp     = new TH1D("hDauPionbarp", "Daughter Pionbar Momentum", 1000, 0., 10.);
	hDauPionbarpt    = new TH1D("hDauPionbarpt", "Daughter Pionbar Transverse Momentum", 1000, 0., 10.);

	hOmegaDauPid = new TH1D("hOmegaDauPid", "Omega Daughter PID", 10000, -4999.5, 5000.5);
	hOmegabarDauPid = new TH1D("hOmegabarDauPid", "Omegabar Daughter PID", 10000, -4999.5, 5000.5);

	// Cut Decider
	hDCAOtoK_signal   = new TH1D("hDCAOtoK_signal", "hDCAOtoK_signal", 1000, 0., 10.);
	hDCAOtoK_sideband = new TH1D("hDCAOtoK_sideband", "hDCAOtoK_sideband", 1000, 0., 10.);
	hDCAOtoL_signal   = new TH1D("hDCAOtoL_signal", "hDCAOtoL_signal", 1000, 0., 10.);
	hDCAOtoL_sideband = new TH1D("hDCAOtoL_sideband", "hDCAOtoL_sideband", 1000, 0., 10.);
	
	// xiatong's analysis
	hCorrKplusO     = new TH1D("hCorrKplusO"    , "K^{+}-#Omega^{-} Correlation"      , 5000, 0.0, 50.0);
    hCorrKplusObar  = new TH1D("hCorrKplusObar" , "K^{+}-#bar{#Omega^{+}} Correlation", 5000, 0.0, 50.0);
    hCorrKminusO    = new TH1D("hCorrKminusO"   , "K^{-}-#Omega^{-} Correlation"      , 5000, 0.0, 50.0);
    hCorrKminusObar = new TH1D("hCorrKminusObar", "K^{-}-#bar{#Omega^{+}} Correlation", 5000, 0.0, 50.0);
	hCorrKplusO_mixed     = new TH1D("hCorrKplusO_mixed"    , "K^{+}-#Omega^{-} Correlation"      , 5000, 0.0, 50.0);
    hCorrKplusObar_mixed  = new TH1D("hCorrKplusObar_mixed" , "K^{+}-#bar{#Omega^{+}} Correlation", 5000, 0.0, 50.0);
    hCorrKminusO_mixed    = new TH1D("hCorrKminusO_mixed"   , "K^{-}-#Omega^{-} Correlation"      , 5000, 0.0, 50.0);
    hCorrKminusObar_mixed = new TH1D("hCorrKminusObar_mixed", "K^{-}-#bar{#Omega^{+}} Correlation", 5000, 0.0, 50.0);

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
	hCorrKplusO_mixed    ->Write(); 
    hCorrKplusObar_mixed ->Write(); 
    hCorrKminusO_mixed   ->Write(); 
    hCorrKminusObar_mixed->Write(); 

	hOmegaM  ->Write();
	hOmegap  ->Write();
	hOmegapt ->Write();
	hOmegay  ->Write();
	hOmegaypt->Write();
	hOmegaphi->Write();
	hOmegaDL ->Write();
	hOmegabarM  ->Write();
	hOmegabarp  ->Write();
	hOmegabarpt ->Write();
	hOmegabary  ->Write();
	hOmegabarypt->Write();
	hOmegabarphi->Write();
	hOmegabarDL ->Write();

	hLambdaM  ->Write();
	hLambdap  ->Write();
	hLambdapt ->Write();
	hLambday  ->Write();
	hLambdaphi->Write();
	hLambdaDL ->Write();
	hLambdabarM  ->Write();
	hLambdabarp  ->Write();
	hLambdabarpt ->Write();
	hLambdabary  ->Write();
	hLambdabarphi->Write();
	hLambdabarDL ->Write();

	hDauKplusp   ->Write();
	hDauKpluspt  ->Write();
	hDauKplusy   ->Write();
	hDauKplusphi ->Write();
	hDauKplusnSigma->Write();
	hDauKminusp  ->Write();
	hDauKminuspt ->Write();
	hDauKminusy  ->Write();
	hDauKminusphi->Write();
	hDauKminusnSigma->Write();

	hDauLambdaM     ->Write();
	hDauLambdap     ->Write();
	hDauLambdapt    ->Write();
	hDauLambday     ->Write();
	hDauLambdaphi   ->Write();
	hDauLambdabarM  ->Write();
	hDauLambdabarp  ->Write();
	hDauLambdabarpt ->Write();
	hDauLambdabary  ->Write();
	hDauLambdabarphi->Write();
	hDauPionp      ->Write();
	hDauPionpt     ->Write();
	hDauPionbarp   ->Write();
	hDauPionbarpt  ->Write();
	hDauProtonp    ->Write();
	hDauProtonpt   ->Write();
	hDauProtonbarp ->Write();
	hDauProtonbarpt->Write();

	hOmegaDauPid->Write();
	hOmegabarDauPid->Write();

	hDCAOtoK_signal->Write();
	hDCAOtoK_sideband->Write();
	hDCAOtoL_signal->Write();
	hDCAOtoL_sideband->Write();

	hgpdEdx      ->Write();
	hgdEdxErr    ->Write();
	hgpinvbeta   ->Write();
	hgm2         ->Write();
	hgpm2        ->Write();
	hgm2nSigmaKaon  ->Write();
	hgm2nSigmaPion  ->Write();
	hgm2nSigmaProton->Write();
	hgp          ->Write();
	hgpT         ->Write();
	hgDCAtoPV    ->Write();
	hgbtofYlocal ->Write();
	hgKpdEdx       ->Write();
	hgKpinvbeta    ->Write();
	hgKm2          ->Write();
	hgKpm2         ->Write();
	hgKp           ->Write();
	hgKpT          ->Write();
	hgKDCAtoPV     ->Write();
	hgKDCAtoO      ->Write();
	hgKpionpdEdx   ->Write();
	hgKptnSigma    ->Write();
	for (int i = 0; i < 15; i++)
	{	
		hgPID2D_pt[i]->Write();
		hgzTPC_pt[i]->Write();
	}

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

	// collect lambdas and omegas
	std::vector<KFParticle> OmegaVec;
	std::vector<int> OmegaDauPionIdVec;
	std::vector<int> OmegaDauProtonIdVec;
	for (int iKFParticle=0; iKFParticle < KFParticlePerformanceInterface->GetNReconstructedParticles(); iKFParticle++)
	{ 
		const KFParticle particle = KFParticleInterface->GetParticles()[iKFParticle]; 
		bool IsOmega = false, IsLambda = false;
		if (fabs(particle.GetPDG()) == OmegaPdg) IsOmega = true; else if (fabs(particle.GetPDG()) == LambdaPdg) IsLambda = true; else continue;
		int upQ; if (particle.GetPDG() > 0) upQ = 1; else if (particle.GetPDG() < 0) upQ = -1; else continue;
		OmegaVec.push_back(particle);

		// Omega/lambda QA
		if (IsOmega)
		{
			if (upQ == 1)
			{
				hOmegaM  ->Fill(particle.GetMass());
				//if (fabs(particle.GetMass()-OmegaPdgMass) > OmegaMassSigma*3) continue; // subject to change
				hOmegap  ->Fill(particle.GetMomentum());
				hOmegapt ->Fill(particle.GetPt());
				hOmegay  ->Fill(particle.GetRapidity());
				hOmegaypt->Fill(particle.GetPt(), particle.GetRapidity());
				hOmegaphi->Fill(particle.GetPhi());
				hOmegaDL ->Fill(particle.GetDecayLength());
				
				// helix
				TVector3 pOmega(particle.GetPx(), particle.GetPy(), particle.GetPz());
				TVector3 xOmega(particle.GetX(), particle.GetY(), particle.GetZ());
				StPicoPhysicalHelix helixOmega(pOmega, xOmega, magnet*kilogauss, particle.GetQ());

				// daughter kaon QA
				for (int iDaughter = 0; iDaughter < particle.NDaughters(); iDaughter++)
				{
					const int daughterId = particle.DaughterIds()[iDaughter];
					const KFParticle daughter = KFParticleInterface->GetParticles()[daughterId];
					const int daughterTrackId = daughter.DaughterIds()[0];
					int trackIndex = trackMap[daughterTrackId];   
					StPicoTrack *daughterTrack = mPicoDst->track(trackIndex);
					hOmegaDauPid->Fill(1.0*daughter.GetPDG());

					// for post-reconstruction DCA cut
					TVector3 pDaughter(daughter.GetPx(), daughter.GetPy(), daughter.GetPz());
					TVector3 xDaughter(daughter.GetX(), daughter.GetY(), daughter.GetZ());
					
					if (daughter.GetPDG() == -KaonPdg)
					{
						hDauKminusp  ->Fill(daughter.GetMomentum());
						hDauKminuspt ->Fill(daughter.GetPt());
						hDauKminusy  ->Fill(daughter.GetRapidity());
						hDauKminusphi->Fill(daughter.GetPhi());
						hDauKminusnSigma->Fill(daughterTrack->nSigmaKaon());

						StPicoPhysicalHelix helixDaughter(pDaughter, xDaughter, magnet*kilogauss, daughter.GetQ());
						pair<double, double> tmps = helixDaughter.pathLengths(helixOmega);
						TVector3 ox1 = helixDaughter.at(tmps.first);
						TVector3 ox2 = helixOmega.at(tmps.second);
						double dca = (ox1 - ox2).Mag();
						CutDecider(particle, hDCAOtoK_signal, hDCAOtoK_sideband, dca);
					}
					else 
					{
						hDauLambdaM  ->Fill(daughter.GetMass());
						hDauLambdap  ->Fill(daughter.GetMomentum());
						hDauLambdapt ->Fill(daughter.GetPt());
						hDauLambdaphi->Fill(daughter.GetPhi());
						hDauLambday  ->Fill(daughter.GetRapidity());

						StPicoPhysicalHelix helixDaughter(pDaughter*(1./pDaughter.Perp())*100, xDaughter, magnet*kilogauss, 1);
						pair<double, double> tmps = helixOmega.pathLengths(helixDaughter);
						TVector3 ox1 = helixOmega.at(tmps.first);
						TVector3 ox2 = helixDaughter.at(tmps.second);
						double dca = (ox1 - ox2).Mag();
						CutDecider(particle, hDCAOtoL_signal, hDCAOtoL_sideband, dca);

						for (int iGrandDaughter = 0; iGrandDaughter < daughter.NDaughters(); iGrandDaughter++)
						{
							const int granddaughterId = daughter.DaughterIds()[iGrandDaughter];
							const KFParticle granddaughter = KFParticleInterface->GetParticles()[granddaughterId];
							if (daughter.GetPDG() == ProtonPdg)
							{
								hDauProtonp  ->Fill(granddaughter.GetMomentum());
								hDauProtonpt ->Fill(granddaughter.GetPt());
							}
							else if (daughter.GetPDG() == PionPdg)
							{
								hDauPionbarp  ->Fill(granddaughter.GetMomentum());
								hDauPionbarpt ->Fill(granddaughter.GetPt());
							}
						}
					}
				}
			}
			else
			{
				hOmegabarM  ->Fill(particle.GetMass());
				//if (fabs(particle.GetMass()-OmegaPdgMass) > OmegaMassSigma*3) continue; // subject to change
				hOmegabarp  ->Fill(particle.GetMomentum());
				hOmegabarpt ->Fill(particle.GetPt());
				hOmegabary  ->Fill(particle.GetRapidity());
				hOmegabarypt->Fill(particle.GetPt(), particle.GetRapidity());
				hOmegabarphi->Fill(particle.GetPhi());
				hOmegabarDL ->Fill(particle.GetDecayLength());

				// helix
				TVector3 pOmega(particle.GetPx(), particle.GetPy(), particle.GetPz());
				TVector3 xOmega(particle.GetX(), particle.GetY(), particle.GetZ());
				StPicoPhysicalHelix helixOmega(pOmega, xOmega, magnet*kilogauss, particle.GetQ());

				// daughter kaon QA
				for (int iDaughter = 0; iDaughter < particle.NDaughters(); iDaughter++)
				{
					const int daughterId = particle.DaughterIds()[iDaughter];
					const KFParticle daughter = KFParticleInterface->GetParticles()[daughterId];
					const int daughterTrackId = daughter.DaughterIds()[0];
					int trackIndex = trackMap[daughterTrackId];   
					StPicoTrack *daughterTrack = mPicoDst->track(trackIndex);					
					hOmegabarDauPid->Fill(1.0*daughter.GetPDG());

					// for post-reconstruction DCA cut
					TVector3 pDaughter(daughter.GetPx(), daughter.GetPy(), daughter.GetPz());
					TVector3 xDaughter(daughter.GetX(), daughter.GetY(), daughter.GetZ());

					if (daughter.GetPDG() ==  KaonPdg)
					{
						hDauKplusp  ->Fill(daughter.GetMomentum());
						hDauKpluspt ->Fill(daughter.GetPt());
						hDauKplusy  ->Fill(daughter.GetRapidity());
						hDauKplusphi->Fill(daughter.GetPhi());
						hDauKplusnSigma->Fill(daughterTrack->nSigmaKaon());

						StPicoPhysicalHelix helixDaughter(pDaughter, xDaughter, magnet*kilogauss, daughter.GetQ());
						pair<double, double> tmps = helixDaughter.pathLengths(helixOmega);
						TVector3 ox1 = helixDaughter.at(tmps.first);
						TVector3 ox2 = helixOmega.at(tmps.second);
						double dca = (ox1 - ox2).Mag();
						CutDecider(particle, hDCAOtoK_signal, hDCAOtoK_sideband, dca);
					}
					else 
					{
						hDauLambdabarM  ->Fill(daughter.GetMass());
						hDauLambdabarp  ->Fill(daughter.GetMomentum());
						hDauLambdabarpt ->Fill(daughter.GetPt());
						hDauLambdabarphi->Fill(daughter.GetPhi());
						hDauLambdabary  ->Fill(daughter.GetRapidity());

						StPicoPhysicalHelix helixDaughter(pDaughter*(1./pDaughter.Perp())*100, xDaughter, magnet*kilogauss, 1);
						pair<double, double> tmps = helixOmega.pathLengths(helixDaughter);
						TVector3 ox1 = helixOmega.at(tmps.first);
						TVector3 ox2 = helixDaughter.at(tmps.second);
						double dca = (ox1 - ox2).Mag();
						CutDecider(particle, hDCAOtoL_signal, hDCAOtoL_sideband, dca);

						for (int iGrandDaughter = 0; iGrandDaughter < daughter.NDaughters(); iGrandDaughter++)
						{
							const int granddaughterId = daughter.DaughterIds()[iGrandDaughter];
							const KFParticle granddaughter = KFParticleInterface->GetParticles()[granddaughterId];
							if (daughter.GetPDG() == -ProtonPdg)
							{
								hDauProtonbarp  ->Fill(granddaughter.GetMomentum());
								hDauProtonbarpt ->Fill(granddaughter.GetPt());
							}
							else if (daughter.GetPDG() == -PionPdg)
							{
								hDauPionp  ->Fill(granddaughter.GetMomentum());
								hDauPionpt ->Fill(granddaughter.GetPt());
							}
						}
					}
				}
			}
		}
		else if (IsLambda)
		{
			if (upQ == 1)
			{
				hLambdaM  ->Fill(particle.GetMass());
				hLambdap  ->Fill(particle.GetMomentum());
				hLambdapt ->Fill(particle.GetPt());
				hLambday  ->Fill(particle.GetRapidity());
				hLambdaphi->Fill(particle.GetPhi());
				hLambdaDL ->Fill(particle.GetDecayLength());
			}
			else
			{
				hLambdabarM  ->Fill(particle.GetMass());
				hLambdabarp  ->Fill(particle.GetMomentum());
				hLambdabarpt ->Fill(particle.GetPt());
				hLambdabary  ->Fill(particle.GetRapidity());
				hLambdabarphi->Fill(particle.GetPhi());
				hLambdabarDL ->Fill(particle.GetDecayLength());
			}
		}

	} // End loop over KFParticles

	// correlation function loop  
  	Int_t nTracks = mPicoDst->numberOfTracks( );
	bool hasOmega = false;
	for (Int_t iTrack = 0; iTrack < nTracks; iTrack++) 
	{
    	StPicoTrack *track = mPicoDst->track(iTrack);
    	if (! track)            continue;
    	if (! track->charge())  continue;
    	if (  track->nHitsFit() < 15) continue;
		if (  track->nHitsDedx() < 15) continue;
		if (  track->dEdxError() < 0.04 || track->dEdxError() > 0.12) continue; // same as kfp

		// TOF Info
		bool hasTOF = false;
		int tofindex = track->bTofPidTraitsIndex();
		float m2 = -999.;
		float beta = -999.;
		if (tofindex >= 0) 
		{
			int tofflag = (mPicoDst->btofPidTraits(tofindex))->btofMatchFlag();
			float tof = (mPicoDst->btofPidTraits(tofindex))->btof();
			float BtofYLocal = (mPicoDst->btofPidTraits(tofindex))->btofYLocal();
			hgbtofYlocal->Fill(BtofYLocal);
			if((tofflag >= 1) && (tof > 0) && (BtofYLocal > -1.8) && (BtofYLocal < 1.8)) hasTOF = true;
		}

		// fill QA
		StPicoPhysicalHelix helix = track->helix(magnet);
		TVector3 pkaon = helix.momentum(magnet*kilogauss);
		hgpdEdx   ->Fill(pkaon.Mag(), track->dEdx());
		hgdEdxErr ->Fill(track->dEdxError());
		hgp       ->Fill(track->gMom().Mag());
		hgpT      ->Fill(track->gMom().Perp());
		hgDCAtoPV ->Fill(track->gDCA(Vertex3D).Mag());
		
		int ptbin = static_cast<int>(floor(track->gMom().Perp()/0.2));
		double zTPC = TMath::Log(track->dEdx() / 1e6*StdEdxPull::EvalPred(pkaon.Mag()/KaonPdgMass,1,1)); 
		hgzTPC_pt[ptbin]->Fill(zTPC);
		if (hasTOF)
		{
			beta = (mPicoDst->btofPidTraits(tofindex))->btofBeta();
			m2 = pkaon.Mag2()*(1.0 / beta / beta - 1.0);
			hgpinvbeta->Fill(pkaon.Mag(), 1./beta);
			hgm2  ->Fill(m2);
			hgpm2 ->Fill(pkaon.Mag(), m2);
			hgm2nSigmaKaon  ->Fill(track->nSigmaKaon(), m2);
			hgm2nSigmaPion  ->Fill(track->nSigmaPion(), m2);
			hgm2nSigmaProton->Fill(track->nSigmaProton(), m2);
			if (m2 < 0.34 && m2 > 0.15) hgKptnSigma->Fill(track->nSigmaKaon(), track->gMom().Perp());
			double zTOF = 1/beta - sqrt(KaonPdgMass*KaonPdgMass/pkaon.Mag2()+1);
			hgPID2D_pt[ptbin]->Fill(zTPC, zTOF);
		}

		// kaon PID cut
		if (track->gMom().Mag() < 0.15 || track->gMom().Mag() > 2) continue;
		if (fabs(track->nSigmaKaon()-0.496) > 2) continue;
		if (!hasTOF && track->gMom().Mag() > 0.6) continue;
		if (track->gMom().Mag() > 0.6 && (m2 > 0.34 || m2 < 0.15)) continue;

		// kaon QA
		hgKpdEdx    ->Fill(pkaon.Mag(), track->dEdx());
		hgKp       ->Fill(track->gMom().Mag());
		hgKpT      ->Fill(track->gMom().Perp());
		hgKDCAtoPV ->Fill(track->gDCA(Vertex3D).Mag());
		if (hasTOF)
		{
			float beta = (mPicoDst->btofPidTraits(tofindex))->btofBeta();
			hgKpinvbeta->Fill(pkaon.Mag(), 1./beta);
			hgKm2  ->Fill(m2);
			hgKpm2 ->Fill(pkaon.Mag(), m2);

			// check kaons misidentified as pions
			if (m2 > -0.06 && m2 < 0.1) hgKpionpdEdx->Fill(pkaon.Mag(), track->dEdx());
		}

		// Omega loop
		my_event current_event;
		const int kaonindex = track->id();
		for (int iOmega=0; iOmega < OmegaVec.size(); iOmega++)
		{ 
			const KFParticle particle = OmegaVec[iOmega]; 

			// Omega cut should be added after this line
			if (fabs(particle.GetMass()-OmegaPdgMass) > 0.0021*3) continue; // subject to change

			current_event.push_back(particle);
			if (IsKaonOmegaDaughter(particle, kaonindex)) continue;
			if (!hasOmega) hasOmega = true;

			// pair-wise should be added after this line
			/* */

			// Omega momentum at DCA to PV
			TVector3 pOmega(particle.GetPx(), particle.GetPy(), particle.GetPz());
			TVector3 xOmega(particle.GetX(), particle.GetY(), particle.GetZ());
			StPicoPhysicalHelix helixOmega(pOmega, xOmega, magnet*kilogauss, particle.GetQ());
            double pathlength = helixOmega.pathLength(Vertex3D, false);
			TVector3 pOmega_tb = helixOmega.momentumAt(pathlength, magnet*kilogauss); 

			// k*
			TLorentzVector lv1(pOmega_tb, OmegaPdgMass);
			TLorentzVector lv2(track->gMom(), KaonPdgMass);
			TLorentzVector P = lv1 + lv2;
			TVector3 pair_beta = P.BoostVector();
			lv1.Boost((-1)*pair_beta); 	
			lv2.Boost((-1)*pair_beta); 		
			if (track->charge() > 0 && particle.GetQ() < 0) hCorrKplusO    ->Fill(0.5*(lv1-lv2).Vect().Mag());
			if (track->charge() > 0 && particle.GetQ() > 0) hCorrKplusObar ->Fill(0.5*(lv1-lv2).Vect().Mag());
			if (track->charge() < 0 && particle.GetQ() < 0) hCorrKminusO   ->Fill(0.5*(lv1-lv2).Vect().Mag());
			if (track->charge() < 0 && particle.GetQ() > 0) hCorrKminusObar->Fill(0.5*(lv1-lv2).Vect().Mag());
		} // End loop over regular Omega

		if (!current_event.IsEmptyEvent()) buffer.Add_Reservoir(current_event, cent, VertexZ, nOmegaEvtProcessed+1);
		
		// mixed event
		std::vector<my_event> mixed_events;
		if (!buffer.IsEmpty(cent, VertexZ)) mixed_events = buffer.Sample_All(cent, VertexZ);
		for (int iMixEvent = 0; iMixEvent < mixed_events.size(); iMixEvent++)
		{
			std::vector<KFParticle> particles = mixed_events[iMixEvent].GetParticles();
			for (int iOmega = 0; iOmega < particles.size(); iOmega++)
			{
				const KFParticle particle = particles[iOmega];
				// Omega momentum at DCA to PV
				TVector3 pOmega(particle.GetPx(), particle.GetPy(), particle.GetPz());
				TVector3 xOmega(particle.GetX(), particle.GetY(), particle.GetZ());
				StPicoPhysicalHelix helixOmega(pOmega, xOmega, magnet*kilogauss, particle.GetQ());
				double pathlength = helixOmega.pathLength(Vertex3D, false);
				TVector3 pOmega_tb = helixOmega.momentumAt(pathlength, magnet*kilogauss); 

				// k*
				TLorentzVector lv1(pOmega_tb, OmegaPdgMass);
				TLorentzVector lv2(track->gMom(), KaonPdgMass);
				TLorentzVector P = lv1 + lv2;
				TVector3 pair_beta = P.BoostVector();
				lv1.Boost((-1)*pair_beta); 	
				lv2.Boost((-1)*pair_beta); 		
				if (track->charge() > 0 && particle.GetQ() < 0) hCorrKplusO_mixed    ->Fill(0.5*(lv1-lv2).Vect().Mag());
				if (track->charge() > 0 && particle.GetQ() > 0) hCorrKplusObar_mixed ->Fill(0.5*(lv1-lv2).Vect().Mag());
				if (track->charge() < 0 && particle.GetQ() < 0) hCorrKminusO_mixed   ->Fill(0.5*(lv1-lv2).Vect().Mag());
				if (track->charge() < 0 && particle.GetQ() > 0) hCorrKminusObar_mixed->Fill(0.5*(lv1-lv2).Vect().Mag());
			}
		}
		
	}
	if (hasOmega) nOmegaEvtProcessed++;

// ======= KFParticle end ======= //

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

bool StKFParticleAnalysisMaker::IsKaonOmegaDaughter(KFParticle particle, int kaonTrackId)
{
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

void StKFParticleAnalysisMaker::CutDecider(KFParticle Omega, TH1D* hist_signal, TH1D* hist_sideband, double value)
{	
	double mass_diff = fabs(Omega.GetMass() - OmegaPdgMass);
	if (mass_diff < 2*OmegaMassSigma) hist_signal->Fill(value);
	if (mass_diff > 4*OmegaMassSigma && mass_diff < 6*OmegaMassSigma) hist_sideband->Fill(value);
}
