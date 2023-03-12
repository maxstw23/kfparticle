#include "StKFParticleAnalysisMaker.h"
#include "PhysicalConstants.h"
#include "phys_constants.h"
#include "StBTofUtil/tofPathLength.hh"
#include "StThreeVectorF.hh"
#include "TVector2.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TTree.h"
#include "TFile.h"
#include "TProfile.h"
#include "TLorentzVector.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
#include "StMessMgr.h"
#include <algorithm>
#include <cmath>
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
//#define OmegaMassSigma     0.0021
#define LambdaPdgMass      1.11568
#define ProtonPdgMass      0.938272
#define PionPdgMass        0.139570
#define KaonPdgMass		   0.493677
#define LambdaPdg          3122
#define OmegaPdg           3334
#define KaonPdg			   321
#define ProtonPdg          2212
#define PionPdg           -211
#define cen_cut            9
#define num_pt_bin         10
#define num_phi_bin        10
#define num_mult_bin       5
#define num_vz_bin         16
#define num_EP_bin         6
#define shift_order_asso   20
#define shift_order_EP     20
#define TPCAssoEtaCut      0.1

const float StKFParticleAnalysisMaker::OmegaMassSigma[]    = {0.00219, 0.00228, 0.00218, 0.00262, 0.00247, 0.00243, 0.00280};
const float StKFParticleAnalysisMaker::OmegabarMassSigma[] = {0.00238, 0.00214, 0.00233, 0.00250, 0.00257, 0.00270, 0.00269};
const float StKFParticleAnalysisMaker::OmegaMassPtLowerBin[] = {0.5, 0.9, 1.2, 1.5, 1.8, 2.1, 2.4, 2.7, 3.0, 3.4, 4.0};


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

	// EPD
	/* Eta Weight */
	TH2D wt ("Order1etaWeight","Order1etaWeight",100,1.5,6.5,9,0,9);
    TH2D wt2("Order2etaWeight","Order2etaWeight",100,1.5,6.5,9,0,9);
	float lin[9] = {-1.950, -1.900, -1.850, -1.706, -1.438, -1.340, -1.045, -0.717, -0.700};
	float cub[9] = {0.1608, 0.1600, 0.1600, 0.1595, 0.1457, 0.1369, 0.1092, 0.0772, 0.0700};
	float par1[9] = {474.651,474.651,474.651,474.651,474.651,3.27243e+02,1.72351,1.72351,1.72351};
	float par2[9] = {3.55515,3.55515,3.55515,3.55515,3.55515,3.56938,-7.69075,-7.69075,-7.69075};
	float par3[9] = {1.80162,1.80162,1.80162,1.80162,1.80162,1.67113,4.72770,4.72770,4.72770};
	for (int iy=1; iy<=9; iy++)
	{
		for (int ix=1; ix<101; ix++)
		{
			double eta = wt.GetXaxis()->GetBinCenter(ix);
			wt.SetBinContent(ix,iy, lin[iy-1]*eta+cub[iy-1]*pow(eta,3));
			wt2.SetBinContent(ix,iy, sqrt(1-1/par1[iy-1]/par1[iy-1]/cosh(eta)/cosh(eta))/(1+exp((abs(eta)-par2[iy-1])/par3[iy-1])));
		}
	}

	/* Set up StEpdEpFinder */
	char fname_in[200]; char fname_out[200];
	sprintf(fname_in,  "cent_EPD_CorrectionInput.root");
	sprintf(fname_out, "cent_EPD_CorrectionOutput_%d.root", mJob);
	//mEpdHits = new TClonesArray("StPicoEpdHit");
	//unsigned int found;
	//chain->SetBranchStatus("EpdHit*",1,&found);
	//cout << "EpdHit Branch returned found= " << found << endl;
	//chain->SetBranchAddress("EpdHit",&mEpdHits);
	mEpFinder = new StEpdEpFinder(9,fname_out,fname_in);
  	mEpFinder->SetnMipThreshold(0.3);    	// recommended by EPD group
  	mEpFinder->SetMaxTileWeight(1.0);     	// recommended by EPD group, 1.0 for low multiplicity (BES)
  	mEpFinder->SetEpdHitFormat(2);         	// 2=pico   
	mEpFinder->SetEtaWeights(1,wt);		// eta weight for 1st-order EP
    //mEpFinder->SetEtaWeights(2,wt2);	// eta weight for 2nd-order EP

	/* TPC weight files */
	fTPCShift = new TFile(Form("TPCShiftInput.root"), "READ");
	if (fTPCShift->IsZombie())
	{
		std::cout << "No TPC shift input! No shift correction used. " << std::endl;
		hTPCAssoShiftInput_cos = 0;
		hTPCAssoShiftInput_sin = 0;
		for (int ewFull=0; ewFull<3; ewFull++) 
		{
			hTPCEPShiftInput_cos[ewFull] = 0;
			hTPCEPShiftInput_sin[ewFull] = 0;
		}
	}
	else
	{
		hTPCAssoShiftInput_cos = (TProfile2D*)fTPCShift->Get(Form("TPCAssoShift_cos"));
		hTPCAssoShiftInput_sin = (TProfile2D*)fTPCShift->Get(Form("TPCAssoShift_sin"));
		for (int ewFull=0; ewFull<3; ewFull++) 
		{
			hTPCEPShiftInput_cos[ewFull] = (TProfile2D*)fTPCShift->Get(Form("TPCEPShiftEW%d_cos", ewFull));
			hTPCEPShiftInput_sin[ewFull] = (TProfile2D*)fTPCShift->Get(Form("TPCEPShiftEW%d_sin", ewFull));
		}
	}

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
	buffer.Init();

	badList.clear();
	runList.clear();

	PerformAnalysis = true;
	PerformMixing = false;
	StoringTree = false;
	CutCent = false;
	PtReweighting = true;
	v2Calculation = true;

	if(!readRunList())return kStFatal;
	if(!readBadList())return kStFatal;

	DeclareHistograms();
	if (StoringTree && PerformMixing) return kStFatal; // cannot perform mixing and storing mixing tree at the same time
	if (StoringTree) DeclareTrees();
	if (PerformMixing) ReadTrees();
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

	if (StoringTree)
	{
		char temp[200];
		if (CutCent) sprintf(temp, "omega_mix_cen%d_%d.root", cen_cut, mJob);
		else    	 sprintf(temp, "omega_mix_%d.root", mJob);
		ftree = new TFile(temp, "RECREATE");
		ftree->cd();
		WriteTrees();
		ftree->Close();
	}

	if (PerformMixing) ftree->Close();

	// EPD
	mEpFinder->Finish();

	// v2
	TFile *fTPCShift_out = new TFile(Form("TPCShiftOutput_%d.root", mJob), "RECREATE");
	fTPCShift_out->cd();
	WriteTPCShift();
	fTPCShift_out->Close();

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
	char temp[200];

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
	for (int i = 0; i < 9; i++)
	{
		sprintf(temp, "hRefMultCorr_cent_%d", i+1);
		hRefMultCorr_cent[i] = new TH1D(temp, temp, 1000, -0.5, 999.5);
	}

	// xiatong's global track QA
	hgpdEdx      = new TH2D("hgpdEdx", "Global dE/dx vs p", 1000, 0., 10., 1000, 0., 10.);
	hgdEdxErr    = new TH1D("hgdEdxErr", "Global dE/dx error", 100, 0., 1.);
	hgpinvbeta   = new TH2D("hgpinvbeta", "Global inverse beta vs p", 1000, 0., 10., 1000, 0., 10.);
	hgm2         = new TH1D("hgm2", "Global m^2", 500, -1., 4.);
	hgpm2        = new TH2D("hgpm2", "Global m^2 vs p", 1000, 0., 10., 500, -1., 4.);
	hgm2nSigmaKaon   = new TH2D("hgm2nSigmaKaon", "Global m2 vs nSigmaKaon", 2000, -10, 10, 1200, -1., 5.);
	hgm2nSigmaPion   = new TH2D("hgm2nSigmaPion", "Global m2 vs nSigmaPion", 2000, -10, 10, 1200, -1., 5.);
	hgm2nSigmaProton = new TH2D("hgm2nSigmaProton", "Global m2 vs nSigmaProton", 2000, -10, 10, 1200, -1., 5.);
	hgptnSigmaKaon    = new TH2D("hgptnSigmaKaon",   "Pt vs nSigmaKaon",   1000, 0., 10., 2000, -10., 10.);
	hgptnSigmaPion    = new TH2D("hgptnSigmaPion",   "Pt vs nSigmaPion",   1000, 0., 10., 2000, -10., 10.);
	hgptnSigmaProton  = new TH2D("hgptnSigmaProton", "Pt vs nSigmaProton", 1000, 0., 10., 2000, -10., 10.);
	hgptm2            = new TH2D("hgptm2", "Pt vs m2", 1000, 0., 10., 500, -1., 4.);
	hgp          = new TH1D("hgp", "Global momentum", 1000, 0., 10.); 
	hgpT         = new TH1D("hgpT", "Global transverse momentum", 1000, 0., 10.);
	hgDCAtoPV    = new TH1D("hgDCAtoPV", "Global DCA to PV", 500, 0., 10.);
	hgbtofYlocal = new TH1D("hgbtofYlocal", "Global K+ BTOF Ylocal", 1000, -5., 5.);

	// v2 and EP 
	hTPCAssoPhi         = new TH1D("hTPCAssoPhi",         "hTPCAssoPhi",         1000, -PI, PI);
	hTPCAssoPhi_shifted = new TH1D("hTPCAssoPhi_shifted", "hTPCAssoPhi_shifted", 1000, -PI, PI);
	hTPCAssoShiftOutput_cos = new TProfile2D(Form("TPCAssoShift_cos"), Form("TPCAssoShift_cos"), 
								                 shift_order_asso, 0.5, 1.0 * shift_order_asso + 0.5, 9, -0.5, 8.5);
	hTPCAssoShiftOutput_sin = new TProfile2D(Form("TPCAssoShift_sin"), Form("TPCAssoShift_sin"), 
								                 shift_order_asso, 0.5, 1.0 * shift_order_asso + 0.5, 9, -0.5, 8.5);
	for (int ewFull = 0; ewFull < 3; ewFull++)
	{
		hTPCEPShiftOutput_cos[ewFull] = new TProfile2D(Form("TPCEPShiftEW%d_cos", ewFull), Form("TPCEPShiftEW%d_cos", ewFull), 
								                 shift_order_EP, 0.5, 1.0 * shift_order_EP + 0.5, 9, -0.5, 8.5);
		hTPCEPShiftOutput_sin[ewFull] = new TProfile2D(Form("TPCEPShiftEW%d_sin", ewFull), Form("TPCEPShiftEW%d_sin", ewFull), 
								                 shift_order_EP, 0.5, 1.0 * shift_order_EP + 0.5, 9, -0.5, 8.5);
		hTPCEP_2[ewFull]       = new TH1D(Form("hTPCEPEW%d_2", ewFull), Form("hTPCEPEW%d_2", ewFull), 1000, 0., 2*PI);
		hTPCEP_2_shifted[ewFull] = new TH1D(Form("hTPCEPEW%d_2_shifted", ewFull), Form("hTPCEPEW%d_2_shifted", ewFull), 1000, 0., 2*PI);										 
	}
	hTPCEP_ew_cos = new TProfile("hTPCEP_ew_cos", "hTPCEP_ew_cos", 2, 0.5, 2.5, -1., 1.);
	hOmega_TPC_v2_pt = new TProfile("hOmega_TPC_v2_pt", "hOmega_TPC_v2_pt", 50, 0., 5., -1., 1.);
	hOmegabar_TPC_v2_pt = new TProfile("hOmegabar_TPC_v2_pt", "hOmegabar_TPC_v2_pt", 50, 0., 5., -1., 1.);
	hOmega_EPD_v2_pt = new TProfile("hOmega_EPD_v2_pt", "hOmega_EPD_v2_pt", 50, 0., 5., -1., 1.);
	hOmegabar_EPD_v2_pt = new TProfile("hOmegabar_EPD_v2_pt", "hOmegabar_EPD_v2_pt", 50, 0., 5., -1., 1.);
	hpiplus_EPD_v2     = new TProfile("hpiplus_EPD_v2",  "hpiplus_EPD_v2",  9, 0.5, 9.5, -1., 1.);
	hpiminus_EPD_v2    = new TProfile("hpiminus_EPD_v2", "hpiminus_EPD_v2", 9, 0.5, 9.5, -1., 1.);
	hproton_EPD_v2     = new TProfile("hproton_EPD_v2",  "hproton_EPD_v2",  9, 0.5, 9.5, -1., 1.);
	hantiproton_EPD_v2 = new TProfile("hantiproton_EPD_v2",  "hantiproton_EPD_v2",  9, 0.5, 9.5, -1., 1.);
	for (int i = 0; i < 9; i++)
	{
		hOmega_EPD_v2[i]    = new TProfile(Form("hOmega_EPD_v2_%d",    i+1), Form("hOmega_EPD_v2_%d",    i+1), 1400, 1, 2.4, -1., 1.);
		hOmegabar_EPD_v2[i] = new TProfile(Form("hOmegabar_EPD_v2_%d", i+1), Form("hOmegabar_EPD_v2_%d", i+1), 1400, 1, 2.4, -1., 1.);
		hOmega_TPC_v2[i]    = new TProfile(Form("hOmega_TPC_v2_%d",    i+1), Form("hOmega_TPC_v2_%d",    i+1), 1400, 1, 2.4, -1., 1.);
		hOmegabar_TPC_v2[i] = new TProfile(Form("hOmegabar_TPC_v2_%d", i+1), Form("hOmegabar_TPC_v2_%d", i+1), 1400, 1, 2.4, -1., 1.);
	}

	hEPD_e_EP_1 = new TH1D("hEPD_e_EP_1", "hEPD_e_EP_1", 1000, 0., 2*PI);
	hEPD_w_EP_1 = new TH1D("hEPD_w_EP_1", "hEPD_w_EP_1", 1000, 0., 2*PI);
    hEPD_e_EP_2 = new TH1D("hEPD_e_EP_2", "hEPD_e_EP_2", 1000, 0., 2*PI);
	hEPD_w_EP_2 = new TH1D("hEPD_w_EP_2", "hEPD_w_EP_2", 1000, 0., 2*PI);
	hEPD_full_EP_1 = new TH1D("hEPD_full_EP_1", "hEPD_full_EP_1", 1000, 0., 2*PI);
	hEPD_full_EP_2 = new TH1D("hEPD_full_EP_2", "hEPD_full_EP_2", 1000, 0., 2*PI);
	hEPD_ew_cos = new TProfile("hEPD_ew_cos", "hEPD_ew_cos", 3, 0.5, 3.5, -1., 1.);

	// 2D pid
	/*
	char temp[200];
	for (int i = 0; i < 15; i++)
	{	
		sprintf(temp, "hgPID2D_pt_%d", i);
		hgPID2D_pt[i] = new TH2D(temp, temp, 20000, -10, 10, 4000, -2, 2);
		sprintf(temp, "hgzTPC_pt_%d", i);
		hgzTPC_pt[i] = new TH1D(temp, temp, 20000, -10, 10);
	}
	*/

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
	// hgptm2_largenSigmaKaon = new TH2D("hgptm2_largenSigmaKaon", "hgptm2_largenSigmaKaon", 2000, 0., 10., 1000, -1., 4.); // m2 vs pt for nSigmaKaon > 6
	// hgptm2_smallnSigmaKaon = new TH2D("hgptm2_smallnSigmaKaon", "hgptm2_smallnSigmaKaon", 2000, 0., 10., 1000, -1., 4.); // m2 vs pt for nSigmaKaon < -6
	hgnSigmaDiff = new TH1D("hgnSigmaDiff", "nSigma difference", 1000, -5., 5.);
	hKaonCt = new TProfile("hKaonCt", "kaon count in events w.o/ and w/ #Omega", 2, -0.5, 1.5, 0, 1000);
	hKpluspt_omega     = new TH1D("hKpluspt_omega",     "Global K transver momentum", 1000, 0., 10.);
	hKpluspt_omegabar  = new TH1D("hKpluspt_omegabar",  "Global K transver momentum", 1000, 0., 10.);
	hKminuspt_omega    = new TH1D("hKminuspt_omega",    "Global K transver momentum", 1000, 0., 10.);
	hKminuspt_omegabar = new TH1D("hKminuspt_omegabar", "Global K transver momentum", 1000, 0., 10.);
	hKpluseta_omega     = new TH1D("hKpluseta_omega",     "Global K transver momentum", 1000, -5., 5.);
	hKpluseta_omegabar  = new TH1D("hKpluseta_omegabar",  "Global K transver momentum", 1000, -5., 5.);
	hKminuseta_omega    = new TH1D("hKminuseta_omega",    "Global K transver momentum", 1000, -5., 5.);
	hKminuseta_omegabar = new TH1D("hKminuseta_omegabar", "Global K transver momentum", 1000, -5., 5.);
	hKplusphi_omega     = new TH1D("hKplusphi_omega",     "Global K transver momentum", 1000, -pi, pi);
	hKplusphi_omegabar  = new TH1D("hKplusphi_omegabar",  "Global K transver momentum", 1000, -pi, pi);
	hKminusphi_omega    = new TH1D("hKminusphi_omega",    "Global K transver momentum", 1000, -pi, pi);
	hKminusphi_omegabar = new TH1D("hKminusphi_omegabar", "Global K transver momentum", 1000, -pi, pi);

	// proton QA
	hProtony     = new TH1D("hProtony", "Proton Rapidity", 1000, -5., 5.);
	hAntiProtony = new TH1D("hAntiProtony", "Anti-proton Rapidity", 1000, -5., 5.);

	// Omega QA
	hOmegaM   = new TH1D("hOmegaM", "Omega Invariant Mass", 1400, 1., 2.4);
	for (int i = 0; i < 9; i++) {sprintf(temp, "hOmegaM_cen_%d", i+1); hOmegaM_cen[i] = new TH1D(temp, temp, 1400, 1., 2.4);}
	hOmegap   = new TH1D("hOmegap", "Omega Momentum", 1000, 0., 10.);
	hOmegapt  = new TH1D("hOmegapt", "Omega Transverse Momentum", 1000, 0., 10.);
	hOmegaeta = new TH1D("hOmegaeta", "Omega Pseudurapidity", 1000, -5., 5.);
	hOmegay   = new TH1D("hOmegay", "Omega Rapidity", 1000, -5., 5.);
	hOmegaypt = new TH2D("hOmegaypt", "Omega Rapidity vs pT", 1000, 0., 10., 1000, -5., 5.);
	hOmegaphi = new TH1D("hOmegaphi", "Omega Phi", 1000, -pi, pi);
	hOmegaDL  = new TH1D("hOmegaDL", "Omega Decay Length", 1000, 0., 10.);
	hOmegaDLminusLambdaDL = new TH1D("hOmegaDLminusLambdaDL", "hOmegaDLminusLambdaDL", 1000, -5., 5.);

	hOmegabarM   = new TH1D("hOmegabarM", "Omegabar Invariant Mass", 1400, 1., 2.4);
	for (int i = 0; i < 9; i++) {sprintf(temp, "hOmegabarM_cen_%d", i+1); hOmegabarM_cen[i] = new TH1D(temp, temp, 1400, 1., 2.4);}
	hOmegabarp   = new TH1D("hOmegabarp", "Omegabar Momentum", 1000, 0., 10.);
	hOmegabarpt  = new TH1D("hOmegabarpt", "Omegabar Transverse Momentum", 1000, 0., 10.);
	hOmegabareta = new TH1D("hOmegabareta", "Omega Pseudurapidity", 1000, -5., 5.);
	hOmegabary   = new TH1D("hOmegabary", "Omegabar Rapidity", 1000, -5., 5.);
	hOmegabarypt = new TH2D("hOmegabarypt", "Omegabar Rapidity vs pT", 1000, 0., 10., 1000, -5., 5.);
	hOmegabarphi = new TH1D("hOmegabarphi", "Omegabar Phi", 1000, -pi, pi);
	hOmegabarDL  = new TH1D("hOmegabarDL", "Omegabar Decay Length", 1000, 0., 10.);
	hOmegabarDLminusLambdaDL = new TH1D("hOmegabarDLminusLambdaDL", "hOmegabarDLminusLambdaDL", 1000, -5., 5.);

	hNumOmega = new TH1D("hNumOmega", "Number of Omega in an event", 10, -0.5, 9.5);
	hOmegaUsed = new TH1D("hOmegaUsed", "Actual Omega Used #Omega/#bar{#Omega}", 4, -0.5, 3.5);
	hOmegaUsed_sideband = new TH1D("hOmegaUsed_sideband", "Actual Omega Used #Omega/#bar{#Omega}", 2, -0.5, 1.5);

	// pT spectrum
	hOmegaM_error_1 = new TH1D("hOmegaM_error_1", "Difference between KFP mass and daughter tracks calculated mass", 1000, -0.1, 0.1);
	hOmegabarM_error_1 = new TH1D("hOmegabarM_error_1", "Difference between KFP mass and daughter tracks calculated mass", 1000, -0.1, 0.1);
	hOmegaM_error_2 = new TH1D("hOmegaM_error_2", "Difference between KFP GetMass and manual four-vector mass", 1000, -0.1, 0.1);
	hOmegabarM_error_2 = new TH1D("hOmegabarM_error_2", "Difference between KFP GetMass and daughter manual four-vector mass", 1000, -0.1, 0.1);
	hDauKaonM_error = new TH1D("hKaonM_error", "Difference between kaon inv mass and PDG mass", 1000, -0.1, 0.1);
	hDauLambdaM_error = new TH1D("hLambdaM_error", "Difference between lambda inv mass and PDG mass", 1000, -0.1, 0.1);

	for (int i = 0; i < 9; i++)
	{
		for (int j = 0; j < num_pt_bin; j++) 
		{
			sprintf(temp, "hOmegaM_pt_%d_%d", i+1, j+1);
			hOmegaM_pt[i][j]    = new TH1D(temp, temp, 1400, 1., 2.4);
			sprintf(temp, "hOmegabarM_pt_%d_%d", i+1, j+1);
			hOmegabarM_pt[i][j] = new TH1D(temp, temp, 1400, 1., 2.4);	
			// sprintf(temp, "hOmegaM_rotbkg_pi_pt_%d", i+1);
			// hOmegaM_rotbkg_pi_pt[i]    = new TH1D(temp, temp, 1400, 1., 2.4);
			// sprintf(temp, "hOmegabarM_rotbkg_pi_pt_%d", i+1);
			// hOmegabarM_rotbkg_pi_pt[i] = new TH1D(temp, temp, 1400, 1., 2.4);
		}
	}

	// Omega v2
	for (int i = 0; i < 9; i++) // for (int i = 0; i < num_phi_bin; i++)
	{
		for (int j = 0; j < num_phi_bin; j++) // for (int j = 0; j < num_pt_bin; j++) 
		{
			sprintf(temp, "hOmegaM_phi_%d_%d", i+1, j+1);
			hOmegaM_phi[i][j]    = new TH1D(temp, temp, 1400, 1., 2.4);
			sprintf(temp, "hOmegabarM_phi_%d_%d", i+1, j+1);
			hOmegabarM_phi[i][j] = new TH1D(temp, temp, 1400, 1., 2.4);	
			// sprintf(temp, "hOmegaM_rotbkg_pi_pt_%d", i+1);
			// hOmegaM_rotbkg_pi_pt[i]    = new TH1D(temp, temp, 1400, 1., 2.4);
			// sprintf(temp, "hOmegabarM_rotbkg_pi_pt_%d", i+1);
			// hOmegabarM_rotbkg_pi_pt[i] = new TH1D(temp, temp, 1400, 1., 2.4);
		}
	}

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
	hDauLambdaDL  = new TH1D("hDauLambdaDL", "Daughter Lambda Decay Length", 1000, 0., 10.);
	hDauLambdabarM   = new TH1D("hDauLambdabarM", "Daughter Lambdabar Invariant Mass", 1200, 0.6, 1.8);
	hDauLambdabarp   = new TH1D("hDauLambdabarp", "Daughter Lambdabar Momentum", 1000, 0., 10.);
	hDauLambdabarpt  = new TH1D("hDauLambdabarpt", "Daughter Lambdabar Transverse Momentum", 1000, 0., 10.);
	hDauLambdabary   = new TH1D("hDauLambdabary", "Daughter Lambdabar Rapidity", 1000, -5., 5.);
	hDauLambdabarphi = new TH1D("hDauLambdabarphi", "Daughter Lambdabar Phi", 1000, -pi, pi);
	hDauLambdabarDL  = new TH1D("hDauLambdabarDL", "Daughter Lambdabar Decay Length", 1000, 0., 10.);

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

	// event mixing QA
	hNumMixedEvent = new TH1D("hNumMixedEvent", "Number of mixed events", 20, -0.5, 19.5);
	hTotalMixedEvent = new TH1D("hTotalMixedEvent", "hTotalMixedEvent", 8000, -0.5, 7999.5);
	
	// xiatong's analysis
	hCorrKplusO        = new TH1D("hCorrKplusO"       , "K^{+}-#Omega^{-} Correlation"      , 5000, 0.0, 50.0);
    hCorrKplusObar     = new TH1D("hCorrKplusObar"    , "K^{+}-#bar{#Omega^{+}} Correlation", 5000, 0.0, 50.0);
    hCorrKminusO       = new TH1D("hCorrKminusO"      , "K^{-}-#Omega^{-} Correlation"      , 5000, 0.0, 50.0);
    hCorrKminusObar    = new TH1D("hCorrKminusObar"   , "K^{-}-#bar{#Omega^{+}} Correlation", 5000, 0.0, 50.0);
	hPtCorrKplusO      = new TH1D("hPtCorrKplusO"     , "K^{+}-#Omega^{-} Correlation"      , 5000, 0.0, 50.0);
    hPtCorrKplusObar   = new TH1D("hPtCorrKplusObar"  , "K^{+}-#bar{#Omega^{+}} Correlation", 5000, 0.0, 50.0);
    hPtCorrKminusO     = new TH1D("hPtCorrKminusO"    , "K^{-}-#Omega^{-} Correlation"      , 5000, 0.0, 50.0);
    hPtCorrKminusObar  = new TH1D("hPtCorrKminusObar" , "K^{-}-#bar{#Omega^{+}} Correlation", 5000, 0.0, 50.0);
	hyCorrKplusO       = new TH1D("hyCorrKplusO"      , "K^{+}-#Omega^{-} Correlation"      , 2000, -10.0, 10.0);
    hyCorrKplusObar    = new TH1D("hyCorrKplusObar"   , "K^{+}-#bar{#Omega^{+}} Correlation", 2000, -10.0, 10.0);
    hyCorrKminusO      = new TH1D("hyCorrKminusO"     , "K^{-}-#Omega^{-} Correlation"      , 2000, -10.0, 10.0);
    hyCorrKminusObar   = new TH1D("hyCorrKminusObar"  , "K^{-}-#bar{#Omega^{+}} Correlation", 2000, -10.0, 10.0);
	hphiCorrKplusO     = new TH1D("hphiCorrKplusO"    , "K^{+}-#Omega^{-} Correlation"      , 1000, -pi, pi);
    hphiCorrKplusObar  = new TH1D("hphiCorrKplusObar" , "K^{+}-#bar{#Omega^{+}} Correlation", 1000, -pi, pi);
    hphiCorrKminusO    = new TH1D("hphiCorrKminusO"   , "K^{-}-#Omega^{-} Correlation"      , 1000, -pi, pi);
    hphiCorrKminusObar = new TH1D("hphiCorrKminusObar", "K^{-}-#bar{#Omega^{+}} Correlation", 1000, -pi, pi);
	hCorrKplusO_mixed     = new TH1D("hCorrKplusO_mixed"    , "K^{+}-#Omega^{-} Correlation"      , 5000, 0.0, 50.0);
    hCorrKplusObar_mixed  = new TH1D("hCorrKplusObar_mixed" , "K^{+}-#bar{#Omega^{+}} Correlation", 5000, 0.0, 50.0);
    hCorrKminusO_mixed    = new TH1D("hCorrKminusO_mixed"   , "K^{-}-#Omega^{-} Correlation"      , 5000, 0.0, 50.0);
    hCorrKminusObar_mixed = new TH1D("hCorrKminusObar_mixed", "K^{-}-#bar{#Omega^{+}} Correlation", 5000, 0.0, 50.0);
	hCorrKplusO_sideband        = new TH1D("hCorrKplusO_sideband"       , "K^{+}-#Omega^{-} Correlation"      , 5000, 0.0, 50.0);
    hCorrKplusObar_sideband     = new TH1D("hCorrKplusObar_sideband"    , "K^{+}-#bar{#Omega^{+}} Correlation", 5000, 0.0, 50.0);
    hCorrKminusO_sideband       = new TH1D("hCorrKminusO_sideband"      , "K^{-}-#Omega^{-} Correlation"      , 5000, 0.0, 50.0);
    hCorrKminusObar_sideband    = new TH1D("hCorrKminusObar_sideband"   , "K^{-}-#bar{#Omega^{+}} Correlation", 5000, 0.0, 50.0);
	hPtCorrKplusO_sideband      = new TH1D("hPtCorrKplusO_sideband"     , "K^{+}-#Omega^{-} Correlation"      , 5000, 0.0, 50.0);
    hPtCorrKplusObar_sideband   = new TH1D("hPtCorrKplusObar_sideband"  , "K^{+}-#bar{#Omega^{+}} Correlation", 5000, 0.0, 50.0);
    hPtCorrKminusO_sideband     = new TH1D("hPtCorrKminusO_sideband"    , "K^{-}-#Omega^{-} Correlation"      , 5000, 0.0, 50.0);
    hPtCorrKminusObar_sideband  = new TH1D("hPtCorrKminusObar_sideband" , "K^{-}-#bar{#Omega^{+}} Correlation", 5000, 0.0, 50.0);
	hyCorrKplusO_sideband       = new TH1D("hyCorrKplusO_sideband"      , "K^{+}-#Omega^{-} Correlation"      , 2000, -10.0, 10.0);
    hyCorrKplusObar_sideband    = new TH1D("hyCorrKplusObar_sideband"   , "K^{+}-#bar{#Omega^{+}} Correlation", 2000, -10.0, 10.0);
    hyCorrKminusO_sideband      = new TH1D("hyCorrKminusO_sideband"     , "K^{-}-#Omega^{-} Correlation"      , 2000, -10.0, 10.0);
    hyCorrKminusObar_sideband   = new TH1D("hyCorrKminusObar_sideband"  , "K^{-}-#bar{#Omega^{+}} Correlation", 2000, -10.0, 10.0);
	hphiCorrKplusO_sideband     = new TH1D("hphiCorrKplusO_sideband"    , "K^{+}-#Omega^{-} Correlation"      , 1000, -pi, pi);
    hphiCorrKplusObar_sideband  = new TH1D("hphiCorrKplusObar_sideband" , "K^{+}-#bar{#Omega^{+}} Correlation", 1000, -pi, pi);
    hphiCorrKminusO_sideband    = new TH1D("hphiCorrKminusO_sideband"   , "K^{-}-#Omega^{-} Correlation"      , 1000, -pi, pi);
    hphiCorrKminusObar_sideband = new TH1D("hphiCorrKminusObar_sideband", "K^{-}-#bar{#Omega^{+}} Correlation", 1000, -pi, pi);

	hNegPtDiff_dphi_KmOb = new TH1D("hNegPtDiff_dphi_KmOb", "#Delta#phi for pairs with negative #Delta p_{T}", 1000, -pi, pi);
	hPosPtDiff_dphi_KmOb = new TH1D("hPosPtDiff_dphi_KmOb", "#Delta#phi for pairs with positive #Delta p_{T}", 1000, -pi, pi);
	hNegPtDiff_dphi_KmO = new TH1D("hNegPtDiff_dphi_KmO", "#Delta#phi for pairs with negative #Delta p_{T}", 1000, -pi, pi);
	hPosPtDiff_dphi_KmO = new TH1D("hPosPtDiff_dphi_KmO", "#Delta#phi for pairs with positive #Delta p_{T}", 1000, -pi, pi);

	hCorrKplusO_y_pT = new TH2D("hCorrKplusO_y_pT", "hCorrKplusO_y_pT", 500, 0.0, 5.0, 200, 0.0, 2.0);
	hCorrKplusObar_y_pT = new TH2D("hCorrKplusObar_y_pT", "hCorrKplusObar_y_pT", 500, 0.0, 5.0, 200, 0.0, 2.0);
	hCorrKminusO_y_pT = new TH2D("hCorrKminusO_y_pT", "hCorrKminusO_y_pT", 500, 0.0, 5.0, 200, 0.0, 2.0); 
	hCorrKminusObar_y_pT = new TH2D("hCorrKminusObar_y_pT", "hCorrKminusObar_y_pT", 500, 0.0, 5.0, 200, 0.0, 2.0);
	hCorrKplusO_y_pT_mixed = new TH2D("hCorrKplusO_y_pT_mixed", "hCorrKplusO_y_pT_mixed", 500, 0.0, 5.0, 200, 0.0, 2.0);
	hCorrKplusObar_y_pT_mixed = new TH2D("hCorrKplusObar_y_pT_mixed", "hCorrKplusObar_y_pT_mixed", 500, 0.0, 5.0, 200, 0.0, 2.0);
	hCorrKminusO_y_pT_mixed = new TH2D("hCorrKminusO_y_pT_mixed", "hCorrKminusO_y_pT_mixed", 500, 0.0, 5.0, 200, 0.0, 2.0); 
	hCorrKminusObar_y_pT_mixed = new TH2D("hCorrKminusObar_y_pT_mixed", "hCorrKminusObar_y_pT_mixed", 500, 0.0, 5.0, 200, 0.0, 2.0);

	hCorrKplusO_y_phi = new TH2D("hCorrKplusO_y_phi", "hCorrKplusO_y_phi", 500, 0.0, pi, 200, 0.0, 2.0);
	hCorrKplusObar_y_phi = new TH2D("hCorrKplusObar_y_phi", "hCorrKplusObar_y_phi", 500, 0.0, pi, 200, 0.0, 2.0);
	hCorrKminusO_y_phi = new TH2D("hCorrKminusO_y_phi", "hCorrKminusO_y_phi", 500, 0.0, pi, 200, 0.0, 2.0); 
	hCorrKminusObar_y_phi = new TH2D("hCorrKminusObar_y_phi", "hCorrKminusObar_y_phi", 500, 0.0, pi, 200, 0.0, 2.0);

	hCorrKplusO_phi_pT = new TH2D("hCorrKplusO_phi_pT", "hCorrKplusO_phi_pT", 500, 0.0, 5.0, 500, 0.0, pi);
	hCorrKplusObar_phi_pT = new TH2D("hCorrKplusObar_phi_pT", "hCorrKplusObar_phi_pT", 500, 0.0, 5.0, 500, 0.0, pi);
	hCorrKminusO_phi_pT = new TH2D("hCorrKminusO_phi_pT", "hCorrKminusO_phi_pT", 500, 0.0, 5.0, 500, 0.0, pi); 
	hCorrKminusObar_phi_pT = new TH2D("hCorrKminusObar_phi_pT", "hCorrKminusObar_phi_pT", 500, 0.0, 5.0, 500, 0.0, pi);

	// a new test observable
    // kaon ratios at different p/pbar bins, in three different scenarios: with one omega, without o/ob, with one omegabar
    hKratio_omega    = new TProfile("hKratio_omega"   , "hKratio_omega"   , 20, 0., 1., 0., 10.);
    hKratio_wo       = new TProfile("hKratio_wo"      , "hKratio_wo"      , 20, 0., 1., 0., 10.);
    hKratio_omegabar = new TProfile("hKratio_omegabar", "hKratio_omegabar", 20, 0., 1., 0., 10.);

	cout << "----------------------------------" << endl;
	cout << "------- histograms claimed -------" << endl;
	cout << "----------------------------------" << endl;

	return;
}

//-----------------------------------------------------------------------------
void StKFParticleAnalysisMaker::DeclareTrees()
{
	char temp[200];
	int split_level = 1;
    int buffer_size = 5000000;

	for (int i = 0; i < num_mult_bin; i++)
	{
		for (int j = 0; j < num_vz_bin; j++)
		{
			for (int k = 0; k < num_EP_bin; k++)
			{
				sprintf(temp, "omega_tree_%d_%d_%d", i, j, k);
				omega_mix[i][j][k] = new TTree(temp, temp);
				gROOT->ProcessLine("#include <vector>");
				omega_mix[i][j][k]->Branch("px", &mix_px, buffer_size, split_level);
				omega_mix[i][j][k]->Branch("py", &mix_py, buffer_size, split_level);
				omega_mix[i][j][k]->Branch("pz", &mix_pz, buffer_size, split_level);
				omega_mix[i][j][k]->Branch("charge", &mix_charge, buffer_size, split_level);
				omega_mix[i][j][k]->Branch("evt_id", &mix_evt_id, "evt_id/I");
				omega_mix[i][j][k]->Branch("run_id", &mix_run_id, "run_id/I");
			}
		}
	}
	cout << "----------------------------------" << endl;
	cout << "--------- Treess claimed ---------" << endl;
	cout << "----------------------------------" << endl;
}

void StKFParticleAnalysisMaker::ReadTrees()
{
	char temp[200];
	sprintf(temp, "./mix/omega_mix_cen_%d.root", cen_cut);
	ftree = new TFile(temp, "READ");

	for (int i = 0; i < num_mult_bin; i++)
	{
		for (int j = 0; j < num_vz_bin; j++)
		{
			for (int k = 0; k < num_EP_bin; k++)
			{
				sprintf(temp, "omega_tree_%d_%d_%d", i, j, k);
				gDirectory->GetObject(temp, omega_mix[i][j][k]); 
			}
		}
	}

	cout << "----------------------------------" << endl;
	cout << "----------- trees read -----------" << endl;
	cout << "----------------------------------" << endl;

	return;
}

void StKFParticleAnalysisMaker::WriteTPCShift()
{
	hTPCAssoShiftOutput_cos->Write();
	hTPCAssoShiftOutput_sin->Write();
	for (int ewFull = 0; ewFull < 3; ewFull++)
	{
		hTPCEPShiftOutput_cos[ewFull]->Write();
		hTPCEPShiftOutput_sin[ewFull]->Write();
	}
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
	for (int i = 0; i < 9; i++) hRefMultCorr_cent[i]->Write();

	hCorrKplusO    ->Write();
    hCorrKplusObar ->Write();
    hCorrKminusO   ->Write();
    hCorrKminusObar->Write();
	hPtCorrKplusO    ->Write();
    hPtCorrKplusObar ->Write();
    hPtCorrKminusO   ->Write();
    hPtCorrKminusObar->Write();
	hyCorrKplusO    ->Write();
    hyCorrKplusObar ->Write();
    hyCorrKminusO   ->Write();
    hyCorrKminusObar->Write();
	hphiCorrKplusO    ->Write();
    hphiCorrKplusObar ->Write();
    hphiCorrKminusO   ->Write();
    hphiCorrKminusObar->Write();
	hCorrKplusO_sideband       ->Write();
    hCorrKplusObar_sideband    ->Write();
    hCorrKminusO_sideband      ->Write();
    hCorrKminusObar_sideband   ->Write();
	hPtCorrKplusO_sideband     ->Write();
    hPtCorrKplusObar_sideband  ->Write();
    hPtCorrKminusO_sideband    ->Write();
    hPtCorrKminusObar_sideband ->Write();
	hyCorrKplusO_sideband      ->Write();
    hyCorrKplusObar_sideband   ->Write();
    hyCorrKminusO_sideband     ->Write();
    hyCorrKminusObar_sideband  ->Write();
	hphiCorrKplusO_sideband    ->Write();
    hphiCorrKplusObar_sideband ->Write();
    hphiCorrKminusO_sideband   ->Write();
    hphiCorrKminusObar_sideband->Write();
	hCorrKplusO_mixed    ->Write(); 
    hCorrKplusObar_mixed ->Write(); 
    hCorrKminusO_mixed   ->Write(); 
    hCorrKminusObar_mixed->Write(); 
	hNegPtDiff_dphi_KmOb->Write(); 
	hPosPtDiff_dphi_KmOb->Write(); 
	hNegPtDiff_dphi_KmO->Write(); 
	hPosPtDiff_dphi_KmO->Write(); 
	hCorrKplusO_y_pT  ->Write(); hCorrKplusObar_y_pT  ->Write(); hCorrKminusO_y_pT  ->Write(); hCorrKminusObar_y_pT  ->Write();
	hCorrKplusO_y_pT_mixed->Write(); hCorrKplusObar_y_pT_mixed->Write(); hCorrKminusO_y_pT_mixed->Write(); hCorrKminusObar_y_pT_mixed->Write();
	hCorrKplusO_y_phi ->Write(); hCorrKplusObar_y_phi ->Write(); hCorrKminusO_y_phi ->Write(); hCorrKminusObar_y_phi ->Write();
	hCorrKplusO_phi_pT->Write(); hCorrKplusObar_phi_pT->Write(); hCorrKminusO_phi_pT->Write(); hCorrKminusObar_phi_pT->Write();

	hKratio_omega   ->Write();
    hKratio_wo      ->Write();
    hKratio_omegabar->Write();

	hProtony->Write();
	hAntiProtony->Write();

	// v2
	for (int ewFull = 0; ewFull < 3; ewFull++)
	{
		hTPCEP_2[ewFull]      ->Write();
		hTPCEP_2_shifted[ewFull]->Write();
	}
	hTPCEP_ew_cos->Write();
	hOmega_TPC_v2_pt->Write();
	hOmegabar_TPC_v2_pt->Write();
	hOmega_EPD_v2_pt->Write();
	hOmegabar_EPD_v2_pt->Write();
	hEPD_e_EP_1->Write();
	hEPD_w_EP_1->Write();
	hEPD_e_EP_2->Write();
	hEPD_w_EP_2->Write();
	hEPD_full_EP_1->Write();
	hEPD_full_EP_2->Write();
	hEPD_ew_cos->Write();
	for (int i = 0; i < 9; i++)
	{
		hOmega_EPD_v2[i]->Write();
		hOmegabar_EPD_v2[i]->Write();
	}
	hpiplus_EPD_v2->Write();
	hpiminus_EPD_v2->Write();
	hproton_EPD_v2->Write();
	hantiproton_EPD_v2->Write();

	hOmegaM  ->Write();
	for (int i = 0; i < 9; i++) hOmegaM_cen[i]->Write();
	hOmegap  ->Write();
	hOmegapt ->Write();
	hOmegaeta->Write();
	hOmegay  ->Write();
	hOmegaypt->Write();
	hOmegaphi->Write();
	hOmegaDL ->Write();
	hOmegaDLminusLambdaDL->Write();
	hOmegabarM  ->Write();
	for (int i = 0; i < 9; i++) hOmegabarM_cen[i]->Write();
	hOmegabarp  ->Write();
	hOmegabarpt ->Write();
	hOmegabareta->Write();
	hOmegabary  ->Write();
	hOmegabarypt->Write();
	hOmegabarphi->Write();
	hOmegabarDL ->Write();
	hOmegabarDLminusLambdaDL->Write();

	hNumOmega->Write();
	hOmegaUsed->Write();
	hOmegaUsed_sideband->Write();

	// pT Spectrum
	hOmegaM_error_1->Write();
	hOmegabarM_error_1->Write();
	hOmegaM_error_2->Write();
	hOmegabarM_error_2->Write();
	hDauKaonM_error->Write();
	hDauLambdaM_error->Write();
	for (int i = 0; i < 9; i++) // centrality
	{
		for (int j = 0; j < num_pt_bin; j++) 
		{
			hOmegaM_pt[i][j]->Write();           hOmegabarM_pt[i][j]->Write();
			//hOmegaM_rotbkg_pi_pt[i]->Write(); hOmegabarM_rotbkg_pi_pt[i]->Write();
		}
	}

	// Omega v2
	for (int i = 0; i < 9; i++) // for (int i = 0; i < num_phi_bin; i++) // phi - psi bins
	{
		for (int j = 0; j < num_phi_bin; j++)
		{
			hOmegaM_phi[i][j]->Write();           hOmegabarM_phi[i][j]->Write();
		}
	}

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
	hDauLambdaDL    ->Write();

	hDauLambdabarM  ->Write();
	hDauLambdabarp  ->Write();
	hDauLambdabarpt ->Write();
	hDauLambdabary  ->Write();
	hDauLambdabarphi->Write();
	hDauLambdabarDL ->Write();

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

	hNumMixedEvent->Write();
	hTotalMixedEvent->Write();

	hgpdEdx      ->Write();
	hgdEdxErr    ->Write();
	hgpinvbeta   ->Write();
	hgm2         ->Write();
	hgpm2        ->Write();
	hgm2nSigmaKaon  ->Write();
	hgm2nSigmaPion  ->Write();
	hgm2nSigmaProton->Write();
	hgptnSigmaKaon->Write();
	hgptnSigmaPion->Write();
	hgptnSigmaProton->Write();
	hgptm2->Write();
	hgp          ->Write();
	hgpT         ->Write();
	hgDCAtoPV    ->Write();
	hgbtofYlocal ->Write();
	// hgKpdEdx       ->Write();
	// hgKpinvbeta    ->Write();
	hgKm2          ->Write();
	hgKpm2         ->Write();
	hgKp           ->Write();
	hgKpT          ->Write();
	hgKDCAtoPV     ->Write();
	hgKDCAtoO      ->Write();
	hgKpionpdEdx   ->Write();
	// hgptm2_largenSigmaKaon->Write();
	// hgptm2_smallnSigmaKaon->Write();
	hgnSigmaDiff->Write();
	hKaonCt->Write();
	hKpluspt_omega     ->Write();  
	hKpluspt_omegabar  ->Write();  
	hKminuspt_omega    ->Write();  
	hKminuspt_omegabar ->Write();
	hKpluseta_omega    ->Write();  
	hKpluseta_omegabar ->Write();  
	hKminuseta_omega   ->Write();  
	hKminuseta_omegabar->Write();
	hKplusphi_omega    ->Write();  
	hKplusphi_omegabar ->Write();  
	hKminusphi_omega   ->Write();  
	hKminusphi_omegabar->Write();
	/*
	for (int i = 0; i < 15; i++)
	{	
		hgPID2D_pt[i]->Write();
		hgzTPC_pt[i]->Write();
	}
	*/

	return;
}

//-----------------------------------------------------------------------------
void StKFParticleAnalysisMaker::WriteTrees()
{	

	for (int i = 0; i < num_mult_bin; i++)
	{
		for (int j = 0; j < num_vz_bin; j++)
		{
			for (int k = 0; k < num_EP_bin; k++)
			{
				omega_mix[i][j][k]->Write();
			}
		}
	}
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
	//time_t time_start;
	//time_t time_now;
	//time(&time_start);

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
	//if(fabs(VertexZ-vpdVz)>3.) return kStOK;       // no vpd cut in low energy? Yes!!

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

	// EPD Event plane
	/*
	for (int i = 0; i < mPicoDst->numberOfEpdHits(); i++) 
	{	
		TClonesArray *middle_step = mPicoDst->picoArray(8); // it is the eigth array, confirmed!!!
		cout << ((StPicoEpdHit*)(*middle_step)[i])->adc() << endl; // now this also works
		//cout << ((StPicoEpdHit*)(mPicoDst->epdHit(i)))->adc() << endl; // this works!!!
		(*mEpdHits)[i] = (StPicoEpdHit*)(mPicoDst->epdHit(i)); 
	}
	*/
	//mEpdHits = mPicoDst->picoArray(8); // grab TClonesArray directly?
	/*
	for (int i = 0; i < mPicoDst->numberOfEpdHits(); i++) 
	{
		cout << "ADC = " << ((StPicoEpdHit*)(*mEpdHits)[i])->adc() << endl;
		cout << "side = " << ((StPicoEpdHit*)(*mEpdHits)[i])->side() << endl;
		cout << "nMIP = " << ((StPicoEpdHit*)(*mEpdHits)[i])->nMIP() << endl;
		cout << "position = " << ((StPicoEpdHit*)(*mEpdHits)[i])->position() << endl;
		cout << "tile = " << ((StPicoEpdHit*)(*mEpdHits)[i])->tile() << endl;
	}
	*/
	StEpdEpInfo result = mEpFinder->Results(mPicoDst->picoArray(StPicoArrays::EpdHit), Vertex3D, cent-1);
	if (cent == cen_cut)
	{
		hEPD_e_EP_1->Fill(result.EastPhiWeightedAndShiftedPsi(1));
		hEPD_w_EP_1->Fill(result.WestPhiWeightedAndShiftedPsi(1));
		hEPD_e_EP_2->Fill(result.EastPhiWeightedAndShiftedPsi(2));
		hEPD_w_EP_2->Fill(result.WestPhiWeightedAndShiftedPsi(2));
		hEPD_full_EP_1->Fill(result.FullPhiWeightedAndShiftedPsi(1));
		hEPD_full_EP_2->Fill(result.FullPhiWeightedAndShiftedPsi(2));
		for (int order = 1; order <= 3; order++) 
			hEPD_ew_cos->Fill(order*1.0, cos(order*1.0*(result.EastPhiWeightedAndShiftedPsi(order)-result.WestPhiWeightedAndShiftedPsi(order))));
	}
	float EP2_EPD_full = result.FullPhiWeightedAndShiftedPsi(2);

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
	hRefMultCorr_cent[cent-1]->Fill(mult_corr);
	///////////////

// ======= KFParticle ======= //
	vector<StLambdaDecayPair> KFParticleLambdaDecayPair;

	// centrality cut and mixed event index
	if (CutCent && cent != cen_cut) return kStOK;

	SetupKFParticle();
	if (InterfaceCantProcessEvent) return kStOK;

	//time(&time_now);
	//int time_diff = (int)difftime(time_now, time_start);
	//cout << "Up to finishing KFP: " << time_diff << endl;

	// collect lambdas and omegas
	std::vector<KFParticle> OmegaVec, OmegaSidebandVec, OmegaVecAll;
	std::vector<int> OmegaDauPionIdVec;
	std::vector<int> OmegaDauProtonIdVec;
	mix_px.resize(0); mix_py.resize(0); mix_pz.resize(0); mix_charge.resize(0); 
	for (int iKFParticle=0; iKFParticle < KFParticlePerformanceInterface->GetNReconstructedParticles(); iKFParticle++)
	{ 
		const KFParticle particle = KFParticleInterface->GetParticles()[iKFParticle]; 
		bool IsOmega = false, IsLambda = false;
		if (fabs(particle.GetPDG()) == OmegaPdg) IsOmega = true; else if (fabs(particle.GetPDG()) == LambdaPdg) IsLambda = true; else continue;
		int upQ; if (particle.GetPDG() > 0) upQ = 1; else if (particle.GetPDG() < 0) upQ = -1; else continue;
		if (IsOmega) OmegaVecAll.push_back(particle);
		if (IsOmega && isGoodOmega(cent, particle)) 
		{
			OmegaVec.push_back(particle);
			if (upQ ==  1) hOmegaUsed->Fill(0.);
			if (upQ == -1) hOmegaUsed->Fill(1.);
		}
		if (IsOmega && isSidebandOmega(cent, particle)) 
		{
			OmegaSidebandVec.push_back(particle);
			if (upQ ==  1) hOmegaUsed_sideband->Fill(0.);
			if (upQ == -1) hOmegaUsed_sideband->Fill(1.);
		}

		// Omega/lambda QA
		if (IsOmega)
		{
			if (upQ == 1)
			{
				hOmegaM  ->Fill(particle.GetMass());
				hOmegaM_cen[cent-1]->Fill(particle.GetMass());

				// pT spectrum
				if (fabs(particle.GetRapidity()) < 0.5)
				{
					for (int i = 0; i < num_pt_bin; i++)
					{
						if (particle.GetPt() >= StKFParticleAnalysisMaker::OmegaMassPtLowerBin[i] && 
							particle.GetPt() <= StKFParticleAnalysisMaker::OmegaMassPtLowerBin[i+1])
							hOmegaM_pt[cent-1][i]->Fill(particle.GetMass());
					}
				}

				// v2 vs pT
				// if (cent == cen_cut)
				// {
				// 	for (int i = 0; i < num_pt_bin; i++)
				// 	{
				// 		float phi_diff = fabs(particle.GetPhi() + pi - EP2_EPD_full);
				// 		if (phi_diff > pi) phi_diff = phi_diff - pi;
				// 		if (phi_diff > pi / 2.) phi_diff = pi - phi_diff;
				// 		int phi_bin = static_cast<int>(floor(phi_diff / (pi / 2. / num_phi_bin)));
				// 		if (particle.GetPt() >= StKFParticleAnalysisMaker::OmegaMassPtLowerBin[i] && 
				// 			particle.GetPt() <= StKFParticleAnalysisMaker::OmegaMassPtLowerBin[i+1])
				// 			hOmegaM_phi[phi_bin][i]->Fill(particle.GetMass());
				// 	}
				// }

				// v2 vs centrality
				float phi_diff = fabs(particle.GetPhi() + pi - EP2_EPD_full);
				if (phi_diff > pi) phi_diff = phi_diff - pi;
				if (phi_diff > pi / 2.) phi_diff = pi - phi_diff;
				int phi_bin = static_cast<int>(floor(phi_diff / (pi / 2. / num_phi_bin)));
				hOmegaM_phi[cent-1][phi_bin]->Fill(particle.GetMass());

				//if (!isGoodOmega(cent, particle)) continue; // subject to change
				hOmegay  ->Fill(particle.GetRapidity());
				hOmegaypt->Fill(particle.GetPt(), particle.GetRapidity());
				// hOmegaDL ->Fill(particle.GetDecayLength());
				
				// helix
				TVector3 pOmega(particle.GetPx(), particle.GetPy(), particle.GetPz());
				TVector3 xOmega(particle.GetX(), particle.GetY(), particle.GetZ());
				TLorentzVector OmegaLorentz(pOmega, particle.GetE());
				StPicoPhysicalHelix helixOmega(pOmega, xOmega, magnet*kilogauss, particle.GetQ());
				double pathlength = helixOmega.pathLength(Vertex3D, false);
				TVector3 pOmega_tb = helixOmega.momentumAt(pathlength, magnet*kilogauss); 
				hOmegap  ->Fill(pOmega_tb.Mag());
				hOmegapt ->Fill(pOmega_tb.Pt());
				hOmegaeta->Fill(pOmega_tb.Eta());
				hOmegaphi->Fill(pOmega_tb.Phi());

				if (StoringTree)
				{
					mix_px.push_back(pOmega_tb.X());
					mix_py.push_back(pOmega_tb.Y());
					mix_pz.push_back(pOmega_tb.Z());
					mix_charge.push_back(particle.GetQ());
					mix_evt_id = evtID;
					mix_run_id = runID;
					// omega_mix[mult_index][vz_index][EP_index]->Fill();
				}

				hOmegaDL->Fill((xOmega-Vertex3D).Mag());

				// daughters QA
				TLorentzVector dauKaon, dauLambda;
				TLorentzVector dauKaonPico, dauLambdaPico;
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
					TVector3 pDaughterPico = daughterTrack->gMom();
					TVector3 xDaughter(daughter.GetX(), daughter.GetY(), daughter.GetZ());
					TVector3 xDaughterPico = daughterTrack->origin();
					
					if (daughter.GetPDG() == -KaonPdg)
					{
						//dauKaon.SetVectM(pDaughter, KaonPdgMass);
						hDauKminusp  ->Fill(daughter.GetMomentum());
						hDauKminuspt ->Fill(daughter.GetPt());
						hDauKminusy  ->Fill(daughter.GetRapidity());
						hDauKminusphi->Fill(daughter.GetPhi());
						hDauKminusnSigma->Fill(daughterTrack->nSigmaKaon());

						StPicoPhysicalHelix helixDaughter(pDaughterPico, xDaughterPico, magnet*kilogauss, daughter.GetQ());
						pair<double, double> tmps = helixDaughter.pathLengths(helixOmega);
						TVector3 ox1 = helixDaughter.at(tmps.first); 
						TVector3 ox2 = helixOmega.at(tmps.second);
						daughter.TransportToProductionVertex();
						dauKaon.SetXYZT(daughter.GetPx(), daughter.GetPy(), daughter.GetPz(), daughter.GetE());
						hDauKaonM_error->Fill(dauKaon.M() - KaonPdgMass);
						dauKaon.SetXYZM(daughter.GetPx(), daughter.GetPy(), daughter.GetPz(), KaonPdgMass);
						dauKaonPico.SetVectM(ox1, KaonPdgMass); // use all StPicoTrack information
						//dauKaon.SetVectM(ox1, KaonPdgMass); // transport daughter to Omega DCA
						double dca = (ox1 - ox2).Mag();
						CutDecider(particle, hDCAOtoK_signal, hDCAOtoK_sideband, dca);
					}
					else 
					{
						//dauLambda.SetVectM(pDaughter, LambdaPdgMass);
						hDauLambdaM  ->Fill(daughter.GetMass());
						hDauLambdap  ->Fill(daughter.GetMomentum());
						hDauLambdapt ->Fill(daughter.GetPt());
						hDauLambdaphi->Fill(daughter.GetPhi());
						hDauLambday  ->Fill(daughter.GetRapidity());
						hDauLambdaDL ->Fill((xDaughter-Vertex3D).Mag());
						hOmegaDLminusLambdaDL->Fill((xDaughter-Vertex3D).Mag()-(xOmega-Vertex3D).Mag());

						StPicoPhysicalHelix helixDaughter(pDaughterPico*(1./pDaughterPico.Perp())*100, xDaughterPico, magnet*kilogauss, 1);
						pair<double, double> tmps = helixOmega.pathLengths(helixDaughter);
						TVector3 ox1 = helixOmega.at(tmps.first);
						TVector3 ox2 = helixDaughter.at(tmps.second);
						daughter.TransportToProductionVertex();
						dauLambda.SetXYZT(daughter.GetPx(), daughter.GetPy(), daughter.GetPz(), daughter.GetE());
						hDauLambdaM_error->Fill(dauLambda.M() - LambdaPdgMass);
						dauLambdaPico.SetVectM(ox2, LambdaPdgMass);
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
				hOmegaM_error_1->Fill(particle.GetMass()-(dauKaon+dauLambda).M());
				hOmegaM_error_2->Fill(particle.GetMass()-OmegaLorentz.M());
				dauKaon.RotateZ(pi);

				// for (int i = 0; i < num_pt_bin; i++)
				// {
				// 	if (particle.GetPt() >= StKFParticleAnalysisMaker::OmegaMassPtLowerBin[i] && 
				// 		particle.GetPt() <= StKFParticleAnalysisMaker::OmegaMassPtLowerBin[i+1])
				// 		hOmegaM_rotbkg_pi_pt[i]->Fill((dauKaon+dauLambda).M());
				// }
			}
			else
			{
				hOmegabarM  ->Fill(particle.GetMass());
				hOmegabarM_cen[cent-1]->Fill(particle.GetMass());

				// pT spectrum
				if (fabs(particle.GetRapidity()) < 0.5)
				{
					for (int i = 0; i < num_pt_bin; i++)
					{
						if (particle.GetPt() >= StKFParticleAnalysisMaker::OmegaMassPtLowerBin[i] && 
							particle.GetPt() <= StKFParticleAnalysisMaker::OmegaMassPtLowerBin[i+1])
							hOmegabarM_pt[cent-1][i]->Fill(particle.GetMass());
					}
				}

				// // v2
				// if (cent == cen_cut) // need to change to a wide centrality range next
				// {
				// 	for (int i = 0; i < num_pt_bin; i++)
				// 	{
				// 		float phi_diff = fabs(particle.GetPhi() + pi - EP2_EPD_full);
				// 		if (phi_diff > pi) phi_diff = phi_diff - pi;
				// 		if (phi_diff > pi / 2.) phi_diff = pi - phi_diff;
				// 		int phi_bin = static_cast<int>(floor(phi_diff / (pi / 2. / num_phi_bin)));
				// 		if (particle.GetPt() >= StKFParticleAnalysisMaker::OmegaMassPtLowerBin[i] && 
				// 			particle.GetPt() <= StKFParticleAnalysisMaker::OmegaMassPtLowerBin[i+1])
				// 			hOmegabarM_phi[phi_bin][i]->Fill(particle.GetMass());
				// 	}
				// }

				// v2 vs centrality
				float phi_diff = fabs(particle.GetPhi() + pi - EP2_EPD_full);
				if (phi_diff > pi) phi_diff = phi_diff - pi;
				if (phi_diff > pi / 2.) phi_diff = pi - phi_diff;
				int phi_bin = static_cast<int>(floor(phi_diff / (pi / 2. / num_phi_bin)));
				hOmegabarM_phi[cent-1][phi_bin]->Fill(particle.GetMass());

				//if (!isGoodOmega(cent, particle)) continue; // subject to change
				hOmegabary  ->Fill(particle.GetRapidity());
				hOmegabarypt->Fill(particle.GetPt(), particle.GetRapidity());
				// hOmegabarDL ->Fill(particle.GetDecayLength());

				// helix
				TVector3 pOmega(particle.GetPx(), particle.GetPy(), particle.GetPz());
				TVector3 xOmega(particle.GetX(), particle.GetY(), particle.GetZ());
				TLorentzVector OmegaLorentz(pOmega, particle.GetE());
				StPicoPhysicalHelix helixOmega(pOmega, xOmega, magnet*kilogauss, particle.GetQ());
				double pathlength = helixOmega.pathLength(Vertex3D, false);
				TVector3 pOmega_tb = helixOmega.momentumAt(pathlength, magnet*kilogauss); 
				hOmegabarp  ->Fill(pOmega_tb.Mag());
				hOmegabarpt ->Fill(pOmega_tb.Pt());
				hOmegabareta->Fill(pOmega_tb.Eta());
				hOmegabarphi->Fill(pOmega_tb.Phi());

				if (StoringTree)
				{
					mix_px.push_back(pOmega_tb.X());
					mix_py.push_back(pOmega_tb.Y());
					mix_pz.push_back(pOmega_tb.Z());
					mix_charge.push_back(particle.GetQ());
					mix_evt_id = evtID;
					mix_run_id = runID;
					// omega_mix[mult_index][vz_index][EP_index]->Fill();
				}

				// daughter kaon QA
				TLorentzVector dauKaon, dauLambda;
				TLorentzVector dauKaonPico, dauLambdaPico;
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
					TVector3 pDaughterPico = daughterTrack->gMom();
					TVector3 xDaughter(daughter.GetX(), daughter.GetY(), daughter.GetZ());
					TVector3 xDaughterPico = daughterTrack->origin();

					if (daughter.GetPDG() ==  KaonPdg)
					{
						//dauKaon.SetVectM(pDaughter, KaonPdgMass);
						hDauKplusp  ->Fill(daughter.GetMomentum());
						hDauKpluspt ->Fill(daughter.GetPt());
						hDauKplusy  ->Fill(daughter.GetRapidity());
						hDauKplusphi->Fill(daughter.GetPhi());
						hDauKplusnSigma->Fill(daughterTrack->nSigmaKaon());

						StPicoPhysicalHelix helixDaughter(pDaughterPico, xDaughterPico, magnet*kilogauss, daughter.GetQ());
						pair<double, double> tmps = helixDaughter.pathLengths(helixOmega);
						TVector3 ox1 = helixDaughter.at(tmps.first);
						TVector3 ox2 = helixOmega.at(tmps.second);
						daughter.TransportToProductionVertex();
						dauKaon.SetXYZT(daughter.GetPx(), daughter.GetPy(), daughter.GetPz(), daughter.GetE());
						hDauKaonM_error->Fill(dauKaon.M() - KaonPdgMass);
						dauKaon.SetXYZM(daughter.GetPx(), daughter.GetPy(), daughter.GetPz(), KaonPdgMass);
						dauKaonPico.SetVectM(ox1, KaonPdgMass);
						//dauKaon.SetVectM(ox1, KaonPdgMass);
						double dca = (ox1 - ox2).Mag();
						CutDecider(particle, hDCAOtoK_signal, hDCAOtoK_sideband, dca);
					}
					else 
					{
						//dauLambda.SetVectM(pDaughter, LambdaPdgMass);
						hDauLambdabarM  ->Fill(daughter.GetMass());
						hDauLambdabarp  ->Fill(daughter.GetMomentum());
						hDauLambdabarpt ->Fill(daughter.GetPt());
						hDauLambdabarphi->Fill(daughter.GetPhi());
						hDauLambdabary  ->Fill(daughter.GetRapidity());
						hDauLambdabarDL ->Fill((xDaughter-Vertex3D).Mag());
						hOmegabarDLminusLambdaDL->Fill((xDaughter-Vertex3D).Mag()-(xOmega-Vertex3D).Mag());

						StPicoPhysicalHelix helixDaughter(pDaughterPico*(1./pDaughterPico.Perp())*100, xDaughterPico, magnet*kilogauss, 1);
						pair<double, double> tmps = helixOmega.pathLengths(helixDaughter);
						TVector3 ox1 = helixOmega.at(tmps.first);
						TVector3 ox2 = helixDaughter.at(tmps.second);
						daughter.TransportToProductionVertex();
						dauLambda.SetXYZT(daughter.GetPx(), daughter.GetPy(), daughter.GetPz(), daughter.GetE());
						hDauLambdaM_error->Fill(dauLambda.M() - LambdaPdgMass);
						dauLambdaPico.SetVectM(ox2, LambdaPdgMass);
						//dauLambda.SetVectM(ox2, LambdaPdgMass);
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
				hOmegabarM_error_1->Fill(particle.GetMass()-(dauKaon+dauLambda).M());
				hOmegabarM_error_2->Fill(particle.GetMass()-OmegaLorentz.M());
				dauKaon.RotateZ(pi);

				// for (int i = 0; i < num_pt_bin; i++)
				// {
				// 	if (particle.GetPt() >= StKFParticleAnalysisMaker::OmegaMassPtLowerBin[i] && 
				// 		particle.GetPt() <= StKFParticleAnalysisMaker::OmegaMassPtLowerBin[i+1])
				// 		hOmegabarM_rotbkg_pi_pt[i]->Fill((dauKaon+dauLambda).M());
				// }
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
	hNumOmega->Fill(OmegaVec.size());

	// correlation function loop  
	if (!PerformAnalysis) return kStOK;
	my_event current_event(OmegaVec);
  	Int_t nTracks = mPicoDst->numberOfTracks();
	std::vector<int> kaon_tracks; kaon_tracks.resize(0);
	float Qx2w = 0, Qy2w = 0, Qx2e = 0, Qy2e = 0; // for EP 
	std::vector<int> track_index;
	int pct = 0; int pbct = 0; // counting protons
	int kpct = 0; int kmct = 0; // counting kaons
	for (Int_t iTrack = 0; iTrack < nTracks; iTrack++) 
	{
    	StPicoTrack *track = mPicoDst->track(iTrack);
    	if (! track)            continue;
    	if (! track->charge())  continue;
    	if (  track->nHitsFit() < 15) continue;
		if (  track->nHitsDedx() < 15) continue;
		if (  track->nHitsFit()*1.0 / track->nHitsMax() < 0.52) continue;
		if (  track->dEdxError() < 0.04 || track->dEdxError() > 0.12) continue; // same as kfp
		if (! track->isPrimary()) continue;
		track_index.push_back(iTrack);

		// track info
		float p = track->gMom().Mag();
		float pt = track->gMom().Perp();
		float phi = track->gMom().Phi();
		float phi_shifted = phi;
		float eta = track->gMom().Eta();
		float dcatopv = track->gDCA(Vertex3D).Mag();
		float nSigmaKaon = track->nSigmaKaon();
		float nSigmaPion = track->nSigmaPion();
		float nSigmaProton = track->nSigmaProton();

		// phi shift for asso
		for (int i = 1; i <= shift_order_asso; i++) // fill correction for output
		{
			hTPCAssoShiftOutput_cos->Fill(i, cent, cos(phi), mWght);
			hTPCAssoShiftOutput_sin->Fill(i, cent, sin(phi), mWght);
		}
		if (hTPCAssoShiftInput_cos != 0) // attempt to read input
		{
			for (int i=1; i <= shift_order_asso; i++)
			{
				float cosAve = hTPCAssoShiftInput_cos->GetBinContent(i, cent);
				float sinAve = hTPCAssoShiftInput_sin->GetBinContent(i, cent);
				phi_shifted += 2.0*(cosAve*sin(phi) - sinAve*cos(phi));
			}
		}
		hTPCAssoPhi->Fill(phi, mWght);
		hTPCAssoPhi_shifted->Fill(phi_shifted, mWght);
		if (eta > TPCAssoEtaCut) // construct east EP
		{
			Qx2e += pt*cos(2*phi_shifted);
			Qy2e += pt*sin(2*phi_shifted);
		}
		else if (eta < -1.0*TPCAssoEtaCut) //construct west EP
		{
			Qx2w += pt*cos(2*phi_shifted);
			Qy2w += pt*sin(2*phi_shifted);
		}

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
		hgp       ->Fill(p);
		hgpT      ->Fill(pt);
		hgDCAtoPV ->Fill(dcatopv);
		
		int ptbin = static_cast<int>(floor(pt/0.2));
		double zTPC = TMath::Log(track->dEdx() / 1e6 / StdEdxPull::EvalPred(pkaon.Mag()/KaonPdgMass,1,1)); 
		hgnSigmaDiff->Fill(zTPC / track->dEdxError() - nSigmaKaon);
		hgptnSigmaKaon->Fill(pt, nSigmaKaon);
		hgptnSigmaPion->Fill(pt, nSigmaPion);
		hgptnSigmaProton->Fill(pt, nSigmaProton);
		
		// if (ptbin >= 0 && ptbin <= 14) hgzTPC_pt[ptbin]->Fill(track->nSigmaKaon());
		if (hasTOF)
		{
			beta = (mPicoDst->btofPidTraits(tofindex))->btofBeta();
			m2 = pkaon.Mag2()*(1.0 / beta / beta - 1.0);
			hgpinvbeta->Fill(pkaon.Mag(), 1./beta);
			hgm2  ->Fill(m2);
			hgpm2 ->Fill(pkaon.Mag(), m2);
			hgptm2->Fill(pt, m2);
			hgm2nSigmaKaon  ->Fill(nSigmaKaon, m2);
			hgm2nSigmaPion  ->Fill(nSigmaPion, m2);
			hgm2nSigmaProton->Fill(nSigmaProton, m2);

			// some kaon QA
			//if (track->nSigmaKaon() >  6) hgptm2_largenSigmaKaon->Fill(track->gMom().Perp(), m2);
			//if (track->nSigmaKaon() < -6) hgptm2_smallnSigmaKaon->Fill(track->gMom().Perp(), m2);
			double zTOF = 1/beta - sqrt(KaonPdgMass*KaonPdgMass/pkaon.Mag2()+1);
			// if (ptbin >= 0 && ptbin <= 14) hgPID2D_pt[ptbin]->Fill(nSigmaKaon, zTOF);
		}

		// primary proton cut for coalescence test
		bool proton_cut = true;
		if (pt < 0.42 || pt > 1.8) proton_cut = false; // use p < 2
		if (!hasTOF && pt > 0.9) proton_cut = false;
		if (pt > 0.9 && (m2 > 1.1 || m2 < 0.75)) proton_cut = false;
		if (fabs(nSigmaProton) > 3) proton_cut = false;
		if (dcatopv > 2) proton_cut = false;
		if (proton_cut)
		{	
			TLorentzVector lv_proton; lv_proton.SetVectM(track->gMom(), ProtonPdgMass);
			if (track->charge() > 0) {hProtony    ->Fill(lv_proton.Rapidity()); if (fabs(lv_proton.Rapidity()) < 0.5) pct++; }
			else					 {hAntiProtony->Fill(lv_proton.Rapidity()); if (fabs(lv_proton.Rapidity()) < 0.5) pbct++;}
			if (track->charge() > 0) hproton_EPD_v2->Fill(cent, cos(2.*(phi-EP2_EPD_full)), pt);
			else 					 hantiproton_EPD_v2->Fill(cent, cos(2.*(phi-EP2_EPD_full)), pt);
		}

		// primary pion cut for coalescence test
		bool pion_cut = true;
		if (pt < 0.28 || pt > 1.2) pion_cut = false; // use p < 2
		if (!hasTOF && pt > 0.6) pion_cut = false;
		if (pt > 0.6 && (m2 > 0.1 || m2 < -0.1)) pion_cut = false;
		if (fabs(nSigmaPion) > 3) pion_cut = false;
		if (dcatopv > 2) pion_cut = false;
		if (pion_cut)
		{
			if (track->charge() > 0) hpiplus_EPD_v2->Fill(cent, cos(2.*(phi-EP2_EPD_full)), pt);
			else 					 hpiminus_EPD_v2->Fill(cent, cos(2.*(phi-EP2_EPD_full)), pt);
		}

		// primary kaon cut
		/******** looser cut ********/
		bool kaon_cut = true;
		if (pt < 0.15 || pt > 1.6) kaon_cut = false; // use p < 1.6
		if (!hasTOF && pt > 0.4) kaon_cut = false;
		if (pt > 0.4 && (m2 > 0.34 || m2 < 0.15)) kaon_cut = false;
		double zTOF = 1/beta - sqrt(KaonPdgMass*KaonPdgMass/pkaon.Mag2()+1);
		KaonPID decider(zTOF, nSigmaKaon, pt);
		if (!decider.IsKaonSimple(3.)) kaon_cut = false;
		if (dcatopv > 2) kaon_cut = false;
		bool isDaughter = false;
		for (int i = 0; i < OmegaVec.size(); i++) if (IsKaonOmegaDaughter(OmegaVec[i], track->id())) isDaughter = true;
		if (isDaughter) kaon_cut = false;

		if (kaon_cut)
		{
			if (track->charge() > 0) kpct++; 
			else                     kmct++;
			kaon_tracks.push_back(iTrack);
		}
		
		/******** stricter cut ********/
		/*
		if (pt < 0.15 || pt> 1.6) continue; // use p < 1.6
		if (!hasTOF) continue;pt
		double zTOF = 1/beta - sqrt(KaonPdgMass*KaonPdgMass/pkaon.Mag2()+1);
		KaonPID decider(zTOF, track->nSigmaKaon(), track->gMom().Perp());
		if (!decider.IsKaon()) continue;
		*/

		// kaon QA and correlation
		if (kaon_cut)
		{
			hgKpdEdx    ->Fill(pkaon.Mag(), track->dEdx());
			hgKp       ->Fill(p);
			hgKpT      ->Fill(pt);
			hgKDCAtoPV ->Fill(dcatopv);
			if (OmegaVec.size() == 1 && OmegaVec[0].GetQ() < 0) 
			{	
				if (track->charge() > 0)
				{
					hKpluspt_omega ->Fill(pt);
					hKpluseta_omega->Fill(eta);
					hKplusphi_omega->Fill(phi);
				}
				else
				{
					hKminuspt_omega ->Fill(pt);
					hKminuseta_omega->Fill(eta);
					hKminusphi_omega->Fill(phi);
				}
			}
			if (OmegaVec.size() == 1 && OmegaVec[0].GetQ() > 0) 
			{
				if (track->charge() > 0)
				{
					hKpluspt_omegabar ->Fill(pt);
					hKpluseta_omegabar->Fill(eta);
					hKplusphi_omegabar->Fill(phi);
				}
				else
				{
					hKminuspt_omegabar ->Fill(pt);
					hKminuseta_omegabar->Fill(eta);
					hKminusphi_omegabar->Fill(phi);
				}
			}
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
			const int kaonindex = track->id();
			for (int iOmega=0; iOmega < OmegaVec.size(); iOmega++)
			{ 
				const KFParticle particle = OmegaVec[iOmega]; 
				if (IsKaonOmegaDaughter(particle, kaonindex)) continue;

				// pair-wise should be added after this line
				/* */

				// Omega momentum at DCA to PV
				TVector3 pOmega(particle.GetPx(), particle.GetPy(), particle.GetPz());
				TVector3 xOmega(particle.GetX(), particle.GetY(), particle.GetZ());
				StPicoPhysicalHelix helixOmega(pOmega, xOmega, magnet*kilogauss, particle.GetQ());
				double pathlength = helixOmega.pathLength(Vertex3D, false);
				TVector3 pOmega_tb = helixOmega.momentumAt(pathlength, magnet*kilogauss); 

				// k*
				TLorentzVector lv1; lv1.SetVectM(pOmega_tb, OmegaPdgMass);
				TLorentzVector lv2; lv2.SetVectM(track->gMom(), KaonPdgMass);
				double dpt = fabs(lv1.Perp()-lv2.Perp());
				double dy  = fabs(lv1.Rapidity() - lv2.Rapidity());
				double dphi = fabs(lv1.Vect().DeltaPhi(lv2.Vect()));
				TLorentzVector P = lv1 + lv2;
				TVector3 pair_beta = P.BoostVector();
				lv1.Boost((-1)*pair_beta); 	
				lv2.Boost((-1)*pair_beta); 	
				double kstar = 0.5*(lv1-lv2).Vect().Mag();

				// pT weight if anti-omega
				float weight = 1.;
				if (PtReweighting) weight = GetPtWeight(particle);
				
				if (track->charge() > 0 && particle.GetQ() < 0) 
				{
					hCorrKplusO   ->Fill(kstar);
					hPtCorrKplusO ->Fill(dpt);
					hyCorrKplusO  ->Fill(dy);
					hphiCorrKplusO->Fill(dphi);

					hCorrKplusO_y_pT  ->Fill(dpt, dy);
					hCorrKplusO_y_phi ->Fill(dphi, dy);
					hCorrKplusO_phi_pT->Fill(dpt, dphi);
				}
				if (track->charge() > 0 && particle.GetQ() > 0) 
				{
					hCorrKplusObar   ->Fill(kstar, weight);
					hPtCorrKplusObar ->Fill(dpt, weight);
					hyCorrKplusObar  ->Fill(dy, weight);
					hphiCorrKplusObar->Fill(dphi, weight);

					hCorrKplusObar_y_pT  ->Fill(dpt, dy, weight);
					hCorrKplusObar_y_phi ->Fill(dphi, dy, weight);
					hCorrKplusObar_phi_pT->Fill(dpt, dphi, weight);
				}
				if (track->charge() < 0 && particle.GetQ() < 0)
				{
					hCorrKminusO   ->Fill(kstar);
					hPtCorrKminusO ->Fill(dpt);
					hyCorrKminusO  ->Fill(dy);
					hphiCorrKminusO->Fill(dphi);
					if (dpt < 0.5) hNegPtDiff_dphi_KmO->Fill(dphi);
					if (dpt > 1.0) hPosPtDiff_dphi_KmO->Fill(dphi);

					hCorrKminusO_y_pT  ->Fill(dpt, dy);
					hCorrKminusO_y_phi ->Fill(dphi, dy);
					hCorrKminusO_phi_pT->Fill(dpt, dphi);
				}
				if (track->charge() < 0 && particle.GetQ() > 0) 
				{
					hCorrKminusObar   ->Fill(kstar, weight);
					hPtCorrKminusObar ->Fill(dpt, weight);
					hyCorrKminusObar  ->Fill(dy, weight);
					hphiCorrKminusObar->Fill(dphi, weight);
					if (dpt < 0.5) hNegPtDiff_dphi_KmOb->Fill(dphi, weight);
					if (dpt > 1.0) hPosPtDiff_dphi_KmOb->Fill(dphi, weight);

					hCorrKminusObar_y_pT  ->Fill(dpt, dy, weight);
					hCorrKminusObar_y_phi ->Fill(dphi, dy, weight);
					hCorrKminusObar_phi_pT->Fill(dpt, dphi, weight);
				}
			} // End loop over regular Omega

			// Omega sideband loop
			for (int iOmega=0; iOmega < OmegaSidebandVec.size(); iOmega++)
			{ 
				const KFParticle particle = OmegaSidebandVec[iOmega]; 
				if (IsKaonOmegaDaughter(particle, kaonindex)) continue;

				// pair-wise should be added after this line
				/* */

				// Omega momentum at DCA to PV
				TVector3 pOmega(particle.GetPx(), particle.GetPy(), particle.GetPz());
				TVector3 xOmega(particle.GetX(), particle.GetY(), particle.GetZ());
				StPicoPhysicalHelix helixOmega(pOmega, xOmega, magnet*kilogauss, particle.GetQ());
				double pathlength = helixOmega.pathLength(Vertex3D, false);
				TVector3 pOmega_tb = helixOmega.momentumAt(pathlength, magnet*kilogauss); 

				// k*
				TLorentzVector lv1; lv1.SetVectM(pOmega_tb, OmegaPdgMass);
				TLorentzVector lv2; lv2.SetVectM(track->gMom(), KaonPdgMass);
				double dpt = fabs(lv1.Perp()-lv2.Perp());
				double dy  = fabs(lv1.Rapidity() - lv2.Rapidity());
				double dphi = fabs(lv1.Vect().DeltaPhi(lv2.Vect()));
				TLorentzVector P = lv1 + lv2;
				TVector3 pair_beta = P.BoostVector();
				lv1.Boost((-1)*pair_beta); 	
				lv2.Boost((-1)*pair_beta); 	
				double kstar = 0.5*(lv1-lv2).Vect().Mag();

				// pT weight if anti-omega
				float weight = 1.;
				if (PtReweighting) weight = GetPtWeight(particle);
				
				if (track->charge() > 0 && particle.GetQ() < 0) 
				{
					hCorrKplusO_sideband   ->Fill(kstar);
					hPtCorrKplusO_sideband ->Fill(dpt);
					hyCorrKplusO_sideband  ->Fill(dy);
					hphiCorrKplusO_sideband->Fill(dphi);

					//hCorrKplusO_y_pT  ->Fill(dpt, dy);
					//hCorrKplusO_y_phi ->Fill(dphi, dy);
					//hCorrKplusO_phi_pT->Fill(dpt, dphi);
				}
				if (track->charge() > 0 && particle.GetQ() > 0) 
				{
					hCorrKplusObar_sideband   ->Fill(kstar, weight);
					hPtCorrKplusObar_sideband ->Fill(dpt, weight);
					hyCorrKplusObar_sideband  ->Fill(dy, weight);
					hphiCorrKplusObar_sideband->Fill(dphi, weight);

					//hCorrKplusObar_y_pT  ->Fill(dpt, dy);
					//hCorrKplusObar_y_phi ->Fill(dphi, dy);
					//hCorrKplusObar_phi_pT->Fill(dpt, dphi);
				}
				if (track->charge() < 0 && particle.GetQ() < 0)
				{
					hCorrKminusO_sideband   ->Fill(kstar);
					hPtCorrKminusO_sideband ->Fill(dpt);
					hyCorrKminusO_sideband  ->Fill(dy);
					hphiCorrKminusO_sideband->Fill(dphi);
					//if (dpt < 0.5) hNegPtDiff_dphi_KmO->Fill(dphi);
					//if (dpt > 1.0) hPosPtDiff_dphi_KmO->Fill(dphi);

					//hCorrKminusO_y_pT  ->Fill(dpt, dy);
					//hCorrKminusO_y_phi ->Fill(dphi, dy);
					//hCorrKminusO_phi_pT->Fill(dpt, dphi);
				}
				if (track->charge() < 0 && particle.GetQ() > 0) 
				{
					hCorrKminusObar_sideband   ->Fill(kstar, weight);
					hPtCorrKminusObar_sideband ->Fill(dpt, weight);
					hyCorrKminusObar_sideband  ->Fill(dy, weight);
					hphiCorrKminusObar_sideband->Fill(dphi, weight);
					//if (dpt < 0.5) hNegPtDiff_dphi_KmOb->Fill(dphi);
					//if (dpt > 1.0) hPosPtDiff_dphi_KmOb->Fill(dphi);

					//hCorrKminusObar_y_pT  ->Fill(dpt, dy);
					//hCorrKminusObar_y_phi ->Fill(dphi, dy);
					//hCorrKminusObar_phi_pT->Fill(dpt, dphi);
				}
			} // End loop over sideband Omega
		}
	}

	// // TPC resolution shuffling
	// float Qx1_e = 0, Qy1_e = 0, Qx2_e = 0, Qy2_e = 0; // for EP res
	// float Qx1_w = 0, Qy1_w = 0, Qx2_w = 0, Qy2_w = 0; // for EP res
	// TVector2 Q1_e, Q2_e, Q1_w, Q2_w; // for EP res
	// std::random_shuffle(track_index.begin(), track_index.end());
	// for (int i = 0; i < track_index.size(); i++) 
	// {
    // 	StPicoTrack *track = mPicoDst->track(track_index[i]);
	// 	float pt = track->gMom().Perp();
	// 	float phi = track->gMom().Phi();

	// 	if (i <= track_index.size() / 2) 
	// 	{
	// 		Qx1_w += pt*cos(phi);
	// 		Qy1_w += pt*sin(phi);
	// 		Qx2_w += pt*cos(2*phi);
	// 		Qy2_w += pt*sin(2*phi);
	// 	}
	// 	else
	// 	{
	// 		Qx1_e += pt*cos(phi);
	// 		Qy1_e += pt*sin(phi);
	// 		Qx2_e += pt*cos(2*phi);
	// 		Qy2_e += pt*sin(2*phi);
	// 	}
	// }

	// event mixing stuff
	int mult_index = MixRefMultBin(cent, mult_corr);
	int vz_index   = static_cast<int>(floor(VertexZ/10.)+8.); if (vz_index == 16) vz_index--;
	int EP_index   = -1;

	// TPC EP
	TVector2 Q2w, Q2e, Q2full; // for EP
	float EP2_TPC_w = 0, EP2_TPC_e = 0, EP2_TPC_full = 0;
	Q2full.Set(Qx2e+Qx2w, Qy2e+Qy2w);
	Q2e.Set(Qx2e, Qy2e), Q2w.Set(Qx2w, Qy2w);
	EP2_TPC_full = Q2full.Phi() / 2.;
	EP2_TPC_e = Q2e.Phi() / 2.;
	EP2_TPC_w = Q2w.Phi() / 2.;

	// shift correction
	float EP2_TPC[3] = {0., 0., 0.}; 
	float EP2_TPC_shifted[3] = {0., 0., 0.};
	EP2_TPC[0] = EP2_TPC_e; EP2_TPC[1] = EP2_TPC_w; EP2_TPC[2] = EP2_TPC_full;
	EP2_TPC_shifted[0] = EP2_TPC_e; EP2_TPC_shifted[1] = EP2_TPC_w; EP2_TPC_shifted[2] = EP2_TPC_full;
	for (int ewFull = 0; ewFull < 3; ewFull++)
	{
		for (int i = 1; i <= shift_order_EP; i++)
		{
			hTPCEPShiftOutput_cos[ewFull]->Fill(i, cent, cos(i*1.0*EP2_TPC[ewFull]), mWght);
			hTPCEPShiftOutput_sin[ewFull]->Fill(i, cent, sin(i*1.0*EP2_TPC[ewFull]), mWght);
		}
	}

	// float EP_1_shift_cos[4] = {-0.191905, 0.00567376, 0.003376, -0.000779899};
	// float EP_1_shift_sin[4] = {0.170762, -0.0436447, 0.00356307, -0.000438169};
	// float EP_2_shift_cos[4] = {-0.0447853, -0.0170472, 0.00125928, 0.000215389};
	// float EP_2_shift_sin[4] = {0.164891, -0.0109107, -0.000975253, 0.000114769};
	EP_index = static_cast<int>(EP2_TPC_full / (PI/6)); if (EP_index == 6) EP_index = 5; // for event mixing index

	// shift EPs
	for (int ewFull = 0; ewFull < 3; ewFull++)
	{
		if (hTPCEPShiftInput_cos[ewFull] != 0)
		{	
			for (int i = 1; i <= shift_order_EP; i++)
			{
				float cosAve = hTPCEPShiftInput_cos[ewFull]->GetBinContent(i, cent);
				float sinAve = hTPCEPShiftInput_sin[ewFull]->GetBinContent(i, cent);
				EP2_TPC_shifted[ewFull] += - 1./i * cos(i*2.*EP2_TPC[ewFull]) * sinAve + 1./i * sin(i*2*EP2_TPC[ewFull]) * cosAve;
			}
		}
		hTPCEP_2[ewFull]->Fill(EP2_TPC[ewFull], mWght);
		hTPCEP_2_shifted[ewFull]->Fill(EP2_TPC_shifted[ewFull], mWght);
	}
	float EP2_TPC_e_shifted = EP2_TPC_shifted[0];
	float EP2_TPC_w_shifted = EP2_TPC_shifted[1];
	float EP2_TPC_full_shifted = EP2_TPC_shifted[2]; 
	hTPCEP_ew_cos->Fill(2., cos(2.*EP2_TPC_e_shifted-2.*EP2_TPC_w_shifted), mWght);// resolution

	// Omega v2
	/* direct v_2, not very applicable */
	for (int i = 0; i < OmegaVec.size(); i++)
	{
		if (OmegaVec[i].GetQ() < 0) 
		{
			float EP2_TPC_Omega = OmegaVec[i].GetEta()>0? EP2_TPC_w_shifted: EP2_TPC_e_shifted; 
			hOmega_TPC_v2_pt->Fill(OmegaVec[i].GetPt(), cos(2*OmegaVec[i].GetPhi() - 2*EP2_TPC_Omega));
			hOmega_EPD_v2_pt->Fill(OmegaVec[i].GetPt(), cos(2*OmegaVec[i].GetPhi() - 2*EP2_EPD_full));
		}
		else if (OmegaVec[i].GetQ() > 0)
		{
			float EP2_TPC_Omega = OmegaVec[i].GetEta()>0? EP2_TPC_w_shifted: EP2_TPC_e_shifted; 
			hOmegabar_TPC_v2_pt->Fill(OmegaVec[i].GetPt(), cos(2*OmegaVec[i].GetPhi() - 2*EP2_TPC_Omega));
			hOmegabar_EPD_v2_pt->Fill(OmegaVec[i].GetPt(), cos(2*OmegaVec[i].GetPhi() - 2*EP2_EPD_full));
		}
	}
	/* v_2 vs M_inv method */
	for (int i = 0; i < OmegaVecAll.size(); i++)
	{
		const KFParticle particle = OmegaVecAll[i];
		if (particle.GetQ() < 0) 
		{
			float EP2_TPC_Omega = particle.GetEta()>0? EP2_TPC_w_shifted: EP2_TPC_e_shifted;
			hOmega_TPC_v2[cent-1]->Fill(particle.GetMass(), cos(2.*particle.GetPhi() - 2.*EP2_TPC_Omega)); 
			hOmega_EPD_v2[cent-1]->Fill(particle.GetMass(), cos(2.*particle.GetPhi() - 2.*EP2_EPD_full));
		}
		else if (particle.GetQ() > 0)
		{
			float EP2_TPC_Omega = particle.GetEta()>0? EP2_TPC_w_shifted: EP2_TPC_e_shifted;
			hOmegabar_TPC_v2[cent-1]->Fill(particle.GetMass(), cos(2.*particle.GetPhi() - 2.*EP2_TPC_Omega)); 
			hOmegabar_EPD_v2[cent-1]->Fill(particle.GetMass(), cos(2.*particle.GetPhi() - 2.*EP2_EPD_full));
		}
	}
	// new observable
	if (OmegaVec.size() == 1 && OmegaVec[0].GetQ() < 0) hKratio_omega   ->Fill(pbct*1.0/pct, kmct*1.0/kpct);
	if (OmegaVec.size() == 0)                           hKratio_wo      ->Fill(pbct*1.0/pct, kmct*1.0/kpct);
	if (OmegaVec.size() == 1 && OmegaVec[0].GetQ() > 0) hKratio_omegabar->Fill(pbct*1.0/pct, kmct*1.0/kpct);

	// filling trees
	if (StoringTree && !OmegaVec.size() == 0) omega_mix[mult_index][vz_index][EP_index]->Fill();
	// if (StoringTree && !OmegaVec.size() == 0) omega_mix[0][0][0]->Fill(); // for testing

	// counting kaon
	if (!current_event.IsEmptyEvent()) hKaonCt->Fill(1.0, kaon_tracks.size());
	else                              {hKaonCt->Fill(0.0, kaon_tracks.size()); return kStOK;}

	// mixed event
	if (!PerformMixing) return kStOK;
	std::vector<float> *px_omega     = 0;
	std::vector<float> *py_omega     = 0;
	std::vector<float> *pz_omega     = 0;
	std::vector<int>   *charge_omega = 0;
	int evtid_omega = 0; int runid_omega = 0;
	omega_mix[mult_index][vz_index][EP_index]->SetBranchAddress("px", &px_omega);
	omega_mix[mult_index][vz_index][EP_index]->SetBranchAddress("py", &py_omega);
	omega_mix[mult_index][vz_index][EP_index]->SetBranchAddress("pz", &pz_omega);
	omega_mix[mult_index][vz_index][EP_index]->SetBranchAddress("charge", &charge_omega);
	omega_mix[mult_index][vz_index][EP_index]->SetBranchAddress("evt_id", &evtid_omega);
	omega_mix[mult_index][vz_index][EP_index]->SetBranchAddress("run_id", &runid_omega);
	for (int iMixEvent = 0; iMixEvent < omega_mix[mult_index][vz_index][EP_index]->GetEntries(); iMixEvent++)
	{	
		omega_mix[mult_index][vz_index][EP_index]->GetEntry(iMixEvent);
		if (runID == runid_omega && evtID == evtid_omega) continue; // no self-correlation
		for (int iOmega = 0; iOmega < px_omega->size(); iOmega++)
		{
			if (charge_omega->at(iOmega) < 0) hOmegaUsed->Fill(2.);
			if (charge_omega->at(iOmega) > 0) hOmegaUsed->Fill(3.);

			TVector3 p_omega(px_omega->at(iOmega), py_omega->at(iOmega), pz_omega->at(iOmega));

			for (int kaon_track = 0; kaon_track < kaon_tracks.size(); kaon_track++) 
			{
				StPicoTrack *track = mPicoDst->track(kaon_tracks[kaon_track]);
				
				// k*
				TLorentzVector lv1; lv1.SetVectM(p_omega,       OmegaPdgMass);
				TLorentzVector lv2; lv2.SetVectM(track->gMom(), KaonPdgMass);
				double dpt = fabs(lv1.Perp()-lv2.Perp());
				double dy  = fabs(lv1.Rapidity() - lv2.Rapidity());
				TLorentzVector P = lv1 + lv2;
				TVector3 pair_beta = P.BoostVector();
				lv1.Boost((-1)*pair_beta); 	
				lv2.Boost((-1)*pair_beta); 		
				if (track->charge() > 0 && charge_omega->at(iOmega) < 0) hCorrKplusO_mixed    ->Fill(0.5*(lv1-lv2).Vect().Mag());
				if (track->charge() > 0 && charge_omega->at(iOmega) > 0) hCorrKplusObar_mixed ->Fill(0.5*(lv1-lv2).Vect().Mag());
				if (track->charge() < 0 && charge_omega->at(iOmega) < 0) hCorrKminusO_mixed   ->Fill(0.5*(lv1-lv2).Vect().Mag());
				if (track->charge() < 0 && charge_omega->at(iOmega) > 0) hCorrKminusObar_mixed->Fill(0.5*(lv1-lv2).Vect().Mag());

				if (track->charge() > 0 && charge_omega->at(iOmega) < 0) hCorrKplusO_y_pT_mixed    ->Fill(dpt, dy);
				if (track->charge() > 0 && charge_omega->at(iOmega) > 0) hCorrKplusObar_y_pT_mixed ->Fill(dpt, dy);
				if (track->charge() < 0 && charge_omega->at(iOmega) < 0) hCorrKminusO_y_pT_mixed   ->Fill(dpt, dy);
				if (track->charge() < 0 && charge_omega->at(iOmega) > 0) hCorrKminusObar_y_pT_mixed->Fill(dpt, dy);
			}

		}
	}


	/*
	std::vector<my_event> mixed_events; mixed_events.resize(0);
	if (!buffer.IsEmpty(mult_index, VertexZ)) mixed_events = buffer.Sample_All(mult_index, VertexZ);
	hNumMixedEvent->Fill(mixed_events.size());
	for (int iMixEvent = 0; iMixEvent < mixed_events.size(); iMixEvent++)
	{
		std::vector<KFParticle> particles = mixed_events[iMixEvent].GetParticles();
		for (int iOmega = 0; iOmega < particles.size(); iOmega++)
		{
			const KFParticle particle = particles[iOmega];
			if (particle.GetQ() < 0) hOmegaUsed->Fill(2.);
			if (particle.GetQ() > 0) hOmegaUsed->Fill(3.);
			// Omega momentum at DCA to PV
			TVector3 pOmega(particle.GetPx(), particle.GetPy(), particle.GetPz());
			TVector3 xOmega(particle.GetX(), particle.GetY(), particle.GetZ());
			StPicoPhysicalHelix helixOmega(pOmega, xOmega, magnet*kilogauss, particle.GetQ());
			double pathlength = helixOmega.pathLength(Vertex3D, false);
			TVector3 pOmega_tb = helixOmega.momentumAt(pathlength, magnet*kilogauss); 

			for (int kaon_track = 0; kaon_track < kaon_tracks.size(); kaon_track++) 
			{
				StPicoTrack *track = mPicoDst->track(kaon_tracks[kaon_track]);
				
				// k*
				TLorentzVector lv1; lv1.SetVectM(pOmega_tb, OmegaPdgMass);
				TLorentzVector lv2; lv2.SetVectM(track->gMom(), KaonPdgMass);
				double dpt = fabs(lv1.Perp()-lv2.Perp());
				double dy  = fabs(lv1.Rapidity() - lv2.Rapidity());
				TLorentzVector P = lv1 + lv2;
				TVector3 pair_beta = P.BoostVector();
				lv1.Boost((-1)*pair_beta); 	
				lv2.Boost((-1)*pair_beta); 		
				if (track->charge() > 0 && particle.GetQ() < 0) hCorrKplusO_mixed    ->Fill(0.5*(lv1-lv2).Vect().Mag());
				if (track->charge() > 0 && particle.GetQ() > 0) hCorrKplusObar_mixed ->Fill(0.5*(lv1-lv2).Vect().Mag());
				if (track->charge() < 0 && particle.GetQ() < 0) hCorrKminusO_mixed   ->Fill(0.5*(lv1-lv2).Vect().Mag());
				if (track->charge() < 0 && particle.GetQ() > 0) hCorrKminusObar_mixed->Fill(0.5*(lv1-lv2).Vect().Mag());

				if (track->charge() > 0 && particle.GetQ() < 0) hCorrKplusO_y_pT_mixed    ->Fill(dpt, dy);
				if (track->charge() > 0 && particle.GetQ() > 0) hCorrKplusObar_y_pT_mixed ->Fill(dpt, dy);
				if (track->charge() < 0 && particle.GetQ() < 0) hCorrKminusO_y_pT_mixed   ->Fill(dpt, dy);
				if (track->charge() < 0 && particle.GetQ() > 0) hCorrKminusObar_y_pT_mixed->Fill(dpt, dy);
			}
		}
	}

	buffer.Add_Reservoir(current_event, mult_index, VertexZ);	
	hTotalMixedEvent->Fill(buffer.TotalStorage());
	*/

// ======= KFParticle end ======= //

	/////////////////////////////////////////////////////////
	return kStOK;

}

bool StKFParticleAnalysisMaker::isGoodOmega(int cent, KFParticle Omega)
{	
	int cent_bin = 0;
	if (cent >= 1 && cent <= 3) cent_bin = 0;
	else cent_bin = cent-3;
	if (Omega.GetQ() > 0) return (fabs(Omega.GetMass() - OmegaPdgMass) < 3*OmegabarMassSigma[cent_bin]);
	else                  return (fabs(Omega.GetMass() - OmegaPdgMass) < 3*OmegaMassSigma   [cent_bin]); 
}

bool StKFParticleAnalysisMaker::isSidebandOmega(int cent, KFParticle Omega)
{	
	int cent_bin = 0;
	if (cent >= 1 && cent <= 3) cent_bin = 0;
	else cent_bin = cent-3;
	if (Omega.GetQ() > 0) return (fabs(Omega.GetMass() - OmegaPdgMass) > 4*OmegabarMassSigma[cent_bin] && fabs(Omega.GetMass() - OmegaPdgMass) < 6*OmegabarMassSigma[cent_bin]);
	else                  return (fabs(Omega.GetMass() - OmegaPdgMass) > 4*OmegaMassSigma   [cent_bin] && fabs(Omega.GetMass() - OmegaPdgMass) < 6*OmegaMassSigma   [cent_bin]); 
}

int StKFParticleAnalysisMaker::MixRefMultBin(int cent, int refmult)
{	
	int bin = -1;
	int refmult_cut[] = {6, 13, 25, 44, 72, 113, 168, 246, 295}; //left inclusive, 1-9
	if (cent == 9) bin = static_cast<int>(floor((refmult-refmult_cut[8]) / 20));
	else bin = static_cast<int>(floor((refmult-refmult_cut[cent])/((refmult_cut[cent+1]-refmult_cut[cent])/5.))); 
	
	if (bin > 4) bin = 4;
	return bin;
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
	if (mass_diff < 2*0.0021) hist_signal->Fill(value);
	if (mass_diff > 4*0.0021 && mass_diff < 6*0.0021) hist_sideband->Fill(value);
}

float StKFParticleAnalysisMaker::GetPtWeight(KFParticle Omega)
{
	double pT = Omega.GetPt();
	int pT_bin = static_cast<int>(floor(pT / 0.05)); 
	float weights[200] = {0, 0, 0, 0.579375, 0, 0, 0, 0.579375, 0, 0, 4.05562, 2.02781, 2.75203, 10.7184, 2.54925, 1.62225, 1.52569, 1.53602, 1.19283, 1.29597, 1.29642, 1.07884, 1.12342, 1.31173, 0.916963, 1.13638, 1.07355, 0.914054, 0.998818, 1.05354, 1.00672, 0.984809, 0.997879, 1.02228, 1.0681, 0.903964, 0.948403, 1.06452, 0.99475, 0.947232, 0.962573, 0.957989, 0.898983, 0.970493, 0.947635, 0.915149, 0.992132, 0.877035, 0.917173, 0.929861, 0.982079, 1.04826, 0.947232, 0.928704, 0.954652, 1.16292, 0.803057, 0.971567, 0.982665, 0.746402, 0.941484, 1.16658, 1.15035, 1.05253, 0.869062, 0.887553, 1.07598, 1.19496, 1.08399, 0.971854, 1.30989, 1.01391, 0.816392, 1.15875, 1.15875, 2.25312, 1.04287, 1.35187, 0.474034, 1.73812, 0.482812, 1.01391, 0.7725, 1.15875, 2.3175, 0.579375, 1.35187, 1.15875, 1.15875, 1.44844, 0, 0, 1.15875, 1.15875, 0, 0, 0, 0, 0, 0, 0, 0.579375, 0, 0.579375, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
	return weights[pT_bin];
}
