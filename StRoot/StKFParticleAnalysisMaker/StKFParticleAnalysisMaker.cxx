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
#include "TEfficiency.h"
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
#define XiPdg              3312
#define phiPdg			   333
#define OmegaPdg           3334
#define KaonPdg			   321
#define ProtonPdg          2212
#define KsPdg			   310
#define Xi1530Pdg		   3324
#define PionPdg           -211
#define cen_cut            9
#define num_pt_bin         10
#define num_phi_bin        10
#define num_mult_bin       5
#define num_vz_bin         16
#define num_EP_bin         6
#define shift_order_asso   24
#define shift_order_POI    24
#define shift_order_EP     5
#define TPCAssoEtaCut      0.5

const float StKFParticleAnalysisMaker::OmegaMassSigma[]    = {0.0028, 0.0035, 0.0037, 0.0031, 0.0035, 0.0041, 0.0039}; // {0.00219, 0.00228, 0.00218, 0.00262, 0.00247, 0.00243, 0.00280};
const float StKFParticleAnalysisMaker::OmegabarMassSigma[] = {0.0035, 0.0038, 0.0032, 0.0032, 0.0041, 0.0052, 0.0055}; // {0.00238, 0.00214, 0.00233, 0.00250, 0.00257, 0.00270, 0.00269};
const float StKFParticleAnalysisMaker::OmegaMassPtLowerBin[] = {0.5, 0.9, 1.2, 1.5, 1.8, 2.1, 2.4, 2.7, 3.0, 3.4, 4.0};
const float StKFParticleAnalysisMaker::LambdaMassSigma[]    = {0.00205, 0.00206, 0.00206, 0.00206, 0.00216, 0.00217, 0.00218, 0.00219, 0.00220};
const float StKFParticleAnalysisMaker::LambdabarMassSigma[] = {0.00202, 0.00203, 0.00203, 0.00204, 0.00204, 0.00205, 0.00207, 0.00207, 0.00208};

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
			wt.SetBinContent(ix,iy, (fabs(eta)>3.8)?lin[iy-1]*eta+cub[iy-1]*pow(eta,3):0);
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
		std::cout << "\n**********************************************" << std::endl;
		std::cout << "No TPC shift input! No shift correction used. " << std::endl;
		std::cout << "**********************************************" << std::endl;
		hTPCAssoShiftInput_cos = 0;
		hTPCAssoShiftInput_sin = 0;
		hTPCPOIShiftInput_cos = 0;
		hTPCPOIShiftInput_sin = 0;
		for (int ewFull=0; ewFull<3; ewFull++) 
		{
			hTPCEPShiftInput_cos[ewFull] = 0;
			hTPCEPShiftInput_sin[ewFull] = 0;
		}
	}
	else
	{	
		std::cout << "\n**********************************************" << std::endl;
		std::cout << "TPC shift input found. " << std::endl;
		std::cout << "**********************************************" << std::endl;
		hTPCAssoShiftInput_cos = (TProfile3D*)fTPCShift->Get(Form("TPCAssoShift_cos"));
		hTPCAssoShiftInput_sin = (TProfile3D*)fTPCShift->Get(Form("TPCAssoShift_sin"));
		hTPCPOIShiftInput_cos = (TProfile3D*)fTPCShift->Get(Form("TPCPOIShift_cos"));
		hTPCPOIShiftInput_sin = (TProfile3D*)fTPCShift->Get(Form("TPCPOIShift_sin"));
		
		for (int ewFull=0; ewFull<3; ewFull++) 
		{
			hTPCEPShiftInput_cos[ewFull] = (TProfile3D*)fTPCShift->Get(Form("TPCEPShiftEW%d_cos", ewFull));
			hTPCEPShiftInput_sin[ewFull] = (TProfile3D*)fTPCShift->Get(Form("TPCEPShiftEW%d_sin", ewFull));
		}
	}

	/* EPD weight files */
	fEPDShift = new TFile(Form("EPDShiftInput.root"), "READ");
	if (fEPDShift->IsZombie())
	{
		std::cout << "\n**********************************************" << std::endl;
		std::cout << "No EPD shift input! No user shift correction used. " << std::endl;
		std::cout << "**********************************************" << std::endl;
		for (int ewFull=0; ewFull<3; ewFull++) 
		{
			hEPDEPShiftInput_1_cos[ewFull] = 0;
			hEPDEPShiftInput_1_sin[ewFull] = 0;
			hEPDEPShiftInput_2_cos[ewFull] = 0;
			hEPDEPShiftInput_2_sin[ewFull] = 0;
		}
		}
		else
		{
		std::cout << "\n**********************************************" << std::endl;
		std::cout << "User EPD shift input found. " << std::endl;
		std::cout << "**********************************************" << std::endl;
		for (int ewFull=0; ewFull<3; ewFull++) 
		{
			hEPDEPShiftInput_1_cos[ewFull] = (TProfile3D*)fEPDShift->Get(Form("EPDEPShiftEW%d_1_cos", ewFull));
			hEPDEPShiftInput_1_sin[ewFull] = (TProfile3D*)fEPDShift->Get(Form("EPDEPShiftEW%d_1_sin", ewFull));
			hEPDEPShiftInput_2_cos[ewFull] = (TProfile3D*)fEPDShift->Get(Form("EPDEPShiftEW%d_2_cos", ewFull));
			hEPDEPShiftInput_2_sin[ewFull] = (TProfile3D*)fEPDShift->Get(Form("EPDEPShiftEW%d_2_sin", ewFull));
		}
		
	}

	std::cout << "Before reading TOF" << std::endl;
	/* TOF Efficiency */
	fTOFEff = new TFile("TOFEfficiency.root", "READ");
	std::cout << "After reading TOF" << std::endl;
	if (fTOFEff->IsZombie())
	{
		std::cout << "\n**********************************************" << std::endl;
		std::cout << "No TOF efficiency input! No TOF efficiency correction used. " << std::endl;
		std::cout << "**********************************************" << std::endl;
		for (int i = 0; i< 9; i++) hTOFEff[i] = 0;
		for (int i = 0; i< 9; i++) hTOFEff_2D[i] = 0;
	}
	else
	{
		std::cout << "\n**********************************************" << std::endl;
		std::cout << "TOF efficiency found. " << std::endl;
		std::cout << "**********************************************" << std::endl;
		for (int i = 0; i < 9; i++) hTOFEff[i] = (TEfficiency*)fTOFEff->Get(Form("hTOFEff_%d", i+1));
		for (int i = 0; i < 9; i++) hTOFEff_2D[i] = (TEfficiency*)fTOFEff->Get(Form("hTOFEff_2D_%d", i+1));
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
	mRunStart = 50000;
	mRunEnd   = 100000;

	PI = M_PI;
	twoPI = 2*M_PI;
	buffer.Init();

	badList.clear();
	runList.clear();

	PerformAnalysis = true;
	PerformMixing = false;
	StoringTree = false;
	CutCent = false;
	PtReweighting = false;
	v2Calculation = true;

	if(!readRunList())return kStFatal;
	//if(!readBadList())return kStFatal;

	dcatoPV_hi = 3.0;
	// pid boundary
	pT_lo = 0.2;
	pT_hi = 2.0;
	
	// for EP
	pT_asso_lo = 0.15;
	pT_asso_hi = 2.0;
	pT_trig_lo = 0.2;
	pT_trig_hi = 2.0;
	eta_trig_cut = 1.0;
	
	// pion cut
	pion_pT_lo = 0.2;
	pion_pT_hi = 1.2;
	pion_pT_TOFth = 0.4;
	pion_m2_lo = -0.1;
	pion_m2_hi = 0.1;

	// proton cut
	proton_pT_lo = 0.3;
	proton_pT_hi = 1.8;
	proton_pT_TOFth = 0.6;
	proton_m2_lo = 0.75;
	proton_m2_hi = 1.1;

	// correlation centrality
	min_cent = 1;
	max_cent = 3;

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
	TFile *fEPDShift_out = new TFile(Form("EPDShiftOutput_%d.root", mJob), "RECREATE");
	fEPDShift_out->cd();
	WriteEPDShift();
	fEPDShift_out->Close();

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
	hEventQA     = new TH1F("hEventQA", "hEventQA", 7, -0.5, 6.5);
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
	for (int i = 0; i < 9; i++)
	{
		hgpT[i]         = new TH1D(Form("hgpT_%d", i+1), "Global transverse momentum", 500, 0., 5.);
		hgpT_TOF[i]     = new TH1D(Form("hgpT_TOF_%d", i+1), "Global transverse momentum with TOF match", 500, 0., 5.);
		hgpTeta[i] 	    = new TH2F(Form("hgpTeta_%d", i+1), "Global transverse momentum vs eta", 500, 0., 5., 400, -2., 2.);
		hgpTeta_TOF[i] 	= new TH2F(Form("hgpTeta_TOF_%d", i+1), "Global transverse momentum vs eta with TOF match", 500, 0., 5., 400, -2., 2.);
	}
	hgDCAtoPV    = new TH1D("hgDCAtoPV", "Global DCA to PV", 500, 0., 10.);
	hgbtofYlocal = new TH1D("hgbtofYlocal", "Global K+ BTOF Ylocal", 1000, -5., 5.);
	hDauProtonFirstPoint_lam_pt = new TProfile("hDauProtonFirstPoint_lam_pt", "hDauProtonFirstPoint_lam_pt", 100, 0., 10., 0., 500.);
	hDauProtonLastPoint_lam_pt  = new TProfile("hDauProtonLastPoint_lam_pt", "hDauProtonLastPoint_lam_pt", 100, 0., 10., 0., 500.);
	hDauPionFirstPoint_lam_pt   = new TProfile("hDauPionFirstPoint_lam_pt", "hDauPionFirstPoint_lam_pt", 100, 0., 10., 0., 500.);
	hDauPionLastPoint_lam_pt    = new TProfile("hDauPionLastPoint_lam_pt", "hDauPionLastPoint_lam_pt", 100, 0., 10., 0., 500.);

	// v2 and EP 
	hTPCAssoPhi         = new TH1D("hTPCAssoPhi",         "hTPCAssoPhi",         1000, -PI, PI);
	hTPCAssoPhi_shifted = new TH1D("hTPCAssoPhi_shifted", "hTPCAssoPhi_shifted", 1000, -PI, PI);
	hTPCPOIPhi         = new TH1D("hTPCPOIPhi",         "hTPCPOIPhi",         1000, -PI, PI);
	hTPCPOIPhi_shifted = new TH1D("hTPCPOIPhi_shifted", "hTPCPOIPhi_shifted", 1000, -PI, PI);
	hTPCAssoPhi_2D = new TH2D("hTPCAssoPhi_2D", "hTPCAssoPhi_2D", 50, -0.5, 49.5, 1000, -PI, PI);
	hTPCAssoPhi_2D_shifted = new TH2D("hTPCAssoPhi_2D_shifted", "hTPCAssoPhi_2D_shifted", 50, -0.5, 49.5, 1000, -PI, PI);
	hTPCPOIPhi_2D = new TH2D("hTPCPOIPhi_2D", "hTPCPOIPhi_2D", 50, -0.5, 49.5, 1000, -PI, PI);
	hTPCPOIPhi_2D_shifted = new TH2D("hTPCPOIPhi_2D_shifted", "hTPCPOIPhi_2D_shifted", 50, -0.5, 49.5, 1000, -PI, PI);
	// TPC weights/shifts
	hTPCAssoShiftOutput_cos = new TProfile3D(Form("TPCAssoShift_cos"), Form("TPCAssoShift_cos"), 
											 50, -0.5, 49.5, shift_order_asso, 0.5, 1.0 * shift_order_asso + 0.5, 9, -0.5, 8.5);
	hTPCAssoShiftOutput_sin = new TProfile3D(Form("TPCAssoShift_sin"), Form("TPCAssoShift_sin"), 
											 50, -0.5, 49.5, shift_order_asso, 0.5, 1.0 * shift_order_asso + 0.5, 9, -0.5, 8.5);
	hTPCPOIShiftOutput_cos = new TProfile3D(Form("TPCPOIShift_cos"), Form("TPCPOIShift_cos"), 
											 50, -0.5, 49.5, shift_order_POI, 0.5, 1.0 * shift_order_POI + 0.5, 9, -0.5, 8.5);
	hTPCPOIShiftOutput_sin = new TProfile3D(Form("TPCPOIShift_sin"), Form("TPCPOIShift_sin"), 
											 50, -0.5, 49.5, shift_order_POI, 0.5, 1.0 * shift_order_POI + 0.5, 9, -0.5, 8.5);
	for (int ewFull = 0; ewFull < 3; ewFull++)
	{
		hTPCEPShiftOutput_cos[ewFull] = new TProfile3D(Form("TPCEPShiftEW%d_cos", ewFull), Form("TPCEPShiftEW%d_cos", ewFull), 
												50, -0.5, 49.5, shift_order_EP, 0.5, 1.0 * shift_order_EP + 0.5, 9, -0.5, 8.5);
		hTPCEPShiftOutput_sin[ewFull] = new TProfile3D(Form("TPCEPShiftEW%d_sin", ewFull), Form("TPCEPShiftEW%d_sin", ewFull), 
												50, -0.5, 49.5, shift_order_EP, 0.5, 1.0 * shift_order_EP + 0.5, 9, -0.5, 8.5);
		for (int i = 0; i < 9; i++)
		{
			hTPCEP_2[i][ewFull]        = new TH1D(Form("hTPCEPEW%d_2_%d", ewFull, i+1),         Form("hTPCEPEW%d_2_%d",         ewFull, i+1), 1000, 0., 2*PI);
			hTPCEP_2_shifted[i][ewFull]= new TH1D(Form("hTPCEPEW%d_2_shifted_%d", ewFull, i+1), Form("hTPCEPEW%d_2_shifted_%d", ewFull, i+1), 1000, 0., 2*PI);
			hTPCEP_2_2D[i][ewFull]         = new TH2D(Form("hTPCEPEW%d_2_2D_%d", ewFull, i+1),         Form("hTPCEPEW%d_2_2D_%d",         ewFull, i+1), 50, -0.5, 49.5, 1000, 0., 2*PI);
			hTPCEP_2_2D_shifted[i][ewFull] = new TH2D(Form("hTPCEPEW%d_2_2D_shifted_%d", ewFull, i+1), Form("hTPCEPEW%d_2_2D_shifted_%d", ewFull, i+1), 50, -0.5, 49.5, 1000, 0., 2*PI);
			hTPCEP_2_shift[i][ewFull] = new TProfile(Form("hTPCEPEW%d_2_shift_%d", ewFull, i+1), Form("hTPCEPEW%d_2_shift_%d", ewFull, i+1), 1000, 0., 2*PI, -PI, PI);
	}
	}
	// EPD user shifts
	for (int ewFull = 0; ewFull < 3; ewFull++)
	{
		hEPDEPShiftOutput_1_cos[ewFull] = new TProfile3D(Form("EPDEPShiftEW%d_1_cos", ewFull), Form("EPDEPShiftEW%d_1_cos", ewFull), 
												50, -0.5, 49.5, shift_order_EP, 0.5, 1.0 * shift_order_EP + 0.5, 9, -0.5, 8.5);
		hEPDEPShiftOutput_1_sin[ewFull] = new TProfile3D(Form("EPDEPShiftEW%d_1_sin", ewFull), Form("EPDEPShiftEW%d_1_sin", ewFull),
												50, -0.5, 49.5, shift_order_EP, 0.5, 1.0 * shift_order_EP + 0.5, 9, -0.5, 8.5);
		hEPDEPShiftOutput_2_cos[ewFull] = new TProfile3D(Form("EPDEPShiftEW%d_2_cos", ewFull), Form("EPDEPShiftEW%d_2_cos", ewFull),
												50, -0.5, 49.5, shift_order_EP, 0.5, 1.0 * shift_order_EP + 0.5, 9, -0.5, 8.5);
		hEPDEPShiftOutput_2_sin[ewFull] = new TProfile3D(Form("EPDEPShiftEW%d_2_sin", ewFull), Form("EPDEPShiftEW%d_2_sin", ewFull),
												50, -0.5, 49.5, shift_order_EP, 0.5, 1.0 * shift_order_EP + 0.5, 9, -0.5, 8.5);
		for (int i = 0; i < 9; i++)
		{
			hEPDEP_1[i][ewFull]         = new TH1D(Form("hEPDEPW%d_1_%d",         ewFull, i+1), Form("hEPDEPW%d_1_%d",         ewFull, i+1), 1000, 0., 2*PI);
			hEPDEP_1_shifted[i][ewFull] = new TH1D(Form("hEPDEPW%d_1_shifted_%d", ewFull, i+1), Form("hEPDEPW%d_1_shifted_%d", ewFull, i+1), 1000, 0., 2*PI);
			hEPDEP_2[i][ewFull]         = new TH1D(Form("hEPDEPW%d_2_%d",         ewFull, i+1), Form("hEPDEPW%d_2_%d",         ewFull, i+1), 1000, 0., 2*PI);
			hEPDEP_2_shifted[i][ewFull] = new TH1D(Form("hEPDEPW%d_2_shifted_%d", ewFull, i+1), Form("hEPDEPW%d_2_shifted_%d", ewFull, i+1), 1000, 0., 2*PI);
			hEPDEP_1_2D[i][ewFull]         = new TH2D(Form("hEPDEPW%d_1_2D_%d",         ewFull, i+1), Form("hEPDEPW%d_1_2D_%d",         ewFull, i+1), 50, -0.5, 49.5, 1000, 0., 2*PI);
			hEPDEP_1_2D_shifted[i][ewFull] = new TH2D(Form("hEPDEPW%d_1_2D_shifted_%d", ewFull, i+1), Form("hEPDEPW%d_1_2D_shifted_%d", ewFull, i+1), 50, -0.5, 49.5, 1000, 0., 2*PI);
			hEPDEP_2_2D[i][ewFull]         = new TH2D(Form("hEPDEPW%d_2_2D_%d",         ewFull, i+1), Form("hEPDEPW%d_2_2D_%d",         ewFull, i+1), 50, -0.5, 49.5, 1000, 0., 2*PI);
			hEPDEP_2_2D_shifted[i][ewFull] = new TH2D(Form("hEPDEPW%d_2_2D_shifted_%d", ewFull, i+1), Form("hEPDEPW%d_2_2D_shifted_%d", ewFull, i+1), 50, -0.5, 49.5, 1000, 0., 2*PI);
			hEPDEP_2_shift[i][ewFull] = new TProfile(Form("hEPDEPW%d_2_shift_%d", ewFull, i+1), Form("hEPDEPW%d_2_shift_%d", ewFull, i+1), 1000, 0., 2*PI, -PI, PI);
		}
	}										
	
	hTPCEP_ew_cos = new TProfile("hTPCEP_ew_cos", "hTPCEP_ew_cos", 9, 0.5, 9.5, -1., 1.);
	hTPCEP_ew_cos->Sumw2();
	hEPDEP_ew_cos_1 = new TProfile("hEPDEP_ew_cos_1", "hEPDEP_ew_cos_1", 9, 0.5, 9.5, -1., 1.);
	hEPDEP_ew_cos_1->Sumw2();
	hEPDEP_ew_cos_2 = new TProfile("hEPDEP_ew_cos_2", "hEPDEP_ew_cos_2", 9, 0.5, 9.5, -1., 1.);
	hEPDEP_ew_cos_2->Sumw2();
	// hEPDEP_ew_cos_1_for_v2 = new TProfile("hEPDEP_ew_cos_1_for_v2", "hEPDEP_ew_cos_1_for_v2", 9, 0.5, 9.5, -1., 1.);
	// hEPDEP_ew_cos_1_for_v2->Sumw2();

	hOmega_TPC_v2_pt = new TProfile("hOmega_TPC_v2_pt", "hOmega_TPC_v2_pt", 50, 0., 5., -1., 1.);
	hOmegabar_TPC_v2_pt = new TProfile("hOmegabar_TPC_v2_pt", "hOmegabar_TPC_v2_pt", 50, 0., 5., -1., 1.);
	hOmega_EPD_v2_pt = new TProfile("hOmega_EPD_v2_pt", "hOmega_EPD_v2_pt", 50, 0., 5., -1., 1.);
	hOmegabar_EPD_v2_pt = new TProfile("hOmegabar_EPD_v2_pt", "hOmegabar_EPD_v2_pt", 50, 0., 5., -1., 1.);

	hpiplus_EPD_v2     = new TProfile("hpiplus_EPD_v2",  "hpiplus_EPD_v2",  9, 0.5, 9.5, -1., 1.);
	hpiplus_EPD_v2->Sumw2();
	hpiminus_EPD_v2    = new TProfile("hpiminus_EPD_v2", "hpiminus_EPD_v2", 9, 0.5, 9.5, -1., 1.);
	hpiminus_EPD_v2->Sumw2();
	hproton_EPD_v2     = new TProfile("hproton_EPD_v2",  "hproton_EPD_v2",  9, 0.5, 9.5, -1., 1.);
	hproton_EPD_v2->Sumw2();
	hantiproton_EPD_v2 = new TProfile("hantiproton_EPD_v2",  "hantiproton_EPD_v2",  9, 0.5, 9.5, -1., 1.);
	hantiproton_EPD_v2->Sumw2();
	hkplus_EPD_v2 = new TProfile("hkplus_EPD_v2",  "hkplus_EPD_v2",  9, 0.5, 9.5, -1., 1.);
	hkplus_EPD_v2->Sumw2();
	hkminus_EPD_v2 = new TProfile("hkminus_EPD_v2",  "hkminus_EPD_v2",  9, 0.5, 9.5, -1., 1.);
	hkminus_EPD_v2->Sumw2();

	hpiplus_TPC_v2     = new TProfile("hpiplus_TPC_v2",  "hpiplus_TPC_v2",  9, 0.5, 9.5, -1., 1.);
	hpiplus_TPC_v2->Sumw2();
	hpiminus_TPC_v2    = new TProfile("hpiminus_TPC_v2", "hpiminus_TPC_v2", 9, 0.5, 9.5, -1., 1.);
	hpiminus_TPC_v2->Sumw2();
	hproton_TPC_v2     = new TProfile("hproton_TPC_v2",  "hproton_TPC_v2",  9, 0.5, 9.5, -1., 1.);
	hproton_TPC_v2->Sumw2();
	hantiproton_TPC_v2 = new TProfile("hantiproton_TPC_v2",  "hantiproton_TPC_v2",  9, 0.5, 9.5, -1., 1.);
	hantiproton_TPC_v2->Sumw2();
	hkplus_TPC_v2 = new TProfile("hkplus_TPC_v2",  "hkplus_TPC_v2",  9, 0.5, 9.5, -1., 1.);
	hkplus_TPC_v2->Sumw2();
	hkminus_TPC_v2 = new TProfile("hkminus_TPC_v2",  "hkminus_TPC_v2",  9, 0.5, 9.5, -1., 1.);
	hkminus_TPC_v2->Sumw2();
	for (int i = 0; i < 9; i++)
	{
		hOmega_EPD_v2[i]    = new TProfile(Form("hOmega_EPD_v2_%d",    i+1), Form("hOmega_EPD_v2_%d",    i+1), 1400, 1, 2.4, -1., 1.);
		hOmega_EPD_v2[i]->Sumw2();
		hOmegabar_EPD_v2[i] = new TProfile(Form("hOmegabar_EPD_v2_%d", i+1), Form("hOmegabar_EPD_v2_%d", i+1), 1400, 1, 2.4, -1., 1.);
		hOmegabar_EPD_v2[i]->Sumw2();
		hOmega_TPC_v2[i]    = new TProfile(Form("hOmega_TPC_v2_%d",    i+1), Form("hOmega_TPC_v2_%d",    i+1), 1400, 1, 2.4, -1., 1.);
		hOmega_TPC_v2[i]->Sumw2();
		hOmegabar_TPC_v2[i] = new TProfile(Form("hOmegabar_TPC_v2_%d", i+1), Form("hOmegabar_TPC_v2_%d", i+1), 1400, 1, 2.4, -1., 1.);
		hOmegabar_TPC_v2[i]->Sumw2();

		hpiplus_EPD_v2_pt[i] = new TProfile(Form("hpiplus_EPD_v2_pt_%d", i+1), Form("hpiplus_EPD_v2_pt_%d", i+1), 1000, 0., 10., -1., 1.);
		hpiplus_EPD_v2_pt[i]->Sumw2();
		hpiminus_EPD_v2_pt[i] = new TProfile(Form("hpiminus_EPD_v2_pt_%d", i+1), Form("hpiminus_EPD_v2_pt_%d", i+1), 1000, 0., 10., -1., 1.);
		hpiminus_EPD_v2_pt[i]->Sumw2();
		hproton_EPD_v2_pt[i] = new TProfile(Form("hproton_EPD_v2_pt_%d", i+1), Form("hproton_EPD_v2_pt_%d", i+1), 1000, 0., 10., -1., 1.);
		hproton_EPD_v2_pt[i]->Sumw2();
		hantiproton_EPD_v2_pt[i] = new TProfile(Form("hantiproton_EPD_v2_pt_%d", i+1), Form("hantiproton_EPD_v2_pt_%d", i+1), 1000, 0., 10., -1., 1.);
		hantiproton_EPD_v2_pt[i]->Sumw2();
		hkplus_EPD_v2_pt[i] = new TProfile(Form("hkplus_EPD_v2_pt_%d", i+1), Form("hkplus_EPD_v2_pt_%d", i+1), 1000, 0., 10., -1., 1.);
		hkplus_EPD_v2_pt[i]->Sumw2();
		hkminus_EPD_v2_pt[i] = new TProfile(Form("hkminus_EPD_v2_pt_%d", i+1), Form("hkminus_EPD_v2_pt_%d", i+1), 1000, 0., 10., -1., 1.);
		hkminus_EPD_v2_pt[i]->Sumw2();
		hpiplus_TPC_v2_pt[i] = new TProfile(Form("hpiplus_TPC_v2_pt_%d", i+1), Form("hpiplus_TPC_v2_pt_%d", i+1), 1000, 0., 10., -1., 1.);
		hpiplus_TPC_v2_pt[i]->Sumw2();
		hpiminus_TPC_v2_pt[i] = new TProfile(Form("hpiminus_TPC_v2_pt_%d", i+1), Form("hpiminus_TPC_v2_pt_%d", i+1), 1000, 0., 10., -1., 1.);
		hpiminus_TPC_v2_pt[i]->Sumw2();
		hproton_TPC_v2_pt[i] = new TProfile(Form("hproton_TPC_v2_pt_%d", i+1), Form("hproton_TPC_v2_pt_%d", i+1), 1000, 0., 10., -1., 1.);
		hproton_TPC_v2_pt[i]->Sumw2();
		hantiproton_TPC_v2_pt[i] = new TProfile(Form("hantiproton_TPC_v2_pt_%d", i+1), Form("hantiproton_TPC_v2_pt_%d", i+1), 1000, 0., 10., -1., 1.);
		hantiproton_TPC_v2_pt[i]->Sumw2();
		hkplus_TPC_v2_pt[i] = new TProfile(Form("hkplus_TPC_v2_pt_%d", i+1), Form("hkplus_TPC_v2_pt_%d", i+1), 1000, 0., 10., -1., 1.);
		hkplus_TPC_v2_pt[i]->Sumw2();
		hkminus_TPC_v2_pt[i] = new TProfile(Form("hkminus_TPC_v2_pt_%d", i+1), Form("hkminus_TPC_v2_pt_%d", i+1), 1000, 0., 10., -1., 1.);
		hkminus_TPC_v2_pt[i]->Sumw2();

		// v1
		hpiplus_EPD_v1_y[i] = new TProfile(Form("hpiplus_EPD_v1_y_%d", i+1), Form("hpiplus_EPD_v1_y_%d", i+1), 10, -1., 1., -1., 1.);
		hpiplus_EPD_v1_y[i]->Sumw2();
		hpiminus_EPD_v1_y[i] = new TProfile(Form("hpiminus_EPD_v1_y_%d", i+1), Form("hpiminus_EPD_v1_y_%d", i+1), 10, -1., 1., -1., 1.);
		hpiminus_EPD_v1_y[i]->Sumw2();
		hproton_EPD_v1_y[i] = new TProfile(Form("hproton_EPD_v1_y_%d", i+1), Form("hproton_EPD_v1_y_%d", i+1), 10, -1., 1., -1., 1.);
		hproton_EPD_v1_y[i]->Sumw2();
		hantiproton_EPD_v1_y[i] = new TProfile(Form("hantiproton_EPD_v1_y_%d", i+1), Form("hantiproton_EPD_v1_y_%d", i+1), 10, -1., 1., -1., 1.);
		hantiproton_EPD_v1_y[i]->Sumw2();
		hkplus_EPD_v1_y[i] = new TProfile(Form("hkplus_EPD_v1_y_%d", i+1), Form("hkplus_EPD_v1_y_%d", i+1), 10, -1., 1., -1., 1.);
		hkplus_EPD_v1_y[i]->Sumw2();
		hkminus_EPD_v1_y[i] = new TProfile(Form("hkminus_EPD_v1_y_%d", i+1), Form("hkminus_EPD_v1_y_%d", i+1), 10, -1., 1., -1., 1.);
		hkminus_EPD_v1_y[i]->Sumw2();
	}

	hEPD_e_EP_1 = new TH1D("hEPD_e_EP_1", "hEPD_e_EP_1", 1000, 0., 2*PI);
	hEPD_w_EP_1 = new TH1D("hEPD_w_EP_1", "hEPD_w_EP_1", 1000, 0., 2*PI);
    hEPD_e_EP_2 = new TH1D("hEPD_e_EP_2", "hEPD_e_EP_2", 1000, 0., 2*PI);
	hEPD_w_EP_2 = new TH1D("hEPD_w_EP_2", "hEPD_w_EP_2", 1000, 0., 2*PI);
	hEPD_full_EP_1 = new TH1D("hEPD_full_EP_1", "hEPD_full_EP_1", 1000, 0., 2*PI);
	hEPD_full_EP_2 = new TH1D("hEPD_full_EP_2", "hEPD_full_EP_2", 1000, 0., 2*PI);
	hEPD_ew_cos = new TProfile("hEPD_ew_cos", "hEPD_ew_cos", 3, 0.5, 3.5, -1., 1.); hEPD_ew_cos->Sumw2();

	// 2D pid
	for (int i = 0; i < 10; i++)
	{	
		hgPID2D_proton_pt[i] = new TH2F(Form("hgPID2D_proton_pt_%d", i), Form("hgPID2D_proton_pt_%d", i), 2000, -10, 10, 400, -2, 2);
		hgPID2D_antiproton_pt[i] = new TH2F(Form("hgPID2D_antiproton_pt_%d", i), Form("hgPID2D_antiproton_pt_%d", i), 2000, -10, 10, 400, -2, 2);	
		hgPID2D_piplus_pt[i]   = new TH2F(Form("hgPID2D_piplus_pt_%d", i), Form("hgPID2D_piplus_pt_%d", i), 2000, -10, 10, 400, -2, 2);
		hgPID2D_piminus_pt[i]  = new TH2F(Form("hgPID2D_piminus_pt_%d", i), Form("hgPID2D_piminus_pt_%d", i), 2000, -10, 10, 400, -2, 2);
		hgPID2D_kplus_pt[i]   = new TH2F(Form("hgPID2D_kplus_pt_%d", i), Form("hgPID2D_kplus_pt_%d", i), 2000, -10, 10, 400, -2, 2);
		hgPID2D_kminus_pt[i]  = new TH2F(Form("hgPID2D_kminus_pt_%d", i), Form("hgPID2D_kminus_pt_%d", i), 2000, -10, 10, 400, -2, 2);
	}

	// testing eTOF
	hgetofcrossingY = new TH1D("hgetofcrossingY", "hgetofcrossingY", 1000, -50., 50.);
	hgm2_etof = new TH1D("hgm2_etof", "hgm2_etof", 500, -1., 4.);
	hgpinvbeta_etof = new TH2F("hgpinvbeta_etof", "hgpinvbeta_etof", 1000, 0., 10., 1000, 0., 10.);
	hgetofmatchflag = new TH1D("hgetofmatchflag", "hgetofmatchflag", 3, -0.5, 2.5);
	hgeta_etof = new TH1D("hgeta_etof", "hgeta_etof", 1000, -5., 5.);

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
	hKplusy_omega = new TH1D("hKplusy_omega", "K+ rapidity", 1000, -5., 5.);
	hKplusy_omegabar = new TH1D("hKplusy_omegabar", "K+ rapidity", 1000, -5., 5.);
	hKminusy_omega = new TH1D("hKminusy_omega", "K- rapidity", 1000, -5., 5.);
	hKminusy_omegabar = new TH1D("hKminusy_omegabar", "K- rapidity", 1000, -5., 5.);
	// sideband
	hKpluspt_omega_sideband      = new TH1D("hKpluspt_omega_sideband",      "Global K transver momentum", 1000, 0., 10.);
	hKpluspt_omegabar_sideband   = new TH1D("hKpluspt_omegabar_sideband",   "Global K transver momentum", 1000, 0., 10.);
	hKminuspt_omega_sideband     = new TH1D("hKminuspt_omega_sideband",     "Global K transver momentum", 1000, 0., 10.);
	hKminuspt_omegabar_sideband  = new TH1D("hKminuspt_omegabar_sideband",  "Global K transver momentum", 1000, 0., 10.);
	hKpluseta_omega_sideband     = new TH1D("hKpluseta_omega_sideband",     "Global K transver momentum", 1000, -5., 5.);
	hKpluseta_omegabar_sideband  = new TH1D("hKpluseta_omegabar_sideband",  "Global K transver momentum", 1000, -5., 5.);
	hKminuseta_omega_sideband    = new TH1D("hKminuseta_omega_sideband",    "Global K transver momentum", 1000, -5., 5.);
	hKminuseta_omegabar_sideband = new TH1D("hKminuseta_omegabar_sideband", "Global K transver momentum", 1000, -5., 5.);
	hKplusphi_omega_sideband     = new TH1D("hKplusphi_omega_sideband",     "Global K transver momentum", 1000, -pi, pi);
	hKplusphi_omegabar_sideband  = new TH1D("hKplusphi_omegabar_sideband",  "Global K transver momentum", 1000, -pi, pi);
	hKminusphi_omega_sideband    = new TH1D("hKminusphi_omega_sideband",    "Global K transver momentum", 1000, -pi, pi);
	hKminusphi_omegabar_sideband = new TH1D("hKminusphi_omegabar_sideband", "Global K transver momentum", 1000, -pi, pi);
	hKplusy_omega_sideband       = new TH1D("hKplusy_omega_sideband",       "K+ rapidity", 1000, -5., 5.);
	hKplusy_omegabar_sideband    = new TH1D("hKplusy_omegabar_sideband",    "K+ rapidity", 1000, -5., 5.);
	hKminusy_omega_sideband      = new TH1D("hKminusy_omega_sideband",      "K- rapidity", 1000, -5., 5.);
	hKminusy_omegabar_sideband   = new TH1D("hKminusy_omegabar_sideband",   "K- rapidity", 1000, -5., 5.);

	// w/ or wo/ lambda
	hKpluspt_wlb_omega    = new TH1D("hKpluspt_wlb_omega", "Global K transver momentum", 1000, 0., 10.);
	hKpluspt_wlb_omegabar = new TH1D("hKpluspt_wlb_omegabar", "Global K transver momentum", 1000, 0., 10.);
	hKminuspt_wl_omega    = new TH1D("hKminuspt_wl_omega", "Global K transver momentum", 1000, 0., 10.);
	hKminuspt_wl_omegabar = new TH1D("hKminuspt_wl_omegabar", "Global K transver momentum", 1000, 0., 10.);
	hKplusy_wlb_omega    = new TH1D("hKplusy_wlb_omega", "K+ rapidity", 1000, -5., 5.);
	hKplusy_wlb_omegabar = new TH1D("hKplusy_wlb_omegabar", "K+ rapidity", 1000, -5., 5.);
	hKminusy_wl_omega    = new TH1D("hKminusy_wl_omega", "K- rapidity", 1000, -5., 5.);
	hKminusy_wl_omegabar = new TH1D("hKminusy_wl_omegabar", "K- rapidity", 1000, -5., 5.);
	hKpluspt_wolb_omega    = new TH1D("hKpluspt_wolb_omega", "Global K transver momentum", 1000, 0., 10.);
	hKpluspt_wolb_omegabar = new TH1D("hKpluspt_wolb_omegabar", "Global K transver momentum", 1000, 0., 10.);
	hKminuspt_wol_omega    = new TH1D("hKminuspt_wol_omega", "Global K transver momentum", 1000, 0., 10.);
	hKminuspt_wol_omegabar = new TH1D("hKminuspt_wol_omegabar", "Global K transver momentum", 1000, 0., 10.);
	hKplusy_wolb_omega    = new TH1D("hKplusy_wolb_omega", "K+ rapidity", 1000, -5., 5.);
	hKplusy_wolb_omegabar = new TH1D("hKplusy_wolb_omegabar", "K+ rapidity", 1000, -5., 5.);
	hKminusy_wol_omega    = new TH1D("hKminusy_wol_omega", "K- rapidity", 1000, -5., 5.);
	hKminusy_wol_omegabar = new TH1D("hKminusy_wol_omegabar", "K- rapidity", 1000, -5., 5.);
	// and their sidebands
	hKpluspt_wlb_omega_sideband    = new TH1D("hKpluspt_wlb_omega_sideband", "Global K transver momentum", 1000, 0., 10.);
	hKpluspt_wlb_omegabar_sideband = new TH1D("hKpluspt_wlb_omegabar_sideband", "Global K transver momentum", 1000, 0., 10.);
	hKminuspt_wl_omega_sideband    = new TH1D("hKminuspt_wl_omega_sideband", "Global K transver momentum", 1000, 0., 10.);
	hKminuspt_wl_omegabar_sideband = new TH1D("hKminuspt_wl_omegabar_sideband", "Global K transver momentum", 1000, 0., 10.);
	hKplusy_wlb_omega_sideband    = new TH1D("hKplusy_wlb_omega_sideband", "K+ rapidity", 1000, -5., 5.);
	hKplusy_wlb_omegabar_sideband = new TH1D("hKplusy_wlb_omegabar_sideband", "K+ rapidity", 1000, -5., 5.);
	hKminusy_wl_omega_sideband    = new TH1D("hKminusy_wl_omega_sideband", "K- rapidity", 1000, -5., 5.);
	hKminusy_wl_omegabar_sideband = new TH1D("hKminusy_wl_omegabar_sideband", "K- rapidity", 1000, -5., 5.);
	hKpluspt_wolb_omega_sideband    = new TH1D("hKpluspt_wolb_omega_sideband", "Global K transver momentum", 1000, 0., 10.);
	hKpluspt_wolb_omegabar_sideband = new TH1D("hKpluspt_wolb_omegabar_sideband", "Global K transver momentum", 1000, 0., 10.);
	hKminuspt_wol_omega_sideband    = new TH1D("hKminuspt_wol_omega_sideband", "Global K transver momentum", 1000, 0., 10.);
	hKminuspt_wol_omegabar_sideband = new TH1D("hKminuspt_wol_omegabar_sideband", "Global K transver momentum", 1000, 0., 10.);
	hKplusy_wolb_omega_sideband    = new TH1D("hKplusy_wolb_omega_sideband", "K+ rapidity", 1000, -5., 5.);
	hKplusy_wolb_omegabar_sideband = new TH1D("hKplusy_wolb_omegabar_sideband", "K+ rapidity", 1000, -5., 5.);
	hKminusy_wol_omega_sideband    = new TH1D("hKminusy_wol_omega_sideband", "K- rapidity", 1000, -5., 5.);
	hKminusy_wol_omegabar_sideband = new TH1D("hKminusy_wol_omegabar_sideband", "K- rapidity", 1000, -5., 5.);

	// for a rapidity check
	hkplus_totaly = new TProfile("hkplus_totaly", "K+ total rapidity in an event", 2, -0.5, 1.5, 0., 500.);
	hkminus_totaly = new TProfile("hkminus_totaly", "K- total rapidity in an event", 2, -0.5, 1.5, 0., 500.);
	hkplus_ntrack = new TProfile("hkplus_ntrack", "K+ number of tracks in an event", 2, -0.5, 1.5, 0., 500.);
	hkminus_ntrack = new TProfile("hkminus_ntrack", "K- number of tracks in an event", 2, -0.5, 1.5, 0., 500.);

	// corresponding Omega inv mass
	hOmegaM_wl = new TH1D("hOmegaM_wl", "Omega inv mass", 1400, 1., 2.4);
	hOmegaM_wlb = new TH1D("hOmegaM_wlb", "Omega inv mass", 1400, 1., 2.4);
	hOmegaM_wol = new TH1D("hOmegaM_wol", "Omega inv mass", 1400, 1., 2.4);
	hOmegaM_wolb = new TH1D("hOmegaM_wolb", "Omega inv mass", 1400, 1., 2.4);
	hOmegabarM_wl = new TH1D("hOmegabarM_wl", "Omegabar inv mass", 1400, 1., 2.4);
	hOmegabarM_wlb = new TH1D("hOmegabarM_wlb", "Omegabar inv mass", 1400, 1., 2.4);
	hOmegabarM_wol = new TH1D("hOmegabarM_wol", "Omegabar inv mass", 1400, 1., 2.4);
	hOmegabarM_wolb = new TH1D("hOmegabarM_wolb", "Omegabar inv mass", 1400, 1., 2.4);

	// proton/pion QA
	hProtony     = new TH1D("hProtony", "Proton Rapidity", 1000, -5., 5.);
	hAntiProtony = new TH1D("hAntiProtony", "Anti-proton Rapidity", 1000, -5., 5.);
	hm2proton_b = new TH1D("hm2proton_b", "m^2 before", 500, -1., 4.);
	hm2proton_r = new TH1D("hm2proton_r", "m^2 regular", 500, -1., 4.);
	hm2proton_a = new TH1D("hm2proton_a", "m^2 after", 500, -1., 4.);
	hm2pion_b = new TH1D("hm2pion_b", "m^2 before", 500, -1., 4.);
	hm2pion_r = new TH1D("hm2pion_r", "m^2 regular", 500, -1., 4.);
	hm2pion_a = new TH1D("hm2pion_a", "m^2 after", 500, -1., 4.);
	hTOFEff_check = new TProfile("hTOFEff_check", "hTOFEff_check", 1000, 0., 10., 0., 1.);
	for (int i = 0; i < 9; i++) hTPCEff_check[i] = new TProfile(Form("hTPCEff_check_%d", i+1), Form("hTPCEff_check_%d", i+1), 1000, 0., 10., 0., 1.);

	// Xi and Lambda v2
	for (int i = 0; i < 9; i++)
	{
		hXiM_cen[i] = new TH1D(Form("hXiM_cen_%d", i+1), Form("hXiM_cen_%d", i+1), 1200, 1., 1.6);
		hXibarM_cen[i] = new TH1D(Form("hXibarM_cen_%d", i+1), Form("hXibarM_cen_%d", i+1), 1200, 1., 1.6);
		hXi_EPD_v2   [i] = new TProfile(Form("hXi_EPD_v2_%d",    i+1), Form("hXi_EPD_v2_%d",    i+1), 1200, 1., 1.6, -1., 1.);
		hXibar_EPD_v2[i] = new TProfile(Form("hXibar_EPD_v2_%d", i+1), Form("hXibar_EPD_v2_%d", i+1), 1200, 1., 1.6, -1., 1.);
		hXi_TPC_v2   [i] = new TProfile(Form("hXi_TPC_v2_%d",    i+1), Form("hXi_TPC_v2_%d",    i+1), 1200, 1., 1.6, -1., 1.);
		hXibar_TPC_v2[i] = new TProfile(Form("hXibar_TPC_v2_%d", i+1), Form("hXibar_TPC_v2_%d", i+1), 1200, 1., 1.6, -1., 1.);
		hphi_EPD_v2[i] = new TProfile(Form("hphi_EPD_v2_%d", i+1), Form("hphi_EPD_v2_%d", i+1), 1200, 0.7, 1.3, -1., 1.);
		hphi_TPC_v2[i] = new TProfile(Form("hphi_TPC_v2_%d", i+1), Form("hphi_TPC_v2_%d", i+1), 1200, 0.7, 1.3, -1., 1.);

		hphiM_cen[i] = new TH1D(Form("hphiM_cen_%d", i+1), Form("hphiM_cen_%d", i+1), 1200, 0.7, 1.3);
		hXi1530M_cen[i] = new TH1D(Form("hXi1530M_cen_%d", i+1), Form("hXi1530M_cen_%d", i+1), 1200, 1.2, 1.8);
		hXi1530barM_cen[i] = new TH1D(Form("hXi1530barM_cen_%d", i+1), Form("hXi1530barM_cen_%d", i+1), 1200, 1.2, 1.8);
		hKsM_cen[i] = new TH1D(Form("hKsM_cen_%d", i+1), Form("hKsM_cen_%d", i+1), 1200, 0.2, 0.8);
		hOmega2012M_cen[i] = new TH1D(Form("hOmega2012M_cen_%d", i+1), Form("hOmega2012M_cen_%d", i+1), 1200, 1.7, 2.3);
		hOmega2012barM_cen[i] = new TH1D(Form("hOmega2012barM_cen_%d", i+1), Form("hOmega2012barM_cen_%d", i+1), 1200, 1.7, 2.3);
		// for (int ptbin=0; ptbin<40; ptbin++)
		// {
		// 	hLambdaM_cen_pt[i][ptbin] = new TH1D(Form("hLambdaM_cen_%d_pt_%d", i+1, ptbin), Form("hLambdaM_cen_%d_pt_%d", i+1, ptbin), 400, 1., 1.2);
		// 	hLambdabarM_cen_pt[i][ptbin] = new TH1D(Form("hLambdabarM_cen_%d_pt_%d", i+1, ptbin), Form("hLambdabarM_cen_%d_pt_%d", i+1, ptbin), 400, 1., 1.2);
		// 	hLambda_EPD_v2_pt[i][ptbin] = new TProfile(Form("hLambda_EPD_v2_%d_pt_%d", i+1, ptbin), Form("hLambda_EPD_v2_%d_pt_%d", i+1, ptbin), 400, 1., 1.2, -1., 1.);
		// 	hLambdabar_EPD_v2_pt[i][ptbin] = new TProfile(Form("hLambdabar_EPD_v2_%d_pt_%d", i+1, ptbin), Form("hLambdabar_EPD_v2_%d_pt_%d", i+1, ptbin), 400, 1., 1.2, -1., 1.);
		// 	hLambda_TPC_v2_pt[i][ptbin] = new TProfile(Form("hLambda_TPC_v2_%d_pt_%d", i+1, ptbin), Form("hLambda_TPC_v2_%d_pt_%d", i+1, ptbin), 400, 1., 1.2, -1., 1.);
		// 	hLambdabar_TPC_v2_pt[i][ptbin] = new TProfile(Form("hLambdabar_TPC_v2_%d_pt_%d", i+1, ptbin), Form("hLambdabar_TPC_v2_%d_pt_%d", i+1, ptbin), 400, 1., 1.2, -1., 1.);
		// }
	}

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
	hOmegaUsed_wlb = new TH1D("hOmegaUsed_wlb", "Actual Omega Used #Omega/#bar{#Omega}", 2, -0.5, 1.5);
	hOmegaUsed_wlb_sideband = new TH1D("hOmegaUsed_wlb_sideband", "Actual Omega Used #Omega/#bar{#Omega}", 2, -0.5, 1.5);
	hOmegaUsed_wolb = new TH1D("hOmegaUsed_wolb", "Actual Omega Used #Omega/#bar{#Omega}", 2, -0.5, 1.5);
	hOmegaUsed_wolb_sideband = new TH1D("hOmegaUsed_wolb_sideband", "Actual Omega Used #Omega/#bar{#Omega}", 2, -0.5, 1.5);
	hOmegaUsed_wl = new TH1D("hOmegaUsed_wl", "Actual Omega Used #Omega/#bar{#Omega}", 2, -0.5, 1.5);
	hOmegaUsed_wl_sideband = new TH1D("hOmegaUsed_wl_sideband", "Actual Omega Used #Omega/#bar{#Omega}", 2, -0.5, 1.5);
	hOmegaUsed_wol = new TH1D("hOmegaUsed_wol", "Actual Omega Used #Omega/#bar{#Omega}", 2, -0.5, 1.5);
	hOmegaUsed_wol_sideband = new TH1D("hOmegaUsed_wol_sideband", "Actual Omega Used #Omega/#bar{#Omega}", 2, -0.5, 1.5);

	// pT spectrum
	hOmegaM_error_1 = new TH1D("hOmegaM_error_1", "Difference between KFP mass and daughter tracks calculated mass", 1000, -0.1, 0.1);
	hOmegabarM_error_1 = new TH1D("hOmegabarM_error_1", "Difference between KFP mass and daughter tracks calculated mass", 1000, -0.1, 0.1);
	hOmegaM_error_2 = new TH1D("hOmegaM_error_2", "Difference between KFP GetMass and manual four-vector mass", 1000, -0.1, 0.1);
	hOmegabarM_error_2 = new TH1D("hOmegabarM_error_2", "Difference between KFP GetMass and daughter manual four-vector mass", 1000, -0.1, 0.1);
	hDauKaonM_error = new TH1D("hKaonM_error", "Difference between kaon inv mass and PDG mass", 1000, -0.1, 0.1);
	hDauLambdaM_error = new TH1D("hLambdaM_error", "Difference between lambda inv mass and PDG mass", 1000, -0.1, 0.1);

	// for (int i = 0; i < 9; i++)
	// {
	// 	for (int j = 0; j < num_pt_bin; j++) 
	// 	{
	// 		sprintf(temp, "hOmegaM_pt_%d_%d", i+1, j+1);
	// 		hOmegaM_pt[i][j]    = new TH1D(temp, temp, 1400, 1., 2.4);
	// 		sprintf(temp, "hOmegabarM_pt_%d_%d", i+1, j+1);
	// 		hOmegabarM_pt[i][j] = new TH1D(temp, temp, 1400, 1., 2.4);	
	// 		// sprintf(temp, "hOmegaM_rotbkg_pi_pt_%d", i+1);
	// 		// hOmegaM_rotbkg_pi_pt[i]    = new TH1D(temp, temp, 1400, 1., 2.4);
	// 		// sprintf(temp, "hOmegabarM_rotbkg_pi_pt_%d", i+1);
	// 		// hOmegabarM_rotbkg_pi_pt[i] = new TH1D(temp, temp, 1400, 1., 2.4);
	// 	}
	// }

	// Omega v2
	// for (int i = 0; i < 9; i++) // for (int i = 0; i < num_phi_bin; i++)
	// {
	// 	for (int j = 0; j < num_phi_bin; j++) // for (int j = 0; j < num_pt_bin; j++) 
	// 	{
	// 		sprintf(temp, "hOmegaM_phi_%d_%d", i+1, j+1);
	// 		hOmegaM_phi[i][j]    = new TH1D(temp, temp, 1400, 1., 2.4);
	// 		sprintf(temp, "hOmegabarM_phi_%d_%d", i+1, j+1);
	// 		hOmegabarM_phi[i][j] = new TH1D(temp, temp, 1400, 1., 2.4);	
	// 		// sprintf(temp, "hOmegaM_rotbkg_pi_pt_%d", i+1);
	// 		// hOmegaM_rotbkg_pi_pt[i]    = new TH1D(temp, temp, 1400, 1., 2.4);
	// 		// sprintf(temp, "hOmegabarM_rotbkg_pi_pt_%d", i+1);
	// 		// hOmegabarM_rotbkg_pi_pt[i] = new TH1D(temp, temp, 1400, 1., 2.4);
	// 	}
	// }

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
	
	// Run-by-Run QA
	hEPD_full_1_runID = new TProfile("hEPD_full_1_runID", "hEPD_full_1_runID", 40000, 19119999.5, 19159999.5, 0., 2*pi);
	hEPD_full_2_runID = new TProfile("hEPD_full_2_runID", "hEPD_full_2_runID", 40000, 19119999.5, 19159999.5, 0., 2*pi);
	hEPD_full_1_day = new TProfile("hEPD_full_1_day", "hEPD_full_1_day", 300, -0.5, 299.5, 0., 2*pi);
	hEPD_full_2_day = new TProfile("hEPD_full_2_day", "hEPD_full_2_day", 300, -0.5, 299.5, 0., 2*pi);

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

	// hCorrKplusO_y_pT = new TH2D("hCorrKplusO_y_pT", "hCorrKplusO_y_pT", 500, 0.0, 5.0, 200, 0.0, 2.0);
	// hCorrKplusObar_y_pT = new TH2D("hCorrKplusObar_y_pT", "hCorrKplusObar_y_pT", 500, 0.0, 5.0, 200, 0.0, 2.0);
	// hCorrKminusO_y_pT = new TH2D("hCorrKminusO_y_pT", "hCorrKminusO_y_pT", 500, 0.0, 5.0, 200, 0.0, 2.0); 
	// hCorrKminusObar_y_pT = new TH2D("hCorrKminusObar_y_pT", "hCorrKminusObar_y_pT", 500, 0.0, 5.0, 200, 0.0, 2.0);
	hCorrKplusO_y_pT_mixed = new TH2D("hCorrKplusO_y_pT_mixed", "hCorrKplusO_y_pT_mixed", 500, 0.0, 5.0, 200, 0.0, 2.0);
	hCorrKplusObar_y_pT_mixed = new TH2D("hCorrKplusObar_y_pT_mixed", "hCorrKplusObar_y_pT_mixed", 500, 0.0, 5.0, 200, 0.0, 2.0);
	hCorrKminusO_y_pT_mixed = new TH2D("hCorrKminusO_y_pT_mixed", "hCorrKminusO_y_pT_mixed", 500, 0.0, 5.0, 200, 0.0, 2.0); 
	hCorrKminusObar_y_pT_mixed = new TH2D("hCorrKminusObar_y_pT_mixed", "hCorrKminusObar_y_pT_mixed", 500, 0.0, 5.0, 200, 0.0, 2.0);

	// hCorrKplusO_y_phi = new TH2D("hCorrKplusO_y_phi", "hCorrKplusO_y_phi", 500, 0.0, pi, 200, 0.0, 2.0);
	// hCorrKplusObar_y_phi = new TH2D("hCorrKplusObar_y_phi", "hCorrKplusObar_y_phi", 500, 0.0, pi, 200, 0.0, 2.0);
	// hCorrKminusO_y_phi = new TH2D("hCorrKminusO_y_phi", "hCorrKminusO_y_phi", 500, 0.0, pi, 200, 0.0, 2.0); 
	// hCorrKminusObar_y_phi = new TH2D("hCorrKminusObar_y_phi", "hCorrKminusObar_y_phi", 500, 0.0, pi, 200, 0.0, 2.0);

	// hCorrKplusO_phi_pT = new TH2D("hCorrKplusO_phi_pT", "hCorrKplusO_phi_pT", 500, 0.0, 5.0, 500, 0.0, pi);
	// hCorrKplusObar_phi_pT = new TH2D("hCorrKplusObar_phi_pT", "hCorrKplusObar_phi_pT", 500, 0.0, 5.0, 500, 0.0, pi);
	// hCorrKminusO_phi_pT = new TH2D("hCorrKminusO_phi_pT", "hCorrKminusO_phi_pT", 500, 0.0, 5.0, 500, 0.0, pi); 
	// hCorrKminusObar_phi_pT = new TH2D("hCorrKminusObar_phi_pT", "hCorrKminusObar_phi_pT", 500, 0.0, 5.0, 500, 0.0, pi);

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
	hTPCPOIShiftOutput_cos->Write();
	hTPCPOIShiftOutput_sin->Write();
	for (int ewFull = 0; ewFull < 3; ewFull++)
	{
		hTPCEPShiftOutput_cos[ewFull]->Write();
		hTPCEPShiftOutput_sin[ewFull]->Write();
	}
}

void StKFParticleAnalysisMaker::WriteEPDShift()
{
	for (int ewFull = 0; ewFull < 3; ewFull++)
	{
		hEPDEPShiftOutput_1_cos[ewFull]->Write();
		hEPDEPShiftOutput_1_sin[ewFull]->Write();
		hEPDEPShiftOutput_2_cos[ewFull]->Write();
		hEPDEPShiftOutput_2_sin[ewFull]->Write();
	}
}

//-----------------------------------------------------------------------------
void StKFParticleAnalysisMaker::WriteHistograms() {

	///////////////////
	hEventQA  ->Write();
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

	hDauProtonFirstPoint_lam_pt->Write();
	hDauProtonLastPoint_lam_pt ->Write();
	hDauPionFirstPoint_lam_pt  ->Write();
	hDauPionLastPoint_lam_pt   ->Write();

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
	// hCorrKplusO_y_pT  ->Write(); hCorrKplusObar_y_pT  ->Write(); hCorrKminusO_y_pT  ->Write(); hCorrKminusObar_y_pT  ->Write();
	hCorrKplusO_y_pT_mixed->Write(); hCorrKplusObar_y_pT_mixed->Write(); hCorrKminusO_y_pT_mixed->Write(); hCorrKminusObar_y_pT_mixed->Write();
	// hCorrKplusO_y_phi ->Write(); hCorrKplusObar_y_phi ->Write(); hCorrKminusO_y_phi ->Write(); hCorrKminusObar_y_phi ->Write();
	// hCorrKplusO_phi_pT->Write(); hCorrKplusObar_phi_pT->Write(); hCorrKminusO_phi_pT->Write(); hCorrKminusObar_phi_pT->Write();

	hKratio_omega   ->Write();
    hKratio_wo      ->Write();
    hKratio_omegabar->Write();

	hProtony->Write();
	hAntiProtony->Write();
	hm2pion_a->Write();
	hm2pion_r->Write();	
	hm2pion_b->Write();
	hm2proton_a->Write();
	hm2proton_r->Write();
	hm2proton_b->Write();
	hTOFEff_check->Write();
	for (int i = 0; i < 9; i++) hTPCEff_check[i]->Write();

	// v2
	for (int ewFull = 0; ewFull < 3; ewFull++)
	{
		for (int i = 0; i < 9; i++)
		{
			hTPCEP_2[i][ewFull]      ->Write();
			hTPCEP_2_shifted[i][ewFull]->Write();
			hTPCEP_2_2D[i][ewFull]   ->Write();
			hTPCEP_2_2D_shifted[i][ewFull]->Write();
			hTPCEP_2_shift[i][ewFull]->Write();

			hEPDEP_1[i][ewFull]      ->Write();
			hEPDEP_1_shifted[i][ewFull]->Write();
			hEPDEP_1_2D[i][ewFull]   ->Write();
			hEPDEP_1_2D_shifted[i][ewFull]->Write();
			hEPDEP_2[i][ewFull]      ->Write();
			hEPDEP_2_shifted[i][ewFull]->Write();
			hEPDEP_2_2D[i][ewFull]   ->Write();
			hEPDEP_2_2D_shifted[i][ewFull]->Write();
			hEPDEP_2_shift[i][ewFull]->Write();
		}
	}
	hTPCAssoPhi->Write();
	hTPCAssoPhi_shifted->Write();
	hTPCPOIPhi->Write();
	hTPCPOIPhi_shifted->Write();
	hTPCAssoPhi_2D->Write();
	hTPCAssoPhi_2D_shifted->Write();
	hTPCPOIPhi_2D->Write();
	hTPCPOIPhi_2D_shifted->Write();
	hTPCEP_ew_cos->Write();
	hEPDEP_ew_cos_1->Write();
	hEPDEP_ew_cos_2->Write();
	// hEPDEP_ew_cos_1_for_v2->Write();

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
		hOmega_TPC_v2[i]->Write();
		hOmegabar_EPD_v2[i]->Write();
		hOmegabar_TPC_v2[i]->Write();

		hpiplus_EPD_v2_pt[i]->Write();
		hpiminus_EPD_v2_pt[i]->Write();
		hproton_EPD_v2_pt[i]->Write();
		hantiproton_EPD_v2_pt[i]->Write();
		hkplus_EPD_v2_pt[i]->Write();
		hkminus_EPD_v2_pt[i]->Write();
		hpiplus_TPC_v2_pt[i]->Write();
		hpiminus_TPC_v2_pt[i]->Write();
		hproton_TPC_v2_pt[i]->Write();
		hantiproton_TPC_v2_pt[i]->Write();
		hkplus_TPC_v2_pt[i]->Write();
		hkminus_TPC_v2_pt[i]->Write();

		// v1
		hpiplus_EPD_v1_y[i]->Write();
		hpiminus_EPD_v1_y[i]->Write();
		hproton_EPD_v1_y[i]->Write();
		hantiproton_EPD_v1_y[i]->Write();
		hkplus_EPD_v1_y[i]->Write();
		hkminus_EPD_v1_y[i]->Write();
	}
	hpiplus_EPD_v2->Write();
	hpiminus_EPD_v2->Write();
	hproton_EPD_v2->Write();
	hantiproton_EPD_v2->Write();
	hkplus_EPD_v2->Write();
	hkminus_EPD_v2->Write();
	hpiplus_TPC_v2->Write();
	hpiminus_TPC_v2->Write();
	hproton_TPC_v2->Write();
	hantiproton_TPC_v2->Write();
	hkplus_TPC_v2->Write();
	hkminus_TPC_v2->Write();

	for (int i = 0; i < 9; i++)
	{
		hXiM_cen[i]->Write();
		hXibarM_cen[i]->Write();
		hXi_EPD_v2[i]->Write();
		hXibar_EPD_v2[i]->Write();
		hXi_TPC_v2[i]->Write();
		hXibar_TPC_v2[i]->Write();
		hphi_EPD_v2[i]->Write();
		hphi_TPC_v2[i]->Write();

		hphiM_cen[i]->Write();
		hXi1530M_cen[i]->Write();
		hXi1530barM_cen[i]->Write();
		hKsM_cen[i]->Write();
		hOmega2012M_cen[i]->Write();
		hOmega2012barM_cen[i]->Write();

		// for (int ptbin=0; ptbin<40; ptbin++)
		// {
		// 	hLambdaM_cen_pt[i][ptbin]->Write();
		// 	hLambdabarM_cen_pt[i][ptbin]->Write();
		// 	hLambda_EPD_v2_pt[i][ptbin]->Write();
		// 	hLambdabar_EPD_v2_pt[i][ptbin]->Write();
		// 	hLambda_TPC_v2_pt[i][ptbin]->Write();
		// 	hLambdabar_TPC_v2_pt[i][ptbin]->Write();
		// }
	}

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
	hOmegaUsed_wlb->Write();
	hOmegaUsed_wlb_sideband->Write();
	hOmegaUsed_wolb->Write();
	hOmegaUsed_wolb_sideband->Write();
	hOmegaUsed_wl->Write();
	hOmegaUsed_wl_sideband->Write();
	hOmegaUsed_wol->Write();
	hOmegaUsed_wol_sideband->Write();

	// corresponding inv mass
	hOmegaM_wlb->Write();
	hOmegaM_wolb->Write();
	hOmegaM_wl->Write();
	hOmegaM_wol->Write();
	hOmegabarM_wlb->Write();
	hOmegabarM_wolb->Write();
	hOmegabarM_wl->Write();
	hOmegabarM_wol->Write();

	// pT Spectrum
	hOmegaM_error_1->Write();
	hOmegabarM_error_1->Write();
	hOmegaM_error_2->Write();
	hOmegabarM_error_2->Write();
	hDauKaonM_error->Write();
	hDauLambdaM_error->Write();
	// for (int i = 0; i < 9; i++) // centrality
	// {
	// 	for (int j = 0; j < num_pt_bin; j++) 
	// 	{
	// 		hOmegaM_pt[i][j]->Write();           hOmegabarM_pt[i][j]->Write();
	// 		//hOmegaM_rotbkg_pi_pt[i]->Write(); hOmegabarM_rotbkg_pi_pt[i]->Write();
	// 	}
	// }

	// // Omega v2
	// for (int i = 0; i < 9; i++) // for (int i = 0; i < num_phi_bin; i++) // phi - psi bins
	// {
	// 	for (int j = 0; j < num_phi_bin; j++)
	// 	{
	// 		hOmegaM_phi[i][j]->Write();           hOmegabarM_phi[i][j]->Write();
	// 	}
	// }

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

	// for eTOF
	hgetofcrossingY ->Write();
	hgm2_etof	->Write();
	hgpinvbeta_etof->Write();
	hgetofmatchflag->Write();
	hgeta_etof	->Write();

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
	for (int i = 0; i < 9; i++)
	{
		hgpT[i]         ->Write();
		hgpT_TOF[i]     ->Write();
		hgpTeta[i]      ->Write();
		hgpTeta_TOF[i]  ->Write();
	}
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
	hKplusy_omega 	->Write();
	hKplusy_omegabar->Write();
	hKminusy_omega 	->Write();
	hKminusy_omegabar->Write();
	
	// sideband
	hKpluspt_omega_sideband	     ->Write();
	hKpluspt_omegabar_sideband	 ->Write();
	hKminuspt_omega_sideband	 ->Write();
	hKminuspt_omegabar_sideband	 ->Write();
	hKpluseta_omega_sideband	 ->Write();
	hKpluseta_omegabar_sideband	 ->Write();
	hKminuseta_omega_sideband	 ->Write();
	hKminuseta_omegabar_sideband->Write();
	hKplusphi_omega_sideband	 ->Write();
	hKplusphi_omegabar_sideband	 ->Write();
	hKminusphi_omega_sideband	 ->Write();
	hKminusphi_omegabar_sideband ->Write();
	hKplusy_omega_sideband	     ->Write();
	hKplusy_omegabar_sideband	 ->Write();
	hKminusy_omega_sideband	     ->Write();
	hKminusy_omegabar_sideband	 ->Write();

	// w/ and wo/ lambda
	hKpluspt_wlb_omega	 ->Write();
	hKpluspt_wlb_omegabar->Write();
	hKminuspt_wl_omega	 ->Write();
	hKminuspt_wl_omegabar->Write();
	hKplusy_wlb_omega	->Write();
	hKplusy_wlb_omegabar->Write();
	hKminusy_wl_omega	->Write();
	hKminusy_wl_omegabar->Write();
	hKpluspt_wolb_omega   ->Write();
	hKpluspt_wolb_omegabar->Write();
	hKminuspt_wol_omega   ->Write();
	hKminuspt_wol_omegabar->Write();
	hKplusy_wolb_omega   ->Write();
	hKplusy_wolb_omegabar->Write();
	hKminusy_wol_omega   ->Write();
	hKminusy_wol_omegabar->Write();
	hKpluspt_wlb_omega_sideband	  ->Write();
	hKpluspt_wlb_omegabar_sideband->Write();
	hKminuspt_wl_omega_sideband	  ->Write();
	hKminuspt_wl_omegabar_sideband->Write();
	hKplusy_wlb_omega_sideband	 ->Write();
	hKplusy_wlb_omegabar_sideband->Write();
	hKminusy_wl_omega_sideband	 ->Write();
	hKminusy_wl_omegabar_sideband->Write();
	hKpluspt_wolb_omega_sideband	->Write();
	hKpluspt_wolb_omegabar_sideband->Write();
	hKminuspt_wol_omega_sideband   ->Write();
	hKminuspt_wol_omegabar_sideband->Write();
	hKplusy_wolb_omega_sideband	  ->Write();
	hKplusy_wolb_omegabar_sideband->Write();
	hKminusy_wol_omega_sideband	  ->Write();
	hKminusy_wol_omegabar_sideband->Write();

	// for the rapidity check
	hkplus_totaly ->Write();
	hkminus_totaly->Write();
	hkplus_ntrack ->Write();
	hkminus_ntrack->Write();

	// Run-by-Run QA
	hEPD_full_1_runID->Write();
	hEPD_full_2_runID->Write();
	hEPD_full_1_day  ->Write();
	hEPD_full_2_day  ->Write();
	
	for (int i = 0; i < 10; i++)
	{	
		hgPID2D_proton_pt[i]->Write();
		hgPID2D_antiproton_pt[i]->Write();
		hgPID2D_piplus_pt[i]->Write();
		hgPID2D_piminus_pt[i]->Write();
		hgPID2D_kplus_pt[i]->Write();
		hgPID2D_kminus_pt[i]->Write();
	}

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
// bool StKFParticleAnalysisMaker::readBadList()
// { 
// 	if(0) return kTRUE;
// 	TString inf=mListDir + "/badList/";
// 	inf += Form("badrun%dList%.1f.list",mRun,mEnergy);
// 	ifstream inrun; 
// 	inrun.open(inf);
// 	if ( inrun.fail() ) {
// 		cout<< "cannot open " << inf.Data() << endl;
// 		return kFALSE;
// 	}
// 	Int_t runid;
// 	while ( inrun >> runid ) { badList.push_back(runid); }
// 	inrun.close();
// 	sort(badList.begin(),badList.end());

// 	vector<int>::iterator it;
// 	it = std::unique (badList.begin(), badList.end());
// 	badList.resize( std::distance(badList.begin(),it) );

// 	cout <<"badrun list :" <<inf.Data() << " loaded." << endl;
// 	cout <<"Total       :" <<badList.size()<< " bad runs. "<< endl;
// 	return kTRUE;
// }

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
	hEventQA->Fill(0.);

    PicoDst = StPicoDst::instance(); 		
	StPicoDst* mPicoDst = PicoDst;
	if(!mPicoDst) {
		LOG_WARN << " No PicoDst! Skip! " << endm;
		return kStOK;
	}

	hEventQA->Fill(1.);
	//     pass event  
	/////////////////////////////////////////////////////////
	StPicoEvent* mEvent= (StPicoEvent*) mPicoDst->event(); 
	if(!mEvent)return kStOK;

	hEventQA->Fill(2.);
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

	hEventQA->Fill(3.);
	const TVector3 Vertex3D=mEvent->primaryVertex();
	const double VertexX = Vertex3D.x(); 
	const double VertexY = Vertex3D.y(); 
	const double VertexZ = Vertex3D.z(); 
	const double VertexR = sqrt(VertexX*VertexX + VertexY*VertexY);
	const double vpdVz   = mEvent->vzVpd();

	//event cut
	//if(refMult <=2 || refMult > 1000) return kStOK;
	mRefMultCorr->init(runID);
	mRefMultCorr->initEvent(refMult,VertexZ,mEvent->ZDCx());
	if(removeBadID(runID)) return kStOK;            
	if(mRefMultCorr->isBadRun(runID)) return kStOK; // reject bad run of StRefMultCorr
	if(!mRefMultCorr->passnTofMatchRefmultCut(1.*refMult, 1.*tofMatch)) return kStOK; // reject pileup of StRefMultCorr

	hEventQA->Fill(4.);
	hVertex2D ->Fill(VertexZ,vpdVz);
	hDiffVz   ->Fill(VertexZ-vpdVz); 

	const double DVz = VertexZ-vpdVz;

	if(fabs(VertexZ) > 70) return kStOK; 
	if(sqrt(pow(VertexX,2.)+pow(VertexY,2.))>2.0) return kStOK; 
	// if(fabs(VertexZ-vpdVz)>4.) return kStOK;       // no vpd cut in low energy? Yes!!

	//check run number
	hEventQA->Fill(5.);
	int runnumberPointer = -999;
	runnumberPointer = CheckrunNumber(runID);
	if(runnumberPointer==-999) return kStOK;

	int dayPointer = (int)((runID)/1000%1000);
	int yearPointer = (int)((runID)/1000000);
	int mRunL=runList.at(0);
	int mDayL=(int) ((mRunL)/1000%1000);
	int mYearL=(int) ((mRunL)/1000000);
	yearPointer -= mYearL;
	dayPointer = yearPointer == 0 ? dayPointer - mDayL : dayPointer + 365 - mDayL;
	int timePointer = dayPointer/mStps;

	//StRefMultCorr
	double refmultWght = mRefMultCorr->getWeight();
	double refmultCorr = mRefMultCorr->getRefMultCorr() ;
	int centrality     = mRefMultCorr->getCentralityBin9();  // 0 - 8  be careful !!!!!!!! 
	//double mWght    = 1;
	//double mult_corr= refMult;
	//int centrality  = findCentrality(mult_corr)-1;  // 0-8 
	if( centrality<0||centrality>=(nCent-1)) return kStOK;
	int cent = centrality+1;  
	
	hEventQA->Fill(6.);
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
	int runID2 = runID%1000000;
	int runID2bin = (runID2-mRunStart)/10+1;
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
	hEPD_full_1_runID->Fill(runID, result.FullPhiWeightedAndShiftedPsi(1));
	hEPD_full_2_runID->Fill(runID, result.FullPhiWeightedAndShiftedPsi(2));
	hEPD_full_1_day->Fill(dayPointer, result.FullPhiWeightedAndShiftedPsi(1));
	hEPD_full_2_day->Fill(dayPointer, result.FullPhiWeightedAndShiftedPsi(2));

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
	std::vector<KFParticle> LambdaVecAll;
	bool hasLambda = false, hasLambdabar = false;
	std::vector<KFParticle> XiVecAll;
	std::vector<KFParticle> phiVecAll;
	std::vector<KFParticle> KsVecAll;
	std::vector<KFParticle> Xi1530VecAll;
	std::vector<int> OmegaDauPionIdVec;
	std::vector<int> OmegaDauProtonIdVec;
	mix_px.resize(0); mix_py.resize(0); mix_pz.resize(0); mix_charge.resize(0); 
	for (int iKFParticle=0; iKFParticle < KFParticlePerformanceInterface->GetNReconstructedParticles(); iKFParticle++)
	{ 
		const KFParticle particle = KFParticleInterface->GetParticles()[iKFParticle]; 
		bool IsOmega = false, IsLambda = false, IsXi = false, Isphi = false, IsKs = false, IsXi1530 = false;
		if      (fabs(particle.GetPDG()) == OmegaPdg)  IsOmega  = true; 
		else if (fabs(particle.GetPDG()) == LambdaPdg) IsLambda = true; 
		else if (fabs(particle.GetPDG()) == XiPdg)     IsXi     = true;
		else if (fabs(particle.GetPDG()) == phiPdg)    Isphi	= true;
		else if (fabs(particle.GetPDG()) == KsPdg)     IsKs     = true;
		else if (fabs(particle.GetPDG()) == Xi1530Pdg) IsXi1530 = true;
		else continue;
		int upQ; if (particle.GetPDG() > 0) upQ = 1; else if (particle.GetPDG() < 0) upQ = -1; else continue;
		if (IsOmega) OmegaVecAll.push_back(particle);
		if (IsXi) XiVecAll.push_back(particle);
		if (IsLambda) LambdaVecAll.push_back(particle);
		if (Isphi) phiVecAll.push_back(particle);
		if (IsKs) KsVecAll.push_back(particle);
		if (IsXi1530) Xi1530VecAll.push_back(particle);
		if (IsOmega && isGoodOmega(cent, particle)) OmegaVec.push_back(particle);
		if (IsLambda && isGoodLambda(cent, particle) && upQ ==  1) hasLambda    = true;
		if (IsLambda && isGoodLambda(cent, particle) && upQ == -1) hasLambdabar = true;
		if (IsOmega && isSidebandOmega(cent, particle)) OmegaSidebandVec.push_back(particle);
			
		// // sec tracks first/last points check
		// if (IsLambda && isGoodLambda(cent, particle)) SetDaughterTrackPointers(iKFParticle);

		// Omega/lambda QA
		if (IsOmega)
		{
			if (upQ == 1)
			{
				hOmegaM  ->Fill(particle.GetMass());
				hOmegaM_cen[cent-1]->Fill(particle.GetMass());

				// pT spectrum
				// if (fabs(particle.GetRapidity()) < 0.5)
				// {
				// 	for (int i = 0; i < num_pt_bin; i++)
				// 	{
				// 		if (particle.GetPt() >= StKFParticleAnalysisMaker::OmegaMassPtLowerBin[i] && 
				// 			particle.GetPt() <= StKFParticleAnalysisMaker::OmegaMassPtLowerBin[i+1])
				// 			hOmegaM_pt[cent-1][i]->Fill(particle.GetMass());
				// 	}
				// }

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
				// hOmegaM_phi[cent-1][phi_bin]->Fill(particle.GetMass());

				if (!isGoodOmega(cent, particle)) continue; // subject to change
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
				// TLorentzVector dauKaon, dauLambda;
				// TLorentzVector dauKaonPico, dauLambdaPico;
				// for (int iDaughter = 0; iDaughter < particle.NDaughters(); iDaughter++)
				// {
				// 	const int daughterId = particle.DaughterIds()[iDaughter];
				// 	const KFParticle daughter = KFParticleInterface->GetParticles()[daughterId];
				// 	const int daughterTrackId = daughter.DaughterIds()[0];
				// 	int trackIndex = trackMap[daughterTrackId];   
				// 	StPicoTrack *daughterTrack = mPicoDst->track(trackIndex);
				// 	if (!daughterTrack) continue;
				// 	hOmegaDauPid->Fill(1.0*daughter.GetPDG());

				// 	// for post-reconstruction DCA cut
				// 	TVector3 pDaughter(daughter.GetPx(), daughter.GetPy(), daughter.GetPz());
				// 	TVector3 pDaughterPico = daughterTrack->gMom();
				// 	TVector3 xDaughter(daughter.GetX(), daughter.GetY(), daughter.GetZ());
				// 	TVector3 xDaughterPico = daughterTrack->origin();
					
				// 	if (daughter.GetPDG() == -KaonPdg)
				// 	{
				// 		//dauKaon.SetVectM(pDaughter, KaonPdgMass);
				// 		hDauKminusp  ->Fill(daughter.GetMomentum());
				// 		hDauKminuspt ->Fill(daughter.GetPt());
				// 		hDauKminusy  ->Fill(daughter.GetRapidity());
				// 		hDauKminusphi->Fill(daughter.GetPhi());
				// 		hDauKminusnSigma->Fill(daughterTrack->nSigmaKaon());

				// 		StPicoPhysicalHelix helixDaughter(pDaughterPico, xDaughterPico, magnet*kilogauss, daughter.GetQ());
				// 		pair<double, double> tmps = helixDaughter.pathLengths(helixOmega);
				// 		TVector3 ox1 = helixDaughter.at(tmps.first); 
				// 		TVector3 ox2 = helixOmega.at(tmps.second);
				// 		daughter.TransportToProductionVertex();
				// 		dauKaon.SetXYZT(daughter.GetPx(), daughter.GetPy(), daughter.GetPz(), daughter.GetE());
				// 		hDauKaonM_error->Fill(dauKaon.M() - KaonPdgMass);
				// 		dauKaon.SetXYZM(daughter.GetPx(), daughter.GetPy(), daughter.GetPz(), KaonPdgMass);
				// 		dauKaonPico.SetVectM(ox1, KaonPdgMass); // use all StPicoTrack information
				// 		//dauKaon.SetVectM(ox1, KaonPdgMass); // transport daughter to Omega DCA
				// 		double dca = (ox1 - ox2).Mag();
				// 		CutDecider(particle, hDCAOtoK_signal, hDCAOtoK_sideband, dca);
				// 	}
				// 	else 
				// 	{
				// 		//dauLambda.SetVectM(pDaughter, LambdaPdgMass);
				// 		hDauLambdaM  ->Fill(daughter.GetMass());
				// 		hDauLambdap  ->Fill(daughter.GetMomentum());
				// 		hDauLambdapt ->Fill(daughter.GetPt());
				// 		hDauLambdaphi->Fill(daughter.GetPhi());
				// 		hDauLambday  ->Fill(daughter.GetRapidity());
				// 		hDauLambdaDL ->Fill((xDaughter-Vertex3D).Mag());
				// 		hOmegaDLminusLambdaDL->Fill((xDaughter-Vertex3D).Mag()-(xOmega-Vertex3D).Mag());

				// 		StPicoPhysicalHelix helixDaughter(pDaughterPico*(1./pDaughterPico.Perp())*100, xDaughterPico, magnet*kilogauss, 1);
				// 		pair<double, double> tmps = helixOmega.pathLengths(helixDaughter);
				// 		TVector3 ox1 = helixOmega.at(tmps.first);
				// 		TVector3 ox2 = helixDaughter.at(tmps.second);
				// 		daughter.TransportToProductionVertex();
				// 		dauLambda.SetXYZT(daughter.GetPx(), daughter.GetPy(), daughter.GetPz(), daughter.GetE());
				// 		hDauLambdaM_error->Fill(dauLambda.M() - LambdaPdgMass);
				// 		dauLambdaPico.SetVectM(ox2, LambdaPdgMass);
				// 		double dca = (ox1 - ox2).Mag();
				// 		CutDecider(particle, hDCAOtoL_signal, hDCAOtoL_sideband, dca);

				// 		for (int iGrandDaughter = 0; iGrandDaughter < daughter.NDaughters(); iGrandDaughter++)
				// 		{
				// 			const int granddaughterId = daughter.DaughterIds()[iGrandDaughter];
				// 			const KFParticle granddaughter = KFParticleInterface->GetParticles()[granddaughterId];
				// 			if (daughter.GetPDG() == ProtonPdg)
				// 			{
				// 				hDauProtonp  ->Fill(granddaughter.GetMomentum());
				// 				hDauProtonpt ->Fill(granddaughter.GetPt());
				// 			}
				// 			else if (daughter.GetPDG() == PionPdg)
				// 			{
				// 				hDauPionbarp  ->Fill(granddaughter.GetMomentum());
				// 				hDauPionbarpt ->Fill(granddaughter.GetPt());
				// 			}
				// 		}
				// 	}
				// }
				// hOmegaM_error_1->Fill(particle.GetMass()-(dauKaon+dauLambda).M());
				// hOmegaM_error_2->Fill(particle.GetMass()-OmegaLorentz.M());
				// dauKaon.RotateZ(pi);

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
				// if (fabs(particle.GetRapidity()) < 0.5)
				// {
				// 	for (int i = 0; i < num_pt_bin; i++)
				// 	{
				// 		if (particle.GetPt() >= StKFParticleAnalysisMaker::OmegaMassPtLowerBin[i] && 
				// 			particle.GetPt() <= StKFParticleAnalysisMaker::OmegaMassPtLowerBin[i+1])
				// 			hOmegabarM_pt[cent-1][i]->Fill(particle.GetMass());
				// 	}
				// }

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
				// hOmegabarM_phi[cent-1][phi_bin]->Fill(particle.GetMass());

				if (!isGoodOmega(cent, particle)) continue; // subject to change
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
				// TLorentzVector dauKaon, dauLambda;
				// TLorentzVector dauKaonPico, dauLambdaPico;
				// for (int iDaughter = 0; iDaughter < particle.NDaughters(); iDaughter++)
				// {
				// 	const int daughterId = particle.DaughterIds()[iDaughter];
				// 	const KFParticle daughter = KFParticleInterface->GetParticles()[daughterId];
				// 	const int daughterTrackId = daughter.DaughterIds()[0];
				// 	int trackIndex = trackMap[daughterTrackId];   
				// 	StPicoTrack *daughterTrack = mPicoDst->track(trackIndex);		
				// 	if (!daughterTrack) continue;			
				// 	hOmegabarDauPid->Fill(1.0*daughter.GetPDG());

				// 	// for post-reconstruction DCA cut
				// 	TVector3 pDaughter(daughter.GetPx(), daughter.GetPy(), daughter.GetPz());
				// 	TVector3 pDaughterPico = daughterTrack->gMom();
				// 	TVector3 xDaughter(daughter.GetX(), daughter.GetY(), daughter.GetZ());
				// 	TVector3 xDaughterPico = daughterTrack->origin();

				// 	if (daughter.GetPDG() ==  KaonPdg)
				// 	{
				// 		//dauKaon.SetVectM(pDaughter, KaonPdgMass);
				// 		hDauKplusp  ->Fill(daughter.GetMomentum());
				// 		hDauKpluspt ->Fill(daughter.GetPt());
				// 		hDauKplusy  ->Fill(daughter.GetRapidity());
				// 		hDauKplusphi->Fill(daughter.GetPhi());
				// 		hDauKplusnSigma->Fill(daughterTrack->nSigmaKaon());

				// 		StPicoPhysicalHelix helixDaughter(pDaughterPico, xDaughterPico, magnet*kilogauss, daughter.GetQ());
				// 		pair<double, double> tmps = helixDaughter.pathLengths(helixOmega);
				// 		TVector3 ox1 = helixDaughter.at(tmps.first);
				// 		TVector3 ox2 = helixOmega.at(tmps.second);
				// 		daughter.TransportToProductionVertex();
				// 		dauKaon.SetXYZT(daughter.GetPx(), daughter.GetPy(), daughter.GetPz(), daughter.GetE());
				// 		hDauKaonM_error->Fill(dauKaon.M() - KaonPdgMass);
				// 		dauKaon.SetXYZM(daughter.GetPx(), daughter.GetPy(), daughter.GetPz(), KaonPdgMass);
				// 		dauKaonPico.SetVectM(ox1, KaonPdgMass);
				// 		//dauKaon.SetVectM(ox1, KaonPdgMass);
				// 		double dca = (ox1 - ox2).Mag();
				// 		CutDecider(particle, hDCAOtoK_signal, hDCAOtoK_sideband, dca);
				// 	}
				// 	else 
				// 	{
				// 		//dauLambda.SetVectM(pDaughter, LambdaPdgMass);
				// 		hDauLambdabarM  ->Fill(daughter.GetMass());
				// 		hDauLambdabarp  ->Fill(daughter.GetMomentum());
				// 		hDauLambdabarpt ->Fill(daughter.GetPt());
				// 		hDauLambdabarphi->Fill(daughter.GetPhi());
				// 		hDauLambdabary  ->Fill(daughter.GetRapidity());
				// 		hDauLambdabarDL ->Fill((xDaughter-Vertex3D).Mag());
				// 		hOmegabarDLminusLambdaDL->Fill((xDaughter-Vertex3D).Mag()-(xOmega-Vertex3D).Mag());

				// 		StPicoPhysicalHelix helixDaughter(pDaughterPico*(1./pDaughterPico.Perp())*100, xDaughterPico, magnet*kilogauss, 1);
				// 		pair<double, double> tmps = helixOmega.pathLengths(helixDaughter);
				// 		TVector3 ox1 = helixOmega.at(tmps.first);
				// 		TVector3 ox2 = helixDaughter.at(tmps.second);
				// 		daughter.TransportToProductionVertex();
				// 		dauLambda.SetXYZT(daughter.GetPx(), daughter.GetPy(), daughter.GetPz(), daughter.GetE());
				// 		hDauLambdaM_error->Fill(dauLambda.M() - LambdaPdgMass);
				// 		dauLambdaPico.SetVectM(ox2, LambdaPdgMass);
				// 		//dauLambda.SetVectM(ox2, LambdaPdgMass);
				// 		double dca = (ox1 - ox2).Mag();
				// 		CutDecider(particle, hDCAOtoL_signal, hDCAOtoL_sideband, dca);

				// 		for (int iGrandDaughter = 0; iGrandDaughter < daughter.NDaughters(); iGrandDaughter++)
				// 		{
				// 			const int granddaughterId = daughter.DaughterIds()[iGrandDaughter];
				// 			const KFParticle granddaughter = KFParticleInterface->GetParticles()[granddaughterId];
				// 			if (daughter.GetPDG() == -ProtonPdg)
				// 			{
				// 				hDauProtonbarp  ->Fill(granddaughter.GetMomentum());
				// 				hDauProtonbarpt ->Fill(granddaughter.GetPt());
				// 			}
				// 			else if (daughter.GetPDG() == -PionPdg)
				// 			{
				// 				hDauPionp  ->Fill(granddaughter.GetMomentum());
				// 				hDauPionpt ->Fill(granddaughter.GetPt());
				// 			}
				// 		}
				// 	}
				// }
				// hOmegabarM_error_1->Fill(particle.GetMass()-(dauKaon+dauLambda).M());
				// hOmegabarM_error_2->Fill(particle.GetMass()-OmegaLorentz.M());
				// dauKaon.RotateZ(pi);

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
		else if (IsXi)
		{
			if (upQ == 1) hXiM_cen   [cent-1]->Fill(particle.GetMass());	
			else          hXibarM_cen[cent-1]->Fill(particle.GetMass());
		}
		else if (IsKs) hKsM_cen[cent-1]->Fill(particle.GetMass());
		else if (Isphi) hphiM_cen[cent-1]->Fill(particle.GetMass());
		else if (IsXi1530)
		{
			if (upQ == 1) hXi1530M_cen   [cent-1]->Fill(particle.GetMass());	
			else          hXi1530barM_cen[cent-1]->Fill(particle.GetMass());
		}
	} // End loop over KFParticles

	// fill Omega inv mass for with Lambda and without Lambda
	for (int iomega=0; iomega < OmegaVecAll.size(); iomega++)
	{
		const KFParticle particle = OmegaVecAll[iomega];
		if (hasLambdabar) 
		{
			if (particle.GetQ() < 0) hOmegaM_wlb   ->Fill(particle.GetMass()); 
			else                     hOmegabarM_wlb->Fill(particle.GetMass());
		}
		else
		{
			if (particle.GetQ() < 0) hOmegaM_wolb   ->Fill(particle.GetMass()); 
			else                     hOmegabarM_wolb->Fill(particle.GetMass());
		}
		if (hasLambda)
		{
			if (particle.GetQ() < 0) hOmegaM_wl   ->Fill(particle.GetMass()); 
			else                     hOmegabarM_wl->Fill(particle.GetMass());
		}
		else
		{
			if (particle.GetQ() < 0) hOmegaM_wol   ->Fill(particle.GetMass()); 
			else                     hOmegabarM_wol->Fill(particle.GetMass());
		}
	}

	// try reconstructing Omega(2012)
	for (int iXi=0; iXi < XiVecAll.size(); iXi++)
	{
		for (int iKs=0; iKs < KsVecAll.size(); iKs++)
		{
			const KFParticle Xi = XiVecAll[iXi];
			const KFParticle Ks = KsVecAll[iKs];
			KFParticle Omega2012;
			Omega2012 += Xi;
			Omega2012 += Ks;
			// if (Omega2012.GetChi2()/Omega2012.GetNDF() > 30) continue;
			if (Xi.GetPDG() > 0) hOmega2012M_cen[cent-1]->Fill(Omega2012.GetMass());
			else 				 hOmega2012barM_cen[cent-1]->Fill(Omega2012.GetMass());
		}
	}
	hNumOmega->Fill(OmegaVec.size());

	// correlation function loop  
	if (!PerformAnalysis) return kStOK;
	my_event current_event(OmegaVec);
	Efficiency eff_finder;
  	Int_t nTracks = mPicoDst->numberOfTracks();
	std::vector<int> kaon_tracks; kaon_tracks.resize(0);
	std::vector<int> pion_tracks; pion_tracks.resize(0);
	std::vector<int> proton_tracks; proton_tracks.resize(0);
	float Qx2w = 0, Qy2w = 0, Qx2e = 0, Qy2e = 0; // for EP 
	std::vector<int> track_index;
	int pct = 0; int pbct = 0; // counting protons
	int kpct = 0; int kmct = 0; // counting kaons
	bool fill_nomega = false;
	bool fill_nomega_sideband = false;
	for (Int_t iTrack = 0; iTrack < nTracks; iTrack++) 
	{
    	StPicoTrack *track = mPicoDst->track(iTrack);
    	if (! track)            continue;
    	if (! track->charge())  continue;
    	if (  track->nHitsFit() < 15) continue;
		if (  track->nHitsDedx() < 15) continue;
		if (  track->nHitsFit()*1.0 / track->nHitsMax() < 0.52 || track->nHitsFit()*1.0 / track->nHitsMax() > 1.05) continue;
		//if (  track->dEdxError() < 0.04 || track->dEdxError() > 0.12) continue; // same as kfp
		if (! track->isPrimary()) continue;
		track_index.push_back(iTrack);

		// track info
		float p = track->gMom().Mag();
		float pt = track->gMom().Perp();
		float phi = track->gMom().Phi();
		float eta = track->gMom().Eta();
		float dcatopv = track->gDCA(Vertex3D).Mag();
		float nSigmaKaon = track->nSigmaKaon();
		float nSigmaPion = track->nSigmaPion();
		float nSigmaProton = track->nSigmaProton();
		if (hTOFEff[8] != 0) hTOFEff_check->Fill(pt, hTOFEff[8]->GetEfficiency(hTOFEff[8]->FindFixBin(pt)));

		float pt_prim = track->pMom().Perp();
		float phi_prim = track->pMom().Phi();
		float eta_prim = track->pMom().Eta();

		// EP
		bool isGoodAsso = true;
		// isGoodAsso &= track->isPrimary();
		isGoodAsso &= pT_asso_lo < pt && pt < pT_asso_hi;
		if (isGoodAsso)
		{
		// fill phi shift for asso
			if (fabs(eta_prim) > TPCAssoEtaCut)
		{
			for (int i = 1; i <= shift_order_asso; i++) // fill correction for output
			{
				hTPCAssoShiftOutput_cos->Fill(dayPointer, i, cent-1, cos(i*1.0*phi_prim), mWght);
				hTPCAssoShiftOutput_sin->Fill(dayPointer, i, cent-1, sin(i*1.0*phi_prim), mWght);
			}
			// shift asso phi
				float phi_shifted_asso = ShiftAssoPhi(phi_prim, dayPointer+1, cent);
				hTPCAssoPhi->Fill(phi_prim, mWght);
				hTPCAssoPhi_shifted->Fill(phi_shifted_asso, mWght);
				hTPCAssoPhi_2D->Fill(dayPointer, phi_prim, mWght);
				hTPCAssoPhi_2D_shifted->Fill(dayPointer, phi_shifted_asso, mWght);
				if (eta_prim > TPCAssoEtaCut) // construct east EP
			{
				Qx2e += pt*cos(2*phi_shifted_asso);
				Qy2e += pt*sin(2*phi_shifted_asso);
			}
				else if (eta_prim < -1.0*TPCAssoEtaCut) //construct west EP
			{
				Qx2w += pt*cos(2*phi_shifted_asso);
				Qy2w += pt*sin(2*phi_shifted_asso);
			}
		}
		}
		bool isGoodPOI = true;
		// isGoodPOI &= track->isPrimary();
		isGoodPOI &= pT_trig_lo < pt && pt < pT_trig_hi;
		isGoodPOI &= fabs(eta_prim) < eta_trig_cut;
		if (isGoodPOI)
		{
		// fill phi shift for POI
		for (int i = 1; i <= shift_order_asso; i++) // fill correction for output
		{
				hTPCPOIShiftOutput_cos->Fill(dayPointer, i, cent-1, cos(i*1.0*phi_prim), mWght);
				hTPCPOIShiftOutput_sin->Fill(dayPointer, i, cent-1, sin(i*1.0*phi_prim), mWght);
		}
		// shift POI phi
			float phi_shifted_POI = ShiftPOIPhi(phi_prim, dayPointer+1, cent);
			hTPCPOIPhi->Fill(phi_prim, mWght);
			hTPCPOIPhi_shifted->Fill(phi_shifted_POI, mWght);
			hTPCPOIPhi_2D->Fill(dayPointer, phi_prim, mWght);
			hTPCPOIPhi_2D_shifted->Fill(dayPointer, phi_shifted_POI, mWght);
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
		hgpT[cent-1]      ->Fill(pt);
		hgpTeta[cent-1]   ->Fill(pt, eta);
		if (hasTOF) 
		{
			hgpT_TOF[cent-1]->Fill(pt);
			hgpTeta_TOF[cent-1]->Fill(pt, eta);
		}
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
			double zTOF_proton = 1/beta - sqrt(ProtonPdgMass*ProtonPdgMass/pkaon.Mag2()+1);
			double zTOF_pion   = 1/beta - sqrt(PionPdgMass*PionPdgMass/pkaon.Mag2()+1);
			double zTOF_kaon   = 1/beta - sqrt(KaonPdgMass*KaonPdgMass/pkaon.Mag2()+1);
			if (ptbin >= 0 && ptbin <= 9 && track->charge() > 0) hgPID2D_proton_pt    [ptbin]->Fill(nSigmaProton, zTOF_proton);
			if (ptbin >= 0 && ptbin <= 9 && track->charge() < 0) hgPID2D_antiproton_pt[ptbin]->Fill(nSigmaProton, zTOF_proton);
			if (ptbin >= 0 && ptbin <= 9 && track->charge() > 0) hgPID2D_piplus_pt    [ptbin]->Fill(nSigmaPion  , zTOF_pion  );
			if (ptbin >= 0 && ptbin <= 9 && track->charge() < 0) hgPID2D_piminus_pt   [ptbin]->Fill(nSigmaPion  , zTOF_pion  );
			if (ptbin >= 0 && ptbin <= 9 && track->charge() > 0) hgPID2D_kplus_pt     [ptbin]->Fill(nSigmaKaon  , zTOF_kaon  );
			if (ptbin >= 0 && ptbin <= 9 && track->charge() < 0) hgPID2D_kminus_pt    [ptbin]->Fill(nSigmaKaon  , zTOF_kaon  );
		}

		// primary proton cut for coalescence test
		bool proton_cut = true;
		if (pt < pT_trig_lo || pt > pT_trig_hi) proton_cut = false; 
		if (fabs(eta_prim) > eta_trig_cut) proton_cut = false;
		ProtonPID proton_pid(0., nSigmaProton, pt); // not using zTOF
		if ((pt > proton_pT_lo && pt < proton_pT_TOFth) && hasTOF) // test efficacy of ProtonPID.h
		{
			hm2proton_b->Fill(m2);
			if (proton_pid.IsProtonSimple(2., track->charge())) hm2proton_a->Fill(m2);
			if (fabs(nSigmaProton) < 2.)       hm2proton_r->Fill(m2);
		}
		if (!hasTOF && pt > proton_pT_TOFth) proton_cut = false;
		if (pt > proton_pT_TOFth && (m2 > proton_m2_hi || m2 < proton_m2_lo)) proton_cut = false;
		if (!proton_pid.IsProtonSimple(2., track->charge())) proton_cut = false; // only 0.2 < pt < 2.0!!!
		// if (fabs(nSigmaProton) > 3) proton_cut = false;
		if (dcatopv > dcatoPV_hi) proton_cut = false;
		if (proton_cut)
		{	
			TLorentzVector lv_proton; lv_proton.SetVectM(track->gMom(), ProtonPdgMass);
			if (track->charge() > 0) {hProtony    ->Fill(lv_proton.Rapidity()); if (fabs(lv_proton.Rapidity()) < 0.5) pct++; }
			else					 {hAntiProtony->Fill(lv_proton.Rapidity()); if (fabs(lv_proton.Rapidity()) < 0.5) pbct++;}

			proton_tracks.push_back(iTrack);
		}

		// primary pion cut for coalescence test
		bool pion_cut = true;
		if (pt < pT_trig_lo || pt > pT_trig_hi) pion_cut = false; // use p < 2
		if (fabs(eta_prim) > eta_trig_cut) pion_cut = false;
		PionPID pion_pid(0., nSigmaPion, pt); // not using zTOF
		if ((pt > pion_pT_lo && pt < pion_pT_TOFth) && hasTOF) // test efficacy of ProtonPID.h
		{
			hm2pion_b->Fill(m2);
			if (pion_pid.IsPionSimple(2., track->charge())) hm2pion_a->Fill(m2);
			if (fabs(nSigmaPion) < 2.)	   hm2pion_r->Fill(m2);				   
		}
		if (!hasTOF && pt > pion_pT_TOFth) pion_cut = false;
		if (pt > pion_pT_TOFth && (m2 > pion_m2_hi || m2 < pion_m2_lo)) pion_cut = false;
		if (!pion_pid.IsPionSimple(2., track->charge())) pion_cut = false; // only up to pt < 1.8!!!
		// if (fabs(nSigmaPion) > 3) pion_cut = false;
		if (dcatopv > dcatoPV_hi) pion_cut = false;
		if (pion_cut) pion_tracks.push_back(iTrack);

		// primary kaon cut
		/******** looser cut ********/
		bool kaon_cut = true;
		if (pt < pT_trig_lo || pt > pT_trig_hi) kaon_cut = false; // use p < 1.6
		if (fabs(eta_prim) > eta_trig_cut) kaon_cut = false;
		if (!hasTOF && pt > 0.4) kaon_cut = false;
		if (pt > 0.4 && (m2 > 0.34 || m2 < 0.15)) kaon_cut = false;
		double zTOF = 1/beta - sqrt(KaonPdgMass*KaonPdgMass/pkaon.Mag2()+1);

		KaonPID decider(zTOF, nSigmaKaon, pt);
		if (!decider.IsKaonSimple(2., track->charge())) kaon_cut = false;
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
			float TOFEff = 1.0;
			float TOFEff2D = 1.0;
			float TPCEff = track->charge() > 0? eff_finder.GetEfficiency1D(pt, cent, "Kp") : eff_finder.GetEfficiency1D(pt, cent, "Km");
			float TPCEff2D = track->charge() > 0? eff_finder.GetEfficiency2D(pt, eta, cent, "Kp") : eff_finder.GetEfficiency2D(pt, eta, cent, "Km");
			if (hTOFEff[cent-1] != 0 && pt > 0.4) TOFEff = hTOFEff[cent-1]->GetEfficiency(hTOFEff[cent-1]->FindFixBin(pt));
			if (hTOFEff_2D[cent-1] != 0 && pt > 0.4) TOFEff2D = hTOFEff_2D[cent-1]->GetEfficiency(hTOFEff_2D[cent-1]->FindFixBin(pt, eta));
			hgKpdEdx    ->Fill(pkaon.Mag(), track->dEdx());
			hgKp       ->Fill(p);
			hgKpT      ->Fill(pt);
			hgKDCAtoPV ->Fill(dcatopv);
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
			if (cent < min_cent || cent > max_cent) continue;
			const int kaonindex = track->id();
			for (int iOmega=0; iOmega < OmegaVec.size(); iOmega++)
			{ 
				const KFParticle particle = OmegaVec[iOmega]; 
				if (IsKaonOmegaDaughter(particle, kaonindex)) continue;

				// couting Omega
				if (!fill_nomega)
				{
					if (particle.GetPDG() > 0) hOmegaUsed->Fill(0.);
					if (particle.GetPDG() < 0) hOmegaUsed->Fill(1.);
					if ( hasLambdabar && particle.GetPDG() > 0) hOmegaUsed_wlb->Fill(0);
					if ( hasLambdabar && particle.GetPDG() < 0) hOmegaUsed_wlb->Fill(1);
					if (!hasLambdabar && particle.GetPDG() > 0) hOmegaUsed_wolb->Fill(0);
					if (!hasLambdabar && particle.GetPDG() < 0) hOmegaUsed_wolb->Fill(1);
					if ( hasLambda && particle.GetPDG() > 0) hOmegaUsed_wl->Fill(0);
					if ( hasLambda && particle.GetPDG() < 0) hOmegaUsed_wl->Fill(1);
					if (!hasLambda && particle.GetPDG() > 0) hOmegaUsed_wol->Fill(0);
					if (!hasLambda && particle.GetPDG() < 0) hOmegaUsed_wol->Fill(1);
				}

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
				TLorentzVector lv2_ori; lv2_ori.SetVectM(track->gMom(), KaonPdgMass);
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
				
				float eff = TPCEff * TOFEff;
				float eff2D = TPCEff2D * TOFEff2D;
				if (track->charge() > 0 && particle.GetQ() < 0) 
				{
					hCorrKplusO   ->Fill(kstar, 1./eff);
					hPtCorrKplusO ->Fill(dpt, 1./eff);
					hyCorrKplusO  ->Fill(dy, 1./eff);
					hphiCorrKplusO->Fill(dphi), 1./eff;

					hKpluspt_omega ->Fill(lv2_ori.Perp(), 1./eff);
					hKpluseta_omega->Fill(lv2_ori.Eta()), 1./eff;
					hKplusphi_omega->Fill(lv2_ori.Phi()), 1./eff;
					hKplusy_omega  ->Fill(lv2_ori.Rapidity(), 1./eff2D);
					// hCorrKplusO_y_pT  ->Fill(dpt, dy);
					// hCorrKplusO_y_phi ->Fill(dphi, dy);
					// hCorrKplusO_phi_pT->Fill(dpt, dphi);

					if (hasLambdabar)
					{
						hKpluspt_wlb_omega ->Fill(lv2_ori.Perp(), 1./eff);
						hKplusy_wlb_omega  ->Fill(lv2_ori.Rapidity(), 1./eff2D);
					}
					else
					{
						hKpluspt_wolb_omega ->Fill(lv2_ori.Perp(), 1./eff);
						hKplusy_wolb_omega  ->Fill(lv2_ori.Rapidity(), 1./eff2D);
					}
				}
				if (track->charge() > 0 && particle.GetQ() > 0) 
				{
					hCorrKplusObar   ->Fill(kstar, weight / eff);
					hPtCorrKplusObar ->Fill(dpt, weight / eff);
					hyCorrKplusObar  ->Fill(dy, weight / eff);
					hphiCorrKplusObar->Fill(dphi, weight/ eff);

					hKpluspt_omegabar ->Fill(lv2_ori.Perp(), 1./eff);
					hKpluseta_omegabar->Fill(lv2_ori.Eta(), 1./eff);
					hKplusphi_omegabar->Fill(lv2_ori.Phi(), 1./eff);
					hKplusy_omegabar  ->Fill(lv2_ori.Rapidity(), 1./eff2D);
					// hCorrKplusObar_y_pT  ->Fill(dpt, dy, weight);
					// hCorrKplusObar_y_phi ->Fill(dphi, dy, weight);
					// hCorrKplusObar_phi_pT->Fill(dpt, dphi, weight);

					if (hasLambdabar)
					{
						hKpluspt_wlb_omegabar ->Fill(lv2_ori.Perp(), 1./eff);
						hKplusy_wlb_omegabar  ->Fill(lv2_ori.Rapidity(), 1./eff2D);
					}
					else
					{
						hKpluspt_wolb_omegabar ->Fill(lv2_ori.Perp(), 1./eff);
						hKplusy_wolb_omegabar  ->Fill(lv2_ori.Rapidity(), 1./eff2D);
					}
				}
				if (track->charge() < 0 && particle.GetQ() < 0)
				{
					hCorrKminusO   ->Fill(kstar, 1./eff);
					hPtCorrKminusO ->Fill(dpt, 1./eff);
					hyCorrKminusO  ->Fill(dy, 1./eff);
					hphiCorrKminusO->Fill(dphi, 1./eff);
					if (dpt < 0.5) hNegPtDiff_dphi_KmO->Fill(dphi);
					if (dpt > 1.0) hPosPtDiff_dphi_KmO->Fill(dphi);
					
					hKminuspt_omega ->Fill(lv2_ori.Perp(), 1./eff);
					hKminuseta_omega->Fill(lv2_ori.Eta(), 1./eff);
					hKminusphi_omega->Fill(lv2_ori.Phi(), 1./eff);
					hKminusy_omega  ->Fill(lv2_ori.Rapidity(), 1./eff2D);
					// hCorrKminusO_y_pT  ->Fill(dpt, dy);
					// hCorrKminusO_y_phi ->Fill(dphi, dy);
					// hCorrKminusO_phi_pT->Fill(dpt, dphi);

					if (hasLambda)
					{
						hKminuspt_wl_omega ->Fill(lv2_ori.Perp(), 1./eff);
						hKminusy_wl_omega  ->Fill(lv2_ori.Rapidity(), 1./eff2D);
					}
					else
					{
						hKminuspt_wol_omega ->Fill(lv2_ori.Perp(), 1./eff);
						hKminusy_wol_omega  ->Fill(lv2_ori.Rapidity(), 1./eff2D);
					}
				}
				if (track->charge() < 0 && particle.GetQ() > 0) 
				{
					hCorrKminusObar   ->Fill(kstar, weight / eff);
					hPtCorrKminusObar ->Fill(dpt, weight / eff);
					hyCorrKminusObar  ->Fill(dy, weight / eff);
					hphiCorrKminusObar->Fill(dphi, weight / eff);
					if (dpt < 0.5) hNegPtDiff_dphi_KmOb->Fill(dphi, weight);
					if (dpt > 1.0) hPosPtDiff_dphi_KmOb->Fill(dphi, weight);

					hKminuspt_omegabar ->Fill(lv2_ori.Perp(), 1./eff);
					hKminuseta_omegabar->Fill(lv2_ori.Eta(), 1./eff);
					hKminusphi_omegabar->Fill(lv2_ori.Phi(), 1./eff);
					hKminusy_omegabar  ->Fill(lv2_ori.Rapidity(), 1./eff2D);
					// hCorrKminusObar_y_pT  ->Fill(dpt, dy, weight);
					// hCorrKminusObar_y_phi ->Fill(dphi, dy, weight);
					// hCorrKminusObar_phi_pT->Fill(dpt, dphi, weight);

					if (hasLambda)
					{
						hKminuspt_wl_omegabar ->Fill(lv2_ori.Perp(), 1./eff);
						hKminusy_wl_omegabar  ->Fill(lv2_ori.Rapidity(), 1./eff2D);
					}
					else
					{
						hKminuspt_wol_omegabar ->Fill(lv2_ori.Perp(), 1./eff);
						hKminusy_wol_omegabar  ->Fill(lv2_ori.Rapidity(), 1./eff2D);
					}
				}
			} // End loop over regular Omega
			fill_nomega = true;

			// Omega sideband loop
			for (int iOmega=0; iOmega < OmegaSidebandVec.size(); iOmega++)
			{ 
				const KFParticle particle = OmegaSidebandVec[iOmega]; 
				if (IsKaonOmegaDaughter(particle, kaonindex)) continue;

				// couting Omega
				if (!fill_nomega_sideband)
				{
					if (particle.GetPDG() > 0) hOmegaUsed_sideband->Fill(0.);
					if (particle.GetPDG() < 0) hOmegaUsed_sideband->Fill(1.);
					if ( hasLambdabar && particle.GetPDG() > 0) hOmegaUsed_wlb_sideband->Fill(0.);
					if ( hasLambdabar && particle.GetPDG() < 0) hOmegaUsed_wlb_sideband->Fill(1.);
					if (!hasLambdabar && particle.GetPDG() > 0) hOmegaUsed_wolb_sideband->Fill(0.);
					if (!hasLambdabar && particle.GetPDG() < 0) hOmegaUsed_wolb_sideband->Fill(1.);
					if ( hasLambda && particle.GetPDG() > 0) hOmegaUsed_wl_sideband->Fill(0.);
					if ( hasLambda && particle.GetPDG() < 0) hOmegaUsed_wl_sideband->Fill(1.);
					if (!hasLambda && particle.GetPDG() > 0) hOmegaUsed_wol_sideband->Fill(0.);
					if (!hasLambda && particle.GetPDG() < 0) hOmegaUsed_wol_sideband->Fill(1.);
				}

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
				TLorentzVector lv2_ori; lv2_ori.SetVectM(track->gMom(), KaonPdgMass);
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
				
				float eff = TPCEff * TOFEff;
				float eff2D = TPCEff2D * TOFEff2D;
				if (track->charge() > 0 && particle.GetQ() < 0) 
				{
					hCorrKplusO_sideband   ->Fill(kstar, 1./eff);
					hPtCorrKplusO_sideband ->Fill(dpt, 1./eff);
					hyCorrKplusO_sideband  ->Fill(dy, 1./eff);
					hphiCorrKplusO_sideband->Fill(dphi, 1./eff);

					hKpluspt_omega_sideband ->Fill(lv2_ori.Perp(), 1./eff);
					hKpluseta_omega_sideband->Fill(lv2_ori.Eta(), 1./eff);
					hKplusphi_omega_sideband->Fill(lv2_ori.Phi(), 1./eff);
					hKplusy_omega_sideband  ->Fill(lv2_ori.Rapidity(), 1./eff2D);
					//hCorrKplusO_y_pT  ->Fill(dpt, dy);
					//hCorrKplusO_y_phi ->Fill(dphi, dy);
					//hCorrKplusO_phi_pT->Fill(dpt, dphi);

					if (hasLambdabar)
					{
						hKpluspt_wlb_omega_sideband ->Fill(lv2_ori.Perp(), 1./eff);
						hKplusy_wlb_omega_sideband  ->Fill(lv2_ori.Rapidity(), 1./eff2D);
					}
					else
					{
						hKpluspt_wolb_omega_sideband ->Fill(lv2_ori.Perp(), 1./eff);
						hKplusy_wolb_omega_sideband  ->Fill(lv2_ori.Rapidity(), 1./eff2D);
					}
				}
				if (track->charge() > 0 && particle.GetQ() > 0) 
				{
					hCorrKplusObar_sideband   ->Fill(kstar, weight/eff);
					hPtCorrKplusObar_sideband ->Fill(dpt, weight/eff);
					hyCorrKplusObar_sideband  ->Fill(dy, weight/eff);
					hphiCorrKplusObar_sideband->Fill(dphi, weight/eff);

					hKpluspt_omegabar_sideband ->Fill(lv2_ori.Perp(), 1./eff);
					hKpluseta_omegabar_sideband->Fill(lv2_ori.Eta(), 1./eff);
					hKplusphi_omegabar_sideband->Fill(lv2_ori.Phi(), 1./eff);
					hKplusy_omegabar_sideband  ->Fill(lv2_ori.Rapidity(), 1./eff2D);
					//hCorrKplusObar_y_pT  ->Fill(dpt, dy);
					//hCorrKplusObar_y_phi ->Fill(dphi, dy);
					//hCorrKplusObar_phi_pT->Fill(dpt, dphi);

					if (hasLambdabar)
					{
						hKpluspt_wlb_omegabar_sideband ->Fill(lv2_ori.Perp(), 1./eff);
						hKplusy_wlb_omegabar_sideband  ->Fill(lv2_ori.Rapidity(), 1./eff2D);
					}
					else
					{
						hKpluspt_wolb_omegabar_sideband ->Fill(lv2_ori.Perp(), 1./eff);
						hKplusy_wolb_omegabar_sideband  ->Fill(lv2_ori.Rapidity(), 1./eff2D);
					}
				}
				if (track->charge() < 0 && particle.GetQ() < 0)
				{
					hCorrKminusO_sideband   ->Fill(kstar, 1./eff);
					hPtCorrKminusO_sideband ->Fill(dpt, 1./eff);
					hyCorrKminusO_sideband  ->Fill(dy, 1./eff);
					hphiCorrKminusO_sideband->Fill(dphi, 1./eff);

					hKminuspt_omega_sideband ->Fill(lv2_ori.Perp(), 1./eff);
					hKminuseta_omega_sideband->Fill(lv2_ori.Eta(), 1./eff);
					hKminusphi_omega_sideband->Fill(lv2_ori.Phi(), 1./eff);
					hKminusy_omega_sideband  ->Fill(lv2_ori.Rapidity(), 1./eff2D);
					//if (dpt < 0.5) hNegPtDiff_dphi_KmO->Fill(dphi);
					//if (dpt > 1.0) hPosPtDiff_dphi_KmO->Fill(dphi);

					//hCorrKminusO_y_pT  ->Fill(dpt, dy);
					//hCorrKminusO_y_phi ->Fill(dphi, dy);
					//hCorrKminusO_phi_pT->Fill(dpt, dphi);

					if (hasLambda)
					{
						hKminuspt_wl_omega_sideband ->Fill(lv2_ori.Perp(), 1./eff);
						hKminusy_wl_omega_sideband  ->Fill(lv2_ori.Rapidity(), 1./eff2D);
					}
					else
					{
						hKminuspt_wol_omega_sideband ->Fill(lv2_ori.Perp(), 1./eff);
						hKminusy_wol_omega_sideband  ->Fill(lv2_ori.Rapidity(), 1./eff2D);
					}

				}
				if (track->charge() < 0 && particle.GetQ() > 0) 
				{
					hCorrKminusObar_sideband   ->Fill(kstar, weight/eff);
					hPtCorrKminusObar_sideband ->Fill(dpt, weight/eff);
					hyCorrKminusObar_sideband  ->Fill(dy, weight/eff);
					hphiCorrKminusObar_sideband->Fill(dphi, weight/eff);

					hKminuspt_omegabar_sideband ->Fill(lv2_ori.Perp(), 1./eff);
					hKminuseta_omegabar_sideband->Fill(lv2_ori.Eta(), 1./eff);
					hKminusphi_omegabar_sideband->Fill(lv2_ori.Phi(), 1./eff);
					hKminusy_omegabar_sideband  ->Fill(lv2_ori.Rapidity(), 1./eff2D);
					//if (dpt < 0.5) hNegPtDiff_dphi_KmOb->Fill(dphi);
					//if (dpt > 1.0) hPosPtDiff_dphi_KmOb->Fill(dphi);

					//hCorrKminusObar_y_pT  ->Fill(dpt, dy);
					//hCorrKminusObar_y_phi ->Fill(dphi, dy);
					//hCorrKminusObar_phi_pT->Fill(dpt, dphi);

					if (hasLambda)
					{
						hKminuspt_wl_omegabar_sideband ->Fill(lv2_ori.Perp(), 1./eff);
						hKminusy_wl_omegabar_sideband  ->Fill(lv2_ori.Rapidity(), 1./eff2D);
					}
					else
					{
						hKminuspt_wol_omegabar_sideband ->Fill(lv2_ori.Perp(), 1./eff);
						hKminusy_wol_omegabar_sideband  ->Fill(lv2_ori.Rapidity(), 1./eff2D);
					}
				}
			} // End loop over sideband Omega
			fill_nomega_sideband = true;
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
	if ((Qx2e == 0. || Qx2w == 0. || Qy2e == 0. || Qy2w == 0.) && !PerformMixing) return kStOK;
	TVector2 Q2w, Q2e, Q2full; // for EP
	float EP2_TPC_w = 0, EP2_TPC_e = 0, EP2_TPC_full = 0;
	Q2full.Set(Qx2e+Qx2w, Qy2e+Qy2w);
	Q2e.Set(Qx2e, Qy2e);
	Q2w.Set(Qx2w, Qy2w);
	EP2_TPC_full = Q2full.Phi() / 2.;
	EP2_TPC_e = Q2e.Phi() / 2.;
	EP2_TPC_w = Q2w.Phi() / 2.;

	// shift correction
	float EP2_TPC[3]         = {EP2_TPC_e, EP2_TPC_w, EP2_TPC_full}; 
	float EP2_TPC_shifted[3] = {EP2_TPC_e, EP2_TPC_w, EP2_TPC_full};
	for (int ewFull = 0; ewFull < 3; ewFull++)
	{
		for (int i = 1; i <= shift_order_EP; i++)
		{
			hTPCEPShiftOutput_cos[ewFull]->Fill(dayPointer, i, cent-1, cos(i*2.*EP2_TPC[ewFull]), mWght);
			hTPCEPShiftOutput_sin[ewFull]->Fill(dayPointer, i, cent-1, sin(i*2.*EP2_TPC[ewFull]), mWght);
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
		
		EP2_TPC_shifted[ewFull] = ShiftTPCEP(EP2_TPC[ewFull], dayPointer+1, cent, ewFull, 2);
		float shift = EP2_TPC_shifted[ewFull] - EP2_TPC[ewFull];

		hTPCEP_2[cent-1][ewFull]->Fill(EP2_TPC[ewFull]);
		hTPCEP_2_shifted[cent-1][ewFull]->Fill(EP2_TPC_shifted[ewFull]);
		hTPCEP_2_2D[cent-1][ewFull]->Fill(dayPointer, EP2_TPC[ewFull]);
		hTPCEP_2_2D_shifted[cent-1][ewFull]->Fill(dayPointer, EP2_TPC_shifted[ewFull]);
		hTPCEP_2_shift[cent-1][ewFull]->Fill(EP2_TPC[ewFull], shift);
	}
	float EP2_TPC_e_shifted = EP2_TPC_shifted[0];
	float EP2_TPC_w_shifted = EP2_TPC_shifted[1];
	float EP2_TPC_full_shifted = EP2_TPC_shifted[2]; 
	hTPCEP_ew_cos->Fill(cent, cos(2.*EP2_TPC_e_shifted-2.*EP2_TPC_w_shifted), mWght);// resolution

	// EPD user shift
	float EP1_EPD[3]         = {result.EastPhiWeightedPsi(1), result.WestPhiWeightedPsi(1), result.FullPhiWeightedPsi(1)}; // EPD EP1
	float EP2_EPD[3]         = {result.EastPhiWeightedPsi(2), result.WestPhiWeightedPsi(2), result.FullPhiWeightedPsi(2)}; // EPD EP2
	float EP1_EPD_shifted[3] = {result.EastPhiWeightedPsi(1), result.WestPhiWeightedPsi(1), result.FullPhiWeightedPsi(1)}; // EPD EP1 initialize
	float EP2_EPD_shifted[3] = {result.EastPhiWeightedPsi(2), result.WestPhiWeightedPsi(2), result.FullPhiWeightedPsi(2)}; // EPD EP2 initialize
	for (int ewFull = 0; ewFull < 3; ewFull++)
		{	
			for (int i = 1; i <= shift_order_EP; i++)
			{
			hEPDEPShiftOutput_1_cos[ewFull]->Fill(dayPointer, i, cent-1, cos(i*1.*EP1_EPD[ewFull]), mWght);
			hEPDEPShiftOutput_1_sin[ewFull]->Fill(dayPointer, i, cent-1, sin(i*1.*EP1_EPD[ewFull]), mWght);
			hEPDEPShiftOutput_2_cos[ewFull]->Fill(dayPointer, i, cent-1, cos(i*2.*EP2_EPD[ewFull]), mWght);
			hEPDEPShiftOutput_2_sin[ewFull]->Fill(dayPointer, i, cent-1, sin(i*2.*EP2_EPD[ewFull]), mWght);
			}
		}
	// shift
	for (int ewFull = 0; ewFull < 3; ewFull++)
	{
		EP1_EPD_shifted[ewFull] = ShiftEPDEP(EP1_EPD[ewFull], dayPointer+1, cent, ewFull, 1);
		EP2_EPD_shifted[ewFull] = ShiftEPDEP(EP2_EPD[ewFull], dayPointer+1, cent, ewFull, 2);
		float shift = EP2_EPD_shifted[ewFull] - EP2_EPD[ewFull];
		
		hEPDEP_1[cent-1][ewFull]->Fill(EP1_EPD[ewFull]);
		hEPDEP_1_shifted[cent-1][ewFull]->Fill(EP1_EPD_shifted[ewFull]);
		hEPDEP_1_2D[cent-1][ewFull]->Fill(dayPointer, EP1_EPD[ewFull]);
		hEPDEP_1_2D_shifted[cent-1][ewFull]->Fill(dayPointer, EP1_EPD_shifted[ewFull]);
		hEPDEP_2[cent-1][ewFull]->Fill(EP2_EPD[ewFull]);
		hEPDEP_2_shifted[cent-1][ewFull]->Fill(EP2_EPD_shifted[ewFull]);
		hEPDEP_2_2D[cent-1][ewFull]->Fill(dayPointer, EP2_EPD[ewFull]);
		hEPDEP_2_2D_shifted[cent-1][ewFull]->Fill(dayPointer, EP2_EPD_shifted[ewFull]);
		hEPDEP_2_shift[cent-1][ewFull]->Fill(EP2_EPD[ewFull], shift);
	}
	float EP1_EPD_e_shifted = EP1_EPD_shifted[0];
	float EP1_EPD_w_shifted = EP1_EPD_shifted[1];
	float EP1_EPD_full_shifted = EP1_EPD_shifted[2];
	float EP2_EPD_e_shifted = EP2_EPD_shifted[0];
	float EP2_EPD_w_shifted = EP2_EPD_shifted[1];
	float EP2_EPD_full_shifted = EP2_EPD_shifted[2];
	hEPDEP_ew_cos_1->Fill(cent, cos(1.*EP1_EPD_shifted[0]-1.*EP1_EPD_shifted[1]), mWght);// resolution
	hEPDEP_ew_cos_2->Fill(cent, cos(2.*EP2_EPD_shifted[0]-2.*EP2_EPD_shifted[1]), mWght);// resolution
	// hEPDEP_ew_cos_1_for_v2->Fill(cent, cos(1.*EP1_EPD_shifted[0]-1.*EP1_EPD_shifted[1]), mWght);// EP resolution of 1st order plane for v_2

	// coalescence v2 and v1
	for (int i = 0; i < pion_tracks.size(); i++) 
	{
		StPicoTrack *track = mPicoDst->track(pion_tracks[i]);
		if (!track) continue;
		float pt = track->gMom().Perp();
		float eta = track->pMom().Eta();
		float phi = track->pMom().Phi();
		float phi_shifted_POI = ShiftPOIPhi(phi, dayPointer+1, cent);
		float y = Eta2y(pt, eta, PionPdgMass);

		float TOFEff = 1.0;
		float TPCEff = 1.0;
		if (hTOFEff[cent-1] != 0 && pt > pion_pT_TOFth) TOFEff = hTOFEff[cent-1]->GetEfficiency(hTOFEff[cent-1]->FindFixBin(pt));
		if (track->charge() > 0) 
		{	
			TPCEff = eff_finder.GetEfficiency1D(pt, cent, "pip");
			hTPCEff_check[cent-1]->Fill(pt, TPCEff);
			hpiplus_EPD_v2->Fill(cent, cos(2.*phi_shifted_POI-EP1_EPD_w_shifted-EP1_EPD_e_shifted), 1./TOFEff/TPCEff);
			hpiplus_EPD_v2_pt[cent-1]->Fill(pt, cos(2.*phi_shifted_POI-EP1_EPD_w_shifted-EP1_EPD_e_shifted));
			if (eta > 0) 
			{
				hpiplus_TPC_v2->Fill(cent, cos(2.*(phi_shifted_POI-EP2_TPC_w_shifted)), 1./TOFEff/TPCEff);
				// hpiplus_EPD_v2->Fill(cent, cos(2.*(phi_shifted_POI-EP2_EPD_w_shifted)), 1./TOFEff/TPCEff);
				hpiplus_TPC_v2_pt[cent-1]->Fill(pt, cos(2.*(phi_shifted_POI-EP2_TPC_w_shifted)));
				// hpiplus_EPD_v2_pt[cent-1]->Fill(pt, cos(2.*(phi_shifted_POI-EP2_EPD_w_shifted)));
				// v1 slope
				hpiplus_EPD_v1_y[cent-1]->Fill(y, cos(phi_shifted_POI-EP1_EPD_w_shifted), 1./TOFEff/TPCEff);
			}
			else        
			{
				hpiplus_TPC_v2->Fill(cent, cos(2.*(phi_shifted_POI-EP2_TPC_e_shifted)), 1./TOFEff/TPCEff);
				// hpiplus_EPD_v2->Fill(cent, cos(2.*(phi_shifted_POI-EP2_EPD_e_shifted)), 1./TOFEff/TPCEff);
				hpiplus_TPC_v2_pt[cent-1]->Fill(pt, cos(2.*(phi_shifted_POI-EP2_TPC_e_shifted)));
				// hpiplus_EPD_v2_pt[cent-1]->Fill(pt, cos(2.*(phi_shifted_POI-EP2_EPD_e_shifted)));
				// v1 slope
				hpiplus_EPD_v1_y[cent-1]->Fill(y, cos(phi_shifted_POI-EP1_EPD_e_shifted), 1./TOFEff/TPCEff);
			} 
		}
		else 					 
		{
			TPCEff = eff_finder.GetEfficiency1D(pt, cent, "pim");
			hpiminus_EPD_v2->Fill(cent, cos(2.*phi_shifted_POI-EP1_EPD_w_shifted-EP1_EPD_e_shifted), 1./TOFEff/TPCEff);
			hpiminus_EPD_v2_pt[cent-1]->Fill(pt, cos(2.*phi_shifted_POI-EP1_EPD_w_shifted-EP1_EPD_e_shifted));
			if (eta > 0)
			{
				hpiminus_TPC_v2->Fill(cent, cos(2.*(phi_shifted_POI-EP2_TPC_w_shifted)), 1./TOFEff/TPCEff);
				// hpiminus_EPD_v2->Fill(cent, cos(2.*(phi_shifted_POI-EP2_EPD_w_shifted)), 1./TOFEff/TPCEff);
				hpiminus_TPC_v2_pt[cent-1]->Fill(pt, cos(2.*(phi_shifted_POI-EP2_TPC_w_shifted)));
				// hpiminus_EPD_v2_pt[cent-1]->Fill(pt, cos(2.*(phi_shifted_POI-EP2_EPD_w_shifted)));
				// v1 slope
				hpiminus_EPD_v1_y[cent-1]->Fill(y, cos(phi_shifted_POI-EP1_EPD_w_shifted), 1./TOFEff/TPCEff);
			}
			else
			{
				hpiminus_TPC_v2->Fill(cent, cos(2.*(phi_shifted_POI-EP2_TPC_e_shifted)), 1./TOFEff/TPCEff);
				// hpiminus_EPD_v2->Fill(cent, cos(2.*(phi_shifted_POI-EP2_EPD_e_shifted)), 1./TOFEff/TPCEff);
				hpiminus_TPC_v2_pt[cent-1]->Fill(pt, cos(2.*(phi_shifted_POI-EP2_TPC_e_shifted)));
				// hpiminus_EPD_v2_pt[cent-1]->Fill(pt, cos(2.*(phi_shifted_POI-EP2_EPD_e_shifted)));
				// v1 slope
				hpiminus_EPD_v1_y[cent-1]->Fill(y, cos(phi_shifted_POI-EP1_EPD_e_shifted), 1./TOFEff/TPCEff);
			}
		}

	}
	for (int i = 0; i < proton_tracks.size(); i++) 
	{
		StPicoTrack *track = mPicoDst->track(proton_tracks[i]);
		if (!track) continue;
		float pt = track->gMom().Perp();
		float eta = track->pMom().Eta();
		float phi = track->pMom().Phi();
		float phi_shifted_POI = ShiftPOIPhi(phi, dayPointer+1, cent);
		float y = Eta2y(pt, eta, ProtonPdgMass);

		float TOFEff = 1.0;
		float TPCEff = 1.0; 
		if (hTOFEff[cent-1] != 0 && pt > proton_pT_TOFth) TOFEff = hTOFEff[cent-1]->GetEfficiency(hTOFEff[cent-1]->FindFixBin(pt));
		if (track->charge() > 0) 
		{
			TPCEff = eff_finder.GetEfficiency1D(pt, cent, "P");
			hproton_EPD_v2->Fill(cent, cos(2.*phi_shifted_POI-EP1_EPD_w_shifted-EP1_EPD_e_shifted), 1./TOFEff/TPCEff);
			hproton_EPD_v2_pt[cent-1]->Fill(pt, cos(2.*phi_shifted_POI-EP1_EPD_w_shifted-EP1_EPD_e_shifted));
			if (eta > 0)
			{
				hproton_TPC_v2->Fill(cent, cos(2.*(phi_shifted_POI-EP2_TPC_w_shifted)), 1./TOFEff/TPCEff);
				// hproton_EPD_v2->Fill(cent, cos(2.*(phi_shifted_POI-EP2_EPD_w_shifted)), 1./TOFEff/TPCEff);
				hproton_TPC_v2_pt[cent-1]->Fill(pt, cos(2.*(phi_shifted_POI-EP2_TPC_w_shifted)));
				// hproton_EPD_v2_pt[cent-1]->Fill(pt, cos(2.*(phi_shifted_POI-EP2_EPD_w_shifted)));
				// v1 slope
				hproton_EPD_v1_y[cent-1]->Fill(y, cos(phi_shifted_POI-EP1_EPD_w_shifted), 1./TOFEff/TPCEff);
			}
			else
			{
				hproton_TPC_v2->Fill(cent, cos(2.*(phi_shifted_POI-EP2_TPC_e_shifted)), 1./TOFEff/TPCEff);
				// hproton_EPD_v2->Fill(cent, cos(2.*(phi_shifted_POI-EP2_EPD_e_shifted)), 1./TOFEff/TPCEff);
				hproton_TPC_v2_pt[cent-1]->Fill(pt, cos(2.*(phi_shifted_POI-EP2_TPC_e_shifted)));
				// hproton_EPD_v2_pt[cent-1]->Fill(pt, cos(2.*(phi_shifted_POI-EP2_EPD_e_shifted)));
				// v1 slope
				hproton_EPD_v1_y[cent-1]->Fill(y, cos(phi_shifted_POI-EP1_EPD_e_shifted), 1./TOFEff/TPCEff);
			}
		}
		else 					 
		{
			TPCEff = eff_finder.GetEfficiency1D(pt, cent, "AP");
			hantiproton_EPD_v2->Fill(cent, cos(2.*phi_shifted_POI-EP1_EPD_w_shifted-EP1_EPD_e_shifted), 1./TOFEff/TPCEff);
			hantiproton_EPD_v2_pt[cent-1]->Fill(pt, cos(2.*phi_shifted_POI-EP1_EPD_w_shifted-EP1_EPD_e_shifted));
			if (eta > 0)
			{
				hantiproton_TPC_v2->Fill(cent, cos(2.*(phi_shifted_POI-EP2_TPC_w_shifted)), 1./TOFEff/TPCEff);
				// hantiproton_EPD_v2->Fill(cent, cos(2.*(phi_shifted_POI-EP2_EPD_w_shifted)), 1./TOFEff/TPCEff);
				hantiproton_TPC_v2_pt[cent-1]->Fill(pt, cos(2.*(phi_shifted_POI-EP2_TPC_w_shifted)));
				// hantiproton_EPD_v2_pt[cent-1]->Fill(pt, cos(2.*(phi_shifted_POI-EP2_EPD_w_shifted)));
				// v1 slope
				hantiproton_EPD_v1_y[cent-1]->Fill(y, cos(phi_shifted_POI-EP1_EPD_w_shifted), 1./TOFEff/TPCEff);
			}
			else
			{
				hantiproton_TPC_v2->Fill(cent, cos(2.*(phi_shifted_POI-EP2_TPC_e_shifted)), 1./TOFEff/TPCEff);
				// hantiproton_EPD_v2->Fill(cent, cos(2.*(phi_shifted_POI-EP2_EPD_e_shifted)), 1./TOFEff/TPCEff);
				hantiproton_TPC_v2_pt[cent-1]->Fill(pt, cos(2.*(phi_shifted_POI-EP2_TPC_e_shifted)));
				// hantiproton_EPD_v2_pt[cent-1]->Fill(pt, cos(2.*(phi_shifted_POI-EP2_EPD_e_shifted)));
				// v1 slope
				hantiproton_EPD_v1_y[cent-1]->Fill(y, cos(phi_shifted_POI-EP1_EPD_e_shifted), 1./TOFEff/TPCEff);
			}
		}
	}

	float kplus_totaly = 0, kminus_totaly = 0;
	float kplus_ntrack = 0, kminus_ntrack = 0;
	for (int i = 0; i < kaon_tracks.size(); i++) 
	{
		StPicoTrack *track = mPicoDst->track(kaon_tracks[i]);
		if (!track) continue;
		float pt = track->gMom().Perp();
		float eta = track->pMom().Eta();
		float phi = track->pMom().Phi();
		float phi_shifted_POI = ShiftPOIPhi(phi, dayPointer+1, cent);
		float y = Eta2y(pt, eta, KaonPdgMass);

		float TOFEff = 1.0, TOFEff_2D = 1.0;
		float TPCEff = 1.0; 
		if (hTOFEff[cent-1] != 0 && pt > 0.4) TOFEff = hTOFEff[cent-1]->GetEfficiency(hTOFEff[cent-1]->FindFixBin(pt));
		if (hTOFEff_2D[cent-1] != 0 && pt > 0.4) TOFEff_2D = hTOFEff_2D[cent-1]->GetEfficiency(hTOFEff_2D[cent-1]->FindFixBin(pt, eta));
		if (track->charge() > 0) 
		{	
			TPCEff = eff_finder.GetEfficiency1D(pt, cent, "Kp");
			if (eta > 0)
			{
				hkplus_TPC_v2->Fill(cent, cos(2.*(phi_shifted_POI-EP2_TPC_w_shifted)), 1./TOFEff/TPCEff);
				hkplus_EPD_v2->Fill(cent, cos(2.*(phi_shifted_POI-EP2_EPD_w_shifted)), 1./TOFEff/TPCEff);
				hkplus_TPC_v2_pt[cent-1]->Fill(pt, cos(2.*(phi_shifted_POI-EP2_TPC_w_shifted)));
				hkplus_EPD_v2_pt[cent-1]->Fill(pt, cos(2.*(phi_shifted_POI-EP2_EPD_w_shifted)));
				// v1 slope
				hkplus_EPD_v1_y[cent-1]->Fill(y, cos(phi_shifted_POI-EP1_EPD_w_shifted), 1./TOFEff/TPCEff);
			}
			else
			{
				hkplus_TPC_v2->Fill(cent, cos(2.*(phi_shifted_POI-EP2_TPC_e_shifted)), 1./TOFEff/TPCEff);
				hkplus_EPD_v2->Fill(cent, cos(2.*(phi_shifted_POI-EP2_EPD_e_shifted)), 1./TOFEff/TPCEff);
				hkplus_TPC_v2_pt[cent-1]->Fill(pt, cos(2.*(phi_shifted_POI-EP2_TPC_e_shifted)));
				hkplus_EPD_v2_pt[cent-1]->Fill(pt, cos(2.*(phi_shifted_POI-EP2_EPD_e_shifted)));
				// v1 slope
				hkplus_EPD_v1_y[cent-1]->Fill(y, cos(phi_shifted_POI-EP1_EPD_e_shifted), 1./TOFEff/TPCEff);
			}

			// mean y
			float TPCEff2D = eff_finder.GetEfficiency2D(pt, eta, cent, "Kp");
			kplus_totaly += Eta2y(pt, eta, KaonPdgMass) * 1.0 / TOFEff_2D / TPCEff2D;
			kplus_ntrack += 1.0 / TOFEff_2D / TPCEff2D;
		}
		else 					 
		{
			TPCEff = eff_finder.GetEfficiency1D(pt, cent, "Km");
			if (eta > 0)
			{
				hkminus_TPC_v2->Fill(cent, cos(2.*(phi_shifted_POI-EP2_TPC_w_shifted)), 1./TOFEff/TPCEff);
				hkminus_EPD_v2->Fill(cent, cos(2.*(phi_shifted_POI-EP2_EPD_w_shifted)), 1./TOFEff/TPCEff);
				hkminus_TPC_v2_pt[cent-1]->Fill(pt, cos(2.*(phi_shifted_POI-EP2_TPC_w_shifted)));
				hkminus_EPD_v2_pt[cent-1]->Fill(pt, cos(2.*(phi_shifted_POI-EP2_EPD_w_shifted)));
				// v1 slope
				hkminus_EPD_v1_y[cent-1]->Fill(y, cos(phi_shifted_POI-EP1_EPD_w_shifted), 1./TOFEff/TPCEff);
			}
			else
			{
				hkminus_TPC_v2->Fill(cent, cos(2.*(phi_shifted_POI-EP2_TPC_e_shifted)), 1./TOFEff/TPCEff);
				hkminus_EPD_v2->Fill(cent, cos(2.*(phi_shifted_POI-EP2_EPD_e_shifted)), 1./TOFEff/TPCEff);
				hkminus_TPC_v2_pt[cent-1]->Fill(pt, cos(2.*(phi_shifted_POI-EP2_TPC_e_shifted)));
				hkminus_EPD_v2_pt[cent-1]->Fill(pt, cos(2.*(phi_shifted_POI-EP2_EPD_e_shifted)));
				// v1 slope
				hkminus_EPD_v1_y[cent-1]->Fill(y, cos(phi_shifted_POI-EP1_EPD_e_shifted), 1./TOFEff/TPCEff);
			}

			// mean y
			float TPCEff2D = eff_finder.GetEfficiency2D(pt, eta, cent, "Km");
			kminus_totaly += Eta2y(pt, eta, KaonPdgMass) * 1.0 / TOFEff_2D / TPCEff2D;
			kminus_ntrack += 1.0 / TOFEff_2D / TPCEff2D;
		}
	}
	if (OmegaVec.size() == 1 && OmegaVec[0].GetQ() < 0) // Omega event
	{
		hkplus_totaly->Fill(0., kplus_totaly);
		hkplus_ntrack->Fill(0., kplus_ntrack);
		hkminus_totaly->Fill(1., kminus_totaly);
		hkminus_ntrack->Fill(1., kminus_ntrack);
	}  
	if (OmegaVec.size() == 1 && OmegaVec[0].GetQ() > 0) // Omegabar event
	{
		hkplus_totaly->Fill(1., kplus_totaly);
		hkplus_ntrack->Fill(1., kplus_ntrack);
		hkminus_totaly->Fill(0., kminus_totaly);
		hkminus_ntrack->Fill(0., kminus_ntrack);
	}

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
		float EP2_TPC_Omega = particle.GetEta()>0? EP2_TPC_w_shifted: EP2_TPC_e_shifted;
		float EP2_EPD_Omega = particle.GetEta()>0? EP2_EPD_w_shifted: EP2_EPD_e_shifted;
		if (particle.GetQ() < 0) 
		{
			hOmega_TPC_v2[cent-1]->Fill(particle.GetMass(), cos(2.*particle.GetPhi() - 2.*EP2_TPC_Omega)); 
			hOmega_EPD_v2[cent-1]->Fill(particle.GetMass(), cos(2.*particle.GetPhi() - 2.*EP2_EPD_Omega));
		}
		else if (particle.GetQ() > 0)
		{
			hOmegabar_TPC_v2[cent-1]->Fill(particle.GetMass(), cos(2.*particle.GetPhi() - 2.*EP2_TPC_Omega)); 
			hOmegabar_EPD_v2[cent-1]->Fill(particle.GetMass(), cos(2.*particle.GetPhi() - 2.*EP2_EPD_Omega));
		}
	}

	// Xi v_2
	for (int i = 0; i < XiVecAll.size(); i++)
	{
		const KFParticle particle = XiVecAll[i];
		float EP2_TPC_Xi = particle.GetEta()>0? EP2_TPC_w_shifted: EP2_TPC_e_shifted;
		float EP2_EPD_Xi = particle.GetEta()>0? EP2_EPD_w_shifted: EP2_EPD_e_shifted;
		if (particle.GetQ() < 0) 
		{
			hXi_TPC_v2[cent-1]->Fill(particle.GetMass(), cos(2.*particle.GetPhi() - 2.*EP2_TPC_Xi)); 
			hXi_EPD_v2[cent-1]->Fill(particle.GetMass(), cos(2.*particle.GetPhi() - 2.*EP2_EPD_Xi));
		}
		else if (particle.GetQ() > 0)
		{
			hXibar_TPC_v2[cent-1]->Fill(particle.GetMass(), cos(2.*particle.GetPhi() - 2.*EP2_TPC_Xi)); 
			hXibar_EPD_v2[cent-1]->Fill(particle.GetMass(), cos(2.*particle.GetPhi() - 2.*EP2_EPD_Xi));
		}
	}

	// phi v_2
	for (int i = 0; i < phiVecAll.size(); i++)
	{
		const KFParticle particle = phiVecAll[i];
		float EP2_TPC_phi = particle.GetEta()>0? EP2_TPC_w_shifted: EP2_TPC_e_shifted;
		float EP2_EPD_phi = particle.GetEta()>0? EP2_EPD_w_shifted: EP2_EPD_e_shifted;
		hphi_TPC_v2[cent-1]->Fill(particle.GetMass(), cos(2.*particle.GetPhi() - 2.*EP2_TPC_phi));
		hphi_EPD_v2[cent-1]->Fill(particle.GetMass(), cos(2.*particle.GetPhi() - 2.*EP2_EPD_phi));
	}
		

	// Lambda v2
	// for (int i = 0; i < LambdaVecAll.size(); i++)
	// {
	// 	const KFParticle particle = LambdaVecAll[i];
	// 	if (particle.GetPDG() > 0) 
	// 	{	
	// 		if (particle.GetPt() >= 4.0) continue;
	// 		int ptbin = static_cast<int>(floor(particle.GetPt()/0.1));
	// 		float EP2_TPC_Lambda = particle.GetEta()>0? EP2_TPC_w_shifted: EP2_TPC_e_shifted;
	// 		hLambda_TPC_v2_pt[cent-1][ptbin]->Fill(particle.GetMass(), cos(2.*particle.GetPhi() - 2.*EP2_TPC_Lambda)); 
	// 		hLambda_EPD_v2_pt[cent-1][ptbin]->Fill(particle.GetMass(), cos(2.*particle.GetPhi() - 2.*EP2_EPD_full));
	// 	}
	// 	else if (particle.GetPDG() < 0)
	// 	{
	// 		if (particle.GetPt() >= 4.0) continue;
	// 		int ptbin = static_cast<int>(floor(particle.GetPt()/0.1));
	// 		float EP2_TPC_Lambda = particle.GetEta()>0? EP2_TPC_w_shifted: EP2_TPC_e_shifted;
	// 		hLambdabar_TPC_v2_pt[cent-1][ptbin]->Fill(particle.GetMass(), cos(2.*particle.GetPhi() - 2.*EP2_TPC_Lambda)); 
	// 		hLambdabar_EPD_v2_pt[cent-1][ptbin]->Fill(particle.GetMass(), cos(2.*particle.GetPhi() - 2.*EP2_EPD_full));
	// 	}
	// }

	// new observable
	if (OmegaVec.size() == 1 && OmegaVec[0].GetQ() < 0) hKratio_omega   ->Fill(pbct*1.0/pct, kmct*1.0/kpct);
	if (OmegaVec.size() == 0)                           hKratio_wo      ->Fill(pbct*1.0/pct, kmct*1.0/kpct);
	if (OmegaVec.size() == 1 && OmegaVec[0].GetQ() > 0) hKratio_omegabar->Fill(pbct*1.0/pct, kmct*1.0/kpct);

	// filling trees
	if (StoringTree && !OmegaVec.size() == 0) omega_mix[mult_index][vz_index][EP_index]->Fill();
	// if (StoringTree && !OmegaVec.size() == 0) omega_mix[0][0][0]->Fill(); // for testing

	// counting kaon, need to be right before event mixing
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
				if (!track) continue;
				
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

bool StKFParticleAnalysisMaker::isGoodLambda(int cent, KFParticle Lambda)
{
	if (Lambda.GetPDG() ==  3122) return (fabs(Lambda.GetMass() - LambdaPdgMass) < 3*LambdaMassSigma   [cent-1]);
	if (Lambda.GetPDG() == -3122) return (fabs(Lambda.GetMass() - LambdaPdgMass) < 3*LambdabarMassSigma[cent-1]);
	return false;
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

float StKFParticleAnalysisMaker::ShiftPOIPhi(float phi, int day, int cent)
{
	float phi_shifted = phi;
	if (hTPCPOIShiftInput_cos != 0) // attempt to read input
	{
		for (int i=1; i <= shift_order_POI; i++)
		{
			float cosAve = hTPCPOIShiftInput_cos->GetBinContent(day, i, cent);
			float sinAve = hTPCPOIShiftInput_sin->GetBinContent(day, i, cent);
			phi_shifted += 2.0/i*(cosAve*sin(i*1.0*phi) - sinAve*cos(i*1.0*phi));
		}
	}
	return phi_shifted;
}

float StKFParticleAnalysisMaker::ShiftAssoPhi(float phi, int day, int cent)
{
	float phi_shifted = phi;
	if (hTPCAssoShiftInput_cos != 0) // attempt to read input
	{
		for (int i=1; i <= shift_order_asso; i++)
		{
			float cosAve = hTPCAssoShiftInput_cos->GetBinContent(day, i, cent);
			float sinAve = hTPCAssoShiftInput_sin->GetBinContent(day, i, cent);
			phi_shifted += 2.0/i*(cosAve*sin(i*1.0*phi) - sinAve*cos(i*1.0*phi));
		}
	}
	return phi_shifted;
}

float StKFParticleAnalysisMaker::ShiftTPCEP(float psi, int day, int cent, int ewFull, int order)
{
	float psi_shifted = psi;
	// if (order == 2) 
	if (hTPCEPShiftInput_cos[ewFull] != 0) // attempt to read input
	{
		for (int i=1; i <= shift_order_EP; i++)
		{
			float cosAve = hTPCEPShiftInput_cos[ewFull]->GetBinContent(day, i, cent);
			float sinAve = hTPCEPShiftInput_sin[ewFull]->GetBinContent(day, i, cent);
			psi_shifted += 2.0/order/i*(cosAve*sin(i*order*psi) - sinAve*cos(i*order*psi));
		}
	}
	return psi_shifted;
}

float StKFParticleAnalysisMaker::ShiftEPDEP(float psi, int day, int cent, int ewFull, int order)
{
	float psi_shifted = psi;
	if (order == 1 && hEPDEPShiftInput_1_cos[ewFull] != 0) // attempt to read input
	{
		for (int i=1; i <= shift_order_EP; i++)
		{
			float cosAve = hEPDEPShiftInput_1_cos[ewFull]->GetBinContent(day, i, cent);
			float sinAve = hEPDEPShiftInput_1_sin[ewFull]->GetBinContent(day, i, cent);
			psi_shifted += 2.0/order/i*(cosAve*sin(i*order*psi) - sinAve*cos(i*order*psi));
		}
	}
	if (order == 2 && hEPDEPShiftInput_2_cos[ewFull] != 0) // attempt to read input
	{
		for (int i=1; i <= shift_order_EP; i++)
		{
			float cosAve = hEPDEPShiftInput_2_cos[ewFull]->GetBinContent(day, i, cent);
			float sinAve = hEPDEPShiftInput_2_sin[ewFull]->GetBinContent(day, i, cent);
			psi_shifted += 2.0/order/i*(cosAve*sin(i*order*psi) - sinAve*cos(i*order*psi));
		}
	}
	return psi_shifted;
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

// rapidity as a function of pseudorapidity
float StKFParticleAnalysisMaker::Eta2y(float pt, float eta, float mass)
{
	float p = pt*cosh(eta);
	float pl = pt*sinh(eta);
	float energy = sqrt(p*p + mass*mass);
	float rapidity = 0.5*log((energy + pl)/(energy - pl));
	return rapidity;
}