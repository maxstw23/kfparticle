#include "StKFParticleAnalysisMaker.h"
#include "PhysicalConstants.h"
#include "phys_constants.h"
#include "StBTofUtil/tofPathLength.hh"
#include "StThreeVectorF.hh"
#include "TVector2.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH2D.h"
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
#include <vector>
#include <set>

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

#define pi TMath::Pi()
#define OmegaPdgMass 1.67245
// #define OmegaMassSigma     0.0021
#define XiPdgMass 1.32171
#define LambdaPdgMass 1.11568
#define ProtonPdgMass 0.938272
#define PionPdgMass 0.139570
#define KaonPdgMass 0.493677
#define LambdaPdg 3122
#define XiPdg 3312
#define phiPdg 333
#define OmegaPdg 3334
#define KaonPdg 321
#define ProtonPdg 2212
#define KsPdg 310
#define Xi1530Pdg 3324
#define PionPdg -211
// #define cen_cut            9
// #define num_pt_bin 10 // previous for Omega
#define num_phi_bin 10
#define num_mult_bin 7
#define num_vz_bin 14
#define num_EP_bin 10
#define num_pt_bin 14 // for Lambda, 0.4 - 1.8
#define shift_order_asso 24
#define shift_order_POI 24
#define shift_order_EP 5
#define TPCAssoEtaCut 0.5
#define USE_P 0
#define num_RecoPar 7
#define num_IdPar 6
#define EPD_part_cut 3.4 // participant is below this value
#define EPD_spec_cut 3.8 // spectator is above this value

#define OmegaMassSigma 0.00308
#define OmegabarMassSigma 0.00314
#define LambdaMassSigma 0.00254
#define LambdabarMassSigma 0.00245

// const float StKFParticleAnalysisMaker::OmegaMassSigma[]    = {0.00272, 0.00314, 0.00269, 0.00311, 0.00286, 0.00319, 0.00296}; // {0.00219, 0.00228, 0.00218, 0.00262, 0.00247, 0.00243, 0.00280};
// const float StKFParticleAnalysisMaker::OmegabarMassSigma[] = {0.00326, 0.00335, 0.00302, 0.00315, 0.00319, 0.00298, 0.00292}; // {0.00238, 0.00214, 0.00233, 0.00250, 0.00257, 0.00270, 0.00269};
// const float StKFParticleAnalysisMaker::OmegaMassPtLowerBin[] = {0.5, 0.9, 1.2, 1.5, 1.8, 2.1, 2.4, 2.7, 3.0, 3.4, 4.0};
// const float StKFParticleAnalysisMaker::LambdaMassSigma[]    = {0.00205, 0.00206, 0.00206, 0.00206, 0.00216, 0.00217, 0.00218, 0.00219, 0.00220};
// const float StKFParticleAnalysisMaker::LambdabarMassSigma[] = {0.00202, 0.00203, 0.00203, 0.00204, 0.00204, 0.00205, 0.00207, 0.00207, 0.00208};

//-----------------------------------------------------------------------------
ClassImp(StKFParticleAnalysisMaker)

	//-----------------------------------------------------------------------------
	StKFParticleAnalysisMaker::StKFParticleAnalysisMaker(const char *name, const char *outName, int sys)
	: StMaker(name)
{
	mOutName = outName;
	sys_tag = sys;
	mRun = 0;
	mEnergy = 0;
	mListDir = "./";
}

//-----------------------------------------------------------------------------
StKFParticleAnalysisMaker::~StKFParticleAnalysisMaker()
{
	SafeDelete(KFParticleInterface);
	SafeDelete(KFParticlePerformanceInterface);
}

//-----------------------------------------------------------------------------
Int_t StKFParticleAnalysisMaker::openFile()
{
	// ======= StRefMultCorr ======= //
	// set up StRefMultCorr for centrality
	mRefMultCorr = CentralityMaker::instance()->getRefMultCorr();
	// ======= StRefMultCorr end ======= //

	// EPD
	/* Eta Weight */
	TH2D wt("Order1etaWeight", "Order1etaWeight", 100, 1.5, 6.5, 9, 0, 9);
	TH2D wt2("Order2etaWeight", "Order2etaWeight", 100, 1.5, 6.5, 9, 0, 9);
	float lin[9] = {-1.950, -1.900, -1.850, -1.706, -1.438, -1.340, -1.045, -0.717, -0.700};
	float cub[9] = {0.1608, 0.1600, 0.1600, 0.1595, 0.1457, 0.1369, 0.1092, 0.0772, 0.0700};
	float par1[9] = {474.651, 474.651, 474.651, 474.651, 474.651, 3.27243e+02, 1.72351, 1.72351, 1.72351};
	float par2[9] = {3.55515, 3.55515, 3.55515, 3.55515, 3.55515, 3.56938, -7.69075, -7.69075, -7.69075};
	float par3[9] = {1.80162, 1.80162, 1.80162, 1.80162, 1.80162, 1.67113, 4.72770, 4.72770, 4.72770};
	for (int iy = 1; iy <= 9; iy++)
	{
		for (int ix = 1; ix < 101; ix++)
		{
			double eta = wt.GetXaxis()->GetBinCenter(ix);
			if (UseParticipant) wt.SetBinContent(ix, iy, (fabs(eta) < EPD_part_cut) ? lin[iy - 1] * eta + cub[iy - 1] * pow(eta, 3) : 0);
			else wt.SetBinContent(ix, iy, (fabs(eta) > EPD_spec_cut) ? lin[iy - 1] * eta + cub[iy - 1] * pow(eta, 3) : 0);
			if (UseParticipant) wt2.SetBinContent(ix, iy, (fabs(eta) > EPD_spec_cut) ? par1[iy - 1] + par2[iy - 1] * pow(eta, 2) + par3[iy - 1] * pow(eta, 4) : 0);
			else wt2.SetBinContent(ix, iy, (fabs(eta) < EPD_part_cut) ? par1[iy - 1] + par2[iy - 1] * pow(eta, 2) + par3[iy - 1] * pow(eta, 4) : 0);
		}
	}

	/* Set up StEpdEpFinder */
	char fname_in[200];
	char fname_out[200];
	if (sys_tag != 1) 
	{
		if (UseParticipant) sprintf(fname_in, "./weight/default/cent_EPD_CorrectionInput_pp.root");
		else sprintf(fname_in, "./weight/default/cent_EPD_CorrectionInput.root"); 
	}
	else 
	{
		if (UseParticipant) sprintf(fname_in, "./weight/default/cent_EPD_CorrectionInput_pp.root");
		else sprintf(fname_in, "./weight/default/cent_EPD_CorrectionInput.root");
	}
	sprintf(fname_out, "cent_EPD_CorrectionOutput_%d.root", mJob);
	// mEpdHits = new TClonesArray("StPicoEpdHit");
	// unsigned int found;
	// chain->SetBranchStatus("EpdHit*",1,&found);
	// cout << "EpdHit Branch returned found= " << found << endl;
	// chain->SetBranchAddress("EpdHit",&mEpdHits);
	mEpFinder = new StEpdEpFinder(9, fname_out, fname_in);
	mEpFinder->SetnMipThreshold(0.3); // recommended by EPD group
	mEpFinder->SetMaxTileWeight(1.0); // recommended by EPD group, 1.0 for low multiplicity (BES)
	mEpFinder->SetEpdHitFormat(2);	  // 2=pico
	mEpFinder->SetEtaWeights(1, wt);  // eta weight for 1st-order EP
	mEpFinder->SetEtaWeights(2,wt2);	// eta weight for 2nd-order EP

	/* TPC weight files */
	if (sys_tag != 1) sprintf(fname_in, "./weight/default/TPCShiftInput.root");
	else sprintf(fname_in, "./weight/default/TPCShiftInput.root");
	fTPCShift = new TFile(fname_in, "READ");
	if (fTPCShift->IsZombie())
	{
		std::cout << "\n**********************************************" << std::endl;
		std::cout << "No TPC shift input! No shift correction used. " << std::endl;
		std::cout << "**********************************************" << std::endl;
		hTPCAssoShiftInput_cos = 0;
		hTPCAssoShiftInput_sin = 0;
		hTPCPOIShiftInput_cos = 0;
		hTPCPOIShiftInput_sin = 0;
		for (int ewFull = 0; ewFull < 3; ewFull++)
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
		hTPCAssoShiftInput_cos = (TProfile3D *)fTPCShift->Get(Form("TPCAssoShift_cos"));
		hTPCAssoShiftInput_sin = (TProfile3D *)fTPCShift->Get(Form("TPCAssoShift_sin"));
		hTPCPOIShiftInput_cos = (TProfile3D *)fTPCShift->Get(Form("TPCPOIShift_cos"));
		hTPCPOIShiftInput_sin = (TProfile3D *)fTPCShift->Get(Form("TPCPOIShift_sin"));

		for (int ewFull = 0; ewFull < 3; ewFull++)
		{
			hTPCEPShiftInput_cos[ewFull] = (TProfile3D *)fTPCShift->Get(Form("TPCEPShiftEW%d_cos", ewFull));
			hTPCEPShiftInput_sin[ewFull] = (TProfile3D *)fTPCShift->Get(Form("TPCEPShiftEW%d_sin", ewFull));
		}
	}

	/* EPD weight files */
	if (sys_tag != 1) 
	{
		if (UseParticipant) sprintf(fname_in, "./weight/default/EPDShiftInput_pp.root");
		else sprintf(fname_in, "./weight/default/EPDShiftInput.root");
	}
	else 
	{
		if (UseParticipant) sprintf(fname_in, "./weight/default/EPDShiftInput_pp.root");
		else sprintf(fname_in, "./weight/default/EPDShiftInput.root");
	}
	fEPDShift = new TFile(fname_in, "READ");
	if (fEPDShift->IsZombie())
	{
		std::cout << "\n**********************************************" << std::endl;
		std::cout << "No EPD shift input! No user shift correction used. " << std::endl;
		std::cout << "**********************************************" << std::endl;
		for (int ewFull = 0; ewFull < 3; ewFull++)
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
		for (int ewFull = 0; ewFull < 3; ewFull++)
		{
			hEPDEPShiftInput_1_cos[ewFull] = (TProfile3D *)fEPDShift->Get(Form("EPDEPShiftEW%d_1_cos", ewFull));
			hEPDEPShiftInput_1_sin[ewFull] = (TProfile3D *)fEPDShift->Get(Form("EPDEPShiftEW%d_1_sin", ewFull));
			hEPDEPShiftInput_2_cos[ewFull] = (TProfile3D *)fEPDShift->Get(Form("EPDEPShiftEW%d_2_cos", ewFull));
			hEPDEPShiftInput_2_sin[ewFull] = (TProfile3D *)fEPDShift->Get(Form("EPDEPShiftEW%d_2_sin", ewFull));
		}
	}

	// obselete, but no need to change
	std::cout << "Before reading TOF" << std::endl;
	/* TOF Efficiency */
	fTOFEff = new TFile("TOFEfficiency.root", "READ");
	std::cout << "After reading TOF" << std::endl;
	if (fTOFEff->IsZombie())
	{
		std::cout << "\n**********************************************" << std::endl;
		std::cout << "No TOF efficiency input! No TOF efficiency correction used. " << std::endl;
		std::cout << "**********************************************" << std::endl;
		for (int i = 0; i < 9; i++)
			hTOFEff[i] = 0;
		for (int i = 0; i < 9; i++)
			hTOFEff_2D[i] = 0;
	}
	else
	{
		std::cout << "\n**********************************************" << std::endl;
		std::cout << "TOF efficiency found. " << std::endl;
		std::cout << "**********************************************" << std::endl;
		for (int i = 0; i < 9; i++)
			hTOFEff[i] = (TEfficiency *)fTOFEff->Get(Form("hTOFEff_%d", i + 1));
		for (int i = 0; i < 9; i++)
			hTOFEff_2D[i] = (TH2D *)fTOFEff->Get(Form("hTOFEff_2D_%d", i + 1));
	}
	if (hTOFEff[0] == 0)
		std::cout << "hTOFEff[0] is null" << std::endl;
	if (hTOFEff_2D[0] == 0)
		std::cout << "hTOFEff_2D[0] is null" << std::endl;

	cout << "-----------------------------------" << endl;
	cout << "------- user's files loaded -------" << endl;
	cout << "-----------------------------------" << endl;
	return kStOK;
}

//-----------------------------------------------------------------------------
Int_t StKFParticleAnalysisMaker::Init()
{

	cout << mOutName << endl;
	const char *sDir = mOutName.Data();
	int lDir = mOutName.Length();
	int iDir = lDir - 6;
	while (!std::isdigit(sDir[iDir]))
		iDir--;
	while (std::isdigit(sDir[iDir]))
		iDir--;
	mJob = std::atoi(&sDir[iDir + 1]);
	cout << "current job id: " << mJob << endl;
	mRunStart = 50000;
	mRunEnd = 100000;
	mDay = 50;

	PI = M_PI;
	twoPI = 2 * M_PI;
	buffer.Init();

	badList.clear();
	runList.clear();

	PerformAnalysis = true;
	PerformMixing = false;
	CheckWeights2D = false;
	UseParticipant = false; // 1st order plane. Also means spectator 2nd
	StoringTree = false;
	CutCent = true; // for omega, should always be true
	PtReweighting = false;
	v2Calculation = true;
	Coal2D = false;
	UsePtForPID = false;
	// v2EPDMethod = "2nd";
	// USE_P = false; // 0 for p, 1 for pt. Note: this applies only for all the v2_vs_pT plots,
	// 				  // not for the inclusive v2 plots (which are always vs pT).
	// 				  // Also, only protons and pions use this flag.

	if (!readRunList())
		return kStFatal;
	// if(!readBadList())return kStFatal;

	dcatoPV_hi = 3.0;
	// pid boundary
	pT_lo = 0.2;
	pT_hi = 2.0;

	// for EP
	pT_asso_lo = 0.15;
	pT_asso_hi = 2.0;
	pT_trig_lo = 0.15;
	pT_trig_hi = 2.0;
	eta_trig_cut = 1.0;

	// coalescence study cut
	y_coal_cut = 0.6; // TOF is flat up to 0.9

	// phi cut
	y_phi_cut = 0.5; // mid-rapidity phi

	// pion cut
	pion_pT_lo = 0.16;
	pion_pT_hi = 1.2;
	pion_pT_TOFth = 0.4;
	pion_m2_lo = -0.01;
	pion_m2_hi = 0.1;

	// proton cut
	proton_pT_lo = 0.24;
	proton_pT_hi = 1.8;
	proton_pT_TOFth = 0.6;
	proton_m2_lo = 0.75;
	proton_m2_hi = 1.1;

	// kaon cut
	kaon_pT_TOFth = 0.4;
	kaon_m2_lo = 0.15;
	kaon_m2_hi = 0.34;

	// correlation centrality
	min_cent = 8;
	max_cent = 9;

	DeclareHistograms();
	if (StoringTree && PerformMixing)
		return kStFatal; // cannot perform mixing and storing mixing tree at the same time
	if (StoringTree)
		DeclareTrees();
	if (PerformMixing)
		ReadTrees();
	Int_t openFileStatus = openFile();
	if (openFileStatus == kStFatal)
		return kStFatal;

	TFile *f = GetTFile(); // These two lines need to be HERE (though I don't know /why/)- don't throw in another function
	if (f)
	{
		f->cd();
		BookVertexPlots();
	}

	return kStOK;
}

//-----------------------------------------------------------------------------
Int_t StKFParticleAnalysisMaker::Finish()
{
	if (mOutName != "")
	{
		TFile *fout = new TFile(mOutName.Data(), "RECREATE");
		fout->cd();
		WriteHistograms();
		fout->Close();
	}

	if (StoringTree)
	{
		char temp[200];
		if (CutCent)
			sprintf(temp, "omega_mix_cen%d%d_%d.root", min_cent, max_cent, mJob);
		else
			sprintf(temp, "omega_mix_%d.root", mJob);
		ftree = new TFile(temp, "RECREATE");
		ftree->cd();
		WriteTrees();
		ftree->Close();
	}

	if (PerformMixing)
		ftree->Close();

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
void StKFParticleAnalysisMaker::DeclareHistograms()
{
	TString RecoPar_names[num_RecoPar] = {"Omega", "Omegabar", "Xi", "Xibar", "Lambda", "Lambdabar", "phi"};
	TString IdPar_names[num_IdPar] = {"piplus", "piminus", "proton", "antiproton", "kplus", "kminus"};
	float RecoPar_invmass_lo[num_RecoPar] = {1., 1., 1., 1., 1., 1., 0.7};
	float RecoPar_invmass_hi[num_RecoPar] = {2.4, 2.4, 1.6, 1.6, 1.2, 1.2, 1.3};

	const int nRuns = runList.size();
	int mRunL = runList.at(0);
	int mRunH = runList.at(nRuns - 1);
	int mDayL = (int)((mRunL) / 1000 % 1000);
	int mDayH = (int)((mRunH) / 1000 % 1000) + 1;
	const int nDays = mDayH - mDayL;
	int nTmps = nDays; //>40?40:nDays;
	mStps = ceil((nDays * 1.0) / (nTmps * 1.0));
	const int nTims = ceil((nDays * 1.0) / (mStps * 1.0));

	cout << nDays << " " << mDayL << " " << mDayH << " " << nTims << " " << mStps << endl;
	char temp[200];

	hNRefMult = new TH1F("RefMult", "Reference Multiplicity", 1000, 0.0, 1000.0);
	hNRefMultA = new TH1F("RefMultA", "Reference MultiplicityA", 1000, 0.0, 1000.0);
	hNRefMultB = new TH1F("RefMultB", "Reference MultiplicityB", 1000, 0.0, 1000.0);
	hVertexXY = new TH2F("VertexXY", "Vertex XY Position", 200, -10.0, 10.0, 200, -10., 10);
	hVertexZ = new TH1F("VertexZ", "Event Vertex Z Position", 300, -150.0, 150.0);
	hVertex2D = new TH2F("Vertex2D", "VertexZ vs VPD Vz", 200, -100.0, 100.0, 200, -100., 100);
	hDiffVz = new TH1F("VertexZdiff", "VertexZ-VPDVz diff", 100, -10.0, 10.0);
	hcent = new TH1F("centrality", "centrality", nCent, 0., nCent);
	hcentw = new TH1F("centralityw", "centralityw", nCent, 0., nCent);

	hcentRefM = new TProfile("hcentRefM", "hcentRefM", nCent, 0., nCent, 0, 1000);
	hcentRefW = new TProfile("hcentRefW", "hcentRefW", nCent, 0., nCent, 0, 1000);
	for (int i = 0; i < 9; i++)
	{
		sprintf(temp, "hRefMultCorr_cent_%d", i + 1);
		hRefMultCorr_cent[i] = new TH1D(temp, temp, 1000, -0.5, 999.5);
	}
	hMultRatioEW_VertexZ = new TProfile("hMultRatioEW_VertexZ", "hMultRatioEW_VertexZ", 300, -150., 150., 0., 500.);

	// xiatong's global track QA
	hEventQA = new TH1F("hEventQA", "hEventQA", 7, -0.5, 6.5);
	hgpdEdx = new TH2D("hgpdEdx", "Global dE/dx vs p", 1000, 0., 10., 1000, 0., 10.);
	hgdEdxErr = new TH1D("hgdEdxErr", "Global dE/dx error", 100, 0., 1.);
	hgpinvbeta = new TH2D("hgpinvbeta", "Global inverse beta vs p", 1000, 0., 10., 1000, 0., 10.);
	hgm2 = new TH1D("hgm2", "Global m^2", 500, -1., 4.);
	hgpm2 = new TH2D("hgpm2", "Global m^2 vs p", 1000, 0., 10., 500, -1., 4.);
	hgm2nSigmaKaon = new TH2D("hgm2nSigmaKaon", "Global m2 vs nSigmaKaon", 2000, -10, 10, 1200, -1., 5.);
	hgm2nSigmaPion = new TH2D("hgm2nSigmaPion", "Global m2 vs nSigmaPion", 2000, -10, 10, 1200, -1., 5.);
	hgm2nSigmaProton = new TH2D("hgm2nSigmaProton", "Global m2 vs nSigmaProton", 2000, -10, 10, 1200, -1., 5.);
	hgptnSigmaKaon = new TH2D("hgptnSigmaKaon", "Pt vs nSigmaKaon", 1000, 0., 10., 2000, -10., 10.);
	hgptnSigmaPion = new TH2D("hgptnSigmaPion", "Pt vs nSigmaPion", 1000, 0., 10., 2000, -10., 10.);
	hgptnSigmaProton = new TH2D("hgptnSigmaProton", "Pt vs nSigmaProton", 1000, 0., 10., 2000, -10., 10.);
	hgptm2 = new TH2D("hgptm2", "Pt vs m2", 1000, 0., 10., 500, -1., 4.);
	// hgp          = new TH1D("hgp", "Global momentum", 1000, 0., 10.);
	for (int i = 0; i < 9; i++)
	{
		hgp[i] = new TH1D(Form("hgp_%d", i + 1), "Global momentum", 500, 0., 5.);
		hgp_TOF[i] = new TH1D(Form("hgp_TOF_%d", i + 1), "Global momentum with TOF match", 500, 0., 5.);
		hgpT[i] = new TH1D(Form("hgpT_%d", i + 1), "Global transverse momentum", 500, 0., 5.);
		hgpT_TOF[i] = new TH1D(Form("hgpT_TOF_%d", i + 1), "Global transverse momentum with TOF match", 500, 0., 5.);
		hgpTeta[i] = new TH2F(Form("hgpTeta_%d", i + 1), "Global transverse momentum vs eta", 500, 0., 5., 400, -2., 2.);
		hgpTeta_TOF[i] = new TH2F(Form("hgpTeta_TOF_%d", i + 1), "Global transverse momentum vs eta with TOF match", 500, 0., 5., 400, -2., 2.);
	}
	hgDCAtoPV = new TH1D("hgDCAtoPV", "Global DCA to PV", 500, 0., 10.);
	hgbtofYlocal = new TH1D("hgbtofYlocal", "Global K+ BTOF Ylocal", 1000, -5., 5.);
	hDauProtonFirstPoint_lam_pt = new TProfile("hDauProtonFirstPoint_lam_pt", "hDauProtonFirstPoint_lam_pt", 100, 0., 10., 0., 500.);
	hDauProtonLastPoint_lam_pt = new TProfile("hDauProtonLastPoint_lam_pt", "hDauProtonLastPoint_lam_pt", 100, 0., 10., 0., 500.);
	hDauPionFirstPoint_lam_pt = new TProfile("hDauPionFirstPoint_lam_pt", "hDauPionFirstPoint_lam_pt", 100, 0., 10., 0., 500.);
	hDauPionLastPoint_lam_pt = new TProfile("hDauPionLastPoint_lam_pt", "hDauPionLastPoint_lam_pt", 100, 0., 10., 0., 500.);

	// y-pT/nq coverage
	for (int i = 0; i < num_RecoPar; i++)
		hRecoPar_y_ptnq[i] = new TH2F(Form("h%s_y_ptnq", RecoPar_names[i].Data()), Form("h%s_y_ptnq", RecoPar_names[i].Data()), 1000, -5., 5., 1000, 0., 10.);
	for (int i = 0; i < num_IdPar; i++)
		hIdPar_y_ptnq[i] = new TH2F(Form("h%s_y_ptnq", IdPar_names[i].Data()), Form("h%s_y_ptnq", IdPar_names[i].Data()), 1000, -5., 5., 1000, 0., 10.);

	// v2 and EP
	hTPCAssoPhi = new TH1D("hTPCAssoPhi", "hTPCAssoPhi", 1000, -PI, PI);
	hTPCAssoPhi_shifted = new TH1D("hTPCAssoPhi_shifted", "hTPCAssoPhi_shifted", 1000, -PI, PI);
	hTPCPOIPhi = new TH1D("hTPCPOIPhi", "hTPCPOIPhi", 1000, -PI, PI);
	hTPCPOIPhi_shifted = new TH1D("hTPCPOIPhi_shifted", "hTPCPOIPhi_shifted", 1000, -PI, PI);
	hTPCAssoPhi_2D = new TH2D("hTPCAssoPhi_2D", "hTPCAssoPhi_2D", mDay, -0.5, mDay - 0.5, 1000, -PI, PI);
	hTPCAssoPhi_2D_shifted = new TH2D("hTPCAssoPhi_2D_shifted", "hTPCAssoPhi_2D_shifted", mDay, -0.5, mDay - 0.5, 1000, -PI, PI);
	hTPCPOIPhi_2D = new TH2D("hTPCPOIPhi_2D", "hTPCPOIPhi_2D", mDay, -0.5, mDay - 0.5, 1000, -PI, PI);
	hTPCPOIPhi_2D_shifted = new TH2D("hTPCPOIPhi_2D_shifted", "hTPCPOIPhi_2D_shifted", mDay, -0.5, mDay - 0.5, 1000, -PI, PI);
	// TPC weights/shifts
	hTPCAssoShiftOutput_cos = new TProfile3D(Form("TPCAssoShift_cos"), Form("TPCAssoShift_cos"),
											 mDay, -0.5, mDay - 0.5, shift_order_asso, 0.5, 1.0 * shift_order_asso + 0.5, 9, -0.5, 8.5);
	hTPCAssoShiftOutput_sin = new TProfile3D(Form("TPCAssoShift_sin"), Form("TPCAssoShift_sin"),
											 mDay, -0.5, mDay - 0.5, shift_order_asso, 0.5, 1.0 * shift_order_asso + 0.5, 9, -0.5, 8.5);
	hTPCPOIShiftOutput_cos = new TProfile3D(Form("TPCPOIShift_cos"), Form("TPCPOIShift_cos"),
											mDay, -0.5, mDay - 0.5, shift_order_POI, 0.5, 1.0 * shift_order_POI + 0.5, 9, -0.5, 8.5);
	hTPCPOIShiftOutput_sin = new TProfile3D(Form("TPCPOIShift_sin"), Form("TPCPOIShift_sin"),
											mDay, -0.5, mDay - 0.5, shift_order_POI, 0.5, 1.0 * shift_order_POI + 0.5, 9, -0.5, 8.5);
	for (int ewFull = 0; ewFull < 3; ewFull++)
	{
		hTPCEPShiftOutput_cos[ewFull] = new TProfile3D(Form("TPCEPShiftEW%d_cos", ewFull), Form("TPCEPShiftEW%d_cos", ewFull),
													   mDay, -0.5, mDay - 0.5, shift_order_EP, 0.5, 1.0 * shift_order_EP + 0.5, 9, -0.5, 8.5);
		hTPCEPShiftOutput_sin[ewFull] = new TProfile3D(Form("TPCEPShiftEW%d_sin", ewFull), Form("TPCEPShiftEW%d_sin", ewFull),
													   mDay, -0.5, mDay - 0.5, shift_order_EP, 0.5, 1.0 * shift_order_EP + 0.5, 9, -0.5, 8.5);
		for (int i = 0; i < 9; i++)
		{
			hTPCEP_2[i][ewFull] = new TH1D(Form("hTPCEPEW%d_2_%d", ewFull, i + 1), Form("hTPCEPEW%d_2_%d", ewFull, i + 1), 1000, 0., 2 * PI);
			hTPCEP_2_shifted[i][ewFull] = new TH1D(Form("hTPCEPEW%d_2_shifted_%d", ewFull, i + 1), Form("hTPCEPEW%d_2_shifted_%d", ewFull, i + 1), 1000, 0., 2 * PI);
			if (CheckWeights2D)
			{
				hTPCEP_2_2D[i][ewFull] = new TH2D(Form("hTPCEPEW%d_2_2D_%d", ewFull, i + 1), Form("hTPCEPEW%d_2_2D_%d", ewFull, i + 1), mDay, -0.5, mDay - 0.5, 1000, 0., 2 * PI);
				hTPCEP_2_2D_shifted[i][ewFull] = new TH2D(Form("hTPCEPEW%d_2_2D_shifted_%d", ewFull, i + 1), Form("hTPCEPEW%d_2_2D_shifted_%d", ewFull, i + 1), mDay, -0.5, mDay - 0.5, 1000, 0., 2 * PI);
			}
			hTPCEP_2_shift[i][ewFull] = new TProfile(Form("hTPCEPEW%d_2_shift_%d", ewFull, i + 1), Form("hTPCEPEW%d_2_shift_%d", ewFull, i + 1), 1000, 0., 2 * PI, -PI, PI);
		}
	}
	// EPD user shifts
	for (int ewFull = 0; ewFull < 3; ewFull++)
	{
		hEPDEPShiftOutput_1_cos[ewFull] = new TProfile3D(Form("EPDEPShiftEW%d_1_cos", ewFull), Form("EPDEPShiftEW%d_1_cos", ewFull),
														 mDay, -0.5, mDay - 0.5, shift_order_EP, 0.5, 1.0 * shift_order_EP + 0.5, 9, -0.5, 8.5);
		hEPDEPShiftOutput_1_sin[ewFull] = new TProfile3D(Form("EPDEPShiftEW%d_1_sin", ewFull), Form("EPDEPShiftEW%d_1_sin", ewFull),
														 mDay, -0.5, mDay - 0.5, shift_order_EP, 0.5, 1.0 * shift_order_EP + 0.5, 9, -0.5, 8.5);
		hEPDEPShiftOutput_2_cos[ewFull] = new TProfile3D(Form("EPDEPShiftEW%d_2_cos", ewFull), Form("EPDEPShiftEW%d_2_cos", ewFull),
														 mDay, -0.5, mDay - 0.5, shift_order_EP, 0.5, 1.0 * shift_order_EP + 0.5, 9, -0.5, 8.5);
		hEPDEPShiftOutput_2_sin[ewFull] = new TProfile3D(Form("EPDEPShiftEW%d_2_sin", ewFull), Form("EPDEPShiftEW%d_2_sin", ewFull),
														 mDay, -0.5, mDay - 0.5, shift_order_EP, 0.5, 1.0 * shift_order_EP + 0.5, 9, -0.5, 8.5);
		for (int i = 0; i < 9; i++)
		{
			hEPDEP_1[i][ewFull] = new TH1D(Form("hEPDEPW%d_1_%d", ewFull, i + 1), Form("hEPDEPW%d_1_%d", ewFull, i + 1), 1000, 0., 2 * PI);
			hEPDEP_1_shifted[i][ewFull] = new TH1D(Form("hEPDEPW%d_1_shifted_%d", ewFull, i + 1), Form("hEPDEPW%d_1_shifted_%d", ewFull, i + 1), 1000, 0., 2 * PI);
			hEPDEP_2[i][ewFull] = new TH1D(Form("hEPDEPW%d_2_%d", ewFull, i + 1), Form("hEPDEPW%d_2_%d", ewFull, i + 1), 1000, 0., 2 * PI);
			hEPDEP_2_shifted[i][ewFull] = new TH1D(Form("hEPDEPW%d_2_shifted_%d", ewFull, i + 1), Form("hEPDEPW%d_2_shifted_%d", ewFull, i + 1), 1000, 0., 2 * PI);
			if (CheckWeights2D)
			{
				hEPDEP_1_2D[i][ewFull] = new TH2D(Form("hEPDEPW%d_1_2D_%d", ewFull, i + 1), Form("hEPDEPW%d_1_2D_%d", ewFull, i + 1), mDay, -0.5, mDay - 0.5, 1000, 0., 2 * PI);
				hEPDEP_1_2D_shifted[i][ewFull] = new TH2D(Form("hEPDEPW%d_1_2D_shifted_%d", ewFull, i + 1), Form("hEPDEPW%d_1_2D_shifted_%d", ewFull, i + 1), mDay, -0.5, mDay - 0.5, 1000, 0., 2 * PI);
				hEPDEP_2_2D[i][ewFull] = new TH2D(Form("hEPDEPW%d_2_2D_%d", ewFull, i + 1), Form("hEPDEPW%d_2_2D_%d", ewFull, i + 1), mDay, -0.5, mDay - 0.5, 1000, 0., 2 * PI);
				hEPDEP_2_2D_shifted[i][ewFull] = new TH2D(Form("hEPDEPW%d_2_2D_shifted_%d", ewFull, i + 1), Form("hEPDEPW%d_2_2D_shifted_%d", ewFull, i + 1), mDay, -0.5, mDay - 0.5, 1000, 0., 2 * PI);
			}
			hEPDEP_2_shift[i][ewFull] = new TProfile(Form("hEPDEPW%d_2_shift_%d", ewFull, i + 1), Form("hEPDEPW%d_2_shift_%d", ewFull, i + 1), 1000, 0., 2 * PI, -PI, PI);
		}
	}

	hTPCEP_ew_cos = new TProfile("hTPCEP_ew_cos", "hTPCEP_ew_cos", 9, 0.5, 9.5, -1., 1.);
	hTPCEP_ew_cos->Sumw2();
	hEPDEP_ew_cos_1 = new TProfile("hEPDEP_ew_cos_1", "hEPDEP_ew_cos_1", 9, 0.5, 9.5, -1., 1.);
	hEPDEP_ew_cos_1->Sumw2();
	hEPDEP_ew_cos_2 = new TProfile("hEPDEP_ew_cos_2", "hEPDEP_ew_cos_2", 9, 0.5, 9.5, -1., 1.);
	hEPDEP_ew_cos_2->Sumw2();
	hEPDEP_ew_cos_3 = new TProfile("hEPDEP_ew_cos_3", "hEPDEP_ew_cos_3", 9, 0.5, 9.5, -1., 1.);
	hEPDEP_ew_cos_3->Sumw2();

	// Reco particles flow
	for (int i = 0; i < num_RecoPar; i++)
	{
		int nBins = static_cast<int>((RecoPar_invmass_hi[i] - RecoPar_invmass_lo[i]) * 1000);
		for (int j = 0; j < 9; j++)
		{
			hRecoParM_cen[i][j] = new TH1D(Form("h%sM_cen_%d", RecoPar_names[i].Data(), j + 1), Form("h%sM_cen_%d", RecoPar_names[i].Data(), j + 1), nBins, RecoPar_invmass_lo[i], RecoPar_invmass_hi[i]);
			hRecoPar_EPD_v2[i][j] = new TProfile(Form("h%s_EPD_v2_%d", RecoPar_names[i].Data(), j + 1), Form("h%s_EPD_v2_%d", RecoPar_names[i].Data(), j + 1), nBins, RecoPar_invmass_lo[i], RecoPar_invmass_hi[i], -1., 1.);
			hRecoPar_EPD_v2[i][j]->Sumw2();
			hRecoPar_TPC_v2[i][j] = new TProfile(Form("h%s_TPC_v2_%d", RecoPar_names[i].Data(), j + 1), Form("h%s_TPC_v2_%d", RecoPar_names[i].Data(), j + 1), nBins, RecoPar_invmass_lo[i], RecoPar_invmass_hi[i], -1., 1.);
			hRecoPar_TPC_v2[i][j]->Sumw2();
		}
	}
	for (int i = 0; i < num_pt_bin; i++)
	{
		int nBins = static_cast<int>((RecoPar_invmass_hi[0] - RecoPar_invmass_lo[0]) * 1000);
		for (int j = 0; j < 9; j++)
		{
			hLambdaM_pt[i][j] = new TH1D(Form("hLambdaM_pt_%d_cen_%d", i + 1, j + 1), Form("hLambdaM_pt_%d_cen_%d", i + 1, j + 1), nBins, RecoPar_invmass_lo[4], RecoPar_invmass_hi[4]);
			hLambdabarM_pt[i][j] = new TH1D(Form("hLambdabarM_pt_%d_cen_%d", i + 1, j + 1), Form("hLambdabarM_pt_%d_cen_%d", i + 1, j + 1), nBins, RecoPar_invmass_lo[5], RecoPar_invmass_hi[5]);
			hLambda_EPD_v2_pt[i][j] = new TProfile(Form("hLambda_EPD_v2_pt_%d_cen_%d", i + 1, j + 1), Form("hLambda_EPD_v2_pt_%d_cen_%d", i + 1, j + 1), nBins, RecoPar_invmass_lo[4], RecoPar_invmass_hi[4], -1., 1.);
			hLambda_EPD_v2_pt[i][j]->Sumw2();
			hLambdabar_EPD_v2_pt[i][j] = new TProfile(Form("hLambdabar_EPD_v2_pt_%d_cen_%d", i + 1, j + 1), Form("hLambdabar_EPD_v2_pt_%d_cen_%d", i + 1, j + 1), nBins, RecoPar_invmass_lo[5], RecoPar_invmass_hi[5], -1., 1.);
			hLambdabar_EPD_v2_pt[i][j]->Sumw2();
			hLambda_TPC_v2_pt[i][j] = new TProfile(Form("hLambda_TPC_v2_pt_%d_cen_%d", i + 1, j + 1), Form("hLambda_TPC_v2_pt_%d_cen_%d", i + 1, j + 1), nBins, RecoPar_invmass_lo[4], RecoPar_invmass_hi[4], -1., 1.);
			hLambda_TPC_v2_pt[i][j]->Sumw2();
			hLambdabar_TPC_v2_pt[i][j] = new TProfile(Form("hLambdabar_TPC_v2_pt_%d_cen_%d", i + 1, j + 1), Form("hLambdabar_TPC_v2_pt_%d_cen_%d", i + 1, j + 1), nBins, RecoPar_invmass_lo[5], RecoPar_invmass_hi[5], -1., 1.);
			hLambdabar_TPC_v2_pt[i][j]->Sumw2();
		}
	}
	// Identified particles flow
	for (int i = 0; i < num_IdPar; i++)
	{
		hIdPar_EPD_v2[i] = new TProfile(Form("h%s_EPD_v2", IdPar_names[i].Data()), Form("h%s_EPD_v2", IdPar_names[i].Data()), 9, 0.5, 9.5, -1., 1.);
		hIdPar_EPD_v2[i]->Sumw2();
		hIdPar_TPC_v2[i] = new TProfile(Form("h%s_TPC_v2", IdPar_names[i].Data()), Form("h%s_TPC_v2", IdPar_names[i].Data()), 9, 0.5, 9.5, -1., 1.);
		hIdPar_TPC_v2[i]->Sumw2();
		for (int j = 0; j < 9; j++)
		{
			hIdPar_EPD_v2_pt_1st[i][j] = new TProfile(Form("h%s_EPD_v2_pt_%d_1st", IdPar_names[i].Data(), j + 1), Form("h%s_EPD_v2_pt_%d_1st", IdPar_names[i].Data(), j + 1), 1000, 0., 10., -1., 1.);
			hIdPar_EPD_v2_pt_1st[i][j]->Sumw2();
			hIdPar_EPD_v2_pt_2nd[i][j] = new TProfile(Form("h%s_EPD_v2_pt_%d_2nd", IdPar_names[i].Data(), j + 1), Form("h%s_EPD_v2_pt_%d_2nd", IdPar_names[i].Data(), j + 1), 1000, 0., 10., -1., 1.);
			hIdPar_EPD_v2_pt_2nd[i][j]->Sumw2();
			// hIdPar_EPD_v1_pt_1st[i][j] = new TProfile(Form("h%s_EPD_v1_pt_%d_1st", IdPar_names[i].Data(), j + 1), Form("h%s_EPD_v1_pt_%d_1st", IdPar_names[i].Data(), j + 1), 1000, 0., 10., -1., 1.);
			// hIdPar_EPD_v1_pt_1st[i][j]->Sumw2();
			// hIdPar_EPD_v3_pt_1st[i][j] = new TProfile(Form("h%s_EPD_v3_pt_%d_1st", IdPar_names[i].Data(), j + 1), Form("h%s_EPD_v3_pt_%d_1st", IdPar_names[i].Data(), j + 1), 1000, 0., 10., -1., 1.);
			// hIdPar_EPD_v3_pt_1st[i][j]->Sumw2();
			hIdPar_TPC_v2_pt[i][j] = new TProfile(Form("h%s_TPC_v2_pt_%d", IdPar_names[i].Data(), j + 1), Form("h%s_TPC_v2_pt_%d", IdPar_names[i].Data(), j + 1), 1000, 0., 10., -1., 1.);
			hIdPar_TPC_v2_pt[i][j]->Sumw2();

			if (Coal2D)
			{
				hIdPar_EPD_v1_eta_pt_1st[i][j] = new TProfile2D(Form("h%s_EPD_v1_eta_pt_%d_1st", IdPar_names[i].Data(), j + 1), Form("h%s_EPD_v1_eta_pt_%d_1st", IdPar_names[i].Data(), j + 1), 300, -1.5, 1.5, 500, 0., 5., -1., 1.);
				hIdPar_EPD_v1_eta_pt_1st[i][j]->Sumw2();
				hIdPar_EPD_v2_eta_pt_1st[i][j] = new TProfile2D(Form("h%s_EPD_v2_eta_pt_%d_1st", IdPar_names[i].Data(), j + 1), Form("h%s_EPD_v2_eta_pt_%d_1st", IdPar_names[i].Data(), j + 1), 300, -1.5, 1.5, 500, 0., 5., -1., 1.);
				hIdPar_EPD_v2_eta_pt_1st[i][j]->Sumw2();
				hIdPar_EPD_v3_eta_pt_1st[i][j] = new TProfile2D(Form("h%s_EPD_v3_eta_pt_%d_1st", IdPar_names[i].Data(), j + 1), Form("h%s_EPD_v3_eta_pt_%d_1st", IdPar_names[i].Data(), j + 1), 300, -1.5, 1.5, 500, 0., 5., -1., 1.);
				hIdPar_EPD_v3_eta_pt_1st[i][j]->Sumw2();
			}

			hIdPar_EPD_v1_y[i][j] = new TProfile(Form("h%s_EPD_v1_y_%d", IdPar_names[i].Data(), j + 1), Form("h%s_EPD_v1_y_%d", IdPar_names[i].Data(), j + 1), 10, -1., 1., -1., 1.);
			hIdPar_EPD_v1_y[i][j]->Sumw2();
		}
	}
	hyLambdaDauProton_east = new TH1D("hyLambdaDauProton_east", "y of Lambda daughter proton in east", 300, -1.5, 1.5);
	hyLambdaDauPion_east = new TH1D("hyLambdaDauPion_east", "y of Lambda daughter pion in east", 300, -1.5, 1.5);
	hyLambdaDauProton_west = new TH1D("hyLambdaDauProton_west", "y of Lambda daughter proton in west", 300, -1.5, 1.5);
	hyLambdaDauPion_west = new TH1D("hyLambdaDauPion_west", "y of Lambda daughter pion in west", 300, -1.5, 1.5);

	hEPD_e_EP_1 = new TH1D("hEPD_e_EP_1", "hEPD_e_EP_1", 1000, 0., 2 * PI);
	hEPD_w_EP_1 = new TH1D("hEPD_w_EP_1", "hEPD_w_EP_1", 1000, 0., 2 * PI);
	hEPD_e_EP_2 = new TH1D("hEPD_e_EP_2", "hEPD_e_EP_2", 1000, 0., 2 * PI);
	hEPD_w_EP_2 = new TH1D("hEPD_w_EP_2", "hEPD_w_EP_2", 1000, 0., 2 * PI);
	hEPD_full_EP_1 = new TH1D("hEPD_full_EP_1", "hEPD_full_EP_1", 1000, 0., 2 * PI);
	hEPD_full_EP_2 = new TH1D("hEPD_full_EP_2", "hEPD_full_EP_2", 1000, 0., 2 * PI);
	hEPD_ew_cos = new TProfile("hEPD_ew_cos", "hEPD_ew_cos", 3, 0.5, 3.5, -1., 1.);
	hEPD_ew_cos->Sumw2();
	hEPDEP_2_Full_Raw_QA_cen4 = new TH1D("hEPDEP_2_Full_Raw_QA_cen4", "hEPDEP_2_Full_Raw_QA_cen4", 1000, 0., 2 * PI);
	hEPDEP_2_Full_PhiWeighted_QA_cen4 = new TH1D("hEPDEP_2_Full_PhiWeighted_QA_cen4", "hEPDEP_2_Full_PhiWeighted_QA_cen4", 1000, 0., 2 * PI);
	hEPDEP_2_Full_PhiWeightedAndShifted_QA_cen4 = new TH1D("hEPDEP_2_Full_PhiWeightedAndShifted_QA_cen4", "hEPDEP_2_Full_PhiWeightedAndShifted_QA_cen4", 1000, 0., 2 * PI);

#if !USE_P
	// 2D pid
	// for (int i = 0; i < 10; i++)
	// {
	// 	hgPID2D_proton_pt[i] = new TH2F(Form("hgPID2D_proton_pt_%d", i), Form("hgPID2D_proton_pt_%d", i), 2000, -10, 10, 200, -0.5, 1.5);
	// 	hgPID2D_antiproton_pt[i] = new TH2F(Form("hgPID2D_antiproton_pt_%d", i), Form("hgPID2D_antiproton_pt_%d", i), 2000, -10, 10, 200, -0.5, 1.5);
	// 	hgPID2D_piplus_pt[i]   = new TH2F(Form("hgPID2D_piplus_pt_%d", i), Form("hgPID2D_piplus_pt_%d", i), 2000, -10, 10, 200, -0.5, 1.5);
	// 	hgPID2D_piminus_pt[i]  = new TH2F(Form("hgPID2D_piminus_pt_%d", i), Form("hgPID2D_piminus_pt_%d", i), 2000, -10, 10, 200, -0.5, 1.5);
	// 	hgPID2D_kplus_pt[i]   = new TH2F(Form("hgPID2D_kplus_pt_%d", i), Form("hgPID2D_kplus_pt_%d", i), 2000, -10, 10, 200, -0.5, 1.5);
	// 	hgPID2D_kminus_pt[i]  = new TH2F(Form("hgPID2D_kminus_pt_%d", i), Form("hgPID2D_kminus_pt_%d", i), 2000, -10, 10, 200, -0.5, 1.5);
	// }

	// for (int i = 0; i < 12; i++)
	// {
	// 	hgPID2D_proton_p[i] = new TH2F(Form("hgPID2D_proton_p_%d", i), Form("hgPID2D_proton_p_%d", i), 2000, -10, 10, 400, -0.5, 1.5);
	// 	hgPID2D_antiproton_p[i] = new TH2F(Form("hgPID2D_antiproton_p_%d", i), Form("hgPID2D_antiproton_p_%d", i), 2000, -10, 10, 400, -0.5, 1.5);
	// 	hgPID2D_piplus_p[i] = new TH2F(Form("hgPID2D_piplus_p_%d", i), Form("hgPID2D_piplus_p_%d", i), 2000, -10, 10, 400, -0.5, 1.5);
	// 	hgPID2D_piminus_p[i] = new TH2F(Form("hgPID2D_piminus_p_%d", i), Form("hgPID2D_piminus_p_%d", i), 2000, -10, 10, 400, -0.5, 1.5);
	// 	hgPID2D_kplus_p[i] = new TH2F(Form("hgPID2D_kplus_p_%d", i), Form("hgPID2D_kplus_p_%d", i), 2000, -10, 10, 400, -0.5, 1.5);
	// 	hgPID2D_kminus_p[i] = new TH2F(Form("hgPID2D_kminus_p_%d", i), Form("hgPID2D_kminus_p_%d", i), 2000, -10, 10, 400, -0.5, 1.5);
	// }
#else
	for (int i = 0; i < 12; i++)
	{
		hgPID2D_proton_p[i] = new TH2F(Form("hgPID2D_proton_p_%d", i), Form("hgPID2D_proton_p_%d", i), 2000, -10, 10, 400, -0.5, 1.5);
		hgPID2D_antiproton_p[i] = new TH2F(Form("hgPID2D_antiproton_p_%d", i), Form("hgPID2D_antiproton_p_%d", i), 2000, -10, 10, 400, -0.5, 1.5);
		hgPID2D_piplus_p[i] = new TH2F(Form("hgPID2D_piplus_p_%d", i), Form("hgPID2D_piplus_p_%d", i), 2000, -10, 10, 400, -0.5, 1.5);
		hgPID2D_piminus_p[i] = new TH2F(Form("hgPID2D_piminus_p_%d", i), Form("hgPID2D_piminus_p_%d", i), 2000, -10, 10, 400, -0.5, 1.5);
		hgPID2D_kplus_p[i] = new TH2F(Form("hgPID2D_kplus_p_%d", i), Form("hgPID2D_kplus_p_%d", i), 2000, -10, 10, 400, -0.5, 1.5);
		hgPID2D_kminus_p[i] = new TH2F(Form("hgPID2D_kminus_p_%d", i), Form("hgPID2D_kminus_p_%d", i), 2000, -10, 10, 400, -0.5, 1.5);
	}
#endif

	// testing eTOF
	hgetofcrossingY = new TH1D("hgetofcrossingY", "hgetofcrossingY", 1000, -50., 50.);
	hgm2_etof = new TH1D("hgm2_etof", "hgm2_etof", 500, -1., 4.);
	hgpinvbeta_etof = new TH2F("hgpinvbeta_etof", "hgpinvbeta_etof", 1000, 0., 10., 1000, 0., 10.);
	hgetofmatchflag = new TH1D("hgetofmatchflag", "hgetofmatchflag", 3, -0.5, 2.5);
	hgeta_etof = new TH1D("hgeta_etof", "hgeta_etof", 1000, -5., 5.);

	// kaon QA
	hgKpdEdx = new TH2D("hgKpdEdx", "Global K dE/dx vs p", 1000, 0., 10., 1000, 0., 10.);
	hgKpinvbeta = new TH2D("hgKpinvbeta", "Global K inverse beta vs p", 1000, 0., 10., 1000, 0., 10.);
	hgKm2 = new TH1D("hgKm2", "Global K m^2", 500, -1., 4.);
	hgKpm2 = new TH2D("hgKpm2", "Global K m^2 vs p", 1000, 0., 10., 500, -1., 4.);
	hgKp = new TH1D("hgKp", "Global K momentum", 1000, 0., 10.);
	hgKpT = new TH1D("hgKpT", "Global K transver momentum", 1000, 0., 10.);
	hgKDCAtoPV = new TH1D("hgKDCAtoPV", "Global K DCA to PV", 500, 0., 10.);
	hgKDCAtoO = new TH1D("hgKDCAtoO", "Global K DCA to Omega", 500, 0., 10.);
	hgKpionpdEdx = new TH2D("hgKpionpdEdx", "Misidentified kaon dEdx", 1000, 0., 10., 1000, 0., 10.);
	// hgptm2_largenSigmaKaon = new TH2D("hgptm2_largenSigmaKaon", "hgptm2_largenSigmaKaon", 2000, 0., 10., 1000, -1., 4.); // m2 vs pt for nSigmaKaon > 6
	// hgptm2_smallnSigmaKaon = new TH2D("hgptm2_smallnSigmaKaon", "hgptm2_smallnSigmaKaon", 2000, 0., 10., 1000, -1., 4.); // m2 vs pt for nSigmaKaon < -6
	hgnSigmaDiff = new TH1D("hgnSigmaDiff", "nSigma difference", 1000, -5., 5.);
	hKaonCt = new TProfile("hKaonCt", "kaon count in events w.o/ and w/ #Omega", 2, -0.5, 1.5, 0, 1000);
	hKpluspt_omega = new TH1D("hKpluspt_omega", "Global K transver momentum", 1000, 0., 10.);
	hKpluspt_omegabar = new TH1D("hKpluspt_omegabar", "Global K transver momentum", 1000, 0., 10.);
	hKminuspt_omega = new TH1D("hKminuspt_omega", "Global K transver momentum", 1000, 0., 10.);
	hKminuspt_omegabar = new TH1D("hKminuspt_omegabar", "Global K transver momentum", 1000, 0., 10.);
	hKpluseta_omega = new TH1D("hKpluseta_omega", "Global K transver momentum", 1000, -5., 5.);
	hKpluseta_omegabar = new TH1D("hKpluseta_omegabar", "Global K transver momentum", 1000, -5., 5.);
	hKminuseta_omega = new TH1D("hKminuseta_omega", "Global K transver momentum", 1000, -5., 5.);
	hKminuseta_omegabar = new TH1D("hKminuseta_omegabar", "Global K transver momentum", 1000, -5., 5.);
	hKplusphi_omega = new TH1D("hKplusphi_omega", "Global K transver momentum", 1000, -pi, pi);
	hKplusphi_omegabar = new TH1D("hKplusphi_omegabar", "Global K transver momentum", 1000, -pi, pi);
	hKminusphi_omega = new TH1D("hKminusphi_omega", "Global K transver momentum", 1000, -pi, pi);
	hKminusphi_omegabar = new TH1D("hKminusphi_omegabar", "Global K transver momentum", 1000, -pi, pi);
	hKplusy_omega = new TH1D("hKplusy_omega", "K+ rapidity", 1000, -5., 5.);
	hKplusy_omegabar = new TH1D("hKplusy_omegabar", "K+ rapidity", 1000, -5., 5.);
	hKminusy_omega = new TH1D("hKminusy_omega", "K- rapidity", 1000, -5., 5.);
	hKminusy_omegabar = new TH1D("hKminusy_omegabar", "K- rapidity", 1000, -5., 5.);
	// sideband
	hKpluspt_omega_sideband = new TH1D("hKpluspt_omega_sideband", "Global K transver momentum", 1000, 0., 10.);
	hKpluspt_omegabar_sideband = new TH1D("hKpluspt_omegabar_sideband", "Global K transver momentum", 1000, 0., 10.);
	hKminuspt_omega_sideband = new TH1D("hKminuspt_omega_sideband", "Global K transver momentum", 1000, 0., 10.);
	hKminuspt_omegabar_sideband = new TH1D("hKminuspt_omegabar_sideband", "Global K transver momentum", 1000, 0., 10.);
	hKpluseta_omega_sideband = new TH1D("hKpluseta_omega_sideband", "Global K transver momentum", 1000, -5., 5.);
	hKpluseta_omegabar_sideband = new TH1D("hKpluseta_omegabar_sideband", "Global K transver momentum", 1000, -5., 5.);
	hKminuseta_omega_sideband = new TH1D("hKminuseta_omega_sideband", "Global K transver momentum", 1000, -5., 5.);
	hKminuseta_omegabar_sideband = new TH1D("hKminuseta_omegabar_sideband", "Global K transver momentum", 1000, -5., 5.);
	hKplusphi_omega_sideband = new TH1D("hKplusphi_omega_sideband", "Global K transver momentum", 1000, -pi, pi);
	hKplusphi_omegabar_sideband = new TH1D("hKplusphi_omegabar_sideband", "Global K transver momentum", 1000, -pi, pi);
	hKminusphi_omega_sideband = new TH1D("hKminusphi_omega_sideband", "Global K transver momentum", 1000, -pi, pi);
	hKminusphi_omegabar_sideband = new TH1D("hKminusphi_omegabar_sideband", "Global K transver momentum", 1000, -pi, pi);
	hKplusy_omega_sideband = new TH1D("hKplusy_omega_sideband", "K+ rapidity", 1000, -5., 5.);
	hKplusy_omegabar_sideband = new TH1D("hKplusy_omegabar_sideband", "K+ rapidity", 1000, -5., 5.);
	hKminusy_omega_sideband = new TH1D("hKminusy_omega_sideband", "K- rapidity", 1000, -5., 5.);
	hKminusy_omegabar_sideband = new TH1D("hKminusy_omegabar_sideband", "K- rapidity", 1000, -5., 5.);

	// w/ or wo/ lambda
	hKpluspt_wlb_omega = new TH1D("hKpluspt_wlb_omega", "Global K transver momentum", 1000, 0., 10.);
	hKpluspt_wlb_omegabar = new TH1D("hKpluspt_wlb_omegabar", "Global K transver momentum", 1000, 0., 10.);
	hKminuspt_wl_omega = new TH1D("hKminuspt_wl_omega", "Global K transver momentum", 1000, 0., 10.);
	hKminuspt_wl_omegabar = new TH1D("hKminuspt_wl_omegabar", "Global K transver momentum", 1000, 0., 10.);
	hKplusy_wlb_omega = new TH1D("hKplusy_wlb_omega", "K+ rapidity", 1000, -5., 5.);
	hKplusy_wlb_omegabar = new TH1D("hKplusy_wlb_omegabar", "K+ rapidity", 1000, -5., 5.);
	hKminusy_wl_omega = new TH1D("hKminusy_wl_omega", "K- rapidity", 1000, -5., 5.);
	hKminusy_wl_omegabar = new TH1D("hKminusy_wl_omegabar", "K- rapidity", 1000, -5., 5.);
	hKpluspt_wolb_omega = new TH1D("hKpluspt_wolb_omega", "Global K transver momentum", 1000, 0., 10.);
	hKpluspt_wolb_omegabar = new TH1D("hKpluspt_wolb_omegabar", "Global K transver momentum", 1000, 0., 10.);
	hKminuspt_wol_omega = new TH1D("hKminuspt_wol_omega", "Global K transver momentum", 1000, 0., 10.);
	hKminuspt_wol_omegabar = new TH1D("hKminuspt_wol_omegabar", "Global K transver momentum", 1000, 0., 10.);
	hKplusy_wolb_omega = new TH1D("hKplusy_wolb_omega", "K+ rapidity", 1000, -5., 5.);
	hKplusy_wolb_omegabar = new TH1D("hKplusy_wolb_omegabar", "K+ rapidity", 1000, -5., 5.);
	hKminusy_wol_omega = new TH1D("hKminusy_wol_omega", "K- rapidity", 1000, -5., 5.);
	hKminusy_wol_omegabar = new TH1D("hKminusy_wol_omegabar", "K- rapidity", 1000, -5., 5.);
	// and their sidebands
	hKpluspt_wlb_omega_sideband = new TH1D("hKpluspt_wlb_omega_sideband", "Global K transver momentum", 1000, 0., 10.);
	hKpluspt_wlb_omegabar_sideband = new TH1D("hKpluspt_wlb_omegabar_sideband", "Global K transver momentum", 1000, 0., 10.);
	hKminuspt_wl_omega_sideband = new TH1D("hKminuspt_wl_omega_sideband", "Global K transver momentum", 1000, 0., 10.);
	hKminuspt_wl_omegabar_sideband = new TH1D("hKminuspt_wl_omegabar_sideband", "Global K transver momentum", 1000, 0., 10.);
	hKplusy_wlb_omega_sideband = new TH1D("hKplusy_wlb_omega_sideband", "K+ rapidity", 1000, -5., 5.);
	hKplusy_wlb_omegabar_sideband = new TH1D("hKplusy_wlb_omegabar_sideband", "K+ rapidity", 1000, -5., 5.);
	hKminusy_wl_omega_sideband = new TH1D("hKminusy_wl_omega_sideband", "K- rapidity", 1000, -5., 5.);
	hKminusy_wl_omegabar_sideband = new TH1D("hKminusy_wl_omegabar_sideband", "K- rapidity", 1000, -5., 5.);
	hKpluspt_wolb_omega_sideband = new TH1D("hKpluspt_wolb_omega_sideband", "Global K transver momentum", 1000, 0., 10.);
	hKpluspt_wolb_omegabar_sideband = new TH1D("hKpluspt_wolb_omegabar_sideband", "Global K transver momentum", 1000, 0., 10.);
	hKminuspt_wol_omega_sideband = new TH1D("hKminuspt_wol_omega_sideband", "Global K transver momentum", 1000, 0., 10.);
	hKminuspt_wol_omegabar_sideband = new TH1D("hKminuspt_wol_omegabar_sideband", "Global K transver momentum", 1000, 0., 10.);
	hKplusy_wolb_omega_sideband = new TH1D("hKplusy_wolb_omega_sideband", "K+ rapidity", 1000, -5., 5.);
	hKplusy_wolb_omegabar_sideband = new TH1D("hKplusy_wolb_omegabar_sideband", "K+ rapidity", 1000, -5., 5.);
	hKminusy_wol_omega_sideband = new TH1D("hKminusy_wol_omega_sideband", "K- rapidity", 1000, -5., 5.);
	hKminusy_wol_omegabar_sideband = new TH1D("hKminusy_wol_omegabar_sideband", "K- rapidity", 1000, -5., 5.);

	// for a rapidity check
	hkplus_totaly = new TProfile("hkplus_totaly", "K+ total rapidity in an event", 2, -0.5, 1.5, 0., 500.);
	hkminus_totaly = new TProfile("hkminus_totaly", "K- total rapidity in an event", 2, -0.5, 1.5, 0., 500.);
	hkplus_ntrack = new TProfile("hkplus_ntrack", "K+ number of tracks in an event", 2, -0.5, 1.5, 0., 500.);
	hkminus_ntrack = new TProfile("hkminus_ntrack", "K- number of tracks in an event", 2, -0.5, 1.5, 0., 500.);

	// checking east west asymmetry
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			for (int k = 0; k < 2; k++)
			{
				hPhi_EW_PN_EWVtx[i][j][k] = new TH1D(Form("hPhi_EW_PN_EWVtx_%d_%d_%d", i, j, k), Form("hPhi_EW_PN_EWVtx_%d_%d_%d", i, j, k), 1000, -pi, pi);
				hPhi_EW_PN_EWVtx[i][j][k]->Sumw2();
				hPhiXi_EW_PN_EWVtx[i][j][k] = new TH1D(Form("hPhiXi_EW_PN_EWVtx_%d_%d_%d", i, j, k), Form("hPhiXi_EW_PN_EWVtx_%d_%d_%d", i, j, k), 1000, -pi, pi);
				hPhiXi_EW_PN_EWVtx[i][j][k]->Sumw2();
			}
		}
	}

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
	hProtony = new TH1D("hProtony", "Proton Rapidity", 1000, -5., 5.);
	hAntiProtony = new TH1D("hAntiProtony", "Anti-proton Rapidity", 1000, -5., 5.);
	hm2proton_b = new TH1D("hm2proton_b", "m^2 before", 500, -1., 4.);
	hm2proton_r = new TH1D("hm2proton_r", "m^2 regular", 500, -1., 4.);
	hm2proton_a = new TH1D("hm2proton_a", "m^2 after", 500, -1., 4.);
	hm2pion_b = new TH1D("hm2pion_b", "m^2 before", 500, -1., 4.);
	hm2pion_r = new TH1D("hm2pion_r", "m^2 regular", 500, -1., 4.);
	hm2pion_a = new TH1D("hm2pion_a", "m^2 after", 500, -1., 4.);
	hTOFEff_check = new TProfile("hTOFEff_check", "hTOFEff_check", 1000, 0., 10., 0., 1.);
	for (int i = 0; i < 9; i++)
		hTPCEff_check[i] = new TProfile(Form("hTPCEff_check_%d", i + 1), Form("hTPCEff_check_%d", i + 1), 1000, 0., 10., 0., 1.);

	// excited states
	for (int i = 0; i < 9; i++)
	{
		hXi1530M_cen[i] = new TH1D(Form("hXi1530M_cen_%d", i + 1), Form("hXi1530M_cen_%d", i + 1), 1200, 1.2, 1.8);
		hXi1530barM_cen[i] = new TH1D(Form("hXi1530barM_cen_%d", i + 1), Form("hXi1530barM_cen_%d", i + 1), 1200, 1.2, 1.8);
		hKsM_cen[i] = new TH1D(Form("hKsM_cen_%d", i + 1), Form("hKsM_cen_%d", i + 1), 1200, 0.2, 0.8);
		hOmega2012M_cen[i] = new TH1D(Form("hOmega2012M_cen_%d", i + 1), Form("hOmega2012M_cen_%d", i + 1), 1200, 1.7, 2.3);
		hOmega2012barM_cen[i] = new TH1D(Form("hOmega2012barM_cen_%d", i + 1), Form("hOmega2012barM_cen_%d", i + 1), 1200, 1.7, 2.3);
	}

	// Omega QA
	hOmegaM = new TH1D("hOmegaM", "Omega Invariant Mass", 1400, 1., 2.4);
	hOmegap = new TH1D("hOmegap", "Omega Momentum", 1000, 0., 10.);
	hOmegapt = new TH1D("hOmegapt", "Omega Transverse Momentum", 1000, 0., 10.);
	hOmegaeta = new TH1D("hOmegaeta", "Omega Pseudurapidity", 1000, -5., 5.);
	hOmegay = new TH1D("hOmegay", "Omega Rapidity", 1000, -5., 5.);
	hOmegaypt = new TH2D("hOmegaypt", "Omega Rapidity vs pT", 1000, 0., 10., 1000, -5., 5.);
	hOmegaphi = new TH1D("hOmegaphi", "Omega Phi", 1000, -pi, pi);
	hOmegaDL = new TH1D("hOmegaDL", "Omega Decay Length", 1000, 0., 10.);
	hOmegaDLminusLambdaDL = new TH1D("hOmegaDLminusLambdaDL", "hOmegaDLminusLambdaDL", 1000, -5., 5.);

	hOmegabarM = new TH1D("hOmegabarM", "Omegabar Invariant Mass", 1400, 1., 2.4);
	hOmegabarp = new TH1D("hOmegabarp", "Omegabar Momentum", 1000, 0., 10.);
	hOmegabarpt = new TH1D("hOmegabarpt", "Omegabar Transverse Momentum", 1000, 0., 10.);
	hOmegabareta = new TH1D("hOmegabareta", "Omega Pseudurapidity", 1000, -5., 5.);
	hOmegabary = new TH1D("hOmegabary", "Omegabar Rapidity", 1000, -5., 5.);
	hOmegabarypt = new TH2D("hOmegabarypt", "Omegabar Rapidity vs pT", 1000, 0., 10., 1000, -5., 5.);
	hOmegabarphi = new TH1D("hOmegabarphi", "Omegabar Phi", 1000, -pi, pi);
	hOmegabarDL = new TH1D("hOmegabarDL", "Omegabar Decay Length", 1000, 0., 10.);
	hOmegabarDLminusLambdaDL = new TH1D("hOmegabarDLminusLambdaDL", "hOmegabarDLminusLambdaDL", 1000, -5., 5.);

	hNumOmega = new TH1D("hNumOmega", "Number of Omega in an event", 10, -0.5, 9.5);
	hOmegaUsed = new TH1D("hOmegaUsed", "Actual Omega Used #Omega/#bar{#Omega}", 4, -0.5, 3.5);
	hOmegaUsed_sideband = new TH1D("hOmegaUsed_sideband", "Actual Omega Used #Omega/#bar{#Omega}", 2, -0.5, 1.5);
	hOmegaEventUsed = new TH1D("hOmegaEventUsed", "Number of events with Omega", 4, -0.5, 3.5);
	hOmegaEventUsed_sideband = new TH1D("hOmegaEventUsed_sideband", "Number of events with Omega", 2, -0.5, 1.5);
	hLambdaUsed = new TH1D("hLambdaUsed", "Actual Lambda Used #Lambda/#bar{#Lambda}", 4, -0.5, 3.5);
	hLambdaUsed_sideband = new TH1D("hLambdaUsed_sideband", "Actual Lambda Used #Lambda/#bar{#Lambda}", 2, -0.5, 1.5);
	hLambdaEventUsed = new TH1D("hLambdaEventUsed", "Number of events with Lambda", 4, -0.5, 3.5);
	hLambdaEventUsed_sideband = new TH1D("hLambdaEventUsed_sideband", "Number of events with Lambda", 2, -0.5, 1.5);

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
	hLambdaM = new TH1D("hLambdaM", "Lambda Invariant Mass", 1200, 0.6, 1.8);
	hLambdap = new TH1D("hLambdap", "Lambda Momentum", 1000, 0., 10.);
	hLambdapt = new TH1D("hLambdapt", "Lambda Transverse Momentum", 1000, 0., 10.);
	hLambday = new TH1D("hLambday", "Lambda Rapidity", 1000, -5., 5.);
	hLambdaphi = new TH1D("hLambdaphi", "Lambda Phi", 1000, -pi, pi);
	hLambdaDL = new TH1D("hLambdaDL", "Lambda Decay Length", 1000, 0., 10.);
	hLambdabarM = new TH1D("hLambdabarM", "Lambdabar Invariant Mass", 1200, 0.6, 1.8);
	hLambdabarp = new TH1D("hLambdabarp", "Lambdabar Momentum", 1000, 0., 10.);
	hLambdabarpt = new TH1D("hLambdabarpt", "Lambdabar Transverse Momentum", 1000, 0., 10.);
	hLambdabary = new TH1D("hLambdabary", "Lambdabar Rapidity", 1000, -5., 5.);
	hLambdabarphi = new TH1D("hLambdabarphi", "Lambdabar Phi", 1000, -pi, pi);
	hLambdabarDL = new TH1D("hLambdabarDL", "Lambdabar Decay Length", 1000, 0., 10.);

	// Daughter kaon QA
	hDauKplusp = new TH1D("hDauKplusp", "Daughter K+ Momentum", 1000, 0., 10.);
	hDauKpluspt = new TH1D("hDauKpluspt", "Daughter K+ Transverse Momentum", 1000, 0., 10.);
	hDauKplusy = new TH1D("hDauKplusy", "Daughter K+ Rapidity", 1000, -5., 5.);
	hDauKplusphi = new TH1D("hDauKplusphi", "Daughter K+ Phi", 1000, -pi, pi);
	hDauKplusnSigma = new TH1D("hDauKplusnSigma", "Daughter K+ nSigma", 1000, -10, 10);
	hDauKminusp = new TH1D("hDauKminusp", "Daughter K- Momentum", 1000, 0., 10.);
	hDauKminuspt = new TH1D("hDauKminuspt", "Daughter K- Transverse Momentum", 1000, 0., 10.);
	hDauKminusy = new TH1D("hDauKminusy", "Daughter K- Rapidity", 1000, -5., 5.);
	hDauKminusphi = new TH1D("hDauKminusphi", "Daughter K- Phi", 1000, -pi, pi);
	hDauKminusnSigma = new TH1D("hDauKminusnSigma", "Daughter K+ nSigma", 1000, -10, 10);

	// Daughter Lambda QA
	hDauLambdaM = new TH1D("hDauLambdaM", "Daughter Lambda Invariant Mass", 1200, 0.6, 1.8);
	hDauLambdap = new TH1D("hDauLambdap", "Daughter Lambda Momentum", 1000, 0., 10.);
	hDauLambdapt = new TH1D("hDauLambdapt", "Daughter Lambda Transverse Momentum", 1000, 0., 10.);
	hDauLambday = new TH1D("hDauLambday", "Daughter Lambda Rapidity", 1000, -5., 5.);
	hDauLambdaphi = new TH1D("hDauLambdaphi", "Daughter Lambda Phi", 1000, -pi, pi);
	hDauLambdaDL = new TH1D("hDauLambdaDL", "Daughter Lambda Decay Length", 1000, 0., 10.);
	hDauLambdabarM = new TH1D("hDauLambdabarM", "Daughter Lambdabar Invariant Mass", 1200, 0.6, 1.8);
	hDauLambdabarp = new TH1D("hDauLambdabarp", "Daughter Lambdabar Momentum", 1000, 0., 10.);
	hDauLambdabarpt = new TH1D("hDauLambdabarpt", "Daughter Lambdabar Transverse Momentum", 1000, 0., 10.);
	hDauLambdabary = new TH1D("hDauLambdabary", "Daughter Lambdabar Rapidity", 1000, -5., 5.);
	hDauLambdabarphi = new TH1D("hDauLambdabarphi", "Daughter Lambdabar Phi", 1000, -pi, pi);
	hDauLambdabarDL = new TH1D("hDauLambdabarDL", "Daughter Lambdabar Decay Length", 1000, 0., 10.);

	// Daughter proton and pion
	hDauProtonp = new TH1D("hDauProtonp", "Daughter Proton Momentum", 1000, 0., 10.);
	hDauProtonpt = new TH1D("hDauProtonpt", "Daughter Proton Transverse Momentum", 1000, 0., 10.);
	hDauPionp = new TH1D("hDauPionp", "Daughter Pion Momentum", 1000, 0., 10.);
	hDauPionpt = new TH1D("hDauPionpt", "Daughter Pion Transverse Momentum", 1000, 0., 10.);
	hDauProtonbarp = new TH1D("hDauProtonbarp", "Daughter Protonbar Momentum", 1000, 0., 10.);
	hDauProtonbarpt = new TH1D("hDauProtonbarpt", "Daughter Protonbar Transverse Momentum", 1000, 0., 10.);
	hDauPionbarp = new TH1D("hDauPionbarp", "Daughter Pionbar Momentum", 1000, 0., 10.);
	hDauPionbarpt = new TH1D("hDauPionbarpt", "Daughter Pionbar Transverse Momentum", 1000, 0., 10.);

	hOmegaDauPid = new TH1D("hOmegaDauPid", "Omega Daughter PID", 10000, -4999.5, 5000.5);
	hOmegabarDauPid = new TH1D("hOmegabarDauPid", "Omegabar Daughter PID", 10000, -4999.5, 5000.5);

	// Cut Decider
	hDCAOtoK_signal = new TH1D("hDCAOtoK_signal", "hDCAOtoK_signal", 1000, 0., 10.);
	hDCAOtoK_sideband = new TH1D("hDCAOtoK_sideband", "hDCAOtoK_sideband", 1000, 0., 10.);
	hDCAOtoL_signal = new TH1D("hDCAOtoL_signal", "hDCAOtoL_signal", 1000, 0., 10.);
	hDCAOtoL_sideband = new TH1D("hDCAOtoL_sideband", "hDCAOtoL_sideband", 1000, 0., 10.);

	// event mixing QA
	hNumMixedEvent = new TH1D("hNumMixedEvent", "Number of mixed events", 20, -0.5, 19.5);
	hTotalMixedEvent = new TH1D("hTotalMixedEvent", "hTotalMixedEvent", 8000, -0.5, 7999.5);

	// Run-by-Run QA
	hEPD_full_1_runID = new TProfile("hEPD_full_1_runID", "hEPD_full_1_runID", 40000, 19119999.5, 19159999.5, 0., 2 * pi);
	hEPD_full_2_runID = new TProfile("hEPD_full_2_runID", "hEPD_full_2_runID", 40000, 19119999.5, 19159999.5, 0., 2 * pi);
	hEPD_full_1_day = new TProfile("hEPD_full_1_day", "hEPD_full_1_day", 300, -0.5, 299.5, 0., 2 * pi);
	hEPD_full_2_day = new TProfile("hEPD_full_2_day", "hEPD_full_2_day", 300, -0.5, 299.5, 0., 2 * pi);

	// xiatong's analysis
	std::string KO_comb_name[4] = {"KplusO", "KplusObar", "KminusO", "KminusObar"};
	std::string LO_comb_name[4] = {"LambdabarO", "LambdabarObar", "LambdaO", "LambdaObar"};
	std::string KL_comb_name[4] = {"KplusLambda", "KplusLambdabar", "KminusLambda", "KminusLambdabar"};
	std::string KP_comb_name[4] = {"KplusProton", "KplusAntiProton", "KminusProton", "KminusAntiProton"};
	for (int i = 0; i < 4; i++)
	{
		hCorrKO[i] = new TH1D(Form("hCorr%s", KO_comb_name[i].c_str()), Form("hCorr%s", KO_comb_name[i].c_str()), 5000, 0.0, 50.0);
		hPtCorrKO[i] = new TH1D(Form("hPtCorr%s", KO_comb_name[i].c_str()), Form("hPtCorr%s", KO_comb_name[i].c_str()), 5000, 0.0, 50.0);
		hyCorrKO[i] = new TH1D(Form("hyCorr%s", KO_comb_name[i].c_str()), Form("hyCorr%s", KO_comb_name[i].c_str()), 2000, -10.0, 10.0);
		hphiCorrKO[i] = new TH1D(Form("hphiCorr%s", KO_comb_name[i].c_str()), Form("hphiCorr%s", KO_comb_name[i].c_str()), 1000, -pi, pi);
		hCorrKO_same[i] = new TH1D(Form("hCorr%s_same", KO_comb_name[i].c_str()), Form("hCorr%s_same", KO_comb_name[i].c_str()), 5000, 0.0, 50.0);
		hCorrKO_mixed[i] = new TH1D(Form("hCorr%s_mixed", KO_comb_name[i].c_str()), Form("hCorr%s_mixed", KO_comb_name[i].c_str()), 5000, 0.0, 50.0);
		hCorrKO_sideband[i] = new TH1D(Form("hCorr%s_sideband", KO_comb_name[i].c_str()), Form("hCorr%s_sideband", KO_comb_name[i].c_str()), 5000, 0.0, 50.0);
		hPtCorrKO_sideband[i] = new TH1D(Form("hPtCorr%s_sideband", KO_comb_name[i].c_str()), Form("hPtCorr%s_sideband", KO_comb_name[i].c_str()), 5000, 0.0, 50.0);
		hyCorrKO_sideband[i] = new TH1D(Form("hyCorr%s_sideband", KO_comb_name[i].c_str()), Form("hyCorr%s_sideband", KO_comb_name[i].c_str()), 2000, -10.0, 10.0);
		hphiCorrKO_sideband[i] = new TH1D(Form("hphiCorr%s_sideband", KO_comb_name[i].c_str()), Form("hphiCorr%s_sideband", KO_comb_name[i].c_str()), 1000, -pi, pi);

		hCorrLO[i] = new TH1D(Form("hCorr%s", LO_comb_name[i].c_str()), Form("hCorr%s", LO_comb_name[i].c_str()), 5000, 0.0, 50.0);
		hyCorrLO[i] = new TH1D(Form("hyCorr%s", LO_comb_name[i].c_str()), Form("hyCorr%s", LO_comb_name[i].c_str()), 2000, -10.0, 10.0);
		hphiCorrLO[i] = new TH1D(Form("hphiCorr%s", LO_comb_name[i].c_str()), Form("hphiCorr%s", LO_comb_name[i].c_str()), 1000, -pi, pi);
		hCorrLO_same[i] = new TH1D(Form("hCorr%s_same", LO_comb_name[i].c_str()), Form("hCorr%s_same", LO_comb_name[i].c_str()), 5000, 0.0, 50.0);
		hCorrLO_mixed[i] = new TH1D(Form("hCorr%s_mixed", LO_comb_name[i].c_str()), Form("hCorr%s_mixed", LO_comb_name[i].c_str()), 5000, 0.0, 50.0);
		hCorrLO_sideband[i] = new TH1D(Form("hCorr%s_sideband", LO_comb_name[i].c_str()), Form("hCorr%s_sideband", LO_comb_name[i].c_str()), 5000, 0.0, 50.0);
		hyCorrLO_sideband[i] = new TH1D(Form("hyCorr%s_sideband", LO_comb_name[i].c_str()), Form("hyCorr%s_sideband", LO_comb_name[i].c_str()), 2000, -10.0, 10.0);
		hphiCorrLO_sideband[i] = new TH1D(Form("hphiCorr%s_sideband", LO_comb_name[i].c_str()), Form("hphiCorr%s_sideband", LO_comb_name[i].c_str()), 1000, -pi, pi);

		hCorrKO_2D_same[i] = new TH2D(Form("hCorr%s_2D_same", KO_comb_name[i].c_str()), Form("hCorr%s_2D_same", KO_comb_name[i].c_str()), 1000, -5.0, 5.0, 1000, -pi, pi);
		hCorrKO_2D_mixed[i] = new TH2D(Form("hCorr%s_2D_mixed", KO_comb_name[i].c_str()), Form("hCorr%s_2D_mixed", KO_comb_name[i].c_str()), 1000, -5.0, 5.0, 1000, -pi, pi);
		hCorrLO_2D_same[i] = new TH2D(Form("hCorr%s_2D_same", LO_comb_name[i].c_str()), Form("hCorr%s_2D_same", LO_comb_name[i].c_str()), 1000, -5.0, 5.0, 1000, -pi, pi);
		hCorrLO_2D_mixed[i] = new TH2D(Form("hCorr%s_2D_mixed", LO_comb_name[i].c_str()), Form("hCorr%s_2D_mixed", LO_comb_name[i].c_str()), 1000, -5.0, 5.0, 1000, -pi, pi);
	
		// Lambda
		hCorrKL[i] = new TH1D(Form("hCorr%s", KL_comb_name[i].c_str()), Form("hCorr%s", KL_comb_name[i].c_str()), 5000, 0.0, 50.0);
		hPtCorrKL[i] = new TH1D(Form("hPtCorr%s", KL_comb_name[i].c_str()), Form("hPtCorr%s", KL_comb_name[i].c_str()), 5000, 0.0, 50.0);
		hyCorrKL[i] = new TH1D(Form("hyCorr%s", KL_comb_name[i].c_str()), Form("hyCorr%s", KL_comb_name[i].c_str()), 2000, -10.0, 10.0);
		hphiCorrKL[i] = new TH1D(Form("hphiCorr%s", KL_comb_name[i].c_str()), Form("hphiCorr%s", KL_comb_name[i].c_str()), 1000, -pi, pi);
		hCorrKL_same[i] = new TH1D(Form("hCorr%s_same", KL_comb_name[i].c_str()), Form("hCorr%s_same", KL_comb_name[i].c_str()), 5000, 0.0, 50.0);
		hCorrKL_mixed[i] = new TH1D(Form("hCorr%s_mixed", KL_comb_name[i].c_str()), Form("hCorr%s_mixed", KL_comb_name[i].c_str()), 5000, 0.0, 50.0);
		hCorrKL_sideband[i] = new TH1D(Form("hCorr%s_sideband", KL_comb_name[i].c_str()), Form("hCorr%s_sideband", KL_comb_name[i].c_str()), 5000, 0.0, 50.0);
		hPtCorrKL_sideband[i] = new TH1D(Form("hPtCorr%s_sideband", KL_comb_name[i].c_str()), Form("hPtCorr%s_sideband", KL_comb_name[i].c_str()), 5000, 0.0, 50.0);
		hyCorrKL_sideband[i] = new TH1D(Form("hyCorr%s_sideband", KL_comb_name[i].c_str()), Form("hyCorr%s_sideband", KL_comb_name[i].c_str()), 2000, -10.0, 10.0);
		hphiCorrKL_sideband[i] = new TH1D(Form("hphiCorr%s_sideband", KL_comb_name[i].c_str()), Form("hphiCorr%s_sideband", KL_comb_name[i].c_str()), 1000, -pi, pi);

		// Proton
		hCorrKP[i] = new TH1D(Form("hCorr%s", KP_comb_name[i].c_str()), Form("hCorr%s", KP_comb_name[i].c_str()), 5000, 0.0, 50.0);
		hPtCorrKP[i] = new TH1D(Form("hPtCorr%s", KP_comb_name[i].c_str()), Form("hPtCorr%s", KP_comb_name[i].c_str()), 5000, 0.0, 50.0);
		hyCorrKP[i] = new TH1D(Form("hyCorr%s", KP_comb_name[i].c_str()), Form("hyCorr%s", KP_comb_name[i].c_str()), 2000, -10.0, 10.0);
		hphiCorrKP[i] = new TH1D(Form("hphiCorr%s", KP_comb_name[i].c_str()), Form("hphiCorr%s", KP_comb_name[i].c_str()), 1000, -pi, pi);
		hCorrKP_same[i] = new TH1D(Form("hCorr%s_same", KP_comb_name[i].c_str()), Form("hCorr%s_same", KP_comb_name[i].c_str()), 5000, 0.0, 50.0);
		hCorrKP_mixed[i] = new TH1D(Form("hCorr%s_mixed", KP_comb_name[i].c_str()), Form("hCorr%s_mixed", KP_comb_name[i].c_str()), 5000, 0.0, 50.0);
	}

	// some QA about tracing back to primary vertex
	hOmegaDCAtoPV = new TH1D("hOmegaDCAtoPV", "hOmegaDCAtoPV", 1000, 0.0, 10.0);
	hOmegabarDCAtoPV = new TH1D("hOmegabarDCAtoPV", "hOmegabarDCAtoPV", 1000, 0.0, 10.0);
	hLambdaDCAtoPV = new TH1D("hLambdaDCAtoPV", "hLambdaDCAtoPV", 1000, 0.0, 10.0);
	hLambdabarDCAtoPV = new TH1D("hLambdabarDCAtoPV", "hLambdabarDCAtoPV", 1000, 0.0, 10.0);
	hOmegaPDiff = new TH1D("hOmegaPDiff", "hOmegaPDiff", 500, 0.0, 5.0);
	hOmegabarPDiff = new TH1D("hOmegabarPDiff", "hOmegabarPDiff", 500, 0.0, 5.0);
	hLambdaPDiff = new TH1D("hLambdaPDiff", "hLambdaPDiff", 500, 0.0, 5.0);
	hLambdabarPDiff = new TH1D("hLambdabarPDiff", "hLambdabarPDiff", 500, 0.0, 5.0);

	hNegPtDiff_dphi_KmOb = new TH1D("hNegPtDiff_dphi_KmOb", "#Delta#phi for pairs with negative #Delta p_{T}", 1000, -pi, pi);
	hPosPtDiff_dphi_KmOb = new TH1D("hPosPtDiff_dphi_KmOb", "#Delta#phi for pairs with positive #Delta p_{T}", 1000, -pi, pi);
	hNegPtDiff_dphi_KmO = new TH1D("hNegPtDiff_dphi_KmO", "#Delta#phi for pairs with negative #Delta p_{T}", 1000, -pi, pi);
	hPosPtDiff_dphi_KmO = new TH1D("hPosPtDiff_dphi_KmO", "#Delta#phi for pairs with positive #Delta p_{T}", 1000, -pi, pi);

	// hCorrKplusO_y_pT = new TH2D("hCorrKplusO_y_pT", "hCorrKplusO_y_pT", 500, 0.0, 5.0, 200, 0.0, 2.0);
	// hCorrKplusObar_y_pT = new TH2D("hCorrKplusObar_y_pT", "hCorrKplusObar_y_pT", 500, 0.0, 5.0, 200, 0.0, 2.0);
	// hCorrKminusO_y_pT = new TH2D("hCorrKminusO_y_pT", "hCorrKminusO_y_pT", 500, 0.0, 5.0, 200, 0.0, 2.0);
	// hCorrKminusObar_y_pT = new TH2D("hCorrKminusObar_y_pT", "hCorrKminusObar_y_pT", 500, 0.0, 5.0, 200, 0.0, 2.0);
	// hCorrKplusO_y_pT_mixed = new TH2D("hCorrKplusO_y_pT_mixed", "hCorrKplusO_y_pT_mixed", 500, 0.0, 5.0, 200, 0.0, 2.0);
	// hCorrKplusObar_y_pT_mixed = new TH2D("hCorrKplusObar_y_pT_mixed", "hCorrKplusObar_y_pT_mixed", 500, 0.0, 5.0, 200, 0.0, 2.0);
	// hCorrKminusO_y_pT_mixed = new TH2D("hCorrKminusO_y_pT_mixed", "hCorrKminusO_y_pT_mixed", 500, 0.0, 5.0, 200, 0.0, 2.0);
	// hCorrKminusObar_y_pT_mixed = new TH2D("hCorrKminusObar_y_pT_mixed", "hCorrKminusObar_y_pT_mixed", 500, 0.0, 5.0, 200, 0.0, 2.0);

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
	hKratio_omega = new TProfile("hKratio_omega", "hKratio_omega", 20, 0., 1., 0., 10.);
	hKratio_wo = new TProfile("hKratio_wo", "hKratio_wo", 20, 0., 1., 0., 10.);
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
	sprintf(temp, "./mix/omega_mix_cen_%d%d.root", min_cent, max_cent);
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
void StKFParticleAnalysisMaker::WriteHistograms()
{

	///////////////////
	hEventQA->Write();
	hNRefMult->Write();
	hNRefMultA->Write();
	hNRefMultB->Write();
	hVertexXY->Write();
	hVertexZ->Write();
	hVertex2D->Write();
	hDiffVz->Write();
	hcent->Write();
	hcentw->Write();

	hcentRefM->Write();
	hcentRefW->Write();
	for (int i = 0; i < 9; i++)
		hRefMultCorr_cent[i]->Write();
	hMultRatioEW_VertexZ->Write();

	hDauProtonFirstPoint_lam_pt->Write();
	hDauProtonLastPoint_lam_pt->Write();
	hDauPionFirstPoint_lam_pt->Write();
	hDauPionLastPoint_lam_pt->Write();

	// y-pt/nq coverage
	for (int i = 0; i < num_RecoPar; i++)
		hRecoPar_y_ptnq[i]->Write();
	for (int i = 0; i < num_IdPar; i++)
		hIdPar_y_ptnq[i]->Write();

	for (int i = 0; i < 4; i++)
	{
		// Omega
		hCorrKO[i]->Write();
		hPtCorrKO[i]->Write();
		hyCorrKO[i]->Write();
		hphiCorrKO[i]->Write();
		hCorrKO_same[i]->Write();
		hCorrKO_mixed[i]->Write();
		hCorrKO_sideband[i]->Write();
		hPtCorrKO_sideband[i]->Write();
		hyCorrKO_sideband[i]->Write();
		hphiCorrKO_sideband[i]->Write();

		hCorrLO[i]->Write();
		hyCorrLO[i]->Write();
		hphiCorrLO[i]->Write();
		hCorrLO_same[i]->Write();
		hCorrLO_mixed[i]->Write();
		hCorrLO_sideband[i]->Write();
		hyCorrLO_sideband[i]->Write();
		hphiCorrLO_sideband[i]->Write();

		hCorrKO_2D_same[i]->Write();
		hCorrKO_2D_mixed[i]->Write();
		hCorrLO_2D_same[i]->Write();
		hCorrLO_2D_mixed[i]->Write();

		// Lambda
		hCorrKL[i]->Write();
		hPtCorrKL[i]->Write();
		hyCorrKL[i]->Write();
		hphiCorrKL[i]->Write();
		hCorrKL_same[i]->Write();
		hCorrKL_mixed[i]->Write();
		hCorrKL_sideband[i]->Write();
		hPtCorrKL_sideband[i]->Write();
		hyCorrKL_sideband[i]->Write();
		hphiCorrKL_sideband[i]->Write();

		// Proton
		hCorrKP[i]->Write();
		hPtCorrKP[i]->Write();
		hyCorrKP[i]->Write();
		hphiCorrKP[i]->Write();
		hCorrKP_same[i]->Write();
		hCorrKP_mixed[i]->Write();
	}

	// some QA about tracing back to primary vertex
	hOmegaDCAtoPV->Write();
	hOmegabarDCAtoPV->Write();
	hLambdaDCAtoPV->Write();
	hLambdabarDCAtoPV->Write();
	hOmegaPDiff->Write();
	hOmegabarPDiff->Write();
	hLambdaPDiff->Write();
	hLambdabarPDiff->Write();

	hNegPtDiff_dphi_KmOb->Write();
	hPosPtDiff_dphi_KmOb->Write();
	hNegPtDiff_dphi_KmO->Write();
	hPosPtDiff_dphi_KmO->Write();
	// hCorrKplusO_y_pT  ->Write(); hCorrKplusObar_y_pT  ->Write(); hCorrKminusO_y_pT  ->Write(); hCorrKminusObar_y_pT  ->Write();
	// hCorrKplusO_y_pT_mixed->Write(); hCorrKplusObar_y_pT_mixed->Write(); hCorrKminusO_y_pT_mixed->Write(); hCorrKminusObar_y_pT_mixed->Write();
	// hCorrKplusO_y_phi ->Write(); hCorrKplusObar_y_phi ->Write(); hCorrKminusO_y_phi ->Write(); hCorrKminusObar_y_phi ->Write();
	// hCorrKplusO_phi_pT->Write(); hCorrKplusObar_phi_pT->Write(); hCorrKminusO_phi_pT->Write(); hCorrKminusObar_phi_pT->Write();

	hKratio_omega->Write();
	hKratio_wo->Write();
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
	for (int i = 0; i < 9; i++)
		hTPCEff_check[i]->Write();

	// v2
	for (int ewFull = 0; ewFull < 3; ewFull++)
	{
		for (int i = 0; i < 9; i++)
		{
			hTPCEP_2[i][ewFull]->Write();
			hTPCEP_2_shifted[i][ewFull]->Write();
			if (CheckWeights2D)
			{
				hTPCEP_2_2D[i][ewFull]->Write();
				hTPCEP_2_2D_shifted[i][ewFull]->Write();
			}
			hTPCEP_2_shift[i][ewFull]->Write();

			hEPDEP_1[i][ewFull]->Write();
			hEPDEP_1_shifted[i][ewFull]->Write();

			hEPDEP_2[i][ewFull]->Write();
			hEPDEP_2_shifted[i][ewFull]->Write();
			if (CheckWeights2D)
			{
				hEPDEP_1_2D[i][ewFull]->Write();
				hEPDEP_1_2D_shifted[i][ewFull]->Write();
				hEPDEP_2_2D[i][ewFull]->Write();
				hEPDEP_2_2D_shifted[i][ewFull]->Write();
			}

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
	hEPDEP_ew_cos_3->Write();

	hEPD_e_EP_1->Write();
	hEPD_w_EP_1->Write();
	hEPD_e_EP_2->Write();
	hEPD_w_EP_2->Write();
	hEPD_full_EP_1->Write();
	hEPD_full_EP_2->Write();
	hEPD_ew_cos->Write();
	hEPDEP_2_Full_Raw_QA_cen4->Write();
	hEPDEP_2_Full_PhiWeighted_QA_cen4->Write();
	hEPDEP_2_Full_PhiWeightedAndShifted_QA_cen4->Write();

	// Reco particle
	for (int i = 0; i < num_RecoPar; i++)
	{
		for (int j = 0; j < 9; j++)
		{
			hRecoParM_cen[i][j]->Write();
			hRecoPar_EPD_v2[i][j]->Write();
			hRecoPar_TPC_v2[i][j]->Write();
		}
	}
	// Lambda
	for (int i = 0; i < num_pt_bin; i++)
	{
		for (int j = 0; j < 9; j++)
		{
			hLambdaM_pt[i][j]->Write();
			hLambdabarM_pt[i][j]->Write();
			hLambda_EPD_v2_pt[i][j]->Write();
			hLambdabar_EPD_v2_pt[i][j]->Write();
			hLambda_TPC_v2_pt[i][j]->Write();
			hLambdabar_TPC_v2_pt[i][j]->Write();
		}
	}
	// Identified particle
	for (int i = 0; i < num_IdPar; i++)
	{
		hIdPar_EPD_v2[i]->Write();
		hIdPar_TPC_v2[i]->Write();
		for (int j = 0; j < 9; j++)
		{
			hIdPar_EPD_v2_pt_1st[i][j]->Write();
			hIdPar_EPD_v2_pt_2nd[i][j]->Write();
			// hIdPar_EPD_v1_pt_1st[i][j]->Write();
			// hIdPar_EPD_v3_pt_1st[i][j]->Write();
			hIdPar_TPC_v2_pt[i][j]->Write();

			if (Coal2D)
			{
				hIdPar_EPD_v1_eta_pt_1st[i][j]->Write();
				hIdPar_EPD_v2_eta_pt_1st[i][j]->Write();
				hIdPar_EPD_v3_eta_pt_1st[i][j]->Write();
			}

			hIdPar_EPD_v1_y[i][j]->Write();
		}
	}
	hyLambdaDauProton_east->Write();
	hyLambdaDauPion_east->Write();
	hyLambdaDauProton_west->Write();
	hyLambdaDauPion_west->Write();

	// excited states
	for (int i = 0; i < 9; i++)
	{
		hXi1530M_cen[i]->Write();
		hXi1530barM_cen[i]->Write();
		hKsM_cen[i]->Write();
		hOmega2012M_cen[i]->Write();
		hOmega2012barM_cen[i]->Write();
	}

	hOmegaM->Write();
	hOmegap->Write();
	hOmegapt->Write();
	hOmegaeta->Write();
	hOmegay->Write();
	hOmegaypt->Write();
	hOmegaphi->Write();
	hOmegaDL->Write();
	hOmegaDLminusLambdaDL->Write();
	hOmegabarM->Write();
	hOmegabarp->Write();
	hOmegabarpt->Write();
	hOmegabareta->Write();
	hOmegabary->Write();
	hOmegabarypt->Write();
	hOmegabarphi->Write();
	hOmegabarDL->Write();
	hOmegabarDLminusLambdaDL->Write();

	hNumOmega->Write();
	hOmegaUsed->Write();
	hOmegaUsed_sideband->Write();
	hOmegaEventUsed->Write();
	hOmegaEventUsed_sideband->Write();
	hLambdaUsed->Write();
	hLambdaUsed_sideband->Write();
	hLambdaEventUsed->Write();
	hLambdaEventUsed_sideband->Write();
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

	hLambdaM->Write();
	hLambdap->Write();
	hLambdapt->Write();
	hLambday->Write();
	hLambdaphi->Write();
	hLambdaDL->Write();
	hLambdabarM->Write();
	hLambdabarp->Write();
	hLambdabarpt->Write();
	hLambdabary->Write();
	hLambdabarphi->Write();
	hLambdabarDL->Write();

	hDauKplusp->Write();
	hDauKpluspt->Write();
	hDauKplusy->Write();
	hDauKplusphi->Write();
	hDauKplusnSigma->Write();
	hDauKminusp->Write();
	hDauKminuspt->Write();
	hDauKminusy->Write();
	hDauKminusphi->Write();
	hDauKminusnSigma->Write();

	hDauLambdaM->Write();
	hDauLambdap->Write();
	hDauLambdapt->Write();
	hDauLambday->Write();
	hDauLambdaphi->Write();
	hDauLambdaDL->Write();

	hDauLambdabarM->Write();
	hDauLambdabarp->Write();
	hDauLambdabarpt->Write();
	hDauLambdabary->Write();
	hDauLambdabarphi->Write();
	hDauLambdabarDL->Write();

	hDauPionp->Write();
	hDauPionpt->Write();
	hDauPionbarp->Write();
	hDauPionbarpt->Write();
	hDauProtonp->Write();
	hDauProtonpt->Write();
	hDauProtonbarp->Write();
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
	hgetofcrossingY->Write();
	hgm2_etof->Write();
	hgpinvbeta_etof->Write();
	hgetofmatchflag->Write();
	hgeta_etof->Write();

	hgpdEdx->Write();
	hgdEdxErr->Write();
	hgpinvbeta->Write();
	hgm2->Write();
	hgpm2->Write();
	hgm2nSigmaKaon->Write();
	hgm2nSigmaPion->Write();
	hgm2nSigmaProton->Write();
	hgptnSigmaKaon->Write();
	hgptnSigmaPion->Write();
	hgptnSigmaProton->Write();
	hgptm2->Write();
	// hgp          ->Write();
	for (int i = 0; i < 9; i++)
	{
		hgp[i]->Write();
		hgp_TOF[i]->Write();
		hgpT[i]->Write();
		hgpT_TOF[i]->Write();
		hgpTeta[i]->Write();
		hgpTeta_TOF[i]->Write();
	}
	hgDCAtoPV->Write();
	hgbtofYlocal->Write();
	// hgKpdEdx       ->Write();
	// hgKpinvbeta    ->Write();
	hgKm2->Write();
	hgKpm2->Write();
	hgKp->Write();
	hgKpT->Write();
	hgKDCAtoPV->Write();
	hgKDCAtoO->Write();
	hgKpionpdEdx->Write();
	// hgptm2_largenSigmaKaon->Write();
	// hgptm2_smallnSigmaKaon->Write();
	hgnSigmaDiff->Write();
	hKaonCt->Write();
	hKpluspt_omega->Write();
	hKpluspt_omegabar->Write();
	hKminuspt_omega->Write();
	hKminuspt_omegabar->Write();
	hKpluseta_omega->Write();
	hKpluseta_omegabar->Write();
	hKminuseta_omega->Write();
	hKminuseta_omegabar->Write();
	hKplusphi_omega->Write();
	hKplusphi_omegabar->Write();
	hKminusphi_omega->Write();
	hKminusphi_omegabar->Write();
	hKplusy_omega->Write();
	hKplusy_omegabar->Write();
	hKminusy_omega->Write();
	hKminusy_omegabar->Write();

	// sideband
	hKpluspt_omega_sideband->Write();
	hKpluspt_omegabar_sideband->Write();
	hKminuspt_omega_sideband->Write();
	hKminuspt_omegabar_sideband->Write();
	hKpluseta_omega_sideband->Write();
	hKpluseta_omegabar_sideband->Write();
	hKminuseta_omega_sideband->Write();
	hKminuseta_omegabar_sideband->Write();
	hKplusphi_omega_sideband->Write();
	hKplusphi_omegabar_sideband->Write();
	hKminusphi_omega_sideband->Write();
	hKminusphi_omegabar_sideband->Write();
	hKplusy_omega_sideband->Write();
	hKplusy_omegabar_sideband->Write();
	hKminusy_omega_sideband->Write();
	hKminusy_omegabar_sideband->Write();

	// w/ and wo/ lambda
	hKpluspt_wlb_omega->Write();
	hKpluspt_wlb_omegabar->Write();
	hKminuspt_wl_omega->Write();
	hKminuspt_wl_omegabar->Write();
	hKplusy_wlb_omega->Write();
	hKplusy_wlb_omegabar->Write();
	hKminusy_wl_omega->Write();
	hKminusy_wl_omegabar->Write();
	hKpluspt_wolb_omega->Write();
	hKpluspt_wolb_omegabar->Write();
	hKminuspt_wol_omega->Write();
	hKminuspt_wol_omegabar->Write();
	hKplusy_wolb_omega->Write();
	hKplusy_wolb_omegabar->Write();
	hKminusy_wol_omega->Write();
	hKminusy_wol_omegabar->Write();
	hKpluspt_wlb_omega_sideband->Write();
	hKpluspt_wlb_omegabar_sideband->Write();
	hKminuspt_wl_omega_sideband->Write();
	hKminuspt_wl_omegabar_sideband->Write();
	hKplusy_wlb_omega_sideband->Write();
	hKplusy_wlb_omegabar_sideband->Write();
	hKminusy_wl_omega_sideband->Write();
	hKminusy_wl_omegabar_sideband->Write();
	hKpluspt_wolb_omega_sideband->Write();
	hKpluspt_wolb_omegabar_sideband->Write();
	hKminuspt_wol_omega_sideband->Write();
	hKminuspt_wol_omegabar_sideband->Write();
	hKplusy_wolb_omega_sideband->Write();
	hKplusy_wolb_omegabar_sideband->Write();
	hKminusy_wol_omega_sideband->Write();
	hKminusy_wol_omegabar_sideband->Write();

	// for the rapidity check
	hkplus_totaly->Write();
	hkminus_totaly->Write();
	hkplus_ntrack->Write();
	hkminus_ntrack->Write();

	// checking east west asymmetry
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			for (int k = 0; k < 2; k++)
			{
				hPhi_EW_PN_EWVtx[i][j][k]->Write();
				hPhiXi_EW_PN_EWVtx[i][j][k]->Write();
			}
		}
	}

	// Run-by-Run QA
	hEPD_full_1_runID->Write();
	hEPD_full_2_runID->Write();
	hEPD_full_1_day->Write();
	hEPD_full_2_day->Write();

#if !USE_P
	// for (int i = 0; i < 10; i++)
	// {
	// 	hgPID2D_proton_pt[i]->Write();
	// 	hgPID2D_antiproton_pt[i]->Write();
	// 	hgPID2D_piplus_pt[i]->Write();
	// 	hgPID2D_piminus_pt[i]->Write();
	// 	hgPID2D_kplus_pt[i]->Write();
	// 	hgPID2D_kminus_pt[i]->Write();
	// }

	// for (int i = 0; i < 12; i++)
	// {
	// 	hgPID2D_proton_p[i]->Write();
	// 	hgPID2D_antiproton_p[i]->Write();
	// 	hgPID2D_piplus_p[i]->Write();
	// 	hgPID2D_piminus_p[i]->Write();
	// 	hgPID2D_kplus_p[i]->Write();
	// 	hgPID2D_kminus_p[i]->Write();
	// }
#else
	for (int i = 0; i < 12; i++)
	{
		hgPID2D_proton_p[i]->Write();
		hgPID2D_antiproton_p[i]->Write();
		hgPID2D_piplus_p[i]->Write();
		hgPID2D_piminus_p[i]->Write();
		hgPID2D_kplus_p[i]->Write();
		hgPID2D_kminus_p[i]->Write();
	}
#endif

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
void StKFParticleAnalysisMaker::setRunEnergyAndListDir(int run, double energy, char ListDir[256])
{
	mRun = run;
	mEnergy = energy;
	mListDir = ListDir;
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
	if (0)
		return kTRUE;
	TString inf = mListDir + "/runList/";
	inf += Form("run%dList%.1f.list", mRun, mEnergy);
	ifstream inrun;
	inrun.open(inf);
	if (inrun.fail())
	{
		cout << "cannot open " << inf.Data() << endl;
		return kFALSE;
	}
	Int_t runid;
	while (inrun >> runid)
	{
		runList.push_back(runid);
	}
	inrun.close();
	sort(runList.begin(), runList.end());

	vector<int>::iterator it;
	it = std::unique(runList.begin(), runList.end());
	runList.resize(std::distance(runList.begin(), it));

	cout << "Run list :" << inf.Data() << " loaded. " << endl;
	cout << "Total    :" << runList.size() << " runs. " << endl;

	if (runList.size() < 1)
	{
		cout << "no run number found!!!" << endl;
		return kFALSE;
	}

	return kTRUE;
}

//-----------------------------------------------------------------------------
bool StKFParticleAnalysisMaker::removeBadID(int runnumber)
{
	for (std::vector<int>::iterator it = badList.begin(); it != badList.end(); ++it)
	{
		if (runnumber == *it)
		{
			return kTRUE;
		}
	}
	return kFALSE;
}

//-----------------------------------------------------------------------------
Int_t StKFParticleAnalysisMaker::CheckrunNumber(int runnumber)
{
	int pointer = -999;
	int id = 0;
	for (std::vector<int>::iterator it = runList.begin(); it != runList.end(); ++it)
	{
		if (runnumber == *it)
			pointer = id;
		id++;
	}

	if (pointer == -999)
		cout << "Run number are not found! " << runnumber << endl;
	return pointer;
}

//-----------------------------------------------------------------------------
Int_t StKFParticleAnalysisMaker::findCentrality(int mult)
{
	// NOT in use
	int centrality = 0;
	float centFull[9] = {30, 46, 65, 89, 118, 154, 200, 262, 305};
	if (mult > centFull[8])
		centrality = 9;
	else if (mult > centFull[7])
		centrality = 8;
	else if (mult > centFull[6])
		centrality = 7;
	else if (mult > centFull[5])
		centrality = 6;
	else if (mult > centFull[4])
		centrality = 5;
	else if (mult > centFull[3])
		centrality = 4;
	else if (mult > centFull[2])
		centrality = 3;
	else if (mult > centFull[1])
		centrality = 2;
	else if (mult > centFull[0])
		centrality = 1;

	return centrality;
}

//-----------------------------------------------------------------------------
void StKFParticleAnalysisMaker::Clear(Option_t *opt)
{
}

//-----------------------------------------------------------------------------
Int_t StKFParticleAnalysisMaker::Make()
{
	// time_t time_start;
	// time_t time_now;
	// time(&time_start);
	hEventQA->Fill(0.);

	PicoDst = StPicoDst::instance();
	StPicoDst *mPicoDst = PicoDst;
	if (!mPicoDst)
	{
		LOG_WARN << " No PicoDst! Skip! " << endm;
		return kStOK;
	}

	hEventQA->Fill(1.);
	//     pass event
	/////////////////////////////////////////////////////////
	StPicoEvent *mEvent = (StPicoEvent *)mPicoDst->event();
	if (!mEvent)
		return kStOK;

	hEventQA->Fill(2.);
	const int runID = mEvent->runId();
	const int evtID = mEvent->eventId();
	const int refMult = mEvent->refMult();
	const int grefMult = mEvent->grefMult();
	const int ranking = mEvent->ranking();
	const int tofMatch = mEvent->nBTOFMatch();

	const double magnet = mEvent->bField();

	if ((!mEvent->isTrigger(610001)) && (!mEvent->isTrigger(610011)) && (!mEvent->isTrigger(610021)) && (!mEvent->isTrigger(610031)) && (!mEvent->isTrigger(610041)) && (!mEvent->isTrigger(610051)))
		return kStOK;

	hEventQA->Fill(3.);
	const TVector3 Vertex3D = mEvent->primaryVertex();
	const double VertexX = Vertex3D.x();
	const double VertexY = Vertex3D.y();
	const double VertexZ = Vertex3D.z();
	const double VertexR = sqrt(VertexX * VertexX + VertexY * VertexY);
	const double vpdVz = mEvent->vzVpd();

	// event cut
	// if(refMult <=2 || refMult > 1000) return kStOK;
	mRefMultCorr->init(runID);
	mRefMultCorr->initEvent(refMult, VertexZ, mEvent->ZDCx());
	if (removeBadID(runID))
		return kStOK;
	if (mRefMultCorr->isBadRun(runID))
		return kStOK; // reject bad run of StRefMultCorr
	if (!mRefMultCorr->passnTofMatchRefmultCut(1. * refMult, 1. * tofMatch))
		return kStOK; // reject pileup of StRefMultCorr

	hEventQA->Fill(4.);
	hVertex2D->Fill(VertexZ, vpdVz);
	hDiffVz->Fill(VertexZ - vpdVz);

	const double DVz = VertexZ - vpdVz;

	if (fabs(VertexZ) > 70)
		return kStOK;
	// if (VertexZ > 0 || VertexZ < -70)
	// 	return kStOK; // for negative vz only
	if (sys_tag == 1 && fabs(VertexZ) > 35) // first systematic
		return kStOK;
	if (sqrt(pow(VertexX, 2.) + pow(VertexY, 2.)) > 2.0)
		return kStOK;
	// if(fabs(VertexZ-vpdVz)>4.) return kStOK;       // no vpd cut in low energy? Yes!!

	// check run number
	hEventQA->Fill(5.);
	int runnumberPointer = -999;
	runnumberPointer = CheckrunNumber(runID);
	if (runnumberPointer == -999)
		return kStOK;

	int dayPointer = (int)((runID) / 1000 % 1000);
	int yearPointer = (int)((runID) / 1000000);
	int mRunL = runList.at(0);
	int mDayL = (int)((mRunL) / 1000 % 1000);
	int mYearL = (int)((mRunL) / 1000000);
	yearPointer -= mYearL;
	dayPointer = yearPointer == 0 ? dayPointer - mDayL : dayPointer + 365 - mDayL;
	int timePointer = dayPointer / mStps;

	// StRefMultCorr
	double refmultWght = mRefMultCorr->getWeight();
	double refmultCorr = mRefMultCorr->getRefMultCorr();
	int centrality = mRefMultCorr->getCentralityBin9(); // 0 - 8  be careful !!!!!!!!
	// double mWght    = 1;
	// double mult_corr= refMult;
	// int centrality  = findCentrality(mult_corr)-1;  // 0-8
	if (centrality < 0 || centrality >= (nCent - 1))
		return kStOK;
	int cent = centrality + 1;

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
	// mEpdHits = mPicoDst->picoArray(8); // grab TClonesArray directly?
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
	int runID2 = runID % 1000000;
	int runID2bin = (runID2 - mRunStart) / 10 + 1;
	StEpdEpInfo result = mEpFinder->Results(mPicoDst->picoArray(StPicoArrays::EpdHit), Vertex3D, cent - 1);

	if (cent >= min_cent && cent <= max_cent)
	{
		hEPD_e_EP_1->Fill(result.EastPhiWeightedAndShiftedPsi(1));
		hEPD_w_EP_1->Fill(result.WestPhiWeightedAndShiftedPsi(1));
		hEPD_e_EP_2->Fill(result.EastPhiWeightedAndShiftedPsi(2));
		hEPD_w_EP_2->Fill(result.WestPhiWeightedAndShiftedPsi(2));
		hEPD_full_EP_1->Fill(result.FullPhiWeightedAndShiftedPsi(1));
		hEPD_full_EP_2->Fill(result.FullPhiWeightedAndShiftedPsi(2));
		for (int order = 1; order <= 3; order++)
			hEPD_ew_cos->Fill(order * 1.0, cos(order * 1.0 * (result.EastPhiWeightedAndShiftedPsi(order) - result.WestPhiWeightedAndShiftedPsi(order))));
	}
	float EP2_EPD_full = result.FullPhiWeightedAndShiftedPsi(2);
	hEPD_full_1_runID->Fill(runID, result.FullPhiWeightedAndShiftedPsi(1));
	hEPD_full_2_runID->Fill(runID, result.FullPhiWeightedAndShiftedPsi(2));
	hEPD_full_1_day->Fill(dayPointer, result.FullPhiWeightedAndShiftedPsi(1));
	hEPD_full_2_day->Fill(dayPointer, result.FullPhiWeightedAndShiftedPsi(2));

	// EPD QA
	if (cent == 4)
	{
		hEPDEP_2_Full_Raw_QA_cen4->Fill(result.FullRawPsi(2));
		hEPDEP_2_Full_PhiWeighted_QA_cen4->Fill(result.FullPhiWeightedPsi(2));
		hEPDEP_2_Full_PhiWeightedAndShifted_QA_cen4->Fill(result.FullPhiWeightedAndShiftedPsi(2));
	}

	///////////////////////////
	hNRefMult->Fill(grefMult);
	hNRefMultA->Fill(mult_corr);
	hNRefMultB->Fill(mult_corr, mWght);
	hVertexXY->Fill(VertexX, VertexY);
	hVertexZ->Fill(VertexZ);
	hcent->Fill(cent);		   // 1-9
	hcentw->Fill(cent, mWght); // 1-9

	///////////////
	hcentRefM->Fill(cent, mult_corr);
	hcentRefW->Fill(cent, mult_corr, mWght);
	hcentRefM->Fill(0., mult_corr);
	hcentRefW->Fill(0., mult_corr, mWght);
	hRefMultCorr_cent[cent - 1]->Fill(mult_corr);
	///////////////

	// ======= KFParticle ======= //
	vector<StLambdaDecayPair> KFParticleLambdaDecayPair;

	SetupKFParticle();
	if (InterfaceCantProcessEvent)
		return kStOK;

	// time(&time_now);
	// int time_diff = (int)difftime(time_now, time_start);
	// cout << "Up to finishing KFP: " << time_diff << endl;

	// collect lambdas and omegas
	std::vector<KFParticle> OmegaVec, OmegaSidebandVec, OmegaVecAll;
	std::vector<KFParticle> LambdaVec, LambdaSidebandVec, LambdaVecAll;
	bool hasLambda = false, hasLambdabar = false;
	std::vector<KFParticle> XiVecAll;
	std::vector<KFParticle> phiVecAll;
	std::vector<KFParticle> KsVecAll;
	std::vector<KFParticle> Xi1530VecAll;
	std::vector<int> OmegaDauPionIdVec;
	std::vector<int> OmegaDauProtonIdVec;
	mix_px.resize(0);
	mix_py.resize(0);
	mix_pz.resize(0);
	mix_charge.resize(0);
	for (int iKFParticle = 0; iKFParticle < KFParticlePerformanceInterface->GetNReconstructedParticles(); iKFParticle++)
	{
		const KFParticle particle = KFParticleInterface->GetParticles()[iKFParticle];
		bool IsOmega = false, IsLambda = false, IsXi = false, Isphi = false, IsKs = false, IsXi1530 = false;
		if (fabs(particle.GetPDG()) == OmegaPdg)
			IsOmega = true;
		else if (fabs(particle.GetPDG()) == LambdaPdg)
			IsLambda = true;
		else if (fabs(particle.GetPDG()) == XiPdg)
			IsXi = true;
		else if (fabs(particle.GetPDG()) == phiPdg)
			Isphi = true;
		else if (fabs(particle.GetPDG()) == KsPdg)
			IsKs = true;
		else if (fabs(particle.GetPDG()) == Xi1530Pdg)
			IsXi1530 = true;
		else
			continue;
		int upQ;
		if (particle.GetPDG() > 0)
			upQ = 1;
		else if (particle.GetPDG() < 0)
			upQ = -1;
		else
			continue;
		if (IsOmega)
			OmegaVecAll.push_back(particle);
		if (IsXi)
			XiVecAll.push_back(particle);
		if (IsLambda)
			LambdaVecAll.push_back(particle);
		if (Isphi)
			phiVecAll.push_back(particle);
		if (IsKs)
			KsVecAll.push_back(particle);
		if (IsXi1530)
			Xi1530VecAll.push_back(particle);
		if (IsOmega && isGoodOmega(cent, particle))
			OmegaVec.push_back(particle);
		if (IsLambda && isGoodLambda(cent, particle))
			LambdaVec.push_back(particle);
		if (IsLambda && isGoodLambda(cent, particle) && upQ == 1)
			hasLambda = true;
		if (IsLambda && isGoodLambda(cent, particle) && upQ == -1)
			hasLambdabar = true;
		if (IsOmega && isSidebandOmega(cent, particle))
			OmegaSidebandVec.push_back(particle);
		if (IsLambda && isSidebandLambda(cent, particle))
			LambdaSidebandVec.push_back(particle);

		// // sec tracks first/last points check
		// if (IsLambda && isGoodLambda(cent, particle)) SetDaughterTrackPointers(iKFParticle);

		// Omega/lambda QA
		if (IsOmega)
		{
			if (upQ == 1)
			{
				if (cent < min_cent || cent > max_cent) hOmegaM->Fill(particle.GetMass());

				// v2 vs centrality
				float phi_diff = fabs(particle.GetPhi() + pi - EP2_EPD_full);
				if (phi_diff > pi)
					phi_diff = phi_diff - pi;
				if (phi_diff > pi / 2.)
					phi_diff = pi - phi_diff;
				int phi_bin = static_cast<int>(floor(phi_diff / (pi / 2. / num_phi_bin)));
				// hOmegaM_phi[cent-1][phi_bin]->Fill(particle.GetMass());

				if (!isGoodOmega(cent, particle))
					continue; // subject to change
				hOmegay->Fill(particle.GetRapidity());
				hOmegaypt->Fill(particle.GetPt(), particle.GetRapidity());
				// hOmegaDL ->Fill(particle.GetDecayLength());

				// helix
				TVector3 pOmega(particle.GetPx(), particle.GetPy(), particle.GetPz());
				TVector3 xOmega(particle.GetX(), particle.GetY(), particle.GetZ());
				TLorentzVector OmegaLorentz(pOmega, particle.GetE());
				StPicoPhysicalHelix helixOmega(pOmega, xOmega, magnet * kilogauss, particle.GetQ());
				double pathlength = helixOmega.pathLength(Vertex3D, false);
				TVector3 pOmega_tb = helixOmega.momentumAt(pathlength, magnet * kilogauss);
				hOmegap->Fill(pOmega_tb.Mag());
				hOmegapt->Fill(pOmega_tb.Pt());
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

				hOmegaDL->Fill((xOmega - Vertex3D).Mag());

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
				// 	TVector3 pDaughterPico = daughterTrack->pMom();
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
				if (cent < min_cent || cent > max_cent) hOmegabarM->Fill(particle.GetMass());

				// v2 vs centrality
				float phi_diff = fabs(particle.GetPhi() + pi - EP2_EPD_full);
				if (phi_diff > pi)
					phi_diff = phi_diff - pi;
				if (phi_diff > pi / 2.)
					phi_diff = pi - phi_diff;
				int phi_bin = static_cast<int>(floor(phi_diff / (pi / 2. / num_phi_bin)));
				// hOmegabarM_phi[cent-1][phi_bin]->Fill(particle.GetMass());

				if (!isGoodOmega(cent, particle))
					continue; // subject to change
				hOmegabary->Fill(particle.GetRapidity());
				hOmegabarypt->Fill(particle.GetPt(), particle.GetRapidity());
				// hOmegabarDL ->Fill(particle.GetDecayLength());

				// helix
				TVector3 pOmega(particle.GetPx(), particle.GetPy(), particle.GetPz());
				TVector3 xOmega(particle.GetX(), particle.GetY(), particle.GetZ());
				TLorentzVector OmegaLorentz(pOmega, particle.GetE());
				StPicoPhysicalHelix helixOmega(pOmega, xOmega, magnet * kilogauss, particle.GetQ());
				double pathlength = helixOmega.pathLength(Vertex3D, false);
				TVector3 pOmega_tb = helixOmega.momentumAt(pathlength, magnet * kilogauss);
				hOmegabarp->Fill(pOmega_tb.Mag());
				hOmegabarpt->Fill(pOmega_tb.Pt());
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
				// 	TVector3 pDaughterPico = daughterTrack->pMom();
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
				if (cent < min_cent || cent > max_cent) hLambdaM->Fill(particle.GetMass());
				hLambdap->Fill(particle.GetMomentum());
				hLambdapt->Fill(particle.GetPt());
				hLambday->Fill(particle.GetRapidity());
				hLambdaphi->Fill(particle.GetPhi());
				hLambdaDL->Fill(particle.GetDecayLength());
			}
			else
			{
				if (cent < min_cent || cent > max_cent) hLambdabarM->Fill(particle.GetMass());
				hLambdabarp->Fill(particle.GetMomentum());
				hLambdabarpt->Fill(particle.GetPt());
				hLambdabary->Fill(particle.GetRapidity());
				hLambdabarphi->Fill(particle.GetPhi());
				hLambdabarDL->Fill(particle.GetDecayLength());
			}
		}
		else if (IsKs)
			hKsM_cen[cent - 1]->Fill(particle.GetMass());
		else if (IsXi1530)
		{
			if (upQ == 1)
				hXi1530M_cen[cent - 1]->Fill(particle.GetMass());
			else
				hXi1530barM_cen[cent - 1]->Fill(particle.GetMass());
		}
	} // End loop over KFParticles

	// fill Omega inv mass for with Lambda and without Lambda
	for (int iomega = 0; iomega < OmegaVecAll.size(); iomega++)
	{
		const KFParticle particle = OmegaVecAll[iomega];
		if (hasLambdabar)
		{
			if (particle.GetQ() < 0)
				hOmegaM_wlb->Fill(particle.GetMass());
			else
				hOmegabarM_wlb->Fill(particle.GetMass());
		}
		else
		{
			if (particle.GetQ() < 0)
				hOmegaM_wolb->Fill(particle.GetMass());
			else
				hOmegabarM_wolb->Fill(particle.GetMass());
		}
		if (hasLambda)
		{
			if (particle.GetQ() < 0)
				hOmegaM_wl->Fill(particle.GetMass());
			else
				hOmegabarM_wl->Fill(particle.GetMass());
		}
		else
		{
			if (particle.GetQ() < 0)
				hOmegaM_wol->Fill(particle.GetMass());
			else
				hOmegabarM_wol->Fill(particle.GetMass());
		}
	}

	// try reconstructing Omega(2012)
	for (int iXi = 0; iXi < XiVecAll.size(); iXi++)
	{
		for (int iKs = 0; iKs < KsVecAll.size(); iKs++)
		{
			const KFParticle Xi = XiVecAll[iXi];
			const KFParticle Ks = KsVecAll[iKs];
			KFParticle Omega2012;
			Omega2012 += Xi;
			Omega2012 += Ks;
			// if (Omega2012.GetChi2()/Omega2012.GetNDF() > 30) continue;
			if (Xi.GetPDG() > 0)
				hOmega2012M_cen[cent - 1]->Fill(Omega2012.GetMass());
			else
				hOmega2012barM_cen[cent - 1]->Fill(Omega2012.GetMass());
		}
	}
	hNumOmega->Fill(OmegaVec.size());

	// correlation function loop
	if (!PerformAnalysis)
		return kStOK;
	my_event current_event(OmegaVec);
	Efficiency eff_finder;
	Int_t nTracks = mPicoDst->numberOfTracks();
	std::vector<int> kaon_tracks;
	kaon_tracks.resize(0);
	std::vector<int> pion_tracks;
	pion_tracks.resize(0);
	std::vector<int> proton_tracks;
	proton_tracks.resize(0);
	float Qx2w = 0, Qy2w = 0, Qx2e = 0, Qy2e = 0; // for EP
	std::vector<int> track_index;
	int pct = 0;
	int pbct = 0; // counting protons
	int kpct = 0;
	int kmct = 0; // counting kaons
	int mult_east = 0;
	int mult_west = 0;
	bool fill_nomega = false;
	bool fill_nomega_sideband = false;
	for (Int_t iTrack = 0; iTrack < nTracks; iTrack++)
	{
		StPicoTrack *track = mPicoDst->track(iTrack);
		if (!track)
			continue;
		if (!track->charge())
			continue;
		if (track->nHitsFit() < 15)
			continue;
		if (track->nHitsDedx() < 15)
			continue;
		if (track->nHitsFit() * 1.0 / track->nHitsMax() < 0.52 || track->nHitsFit() * 1.0 / track->nHitsMax() > 1.05)
			continue;
		// if (  track->dEdxError() < 0.04 || track->dEdxError() > 0.12) continue; // same as kfp
		if (!track->isPrimary())
			continue;
		track_index.push_back(iTrack);

		// track info
		float p = track->pMom().Mag();
		float pt = track->pMom().Perp();
		float phi = track->pMom().Phi();
		float eta = track->pMom().Eta();
		float dcatopv = track->gDCA(Vertex3D).Mag();
		float nSigmaKaon = track->nSigmaKaon();
		float nSigmaPion = track->nSigmaPion();
		float nSigmaProton = track->nSigmaProton();
		if (hTOFEff[8] != 0)
			hTOFEff_check->Fill(pt, hTOFEff[8]->GetEfficiency(hTOFEff[8]->FindFixBin(pt)));

		// checking east west asymmetry (East=0 (y<0), West=1 (y>0), Pos=0, Neg=1)
		int eastwest = eta > 0 ? 1 : 0;
		int posneg = track->charge() > 0 ? 0 : 1;
		int eastwestvertex = VertexZ > 0 ? 1 : 0;
		hPhi_EW_PN_EWVtx[eastwest][posneg][eastwestvertex]->Fill(phi, mWght);
		for (int idXi = 0; idXi < XiVecAll.size(); idXi++)
		{
			const KFParticle particle = XiVecAll[idXi];
			// check if mass is within 3 sigma
			if (fabs(particle.GetMass() - XiPdgMass) > 3 * OmegaMassSigma) // can just use omega sigma for rough estimation
				continue;
			eastwest = particle.GetEta() > 0 ? 1 : 0;
			posneg = particle.GetQ() > 0 ? 0 : 1;
			hPhiXi_EW_PN_EWVtx[eastwest][posneg][eastwestvertex]->Fill(particle.GetPhi(), mWght);
		}
		if (eta < 0) mult_east++;
		else mult_west++;

		// EP
		bool isGoodAsso = true;
		// isGoodAsso &= track->isPrimary();
		isGoodAsso &= pT_asso_lo < pt && pt < pT_asso_hi;
		if (isGoodAsso)
		{
			// fill phi shift for asso
			if (fabs(eta) > TPCAssoEtaCut)
			{
				for (int i = 1; i <= shift_order_asso; i++) // fill correction for output
				{
					hTPCAssoShiftOutput_cos->Fill(dayPointer, i, cent - 1, cos(i * 1.0 * phi), mWght);
					hTPCAssoShiftOutput_sin->Fill(dayPointer, i, cent - 1, sin(i * 1.0 * phi), mWght);
				}
				// shift asso phi
				float phi_shifted_asso = ShiftAssoPhi(phi, dayPointer + 1, cent);
				hTPCAssoPhi->Fill(phi, mWght);
				hTPCAssoPhi_shifted->Fill(phi_shifted_asso, mWght);
				hTPCAssoPhi_2D->Fill(dayPointer, phi, mWght);
				hTPCAssoPhi_2D_shifted->Fill(dayPointer, phi_shifted_asso, mWght);
				if (eta > TPCAssoEtaCut) // construct east EP
				{
					Qx2e += pt * cos(2 * phi_shifted_asso);
					Qy2e += pt * sin(2 * phi_shifted_asso);
				}
				else if (eta < -1.0 * TPCAssoEtaCut) // construct west EP
				{
					Qx2w += pt * cos(2 * phi_shifted_asso);
					Qy2w += pt * sin(2 * phi_shifted_asso);
				}
			}
		}
		bool isGoodPOI = true;
		// isGoodPOI &= track->isPrimary();
		isGoodPOI &= pT_trig_lo < pt && pt < pT_trig_hi;
		isGoodPOI &= fabs(eta) < eta_trig_cut;
		if (isGoodPOI)
		{
			// fill phi shift for POI
			for (int i = 1; i <= shift_order_asso; i++) // fill correction for output
			{
				hTPCPOIShiftOutput_cos->Fill(dayPointer, i, cent - 1, cos(i * 1.0 * phi), mWght);
				hTPCPOIShiftOutput_sin->Fill(dayPointer, i, cent - 1, sin(i * 1.0 * phi), mWght);
			}
			// shift POI phi
			float phi_shifted_POI = ShiftPOIPhi(phi, dayPointer + 1, cent);
			hTPCPOIPhi->Fill(phi, mWght);
			hTPCPOIPhi_shifted->Fill(phi_shifted_POI, mWght);
			hTPCPOIPhi_2D->Fill(dayPointer, phi, mWght);
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
			if ((tofflag >= 1) && (tof > 0) && (BtofYLocal > -1.8) && (BtofYLocal < 1.8))
				hasTOF = true;
		}

		// fill QA
		StPicoPhysicalHelix helix = track->helix(magnet);
		TVector3 pkaon = helix.momentum(magnet * kilogauss);
		hgpdEdx->Fill(pkaon.Mag(), track->dEdx());
		hgdEdxErr->Fill(track->dEdxError());
		// hgp       ->Fill(p);
		hgpTeta[cent - 1]->Fill(pt, eta);
		if (hasTOF)
			hgpTeta_TOF[cent - 1]->Fill(pt, eta);
		if (fabs(eta) < y_coal_cut)
		{
			hgp[cent - 1]->Fill(p);
			hgpT[cent - 1]->Fill(pt);
			if (hasTOF)
			{
				hgp_TOF[cent - 1]->Fill(p);
				hgpT_TOF[cent - 1]->Fill(pt);
			}
		}
		hgDCAtoPV->Fill(dcatopv);

		int ptbin = static_cast<int>(floor(pt / 0.2));
		int pbin = static_cast<int>(floor(p / 0.2));
		double zTPC = TMath::Log(track->dEdx() / 1e6 / StdEdxPull::EvalPred(pkaon.Mag() / KaonPdgMass, 1, 1));
		hgnSigmaDiff->Fill(zTPC / track->dEdxError() - nSigmaKaon);
		hgptnSigmaKaon->Fill(pt, nSigmaKaon);
		hgptnSigmaPion->Fill(pt, nSigmaPion);
		hgptnSigmaProton->Fill(pt, nSigmaProton);

		// if (ptbin >= 0 && ptbin <= 14) hgzTPC_pt[ptbin]->Fill(track->nSigmaKaon());
		if (hasTOF)
		{
			beta = (mPicoDst->btofPidTraits(tofindex))->btofBeta();
			m2 = pkaon.Mag2() * (1.0 / beta / beta - 1.0);
			hgpinvbeta->Fill(pkaon.Mag(), 1. / beta);
			hgm2->Fill(m2);
			hgpm2->Fill(pkaon.Mag(), m2);
			hgptm2->Fill(pt, m2);
			hgm2nSigmaKaon->Fill(nSigmaKaon, m2);
			hgm2nSigmaPion->Fill(nSigmaPion, m2);
			hgm2nSigmaProton->Fill(nSigmaProton, m2);

			// some kaon QA
			// if (track->nSigmaKaon() >  6) hgptm2_largenSigmaKaon->Fill(track->pMom().Perp(), m2);
			// if (track->nSigmaKaon() < -6) hgptm2_smallnSigmaKaon->Fill(track->pMom().Perp(), m2);
			double zTOF_proton = 1 / beta - sqrt(ProtonPdgMass * ProtonPdgMass / pkaon.Mag2() + 1);
			double zTOF_pion = 1 / beta - sqrt(PionPdgMass * PionPdgMass / pkaon.Mag2() + 1);
			double zTOF_kaon = 1 / beta - sqrt(KaonPdgMass * KaonPdgMass / pkaon.Mag2() + 1);

#if !USE_P
			// using pT
			// if (ptbin >= 0 && ptbin <= 9 && track->charge() > 0) hgPID2D_proton_pt    [ptbin]->Fill(nSigmaProton, m2); // zTOF_proton);
			// if (ptbin >= 0 && ptbin <= 9 && track->charge() < 0) hgPID2D_antiproton_pt[ptbin]->Fill(nSigmaProton, m2); // zTOF_proton);
			// if (ptbin >= 0 && ptbin <= 9 && track->charge() > 0) hgPID2D_piplus_pt    [ptbin]->Fill(nSigmaPion  , m2); // zTOF_pion  );
			// if (ptbin >= 0 && ptbin <= 9 && track->charge() < 0) hgPID2D_piminus_pt   [ptbin]->Fill(nSigmaPion  , m2); // zTOF_pion  );
			// if (ptbin >= 0 && ptbin <= 9 && track->charge() > 0) hgPID2D_kplus_pt     [ptbin]->Fill(nSigmaKaon  , m2); // zTOF_kaon  );
			// if (ptbin >= 0 && ptbin <= 9 && track->charge() < 0) hgPID2D_kminus_pt    [ptbin]->Fill(nSigmaKaon  , m2); // zTOF_kaon  );

			// using p
			// if (pbin >= 0 && pbin <= 11 && track->charge() > 0)
			// 	hgPID2D_proton_p[pbin]->Fill(nSigmaProton, m2);
			// if (pbin >= 0 && pbin <= 11 && track->charge() < 0)
			// 	hgPID2D_antiproton_p[pbin]->Fill(nSigmaProton, m2);
			// if (pbin >= 0 && pbin <= 11 && track->charge() > 0)
			// 	hgPID2D_piplus_p[pbin]->Fill(nSigmaPion, m2);
			// if (pbin >= 0 && pbin <= 11 && track->charge() < 0)
			// 	hgPID2D_piminus_p[pbin]->Fill(nSigmaPion, m2);
			// if (pbin >= 0 && pbin <= 11 && track->charge() > 0)
			// 	hgPID2D_kplus_p[pbin]->Fill(nSigmaKaon, m2);
			// if (pbin >= 0 && pbin <= 11 && track->charge() < 0)
			// 	hgPID2D_kminus_p[pbin]->Fill(nSigmaKaon, m2);
#else
			// using p
			if (pbin >= 0 && pbin <= 11 && track->charge() > 0)
				hgPID2D_proton_p[pbin]->Fill(nSigmaProton, m2);
			if (pbin >= 0 && pbin <= 11 && track->charge() < 0)
				hgPID2D_antiproton_p[pbin]->Fill(nSigmaProton, m2);
			if (pbin >= 0 && pbin <= 11 && track->charge() > 0)
				hgPID2D_piplus_p[pbin]->Fill(nSigmaPion, m2);
			if (pbin >= 0 && pbin <= 11 && track->charge() < 0)
				hgPID2D_piminus_p[pbin]->Fill(nSigmaPion, m2);
			if (pbin >= 0 && pbin <= 11 && track->charge() > 0)
				hgPID2D_kplus_p[pbin]->Fill(nSigmaKaon, m2);
			if (pbin >= 0 && pbin <= 11 && track->charge() < 0)
				hgPID2D_kminus_p[pbin]->Fill(nSigmaKaon, m2);
#endif
		}

		// change pt to p
		if (USE_P)
			pt = p;

		// primary proton cut for coalescence test
		bool proton_cut = true;
		if (pt < pT_trig_lo || pt > pT_trig_hi)
			proton_cut = false;
		if (fabs(eta) > eta_trig_cut)
			proton_cut = false;
		ProtonPID proton_pid(0., nSigmaProton, pt, p, UsePtForPID);				   // not using zTOF
		if ((pt > proton_pT_lo && pt < proton_pT_TOFth) && hasTOF) // test efficacy of ProtonPID.h
		{
			hm2proton_b->Fill(m2);
			if (proton_pid.IsProtonSimple(2., track->charge()))
				hm2proton_a->Fill(m2);
			if (fabs(nSigmaProton) < 2.)
				hm2proton_r->Fill(m2);
		}
		if (!hasTOF && pt > proton_pT_TOFth)
			proton_cut = false;
		if (pt > proton_pT_TOFth && (m2 > proton_m2_hi || m2 < proton_m2_lo))
			proton_cut = false;
		if (!proton_pid.IsProtonSimple(2., track->charge()))
			proton_cut = false; // only 0.2 < pt < 2.0!!!
		// if (fabs(nSigmaProton) > 3) proton_cut = false;
		if (dcatopv > dcatoPV_hi)
			proton_cut = false;
		if (proton_cut)
		{
			TLorentzVector lv_proton;
			lv_proton.SetVectM(track->pMom(), ProtonPdgMass);
			if (track->charge() > 0)
			{
				hProtony->Fill(lv_proton.Rapidity());
				if (fabs(lv_proton.Rapidity()) < 0.5)
					pct++;
			}
			else
			{
				hAntiProtony->Fill(lv_proton.Rapidity());
				if (fabs(lv_proton.Rapidity()) < 0.5)
					pbct++;
			}

			proton_tracks.push_back(iTrack);
		}

		// primary pion cut for coalescence test
		bool pion_cut = true;
		if (pt < pT_trig_lo || pt > pT_trig_hi)
			pion_cut = false; // use p < 2
		if (fabs(eta) > eta_trig_cut)
			pion_cut = false;
		PionPID pion_pid(0., nSigmaPion, pt, p, UsePtForPID);				   // not using zTOF
		if ((pt > pion_pT_lo && pt < pion_pT_TOFth) && hasTOF) // test efficacy of ProtonPID.h
		{
			hm2pion_b->Fill(m2);
			if (pion_pid.IsPionSimple(2., track->charge()))
				hm2pion_a->Fill(m2);
			if (fabs(nSigmaPion) < 2.)
				hm2pion_r->Fill(m2);
		}
		if (!hasTOF && pt > pion_pT_TOFth)
			pion_cut = false;
		if (pt > pion_pT_TOFth && (m2 > pion_m2_hi || m2 < pion_m2_lo))
			pion_cut = false;
		if (!pion_pid.IsPionSimple(2., track->charge()))
			pion_cut = false; // only up to pt < 1.8!!!
		// if (fabs(nSigmaPion) > 3) pion_cut = false;
		if (dcatopv > dcatoPV_hi)
			pion_cut = false;
		if (pion_cut)
			pion_tracks.push_back(iTrack);

		// change pt back to pt
		if (USE_P)
			pt = track->pMom().Perp();

		// primary kaon cut
		/******** looser cut ********/
		bool kaon_cut = true;
		if (pt < pT_trig_lo || pt > pT_trig_hi)
			kaon_cut = false; // use p < 1.6
		if (fabs(eta) > eta_trig_cut)
			kaon_cut = false;
		if (!hasTOF && pt > kaon_pT_TOFth)
			kaon_cut = false;
		if (pt > kaon_pT_TOFth && (m2 > kaon_m2_hi || m2 < kaon_m2_lo))
			kaon_cut = false;
		double zTOF = 1 / beta - sqrt(KaonPdgMass * KaonPdgMass / pkaon.Mag2() + 1);

		KaonPID decider(zTOF, nSigmaKaon, pt, p, UsePtForPID);
		if (!decider.IsKaonSimple(2., track->charge()))
			kaon_cut = false;
		if (dcatopv > 2.)
			kaon_cut = false;
		bool isDaughter = false;
		for (int i = 0; i < OmegaVec.size(); i++)
			if (IsKaonOmegaDaughter(OmegaVec[i], track->id()))
				isDaughter = true;
		if (isDaughter)
			kaon_cut = false;
		if (kaon_cut)
		{
			if (track->charge() > 0)
				kpct++;
			else
				kmct++;
			kaon_tracks.push_back(iTrack);
		}

		/******** stricter cut ********/
		/*
		if (pt < 0.15 || pt> 1.6) continue; // use p < 1.6
		if (!hasTOF) continue;pt
		double zTOF = 1/beta - sqrt(KaonPdgMass*KaonPdgMass/pkaon.Mag2()+1);
		KaonPID decider(zTOF, track->nSigmaKaon(), track->pMom().Perp());
		if (!decider.IsKaon()) continue;
		*/

		// kaon QA and correlation
		if (kaon_cut)
		{
			float TOFEff = 1.0;
			float TOFEff2D = 1.0;
			float TPCEff = track->charge() > 0 ? eff_finder.GetEfficiency1D(pt, cent, "Kp") : eff_finder.GetEfficiency1D(pt, cent, "Km");
			float TPCEff2D = track->charge() > 0 ? eff_finder.GetEfficiency2D(pt, eta, cent, "Kp") : eff_finder.GetEfficiency2D(pt, eta, cent, "Km");
			if (hTOFEff[cent - 1] != 0 && pt > 0.4)
				TOFEff = hTOFEff[cent - 1]->GetEfficiency(hTOFEff[cent - 1]->FindFixBin(pt));
			if (hTOFEff_2D[cent - 1] != 0 && pt > 0.4)
				TOFEff2D = hTOFEff_2D[cent - 1]->GetBinContent(hTOFEff_2D[cent - 1]->FindFixBin(pt, eta));
			hgKpdEdx->Fill(pkaon.Mag(), track->dEdx());
			hgKp->Fill(p);
			hgKpT->Fill(pt);
			hgKDCAtoPV->Fill(dcatopv);
			if (hasTOF)
			{
				float beta = (mPicoDst->btofPidTraits(tofindex))->btofBeta();
				hgKpinvbeta->Fill(pkaon.Mag(), 1. / beta);
				hgKm2->Fill(m2);
				hgKpm2->Fill(pkaon.Mag(), m2);

				// check kaons misidentified as pions
				if (m2 > -0.06 && m2 < 0.1)
					hgKpionpdEdx->Fill(pkaon.Mag(), track->dEdx());
			}
		}
	}
	hMultRatioEW_VertexZ->Fill(VertexZ, mult_east * 1.0 / mult_west);

	// // TPC resolution shuffling
	// float Qx1_e = 0, Qy1_e = 0, Qx2_e = 0, Qy2_e = 0; // for EP res
	// float Qx1_w = 0, Qy1_w = 0, Qx2_w = 0, Qy2_w = 0; // for EP res
	// TVector2 Q1_e, Q2_e, Q1_w, Q2_w; // for EP res
	// std::random_shuffle(track_index.begin(), track_index.end());
	// for (int i = 0; i < track_index.size(); i++)
	// {
	// 	StPicoTrack *track = mPicoDst->track(track_index[i]);
	// 	float pt = track->pMom().Perp();
	// 	float phi = track->pMom().Phi();

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
	bool mult_index_bad_flag = false;
	bool vz_index_bad_flag = false;
	bool EP_index_bad_flag = false;
	int mult_index = MixRefMultBin(cent, mult_corr);
	if (mult_index < 0 || mult_index >= num_mult_bin)
		mult_index_bad_flag = true;
	float vz_cuts[] = {-70., -60., -50., -40., -30., -20., -10., 0., 10., 20., 30., 40., 50., 60., 70.};
	int vz_index = -1;
	for (int i = 0; i < num_vz_bin; i++)
	{
		if (VertexZ >= vz_cuts[i] && VertexZ < vz_cuts[i + 1])
			vz_index = i;
	}
	if (vz_index < 0 || vz_index >= num_vz_bin)
		vz_index_bad_flag = true;
	int EP_index = -1;

	// TPC EP
	float EP2_TPC_w = -999, EP2_TPC_e = -999, EP2_TPC_full = -999;
	bool goodTPCEP = false;
	if (!(Qx2e == 0. || Qx2w == 0. || Qy2e == 0. || Qy2w == 0.))
	{
		goodTPCEP = true;
		TVector2 Q2w, Q2e, Q2full; // for EP
		Q2full.Set(Qx2e + Qx2w, Qy2e + Qy2w);
		Q2e.Set(Qx2e, Qy2e);
		Q2w.Set(Qx2w, Qy2w);
		EP2_TPC_full = Q2full.Phi() / 2.;
		EP2_TPC_e = Q2e.Phi() / 2.;
		EP2_TPC_w = Q2w.Phi() / 2.;
	}

	// shift correction
	float EP2_TPC[3] = {EP2_TPC_e, EP2_TPC_w, EP2_TPC_full};
	float EP2_TPC_shifted[3] = {EP2_TPC_e, EP2_TPC_w, EP2_TPC_full};
	for (int ewFull = 0; ewFull < 3; ewFull++)
	{
		for (int i = 1; i <= shift_order_EP; i++)
		{
			if (!goodTPCEP)
				continue;
			hTPCEPShiftOutput_cos[ewFull]->Fill(dayPointer, i, cent - 1, cos(i * 2. * EP2_TPC[ewFull]), mWght);
			hTPCEPShiftOutput_sin[ewFull]->Fill(dayPointer, i, cent - 1, sin(i * 2. * EP2_TPC[ewFull]), mWght);
		}
	}

	// float EP_1_shift_cos[4] = {-0.191905, 0.00567376, 0.003376, -0.000779899};
	// float EP_1_shift_sin[4] = {0.170762, -0.0436447, 0.00356307, -0.000438169};
	// float EP_2_shift_cos[4] = {-0.0447853, -0.0170472, 0.00125928, 0.000215389};
	// float EP_2_shift_sin[4] = {0.164891, -0.0109107, -0.000975253, 0.000114769};
	EP_index = static_cast<int>(EP2_TPC_full / (PI / num_EP_bin * 1.0));
	if (EP_index == num_EP_bin)
		EP_index = num_EP_bin - 1; // for event mixing index
	if (!goodTPCEP)
		EP_index = -1;
	if (EP_index < 0 || EP_index >= num_EP_bin)
		EP_index_bad_flag = true;

	// shift EPs
	float EP2_TPC_e_shifted = EP2_TPC_shifted[0];
	float EP2_TPC_w_shifted = EP2_TPC_shifted[1];
	float EP2_TPC_full_shifted = EP2_TPC_shifted[2];
	if (goodTPCEP)
	{
		for (int ewFull = 0; ewFull < 3; ewFull++)
		{
			EP2_TPC_shifted[ewFull] = ShiftTPCEP(EP2_TPC[ewFull], dayPointer + 1, cent, ewFull, 2);
			float shift = EP2_TPC_shifted[ewFull] - EP2_TPC[ewFull];

			hTPCEP_2[cent - 1][ewFull]->Fill(EP2_TPC[ewFull]);
			hTPCEP_2_shifted[cent - 1][ewFull]->Fill(EP2_TPC_shifted[ewFull]);
			if (CheckWeights2D)
			{
				hTPCEP_2_2D[cent - 1][ewFull]->Fill(dayPointer, EP2_TPC[ewFull]);
				hTPCEP_2_2D_shifted[cent - 1][ewFull]->Fill(dayPointer, EP2_TPC_shifted[ewFull]);
			}
			hTPCEP_2_shift[cent - 1][ewFull]->Fill(EP2_TPC[ewFull], shift);
		}
		EP2_TPC_e_shifted = EP2_TPC_shifted[0];
		EP2_TPC_w_shifted = EP2_TPC_shifted[1];
		EP2_TPC_full_shifted = EP2_TPC_shifted[2];
		hTPCEP_ew_cos->Fill(cent, cos(2. * EP2_TPC_e_shifted - 2. * EP2_TPC_w_shifted), mWght); // resolution
	}

	// EPD user shift
	bool goodEPDEP = true;
	// DEBUGGING EPD
	if (result.EastRawPsi(1) > 0.026463 && result.EastRawPsi(1) < 0.026464) goodEPDEP = false;
	if (result.WestRawPsi(1) > 0.026463 && result.WestRawPsi(1) < 0.026464) goodEPDEP = false;
	if (result.EastRawPsi(2) > 0.026463 && result.EastRawPsi(2) < 0.026464) goodEPDEP = false;
	if (result.WestRawPsi(2) > 0.026463 && result.WestRawPsi(2) < 0.026464) goodEPDEP = false;
	float EP1_EPD[3] = {result.EastPhiWeightedPsi(1), result.WestPhiWeightedPsi(1), result.FullPhiWeightedPsi(1)};		   // EPD EP1
	float EP2_EPD[3] = {result.EastPhiWeightedPsi(2), result.WestPhiWeightedPsi(2), result.FullPhiWeightedPsi(2)};		   // EPD EP2
	float EP1_EPD_shifted[3] = {result.EastPhiWeightedPsi(1), result.WestPhiWeightedPsi(1), result.FullPhiWeightedPsi(1)}; // EPD EP1 initialize
	float EP2_EPD_shifted[3] = {result.EastPhiWeightedPsi(2), result.WestPhiWeightedPsi(2), result.FullPhiWeightedPsi(2)}; // EPD EP2 initialize
	float EP1_EPD_e_shifted = EP1_EPD_shifted[0];
	float EP1_EPD_w_shifted = EP1_EPD_shifted[1];
	float EP1_EPD_full_shifted = EP1_EPD_shifted[2];
	float EP2_EPD_e_shifted = EP2_EPD_shifted[0];
	float EP2_EPD_w_shifted = EP2_EPD_shifted[1];
	float EP2_EPD_full_shifted = EP2_EPD_shifted[2];
	if (goodEPDEP)
	{
		for (int ewFull = 0; ewFull < 3; ewFull++)
		{
			for (int i = 1; i <= shift_order_EP; i++)
			{
				hEPDEPShiftOutput_1_cos[ewFull]->Fill(dayPointer, i, cent - 1, cos(i * 1. * EP1_EPD[ewFull]), mWght);
				hEPDEPShiftOutput_1_sin[ewFull]->Fill(dayPointer, i, cent - 1, sin(i * 1. * EP1_EPD[ewFull]), mWght);
				hEPDEPShiftOutput_2_cos[ewFull]->Fill(dayPointer, i, cent - 1, cos(i * 2. * EP2_EPD[ewFull]), mWght);
				hEPDEPShiftOutput_2_sin[ewFull]->Fill(dayPointer, i, cent - 1, sin(i * 2. * EP2_EPD[ewFull]), mWght);
			}
		}
		// shift
		for (int ewFull = 0; ewFull < 3; ewFull++)
		{
			EP1_EPD_shifted[ewFull] = ShiftEPDEP(EP1_EPD[ewFull], dayPointer + 1, cent, ewFull, 1);
			EP2_EPD_shifted[ewFull] = ShiftEPDEP(EP2_EPD[ewFull], dayPointer + 1, cent, ewFull, 2);
			float shift = EP2_EPD_shifted[ewFull] - EP2_EPD[ewFull];

			hEPDEP_1[cent - 1][ewFull]->Fill(EP1_EPD[ewFull]);
			hEPDEP_1_shifted[cent - 1][ewFull]->Fill(EP1_EPD_shifted[ewFull]);
			hEPDEP_2[cent - 1][ewFull]->Fill(EP2_EPD[ewFull]);
			hEPDEP_2_shifted[cent - 1][ewFull]->Fill(EP2_EPD_shifted[ewFull]);
			if (CheckWeights2D)
			{
				hEPDEP_1_2D[cent - 1][ewFull]->Fill(dayPointer, EP1_EPD[ewFull]);
				hEPDEP_1_2D_shifted[cent - 1][ewFull]->Fill(dayPointer, EP1_EPD_shifted[ewFull]);
				hEPDEP_2_2D[cent - 1][ewFull]->Fill(dayPointer, EP2_EPD[ewFull]);
				hEPDEP_2_2D_shifted[cent - 1][ewFull]->Fill(dayPointer, EP2_EPD_shifted[ewFull]);
			}
			hEPDEP_2_shift[cent - 1][ewFull]->Fill(EP2_EPD[ewFull], shift);
		}
		EP1_EPD_e_shifted = EP1_EPD_shifted[0];
		EP1_EPD_w_shifted = EP1_EPD_shifted[1];
		EP1_EPD_full_shifted = EP1_EPD_shifted[2];
		EP2_EPD_e_shifted = EP2_EPD_shifted[0];
		EP2_EPD_w_shifted = EP2_EPD_shifted[1];
		EP2_EPD_full_shifted = EP2_EPD_shifted[2];
		hEPDEP_ew_cos_1->Fill(cent, cos(1. * EP1_EPD_shifted[0] - 1. * EP1_EPD_shifted[1]), mWght); // resolution
		hEPDEP_ew_cos_2->Fill(cent, cos(2. * EP2_EPD_shifted[0] - 2. * EP2_EPD_shifted[1]), mWght); // resolution
		hEPDEP_ew_cos_3->Fill(cent, cos(3. * EP2_EPD_shifted[0] - 3. * EP2_EPD_shifted[1]), mWght); // resolution
	}
	// coalescence v2 and v1
	std::vector<std::vector<int>> IdPar_tracks;
	IdPar_tracks.push_back(pion_tracks);
	IdPar_tracks.push_back(proton_tracks);
	IdPar_tracks.push_back(kaon_tracks);
	int IdPar_nq[num_IdPar] = {2, 2, 3, 3, 2, 2};
	float IdPar_PdgMass[num_IdPar] = {PionPdgMass, PionPdgMass, ProtonPdgMass, ProtonPdgMass, KaonPdgMass, KaonPdgMass};
	float IdPar_pT_lo[num_IdPar] = {pion_pT_lo, pion_pT_lo, proton_pT_lo, proton_pT_lo, pion_pT_lo, pion_pT_lo}; // kaon uses same cut as pion
	float IdPar_pT_hi[num_IdPar] = {pion_pT_hi, pion_pT_hi, proton_pT_hi, proton_pT_hi, pion_pT_hi, pion_pT_hi};
	for (int idpar = 0; idpar < IdPar_tracks.size(); idpar++)
	{
		std::vector<int> tracks = IdPar_tracks[idpar];
		for (int i = 0; i < tracks.size(); i++)
		{
			StPicoTrack *track = mPicoDst->track(tracks[i]);
			if (!track)
				continue;
			int charge_index = track->charge() > 0 ? 0 : 1;
			float p = track->pMom().Mag();
			float pt = track->pMom().Perp();
			float eta = track->pMom().Eta();
			float phi = track->pMom().Phi();
			float phi_shifted_POI = ShiftPOIPhi(phi, dayPointer + 1, cent);
			float y = Eta2y(pt, eta, IdPar_PdgMass[2 * idpar + charge_index]);

			hIdPar_y_ptnq[2 * idpar + charge_index]->Fill(y, pt / IdPar_nq[2 * idpar + charge_index]);
			if (goodEPDEP)
			{
				if (fabs(y) <= y_coal_cut)
				{
					if (pt > IdPar_pT_lo[2 * idpar + charge_index] && pt < IdPar_pT_hi[2 * idpar + charge_index])
						hIdPar_EPD_v2[2 * idpar + charge_index]->Fill(cent, cos(2. * phi_shifted_POI - EP1_EPD_w_shifted - EP1_EPD_e_shifted));
					
					hIdPar_EPD_v2_pt_1st[2 * idpar + charge_index][cent - 1]->Fill(pt, cos(2. * phi_shifted_POI - EP1_EPD_w_shifted - EP1_EPD_e_shifted));
					hIdPar_EPD_v2_pt_2nd[2 * idpar + charge_index][cent - 1]->Fill(pt, cos(2. * phi_shifted_POI - 2. * EP2_EPD_w_shifted));
					hIdPar_EPD_v2_pt_2nd[2 * idpar + charge_index][cent - 1]->Fill(pt, cos(2. * phi_shifted_POI - 2. * EP2_EPD_e_shifted));

					if (pt > IdPar_pT_lo[2 * idpar + charge_index] && pt < IdPar_pT_hi[2 * idpar + charge_index])
						hIdPar_EPD_v1_y[2 * idpar + charge_index][cent - 1]->Fill(y, cos(phi_shifted_POI - EP1_EPD_w_shifted));
					if (pt > IdPar_pT_lo[2 * idpar + charge_index] && pt < IdPar_pT_hi[2 * idpar + charge_index])
						hIdPar_EPD_v1_y[2 * idpar + charge_index][cent - 1]->Fill(y, cos(phi_shifted_POI - EP1_EPD_e_shifted));
				}
				if (Coal2D)
				{
					hIdPar_EPD_v1_eta_pt_1st[2 * idpar + charge_index][cent - 1]->Fill(eta, pt, cos(phi_shifted_POI - EP1_EPD_w_shifted));
					hIdPar_EPD_v1_eta_pt_1st[2 * idpar + charge_index][cent - 1]->Fill(eta, pt, cos(phi_shifted_POI - EP1_EPD_e_shifted));
					hIdPar_EPD_v2_eta_pt_1st[2 * idpar + charge_index][cent - 1]->Fill(eta, pt, cos(2. * phi_shifted_POI - EP1_EPD_w_shifted - EP1_EPD_e_shifted));
					hIdPar_EPD_v3_eta_pt_1st[2 * idpar + charge_index][cent - 1]->Fill(eta, pt, cos(3. * phi_shifted_POI - 3. * EP1_EPD_w_shifted));
					hIdPar_EPD_v3_eta_pt_1st[2 * idpar + charge_index][cent - 1]->Fill(eta, pt, cos(3. * phi_shifted_POI - 3. * EP1_EPD_e_shifted));
				}
				// if (Coal2D)
				// {
				// 	hIdPar_EPD_v2_y_pt[2 * idpar + charge_index][cent - 1]->Fill(y, pt, cos(2. * phi_shifted_POI - 2. * EP2_EPD_w_shifted));
				// 	hIdPar_EPD_v2_y_pt[2 * idpar + charge_index][cent - 1]->Fill(y, pt, cos(2. * phi_shifted_POI - 2. * EP2_EPD_e_shifted));
				// }
			}
			if (goodTPCEP)
			{
				if (fabs(y) <= y_coal_cut)
				{
					if (eta > 0)
					{
						if (pt > IdPar_pT_lo[2 * idpar + charge_index] && pt < IdPar_pT_hi[2 * idpar + charge_index])
							hIdPar_TPC_v2[2 * idpar + charge_index]->Fill(cent, cos(2. * (phi_shifted_POI - EP2_TPC_w_shifted)));
						hIdPar_TPC_v2_pt[2 * idpar + charge_index][cent - 1]->Fill(pt, cos(2. * (phi_shifted_POI - EP2_TPC_w_shifted)));
						// if (Coal2D)
						// 	hIdPar_TPC_v2_y_pt[2 * idpar + charge_index][cent - 1]->Fill(y, pt, cos(2. * (phi_shifted_POI - EP2_TPC_w_shifted)));
					}
					else
					{
						if (pt > IdPar_pT_lo[2 * idpar + charge_index] && pt < IdPar_pT_hi[2 * idpar + charge_index])
							hIdPar_TPC_v2[2 * idpar + charge_index]->Fill(cent, cos(2. * (phi_shifted_POI - EP2_TPC_e_shifted)));
						hIdPar_TPC_v2_pt[2 * idpar + charge_index][cent - 1]->Fill(pt, cos(2. * (phi_shifted_POI - EP2_TPC_e_shifted)));
						// if (Coal2D)
						// 	hIdPar_TPC_v2_y_pt[2 * idpar + charge_index][cent - 1]->Fill(y, pt, cos(2. * (phi_shifted_POI - EP2_TPC_e_shifted)));
					}
				}
			}
		}
	}

	// if (OmegaVec.size() == 1 && OmegaVec[0].GetQ() < 0) // Omega event
	// {
	// 	hkplus_totaly->Fill(0., kplus_totaly);
	// 	hkplus_ntrack->Fill(0., kplus_ntrack);
	// 	hkminus_totaly->Fill(1., kminus_totaly);
	// 	hkminus_ntrack->Fill(1., kminus_ntrack);
	// }
	// if (OmegaVec.size() == 1 && OmegaVec[0].GetQ() > 0) // Omegabar event
	// {
	// 	hkplus_totaly->Fill(1., kplus_totaly);
	// 	hkplus_ntrack->Fill(1., kplus_ntrack);
	// 	hkminus_totaly->Fill(0., kminus_totaly);
	// 	hkminus_ntrack->Fill(0., kminus_ntrack);
	// }

	// RecoPar v2
	std::vector<std::vector<KFParticle>> RecoParVecAll;
	RecoParVecAll.push_back(OmegaVecAll);
	RecoParVecAll.push_back(XiVecAll);
	RecoParVecAll.push_back(LambdaVecAll);
	RecoParVecAll.push_back(phiVecAll);
	int RecoPar_nq[num_RecoPar] = {3, 3, 3, 3, 3, 3, 2};
	for (int recopar = 0; recopar < RecoParVecAll.size(); recopar++)
	{
		std::vector<KFParticle> ParVecAll = RecoParVecAll[recopar];
		for (int i = 0; i < ParVecAll.size(); i++)
		{
			const KFParticle particle = ParVecAll[i];
			float EP2_TPC_RecoPar = particle.GetEta() > 0 ? EP2_TPC_w_shifted : EP2_TPC_e_shifted;
			// float EP2_EPD_RecoPar = particle.GetEta() > 0 ? EP2_EPD_w_shifted : EP2_EPD_e_shifted;
			int charge_index = particle.GetPDG() > 0 ? 0 : 1;
			hRecoPar_y_ptnq[2 * recopar + charge_index]->Fill(particle.GetRapidity(), particle.GetPt() / RecoPar_nq[2 * recopar + charge_index]);

			// pt bin for lambda and lambdabar
			if (recopar == 2)
			{
				int ptbin = static_cast<int>(particle.GetPt() / 0.1);
				if (ptbin <= 4 + num_pt_bin - 1 && ptbin >=4 && fabs(particle.GetRapidity()) > y_coal_cut)
				{
					if (particle.GetPDG() > 0)
					{
						hLambdaM_pt[ptbin - 4][cent - 1]->Fill(particle.GetMass());
						if (goodTPCEP)
							hLambda_TPC_v2_pt[ptbin - 4][cent - 1]->Fill(particle.GetMass(), cos(2. * particle.GetPhi() - 2. * EP2_TPC_RecoPar));
						if (goodEPDEP)
						{
							hLambda_EPD_v2_pt[ptbin - 4][cent - 1]->Fill(particle.GetMass(), cos(2. * particle.GetPhi() - 2. * EP2_EPD_w_shifted));
							hLambda_EPD_v2_pt[ptbin - 4][cent - 1]->Fill(particle.GetMass(), cos(2. * particle.GetPhi() - 2. * EP2_EPD_e_shifted));
						}
					}
					if (particle.GetPDG() < 0)
					{
						hLambdabarM_pt[ptbin - 4][cent - 1]->Fill(particle.GetMass());
						if (goodTPCEP)
							hLambdabar_TPC_v2_pt[ptbin - 4][cent - 1]->Fill(particle.GetMass(), cos(2. * particle.GetPhi() - 2. * EP2_TPC_RecoPar));
						if (goodEPDEP)
						{
							hLambdabar_EPD_v2_pt[ptbin - 4][cent - 1]->Fill(particle.GetMass(), cos(2. * particle.GetPhi() - 2. * EP2_EPD_w_shifted));
							hLambdabar_EPD_v2_pt[ptbin - 4][cent - 1]->Fill(particle.GetMass(), cos(2. * particle.GetPhi() - 2. * EP2_EPD_e_shifted));
						}
					}
				}
			}

			// general cuts for all RecoPar
			if (particle.GetPt() < 0.14 * RecoPar_nq[2 * recopar + charge_index] || particle.GetPt() > 0.6 * RecoPar_nq[2 * recopar + charge_index])
				continue; // pt cut
			if (fabs(particle.GetRapidity()) > y_coal_cut)
				continue; // rapidity cut
			hRecoParM_cen[2 * recopar + charge_index][cent - 1]->Fill(particle.GetMass());
			if (goodTPCEP)
				hRecoPar_TPC_v2[2 * recopar + charge_index][cent - 1]->Fill(particle.GetMass(), cos(2. * particle.GetPhi() - 2. * EP2_TPC_RecoPar));
			if (goodEPDEP)
			{
				hRecoPar_EPD_v2[2 * recopar + charge_index][cent - 1]->Fill(particle.GetMass(), cos(2. * particle.GetPhi() - 2. * EP2_EPD_w_shifted));
				hRecoPar_EPD_v2[2 * recopar + charge_index][cent - 1]->Fill(particle.GetMass(), cos(2. * particle.GetPhi() - 2. * EP2_EPD_e_shifted));
			}
			
			if (recopar == 2) // Lambda
			{
				// check if good Lambda
				if (!isGoodLambda(cent, particle))
					continue;
				for (int iDaughter = 0; iDaughter < particle.NDaughters(); iDaughter++)
				{
					const int daughterId = particle.DaughterIds()[iDaughter];
					const KFParticle daughter = KFParticleInterface->GetParticles()[daughterId];
					if (fabs(daughter.GetPDG()) == ProtonPdg)
					{
						const int globalTrackId = daughter.DaughterIds()[0];
						int picoTrackIndex = trackMap[globalTrackId];
						StPicoTrack *dauTrack = mPicoDst->track(picoTrackIndex);
						if (!dauTrack)
							continue;
						float dauRap = Eta2y(dauTrack->pMom().Perp(), dauTrack->pMom().Eta(), ProtonPdgMass);
						if (particle.GetRapidity() < 0) hyLambdaDauProton_east->Fill(dauRap);
						if (particle.GetRapidity() > 0) hyLambdaDauProton_west->Fill(dauRap);
					}
					if (fabs(daughter.GetPDG()) == -PionPdg) // PionPdg is -211
					{
						const int globalTrackId = daughter.DaughterIds()[0];
						int picoTrackIndex = trackMap[globalTrackId];
						StPicoTrack *dauTrack = mPicoDst->track(picoTrackIndex);
						if (!dauTrack)
							continue;
						float dauRap = Eta2y(dauTrack->pMom().Perp(), dauTrack->pMom().Eta(), PionPdgMass);
						if (particle.GetRapidity() < 0) hyLambdaDauPion_east->Fill(dauRap);
						if (particle.GetRapidity() > 0) hyLambdaDauPion_west->Fill(dauRap);
					}
				} // iDaughter
			}
		}
	}

	// if (mult_index_bad_flag)
	// 	cout << "mult_index_bad_flag" << endl;
	// if (vz_index_bad_flag)
	// 	cout << "vz_index_bad_flag" << endl;
	// if (EP_index_bad_flag)
	// 	cout << "EP_index_bad_flag" << endl;
	bool index_bad_flag = mult_index_bad_flag || vz_index_bad_flag || EP_index_bad_flag;
	bool event_bad_flag = OmegaVec.size() == 0 || index_bad_flag;

	// Omega loop
	if (cent < min_cent || cent > max_cent)
		return kStOK;

	// remove all daughter kaons
	std::set<int> daughter_kaon_tracks;
	for (int iOmega = 0; iOmega < OmegaVec.size(); iOmega++)
	{
		const KFParticle particle = OmegaVec[iOmega];
		for (int kaon_track = 0; kaon_track < kaon_tracks.size(); kaon_track++)
		{
			StPicoTrack *track = mPicoDst->track(kaon_tracks[kaon_track]);
			if (IsKaonOmegaDaughter(particle, track->id()))
				daughter_kaon_tracks.insert(kaon_tracks[kaon_track]);
		}
	}
	bool filled = false;
	for (int iOmega = 0; iOmega < OmegaVec.size(); iOmega++)
	{
		const KFParticle particle = OmegaVec[iOmega];

		// couting Omega
		if (!filled)
		{
			if (particle.GetPDG() > 0) hOmegaEventUsed->Fill(0.);
			if (particle.GetPDG() < 0) hOmegaEventUsed->Fill(1.);
			filled = true;
		} 
		if (particle.GetPDG() > 0)
			hOmegaUsed->Fill(0.);
		if (particle.GetPDG() < 0)
			hOmegaUsed->Fill(1.);
		if (hasLambdabar && particle.GetPDG() > 0)
			hOmegaUsed_wlb->Fill(0);
		if (hasLambdabar && particle.GetPDG() < 0)
			hOmegaUsed_wlb->Fill(1);
		if (!hasLambdabar && particle.GetPDG() > 0)
			hOmegaUsed_wolb->Fill(0);
		if (!hasLambdabar && particle.GetPDG() < 0)
			hOmegaUsed_wolb->Fill(1);
		if (hasLambda && particle.GetPDG() > 0)
			hOmegaUsed_wl->Fill(0);
		if (hasLambda && particle.GetPDG() < 0)
			hOmegaUsed_wl->Fill(1);
		if (!hasLambda && particle.GetPDG() > 0)
			hOmegaUsed_wol->Fill(0);
		if (!hasLambda && particle.GetPDG() < 0)
			hOmegaUsed_wol->Fill(1);

		// Omega momentum at DCA to PV
		TVector3 pOmega(particle.GetPx(), particle.GetPy(), particle.GetPz());
		TVector3 xOmega(particle.GetX(), particle.GetY(), particle.GetZ());
		StPicoPhysicalHelix helixOmega(pOmega, xOmega, magnet * kilogauss, particle.GetQ());
		double pathlength = helixOmega.pathLength(Vertex3D, false);
		TVector3 pOmega_tb = helixOmega.momentumAt(pathlength, magnet * kilogauss);
		TVector3 dcaOmega = helixOmega.at(pathlength) - Vertex3D;
		if (particle.GetPDG() > 0)
			hOmegaDCAtoPV->Fill(dcaOmega.Mag());
		if (particle.GetPDG() < 0)
			hOmegabarDCAtoPV->Fill(dcaOmega.Mag());
		if (particle.GetPDG() > 0)
			hOmegaPDiff->Fill((pOmega - pOmega_tb).Mag());
		if (particle.GetPDG() < 0)
			hOmegabarPDiff->Fill((pOmega - pOmega_tb).Mag());

		// kaon loop
		for (int kaon_track = 0; kaon_track < kaon_tracks.size(); kaon_track++)
		{
			StPicoTrack *track = mPicoDst->track(kaon_tracks[kaon_track]);
			// TOF rapidity cut
			if (fabs(Eta2y(track->pMom().Perp(), track->pMom().Eta(), KaonPdgMass)) > y_coal_cut)
				continue;
			// Omega daughter removal
			if (daughter_kaon_tracks.find(kaon_tracks[kaon_track]) != daughter_kaon_tracks.end())
				continue;

			float pt = track->pMom().Perp();
			float eta = track->pMom().Eta();
			float TOFEff = 1.0;
			float TOFEff2D = 1.0;
			float TPCEff = track->charge() > 0 ? eff_finder.GetEfficiency1D(pt, cent, "Kp") : eff_finder.GetEfficiency1D(pt, cent, "Km");
			float TPCEff2D = track->charge() > 0 ? eff_finder.GetEfficiency2D(pt, eta, cent, "Kp") : eff_finder.GetEfficiency2D(pt, eta, cent, "Km");
			if (hTOFEff[cent - 1] != 0 && pt > 0.4)
				TOFEff = hTOFEff[cent - 1]->GetEfficiency(hTOFEff[cent - 1]->FindFixBin(pt));
			if (hTOFEff_2D[cent - 1] != 0 && pt > 0.4)
				TOFEff2D = hTOFEff_2D[cent - 1]->GetBinContent(hTOFEff_2D[cent - 1]->FindFixBin(pt, eta));

			// pair-wise should be added after this line
			/* */

			// k*
			TLorentzVector lv1;
			lv1.SetVectM(pOmega_tb, OmegaPdgMass);
			TLorentzVector lv2;
			lv2.SetVectM(track->pMom(), KaonPdgMass);
			TLorentzVector lv2_ori;
			lv2_ori.SetVectM(track->pMom(), KaonPdgMass);
			double dpt = fabs(lv1.Perp() - lv2.Perp());
			double dy = fabs(lv1.Rapidity() - lv2.Rapidity());
			double dphi = fabs(lv1.Vect().DeltaPhi(lv2.Vect()));
			TLorentzVector P = lv1 + lv2;
			TVector3 pair_beta = P.BoostVector();
			lv1.Boost((-1) * pair_beta);
			lv2.Boost((-1) * pair_beta);
			double kstar = 0.5 * (lv1 - lv2).Vect().Mag();

			float eff = TPCEff * TOFEff;
			float eff2D = TPCEff2D * TOFEff2D;

			// "KplusO", "KplusObar", "KminusO", "KminusObar"
			int corr_index = (track->charge() < 0) * 2 + (particle.GetQ() > 0);
			hCorrKO[corr_index]->Fill(kstar, 1. / eff);
			hPtCorrKO[corr_index]->Fill(dpt, 1. / eff);
			hyCorrKO[corr_index]->Fill(dy, 1. / eff);
			hphiCorrKO[corr_index]->Fill(dphi, 1. / eff);
			hCorrKO_same[corr_index]->Fill(kstar);
			hCorrKO_2D_same[corr_index]->Fill(dy, dphi);
		}

		// lambda loop
		for (int iLambda = 0; iLambda < LambdaVec.size(); iLambda++)
		{
			KFParticle Lambda = LambdaVec[iLambda];

			// do we need to worry about efficiency?

			TVector3 pLambda(Lambda.GetPx(), Lambda.GetPy(), Lambda.GetPz());
			TVector3 xLambda(Lambda.GetX(), Lambda.GetY(), Lambda.GetZ());
			StPicoPhysicalHelix helixLambda(pLambda, xLambda, 0, 0);
			double pathlength = helixLambda.pathLength(Vertex3D, false);
			TVector3 pLambda_tb = helixLambda.momentumAt(pathlength, 0);
			TVector3 posLambda_to_PV = xLambda - Vertex3D;
			double rdotp = posLambda_to_PV.Dot(pLambda);
			double dcaLambda_to_PV = sqrt(posLambda_to_PV.Mag2() - rdotp * rdotp / pLambda.Mag2());
			if (iOmega == 0)
			{
				if (Lambda.GetPDG() > 0)
					hLambdaDCAtoPV->Fill(dcaLambda_to_PV);
				if (Lambda.GetPDG() < 0)
					hLambdabarDCAtoPV->Fill(dcaLambda_to_PV);
				if (Lambda.GetPDG() > 0)
					hLambdaPDiff->Fill((pLambda - pLambda_tb).Mag());
				if (Lambda.GetPDG() < 0)
					hLambdabarPDiff->Fill((pLambda - pLambda_tb).Mag());
			}
			if (dcaLambda_to_PV > 0.4)
				continue;

			// k* for Lambda
			TLorentzVector lv1;
			lv1.SetVectM(pOmega_tb, OmegaPdgMass);
			TLorentzVector lv2;
			lv2.SetVectM(pLambda_tb, LambdaPdgMass);
			double dpt = fabs(lv1.Perp() - lv2.Perp());
			double dy = fabs(lv1.Rapidity() - lv2.Rapidity());
			double dphi = fabs(lv1.Vect().DeltaPhi(lv2.Vect()));
			TLorentzVector P = lv1 + lv2;
			TVector3 pair_beta = P.BoostVector();
			lv1.Boost((-1) * pair_beta);
			lv2.Boost((-1) * pair_beta);
			double kstar = 0.5 * (lv1 - lv2).Vect().Mag();

			// "LambdabarO"
			// "LambdabarO", "LambdabarObar", "LambdaO", "LambdaObar"
			int corr_index = (Lambda.GetPDG() > 0) * 2 + (particle.GetQ() > 0);
			hCorrLO[corr_index]->Fill(kstar);
			hyCorrLO[corr_index]->Fill(dy);
			hphiCorrLO[corr_index]->Fill(dphi);
			hCorrLO_same[corr_index]->Fill(kstar); // here this is kind of redundant
			hCorrLO_2D_same[corr_index]->Fill(dy, dphi);
		}
	} // End loop over regular Omega

	// Omega sideband loop
	filled = false;
	for (int iOmega = 0; iOmega < OmegaSidebandVec.size(); iOmega++)
	{
		const KFParticle particle = OmegaSidebandVec[iOmega];

		// couting Omega
		if (!filled)
		{
			if (particle.GetPDG() > 0) hOmegaEventUsed_sideband->Fill(0.);
			if (particle.GetPDG() < 0) hOmegaEventUsed_sideband->Fill(1.);
			filled = true;
		} 
		if (particle.GetPDG() > 0)
			hOmegaUsed_sideband->Fill(0.);
		if (particle.GetPDG() < 0)
			hOmegaUsed_sideband->Fill(1.);
		if (hasLambdabar && particle.GetPDG() > 0)
			hOmegaUsed_wlb_sideband->Fill(0.);
		if (hasLambdabar && particle.GetPDG() < 0)
			hOmegaUsed_wlb_sideband->Fill(1.);
		if (!hasLambdabar && particle.GetPDG() > 0)
			hOmegaUsed_wolb_sideband->Fill(0.);
		if (!hasLambdabar && particle.GetPDG() < 0)
			hOmegaUsed_wolb_sideband->Fill(1.);
		if (hasLambda && particle.GetPDG() > 0)
			hOmegaUsed_wl_sideband->Fill(0.);
		if (hasLambda && particle.GetPDG() < 0)
			hOmegaUsed_wl_sideband->Fill(1.);
		if (!hasLambda && particle.GetPDG() > 0)
			hOmegaUsed_wol_sideband->Fill(0.);
		if (!hasLambda && particle.GetPDG() < 0)
			hOmegaUsed_wol_sideband->Fill(1.);

		// Omega momentum at DCA to PV
		TVector3 pOmega(particle.GetPx(), particle.GetPy(), particle.GetPz());
		TVector3 xOmega(particle.GetX(), particle.GetY(), particle.GetZ());
		StPicoPhysicalHelix helixOmega(pOmega, xOmega, magnet * kilogauss, particle.GetQ());
		double pathlength = helixOmega.pathLength(Vertex3D, false);
		TVector3 pOmega_tb = helixOmega.momentumAt(pathlength, magnet * kilogauss);
		TVector3 dcaOmega = helixOmega.at(pathlength) - Vertex3D;

		// kaon loop
		for (int kaon_track = 0; kaon_track < kaon_tracks.size(); kaon_track++)
		{
			StPicoTrack *track = mPicoDst->track(kaon_tracks[kaon_track]);
			// TOF rapidity cut
			if (fabs(Eta2y(track->pMom().Perp(), track->pMom().Eta(), KaonPdgMass)) > y_coal_cut)
				continue;
			// Omega daughter removal
			if (daughter_kaon_tracks.find(kaon_tracks[kaon_track]) != daughter_kaon_tracks.end())
				continue;

			float pt = track->pMom().Perp();
			float eta = track->pMom().Eta();
			float TOFEff = 1.0;
			float TOFEff2D = 1.0;
			float TPCEff = track->charge() > 0 ? eff_finder.GetEfficiency1D(pt, cent, "Kp") : eff_finder.GetEfficiency1D(pt, cent, "Km");
			float TPCEff2D = track->charge() > 0 ? eff_finder.GetEfficiency2D(pt, eta, cent, "Kp") : eff_finder.GetEfficiency2D(pt, eta, cent, "Km");
			if (hTOFEff[cent - 1] != 0 && pt > 0.4)
				TOFEff = hTOFEff[cent - 1]->GetEfficiency(hTOFEff[cent - 1]->FindFixBin(pt));
			if (hTOFEff_2D[cent - 1] != 0 && pt > 0.4)
				TOFEff2D = hTOFEff_2D[cent - 1]->GetBinContent(hTOFEff_2D[cent - 1]->FindFixBin(pt, eta));

			// pair-wise should be added after this line
			/* */

			// k*
			TLorentzVector lv1;
			lv1.SetVectM(pOmega_tb, OmegaPdgMass);
			TLorentzVector lv2;
			lv2.SetVectM(track->pMom(), KaonPdgMass);
			TLorentzVector lv2_ori;
			lv2_ori.SetVectM(track->pMom(), KaonPdgMass);
			double dpt = fabs(lv1.Perp() - lv2.Perp());
			double dy = fabs(lv1.Rapidity() - lv2.Rapidity());
			double dphi = fabs(lv1.Vect().DeltaPhi(lv2.Vect()));
			TLorentzVector P = lv1 + lv2;
			TVector3 pair_beta = P.BoostVector();
			lv1.Boost((-1) * pair_beta);
			lv2.Boost((-1) * pair_beta);
			double kstar = 0.5 * (lv1 - lv2).Vect().Mag();

			float eff = TPCEff * TOFEff;
			float eff2D = TPCEff2D * TOFEff2D;

			// "KplusO", "KplusObar", "KminusO", "KminusObar"
			int corr_index = (track->charge() < 0) * 2 + (particle.GetQ() > 0);
			hCorrKO_sideband[corr_index]->Fill(kstar, 1. / eff);
			hPtCorrKO_sideband[corr_index]->Fill(dpt, 1. / eff);
			hyCorrKO_sideband[corr_index]->Fill(dy, 1. / eff);
			hphiCorrKO_sideband[corr_index]->Fill(dphi, 1. / eff);
		}

		// lambda loop
		for (int iLambda = 0; iLambda < LambdaVec.size(); iLambda++)
		{
			KFParticle Lambda = LambdaVec[iLambda];

			// do we need to worry about efficiency?

			TVector3 pLambda(Lambda.GetPx(), Lambda.GetPy(), Lambda.GetPz());
			TVector3 xLambda(Lambda.GetX(), Lambda.GetY(), Lambda.GetZ());
			StPicoPhysicalHelix helixLambda(pLambda, xLambda, 0, 0);
			double pathlength = helixLambda.pathLength(Vertex3D, false);
			TVector3 pLambda_tb = helixLambda.momentumAt(pathlength, 0);
			TVector3 posLambda_to_PV = xLambda - Vertex3D;
			double rdotp = posLambda_to_PV.Dot(pLambda);
			double dcaLambda_to_PV = sqrt(posLambda_to_PV.Mag2() - rdotp * rdotp / pLambda.Mag2());
			if (dcaLambda_to_PV > 0.4)
				continue;

			// k* for Lambda
			TLorentzVector lv1;
			lv1.SetVectM(pOmega_tb, OmegaPdgMass);
			TLorentzVector lv2;
			lv2.SetVectM(pLambda_tb, LambdaPdgMass);
			double dpt = fabs(lv1.Perp() - lv2.Perp());
			double dy = fabs(lv1.Rapidity() - lv2.Rapidity());
			double dphi = fabs(lv1.Vect().DeltaPhi(lv2.Vect()));
			TLorentzVector P = lv1 + lv2;
			TVector3 pair_beta = P.BoostVector();
			lv1.Boost((-1) * pair_beta);
			lv2.Boost((-1) * pair_beta);
			double kstar = 0.5 * (lv1 - lv2).Vect().Mag();

			// "LambdabarO"
			// "LambdabarO", "LambdabarObar", "LambdaO", "LambdaObar"
			int corr_index = (Lambda.GetPDG() > 0) * 2 + (particle.GetQ() > 0);
			hCorrLO_sideband[corr_index]->Fill(kstar);
			hyCorrLO_sideband[corr_index]->Fill(dy);
			hphiCorrLO_sideband[corr_index]->Fill(dphi);
		}
	} // End loop over sideband Omega

	// Lambda loop
	filled = false;
	for (int iOmega = 0; iOmega < LambdaVec.size(); iOmega++)
	{
		const KFParticle particle = LambdaVec[iOmega];

		// couting Lambda
		if (!filled)
		{
			if (particle.GetPDG() > 0) hLambdaEventUsed->Fill(0.);
			if (particle.GetPDG() < 0) hLambdaEventUsed->Fill(1.);
			filled = true;
		}

		if (particle.GetPDG() > 0)
			hLambdaUsed->Fill(0.);
		if (particle.GetPDG() < 0)
			hLambdaUsed->Fill(1.);

		// Omega momentum at DCA to PV
		TVector3 pOmega(particle.GetPx(), particle.GetPy(), particle.GetPz());
		TVector3 xOmega(particle.GetX(), particle.GetY(), particle.GetZ());
		StPicoPhysicalHelix helixOmega(pOmega, xOmega, 0, 0);
		double pathlength = helixOmega.pathLength(Vertex3D, false);
		TVector3 pOmega_tb = helixOmega.momentumAt(pathlength, 0);
		TVector3 dcaOmega = helixOmega.at(pathlength) - Vertex3D;

		// kaon loop
		for (int kaon_track = 0; kaon_track < kaon_tracks.size(); kaon_track++)
		{
			StPicoTrack *track = mPicoDst->track(kaon_tracks[kaon_track]);
			// TOF rapidity cut
			if (fabs(Eta2y(track->pMom().Perp(), track->pMom().Eta(), KaonPdgMass)) > y_coal_cut)
				continue;
			// Omega daughter removal
			if (daughter_kaon_tracks.find(kaon_tracks[kaon_track]) != daughter_kaon_tracks.end())
				continue;

			float pt = track->pMom().Perp();
			float eta = track->pMom().Eta();
			float TOFEff = 1.0;
			float TOFEff2D = 1.0;
			float TPCEff = track->charge() > 0 ? eff_finder.GetEfficiency1D(pt, cent, "Kp") : eff_finder.GetEfficiency1D(pt, cent, "Km");
			float TPCEff2D = track->charge() > 0 ? eff_finder.GetEfficiency2D(pt, eta, cent, "Kp") : eff_finder.GetEfficiency2D(pt, eta, cent, "Km");
			if (hTOFEff[cent - 1] != 0 && pt > 0.4)
				TOFEff = hTOFEff[cent - 1]->GetEfficiency(hTOFEff[cent - 1]->FindFixBin(pt));
			if (hTOFEff_2D[cent - 1] != 0 && pt > 0.4)
				TOFEff2D = hTOFEff_2D[cent - 1]->GetBinContent(hTOFEff_2D[cent - 1]->FindFixBin(pt, eta));

			// pair-wise should be added after this line
			/* */

			// k*
			TLorentzVector lv1;
			lv1.SetVectM(pOmega_tb, LambdaPdgMass);
			TLorentzVector lv2;
			lv2.SetVectM(track->pMom(), KaonPdgMass);
			TLorentzVector lv2_ori;
			lv2_ori.SetVectM(track->pMom(), KaonPdgMass);
			double dpt = fabs(lv1.Perp() - lv2.Perp());
			double dy = fabs(lv1.Rapidity() - lv2.Rapidity());
			double dphi = fabs(lv1.Vect().DeltaPhi(lv2.Vect()));
			TLorentzVector P = lv1 + lv2;
			TVector3 pair_beta = P.BoostVector();
			lv1.Boost((-1) * pair_beta);
			lv2.Boost((-1) * pair_beta);
			double kstar = 0.5 * (lv1 - lv2).Vect().Mag();

			float eff = TPCEff * TOFEff;
			float eff2D = TPCEff2D * TOFEff2D;

			// "KplusO", "KplusObar", "KminusO", "KminusObar"
			int corr_index = (track->charge() < 0) * 2 + (particle.GetPDG() < 0);
			hCorrKL[corr_index]->Fill(kstar, 1. / eff);
			hPtCorrKL[corr_index]->Fill(dpt, 1. / eff);
			hyCorrKL[corr_index]->Fill(dy, 1. / eff);
			hphiCorrKL[corr_index]->Fill(dphi, 1. / eff);
			hCorrKL_same[corr_index]->Fill(kstar);
		}
	} // End loop over regular Lambda

	// Lambda sideband loop
	for (int iOmega = 0; iOmega < LambdaSidebandVec.size(); iOmega++)
	{
		const KFParticle particle = LambdaSidebandVec[iOmega];

		// couting Lambda
		if (!filled)
		{
			if (particle.GetPDG() > 0) hLambdaEventUsed_sideband->Fill(0.);
			if (particle.GetPDG() < 0) hLambdaEventUsed_sideband->Fill(1.);
			filled = true;
		}

		if (particle.GetPDG() > 0)
			hLambdaUsed_sideband->Fill(0.);
		if (particle.GetPDG() < 0)
			hLambdaUsed_sideband->Fill(1.);

		// Omega momentum at DCA to PV
		TVector3 pOmega(particle.GetPx(), particle.GetPy(), particle.GetPz());
		TVector3 xOmega(particle.GetX(), particle.GetY(), particle.GetZ());
		StPicoPhysicalHelix helixOmega(pOmega, xOmega, 0, 0);
		double pathlength = helixOmega.pathLength(Vertex3D, false);
		TVector3 pOmega_tb = helixOmega.momentumAt(pathlength, 0);
		TVector3 dcaOmega = helixOmega.at(pathlength) - Vertex3D;

		// kaon loop
		for (int kaon_track = 0; kaon_track < kaon_tracks.size(); kaon_track++)
		{
			StPicoTrack *track = mPicoDst->track(kaon_tracks[kaon_track]);
			// TOF rapidity cut
			if (fabs(Eta2y(track->pMom().Perp(), track->pMom().Eta(), KaonPdgMass)) > y_coal_cut)
				continue;
			// Omega daughter removal
			if (daughter_kaon_tracks.find(kaon_tracks[kaon_track]) != daughter_kaon_tracks.end())
				continue;

			float pt = track->pMom().Perp();
			float eta = track->pMom().Eta();
			float TOFEff = 1.0;
			float TOFEff2D = 1.0;
			float TPCEff = track->charge() > 0 ? eff_finder.GetEfficiency1D(pt, cent, "Kp") : eff_finder.GetEfficiency1D(pt, cent, "Km");
			float TPCEff2D = track->charge() > 0 ? eff_finder.GetEfficiency2D(pt, eta, cent, "Kp") : eff_finder.GetEfficiency2D(pt, eta, cent, "Km");
			if (hTOFEff[cent - 1] != 0 && pt > 0.4)
				TOFEff = hTOFEff[cent - 1]->GetEfficiency(hTOFEff[cent - 1]->FindFixBin(pt));
			if (hTOFEff_2D[cent - 1] != 0 && pt > 0.4)
				TOFEff2D = hTOFEff_2D[cent - 1]->GetBinContent(hTOFEff_2D[cent - 1]->FindFixBin(pt, eta));

			// pair-wise should be added after this line
			/* */

			// k*
			TLorentzVector lv1;
			lv1.SetVectM(pOmega_tb, LambdaPdgMass);
			TLorentzVector lv2;
			lv2.SetVectM(track->pMom(), KaonPdgMass);
			TLorentzVector lv2_ori;
			lv2_ori.SetVectM(track->pMom(), KaonPdgMass);
			double dpt = fabs(lv1.Perp() - lv2.Perp());
			double dy = fabs(lv1.Rapidity() - lv2.Rapidity());
			double dphi = fabs(lv1.Vect().DeltaPhi(lv2.Vect()));
			TLorentzVector P = lv1 + lv2;
			TVector3 pair_beta = P.BoostVector();
			lv1.Boost((-1) * pair_beta);
			lv2.Boost((-1) * pair_beta);
			double kstar = 0.5 * (lv1 - lv2).Vect().Mag();

			float eff = TPCEff * TOFEff;
			float eff2D = TPCEff2D * TOFEff2D;

			// "KplusO", "KplusObar", "KminusO", "KminusObar"
			int corr_index = (track->charge() < 0) * 2 + (particle.GetPDG() < 0);
			hCorrKL_sideband[corr_index]->Fill(kstar, 1. / eff);
			hPtCorrKL_sideband[corr_index]->Fill(dpt, 1. / eff);
			hyCorrKL_sideband[corr_index]->Fill(dy, 1. / eff);
			hphiCorrKL_sideband[corr_index]->Fill(dphi, 1. / eff);
		}
	} // End loop over sideband Lambda

	// new observable
	if (OmegaVec.size() == 1 && OmegaVec[0].GetQ() < 0)
		hKratio_omega->Fill(pbct * 1.0 / pct, kmct * 1.0 / kpct);
	if (OmegaVec.size() == 0)
		hKratio_wo->Fill(pbct * 1.0 / pct, kmct * 1.0 / kpct);
	if (OmegaVec.size() == 1 && OmegaVec[0].GetQ() > 0)
		hKratio_omegabar->Fill(pbct * 1.0 / pct, kmct * 1.0 / kpct);

	// filling trees
	if (StoringTree && !event_bad_flag)
		omega_mix[mult_index][vz_index][EP_index]->Fill();
	// if (StoringTree && !OmegaVec.size() == 0) omega_mix[0][0][0]->Fill(); // for testing

	// counting kaon, need to be right before event mixing
	if (!current_event.IsEmptyEvent())
		hKaonCt->Fill(1.0, kaon_tracks.size());
	else
	{
		hKaonCt->Fill(0.0, kaon_tracks.size());
		return kStOK;
	}

	// mixed event
	if (!PerformMixing)
		return kStOK;
	if (OmegaVec.size() == 0)
	{
		std::cout << "OmegaVec.size() == 0" << std::endl;
		std::cout << "Will not mix with events with no Omega" << std::endl;
		return kStOK;
	}
	std::vector<float> *px_omega = 0;
	std::vector<float> *py_omega = 0;
	std::vector<float> *pz_omega = 0;
	std::vector<int> *charge_omega = 0;
	int evtid_omega = 0;
	int runid_omega = 0;
	omega_mix[mult_index][vz_index][EP_index]->SetBranchAddress("px", &px_omega);
	omega_mix[mult_index][vz_index][EP_index]->SetBranchAddress("py", &py_omega);
	omega_mix[mult_index][vz_index][EP_index]->SetBranchAddress("pz", &pz_omega);
	omega_mix[mult_index][vz_index][EP_index]->SetBranchAddress("charge", &charge_omega);
	omega_mix[mult_index][vz_index][EP_index]->SetBranchAddress("evt_id", &evtid_omega);
	omega_mix[mult_index][vz_index][EP_index]->SetBranchAddress("run_id", &runid_omega);
	for (int iMixEvent = 0; iMixEvent < omega_mix[mult_index][vz_index][EP_index]->GetEntries(); iMixEvent++)
	{
		omega_mix[mult_index][vz_index][EP_index]->GetEntry(iMixEvent);
		if (runID == runid_omega && evtID == evtid_omega)
			continue; // no self-correlation
		for (int iOmega = 0; iOmega < px_omega->size(); iOmega++)
		{
			if (charge_omega->at(iOmega) < 0)
				hOmegaUsed->Fill(2.);
			if (charge_omega->at(iOmega) > 0)
				hOmegaUsed->Fill(3.);

			TVector3 p_omega(px_omega->at(iOmega), py_omega->at(iOmega), pz_omega->at(iOmega));

			for (int kaon_track = 0; kaon_track < kaon_tracks.size(); kaon_track++)
			{
				StPicoTrack *track = mPicoDst->track(kaon_tracks[kaon_track]);
				if (fabs(Eta2y(track->pMom().Perp(), track->pMom().Eta(), KaonPdgMass)) > y_coal_cut)
					continue;

				// k*
				TLorentzVector lv1;
				lv1.SetVectM(p_omega, OmegaPdgMass);
				TLorentzVector lv2;
				lv2.SetVectM(track->pMom(), KaonPdgMass);
				double dpt = fabs(lv1.Perp() - lv2.Perp());
				double dy = fabs(lv1.Rapidity() - lv2.Rapidity());
				double dphi = fabs(lv1.Vect().DeltaPhi(lv2.Vect()));
				TLorentzVector P = lv1 + lv2;
				TVector3 pair_beta = P.BoostVector();
				lv1.Boost((-1) * pair_beta);
				lv2.Boost((-1) * pair_beta);
				double kstar = 0.5 * (lv1 - lv2).Vect().Mag();

				int corr_index = (track->charge() < 0) * 2 + (charge_omega->at(iOmega) > 0);
				hCorrKO_mixed[corr_index]->Fill(kstar);
				hCorrKO_2D_mixed[corr_index]->Fill(dy, dphi);
			}

			for (int iLambda = 0; iLambda < LambdaVec.size(); iLambda++)
			{
				KFParticle Lambda = LambdaVec[iLambda];

				// do we need to worry about efficiency?

				TVector3 pLambda(Lambda.GetPx(), Lambda.GetPy(), Lambda.GetPz());
				TVector3 xLambda(Lambda.GetX(), Lambda.GetY(), Lambda.GetZ());
				StPicoPhysicalHelix helixLambda(pLambda, xLambda, 0, 0);
				double pathlength = helixLambda.pathLength(Vertex3D, false);
				TVector3 pLambda_tb = helixLambda.momentumAt(pathlength, 0);
				TVector3 posLambda_to_PV = xLambda - Vertex3D;
				double rdotp = posLambda_to_PV.Dot(pLambda);
				double dcaLambda_to_PV = sqrt(posLambda_to_PV.Mag2() - rdotp * rdotp / pLambda.Mag2());
				if (dcaLambda_to_PV > 0.4)
					continue;

				// k* for Lambda
				TLorentzVector lv1;
				lv1.SetVectM(p_omega, OmegaPdgMass);
				TLorentzVector lv2;
				lv2.SetVectM(pLambda_tb, LambdaPdgMass);
				double dpt = fabs(lv1.Perp() - lv2.Perp());
				double dy = fabs(lv1.Rapidity() - lv2.Rapidity());
				double dphi = fabs(lv1.Vect().DeltaPhi(lv2.Vect()));
				TLorentzVector P = lv1 + lv2;
				TVector3 pair_beta = P.BoostVector();
				lv1.Boost((-1) * pair_beta);
				lv2.Boost((-1) * pair_beta);
				double kstar = 0.5 * (lv1 - lv2).Vect().Mag();

				// "LambdabarO"
				// "LambdabarO", "LambdabarObar", "LambdaO", "LambdaObar"
				int corr_index = (Lambda.GetPDG() > 0) * 2 + (charge_omega->at(iOmega) > 0);
				hCorrLO_mixed[corr_index]->Fill(kstar);
				hCorrLO_2D_mixed[corr_index]->Fill(dy, dphi);
			}
		}
	}
	// ======= KFParticle end ======= //

	/////////////////////////////////////////////////////////
	return kStOK;
}

bool StKFParticleAnalysisMaker::isGoodOmega(int cent, KFParticle Omega)
{
	int cent_bin = 0;
	if (cent >= 1 && cent <= 3)
		cent_bin = 0;
	else
		cent_bin = cent - 3;
	if (Omega.GetQ() > 0)
		return (fabs(Omega.GetMass() - OmegaPdgMass) < 3 * OmegabarMassSigma);
	else
		return (fabs(Omega.GetMass() - OmegaPdgMass) < 3 * OmegaMassSigma);
}

bool StKFParticleAnalysisMaker::isGoodLambda(int cent, KFParticle Lambda)
{
	if (Lambda.GetPDG() == 3122)
		return (fabs(Lambda.GetMass() - LambdaPdgMass) < 3 * LambdaMassSigma);
	if (Lambda.GetPDG() == -3122)
		return (fabs(Lambda.GetMass() - LambdaPdgMass) < 3 * LambdabarMassSigma);
	return false;
}

bool StKFParticleAnalysisMaker::isSidebandOmega(int cent, KFParticle Omega)
{
	int cent_bin = 0;
	if (cent >= 1 && cent <= 3)
		cent_bin = 0;
	else
		cent_bin = cent - 3;
	if (Omega.GetQ() > 0)
		return (fabs(Omega.GetMass() - OmegaPdgMass) > 4 * OmegabarMassSigma && fabs(Omega.GetMass() - OmegaPdgMass) < 6 * OmegabarMassSigma);
	else
		return (fabs(Omega.GetMass() - OmegaPdgMass) > 4 * OmegaMassSigma && fabs(Omega.GetMass() - OmegaPdgMass) < 6 * OmegaMassSigma);
}

bool StKFParticleAnalysisMaker::isSidebandLambda(int cent, KFParticle Lambda)
{
	if (Lambda.GetPDG() == 3122)
		return (fabs(Lambda.GetMass() - LambdaPdgMass) > 4 * LambdaMassSigma && fabs(Lambda.GetMass() - LambdaPdgMass) < 6 * LambdaMassSigma);
	if (Lambda.GetPDG() == -3122)
		return (fabs(Lambda.GetMass() - LambdaPdgMass) > 4 * LambdabarMassSigma && fabs(Lambda.GetMass() - LambdaPdgMass) < 6 * LambdabarMassSigma);
	return false;
}

int StKFParticleAnalysisMaker::MixRefMultBin(int cent, int refmult)
{
	if (cent != 8 && cent != 9)
		return -1;
	int bin = -1;
	int refmult_cut[] = {6, 9, 13, 18, 25, 33, 44, 57, 72, 91, 113, 138, 168, 204, 246, 295}; // left inclusive, 27 GeV

	bin = static_cast<int>(floor((refmult - refmult_cut[14]) / 20));
	if (bin > 6)
		bin = 6;
	return bin;
}

//------------------------------------------------------------
bool StKFParticleAnalysisMaker::isGoodObs(double Obs)
{

	bool mGoodObs = true;
	mGoodObs = mGoodObs && !std::isnan(Obs);
	mGoodObs = mGoodObs && !std::isinf(Obs);
	// mGoodObs = mGoodObs && fabs(Obs)<2.1;

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
void StKFParticleAnalysisMaker::SetupKFParticle()
{
	int maxGBTrackIndex = -1; // find max global track index
	for (unsigned int iTrack = 0; iTrack < PicoDst->numberOfTracks(); iTrack++)
	{
		StPicoTrack *track = PicoDst->track(iTrack);
		if (!track)
			continue;
		if (track->id() > maxGBTrackIndex)
			maxGBTrackIndex = track->id();
	}
	vector<KFMCTrack> mcTracks(0);
	vector<int> triggeredTracks;
	vector<int> mcIndices(maxGBTrackIndex + 1);
	for (unsigned int iIndex = 0; iIndex < mcIndices.size(); iIndex++)
		mcIndices[iIndex] = -1;
	if (maxGBTrackIndex > 0)
		KFParticleInterface->ResizeTrackPidVectors(maxGBTrackIndex + 1);

	if (!KFParticleInterface->ProcessEvent(PicoDst, triggeredTracks, sys_tag))
		InterfaceCantProcessEvent = true;
	else
		InterfaceCantProcessEvent = false;

	trackMap.resize(maxGBTrackIndex + 1, -1); // make a map from trackID to track index in global track array
	for (unsigned int iTrack = 0; iTrack < PicoDst->numberOfTracks(); iTrack++)
	{
		StPicoTrack *track = PicoDst->track(iTrack);
		if (!track)
			continue;
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
void StKFParticleAnalysisMaker::SetDaughterTrackPointers(int iKFParticle)
{ // Get Daughter Tracks to calculate decay variables manually
	const KFParticle particle = KFParticleInterface->GetParticles()[iKFParticle];
	int upQ;
	if (particle.GetPDG() == LambdaPdg)
		upQ = 1;
	else if (particle.GetPDG() == -1 * LambdaPdg)
		upQ = -1;
	else
		upQ = 0;
	ProtonTrackIndex = -99999, PionTrackIndex = -99999;
	for (int iDaughter = 0; iDaughter < particle.NDaughters(); iDaughter++)
	{
		const int daughterId = particle.DaughterIds()[iDaughter];
		const KFParticle daughter = KFParticleInterface->GetParticles()[daughterId];
		const int globalTrackId = daughter.DaughterIds()[0];
		int trackIndex = trackMap[globalTrackId];
		if (daughter.GetPDG() == upQ * ProtonPdg)
			ProtonTrackIndex = trackIndex;
		else if (daughter.GetPDG() == upQ * PionPdg)
			PionTrackIndex = trackIndex;
	} // iDaughter
	ProtonTrack = PicoDst->track(ProtonTrackIndex);
	PionTrack = PicoDst->track(PionTrackIndex);
} // void SetDaughterTrackPointers

bool StKFParticleAnalysisMaker::IsKaonOmegaDaughter(KFParticle particle, int kaonTrackId)
{
	if (fabs(particle.GetPDG()) != OmegaPdg)
		return false;
	for (int iDaughter = 0; iDaughter < particle.NDaughters(); iDaughter++)
	{
		const int daughterId = particle.DaughterIds()[iDaughter];
		const KFParticle daughter = KFParticleInterface->GetParticles()[daughterId];
		if (fabs(daughter.GetPDG()) != KaonPdg)
			continue;
		const int globalTrackId = daughter.DaughterIds()[0];

		if (globalTrackId == kaonTrackId)
			return true;
	} // iDaughter
	return false;
}

float StKFParticleAnalysisMaker::ShiftPOIPhi(float phi, int day, int cent)
{
	float phi_shifted = phi;
	if (hTPCPOIShiftInput_cos != 0) // attempt to read input
	{
		for (int i = 1; i <= shift_order_POI; i++)
		{
			float cosAve = hTPCPOIShiftInput_cos->GetBinContent(day, i, cent);
			float sinAve = hTPCPOIShiftInput_sin->GetBinContent(day, i, cent);
			phi_shifted += 2.0 / i * (cosAve * sin(i * 1.0 * phi) - sinAve * cos(i * 1.0 * phi));
		}
	}
	return phi_shifted;
}

float StKFParticleAnalysisMaker::ShiftAssoPhi(float phi, int day, int cent)
{
	float phi_shifted = phi;
	if (hTPCAssoShiftInput_cos != 0) // attempt to read input
	{
		for (int i = 1; i <= shift_order_asso; i++)
		{
			float cosAve = hTPCAssoShiftInput_cos->GetBinContent(day, i, cent);
			float sinAve = hTPCAssoShiftInput_sin->GetBinContent(day, i, cent);
			phi_shifted += 2.0 / i * (cosAve * sin(i * 1.0 * phi) - sinAve * cos(i * 1.0 * phi));
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
		for (int i = 1; i <= shift_order_EP; i++)
		{
			float cosAve = hTPCEPShiftInput_cos[ewFull]->GetBinContent(day, i, cent);
			float sinAve = hTPCEPShiftInput_sin[ewFull]->GetBinContent(day, i, cent);
			psi_shifted += 2.0 / order / i * (cosAve * sin(i * order * psi) - sinAve * cos(i * order * psi));
		}
	}
	// folding
	float wraparound = 2.0 * TMath::Pi() / order;
	if (psi_shifted < 0)
		psi_shifted += wraparound;
	if (psi_shifted > wraparound)
		psi_shifted -= wraparound;
	return psi_shifted;
}

float StKFParticleAnalysisMaker::ShiftEPDEP(float psi, int day, int cent, int ewFull, int order)
{
	float psi_shifted = psi;
	if (order == 1 && hEPDEPShiftInput_1_cos[ewFull] != 0) // attempt to read input
	{
		for (int i = 1; i <= shift_order_EP; i++)
		{
			float cosAve = hEPDEPShiftInput_1_cos[ewFull]->GetBinContent(day, i, cent);
			float sinAve = hEPDEPShiftInput_1_sin[ewFull]->GetBinContent(day, i, cent);
			psi_shifted += 2.0 / order / i * (cosAve * sin(i * order * psi) - sinAve * cos(i * order * psi));
		}
	}
	if (order == 2 && hEPDEPShiftInput_2_cos[ewFull] != 0) // attempt to read input
	{
		for (int i = 1; i <= shift_order_EP; i++)
		{
			float cosAve = hEPDEPShiftInput_2_cos[ewFull]->GetBinContent(day, i, cent);
			float sinAve = hEPDEPShiftInput_2_sin[ewFull]->GetBinContent(day, i, cent);
			psi_shifted += 2.0 / order / i * (cosAve * sin(i * order * psi) - sinAve * cos(i * order * psi));
		}
	}
	// folding
	float wraparound = 2.0 * TMath::Pi() / order;
	if (psi_shifted < 0)
		psi_shifted += wraparound;
	if (psi_shifted > wraparound)
		psi_shifted -= wraparound;
	return psi_shifted;
}

void StKFParticleAnalysisMaker::CutDecider(KFParticle Omega, TH1D *hist_signal, TH1D *hist_sideband, double value)
{
	double mass_diff = fabs(Omega.GetMass() - OmegaPdgMass);
	if (mass_diff < 2 * 0.0021)
		hist_signal->Fill(value);
	if (mass_diff > 4 * 0.0021 && mass_diff < 6 * 0.0021)
		hist_sideband->Fill(value);
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
	float p = pt * cosh(eta);
	float pl = pt * sinh(eta);
	float energy = sqrt(p * p + mass * mass);
	float rapidity = 0.5 * log((energy + pl) / (energy - pl));
	return rapidity;
}