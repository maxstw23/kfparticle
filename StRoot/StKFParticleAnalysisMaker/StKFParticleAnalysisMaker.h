#ifndef StKFParticleAnalysisMaker_h
#define StKFParticleAnalysisMaker_h

#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"
#include "StRoot/StEpdUtil/StEpdEpFinder.h"
#include "StRoot/StEpdUtil/StEpdEpInfo.h"

#include "StKFParticleInterface.h"
#include "StKFParticlePerformanceInterface.h"

#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoEpdHit.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"
#include "StPicoEvent/StPicoPhysicalHelix.h"
#include "TChain.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "StMaker.h"
#include "TString.h"
#include "TObject.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
#include "TEfficiency.h"
#include <cmath>
#include <vector>
#include <algorithm>
#include <iostream>
#include "TF1.h"
#include "StTrackHelix.h"
#include "my_event.h"
#include "MixedBuffer.h"
#include "KaonPID.h"
#include "ProtonPID.h"
#include "PionPID.h"
#include "Efficiency.h"
#include "MyConstant.h" // must include

class StPicoDst;
class StChain;
class StPicoDstMaker;
class TString;
class KFParticle;
class StKFParticleInterface;
class StKFParticlePerformanceInterface;
class TH1F;
class TH2F;
class TH2D;
class TH3F;
class TH3D;
class TF1;
class TProfile;
class TProfile2D;
class TProfile3D;
class TEfficiency;
class CentralityMaker;
class StRefMultCorr;
class my_event;
class MixedBuffer;

#define USE_P 0
#define num_mult_bin 7
#define num_vz_bin 14
#define num_EP_bin 10
#define num_RecoPar 7
#define num_IdPar 6

class StKFParticleAnalysisMaker : public StMaker
{
public:
	StKFParticleAnalysisMaker(const char *name, const char *outName, int sys);
	virtual ~StKFParticleAnalysisMaker();

	virtual Int_t Init();
	virtual Int_t Make();
	virtual void Clear(Option_t *opt = "");
	virtual Int_t Finish();

	void setRunEnergyAndListDir(int run, double energy, char ListDir[256]);

private:
	// KFParticle
	StKFParticleInterface *KFParticleInterface;
	StKFParticlePerformanceInterface *KFParticlePerformanceInterface;
	void SetupKFParticle();
	void SetDaughterTrackPointers(int iKFParticle);
	bool IsKaonOmegaDaughter(KFParticle particle, int kaonTrackId);
	void CutDecider(KFParticle Omega, TH1D *hist_signal, TH1D *hist_sideband, double value);
	bool InterfaceCantProcessEvent;
	int ProtonTrackIndex, PionTrackIndex;
	vector<int> trackMap;
	StPicoDst *PicoDst;
	StPicoTrack *ProtonTrack, *PionTrack;
	void BookVertexPlots();

	StPicoDstMaker *mPicoDstMaker;
	StRefMultCorr *mRefMultCorr;

	int sys_tag;
	bool UseParticipant; // for first-order EPD plane, 2nd order always participant
	int mRun;
	int mRunStart;
	int mRunEnd;
	int mDay; // # of days used in EP correction
	double mEnergy;
	TString mListDir;

	TString mOutName;
	double PI;
	double twoPI;

	int mJob;
	// static const float OmegaMassSigma[7];
	// static const float OmegabarMassSigma[7];
	// static const float OmegaMassPtLowerBin[11];
	// static const float LambdaMassSigma[9];
	// static const float LambdabarMassSigma[9];

	// cut params for coalescence
	float pT_lo, pT_hi;
	float pT_asso_lo, pT_asso_hi;
	float pT_trig_lo, pT_trig_hi;
	float eta_trig_cut;
	float y_coal_cut;
	float y_phi_cut;
	float pion_pT_lo, pion_pT_hi;
	float proton_pT_lo, proton_pT_hi;
	float pion_pT_TOFth; // threshold above which TOF becomes required
	float proton_pT_TOFth;
	float kaon_pT_TOFth;
	float pion_m2_lo, pion_m2_hi;
	float proton_m2_lo, proton_m2_hi;
	float kaon_m2_lo, kaon_m2_hi;
	float dcatoPV_hi;

	// cut for correlation
	float lambda_dca_hi; // DCA between Lambda and PV

	// for correlation
	int min_cent;
	int max_cent;

	// for event mixing
	std::vector<float> mix_px;
	std::vector<float> mix_py;
	std::vector<float> mix_pz;
	std::vector<int> mix_charge;
	int mix_evt_id, mix_run_id;
	bool PerformMixing;
	bool CheckWeights2D; // whether to check EPD and TPC weights
	bool StoringTree;
	bool CutCent;
	bool PtReweighting;

	// v2, EPD stuff
	TString v2EPDMethod; // "Gang" or "2nd"
	bool v2Calculation;
	bool UsePtForPID;
	bool Coal2D; // if true, calculate v2 for coalescence in 2D
	// bool usePorPt; // 0 for p, 1 for pT
	/* Set up StEpdEpFinder */
	TClonesArray *mEpdHits;
	TChain *chain;
	StEpdEpFinder *mEpFinder;

	////////////////
	TH1F *hNRefMult;
	TH1F *hNRefMultA;
	TH1F *hNRefMultB;
	TH2F *hVertexXY;
	TH1F *hVertexZ;
	TH2F *hVertex2D;
	TH1F *hDiffVz;
	TH1F *hcent;
	TH1F *hcentw;

	TProfile *hcentRefM;
	TProfile *hcentRefW;
	TH1D *hRefMultCorr_cent[9];
	TProfile *hMultRatioEW_VertexZ;

	// xiatong's analysis
	bool PerformAnalysis;
	// Omega
	TH1D *hCorrKO[4];
	TH1D *hPtCorrKO[4];
	TH1D *hyCorrKO[4];
	TH1D *hphiCorrKO[4];

	TH1D *hCorrKO_same[4]; // for event-mixing, no efficiency correction
	TH1D *hCorrKO_mixed[4];

	TH1D *hCorrKO_sideband[4];
	TH1D *hPtCorrKO_sideband[4];
	TH1D *hyCorrKO_sideband[4];
	TH1D *hphiCorrKO_sideband[4];

	TH2D *hCorrKO_2D_same[4];
	TH2D *hCorrKO_2D_mixed[4];
	TH2D *hCorrLO_2D_same[4];
	TH2D *hCorrLO_2D_mixed[4];

	// Lambda
	TH1D *hCorrKL[4];
	TH1D *hPtCorrKL[4];
	TH1D *hyCorrKL[4];
	TH1D *hphiCorrKL[4];

	TH1D *hCorrKL_same[4]; // for event-mixing, no efficiency correction
	TH1D *hCorrKL_mixed[4];

	TH1D *hCorrKL_sideband[4];
	TH1D *hPtCorrKL_sideband[4];
	TH1D *hyCorrKL_sideband[4];
	TH1D *hphiCorrKL_sideband[4];

	// Proton
	TH1D *hCorrKP[4];
	TH1D *hPtCorrKP[4];
	TH1D *hyCorrKP[4];
	TH1D *hphiCorrKP[4];

	TH1D *hCorrKP_same[4]; // for event-mixing, no efficiency correction
	TH1D *hCorrKP_mixed[4];

	// some QA about tracing back to primary vertex
	TH1D *hOmegaDCAtoPV;
	TH1D *hOmegabarDCAtoPV;
	TH1D *hLambdaDCAtoPV;
	TH1D *hLambdabarDCAtoPV;
	TH1D *hOmegaPDiff; // difference between momentum at decay vertex and DCA to PV (with or without tracing back)
	TH1D *hOmegabarPDiff;
	TH1D *hLambdaPDiff;
	TH1D *hLambdabarPDiff;

	TH1D *hNegPtDiff_dphi_KmOb, *hPosPtDiff_dphi_KmOb, *hNegPtDiff_dphi_KmO, *hPosPtDiff_dphi_KmO;
	// TH2D *hCorrKplusO_y_pT, *hCorrKplusObar_y_pT, *hCorrKminusO_y_pT, *hCorrKminusObar_y_pT;
	// TH2D *hCorrKplusO_y_pT_mixed, *hCorrKplusObar_y_pT_mixed, *hCorrKminusO_y_pT_mixed, *hCorrKminusObar_y_pT_mixed;
	// TH2D *hCorrKplusO_y_phi, *hCorrKplusObar_y_phi, *hCorrKminusO_y_phi, *hCorrKminusObar_y_phi;
	// TH2D *hCorrKplusO_phi_pT, *hCorrKplusObar_phi_pT, *hCorrKminusO_phi_pT, *hCorrKminusObar_phi_pT;

	// lambda
	TH1D *hCorrLO[4];
	TH1D *hyCorrLO[4];
	TH1D *hphiCorrLO[4];
	TH1D *hCorrLO_sideband[4];
	TH1D *hyCorrLO_sideband[4];
	TH1D *hphiCorrLO_sideband[4];
	TH1D *hCorrLO_same[4];
	TH1D *hCorrLO_mixed[4];

	// a new test observable
	// kaon ratios at different p/pbar bins, in three different scenarios: with one omega, without o/ob, with one omegabar
	TProfile *hKratio_omega;
	TProfile *hKratio_wo;
	TProfile *hKratio_omegabar;

	// xiatong's QA
	TH1F *hEventQA;
	TH2D *hgpdEdx;
	TH1D *hgdEdxErr;
	TH2D *hgpinvbeta;
	TH1D *hgm2;
	TH2D *hgpm2;
	TH2D *hgptm2;
	TH2D *hgm2nSigmaKaon;
	TH2D *hgm2nSigmaPion;
	TH2D *hgm2nSigmaProton;
	TH2D *hgptnSigmaKaon;
	TH2D *hgptnSigmaPion;
	TH2D *hgptnSigmaProton;
	TH1D *hgp[9];
	TH1D *hgp_TOF[9];
	TH1D *hgpT[9];
	TH1D *hgpT_TOF[9];
	TH2F *hgpTeta[9];
	TH2F *hgpTeta_TOF[9];
	TProfile *hDauProtonFirstPoint_lam_pt;
	TProfile *hDauProtonLastPoint_lam_pt;
	TProfile *hDauPionFirstPoint_lam_pt;
	TProfile *hDauPionLastPoint_lam_pt;

	// y-pT/nq coverage
	TH2F *hIdPar_y_ptnq[num_IdPar];
	TH2F *hRecoPar_y_ptnq[num_RecoPar];

	// eTOF test
	TH1D *hgetofcrossingY;
	TH1D *hgm2_etof;
	TH2F *hgpinvbeta_etof;
	TH1D *hgetofmatchflag;
	TH1D *hgeta_etof;

	// TOF eff
	TH1D *hgDCAtoPV;
	TH1D *hgbtofYlocal;
	TH2D *hgKpdEdx;
	TH2D *hgKpinvbeta;
	TH1D *hgKm2;
	TH2D *hgKpm2;
	TH1D *hgKp;
	TH1D *hgKpT;
	TH1D *hgKDCAtoPV;
	TH1D *hgKDCAtoO;
	TH2D *hgKpionpdEdx;
	TH2D *hgKptnSigma;
	// TH2D *hgptm2_largenSigmaKaon;
	// TH2D *hgptm2_smallnSigmaKaon;
	TH1D *hgnSigmaDiff;
#if !USE_P
	// TH2F* hgPID2D_proton_pt[10];
	// TH2F* hgPID2D_antiproton_pt[10];
	// TH2F* hgPID2D_piplus_pt[10];
	// TH2F* hgPID2D_piminus_pt[10];
	// TH2F* hgPID2D_kplus_pt[10];
	// TH2F* hgPID2D_kminus_pt[10];

	TH2F *hgPID2D_proton_p[12];
	TH2F *hgPID2D_antiproton_p[12];
	TH2F *hgPID2D_piplus_p[12];
	TH2F *hgPID2D_piminus_p[12];
	TH2F *hgPID2D_kplus_p[12];
	TH2F *hgPID2D_kminus_p[12];
#else
	TH2F *hgPID2D_proton_p[12];
	TH2F *hgPID2D_antiproton_p[12];
	TH2F *hgPID2D_piplus_p[12];
	TH2F *hgPID2D_piminus_p[12];
	TH2F *hgPID2D_kplus_p[12];
	TH2F *hgPID2D_kminus_p[12];
#endif

	// TH1D* hgzTPC_pt[15];
	TProfile *hKaonCt;
	TH1D *hKpluspt_omega;
	TH1D *hKpluspt_omegabar;
	TH1D *hKminuspt_omega;
	TH1D *hKminuspt_omegabar;
	TH1D *hKpluseta_omega;
	TH1D *hKpluseta_omegabar;
	TH1D *hKminuseta_omega;
	TH1D *hKminuseta_omegabar;
	TH1D *hKplusphi_omega;
	TH1D *hKplusphi_omegabar;
	TH1D *hKminusphi_omega;
	TH1D *hKminusphi_omegabar;
	TH1D *hKplusy_omega;
	TH1D *hKplusy_omegabar;
	TH1D *hKminusy_omega;
	TH1D *hKminusy_omegabar;

	TH1D *hKpluspt_omega_sideband;
	TH1D *hKpluspt_omegabar_sideband;
	TH1D *hKminuspt_omega_sideband;
	TH1D *hKminuspt_omegabar_sideband;
	TH1D *hKpluseta_omega_sideband;
	TH1D *hKpluseta_omegabar_sideband;
	TH1D *hKminuseta_omega_sideband;
	TH1D *hKminuseta_omegabar_sideband;
	TH1D *hKplusphi_omega_sideband;
	TH1D *hKplusphi_omegabar_sideband;
	TH1D *hKminusphi_omega_sideband;
	TH1D *hKminusphi_omegabar_sideband;
	TH1D *hKplusy_omega_sideband;
	TH1D *hKplusy_omegabar_sideband;
	TH1D *hKminusy_omega_sideband;
	TH1D *hKminusy_omegabar_sideband;

	TH1D *hKpluspt_wlb_omega;
	TH1D *hKpluspt_wlb_omegabar;
	TH1D *hKminuspt_wl_omega;
	TH1D *hKminuspt_wl_omegabar;
	TH1D *hKplusy_wlb_omega;
	TH1D *hKplusy_wlb_omegabar;
	TH1D *hKminusy_wl_omega;
	TH1D *hKminusy_wl_omegabar;
	TH1D *hKpluspt_wolb_omega;
	TH1D *hKpluspt_wolb_omegabar;
	TH1D *hKminuspt_wol_omega;
	TH1D *hKminuspt_wol_omegabar;
	TH1D *hKplusy_wolb_omega;
	TH1D *hKplusy_wolb_omegabar;
	TH1D *hKminusy_wol_omega;
	TH1D *hKminusy_wol_omegabar;
	TH1D *hKpluspt_wlb_omega_sideband;
	TH1D *hKpluspt_wlb_omegabar_sideband;
	TH1D *hKminuspt_wl_omega_sideband;
	TH1D *hKminuspt_wl_omegabar_sideband;
	TH1D *hKplusy_wlb_omega_sideband;
	TH1D *hKplusy_wlb_omegabar_sideband;
	TH1D *hKminusy_wl_omega_sideband;
	TH1D *hKminusy_wl_omegabar_sideband;
	TH1D *hKpluspt_wolb_omega_sideband;
	TH1D *hKpluspt_wolb_omegabar_sideband;
	TH1D *hKminuspt_wol_omega_sideband;
	TH1D *hKminuspt_wol_omegabar_sideband;
	TH1D *hKplusy_wolb_omega_sideband;
	TH1D *hKplusy_wolb_omegabar_sideband;
	TH1D *hKminusy_wol_omega_sideband;
	TH1D *hKminusy_wol_omegabar_sideband;

	// for a rapidity check
	TProfile *hkplus_totaly;
	TProfile *hkminus_totaly;
	TProfile *hkplus_ntrack;
	TProfile *hkminus_ntrack;

	// these are for bTOF
	TH1D *hm2proton_b;		 // before ProtonPID.h cut
	TH1D *hm2proton_r;		 // after regular cut
	TH1D *hm2proton_a;		 // after
	TH1D *hm2pion_b;		 // before PionPID.h cut
	TH1D *hm2pion_r;		 // after regular cut
	TH1D *hm2pion_a;		 // after
	TProfile *hTOFEff_check; // confirming TOF eff correctness
	TProfile *hTPCEff_check[9];

	// v2 and EP
	TH1D *hEPD_e_EP_1, *hEPD_w_EP_1;
	TH1D *hEPD_e_EP_2, *hEPD_w_EP_2;
	TH1D *hEPD_full_EP_1, *hEPD_full_EP_2;
	TProfile *hEPD_ew_cos;
	// 7 reco particles
	// Omega, Omegabar, Xi, Xibar, Lambda, Lambdabar, phi
	TH1D *hRecoParM_cen[num_RecoPar][9]; // inv mass
	TProfile *hRecoPar_EPD_v2[num_RecoPar][9]; // v2 vs. inv mass
	TProfile *hRecoPar_TPC_v2[num_RecoPar][9];
	// 6 identified particles
	// piplus, piminus, proton, antiproton, kplus, kminus
	TProfile *hIdPar_EPD_v2[num_IdPar]; // v2 vs. centrality
	TProfile *hIdPar_TPC_v2[num_IdPar];
	TProfile *hIdPar_EPD_v2_pt_1st[num_IdPar][9]; // 1st order EP
	TProfile *hIdPar_EPD_v1_pt_1st[num_IdPar][9]; 
	TProfile *hIdPar_EPD_v2_pt_2nd[num_IdPar][9]; // 2nd order EP
	TProfile *hIdPar_EPD_v3_pt_1st[num_IdPar][9];
	TProfile *hIdPar_TPC_v2_pt[num_IdPar][9];
	// TProfile2D *hIdPar_EPD_v2_y_pt[num_IdPar][9]; // v2 vs. y and pT
	// TProfile2D *hIdPar_TPC_v2_y_pt[num_IdPar][9];
	TProfile2D *hIdPar_EPD_v1_eta_pt_1st[num_IdPar][9];
	TProfile2D *hIdPar_EPD_v2_eta_pt_1st[num_IdPar][9];
	TProfile2D *hIdPar_EPD_v3_eta_pt_1st[num_IdPar][9]; // for asynchronous analysis
	TProfile *hIdPar_EPD_v1_y[num_IdPar][9]; // v1 vs. y

	// QA: Lambda daughters rapidity spread, including both Lambda and Lambdabar
	// (for east Lambda, what is the rapidity spread, etc)
	TH1D *hyLambdaDauProton_east;
	TH1D *hyLambdaDauPion_east;
	TH1D *hyLambdaDauProton_west;
	TH1D *hyLambdaDauPion_west;

	// excited states
	TH1D *hXi1530M_cen[9];
	TH1D *hXi1530barM_cen[9];
	TH1D *hKsM_cen[9];
	TH1D *hOmega2012M_cen[9];
	TH1D *hOmega2012barM_cen[9];

	// mixed QA
	TH1D *hNumMixedEvent;
	TH1D *hTotalMixedEvent;

	TH1D *hProtony;
	TH1D *hAntiProtony;
	TH1D *hOmegaM;
	TH1D *hOmegap;
	TH1D *hOmegapt;
	TH1D *hOmegaeta;
	TH1D *hOmegay;
	TH2D *hOmegaypt;
	TH1D *hOmegaphi;
	TH1D *hOmegaDL;
	TH1D *hOmegaDLminusLambdaDL;

	TH1D *hOmegabarM;
	TH1D *hOmegabarp;
	TH1D *hOmegabarpt;
	TH1D *hOmegabareta;
	TH1D *hOmegabary;
	TH2D *hOmegabarypt;
	TH1D *hOmegabarphi;
	TH1D *hOmegabarDL;
	TH1D *hOmegabarDLminusLambdaDL;

	TH1D *hNumOmega;
	TH1D *hOmegaUsed;
	TH1D *hOmegaUsed_sideband;
	TH1D *hOmegaEventUsed;
	TH1D *hOmegaEventUsed_sideband;
	TH1D *hOmegaUsed_wlb;
	TH1D *hOmegaUsed_wlb_sideband;
	TH1D *hOmegaUsed_wolb;
	TH1D *hOmegaUsed_wolb_sideband;
	TH1D *hOmegaUsed_wl;
	TH1D *hOmegaUsed_wl_sideband;
	TH1D *hOmegaUsed_wol;
	TH1D *hOmegaUsed_wol_sideband;

	// for Lambda
	TH1D *hLambdaUsed;
	TH1D *hLambdaUsed_sideband;
	TH1D *hLambdaEventUsed;
	TH1D *hLambdaEventUsed_sideband;
	// for Proton
	TH1D *hProtonEventUsed;
	TH1D *hProtonEventUsed_sideband;

	// corresponding Omega inv mass
	TH1D *hOmegaM_wlb;
	TH1D *hOmegaM_wolb;
	TH1D *hOmegaM_wl;
	TH1D *hOmegaM_wol;
	TH1D *hOmegabarM_wlb;
	TH1D *hOmegabarM_wolb;
	TH1D *hOmegabarM_wl;
	TH1D *hOmegabarM_wol;

	// pT spectrum
	TH1D *hOmegaM_error_1;
	TH1D *hOmegabarM_error_1;
	TH1D *hOmegaM_error_2;
	TH1D *hOmegabarM_error_2;
	TH1D *hDauKaonM_error;
	TH1D *hDauLambdaM_error;
	// TH1D* hOmegaM_pt[9][10];
	// TH1D* hOmegabarM_pt[9][10];
	// TH1D* hOmegaM_phi[9][10]; // if testing v2 vs pT, change centrality bins (9) to pT bins (10)
	// TH1D* hOmegabarM_phi[9][10];
	// TH1D* hOmegaM_rotbkg_pi_pt[10];
	// TH1D* hOmegabarM_rotbkg_pi_pt[10];

	TH1D *hLambdaM;
	TH1D *hLambdap;
	TH1D *hLambdapt;
	TH1D *hLambday;
	TH1D *hLambdaphi;
	TH1D *hLambdaDL;

	TH1D *hLambdabarM;
	TH1D *hLambdabarp;
	TH1D *hLambdabarpt;
	TH1D *hLambdabary;
	TH1D *hLambdabarphi;
	TH1D *hLambdabarDL;

	TH1D *hDauKplusp;
	TH1D *hDauKpluspt;
	TH1D *hDauKplusy;
	TH1D *hDauKplusphi;
	TH1D *hDauKplusnSigma;
	TH1D *hDauKminusp;
	TH1D *hDauKminuspt;
	TH1D *hDauKminusy;
	TH1D *hDauKminusphi;
	TH1D *hDauKminusnSigma;

	TH1D *hDauLambdaM;
	TH1D *hDauLambdap;
	TH1D *hDauLambdapt;
	TH1D *hDauLambday;
	TH1D *hDauLambdaphi;
	TH1D *hDauLambdaDL;
	TH1D *hDauLambdabarM;
	TH1D *hDauLambdabarp;
	TH1D *hDauLambdabarpt;
	TH1D *hDauLambdabary;
	TH1D *hDauLambdabarphi;
	TH1D *hDauLambdabarDL;

	TH1D *hDauProtonpt;
	TH1D *hDauProtonp;
	TH1D *hDauPionpt;
	TH1D *hDauPionp;
	TH1D *hDauProtonbarpt;
	TH1D *hDauProtonbarp;
	TH1D *hDauPionbarpt;
	TH1D *hDauPionbarp;

	TH1D *hOmegaDauPid;
	TH1D *hOmegabarDauPid;

	TH1D *hDCAOtoK_signal;
	TH1D *hDCAOtoK_sideband;
	TH1D *hDCAOtoL_signal;
	TH1D *hDCAOtoL_sideband;

	// Run-by-run QA
	TProfile *hEPD_full_1_runID;
	TProfile *hEPD_full_2_runID;
	TProfile *hEPD_full_1_day;
	TProfile *hEPD_full_2_day;

	// mixed event buffer
	MixedBuffer buffer;
	TFile *ftree;
	TTree *omega_mix[num_mult_bin][num_vz_bin][num_EP_bin];

	// TPC weights
	TFile *fTPCShift;
	TH1D *hTPCAssoPhi;
	TH1D *hTPCAssoPhi_shifted;
	TH2D *hTPCAssoPhi_2D;
	TH2D *hTPCAssoPhi_2D_shifted;
	TH1D *hTPCPOIPhi;
	TH1D *hTPCPOIPhi_shifted;
	TH2D *hTPCPOIPhi_2D;
	TH2D *hTPCPOIPhi_2D_shifted;
	TProfile3D *hTPCAssoShiftInput_sin;
	TProfile3D *hTPCAssoShiftInput_cos;
	TProfile3D *hTPCPOIShiftInput_sin;
	TProfile3D *hTPCPOIShiftInput_cos;
	TProfile3D *hTPCEPShiftInput_cos[3];
	TProfile3D *hTPCEPShiftInput_sin[3];
	TProfile3D *hTPCAssoShiftOutput_sin;
	TProfile3D *hTPCAssoShiftOutput_cos;
	TProfile3D *hTPCPOIShiftOutput_sin;
	TProfile3D *hTPCPOIShiftOutput_cos;
	TProfile3D *hTPCEPShiftOutput_cos[3];
	TProfile3D *hTPCEPShiftOutput_sin[3];
	TH1D *hTPCEP_2[9][3];
	TH1D *hTPCEP_2_shifted[9][3];
	TH2D *hTPCEP_2_2D[9][3];
	TH2D *hTPCEP_2_2D_shifted[9][3];
	TProfile *hTPCEP_2_shift[9][3];
	TProfile *hTPCEP_ew_cos;

	// EPD weights
	TFile *fEPDShift;
	TProfile3D *hEPDEPShiftInput_1_cos[3];
	TProfile3D *hEPDEPShiftInput_1_sin[3];
	TProfile3D *hEPDEPShiftOutput_1_cos[3];
	TProfile3D *hEPDEPShiftOutput_1_sin[3];
	TProfile3D *hEPDEPShiftInput_2_cos[3];
	TProfile3D *hEPDEPShiftInput_2_sin[3];
	TProfile3D *hEPDEPShiftOutput_2_cos[3];
	TProfile3D *hEPDEPShiftOutput_2_sin[3];
	TH1D *hEPDEP_1[9][3];
	TH1D *hEPDEP_2[9][3];
	TH1D *hEPDEP_1_shifted[9][3];
	TH1D *hEPDEP_2_shifted[9][3];
	TH2D *hEPDEP_1_2D[9][3];
	TH2D *hEPDEP_2_2D[9][3];
	TH2D *hEPDEP_1_2D_shifted[9][3];
	TH2D *hEPDEP_2_2D_shifted[9][3];
	TProfile *hEPDEP_2_shift[9][3];
	TProfile *hEPDEP_ew_cos_1;
	TProfile *hEPDEP_ew_cos_2;
	TProfile *hEPDEP_ew_cos_3;
	TH1D *hEPDEP_2_Full_Raw_QA_cen4;
	TH1D *hEPDEP_2_Full_PhiWeighted_QA_cen4;
	TH1D *hEPDEP_2_Full_PhiWeightedAndShifted_QA_cen4;

	// Checking east west asymmetry
	TH1D *hPhi_EW_PN_EWVtx[2][2][2];
	TH1D *hPhiXi_EW_PN_EWVtx[2][2][2];

	// TOF Efficiency
	TFile *fTOFEff;
	TEfficiency *hTOFEff[9];
	TH2D *hTOFEff_2D[9];

	/////////////////////////////////////
	int mStps;

	std::vector<Int_t> badList;
	std::vector<Int_t> runList;

	Int_t findCentrality(int mult);
	Int_t CheckrunNumber(int runnumber);
	bool readRunList();
	bool readBadList();
	bool removeBadID(int runnumber);
	Int_t openFile();

	void DeclareHistograms();
	void DeclareTrees();
	void WriteHistograms();
	void WriteTPCShift();
	void WriteEPDShift();
	void WriteTrees();
	void ReadTrees();

	int MixRefMultBin(int cent, int refmult);
	bool isGoodOmega(int cent, KFParticle Omega);
	bool isSidebandOmega(int cent, KFParticle Omega);
	bool isGoodLambda(int cent, KFParticle Lambda);
	bool isSidebandLambda(int cent, KFParticle Lambda);
	bool isGoodObs(double obs);
	float GetPtWeight(KFParticle Omega);
	float Eta2y(float pt, float eta, float mass);
	float ShiftPOIPhi(float phi, int day, int cent);
	float ShiftAssoPhi(float phi, int day, int cent);
	float ShiftTPCEP(float psi, int day, int cent, int ewFull, int order);
	float ShiftEPDEP(float psi, int day, int cent, int ewFull, int order);

	ClassDef(StKFParticleAnalysisMaker, 1)
};

#endif
