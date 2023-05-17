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

class StKFParticleAnalysisMaker : public StMaker 
{
public:
	StKFParticleAnalysisMaker(const char *name, const char *outName);
	virtual ~StKFParticleAnalysisMaker();

	virtual Int_t Init();
	virtual Int_t Make();
	virtual void  Clear(Option_t *opt="");
	virtual Int_t Finish();

	void    setRunEnergyAndListDir(int run,double energy,char ListDir[256]);            

private:
	// efficiency
	static const float P0_P[9];
	static const float P1_P[9];
	static const float P2_P[9];
	static const float P0_AP[9];
	static const float P1_AP[9];
	static const float P2_AP[9];
	static const float P0_pip[9];
	static const float P1_pip[9];
	static const float P2_pip[9];
	static const float P0_pim[9];
	static const float P1_pim[9];
	static const float P2_pim[9];
	static const float P0_Kp[9];
	static const float P1_Kp[9];
	static const float P2_Kp[9];
	static const float P0_Km[9];
	static const float P1_Km[9];
	static const float P2_Km[9];

	// KFParticle
	StKFParticleInterface *KFParticleInterface;
	StKFParticlePerformanceInterface *KFParticlePerformanceInterface;
	void SetupKFParticle();
	void SetDaughterTrackPointers(int iKFParticle);
	bool IsKaonOmegaDaughter(KFParticle particle, int kaonTrackId);
	void CutDecider(KFParticle Omega, TH1D* hist_signal, TH1D* hist_sideband, double value);
	bool InterfaceCantProcessEvent;
	int ProtonTrackIndex, PionTrackIndex;
	vector<int> trackMap;
	StPicoDst *PicoDst;
	StPicoTrack *ProtonTrack, *PionTrack;
	void BookVertexPlots();

	StPicoDstMaker *mPicoDstMaker;
	StRefMultCorr *mRefMultCorr;

	int        mRun;            
	double     mEnergy;            
	TString    mListDir;            

	TString    mOutName;
	double     PI;
	double     twoPI;

	int        mJob;
	static const float OmegaMassSigma[7];
	static const float OmegabarMassSigma[7];
	static const float OmegaMassPtLowerBin[11];
	static const float LambdaMassSigma[9];
	static const float LambdabarMassSigma[9];	

	// cut params for coalescence
	float pT_lo, pT_hi;
	float pion_pT_lo, pion_pT_hi;
	float proton_pT_lo, proton_pT_hi;
	float pion_pT_TOFth; // threshold above which TOF becomes required
	float proton_pT_TOFth;
	float pion_m2_lo, pion_m2_hi;
	float proton_m2_lo, proton_m2_hi;
	float dcatoPV_hi;

	// for event mixing
	std::vector<float> mix_px;
	std::vector<float> mix_py;
	std::vector<float> mix_pz;
	std::vector<int>   mix_charge;
	int mix_evt_id, mix_run_id;
	bool PerformMixing;
	bool StoringTree;
	bool CutCent;
	bool PtReweighting;

	// v2, EPD stuff
	bool v2Calculation;
	/* Set up StEpdEpFinder */
	TClonesArray* mEpdHits;
	TChain* chain;
	StEpdEpFinder* mEpFinder;

	////////////////
	TH1F *hNRefMult;
	TH1F *hNRefMultA;
	TH1F *hNRefMultB;
	TH2F *hVertexXY;
	TH1F *hVertexZ;
	TH2F *hVertex2D; 
	TH1F *hDiffVz  ; 
	TH1F *hcent;
	TH1F *hcentw;

	TProfile *hcentRefM ; 
	TProfile *hcentRefW ; 
	TH1D *hRefMultCorr_cent[9];

	// xiatong's analysis
	bool PerformAnalysis;
	TH1D *hCorrKplusO, *hCorrKplusObar, *hCorrKminusO, *hCorrKminusObar;
	TH1D *hPtCorrKplusO, *hPtCorrKplusObar, *hPtCorrKminusO, *hPtCorrKminusObar;
	TH1D *hyCorrKplusO, *hyCorrKplusObar, *hyCorrKminusO, *hyCorrKminusObar;
	TH1D *hphiCorrKplusO, *hphiCorrKplusObar, *hphiCorrKminusO, *hphiCorrKminusObar;

	TH1D *hCorrKplusO_mixed, *hCorrKplusObar_mixed, *hCorrKminusO_mixed, *hCorrKminusObar_mixed; 

	TH1D *hCorrKplusO_sideband, *hCorrKplusObar_sideband, *hCorrKminusO_sideband, *hCorrKminusObar_sideband;
	TH1D *hPtCorrKplusO_sideband, *hPtCorrKplusObar_sideband, *hPtCorrKminusO_sideband, *hPtCorrKminusObar_sideband;
	TH1D *hyCorrKplusO_sideband, *hyCorrKplusObar_sideband, *hyCorrKminusO_sideband, *hyCorrKminusObar_sideband;
	TH1D *hphiCorrKplusO_sideband, *hphiCorrKplusObar_sideband, *hphiCorrKminusO_sideband, *hphiCorrKminusObar_sideband;

	TH1D *hNegPtDiff_dphi_KmOb, *hPosPtDiff_dphi_KmOb, *hNegPtDiff_dphi_KmO, *hPosPtDiff_dphi_KmO;
	// TH2D *hCorrKplusO_y_pT, *hCorrKplusObar_y_pT, *hCorrKminusO_y_pT, *hCorrKminusObar_y_pT;
	TH2D *hCorrKplusO_y_pT_mixed, *hCorrKplusObar_y_pT_mixed, *hCorrKminusO_y_pT_mixed, *hCorrKminusObar_y_pT_mixed;
	// TH2D *hCorrKplusO_y_phi, *hCorrKplusObar_y_phi, *hCorrKminusO_y_phi, *hCorrKminusObar_y_phi;
	// TH2D *hCorrKplusO_phi_pT, *hCorrKplusObar_phi_pT, *hCorrKminusO_phi_pT, *hCorrKminusObar_phi_pT;

	// a new test observable
    // kaon ratios at different p/pbar bins, in three different scenarios: with one omega, without o/ob, with one omegabar
    TProfile* hKratio_omega   ;  
    TProfile* hKratio_wo      ;  
    TProfile* hKratio_omegabar; 

	// xiatong's QA
	TH1F *hEventQA    ;
	TH2D *hgpdEdx     ;
	TH1D *hgdEdxErr   ;
	TH2D *hgpinvbeta  ;
	TH1D *hgm2        ; 
	TH2D *hgpm2       ;
	TH2D *hgptm2      ;
	TH2D *hgm2nSigmaKaon;
	TH2D *hgm2nSigmaPion;
	TH2D *hgm2nSigmaProton;
	TH2D *hgptnSigmaKaon;
	TH2D *hgptnSigmaPion;
	TH2D *hgptnSigmaProton;
	TH1D *hgp         ; 
	TH1D *hgpT[9]     ;
	TH1D *hgpT_TOF[9] ;

	// TOF eff
	TH1D *hgDCAtoPV   ;
	TH1D *hgbtofYlocal;
	TH2D *hgKpdEdx    ;
	TH2D *hgKpinvbeta ;
	TH1D *hgKm2       ;   
	TH2D *hgKpm2      ;   
	TH1D *hgKp        ;   
	TH1D *hgKpT       ;
	TH1D *hgKDCAtoPV  ;   
	TH1D *hgKDCAtoO   ; 
	TH2D *hgKpionpdEdx;
	TH2D *hgKptnSigma ;
	// TH2D *hgptm2_largenSigmaKaon;
	// TH2D *hgptm2_smallnSigmaKaon;
	TH1D *hgnSigmaDiff;
	// TH2D* hgPID2D_proton_pt[10];
	// TH2D* hgPID2D_pion_pt[10];
	//TH1D* hgzTPC_pt[15];
	TProfile* hKaonCt;
	TH1D *hKpluspt_omega    ;  
	TH1D *hKpluspt_omegabar ;  
	TH1D *hKminuspt_omega   ;  
	TH1D *hKminuspt_omegabar;
	TH1D *hKpluseta_omega    ;  
	TH1D *hKpluseta_omegabar ;  
	TH1D *hKminuseta_omega   ;  
	TH1D *hKminuseta_omegabar;
	TH1D *hKplusphi_omega    ;  
	TH1D *hKplusphi_omegabar ;  
	TH1D *hKminusphi_omega   ;  
	TH1D *hKminusphi_omegabar;
	TH1D *hKplusy_omega    ;
	TH1D *hKplusy_omegabar ;
	TH1D *hKminusy_omega   ;
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

	TH1D *hm2proton_b; // before ProtonPID.h cut
	TH1D *hm2proton_r; // after regular cut
	TH1D *hm2proton_a; // after
	TH1D *hm2pion_b; // before PionPID.h cut
	TH1D *hm2pion_r; // after regular cut
	TH1D *hm2pion_a; // after
	TProfile *hTOFEff_check; // confirming TOF eff correctness
	TProfile *hTPCEff_check[9];

	// v2 and EP
	TH1D *hTPCEP_2[3]      ;
	TH1D *hTPCEP_2_shifted[3];
	TProfile *hTPCEP_ew_cos;
	TProfile *hOmega_TPC_v2_pt, *hOmegabar_TPC_v2_pt;
	TProfile *hOmega_EPD_v2_pt, *hOmegabar_EPD_v2_pt;
	TH1D *hEPD_e_EP_1, *hEPD_w_EP_1;
	TH1D *hEPD_e_EP_2, *hEPD_w_EP_2;
	TH1D *hEPD_full_EP_1, *hEPD_full_EP_2;
	TProfile *hEPD_ew_cos;
	TProfile *hOmega_EPD_v2[9];
	TProfile *hOmega_TPC_v2[9];
	TProfile *hOmegabar_EPD_v2[9];
	TProfile *hOmegabar_TPC_v2[9];
	TProfile *hpiplus_EPD_v2;
	TProfile *hpiminus_EPD_v2;
	TProfile *hproton_EPD_v2;
	TProfile *hantiproton_EPD_v2;
	TProfile *hkplus_EPD_v2;
	TProfile *hkminus_EPD_v2;
	TProfile *hpiplus_TPC_v2;
	TProfile *hpiminus_TPC_v2;
	TProfile *hproton_TPC_v2;
	TProfile *hantiproton_TPC_v2;
	TProfile *hkplus_TPC_v2;
	TProfile *hkminus_TPC_v2;
	TProfile *hpiplus_EPD_v2_pt[9];
	TProfile *hpiminus_EPD_v2_pt[9];
	TProfile *hproton_EPD_v2_pt[9];
	TProfile *hantiproton_EPD_v2_pt[9];
	TProfile *hkplus_EPD_v2_pt[9];
	TProfile *hkminus_EPD_v2_pt[9];
	TProfile *hpiplus_TPC_v2_pt[9];
	TProfile *hpiminus_TPC_v2_pt[9];
	TProfile *hproton_TPC_v2_pt[9];
	TProfile *hantiproton_TPC_v2_pt[9];
	TProfile *hkplus_TPC_v2_pt[9];
	TProfile *hkminus_TPC_v2_pt[9];
	
	// Xi
	TH1D *hXiM_cen[9];
	TH1D *hXibarM_cen[9];
	TProfile *hXi_EPD_v2[9];
	TProfile *hXi_TPC_v2[9];
	TProfile *hXibar_EPD_v2[9];
	TProfile *hXibar_TPC_v2[9];

	// Lambda
	TH1D *hLambdaM_cen_pt[9][40];
	TH1D *hLambdabarM_cen_pt[9][40];
	TProfile *hLambda_EPD_v2_pt[9][40];
	TProfile *hLambda_TPC_v2_pt[9][40];
	TProfile *hLambdabar_EPD_v2_pt[9][40];
	TProfile *hLambdabar_TPC_v2_pt[9][40];

	// mixed QA
	TH1D* hNumMixedEvent;
	TH1D* hTotalMixedEvent;

	TH1D* hProtony;
	TH1D* hAntiProtony;
	TH1D *hOmegaM     ;
	TH1D *hOmegaM_cen[9];
	TH1D *hOmegap     ;
	TH1D *hOmegapt    ;
	TH1D *hOmegaeta   ;
	TH1D *hOmegay     ;
	TH2D *hOmegaypt   ;
	TH1D *hOmegaphi   ;
	TH1D *hOmegaDL    ;
	TH1D *hOmegaDLminusLambdaDL;

	TH1D *hOmegabarM  ;
	TH1D *hOmegabarM_cen[9];
	TH1D *hOmegabarp  ;
	TH1D *hOmegabarpt ;
	TH1D *hOmegabareta;
	TH1D *hOmegabary  ;
	TH2D *hOmegabarypt;
	TH1D *hOmegabarphi;
	TH1D *hOmegabarDL ;
	TH1D *hOmegabarDLminusLambdaDL;

	TH1D *hNumOmega;
	TH1D *hOmegaUsed;
	TH1D *hOmegaUsed_sideband;
	TH1D *hOmegaUsed_wlb;
	TH1D *hOmegaUsed_wlb_sideband;
	TH1D *hOmegaUsed_wolb;
	TH1D *hOmegaUsed_wolb_sideband;
	TH1D *hOmegaUsed_wl;
	TH1D *hOmegaUsed_wl_sideband;
	TH1D *hOmegaUsed_wol;
	TH1D *hOmegaUsed_wol_sideband;

	// pT spectrum
	TH1D* hOmegaM_error_1;
	TH1D* hOmegabarM_error_1;
	TH1D* hOmegaM_error_2;
	TH1D* hOmegabarM_error_2;
	TH1D* hDauKaonM_error;
	TH1D* hDauLambdaM_error;
	// TH1D* hOmegaM_pt[9][10]; 
	// TH1D* hOmegabarM_pt[9][10];
	TH1D* hOmegaM_phi[9][10]; // if testing v2 vs pT, change centrality bins (9) to pT bins (10)
	TH1D* hOmegabarM_phi[9][10];
	//TH1D* hOmegaM_rotbkg_pi_pt[10];
	//TH1D* hOmegabarM_rotbkg_pi_pt[10];

	TH1D *hLambdaM     ;
	TH1D *hLambdap     ;
	TH1D *hLambdapt    ;
	TH1D *hLambday     ;
	TH1D *hLambdaphi   ;
	TH1D *hLambdaDL    ;

	TH1D *hLambdabarM  ;
	TH1D *hLambdabarp  ;
	TH1D *hLambdabarpt ;
	TH1D *hLambdabary  ;
	TH1D *hLambdabarphi;
	TH1D *hLambdabarDL ;

	TH1D *hDauKplusp     ;
	TH1D *hDauKpluspt    ;
	TH1D *hDauKplusy     ;
	TH1D *hDauKplusphi   ;
	TH1D *hDauKplusnSigma;
	TH1D *hDauKminusp     ;
	TH1D *hDauKminuspt    ;
	TH1D *hDauKminusy     ;
	TH1D *hDauKminusphi   ;
	TH1D *hDauKminusnSigma;

	TH1D *hDauLambdaM     ;
	TH1D *hDauLambdap     ;
	TH1D *hDauLambdapt    ;
	TH1D *hDauLambday     ;
	TH1D *hDauLambdaphi   ;
	TH1D *hDauLambdaDL    ;
	TH1D *hDauLambdabarM  ;
	TH1D *hDauLambdabarp  ;
	TH1D *hDauLambdabarpt ;
	TH1D *hDauLambdabary  ;
	TH1D *hDauLambdabarphi;
	TH1D *hDauLambdabarDL ;

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

	TH1D* hDCAOtoK_signal;
	TH1D* hDCAOtoK_sideband;
	TH1D* hDCAOtoL_signal;
	TH1D* hDCAOtoL_sideband;


	// mixed event buffer
	MixedBuffer buffer;
	TFile *ftree;
	TTree *omega_mix[5][16][6]; 

	// TPC weights
	TFile *fTPCShift;
	TH1D *hTPCAssoPhi;
	TH1D *hTPCAssoPhi_shifted;
	TH1D *hTPCPOIPhi;
	TH1D *hTPCPOIPhi_shifted;
	TProfile2D *hTPCAssoShiftInput_sin;
	TProfile2D *hTPCAssoShiftInput_cos;
	TProfile2D *hTPCPOIShiftInput_sin;
	TProfile2D *hTPCPOIShiftInput_cos;
	TProfile2D *hTPCEPShiftInput_cos[3];
	TProfile2D *hTPCEPShiftInput_sin[3];
	TProfile2D *hTPCAssoShiftOutput_sin;
	TProfile2D *hTPCAssoShiftOutput_cos;
	TProfile2D *hTPCPOIShiftOutput_sin;
	TProfile2D *hTPCPOIShiftOutput_cos;
	TProfile2D *hTPCEPShiftOutput_cos[3];
	TProfile2D *hTPCEPShiftOutput_sin[3];

	// TOF Efficiency
	TFile *fTOFEff;
	TEfficiency *hTOFEff[9];

	/////////////////////////////////////
	int mStps;  

	std::vector<Int_t> badList;
	std::vector<Int_t> runList;

	Int_t findCentrality(int mult);
	Int_t CheckrunNumber(int runnumber);            
	bool  readRunList();            
	bool  readBadList();            
	bool  removeBadID(int runnumber);            
	Int_t openFile();

	void  DeclareHistograms();
	void  DeclareTrees();
	void  WriteHistograms();
	void  WriteTPCShift();
	void  WriteTrees();
	void  ReadTrees();

	int MixRefMultBin(int cent, int refmult);
	bool isGoodOmega(int cent, KFParticle Omega);
	bool isSidebandOmega(int cent, KFParticle Omega);
	bool isGoodLambda(int cent, KFParticle Lambda);
	bool isGoodObs(double obs);
	float GetPtWeight(KFParticle Omega);
	float Eta2y(float pt, float eta, float mass);
	float ShiftPOIPhi(float phi, int cent);
	float ShiftAssoPhi(float phi, int cent);
		
	ClassDef(StKFParticleAnalysisMaker, 1)
};

#endif


