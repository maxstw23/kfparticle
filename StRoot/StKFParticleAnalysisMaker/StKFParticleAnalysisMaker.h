#ifndef StKFParticleAnalysisMaker_h
#define StKFParticleAnalysisMaker_h

#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"
#include "StRoot/StEpdUtil/StEpdEpFinder.h"

#include "StKFParticleInterface.h"
#include "StKFParticlePerformanceInterface.h"

#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"
#include "StPicoEvent/StPicoPhysicalHelix.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "StMaker.h"
#include "TString.h"
#include "TObject.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
#include <cmath>
#include <vector>
#include <algorithm>
#include <iostream>
#include "TF1.h"
#include "StTrackHelix.h"
#include "my_event.h"
#include "MixedBuffer.h"
#include "MyConstant.h" // must include

class StPicoDst;
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
	// KFParticle
	StKFParticleInterface *KFParticleInterface;
	StKFParticlePerformanceInterface *KFParticlePerformanceInterface;
	void SetupKFParticle();
	void SetDaughterTrackPointers(int iKFParticle);
	bool IsKaonOmegaDaughter(KFParticle particle, int kaonTrackId);
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
	int 	   nOmegaEvtProcessed;

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

	// xiatong's analysis
	TH1D *hCorrKplusO, *hCorrKplusObar, *hCorrKminusO, *hCorrKminusObar;
	TH1D *hCorrKplusO_mixed, *hCorrKplusObar_mixed, *hCorrKminusO_mixed, *hCorrKminusObar_mixed; 

	// xiatong's QA
	TH2D *hgpdEdx     ;
	TH1D *hgdEdxErr   ;
	TH2D *hgpinvbeta  ;
	TH1D *hgm2        ; 
	TH2D *hgpm2       ;
	TH2D *hgm2nSigmaKaon;
	TH2D *hgm2nSigmaPion;
	TH2D *hgm2nSigmaProton;
	TH1D *hgp         ; 
	TH1D *hgpT        ;
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

	TH1D *hOmegaM     ;
	TH1D *hOmegap     ;
	TH1D *hOmegapt    ;
	TH1D *hOmegay     ;
	TH2D *hOmegaypt   ;
	TH1D *hOmegaphi   ;
	TH1D *hOmegaDL    ;

	TH1D *hOmegabarM  ;
	TH1D *hOmegabarp  ;
	TH1D *hOmegabarpt ;
	TH1D *hOmegabary  ;
	TH2D *hOmegabarypt;
	TH1D *hOmegabarphi;
	TH1D *hOmegabarDL ;

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
	TH1D *hDauLambdabarM  ;
	TH1D *hDauLambdabarp  ;
	TH1D *hDauLambdabarpt ;
	TH1D *hDauLambdabary  ;
	TH1D *hDauLambdabarphi;

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
	void  WriteHistograms();

	bool isGoodObs(double obs);
		
	ClassDef(StKFParticleAnalysisMaker, 1)
};

#endif


