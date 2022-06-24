#ifndef MYTOOLKIT_H
#define MYTOOLKIT_H

#include "TVector3.h"
#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"
#include "StTrackHelix.h"


//------------------------------------------------------------
double getDcaToPV(StPicoTrack * track, const TVector3& pv, double magnet)
{
    // calculate the Dca to PV of a global track in Hui LONG's way.
    // first prepare the parameter of the helix.
    //StPhysicalHelixD helix = track->helix(magnet); //inner helix. good for dca to PV.
    StPicoPhysicalHelix helix = track->helix(magnet); //inner helix. good for dca to PV.
    int hrot = helix.h();
    float psic = helix.phase();
    TVector3 origin = helix.origin();
    float Xc = helix.xcenter() - pv.x();
    float Yc = helix.ycenter() - pv.y();
    float Zc = origin.z() - pv.z();

    float r  = 1./helix.curvature();
    float Vz = r*tan(helix.dipAngle());

    double PI = TMath::Pi();

    //then call Hui Long's func

    float f;
    float t1,t;
    float t2,a,b,c;
    if(fabs(Vz/r)<0.2){
	f=sqrt(Xc*Xc+Yc*Yc);
	t1=-Xc*r/f;
	t2=-Yc*r/f;
	a=r*cos(psic);
	b=r*sin(psic);
	c=acos((a*t1+b*t2)/r/r);
	t1=sqrt((f-r)*(f-r)+(Zc-c*Vz)*(Zc-c*Vz));
	c-=2*PI;
	t2=sqrt((f-r)*(f-r)+(Zc+c*Vz)*(Zc+c*Vz));
	if(t1<t2) return t1;
	else return t2;
    }

    //cout<<"test here: "<<endl;
    a=sqrt(Xc*Xc+Yc*Yc)*r/Vz/Vz;
    c=hrot*(psic-atan2(Yc,Xc));
    b=-Zc/Vz+c;
    int n1,n2;

    if(a>0.){
	if(b-a<0.) n1=(int)( (b-a-PI/2.)/PI);
	else n1=(int)( (b-a+PI/2.)/PI);
	if(b+a<0.) n2=(int)( (b+a-PI/2.)/PI);
	else n2=(int)( (b+a+PI/2.)/PI);
    }
    else{
	if(b+a<0.) n1=(int)( (b+a-PI/2.)/PI);
	else n1=(int)( (b+a+PI/2.)/PI);
	if(b-a<0.) n2=(int)( (b-a-PI/2.)/PI);
	else n2=(int)( (b-a+PI/2.)/PI);
    }
    float bound1,bound2,fb1,fb2,gb1,gb2;
    float save =1.e+10;

    for(int n=n1;n<=n2;n++){
	int i=0;
	bound1=n*PI-PI/2.0;
	bound2=n*PI+PI/2.0;
	fb1=a*sin(bound1);
	fb2=a*sin(bound2);
	gb1=bound1-b;
	gb2=bound2-b;
	if((gb1>fb1&&gb2>fb2)||(gb1<fb1&&gb2<fb2))continue;
	if(gb1>fb1&&gb2<fb2){
	    t1=bound1;
	    t2=bound2;
	}
	else{
	    t1=bound2;
	    t2=bound1;
	}
	while(fabs(t2-t1)>0.000001){
	    i++;
	    t=(t1+t2)/2.0;
	    if(a*sin(t)+b-t>0.) t2=t;
	    else t1=t;
	    if(i>1000)break;
	};
	f=t1-c;
	a=sqrt((Xc+r*cos(psic+f*hrot))*(Xc+r*cos(psic+f*hrot))+(Yc+r*sin(psic+f*hrot))*(Yc+r*sin(psic+f*hrot))+(Zc+f*Vz)*(Zc+f*Vz));
	if(a<save) save=a;
    }
    return save;
}


//------------------------------------------------------------
void  dcaPToPi(float *fi_p,float *fi_pi, StTrackHelix* proton,StTrackHelix* pion,float *d_root,float *v0,float alfa)
{

    float R2_TPC = 30000.; //hard-coded cuts
    float Z_TPC  = 180.  ;
    float t1,t2;

    float Co_a,Co_b,Co_c,Co_e,Co_f,dfi_p,dfi_pi;
    float new_d,x1,y1,z1,x2,y2,z2;


    Co_a=proton->r*proton->r+proton->Vz*proton->Vz;
    Co_b=proton->r*pion->r*cos(alfa)-proton->Vz*pion->Vz;  //XZHU: note the assumption here: helicities of two Helixes should be opposite
    Co_c=pion->r*pion->r+pion->Vz*pion->Vz;
    Co_e=(*d_root)*proton->Vz;
    Co_f=-(*d_root)*pion->Vz;
    dfi_p=0.;
    dfi_pi=0.;
    if(fabs(Co_a*Co_c-Co_b*Co_b)>0.){
	dfi_p=(Co_e*Co_c-Co_b*Co_f)/(Co_a*Co_c-Co_b*Co_b);

	dfi_pi=(Co_a*Co_f-Co_b*Co_e)/(Co_a*Co_c-Co_b*Co_b);
    }

    t1=*fi_p+dfi_p;
    t2=*fi_pi+dfi_pi;

    x1=proton->Xc+proton->r*cos(t1*proton->h+proton->theta);
    y1=proton->Yc+proton->r*sin(t1*proton->h+proton->theta);
    z1=proton->Zc+proton->Vz*t1;

    x2=pion->Xc+pion->r*cos(t2*pion->h+pion->theta);
    y2=pion->Yc+pion->r*sin(t2*pion->h+pion->theta);
    z2=pion->Zc+pion->Vz*t2;
    new_d=sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2));
    if(fabs(z1)>Z_TPC||fabs(z2)>Z_TPC)new_d=500.;
    if((x1*x1+y1*y1)>R2_TPC||(x2*x2+y2*y2)>R2_TPC)new_d=500.;

    if(new_d<fabs(*d_root)){
	*fi_p=t1*proton->h+proton->theta;
	*fi_pi=t2*pion->h+pion->theta;

	*d_root=new_d;

	*v0=(x2+x1)/2;
	*(v0+1)=(y2+y1)/2;
	*(v0+2)=(z2+z1)/2;
    }
    else{
	*fi_p=(t1-dfi_p)*proton->h+proton->theta;
	*fi_pi=(t2-dfi_pi)*pion->h+pion->theta;

	*d_root = fabs(*d_root); //XZHU: remove the negative dca bug

	*(v0+2)=((pion->Vz*(t2-dfi_pi)+pion->Zc)+(proton->Vz*(t1-dfi_p)+proton->Zc))/2.;
    }
    return;
}


//-----------------------------------------------------
double closestDistance(StPicoTrack * track1, StPicoTrack * track2, double magnet, const TVector3& pv, TVector3& xv0, TVector3& op1, TVector3& op2 )
{
    double PI = TMath::Pi();
    float x[3],p1[3],p2[3], d;
    StTrackHelix protonobject, pionobject;
    StTrackHelix* proton = &protonobject;
    StTrackHelix* pion = &pionobject;

    //fill StTrackHelix
    //StPhysicalHelixD helix = track1->helix(magnet); //inner helix. good for dca to PV.
    StPicoPhysicalHelix helix = track1->helix(magnet); //inner helix. good for dca to PV.
    int hrot = helix.h();
    float ttheta = helix.phase();
    TVector3 origin = helix.origin();
    float Xc = helix.xcenter();
    float Yc = helix.ycenter();
    float Zc = origin.z();
    float r  = 1./helix.curvature();
    float Vz = r*tan(helix.dipAngle());
    float pt = (helix.momentum(magnet*kilogauss)).Perp();
    float Pz = (helix.momentum(magnet*kilogauss)).Z();

    proton -> Xc = Xc;
    proton -> Yc = Yc;
    proton -> Zc = Zc;
    proton -> r  = r;
    proton -> theta = ttheta;
    proton -> h  = hrot;
    proton -> Vz = Vz;
    proton -> pt = pt;
    proton -> Pz = Pz;
    proton -> Flag = 3;

    helix = track2->helix(magnet); //inner helix. good for dca to PV.
    hrot = helix.h();
    ttheta = helix.phase();
    origin = helix.origin();
    Xc = helix.xcenter();
    Yc = helix.ycenter();
    Zc = origin.z();
    r  = 1./helix.curvature();
    Vz = r*tan(helix.dipAngle());
    pt = (helix.momentum(magnet*kilogauss)).Perp();
    Pz = (helix.momentum(magnet*kilogauss)).Z();

    pion -> Xc = Xc;
    pion -> Yc = Yc;
    pion -> Zc = Zc;
    pion -> r  = r;
    pion -> theta = ttheta;
    pion -> h  = hrot;
    pion -> Vz = Vz;
    pion -> pt = pt;
    pion -> Pz = Pz;
    pion -> Flag = 3;

    //copy Hui Long's 
    float d_root1,d_root2,d_p_pi,theta,r_p,r_pi,zij,rec_zij,rec_zij2;
    float fi_root1_p,fi_root2_p,fi_root1_pi,fi_root2_pi;
    float reci,reci2,recj,recj2;
    float v0_root1[3],v0_root2[3];
    float fi_p,fi_pi,ti,tj;
    int n1,n2,m1,m2,i,j;
    float alfa=0;

    float R2_TPC = 30000.; //hard-coded cut

    d_root2=500.;
    d_root1=500.;

    d_p_pi=sqrt((pion->Xc-proton->Xc)*(pion->Xc-proton->Xc)+(pion->Yc-proton->Yc)*(pion->Yc-proton->Yc));

    theta=atan2((pion->Yc-proton->Yc),(pion->Xc-proton->Xc));

    r_p=proton->r;   //radium of proton helix
    r_pi=pion->r;    //radium of pion helix
    if(d_p_pi>=(r_p+r_pi)||d_p_pi<fabs(r_p-r_pi))return 500.;
    if(d_p_pi>=(r_p+r_pi)){           //XZHU: will not come here, it will be rare that two tracks with d_p_pi == r_p+r_pi come from a single particle decay. that means their momenta should be in the same direction or in opposite direction.
	fi_root1_p=0.;
	fi_root2_p=0.;
	fi_root1_pi=PI;
	fi_root2_pi=PI;
    }
    else {
	alfa=acos((+r_pi*r_pi-d_p_pi*d_p_pi+r_p*r_p)/2./r_p/r_pi);

	fi_root1_p=acos((-r_pi*r_pi+d_p_pi*d_p_pi+r_p*r_p)/2./r_p/d_p_pi);
	fi_root2_p=-fi_root1_p;
	fi_root1_pi=alfa+fi_root1_p;
	fi_root2_pi=- fi_root1_pi;

    }

    v0_root1[2]=-1000.;
    v0_root2[2]=1000.;

    fi_root1_p=fi_root1_p+theta;
    fi_root2_p=fi_root2_p+theta;
    fi_root1_pi=fi_root1_pi+theta;
    fi_root2_pi=fi_root2_pi+theta;

    v0_root1[0]=proton->Xc+proton->r*cos(fi_root1_p);
    v0_root2[0]=proton->Xc+proton->r*cos(fi_root2_p);
    v0_root1[1]=proton->Yc+proton->r*sin(fi_root1_p);
    v0_root2[1]=proton->Yc+proton->r*sin(fi_root2_p);
    fi_root1_p=(fi_root1_p-proton->theta)/proton->h;
    fi_root2_p=(fi_root2_p-proton->theta)/proton->h;
    fi_root1_pi=(fi_root1_pi-pion->theta)/pion->h;
    fi_root2_pi=(fi_root2_pi-pion->theta)/pion->h;
    fi_root1_p=fi_root1_p-(int)(fi_root1_p/2/PI)*2*PI; //XZHU: try floor() here, instead of (int). since a wide range of period is tried below, the result should be the same. in fact, fi_root1_p should always be less than 0 (0,-2pi for 1 period), since we are looking backward for the position where the helix originates.
    fi_root2_p=fi_root2_p-(int)(fi_root2_p/2/PI)*2*PI;
    fi_root1_pi=fi_root1_pi-(int)(fi_root1_pi/2/PI)*2*PI;
    fi_root2_pi=fi_root2_pi-(int)(fi_root2_pi/2/PI)*2*PI;

    if((v0_root1[0]*v0_root1[0]+v0_root1[1]*v0_root1[1])<R2_TPC){ //XZHU: if the overlapping position is close (25 cm) to the outer layer of TPC. Even if the dca between two tracks is small, it is not interesting to us anymore. We find the right periods here, which is not necessary for particles with pt > 0.15GeV/c. in recHyperon, the period finding is not used.
	rec_zij=1000.;
	n1=-((int)(fabs((proton->Zc-pv.z())/proton->Vz/2./PI))+1);
	n2=proton->Flag;  //XZHU: it is meanlingless to go forward in a helix to get the dca point where the particle (or helix) should orginates from! so i suggest to decrease this flag to 0, since fi_root1_pi is [-2pi,2pi]
	m1=-((int)(fabs((pion->Zc-pv.z())/pion->Vz/2./PI))+1);;
	m2=pion->Flag;
	if(n1<-10)n1=-10;
	if(m1<-10)m1=-10;

	reci=fi_root1_p;
	recj=fi_root1_pi;
	for(i=n1;i<=n2;i++) {
	    for(j=m1;j<=m2;j++) {

		ti=fi_root1_p+i*2*PI;
		tj=fi_root1_pi+j*2*PI;
		zij=-proton->Zc-proton->Vz*ti+pion->Zc+pion->Vz*tj;
		if(fabs(zij)<fabs(rec_zij)){
		    rec_zij=zij;
		    reci=ti;
		    recj=tj;
		}
	    }
	}

	fi_root1_p=reci;
	fi_root1_pi=recj;
	d_root1=rec_zij;

	dcaPToPi(&fi_root1_p,&fi_root1_pi,proton,pion,&d_root1,&v0_root1[0],alfa);

    }

    if((v0_root2[0]*v0_root2[0]+v0_root2[1]*v0_root2[1])<R2_TPC){
	rec_zij2=1000.;
	n1=-((int)(fabs((proton->Zc-pv.z())/proton->Vz/2./PI))+1);
	n2=proton->Flag;
	m1=-((int)(fabs((pion->Zc-pv.z())/pion->Vz/2./PI))+1);
	m2=pion->Flag;
	if(n1<-10)n1=-10;
	if(m1<-10)m1=-10;

	reci2=fi_root2_p;
	recj2=fi_root2_pi;
	for(i=n1;i<=n2;i++){
	    for(j=m1;j<=m2;j++){
		ti=fi_root2_p+i*2*PI;
		tj=fi_root2_pi+j*2*PI;
		zij=-proton->Zc-proton->Vz*ti+pion->Zc+pion->Vz*tj;
		if(fabs(zij)<fabs(rec_zij2)){
		    rec_zij2=zij;
		    reci2=ti;
		    recj2=tj;
		}
	    }
	}
	fi_root2_p=reci2;
	fi_root2_pi=recj2;
	d_root2=rec_zij2;

	dcaPToPi(&fi_root2_p,&fi_root2_pi,proton,pion,&d_root2,&v0_root2[0],alfa);
    }


    //if(fabs(d_root2)>d_tolerance&&fabs(d_root1)>d_tolerance)return 1; //XZHU: apply the dca cut here
    if(fabs(d_root1)<fabs(d_root2)){
	fi_p=fi_root1_p;
	fi_pi=fi_root1_pi;
	*x=v0_root1[0];
	*(x+1)=v0_root1[1];
	*(x+2)=v0_root1[2];
	d=d_root1;
    }
    else{
	fi_p=fi_root2_p;
	fi_pi=fi_root2_pi;
	*x=v0_root2[0];
	*(x+1)=v0_root2[1];
	*(x+2)=v0_root2[2];
	d=d_root2;
    }


    *p1=proton->pt*cos(fi_p+proton->h*PI/2.);
    *p2=pion->pt*cos(fi_pi+pion->h*PI/2.);

    *(p1+1)=proton->pt*sin(fi_p+proton->h*PI/2.);
    *(p2+1)=pion->pt*sin(fi_pi+pion->h*PI/2.);
    *(p1+2)=proton->Pz;
    *(p2+2)=pion->Pz;

    //set the StThreeVectors
    xv0.SetX(x[0]);
    xv0.SetY(x[1]);
    xv0.SetZ(x[2]);

    op1.SetX(p1[0]);
    op1.SetY(p1[1]);
    op1.SetZ(p1[2]);

    op2.SetX(p2[0]);
    op2.SetY(p2[1]);
    op2.SetZ(p2[2]);

    return d;
}


//------------------------------------------------------------
double CalcChiFromReso(double res, double (*CalcResoFromChi)(double)) {
	// Calculates chi from the event plane resolution
	double chi   = 2.0;
	double delta = 1.0;

	for (int i = 0; i < 15; i++) {
		chi   = (CalcResoFromChi(chi) < res) ? chi + delta : chi - delta;
		delta = delta / 2.;
	}

	return chi;
}

//------------------------------------------------------------
double CalcResoPsiK2(double chi) {
	// Calculates the event plane resolution as a function of chi
	//  for the case k=2.

	double con = sqrt(M_PI/2)/2; //0.626657;
	double arg = chi*chi/4.;
	double halfpi = M_PI/2; //1.570796;
	double besselOneHalf = sqrt(arg/halfpi) * sinh(arg)/arg;
	double besselThreeHalfs = sqrt(arg/halfpi) * (cosh(arg)/arg - sinh(arg)/(arg*arg));
	double res = con * chi * exp(-arg) * (besselOneHalf + besselThreeHalfs);

	// Approximations.
	//res = 0.25*chi*chi-0.01141*chi*chi*chi-0.034726*chi*chi*chi*chi+0.006815*chi*chi*chi*chi*chi;

	return res;
}


//------------------------------------------------------------
double CalcResoPsiK1(double chi) {
	// Calculates the event plane resolution as a function of chi

	double con = sqrt(M_PI/2)/2; //0.626657;
	double arg = chi*chi/4.;

	double res = con * chi * exp(-arg) * (TMath::BesselI0(arg) + TMath::BesselI1(arg));

	return res;
}


#endif
