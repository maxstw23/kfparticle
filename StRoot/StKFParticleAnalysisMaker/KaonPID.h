#ifndef KaonPID_hh
#define KaonPID_hh

class KaonPID
{
private:
    // PID parameter
    float zTOF;
    float nSigma;
    float pT;
    float p;
    bool usingPt;

    // cut parameters
    static const float num_std;
    static const float nSigma_mean_kplus[9];
    static const float nSigma_std_kplus[9];
    static const float nSigma_mean_kminus[9];
    static const float nSigma_std_kminus[9];

    static const float nSigma_mean_kplus_p[11];
    static const float nSigma_std_kplus_p[11];
    static const float nSigma_mean_kminus_p[11];
    static const float nSigma_std_kminus_p[11];

    bool IsKaonSimpleUsingPt(float nSigmaCut, int charge);
    bool IsKaonSimpleUsingP(float nSigmaCut, int charge);
    
public:
    KaonPID(float _zTOF, float _nSigma, float _pT, float _p, bool usingPt): zTOF(_zTOF), nSigma(_nSigma), pT(_pT), p(_p), usingPt(usingPt) {}
    // bool IsKaon();
    bool IsKaonSimple(float nSigmaCut, int charge);
};
#endif