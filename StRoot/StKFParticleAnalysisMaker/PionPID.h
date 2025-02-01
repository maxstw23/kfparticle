#ifndef PionPID_hh
#define PionPID_hh

class PionPID
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
    static const float nSigma_mean_piplus[10];
    static const float nSigma_std_piplus[10];
    static const float nSigma_mean_piminus[10];
    static const float nSigma_std_piminus[10];

    static const float nSigma_mean_piplus_p[12];
    static const float nSigma_std_piplus_p[12];
    static const float nSigma_mean_piminus_p[12];
    static const float nSigma_std_piminus_p[12];

    bool IsPionSimpleUsingPt(float nSigmaCut, int charge);
    bool IsPionSimpleUsingP(float nSigmaCut, int charge);
    
public:
    PionPID(float _zTOF, float _nSigma, float _pT, float _p, bool _usingPt): zTOF(_zTOF), nSigma(_nSigma), pT(_pT), p(_p), usingPt(_usingPt) {}
    // bool IsPion();
    bool IsPionSimple(float nSigmaCut, int charge);
};
#endif