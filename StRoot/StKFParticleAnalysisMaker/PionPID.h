#ifndef PionPID_hh
#define PionPID_hh

class PionPID
{
private:
    // PID parameter
    float zTOF;
    float nSigma;
    float pT;

    // cut parameters
    static const float num_std;
    static const float nSigma_mean[8];
    static const float nSigma_std[8];
    static const float zTOF_mean[8];
    static const float zTOF_std[8];
    
public:
    PionPID(float _zTOF, float _nSigma, float _pT): zTOF(_zTOF), nSigma(_nSigma), pT(_pT) {}
    bool IsPion();
    bool IsPionSimple(float nSigmaCut);
};
#endif