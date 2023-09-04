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
    static const float nSigma_mean_piplus[9];
    static const float nSigma_std_piplus[9];
    static const float zTOF_mean_piplus[9];
    static const float zTOF_std_piplus[9];
    static const float nSigma_mean_piminus[9];
    static const float nSigma_std_piminus[9];
    static const float zTOF_mean_piminus[9];
    static const float zTOF_std_piminus[9];
    
public:
    PionPID(float _zTOF, float _nSigma, float _pT): zTOF(_zTOF), nSigma(_nSigma), pT(_pT) {}
    // bool IsPion();
    bool IsPionSimple(float nSigmaCut, int charge);
};
#endif