#ifndef ProtonPID_hh
#define ProtonPID_hh

class ProtonPID
{
private:
    // PID parameter
    float zTOF;
    float nSigma;
    float pT;
    float p; 
    bool usingPt;

    // cut parameters using pT
    static const float num_std;
    static const float nSigma_mean_proton[9];
    static const float nSigma_std_proton[9];
    static const float nSigma_mean_antiproton[9];
    static const float nSigma_std_antiproton[9];

    // cut parameters using p
    static const float nSigma_mean_proton_p[11];
    static const float nSigma_std_proton_p[11];
    static const float nSigma_mean_antiproton_p[11];
    static const float nSigma_std_antiproton_p[11];

    bool IsProtonSimpleUsingPt(float nSigmaCut, int charge);
    bool IsProtonSimpleUsingP(float nSigmaCut, int charge);
    
public:
    ProtonPID(float _zTOF, float _nSigma, float _pT, float _p, bool _usingPt): zTOF(_zTOF), nSigma(_nSigma), pT(_pT), p(_p), usingPt(_usingPt) {}
    // bool IsProton();
    bool IsProtonSimple(float nSigmaCut, int charge);
};
#endif