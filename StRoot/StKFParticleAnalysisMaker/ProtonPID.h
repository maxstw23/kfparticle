#ifndef ProtonPID_hh
#define ProtonPID_hh

class ProtonPID
{
private:
    // PID parameter
    float zTOF;
    float nSigma;
    float pT;

    // cut parameters
    static const float num_std;
    static const float nSigma_mean[9];
    static const float nSigma_std[9];
    static const float zTOF_mean[9];
    static const float zTOF_std[9];
    
public:
    ProtonPID(float _zTOF, float _nSigma, float _pT): zTOF(_zTOF), nSigma(_nSigma), pT(_pT) {}
    bool IsProton();
    bool IsProtonSimple(float nSigmaCut);
};
#endif