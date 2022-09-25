#ifndef KaonPID_hh
#define KaonPID_hh

class KaonPID
{
private:
    // PID parameter
    float zTOF;
    float nSigma;
    float pT;

    // cut parameters
    static const float num_std = 3.5;
    static const float nSigma_mean[11];
    static const float nSigma_std[11];
    static const float zTOF_mean[11];
    static const float zTOF_std[11];
    
public:
    KaonPID(float _zTOF, float _nSigma, float _pT): zTOF(_zTOF), nSigma(_nSigma), pT(_pT) {}
    bool IsKaon();
}
#endif