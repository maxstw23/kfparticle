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
    static const float num_std;
    static const float nSigma_mean[13];
    static const float nSigma_std[13];
    static const float zTOF_mean[13];
    static const float zTOF_std[13];
    
public:
    KaonPID(float _zTOF, float _nSigma, float _pT): zTOF(_zTOF), nSigma(_nSigma), pT(_pT) {}
    bool IsKaon();
    bool IsKaonSimple(float nSigmaCut);
};
#endif