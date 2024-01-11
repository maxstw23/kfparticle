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
    static const float nSigma_mean_kplus[9];
    static const float nSigma_std_kplus[9];
    // static const float zTOF_mean_kplus[9];
    // static const float zTOF_std_kplus[9];
    static const float nSigma_mean_kminus[9];
    static const float nSigma_std_kminus[9];
    // static const float zTOF_mean_kminus[9];
    // static const float zTOF_std_kminus[9];
    
public:
    KaonPID(float _zTOF, float _nSigma, float _pT): zTOF(_zTOF), nSigma(_nSigma), pT(_pT) {}
    // bool IsKaon();
    bool IsKaonSimple(float nSigmaCut, int charge);
};
#endif