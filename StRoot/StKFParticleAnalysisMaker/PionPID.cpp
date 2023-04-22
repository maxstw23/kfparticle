#include "PionPID.h"
#include <cmath>

const float PionPID::num_std = 3.5;
const float PionPID::nSigma_mean[] = {0.658390169563609, 0.6838122718403846, 0.6164770413749695, 0.5716575892875071, 0.5327371322112906, 0.5070971873729582, 0.48019457560336803, 0.47204955737480003, 0.4749610090155105};
const float PionPID::nSigma_std[] = {1.1979940379121445, 1.1379816017956788, 1.0598324166186672, 1.0133290967442161, 0.9736662565846723, 0.9518248857673615, 0.9286219803234146, 0.8934160270420024, 0.8386307823608126};
const float PionPID::zTOF_mean[] = {-0.002120831462857123, 3.3437185297351714e-05, -0.0018919150235819188, -0.0023316985562190154, -0.0026139302055204297, -0.002836828336971697, -0.0029467928940488572, -0.003014613235037619, -0.003027950674638674};
const float PionPID::zTOF_std[] = {0.020545365134463627, 0.02064389791536253, 0.011651848134411393, 0.010898485309285019, 0.010480829619920267, 0.010019991056846445, 0.009787114846263463, 0.0092207778642307, 0.008996853589030239};

bool PionPID::IsPion()
// deprecated
{
    int pTbin = static_cast<int>(floor(pT / 0.2));
    if (pTbin < 1 || pTbin > 8) return false;

    // rectangular 2D cut
    if (zTOF > zTOF_mean[pTbin-1] + num_std * zTOF_std[pTbin-1] || zTOF < zTOF_mean[pTbin-1] - num_std * zTOF_std[pTbin-1]) return false;
    if (nSigma > nSigma_mean[pTbin-1] + num_std * nSigma_std[pTbin-1] || nSigma < nSigma_mean[pTbin-1] - num_std * nSigma_std[pTbin-1]) return false;

    // decision boundary cut
    float x = nSigma;
    if (pTbin == 2)
    {
        if (zTOF < -0.00740428142528834*x - 0.309396675573395*sqrt(4.89619390913469e-5*x*x - 0.0608378198738696*x + 1) + 0.237735090905258) return false;
    }
    if (pTbin == 3)
    {
        if (zTOF < -0.0281648248295252*x - 0.962886729723417*sqrt(0.000443038401028118*x*x - 0.0733453661219122*x + 1) + 0.900787586767528) return false;
    } 
    if (pTbin == 4)
    {
        if (zTOF > -0.0032406934598125*x + 0.13698443471894*sqrt(0.00280927391482182*x*x + 0.00221847237193595*x + 1) - 0.0716263807099521) return false;
        if (zTOF < -0.00336998016087225*x - 0.068836431881989*sqrt(-0.00207507565876224*x*x - 0.209026087342191*x + 1) + 0.036153467837767) return false;
    }
    if (pTbin == 5)
    {
        if (zTOF > -0.00213967993133456*x + 0.0676440825750894*sqrt(0.00250094845330029*x*x - 0.00751326905926489*x + 1) - 0.0243909416305087) return false;
        if (zTOF < -0.0038428148130518*x - 0.0747002505703554*sqrt(-0.000587438544893784*x*x - 0.207476923610012*x + 1) + 0.0460501491083069) return false;
    }
    if (pTbin == 6)
    {
        if (zTOF > -0.00175411118868176*x + 0.045714123736958*sqrt(0.000712016627218395*x*x + 0.0104467802787154*x + 1) - 0.0156521506962398) return false;
        if (zTOF < -0.00388033959336782*x - 0.0748278805316802*sqrt(0.000117474696868625*x*x - 0.213332869998102*x + 1) + 0.0499585395905155) return false;
    }
    if (pTbin == 7)
    {
        if (zTOF > -0.00135823151446061*x + 0.0393009793559247*sqrt(-8.1076042377777e-6*x*x + 0.0375948173484308*x + 1) - 0.0158373005570447) return false;
        if (zTOF < -0.00311068719528727*x - 0.0700715902989559*sqrt(-0.000765659313779474*x*x - 0.217735809837499*x + 1) + 0.0482913994238915) return false;
    }
    if (pTbin == 8)
    {
        if (zTOF > -0.000810505839619179*x + 0.035206338530538*sqrt(-0.00195961733178985*x*x + 0.0541320905456119*x + 1) - 0.0164687958929413) return false;
        if (zTOF < -0.00191783901912651*x - 0.0639157056402445*sqrt(-0.00304299323779905*x*x - 0.221235277690126*x + 1) + 0.0443158469801602) return false;
    }

    return true;
}

bool PionPID::IsPionSimple(float nSigmaCut)
{
    int pTbin = static_cast<int>(floor(pT / 0.2));
    if (pTbin < 1 || pTbin > 9) return false;

    // loose nSigma cut
    if (fabs(nSigma-nSigma_mean[pTbin-1])*1.0 / nSigma_std[pTbin-1] > nSigmaCut) return false;

    return true;
}