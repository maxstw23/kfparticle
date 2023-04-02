#include "PionPID.h"
#include <cmath>

const float PionPID::num_std = 3.5;
const float PionPID::nSigma_mean[] = {0.658390169563609, 0.6838122718403846, 0.6164769740917874, 0.5716576204508869, 0.532737224968403, 0.5070869155332025, 0.4801548203441179, 0.47136223233986124};
const float PionPID::nSigma_std[] = {1.1979940379121445, 1.1379816017956788, 1.0598323961061005, 1.0133290940255493, 0.9736662352405308, 0.9518275871331641, 0.9286337304937482, 0.8936356125393272};
const float PionPID::zTOF_mean[] = {-0.002120831462857123, 3.3437185297351714e-05, -0.0018919059822103306, -0.0023316996267758396, -0.002613932781838275, -0.002836614066260732, -0.0029460687860616342, -0.003007975815123568};
const float PionPID::zTOF_std[] = {0.020545365134463627, 0.02064389791536253, 0.011651862698708196, 0.010898483329609066, 0.010480827167512245, 0.010020151111511631, 0.00978757780507371, 0.009223919424297105};

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
    if (pTbin < 1 || pTbin > 8) return false;

    // loose nSigma cut
    if (fabs(nSigma-nSigma_mean[pTbin-1]) > nSigmaCut) return false;

    return true;
}