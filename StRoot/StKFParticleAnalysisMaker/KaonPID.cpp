#include "KaonPID.h"
#include <cmath>

const float KaonPID::num_std = 3.5;
const float KaonPID::nSigma_mean_kplus[] = {0.5002158135362356, 0.4182163349133706, 0.48130210103327153, 0.5531276662139721, 0.6000251832114242, 0.6361808140088189, 0.6664394949072706, 0.6955447308554751, 0.7319082876913042};
const float KaonPID::nSigma_std_kplus[] = {1.0074862078322409, 1.0398423167015347, 1.046268368296106, 1.0356657399676232, 1.0314912742131301, 1.0365436773428698, 1.0532232444704066, 1.0818484732311449, 1.1239060512485277};
const float KaonPID::nSigma_mean_kminus[] = {0.4603743421422999, 0.39048514414138397, 0.45675300060471397, 0.5352877387650564, 0.5909949095530086, 0.6405032626406305, 0.6909996433334052, 0.7520789769114872, 0.8358994928316273};
const float KaonPID::nSigma_std_kminus[] = {1.0082989744309434, 1.0411282202610077, 1.0474884913921618, 1.0367352596171975, 1.033641835972795, 1.0449726006752518, 1.0759432707172976, 1.1280692090651037, 1.203452346603288};

const float KaonPID::nSigma_mean_kplus_p[] = {0.3411630321456646, 0.37614693269549027, 0.4527071749870946, 0.5517597129577906, 0.6092439348954245, 0.6531890747394016, 0.6892532863649086, 0.7214791027032803, 0.7529438251293079, 0.7907644936660417, 0.834242912581453};
const float KaonPID::nSigma_std_kplus_p[] = {0.9862859189636303, 1.0349729119274007, 1.0421608419111903, 1.0352786846066142, 1.0267692482316846, 1.0251578745514163, 1.0326039213945157, 1.046253082448786, 1.0668808372566327, 1.0952903649441064, 1.131671442503555};
const float KaonPID::nSigma_mean_kminus_p[] = {0.2990107698406812, 0.345205430247423, 0.4264781879876914, 0.5300847786306762, 0.5942839000297531, 0.6470026496240403, 0.6960391395497489, 0.7455053913977175, 0.7996793611012166, 0.8701304397529926, 0.9617472493460847};
const float KaonPID::nSigma_std_kminus_p[] = {0.9859230982444928, 1.0355911815032441, 1.0430912592514512, 1.035963938421809, 1.0264824786765663, 1.025615867474706, 1.0383492697231576, 1.0630332033554029, 1.10061946788938, 1.1540447493766162, 1.2214916062171755};

bool KaonPID::IsKaonSimpleUsingPt(float nSigmaCut, int charge)
{
    int pTbin = static_cast<int>(floor(pT / 0.2));
    if (pTbin < 1 || pTbin > 9) return false;
    if (std::abs(charge) != 1) return false;

    // loose nSigma cut
    if (charge > 0)
    {
        if (fabs(nSigma - nSigma_mean_kplus[pTbin-1])*1.0 / nSigma_std_kplus[pTbin-1] > nSigmaCut) return false;
        return true;
    }
    else
    {
        if (fabs(nSigma - nSigma_mean_kminus[pTbin-1])*1.0 / nSigma_std_kminus[pTbin-1]> nSigmaCut) return false;
        return true;
    }
}

bool KaonPID::IsKaonSimpleUsingP(float nSigmaCut, int charge)
{
    int pbin = static_cast<int>(floor(p / 0.2));
    if (pbin < 1 || pbin > 11) return false;
    if (std::abs(charge) != 1) return false;

    // loose nSigma cut
    if (charge > 0)
    {
        if (fabs(nSigma - nSigma_mean_kplus_p[pbin-1])*1.0 / nSigma_std_kplus_p[pbin-1] > nSigmaCut) return false;
        return true;
    }
    else
    {
        if (fabs(nSigma - nSigma_mean_kminus_p[pbin-1])*1.0 / nSigma_std_kminus_p[pbin-1]> nSigmaCut) return false;
        return true;
    }
}

bool KaonPID::IsKaonSimple(float nSigmaCut, int charge)
{
    if (usingPt) return IsKaonSimpleUsingPt(nSigmaCut, charge);
    else return IsKaonSimpleUsingP(nSigmaCut, charge);
}