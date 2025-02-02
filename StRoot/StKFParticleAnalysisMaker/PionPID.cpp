#include "PionPID.h"
#include <cmath>

const float PionPID::num_std = 3.5;
const float PionPID::nSigma_mean_piplus[] = {0.5172574838601324, 0.5406139994545175, 0.5685486670476058, 0.5344733663771382, 0.49436089694291246, 0.4603414809242065, 0.42720096570839733, 0.3901453906880524, 0.34239177550148964, 0.27799814684995383};
const float PionPID::nSigma_std_piplus[] = {1.0785663045184555, 1.0498564493999676, 1.0214694978206957, 1.0062428420212832, 0.9969307745273558, 0.9914228113672997, 0.988554424539734, 0.9918256410085884, 1.0083336866249117, 1.0432030568700255};
const float PionPID::nSigma_mean_piminus[] = {0.4961689158439063, 0.529441156597391, 0.5664168427060067, 0.5355380540169318, 0.4961384249576719, 0.4625360715264887, 0.43109723945283296, 0.3998787019053495, 0.36377286926802105, 0.319541912172933};
const float PionPID::nSigma_std_piminus[] = {1.0804506814948749, 1.049048291950112, 1.0206250120643778, 1.0055839710166967, 0.996588029744264, 0.990551389241719, 0.9867527998943671, 0.9876729273103109, 0.9968284189373825, 1.0164063359251159};

const float PionPID::nSigma_mean_piplus_p[] = {0.4077357290429908, 0.5148235716230026, 0.5714997885232493, 0.5583536843914758, 0.5393770178733244, 0.5297366239737236, 0.5268395050550745, 0.5250489804750966, 0.5210443982761432, 0.5116120710692182, 0.49420665694346766, 0.4640701026167111};
const float PionPID::nSigma_std_piplus_p[] = {1.085731220246315, 1.0556653435823597, 1.0280733571031964, 1.0133440754325662, 1.0027893488373891, 0.9943464293418831, 0.98843833764244, 0.984092944177352, 0.9831085959482985, 0.9903835837612632, 1.0095964408559732, 1.0444220522648946};
const float PionPID::nSigma_mean_piminus_p[] = {0.37666939207304817, 0.5005848330640008, 0.5685627437194825, 0.5599042341984612, 0.5427193766972309, 0.5336230808888304, 0.5307827968302206, 0.5298701882257401, 0.5304160181185413, 0.5291658955328298, 0.5242097295305626, 0.5143501397838433};
const float PionPID::nSigma_std_piminus_p[] = {1.0889318231391871, 1.0548886161617546, 1.0270834475512338, 1.0124720971949202, 1.0021686430086525, 0.9938939236510239, 0.9870980178508275, 0.982148573528944, 0.9797501405072274, 0.9821963552157307, 0.9904646493553048, 1.0062328015186137};
bool PionPID::IsPionSimpleUsingPt(float nSigmaCut, int charge)
{
    int pTbin = static_cast<int>(floor(pT / 0.2));
    if (pTbin > 9) return false;
    if (std::abs(charge) != 1) return false;

    // loose nSigma cut
    if (charge > 0)
    {
        if (fabs(nSigma-nSigma_mean_piplus[pTbin])*1.0 / nSigma_std_piplus[pTbin] > nSigmaCut) return false;
        return true;
    }
    else
    {
        
        if (fabs(nSigma-nSigma_mean_piminus[pTbin])*1.0 / nSigma_std_piminus[pTbin] > nSigmaCut) return false;
        return true;
    }
}

bool PionPID::IsPionSimpleUsingP(float nSigmaCut, int charge)
{
    int pbin = static_cast<int>(floor(p / 0.2));
    if (pbin > 11) return false;
    if (std::abs(charge) != 1) return false;

    // loose nSigma cut
    if (charge > 0)
    {
        if (fabs(nSigma-nSigma_mean_piplus_p[pbin])*1.0 / nSigma_std_piplus_p[pbin] > nSigmaCut) return false;
        return true;
    }
    else
    {
        if (fabs(nSigma-nSigma_mean_piminus_p[pbin])*1.0 / nSigma_std_piminus_p[pbin] > nSigmaCut) return false;
        return true;
    }
}

bool PionPID::IsPionSimple(float nSigmaCut, int charge)
{
    if (usingPt) return IsPionSimpleUsingPt(nSigmaCut, charge);
    else return IsPionSimpleUsingP(nSigmaCut, charge);
}

