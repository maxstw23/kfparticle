#include "ProtonPID.h"
#include <cmath>

const float ProtonPID::num_std = 3.5;
const float ProtonPID::nSigma_mean_proton[] = {0.32856898912065524, 0.4098843660632264, 0.4600929184220844, 0.4266253776559791, 0.45284868030347997, 0.5095599225084195, 0.5529033686842685, 0.5834227289551327, 0.6024278635830862};
const float ProtonPID::nSigma_std_proton[] = {1.0325488119903803, 1.019138152033574, 1.0363375865903763, 1.042202218162307, 1.0544129337920174, 1.0453919754206316, 1.0367009805622085, 1.0295642209306577, 1.02436372771767};
const float ProtonPID::nSigma_mean_antiproton[] = {0.23958830290854294, 0.36272702881341734, 0.40074475979717916, 0.3652481279130016, 0.3919338473398265, 0.45447937632791846, 0.5048410121895403, 0.5416066733748667, 0.5669124910024287};
const float ProtonPID::nSigma_std_antiproton[] = {1.0440974362136834, 1.0192500504204944, 1.0371068367543788, 1.0418732888278077, 1.055869165454058, 1.0465073389441732, 1.0380730977497155, 1.0314764267490484, 1.0277298808381354};

const float ProtonPID::nSigma_mean_proton_p[] = {-0.1550634774780441, 0.24924573955397406, 0.3793280268387734, 0.4684237826084941, 0.43468021245284333, 0.5262249707145028, 0.591705622429231, 0.6421931481973874, 0.6761232766031283, 0.7002805282104168, 0.7172769101663811};
const float ProtonPID::nSigma_std_proton_p[] = {0.9643696817681089, 0.9971771227158162, 1.023114280268588, 1.0428492811822598, 1.0433730036271651, 1.0384517560166173, 1.0325027128278033, 1.0267308593983375, 1.0216732817625107, 1.0175972659489643, 1.0142656390934357};
const float ProtonPID::nSigma_mean_antiproton_p[] = {-0.1795547981328962, 0.19226420926466237, 0.31487666997733277, 0.40750921339169044, 0.37381958740327537, 0.47009725780952677, 0.54019083534965, 0.5964253956693281, 0.6357168031818341, 0.6649653130519547, 0.6862084750978683};
const float ProtonPID::nSigma_std_antiproton_p[] = {0.9924728608957202, 0.9965379652685963, 1.022582641393169, 1.0427466174430926, 1.044210729492363, 1.0401612570625598, 1.0339334732185799, 1.0276356555162363, 1.0237829631958855, 1.0201133204862405, 1.0178713576486667};
bool ProtonPID::IsProtonSimpleUsingPt(float nSigmaCut, int charge)
{
    int pTbin = static_cast<int>(floor(pT / 0.2));
    if (pTbin < 1 || pTbin > 9) return false;
    if (std::abs(charge) != 1) return false;

    // loose nSigma cut
    if (charge > 0)
    {
        if (fabs(nSigma-nSigma_mean_proton[pTbin-1])*1.0 / nSigma_std_proton[pTbin-1] > nSigmaCut) return false;
        return true;
    }
    else
    {
        if (fabs(nSigma-nSigma_mean_antiproton[pTbin-1])*1.0 / nSigma_std_antiproton[pTbin-1] > nSigmaCut) return false;
        return true;
    }
}

bool ProtonPID::IsProtonSimpleUsingP(float nSigmaCut, int charge)
{
    int pbin = static_cast<int>(floor(p / 0.2));
    if (pbin < 1 || pbin > 11) return false;
    if (std::abs(charge) != 1) return false;

    // loose nSigma cut
    if (charge > 0)
    {
        if (fabs(nSigma-nSigma_mean_proton_p[pbin-1])*1.0 / nSigma_std_proton_p[pbin-1] > nSigmaCut) return false;
        return true;
    }
    else
    {
        if (fabs(nSigma-nSigma_mean_antiproton_p[pbin-1])*1.0 / nSigma_std_antiproton_p[pbin-1] > nSigmaCut) return false;
        return true;
    }
}

bool ProtonPID::IsProtonSimple(float nSigmaCut, int charge)
{
    if (usingPt) return IsProtonSimpleUsingPt(nSigmaCut, charge);
    else return IsProtonSimpleUsingP(nSigmaCut, charge);
}