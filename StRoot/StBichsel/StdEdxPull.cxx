#include "Bichsel.h"
#include "StdEdxModel.h"
#include "StdEdxPull.h"
//________________________________________________________________________________
Double_t StdEdxPull::EvalPred(Double_t betagamma, UChar_t fit, Int_t charge) {
  Double_t dedx_expected;
  Double_t dx2 = 1;
  if (TMath::Abs(charge) > 1) dx2 = TMath::Log2(5.);
  if (! fit) { // I70
    dedx_expected = 1.e-6*charge*charge*Bichsel::Instance()->GetI70M(TMath::Log10(betagamma),dx2); 
  } else if ( fit == 1) {     // Ifit
    dedx_expected = 1.e-6*charge*charge*TMath::Exp(Bichsel::Instance()->GetMostProbableZ(TMath::Log10(betagamma),dx2));
  } else {     // dNdx
    dedx_expected = StdEdxModel::instance()->dNdx(betagamma,charge);
  }
  return dedx_expected;
}
//________________________________________________________________________________
Double_t StdEdxPull::EvalDeV(Double_t dEdx, Double_t betagamma, UChar_t fit, Int_t charge) {
  return TMath::Log(dEdx/EvalPred(betagamma, fit, charge));
}
//________________________________________________________________________________
Double_t StdEdxPull::Eval(Double_t dEdx, Double_t dEdxError, Double_t betagamma, UChar_t fit, Int_t charge) {
  return (dEdxError > 0) ? EvalDeV(dEdx, betagamma, fit, charge)/dEdxError : -999;
}
//________________________________________________________________________________
Double_t StdEdxPull::EvalPred2(Double_t betagamma, Double_t dx2, UChar_t fit, Int_t charge) {
  Double_t dedx_expected;
  if (TMath::Abs(charge) > 1) dx2 = TMath::Log2(5.);
  if (! fit) { // I70
    dedx_expected = 1.e-6*charge*charge*Bichsel::Instance()->GetI70M(TMath::Log10(betagamma),dx2); 
  } else if ( fit == 1) {     // Ifit
    dedx_expected = 1.e-6*charge*charge*TMath::Exp(Bichsel::Instance()->GetMostProbableZ(TMath::Log10(betagamma),dx2));
  } else {     // dNdx
    dedx_expected = StdEdxModel::instance()->dNdx(betagamma,charge);
  }
  return dedx_expected;
}
//________________________________________________________________________________
Double_t StdEdxPull::EvalDeV2(Double_t dEdx, Double_t betagamma, Double_t dx2, UChar_t fit, Int_t charge) {
  return TMath::Log(dEdx/EvalPred2(betagamma, dx2, fit, charge));
}
//________________________________________________________________________________
Double_t StdEdxPull::Eval2(Double_t dEdx, Double_t dEdxError, Double_t betagamma, Double_t dx2, UChar_t fit, Int_t charge) {
  return (dEdxError > 0) ? EvalDeV2(dEdx, betagamma, dx2, fit, charge)/dEdxError : -999;
}
