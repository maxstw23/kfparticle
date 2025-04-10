class Bichsel;
//Bichsel *m_Bichsel = Bichsel::instance();
Bichsel *m_Bichsel = 0;
TF1 *sipi = 0, *sie = 0, *sip = 0, *siK = 0, *sid = 0, *sit = 0, *sihe3 = 0;
//________________________________________________________________________________
Double_t sifunc(Double_t *x,Double_t *par) {
  Double_t pove   = pow(10.,x[0]);
  Double_t ppion  = pove*0.13956995;
  Double_t poverm = ppion/par[0];
  Int_t    k      = par[1];
  Double_t charge2     = 1;
  if (k == 1) {
    charge2 = 4;
    poverm *= 2;
  }
  return TMath::Log(charge2*BetheBloch::Sirrf(poverm,60,k==3)/BetheBloch::Sirrf(pove,60,0));
}
//________________________________________________________________________________
Double_t bichselZ(Double_t *x,Double_t *par) {
  Double_t ppion  = pow(10.,x[0]);
  Double_t pove   = ppion;///0.13956995;
  Double_t poverm = ppion/par[0];
  Int_t    k      = par[1];
  Double_t charge2     = 1;
  if (k == 1) {
    charge2 = 4;
    poverm *= 2;
  }
  Double_t val1 = m_Bichsel->GetMostProbableZ(TMath::Log10(poverm),1.);
  Double_t val2 = m_Bichsel->GetMostProbableZ(TMath::Log10(pove),1.);
  Double_t val  = val1 - val2;
  //  cout << "pove\t" << pove << "\t" << val1 << "\tpoverm\t" << val2 << "\t" << val << endl;
  return val;
}
//________________________________________________________________________________

void Bichsel(const Char_t *tag = "") {
#if 1
  if (gClassTable->GetID("StBichsel")) {
    gSystem->Load("libTable");
    gSystem->Load("St_base");
    gSystem->Load("StarClassLibrary");
    gSystem->Load("StBichsel");
  }
#endif
  if (!m_Bichsel) m_Bichsel = Bichsel::Instance();
  TLegend *leg = new TLegend(0.72,0.7,0.9,0.9,"");//TLegend(0.79,0.91,0.89,0.89,"");
  sipi = new TF1(Form("%sbipi",tag),bichselZ,-1.,4.0,2);
  sipi->SetParameter(0,0.13956995);
  sipi->SetParameter(1,0);
  sie = new TF1(Form("%sbie",tag),bichselZ,-1.,4.0,2);
  sie->SetParameter(0,0.51099907e-3);
  sie->SetParameter(1,3);
  sie->SetLineColor(2);
  leg->AddEntry(sie,"e-#pi","L");
//   TF1 *sie1 = new TF1(Form("%sbie1",tag),bichselZ,-1.,4.0,2);
//   sie1->SetParameter(0,0.51099907e-3);
//   sie1->SetParameter(1,0);
//   sie1->SetLineColor(2);
  sip = new TF1(Form("%sbi",tag),bichselZ,-1.,4.0,2);
  sip->SetParameter(0,0.93827231);
  sip->SetParameter(1,0);
  sip->SetLineColor(3);
  leg->AddEntry(sip,"p-#pi","L");
  sim = new TF1(Form("%sbim",tag),bichselZ,-1.,4.0,2);
  sim->SetParameter(0,0.1056584);
  sim->SetParameter(1,0);
  sim->SetLineColor(8);
  leg->AddEntry(sim,"#mu-#pi","L");
  siK = new TF1(Form("%sbiK",tag),bichselZ,-1.,4.0,2);
  siK->SetParameter(0,0.493677);
  siK->SetParameter(1,0);
  siK->SetLineColor(4);
  leg->AddEntry(siK,"K-#pi","L");
  sid = new TF1(Form("%sbid",tag),bichselZ,-1.,4.0,2);
  sid->SetParameter(0,0.1876E+01);
  sid->SetParameter(1,0);
  sid->SetLineColor(6);
  leg->AddEntry(sid,"d-#pi","L");
  TF1 *sit = new TF1(Form("%sbit",tag),bichselZ,-1.,4.0,2);
  sit->SetParameter(0,0.2809E+01);
  sit->SetParameter(1,0);
  sit->SetLineColor(7);
  leg->AddEntry(sit,"t-#pi","L");

//   sihe3 = new TF1(Form("%sbihe3",tag),bichselZ,-1.,4.0,2);
//   sihe3->SetParameter(0,0.2809E+01);
//   sihe3->SetParameter(1,1);
//   sihe3->SetLineColor(1);
//   leg->AddEntry(sihe3,"He3-#pi","L");
  sipi->Draw("same");
  sim->Draw("same");
  sie->Draw("same");
  //  sie1->Draw("same");
  siK->Draw("same");
  sip->Draw("same");
  sid->Draw("same");
  sit->Draw("same");
//   sihe3->Draw("same");
  leg->Draw();
}
Double_t Divation(Double_t *x,Double_t *par) {
  Double_t val[7];
  //  val[0]  = sipi->Eval(x);
  TF1 *f;
  f = (TF1 *) gROOT->GetListOfFunctions()->FindObject("bie");  if (f) val[1] = f->Eval(x);
  f = (TF1 *) gROOT->GetListOfFunctions()->FindObject("biK");  if (f) val[2] = f->Eval(x);
  f = (TF1 *) gROOT->GetListOfFunctions()->FindObject("bip");  if (f) val[3] = f->Eval(x);
  f = (TF1 *) gROOT->GetListOfFunctions()->FindObject("bid");  if (f) val[4] = f->Eval(x);
  f = (TF1 *) gROOT->GetListOfFunctions()->FindObject("bit");  if (f) val[5] = f->Eval(x);
  f = (TF1 *) gROOT->GetListOfFunctions()->FindObject("bihe3");if (f) val[6] = f->Eval(x);
  Double_t min = 9999.;
  for (int i = 1; i <7; i++) {
    Double_t dev = TMath::Abs(val[i]);
    if (dev < min) min = dev;
  }
  return min;
}
