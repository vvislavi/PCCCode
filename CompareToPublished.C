#include "Scripts/CommonFunctions.C"
#include "Scripts/BuildCanvas.C"
#include "AliGFWFlowContainer.h"
#include "Scripts/ProcessAxii.C"

TFile *tfP=0;
AliGFWFlowContainer *fc[3]={0,0,0};
TString fiNas[] = {"LHC15o.root","LHC15o_pass1.root","LHC18qr.root"};//"LHC18qr.root","Merged.root"};
TString dsNam[] = {"LHC15o", "LHC15o_pass1","LHC18qr"};//, "LHC15o18qr"};
Int_t cols[] = {kRed+1,kBlue+1,kBlack};
Int_t linStyle[] = {2,2,1};

void ShapeHist(TH1 *inh) {
  ProcessAxisPx(inh->GetYaxis(),25,0.01,30,1.75,507,0.03);
  ProcessAxisPx(inh->GetXaxis(),25,0,30,3,507,0.03);
}
void Init() {
  if(!tfP) tfP = getFile("Published/HEPData-ins1666817-v1-root.root");
  for(Int_t i=0;i<3;i++) {
    if(fc[i]) continue;
    auto tf = getFile(fiNas[i]);
    fc[i] = (AliGFWFlowContainer*)getObj(tf,"FlowContainer_0");
    tf->Close();
  };
}
void SetIDName(TString newID) {
  for(Int_t i=0;i<3;i++) fc[i]->SetIDName(newID);
}
TH1 **BuildHists(Int_t nHar, Int_t nCor=2, TString idName="ChFull") {
  SetIDName(idName);
  TH1 **rh = new TH1*[3];
  for(Int_t i=0;i<3;i++) {
    rh[i] = (nCor==2)?fc[i]->GetVN2VsMulti(nHar):fc[i]->GetVN4VsMulti(nHar);
    rh[i]->SetLineColor(cols[i]);
    rh[i]->SetLineWidth(2);
    rh[i]->SetLineStyle(linStyle[i]);
    ShapeHist(rh[i]);
  };
  return rh;
}
TH1 **makeRatios(TH1 **num, TGraph *den) {
  TH1 **rh = new TH1*[3];
  Double_t *xv = den->GetX();
  Double_t *yv = den->GetY();
  for(Int_t i=0;i<3;i++) {
    rh[i] = (TH1*)num[i]->Clone(Form("%s_ratio",num[i]->GetName()));
    rh[i]->Reset();
    for(Int_t j=0;j<den->GetN();j++) {
      Int_t bi = rh[i]->FindBin(xv[j]);
      if(bi<1 || bi>rh[i]->GetNbinsX()) continue;
      rh[i]->SetBinContent(bi,num[i]->GetBinContent(bi)/yv[j]);
    }
  }
  return rh;
}
TCanvas *CompareToPublished(Int_t tNo, Int_t nHar, Int_t nCor, TString yLabel="v_{2}{2}", TString idName="ChGap") {
  Init();
  TCanvas *c = BuildCanvasForRatios();
  TPad *top = (TPad*)c->FindObject("top");
  TPad *bot = (TPad*)c->FindObject("bot");
  auto tdf = getObj(tfP,Form("Table %i",tNo));
  auto gre = (TGraphAsymmErrors*)getObj(tdf,"Graph1D_y1");
  gre->SetTitle("");
  gre->SetLineColor(kGray+1);
  top->cd();
  TH1 **h = BuildHists(nHar,nCor,idName);
  h[0]->GetYaxis()->SetTitle(yLabel.Data());
  h[0]->GetYaxis()->SetRangeUser(0.0001,0.12);
  for(Int_t i=0;i<3;i++) h[i]->Draw(i?"SAME":"");
  gre->Draw("SAME P");
  TLegend *tleg = Legend(0.53457,0.0273881,0.834739,0.326497);
  tleg->SetTextSize(25);
  tleg->AddEntry(gre,"JHEP09 (2018) 06","L");
  for(Int_t i=0;i<3;i++) tleg->AddEntry(h[i],dsNam[i].Data(),"L");
  tleg->Draw();
  bot->cd();
  TH1 **hr = makeRatios(h,gre);
  for(Int_t i=0;i<3;i++) {
    hr[i]->GetYaxis()->SetRangeUser(0.81,1.19);
    hr[i]->SetTitle("; V0M %; Ratio to pub.");
    hr[i]->Draw(i?"SAME":"");
  };
  return c;
}
void CompareToPublished() {
  TCanvas *c = CompareToPublished(1,2,2,"v_{2}{2, |#Delta#eta|>X}","ChGap");
  c->Print("Plots20210816/vnComp.pdf(");
  delete c;
  c = CompareToPublished(193,2,2,"v_{2}{2}","ChFull");
  c->Print("Plots20210816/vnComp.pdf");
  delete c;
  c = CompareToPublished(194,2,4,"v_{2}{4}","ChFull");
  c->Print("Plots20210816/vnComp.pdf");
  delete c;
  c = CompareToPublished(3,3,2,"v_{3}{2, |#Delta#eta|>X}","ChGap");
  c->Print("Plots20210816/vnComp.pdf");
  delete c;
  c = CompareToPublished(4,4,2,"v_{4}{2, |#Delta#eta|>X}","ChGap");
  c->Print("Plots20210816/vnComp.pdf)");
  delete c;

}
