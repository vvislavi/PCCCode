#include "Scripts/BuildCanvas.C"
#include "Scripts/PCCContainer.h"
#include "TLegend.h"
#include "TString.h"
#include "Scripts/CommonFunctions.C"
#include "Scripts/ProcessAxii.C"
using namespace PCCSpace;

PCCContainer *pc18=0;
PCCContainer *pc15=0;
Int_t nrb=1;
void ShapeHist(TH1 *inh) {
  ProcessAxisPx(inh->GetYaxis(),25,0.01,30,1.75,507,0.03);
  ProcessAxisPx(inh->GetXaxis(),25,0,30,3,507,0.03);
}
TCanvas *DrawComparison(PCCContainer::lFunc sf, Double_t ymin=0.001, Double_t ymax=0.3, TString exName="") {
  if(!pc15) {
    pc15 = new PCCContainer("LHC15o.root",1);
    pc18 = new PCCContainer("LHC15o_pass1.root",1);//LHC18qr.root",1);
  };
  pc15->RebinMulti(nrb,0);
  pc18->RebinMulti(nrb,0);
  TCanvas *rc = BuildCanvasForRatios(500,0.33,0.02);
  TH1 *h18 = pc18->Var(sf);
  TH1 *h15 = pc15->Var(sf);
  h18->GetXaxis()->SetTitle("V0M %");
  h15->GetXaxis()->SetTitle("V0M %");
  TPad *tp = (TPad*)rc->FindObject("top");
  tp->cd();
  TLegend *tleg = 0;
  if(sf==kPCC3) tleg = Legend(0.193333,0.671346,0.493333,0.971772);
  else tleg = Legend(0.183811,0.0436661,0.48398,0.342775);
  if(!exName.IsNull()) tleg->AddEntry((TObject*)0x0,exName.Data(),"");
  tleg->AddEntry(h15,"LHC15o pass2","L");
  tleg->AddEntry(h18,"LHC15o pass1", "L");//"LHC18qr pass3","L");
  h15->SetLineColor(kRed+1);
  h18->SetLineColor(kBlue+1);
  h15->SetLineWidth(2);
  h18->SetLineWidth(2);
  h18->GetYaxis()->SetRangeUser(ymin,ymax);
  ShapeHist(h18);
  ShapeHist(h15);
  h18->Draw();
  h15->Draw("SAME");
  tleg->Draw();
  TH1 *hrat = (TH1*)h15->Clone(Form("%s_ratio",h15->GetName()));
  hrat->Divide(h18);
  hrat->GetYaxis()->SetRangeUser(0.49,1.51);
  hrat->GetYaxis()->SetTitle("Ratio, pass2/pass1");
  TPad *bt = (TPad*)rc->FindObject("bot");
  bt->cd();
  hrat->Draw();
  return rc;
}
void DrawComparisons() {
  TCanvas *c;
  c = DrawComparison(kPCC2);
  c->Print("Plots20210816/PCC.pdf(");
  delete c;
  c = DrawComparison(kPCC3Sub,0.001,0.3,"3-sub");
  c->Print("Plots20210816/PCC.pdf");
  delete c;
  c = DrawComparison(kPCC3);
  c->Print("Plots20210816/PCC.pdf");
  delete c;
  c = DrawComparison(kMultiHar,-0.99,0.99);
  c->Print("Plots20210816/PCC.pdf)");
  delete c;
}
