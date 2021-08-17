#ifndef BUILDCANVAS__C
#define BUILDCANVAS__C
#include "TCanvas.h"
#include "TPad.h"
TCanvas *BuildCanvasForRatios(Double_t ph=500, Double_t rr=0.33, Double_t tm=0.01, Double_t bm=0.12, Double_t lm=0.15) {
  Double_t cheight = ph*(1.+rr)/(1.-tm-bm);
  Double_t cwidth = ph/(1.-tm-lm);
  TCanvas *lc = new TCanvas("c","c",cwidth,cheight);
  lc->SetMargin(0,0,0,0);
  Double_t bheightabs = ph*rr+cheight*bm;
  Double_t bheightfrac = bheightabs/cheight;
  Double_t theightfrac = 1-bheightfrac;
  Double_t theightabs = cheight*theightfrac;
  
  TPad *top = new TPad("top","top",0,bheightfrac,1,1);
  top->SetMargin(lm,tm,0,(tm*cheight)/theightabs);
  lc->cd();
  top->Draw();
  lc->cd();
  TPad *bot = new TPad("bot","bot",0,0,1,bheightfrac);
  bot->SetMargin(lm,tm,(bm*cheight)/bheightabs,0);
  bot->Draw();
  lc->cd();
  return lc;
};

TCanvas *BuildCanvasForSpectra(Double_t ph=600, Double_t pw=500, Double_t rr=0.33, Double_t tm=0.02, Double_t bm=0.1, Double_t lm=0.18) {
  Double_t cheight = ph*(1.+rr)/(1.-tm-bm);
  Double_t cwidth = pw/(1.-tm-lm);
  TCanvas *lc = new TCanvas("c","c",cwidth,cheight);
  lc->SetMargin(0,0,0,0);
  Double_t bheightabs = ph*rr+cheight*bm;
  Double_t bheightfrac = bheightabs/cheight;
  Double_t theightfrac = 1-bheightfrac;
  Double_t theightabs = cheight*theightfrac;
  
  TPad *top = new TPad("top","top",0,bheightfrac,1,1);
  top->SetMargin(lm,tm,0,(tm*cheight)/theightabs);
  lc->cd();
  top->Draw();
  lc->cd();
  TPad *bot = new TPad("bot","bot",0,0,1,bheightfrac);
  bot->SetMargin(lm,tm,(bm*cheight)/bheightabs,0);
  bot->Draw();
  lc->cd();
  // TPad *mask = new TPad("mask","mask",0,0,1,1);
  // mask->SetFillStyle(0);
  // mask->Draw();
  return lc;
};
TCanvas *BuildSingleCanvas(Double_t ph=500, Double_t pw=500, Double_t lm=0.15, Double_t bm=0.13, Double_t tm=0.04, Double_t rm=-1) {
  if(rm<0) rm = tm;
  Double_t cwidth = pw/(1-lm-tm);
  Double_t cheight = ph/(1-bm-tm);
  TCanvas *c = new TCanvas("c","c",cwidth,cheight);
  c->SetMargin(lm,rm,bm,tm);
  return c;
};
TCanvas *BuildParallelCanvases(Double_t ph=500, Double_t pw=500, Double_t lm=0.15, Double_t bm=0.13, Double_t tm=0.04) {
  Double_t cwidth = 2*pw/(1-lm-tm);
  Double_t cheight = ph/(1-bm-tm);
  TCanvas *c = new TCanvas("c","c",cwidth,cheight);
  c->SetMargin(0,0,0,0);
  c->cd();
  TPad *tpl = new TPad("tpl","tpl",0,0,0.5,1);
  tpl->SetMargin(lm,tm,bm,tm);
  tpl->Draw();
  c->cd();
  TPad *tpr = new TPad("tpr","tpr",0.5,0,1,1);
  tpr->SetMargin(lm,tm,bm,tm);
  tpr->Draw();
  return c;
};
TCanvas *BuildParallelCanvasesMerged(Double_t ph=500, Double_t pw=500, Double_t lm=0.15, Double_t bm=0.13, Double_t tm=0.04, Double_t rm=0.04) {
  Double_t p1w = pw/(1-lm);
  Double_t p2w = pw/(1-rm);
  Double_t cwidth = p1w+p2w;
  Double_t DivPoint = p1w/cwidth;
  Double_t cheight = ph/(1-bm-tm);
  TCanvas *c = new TCanvas("c","c",cwidth,cheight);
  c->SetMargin(0,0,0,0);
  c->cd();
  TPad *tpl = new TPad("tpl","tpl",0,0,DivPoint,1);
  tpl->SetMargin(lm,0.,bm,tm);
  tpl->Draw();
  c->cd();
  TPad *tpr = new TPad("tpr","tpr",DivPoint,0,1,1);
  tpr->SetMargin(0.,rm,bm,tm);
  tpr->Draw();
  return c;
};

TCanvas *BuildThreeInFour(Double_t ph=500, Double_t pw=500,Double_t lm=0.2,Double_t bm=0.13,Double_t tm=0.04,Double_t rm=0.04) {
  //Calculating dimensions
  Double_t lpw = pw/(1-lm);
  Double_t rpw = pw/(1-rm);
  Double_t tlph = ph/(1-tm);
  Double_t blph = ph/(1-bm);
  Double_t canw = lpw+rpw;
  Double_t canh = tlph+blph;
  
  Double_t MpxT = tm*tlph;
  Double_t MpxB = bm*blph;
  Double_t trph = ph+MpxT+MpxB;
  Double_t brph = canh-trph;
  Double_t MfrT = MpxT/trph;
  Double_t MfrB = MpxB/trph;

  //Building canvas with pads
  TCanvas *c = new TCanvas("can","can",canw,canh);
  c->SetMargin(0,0,0,0);
  //Top right pad
  TPad *tptl = new TPad("tl","tl",0,blph/canh,lpw/canw,1);
  c->cd();
  tptl->SetMargin(lm,0,0,tm);
  tptl->Draw();
  c->cd();
  TPad *tpbl = new TPad("bl","bl",0,0,lpw/canw,blph/canh);
  tpbl->Draw();
  tpbl->SetMargin(lm,0,bm,0);
  c->cd();
  TPad *tptr = new TPad("tr","tr",lpw/canw,brph/canh,1,1);
  tptr->Draw();
  tptr->SetMargin(0,rm,MfrB,MfrT);
  c->cd();
  TPad *tpbr = new TPad("br","br",lpw/canw,0,1,brph/canh);
  tpbr->Draw();
  c->cd();
  TPad *mask = new TPad("mask","mask",0,0,1,1);
  mask->Draw();
  mask->SetFillStyle(0);
  return c;
};

TCanvas *BuildCanvasMatrix(Int_t Nx=2, Int_t Ny=2, Double_t ph=500, Double_t pw=500, Double_t lm=0.2, Double_t bm=0.2,Double_t tm=0.04, Double_t rm=0.04) {
  Double_t *wids = new Double_t[Nx];
  Double_t *higs = new Double_t[Ny];
  Double_t cwid=0;
  Double_t chig=0;
  for(Int_t i=0;i<Nx;i++) {
    Double_t oamar = 0;
    if(!i) oamar+=lm;
    if(i==Nx-1) oamar+=rm;
    wids[i] = pw/(1-oamar);
    cwid+=wids[i];
  };
  for(Int_t i=0;i<Ny;i++) {
    Double_t oamar = 0;
    if(!i) oamar+=tm;
    if(i==Ny-1) oamar+=bm;
    higs[i] = ph/(1-oamar);
    chig+=higs[i];
  };
  TCanvas *c = new TCanvas("can","can",cwid,chig);
  c->SetMargin(0,0,0,0);
  for(Int_t i=0;i<Nx;i++) {
    Double_t xof=0;
    for(Int_t k=0;k<i;k++) xof+=wids[k];
    for(Int_t j=0;j<Ny;j++) {
      Double_t yof=chig;
      for(Int_t k=0;k<=j;k++) yof-=higs[k];
      if(yof<0) yof=0;
      Double_t x1 = xof/cwid;
      Double_t y1 = yof/chig;
      Double_t x2 = (xof+wids[i])/cwid;
      Double_t y2 = (yof+higs[j])/chig;
      TPad *tp = new TPad(Form("tp_%i_%i",i,j),Form("tp_%i_%i",i,j),x1,y1,x2,y2);
      tp->SetMargin(0,0,0,0);
      if(!i) tp->SetLeftMargin(lm);
      if(i==Nx-1) tp->SetRightMargin(rm);
      if(!j) tp->SetTopMargin(tm);
      if(j==Ny-1) tp->SetBottomMargin(bm);
      c->cd();
      tp->Draw();
    };
  };
  TPad *mask = new TPad("mask","mask",0,0,1,1);
  c->cd();
  mask->SetFillStyle(0);
  mask->Draw();
  return c;
};
TCanvas *BuildSpagethiCanvas(Int_t ny, Double_t pw=500, Double_t th=500, Double_t bh=200, Double_t tm=0.04, Double_t bm=0.15, Double_t lm=0.15, Double_t rm=0.04) {
  //  Double_t tph = th/(1-tm); //top pad height
  //  Double_t bph = bh/(1-bm); //bottom pad height
  Double_t canHeight = th;
  for(Int_t i=0;i<ny-1;i++) canHeight+=bh;
  canHeight = canHeight/(1-tm-bm);
  Double_t tmabs = canHeight*tm;
  Double_t tph = th+tmabs; //Top pad height
  Double_t tmfrac = tmabs/tph; //Fraction for the top margin
  Double_t bmabs = canHeight*bm;
  Double_t bph = bh+bmabs; //Bottom pad height
  Double_t bmfrac = bmabs/bph; //Fraction for the bottom margin

  Double_t canWidth = (pw)/(1-lm-rm);
  TCanvas *c = new TCanvas("can","can",canWidth,canHeight);
  Double_t ypos[ny+1];
  ypos[0] = 1;
  ypos[1] = 1-(tph/canHeight);
  for(Int_t i=1; i<ny-1;i++) 
    ypos[i+1] = ypos[i] - bh/canHeight;
  ypos[ny] = 0;
  for(Int_t i=0;i<ny;i++) {
    c->cd();
    TPad *tp = new TPad(Form("tp_%i",i),Form("tp_%i",i),0,ypos[i+1],1,ypos[i]);
    tp->SetMargin(lm,rm,0,0);
    if(!i) tp->SetTopMargin(tmfrac);
    if(i==ny-1) tp->SetBottomMargin(bmfrac);
    tp->Draw();
  };
  c->cd();
  TPad *mask = new TPad("mask","mask",0,0,1,1);
  mask->SetFillStyle(0);
  mask->Draw();
  return c;
};
TCanvas *BuildParallelWithRatios(Double_t ph=500, Double_t pw=500, Double_t rr=0.33, Double_t lm=0.15, Double_t bm=0.13, Double_t tm=0.04, Double_t rm=0.04) {
  Double_t p1w = pw/(1-lm);
  Double_t p2w = pw;
  Double_t p3w = pw/(1-rm);
  Double_t cwidth = p1w+p2w+p3w;
  Double_t toph = ph/(1-tm);
  Double_t both = (ph*rr)/(1-bm);
  Double_t cheight = toph+both;
  Double_t DivPointX1 = p1w/cwidth;
  Double_t DivPointX2 = (p1w+p2w)/cwidth;
  Double_t DivPointY = both/cheight;

  TCanvas *c = new TCanvas("c","c",cwidth,cheight);
  c->SetMargin(0,0,0,0);
  c->cd();
  TPad *tlp = new TPad("tl","tl",0,DivPointY,DivPointX1,1);
  tlp->SetMargin(lm,0.,0,tm);
  tlp->Draw();
  c->cd();
  TPad *tmp = new TPad("tm","tm",DivPointX1,DivPointY,DivPointX2,1);
  tmp->SetMargin(0.,0,0,tm);
  tmp->Draw();
  c->cd();
  TPad *trp = new TPad("tr","tr",DivPointX2,DivPointY,1,1);
  trp->SetMargin(0.,rm,0,tm);
  trp->Draw();
  c->cd();
  TPad *blp = new TPad("bl","bl",0,0,DivPointX1,DivPointY);
  blp->SetMargin(lm,0,bm,0);
  blp->Draw();
  c->cd();
  TPad *bmp = new TPad("bm","bm",DivPointX1,0,DivPointX2,DivPointY);
  bmp->SetMargin(0,0,bm,0);
  bmp->Draw();
  c->cd();
  TPad *brp = new TPad("br","br",DivPointX2,0,1,DivPointY);
  brp->SetMargin(0,rm,bm,0);
  brp->Draw();
  c->cd();
  TPad *mask = new TPad("mask","mask",0,0,1,1);
  mask->SetFillStyle(0);
  mask->Draw();
  return c;
};
#endif
