/*Notes to myself:
  Removed the l_covConst based rebinning from CalculateCovariance in favor of a global PresetWeights function.
  This is b/c I'm not sure whether this should be used for MH or not (In particular, cov23 rebinning is done using w^3, PCC w/ w^{2} and other terms with w)
*/
#include "PCCContainer.h"
using namespace PCCSpace;
PCCContainer::PCCContainer():
  TNamed("",""),
  fInitialized(kFALSE),
  fObjs(0),
  fSysts(0),
  fBootstrapMean(kFALSE),
  fUseFSRebin(kFALSE),
  fNrb(0),
  fBrb(0),
  fSystAvgRange(0),
  fMaxAllowedSyst(1),
  fMaxSyst(0),
  fSystName(""),
  fMultiDist(0),
  fMeanNch(0),
  fFC(0),
  fck(0),
  fckbs(0),
  fcov2(0),
  fcov3(0),
  fcov23(0),
  fcov2nopt(0),
  fcov3nopt(0),
  fcov23nopt(0),
  fmpt(0),
  fNRebinForSyst(0),
  fRecGen(0),
  fNchMC(0),
  fNchCorrFunc(0),
  f_mpt_tl(0),
  f_cov_tl(0),
  fxTitle("V0M")//"#it{N}_{ch} (|#eta|<0.8)")
{
};
PCCContainer::PCCContainer(TString nomFile, Int_t MaxSyst): //Constructor to read file (nominal + systematics)
  TNamed("",""),
  fInitialized(kFALSE),
  fObjs(0),
  fSysts(0),
  fBootstrapMean(kFALSE),
  fUseFSRebin(kFALSE),
  fNrb(0),
  fBrb(0),
  fSystAvgRange(0),
  fMaxAllowedSyst(1),
  fMaxSyst(MaxSyst),
  fSystName(""),
  fMultiDist(0),
  fMeanNch(0),
  fFC(0),
  fck(0),
  fckbs(0),
  fcov2(0),
  fcov3(0),
  fcov23(0),
  fcov2nopt(0),
  fcov3nopt(0),
  fcov23nopt(0),
  fmpt(0),
  fNRebinForSyst(0),
  fRecGen(0),
  fNchMC(0),
  fNchCorrFunc(0),
  f_mpt_tl(0),
  f_cov_tl(0),
  fxTitle("V0M")//"#it{N}_{ch} (|#eta|<0.8)")
{
  ConstructContainer(nomFile);
  // if(!systFile.IsNull()) ReadSystematics(systFile);
};
PCCContainer::~PCCContainer()
{
  delete fObjs;
  delete fSysts;
  delete fRecGen;
  delete fNchMC;
  delete fNchCorrFunc;
  delete fMultiDist;
  delete f_mpt_tl;
  delete f_cov_tl;
  if(fFC) delete fFC;
};
TFile *PCCContainer::getFile(TString fina) {
  if(fina.IsNull()) { printf("File name not specified\n"); return 0; };
  TFile *retf = new TFile(fina.Data());
  if(!retf || retf->IsZombie()) { delete fina; printf("File %s not opened properly!\n",fina.Data()); return 0; };
  return retf;
}
TObject *PCCContainer::getObj(TFile *infi, TString objName) {
  if(!infi) {printf("File not provided!"); return 0; };
  TObject *retObj = (TObject*)infi->Get(objName.Data());
  if(!retObj) {printf("Could not find %s in file %s\n",objName.Data(),infi->GetName()); return 0; };
  return retObj;

};
void PCCContainer::HistOper(TH1 *reth, Double_t cvl, Double_t lpower, Bool_t Abs, Bool_t RemoveErrors) {
  for(Int_t i=1;i<=reth->GetNbinsX();i++) {
    Double_t bc = reth->GetBinContent(i);
    if(!bc) continue;
    bc+=cvl;
    Double_t bcn = TMath::Power(bc,lpower);
    if(Abs) bcn=TMath::Abs(bcn);
    reth->SetBinContent(i,bcn);
    if(RemoveErrors) { reth->SetBinError(i,0); }
    else {
      Double_t err = reth->GetBinError(i);
      err = err*lpower*TMath::Power(bc,lpower-1);
      reth->SetBinError(i,err);
    }
  }
}
TH1 *PCCContainer::HistSqrt(TH1 *inh, Bool_t RemoveErrors) {
  TH1 *reth = (TH1*)inh->Clone(Form("%s_sqrt",inh->GetName()));
  reth->Reset();
  for(Int_t i=1;i<=inh->GetNbinsX();i++) {
    Double_t bc = inh->GetBinContent(i);
    if(!bc || bc<0) continue;
    Double_t be = inh->GetBinError(i);
    if(RemoveErrors) be=0;
    bc = TMath::Sqrt(bc);
    be = be/(2*bc);
    reth->SetBinContent(i,bc);
    reth->SetBinError(i,be);
  }
  return reth;
}
void PCCContainer::ApplyBarlow(TH1 *target, TH1 *barlow) {
  if(!target || !barlow) {printf("ApplyBarlow: missing a histogram...\n"); return; };
  if(target->GetNbinsX()!=barlow->GetNbinsX()) { printf("ApplyBarlow: different number of bins\n"); return; };
  for(Int_t i=1;i<=target->GetNbinsX();i++) {
    if(barlow->GetBinContent(i)>1) continue;
    target->SetBinContent(i,0);
  };
}
void PCCContainer::AddInQuad(TH1 *h1, TH1 *h2) {
  if(!h1 || !h2) {printf("AddInQuad: One of the histograms is missing!\n"); return; };
  if(h1->GetNbinsX()!=h2->GetNbinsX()) { printf("AddInQuad: number of bins in provided histograms are not the same!\n"); return; };
  for(Int_t i=1;i<=h1->GetNbinsX();i++) {
    Double_t bc1 = h1->GetBinContent(i);
    Double_t bc2 = h2->GetBinContent(i);
    Double_t be1 = h1->GetBinError(i);
    Double_t be2 = h2->GetBinError(i);
    if(bc1==0 && bc2==0) continue;
    h1->SetBinContent(i,TMath::Sqrt(bc1*bc1+bc2*bc2));
    h1->SetBinError(i,TMath::Sqrt(be1*be1+be2*be2));
  }
}
void PCCContainer::ApplyErrors(TH1 *inh, TH1 *errors, Bool_t relative) {
  if(!inh || !errors) {printf("ApplyErrors: One of the histograms is missing!\n"); return; };
  if(inh->GetNbinsX()!=errors->GetNbinsX()) { printf("ApplyErrors: number of bins in provided histograms are not the same (%i vs %i)!\n",inh->GetNbinsX(),errors->GetNbinsX()); return; }; // && fNRebinForSyst<2
  for(Int_t i=1;i<=inh->GetNbinsX();i++) {
    Double_t bc1 = inh->GetBinContent(i);
    Double_t bce = inh->GetBinError(i);
    if(bc1==0 && bce==0) continue;
    Double_t bc2;
    if(fNRebinForSyst>1) {
      Double_t bcent = inh->GetBinCenter(i);
      Int_t tind = errors->FindBin(bcent);
      bc2 = errors->GetBinContent(tind);
    } else bc2 = errors->GetBinContent(i);
    bce = relative?(bc1*bc2):bc2;
    inh->SetBinError(i,bce);
  }
}
void PCCContainer::ApplyErrors(TH1 *inh, TF1 *errors, Bool_t relative) {
  if(!inh || !errors) {printf("ApplyErrors: Input histogram or TF1 is missing!\n"); return; };
  for(Int_t i=1;i<=inh->GetNbinsX();i++) {
    Double_t bc1 = inh->GetBinContent(i);
    Double_t bce = inh->GetBinError(i);
    if(bc1==0 && bce==0) continue;
    Double_t bc2 = errors->Eval(inh->GetBinCenter(i));
    bce = relative?(bc1*bc2):bc2;
    inh->SetBinError(i,bce);
  }
}

Bool_t PCCContainer::ConstructContainer(TString fina) {
  fInitialized = kFALSE;
  TFile *infi = getFile(fina);
  if(!infi) return kFALSE;
  ConstructContainer(infi,0);
  for(Int_t i=0;i<fMaxSyst;i++) {
    PCCContainer *systCont = new PCCContainer();
    Bool_t c_found = systCont->ConstructContainer(infi,i);
    if(!c_found) { delete systCont; continue; };//If not found, then delete to not keep trash around
    systCont->SetName(Form("PCCCont_Syst%i",i));
    if(!fSysts) { fSysts = new TList(); fSysts->SetOwner(kTRUE); };
    fSysts->Add(systCont);
  }
  infi->Close();
  delete infi;
  fInitialized=kTRUE;
  return kTRUE;
}
Bool_t PCCContainer::ConstructContainer(TFile *infi, Int_t lInd) {
  // printf("Constructor called!\n");
  if(!infi) return 0;
  fFC = (AliGFWFlowContainer*)getObj(infi,Form("FlowContainer__%i",lInd));
  if(!fFC) fFC = (AliGFWFlowContainer*)getObj(infi,Form("FlowContainer_%i",lInd));
  if(!fFC) fFC = (AliGFWFlowContainer*)getObj(infi,"FlowContainer");
  if(!fFC) return 0;
  f_mpt_tl = (TList*)getObj(infi,Form("MPTDiff_%i",lInd));
  if(!f_mpt_tl) f_mpt_tl = (TList*)getObj(infi,"MPTDiff");
  if(!f_mpt_tl) {delete fFC; return 0; };
  if(f_mpt_tl->FindObject("ckcont_ch")) fck = (AliCkContainer*)f_mpt_tl->FindObject("ckcont_ch");
                             else fckbs = (AliProfileBS*)f_mpt_tl->FindObject("varpt_ch");
  if(!fck && !fckbs) { delete fFC; return 0; };
  if(fckbs) fmpt = (AliProfileBS*)f_mpt_tl->FindObject("MeanPtClosure_ch");
       else fmpt = (AliProfileBS*)fck->getObsList()->At(2);
  if(fck) fck->SetTitle(Form(";%s; #it{c}_{K}",fxTitle.Data()));
  else fckbs->SetTitle(Form(";%s; #it{c}_{K}",fxTitle.Data()));

  f_cov_tl = (TList*)getObj(infi,Form("Covariance_%i",lInd));
  if(!f_cov_tl) f_cov_tl = (TList*)getObj(infi,"Covariance");
  if(!f_cov_tl) { delete fFC; delete fck; delete fckbs; delete f_mpt_tl; return 0; };
  if(fckbs) { //If we are working with the "old" setup
    fcov2 = (AliProfileBS*)f_cov_tl->FindObject("cov_ch");
    fcov3 = (AliProfileBS*)f_cov_tl->FindObject("cov_v3_ch");
    fcov23 = (AliProfileBS*)f_cov_tl->FindObject("cov_v23_ch");
    if(!fcov2 || !fcov3 || !fcov23) { delete fcov2; delete fcov3; delete fcov23; };
  } else { //For newer setup -- would have been so much more simple if I kept the naming didn't change!
    fcov2 = (AliProfileBS*)f_cov_tl->FindObject("covmpt_ch");
    fcov3 = (AliProfileBS*)f_cov_tl->FindObject("covmpt_v3_ch");
    fcov23 = (AliProfileBS*)f_cov_tl->FindObject("covmpt_v23_ch");
    fcov2nopt = (AliProfileBS*)f_cov_tl->FindObject("covnopt_ch");
    fcov3nopt = (AliProfileBS*)f_cov_tl->FindObject("covnopt_v3_ch");
    fcov23nopt = (AliProfileBS*)f_cov_tl->FindObject("covnopt_v23_ch");
    if(!fcov2 || !fcov3 || !fcov23 || !fcov2nopt || !fcov3nopt || !fcov23nopt ) { delete fcov2; delete fcov3; delete fcov23; delete fcov2nopt; delete fcov3nopt; delete fcov23nopt; };
  }
  if(!fcov2) { delete fFC; delete fck; delete fckbs; return 0;};
  fInitialized=kTRUE;
  return kTRUE;
}
void PCCContainer::PresetWeights(AliProfileBS *bs) {
  if(fcov2) fcov2->PresetWeights(bs);
  if(fcov23) fcov23->PresetWeights(bs);
  if(fcov3) fcov3->PresetWeights(bs);
  if(fcov2nopt) fcov2nopt->PresetWeights(bs);
  if(fcov23nopt) fcov23nopt->PresetWeights(bs);
  if(fcov3nopt) fcov3nopt->PresetWeights(bs);
  if(fmpt) fmpt->PresetWeights(bs);
}
Bool_t PCCContainer::OverrideMeanPt(TString infi, Int_t lInd) {
  //Obviously, obsolete, when working with the AliCkContainer
  // TFile *tf = new TFile(infi.Data());
  // TList *tl = (TList*)getObj(tf,Form("MPTProfileList_%i",lInd));
  // fmpt = (TProfile*)tl->At(0);
  // fmpt = (TProfile*)fmpt->Clone("MeanPt");
  // fmpt->SetDirectory(0);
  // tf->Close();
  return kTRUE;
}
Ssiz_t PCCContainer::findStrInd(TString instr, TString delim) {
  Ssiz_t ind1=-1;
  Ssiz_t ind2=-1;
  Bool_t found=kFALSE;
  do {
    ind1=ind2+1;
    ind2 = instr.Index(delim.Data(),ind1);
    if(ind2>-1) found = kTRUE;
  } while (ind2>=0);
  if(!found) return -1;
  return ind1;
}
/*void PCCContainer::RebinMulti(Int_t nrb) { //with FS rebinning, this is not relevant
  if(!fFC || !fck || !fcov2) return;
  fFC->RebinMulti(nrb);
  fck->RebinMulti(nrb);
  // fck->CalculateCks();
  fcov2->RebinMulti(nrb);
  fcov3->RebinMulti(nrb);
  fcov23->RebinMulti(nrb);
  fcov2nopt->RebinMulti(nrb);
  fcov3nopt->RebinMulti(nrb);
  fcov23nopt->RebinMulti(nrb);
  // fmpt->RebinMulti(nrb); //Rebinned during the fck->Rebin
  TObjArray *oba = fFC->GetSubProfiles();
  if(oba)
    for(Int_t i=0;i<oba->GetEntries();i++)
      ((TProfile2D*)oba->At(i))->RebinX(nrb);
  if(fSysts) {
    for(Int_t i=0;i<fSysts->GetEntries();i++) ((PCCContainer*)fSysts->At(i))->RebinMulti(nrb);
  }
  MakeNchMC(); //Only done if fRecGen is set
}*/
void PCCContainer::RebinMulti(Int_t nrb, Double_t *chbins) {
  //New approach:
  fNrb = nrb;
  fBrb = chbins;
  MakeNchMC();
  if(fMultiDist && chbins) {
    if(fMeanNch) delete fMeanNch;
    fMeanNch = new TProfile("Mean_Nch_in_bin",";#it{N}_{ch}; #LT#it{N}_{ch}#GT",nrb,chbins);
    fMeanNch->SetDirectory(0);
    for(Int_t i=1;i<=fMultiDist->GetNbinsX();i++) fMeanNch->Fill(fMultiDist->GetBinCenter(i),fMultiDist->GetBinCenter(i),fMultiDist->GetBinContent(i));
  }
  if(fSysts) {
    for(Int_t i=0;i<fSysts->GetEntries();i++) ((PCCContainer*)fSysts->At(i))->RebinMulti(nrb,chbins);
  }
  if(fUseFSRebin) return;
  if(!fFC || !fck || !fcov2) return;
  if(fFC) fFC->SetMultiRebin(nrb,chbins);
  if(fck) {
    fck->RebinMulti(nrb,chbins);
    fck->CalculateCks();
  };
  if(fcov2) fcov2->RebinMulti(nrb,chbins);
  if(fcov23) fcov23->RebinMulti(nrb,chbins);
  if(fcov3) fcov3->RebinMulti(nrb,chbins);
  if(fcov2nopt) fcov2nopt->RebinMulti(nrb,chbins);
  if(fcov23nopt) fcov23nopt->RebinMulti(nrb,chbins);
  if(fcov3nopt) fcov3nopt->RebinMulti(nrb,chbins);
  if(fmpt) fmpt->RebinMulti(nrb,chbins);
  // MakeNchMC(); //Only done if fRecGen is set
  // if(fMultiDist) {
  //   if(fMeanNch) delete fMeanNch;
  //   fMeanNch = new TProfile("Mean_Nch_in_bin",";#it{N}_{ch}; #LT#it{N}_{ch}#GT",nrb,chbins);
  //   fMeanNch->SetDirectory(0);
  //   for(Int_t i=1;i<=fMultiDist->GetNbinsX();i++) fMeanNch->Fill(fMultiDist->GetBinCenter(i),fMultiDist->GetBinCenter(i),fMultiDist->GetBinContent(i));
  // }
}
TH1 *PCCContainer::CalculateCovariance(const Int_t &nrb, AliProfileBS* l_covConst, AliProfileBS* l_covLin, AliProfileBS* l_mpt, AliProfileBS *rbWeights) {
  if(!rbWeights) rbWeights = l_covConst;
  if(!fUseFSRebin) {
    l_mpt->PresetWeights(rbWeights);
  //   l_covConst->PresetWeights(rbWeights);
  //   l_covLin->PresetWeights(rbWeights);
  };
  TH1 *reth = l_covConst->getHist(nrb);
  if(!l_covLin) return reth;
  TH1 *rethnopt = l_covLin->getHist(nrb);
  TH1 *mpt = l_mpt->getHist(nrb);
  if(!fUseFSRebin) {
    l_mpt->PresetWeights(0);
  //   l_covConst->PresetWeights(0);
  //   l_covLin->PresetWeights(0);
  };
  rethnopt->Multiply(mpt);
  reth->Add(rethnopt,-1);
  delete rethnopt;
  delete mpt;
  return reth;
}
TH1 *PCCContainer::getCov2(Int_t nrb, Bool_t bootstrap) {
  TH1 *reth = CalculateCovariance(nrb,fcov2,fcov2nopt,fmpt);
  if(!bootstrap) return reth;
  return getBootstrapped(reth,kCov2);
}
TH1 *PCCContainer::getCov3(Int_t nrb, Bool_t bootstrap) {
  TH1 *reth = CalculateCovariance(nrb,fcov3,fcov3nopt,fmpt);
  if(!bootstrap) return reth;
  return getBootstrapped(reth,kCov3);
}
TH1 *PCCContainer::getCov23(Int_t nrb, Bool_t bootstrap) {
  TH1 *reth = CalculateCovariance(nrb,fcov23,fcov23nopt,fmpt);
  if(!bootstrap) return reth;
  return getBootstrapped(reth,kCov23);
}
TH1 *PCCContainer::getCk(Int_t nrb, Bool_t bootstrap) {
  if(!fck && !fckbs) {printf("Ck could not be read. Try reinitializing the PCCContainer!\n"); return 0; };
  TH1 *reth = fck?fck->getHist(nrb):fckbs->getHist(nrb);
  //Need to clone, b/c this is not calculated on the fly, and then gets deleted in subsequent methods. Should be implemented in AliCkContainer instead, but a quickfix for now
  reth = (TH1*)reth->Clone(Form("retCk_%s",reth->GetName()));
  if(!bootstrap) return reth;
  return getBootstrapped(reth,kCk);
}
TH1 *PCCContainer::getVar2(Int_t nrb, Bool_t bootstrap) {
  TH1 *c24 = getVarNHar(nrb,kC22,kC24);
  c24->SetTitle(Form(";%s; Var(#it{v}_{2}^{2})",fxTitle.Data()));
  if(!bootstrap) return c24;
  return getBootstrapped(c24,kVar2);
}
TH1 *PCCContainer::getVar3(Int_t nrb, Bool_t bootstrap) {
  TH1 *c34 = getVarNHar(nrb,kC32,kC34);
  c34->SetTitle(Form(";%s; Var(#it{v}_{3}^{2})",fxTitle.Data()));
  if(!bootstrap) return c34;
  return getBootstrapped(c34,kVar3);
}

TH1 *PCCContainer::getVarNHar(Int_t nrb, lFunc f_cn2, lFunc f_cn4) {
  TH1 *cn4 = (this->*f_cn4)(nrb,kFALSE);
  TH1 *cn2 = (this->*f_cn2)(nrb,kFALSE);
  cn2->Multiply(cn2);
  cn4->Add(cn2);
  cn4->SetName("VarV2");
  cn4->SetTitle(Form(";%s; Var(#it{v}_{n}^{2})",fxTitle.Data()));
  delete cn2;
  return cn4;
}

TH1 *PCCContainer::getVar3Sub(Int_t nrb, Bool_t bootstrap) {
  TH1 *c24 = getC243Sub(nrb,bootstrap);
  // TH1 *c24 = getC24(nrb,bootstrap);
  TH1 *c22 = getC22(nrb,bootstrap);
  c22->Multiply(c22);
  c24->Add(c22);
  c24->SetName("VarV2");
  c24->SetTitle(Form(";%s; Var(#it{v}_{2}^{2})",fxTitle.Data()));
  delete c22;
  if(!bootstrap) return c24;
  return getBootstrapped(c24,kVar3Sub);
}
TH1 *PCCContainer::getCN4(Int_t nrb, Int_t nHarmonic) {
  AliGFWFlowContainer *fClone = (nrb<0)?fFC:(AliGFWFlowContainer*)fFC->Clone("TempClone");
  if(nrb>=0) {
    fClone->OverrideMainWithSub(nrb,kFALSE);
    Int_t nxb=0;
    Double_t *xbs = fFC->GetMultiRebin(nxb);
    fClone->SetMultiRebin(nxb,xbs);
  };
  fClone->SetIDName("ChFull");
  TH1 *reth = fClone->GetCN4VsMulti(nHarmonic,0);
  if(nrb>=0) delete fClone;
  return reth;
}
TH1 *PCCContainer::getCN2(Int_t nrb, Int_t nHarmonic) {
  AliGFWFlowContainer *fClone = (nrb<0)?fFC:(AliGFWFlowContainer*)fFC->Clone("TempClone");
  if(nrb>=0) {
    fClone->OverrideMainWithSub(nrb,kFALSE);
    Int_t nxb;
    Double_t *xbs = fFC->GetMultiRebin(nxb);
    fClone->SetMultiRebin(nxb,xbs);
  };
  fClone->SetIDName("ChGap");
  TH1 *reth = fClone->GetCN2VsX(nHarmonic,kFALSE,0);
  if(nrb>=0) delete fClone;
  return reth;
}
TH1 *PCCContainer::getC22(Int_t nrb, Bool_t bootstrap) {
  TH1 *reth = getCN2(nrb,2);
  if(!bootstrap) return reth;
  return getBootstrapped(reth,kC22);
}
TH1 *PCCContainer::getC32(Int_t nrb, Bool_t bootstrap) {
  TH1 *reth = getCN2(nrb,3);
  if(!bootstrap) return reth;
  return getBootstrapped(reth,kC32);
}
TH1 *PCCContainer::getC24(Int_t nrb, Bool_t bootstrap) {
  TH1 *reth = getCN4(nrb,2);
  if(!bootstrap) return reth;
  return getBootstrapped(reth,kC24);
}
TH1 *PCCContainer::getC34(Int_t nrb, Bool_t bootstrap) {
  TH1 *reth = getCN4(nrb,3);
  if(!bootstrap) return reth;
  return getBootstrapped(reth,kC34);
}

TH1 *PCCContainer::getC243Sub(AliGFWFlowContainer *infc, Int_t nrb, Bool_t bootstrap, TString ms, TString conf1, TString conf2) {
  infc->SetIDName(ms.Data());
  TProfile *mainconf = infc->GetCorrXXVsMulti("24",0);
  infc->SetIDName(conf1.Data());
  TProfile *c1 = infc->GetCorrXXVsMulti("22",0);
  c1->SetName("pf221");
  infc->SetIDName(conf2.Data());
  TProfile *c2 = infc->GetCorrXXVsMulti("22",0);
  c1->SetName("pf222");
  TH1 *m4 = ProfToHist(mainconf);
  delete mainconf;
  TH1 *m21 = ProfToHist(c1);
  delete c1;
  TH1 *m22 = ProfToHist(c2);
  delete c2;
  m21->Multiply(m22);
  m4->Add(m21,-2);
  delete m21;
  delete m22;
  m4->SetName(Form("c24_%s",ms.Data()));
  return m4;
}
TH1 *PCCContainer::getC243Sub(Int_t nrb, Bool_t bootstrap) {
  AliGFWFlowContainer *fClone = (nrb<0)?fFC:(AliGFWFlowContainer*)fFC->Clone("TempClone");
  if(nrb>=0) {
    fClone->OverrideMainWithSub(nrb,kFALSE);
    Int_t nxb;
    Double_t *xbs = fFC->GetMultiRebin(nxb);
    fClone->SetMultiRebin(nxb,xbs);
  };
  TH1 *conf1 = getC243Sub(fClone,nrb,bootstrap,"LLMR","LM","LR");
  TH1 *conf2 = getC243Sub(fClone,nrb,bootstrap,"LMMR","LM","MR");
  TH1 *conf3 = getC243Sub(fClone,nrb,bootstrap,"LMRR","LR","MR");
  if(nrb>=0) delete fClone;
  conf1->SetBit(TH1::kIsAverage);
  conf2->SetBit(TH1::kIsAverage);
  conf3->SetBit(TH1::kIsAverage);
  conf1->Add(conf2);
  conf1->Add(conf3);
  delete conf2;
  delete conf3;
  conf1->SetName("Var24_3Sub");
  if(!bootstrap) return conf1;
  return getBootstrapped(conf1,kC243Sub);
}
TH1 *PCCContainer::getSigmaNHar(Int_t nrb, Bool_t bootstrap, lFunc varFunc) {
  TH1 *varv2 = (this->*varFunc)(nrb,bootstrap);
  TH1 *reth = HistSqrt(varv2,kFALSE);
  delete varv2;
  reth->SetTitle(Form(";%s; #sigma(#it{v}_{n}^{2})",fxTitle.Data()));
  return reth;
};
TH1 *PCCContainer::getSigma2(Int_t nrb, Bool_t bootstrap) {
  TH1 *reth = getSigmaNHar(nrb,bootstrap,kVar2);
  reth->SetTitle(Form(";%s; #sigma(#it{v}_{2}^{2})",fxTitle.Data()));
  return reth;
};
TH1 *PCCContainer::getSigma3(Int_t nrb, Bool_t bootstrap) {
  TH1 *reth = getSigmaNHar(nrb,bootstrap,kVar3);
  reth->SetTitle(Form(";%s; #sigma(#it{v}_{3}^{2})",fxTitle.Data()));
  return reth;
};
TH1 *PCCContainer::getSigma3Sub(Int_t nrb, Bool_t bootstrap) {
  TH1 *varv2 = getVar3Sub(nrb,bootstrap);
  TH1 *reth = HistSqrt(varv2,kFALSE);
  delete varv2;
  reth->SetTitle(Form(";%s; #sigma(#it{v}_{2}^{2})",fxTitle.Data()));
  reth->SetBit(TH1::kIsAverage, kFALSE);
  return reth;
};
TH1 *PCCContainer::getPCC2(Int_t nrb, Bool_t bootstrap) {
  TH1 *reth = getPCCNHar(nrb,kCov2,kSigma2);
  reth->SetTitle(Form(";%s; #rho(v_{2}^{2}, [#it{p}_{T}])",fxTitle.Data()));
  if(!bootstrap) return reth;
  return getBootstrapped(reth,kPCC2);
}
TH1 *PCCContainer::getPCC3(Int_t nrb, Bool_t bootstrap) {
  TH1 *reth = getPCCNHar(nrb,kCov3,kSigma3);
  reth->SetTitle(Form(";%s; #rho(v_{3}^{2}, [#it{p}_{T}])",fxTitle.Data()));
  if(!bootstrap) return reth;
  return getBootstrapped(reth,kPCC3);
}
TH1 *PCCContainer::getPCCNHar(Int_t nrb, lFunc covFunc, lFunc sigmaFunc) {
  TH1 *reth = (this->*covFunc)(nrb,0);
  TH1 *sigmavn = (this->*sigmaFunc)(nrb,0);
  TH1 *ck = (this->*kCk)(nrb,0);
  TH1 *sqrtCk = HistSqrt(ck);
  reth->Divide(sigmavn);
  reth->Divide(sqrtCk);
  delete sigmavn;
  delete ck;
  delete sqrtCk;
  reth->SetName("PCC");
  reth->SetTitle(Form(";%s; #rho(v_{n}^{2}, [#it{p}_{T}])",fxTitle.Data()));
  return reth;
}

TH1 *PCCContainer::getSqCkOverMpt(Int_t nrb, Bool_t bootstrap) {
  TH1 *lck = getCk(nrb,bootstrap);
  TH1 *cksq = HistSqrt(lck,kFALSE);
  delete lck;
  cksq->SetName("SqrtCkOverMpt");
  TH1 *mpt = getMpt();
  cksq->Divide(mpt);
  delete mpt;
  return cksq;
}
TH1 *PCCContainer::getPCC3Sub(Int_t nrb, Bool_t bootstrap) {
  TH1 *reth = (this->*kCov2)(nrb,0);
  TH1 *sigmav2 = (this->*kSigma3Sub)(nrb,0);
  TH1 *ck = (this->*kCk)(nrb,0);
  TH1 *sqrtCk = HistSqrt(ck);
  reth->Divide(sigmav2);
  reth->Divide(sqrtCk);
  delete sigmav2;
  delete ck;
  delete sqrtCk;
  reth->SetName("PCC");
  reth->SetTitle(Form(";%s; #rho(v_{2}^{2}, [#it{p}_{T}])",fxTitle.Data()));
  if(!bootstrap) return reth;
  return getBootstrapped(reth,kPCC3Sub);
}
TProfile *PCCContainer::getStatistics(lFunc sf) {
  Int_t nSubs = fcov2->getNSubs();
  TH1 **htmp = new TH1*[nSubs];
  for(Int_t i=0;i<nSubs;i++) htmp[i] = (this->*sf)(i,kFALSE);
  Int_t nxbins = htmp[0]->GetNbinsX();
  Double_t *xbins = new Double_t[nxbins+1];
  htmp[0]->GetXaxis()->GetLowEdge(xbins);
  xbins[nxbins] = htmp[0]->GetXaxis()->GetBinUpEdge(nxbins);
  TProfile *retpf = new TProfile("Statistics_Profile","StatProf",nxbins,xbins);
  for(Int_t i=1;i<=htmp[0]->GetNbinsX();i++)
    for(Int_t j=0;j<nSubs;j++)
      retpf->Fill(htmp[0]->GetBinCenter(i),htmp[j]->GetBinContent(i));
  for(Int_t i=0; i<nSubs;i++) delete htmp[i];
  delete [] htmp;
  delete [] xbins;
  return retpf;
}
TProfile *PCCContainer::getStatistics(lFunc sf, Int_t lSyst) {
  //Bootstrapping relative differences between nominal and systematic check
  Int_t nSubs = fcov2->getNSubs();
  TH1 **htmp = new TH1*[nSubs];
  for(Int_t i=0;i<nSubs;i++) { htmp[i] = (this->*sf)(i,kFALSE); htmp[i]->SetName(Form("Nom_%i",i)); };
  PCCContainer *trg = (PCCContainer*)fSysts->FindObject(Form("PCCCont_Syst%i",lSyst));
  TH1 **hsys = new TH1*[nSubs];
  for(Int_t i=0;i<nSubs;i++) hsys[i] = (trg->*sf)(i,kFALSE);
  for(Int_t i=0;i<nSubs;i++) { hsys[i]->Add(htmp[i],-1); hsys[i]->Divide(htmp[i]); };
  Int_t nxbins = htmp[0]->GetNbinsX();
  Double_t *xbins = new Double_t[nxbins+1];
  htmp[0]->GetXaxis()->GetLowEdge(xbins);
  xbins[nxbins] = htmp[0]->GetXaxis()->GetBinUpEdge(nxbins);
  TProfile *retpf = new TProfile("Statistics_Profile","StatProf",nxbins,xbins);
  for(Int_t i=1;i<=hsys[0]->GetNbinsX();i++)
    for(Int_t j=0;j<nSubs;j++)
      retpf->Fill(htmp[0]->GetBinCenter(i),hsys[j]->GetBinContent(i));
  for(Int_t i=0; i<nSubs;i++) { delete htmp[i]; delete hsys[i]; };
  delete [] htmp;
  delete [] hsys;
  delete [] xbins;
  return retpf;
}
TH2 *PCCContainer::getStatisticsDist(lFunc sf,Int_t nybins, Double_t ymin, Double_t ymax) {
  Int_t nSubs = fcov2->getNSubs();
  TH1 **htmp = new TH1*[nSubs];
  TH1 *hCent = (this->*sf)(-1,kFALSE);
  for(Int_t i=0;i<nSubs;i++) htmp[i] = (this->*sf)(i,kFALSE);
  Int_t nxbins = htmp[0]->GetNbinsX();
  Double_t *xbins = new Double_t[nxbins+1];
  htmp[0]->GetXaxis()->GetLowEdge(xbins);
  xbins[nxbins] = htmp[0]->GetXaxis()->GetBinUpEdge(nxbins);
  TH2 *reth = new TH2D("statDist","statDist",nxbins,xbins[0],xbins[nxbins],nybins,ymin,ymax);
  reth->GetXaxis()->Set(nxbins,xbins);
  for(Int_t i=1;i<=htmp[0]->GetNbinsX();i++)
    for(Int_t j=0;j<nSubs;j++)
      reth->Fill(htmp[0]->GetBinCenter(i),htmp[j]->GetBinContent(i));//-hCent->GetBinContent(i));
  for(Int_t i=0; i<nSubs;i++) delete htmp[i];
  delete [] htmp;
  delete [] xbins;
  delete hCent;
  return reth;
}

TH1 *PCCContainer::getSystematicsObs(lFunc sf, Int_t lSyst, Bool_t bootstrap) {
  if(sf==kDisabled) return 0;
  if(fNrb>0) bootstrap=kTRUE;
  //Rewriting so that I can add final state rebinning at the end
  TH1 *reth = 0;
  if(sf == kMpt || sf == kSQCK) reth = (this->*sf)(-1,0);
  else if(!lSyst && !fSysts) reth = (this->*sf)(-1,bootstrap);
  else {
    if(!fSysts) {printf("Systematics were not read!\n"); return 0; };
    PCCContainer *trg = (PCCContainer*)fSysts->FindObject(Form("PCCCont_Syst%i",lSyst));
    if(!trg) { printf("Systematics index %i not found!\n",lSyst); return 0; };
    reth = (trg->*sf)(-1,bootstrap);
  }
  if(!reth) { printf("Ops! Went through the whole getSystematicsObs and did not fetch what you've asked me to :(\n"); return 0; };
  FSRebin(&reth);
  return reth;
  //Old method with no FS rebinning
  // if(sf== kMpt || sf == kSQCK) return (this->*sf)(-1,0);
  // // printf("Picking up systematics with index %i\n",lSyst);
  // if(!lSyst && !fSysts) return (this->*sf)(-1,bootstrap);
  // if(!fSysts) {printf("Systematics were not read!\n"); return 0; };
  // PCCContainer *trg = (PCCContainer*)fSysts->FindObject(Form("PCCCont_Syst%i",lSyst));
  // if(!trg) { printf("Systematics index %i not found!\n",lSyst); return 0; };
  // return (trg->*sf)(-1,bootstrap);
}
TH1 *PCCContainer::getSystematics(lFunc sf, Int_t lSyst, Bool_t rel, Bool_t bootstrap, Bool_t applyBarlow, Bool_t takeAbs) {
  if(sf==kDisabled) return 0;
  TH1 *mod = getSystematicsObs(sf,lSyst,bootstrap);//(trg->*sf)(-1,bootstrap);
  Double_t *relErrs = 0;
  if(bootstrap) {
    relErrs = new Double_t[mod->GetNbinsX()];
    for(Int_t i=1;i<=mod->GetNbinsX();i++) if(mod->GetBinContent(i)==0) continue; else relErrs[i-1] = mod->GetBinError(i)/mod->GetBinContent(i);
  }
  if(!mod) {printf("Could not fetch systematic index %i\n", lSyst); return 0; };
  mod->SetName("Modded");
  //If we have final state rebinning, then need bootstrapping on the nominal here. Otherwise, final state rebinning will be bogus
  TH1 *nom = 0;
  if(fUseFSRebin) {
    nom = (this->*sf)(-1,fNrb>0?kTRUE:bootstrap);
    FSRebin(&nom);
  } else nom = (this->*sf)(-1,bootstrap);
  mod->SetName(Form("RelSyst_%s_Ind%i",nom->GetName(),lSyst));
  mod->Add(nom,-1);
  if(rel)
    mod->Divide(nom);
  HistOper(mod,0,1,kTRUE,kTRUE);
  if(bootstrap) { for(Int_t i=1;i<=mod->GetNbinsX();i++) if(relErrs[i-1]==0) continue; else mod->SetBinError(i,relErrs[i-1]*(rel?1:TMath::Abs(mod->GetBinContent(i)))); };  //):(mod->GetBinContent(i)*relErrs[i-1])); };
  delete nom;
  if(!applyBarlow)
    return mod;
  TH1 *barlow = getBarlowTest(sf,lSyst,getCorrelationError(lSyst));
  if(!barlow) return mod;
  ApplyBarlow(mod,barlow);
  delete barlow;
  return mod;
}
TH1 *PCCContainer::getBarlowTest(lFunc sf, Int_t lSyst, Double_t corrFactor) {
  if(sf==kDisabled) return 0;
  TH1 *mod = getSystematicsObs(sf,lSyst,kTRUE);
  if(!mod) return 0;
  mod->SetName("ModdedBS");
  TH1 *nom = (this->*sf)(-1,kTRUE);
  nom->SetName("NominalBS");
  TH1 *absDif = getSystematics(sf,lSyst,kFALSE,kFALSE,kFALSE);
  if(!absDif) return 0;
  for(Int_t i=1;i<=absDif->GetNbinsX();i++) {
    Double_t bc = absDif->GetBinContent(i);
    if(!bc) continue;
    Double_t s1 = mod->GetBinError(i);
    Double_t s2 = nom->GetBinError(i);
    Double_t denom = TMath::Sqrt(TMath::Abs((s1*s1+corrFactor*s2*s2)));
    if(!denom) { absDif->SetBinError(i,0); continue; };
    absDif->SetBinContent(i,bc/denom);
    absDif->SetBinError(i,0);
  }
  delete nom;
  delete mod;
  return absDif;
}
Bool_t PCCContainer::BuildIndexMap(Bool_t force) {
  if(!fSysts) {printf("Systematics not read!\n"); return 0; };
  if(!fSysts->GetEntries()) {printf("No systematics available!\n"); return 0; };
  if(force && !((int)fIndMap.size())) fIndMap.clear();
  for(Int_t i=0;i<fMaxSyst;i++) {
    Int_t tval = getTarget(i);
    if(!fIndMap.count(tval)) fIndMap.insert(make_pair(tval,vector<int>{}));
    fIndMap.at(tval).push_back(i);
  };
  return kTRUE;
};
vector<TH1*> PCCContainer::getSystSubset(lFunc sf, Int_t KeyInd, Bool_t relative, Bool_t inclErrors, Bool_t applyBarlow, Bool_t takeAbs) {
  if(!(int)fIndMap.size()) if(!BuildIndexMap(kTRUE)) { printf("getSystSubset: could not build map!\n"); return nullvec; };
  if(!fIndMap.count(KeyInd)) { printf("getSystSubset: could not find key %i in the map!\n",KeyInd); return nullvec; };
  vector<int> inds = fIndMap.at(KeyInd);
  vector<TH1*> retvec;
  for(Int_t j=0;j<(int)inds.size();j++) {
    retvec.push_back(getSystematics(sf,inds.at(j),relative,inclErrors,applyBarlow,takeAbs));
  };
  return retvec;
}
vector<TH1*> PCCContainer::getMergedErrors(lFunc sf, Bool_t relative, Bool_t inclErrors, Bool_t applyBarlow, Bool_t weightByErrors) {
  if(!(int)fIndMap.size()) if(!BuildIndexMap(kTRUE)) { printf("getMergedErrors: could not build map!\n"); return nullvec; };
  vector<TH1*> reth;
  for(auto errcind=fIndMap.begin(); errcind!=fIndMap.end(); errcind++) {
    Int_t targetInd = errcind->first;
    if(targetInd==0) continue; //0 contains all those that are not relevant
    vector<TH1*> lSubsets = getSystSubset(sf,targetInd,relative,inclErrors,applyBarlow,!weightByErrors); //if we want to weight by errors, then we shouldn't take the absolute value here
    if((int)lSubsets.size()==0) continue; //if no subsets, then continue
    TH1 *rh = (TH1*)lSubsets.at(0)->Clone(Form("MergedErrors_%i",targetInd));
    if(!rh) continue;
    if(!inclErrors) RemoveErrors(rh); //To remove errors
    reth.push_back(rh);
    if(fSystOverride.count(targetInd)) {
      Double_t newerr = fSystOverride.at(targetInd);
      printf("Found a bins with key index %i, value %f\n",targetInd,newerr);
      for(Int_t i=1;i<=rh->GetNbinsX();i++) rh->SetBinContent(i,newerr);
      continue;
    }

    if((int)lSubsets.size()==1) continue; //if only one, then no need to continue;
    //loop over subsets and pick max
    for(auto sh = lSubsets.begin()+1; sh!=lSubsets.end(); sh++) {
      //Looping over bins
      for(Int_t i=1;i<=rh->GetNbinsX();i++) {
        if(inclErrors) {
          Double_t bc1 = rh->GetBinContent(i);
          Double_t bc2 = (*sh)->GetBinContent(i);
          Double_t be1 = rh->GetBinError(i);
          Double_t be2 = (*sh)->GetBinError(i);
          if(weightByErrors) {
            be1*=be1;
            be2*=be2;
            Double_t val1 = be1>0?bc1/be1:0;
            Double_t val2 = be2>0?bc2/be2:0;
            Double_t err1 = be1>0?1./be1:0;
            Double_t err2 = be2>0?1./be2:0;
            rh->SetBinContent(i,val1+val2);
            rh->SetBinError(i,err1+err2);
          } else {
            rh->SetBinContent(i,TMath::Max(bc1,bc2));
            rh->SetBinError(i,bc1>bc2?be1:be2);
          }
        } else {
          if(weightByErrors) {
            Double_t bc1 = rh->GetBinContent(i);
            Double_t bc2 = (*sh)->GetBinContent(i);
            Double_t be1 = rh->GetBinError(i);
            Double_t be2 = (*sh)->GetBinError(i);
            be1*=be1;
            be2*=be2;
            Double_t val1 = be1>0?bc1/be1:0;
            Double_t val2 = be2>0?bc2/be2:0;
            Double_t err1 = be1>0?1./be1:0;
            Double_t err2 = be2>0?1./be2:0;
            rh->SetBinContent(i,val1+val2);
            rh->SetBinError(i,err1+err2);
          } else {
            Double_t vmax = TMath::Max(rh->GetBinContent(i),(*sh)->GetBinContent(i));
            rh->SetBinContent(i,vmax);
          };
        };
      }
      (*sh)->Delete(); //Clean after myself
      if(weightByErrors) {
        for(Int_t i=1;i<=rh->GetNbinsX();i++) {
          Double_t bc = rh->GetBinContent(i);
          Double_t be = rh->GetBinError(i);
          if(be==0) continue;
          rh->SetBinContent(i,TMath::Abs(bc/be));
          rh->SetBinError(i,inclErrors?(1./TMath::Sqrt(be)):0);
        }
      }
    }
  }
  return reth;
}
TH1 *PCCContainer::getSystematicsSummed(lFunc sf, Bool_t relative, Bool_t inclErrors, Bool_t applyBarlow) {
  vector<TH1*> mergedErrors;
  if(fNRebinForSyst>1) {
    printf("Rebinning for systematics is not implemented yet. The standard rebinning is not valid here, and the FS rebinning isn't implemented yet.\n");
    printf("Will be crashing now... Sorry!\n");
    return 0;
/*    PCCContainer *tmpc = (PCCContainer*)this->Clone("TempPCCCont");
    tmpc->RebinMulti(fNRebinForSyst);
    mergedErrors = tmpc->getMergedErrors(sf,relative,applyBarlow);
    delete tmpc;*/
  } else mergedErrors = getMergedErrors(sf,relative,inclErrors,applyBarlow);
  if((int)mergedErrors.size()==0) {printf("getSystematicsSummed: could not find any errors to sum\n"); return 0; };
  TH1 *reth = (TH1*)mergedErrors.at(0)->Clone("MergedErrors_Summed");
  reth->Reset();
  for(auto ht = mergedErrors.begin(); ht!=mergedErrors.end();ht++) {
    AddInQuad(reth,*ht);
    (*ht)->Delete();
  };
  if(fSystAvgRange>0) {
    Int_t nbins = reth->GetNbinsX();
    Int_t nSteps = (Int_t)(nbins-(nbins%fSystAvgRange))/fSystAvgRange;
    if(nbins%fSystAvgRange) nSteps++;
    for(Int_t i=0;i<nSteps;i++) {
      Int_t iStart = reth->FindBin(i*fSystAvgRange);
      Int_t iStop = reth->FindBin((i+1)*fSystAvgRange);
      if(iStart==0) iStart++;
      if(iStop==nbins) iStop--;
      Int_t bc=0;
      Double_t sum=0;
      for(Int_t j=iStart;j<iStop;j++) {
        if(reth->GetBinContent(j)==0) continue;
        if(reth->GetBinContent(j)>fMaxAllowedSyst) continue;
        sum+=reth->GetBinContent(j);
        bc++;
      }
      if(!bc) continue;
      sum/=bc;
      for(Int_t j=iStart; j<iStop; j++) reth->SetBinContent(j,sum);
    }
  }
  return reth;
};
TH1 *PCCContainer::VarSyst(lFunc sf, Bool_t relative, Bool_t applyBarlow, TObject *inErr) {
  if(sf==kDisabled) return 0;
  TH1 *cvals = (this->*sf)(-1,kFALSE);
  FSRebin(&cvals);
  TH1 *hErr = dynamic_cast<TH1*>(inErr);
  TF1 *fErr = dynamic_cast<TF1*>(inErr);
  if(!inErr) {
    TH1 *systs = getSystematicsSummed(sf,relative,applyBarlow);
    systs->SetName("Systematics_Summed");
    ApplyErrors(cvals,systs,relative);
    delete systs;
  } else if(hErr) ApplyErrors(cvals,hErr,relative);
          else ApplyErrors(cvals,fErr,relative);
  return cvals;
}

void PCCContainer::ApplyBootstrapErrors(TH1 *inh, TProfile *inpf) {
  if(inh->GetNbinsX()!=inpf->GetNbinsX()) {printf("Bootstrap warning: number of bins in histogram and profile is not the same! Not applying...\n"); return; };
  for(Int_t i=1;i<=inh->GetNbinsX();i++) {
    if(inh->GetBinContent(i)==0) continue;
    inh->SetBinError(i,inpf->GetBinError(i));
    if(fBootstrapMean) inh->SetBinContent(i,inpf->GetBinContent(i));
  };
}
void PCCContainer::SetRecGenMatrix(TString infi, TString hname) {
  TFile *tf = new TFile(infi.Data());
  if(!tf) {printf("SetRecGenMatrix: Could not open file %s\n",infi.Data()); return; };
  if(tf->IsZombie()) {printf("SetRecGenMatrix: File %s is a zombie\n",infi.Data()); return; };
  TH2 *htemp = (TH2*)tf->Get(hname.Data());
  if(!htemp) {printf("SetRecGenMatrix: Could not find matrix %s in file %s\n", hname.Data(), infi.Data()); return; };
  SetRecGenMatrix(htemp);
  tf->Close();
}
void PCCContainer::SetNchCorrPar(Double_t offset, Double_t slope) {
  if(!fNchCorrFunc) fNchCorrFunc = new TF1("NchCorrFunction","pol1",0,3000);
  fNchCorrFunc->SetParameters(offset, slope);
}
void PCCContainer::MakeNchMC() {
  if(!fRecGen) return; //Avoid printout here, otherwise it will flood the output
  if(!fck) {printf("MakeNchMC: fCk does not exist; not sure where to pick up binning from\n"); return; };
  //New method: fetch Nch_rec bins straight from the fck output. This way, we don't differentiate between dif. rebinning methods
  //fetching bins:
  TH1 *hTemp = fck->getHist(-1);
  hTemp->SetName("TempHist");
  Int_t Nxbins = hTemp->GetNbinsX();
  Double_t *xbins = new Double_t[Nxbins+1];
  hTemp->GetXaxis()->GetLowEdge(xbins);
  xbins[Nxbins] = hTemp->GetXaxis()->GetBinUpEdge(Nxbins);
  delete hTemp;
  //getting profile with initial binning
  if(fNchMC) delete fNchMC;
  TProfile *tempProf1 = fRecGen->ProfileY("tempPF",1,-1,"s"); //mean and RMS as an error
  TProfile *tppTemp = new TProfile("MeanNch_prof","MeanNch_prof",Nxbins,xbins);
  //fill Nch_gen
  for(Int_t i=1;i<=fMultiDist->GetNbinsX();i++)
    tppTemp->Fill(tempProf1->GetBinCenter(i), tempProf1->GetBinContent(i), fMultiDist->GetBinContent(i));
  fNchMC = new TH1D("Nch_Rec_vs_Gen","Nch_Rec_vs_Gen",Nxbins,xbins);
  fNchMC->SetDirectory(0);
  for(Int_t i=1;i<=fNchMC->GetNbinsX();i++) {
    fNchMC->SetBinContent(i,tppTemp->GetBinContent(i));
    fNchMC->SetBinError(i,(tppTemp->GetBinEffectiveEntries(i)>0)?(tppTemp->GetBinError(i)/tppTemp->GetBinEffectiveEntries(i)):1);
  }
  delete tempProf1;
  delete tppTemp;
  delete [] xbins;
  return;

  //Old method:
  /*
  //rebinning and getting projections
  TProfile *tempProf2 = (TProfile*)tempProf1->Rebin(Nxbins,"tempPF2",xbins);
  fNchMC = (TH1*)tempProf2->ProjectionX("Nch_Rec_vs_Gen");
  //cleaning up
  delete tempProf1;
  delete tempProf2;
  delete [] xbins;
  return;*/
}
TGraphErrors *PCCContainer::NchRecToGen(TH1 *inh, Bool_t includeXError, Double_t xErrorSize) {
  if(!fRecGen && !fNchCorrFunc) {
    printf("NchRecToGen: Neither fRecGen nor fNchCorrFunc exist. Returning uncorrected graph...\n");
    TGraphErrors *retgr = new TGraphErrors();
    for(Int_t i=1;i<=inh->GetNbinsX();i++) {
      Double_t yval = inh->GetBinContent(i);
      Double_t yerr = inh->GetBinError(i);
      Double_t xval = inh->GetBinCenter(i);
      Double_t xerr = includeXError?(inh->GetBinWidth(i)/2.):xErrorSize;
      Int_t grsp = retgr->GetN();
      retgr->SetPoint(grsp,xval,yval);
      retgr->SetPointError(grsp,xerr,yerr);
    }
    return retgr;
  };
  TGraphErrors *retgr = new TGraphErrors();
  if(!fNchCorrFunc) { //If TF1 exists, prefer that
    MakeNchMC();
    if(!fNchMC) return 0;
    for(Int_t i=1;i<=inh->GetNbinsX();i++) {
      if(!fNchMC->GetBinContent(i) || !inh->GetBinContent(i)) continue; //Either Nch Gen ain't there, or the contents are empty
      Double_t yval = inh->GetBinContent(i);
      Double_t yerr = inh->GetBinError(i);
      Double_t xval = fNchMC->GetBinContent(i);
      Double_t xerr = includeXError?fNchMC->GetBinError(i):xErrorSize;
      Int_t grsp = retgr->GetN();
      retgr->SetPoint(grsp,xval,yval);
      retgr->SetPointError(grsp,xerr,yerr);
    };
  } else {
    for(Int_t i=1;i<=inh->GetNbinsX();i++) {
      if(!inh->GetBinContent(i)) continue; //Either Nch Gen ain't there, or the contents are empty
      //Old method -> Use bin center
      // Double_t bcent = inh->GetBinCenter(i);
      // Double_t bcvar = bcent+inh->GetBinWidth(i)/2;
      // Double_t xval = fNchCorrFunc->Eval(bcent);
      // Double_t xerr = includeXError?(fNchCorrFunc->Eval(bcvar)-xval):xErrorSize;
      //New method: use weighted Nch:
      Double_t bcent, berr;
      if(!fMeanNch) { //If MeanNch is not calculated, then we have unit bins and no need to worry
        bcent = inh->GetBinCenter(i);
        berr  = inh->GetBinWidth(i)/2;
      } else {
        // Double_t tx = inh->GetBinCenter(i);
        // Int_t ix = fMeanNch->FindBin(tx);
        //In principle, ix == i, or should be. check this:
        bcent = fMeanNch->GetBinContent(i);
        berr  = (fMeanNch->GetBinEffectiveEntries(i)>0)?(fMeanNch->GetBinError(i)/fMeanNch->GetBinEffectiveEntries(i)):1;
      }
      Double_t xval = fNchCorrFunc->Eval(bcent);
      Double_t xerr = includeXError?(fNchCorrFunc->Eval(bcent+berr)-xval):xErrorSize;
      Double_t yval = inh->GetBinContent(i);
      Double_t yerr = inh->GetBinError(i);
      Int_t grsp = retgr->GetN();
      retgr->SetPoint(grsp,xval,yval);
      retgr->SetPointError(grsp,xerr,yerr);
    };
  }
  return retgr;
}
TH1 *PCCContainer::getMultiHarmonic(Int_t nrb, Bool_t bootstrap) {
  //Calculating numerator
  if(!fcov23) return 0;
  TH1 *hCovV2V3Pt = (this->*kCov23)(nrb,kFALSE);
  TH1 *hCovV2Pt   = (this->*kCov2)(nrb,kFALSE);
  TH1 *hCovV3Pt   = (this->*kCov3)(nrb,kFALSE);
  TH1 *hC22sq     = (this->*kC22)(nrb,kFALSE);
  TH1 *hC32sq     = (this->*kC32)(nrb,kFALSE);
  hCovV2Pt->Multiply(hC32sq);
  hCovV3Pt->Multiply(hC22sq);
  hCovV2V3Pt->Add(hCovV2Pt,-1);
  hCovV2V3Pt->Add(hCovV3Pt,-1);
  delete hC32sq;
  delete hC22sq;
  delete hCovV2Pt;
  delete hCovV3Pt;
  // printf("Calculating denominator...\n");
  //Calculating denominator
  TH1 *sigmaV2    = (this->*kSigma2)(nrb,kFALSE);
  TH1 *sigmaV3    = (this->*kSigma3)(nrb,kFALSE);
  TH1 *ck         = (this->*kCk)(nrb,kFALSE);
  TH1 *hcksqrt    = HistSqrt(ck,kFALSE);
  hcksqrt->Multiply(sigmaV2);
  hcksqrt->Multiply(sigmaV3);
  delete sigmaV2;
  delete sigmaV3;
  delete ck;
  //Ratio and cleaning up
  // printf("Dividing...\n");
  hCovV2V3Pt->Divide(hcksqrt);
  hCovV2V3Pt->SetTitle(Form(";%s; #rho(v_{2}^{2}, v_{3}^{2}, [#it{p}_{T}])",fxTitle.Data()));
  if(!bootstrap) return hCovV2V3Pt;
  return getBootstrapped(hCovV2V3Pt,kMultiHar);
}
TH1 *PCCContainer::getMHTerm1(Int_t nrb, Bool_t bootstrap) {
  //Calculating numerator
  if(!fcov23) return 0;
  TH1 *hCovV2V3Pt = (this->*kCov23)(nrb,kFALSE);
  // printf("Calculating denominator...\n");
  //Calculating denominator
  TH1 *sigmaV2    = (this->*kSigma2)(nrb,kFALSE);
  TH1 *sigmaV3    = (this->*kSigma3)(nrb,kFALSE);
  TH1 *ck         = (this->*kCk)(nrb,kFALSE);
  TH1 *hcksqrt    = HistSqrt(ck,kFALSE);
  hcksqrt->Multiply(sigmaV2);
  hcksqrt->Multiply(sigmaV3);
  delete sigmaV2;
  delete sigmaV3;
  delete ck;
  //Ratio and cleaning up
  // printf("Dividing...\n");
  hCovV2V3Pt->Divide(hcksqrt);
  hCovV2V3Pt->SetTitle(Form(";%s; #rho(v_{2}^{2}, v_{3}^{2}, [#it{p}_{T}])",fxTitle.Data()));
  if(!bootstrap) return hCovV2V3Pt;
  return getBootstrapped(hCovV2V3Pt,kMH1);
}
TH1 *PCCContainer::getMHTerm2(Int_t nrb, Bool_t bootstrap) {
  //Calculating numerator
  TH1 *hCovV2Pt   = (this->*kCov2)(nrb,kFALSE);
  TH1 *hC32sq     = (this->*kC32)(nrb,kFALSE);
  hCovV2Pt->Multiply(hC32sq);
  delete hC32sq;
  // printf("Calculating denominator...\n");
  //Calculating denominator
  TH1 *sigmaV2    = (this->*kSigma2)(nrb,kFALSE);
  TH1 *sigmaV3    = (this->*kSigma3)(nrb,kFALSE);
  TH1 *ck         = (this->*kCk)(nrb,kFALSE);
  TH1 *hcksqrt    = HistSqrt(ck,kFALSE);
  hcksqrt->Multiply(sigmaV2);
  hcksqrt->Multiply(sigmaV3);
  delete sigmaV2;
  delete sigmaV3;
  delete ck;
  //Ratio and cleaning up
  // printf("Dividing...\n");
  hCovV2Pt->Divide(hcksqrt);
  hCovV2Pt->SetTitle(Form(";%s; #rho(v_{2}^{2}, v_{3}^{2}, [#it{p}_{T}])",fxTitle.Data()));
  if(!bootstrap) return hCovV2Pt;
  return getBootstrapped(hCovV2Pt,kMH2);
}
TH1 *PCCContainer::getMHTerm3(Int_t nrb, Bool_t bootstrap) {
  //Calculating numerator
  TH1 *hCovV3Pt   = (this->*kCov3)(nrb,kFALSE);
  TH1 *hC22sq     = (this->*kC22)(nrb,kFALSE);
  hCovV3Pt->Multiply(hC22sq);
  delete hC22sq;
  // printf("Calculating denominator...\n");
  //Calculating denominator
  TH1 *sigmaV2    = (this->*kSigma2)(nrb,kFALSE);
  TH1 *sigmaV3    = (this->*kSigma3)(nrb,kFALSE);
  TH1 *ck         = (this->*kCk)(nrb,kFALSE);
  TH1 *hcksqrt    = HistSqrt(ck,kFALSE);
  hcksqrt->Multiply(sigmaV2);
  hcksqrt->Multiply(sigmaV3);
  delete sigmaV2;
  delete sigmaV3;
  delete ck;
  //Ratio and cleaning up
  // printf("Dividing...\n");
  hCovV3Pt->Divide(hcksqrt);
  hCovV3Pt->SetTitle(Form(";%s; #rho(v_{2}^{2}, v_{3}^{2}, [#it{p}_{T}])",fxTitle.Data()));
  if(!bootstrap) return hCovV3Pt;
  return getBootstrapped(hCovV3Pt,kMH3);
}
TH1 *PCCContainer::getMHTerm23(Int_t nrb, Bool_t bootstrap) {
  //Calculating numerator
  TH1 *hCovV2Pt   = (this->*kCov2)(nrb,kFALSE);
  TH1 *hC32sq     = (this->*kC32)(nrb,kFALSE);
  TH1 *hCovV3Pt   = (this->*kCov3)(nrb,kFALSE);
  TH1 *hC22sq     = (this->*kC22)(nrb,kFALSE);
  hCovV3Pt->Multiply(hC22sq);
  hCovV2Pt->Multiply(hC32sq);
  hCovV2Pt->Add(hCovV3Pt);
  delete hCovV3Pt;
  delete hC22sq;
  delete hC32sq;
  // printf("Calculating denominator...\n");
  //Calculating denominator
  TH1 *sigmaV2    = (this->*kSigma2)(nrb,kFALSE);
  TH1 *sigmaV3    = (this->*kSigma3)(nrb,kFALSE);
  TH1 *ck         = (this->*kCk)(nrb,kFALSE);
  TH1 *hcksqrt    = HistSqrt(ck,kFALSE);
  hcksqrt->Multiply(sigmaV2);
  hcksqrt->Multiply(sigmaV3);
  delete sigmaV2;
  delete sigmaV3;
  delete ck;
  //Ratio and cleaning up
  // printf("Dividing...\n");
  hCovV2Pt->Divide(hcksqrt);
  hCovV2Pt->SetTitle(Form(";%s; #rho(v_{2}^{2}, v_{3}^{2}, [#it{p}_{T}])",fxTitle.Data()));
  if(!bootstrap) return hCovV2Pt;
  return getBootstrapped(hCovV2Pt,kMH23);
}


TH1 *PCCContainer::GetTerms() {
  // fFC->SetIDName("ChSC");
  // TProfile *tpf1 = fFC->GetCorrXXVsMulti("234");
  fFC->SetIDName("ChGap");
  TProfile *tpfv3 = fFC->GetCorrXXVsMulti("32");
  TProfile *tpfv2 = fFC->GetCorrXXVsMulti("22");

  TH1 *h1 = fcov23->getHist(-1);//ProfToHist(tpf1);
  TH1 *h2 = ProfToHist(tpfv3);
  TH1 *h3 = ProfToHist(tpfv2);

  h3->Multiply(fcov3->getHist(-1));
  h2->Multiply(fcov2->getHist(-1));

  h1->Add(h2,-1);
  h1->Add(h3,-1);
  return h1;
}
TString PCCContainer::GetSFDescription(lFunc sf) {
  if(sf==kCk)  return "kCk";
  if(sf==kSigma2) return "kSigma2";
  if(sf==kSigma3Sub) return "kSigma3Sub";
  if(sf==kCov2) return "kCov2";
  if(sf==kPCC2) return "kPCC2";
  if(sf==kPCC3Sub) return "kPCC3Sub";
  return "Undefined";
}
void PCCContainer::ApplyClosureCorrection(lFunc sf, TGraphErrors *ingr) {
  if(fCorrFile.IsNull()) { printf("PCCContainer::ApplyClosureCorrection: correction file not defined!\n"); return; };
  TFile *tempf = new TFile(fCorrFile.Data(),"READ");
  if(!tempf) { printf("PCCContainer::ApplyClosureCorrection: file %s not found!\n", fCorrFile.Data()); return; };
  if(tempf->IsZombie()) {  printf("PCCContainer::ApplyClosureCorrection: file %s is a zombie!\n", fCorrFile.Data()); return;  };
  TString funcName = GetSFDescription(sf);
  TF1 *tfunc = (TF1*)tempf->Get(funcName.Data());
  if(!tfunc) { printf("PCCContainer::ApplyClosureCorrection: function %s not found in file %s. Not applying...\n",funcName.Data(),fCorrFile.Data()); tempf->Close(); return; };
  Double_t *xv = ingr->GetX();
  Double_t *yv = ingr->GetY();
  Double_t *ye = ingr->GetEY();
  for(Int_t i=0;i<ingr->GetN();i++) {
    Double_t sf = tfunc->Eval(xv[i]);
    if(sf==0) { printf("PCCContainer::ApplyClosureCorrection: could not evaluate function %s at point %f!\n",funcName.Data(),xv[i]); continue; };
    yv[i] = yv[i]/sf;
    ye[i] = ye[i]/sf;
  };
  delete tfunc;
  tempf->Close();
}
void PCCContainer::f_ResetBin(Int_t nbin) {
    if(fcov2) fcov2->ResetBin(nbin);
    if(fcov3) fcov3->ResetBin(nbin);
    if(fcov23) fcov23->ResetBin(nbin);
    if(fcov2nopt) fcov2nopt->ResetBin(nbin);
    if(fcov3nopt) fcov3nopt->ResetBin(nbin);
    if(fcov23nopt) fcov23nopt->ResetBin(nbin);
}
void PCCContainer::ResetBin(Int_t nbin, Int_t lSyst) {
  if(!lSyst || lSyst<0) f_ResetBin(nbin);
  if(lSyst>0) if(fSysts->FindObject(Form("PCCCont_Syst%i",lSyst))) ((PCCContainer*)fSysts->FindObject(Form("PCCCont_Syst%i",lSyst)))->ResetBin(nbin,0);
}
void PCCContainer::ResetBin(Int_t nbin, vector<Int_t> lSInd) {
  for(auto i:lSInd) ResetBin(nbin,i);
}
void PCCContainer::OverrideSystematics(Int_t keyInd, Double_t newval) {
  if(fSystOverride.count(keyInd)) fSystOverride.at(keyInd) = newval;
  else fSystOverride.insert(make_pair(keyInd,newval));
}
void PCCContainer::FSRebin(TH1 **inh) {
  if(!fNrb) return; //Rebin only when needed
  if(!fUseFSRebin) return;
  TH1 *sH = *inh;
  TString hName(sH->GetName());
  sH->SetName(Form("%s_BU",hName.Data()));
  TH1 *tH = sH->Rebin(fNrb,hName.Data(),fBrb);
  tH->Reset();
  for(Int_t i=1;i<=sH->GetNbinsX();i++) {
    Double_t bc = sH->GetBinContent(i);
    Double_t be = sH->GetBinError(i);
    if(be==0) continue;
    be = (1./(be*be));
    Double_t bind = tH->FindBin(sH->GetBinCenter(i));
    tH->SetBinContent(bind,tH->GetBinContent(bind)+bc*be);
    tH->SetBinError(bind,tH->GetBinError(bind)+be);
  }
  for(Int_t i=1;i<=tH->GetNbinsX();i++) {
    Double_t bc = tH->GetBinContent(i);
    Double_t be = tH->GetBinError(i);
    if(be==0) continue;
    tH->SetBinContent(i,bc/be);
    tH->SetBinError(i,TMath::Sqrt(1./be));
  }
  delete sH;
  (*inh) = tH;
};
TGraphErrors *PCCContainer::f_HtoGr(TH1 *inh, Double_t xOffset, Double_t xError) {
  if(!inh) {printf("Input histogram is a NULL!\n"); return 0; };
  TGraphErrors *retgr = new TGraphErrors();
  for(Int_t i=1;i<=inh->GetNbinsX();i++) {
    Double_t v_x = inh->GetBinCenter(i)+xOffset;
    Double_t v_y = inh->GetBinContent(i);
    Double_t e_x = inh->GetBinWidth(i)/2;
    if(xError>=0) e_x=xError;
    Double_t e_y = inh->GetBinError(i);
    retgr->SetPoint(i-1,v_x,v_y);
    retgr->SetPointError(i-1,e_x,e_y);
  }
  retgr->GetXaxis()->SetTitle(inh->GetXaxis()->GetTitle());
  retgr->GetYaxis()->SetTitle(inh->GetYaxis()->GetTitle());
  return retgr;
}
