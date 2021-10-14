#ifndef PCCContainer__h
#define PCCContainer__h

#include "TFile.h"
#include "AliGFWFlowContainer.h"
#include "AliCkContainer.h"
#include "TString.h"
#include "TList.h"
#include "TProfile.h"
#include "TH1D.h"
#include "TNamed.h"
#include "TRegexp.h"
#include "AliProfileBS.h"
#include "TObjArray.h"
#include "TProfile2D.h"
#include "SystConfig.C"
#include "TH2D.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include <vector>
class PCCContainer:public TNamed {
public:
  typedef TH1* (PCCContainer::*lFunc)(Int_t, Bool_t);
  PCCContainer();
  PCCContainer(TString, Int_t MaxSyst=0);
  ~PCCContainer();
  static TString GetSFDescription(lFunc sf);
  Bool_t ConstructContainer(TString infi);
  Bool_t ConstructContainer(TFile *infi, Int_t lInd);
  Bool_t OverrideMeanPt(TString infi, Int_t lInd); //Temporary method to read mean pt
  // void RebinMulti(Int_t nrb);
  void RebinMulti(Int_t nrb, Double_t *chbins);
  TH1 *Var(lFunc sf, Int_t ind=0) { return getSystematicsObs(sf,ind,kTRUE); };//(this->*sf)(ind,kTRUE); };
  TH1 *VarSyst(lFunc, Bool_t relative=kTRUE, Bool_t applyBarlow=kTRUE, TObject *inErr=0);
  //Simple graph with statistical uncertainties
  TGraphErrors *Gr(lFunc sf, Double_t xOffset=0, Double_t xError=-1) { TH1 *ht = Var(sf); TGraphErrors *gr = f_HtoGr(ht,xOffset,xError); delete ht; return gr; };
  //Simple graph with syst. uncertainties
  TGraphErrors *GrSyst(lFunc sf, Bool_t relative=kTRUE, Bool_t applyBarlow=kTRUE, TObject *inErr=0, Double_t xOffset=0, Double_t xError=-1) { TH1 *ht = VarSyst(sf,relative,applyBarlow,inErr); TGraphErrors *gr = f_HtoGr(ht,xOffset,xError); delete ht; return gr; };
  //Graph with stat + Nch unfolding and closure correction
  TGraphErrors *GrVar(lFunc sf, Bool_t includeXError=kFALSE, Double_t xErrorSize=0.5) { TH1 *htemp = Var(sf,0); TGraphErrors *retgr = NchRecToGen(htemp,includeXError,xErrorSize); delete htemp; ApplyClosureCorrection(sf, retgr); return retgr; };
  //Graph with syst + Nch unfolding and closure correction
  TGraphErrors *GrVarSyst(lFunc sf, Bool_t relative=kTRUE, Bool_t applyBarlow=kTRUE, Bool_t includeXError=kFALSE, Double_t xErrorSize=0.5) { TH1 *htemp = VarSyst(sf,relative,applyBarlow); TGraphErrors *retgr = NchRecToGen(htemp,includeXError,xErrorSize); ApplyClosureCorrection(sf, retgr); delete htemp; return retgr; }
  TH1 *Syst(lFunc sf, Int_t ind, Bool_t rel=kTRUE, Bool_t applyBarlow=kTRUE) { return getSystematics(sf,ind,rel,kFALSE,applyBarlow); };
  TH1 *Barlow(lFunc sf, Int_t ind) { return  getBarlowTest(sf, ind, getCorrelationError(ind)); };
  Bool_t isRelevant(Int_t ind) { return getCorrelationError(ind); };
  void SetRecGenMatrix(TString, TString hname="RecoGen");
  void SetNchCorrPar(Double_t offset, Double_t slope);
  void SetCorrectionFile(TString infile) { fCorrFile = infile; };
  void ApplyClosureCorrection(lFunc, TGraphErrors *ingr);
  void OverrideSystematics(Int_t keyInd, Double_t newval);
  void ClearSystematicOverride() { fSystOverride.clear(); };
  Bool_t isInitialized() { return fInitialized; };
  void SetBootstrapMean(Bool_t newval) { fBootstrapMean = newval; };
  void SetUseFSRebin(Bool_t newval) { fUseFSRebin=newval; if(newval && fBrb && fNrb>0) RebinMulti(fNrb,fBrb); };
  void ResetBin(Int_t nbin, Int_t lSyst=-1);
  void ResetBin(Int_t nbin, vector<Int_t> lSInd);
  void PresetWeights(AliProfileBS *bs);
  Bool_t fInitialized;
  TList *fObjs;
  TList *fSysts;
  Bool_t fBootstrapMean;
  Bool_t fUseFSRebin;
  Int_t fNrb;
  Double_t *fBrb;
  Int_t fSystAvgRange;
  Double_t fMaxAllowedSyst;
  Int_t fMaxSyst;
  TString fSystName;
  TH1 *fMultiDist;
  TString fCorrFile;
  TProfile *fMeanNch;
  AliGFWFlowContainer *fFC;
  AliCkContainer *fck;
  AliProfileBS *fckbs;
  AliProfileBS *fcov2;
  AliProfileBS *fcov3;
  AliProfileBS *fcov23;
  AliProfileBS *fcov2nopt;
  AliProfileBS *fcov3nopt;
  AliProfileBS *fcov23nopt;
  AliProfileBS *fmpt;
  map<int,vector<int>> fIndMap;
  map<int,double> fSystOverride;
  Int_t fNRebinForSyst;
  TH2 *fRecGen;
  TH1 *fNchMC;
  TF1 *fNchCorrFunc;
  TList *f_mpt_tl;
  TList *f_cov_tl;
  TString fxTitle;
  void ApplyBootstrapErrors(TH1 *inh, TProfile *inpf);
  Bool_t BuildIndexMap(Bool_t force=kFALSE);
  void SetRebinForSyst(Int_t newval) { fNRebinForSyst=newval; };
  void SetSystAverageRange(Int_t lrange, Double_t lMaxAllowed=1) {fSystAvgRange=lrange; fMaxAllowedSyst=lMaxAllowed; };
  TH1 *CalculateCovariance(const Int_t &nrb, AliProfileBS* l_covConst, AliProfileBS* l_covLin, AliProfileBS* l_mpt, AliProfileBS *rbWeights=0);
  vector<TH1*> nullvec;
  //basic calls:
  TH1 *getCov2(Int_t nrb=-1, Bool_t bootstrap=kTRUE);
  TH1 *getCov3(Int_t nrb=-1, Bool_t bootstrap=kTRUE);
  TH1 *getCov23(Int_t nrb=-1, Bool_t bootstrap=kTRUE);
  TH1 *getCk(Int_t nrb=-1, Bool_t bootstrap=kTRUE);
  TH1 *getVar2(Int_t nrb=-1, Bool_t bootstrap=kTRUE);
  TH1 *getVar3(Int_t nrb=-1, Bool_t bootstrap=kTRUE);
  TH1 *getSigma2(Int_t nrb=-1, Bool_t bootstrap=kTRUE);
  TH1 *getSigma3(Int_t nrb=-1, Bool_t bootstrap=kTRUE);
  TH1 *getPCC2(Int_t nrb=-1, Bool_t bootstrap=kTRUE);
  TH1 *getPCC3(Int_t nrb=-1, Bool_t bootstrap=kTRUE);
  TH1 *getMpt(Int_t nrb=-1, Bool_t bootstrap=kFALSE) { return fmpt->ProjectionX("Mpt"); }; //No bootstrapping, etc., but need arguments for compatibility
  TH1 *getSqCkOverMpt(Int_t nrb=-1, Bool_t bootstrap=kTRUE);
  TH1 *getVar3Sub(Int_t nrb=-1, Bool_t bootstrap=kTRUE);
  TH1 *getSigma3Sub(Int_t nrb=-1, Bool_t bootstrap=kTRUE);
  TH1 *getPCC3Sub(Int_t nrb=-1, Bool_t bootstrap=kTRUE);
  TH1 *getMultiHarmonic(Int_t nrb=-1, Bool_t bootstrap=kTRUE);

  TH1 *getC22(Int_t nrb, Bool_t bootstrap);
  TH1 *getC24(Int_t nrb=-1, Bool_t bootstrap=kTRUE);
  TH1 *getC32(Int_t nrb, Bool_t bootstrap);
  TH1 *getC34(Int_t nrb=-1, Bool_t bootstrap=kTRUE);
  TH1 *getC243Sub(Int_t nrb=-1, Bool_t bootstrap=kTRUE);
  TH1 *getC243Sub(AliGFWFlowContainer *infc, Int_t nrb, Bool_t bootstrap, TString ms, TString conf1, TString conf2);

  //Debugging
  TH1 *getMHTerm1(Int_t nrb=-1, Bool_t bootstrap=kTRUE);
  TH1 *getMHTerm2(Int_t nrb=-1, Bool_t bootstrap=kTRUE);
  TH1 *getMHTerm3(Int_t nrb=-1, Bool_t bootstrap=kTRUE);
  TH1 *getMHTerm23(Int_t nrb=-1, Bool_t bootstrap=kTRUE);

// private:
  TFile *getFile(TString fina);
  TObject *getObj(TFile *infi, TString objName);
  Ssiz_t findStrInd(TString, TString);
  TH1 *ProfToHist(TProfile *inpf) { if(!inpf) {printf("ProfToHist: null profile provided!\n"); return 0; }; return (TH1*)inpf->ProjectionX(Form("%s_hist",inpf->GetName())); };
  void HistOper(TH1 *reth, Double_t cvl=0, Double_t lpower=1, Bool_t Abs=kFALSE, Bool_t RemoveErrors=kFALSE);
  void AddInQuad(TH1 *h1, TH1 *h2);
  void SetRecGenMatrix(TH2 *inh) { if(fRecGen) delete fRecGen; if(!inh) { printf("Unsetting NchRecGen matrix!\n"); return; }fRecGen = (TH2*)inh->Clone("RecGenMatrix"); fRecGen->SetDirectory(0); };
  void MakeNchMC();
  void SetxTitle(TString newval) { fxTitle = newval; };
  TGraphErrors *NchRecToGen(TH1 *inh, Bool_t includeXError=kFALSE, Double_t xErrorSize=0.5); //if includeXError is true, MC Nch spread as x-error; otherwise, use xErrorSize

  //Kind-of template for bootstrapping:
  TH1 *getBootstrapped(TH1 *centVal, lFunc sf) {
    TProfile *statProf = getStatistics(sf);
    ApplyBootstrapErrors(centVal,statProf);
    delete statProf;
    return centVal;
  };
  //Envelope classes. I think this could be redone in templates to make things more brief, since the callers are somewhat repetetive
  TH1 *getCN4(Int_t nrb, Int_t nHarmonic);
  TH1 *getCN2(Int_t nrb, Int_t nHarmonic);
  TH1 *getVarNHar(Int_t nrb, lFunc f_cn2, lFunc f_cn4);
  TH1 *getSigmaNHar(Int_t nrb, Bool_t bootstrap, lFunc varFunc);
  TH1 *getPCCNHar(Int_t nrb, lFunc covFunc, lFunc sigmaFunc);

  // TH1 *getPCCNhar(Int_t nrb=-1, Bool_t bootstrap=kTRUE);


  void ApplyBarlow(TH1 *target, TH1 *barlow);
  void ApplyErrors(TH1 *inh, TH1 *errors, Bool_t relative=kTRUE);
  void ApplyErrors(TH1 *inh, TF1 *errors, Bool_t relative=kTRUE);
  void RemoveErrors(TH1 *inh) {for(Int_t i=1;i<=inh->GetNbinsX();i++) if(inh->GetBinError(i)) inh->SetBinError(i,0); };
  TH1 *HistSqrt(TH1 *inh, Bool_t removeErrors=kFALSE);
  TProfile *getStatistics(lFunc);
  TH2 *getStatisticsDist(lFunc,Int_t nbins=100, Double_t ymin=-1, Double_t ymax=1);
  TProfile *getStatistics(lFunc, Int_t);
  TH1 *getSystematicsObs(lFunc, Int_t lSyst, Bool_t bootstrap=kFALSE);
  TH1 *getBarlowTest(lFunc, Int_t lSyst, Double_t corrFactor=1);
  TH1 *getSystematics(lFunc, Int_t lSyst, Bool_t relative=kTRUE, Bool_t bootstrap=kFALSE, Bool_t applyBarlow=kFALSE, Bool_t takeAbs=kTRUE);
  TH1* getSystematicsSummed(lFunc, Bool_t relative=kTRUE, Bool_t inclErrors=kFALSE, Bool_t applyBarlow=kFALSE);
  vector<TH1*> getMergedErrors(lFunc, Bool_t relative=kTRUE, Bool_t inclErrors=kFALSE, Bool_t applyBarlow=kFALSE, Bool_t weightByErrors=kFALSE);
  vector<TH1*> getSystSubset(lFunc, Int_t KeyInd, Bool_t relative=kTRUE, Bool_t inclErrors=kFALSE, Bool_t ApplyBarlow=kFALSE, Bool_t takeAbs=kTRUE);
  TH1 *GetTerms();
  void FSRebin(TH1 **inh);
private:
  TGraphErrors *f_HtoGr(TH1 *inh, Double_t xOffset=0, Double_t xError=-1);
  void f_ResetBin(Int_t nbin);

  ClassDef(PCCContainer, 1);
};
namespace PCCSpace {
  PCCContainer::lFunc kCov2 = &PCCContainer::getCov2; //Cov(v2,pt)
  PCCContainer::lFunc kCov3 = &PCCContainer::getCov3; //Cov(v3,pt)
  PCCContainer::lFunc kCov23 = &PCCContainer::getCov23; //Cov(v2,v3,pt)
  PCCContainer::lFunc kCk = &PCCContainer::getCk; //cK
  PCCContainer::lFunc kVar2 = &PCCContainer::getVar2; //Var(v_2^2)
  PCCContainer::lFunc kVar3 = &PCCContainer::getVar3; //Var(v_3^2)
  PCCContainer::lFunc kSigma2 = &PCCContainer::getSigma2; //sigma(v_2^2)
  PCCContainer::lFunc kSigma3 = &PCCContainer::getSigma3; //sigma(v_3^2)
  PCCContainer::lFunc kPCC2 = &PCCContainer::getPCC2; //PCC(v2,pt)
  PCCContainer::lFunc kPCC3 = &PCCContainer::getPCC3; //PCC(v3,pt)
  PCCContainer::lFunc kMpt = &PCCContainer::getMpt; //[pt]
  PCCContainer::lFunc kSQCK = &PCCContainer::getSqCkOverMpt; //ck/sqrt{pt}
  PCCContainer::lFunc kVar3Sub = &PCCContainer::getVar3Sub; //Var(v_2^2,3-sub)
  PCCContainer::lFunc kSigma3Sub = &PCCContainer::getSigma3Sub; //Sigma(v_2^2,3-sub)
  PCCContainer::lFunc kPCC3Sub = &PCCContainer::getPCC3Sub; //PCC(v2, pt), 3-sub
  PCCContainer::lFunc kMultiHar= &PCCContainer::getMultiHarmonic; //PCC(v_2^2,v_3^2,pt)
  PCCContainer::lFunc kC22 = &PCCContainer::getC22; //c_2{2}
  PCCContainer::lFunc kC24 = &PCCContainer::getC24; //c_2{4}
  PCCContainer::lFunc kC32 = &PCCContainer::getC32; //c_3{2}
  PCCContainer::lFunc kC34 = &PCCContainer::getC34; //c_3{4}
  PCCContainer::lFunc kC243Sub = &PCCContainer::getC243Sub; //SC{2,4}, 3-sub
  PCCContainer::lFunc kSame = 0; //NA
  PCCContainer::lFunc kDisabled = 0; //NA
  PCCContainer::lFunc kMH1 = &PCCContainer::getMHTerm1; //SC{2,4}, 3-sub
  PCCContainer::lFunc kMH2 = &PCCContainer::getMHTerm2; //SC{2,4}, 3-sub
  PCCContainer::lFunc kMH3 = &PCCContainer::getMHTerm3; //SC{2,4}, 3-sub
  PCCContainer::lFunc kMH23 = &PCCContainer::getMHTerm23; //SC{2,4}, 3-sub

}
#endif
