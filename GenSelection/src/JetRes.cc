#include "../include/JetRes.hh"
#include "TFile.h"
#include "TList.h"
#include "TRandom.h"

JetRes::JetRes(bool iData) { 
  std::string lRes = "/afs/cern.ch/user/p/pharris/pharris/public/bacon/prod/CMSSW_7_4_7/src/BaconAnalyzer/MJSelection/Json/fitgraph_v4.root";
  lRes = "/afs/cern.ch/user/p/pharris/pharris/public/bacon/prod/CMSSW_7_4_7/src/BaconAnalyzer/MJSelection/Json/";
  if(iData)  lRes += "fitgraph_v6_vtx.root";
  if(!iData) lRes += "fitgraph_v4.root";
  fData = iData;
  TFile *lFile = new TFile(lRes.c_str()); 
  TList *lList = (TList*) lFile->GetListOfKeys(); 
  for(int i0 = 0; i0 < lList->GetEntries(); i0++) { 
    TObject* pKey  = (TObject*) lList->At(i0); 
    TGraph * pFunc = (TGraph*) lFile->FindObjectAny(pKey->GetName()); 
    std::string pName = pFunc->GetName(); 
    double *lId = getId(pName+"_");
    TRes* pEtaFunc = new TRes(pFunc,lId[2],lId[3]);
    if(lId[0] == 0) { 
      if(lId[1] == 0) fBase .push_back(pEtaFunc);
      if(lId[1] == 1) fHBase.push_back(pEtaFunc);
      if(lId[1] == 2) fRec  .push_back(pEtaFunc);
    }
    if(lId[0] == 1) { 
      if(lId[1] == 0) fS1Base .push_back(pEtaFunc);
      if(lId[1] == 1) fS1HBase.push_back(pEtaFunc);
      if(lId[1] == 2) fS1Rec  .push_back(pEtaFunc);
    }
    if(lId[0] == 2) { 
      if(lId[1] == 0) fS2Base .push_back(pEtaFunc);
      if(lId[1] == 1) fS2HBase.push_back(pEtaFunc);
      if(lId[1] == 2) fS2Rec  .push_back(pEtaFunc);
    }
  }
  fRandom = new TRandom(0xDEADBEEF);
  std::string lName = "/afs/cern.ch/user/p/pharris/pharris/public/bacon/prod/CMSSW_7_4_7/src/BaconAnalyzer/MJSelection/prod/files/";
  if( iData) lName += "recoilfit_zmmData_pf_v2.root";
  if(!iData) lName += "recoilfit_zmmMC_pf_v2.root";
  fCorrector = new RecoilCorrector(lName);
}
std::vector<TLorentzVector> JetRes::smearEvent(double iGMet,double iGMetPhi,double iHT,std::vector<baconhep::TGenJet*> &iJets) {
  std::vector<TLorentzVector> lSmeared;
  TLorentzVector lMet; lMet.SetPtEtaPhiM(iGMet,0,iGMetPhi,0);
  double lHT = iHT;
  for(unsigned int i0 = 0; i0 < iJets.size(); i0++) {
    TLorentzVector pVec = smearedJet(iJets[i0],0.,0.);
    lSmeared.push_back(pVec);
    lMet += pVec;
    lHT += -iJets[i0]->pt+pVec.Pt();
  }
  double lDMet    = lMet.Pt();
  double lDMetPhi = lMet.Phi();
  double lU1 =0,lU2=0;
  int lNCand = 0.;
  //std::cout << "===> PreMET : " << lDMet << " - " << lDMetPhi << " ";
  fCorrector->CorrectAll(lDMet,lDMetPhi,lMet.Pt(),lMet.Phi(),0,0,lU1,lU2,0,0,lNCand);
  //std::cout << "===> PostMET : " << lDMet << " - " << lDMetPhi << " Base : " << iGMet << " - " << iGMet;
  TLorentzVector lVec; lVec.SetPtEtaPhiM(lDMet,0,lDMetPhi,0); //lVec.RotateZ(TMath::Pi());                                                                                                       
  for(unsigned int i0 = 0; i0 < lSmeared.size(); i0++) lVec -= lSmeared[i0]; 
  //std::cout << " ---> " << lVec.Pt() << " -- " << lVec.Phi() << std::endl;
  lSmeared.push_back(lVec);
  TLorentzVector lVHT; lVHT.SetPtEtaPhiM(lHT,0,0,0); //lVec.RotateZ(TMath::Pi());                                                                                                       
  lSmeared.push_back(lVHT);
  return lSmeared;
}
double* JetRes::getId(std::string iName) { 
  std::string delimiter = "_";
  double *lId = new double[4];
  std::string acc = "";
  int lCount = 0; 
  for(unsigned int i0 = 0; i0 < iName.size(); i0++) {
     if(iName[i0] == delimiter) {
       if(lCount == 1 && acc == "sigm") lId[0] = 0;
       if(lCount == 1 && acc == "sig1") lId[0] = 1;
       if(lCount == 1 && acc == "sig2") lId[0] = 2;
       if(lCount > 1) { 
	 double pVal = atof(acc.c_str()); 
	 lId[lCount-1] = pVal;
       }
       acc = "";
       lCount++; 
       continue;
     }
     acc += iName[i0];
  }
  return lId;
}
double JetRes::res(double iCSV,double iPt,double iEta,double iSig) { 
  double lCorr = 1. + iSig*0.1;
  if(fData) lCorr *= 0.95;
  TRes *pRes = 0; 
  double lEta = iEta; if(iCSV > 0.5 && fabs(iEta) > 2.5) lEta = 2.45*lEta/fabs(lEta); 
  if(fabs(lEta) > 5) lEta = 4.9*lEta/fabs(lEta);
  if(iCSV  > 0.5 && iCSV < 100000)  pRes = find(fHBase,lEta);
  if(iCSV <= 0.5 && iCSV > -15)     pRes = find(fBase ,lEta);
  if(iCSV < -15 && iCSV > -1000)   pRes = find(fRec  ,iEta);
  return TMath::Max(pRes->func->Eval(iPt)*lCorr,0.1);
}
double JetRes::smear(double iCSV,double iPt,double iEta,double iSig) { 
  double lCorr = 1. + iSig*0.1;
  if(fData) lCorr *= 0.95;
  double lSigma = 0.1;
  double lEta = iEta; if(iCSV > 0.5 && fabs(iEta) > 2.5) lEta = 2.45*lEta/fabs(lEta); 
  if(fabs(lEta) > 5) lEta = 4.9*lEta/fabs(lEta);
  if(iCSV > 0.5 && iCSV < 100000)  lSigma = findSig(fHBase,fS1HBase,fS2HBase,iPt,lEta);
  if(iCSV < 0.5 && iCSV > -15)     lSigma = findSig(fBase ,fS1Base ,fS2Base ,iPt,lEta);
  if(iCSV < -15 && iCSV > -1000)   lSigma = findSig(fRec  ,fS1Rec  ,fS2Rec  ,iPt,iEta);
  return fRandom->Gaus(iPt,lSigma*lCorr);
}
double JetRes::findSig(std::vector<TRes*> &iBase,std::vector<TRes*> &iSig1,std::vector<TRes*> &iSig2,double iPt,double iEta) { 
  TRes* lSigM = find(iBase,iEta);
  TRes* lSig1 = find(iSig1,iEta);
  TRes* lSig2 = find(iSig2,iEta);
  double lM  = TMath::Max(lSigM->func->Eval(iPt),0.1);
  double lS1 = TMath::Max(lSig1->func->Eval(iPt),0.05);
  double lS2 = TMath::Max(lSig2->func->Eval(iPt),lM);
  //lM = f*lS1 + (1-f)*lS2 => lM -lS2 = f*(lS1-lS2) => 
  double lF = (lM - lS2)/(lS1-lS2);
  //lF = 1-(1-lF)/2.;
  double lRand = fRandom->Uniform(0,1);
  if(lRand < lF) return lS1;
  return lS2;
}
JetRes::TRes* JetRes::find(std::vector<TRes*> &iBase,double iEta) { 
  for(unsigned int i0 = 0; i0 < iBase.size(); i0++) { 
    if(iEta < iBase[i0]->etaMin) continue;
    if(iEta > iBase[i0]->etaMax) continue;
    return iBase[i0];
  }
  std::cout << " ---> Should not be here " << iEta << " -- " << iBase.size() << std::endl;
  return iBase[0];
}
TLorentzVector JetRes::smearedJet(baconhep::TGenJet *iJet,double iCSV,double iSig) { 
  double lOldJPt   = iJet->pt;
  double lJPt      = smear(iCSV,iJet->pt,iJet->eta,iSig);//fRandom->Gaus(iVec.Pt() , lSigma);
  double lJEta     = iJet->eta;//fRandom->Gaus(iVec.Eta(),lSigmaEta);
  double lJPhi     = iJet->phi;//fRandom->Gaus(iVec.Phi(),lSigmaPhi);
  double lJM       = (iJet->pt/lOldJPt) * iJet->mass;
  TLorentzVector lVec;
  lVec.SetPtEtaPhiM(lJPt,lJEta,lJPhi,lJM);
  return lVec;
}
