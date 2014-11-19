#include "../include/LeptonsWW.hh"

LeptonsWW::LeptonsWW(TTree *iTree) { 
  setupSelector();
  fEffWeight = 0; 
  fRandom    = new TRandom2(0xDEADBEEF);
  iTree->Branch("effweight",&fEffWeight,"fEffWeight/F");
}
LeptonsWW::~LeptonsWW() {
  delete fFakeRate;
}
void LeptonsWW::setupSelector() { 
  fFakeRate = new TF1("Fakes","TMath::Landau(x,[0],[1])+[2]");
  fFakeRate->SetParameter(0,-3.22);
  fFakeRate->SetParameter(1,3.6);
  fFakeRate->SetParameter(2,0.0026);
}
float LeptonsWW::smearLep(int iId, float &iPt) { 
  //Using common sense numbers for now
  float lPt = iPt;
  if(fabs(iId) == 11) lPt = fRandom->Gaus(iPt,2.); // 2GeV for muons
  if(fabs(iId) == 13) lPt = fRandom->Gaus(iPt,3.); // 3GeV for electrons
  if(fabs(iId) == 17) lPt = fRandom->Gaus(iPt,2.); // 2GeV for muons
  if(fabs(iId) == 18) lPt = fRandom->Gaus(iPt,3.); // 3GeV for electrons
  if(fabs(iId) == 19) lPt = fRandom->Gaus(iPt,4.); // 4GeV for taus 
  return lPt;
}
double LeptonsWW::effWeight(double iPt,double iEff,double iSigma,double iMed) { 
  double lX   = (iPt-iMed)/iSigma;
  double lVal = iEff*(TMath::Erf(lX) + 1.)/2.;
  return lVal;
}
void LeptonsWW::fillWeight(int iId,double iPt) { 
  //Compute Lepton Smearing + Eff Weights
  fEffWeight = 1.;
  if(fabs(iId) == 11 ) fEffWeight *= effWeight(iPt,0.8,15);  //=> Flat 90% eff for muons (we can do better)  => at 10 GeV this drops to  50% 
  if(fabs(iId) == 13 ) fEffWeight *= effWeight(iPt,0.9,8);   //=> Flat 90% eff for muons (we can do better)  => at 10 GeV this drops to  50% 
  if(fabs(iId) == 17 ) fEffWeight *= effWeight(iPt,0.9,8);   //=> Flat 90% eff for muons (we can do better)  => at 10 GeV this drops to  50% 
  if(fabs(iId) == 18 ) fEffWeight *= effWeight(iPt,0.8,15);  //=> Flat 80& eff for electrons                 => at 10 GeV this drops to 30% 
  if(fabs(iId) == 19 ) fEffWeight *= effWeight(iPt,0.6);     //=> Flat 60& eff for electrons                 => at 10 GeV this drops to 30% 
  if(fabs(iId) == 20 ) fEffWeight *= fFakeRate->Eval(iPt);
  if(fabs(iId) == 21 ) fEffWeight *= 0.01 * 2.;//See discussion below + x2 for degeneracy (either can fake)
  //if(fabs(fId2) == 21) fEffWeight *= 0.008;//See discussion below
  return;
}
LeptonsWW::lepton LeptonsWW::convert(baconhep::TGenParticle *iPart) { 
  TLorentzVector lVec;
  double pPt = iPart->pt;//smearLep(iPart->pdgId,iPart->pt);
  lVec.SetPtEtaPhiM(pPt,iPart->eta,iPart->phi,iPart->mass); 
  lepton lLepton;
  lLepton.first  = lVec;
  lLepton.second = iPart->pdgId;
  return lLepton;
}
int LeptonsWW::parent(baconhep::TGenParticle *iGen, TClonesArray *iGenParts) { 
  int pParentPdgId = iGen->pdgId;
  baconhep::TGenParticle *pGen = iGen;
  while(pParentPdgId == iGen->pdgId) { 
    pParentPdgId = pGen->parent > -1 ? ((baconhep::TGenParticle*) iGenParts->At(pGen->parent))->pdgId : -1;
    if(pGen->parent > -1) pGen = (baconhep::TGenParticle*) iGenParts->At(pGen->parent);
  } 
  return pParentPdgId;
}
void LeptonsWW::getAllLeps(std::vector<lepton> &iLeptons,TClonesArray *iGenParts,bool iFromBoson) { 
  for(int i0 = 0; i0 < iGenParts->GetEntriesFast(); i0++) { 
    baconhep::TGenParticle *pGen = (baconhep::TGenParticle*) iGenParts->At(i0);
    if(fabs(pGen->status) != 1) continue;
    if(fabs(pGen->pdgId) != 11 && fabs(pGen->pdgId) != 13 && fabs(pGen->pdgId) != 15) continue; 
    if(fabs(pGen->eta)    > 2.4 || fabs(pGen->pt) < 10) continue;
    int pParentPdgId = parent(pGen,iGenParts);// pGen->parent > -1 ? ((baconhep::TGenParticle*) iGenParts->At(pGen->parent))->pdgId : -1;
    if(iFromBoson && !(fabs(pParentPdgId)  == 23 ||   fabs(pParentPdgId)  == 24)) continue; //!!!! Leptons can descend form leptons
    lepton pLepton = convert(pGen); 
    bool pInsert = false;
    for(std::vector<lepton>::iterator pIter =  iLeptons.begin(); pIter != iLeptons.end(); pIter++) {
      if(pLepton.first.Pt() > pIter->first.Pt()) { 
	iLeptons.insert(pIter,pLepton);
	pInsert = true;
      }
      if(pInsert) break;
    }
    if(!pInsert) iLeptons.push_back(pLepton); 
  }
}
void LeptonsWW::reorderZ(std::vector<lepton> &iLeptons,std::vector<int> iOrder) { 
  //Find the best Z candidate
  double lDiff    = -1000; 
  int lId0 = -1; 
  int lId1 = -1;
  for(unsigned int i0 = 0; i0 < iLeptons.size(); i0++) { 
    for(unsigned int i1 = i0; i1 < iLeptons.size(); i1++) { 
      if(fabs(iLeptons[i0].second) != fabs(iLeptons[i1].second)) continue;
      double pMass = (iLeptons[i0].first + iLeptons[i1].first).M();
      double pDiff =  91.3-pMass;
      if(fabs(lDiff) > fabs(pDiff)) { 
	lDiff = pDiff;
	lId0 = i0;
	lId1 = i1;
      }
    }
  }
  if(lId0 == -1) return;
  if(lDiff+91.3 > 30  && lDiff+91.3 < 120 ) {
    iOrder[0] = lId0; 
    iOrder[1] = lId1; 
  } else { return;}
  //Find Second Best Z candidate
  lDiff = -1000;
  lId0  = -1;
  lId1  = -1;
  for(int i0 = 0; i0 < int(iLeptons.size()); i0++) { 
    if(i0 == iOrder[0] || i0 == iOrder[1]) continue;
    for(int i1 = i0; i1 < int(iLeptons.size()); i1++) { 
      if(i1 == iOrder[0] || i1 == iOrder[1]) continue;
      if(fabs(iLeptons[i0].second) != fabs(iLeptons[i1].second)) continue;
      double pMass = (iLeptons[i0].first + iLeptons[i1].first).M();
      double pDiff =  91.3-pMass;
      if(fabs(lDiff) > fabs(pDiff)) { 
	lDiff = pDiff;
	lId0 = i0;
	lId1 = i1;
      }
    }
  }
  if(lDiff+91.3 > 10  && lDiff+91.3 < 120 ) {
    iOrder[2] = lId0; 
    iOrder[3] = lId1; 
    return;
  }
  int iFill = 2;
  for(int i0 = 0; i0 < int(iLeptons.size()); i0++) { 
    if(i0 == iOrder[0] || i0 == iOrder[1]) continue;
    iOrder[iFill] = i0;
    iFill++;
  }
}
void LeptonsWW::fill(float &iPt,float &iEta,float &iPhi,float &iM,int &iId,lepton &iLepton)  { 
  iPt  = iLepton.first.Pt();
  iPhi = iLepton.first.Phi();
  iEta = iLepton.first.Eta();  
  iM   = iLepton.first.M();
  iId  = iLepton.second;
}
bool LeptonsWW::findRealLep(Gen &iGen) { 
  std::vector<lepton> lLeptons;
  getAllLeps(lLeptons,iGen.fGenParts,true);//false);
  std::vector<int> lOrder; for(int i0 = 0; i0 < 4; i0++) lOrder.push_back(i0);
  reorderZ(lLeptons,lOrder);
  //Do we wnat to require closes to the 
  if(lLeptons.size() > 0) fill(iGen.fPt1,iGen.fEta1,iGen.fPhi1,iGen.fM1,iGen.fId1,lLeptons[lOrder[0]]);
  if(lLeptons.size() > 1) fill(iGen.fPt2,iGen.fEta2,iGen.fPhi2,iGen.fM2,iGen.fId2,lLeptons[lOrder[1]]);
  if(lLeptons.size() > 2) fill(iGen.fPt3,iGen.fEta3,iGen.fPhi3,iGen.fM3,iGen.fId3,lLeptons[lOrder[2]]);
  if(lLeptons.size() > 3) fill(iGen.fPt4,iGen.fEta4,iGen.fPhi4,iGen.fM4,iGen.fId4,lLeptons[lOrder[3]]);
  for(int i0 = 0; i0 < TMath::Min(int(lLeptons.size()),4); i0++) { 
    fillWeight(lLeptons[i0].second,lLeptons[i0].first.Pt());
  }
  return (lLeptons.size() > 1);
}
