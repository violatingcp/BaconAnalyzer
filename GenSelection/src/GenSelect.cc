#include "../include/GenSelect.hh"
#include "../include/Tools.hh"
#include "TMath.h"
#include "TLorentzVector.h"
#include <iostream>

GenSelect::GenSelect(TTree *iTree,TTree *iOTree) { 
  loadGen(iTree);
  setupGenTree(iOTree);
}
GenSelect::~GenSelect() { 
}
void GenSelect::loadGen(TTree *iTree) { 
  baconhep::TGenEventInfo::Class() ->IgnoreTObjectStreamer();
  baconhep::TGenParticle ::Class() ->IgnoreTObjectStreamer();
  fGen               = new baconhep::TGenEventInfo();
  fGenParts          = new TClonesArray("baconhep::TGenParticle");
  iTree->SetBranchAddress("GenEvtInfo"   ,&fGen);         fGenBr         = iTree->GetBranch("GenEvtInfo"  ); 
  iTree->SetBranchAddress("GenParticle"  ,&fGenParts);    fGenPartBr     = iTree->GetBranch("GenParticle" ); 
}
void GenSelect::setupGenTree(TTree *iTree) {
  iTree->Branch("decaymode",&fDecay     ,"fDecay/I");
  iTree->Branch("evtweight" ,&fEvtWeight  ,"fEvtWeight/F");
  iTree->Branch("mcweight" ,&fMCWeight  ,"fMCWeight/F");
  iTree->Branch("mcweight2",&fMCWeight2 ,"fMCWeight2/F");
  iTree->Branch("xs"       ,&fXS        ,"fXS/F");
  iTree->Branch("xs2"      ,&fXS2       ,"fXS2/F");
  iTree->Branch("v_pt"    ,&fVPt,   "fVPt/F");
  iTree->Branch("v_y"     ,&fVEta,  "fVEta/F");
  iTree->Branch("v_phi"   ,&fVPhi,  "fVPhi/F");
  iTree->Branch("v_m"     ,&fVM,    "fVM/F");
  iTree->Branch("v_id"    ,&fVId,   "fVId/I");

  iTree->Branch("pt_1"    ,&fPt1,   "fPt1/F");
  iTree->Branch("eta_1"   ,&fEta1,  "fEta1/F");
  iTree->Branch("phi_1"   ,&fPhi1,  "fPhi1/F");
  iTree->Branch("m_1"     ,&fM1,    "fM1/F");
  iTree->Branch("cosphi_1",&fCosPhi1,   "fCosPhi1/F");
  iTree->Branch("id_1"    ,&fId1,   "fId1/I");

  iTree->Branch("pt_2"    ,&fPt2,   "fPt2/F");
  iTree->Branch("eta_2"   ,&fEta2,  "fEta2/F");
  iTree->Branch("phi_2"   ,&fPhi2,  "fPhi2/F");
  iTree->Branch("m_2"     ,&fM2,    "fM2/F");
  iTree->Branch("cosphi_2",&fCosPhi2,   "fCosPhi2/F");
  iTree->Branch("id_2"    ,&fId2,   "fId2/I");
  iTree->Branch("npartons",&fNPartons,"fNPartons/I");
}
void GenSelect::clearGen(int i0){
  fGenBr->GetEntry(i0);
  fGenParts->Clear(); fGenPartBr  ->GetEntry(i0);  
  fDecay = -1;
  fEvtWeight    = 1.;
  fVPt       = 0;
  fVEta      = 0;
  fVPhi      = 0;
  fVM        = 0;
  fVId       = 0; 
 
  fPt1       = 0;
  fEta1      = 0;
  fPhi1      = 0;
  fM1        = 0;
  fCosPhi1   = 0; 
  fId1       = 0;
 
  fPt2       = 0;
  fEta2      = 0;
  fPhi2      = 0;
  fM2        = 0;
  fCosPhi2   = 0; 
  fId2       = 0;
  fNPartons  = 0; 
}
void GenSelect::replace(int iId,baconhep::TGenParticle *iGen,std::vector<baconhep::TGenParticle*> &iBoson,bool iRemove) { 
  int pId = 0; 
  for(std::vector<baconhep::TGenParticle*>::iterator pIter = iBoson.begin();  pIter != iBoson.end(); pIter++) { 
    if(iId == pId) { 
      //std::cout << " Check - " << (*pIter)->pt << " -- " << (*pIter)->status << " -- " << iBoson.size() << std::endl;
      iBoson.insert(pIter,iGen);
      if(iRemove) iBoson.erase(iBoson.begin()+pId+1); 
      return;
    }
    pId++;
  }
}
void GenSelect::findBosons(std::vector<baconhep::TGenParticle*> & iBoson,int iIdMin,int iIdMax,int iStatusMin,int iStatusMax) { 
  for(int i0 = 0; i0 < fGenParts->GetEntriesFast(); i0++) { 
    baconhep::TGenParticle *pGen = (baconhep::TGenParticle*) fGenParts->At(i0);
    //std::cout <<" ===> " << pGen->pdgId << " -- " << pGen->status << std::endl;
    if(fabs(pGen->pdgId)  <  iIdMin     ) continue;
    if(fabs(pGen->pdgId)  >  iIdMax     ) continue;
    if(fabs(pGen->status) <  iStatusMin ) continue;
    if(fabs(pGen->status) >  iStatusMax ) continue;
    if(pGen->pt == 0) continue;
    int pBoson = -1;
    for(unsigned int i1 = 0; i1 < iBoson.size(); i1++) { 
      if(Tools::deltaR(pGen->eta,pGen->phi,iBoson[i1]->eta,iBoson[i1]->phi) < 1.2 ) {
	if(pGen->status >= iBoson[i1]->status  ) pBoson = i1; 
	continue;
      }
      if(pGen->pt > iBoson[i1]->pt && (pBoson == -1 && iBoson.size() > 1) && i1 > 0) pBoson = 10 + i1;
    }
    if(pBoson <  10 && pBoson > -1) replace(pBoson,   pGen,iBoson,true);
    if(pBoson >   9) replace(pBoson-10,pGen,iBoson,false);
    if(pBoson == -1) iBoson.push_back(pGen);
  }
}
void GenSelect::npartons() { 
  fNPartons = 0; 
  for(int i0 = 0; i0 < fGenParts->GetEntriesFast(); i0++) { 
    baconhep::TGenParticle *pGen = (baconhep::TGenParticle*) fGenParts->At(i0);
    if(pGen->status != 3) continue;
    if(pGen->parent < 0) continue;
    if(pGen->pt     < 30) continue;
    if((fabs(pGen->pdgId) > -1 && fabs(pGen->pdgId) < 7) || fabs(pGen->pdgId) == 21) fNPartons++;
  }
}
void GenSelect::fillGen(double iTot,std::vector<TLorentzVector> &iVec) {
  fXS2       = fGen->xs;
  fMCWeight2 = fXS2/iTot;
  fEvtWeight = fGen->weight;
  std::vector<baconhep::TGenParticle*> lBosons;
  std::vector<baconhep::TGenParticle*> lDaughters;
  findBosons(lBosons,23,23,20,70);
  findBosons(lBosons,22,22,20,70);
  findBosons(lBosons,24,24,20,70);
  if(lBosons.size() == 0) findBosons(lBosons,23,23,-1,70);
  if(lBosons.size() == 0) findBosons(lBosons,24,24,-1,70);
  if(lBosons.size() == 0) findBosons(lBosons,25,25,-1,70);
  if(lBosons.size() == 0) findBosons(lBosons,32,32,-1,70);
  findBosons(lDaughters,1000022,1000022,-1,70);
  findBosons(lDaughters,12,12,0,2);
  findBosons(lDaughters,14,14,0,2);
  findBosons(lDaughters,16,16,0,2);
  findBosons(lDaughters,11,11,0,2);
  findBosons(lDaughters,13,13,0,2);
  findBosons(lDaughters,15,15,0,2);
  TLorentzVector lLep1,lLep2;
  if(lDaughters.size() > 0) { 
    fPt1  = lDaughters[0]->pt; 
    fEta1 = lDaughters[0]->eta; 
    fPhi1 = lDaughters[0]->phi; 
    fM1   = lDaughters[0]->mass; 
    fId1  = lDaughters[0]->pdgId;
    lLep2.SetPtEtaPhiM(fPt1,fEta1,fPhi1,fM1);
  }
  if(lDaughters.size() > 1) { 
    fPt2  = lDaughters[1]->pt; 
    fEta2 = lDaughters[1]->eta; 
    fPhi2 = lDaughters[1]->phi; 
    fM2   = lDaughters[1]->mass;
    fId2  = lDaughters[1]->pdgId;
    lLep2.SetPtEtaPhiM(fPt2,fEta2,fPhi2,fM2);
  }
  if(abs(fId2) %2 == 0 || abs(fId1) % 2 == 0) fDecay = 1.;
  if(abs(fId2) %2 == 0 && abs(fId1) % 2 == 0) fDecay = 2.;
  TLorentzVector lVMass;
  for(int i0 = 0; i0 < TMath::Min(int(lBosons.size()),1); i0++) { 
    TLorentzVector pVec;
    pVec.SetPtEtaPhiM(lBosons[i0]->pt,lBosons[i0]->eta,lBosons[i0]->phi,lBosons[i0]->mass);
    if(fId1 < 0 && lBosons[i0]->pdgId > 0) { lLep1.Boost(-pVec.BoostVector()); fCosPhi1 = cos(lLep1.Angle(pVec.Vect())); }
    if(fId1 > 0 && lBosons[i0]->pdgId < 0) { lLep1.Boost(-pVec.BoostVector()); fCosPhi1 = cos(lLep1.Angle(pVec.Vect())); }
    if(fId2 < 0 && lBosons[i0]->pdgId > 0) { lLep2.Boost(-pVec.BoostVector()); fCosPhi2 = cos(lLep2.Angle(pVec.Vect())); }
    if(fId2 > 0 && lBosons[i0]->pdgId < 0) { lLep2.Boost(-pVec.BoostVector()); fCosPhi2 = cos(lLep2.Angle(pVec.Vect())); }
    fVId   = lBosons[i0]->pdgId;
    lVMass += pVec;
  }
  //Should be checked with Gen Met
  fVPt     = lVMass.Pt();
  fVEta    = lVMass.Eta();
  fVPhi    = lVMass.Phi();
  fVM      = lVMass.M();
  npartons();
  for(unsigned int i0 = 0; i0 < lDaughters.size(); i0++) { 
    if(fabs(lDaughters[i0]->pdgId) == 11 && fabs(lDaughters[i0]->pdgId) ==13) continue;
    TLorentzVector pVec(0,0,0,0);
    pVec.SetPtEtaPhiM(lDaughters[i0]->pt,lDaughters[i0]->eta,lDaughters[i0]->phi,lDaughters[i0]->mass);
    iVec.push_back(pVec);
  } 
}
