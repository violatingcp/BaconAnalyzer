#include "../include/GenV.hh"
#include "../include/Tools.hh"
#include "TMath.h"
#include "TLorentzVector.h"
#include <iostream>

GenV::GenV(TTree *iTree,TTree *iOTree) { 
  loadGen(iTree);
  setupGenTree(iOTree);
}
GenV::~GenV() { }
void GenV::loadGen(TTree *iTree) { 
  baconhep::TGenEventInfo::Class() ->IgnoreTObjectStreamer();
  baconhep::TGenParticle ::Class() ->IgnoreTObjectStreamer();
  baconhep::TGenJet      ::Class() ->IgnoreTObjectStreamer();
  fGenParts          = new TClonesArray("baconhep::TGenParticle");
  fGenJets           = new TClonesArray("baconhep::TGenJet");
  fGenFatJets        = new TClonesArray("baconhep::TGenJet");
  iTree->SetBranchAddress("GenEvtInfo"   ,&fGen);         fGenBr         = iTree->GetBranch("GenEvtInfo"  ); 
  iTree->SetBranchAddress("GenParticle"  ,&fGenParts);    fGenPartBr     = iTree->GetBranch("GenParticle" ); 
  iTree->SetBranchAddress("GenJet"       ,&fGenJets);     fGenJetBr      = iTree->GetBranch("GenJet"      ); 
  iTree->SetBranchAddress("GenFatJet"    ,&fGenFatJets);  fGenFatJetBr   = iTree->GetBranch("GenFatJet"   ); 
}
void GenV::setupGenTree(TTree *iTree) {
  iTree->Branch("decaymode",&fDecay    ,"fDecay/I");
  iTree->Branch("mcweight" ,&fMCWeight ,"fMCWeight/F");
  //iTree->Branch("effweight",&fEffWeight,"fEffWeight/F");

  iTree->Branch("dm_pt"    ,&fDPt,   "fDPt/F");
  iTree->Branch("dm_eta"   ,&fDEta,  "fDEta/F");
  iTree->Branch("dm_phi"   ,&fDPhi,  "fDPhi/F");
  iTree->Branch("dm_m"     ,&fDM,    "fDM/F");
  iTree->Branch("dm_id"    ,&fDId,   "fDId/I");
  //  iTree->Branch("v_id"    ,&fVId,   "fVId/I");
  iTree->Branch("d1_pt"    ,&fD1Pt,   "fD1Pt/F");
  iTree->Branch("d1_eta"   ,&fD1Eta,  "fD1Eta/F");
  iTree->Branch("d1_phi"   ,&fD1Phi,  "fD1Phi/F");
  iTree->Branch("d1_m"     ,&fD1M,    "fD1M/F");
  iTree->Branch("d1_id"    ,&fD1Id,   "fD1Id/I");

  iTree->Branch("d2_pt"    ,&fD2Pt,   "fD2Pt/F");
  iTree->Branch("d2_eta"   ,&fD2Eta,  "fD2Eta/F");
  iTree->Branch("d2_phi"   ,&fD2Phi,  "fD2Phi/F");
  iTree->Branch("d2_m"     ,&fD2M,    "fD2M/F");
  iTree->Branch("d2_id"    ,&fD2Id,   "fD2Id/I");

  iTree->Branch("q1_pt"    ,&fQ1Pt,   "fQ1Pt/F");
  iTree->Branch("q1_eta"   ,&fQ1Eta,  "fQ1Eta/F");
  iTree->Branch("q1_phi"   ,&fQ1Phi,  "fQ1Phi/F");
  iTree->Branch("q1_m"     ,&fQ1M,    "fQ1M/F");
  iTree->Branch("q1_id"    ,&fQ1Id,   "fQ1Id/I");

  iTree->Branch("q2_pt"    ,&fQ2Pt,   "fQ2Pt/F");
  iTree->Branch("q2_eta"   ,&fQ2Eta,  "fQ2Eta/F");
  iTree->Branch("q2_phi"   ,&fQ2Phi,  "fQ2Phi/F");
  iTree->Branch("q2_m"     ,&fQ2M,    "fQ2M/F");
  iTree->Branch("q2_id"    ,&fQ2Id,   "fQ2Id/I");

  iTree->Branch("vj_pt"    ,&fBJPt,   "fBJPt/F");
  iTree->Branch("vj_eta"   ,&fBJEta,  "fBJEta/F");
  iTree->Branch("vj_phi"   ,&fBJPhi,  "fBJPhi/F");
  iTree->Branch("vj_m"     ,&fBJM,    "fBJM/F");
  iTree->Branch("vj_id"    ,&fBJId,   "fBJId/I");

  iTree->Branch("v_pt"    ,&fBPt,   "fBPt/F");
  iTree->Branch("v_eta"   ,&fBEta,  "fBEta/F");
  iTree->Branch("v_phi"   ,&fBPhi,  "fBPhi/F");
  iTree->Branch("v_m"     ,&fBM,    "fBM/F");
  iTree->Branch("v_id"    ,&fBId,   "fBId/I");

  iTree->Branch("x_pt"    ,&fV1Pt,   "fV1Pt/F");
  iTree->Branch("x_eta"   ,&fV1Eta,  "fV1Eta/F");
  iTree->Branch("x_phi"   ,&fV1Phi,  "fV1Phi/F");
  iTree->Branch("x_m"     ,&fV1M,    "fV1M/F");
  iTree->Branch("x_id"    ,&fV1Id,   "fV1Id/I");
  /*
  iTree->Branch("j1_pt"    ,&fJ1Pt,   "fV1Pt/F");
  iTree->Branch("j1_eta"   ,&fJ1Eta,  "fV1Eta/F");
  iTree->Branch("j1_phi"   ,&fJ1Phi,  "fV1Phi/F");
  iTree->Branch("j1_m"     ,&fJ1M,    "fV1M/F");
  iTree->Branch("j1_id"    ,&fJ1Id,   "fV1Id/I");
  */
}
void GenV::clearGen(int i0){
  fGenBr->GetEntry(i0);
  fGenParts   ->Clear(); fGenPartBr  ->GetEntry(i0);  
  fGenJets    ->Clear(); fGenJetBr   ->GetEntry(i0);  
  fGenFatJets ->Clear(); fGenFatJetBr->GetEntry(i0);  
  fDecay     = -1;
  fMCWeight  = fGen->xs;
  fEffWeight = 1.;

  fDPt   = 0;
  fDEta  = 0;
  fDPhi  = 0;
  fDM    = 0;
  fDId   = 0; 

  fD1Pt  = 0;
  fD1Eta = 0;
  fD1Phi = 0;
  fD1M   = 0;
  fD1Id  = 0;

  fD2Pt  = 0;
  fD2Eta = 0;
  fD2Phi = 0;
  fD2M   = 0;
  fD2Id  = 0;

  fQ1Pt  = 0;
  fQ1Eta = 0;
  fQ1Phi = 0;
  fQ1M   = 0;
  fQ1Id  = 0;

  fQ2Pt  = 0;
  fQ2Eta = 0;
  fQ2Phi = 0;
  fQ2M   = 0;
  fQ2Id  = 0;

  fBJPt  = 0;
  fBJEta = 0;
  fBJPhi = 0;
  fBJM   = 0;

  fBPt  = 0;
  fBEta = 0;
  fBPhi = 0;
  fBM   = 0;
  fBId  = 0;

  fV1Pt  = 0;
  fV1Eta = 0;
  fV1Phi = 0;
  fV1M   = 0;
  fV1Id  = 0; 
}
void GenV::replace(int iId,baconhep::TGenParticle *iGen,std::vector<baconhep::TGenParticle*> &iBoson,bool iRemove) { 
  int pId = 0; 
  for(std::vector<baconhep::TGenParticle*>::iterator pIter = iBoson.begin();  pIter != iBoson.end(); pIter++) { 
    if(iId == pId) { 
      iBoson.insert(pIter,iGen);
      if(iRemove) iBoson.erase(iBoson.begin()+pId+1); 
      return;
    }
    pId++;
  }
}
void GenV::findBosons(std::vector<baconhep::TGenParticle*> & iBoson,int iIdMin,int iIdMax,int iStatus) { 
  for(int i0 = 0; i0 < fGenParts->GetEntriesFast(); i0++) { 
    baconhep::TGenParticle *pGen = (baconhep::TGenParticle*) fGenParts->At(i0);
    if(fabs(pGen->pdgId)  <  iIdMin  ) continue;
    if(fabs(pGen->pdgId)  >  iIdMax  ) continue;
    if(fabs(pGen->status) != iStatus ) continue;
    if(pGen->pt == 0) continue;
    int pBoson = -1;
    for(unsigned int i1 = 0; i1 < iBoson.size(); i1++) { 
      if(Tools::deltaR(pGen->eta,pGen->phi,iBoson[i1]->eta,iBoson[i1]->phi) < 0.5 || (pGen->pdgId == 32 && iBoson[i1]->pdgId == 32) ) {
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
void GenV::findDaughters(baconhep::TGenParticle *iParticle,std::vector<baconhep::TGenParticle*> & iDaughter,int iIdMin,int iIdMax) { 
  int lId = -1;
  for(int i0 = 0; i0 < fGenParts->GetEntriesFast(); i0++) { 
    baconhep::TGenParticle *pGen = (baconhep::TGenParticle*) fGenParts->At(i0);
    if(pGen == iParticle) lId = i0;
  }
  for(int i0 = 0; i0 < fGenParts->GetEntriesFast(); i0++) { 
    baconhep::TGenParticle *pGen = (baconhep::TGenParticle*) fGenParts->At(i0);
    if(pGen->parent != lId) continue;
    if(fabs(pGen->pdgId)  <  iIdMin  ) continue;
    if(fabs(pGen->pdgId)  >  iIdMax  ) continue;
    if(pGen->pt == 0) continue;
    int pBoson = -1;
    for(unsigned int i1 = 0; i1 < iDaughter.size(); i1++) { 
      if(Tools::deltaR(pGen->eta,pGen->phi,iDaughter[i1]->eta,iDaughter[i1]->phi) < 0.5 || (pGen->pdgId == 32 && iDaughter[i1]->pdgId == 32) ) {
	if(pGen->status >= iDaughter[i1]->status  ) pBoson = i1; 
	continue;
      }
      if(pGen->pt > iDaughter[i1]->pt && (pBoson == -1 && iDaughter.size() > 1) && i1 > 0) pBoson = 10 + i1;
    }
    if(pBoson <  10 && pBoson > -1) replace(pBoson,   pGen,iDaughter,true);
    if(pBoson >   9)                replace(pBoson-10,pGen,iDaughter,false);
    if(pBoson == -1) iDaughter.push_back(pGen);
  }
}
void GenV::fillGen(int iPy6) {
  std::vector<baconhep::TGenParticle*> lBosons;
  std::vector<baconhep::TGenParticle*> lDM;
  std::vector<baconhep::TGenParticle*> lDaughters;
  if(iPy6 == 0) findBosons(lDM,9100012,9100012,1);
  if(iPy6 == 1) findBosons(lDM,1000022,1000022,1);
  if(iPy6 == 2) findBosons(lDM,25,25,1);
  TLorentzVector lVMass; 
  lVMass .SetPtEtaPhiM(0.0,0,0,0);
  for(int i0 = 0; i0 < TMath::Min(int(lDM.size()),2); i0++) { 
    TLorentzVector pVec;
    pVec.SetPtEtaPhiM(lDM[i0]->pt,lDM[i0]->eta,lDM[i0]->phi,lDM[i0]->mass);
    if(i0 == 0) { 
      fD1Pt  = lDM[i0]->pt;
      fD1Eta = lDM[i0]->eta;
      fD1Phi = lDM[i0]->phi;
      fD1M   = lDM[i0]->mass;
      fD1Id  = lDM[i0]->pdgId;
    }
    if(i0 == 1) { 
      fD2Pt  = lDM[i0]->pt;
      fD2Eta = lDM[i0]->eta;
      fD2Phi = lDM[i0]->phi;
      fD2M   = lDM[i0]->mass;
      fD2Id  = lDM[i0]->pdgId;
    }
    lVMass += pVec;
  }
  if(lVMass.Pt() == 0)  return;
  //Should be checked with Gen Met
  fDPt     = lVMass.Pt();
  fDEta    = lVMass.Eta();
  fDPhi    = lVMass.Phi();
  fDM      = lVMass.M();

  fV1Pt     = lVMass.Pt();
  fV1Eta    = lVMass.Eta();
  fV1Phi    = lVMass.Phi();
  fV1M      = lVMass.M();

  if(iPy6 != 0) {  
    findBosons(lBosons,23,24,3); 
  } else {
    findBosons(lBosons,23,24,62); 
  }
  findDaughters(lBosons[0],lDaughters,0,20);
  lVMass .SetPtEtaPhiM(0.0,0,0,0);  
  for(int i0 = 0; i0 < TMath::Min(int(lDaughters.size()),2); i0++) { 
    TLorentzVector pVec;
    pVec.SetPtEtaPhiM(lDaughters[i0]->pt,lDaughters[i0]->eta,lDaughters[i0]->phi,lDaughters[i0]->mass);
    if(i0 == 0) { 
      fQ1Pt  = lDaughters[i0]->pt;
      fQ1Eta = lDaughters[i0]->eta;
      fQ1Phi = lDaughters[i0]->phi;
      fQ1M   = lDaughters[i0]->mass;
      fQ1Id  = lDaughters[i0]->pdgId;
    }
    if(i0 == 1) { 
      fQ2Pt  = lDaughters[i0]->pt;
      fQ2Eta = lDaughters[i0]->eta;
      fQ2Phi = lDaughters[i0]->phi;
      fQ2M   = lDaughters[i0]->mass;
      fQ2Id  = lDaughters[i0]->pdgId;
    }
    lVMass += pVec;
  }
  if(lVMass.Pt() == 0)  return;
  //std::cout << " ===> " << iPy6 << " -- " << lDM.size() << " -- " << lBosons.size() << " -- " << lDaughters.size() << " -- j1: " << fQ1Pt << " -- j2: " << fQ2Pt << " -- " << lVMass.Pt() << " -- " << lVMass.Eta() << std::endl;
  //Should be checked with Gen Met
  fBJPt     = lVMass.Pt();
  fBJEta    = lVMass.Eta();
  fBJPhi    = lVMass.Phi();
  fBJM      = lVMass.M();

  fBPt      = lBosons[0]->pt;
  fBEta     = lBosons[0]->eta;
  fBPhi     = lBosons[0]->phi;
  fBM       = lBosons[0]->mass;
  fBId      = lBosons[0]->pdgId;
}
