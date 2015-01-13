#include "../include/GenJet.hh"
#include "../include/Tools.hh"
#include "TMath.h"
#include "TLorentzVector.h"
#include <iostream>

GenJet::GenJet(TTree *iTree,TTree *iOTree) { 
  loadGen(iTree);
  setupGenTree(iOTree);
}
GenJet::~GenJet() { }
void GenJet::loadGen(TTree *iTree) { 
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
void GenJet::setupGenTree(TTree *iTree) {
  iTree->Branch("decaymode",&fDecay    ,"fDecay/I");
  iTree->Branch("mcweight" ,&fMCWeight ,"fMCWeight/F");
  //iTree->Branch("effweight",&fEffWeight,"fEffWeight/F");

  iTree->Branch("v_pt"    ,&fVPt,   "fVPt/F");
  iTree->Branch("v_y"     ,&fVEta,  "fVEta/F");
  iTree->Branch("v_phi"   ,&fVPhi,  "fVPhi/F");
  iTree->Branch("v_m"     ,&fVM,    "fVM/F");
  //  iTree->Branch("v_id"    ,&fVId,   "fVId/I");

  iTree->Branch("dm_pt"    ,&fV1Pt,   "fV1Pt/F");
  iTree->Branch("dm_eta"   ,&fV1Eta,  "fV1Eta/F");
  iTree->Branch("dm_phi"   ,&fV1Phi,  "fV1Phi/F");
  iTree->Branch("dm_m"     ,&fV1M,    "fV1M/F");
  iTree->Branch("dm_id"    ,&fV1Id,   "fV1Id/I");

  iTree->Branch("j1_pt"    ,&fJ1Pt,   "fV1Pt/F");
  iTree->Branch("j1_eta"   ,&fJ1Eta,  "fV1Eta/F");
  iTree->Branch("j1_phi"   ,&fJ1Phi,  "fV1Phi/F");
  iTree->Branch("j1_m"     ,&fJ1M,    "fV1M/F");
  iTree->Branch("j1_id"    ,&fJ1Id,   "fV1Id/I");

}
void GenJet::clearGen(int i0){
  fGenBr->GetEntry(i0);
  fGenParts->Clear(); fGenPartBr  ->GetEntry(i0);  
  fGenJets ->Clear(); fGenJetBr   ->GetEntry(i0);  
  fDecay = -1;
  fMCWeight  = fGen->xs;
  fEffWeight = 1.;

  fVPt   = 0;
  fVEta  = 0;
  fVPhi  = 0;
  fVM    = 0;
  //fVId   = 0; 

  fV1Pt  = 0;
  fV1Eta = 0;
  fV1Phi = 0;
  fV1M   = 0;
  fV1Id  = 0; 

  fJ1Pt  = 0;
  fJ1Eta = 0;
  fJ1Phi = 0;
  fJ1M   = 0;
  fJ1Id  = 0; 
}
void GenJet::replace(int iId,baconhep::TGenParticle *iGen,std::vector<baconhep::TGenParticle*> &iBoson,bool iRemove) { 
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
void GenJet::findBosons(std::vector<baconhep::TGenParticle*> & iBoson,int iIdMin,int iIdMax) { 
  for(int i0 = 0; i0 < fGenParts->GetEntriesFast(); i0++) { 
    baconhep::TGenParticle *pGen = (baconhep::TGenParticle*) fGenParts->At(i0);
    if(fabs(pGen->pdgId)  <  iIdMin ) continue;
    if(fabs(pGen->pdgId)  >  iIdMax ) continue;
    if(fabs(pGen->status) <  20     ) continue;
    if(fabs(pGen->status) >  30      && iIdMin != 32) continue;
    if(fabs(pGen->status) <  40      && iIdMin == 32) continue;
    if(pGen->pt == 0) continue;
    int pBoson = -1;
    for(unsigned int i1 = 0; i1 < iBoson.size(); i1++) { 
      if(Tools::deltaR(pGen->eta,pGen->phi,iBoson[i1]->eta,iBoson[i1]->phi) < 0.5 || (pGen->pdgId == 32 && iBoson[i1]->pdgId == 32) ) {
	if(pGen->status >= iBoson[i1]->status  ) pBoson = i1; 
	continue;
      }
      if(pGen->pt > iBoson[i1]->pt && (pBoson == -1 && iBoson.size() > 1) && i1 > 0) pBoson = 10 + i1;
    }
    //if(iIdMin == 24) std::cout <<" 1---> Size " << iBoson.size() << " -- " << pBoson << " -- " << pGen->status << " -- " << pGen->pdgId << " -- " << i0 << " -- " << pGen->pt << " -- " << pGen->phi << std::endl;
    if(pBoson <  10 && pBoson > -1) replace(pBoson,   pGen,iBoson,true);
    //if(iIdMin == 32) std::cout <<" 2---> Size " << iBoson.size() << " -- " << pBoson << std::endl;
    if(pBoson >   9) replace(pBoson-10,pGen,iBoson,false);
    if(pBoson == -1) iBoson.push_back(pGen);
  }
}
void GenJet::findBosonsPythia6(std::vector<baconhep::TGenParticle*> & iBoson,int iIdMin,int iIdMax,int iStatus) { 
  for(int i0 = 0; i0 < fGenParts->GetEntriesFast(); i0++) { 
    baconhep::TGenParticle *pGen = (baconhep::TGenParticle*) fGenParts->At(i0);
    if(fabs(pGen->pdgId)  <  iIdMin ) continue;
    if(fabs(pGen->pdgId)  >  iIdMax ) continue;
    if(fabs(pGen->status) != iStatus) continue;
    if(pGen->pt == 0) continue;
    int pBoson = -1;
    for(unsigned int i1 = 0; i1 < iBoson.size(); i1++) { 
      if(Tools::deltaR(pGen->eta,pGen->phi,iBoson[i1]->eta,iBoson[i1]->phi) < 0.5 || (pGen->pdgId == 32 && iBoson[i1]->pdgId == 32) ) {
	if(pGen->status < iBoson[i1]->status  ) pBoson = i1; 
	continue;
      }
      if(pGen->pt > iBoson[i1]->pt && (pBoson == -1 && iBoson.size() > 1) && i1 > 0) pBoson = 10 + i1;
    }
    if(pBoson <  10 && pBoson > -1) replace(pBoson,   pGen,iBoson,true);
    if(pBoson >   9)                replace(pBoson-10,pGen,iBoson,false);
    if(pBoson == -1) iBoson.push_back(pGen);
  }
}
void GenJet::fillGen(bool iPy6) {
  std::vector<baconhep::TGenParticle*> lBosons;
  std::vector<baconhep::TGenParticle*> lHiggs;
  std::vector<baconhep::TGenParticle*> lJets;
  //findBosons(lBosons,24,24); 
  if(iPy6) { 
    //findBosonsPythia6(lBosons,18,18,1); 
    //findBosonsPythia6(lBosons,16,16,1); 
    //findBosonsPythia6(lBosons,14,14,1); 
    //findBosonsPythia6(lBosons,12,12,1); 
    findBosonsPythia6(lBosons,25,25,2); 
    findBosonsPythia6(lHiggs ,25,25,3); 
    //findBosons(lBosons, 0,10); 
    //findBosons(lBosons,21,21);
  } else {
    findBosons(lBosons,32,32); 
    findBosons(lBosons, 0,10); 
    findBosons(lBosons,21,21); 
  }
  TLorentzVector lVMass; 
  lVMass .SetPtEtaPhiM(0.0,0,0,0);
  //if( lBosons.size() > 2) std::cout << "Xsoze---> " << lBosons.size() << std::endl;
  //if( lBosons.size() > 2) std::cout <<" -XXX-> " << lBosons[1]->pdgId << std::endl;
  int lMin = 2; 
  if(iPy6) lMin = 1000;
  for(int i0 = 0; i0 < TMath::Min(int(lBosons.size()),lMin); i0++) { 
    TLorentzVector pVec;
    pVec.SetPtEtaPhiM(lBosons[i0]->pt,lBosons[i0]->eta,lBosons[i0]->phi,lBosons[i0]->mass);
    if(i0 == 0) { 
      fV1Pt  = lBosons[i0]->pt;
      fV1Eta = lBosons[i0]->eta;
      fV1Phi = lBosons[i0]->phi;
      fV1M   = lBosons[i0]->mass;
      fV1Id  = lBosons[i0]->pdgId;
    }
    if(i0 == 1) { 
      fJ1Pt  = lBosons[i0]->pt;
      fJ1Eta = lBosons[i0]->eta;
      fJ1Phi = lBosons[i0]->phi;
      fJ1M   = lBosons[i0]->mass;
      fJ1Id  = lBosons[i0]->pdgId;
    }
    lVMass += pVec;
  }
  if(lVMass.Pt() == 0)  return;
  //Should be checked with Gen Met
  fVPt     = lVMass.Pt();
  fVEta    = lVMass.Eta();
  fVPhi    = lVMass.Phi();
  fVM      = lVMass.M();
  if(lHiggs.size() > 0)  { 
    fVPt     = lHiggs[0]->pt;
    fVEta    = lHiggs[0]->eta;
    fVPhi    = lHiggs[0]->phi;
    fVM      = lHiggs[0]->mass;
  }
  if(iPy6) { 
    fV1Pt  = lVMass.Pt();
    fV1Eta = lVMass.Eta();
    fV1Phi = lVMass.Phi();
    fV1M   = lVMass.M();
    fV1Id  = 18;
  }
}
