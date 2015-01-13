#include "../include/Gen.hh"
#include "../include/Tools.hh"
#include "TMath.h"
#include "TLorentzVector.h"
#include <iostream>

Gen::Gen(TTree *iTree,TTree *iOTree) { 
  loadGen(iTree);
  setupGenTree(iOTree);
}
Gen::~Gen() { 
}
void Gen::loadGen(TTree *iTree) { 
  baconhep::TGenEventInfo::Class() ->IgnoreTObjectStreamer();
  baconhep::TGenParticle ::Class() ->IgnoreTObjectStreamer();
  baconhep::TGenJet      ::Class() ->IgnoreTObjectStreamer();
  fGenParts          = new TClonesArray("baconhep::TGenParticle");
  fGenJets           = new TClonesArray("baconhep::TGenJet");
  iTree->SetBranchAddress("GenEvtInfo"   ,&fGen);         fGenBr         = iTree->GetBranch("GenEvtInfo"  ); 
  iTree->SetBranchAddress("GenParticle"  ,&fGenParts);    fGenPartBr     = iTree->GetBranch("GenParticle" ); 
  iTree->SetBranchAddress("GenJet"       ,&fGenJets);     fGenJetBr      = iTree->GetBranch("GenJet"      ); 
}
void Gen::setupGenTree(TTree *iTree) {
  iTree->Branch("decaymode",&fDecay    ,"fDecay/I");
  iTree->Branch("mcweight" ,&fMCWeight ,"fMCWeight/F");
  //iTree->Branch("effweight",&fEffWeight,"fEffWeight/F");

  iTree->Branch("v_pt"    ,&fVPt,   "fVPt/F");
  iTree->Branch("v_y"     ,&fVEta,  "fVEta/F");
  iTree->Branch("v_phi"   ,&fVPhi,  "fVPhi/F");
  iTree->Branch("v_m"     ,&fVM,    "fVM/F");
  //  iTree->Branch("v_id"    ,&fVId,   "fVId/I");

  iTree->Branch("v1_pt"    ,&fV1Pt,   "fV1Pt/F");
  iTree->Branch("v1_y"     ,&fV1Eta,  "fV1Eta/F");
  iTree->Branch("v1_phi"   ,&fV1Phi,  "fV1Phi/F");
  iTree->Branch("v1_m"     ,&fV1M,    "fV1M/F");
  iTree->Branch("v1_id"    ,&fV1Id,   "fV1Id/I");

  iTree->Branch("v2_pt"    ,&fV2Pt,   "fV2Pt/F");
  iTree->Branch("v2_y"     ,&fV2Eta,  "fV2Eta/F");
  iTree->Branch("v2_phi"   ,&fV2Phi,  "fV2Phi/F");
  iTree->Branch("v2_m"     ,&fV2M,    "fV2M/F");
  iTree->Branch("v2_id"    ,&fV2Id,   "fV2Id/I");

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
  iTree->Branch("cosphi_2"    ,&fCosPhi2,   "fCosPhi2/F");
  iTree->Branch("id_2"    ,&fId2,   "fId2/I");

  iTree->Branch("pt_3"    ,&fPt3,   "fPt3/F");
  iTree->Branch("eta_3"   ,&fEta3,  "fEta3/F");
  iTree->Branch("phi_3"   ,&fPhi3,  "fPhi3/F");
  iTree->Branch("m_3"     ,&fM3,    "fM3/F");
  iTree->Branch("id_3"    ,&fId3,   "fId3/I");

  iTree->Branch("pt_4"    ,&fPt4,   "fPt4/F");
  iTree->Branch("eta_4"   ,&fEta4,  "fEta4/F");
  iTree->Branch("phi_4"   ,&fPhi4,  "fPhi4/F");
  iTree->Branch("m_4"     ,&fM4,    "fM4/F");
  iTree->Branch("id_4"    ,&fId4,   "fId4/I");
}
void Gen::clearGen(int i0){
  fGenBr->GetEntry(i0);
  fGenParts->Clear(); fGenPartBr  ->GetEntry(i0);  
  fGenJets ->Clear(); fGenJetBr   ->GetEntry(i0);  
  fDecay = -1;
  //fMCWeight  = 1.;
  fEffWeight = 1.;
  fV1Pt  = 0;
  fV1Eta = 0;
  fV1Phi = 0;
  fV1M   = 0;
  fV1Id  = 0; 

  fV2Pt  = 0;
  fV2Eta = 0;
  fV2Phi = 0;
  fV2M   = 0;
  fV2Id  = 0; 

  fPt1  = 0;
  fEta1 = 0;
  fPhi1 = 0;
  fM1   = 0;
  fCosPhi1 = 0; 
  fId1  = 0;
 
  fPt2  = 0;
  fEta2 = 0;
  fPhi2 = 0;
  fM2   = 0;
  fCosPhi2 = 0; 
  fId2     = 0;

  fPt3  = 0;
  fEta3 = 0;
  fPhi3 = 0;
  fM3   = 0;
  fId3  = 0;
 
  fPt4  = 0;
  fEta4 = 0;
  fPhi4 = 0;
  fM4   = 0;
  fId4  = 0;
}
void Gen::replace(int iId,baconhep::TGenParticle *iGen,std::vector<baconhep::TGenParticle*> &iBoson,bool iRemove) { 
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
void Gen::findBosons(std::vector<baconhep::TGenParticle*> & iBoson) { 
  for(int i0 = 0; i0 < fGenParts->GetEntriesFast(); i0++) { 
    baconhep::TGenParticle *pGen = (baconhep::TGenParticle*) fGenParts->At(i0);
    if(fabs(pGen->pdgId) != 23 && fabs(pGen->pdgId) != 24) continue;
    int pBoson = -1;
    for(unsigned int i1 = 0; i1 < iBoson.size(); i1++) { 
      if(Tools::deltaR(pGen->eta,pGen->phi,iBoson[i1]->eta,iBoson[i1]->phi) < 0.5) {
	if(pGen->status < iBoson[i1]->status ) pBoson = i1; 
	continue;
      }
      if(pGen->pt > iBoson[i1]->pt && (pBoson == -1 && iBoson.size() > 1)) pBoson = 10 + i1;
    }
    if(pBoson <  10) replace(pBoson,   pGen,iBoson,true);
    if(pBoson >   9) replace(pBoson-10,pGen,iBoson,false);
    if(pBoson == -1) iBoson.push_back(pGen);
  }
}
void Gen::findBosonsPythia6(std::vector<baconhep::TGenParticle*> & iBoson) { 
  for(int i0 = 0; i0 < fGenParts->GetEntriesFast(); i0++) { 
    baconhep::TGenParticle *pGen = (baconhep::TGenParticle*) fGenParts->At(i0);
    if(fabs(pGen->pdgId) != 23 && fabs(pGen->pdgId) != 24) continue;
    int pBoson = -1;
    for(unsigned int i1 = 0; i1 < iBoson.size(); i1++) { 
      if(Tools::deltaR(pGen->eta,pGen->phi,iBoson[i1]->eta,iBoson[i1]->phi) < 0.5) {
	if(pGen->status < iBoson[i1]->status ) pBoson = i1; 
	continue;
      }
      if(pGen->pt > iBoson[i1]->pt && (pBoson == -1 && iBoson.size() > 1)) pBoson = 10 + i1;
    }
    if(pBoson <  10) replace(pBoson,   pGen,iBoson,true);
    if(pBoson >   9) replace(pBoson-10,pGen,iBoson,false);
    if(pBoson == -1) iBoson.push_back(pGen);
  }
}
//void bosonPtCheck()
void Gen::fillGen(bool iPythia6) {
  std::vector<baconhep::TGenParticle*> lBosons;
  if(!iPythia6) findBosons(lBosons); 
  if(iPythia6 ) findBosonsPythia6(lBosons); 
  TLorentzVector lVMass; 
  lVMass .SetPtEtaPhiM(0,0,0,0);
  TLorentzVector lLep1; lLep1.SetPtEtaPhiM(fPt1,fEta1,fPhi1,fM1);
  TLorentzVector lLep2; lLep2.SetPtEtaPhiM(fPt2,fEta2,fPhi2,fM2);
  for(int i0 = 0; i0 < TMath::Min(int(lBosons.size()),2); i0++) { 
    TLorentzVector pVec;
    pVec.SetPtEtaPhiM(lBosons[i0]->pt,lBosons[i0]->eta,lBosons[i0]->phi,lBosons[i0]->mass);
    if(fId1 < 0 && lBosons[i0]->pdgId > 0) { lLep1.Boost(-pVec.BoostVector()); fCosPhi1 = cos(lLep1.Angle(pVec.Vect())); }//std::cout << " ==> 1 " << lLep1.E() << std::endl;}
    if(fId1 > 0 && lBosons[i0]->pdgId < 0) { lLep1.Boost(-pVec.BoostVector()); fCosPhi1 = cos(lLep1.Angle(pVec.Vect())); }//std::cout << " ==> 1 " << lLep1.E() << std::endl;}
    if(fId2 < 0 && lBosons[i0]->pdgId > 0) { lLep2.Boost(-pVec.BoostVector()); fCosPhi2 = cos(lLep2.Angle(pVec.Vect())); }//std::cout << " ==> 2 " << lLep2.E() << std::endl;}
    if(fId2 > 0 && lBosons[i0]->pdgId < 0) { lLep2.Boost(-pVec.BoostVector()); fCosPhi2 = cos(lLep2.Angle(pVec.Vect())); }//std::cout << " ==> 2 " << lLep2.E() << std::endl;}
    lVMass += pVec;
  }
  //std::cout << fId1 << " -- " << fId2 << " -- " << lBosons[0]->pdgId << " -- " << lBosons[1]->pdgId << " -- " << fCosPhi1 << " -- " << fCosPhi2  << std::endl;
  if(lBosons.size() > 0) { 
    fV1Pt  = lBosons[0]->pt;
    fV1Eta = lBosons[0]->eta;
    fV1Phi = lBosons[0]->phi;
    fV1M   = lBosons[0]->mass;
    fV1Id  = lBosons[0]->pdgId;
  }
  if(lBosons.size() > 1) { 
    fV2Pt  = lBosons[1]->pt;
    fV2Eta = lBosons[1]->eta;
    fV2Phi = lBosons[1]->phi;
    fV2M   = lBosons[1]->mass;
    fV2Id  = lBosons[1]->pdgId;
    
  }
  //Should be checked with Gen Met
  fVPt     = lVMass.Pt();
  fVEta    = lVMass.Eta();
  fVPhi    = lVMass.Phi();
  fVM      = lVMass.M();
}
