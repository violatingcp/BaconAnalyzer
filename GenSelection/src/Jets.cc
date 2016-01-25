#include "../include/Jets.hh"
#include "../include/Tools.hh"

Jets::Jets(TTree *iTree,TTree *iOTree) { 
  setupJetTree(iOTree);
  baconhep::TGenJet::Class() ->IgnoreTObjectStreamer();
  fGenJets          = new TClonesArray("baconhep::TGenJet");
  fGenFatJets       = new TClonesArray("baconhep::TGenJet");
  iTree->SetBranchAddress("GenJet"       ,&fGenJets);     fGenJetBr      = iTree->GetBranch("GenJet"      ); 
  iTree->SetBranchAddress("GenFatJet"    ,&fGenFatJets);  fGenFatJetBr   = iTree->GetBranch("GenFatJet"   );   
}
void Jets::setupJetTree(TTree *iTree) { 
  iTree->Branch("jpt_1"    ,&fGJetPt1,  "fGJetPt1/F");
  iTree->Branch("jeta_1"   ,&fGJetEta1, "fGJetEta1/F");
  iTree->Branch("jphi_1"   ,&fGJetPhi1, "fGJetPhi1/F");
  iTree->Branch("jm_1"     ,&fGJetM1,   "fGJetM1/F");
  iTree->Branch("jid_1"    ,&fGJetId1,  "fGJetId1/F");

  iTree->Branch("jpt_2"    ,&fGJetPt2,  "fGJetPt2/F");
  iTree->Branch("jeta_2"   ,&fGJetEta2, "fGJetEta2/F");
  iTree->Branch("jphi_2"   ,&fGJetPhi2, "fGJetPhi2/F");
  iTree->Branch("jm_2"     ,&fGJetM2,   "fGJetM2/F");
  iTree->Branch("jid_2"    ,&fGJetId2,  "fGJetId2/F");

  iTree->Branch("jpt_3"    ,&fGJCPt,     "fGJCPt/F");
  iTree->Branch("jeta_3"   ,&fGJCEta,    "fGJCEta/F");
  iTree->Branch("jphi_3"   ,&fGJCPhi,    "fGJCPhi/F");
  iTree->Branch("jm_3"     ,&fGJCM,      "fGJCM/F");
  iTree->Branch("jid_3"    ,&fGJCId,     "fGJCId/F");

  iTree->Branch("fjpt"     ,&fGFJPt,     "fGFJPt/F");
  iTree->Branch("fjeta"    ,&fGFJEta,    "fGFJEta/F");
  iTree->Branch("fjphi"    ,&fGFJPhi,    "fGFJPhi/F");
  iTree->Branch("fjm"      ,&fGFJM,      "fGFJM/F");
  iTree->Branch("fjmtrim"  ,&fGFJMTrim,  "fGFJMTrim/F");

  iTree->Branch("mjj"      ,&fMJJ,      "fMJJ/F");
  iTree->Branch("ptjj"     ,&fPtJJ,     "fPtJJ/F");
  iTree->Branch("phijj"    ,&fPhiJJ,    "fPhiJJ/F");
  iTree->Branch("etajj"    ,&fEtaJJ,    "fEtaJJ/F");
  iTree->Branch("jdeta"    ,&fJDEta,    "fJDEta/F");
  iTree->Branch("jdphi"    ,&fJDPhi,    "fJDPhi/F");
  iTree->Branch("njets"    ,&fNJets,    "fNJets/F");
  iTree->Branch("njets50"  ,&fNJets50,  "fNJets50/F");
}
void Jets::clearJets(int i0){
  fGenJets   ->Clear(); fGenJetBr   ->GetEntry(i0);  
  fGenFatJets->Clear(); fGenFatJetBr->GetEntry(i0);  
  fGJetPt1=0;
  fGJetEta1=0;
  fGJetPhi1=0;
  fGJetM1=0;
  fGJetId1=0;
  fGJetPt2=0;
  fGJetEta2=0;
  fGJetPhi2=0;
  fGJetId2=0;
  fGJetM2=0;
  fGJCPt=0;
  fGJCEta=0;
  fGJCPhi=0;
  fGJCM=0;
  fGJCId=0;

  fGFJPt=0;
  fGFJEta=0;
  fGFJPhi=0;
  fGFJM=0;
  fGFJMTrim=0;
   
  fMJJ=0;
  fPtJJ=0;
  fPhiJJ=0;
  fEtaJJ=0;
  fJDEta=0;
  fJDPhi=0;
  fNJets=0;
  fNJets50=0;
}
void Jets::insert(baconhep::TGenJet *iJet,std::vector<baconhep::TGenJet*> &lJets) { 
  for(std::vector<baconhep::TGenJet*>::iterator pIter = lJets.begin(); pIter != lJets.end(); pIter++) { 
    if(iJet->pt > (*pIter)->pt) { 
      lJets.insert(pIter,iJet);
      return;
    }
  }
  lJets.push_back(iJet);
}

void Jets::selectJets(std::vector<TLorentzVector> &iVec) { 
  std::vector<baconhep::TGenJet*> lJets;
  for(int i0 = 0; i0 < fGenJets->GetEntriesFast(); i0++) {
    baconhep::TGenJet *pJet = (baconhep::TGenJet*) fGenJets->At(i0);
    bool pMatch = false;
    for(int i1 = 0; i1 < int(iVec.size()); i1++) { 
      if(Tools::deltaR(iVec[i1].Eta(),iVec[i1].Phi(),pJet->eta,pJet->phi) < 0.5) pMatch = true;
    }
    if(pMatch) continue;
    if(pJet->pt > 30) fNJets++;
    if(pJet->pt > 50) fNJets50++;
    if(pJet->pt < 15) continue;
    insert(pJet,lJets);
  }
  if(lJets.size() >  0) fGJetPt1  = lJets[0]->pt;
  if(lJets.size() >  0) fGJetEta1 = lJets[0]->eta;
  if(lJets.size() >  0) fGJetPhi1 = lJets[0]->phi;
  if(lJets.size() >  0) fGJetM1   = lJets[0]->mass;
  if(lJets.size() >  0) fGJetId1  = lJets[0]->pdgId;
  if(lJets.size() >  1) fGJetPt2  = lJets[1]->pt;
  if(lJets.size() >  1) fGJetEta2 = lJets[1]->eta;
  if(lJets.size() >  1) fGJetPhi2 = lJets[1]->phi;
  if(lJets.size() >  1) fGJetM2   = lJets[1]->mass;
  if(lJets.size() >  1) fGJetId2  = lJets[1]->pdgId;
  if(lJets.size() >  2) fGJCPt    = lJets[2]->pt;
  if(lJets.size() >  2) fGJCEta   = lJets[2]->eta;;
  if(lJets.size() >  2) fGJCPhi   = lJets[2]->phi;
  if(lJets.size() >  2) fGJCM     = lJets[2]->mass;
  if(lJets.size() >  2) fGJCId    = lJets[2]->pdgId;

  std::vector<baconhep::TGenJet*> lVJets;
  for(int i0 = 0; i0 < fGenFatJets->GetEntriesFast(); i0++) {
    baconhep::TGenJet *pJet = (baconhep::TGenJet*) fGenFatJets->At(i0);
    bool pMatch = false;
    for(int i1 = 0; i1 < int(iVec.size()); i1++) { 
      if(Tools::deltaR(iVec[i1].Eta(),iVec[i1].Phi(),pJet->eta,pJet->phi) < 0.5) pMatch = true;
    }
    if(pMatch) continue;
    if(pJet->pt < 15) continue;
    insert(pJet,lVJets);
  }
  if(lVJets.size() >  0) fGFJPt    = lVJets[0]->pt;
  if(lVJets.size() >  0) fGFJEta   = lVJets[0]->eta;
  if(lVJets.size() >  0) fGFJPhi   = lVJets[0]->phi;
  if(lVJets.size() >  0) fGFJM     = lVJets[0]->mass;
  if(lVJets.size() >  0) fGFJMTrim = lVJets[0]->mtrim;
  /*
  baconhep::TGenJet *lCentral = 0;
  for(int i0 = 2; i0 < int(lJets.size()); i0++) {
    baconhep::TGenJet *pJet = (baconhep::TGenJet*) lJets[i0];
    if(!((pJet->eta < fJetEta1 && pJet->eta > fJetEta2) || (pJet->eta > fJetEta1 && pJet->eta < fJetEta2))) continue;  
    if(lCentral == 0) lCentral = pJet;
    if(pJet->pt > lCentral->pt) pJet = lCentral;
  }
  if(lCentral != 0) fJCPt  = lCentral->pt;
  if(lCentral != 0) fJCEta = lCentral->eta;;
  if(lCentral != 0) fJCPhi = lCentral->phi;
  if(lCentral != 0) fJCM   = lCentral->mass;
  */
  return;
}
void Jets::computeFinal() { 
  if(fGJetPt1 == 0 || fGJetPt2 == 0) return;
  TLorentzVector lVec1,lVec2;
  lVec1.SetPtEtaPhiM(fGJetPt1,fGJetEta1,fGJetPhi1,fGJetM1);
  lVec2.SetPtEtaPhiM(fGJetPt2,fGJetEta2,fGJetPhi2,fGJetM2);
  fMJJ   = (lVec1+lVec2).M();
  fPtJJ  = (lVec1+lVec2).Pt();
  fPhiJJ = (lVec1+lVec2).Phi();
  fEtaJJ = (lVec1+lVec2).Eta(); 
  fJDEta =  fGJetEta1-fGJetEta2;
  fJDPhi = (fGJetPhi1-fGJetPhi2);
  if(fabs(fJDPhi) > 2.*TMath::Pi()-fabs(fJDPhi) && fJDPhi > 0)  fJDPhi -= 2.*TMath::Pi();
  if(fabs(fJDPhi) > 2.*TMath::Pi()-fabs(fJDPhi) && fJDPhi < 0)  fJDPhi += 2.*TMath::Pi();
}
void Jets::fillJets(std::vector<TLorentzVector> &iVec) { 
  selectJets(iVec);
  computeFinal();
}
