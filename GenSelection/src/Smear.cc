#include "../include/Smear.hh"
#include "../include/JetRes.hh"
#include "../include/Tools.hh"

Smear::Smear(TTree *iTree,TTree *iOTree) { 
  setupTree(iOTree);
  fJetRes = new JetRes(false);
}
void Smear::setupTree(TTree *iTree) { 
  iTree->Branch("jpt_1r"    ,&fGJetPt1,  "fGJetPt1/F");
  iTree->Branch("jeta_1r"   ,&fGJetEta1, "fGJetEta1/F");
  iTree->Branch("jphi_1r"   ,&fGJetPhi1, "fGJetPhi1/F");
  iTree->Branch("jm_1r"     ,&fGJetM1,   "fGJetM1/F");

  iTree->Branch("jpt_2r"    ,&fGJetPt2,  "fGJetPt2/F");
  iTree->Branch("jeta_2r"   ,&fGJetEta2, "fGJetEta2/F");
  iTree->Branch("jphi_2r"   ,&fGJetPhi2, "fGJetPhi2/F");
  iTree->Branch("jm_2r"     ,&fGJetM2,   "fGJetM2/F");
  iTree->Branch("jid_r2"    ,&fGJetId2,  "fGJetId2/F");

  iTree->Branch("jpt_3r"    ,&fGJCPt,     "fGJCPt/F");
  iTree->Branch("jeta_3r"   ,&fGJCEta,    "fGJCEta/F");
  iTree->Branch("jphi_3r"   ,&fGJCPhi,    "fGJCPhi/F");
  iTree->Branch("jm_3r"     ,&fGJCM,      "fGJCM/F");
  iTree->Branch("jid_3r"    ,&fGJCId,     "fGJCId/F");

  iTree->Branch("met"       ,&fMet,       "fMet/F");
  iTree->Branch("metphi"    ,&fMetPhi,    "fMetPhi/F");
  iTree->Branch("htr"       ,&fHT,        "fHT/F");
}
void Smear::clear(int i0){
  fGJetPt1=0;
  fGJetEta1=0;
  fGJetPhi1=0;
  fGJetM1=0;
  fGJetPt2=0;
  fGJetEta2=0;
  fGJetPhi2=0;
  fGJetM2=0;
  fGJCPt=0;
  fGJCEta=0;
  fGJCPhi=0;
  fGJCM=0;
  fHT=0;
  fMet=0;
  fMetPhi=0;
}
void Smear::smear(double iMet,double iMetPhi,double iHT,std::vector<baconhep::TGenJet*> &iJets) { 
  std::vector<TLorentzVector> lVec = fJetRes->smearEvent(iMet,iMetPhi,iHT,iJets);
  if(lVec.size() >  2) fGJetPt1  = lVec[0].Pt();
  if(lVec.size() >  2) fGJetEta1 = lVec[0].Eta();
  if(lVec.size() >  2) fGJetPhi1 = lVec[0].Phi();
  if(lVec.size() >  2) fGJetM1   = lVec[0].M();
  if(lVec.size() >  3) fGJetPt2  = lVec[1].Pt();
  if(lVec.size() >  3) fGJetEta2 = lVec[1].Eta();
  if(lVec.size() >  3) fGJetPhi2 = lVec[1].Phi();
  if(lVec.size() >  3) fGJetM2   = lVec[1].M();
  if(lVec.size() >  4) fGJCPt    = lVec[2].Pt();
  if(lVec.size() >  4) fGJCEta   = lVec[2].Eta();
  if(lVec.size() >  4) fGJCPhi   = lVec[2].Phi();
  if(lVec.size() >  4) fGJCM     = lVec[2].M();
  fMet    = lVec[lVec.size()-2].Pt();
  fMetPhi = lVec[lVec.size()-2].Phi();
  fHT     = lVec[lVec.size()-1].Pt();
  return;
}

