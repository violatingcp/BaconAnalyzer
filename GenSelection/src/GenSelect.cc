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
  iTree->Branch("q2" ,&fQ2  ,"fQ2/F");
  iTree->Branch("evtweight" ,&fEvtWeight  ,"fEvtWeight/F");
  iTree->Branch("mcweight" ,&fMCWeight  ,"fMCWeight/F");
  iTree->Branch("mcweight2",&fMCWeight2 ,"fMCWeight2/F");
  iTree->Branch("xs"       ,&fXS        ,"fXS/F");
  iTree->Branch("xs2"      ,&fXS2       ,"fXS2/F");

  iTree->Branch("med_pt"    ,&fMPt,   "fMPt/F");
  iTree->Branch("med_y"     ,&fMEta,  "fMEta/F");
  iTree->Branch("med_phi"   ,&fMPhi,  "fMPhi/F");
  iTree->Branch("med_m"     ,&fMM,    "fMM/F");
  iTree->Branch("med_id"    ,&fMId,   "fMId/I");
  iTree->Branch("med_status",&fMStatus,"fMStatus/I");

  iTree->Branch("med2_pt"    ,&fM2Pt,   "fM2Pt/F");
  iTree->Branch("med2_y"     ,&fM2Eta,  "fM2Eta/F");
  iTree->Branch("med2_phi"   ,&fM2Phi,  "fM2Phi/F");
  iTree->Branch("med2_m"     ,&fM2M,    "fM2M/F");
  iTree->Branch("med2_id"    ,&fM2Id,   "fM2Id/I");
  iTree->Branch("med2_status"    ,&fM2Status,   "fM2Status/I");

  iTree->Branch("v_pt"    ,&fVPt,    "fVPt/F");
  iTree->Branch("v_y"     ,&fVEta,   "fVEta/F");
  iTree->Branch("v_phi"   ,&fVPhi,   "fVPhi/F");
  iTree->Branch("v_m"     ,&fVM,     "fVM/F");
  iTree->Branch("v_id"    ,&fVId,    "fVId/I");
  iTree->Branch("v_status",&fVStatus,"fVStatus/I");
  iTree->Branch("v_iso04" ,&fVIso,   "fVId/F");
  iTree->Branch("v_dyniso",&fDVIso,  "fDVIso/F");

  iTree->Branch("vme_pt"    ,&fVMEPt,    "fVMEPt/F");
  iTree->Branch("vme_y"     ,&fVMEEta,   "fVMEEta/F");
  iTree->Branch("vme_phi"   ,&fVMEPhi,   "fVMEPhi/F");
  iTree->Branch("vme_m"     ,&fVMEM,     "fVMEM/F");
  iTree->Branch("vme_id"    ,&fVMEId,    "fVMEId/I");
  iTree->Branch("vme_status",&fVMEStatus,"fVMEStatus/I");
  iTree->Branch("vme_iso04" ,&fVMEIso,   "fVMEIso/F");
  iTree->Branch("vme_dyniso",&fDVMEIso,  "fDVMEIso/F");

  iTree->Branch("pt_1"    ,&fPt1,   "fPt1/F");
  iTree->Branch("eta_1"   ,&fEta1,  "fEta1/F");
  iTree->Branch("phi_1"   ,&fPhi1,  "fPhi1/F");
  iTree->Branch("m_1"     ,&fM1,    "fM1/F");
  iTree->Branch("cosphi_1",&fCosPhi1,   "fCosPhi1/F");
  iTree->Branch("id_1"    ,&fId1,   "fId1/I");
  iTree->Branch("id_1status",&fStatusId1,   "fId1/I");

  iTree->Branch("pt_2"    ,&fPt2,   "fPt2/F");
  iTree->Branch("eta_2"   ,&fEta2,  "fEta2/F");
  iTree->Branch("phi_2"   ,&fPhi2,  "fPhi2/F");
  iTree->Branch("m_2"     ,&fM2,    "fM2/F");
  iTree->Branch("cosphi_2",&fCosPhi2,   "fCosPhi2/F");
  iTree->Branch("id_2"    ,&fId2,   "fId2/I");
  iTree->Branch("id_2status"    ,&fStatusId2,   "fStatusId2/I");

  iTree->Branch("toppt"    ,&fTPt,   "fTPt/F");
  iTree->Branch("topeta"   ,&fTEta,  "fTEta/F");
  iTree->Branch("topphi"   ,&fTPhi,  "fTPhi/F");
  iTree->Branch("topm"     ,&fTM,    "fTM/F");

  iTree->Branch("top1pt"    ,&fT1Pt,   "fT1Pt/F");
  iTree->Branch("top1eta"   ,&fT1Eta,  "fT1Eta/F");
  iTree->Branch("top1phi"   ,&fT1Phi,  "fT1Phi/F");
  iTree->Branch("top1m"     ,&fT1M,    "fT1M/F");

  iTree->Branch("npartons",&fNPartons,"fNPartons/I");
  iTree->Branch("genmet"     ,&fMet,     "fMet/F");
  iTree->Branch("genmetphi"  ,&fMetPhi,  "fMetPhi/F");
}
void GenSelect::clearGen(int i0){
  fGenBr->GetEntry(i0);
  fGenParts->Clear(); fGenPartBr  ->GetEntry(i0);  
  fQ2 = -1;
  fDecay = -1;
  fEvtWeight    = 1.;
  fMPt       = 0;
  fMEta      = 0;
  fMPhi      = 0;
  fMM        = 0;
  fMId       = 0; 
  fMStatus   = 0; 

  fM2Pt       = 0;
  fM2Eta      = 0;
  fM2Phi      = 0;
  fM2M        = 0;
  fM2Id       = 0; 
  fM2Status   = 0; 

  fVPt       = 0;
  fVEta      = 0;
  fVPhi      = 0;
  fVM        = 0;
  fVId       = 0; 
  fVStatus   = 0; 
  fVIso      = 0;  
  fDVIso     = 0;  

  fVMEPt       = 0;
  fVMEEta      = 0;
  fVMEPhi      = 0;
  fVMEM        = 0;
  fVMEId       = 0; 
  fVMEStatus   = 0; 
  fVMEIso      = 0;
  fDVMEIso     = 0;  

  fPt1       = 0;
  fEta1      = 0;
  fPhi1      = 0;
  fM1        = 0;
  fCosPhi1   = 0; 
  fId1       = 0;
  fStatusId1 = 0;
 
  fPt2       = 0;
  fEta2      = 0;
  fPhi2      = 0;
  fM2        = 0;
  fCosPhi2   = 0; 
  fId2       = 0;
  fStatusId2 = 0;

  fTPt       = 0;
  fTEta      = 0;
  fTPhi      = 0;
  fTM        = 0;

  fNPartons  = 0; 
  fMet       = 0;
  fMetPhi    = 0; 
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
void GenSelect::add(baconhep::TGenParticle *iGen,std::vector<baconhep::TGenParticle*> &iBoson) { 
  for(std::vector<baconhep::TGenParticle*>::iterator pIter = iBoson.begin();  pIter != iBoson.end(); pIter++) { 
    if((*pIter)->pt > iGen->pt) continue;
    iBoson.insert(pIter,iGen);
    return;
  }
  iBoson.push_back(iGen);
}
void GenSelect::findBosons(std::vector<baconhep::TGenParticle*> & iBoson,int iIdMin,int iIdMax,int iStatusMin,int iStatusMax,bool iLast) { 
  baconhep::TGenParticle *lGen = 0;
  for(int i0 = 0; i0 < fGenParts->GetEntriesFast(); i0++) { 
    baconhep::TGenParticle *pGen = (baconhep::TGenParticle*) fGenParts->At(i0);
    //if(pGen->pt < 15) continue; 
    if(pGen->pt < 100 && (iIdMin > 20 && iIdMin < 25)) continue; 
    if(fabs(pGen->pdgId)  <  iIdMin     ) continue;
    if(fabs(pGen->pdgId)  >  iIdMax     ) continue;
    if(fabs(pGen->status) <  iStatusMin ) continue;
    if(fabs(pGen->status) >  iStatusMax ) continue;
    if(pGen->pt == 0) continue;
    int pBoson = -1;
    if(iStatusMax != 2) { 
      for(unsigned int i1 = 0; i1 < iBoson.size(); i1++) { 
	if(Tools::deltaR(pGen->eta,pGen->phi,iBoson[i1]->eta,iBoson[i1]->phi) < 1.0 ) {
	  if(pGen->status >= iBoson[i1]->status ) pBoson = i1; 
	  continue;
	}
	if(pGen->pt > iBoson[i1]->pt && (pBoson == -1 && iBoson.size() > 1) && i1 > 0) pBoson = 10 + i1;
      }
    }
    //if(pBoson <  10 && pBoson > -1) replace(pBoson,   pGen,iBoson,true);
    //if(pBoson >   9) replace(pBoson-10,pGen,iBoson,false);
    //if(pBoson == -1) add(pGen,iBoson);
    //if(iBoson.size() > 0 && !iLast) break;
    if(lGen == 0) lGen = pGen;
    if(lGen->pt < pGen->pt) lGen = pGen;
  }
  if(lGen != 0) add(lGen,iBoson);
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
float GenSelect::iso(double iEta,double iPhi,int iIdMin,int iIdMax,double iStatusMin,double iStatusMax,double iDRMax) { 
  double lPtIso = 0; 
  //if(iStatusMin == 1) std::cout << "================================================================> " << std::endl;
  for(int i0 = 0; i0 < fGenParts->GetEntriesFast(); i0++) { 
    baconhep::TGenParticle *pGen = (baconhep::TGenParticle*) fGenParts->At(i0);
    if(double(fabs(pGen->status)) <  iStatusMin ) continue;
    if(double(fabs(pGen->status)) >  iStatusMax ) continue;
    if(fabs(pGen->pdgId)          <  iIdMin     ) continue;
    if(fabs(pGen->pdgId)          >  iIdMax     ) continue;
    double pDR = Tools::deltaR(pGen->eta,pGen->phi,iEta,iPhi);
    if(pDR < 0.005  ) continue;
    if(pDR > iDRMax ) continue;
    lPtIso += (pGen->pt)*(1./((1-cos(pDR))/(1-cos(iDRMax))));
    //if(iStatusMin == 1 && fabs(pGen->status) < 40) std::cout << "----> " << pGen->pdgId << " -- " << pGen->status << " -- " << pGen->pt << " -- " << pGen->eta << std::endl;
  }
  //if(iStatusMin == 1) std::cout << "================================================================> " << std::endl;
  return lPtIso;
}
void GenSelect::met() { 
  fMet = 0; 
  fMetPhi = 0; 
  TLorentzVector lVec; lVec.SetPtEtaPhiM(0,0,0,0);
  for(int i0 = 0; i0 < fGenParts->GetEntriesFast(); i0++) { 
    baconhep::TGenParticle *pGen = (baconhep::TGenParticle*) fGenParts->At(i0);
    if(fabs(pGen->pdgId) < 11 || (fabs(pGen->pdgId) > 19 && fabs(pGen->pdgId) < 100000)) continue;
    if( abs(pGen->pdgId) % 2 == 1) continue;
    if(fabs(pGen->status) != 1) continue;
    TLorentzVector pVec; pVec.SetPtEtaPhiM(pGen->pt,0,pGen->phi,0);
    lVec += pVec;
  }
  fMet    = lVec.Pt();
  fMetPhi = lVec.Phi();
}
void GenSelect::fillGen(double iTot,std::vector<TLorentzVector> &iVec) {
  fXS2       = fGen->xs;
  fMCWeight2 = fXS2/iTot;
  fEvtWeight = fGen->weight;
  fQ2        = fGen->scalePDF;
  std::vector<baconhep::TGenParticle*> lMeds;
  std::vector<baconhep::TGenParticle*> lBosons;
  std::vector<baconhep::TGenParticle*> lDaughters;
  std::vector<baconhep::TGenParticle*> lTops;
  std::vector<baconhep::TGenParticle*> lFirstTops;
  std::vector<baconhep::TGenParticle*> lMEBosons;
  findBosons(lMeds,55,55,-1,70);
  findBosons(lMeds,54,54,-1,70);
  findBosons(lMeds,25,25,-1,70);
  findBosons(lMeds,32,32,-1,70);
  findBosons(lMeds,9900032,9900032,-1,70);
  findBosons(lMeds,22,22,-1,70);
  findBosons(lBosons,32,32,-1,70);
  findBosons(lBosons,25,25,-1,70);
  findBosons(lBosons,9900032,9900032,-1,70);
  //findBosons(lBosons,9100000,9100000,20,70);
  findBosons(lBosons,23,23,20,70);
  findBosons(lBosons,22,22,20,70);
  findBosons(lBosons,24,24,20,70);
  if(lBosons.size() == 0) findBosons(lBosons,23,23,-1,70);
  if(lBosons.size() == 0) findBosons(lBosons,24,24,-1,70);
  if(lBosons.size() == 0) findBosons(lBosons,22,22,-1,70);
  findBosons(lMEBosons,23,23,20,30);
  findBosons(lMEBosons,22,22,20,30);
  findBosons(lMEBosons,24,24,20,30);
  //if(lBosons.size() == 0) findBosons(lBosons,25,25,-1,70);
  //if(lBosons.size() == 0) findBosons(lBosons,32,32,50,70);
  findBosons(lDaughters,1000022,1000022,-1,70);
  findBosons(lDaughters,91000022,91000022,-1,70);
  findBosons(lDaughters,12,12,0,2);
  findBosons(lDaughters,14,14,0,2);
  findBosons(lDaughters,16,16,0,2);
  findBosons(lDaughters,11,11,0,2);
  findBosons(lDaughters,13,13,0,2);
  findBosons(lDaughters,15,15,0,2);
  if(lMeds.size() > 0) { 
    fMPt  = lMeds[0]->pt; 
    fMEta = lMeds[0]->eta; 
    fMPhi = lMeds[0]->phi; 
    fMM   = lMeds[0]->mass; 
    fMId  = lMeds[0]->pdgId;
    fMStatus  = lMeds[0]->status;
  }
  if(lMeds.size() > 1) { 
    fM2Pt  = lMeds[1]->pt; 
    fM2Eta = lMeds[1]->eta; 
    fM2Phi = lMeds[1]->phi; 
    fM2M   = lMeds[1]->mass; 
    fM2Id  = lMeds[1]->pdgId;
    fM2Status  = lMeds[1]->status;
  }
  TLorentzVector lLep1,lLep2;
  if(lDaughters.size() > 0) { 
    fPt1  = lDaughters[0]->pt; 
    fEta1 = lDaughters[0]->eta; 
    fPhi1 = lDaughters[0]->phi; 
    fM1   = lDaughters[0]->mass; 
    fId1  = lDaughters[0]->pdgId;
    fStatusId1  = lDaughters[0]->status;
    lLep2.SetPtEtaPhiM(fPt1,fEta1,fPhi1,fM1);
  }
  if(lDaughters.size() > 1) { 
    fPt2  = lDaughters[1]->pt; 
    fEta2 = lDaughters[1]->eta; 
    fPhi2 = lDaughters[1]->phi; 
    fM2   = lDaughters[1]->mass;
    fId2  = lDaughters[1]->pdgId;
    fStatusId2  = lDaughters[1]->status;
    lLep2.SetPtEtaPhiM(fPt2,fEta2,fPhi2,fM2);
  }
  if(abs(fId2) %2 == 0 || abs(fId1) % 2 == 0) fDecay = 1.;
  if(abs(fId2) %2 == 0 && abs(fId1) % 2 == 0) fDecay = 2.;
  findBosons(lTops,6,6,-1,70);
  findBosons(lFirstTops,6,6,-1,70,false);
  bool lSkip = false; //if(lTops.size() > 0 && fabs(lBosons[lBosons.size()-1]->pdgId) == 24) lSkip=true;
  TLorentzVector lVMass;
  for(int i0 = 0; i0 < TMath::Min(int(lBosons.size()),1); i0++) { 
    //if(lBosons.size() > 1 && lBosons[i0]->status < 60) continue;
    if(lBosons[i0]->pdgId != 32 && lBosons[i0]->pdgId != 23 && fabs(lBosons[i0]->pdgId) != 24 && fabs(lBosons[i0]->pdgId) != 22 && lBosons[i0]->pdgId != 9900032) continue;
    TLorentzVector pVec;
    pVec.SetPtEtaPhiM(lBosons[i0]->pt,lBosons[i0]->eta,lBosons[i0]->phi,lBosons[i0]->mass);
    if(fId1 < 0 && lBosons[i0]->pdgId > 0) { lLep1.Boost(-pVec.BoostVector()); fCosPhi1 = cos(lLep1.Angle(pVec.Vect())); }
    if(fId1 > 0 && lBosons[i0]->pdgId < 0) { lLep1.Boost(-pVec.BoostVector()); fCosPhi1 = cos(lLep1.Angle(pVec.Vect())); }
    if(fId2 < 0 && lBosons[i0]->pdgId > 0) { lLep2.Boost(-pVec.BoostVector()); fCosPhi2 = cos(lLep2.Angle(pVec.Vect())); }
    if(fId2 > 0 && lBosons[i0]->pdgId < 0) { lLep2.Boost(-pVec.BoostVector()); fCosPhi2 = cos(lLep2.Angle(pVec.Vect())); }
    fVId      = lBosons[i0]->pdgId;
    fVStatus  = lBosons[i0]->status;
    lVMass = pVec;
    if(lSkip) break;
  }
  //Should be checked with Gen Met
  fVPt     = lVMass.Pt();
  fVEta    = lVMass.Rapidity();
  fVPhi    = lVMass.Phi();
  fVM      = lVMass.M();
  fVIso    = iso(lVMass.Eta(),  fVPhi  , -1,10,0, 2,0.4);
  fDVIso   = iso(lVMass.Eta(),  fVPhi  , -1,10,0, 2,TMath::Min(91.18/fVPt/sqrt(0.1),1.0));  //Dynamic cone
  if(lMEBosons.size() > 0) { 
    fVMEPt     = lMEBosons[0]->pt;
    fVMEEta    = lMEBosons[0]->eta;
    fVMEPhi    = lMEBosons[0]->phi;
    fVMEM      = lMEBosons[0]->mass;
    fVMEStatus = lMEBosons[0]->status;
    fVMEId     = lMEBosons[0]->pdgId;
    fVMEIso  = iso(fVMEEta,fVMEPhi,-1,10,20,30,0.4);
    fDVMEIso = iso(fVMEEta,fVMEPhi,-1,10,20,30,TMath::Min(91.18/fVPt/sqrt(0.1),1.0));  //Dynamic cone
  }
  if(lTops.size() > 0) { 
    fTPt    = lTops[0]->pt;
    fTEta   = lTops[0]->eta;
    fTPhi   = lTops[0]->phi;
    fTM     = lTops[0]->mass;
  }
  if(lFirstTops.size() > 0) { 
    fT1Pt    = lFirstTops[0]->pt;
    fT1Eta   = lFirstTops[0]->eta;
    fT1Phi   = lFirstTops[0]->phi;
    fT1M     = lFirstTops[0]->mass;
  }
  npartons();
  for(unsigned int i0 = 0; i0 < lDaughters.size(); i0++) { 
    if(fabs(lDaughters[i0]->pdgId) != 11 && fabs(lDaughters[i0]->pdgId) !=13 && fabs(lDaughters[i0]->pdgId) !=15) continue;
    TLorentzVector pVec(0,0,0,0);
    pVec.SetPtEtaPhiM(lDaughters[i0]->pt,lDaughters[i0]->eta,lDaughters[i0]->phi,lDaughters[i0]->mass);
    iVec.push_back(pVec);
  } 
  if(fVId == 22)   iVec.push_back(lVMass);
  met();
}
void GenSelect::getTop(std::vector<baconhep::TGenParticle*> &iParts) {
  std::vector<int> lTops;
  for(int i0 = 0; i0 < fGenParts->GetEntriesFast(); i0++) { 
    baconhep::TGenParticle *pGen = (baconhep::TGenParticle*) fGenParts->At(i0);
    if(fabs(pGen->pdgId) == 6) lTops.push_back(i0);
    for(unsigned int i1 = 0; i1 < lTops.size(); i1++) { 
      if(pGen->parent == lTops[i1]) iParts.push_back(pGen);
      if(pGen->parent == lTops[i1] && fabs(pGen->pdgId) == 24) {lTops.push_back(i0); continue;}
    }
  }
  /*
  std::cout << "=================================" << std::endl;
  for(unsigned int i0 = 0; i0 < iParts.size(); i0++) { 
    std::cout << "===> " << iParts[i0]->pdgId << " -- " << iParts[i0]->eta << " -- " << iParts[i0]->pt << " -- " << iParts[i0]->status << std::endl;
  }
  */
}
