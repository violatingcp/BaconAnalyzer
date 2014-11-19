#include "../include/Tools.hh"
#include "TH1F.h"

Tools::Tools() {}
Tools::~Tools() {}

double Tools::deltaPhi(double iPhi1,double iPhi2) { 
  double lDPhi = fabs(iPhi1-iPhi2); 
  if(lDPhi > 2*TMath::Pi()-lDPhi) lDPhi = 2*TMath::Pi()-lDPhi;
  return lDPhi;
}
double Tools::deltaR  (double iEta1,double iPhi1,double iEta2,double iPhi2) { 
  double lDEta = fabs(iEta1-iEta2);
  double lDPhi = deltaPhi(iPhi1,iPhi2);
  return sqrt(lDEta*lDEta+lDPhi*lDPhi);
}
TTree* Tools::load(std::string iName,TFile *iFile,int &iNMCEvents) { 
  //iFile = new TFile(iName.c_str());
  iFile = TFile::Open(iName.c_str());
  TTree *lTree = (TTree*)  iFile->FindObjectAny("Events");
  TH1F  *lHist = (TH1F*)   iFile->FindObjectAny("TotalEvents");
  //if(lHist == 0) lHist = (TH1F*)  iFile->FindObjectAny("TotalEvants");
  iNMCEvents = lHist->GetEntries();
  return lTree;
  iNMCEvents = 1.;
}
double Tools::mT(double iE1,double iPhi1,double iE2,double iPhi2) { 
  double lEx = iE1*TMath::Cos(iPhi1) + iE2*TMath::Cos(iPhi2);
  double lEy = iE1*TMath::Sin(iPhi1) + iE2*TMath::Sin(iPhi2);
  double lEt = sqrt((iE1+iE2)*(iE1+iE2) - lEx*lEx -lEy*lEy);
  return lEt;
}
double Tools::pZeta(int iId,double iPt1,double iEta1,double iPhi1,double iPt2,double iEta2,double iPhi2,double iMet,double iMetPhi,double iM2) { 
  double lTM = iM2;
  double lMM = 0.1056;
  double lEM = 0.0005;
  double lM1 = lMM;
  double lM2 = lTM;
  if(iId != 3)  lM1 = lTM;
  if(iId == 1)  lM2 = lMM;  
  if(iId == 2 || iId == 3)  lM2 = lEM;  
  TLorentzVector lLep1; lLep1.SetPtEtaPhiM(iPt1    ,iEta1    ,iPhi1     ,lM1);
  TLorentzVector lLep2; lLep2.SetPtEtaPhiM(iPt2    ,iEta2    ,iPhi2     ,lM2);
  double         lPPhi = (iPhi1+iPhi2)/2.; if(fabs(iPhi1-iPhi2) > TMath::Pi()) lPPhi+=TMath::Pi();
  double         lPX   = cos(lPPhi);
  double         lPY   = sin(lPPhi);
  double         lPMet = iMet*cos(iMetPhi)*lPX + iMet*sin(iMetPhi)*lPY;
  double         lPVis = (lLep1+lLep2).Px()*lPX + (lLep1+lLep2).Py()*lPY;
  return (0.85*lPVis-lPMet);
}
double Tools::pVis(int iId,double iPt1,double iEta1,double iPhi1,double iPt2,double iEta2,double iPhi2,double iM2) { 
  double lTM = iM2;
  double lMM = 0.1056;
  double lEM = 0.0005;
  double lM1 = lMM;
  double lM2 = lTM;
  if(iId != 3)  lM1 = lTM;
  if(iId == 1)  lM2 = lMM;  
  if(iId == 2 || iId == 3)  lM2 = lEM;  
  TLorentzVector lLep1; lLep1.SetPtEtaPhiM(iPt1    ,iEta1    ,iPhi1     ,lM1);
  TLorentzVector lLep2; lLep2.SetPtEtaPhiM(iPt2    ,iEta2    ,iPhi2     ,lM2);
  double         lPPhi = (iPhi1+iPhi2)/2.; if(fabs(iPhi1-iPhi2) > TMath::Pi()) lPPhi+=TMath::Pi();
  double         lPX   = cos(lPPhi);
  double         lPY   = sin(lPPhi);
  double         lPVis = (lLep1+lLep2).Px()*lPX + (lLep1+lLep2).Py()*lPY;
  return lPVis;
}
double Tools::pMet(int iId,double iPt1,double iEta1,double iPhi1,double iPt2,double iEta2,double iPhi2,double iMet,double iMetPhi,double iM2) { 
  double lTM = iM2;
  double lMM = 0.1056;
  double lEM = 0.0005;
  double lM1 = lMM;
  double lM2 = lTM;
  if(iId != 3)  lM1 = lTM;
  if(iId == 1)  lM2 = lMM;  
  if(iId == 2 || iId == 3)  lM2 = lEM;  
  TLorentzVector lLep1; lLep1.SetPtEtaPhiM(iPt1    ,iEta1    ,iPhi1     ,lM1);
  TLorentzVector lLep2; lLep2.SetPtEtaPhiM(iPt2    ,iEta2    ,iPhi2     ,lM2);
  double         lPPhi = (iPhi1+iPhi2)/2.; if(fabs(iPhi1-iPhi2) > TMath::Pi()) lPPhi+=TMath::Pi();
  double         lPX   = cos(lPPhi);
  double         lPY   = sin(lPPhi);
  double         lPMet = iMet*cos(iMetPhi)*lPX + iMet*sin(iMetPhi)*lPY;
  return lPMet;
}
TLorentzVector Tools::system(int iId,double iPt1,double iEta1,double iPhi1,double iPt2,double iEta2,double iPhi2,double iMet,double iMetPhi,double iM2) { 
  double lTM = iM2;
  double lMM = 0.1056;
  double lEM = 0.0005;
  double lM1 = lMM;
  double lM2 = lTM;
  if(iId != 3)  lM1 = lTM;
  if(iId == 1)  lM2 = lMM;  
  if(iId == 2 || iId == 3)  lM2 = lEM;  
  TLorentzVector lLep1; lLep1.SetPtEtaPhiM(iPt1    ,iEta1    ,iPhi1     ,lM1);
  TLorentzVector lLep2; lLep2.SetPtEtaPhiM(iPt2    ,iEta2    ,iPhi2     ,lM2);
  TLorentzVector lMet;  lMet .SetPtEtaPhiM(iMet    ,0.       ,iMetPhi   ,0);
  return (lLep1+lLep2+lMet);
}
double  Tools::collinear(int iId,double iPt1,double iEta1,double iPhi1,double iPt2,double iEta2,double iPhi2,double iMet,double iMetPhi,double iM2) { 
  double metx = iMet*cos(iMetPhi);
  double mety = iMet*sin(iMetPhi);

  double lTM = iM2;
  double lMM = 0.1056;
  double lEM = 0.0005;
  double lM1 = lMM;
  double lM2 = lTM;
  if(iId != 3)  lM1 = lTM;
  if(iId == 1)  lM2 = lMM;  
  if(iId == 2 || iId == 3)  lM2 = lEM;  

  TLorentzVector lep1, lep2;
  lep1.SetPtEtaPhiM(iPt1,iEta1,iPhi1,lM1);
  lep2.SetPtEtaPhiM(iPt2,iEta2,iPhi2,lM2);
  double pm1  = (metx*sin(lep2.Phi())-mety*cos(lep2.Phi()))/(sin(lep1.Theta())*sin(lep2.Phi()-lep1.Phi()));
  double pm2  = (metx*sin(lep1.Phi())-mety*cos(lep1.Phi()))/(sin(lep2.Theta())*sin(lep1.Phi()-lep2.Phi()));
  double camass = 0.;
  if(pm1 > 0 && pm2 > 0)
    {
      double x1   = lep1.Vect().Mag() / (lep1.Vect().Mag() + pm1);
      double x2   = lep2.Vect().Mag() / (lep2.Vect().Mag() + pm2);
      camass = (lep1+lep2).M()/sqrt(x1*x2);
    }
  return camass;
}
TLorentzVector getVisMass(int iId,double iPt1,double iEta1,double iPhi1,double iPt2,double iEta2,double iPhi2,double iM2) { 
  double lTM = iM2;
  double lMM = 0.1056;
  double lEM = 0.0005;
  TLorentzVector lTau; 
  TLorentzVector lOth;
  double lM1 = lMM;
  double lM2 = lTM;
  if(iId == 1)  lM1 = lEM;  
  if(iId == 3)  lM2 = lEM;  
  lTau.SetPtEtaPhiM(iPt1    ,iEta1    ,iPhi1     ,lM1);
  lOth.SetPtEtaPhiM(iPt2    ,iEta2    ,iPhi2     ,lM2);
  lTau += lOth;
  return lTau;
}
double u1(double iPt,double iPhi,double iMet,double iMPhi,double iTruePhi) { 
  double lUX    = iMet*cos(iMPhi) + iPt*cos(iPhi);
  double lUY    = iMet*sin(iMPhi) + iPt*sin(iPhi);
  double lU     = sqrt(lUX*lUX + lUY*lUY);
  double lCos   = -(lUX *cos(iTruePhi) + lUY*sin(iTruePhi))/lU;
  //double lSin   =  (-lUY*cos(iPhi) + lUX*sin(iTruePhi))/lU;
  return lU*lCos;
}
double u2(double iPt,double iPhi,double iMet,double iMPhi,double iTruePhi) { 
  double lUX    = iMet*cos(iMPhi) + iPt*cos(iPhi);
  double lUY    = iMet*sin(iMPhi) + iPt*sin(iPhi);
  double lU     = sqrt(lUX*lUX + lUY*lUY);
  //double lCos   = -(lUX *cos(iTruePhi) + lUY*sin(iTruePhi))/lU;
  double lSin   =  (-lUY*cos(iTruePhi) + lUX*sin(iTruePhi))/lU;
  return lU*lSin;
}
double uphi(double iPt,double iPhi,double iMet,double iMPhi){//,double iTruePhi) { 
  double lUX    = iMet*cos(iMPhi) + iPt*cos(iPhi);
  double lUY    = iMet*sin(iMPhi) + iPt*sin(iPhi);
  return atan2(lUY,lUX);
  ///TVector3 lU(lUX,lUY,0);
  //return lU.Phi();
}
/*
TLorentzVector mass(TMuon *iMuon1,TMuon *iMuon2) { 
  TLorentzVector lL1(0,0,0,0); 
  TLorentzVector lL2(0,0,0,0);
  lL1.SetPtEtaPhiM(iMuon1->pt,iMuon1->eta,iMuon1->phi,0.1054);
  lL2.SetPtEtaPhiM(iMuon2->pt,iMuon2->eta,iMuon2->phi,0.1054);
  lL1 = lL1 + lL2;
  return lL1;  
}
*/

/*
void loadZMuMu(TMuon *iMuon1,TMuon *iMuon2) { 
  for(int i0 = 0; i0 < fMuons->GetEntriesFast(); i0++) { 
    TMuon *p0Muon = (TMuon*) fMuons->At(i0);
    for(int i1 = 0; i1 < i0; i1++) { 
      TMuon *p1Muon = (TMuon*) fMuons->At(i1);
      TLorentzVector lM = mass(p0Muon,p1Muon);
      if(lM.M() > 60 && lM.M() < 120) { 
	fZM   = lM.M();
	fZPt  = lM.Pt();
	fZY   = lM.Rapidity();
	fZPz  = lM.Pz();
	fZPhi = lM.Phi();
	iMuon1 = p0Muon;
	iMuon2 = p1Muon;
	return;
      }
    }
  }
}
*/
