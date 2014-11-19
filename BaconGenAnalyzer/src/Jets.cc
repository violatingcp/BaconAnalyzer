#include "../include/Jets.hh"
#include "../include/Tools.hh"

Jets::Jets(TTree *iOTree) { 
  setupJetTree(iOTree);
  fRandom = new TRandom2(0xDEADBEEF);
}
void Jets::setupJetTree(TTree *iTree) { 
  iTree->Branch("jpt_1"    ,&fJetPt1,  "fJetPt1/F");
  iTree->Branch("jeta_1"   ,&fJetEta1, "fJetEta1/F");
  iTree->Branch("jphi_1"   ,&fJetPhi1, "fJetPhi1/F");
  iTree->Branch("jm_1"     ,&fJetM1,   "fJetM1/F");
  iTree->Branch("jid_1"    ,&fJetId1,  "fJetId1/I");

  iTree->Branch("jpt_2"    ,&fJetPt2,  "fJetPt2/F");
  iTree->Branch("jeta_2"   ,&fJetEta2, "fJetEta2/F");
  iTree->Branch("jphi_2"   ,&fJetPhi2, "fJetPhi2/F");
  iTree->Branch("jm_2"     ,&fJetM2,   "fJetM2/F");
  iTree->Branch("jid_2"    ,&fJetId2,  "fJetId2/I");

  iTree->Branch("jcpt"    ,&fJCPt,     "fJCPt/F");
  iTree->Branch("jceta"   ,&fJCEta,    "fJCEta/F");
  iTree->Branch("jcphi"   ,&fJCPhi,    "fJCPhi/F");
  iTree->Branch("jcm"     ,&fJCM,      "fJCM/F");
  iTree->Branch("jcid"    ,&fJCId,     "fJCId/I");


  iTree->Branch("jgpt_1"    ,&fGJetPt1,  "fJetPt1/F");
  iTree->Branch("jgeta_1"   ,&fGJetEta1, "fJetEta1/F");
  iTree->Branch("jgphi_1"   ,&fGJetPhi1, "fJetPhi1/F");
  iTree->Branch("jgm_1"     ,&fGJetM1,   "fJetM1/F");

  iTree->Branch("jgpt_2"    ,&fGJetPt2,  "fJetPt2/F");
  iTree->Branch("jgeta_2"   ,&fGJetEta2, "fJetEta2/F");
  iTree->Branch("jgphi_2"   ,&fGJetPhi2, "fJetPhi2/F");
  iTree->Branch("jgm_2"     ,&fGJetM2,   "fJetM2/F");

  iTree->Branch("jgcpt"    ,&fGJCPt,     "fJCPt/F");
  iTree->Branch("jgceta"   ,&fGJCEta,    "fJCEta/F");
  iTree->Branch("jgcphi"   ,&fGJCPhi,    "fJCPhi/F");
  iTree->Branch("jgcm"     ,&fGJCM,      "fJCM/F");

  iTree->Branch("mjj"     ,&fMJJ,      "fMJJ/F");
  iTree->Branch("jdeta"   ,&fJDEta,    "fJDEta/F");
  iTree->Branch("jdphi"   ,&fJDPhi,    "fJDPhi/F");
  iTree->Branch("njets"   ,&fNJets,    "fNJets/F");
  
  //Smear in Pt
  fPtSmear = new TF1("JetPt","sqrt(((TMath::Sign(1,[0])*sq([0]/x))+(sq([1])*(x^([3]-1))))+sq([2]))");
  fPtPar0 = new double[10];
  fPtPar1 = new double[10];
  fPtPar2 = new double[10];
  fPtPar3 = new double[10];
  //9 Eta Bins
  fPtPar0[0] =  1.41584;  fPtPar1[0] = 0.209477; fPtPar2[0] = 0; fPtPar3[0] =  0.588872;
  fPtPar0[1] =  1.41584;  fPtPar1[1] = 0.209477; fPtPar2[1] = 0; fPtPar3[1] =  0.588872;
  fPtPar0[2] =  1.65966;  fPtPar1[2] = 0.223683; fPtPar2[2] = 0; fPtPar3[2] =  0.60873;
  fPtPar0[3] =  2.81978;  fPtPar1[3] = 0.272373; fPtPar2[3] = 0; fPtPar3[3] =  0.579396;
  fPtPar0[4] =  2.56933;  fPtPar1[4] = 0.305802; fPtPar2[4] = 0; fPtPar3[4] =  0.398929;
  fPtPar0[5] =  1.04792;  fPtPar1[5] = 0.466763; fPtPar2[5] = 0; fPtPar3[5] =  0.193137;
  fPtPar0[6] = -1.12329;  fPtPar1[6] = 0.657891; fPtPar2[6] = 0; fPtPar3[6] =  0.139595;
  fPtPar0[7] = -0.561649; fPtPar1[7] = 0.420293; fPtPar2[7] = 0; fPtPar3[7] =  0.392398;
  fPtPar0[8] = -0.499735; fPtPar1[8] = 0.336391; fPtPar2[8] = 0; fPtPar3[8] =  0.430689;
  fPtPar0[9] = -0.349206; fPtPar1[9] = 0.297831; fPtPar2[9] = 0; fPtPar3[9] =  0.471121;
  
  //Smear in Eta
  fEtaSmear = new TF1("JetEta","(sqrt((sq([0]/x)+(sq([1])/x))+sq([2]))+([3]/x))+(([4]/x)/sqrt(x))");
  fEtaPar0 = new double[10];
  fEtaPar1 = new double[10];
  fEtaPar2 = new double[10];
  fEtaPar3 = new double[10];
  fEtaPar4 = new double[10];
  //9 bins again
  fEtaPar0[0] = 0.710478; fEtaPar1[0] = 0.143847;     fEtaPar2[0] = -2.52411e-05; fEtaPar3[0] = 0.77394;  fEtaPar4[0] = -1.89622;
  fEtaPar0[1] = 527.518;  fEtaPar1[1] = -0.0143625;   fEtaPar2[1] = 0.316441;     fEtaPar3[1] =-526.599;  fEtaPar4[1] =  0.244142;
  fEtaPar0[2] = 0.494977; fEtaPar1[2] =  1.48277e-06; fEtaPar2[2] = 0.0171135;    fEtaPar3[2] = 0.354901; fEtaPar4[2] =  0.554065;
  fEtaPar0[3] =  2.88983; fEtaPar1[3] = -1.71576e-06; fEtaPar2[3] = 0.0402796;    fEtaPar3[3] = -1.72649; fEtaPar4[3] = -0.124126;
  fEtaPar0[4] = 0.851656; fEtaPar1[4] =  4.56968e-07; fEtaPar2[4] = 0.0441276;    fEtaPar3[4] =-0.101993; fEtaPar4[4] =  0.773812;
  fEtaPar0[5] =  9.64435; fEtaPar1[5] =  0.458594 ;   fEtaPar2[5] = 1.92485e-08;  fEtaPar3[5] = -8.67043; fEtaPar4[5] = -0.0541106;
  fEtaPar0[6] = 0.343262; fEtaPar1[6] = -3.39452e-07; fEtaPar2[6] = 0.00849674;   fEtaPar3[6] =  1.05358; fEtaPar4[6] = -1.24367;
  fEtaPar0[7] = 0.572596; fEtaPar1[7] = -1.09687e-07; fEtaPar2[7] = 0.0094876;    fEtaPar3[7] = 0.799819; fEtaPar4[7] = -1.23444;
  fEtaPar0[8] = 0.622981; fEtaPar1[8] =  0.100943;    fEtaPar2[8] = 0.00744374;   fEtaPar3[8] = 0.317454; fEtaPar4[8] = -0.324557;
  fEtaPar0[9] =  355.708; fEtaPar1[9] =  2.20794;     fEtaPar2[9] = 0.032666;     fEtaPar3[9] = -354.691; fEtaPar4[9] = -0.857295;

  //No I'm going to kill myself
  fPhiSmear = new TF1("JetPhi","(sqrt((sq([0]/x)+(sq([1])/x))+sq([2]))+([3]/x))+(([4]/x)/sqrt(x))");
  fPhiPar0  = new double[10];
  fPhiPar1  = new double[10];
  fPhiPar2  = new double[10];
  fPhiPar3  = new double[10];
  fPhiPar4  = new double[10];
  //10 bins again
  fPhiPar0[0] = 259.189 ;    fPhiPar1[0] =  0.00132792;  fPhiPar2[0] = -0.311411;     fPhiPar3[0] = -258.647; fPhiPar4[0] =   0;
  fPhiPar0[1] = 0.765787;    fPhiPar1[1] = -3.90638e-06; fPhiPar2[1] = -4.70224e-08;  fPhiPar3[1] = 0.11831;  fPhiPar4[1] =  -1.4675;
  fPhiPar0[2] = 2.11446;     fPhiPar1[2] =  0.203329;    fPhiPar2[2] = -0.0175832;    fPhiPar3[2] = -1.67946; fPhiPar4[2] =  -0.00853474;
  fPhiPar0[3] = 1.9027;      fPhiPar1[3] = -4.56986e-06; fPhiPar2[3] = 0.0304793;     fPhiPar3[3] = -1.09124; fPhiPar4[3] =  -0.136402;
  fPhiPar0[4] = 11.1957;     fPhiPar1[4] =  0.643236;    fPhiPar2[4] = 0.00711422;    fPhiPar3[4] = -10.7613; fPhiPar4[4] =   0.280927;
  fPhiPar0[5] = 0.00336639;  fPhiPar1[5] =  0.0880209;   fPhiPar2[5] = -0.0023084;    fPhiPar3[5] = 0.214304; fPhiPar4[5] =  -0.416353;
  fPhiPar0[6] = 2.92001e-07; fPhiPar1[6] =  0.0718389;   fPhiPar2[6] = -0.00385579;   fPhiPar3[6] = 0.403668; fPhiPar4[6] =  -0.62698;
  fPhiPar0[7] = 0.38469 ;    fPhiPar1[7] =  0.0755727;   fPhiPar2[7] = -0.0044353;    fPhiPar3[7] = 0.453887; fPhiPar4[7] =  -1.8947;
  fPhiPar0[8] = 3.32512e-06; fPhiPar1[8] =  0.063941;    fPhiPar2[8] = -0.00387593;   fPhiPar3[8] = 0.301932; fPhiPar4[8] =  -0.825352;
  fPhiPar0[9] = 926.978;     fPhiPar1[9] =  2.52747;     fPhiPar2[9] = 0.0304001;     fPhiPar3[9] = -926.224; fPhiPar4[9] =  -1.94117;
}
void Jets::clearJets(){
  fJetPt1=0;
  fJetEta1=0;
  fJetPhi1=0;
  fJetM1=0;
  fJetId1=0;
  
  fJetPt2=0;
  fJetEta2=0;
  fJetPhi2=0;
  fJetM2=0;
  fJetId2=0;
  
  fJCPt=0;
  fJCEta=0;
  fJCPhi=0;
  fJCM=0;
  fJCId=0;

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
   
  fMJJ=0;
  fJDEta=0;
  fJDPhi=0;
  fNJets=0;
  fLeptons.clear();
}
void Jets::buildLeptonVec(Gen &iGen) { 
  if(iGen.fPt1 > 0) {
    TLorentzVector lVec; lVec.SetPtEtaPhiM(iGen.fPt1,iGen.fEta1,iGen.fPhi1,iGen.fM1);
    fLeptons.push_back(lVec);
  }
  if(iGen.fPt2 > 0) {
    TLorentzVector lVec; lVec.SetPtEtaPhiM(iGen.fPt1,iGen.fEta1,iGen.fPhi1,iGen.fM1);
    fLeptons.push_back(lVec);
  }
  if(iGen.fPt3 > 0) {
    TLorentzVector lVec; lVec.SetPtEtaPhiM(iGen.fPt1,iGen.fEta1,iGen.fPhi1,iGen.fM1);
    fLeptons.push_back(lVec);
  }
  if(iGen.fPt4 > 0) {
    TLorentzVector lVec; lVec.SetPtEtaPhiM(iGen.fPt1,iGen.fEta1,iGen.fPhi1,iGen.fM1);
    fLeptons.push_back(lVec);
  }
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
void Jets::selectJets(Gen &iGen) { 
  clearJets();
  buildLeptonVec(iGen);
  std::vector<baconhep::TGenJet*> lJets;
  for(int i0 = 0; i0 < iGen.fGenJets->GetEntriesFast(); i0++) {
    baconhep::TGenJet *pJet = (baconhep::TGenJet*) iGen.fGenJets->At(i0);
    bool pMatch = false;
    for(int i1 = 0; i1 < int(fLeptons.size()); i1++) { 
      if(Tools::deltaR(fLeptons[i1].Eta(),fLeptons[i1].Phi(),pJet->eta,pJet->phi) < 0.5) pMatch = true;
    }
    if(pMatch) continue;
    if(pJet->pt > 30) fNJets++;
    if(pJet->pt < 15) continue;
    insert(pJet,lJets);
  }
  if(lJets.size() >  0) fJetPt1  = lJets[0]->pt;
  if(lJets.size() >  0) fJetEta1 = lJets[0]->eta;
  if(lJets.size() >  0) fJetPhi1 = lJets[0]->phi;
  if(lJets.size() >  0) fJetM1   = lJets[0]->mass;
  if(lJets.size() >  1) fJetPt2  = lJets[1]->pt;
  if(lJets.size() >  1) fJetEta2 = lJets[1]->eta;
  if(lJets.size() >  1) fJetPhi2 = lJets[1]->phi;
  if(lJets.size() >  1) fJetM2   = lJets[1]->mass;
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
  return;
}
void Jets::selectJets(GenJet &iGen) { 
  clearJets();
  std::vector<baconhep::TGenJet*> lJets;
  for(int i0 = 0; i0 < iGen.fGenJets->GetEntriesFast(); i0++) {
    baconhep::TGenJet *pJet = (baconhep::TGenJet*) iGen.fGenJets->At(i0);
    bool pMatch = false;
    for(int i1 = 0; i1 < int(fLeptons.size()); i1++) { 
      if(Tools::deltaR(fLeptons[i1].Eta(),fLeptons[i1].Phi(),pJet->eta,pJet->phi) < 0.5) pMatch = true;
    }
    if(pMatch) continue;
    if(pJet->pt > 30) fNJets++;
    if(pJet->pt < 15) continue;
    insert(pJet,lJets);
  }
  if(lJets.size() >  0) fJetPt1  = lJets[0]->pt;
  if(lJets.size() >  0) fJetEta1 = lJets[0]->eta;
  if(lJets.size() >  0) fJetPhi1 = lJets[0]->phi;
  if(lJets.size() >  0) fJetM1   = lJets[0]->mass;
  if(lJets.size() >  1) fJetPt2  = lJets[1]->pt;
  if(lJets.size() >  1) fJetEta2 = lJets[1]->eta;
  if(lJets.size() >  1) fJetPhi2 = lJets[1]->phi;
  if(lJets.size() >  1) fJetM2   = lJets[1]->mass;
  if(lJets.size() >  2) fJCPt    = lJets[2]->pt;
  if(lJets.size() >  2) fJCEta   = lJets[2]->eta;;
  if(lJets.size() >  2) fJCPhi   = lJets[2]->phi;
  if(lJets.size() >  2) fJCM     = lJets[2]->mass;

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
  if(fJetPt1 == 0 || fJetPt2 == 0) return;
  TLorentzVector lVec1,lVec2;
  lVec1.SetPtEtaPhiM(fJetPt1,fJetEta1,fJetPhi1,fJetM1);
  lVec2.SetPtEtaPhiM(fJetPt2,fJetEta2,fJetPhi2,fJetM2);
  fMJJ   = (lVec1+lVec2).M();
  fJDEta =  fJetEta1-fJetEta2;
  fJDPhi = (fJetPhi1-fJetPhi2);
  if(fabs(fJDPhi) > 2.*TMath::Pi()-fabs(fJDPhi) && fJDPhi > 0)  fJDPhi -= 2.*TMath::Pi();
  if(fabs(fJDPhi) > 2.*TMath::Pi()-fabs(fJDPhi) && fJDPhi < 0)  fJDPhi += 2.*TMath::Pi();
}
void Jets::setSmear(int iBin) { 
  fPtSmear->SetParameter(0,fPtPar0[iBin]);
  fPtSmear->SetParameter(1,fPtPar1[iBin]);
  fPtSmear->SetParameter(2,fPtPar2[iBin]);
  fPtSmear->SetParameter(3,fPtPar3[iBin]);

  fPhiSmear->SetParameter(0,fPhiPar0[iBin]);
  fPhiSmear->SetParameter(1,fPhiPar1[iBin]);
  fPhiSmear->SetParameter(2,fPhiPar2[iBin]);
  fPhiSmear->SetParameter(3,fPhiPar3[iBin]);
  fPhiSmear->SetParameter(4,fPhiPar4[iBin]);

  fEtaSmear->SetParameter(0,fEtaPar0[iBin]);
  fEtaSmear->SetParameter(1,fEtaPar1[iBin]);
  fEtaSmear->SetParameter(2,fEtaPar2[iBin]);
  fEtaSmear->SetParameter(3,fEtaPar3[iBin]);
  fEtaSmear->SetParameter(4,fEtaPar4[iBin]);
}
void Jets::smear(float &iJPt,float &iJEta,float &iJPhi,float &iJM) { 
  //PAS JME-10-014. + http://arxiv.org/pdf/1107.4277v1.pdf
  //http://cmslxr.fnal.gov/lxr/source/CondFormats/JetMETObjects/data/Spring10_EtaResolution_AK5PF.txt
  //http://cmslxr.fnal.gov/lxr/source/CondFormats/JetMETObjects/data/Spring10_PhiResolution_AK5PF.txt
  //http://cmslxr.fnal.gov/lxr/source/CondFormats/JetMETObjects/data/Spring10_PtResolution_AK5PF.txt
  int iBin = 9-int(fabs(iJEta)*2.);
  if(iBin < 0) iBin = 0;
  setSmear(iBin);
  double lSigma    = fPtSmear ->Eval(iJPt)*iJPt;
  double lSigmaEta = fEtaSmear->Eval(iJPt);//*iJPt;
  double lSigmaPhi = fPhiSmear->Eval(iJPt);//*iJPt;
  double lOldJPt   = iJPt;
  iJPt  = fRandom->Gaus(iJPt, lSigma*2.);
  iJEta = fRandom->Gaus(iJEta,lSigmaEta);
  iJPhi = fRandom->Gaus(iJPhi,lSigmaPhi);
  iJM   = (iJPt/lOldJPt) * iJM;
}
void Jets::smear() { 
  fGJetPt1 = fJetPt1;
  fGJetEta1= fJetEta1;
  fGJetPhi1= fJetPhi1;
  fGJetM1  = fJetM1;

  fGJetPt2  = fJetPt2;
  fGJetEta2 = fJetEta2;
  fGJetPhi2 = fJetPhi2;
  fGJetM2   = fJetM2;

  fGJCPt    = fJCPt;
  fGJCEta   = fJCEta;
  fGJCPhi   = fJCPhi;
  fGJCM     = fJCM;
 
  smear(fJetPt1,fJetEta1,fJetPhi1,fJetM1);
  smear(fJetPt2,fJetEta2,fJetPhi2,fJetM2);
  smear(fJCPt  ,fJCEta  ,fJCPhi  ,fJCM  );
}
void Jets::fillJets(Gen &iGen) { 
  selectJets(iGen);
  smear(); //=> *add smearing functions !!! after Selection
  computeFinal();
}
void Jets::fillJets(GenJet &iGen) { 
  selectJets(iGen);
  smear(); 
  computeFinal();
}
