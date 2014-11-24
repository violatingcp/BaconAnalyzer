#include "../include/JetLoader.hh"
#include <cmath>
#include <iostream>
#include "TVector2.h"

using namespace baconhep;

JetLoader::JetLoader(TTree *iTree,bool iData,std::string iHLTFile) { 
  fJets     = new TClonesArray("baconhep::TJet");
  fAddJets  = new TClonesArray("baconhep::TAddJet");
  iTree->SetBranchAddress("Jet05",       &fJets);
  iTree->SetBranchAddress("AddJet05",    &fAddJets);
  fJetBr         = iTree->GetBranch("Jet05");
  fAddJetBr      = iTree->GetBranch("AddJet05");

  fFatJets     = new TClonesArray("baconhep::TJet");
  fFatAddJets  = new TClonesArray("baconhep::TAddJet");
  iTree->SetBranchAddress("Jet08CHS",       &fFatJets);
  iTree->SetBranchAddress("AddJet08CHS",    &fFatAddJets);
  fFatJetBr      = iTree->GetBranch("Jet08CHS");
  fFatAddJetBr   = iTree->GetBranch("AddJet08CHS");

  fTrigger = new TTrigger(iHLTFile);

  fFPtr1   = new TLorentzVector();
  fFPtr2   = new TLorentzVector();
  fPtr1    = new TLorentzVector();
  fPtr2    = new TLorentzVector();
  fPtr3    = new TLorentzVector();
  fPtr4    = new TLorentzVector();
  fDiJet   = new TLorentzVector();

  fReader = new TMVA::Reader( "!Color:!Silent" );    
  std::string lJetName    = "BDT";
  std::string lWeightFile = "/afs/cern.ch/work/p/pharris/public/Dan/Some/weights/TMVAClassificationCategory_BDT_simple_alpha.weights.xml";
  fQG1 = 0; fReader->AddVariable("jet1csv1",   &fQG1);
  fQG2 = 0; fReader->AddVariable("jet1csv2",   &fQG2);
  fJT1 = 0; fReader->AddVariable("jet1tau1",   &fJT1);
  fJT2 = 0; fReader->AddVariable("jet1tau2",   &fJT2);
  fJM1 = 0; fReader->AddVariable("jet1mass_m1",&fJM1);
  fJM2 = 0; fReader->AddVariable("jet1mass_t1",&fJM2);
  fJC2 = 0; fReader->AddVariable("jet1c2_2P0" ,&fJC2);
  fJQ  = 0; fReader->AddVariable("jet1qjet"   ,&fJQ );
  fReader->BookMVA(lJetName .c_str(),lWeightFile.c_str());

  std::vector<JetCorrectorParameters> corrParams;
  corrParams.push_back(JetCorrectorParameters("/afs/cern.ch/work/p/pharris/public/JEC/FT53_V21A_AN6_L1FastJet_AK7PFchs.txt"));
  corrParams.push_back(JetCorrectorParameters("/afs/cern.ch/work/p/pharris/public/JEC/FT53_V21A_AN6_L2Relative_AK7PFchs.txt"));
  corrParams.push_back(JetCorrectorParameters("/afs/cern.ch/work/p/pharris/public/JEC/FT53_V21A_AN6_L3Absolute_AK7PFchs.txt"));
  if(iData) corrParams.push_back(JetCorrectorParameters("/afs/cern.ch/work/p/pharris/public/JEC/FT53_V21A_AN6_L2L3Residual_AK7PFchs.txt"));
  //JetCorrectorParameters param("/afs/cern.ch/work/p/pharris/public/JEC/FT53_V21A_AN6_L3Absolute_AK7PFunc.txt");
  fJetCorr = new FactorizedJetCorrector(corrParams);
  //fJetUnc = new JetCorrectionUncertainty(param);
}
JetLoader::~JetLoader() { 
  delete fJets;
  delete fFatJets;
  delete fJetBr;
  delete fFatJetBr;
  delete fAddJets;
  delete fFatAddJets;
  delete fAddJetBr;
  delete fFatAddJetBr;

  delete fPtr1;
  delete fPtr2; 
  delete fPtr3;
  delete fPtr4;
  delete fFPtr1sj1;
  delete fFPtr2sj1; 
  delete fFPtr1sj2;
  delete fFPtr2sj2; 
  delete fDiJet;
}
void JetLoader::reset(TJet &iJet,TAddJet &iAddJet) { 
  //Jet Properties => Probably a better way (but who cares for now)
  iJet.pt   = 0; 
  iJet.eta  = 0; 
  iJet.phi  = 0; 
  iJet.csv  = 0; 
  iJet.csv1 = 0; 
  iJet.csv2 = 0; 
  iJet.mva  = 0; 
  iJet.qgid = 0; 
  iJet.qg1  = 0; 
  iJet.qg2  = 0; 
  iJet.tau1 = 0; 
  iJet.tau2 = 0;  
  iJet.tau3 = 0; 
  iJet.tau4 = 0; 
  iJet.prunedm    = 0; 
  iJet.nCharged   = 0; 
  iJet.nNeutrals  = 0; 
  iJet.nParticles = 0; 
  iJet.beta       = 0; 
  iJet.betaStar   = 0; 
  iJet.dR2Mean    = 0; 
  iJet.ptD        = 0; 
  iJet.q          = 0; 
  iJet.pullAngle  = 0; 
  iJet.pullY      = 0; 
  iJet.pullPhi    = 0; 
  iJet.chEmFrac   = 0; 
  iJet.neuEmFrac  = 0; 
  iJet.chHadFrac  = 0; 
  iJet.neuHadFrac = 0; 
  iJet.mcFlavor   = 0; 
  iJet.mcFlavorPhys = 0; 
  iJet.genpt      = 0; 
  iJet.geneta     = 0; 
  iJet.genphi     = 0; 
  iJet.genm       = 0; 
  iAddJet.pt_p1     = 0; 
  iAddJet.phi_p1    = 0; 
  iAddJet.eta_p1    = 0; 
  iAddJet.mass_p1   = 0; 
  iAddJet.pt_t1     = 0; 
  iAddJet.phi_t1    = 0; 
  iAddJet.eta_t1    = 0; 
  iAddJet.mass_t1   = 0; 
  iAddJet.qjet      = 0; 
  iAddJet.c2_0      = 0; 
  iAddJet.c2_1P0    = 0; 
  iAddJet.c2_2P0    = 0; 
}
void JetLoader::reset() { 
  fNJets      = 0;
  fNoiseClean = 0; 
  fNBTags     = 0; 
  fNBTags10   = 0; 
  fNQTags     = 0; 
  fJDEta      = 0; 
  fJDPhi      = 0; 
  fJDPullY    = 0; 
  fJDPullPhi  = 0; 

  fPtr1->SetPtEtaPhiM(1e-9,0,0,0);
  fPtr2->SetPtEtaPhiM(1e-9,0,0,0);
  fPtr3->SetPtEtaPhiM(1e-9,0,0,0);
  fPtr4->SetPtEtaPhiM(1e-9,0,0,0);
  fDiJet->SetPtEtaPhiM(1e-9,0,0,0);
  fFPtr1->SetPtEtaPhiM(1e-9,0,0,0);
  fFPtr2->SetPtEtaPhiM(1e-9,0,0,0);

  fFPtr1sj1->SetPtEtaPhiM(1e-9,0,0,0);
  fFPtr2sj1->SetPtEtaPhiM(1e-9,0,0,0);
  fFPtr1sj2->SetPtEtaPhiM(1e-9,0,0,0);
  fFPtr2sj2->SetPtEtaPhiM(1e-9,0,0,0);

  reset(fJet1,fAJet1); 
  reset(fJet2,fAJet2);
  reset(fJet3,fAJet3);
  reset(fJet4,fAJet4);

  reset(fFJet1,fFAJet1); 
  reset(fFJet2,fFAJet2);
}
void JetLoader::setupTree(TTree *iTree) { 
  reset();
  fTree = iTree;
  fTree->Branch("njets"        , &fNJets     , "fNJets/i");
  fTree->Branch("noiseCleaning", &fNoiseClean, "fNoiseClean/i");
  fTree->Branch("nbtags"       , &fNBTags    , "fNBTags/i");
  fTree->Branch("nbtags10"     , &fNBTags10  , "fNBTags10/i");
  fTree->Branch("nqtags"       , &fNQTags    , "fNQTags/i");
  fTree->Branch("jdeta"        , &fJDEta     , "fJDEta/F");
  fTree->Branch("jdphi"        , &fJDPhi     , "fJDPhi/F");
  fTree->Branch("jdpullY"      , &fJDPullY   , "fJDPullY/F");
  fTree->Branch("jdpullPhi"    , &fJDPullPhi , "fJDPullPhi/F");
    
  fTree->Branch("jet1"         , "TLorentzVector", &fPtr1);
  fTree->Branch("jet2"         , "TLorentzVector", &fPtr2);
  fTree->Branch("jet3"         , "TLorentzVector", &fPtr3);
  fTree->Branch("jet4"         , "TLorentzVector", &fPtr4);

  fTree->Branch("fjet1"         , "TLorentzVector", &fFPtr1);
  fTree->Branch("fjet2"         , "TLorentzVector", &fFPtr2);
  //Subjet
  fTree->Branch("fjet1sj1"      , "TLorentzVector", &fFPtr1sj1);
  fTree->Branch("fjet2sj1"      , "TLorentzVector", &fFPtr2sj1);
  fTree->Branch("fjet1sj2"      , "TLorentzVector", &fFPtr1sj2);
  fTree->Branch("fjet2sj2"      , "TLorentzVector", &fFPtr2sj2);
  fTree->Branch("dijet"         , "TLorentzVector", &fDiJet);
  //Jet Properites
  fTree->Branch("jet1CHF"      , &fJet1.chHadFrac    , "fJet1CHF/F");
  fTree->Branch("jet1NHF"      , &fJet1.neuHadFrac   , "fJet1NHF/F");
  fTree->Branch("jet1NEMF"     , &fJet1.neuEmFrac    , "fJet1NEMF/F");
  fTree->Branch("jet1Unc"      , &fJet1.unc          , "fJet1Unc/F");
  fTree->Branch("jet1Btag"     , &fJet1.csv          , "fJet1Btag/F");
  fTree->Branch("jet1QGtag"    , &fJet1.qgid         , "fJet1QGtag/F");
  fTree->Branch("jet1PullY"    , &fJet1.pullY        , "fJet1Pull/F");
  fTree->Branch("jet1PullPhi"  , &fJet1.pullPhi      , "fJet1PullPhi/F");
  fTree->Branch("jet1q"        , &fJet1.q            , "fJet1q/F");
  fTree->Branch("jet1PartonId" , &fJet1.mcFlavorPhys, "fJet1PartonId/I");

  fTree->Branch("fjet1CHF"      , &fFJet1.chHadFrac     , "fFJet1CHF/F");
  fTree->Branch("fjet1NHF"      , &fFJet1.neuHadFrac    , "fFJet1NHF/F");
  fTree->Branch("fjet1NEMF"     , &fFJet1.neuEmFrac     , "fFJet1NEMF/F");
  fTree->Branch("fjet1Unc"      , &fFJet1.unc           , "fFJet1Unc/F");
  fTree->Branch("fjet1Btag"     , &fFJet1.csv           , "fFJet1Btag/F");
  fTree->Branch("fjet1QGtag"    , &fFJet1.qgid          , "fFJet1QGtag/F");
  fTree->Branch("fjet1PullY"    , &fFJet1.pullY         , "fFJet1Pull/F");
  fTree->Branch("fjet1PullPhi"  , &fFJet1.pullPhi       , "fFJet1PullPhi/F");
  fTree->Branch("fjet1PullAngle", &fFJet1.pullAngle     , "fFJet1PullAngle/F");
  fTree->Branch("fjet1PartonId" , &fFJet1.mcFlavorPhys  , "fFJet1PartonId/I");
  fTree->Branch("fjet1tau1"     , &fFJet1.tau1         , "fFJet1Tau1/F");
  fTree->Branch("fjet1tau2"     , &fFJet1.tau2         , "fFJet1Tau2/F");
  fTree->Branch("fjet1tau3"     , &fFJet1.tau3         , "fFJet1Tau3/F");
  fTree->Branch("fjet1tau4"     , &fFJet1.tau4         , "fFJet1Tau4/F");
  fTree->Branch("fjet1mtrim"    , &fFAJet1.mass_t1     , "fFJet1MTrim/F");
  fTree->Branch("fjet1mprune"   , &fFJet1.prunedm      , "fFJet1MPrune/F");
  fTree->Branch("fjet1npart"    , &fFJet1.nParticles   , "fFJet1NPart/F");
  fTree->Branch("fjet1qg1"      , &fFJet1.csv1         , "fFJet1csv1/F");
  fTree->Branch("fjet1qg2"      , &fFJet1.csv2         , "fFJet1csv2/F");
  fTree->Branch("fjet1csv1"     , &fFJet1.qg1          , "fFJet1qg1/F");
  fTree->Branch("fjet1csv2"     , &fFJet1.qg2          , "fFJet1qg2/F");
  fTree->Branch("fjet1q"        , &fFJet1.q            , "fFJet1q/F");
  fTree->Branch("fjet1qjet"     , &fFAJet1.qjet        , "fFJet1qjet/F");
  fTree->Branch("fjet1c2b0"     , &fFAJet1.c2_0        , "fFJet1c2_0/F");
  fTree->Branch("fjet1c2b1"     , &fFAJet1.c2_1P0      , "fFJet1c2_1P0/F");
  fTree->Branch("fjet1c2b2"     , &fFAJet1.c2_2P0      , "fFJet1c2_2P0/F");
  fTree->Branch("fjet1Wmva"     , &fFJet1.mva          , "fFJet1mva/F");

  fTree->Branch("jet2CHF"     , &fJet2.chHadFrac    , "fJet2CHF/F");
  fTree->Branch("jet2NHF"     , &fJet2.neuHadFrac   , "fJet2NHF/F");
  fTree->Branch("jet2NEMF"    , &fJet2.neuEmFrac    , "fJet2NEMF/F");
  fTree->Branch("jet2Unc"     , &fJet2.unc          , "fJet2Unc/F");
  fTree->Branch("jet2Btag"    , &fJet2.csv          , "fJet2Btag/F");
  fTree->Branch("jet2QGtag"   , &fJet2.qgid         , "fJet2QGtag/F");
  fTree->Branch("jet2q"       , &fJet2.q            , "fFJet2q/F");
  fTree->Branch("jet2PartonId", &fJet2.mcFlavorPhys , "fJet2PartonId/I");

  fTree->Branch("fjet2CHF"      , &fFJet2.chHadFrac    , "fJet2CHF/F");
  fTree->Branch("fjet2NHF"      , &fFJet2.neuHadFrac   , "fJet2NHF/F");
  fTree->Branch("fjet2NEMF"     , &fFJet2.neuEmFrac    , "fJet2NEMF/F");
  fTree->Branch("fjet2Unc"      , &fFJet2.unc          , "fJet2Unc/F");
  fTree->Branch("fjet2Btag"     , &fFJet2.csv          , "fJet2Btag/F");
  fTree->Branch("fjet2QGtag"    , &fFJet2.qgid         , "fJet2QGtag/F");
  fTree->Branch("fjet2pullY"    , &fFJet2.pullY        , "fFJet2PullEta/F");
  fTree->Branch("fjet2pullPhi"  , &fFJet2.pullPhi      , "fFJet2PullPhi/F");
  fTree->Branch("fjet2PullAngle", &fFJet2.pullAngle    , "fFJet2PullAngle/F");
  fTree->Branch("fjet2PartonId" , &fFJet2.mcFlavorPhys , "fJet2PartonId/I");

  fTree->Branch("fjet2tau1"    , &fFJet2.tau1         , "fFJet2Tau1/F");
  fTree->Branch("fjet2tau2"    , &fFJet2.tau2         , "fFJet2Tau2/F");
  fTree->Branch("fjet2tau3"    , &fFJet2.tau3         , "fFJet2Tau3/F");
  fTree->Branch("fjet2tau4"    , &fFJet2.tau4         , "fFJet2Tau4/F");
  fTree->Branch("fjet2mtrim"   , &fFAJet2.mass_t1     , "fFJet2MTrim/F");
  fTree->Branch("fjet2mprune"  , &fFJet2.prunedm      , "fFJet2MPrune/F");
  fTree->Branch("fjet2npart"   , &fFJet2.nParticles   , "fFJet2NPart/F");
  fTree->Branch("fjet2qg1"     , &fFJet2.csv1         , "fFJet2csv1/F");
  fTree->Branch("fjet2qg2"     , &fFJet2.csv2         , "fFJet2csv2/F");
  fTree->Branch("fjet2csv1"    , &fFJet2.qg1          , "fFJet2qg1/F");
  fTree->Branch("fjet2csv2"    , &fFJet2.qg2          , "fFJet2qg2/F");
  fTree->Branch("fjet2q"       , &fFJet2.q            , "fFJet2q/F");
  fTree->Branch("fjet2qjet"    , &fFAJet2.qjet        , "fFJet2qjet/F");
  fTree->Branch("fjet2c2b0"    , &fFAJet2.c2_0        , "fFJet1c2_0/F");
  fTree->Branch("fjet2c2b1"    , &fFAJet2.c2_1P0      , "fFJet1c2_1P0/F");
  fTree->Branch("fjet2c2b2"    , &fFAJet2.c2_2P0      , "fFJet1c2_2P0/F");
  fTree->Branch("fjet2Wmva"    , &fFJet2.mva          , "fFJet2mva/F");
  
  fTree->Branch("jet3CHF"     , &fJet3.chHadFrac    , "fJet3CHF/F");
  fTree->Branch("jet3NHF"     , &fJet3.neuHadFrac   , "fJet3NHF/F");
  fTree->Branch("jet3NEMF"    , &fJet3.neuEmFrac    , "fJet3NEMF/F");
  fTree->Branch("jet3Unc"     , &fJet3.unc          , "fJet3Unc/F");
  fTree->Branch("jet3Btag"    , &fJet3.csv          , "fJet3Btag/F");
  fTree->Branch("jet3QGtag"   , &fJet3.qgid         , "fJet3QGtag/F");
  fTree->Branch("jet3PartonId", &fJet3.mcFlavorPhys , "fJet3PartonId/I");
  fTree->Branch("jet3tau1"    , &fJet3.tau1         , "fJet3Tau1/F");
  fTree->Branch("jet3tau2"    , &fJet3.tau2         , "fJet3Tau2/F");
  fTree->Branch("jet3tau3"    , &fJet3.tau3         , "fJet3Tau3/F");
  fTree->Branch("jet3tau4"    , &fJet3.tau4         , "fJet3Tau4/F");
  fTree->Branch("jet3mtrim"   , &fAJet3.mass_t1     , "fJet3MTrim/F");
  fTree->Branch("jet3mprune"  , &fJet3.prunedm      , "fJet3MPrune/F");
  fTree->Branch("jet3npart"   , &fJet3.nParticles   , "fJet3NPart/F");
  fTree->Branch("jet3qg1"     , &fJet3.csv1         , "fJet3csv1/F");
  fTree->Branch("jet3qg2"     , &fJet3.csv2         , "fJet3csv2/F");
  fTree->Branch("jet3csv1"    , &fJet3.qg1          , "fJet3qg1/F");
  fTree->Branch("jet3csv2"    , &fJet3.qg2          , "fJet3qg2/F");
  fTree->Branch("jet3q"       , &fJet3.q            , "fJet3q/F");
  //fTree->Branch("jet3pull"    , &fJet3.pull         , "fJet3pull/F");

  fTree->Branch("jet4Btag"    , &fJet4.csv          , "fJet4Btag/F");
  fTree->Branch("jet4QGtag"   , &fJet4.qgid         , "fJet4QGtag/F");

  //Redundant Info (for TTree draw => maybe I get rid of this)
  /*
  fTree->Branch("jpt_1"  ,&fJet1.pt  ,"fPt1/F");
  fTree->Branch("jeta_1" ,&fJet1.eta ,"fEta1/F");
  fTree->Branch("jphi_1" ,&fJet1.phi ,"fPhi1/F");
  fTree->Branch("jm_1"   ,&fJet1.mass,"fM1/F");

  fTree->Branch("jpt_2"  ,&fJet2.pt  ,"fPt2/F");
  fTree->Branch("jeta_2" ,&fJet2.eta ,"fEta2/F");
  fTree->Branch("jphi_2" ,&fJet2.phi ,"fPhi2/F");
  fTree->Branch("jm_2"   ,&fJet2.mass,"fM2/F");

  fTree->Branch("jpt_3"  ,&fJet3.pt  ,"fPt3/F");
  fTree->Branch("jeta_3" ,&fJet3.eta ,"fEta3/F");
  fTree->Branch("jphi_3" ,&fJet3.phi ,"fPhi3/F");
  fTree->Branch("jm_3"   ,&fJet3.mass,"fM3/F");

  fTree->Branch("jpt_4"  ,&fJet4.pt  ,"fPt4/F");
  fTree->Branch("jeta_4" ,&fJet4.eta ,"fEta4/F");
  fTree->Branch("jphi_4" ,&fJet4.phi ,"fPhi4/F");
  fTree->Branch("jm_4"   ,&fJet4.mass,"fM4/F");
  */
}
void JetLoader::load(int iEvent) { 
  fJets        ->Clear();
  fAddJets     ->Clear();
  fFatJets     ->Clear();
  fFatAddJets  ->Clear();
  fJetBr       ->GetEntry(iEvent);
  fAddJetBr    ->GetEntry(iEvent);
  fFatJetBr    ->GetEntry(iEvent);
  fFatAddJetBr ->GetEntry(iEvent);
}
double JetLoader::correction(TJet &iJet,double iRho) { 
  TLorentzVector lVec; lVec.SetPtEtaPhiM(iJet.ptRaw,iJet.eta,iJet.phi,iJet.mass);
  fJetCorr->setJetEta(iJet.eta);
  fJetCorr->setJetPt (iJet.ptRaw);
  fJetCorr->setJetPhi(iJet.phi);
  fJetCorr->setJetE  (lVec.E());
  fJetCorr->setRho   (iRho);
  fJetCorr->setJetA  (iJet.area);
  fJetCorr->setJetEMF(-99.0);     
  return ((fJetCorr->getCorrection())*iJet.ptRaw);
}
bool JetLoader::selectJets(std::vector<TLorentzVector> &iVetoes,double iRho) {
  reset(); 
  std::vector<TJet*> lJets;
  for  (int i0 = 0; i0 < fJets->GetEntriesFast(); i0++) { 
    TJet *pJet = (TJet*)((*fJets)[i0]);
    if(pJet->pt < 30 )  continue;
    if(!passPUId(pJet) && pJet->pt < 50 ) continue;
    //Veto
    bool pMatch = false;
    for(unsigned int i1 = 0; i1 < iVetoes.size(); i1++) {
      double pDEta = pJet->eta      - iVetoes[i1].Eta();
      double pDPhi = fabs(pJet->phi - iVetoes[i1].Phi());
      if(fabs(pDPhi) > 2.*TMath::Pi()-fabs(pDPhi)) pDPhi =  2.*TMath::Pi()-fabs(pDPhi);
      if(sqrt(pDPhi*pDPhi+pDEta*pDEta) > 0.5) continue;
      pMatch = true;
    }
    if(pMatch) continue;
    if(passLoose(pJet)    && pJet->pt > 30 && passPUId(pJet)) fNJets++;
    if(pJet->csv  > 0.25  && pJet->pt > 20 && passPUId(pJet)) fNBTags++;    
    if(pJet->csv  > 0.25  && pJet->pt > 10 && passPUId(pJet)) fNBTags10++;    
    if(pJet->qgid > 0.4   && pJet->pt > 20 && passPUId(pJet)) fNQTags++;
    //Start filling the collection
    if(lJets.size() == 0) {lJets.push_back(pJet); continue;}
    //Vector Manipulation to get sort the jets in corrected pT (collection is sorted in raw pt)
    bool pFill = false;
    for( std::vector<TJet*>::iterator pJetIter = lJets.begin(); pJetIter != lJets.end(); pJetIter++) { 
      if((*pJetIter)->pt > pJet->pt) continue;
      lJets.insert(pJetIter,pJet);
      pFill = true;
      break;
    } 
    if(!pFill) lJets.push_back(pJet);
    //Limit this to the top 4 Jets
    //if(lJets.size() > 4) lJets.pop_back();
  }
  //-1 means no corr
  if(lJets.size() > 0)   fillVars(lJets[0],-1,fPtr1,fJet1,fAJet1,0,0);
  if(lJets.size() > 1)   fillVars(lJets[1],-1,fPtr2,fJet2,fAJet2,0,0);
  if(lJets.size() > 2)   fillVars(lJets[2],-1,fPtr3,fJet3,fAJet3,0,0);
  if(lJets.size() > 3)   fillVars(lJets[3],-1,fPtr4,fJet4,fAJet4,0,0);

  TJet *lFJet1 = 0; 
  TJet *lFJet2 = 0; 
  if(lJets.size() > 0)   lFJet1  = fatJet(lJets[0]);
  if(lJets.size() > 1)   lFJet2  = fatJet(lJets[1]);
  if(lFJet1 != 0     )   fillVars(lFJet1,iRho,fFPtr1,fFJet1,fFAJet1,fFPtr1sj1,fFPtr1sj2);
  if(lFJet2 != 0     )   fillVars(lFJet2,iRho,fFPtr2,fFJet2,fFAJet2,fFPtr2sj1,fFPtr2sj2);

  if(lJets.size() > 1)   {
    TLorentzVector lDiJet = *fPtr1 + *fPtr2;
    fDiJet->SetPtEtaPhiM(lDiJet.Pt(),lDiJet.Eta(),lDiJet.Phi(),lDiJet.M());
    fJDPhi = fPtr1->DeltaPhi(*fPtr2);
    fJDEta = fabs(fPtr1->Eta() - fPtr2->Eta());
    float lJDY = (fPtr1->Rapidity() - fPtr2->Rapidity());
    //double lJMag = sqrt(lJets[0]->pullY*lJets[0]->pullY+lJets[0]->pullPhi*lJets[0]->pullPhi);
    //double l2Mag = sqrt(fPtr2->Rapidity()*fPtr2->Rapidity()+fPtr2->Phi()*fPtr2->Phi());
    TVector2 lJ1Vec; lJ1Vec.Set(lJets[0]->pullY,lJets[0]->pullPhi);
    TVector2 lJ2Vec; lJ2Vec.Set(lJDY,fJDPhi);
    lJ1Vec.Rotate(-lJ2Vec.Phi());
    fJDPullY   = lJ1Vec.Y();
    fJDPullPhi = lJ1Vec.X();
  }
  if(lJets.size() > 0) fNoiseClean |= int(fJet1.chHadFrac >0.2) << 0;
  if(lJets.size() > 0) fNoiseClean |= int(fJet1.neuHadFrac>0.7) << 1;
  if(lJets.size() > 0) fNoiseClean |= int(fJet1.neuEmFrac >0.7) << 2;
  if(lJets.size() > 1) fNoiseClean |= int(fJet2.chHadFrac >0.2) << 3;
  if(lJets.size() > 1) fNoiseClean |= int(fJet2.neuHadFrac>0.7) << 4;
  if(lJets.size() > 1) fNoiseClean |= int(fJet2.neuEmFrac >0.7) << 5;
  fHLTMatch   = passTrigObj(&fJet1,0) << 0;
  if(lJets.size() > 1) fHLTMatch   = ((passTrigObj(&fJet1,1) && passTrigObj(&fJet2,1))) << 1;// ||  (passTrigObj(&fJet1,2) && passTrigObj(&fJet2,1)))  << 1;
  return true;
}
void JetLoader::fillVars(TJet *iJet,double iRho,TLorentzVector *iPtr,TJet &iSaveJet,TAddJet &iASaveJet,TLorentzVector *iPtrsj1,TLorentzVector *iPtrsj2) { 
  double lPt =  iJet->pt;
  if(iRho  < 0) lPt = correction(*iJet,iRho);
  iSaveJet.pt   = lPt;
  iSaveJet.eta  = iJet->eta;
  iSaveJet.phi  = iJet->phi;
  iSaveJet.mass = iJet->mass;
  if(iJet->pt > 0) iPtr->SetPtEtaPhiM(lPt,iJet->eta,iJet->phi,iJet->mass);
  //if(iJet->pt > 0) iGPtr->SetPtEtaPhiM(iJet->pt,iJet->eta,iJet->phi,iJet->mass);
  iSaveJet.chHadFrac    = iJet->chHadFrac;
  iSaveJet.chEmFrac     = iJet->chEmFrac;
  iSaveJet.neuHadFrac   = iJet->neuHadFrac;
  iSaveJet.neuEmFrac    = iJet->neuEmFrac;
  iSaveJet.unc          = iJet->unc;
  iSaveJet.csv          = iJet->csv;
  iSaveJet.qgid         = iJet->qgid;
  iSaveJet.mcFlavorPhys = iJet->mcFlavorPhys;
  iSaveJet.pullY        = iJet->pullY;
  iSaveJet.pullAngle    = iJet->pullAngle;
  iSaveJet.q            = iJet->q;
  TAddJet *lAJet  = addJet(iJet);
  if(lAJet == 0) return;

  iSaveJet.tau1         = iJet->tau1;
  iSaveJet.tau2         = iJet->tau2;
  iSaveJet.tau3         = iJet->tau3;
  iSaveJet.tau4         = iJet->tau4;
  iSaveJet.mass         = iJet->mass;
  iSaveJet.prunedm      = iJet->prunedm;
  iSaveJet.nParticles   = iJet->nParticles;
  iSaveJet.csv1         = iJet->csv1;
  iSaveJet.csv2         = iJet->csv2;
  iSaveJet.qg1          = iJet->qg1;
  iSaveJet.qg2          = iJet->qg2;
  iSaveJet.pullAngle    = iJet->pullAngle;
  iSaveJet.hltMatchBits = iJet->hltMatchBits;
  iASaveJet.pt_p1       = lAJet->pt_p1    ;
  iASaveJet.phi_p1      = lAJet->phi_p1   ;
  iASaveJet.eta_p1      = lAJet->eta_p1   ;
  iASaveJet.mass_p1     = lAJet->mass_p1  ;
  iASaveJet.pt_t1       = lAJet->pt_t1    ;
  iASaveJet.phi_t1      = lAJet->phi_t1   ;
  iASaveJet.eta_t1      = lAJet->eta_t1   ;
  iASaveJet.mass_t1     = lAJet->mass_t1  ;
  iASaveJet.qjet        = lAJet->qjet     ;
  iASaveJet.c2_0        = lAJet->c2_0     ;
  iASaveJet.c2_1P0      = lAJet->c2_1P0   ;
  iASaveJet.c2_2P0      = lAJet->c2_2P0   ;
  iSaveJet .mva         = mva(iJet,lAJet);
  if(lAJet->sj1_pt > 0 && iPtrsj1 != 0) iPtrsj1->SetPtEtaPhiM(lAJet->sj1_pt,lAJet->sj1_eta,lAJet->sj1_phi,lAJet->sj1_m);
  if(lAJet->sj2_pt > 0 && iPtrsj2 != 0) iPtrsj2->SetPtEtaPhiM(lAJet->sj2_pt,lAJet->sj2_eta,lAJet->sj2_phi,lAJet->sj2_m);
}
bool JetLoader::vetoJet() {
  for  (int i0 = 0; i0 < fJets->GetEntriesFast(); i0++) { 
    TJet *pJet = (TJet*)((*fJets)[i0]);
    if(passVeto(pJet)) return true;
  }
  return false;
}
//PFJet Ids
bool JetLoader::passLoose(TJet *iJet) { 
  if(iJet->neuEmFrac        >  0.99)                         return false;
  if(iJet->neuHadFrac       >  0.99)                         return false;
  if(iJet->nParticles       <  2)                            return false;
  if(iJet->chHadFrac        <= 0     && fabs(iJet->eta) < 2.4 )    return false;
  if(iJet->chEmFrac         >  0.99  && fabs(iJet->eta) < 2.4 )    return false;
  if(iJet->nCharged         < 1      && fabs(iJet->eta) < 2.4 )    return false;
  return true;
}
bool JetLoader::passTight(TJet *iJet) { 
  if(iJet->neuEmFrac        >  0.9)                          return false;
  if(iJet->neuHadFrac       >  0.9)                          return false;
  if(iJet->nParticles       <  2)                            return false;
  if(iJet->chHadFrac        <= 0     && fabs(iJet->eta) < 2.4 )    return false;
  if(iJet->chEmFrac         >  0.9   && fabs(iJet->eta) < 2.4 )    return false;
  if(iJet->nCharged         < 1      && fabs(iJet->eta) < 2.4 )    return false;
  return true;
}
bool JetLoader::passVeto(TJet *iJet) { 
  return passLoose(iJet);
}
bool JetLoader::passPUId(TJet *iJet) { 
  if(iJet->pt >  0 && iJet->pt < 20    && fabs(iJet->eta) < 2.5)                            return (iJet->mva > -0.95); 
  if(iJet->pt >  0 && iJet->pt < 20    && fabs(iJet->eta) > 2.5  && fabs(iJet->eta) < 2.75) return (iJet->mva > -0.96);
  if(iJet->pt >  0 && iJet->pt < 20    && fabs(iJet->eta) > 2.75 && fabs(iJet->eta) < 3.0 ) return (iJet->mva > -0.94);
  if(iJet->pt >  0 && iJet->pt < 20    && fabs(iJet->eta) > 3.0  )                          return (iJet->mva > -0.95);
  if(iJet->pt > 20 && iJet->pt < 20000 && fabs(iJet->eta) < 2.5)                            return (iJet->mva > -0.63); 
  if(iJet->pt > 20 && iJet->pt < 20000 && fabs(iJet->eta) > 2.5  && fabs(iJet->eta) < 2.75) return (iJet->mva > -0.60);
  if(iJet->pt > 20 && iJet->pt < 20000 && fabs(iJet->eta) > 2.75 && fabs(iJet->eta) < 3.0 ) return (iJet->mva > -0.55);
  if(iJet->pt > 20 && iJet->pt < 20000 && fabs(iJet->eta) > 3.0  )                          return (iJet->mva > -0.45);
  return false;
}
TAddJet *JetLoader::addJet(TJet *iJet) { 
  int lIndex = -1;
  TAddJet *lJet = 0; 
  for(int i0 = 0; i0 < fFatJets->GetEntriesFast(); i0++) { 
    if((*fFatJets)[i0] == iJet) { lIndex = i0; break;}
  }
  if(lIndex == -1) return 0;
  for  (int i0 = 0; i0 < fFatAddJets->GetEntriesFast(); i0++) { 
    TAddJet *pJet = (TAddJet*)((*fFatAddJets)[i0]);
    if(pJet->index == fabs(lIndex)) { lJet = pJet; break;}
  }
  return lJet;
}
//Two triggers to be added
//MonoCentralPFJet80_PFMETnoMu
//HLT_DiPFJet40_PFMETnoMu65_MJJ800VBF_AllJets_v"
void JetLoader::addTrigger(std::string iName) { 
  fTrigString.push_back(iName);
}
bool JetLoader::passTrigObj(TJet *iJet,int iId) {
  bool lPass = false;
  if(fTrigger->passObj(fTrigString[iId],iJet->hltMatchBits))  lPass = true;
  return lPass;
}
TJet *JetLoader::fatJet(TJet *iJet) { 
  TJet *lJet = 0; 
  for(int i0 = 0; i0 < fFatJets->GetEntriesFast(); i0++) { 
    TJet *pJet   = (TJet*) (*fFatJets)[i0];
    double pDEta = pJet->eta-iJet->eta; 
    double pDPhi = pJet->phi-iJet->phi; 
    if(fabs(pDPhi) > TMath::Pi()*2.-fabs(pDPhi)) pDPhi = 2.*TMath::Pi()-fabs(pDPhi);
    if(sqrt(pDEta*pDEta+pDPhi*pDPhi) < 0.4) lJet = pJet;
    std::cout << "==> " << sqrt(pDEta*pDEta+pDPhi*pDPhi) << std::endl;
  }
  return lJet;
}
float JetLoader::mva(TJet *iJet,TAddJet *iAJet) { 
  fQG1 = iJet->qg1; 
  fQG2 = iJet->qg2; 
  fJT1 = iJet->tau1;
  fJT2 = iJet->tau2;
  fJM1 = iAJet->mass_m1;
  fJM2 = iAJet->mass_t1;
  fJC2 = iAJet->c2_2P0;
  fJQ  = iAJet->qjet;
  float lMVA      = float(fReader->EvaluateMVA("BDT"));
  return lMVA;
}
