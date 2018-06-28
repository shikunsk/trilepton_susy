//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////

#ifndef FinalYield_h
#define FinalYield_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "vector"
//Include MT2 Calculator
#include "mt2_bisect.h"

using namespace std;

class FinalYield {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   ULong64_t       _eventNb;
   ULong64_t       _runNb;
   ULong64_t       _lumiBlock;
   Bool_t          passmetfilters;
   Bool_t          passHLT_DoubleMu8_Mass8_PFHT300;
   Bool_t          passHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL;
   Bool_t          passHLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL;
   Bool_t          passHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ;
   Bool_t          passHLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ;
   Bool_t          passHLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT300;
   Bool_t          passHLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL;
   Bool_t          passHLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ;
   Bool_t          passHLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT300;
   Bool_t          passHLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL;
   Bool_t          passHLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL;
   Bool_t          passHLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL;
   Bool_t          passHLT_IsoMu24_eta2p1;
   Bool_t          passHLT_Ele27_eta2p1_WPLoose_Gsf;
   Bool_t          passHLT_Ele32_eta2p1_WPLoose_Gsf;
   Bool_t          passHLT_Mu8;
   Bool_t          passHLT_Mu17;
   Bool_t          passHLT_Mu8_TrkIsoVVL;
   Bool_t          passHLT_Mu17_TrkIsoVVL;
   Bool_t          passHLT_Ele12_CaloIdM_TrackIdM_PFJet30;
   Bool_t          passHLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30;
   Bool_t          passHLT_Ele17_CaloIdM_TrackIdM_PFJet30;
   Bool_t          passHLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30;
   Bool_t          passHLT_Ele27_WPTight_Gsf;
   Bool_t          passHLT_IsoMu24;
   Bool_t          passHLT_IsoTkMu24;
   Bool_t          passHLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ;
   Bool_t          passHLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ;
   Bool_t          passHLT_Ele22_eta2p1_WPLoose_GSF_LooseIsoPFtau20_SingleL1;
   Bool_t          passHLT_Ele24_eta2p1_WPLoose_GSF_LooseIsoPFtau30;
   Bool_t          passHLT_Ele36_eta2p1_WPLoose_GSF_LooseIsoPFtau20_SingleL1;
   Bool_t          passHLT_IsoMu19_eta2p1_LooseIsoPFTau20_SingleL1_v;
   Double_t        Prescale_passHLT_Mu8;
   Double_t        Prescale_passHLT_Mu17;
   Double_t        Prescale_passHLT_Mu8_TrkIsoVVL;
   Double_t        Prescale_passHLT_Mu17_TrkIsoVVL;
   Double_t        Prescale_passHLT_Ele12_CaloIdM_TrackIdM_PFJet30;
   Double_t        Prescale_passHLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30;
   Double_t        Prescale_passHLT_Ele17_CaloIdM_TrackIdM_PFJet30;
   Double_t        Prescale_passHLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30;
   Double_t        _weight;
   Double_t        _genHT;
   Int_t           _nLeptons;
   Int_t           _index1;
   Int_t           _index2;
   Bool_t          _sb;
   Bool_t          _doubleF;
   Double_t        _PVchi2;
   Double_t        _PVerr[3];
   Int_t           _n_PV;
   Double_t        _met;
   Double_t        _met_phi;
   Double_t        _met_JECup;
   Double_t        _met_JECdown;
   Double_t        _met_phi_JECup;
   Double_t        _met_phi_JECdown;
   Double_t        HT;
   Double_t        _genmet;
   Double_t        _genmet_phi;
   Int_t           _n_bJets;
   Int_t           _n_Jets;
   Int_t           _n_Jets40;
   Int_t           trueNVtx;
   vector<int>     *_closeIndex;
   vector<int>     *_indeces;
   vector<int>     *_flavors;
   vector<int>     *_charges;
   vector<int>     *_chargesMC;
   vector<double>  *_isolation;
   vector<double>  *_miniisolation_0p2;
   vector<double>  *_miniisolation_0p3;
   vector<double>  *_miniisoneutral;
   vector<double>  *_miniisolationcharged_0p2;
   vector<int>     *_trackSelectionMultiplicity;
   vector<double>  *_muonSegmentComp;
   vector<bool>    *_multiisolation_T;
   vector<bool>    *_multiisolation_M;
   vector<bool>    *_multiisolation_L;
   vector<double>  *_pfisocharged;
   vector<double>  *_pfisoneutral;
   vector<double>  *_pfisophoton;
   vector<double>  *_isolationMCraw;
   vector<double>  *_isolationMCnonu;
   vector<double>  *_isolationMCdr03;
   vector<double>  *_isolationMCdr03nonu;
   vector<double>  *_ptrel;
   vector<double>  *_ptrel2;
   vector<double>  *_ptratio;
   vector<double>  *_mvaValue;
   vector<double>  *_mt;
   vector<double>  *_mllZ;
   vector<double>  *_mllG;
   vector<double>  *_mll;
   vector<vector<double> > *_mllv;
   vector<int>     *_origin;
   vector<int>     *_originReduced;
   vector<double>  *_ipPV;
   vector<double>  *_ipPVerr;
   vector<double>  *_ipPVmc;
   vector<double>  *_ipZPV;
   vector<double>  *_ipZPVerr;
   vector<double>  *_3dIP;
   vector<double>  *_3dIPerr;
   vector<double>  *_3dIPsig;
   vector<double>  *_closeJetPtAll;
   vector<double>  *_closeJetEtaAll;
   vector<double>  *_closeJetPhiAll;
   vector<double>  *_closeJetMAll;
   vector<double>  *_closeJetCSVAll;
   vector<int>     *_closeJetNconstAll;
   vector<double>  *_closeJetAngAll;
   vector<double>  *_ptRelAll;
   vector<double>  *_closeJetPtAllMC;
   vector<double>  *_closeJetPtAllstatus;
   vector<int>     *_partonIdMatched;
   vector<bool>    *_sameParton;
   vector<bool>    *_isloose;
   vector<bool>    *_istight;
   vector<bool>    *_isMediumMuon;
   vector<bool>    *_istightIso;
   vector<bool>    *_istightID;
   vector<bool>    *_istightIDWP2016_RA7;
   vector<bool>    *_isFOIDWP2016_RA7;
   vector<bool>    *_istightIDWP2016;
   vector<bool>    *_isFOIDWP2016;
   vector<bool>    *_isVetoIDWP2016;
   vector<bool>    *_istightIDWP2016_EWK;
   vector<bool>    *_isFOIDWP2016_EWK;
   vector<bool>    *_isVetoIDWP2016_EWK;
   vector<bool>    *_istrigemulID;
   vector<bool>    *_istrigemulISO;
   vector<double>  *_mompt;
   vector<double>  *_momphi;
   vector<double>  *_mometa;
   vector<int>     *_mompdg;
   vector<double>  *_lmva;
   vector<double>  *_lPt;
   vector<double>  *_lEta;
   vector<double>  *_lPhi;
   vector<double>  *_lE;
   vector<double>  *_lPtmc;
   vector<double>  *_lEtamc;
   vector<double>  *_lPhimc;
   vector<double>  *_lEmc;
   vector<double>  *_nuPtmc;
   vector<double>  *_nuEtamc;
   vector<double>  *_nuPhimc;
   vector<double>  *_nuEmc;
   vector<double>  *_mtmc;
   vector<double>  *_genWPhi;
   vector<double>  *_genZPhi;
   vector<double>  *_genWEta;
   vector<double>  *_genZEta;
   vector<double>  *_genWPt;
   vector<double>  *_genZPt;
   vector<double>  *_genWE;
   vector<double>  *_genZE;
   vector<int>     *_lnmisshits;
   vector<int>     *_lchargeGSF;
   vector<int>     *_lchargeCTF;
   vector<int>     *_lchargePixSc;
   vector<double>  *_bTagged;
   vector<double>  *_jetEta;
   vector<double>  *_jetPhi;
   vector<double>  *_jetPt;
   vector<double>  *_jetM;
   vector<double>  *_csv;
   vector<double>  *_jetDeltaRloose;

    vector<double>  *_tau_dz;
    vector <bool> *_TauIs_againstMuonLoose3;
    vector <bool> *_TauIs_byLooseCombinedIsolationDeltaBetaCorr3Hits;
    vector <bool> *_TauIs_decayModeFindingNewDMs;
    vector <bool> *_TauIs_decayModeFinding;
    vector <bool> *_TauIs_byLooseIsolationMVArun2v1DBdR03oldDMwLT;
    vector <bool> *_TauIs_byMediumIsolationMVArun2v1DBdR03oldDMwLT;
    vector <bool> *_TauIs_byVLooseIsolationMVArun2v1DBoldDMwLT;
    vector <bool> *_TauIs_byTightIsolationMVArun2v1DBoldDMwLT;
    
    vector<double> *_jetJECuncty;
    vector<double> *_jetbtagSF;
    vector<double> *_jetbtagSF_up;
    vector<double> *_jetbtagSF_down;
    vector<double> *_jetbtagEff;

   // List of branches
   TBranch        *b__eventNb;   //!
   TBranch        *b__runNb;   //!
   TBranch        *b__lumiBlock;   //!
   TBranch        *b_passmetfilters;
   TBranch        *b_passHLT_DoubleMu8_Mass8_PFHT300;   //!
   TBranch        *b_passHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL;   //!
   TBranch        *b_passHLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL;   //!
   TBranch        *b_passHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ;   //!
   TBranch        *b_passHLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ;   //!
   TBranch        *b_passHLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT300;   //!
   TBranch        *b_passHLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL;   //!
   TBranch        *b_passHLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ;   //!
   TBranch        *b_passHLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT300;   //!
   TBranch        *b_passHLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL;   //!
   TBranch        *b_passHLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL;   //!
   TBranch        *b_passHLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL;   //!
   TBranch        *b_passHLT_IsoMu24_eta2p1;   //!
   TBranch        *b_passHLT_Ele27_eta2p1_WPLoose_Gsf;   //!
   TBranch        *b_passHLT_Ele32_eta2p1_WPLoose_Gsf;   //!
   TBranch        *b_passHLT_Mu8;
   TBranch        *b_passHLT_Mu17;
   TBranch        *b_passHLT_Mu8_TrkIsoVVL;
   TBranch        *b_passHLT_Mu17_TrkIsoVVL;
   TBranch        *b_passHLT_Ele12_CaloIdM_TrackIdM_PFJet30;
   TBranch        *b_passHLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30;
   TBranch        *b_passHLT_Ele17_CaloIdM_TrackIdM_PFJet30;
   TBranch        *b_passHLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30;
   TBranch        *b_passHLT_Ele27_WPTight_Gsf;
   TBranch        *b_passHLT_IsoMu24;
   TBranch        *b_passHLT_IsoTkMu24;
   TBranch        *b_passHLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ;
   TBranch        *b_passHLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ;
   TBranch        *b_passHLT_Ele22_eta2p1_WPLoose_GSF_LooseIsoPFtau20_SingleL1;
   TBranch        *b_passHLT_Ele24_eta2p1_WPLoose_GSF_LooseIsoPFtau30;
   TBranch        *b_passHLT_Ele36_eta2p1_WPLoose_GSF_LooseIsoPFtau20_SingleL1;
   TBranch        *b_passHLT_IsoMu19_eta2p1_LooseIsoPFTau20_SingleL1_v;
   TBranch        *b_Prescale_passHLT_Mu8;
   TBranch        *b_Prescale_passHLT_Mu17;
   TBranch        *b_Prescale_passHLT_Mu8_TrkIsoVVL;
   TBranch        *b_Prescale_passHLT_Mu17_TrkIsoVVL;
   TBranch        *b_Prescale_passHLT_Ele12_CaloIdM_TrackIdM_PFJet30;
   TBranch        *b_Prescale_passHLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30;
   TBranch        *b_Prescale_passHLT_Ele17_CaloIdM_TrackIdM_PFJet30;
   TBranch        *b_Prescale_passHLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30;
   TBranch        *b__weight;   //!
   TBranch        *b__genHT;   //!
   TBranch        *b__nLeptons;   //!
   TBranch        *b__index1;   //!
   TBranch        *b__index2;   //!
   TBranch        *b__sb;   //!
   TBranch        *b__doubleF;   //!
   TBranch        *b__PVchi2;   //!
   TBranch        *b__PVerr;   //!
   TBranch        *b__n_PV;   //!
   TBranch        *b__met;   //!
   TBranch        *b__met_phi;   //!
   TBranch        *b__met_JECup;   //!
   TBranch        *b__met_JECdown;   //!
   TBranch        *b__met_phi_JECup;   //!
   TBranch        *b__met_phi_JECdown;   //!
   TBranch        *b_HT;   //!
   TBranch        *b__genmet;   //!
   TBranch        *b__genmet_phi;   //!
   TBranch        *b__n_bJets;   //!
   TBranch        *b__n_Jets;   //!
   TBranch        *b__n_Jets40;//!
   TBranch        *b_trueNVtx;
   TBranch        *b__closeIndex;   //!
   TBranch        *b__indeces;   //!
   TBranch        *b__flavors;   //!
   TBranch        *b__charges;   //!
   TBranch        *b__chargesMC;   //!
   TBranch        *b__isolation;   //!
   TBranch        *b__miniisolation_0p2;   //!
   TBranch        *b__miniisolation_0p3;   //!
   TBranch        *b__miniisolationcharged_0p2;   //!
   TBranch        *b__miniisoneutral;   //!
   TBranch        *b__trackSelectionMultiplicity; //!
   TBranch        *b__muonSegmentComp;
   TBranch        *b__multiisolation_T;   //!
   TBranch        *b__multiisolation_M;   //!
   TBranch        *b__multiisolation_L;   //!
   TBranch        *b__pfisocharged;   //!
   TBranch        *b__pfisoneutral;   //!
   TBranch        *b__pfisophoton;   //!
   TBranch        *b__isolationMCraw;   //!
   TBranch        *b__isolationMCnonu;   //!
   TBranch        *b__isolationMCdr03;   //!
   TBranch        *b__isolationMCdr03nonu;   //!
   TBranch        *b__ptrel;   //!
   TBranch        *b__ptrel2;   //!
   TBranch        *b__ptratio;   //!
   TBranch        *b__mvaValue;   //!
   TBranch        *b__mt;   //!
   TBranch        *b__mllZ;   //!
   TBranch        *b__mllG;   //!
   TBranch        *b__mll;   //!
   TBranch        *b__mllv;   //!
   TBranch        *b__origin;   //!
   TBranch        *b__originReduced;   //!
   TBranch        *b__ipPV;   //!
   TBranch        *b__ipPVerr;   //!
   TBranch        *b__ipPVmc;   //!
   TBranch        *b__ipZPV;   //!
   TBranch        *b__ipZPVerr;   //!
   TBranch        *b__3dIP;   //!
   TBranch        *b__3dIPerr;   //!
   TBranch        *b__3dIPsig;   //!
   TBranch        *b__closeJetPtAll;   //!
   TBranch        *b__closeJetEtaAll;   //!
   TBranch        *b__closeJetPhiAll;   //!
   TBranch        *b__closeJetMAll;   //!
   TBranch        *b__closeJetCSVAll;   //!
   TBranch        *b__closeJetNconstAll;   //!
   TBranch        *b__closeJetAngAll;   //!
   TBranch        *b__ptRelAll;   //!
   TBranch        *b__closeJetPtAllMC;   //!
   TBranch        *b__closeJetPtAllstatus;   //!
   TBranch        *b__partonIdMatched;   //!
   TBranch        *b__sameParton;   //!
   TBranch        *b__isloose;   //!
   TBranch        *b__istight;   //!
   TBranch        *b__isMediumMuon;
   TBranch        *b__istightIso;   //!
   TBranch        *b__istightID;   //!
   TBranch        *b__istightIDWP2016_RA7; //!
   TBranch        *b__isFOIDWP2016_RA7;//!
   TBranch        *b__istightIDWP2016; //!
   TBranch        *b__isFOIDWP2016;//!
   TBranch        *b__isVetoIDWP2016;//!
   TBranch        *b__isVetoIDWP2016_EWK;
   TBranch        *b__isFOIDWP2016_EWK;
   TBranch        *b__istightIDWP2016_EWK;
   TBranch        *b__istrigemulID;//!
   TBranch        *b__istrigemulISO;//!
   TBranch        *b__mompt;   //!
   TBranch        *b__momphi;   //!
   TBranch        *b__mometa;   //!
   TBranch        *b__mompdg;   //!
   TBranch        *b__lmva;   //!
   TBranch        *b__lPt;   //!
   TBranch        *b__lEta;   //!
   TBranch        *b__lPhi;   //!
   TBranch        *b__lE;   //!
   TBranch        *b__lPtmc;   //!
   TBranch        *b__lEtamc;   //!
   TBranch        *b__lPhimc;   //!
   TBranch        *b__lEmc;   //!
   TBranch        *b__nuPtmc;   //!
   TBranch        *b__nuEtamc;   //!
   TBranch        *b__nuPhimc;   //!
   TBranch        *b__nuEmc;   //!
   TBranch        *b__mtmc;   //!
   TBranch        *b__genWPhi;
   TBranch        *b__genZPhi;
   TBranch        *b__genWEta;
   TBranch        *b__genZEta;
   TBranch        *b__genWPt;
   TBranch        *b__genZPt;
   TBranch        *b__genWE;
   TBranch        *b__genZE;
   TBranch        *b__lnmisshits;   //!
   TBranch        *b__lchargeGSF;   //!
   TBranch        *b__lchargeCTF;   //!
   TBranch        *b__lchargePixSc;   //!
   TBranch        *b__bTagged;   //!
   TBranch        *b__jetEta;   //!
   TBranch        *b__jetPhi;   //!
   TBranch        *b__jetPt;   //!
   TBranch        *b__jetM;   //!
   TBranch        *b__csv;   //!
   TBranch        *b__jetDeltaRloose;   //!
    
    TBranch        *b__tau_dz; //!
    //TBranch       *b__TauIs_againstElectronLooseMVA5;
    TBranch       *b__TauIs_againstMuonLoose3;
    TBranch       *b__TauIs_byLooseCombinedIsolationDeltaBetaCorr3Hits;
    TBranch       *b__TauIs_decayModeFindingNewDMs;
    TBranch       *b__TauIs_decayModeFinding;
    TBranch       *b__TauIs_byLooseIsolationMVArun2v1DBdR03oldDMwLT;
    TBranch       *b__TauIs_byMediumIsolationMVArun2v1DBdR03oldDMwLT;
    TBranch       *b__TauIs_byVLooseIsolationMVArun2v1DBoldDMwLT;
    TBranch       *b__TauIs_byTightIsolationMVArun2v1DBoldDMwLT;
    
    TBranch       *b__jetJECuncty;
    TBranch       *b__jetbtagSF;
    TBranch       *b__jetbtagSF_up;
    TBranch       *b__jetbtagSF_down;
    TBranch       *b__jetbtagEff;

   FinalYield(TTree *tree=0);
   virtual ~FinalYield();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
    
    
    /*
     *Define Functions
     */
    ////////////////////////////////////////////
    ////////////////////////////////////////////
    bool TriggerSelectionData(TString& Dataset);
    bool TriggerSelectionMC(TString& MCsample);
    bool bJetandMETcut(double& MET, const double& VarCoef);
    bool ObjectSelection(int& n);
    bool FOSelection(int& n);
    bool DeltaRSelection(int& n, int& l);
    bool RemoveJet(int& n, int& l1,int& l2);
    
    double CalcInvariantMass(int& n, int& l);
    double CalcTrilepInvariantMass(int& n,int& l, int& m);
    double CalcTrileppz(int& n,int& l, int& m);
    double CalcTransMass(int& n, double& MET, double& METPhi);
    double CalcMt2(double* pa, double* pb, double* pmiss);
    double CalcdR(int& n, int& l);
    double CalcAngZW(int& l1, int& l2, int& l3);
    double CalcAngZMET(int& l1, int& l2);
    double CalcAngZl1l2(int& l1, int& l2);
    double CalcZpt(int& l1, int& l2);
    double CalcZpz(int& l1, int& l2);
    double Calcu1(double& uphi, double& bosonphi, double& umag);
    double Calcu2(double& uphi, double& bosonphi, double& umag);
    double CalcFR(int& n, const double& VarCoef);
    double CalcBTagSF(const double& VarCoef);
    double CalcLepSF_ele_FullSim(TH2D* SF, TH2D* RecoSF, int& n, const double& VarCoef);
    double CalcLepSF_ele_FastSim(TH2D* SF, int& n, const double& VarCoef);
    double CalcLepSF_mu_FullSim(TH2D* SF1, TH2D* SF2, TH2D* SF3, TH2D* SF4, TGraphAsymmErrors* RecoSF, int& n, const double& VarCoef);
    double CalcLepSF_mu_FastSim(TH2D* SF1, TH2D* SF2, TH2D* SF3, TH2D* SF4, int& n, const double& VarCoef);
    double CalcLepSF_tau_FullSim(const double& VarCoef);
    double CalcLepSF_tau_FastSim(TH2D* SF, int& n, const double& VarCoef);
    
    int CalcSRIndex_A(double& MET, double& Mll, double& MT);
    int CalcSRIndex_B(double& MET, double& Mll, double& MT);
    int CalcSRIndex_C(double& MET, double& Mll, double& MT);
    int CalcSRIndex_D(double& MET, double& Mll, double& MT);
    int CalcSRIndex_E(double& MET, double& Mll, double& MT);
    int CalcSRIndex_F(double& MET, double& Mll, double& MT);
    
    void ShapeForLimitSetting_WZ(TH1F* Hist, int& nbin, const double& VarCoef);
    void AssignUncForPlot_SRBG_WZ(TH1F* Hist, int& nbin, TH1F* HistUp1, TH1F* HistDown1, TH1F* HistUp2, TH1F* HistDown2, TH1F* HistUp3, TH1F* HistDown3, TH1F* HistUp4, TH1F* HistDown4, TH1F* HistStatUnc, TH1F* HistTotalUnc);
    void AssignUncForPlot_SRBG_NonWZ(TH1F* Hist, double uncertainty, int& nbin, TH1F* HistUp1, TH1F* HistDown1, TH1F* HistUp2, TH1F* HistDown2, TH1F* HistUp3, TH1F* HistDown3, TH1F* HistUp4, TH1F* HistDown4, TH1F* HistStatUnc, TH1F* HistTotalUnc);
    void AssignUncForPlot_SRBG_Nonprompt(TH1F* Hist, double uncertainty, int& nbin, TH1F* HistUp, TH1F* HistDown, TH1F* HistStatUnc, TH1F* HistTotalUnc);
    void AssignUncForPlot_VarBG_WZ(TH1F* Hist, double uncertainty, int& nbin, TH1F* HistStatUnc, TH1F* HistTotalUnc);
    void AssignUncForPlot_VarBG_NonWZ(TH1F* Hist, double uncertainty, int& nbin, TH1F* HistStatUnc, TH1F* HistTotalUnc);
    void AssignUncForPlot_VarBG_Nonprompt(TH1F* Hist, double uncertainty, int& nbin, TH1F* HistStatUnc, TH1F* HistTotalUnc);
    ////////////////////////////////////////////
    ////////////////////////////////////////////
    /*
     *END OF Define Functions
     */
    /*
     *MT2 Calculator
     */
    mt2 mt2_calculator;
    /*
     *MT2 Calculator
     */
};

#endif

#ifdef FinalYield_cxx
FinalYield::FinalYield(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
  Init(tree);
}

FinalYield::~FinalYield()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t FinalYield::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t FinalYield::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void FinalYield::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   _closeIndex = 0;
   _indeces = 0;
   _flavors = 0;
   _charges = 0;
   _chargesMC = 0;
   _isolation = 0;
   _miniisolation_0p2 = 0;
   _miniisolation_0p3 = 0;
    _miniisolationcharged_0p2=0;
    _miniisoneutral=0;
    _trackSelectionMultiplicity=0;
    _muonSegmentComp=0;
   _multiisolation_T = 0;
   _multiisolation_M = 0;
   _multiisolation_L = 0;
   _pfisocharged = 0;
   _pfisoneutral = 0;
   _pfisophoton = 0;
   _isolationMCraw = 0;
   _isolationMCnonu = 0;
   _isolationMCdr03 = 0;
   _isolationMCdr03nonu = 0;
   _ptrel = 0;
   _ptrel2 = 0;
   _ptratio = 0;
   _mvaValue = 0;
   _mt = 0;
   _mllZ = 0;
   _mllG = 0;
   _mll = 0;
   _mllv = 0;
   _origin = 0;
   _originReduced = 0;
   _ipPV = 0;
   _ipPVerr = 0;
   _ipPVmc = 0;
   _ipZPV = 0;
   _ipZPVerr = 0;
   _3dIP = 0;
   _3dIPerr = 0;
   _3dIPsig = 0;
   _closeJetPtAll = 0;
   _closeJetEtaAll = 0;
   _closeJetPhiAll = 0;
   _closeJetMAll = 0;
   _closeJetCSVAll = 0;
   _closeJetNconstAll = 0;
   _closeJetAngAll = 0;
   _ptRelAll = 0;
   _closeJetPtAllMC = 0;
   _closeJetPtAllstatus = 0;
   _partonIdMatched = 0;
   _sameParton = 0;
   _isloose = 0;
   _istightIDWP2016_RA7 = 0;
   _isFOIDWP2016_RA7 = 0;
   _istightIDWP2016 = 0;
   _isFOIDWP2016 = 0;
   _isVetoIDWP2016 = 0;
   _isVetoIDWP2016_EWK = 0;
   _isFOIDWP2016_EWK = 0;
   _istightIDWP2016_EWK = 0;
   _istrigemulID=0;
   _istrigemulISO=0;
   _istight = 0;
   _isMediumMuon = 0;
   _istightIso = 0;
   _istightID = 0;
   _mompt = 0;
   _momphi = 0;
   _mometa = 0;
   _mompdg = 0;
    _lmva = 0;
   _lPt = 0;
   _lEta = 0;
   _lPhi = 0;
   _lE = 0;
   _lPtmc = 0;
   _lEtamc = 0;
   _lPhimc = 0;
   _lEmc = 0;
   _nuPtmc = 0;
   _nuEtamc = 0;
   _nuPhimc = 0;
   _nuEmc = 0;
   _mtmc = 0;
    _genWPhi=0;
    _genZPhi=0;
    _genWEta=0;
    _genZEta=0;
    _genWPt=0;
    _genZPt=0;
    _genWE=0;
    _genZE=0;
   _lnmisshits = 0;
   _lchargeGSF = 0;
   _lchargeCTF = 0;
   _lchargePixSc = 0;
   _bTagged = 0;
   _jetEta = 0;
   _jetPhi = 0;
   _jetPt = 0;
   _jetM = 0;
   _csv = 0;
   _jetDeltaRloose = 0;
    
    Prescale_passHLT_Mu8=0;
    Prescale_passHLT_Mu17=0;
    Prescale_passHLT_Mu8_TrkIsoVVL=0;
    Prescale_passHLT_Mu17_TrkIsoVVL=0;
    Prescale_passHLT_Ele12_CaloIdM_TrackIdM_PFJet30=0;
    Prescale_passHLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30=0;
    Prescale_passHLT_Ele17_CaloIdM_TrackIdM_PFJet30=0;
    Prescale_passHLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30=0;
    
    _tau_dz =0;
    _TauIs_againstMuonLoose3=0;
    _TauIs_byLooseCombinedIsolationDeltaBetaCorr3Hits=0;
    _TauIs_decayModeFindingNewDMs=0;
    _TauIs_decayModeFinding=0;
    _TauIs_byLooseIsolationMVArun2v1DBdR03oldDMwLT=0;
    _TauIs_byMediumIsolationMVArun2v1DBdR03oldDMwLT=0;
    _TauIs_byVLooseIsolationMVArun2v1DBoldDMwLT=0;
    _TauIs_byTightIsolationMVArun2v1DBoldDMwLT=0;
    
    _jetJECuncty=0;
    _jetbtagSF=0;
    _jetbtagSF_up=0;
    _jetbtagSF_down=0;
    _jetbtagEff=0;
    
    
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("_eventNb", &_eventNb, &b__eventNb);
   fChain->SetBranchAddress("_runNb", &_runNb, &b__runNb);
   fChain->SetBranchAddress("_lumiBlock", &_lumiBlock, &b__lumiBlock);
   fChain->SetBranchAddress("passmetfilters", &passmetfilters, &b_passmetfilters);
   fChain->SetBranchAddress("passHLT_DoubleMu8_Mass8_PFHT300", &passHLT_DoubleMu8_Mass8_PFHT300, &b_passHLT_DoubleMu8_Mass8_PFHT300);
   fChain->SetBranchAddress("passHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL", &passHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL, &b_passHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL);
   fChain->SetBranchAddress("passHLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL", &passHLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL, &b_passHLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL);
   fChain->SetBranchAddress("passHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ", &passHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ, &b_passHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ);
   fChain->SetBranchAddress("passHLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ", &passHLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ, &b_passHLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ);
   fChain->SetBranchAddress("passHLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT300", &passHLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT300, &b_passHLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT300);
   fChain->SetBranchAddress("passHLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL", &passHLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL, &b_passHLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL);
   fChain->SetBranchAddress("passHLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", &passHLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ, &b_passHLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ);
   fChain->SetBranchAddress("passHLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT300", &passHLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT300, &b_passHLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT300);
   fChain->SetBranchAddress("passHLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL", &passHLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL, &b_passHLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL);
   fChain->SetBranchAddress("passHLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL", &passHLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL, &b_passHLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL);
   fChain->SetBranchAddress("passHLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL", &passHLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL, &b_passHLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL);
   fChain->SetBranchAddress("passHLT_IsoMu24_eta2p1", &passHLT_IsoMu24_eta2p1, &b_passHLT_IsoMu24_eta2p1);
   fChain->SetBranchAddress("passHLT_Ele27_eta2p1_WPLoose_Gsf", &passHLT_Ele27_eta2p1_WPLoose_Gsf, &b_passHLT_Ele27_eta2p1_WPLoose_Gsf);
   fChain->SetBranchAddress("passHLT_Ele32_eta2p1_WPLoose_Gsf", &passHLT_Ele32_eta2p1_WPLoose_Gsf, &b_passHLT_Ele32_eta2p1_WPLoose_Gsf);
   fChain->SetBranchAddress("passHLT_Mu8", &passHLT_Mu8, &b_passHLT_Mu8);
   fChain->SetBranchAddress("passHLT_Mu17", &passHLT_Mu17, &b_passHLT_Mu17);
   fChain->SetBranchAddress("passHLT_Mu8_TrkIsoVVL", &passHLT_Mu8_TrkIsoVVL, &b_passHLT_Mu8_TrkIsoVVL);
   fChain->SetBranchAddress("passHLT_Mu17_TrkIsoVVL", &passHLT_Mu17_TrkIsoVVL, &b_passHLT_Mu17_TrkIsoVVL);
   fChain->SetBranchAddress("passHLT_Ele12_CaloIdM_TrackIdM_PFJet30", &passHLT_Ele12_CaloIdM_TrackIdM_PFJet30, &b_passHLT_Ele12_CaloIdM_TrackIdM_PFJet30);
   fChain->SetBranchAddress("passHLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30", &passHLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30, &b_passHLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30);
    fChain->SetBranchAddress("passHLT_Ele17_CaloIdM_TrackIdM_PFJet30", &passHLT_Ele17_CaloIdM_TrackIdM_PFJet30, &b_passHLT_Ele17_CaloIdM_TrackIdM_PFJet30);
    fChain->SetBranchAddress("passHLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30", &passHLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30, &b_passHLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30);
    fChain->SetBranchAddress("passHLT_Ele27_WPTight_Gsf", &passHLT_Ele27_WPTight_Gsf, &b_passHLT_Ele27_WPTight_Gsf);
    fChain->SetBranchAddress("passHLT_IsoMu24", &passHLT_IsoMu24, &b_passHLT_IsoMu24);
    fChain->SetBranchAddress("passHLT_IsoTkMu24", &passHLT_IsoTkMu24, &b_passHLT_IsoTkMu24);
    fChain->SetBranchAddress("passHLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ", &passHLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ, &b_passHLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ);
    fChain->SetBranchAddress("passHLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ", &passHLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ, &b_passHLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ);
    fChain->SetBranchAddress("passHLT_Ele22_eta2p1_WPLoose_GSF_LooseIsoPFtau20_SingleL1", &passHLT_Ele22_eta2p1_WPLoose_GSF_LooseIsoPFtau20_SingleL1, &b_passHLT_Ele22_eta2p1_WPLoose_GSF_LooseIsoPFtau20_SingleL1);
    fChain->SetBranchAddress("passHLT_Ele24_eta2p1_WPLoose_GSF_LooseIsoPFtau30", &passHLT_Ele24_eta2p1_WPLoose_GSF_LooseIsoPFtau30, &b_passHLT_Ele24_eta2p1_WPLoose_GSF_LooseIsoPFtau30);
    fChain->SetBranchAddress("passHLT_Ele36_eta2p1_WPLoose_GSF_LooseIsoPFtau20_SingleL1", &passHLT_Ele36_eta2p1_WPLoose_GSF_LooseIsoPFtau20_SingleL1, &b_passHLT_Ele36_eta2p1_WPLoose_GSF_LooseIsoPFtau20_SingleL1);
    fChain->SetBranchAddress("passHLT_IsoMu19_eta2p1_LooseIsoPFTau20_SingleL1_v", &passHLT_IsoMu19_eta2p1_LooseIsoPFTau20_SingleL1_v, &b_passHLT_IsoMu19_eta2p1_LooseIsoPFTau20_SingleL1_v);
    fChain->SetBranchAddress("Prescale_passHLT_Mu8", &Prescale_passHLT_Mu8, &b_Prescale_passHLT_Mu8);
    fChain->SetBranchAddress("Prescale_passHLT_Mu17", &Prescale_passHLT_Mu17, &b_Prescale_passHLT_Mu17);
    fChain->SetBranchAddress("Prescale_passHLT_Mu8_TrkIsoVVL", &Prescale_passHLT_Mu8_TrkIsoVVL, &b_Prescale_passHLT_Mu8_TrkIsoVVL);
    fChain->SetBranchAddress("Prescale_passHLT_Mu17_TrkIsoVVL", &Prescale_passHLT_Mu17_TrkIsoVVL, &b_Prescale_passHLT_Mu17_TrkIsoVVL);
    fChain->SetBranchAddress("Prescale_passHLT_Ele12_CaloIdM_TrackIdM_PFJet30", &Prescale_passHLT_Ele12_CaloIdM_TrackIdM_PFJet30, &b_Prescale_passHLT_Ele12_CaloIdM_TrackIdM_PFJet30);
    fChain->SetBranchAddress("Prescale_passHLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30", &Prescale_passHLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30, &b_Prescale_passHLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30);
    fChain->SetBranchAddress("Prescale_passHLT_Ele17_CaloIdM_TrackIdM_PFJet30", &Prescale_passHLT_Ele17_CaloIdM_TrackIdM_PFJet30, &b_Prescale_passHLT_Ele17_CaloIdM_TrackIdM_PFJet30);
    fChain->SetBranchAddress("Prescale_passHLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30", &Prescale_passHLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30, &b_Prescale_passHLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30);
   fChain->SetBranchAddress("_weight", &_weight, &b__weight);
   fChain->SetBranchAddress("_genHT", &_genHT, &b__genHT);
   fChain->SetBranchAddress("_nLeptons", &_nLeptons, &b__nLeptons);
   fChain->SetBranchAddress("_index1", &_index1, &b__index1);
   fChain->SetBranchAddress("_index2", &_index2, &b__index2);
   fChain->SetBranchAddress("_sb", &_sb, &b__sb);
   fChain->SetBranchAddress("_doubleF", &_doubleF, &b__doubleF);
   fChain->SetBranchAddress("_PVchi2", &_PVchi2, &b__PVchi2);
   fChain->SetBranchAddress("_PVerr", _PVerr, &b__PVerr);
   fChain->SetBranchAddress("_n_PV", &_n_PV, &b__n_PV);
   fChain->SetBranchAddress("_met", &_met, &b__met);
   fChain->SetBranchAddress("_met_phi", &_met_phi, &b__met_phi);
   fChain->SetBranchAddress("_met_JECup", &_met_JECup, &b__met_JECup);
   fChain->SetBranchAddress("_met_JECdown", &_met_JECdown, &b__met_JECdown);
   fChain->SetBranchAddress("_met_phi_JECup", &_met_phi_JECup, &b__met_phi_JECup);
   fChain->SetBranchAddress("_met_phi_JECdown", &_met_phi_JECdown, &b__met_phi_JECdown);
   fChain->SetBranchAddress("HT", &HT, &b_HT);
   fChain->SetBranchAddress("_genmet", &_genmet, &b__genmet);
   fChain->SetBranchAddress("_genmet_phi", &_genmet_phi, &b__genmet_phi);
   fChain->SetBranchAddress("_n_bJets", &_n_bJets, &b__n_bJets);
   fChain->SetBranchAddress("_n_Jets", &_n_Jets, &b__n_Jets);
   fChain->SetBranchAddress("_n_Jets40", &_n_Jets40, &b__n_Jets40);
   fChain->SetBranchAddress("trueNVtx", &trueNVtx, &b_trueNVtx);
   fChain->SetBranchAddress("_closeIndex", &_closeIndex, &b__closeIndex);
   fChain->SetBranchAddress("_indeces", &_indeces, &b__indeces);
   fChain->SetBranchAddress("_flavors", &_flavors, &b__flavors);
   fChain->SetBranchAddress("_charges", &_charges, &b__charges);
   fChain->SetBranchAddress("_chargesMC", &_chargesMC, &b__chargesMC);
   fChain->SetBranchAddress("_isolation", &_isolation, &b__isolation);
   fChain->SetBranchAddress("_miniisolation_0p2", &_miniisolation_0p2, &b__miniisolation_0p2);
   fChain->SetBranchAddress("_miniisolation_0p3", &_miniisolation_0p3, &b__miniisolation_0p3);
   fChain->SetBranchAddress("_miniisolationcharged_0p2", &_miniisolationcharged_0p2, &b__miniisolationcharged_0p2);
   fChain->SetBranchAddress("_miniisoneutral", &_miniisoneutral, &b__miniisoneutral);
   fChain->SetBranchAddress("_trackSelectionMultiplicity", &_trackSelectionMultiplicity, &b__trackSelectionMultiplicity);
   fChain->SetBranchAddress("_muonSegmentComp", &_muonSegmentComp, &b__muonSegmentComp);
   fChain->SetBranchAddress("_multiisolation_T", &_multiisolation_T, &b__multiisolation_T);
   fChain->SetBranchAddress("_multiisolation_M", &_multiisolation_M, &b__multiisolation_M);
   fChain->SetBranchAddress("_multiisolation_L", &_multiisolation_L, &b__multiisolation_L);
   fChain->SetBranchAddress("_pfisocharged", &_pfisocharged, &b__pfisocharged);
   fChain->SetBranchAddress("_pfisoneutral", &_pfisoneutral, &b__pfisoneutral);
   fChain->SetBranchAddress("_pfisophoton", &_pfisophoton, &b__pfisophoton);
   fChain->SetBranchAddress("_isolationMCraw", &_isolationMCraw, &b__isolationMCraw);
   fChain->SetBranchAddress("_isolationMCnonu", &_isolationMCnonu, &b__isolationMCnonu);
   fChain->SetBranchAddress("_isolationMCdr03", &_isolationMCdr03, &b__isolationMCdr03);
   fChain->SetBranchAddress("_isolationMCdr03nonu", &_isolationMCdr03nonu, &b__isolationMCdr03nonu);
   fChain->SetBranchAddress("_ptrel", &_ptrel, &b__ptrel);
   fChain->SetBranchAddress("_ptrel2", &_ptrel2, &b__ptrel2);
   fChain->SetBranchAddress("_ptratio", &_ptratio, &b__ptratio);
   fChain->SetBranchAddress("_mvaValue", &_mvaValue, &b__mvaValue);
   fChain->SetBranchAddress("_mt", &_mt, &b__mt);
   fChain->SetBranchAddress("_mllZ", &_mllZ, &b__mllZ);
   fChain->SetBranchAddress("_mllG", &_mllG, &b__mllG);
   fChain->SetBranchAddress("_mll", &_mll, &b__mll);
   fChain->SetBranchAddress("_mllv", &_mllv, &b__mllv);
   fChain->SetBranchAddress("_origin", &_origin, &b__origin);
   fChain->SetBranchAddress("_originReduced", &_originReduced, &b__originReduced);
   fChain->SetBranchAddress("_ipPV", &_ipPV, &b__ipPV);
   fChain->SetBranchAddress("_ipPVerr", &_ipPVerr, &b__ipPVerr);
   fChain->SetBranchAddress("_ipPVmc", &_ipPVmc, &b__ipPVmc);
   fChain->SetBranchAddress("_ipZPV", &_ipZPV, &b__ipZPV);
   fChain->SetBranchAddress("_ipZPVerr", &_ipZPVerr, &b__ipZPVerr);
   fChain->SetBranchAddress("_3dIP", &_3dIP, &b__3dIP);
   fChain->SetBranchAddress("_3dIPerr", &_3dIPerr, &b__3dIPerr);
   fChain->SetBranchAddress("_3dIPsig", &_3dIPsig, &b__3dIPsig);
   fChain->SetBranchAddress("_closeJetPtAll", &_closeJetPtAll, &b__closeJetPtAll);
   fChain->SetBranchAddress("_closeJetEtaAll", &_closeJetEtaAll, &b__closeJetEtaAll);
   fChain->SetBranchAddress("_closeJetPhiAll", &_closeJetPhiAll, &b__closeJetPhiAll);
   fChain->SetBranchAddress("_closeJetMAll", &_closeJetMAll, &b__closeJetMAll);
   fChain->SetBranchAddress("_closeJetCSVAll", &_closeJetCSVAll, &b__closeJetCSVAll);
   fChain->SetBranchAddress("_closeJetNconstAll", &_closeJetNconstAll, &b__closeJetNconstAll);
   fChain->SetBranchAddress("_closeJetAngAll", &_closeJetAngAll, &b__closeJetAngAll);
   fChain->SetBranchAddress("_ptRelAll", &_ptRelAll, &b__ptRelAll);
   fChain->SetBranchAddress("_closeJetPtAllMC", &_closeJetPtAllMC, &b__closeJetPtAllMC);
   fChain->SetBranchAddress("_closeJetPtAllstatus", &_closeJetPtAllstatus, &b__closeJetPtAllstatus);
   fChain->SetBranchAddress("_partonIdMatched", &_partonIdMatched, &b__partonIdMatched);
   fChain->SetBranchAddress("_sameParton", &_sameParton, &b__sameParton);
   fChain->SetBranchAddress("_isloose", &_isloose, &b__isloose);
   fChain->SetBranchAddress("_istight", &_istight, &b__istight);
   fChain->SetBranchAddress("_isMediumMuon", &_isMediumMuon, &b__isMediumMuon);
   fChain->SetBranchAddress("_istightIso", &_istightIso, &b__istightIso);
   fChain->SetBranchAddress("_istightID", &_istightID, &b__istightID);
   fChain->SetBranchAddress("_istightIDWP2016_RA7", &_istightIDWP2016_RA7, &b__istightIDWP2016_RA7);
   fChain->SetBranchAddress("_isFOIDWP2016_RA7", &_isFOIDWP2016_RA7, &b__isFOIDWP2016_RA7);
   fChain->SetBranchAddress("_isVetoIDWP2016", &_isVetoIDWP2016, &b__isVetoIDWP2016);
   fChain->SetBranchAddress("_istightIDWP2016", &_istightIDWP2016, &b__istightIDWP2016);
   fChain->SetBranchAddress("_isFOIDWP2016", &_isFOIDWP2016, &b__isFOIDWP2016);
   fChain->SetBranchAddress("_isVetoIDWP2016_EWK", &_isVetoIDWP2016_EWK, &b__isVetoIDWP2016_EWK);
   fChain->SetBranchAddress("_isFOIDWP2016_EWK", &_isFOIDWP2016_EWK, &b__isFOIDWP2016_EWK);
   fChain->SetBranchAddress("_istightIDWP2016_EWK", &_istightIDWP2016_EWK, &b__istightIDWP2016_EWK);
   fChain->SetBranchAddress("_istrigemulID", &_istrigemulID, &b__istrigemulID);
   fChain->SetBranchAddress("_istrigemulISO", &_istrigemulISO, &b__istrigemulISO);
   fChain->SetBranchAddress("_mompt", &_mompt, &b__mompt);
   fChain->SetBranchAddress("_momphi", &_momphi, &b__momphi);
   fChain->SetBranchAddress("_mometa", &_mometa, &b__mometa);
   fChain->SetBranchAddress("_mompdg", &_mompdg, &b__mompdg);
   fChain->SetBranchAddress("_lmva", &_lmva, &b__lmva);
   fChain->SetBranchAddress("_lPt", &_lPt, &b__lPt);
   fChain->SetBranchAddress("_lEta", &_lEta, &b__lEta);
   fChain->SetBranchAddress("_lPhi", &_lPhi, &b__lPhi);
   fChain->SetBranchAddress("_lE", &_lE, &b__lE);
   fChain->SetBranchAddress("_lPtmc", &_lPtmc, &b__lPtmc);
   fChain->SetBranchAddress("_lEtamc", &_lEtamc, &b__lEtamc);
   fChain->SetBranchAddress("_lPhimc", &_lPhimc, &b__lPhimc);
   fChain->SetBranchAddress("_lEmc", &_lEmc, &b__lEmc);
   fChain->SetBranchAddress("_nuPtmc", &_nuPtmc, &b__nuPtmc);
   fChain->SetBranchAddress("_nuEtamc", &_nuEtamc, &b__nuEtamc);
   fChain->SetBranchAddress("_nuPhimc", &_nuPhimc, &b__nuPhimc);
   fChain->SetBranchAddress("_nuEmc", &_nuEmc, &b__nuEmc);
   fChain->SetBranchAddress("_mtmc", &_mtmc, &b__mtmc);
    fChain->SetBranchAddress("_genWPhi", &_genWPhi, &b__genWPhi);
    fChain->SetBranchAddress("_genZPhi", &_genZPhi, &b__genZPhi);
    fChain->SetBranchAddress("_genWEta", &_genWEta, &b__genWEta);
    fChain->SetBranchAddress("_genZEta", &_genZEta, &b__genZEta);
    fChain->SetBranchAddress("_genWPt", &_genWPt, &b__genWPt);
    fChain->SetBranchAddress("_genZPt", &_genZPt, &b__genZPt);
    fChain->SetBranchAddress("_genWE", &_genWE, &b__genWE);
    fChain->SetBranchAddress("_genZE", &_genZE, &b__genZE);
   fChain->SetBranchAddress("_lnmisshits", &_lnmisshits, &b__lnmisshits);
   fChain->SetBranchAddress("_lchargeGSF", &_lchargeGSF, &b__lchargeGSF);
   fChain->SetBranchAddress("_lchargeCTF", &_lchargeCTF, &b__lchargeCTF);
   fChain->SetBranchAddress("_lchargePixSc", &_lchargePixSc, &b__lchargePixSc);
   fChain->SetBranchAddress("_bTagged", &_bTagged, &b__bTagged);
   fChain->SetBranchAddress("_jetEta", &_jetEta, &b__jetEta);
   fChain->SetBranchAddress("_jetPhi", &_jetPhi, &b__jetPhi);
   fChain->SetBranchAddress("_jetPt", &_jetPt, &b__jetPt);
   fChain->SetBranchAddress("_jetM", &_jetM, &b__jetM);
   fChain->SetBranchAddress("_csv", &_csv, &b__csv);
   fChain->SetBranchAddress("_jetDeltaRloose", &_jetDeltaRloose, &b__jetDeltaRloose);
    
    fChain->SetBranchAddress("_tau_dz", &_tau_dz, &b__tau_dz);
    fChain->SetBranchAddress("_TauIs_againstMuonLoose3", &_TauIs_againstMuonLoose3, &b__TauIs_againstMuonLoose3);
    fChain->SetBranchAddress("_TauIs_byLooseCombinedIsolationDeltaBetaCorr3Hits", &_TauIs_byLooseCombinedIsolationDeltaBetaCorr3Hits, &b__TauIs_byLooseCombinedIsolationDeltaBetaCorr3Hits);
    fChain->SetBranchAddress("_TauIs_decayModeFindingNewDMs", &_TauIs_decayModeFindingNewDMs, &b__TauIs_decayModeFindingNewDMs);
    fChain->SetBranchAddress("_TauIs_decayModeFinding", &_TauIs_decayModeFinding, &b__TauIs_decayModeFinding);
    fChain->SetBranchAddress("_TauIs_byLooseIsolationMVArun2v1DBdR03oldDMwLT", &_TauIs_byLooseIsolationMVArun2v1DBdR03oldDMwLT, &b__TauIs_byLooseIsolationMVArun2v1DBdR03oldDMwLT);
    fChain->SetBranchAddress("_TauIs_byMediumIsolationMVArun2v1DBdR03oldDMwLT", &_TauIs_byMediumIsolationMVArun2v1DBdR03oldDMwLT, &b__TauIs_byMediumIsolationMVArun2v1DBdR03oldDMwLT);
    fChain->SetBranchAddress("_TauIs_byVLooseIsolationMVArun2v1DBoldDMwLT", &_TauIs_byVLooseIsolationMVArun2v1DBoldDMwLT, &b__TauIs_byVLooseIsolationMVArun2v1DBoldDMwLT);
    fChain->SetBranchAddress("_TauIs_byTightIsolationMVArun2v1DBoldDMwLT", &_TauIs_byTightIsolationMVArun2v1DBoldDMwLT, &b__TauIs_byTightIsolationMVArun2v1DBoldDMwLT);
    
    fChain->SetBranchAddress("_jetJECuncty", &_jetJECuncty, &b__jetJECuncty);
    fChain->SetBranchAddress("_jetbtagSF", &_jetbtagSF, &b__jetbtagSF);
    fChain->SetBranchAddress("_jetbtagSF_up", &_jetbtagSF_up, &b__jetbtagSF_up);
    fChain->SetBranchAddress("_jetbtagSF_down", &_jetbtagSF_down, &b__jetbtagSF_down);
    fChain->SetBranchAddress("_jetbtagEff", &_jetbtagEff, &b__jetbtagEff);
    
   Notify();
}

Bool_t FinalYield::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void FinalYield::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t FinalYield::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef FinalYield_cxx
