#ifndef SSb13_takeimai_H
#define SSb13_takeimai_H

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Conversion.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/RecoCandidate/interface/IsoDepositFwd.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"

#include "SUSYAnalyzer/PatAnalyzer/interface/GenParticleManager.h"
#include "SUSYAnalyzer/PatAnalyzer/interface/Statistics.h"
#include "SUSYAnalyzer/PatAnalyzer/interface/Tools.h"
#include "SUSYAnalyzer/PatAnalyzer/interface/OnTheFlyCorrections.hh"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "EgammaAnalysis/ElectronTools/interface/EGammaMvaEleEstimatorCSA14.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "RecoTauTag/RecoTau/interface/RecoTauVertexAssociator.h"
//#include "DataFormats/Common/interface/TriggerResults.h"
#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondTools/BTau/interface/BTagCalibrationReader.h"

#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
//Root Classes
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TH1I.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TTree.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TString.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TLegend.h"
#include "TClonesArray.h"

//Standard C++ classes
#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <utility>
#include <ostream>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <memory>
#include <iomanip>

using namespace std;

const int nLeptonsMax = 10;
const int nJetsMax = 30;

class SSb13_takeimai : public edm::EDAnalyzer {
public:
    
    explicit SSb13_takeimai(const edm::ParameterSet & iConfig);
    ~SSb13_takeimai(){};
    
private:
    
    //virtual void analyze(edm::Event & iEvent, const edm::EventSetup & iSetup);
    virtual void analyze(const edm::Event&, const edm::EventSetup&);
    virtual void beginJob();
    virtual void endJob(void);
    
    void fillRegVars(const pat::Jet *jet, double genpt, const pat::Muon* mu);
    void fillRegVars(const pat::Jet *jet, double genpt, const pat::Electron* el);
    
    void fillMCVars(const GenParticle* mc, const int leptonCounter);
    void fillCloseJetVars(const int leptonCounter, Vertex::Point PV);
    void matchCloseJet(const int leptonCounter);
    void fillIsoMCVars(const int leptonCounter);
    
    bool PassTriggerLeg(std::string triggerlegstring, std::string triggerlegstringalt,const pat::Muon *muonit, const edm::Event&);
    bool PassTriggerLeg(std::string triggerlegstring, std::string triggerlegstringalt,const pat::Electron *eleit, const edm::Event&);
    bool PassTriggerLeg(std::string triggerlegstring, const pat::Muon *muonit, const edm::Event& theevent){return PassTriggerLeg(triggerlegstring,"Noalttrigger",muonit,theevent);};
    bool PassTriggerLeg(std::string triggerlegstring, const pat::Electron *eleit, const edm::Event& theevent){ return PassTriggerLeg(triggerlegstring,"Noalttrigger",eleit,theevent);};
    bool ApplySkim(const edm::Event& ,std::string skim);

    int GenParticleInfo(const edm::Event& iEvent);
    void MatchToGenParticle(const edm::Event& iEvent, const int & leptonctr, const double &recolpt, const double &recoleta, const double &recolphi, const int &recolpdgid, vector <bool> & checkusedpar);
    void bookTree(), ClearStuff(),InitLepton();
    
    std::vector<const pat::Jet* > SelectedJetsAll;
    edm::Handle<GenParticleCollection> TheGenParticles;

    Vertex::Point PVmc;
    
    bool IsMC, IsFS;
    bool isgetbylabel,usemva;
    std::string Sample;
    std::string Skim;

    edm::EDGetTokenT<vector<PileupSummaryInfo> > puInfo_token; 
    edm::EDGetTokenT<GenEventInfoProduct> geninfo_token; 
    edm::EDGetTokenT<GenParticleCollection> genpart_token; 
    edm::EDGetTokenT<edm::TriggerResults> trgresults_token; 
    edm::EDGetTokenT<reco::BeamSpot> beamspot_token; 
    edm::EDGetTokenT<std::vector<Vertex> > goodpv_token; 
    edm::EDGetTokenT<std::vector<pat::MET> > met_token; 
    edm::EDGetTokenT<pat::PackedCandidateCollection> pfcands_token; 
    edm::EDGetTokenT<std::vector<pat::Muon> > patmuon_token; 
    edm::EDGetTokenT<std::vector<pat::Electron> > patelectron_token; 
    edm::EDGetTokenT<edm::View<pat::Electron> > patelectronview_token;
    edm::EDGetTokenT<std::vector<pat::Tau> > pattau_token;
    edm::EDGetTokenT<std::vector<reco::Conversion> > conv_token; 
    edm::EDGetTokenT<std::vector< pat::Jet> > jet_token; 
    edm::EDGetTokenT<double> rhoJets_token; 
    edm::EDGetTokenT<double> rhoJetsNC_token;
    edm::EDGetTokenT<l1extra::L1EmParticleCollection> emIsolColl_token; 
    edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> trigobject_token; 
    edm::EDGetTokenT<pat::PackedTriggerPrescales> trgpresc_token;
    edm::EDGetTokenT<pat::PackedTriggerPrescales> L1presc_token;
    edm::EDGetTokenT<edm::TriggerResults> metfilters_token ; 
    edm::EDGetTokenT<bool>filterbadChCand_token ;
    edm::EDGetTokenT<bool>filterbadPFMuon_token ;

    
    
    edm::InputTag IT_muon;
    edm::InputTag IT_electron;
    edm::InputTag IT_tau;
    edm::InputTag IT_tauDiscriminator;
    edm::InputTag IT_jet;
    edm::InputTag IT_pfmet;
    edm::InputTag IT_pileup;
    edm::InputTag IT_beamspot;
    edm::InputTag IT_hltresults;
    edm::InputTag IT_METFilters;
    edm::InputTag IT_gen;
    edm::InputTag IT_genpart;
    edm::InputTag IT_goodVtx;
    edm::InputTag IT_pfcands;
    edm::InputTag IT_conv;
    edm::InputTag IT_triggerObjects;
    edm::InputTag IT_triggerPrescales;
    edm::InputTag IT_L1Prescales;


    edm::Service<TFileService> fs;
    FILE *outfile;
    
    TH1F *Nvtx;

    BTagCalibration *calib;
    BTagCalibrationReader *reader;
    BTagCalibrationReader *reader_up;
    BTagCalibrationReader *reader_do;
    TH2F *btagEffs2D[3];

    TMVA::Reader *tmvareaderele, *tmvareadermu;
    float LepGood_pt_forMVA, LepGood_eta_forMVA;
    float LepGood_jetNDauChargedMVASel_forMVA; 
    float LepGood_miniRelIsoCharged_forMVA, LepGood_miniRelIsoNeutral_forMVA, LepGood_jetPtRelv2_forMVA,LepGood_jetPtRatio_forMVA, LepGood_jetBTagCSV_forMVA;
    float LepGood_sip3d_forMVA,LepGood_dxy_forMVA,LepGood_dz_forMVA,LepGood_mvaIdSpring16GP_forMVA;
    float LepGood_segmentCompatibility_forMVA;
    //desired output variables
    TTree* outputTree;
    
    string _corrLevel;
    
    
    double _relIsoCutE;
    double _relIsoCutMu;
    double _relIsoCutEloose;
    double _relIsoCutMuloose;
    
    bool _chargeConsistency;
    bool _useconversions;
    
    double _tightD0Mu;
    double _tightD0E;
    double _looseD0Mu;
    double _looseD0E;
    
    double _jetPtCut;
    double _jetEtaCut;

    double _eleMinPt;
    double _muonMinPt;
    double _tauMinPt;
    
    double _tauPt;
    double _tauEta;
    
    bool _regression;
    
    // MVA values and categories (optional)
    edm::EDGetTokenT<edm::ValueMap<float> > mvaValuesMapToken_;
    edm::EDGetTokenT<edm::ValueMap<int> > mvaCategoriesMapToken_;
    double _miniisocut[2];
    double _ptratiocut[2];
    double _ptrelcut[2];
    double _multiConst[3][3];


    
    /*    std::vector<std::string> myManualCatWeigths;
    vector<string> myManualCatWeigthsTrig;
    EGammaMvaEleEstimatorCSA14* myMVATrig;*/
    //double looseMVA[3][2][3]; //Eta:0.8-1.479-2.5//Loose/Tight//Pt:0-10-20
    double looseMVA[3][2];

    double myRhoJetsNC;
    double myRhoJets;

    //genlevel particles
    GenParticleManager GPM;
    OnTheFlyCorrections* fMetCorrector;
    
    int _n_bJets;
    int _n_Jets;
    int _n_Jets40;
    
    
    vector<double> _jetEta;
    vector<double>  _jetPhi;
    vector<double>  _jetPt;
    vector<double>  _jetM;
    vector<double>  _bTagged;
    vector<double>  _csv;
    vector<double>  _jetDeltaRloose;
    vector<int> _jetpFlav;
    vector<int> _jethFlav;
    vector<double>  _jetJECuncty ;
    vector<double> _jetbtagSF;
    vector<double> _jetbtagSF_up;
    vector<double> _jetbtagSF_down;
    vector<double> _jetbtagEff;


    //    vector<double>  ;
    /*
    double _jetEta[nJetsMax];
    double _jetPhi[nJetsMax];
    double _jetPt[nJetsMax];
    double _jetM[nJetsMax];
    bool _bTagged[nJetsMax];
    double _csv[nJetsMax];
    double _jetDeltaR[nJetsMax][nLeptonsMax];
    double _jetDeltaRloose[nJetsMax];
    */
    int _leptonIndex;

    
    
    TClonesArray* _leptonP4;
    TClonesArray* _jetP4;
    TClonesArray* _jetAllP4;
    
    int _nLeptons;
    int _nEle;
    int _nMu;
    int _nTau;
    double _weight;
    double _genHT;
    
    int _eventType; //ee,mm,em
    bool _sb;
    bool _doubleF;
    int _index1 = -1;
    int _index2 = -1;

    int nbofL1EG =0;
    /*   int _closeIndex[nLeptonsMax];    
     int _indeces[nLeptonsMax];
    int _flavors[nLeptonsMax];
    double _charges[nLeptonsMax];
    double _chargesMC[nLeptonsMax];
    double _isolation[nLeptonsMax];
    double _miniisolation[nLeptonsMax][2];
    bool _multiisolation[nLeptonsMax][3];//3 WP
    double _isolationComponents[nLeptonsMax][4];
    double _isolationMC[nLeptonsMax][4]; //all daughters; all visible daughters; all daughters in the cone; all visible daughters in the cone
    double _ptrel[nLeptonsMax];
    double _ptrel2[nLeptonsMax];
    double _ptratio[nLeptonsMax];
    double _mvaValue[nLeptonsMax];
    double _mt[nLeptonsMax];
    double _mllZ[nLeptonsMax];
    double _mllG[nLeptonsMax];
    double _mll[nLeptonsMax][nLeptonsMax];
    vector <vector<double>> _mllv ;
    int _origin[nLeptonsMax];
    int _originReduced[nLeptonsMax];
    double _PVchi2;
    double _PVerr[3];
    double _ipPV[nLeptonsMax];
    double _ipPVerr[nLeptonsMax];
    double _ipPVmc[nLeptonsMax];
    double _ipZPV[nLeptonsMax];
    double _ipZPVerr[nLeptonsMax];
    double _3dIP[nLeptonsMax];
    double _3dIPerr[nLeptonsMax];
    double _3dIPsig[nLeptonsMax];
    double _closeJetPtAll[nLeptonsMax];
    double _closeJetEtaAll[nLeptonsMax];
    double _closeJetPhiAll[nLeptonsMax];
    double _closeJetMAll[nLeptonsMax];
    double _closeJetCSVAll[nLeptonsMax];
    int _closeJetNconstAll[nLeptonsMax];
    double _closeJetAngAll[nLeptonsMax];
    double _ptRelAll[nLeptonsMax];
    double _closeJetPtAllMC[nLeptonsMax];
    double _closeJetPtAllstatus[nLeptonsMax];
    int _partonIdMatched[nLeptonsMax];
    bool _sameParton[nLeptonsMax];
    bool _isloose[nLeptonsMax];
    bool _istight[nLeptonsMax];
    bool _istightIso[nLeptonsMax];
    bool _istightID[nLeptonsMax];
    double _lPt[nLeptonsMax], _lEta[nLeptonsMax], _lPhi[nLeptonsMax], _lE[nLeptonsMax];
    double _lPtmc[nLeptonsMax], _lEtamc[nLeptonsMax], _lPhimc[nLeptonsMax], _lEmc[nLeptonsMax];
    double _nuPtmc[nLeptonsMax], _nuEtamc[nLeptonsMax], _nuPhimc[nLeptonsMax], _nuEmc[nLeptonsMax];
    double _mtmc[nLeptonsMax];
    double _mompt[nLeptonsMax];
    double _momphi[nLeptonsMax];
    double _mometa[nLeptonsMax];
    int _mompdg[nLeptonsMax];
    int  _lnmisshits[nLeptonsMax];
    int _lchargeGSF[nLeptonsMax], _lchargeCTF[nLeptonsMax], _lchargePixSc[nLeptonsMax];
    //double charges_ScPix[nLeptonsMax], charges_Ctf[nLeptonsMax],charges_Gsf[nLeptonsMax];
    */












    
    int _n_PV;
    int trueNVtx,nVtxNow,nVtxBefore,nVtxAfter;
    int _n_electrons;
    int _n_muons;
    int _n_taus;
    
    unsigned long _eventNb;
    unsigned long _runNb;
    unsigned long _lumiBlock;
    

    double _PVchi2;
    double _PVerr[3];

    bool passHLT_PFMET170_NotCleaned, passHLT_PFMET170_HBHECleaned, passHLT_PFMET170_BeamHaloCleaned, passHLT_PFMET120_Mu5,passHLT_Ele15_IsoVVVL_PFHT350_PFMET50;
    bool passHLT_IsoMu16_eta2p1_MET30,passHLT_IsoMu19_eta2p1_LooseIsoPFTau20_v,passHLT_Ele22_eta2p1_WPLoose_GSF_LooseIsoPFtau20_SingleL1_v,passHLT_Mu17_Mu8_SameSign_DZ,passHLT_DoubleMu3_PFMET50,passHLT_TripleMu_5_3_3;
    bool passHLT_Ele8_CaloIdM_TrackIdM_PFJet30,passHLT_Ele12_CaloIdM_TrackIdM_PFJet30,passHLT_Ele17_CaloIdM_TrackIdM_PFJet30,passHLT_Ele23_CaloIdM_TrackIdM_PFJet30;
    bool passHLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30,passHLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30,passHLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30,passHLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30;
    bool passHLT_Ele15_IsoVVVL_PFHT350;
    bool passHLT_Ele10_CaloIdM_TrackIdM_CentralPFJet30_BTagCSV_p13;
    bool passHLT_Mu10_CentralPFJet30_BTagCSV_p13 ;
    bool passHLT_Mu15_IsoVVVL_PFHT350;
    bool passHLT_Mu8, passHLT_Mu17, passHLT_Mu20, passHLT_Mu27;
    bool passHLT_Mu8_TrkIsoVVL,passHLT_Mu17_TrkIsoVVL;
    bool passcsctighthalo2015, passglobaltighthalo2016, passglobalsupertighthalo2016;
    bool passmetfilters;

    bool passHLT_TripleMu_12_10_5;
    bool passHLT_DiMu9_Ele9_CaloIdL_TrackIdL;
    bool passHLT_Mu8_DiEle12_CaloIdL_TrackIdL;
    bool passHLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL;

    bool passHLT_DoubleMu8_Mass8_PFHT300;
    bool passHLT_DoubleMu8_Mass8_PFHT250;
    bool passHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL;
    bool passHLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL;
    bool passHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ;
    bool passHLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ;

    bool passHLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT300;
    bool passHLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT250;
    bool passHLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL;
    bool passHLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ;
    bool passHLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL;
    bool passHLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ;

    bool passHLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT300;
    bool passHLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT250;
    bool passHLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL;
    bool passHLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL;
    bool passHLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL;
    bool passHLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL;
    bool passHLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL;
  
    bool passHLT_IsoMu24_eta2p1;
    bool passHLT_IsoMu18,passHLT_IsoMu20,passHLT_IsoMu22;
    bool passHLT_IsoTkMu18,passHLT_IsoTkMu20,passHLT_IsoTkMu22;

    bool passHLT_Ele23_WPLoose_Gsf,passHLT_Ele27_WPLoose_Gsf;
    bool passHLT_Ele27_eta2p1_WPLoose_Gsf;
    bool passHLT_Ele32_eta2p1_WPLoose_Gsf;

    bool passHLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW;
    bool passHLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL;
    
    bool passHLT_Ele23_CaloIdL_TrackIdL_IsoVL;


    /*
    bool hltL1sL1Mu6HTT150ORL1Mu8HTT125ORL1HTT125ORL1HTT150ORL1HTT175[nLeptonsMax];
    bool  hltDoubleMu8Mass8L3Filtered[nLeptonsMax];


    bool  passleg1L[nLeptonsMax];


    //bool passleg
    bool  passleg2LSFandHT[nLeptonsMax];
    bool  passleg2LSFIso[nLeptonsMax];
    bool  passleg2LOFandHT[nLeptonsMax];
    bool  passlegMu23Ele12Iso[nLeptonsMax];
    bool  passlegEle23Mu8Iso[nLeptonsMax];
    */
    
    //Prescale Factors (HLT*L1)
    double Prescale_passHLT_Mu8;
    double Prescale_passHLT_Mu17;
    double Prescale_passHLT_Mu8_TrkIsoVVL;
    double Prescale_passHLT_Mu17_TrkIsoVVL;
    double Prescale_passHLT_Ele12_CaloIdM_TrackIdM_PFJet30;
    double Prescale_passHLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30;
    double Prescale_passHLT_Ele17_CaloIdM_TrackIdM_PFJet30;
    double Prescale_passHLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30;

    
    vector < int >_closeIndex;
    vector < int >_indeces;
    vector < int >_flavors;
    vector < int >_charges;
    vector < int >_chargesMC;
    vector < double >_isolation;
    vector < double >_miniisolation_0p2;
    vector < double >_miniisolationcharged_0p2;
    vector < double >_miniisolation_0p3;
    vector < double > _miniisoneutral;
    vector < int > _trackSelectionMultiplicity;
    vector < double > _muonSegmentComp;
    vector < bool > _multiisolation_T;
    vector < bool > _multiisolation_M;        
    vector < bool > _multiisolation_L;
    
    vector < double >_lmvawithiso;
    vector < double >_lmva;
    vector < double >_pfisocharged;
    vector < double >_pfisoneutral;
    vector < double >_pfisophoton;
    vector < double >_isolationMCraw;
    vector < double >_isolationMCnonu;
    vector < double >_isolationMCdr03;
    vector < double >_isolationMCdr03nonu;
    vector < double >_ptrel;
    vector < double >_ptrel2;
    vector < double >_ptratio;
    vector < double >_ljetmultiplicity;
    vector < double >_mvaValue;
    vector < double >_mt;
    vector < double >_mllZ;
    vector < double >_mllG;
    vector < double >_mll;
    vector <vector<double>> _mllv ;
    vector < int >_origin;
    vector < int >_originReduced;
    vector < double >_ipPV;
    vector < double >_ipPVerr;
    vector < double >_ipPVmc;
    vector < double >_ipZPV;
    vector < double >_ipZPVerr;
    vector < double >_3dIP;
    vector < double >_3dIPerr;
    vector < double >_3dIPsig;
    vector < double >_closeJetPtAll;
    vector < double >_closeJetEtaAll;
    vector < double >_closeJetPhiAll;
    vector < double >_closeJetEAll;
    vector < double >_closeJetCSVAll;
    vector < int >_closeJetNconstAll;
    vector < double >_closeJetAngAll;
    vector < double >_ptRelAll;
    vector < double >_closeJetPtAllMC;
    vector < double >_closeJetPtAllstatus;
    vector < int >_partonIdMatched;
    vector < bool > _sameParton;
    vector < bool > _isloose;
    vector < bool > _issoftmuon;
    vector < bool > _issoftmuonwithiso;

    vector < bool > _istight;
    vector < bool > _isMediumMuon;
    vector < bool > _istrigemulID;
    vector < bool > _istrigemulISO;
    vector < bool > _istightIso;
    vector < bool > _istightID;
    vector < bool > _istightIDWP2016_noIP3D;
    vector < bool > _istightIDWP2016_IsoEmul;
    vector < bool > _istightIDWP2016;
    vector < bool > _istightIDWP2016_RA7;
    vector < bool > _istightIDWP2016_EWK;
    vector < bool > _isFOIDWP2016_noIP3D;
    vector < bool > _isFOIDWP2016_IsoEmul;
    vector < bool > _isFOIDWP2016;
    vector < bool > _isFOIDWP2016_RA7;
    vector < bool > _isFOIDWP2016_EWK;
    vector < bool > _isVetoIDWP2016_RA7;
    vector < bool > _isVetoIDWP2016_EWK;
    
    vector < bool > _islooseMT2;
    vector < bool > _iscutbasedWPtight;
    vector < double >_mompt;
    vector < double >_momphi;
    vector < double >_mometa;
    vector < int >_mompdg;
    vector <  double >_lPt, _lEta, _lPhi, _lE;
    vector <  double >_lisocorrPt;
    vector < double >_lPtmc, _lEtamc, _lPhimc, _lEmc;
    vector < double >_nuPtmc, _nuEtamc, _nuPhimc, _nuEmc;
    vector < double > _genWPhi, _genZPhi, _genWEta, _genZEta, _genWPt, _genZPt, _genWE, _genZE;
    vector <  double > _lgenPt, _lgenEta,_lgenPhi;
    vector < int> _lgenpdgId;
    vector < double > _mtmc;
    vector < int > _lnmisshits;
    vector < bool > _passconvveto;
    
    vector <double> _tau_dz;
    vector <bool> _TauIs_againstMuonLoose3;
    vector <bool> _TauIs_byLooseCombinedIsolationDeltaBetaCorr3Hits;
    vector <bool> _TauIs_decayModeFindingNewDMs;
    vector <bool> _TauIs_decayModeFinding;
    vector <bool> _TauIs_byLooseIsolationMVArun2v1DBdR03oldDMwLT;
    vector <bool> _TauIs_byMediumIsolationMVArun2v1DBdR03oldDMwLT;
    vector <bool> _TauIs_byVLooseIsolationMVArun2v1DBoldDMwLT;
    vector <bool> _TauIs_byTightIsolationMVArun2v1DBoldDMwLT;
    
    //    triggers https://docs.google.com/spreadsheets/d/1y0E_3zNfBo_tEW5EyDbghwn8_cGa28KOVrQ3DHddh18/edit#gid=1670905758
    vector < int >_lchargeGSF, _lchargeCTF, _lchargePixSc;
    vector < bool > passleg1L;
    vector <  double > _lsietaieta, _lhovere, _lfull5x5sietaieta,_ldphiin,_ldetain,_l1oemin1op; 
    //HLT_TripleMu_12_10_5
    vector <bool>      hltL1sL1TripleMu553;
    vector <bool>      hltL1TripleMu553L2TriMuFiltered3L3TriMuFiltered12105;
    vector <bool>      hltL1TripleMu553L2TriMuFiltered3L3TriMuFiltered10105;
    vector <bool>      hltL1TripleMu553L2TriMuFiltered3L3TriMuFiltered5;

    //HLT_DiMu9_Ele9_CaloIdL_TrackIdL
    vector <bool>      hltDiMu9Ele9CaloIdLTrackIdLMuonlegL3Filtered9;
    vector <bool>      hltDiMu9Ele9CaloIdLTrackIdLElectronlegDphiFilter;

    //HLT_Mu8_DiEle12_CaloIdL_TrackIdL
    vector <bool>      hltMu8DiEle12CaloIdLTrackIdLMuonlegL3Filtered8;
    vector <bool>      hltMu8DiEle12CaloIdLTrackIdLElectronlegDphiFilter;


    //HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL
    vector <bool>      hltL1sL1TripleEG14108;
    vector <bool>      hltEle16Ele12Ele8CaloIdLTrackIdLDphiLeg1Filter;
    vector <bool>      hltEle16Ele12Ele8CaloIdLTrackIdLDphiLeg2Filter;
    vector <bool>      hltEle16Ele12Ele8CaloIdLTrackIdLDphiLeg3Filter;

    //HLT_DoubleMu8_Mass8_PFHT300
    vector <bool>      hltL1sL1Mu6HTT150ORL1Mu8HTT125ORL1HTT175;
    vector <bool>      hltDoubleMu8Mass8L3Filtered;


    //HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT300
    vector <bool>      hltL1sL1DoubleEG6HTT150orL1HTT175;
    vector <bool>      hltDoubleEle8CaloIdMGsfTrackIdMDphiFilter;


    //HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT300
    vector <bool> hltMu8Ele8CaloIdMGsfTrackIdMDphiFilter;
    vector <bool> hltMuon8L3Filtered0 ;

    //HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL

    vector <bool> hltL1sL1DoubleMu103p5ORDoubleMu125;
    vector <bool> hltDiMuonGlb17Glb8RelTrkIsoFiltered0p4;
    vector <bool> hltL3fL1sDoubleMu103p5L1f0L2f10OneMuL3Filtered17;
    vector <bool> hltL3pfL1sDoubleMu103p5L1f0L2pf0L3PreFiltered8;



    //HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL(_DZ)

    vector <bool> hltL3fL1sDoubleMu103p5L1f0L2f10L3Filtered17;
    vector <bool> hltDiMuonGlbFiltered17TrkFiltered8;
    vector <bool> hltDiMuonGlb17Trk8RelTrkIsoFiltered0p4;
    vector <bool> hltDiMuonGlb17Glb8RelTrkIsoFiltered0p4DzFiltered0p2;
    vector <bool> hltDiMuonGlb17Trk8RelTrkIsoFiltered0p4DzFiltered0p2;

    //HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL(_DZ)
    vector <bool> hltL1sL1DoubleEG2210;
    vector <bool> hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter;
    vector <bool> hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter;
    vector <bool> hltEle23Ele12CaloIdLTrackIdLIsoVLDZFilter ; 

    //HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL(_DZ)
    vector <bool> hltL1sL1DoubleEG1510;
    vector <bool> hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter;
    vector <bool> hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter;
    vector <bool> hltEle17Ele12CaloIdLTrackIdLIsoVLDZFilter ; 


    //HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL
    vector <bool> hltL1sL1Mu20EG10; //L1 seed for Mu23 ele12 
    vector <bool> hltL1sSingleMu20erlorSingleMu22lorSingleMu25; //L1 seed for Mu23 ele8
    //    vector <bool> hltMu23TrkIsoVVLEle12CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered23; 
    vector <bool> hltMu23TrkIsoVVLEle8CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered23;
    vector <bool> hltMu23TrkIsoVVLEle8CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter;
    //HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL
    vector <bool> hltL1sL1Mu5EG20; 
    vector <bool> hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered8; 
    vector <bool> hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter; 

    //HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL
    vector <bool> hltL1sL1Mu12EG10; 
    vector <bool> hltMu17TrkIsoVVLEle12CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered17; 
    vector <bool> hltMu17TrkIsoVVLEle12CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter; 

    //HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL
    vector <bool> hltL1sL1Mu5EG15; 
    vector <bool> hltMu8TrkIsoVVLEle17CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered8; 
    vector <bool> hltMu8TrkIsoVVLEle17CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter; 


    //HLT_IsoMu24_eta2p1
    vector <bool> hltL3crIsoL1sMu20Eta2p1L1f0L2f10QL3f24QL3trkIsoFiltered0p09; 

    
    //HLT_IsoMu20

    vector<bool> hltL1sSingleMu18;
    vector <bool> hltL3crIsoL1sMu16L1f0L2f10QL3f20QL3trkIsoFiltered0p09;

    //HLT_Iso(Tk)Mu18
    vector<bool> hltL1sSingleMu16;
    vector <bool> hltL3crIsoL1sMu16L1f0L2f10QL3f18QL3trkIsoFiltered0p09;
    vector <bool> hltL3fL1sMu16L1f0Tkf18QL3trkIsoFiltered0p09;

    //HLT_Iso(Tk)Mu22
    vector<bool> hltL1sSingleMu20;
    vector <bool> hltL3crIsoL1sMu20L1f0L2f10QL3f22QL3trkIsoFiltered0p09; 
    vector <bool> hltL3fL1sMu20L1f0Tkf22QL3trkIsoFiltered0p09;

    //HLT_Ele23_WPLoose_Gsf
    vector <bool> hltL1sSingleEG40IorSingleIsoEG22erIorSingleIsoEG24erIorSingleIsoEG24IorSingleIsoEG26;
    vector <bool> hltEle23WPLooseGsfTrackIsoFilter;
    vector <bool> hltEle23CaloIdLTrackIdLIsoVLTrackIsoFilter;
                  
    //HLT_Ele27_eta2p1_WPLoose_Gsf
    vector <bool> hltEle27WPLooseGsfTrackIsoFilter; 

    //HLT_Ele27_WPLoose_Gsf
    vector <bool> hltEle27noerWPLooseGsfTrackIsoFilter;

    //HLT_Ele32_eta2p1_WPLoose_Gsf
    vector <bool> hltEle32WPLooseGsfTrackIsoFilter; 

    //HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW
    vector <bool> hltDiEle33CaloIdLNewPixelMatchUnseededFilter;

    //HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL
    vector <bool> hltL3fL1sMu22orMu25orMu20EG15orMu5EG20L1f0L2f10QL3Filtered30Q;
    vector <bool>   hltEle30CaloIdLGsfTrkIdVLDPhiUnseededFilter;

    vector <bool>  hltL1sL1SingleMu16ORSingleMu25;
    vector <bool> hltDiMuonGlb27Trk8DzFiltered0p2;
    
    vector <bool> hltL3fL1sMu16orMu25L1f0L2f25L3Filtered27;
    vector <bool> hltL3fL1sMu16orMu25L1f0L2f10QL3Filtered30Q;
  
    //HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW
    vector <bool> hltDiEle33CaloIdLGsfTrkIdVLDPhiUnseededFilter;
    //HLT_DoubleMu3_PFMET50
    vector <bool> hltL3fL1sL1DoubleMu0ETM40lorDoubleMu0ETM55 ;

    //HLT_IsoMu19_eta2p1_LooseIsoPFTau20_SingleL1
    vector <bool> hltL3crIsoL1sSingleMu18erIorSingleMu20erL1f0L2f10QL3f19QL3trkIsoFiltered0p09;
    //HLT_Ele22_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1
    vector <bool> hltEle22WPLooseL1SingleIsoEG20erGsfTrackIsoFilter;
    double _met;
    double _met_phi;
    double _met_JECup;
    double _met_phi_JECup;
    double _met_JECdown;
    double _met_phi_JECdown;


    
    double HT;
    
    double _genmet;
    double _genmet_phi;
    double _mgluino,_mLSP, _mchargino , _mstop;



    long _nEventsTotal;
    long _nEventsTotalCounted;
    long _nEventsFiltered;
    
    TH1D* _hCounter;
    
    double _regVars[15];
    double hJet_ptRaw;
    double hJet_genPt;
    double hJet_pt;
    double hJet_phi;
    double hJet_eta;
    double hJet_e;
    
    double hJet_ptLeadTrack;
    
    double hJet_vtx3dL;
    double hJet_vtx3deL;
    double hJet_vtxMass;
    double hJet_vtxPt;
    
    double hJet_cef;
    double hJet_nconstituents;
    double hJet_JECUnc;
    
    double hJet_SoftLeptptRel;
    double hJet_SoftLeptPt;
    double hJet_SoftLeptdR;
    
    double hJet_SoftLeptIdlooseMu;
    double hJet_SoftLeptId95;
};

#endif
