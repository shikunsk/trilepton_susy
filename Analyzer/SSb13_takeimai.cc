#include "SSb13_takeimai.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/View.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "CommonTools/CandAlgos/interface/CandMatcher.h"
#include "PhysicsTools/HepMCCandAlgos/interface/MCTruthPairSelector.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Common/interface/MergeableCounter.h"

#include "PhysicsTools/CandUtils/interface/CandMatcherNew.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include <DataFormats/L1Trigger/interface/L1EtMissParticle.h>

using namespace std;
using namespace edm;
using namespace reco;
using namespace tools;
using namespace math;
using namespace reco::tau;
const int statusgenparticles = 23;
const double btag_cut = 0.800;
SSb13_takeimai::SSb13_takeimai(const edm::ParameterSet & iConfig) :
_relIsoCutE(0.10),
_relIsoCutMu(0.10),
//_relIsoCutEloose(0.5), //0.6
//_relIsoCutMuloose(0.5), //1.0
_relIsoCutEloose(999999.), //0.6
_relIsoCutMuloose(999999.), //1.0
_chargeConsistency(false),
_useconversions(false),
_tightD0Mu(0.02),
_tightD0E(0.05),
_looseD0Mu(10.05),
_looseD0E(10.05),
//_looseD0Mu(0.2),
//_looseD0E(9999999.),
_jetPtCut(20.), //20
_jetEtaCut(2.4),
_eleMinPt(5.),//Ele Loose Pt 5
_muonMinPt(5.), //Mu Loose Pt 5
//_tauMinPt(1800000.), //15
_tauMinPt(20.), //20
_tauPt(1800000.),
_tauEta(2.3),
_regression(false),
mvaValuesMapToken_(consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("mvaValuesMap"))),
mvaCategoriesMapToken_(consumes<edm::ValueMap<int> >(iConfig.getParameter<edm::InputTag>("mvaCategoriesMap")))
{
    Sample              = iConfig.getUntrackedParameter<std::string>("SampleLabel") ;
    Skim                = iConfig.getUntrackedParameter<std::string>("theskim");
    IT_muon             = iConfig.getParameter<edm::InputTag>("MuonLabel") ;
    IT_electron         = iConfig.getParameter<edm::InputTag>("ElectronLabel") ;
    IT_tau              = iConfig.getParameter<edm::InputTag>("TauLabel") ;
    //IT_tauDiscriminator = iConfig.getParameter<edm::InputTag>("TauDiscriminatorLabel") ;
    IT_jet              = iConfig.getParameter<edm::InputTag>("JetLabel");
    IT_pfmet            = iConfig.getParameter<edm::InputTag>("METLabel")  ;
    IT_pileup            = iConfig.getParameter<edm::InputTag>("PULabel")  ;
    IT_beamspot         = iConfig.getParameter<edm::InputTag>("BeamSpotLabel");
    IT_hltresults       = iConfig.getParameter<edm::InputTag>("HLTResultsLabel");
    IsMC                = iConfig.getUntrackedParameter<bool>("isMCLabel");
    IsFS                = iConfig.getUntrackedParameter<bool>("isFSLabel");
    IT_METFilters       = iConfig.getParameter<edm::InputTag>("METFilter");

    puInfo_token = consumes<std::vector<PileupSummaryInfo> >(iConfig.getParameter<edm::InputTag>("PULabel"));     

    metfilters_token = consumes<TriggerResults>(IT_METFilters);
    filterbadChCand_token  = consumes<bool>(edm::InputTag("BadChargedCandidateFilter"));
    filterbadPFMuon_token = consumes<bool>(edm::InputTag("BadPFMuonFilter"));
    //    filterbadChCand_token  = consumes<bool>(edm::InputTag("BadChCandFilter"));
    //filterbadPFMuon_token = consumes<bool>(edm::InputTag("BadPFMuon"));
    
    IT_gen = edm::InputTag("generator");
    geninfo_token = consumes<GenEventInfoProduct>(IT_gen); 

    IT_genpart= edm::InputTag("prunedGenParticles");
    genpart_token  = consumes<GenParticleCollection>(IT_genpart);
    
    trgresults_token= consumes<TriggerResults>(IT_hltresults);
    
    beamspot_token = consumes<reco::BeamSpot>( IT_beamspot );
    
    IT_goodVtx = edm::InputTag("offlineSlimmedPrimaryVertices");
    goodpv_token = consumes<std::vector<Vertex> > (IT_goodVtx);

    met_token = consumes<std::vector<pat::MET> > (IT_pfmet);
    
    IT_pfcands = edm::InputTag("packedPFCandidates");
    pfcands_token  = consumes<pat::PackedCandidateCollection>(IT_pfcands);

    patmuon_token = consumes<std::vector<pat::Muon> >(IT_muon);

    patelectron_token = consumes<std::vector<pat::Electron> >(IT_electron);

    patelectronview_token= consumes< edm::View<pat::Electron> >(IT_electron);
    
    pattau_token =consumes<std::vector<pat::Tau> >(IT_tau);
    
    IT_conv = edm::InputTag("reducedEgamma","reducedConversions");
    conv_token = consumes< std::vector<reco::Conversion> >(IT_conv);
    
    jet_token = consumes< std::vector< pat::Jet> >(IT_jet);

    rhoJets_token = consumes<double>(edm::InputTag("fixedGridRhoFastjetAll",""));

    rhoJetsNC_token = consumes<double>(edm::InputTag("fixedGridRhoFastjetCentralNeutral",""));

    emIsolColl_token = consumes< l1extra::L1EmParticleCollection > (edm::InputTag("l1extraParticles","Isolated"));

    IT_triggerObjects= edm::InputTag("selectedPatTrigger");
    trigobject_token = consumes<pat::TriggerObjectStandAloneCollection>(IT_triggerObjects);
    
    IT_triggerPrescales= edm::InputTag("patTrigger");
    IT_L1Prescales= edm::InputTag("patTrigger","l1min");
    trgpresc_token = consumes<pat::PackedTriggerPrescales>(IT_triggerPrescales);
    L1presc_token = consumes<pat::PackedTriggerPrescales>(IT_L1Prescales);
    
isgetbylabel = false;
    //outfile = fopen("FakeSync.txt", "w");

 usemva = true;
 if(usemva){
 tmvareaderele = new TMVA::Reader( "!Color:!Silent" );
 tmvareaderele->AddVariable("LepGood_pt",&LepGood_pt_forMVA);
 tmvareaderele->AddVariable("LepGood_eta",&LepGood_eta_forMVA);
 tmvareaderele->AddVariable("LepGood_jetNDauChargedMVASel",&LepGood_jetNDauChargedMVASel_forMVA);
 tmvareaderele->AddVariable("LepGood_miniRelIsoCharged",&LepGood_miniRelIsoCharged_forMVA);
 tmvareaderele->AddVariable("LepGood_miniRelIsoNeutral",&LepGood_miniRelIsoNeutral_forMVA);
 tmvareaderele->AddVariable("LepGood_jetPtRelv2",&LepGood_jetPtRelv2_forMVA);
 tmvareaderele->AddVariable("min(LepGood_jetPtRatiov2,1.5)",&LepGood_jetPtRatio_forMVA);
 tmvareaderele->AddVariable("max(LepGood_jetBTagCSV,0)",&LepGood_jetBTagCSV_forMVA);
 tmvareaderele->AddVariable("LepGood_sip3d",&LepGood_sip3d_forMVA);
 tmvareaderele->AddVariable("log(abs(LepGood_dxy))",&LepGood_dxy_forMVA);
 tmvareaderele->AddVariable("log(abs(LepGood_dz))",&LepGood_dz_forMVA);
 tmvareaderele->AddVariable("LepGood_mvaIdSpring16GP",&LepGood_mvaIdSpring16GP_forMVA);
 //tmvareaderele->BookMVA("BDTG method","/cms/data/store/user/t2/users/takeimai/CMSSW_8_0_13/src/SUSYAnalyzer/PatAnalyzer/test/el_BDTG.weights.xml");
 tmvareaderele->BookMVA("BDTG method","forMoriond16_el_sigTTZ_bkgTT_BDTG.weights.xml");


 tmvareadermu = new TMVA::Reader( "!Color:!Silent");
 tmvareadermu->AddVariable("LepGood_pt",&LepGood_pt_forMVA);
 tmvareadermu->AddVariable("LepGood_eta",&LepGood_eta_forMVA);
 tmvareadermu->AddVariable("LepGood_jetNDauChargedMVASel",&LepGood_jetNDauChargedMVASel_forMVA);
 tmvareadermu->AddVariable("LepGood_miniRelIsoCharged",&LepGood_miniRelIsoCharged_forMVA);
 tmvareadermu->AddVariable("LepGood_miniRelIsoNeutral",&LepGood_miniRelIsoNeutral_forMVA);
 tmvareadermu->AddVariable("LepGood_jetPtRelv2",&LepGood_jetPtRelv2_forMVA);
 tmvareadermu->AddVariable("min(LepGood_jetPtRatiov2,1.5)",&LepGood_jetPtRatio_forMVA);
 tmvareadermu->AddVariable("max(LepGood_jetBTagCSV,0)",&LepGood_jetBTagCSV_forMVA);
 tmvareadermu->AddVariable("LepGood_sip3d",&LepGood_sip3d_forMVA);
 tmvareadermu->AddVariable("log(abs(LepGood_dxy))",&LepGood_dxy_forMVA);
 tmvareadermu->AddVariable("log(abs(LepGood_dz))",&LepGood_dz_forMVA);
 tmvareadermu->AddVariable("LepGood_segmentCompatibility",&LepGood_segmentCompatibility_forMVA);
 //tmvareadermu->BookMVA("BDTG method","/cms/data/store/user/t2/users/takeimai/CMSSW_8_0_13/src/SUSYAnalyzer/PatAnalyzer/test/mu_BDTG.weights.xml");
 tmvareadermu->BookMVA("BDTG method","forMoriond16_mu_sigTTZ_bkgTT_BDTG.weights.xml");
 }
}


void SSb13_takeimai::beginJob()
{
  cout << "deb0 "<< endl;
    Nvtx           = fs->make<TH1F>("N_{vtx}"        , "Number of vertices;N_{vtx};events / 1"  ,    50, 0., 50.);
    cout << "deb0.1 "<< endl;
    _hCounter = fs->make<TH1D>("hCounter", "Events counter", 5,0,5);
    cout << "deb0.2"<< endl;
    //   outputTree = new TTree("fakeTree","fakeTree");
    outputTree = fs->make<TTree>("fakeTree","fakeTree");
    cout << "deb0.3 "<< endl;
    bookTree();
    cout << "deb1 "<< endl; 
    _leptonP4 = new TClonesArray("TLorentzVector", nLeptonsMax);
    for (int i=0; i!=nLeptonsMax; ++i) {
        new ( (*_leptonP4)[i] ) TLorentzVector();
    }
    
    _jetP4 = new TClonesArray("TLorentzVector", nJetsMax);
    for (int i=0; i!=nJetsMax; ++i) {
        new ( (*_jetP4)[i] ) TLorentzVector();
    }
        cout << "deb2 "<< endl; 
    GPM = GenParticleManager();
    cout << "deb3 "<< endl; 

    //    bool isData = !(Sample=="ElectronsMC");
    if (!IsMC)
      fMetCorrector = new OnTheFlyCorrections("Spring16_25nsV6_DATA", !IsMC);  
      // fMetCorrector = new OnTheFlyCorrections("Summer15_25nsV7_DATA", !IsMC);
    else
      fMetCorrector = new OnTheFlyCorrections("Spring16_25nsV6_MC", !IsMC); 
    _corrLevel = "L3Absolute";
    if (!IsMC) _corrLevel = "L2L3Residual";

    
    //calib = new BTagCalibration("csvv2","/cms/data/store/user/t2/users/takeimai/CMSSW_8_0_13/src/SUSYAnalyzer/PatAnalyzer/btagSFandEff/CSVv2_ichep.csv");
    calib = new BTagCalibration("csvv2","btagSFandEff/CSVv2_ichep.csv");
    //BTagCalibration calibbis ("csvv2","/cms/data/store/user/t2/users/takeimai/CMSSW_8_0_13/src/SUSYAnalyzer/PatAnalyzer/btagSFandEff/CSVv2_ichep.csv");
    BTagCalibration calibbis ("csvv2","btagSFandEff/CSVv2_ichep.csv");


    reader = new BTagCalibrationReader(BTagEntry::OP_MEDIUM,"central");
    reader->load(calibbis,                // calibration instance
		BTagEntry::FLAV_B,    // btag flavour
		"mujets") ;              // measurement type

    reader_up = new BTagCalibrationReader(BTagEntry::OP_MEDIUM,"up");
    reader_up->load(calibbis,                // calibration instance
		BTagEntry::FLAV_B,    // btag flavour
		"mujets") ;              // measurement type

    reader_do = new BTagCalibrationReader(BTagEntry::OP_MEDIUM,"down");
    reader_do->load(calibbis,                // calibration instance
		BTagEntry::FLAV_B,    // btag flavour
		"mujets") ;              // measurement type


    
    //TFile *fileforbtageff = new TFile("/cms/data/store/user/t2/users/takeimai/CMSSW_8_0_13/src/SUSYAnalyzer/PatAnalyzer/btagSFandEff/bTagEffs.root","read");
    TFile *fileforbtageff = new TFile("btagSFandEff/bTagEffs.root","read");
    //    const char *fl[3] = {"b","c","udsg"};

    TString flstring[3] = {"b","c","udsg"};
    double ptbinsB[18]={20,25,30,35,40,45,50,60,70,80,90,100,120,150,200,300,400,600};
    double etabinsB[18]={0,0.4,0.8,1.2,1.6,2.0,2.4,2.8};
    for (int i=0; i!=3;++i) {

      btagEffs2D[i] = new TH2F("eff_total_M_"+flstring[i],"",17,ptbinsB,7,etabinsB);
      btagEffs2D[i]->Read("eff_total_M_"+flstring[i]);
      std::cout<<"Read "<<"eff_total_M_"+flstring[i]<<std::endl;

    }

    fileforbtageff->Close();


    cout << "deb5 "<< endl;

    //Spring 2016 NEW WP
    looseMVA[0][0] = -0.70;
    looseMVA[1][0] = -0.83;
    looseMVA[2][0] = -0.92;

    looseMVA[0][1] = 0.87;
    looseMVA[1][1] = 0.60;
    looseMVA[2][1] = 0.17;
    
    /*
    //Loose//Pt<10
    looseMVA[0][0][0] = 0.46;
    looseMVA[1][0][0] = -0.03;
    looseMVA[2][0][0] = 0.06;
    //Loose//10<Pt<20
    looseMVA[0][0][1] = -0.86;
    looseMVA[1][0][1] = -0.85;
    looseMVA[2][0][1] = -0.81;
    //Loose//Pt>20
    looseMVA[0][0][2] = -0.96;
    looseMVA[1][0][2] = -0.96;
    looseMVA[2][0][2] = -0.95;
    
    //Tight//Pt<10
    looseMVA[0][1][0] = 0.77;
    looseMVA[1][1][0] = 0.56;
    looseMVA[2][1][0] = 0.48;
    //Tight//10<Pt<20
    looseMVA[0][1][1] = 0.77;
    looseMVA[1][1][1] = 0.56;
    looseMVA[2][1][1] = 0.48;
    //Tight//Pt>20
    looseMVA[0][1][2] = 0.52;
    looseMVA[1][1][2] = 0.11;
    looseMVA[2][1][2] = -0.01;
    */

    _nEventsTotal = 0;
    _nEventsFiltered = 0;
    _nEventsTotalCounted = 0;
    
    /*    _miniisocut[0] = 0.1;
    _ptratiocut[0] = 0.7;
    _ptrelcut[0] = 7;*/

    //Electron values for multi iso
    //New wp?  https://indico.cern.ch/event/450052/contribution/0/attachments/1162935/1675309/lepawareJECv2_bkg_wp_300915.pdf
    _miniisocut[0] = 0.12;//0.13
    _ptratiocut[0] = 0.80;//0.81
    _ptrelcut[0] = 7.2;//7.2
    
    _miniisocut[1] = 0.16;//0.21
    _ptratiocut[1] = 0.76;//0.8
    _ptrelcut[1] = 7.2;//6.9
    
    _multiConst[0][0] = 0.22;
    _multiConst[0][1] = 0.63;
    _multiConst[0][2] = 6;
    
    _multiConst[1][0] = 0.14;
    _multiConst[1][1] = 0.68;
    _multiConst[1][2] = 6.7;
    
    _multiConst[2][0] = 0.10;
    _multiConst[2][1] = 0.70;
    _multiConst[2][2] = 7;
    
}

void SSb13_takeimai::endJob() {
    //outputTree -> Write();
    // store nEventsTotal and nEventsFiltered in preferred way
    std::cout<<_nEventsTotal<<std::endl;
    std::cout<<_nEventsFiltered<<std::endl;
    
    delete fMetCorrector;
}

void SSb13_takeimai::analyze(const edm::Event& iEvent, const edm::EventSetup& iEventSetup)
{




  //bool islepton;

  double jetptcut =20;
  bool debug =false;
  if(debug) cout << "beginning of analyze method"<<endl;
  ClearStuff();
  
  const int genpartsize = GenParticleInfo(iEvent);
  vector<bool> genpartalreadymatched;
  for( int i = 0; i < genpartsize;i++){
    genpartalreadymatched.push_back(false);
  }

  //============ Total number of events is the sum of the events ============
  //============ in each of these luminosity blocks ============

  _weight =1;

  if(IsMC){ 
    edm::Handle<GenEventInfoProduct> pdfvariables;
    if(isgetbylabel)    iEvent.getByLabel("generator", pdfvariables);
    else iEvent.getByToken(geninfo_token, pdfvariables); 
    _weight=pdfvariables->weight();

  }
  //  if (_weight > 0) _weight = 1;
  // else _weight = -1;
  _nEventsTotalCounted+=_weight;
  _hCounter->Fill(0.,_weight);
  //============ Counter done ============
  
  //Now applying skim :-) 
  
 
  if (!ApplySkim(iEvent, Skim)) return;
  //  if(!IsMC) {if(!PassMETFilters(iEvent,metfilters_token,filterbadChCand_token,filterbadPFMuon_token) && !IsMC) return;}


 
 



  if(debug) cout << "Debug 0 " << endl; 
  if (IsMC) {
        _genHT = 0;
        //******************************************************************************************************************
        // Gen level particles                  ****************************************************************************
        //******************************************************************************************************************
	if(isgetbylabel) iEvent.getByLabel("prunedGenParticles", TheGenParticles);
	else iEvent.getByToken(genpart_token, TheGenParticles); 
        std::vector<const GenParticle*> vGenElectrons, vGenMuons, vGenNPElectrons, vGenNPMuons, vGenW;
        if( TheGenParticles.isValid() ) {
            GPM.SetCollection(TheGenParticles);
            GPM.Classify();
            vGenMuons = GPM.filterByStatus(GPM.getPromptMuons(),statusgenparticles);
            vGenElectrons = GPM.filterByStatus(GPM.getPromptElectrons(),statusgenparticles);
            vGenNPMuons = GPM.filterByStatus(GPM.getNonPromptMuons(),statusgenparticles);
            vGenNPElectrons = GPM.filterByStatus(GPM.getNonPromptElectrons(),statusgenparticles);
	    
	    
	    //std::cout<<"** vGenElectronssize "<<vGenElectrons.size()<<std::endl;
	    // std::cout<<"** vGenMuonssize "<<vGenMuons.size()<<std::endl;
	    //std::cout<<"** vGenNPElectronssize "<<vGenNPElectrons.size()<<std::endl;
            //std::cout<<"** vGenNPMuonssize "<<vGenNPMuons.size()<<std::endl;

	    //if( vGenMuons.size() + vGenElectrons.size()<2 ) cout << "failed lepton size" << endl;
	    //            if( vGenMuons.size() + vGenElectrons.size()<2 &&Twolatgenlevel ) return;
	    if(debug) cout << "Debug 1 " <<endl;
            TLorentzVector Gen0;
            Gen0.SetPtEtaPhiE( 0, 0, 0, 0);
            for(GenParticleCollection::const_reverse_iterator p = TheGenParticles->rbegin() ; p != TheGenParticles->rend() ; p++ ) {
	      
                int id = TMath::Abs(p->pdgId());
                
                //Find the W boson and Z boson at the gen level
                if (id==24) {//W boson
                    _genWPhi.push_back(p->phi());
                    _genWEta.push_back(p->eta());
                    _genWPt.push_back(p->pt());
                    _genWE.push_back(p->energy());
                }
                if (id==23) {//Z boson
                    _genZPhi.push_back(p->phi());
                    _genZEta.push_back(p->eta());
                    _genZPt.push_back(p->pt());
                    _genZE.push_back(p->energy());
                }
                
		/*		if(fabs(id) ==11||fabs(id) ==13 || fabs(id) ==15 ){
		  cout << "found lepton: id, status, pt, eta phi " <<p->pdgId()<<","<< p->status()<<","<<p->pt()<< ","<<p->eta()<<","<< p->phi() <<endl;
		  }*/
                if ( (id == 12 || id == 14 || id == 16 ) && (p->status() == 1) ) {
                    TLorentzVector Gen;
                    Gen.SetPtEtaPhiE( p->pt(), p->eta(), p->phi(), p->energy() );
                    Gen0 += Gen;
                }
                if ((id == 1 || id == 2 || id == 3 || id == 4 || id == 5 || id == 21 || id == 22 ) && (p->status() == 23)){
                    _genHT += p->pt();
                }
            }
            if (Gen0.E()!=0) {
                _genmet = Gen0.Pt();
                _genmet_phi = Gen0.Phi();
            } else {
                _genmet = 0;
                _genmet_phi = 0;
            }
        }
        //******************************************************************************************************************
        //******************************************************************************************************************
        //******************************************************************************************************************
        
        //**************************************************************************************
        // MC
        //**************************************************************************************
    }
    
    _runNb = iEvent.id().run();
    _eventNb = iEvent.id().event();
    _lumiBlock = iEvent.luminosityBlock();
    //    cout << "runnb, ls, eventnb " << _runNb << ", " <<_lumiBlock<<", "<< _eventNb<<endl;
    //    if ( _eventNb!=204673540) return;
    if(debug) cout << "Debug 2 " <<endl;
    nbofL1EG =0;

    passmetfilters =true;
    passcsctighthalo2015=true; 
    passglobaltighthalo2016=true; 
    passglobalsupertighthalo2016=true;

    passHLT_PFMET170_NotCleaned=IsFS;
    passHLT_PFMET170_HBHECleaned=IsFS;
    passHLT_PFMET170_BeamHaloCleaned=IsFS;
    passHLT_PFMET120_Mu5=IsFS;
    passHLT_Ele15_IsoVVVL_PFHT350_PFMET50=IsFS;
    passHLT_IsoMu16_eta2p1_MET30=IsFS;
    passHLT_IsoMu19_eta2p1_LooseIsoPFTau20_v=IsFS;
    passHLT_Ele22_eta2p1_WPLoose_GSF_LooseIsoPFtau20_SingleL1_v=IsFS;
    passHLT_Mu17_Mu8_SameSign_DZ=IsFS;
    passHLT_DoubleMu3_PFMET50=IsFS;
    passHLT_TripleMu_5_3_3=IsFS;
    passHLT_Ele8_CaloIdM_TrackIdM_PFJet30=IsFS;
    passHLT_Ele12_CaloIdM_TrackIdM_PFJet30=IsFS;
    passHLT_Ele17_CaloIdM_TrackIdM_PFJet30=IsFS;
    passHLT_Ele23_CaloIdM_TrackIdM_PFJet30=IsFS;
    passHLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30=IsFS;
    passHLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30=IsFS;
    passHLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30=IsFS;
    passHLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30=IsFS;
    passHLT_Ele15_IsoVVVL_PFHT350=IsFS;
    passHLT_Ele10_CaloIdM_TrackIdM_CentralPFJet30_BTagCSV_p13=IsFS;
    passHLT_Mu10_CentralPFJet30_BTagCSV_p13=IsFS;
    passHLT_Mu8=IsFS;
    passHLT_Mu17=IsFS;
    passHLT_Mu20=IsFS;
    passHLT_Mu27=IsFS;
    passHLT_Mu8_TrkIsoVVL=IsFS;
    passHLT_Mu17_TrkIsoVVL=IsFS;
    passHLT_Mu15_IsoVVVL_PFHT350=IsFS;

    



    passHLT_DoubleMu8_Mass8_PFHT300 = IsFS;
    passHLT_DoubleMu8_Mass8_PFHT250 = IsFS;

    passHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL= IsFS;
    passHLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL= IsFS;
    passHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ= IsFS;
    passHLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ= IsFS;

    passHLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT300= IsFS;
    passHLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT250= IsFS;
    passHLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL= IsFS;
    passHLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ= IsFS;
    passHLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL= IsFS;
    passHLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ= IsFS;

    passHLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT300= IsFS;
    passHLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT250= IsFS;
    passHLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL= IsFS;
    passHLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL= IsFS;
    passHLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL= IsFS;
    passHLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL= IsFS;
    passHLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL= IsFS;
    passHLT_IsoMu24_eta2p1= IsFS;
    passHLT_IsoMu20= IsFS;
    passHLT_IsoTkMu20= IsFS;
    passHLT_IsoMu18= IsFS;
    passHLT_IsoTkMu18= IsFS;
    passHLT_IsoMu22= IsFS;
    passHLT_IsoTkMu22= IsFS;


    passHLT_Ele23_WPLoose_Gsf= IsFS;
    passHLT_Ele27_eta2p1_WPLoose_Gsf= IsFS;
    passHLT_Ele27_WPLoose_Gsf= IsFS;
    passHLT_Ele32_eta2p1_WPLoose_Gsf= IsFS;  
    passHLT_TripleMu_12_10_5= IsFS;
    passHLT_DiMu9_Ele9_CaloIdL_TrackIdL= IsFS;
    passHLT_Mu8_DiEle12_CaloIdL_TrackIdL= IsFS;
    passHLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL= IsFS;
    
    passHLT_Ele23_CaloIdL_TrackIdL_IsoVL =IsFS;
    passHLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW= IsFS;
    
    //Prescales
    Prescale_passHLT_Mu8=-99;
    Prescale_passHLT_Mu17=-99;
    Prescale_passHLT_Mu8_TrkIsoVVL=-99;
    Prescale_passHLT_Mu17_TrkIsoVVL=-99;
    Prescale_passHLT_Ele12_CaloIdM_TrackIdM_PFJet30=-99;
    Prescale_passHLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30=-99;
    Prescale_passHLT_Ele17_CaloIdM_TrackIdM_PFJet30=-99;
    Prescale_passHLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30=-99;



    edm::Handle<TriggerResults> METFilterResults;

    iEvent.getByToken(metfilters_token, METFilterResults);
    if(!PassMETFilters(iEvent,metfilters_token,filterbadChCand_token,filterbadPFMuon_token)) passmetfilters=false; 
    if( !METFilterResults.failedToGet() ) {
        int N_MetFilters = METFilterResults->size();
        const edm::TriggerNames & metfilterName = iEvent.triggerNames(*METFilterResults);
        
        for( int i_Metfilter = 0; i_Metfilter < N_MetFilters; ++i_Metfilter ) {
            //      cout <<  trigName.triggerName(i_Trig) << endl;
            if (!METFilterResults.product()->accept(i_Metfilter)) {
                TString MetfilterPath =metfilterName.triggerName(i_Metfilter);
                if(MetfilterPath.Index("Flag_CSCTightHalo2015Filter")>=0) passcsctighthalo2015 =false;
                if(MetfilterPath.Index("Flag_globalTightHalo2016Filter")>=0) passglobaltighthalo2016 =false;
                if(MetfilterPath.Index("Flag_globalSuperTightHalo2016Filter")>=0) passglobalsupertighthalo2016 =false;
            }
        }
    }



    edm::Handle<TriggerResults> trigResults;
    if(isgetbylabel)    iEvent.getByLabel(IT_hltresults, trigResults);
    else iEvent.getByToken(trgresults_token, trigResults);

    edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
    edm::Handle<pat::PackedTriggerPrescales> L1Prescales;
    //edm::InputTag IT_triggerPrescales("patTrigger");
    if(isgetbylabel)  iEvent.getByLabel(IT_triggerPrescales, triggerPrescales); 
    else iEvent.getByToken(trgpresc_token, triggerPrescales);
    if(isgetbylabel)  iEvent.getByLabel(IT_L1Prescales, L1Prescales);
    else iEvent.getByToken(L1presc_token, L1Prescales);


    if(debug) cout << "Debug 3 " <<endl;
    //    if( trigResults.failedToGet() ) cout << "--- NO TRIGGER RESULTS !! ---" << endl;
    if( !trigResults.failedToGet() ) {
        
     	int N_Triggers = trigResults->size();
        
     	const edm::TriggerNames & trigName = iEvent.triggerNames(*trigResults);


     	for( int i_Trig = 0; i_Trig < N_Triggers; ++i_Trig ) {
	           
            if (trigResults.product()->accept(i_Trig)) {
                
                TString TrigPath =trigName.triggerName(i_Trig);
                if(TrigPath.Contains("HLT_PFMET170_NotCleaned"))passHLT_PFMET170_NotCleaned =true;
                if(TrigPath.Contains("HLT_PFMET170_HBHECleaned"))passHLT_PFMET170_HBHECleaned =true;
                if(TrigPath.Contains("HLT_PFMET170_BeamHaloCleaned"))passHLT_PFMET170_BeamHaloCleaned =true;
                if(TrigPath.Contains("HLT_PFMET120_Mu5"))passHLT_PFMET120_Mu5 =true;
                if(TrigPath.Contains("HLT_Ele15_IsoVVVL_PFHT350_PFMET50"))passHLT_Ele15_IsoVVVL_PFHT350_PFMET50 =true;
                if(TrigPath.Contains("HLT_IsoMu16_eta2p1_MET30"))passHLT_IsoMu16_eta2p1_MET30 =true;
                if(TrigPath.Contains("HLT_IsoMu19_eta2p1_LooseIsoPFTau20_v"))passHLT_IsoMu19_eta2p1_LooseIsoPFTau20_v =true;
                if(TrigPath.Contains("HLT_Ele22_eta2p1_WPLoose_GSF_LooseIsoPFtau20_SingleL1_v"))passHLT_Ele22_eta2p1_WPLoose_GSF_LooseIsoPFtau20_SingleL1_v =true;
                if(TrigPath.Contains("HLT_Mu17_Mu8_SameSign_DZ") || TrigPath.Contains("HLT_Mu20_Mu10_SameSign_DZ") )passHLT_Mu17_Mu8_SameSign_DZ =true;
                if(TrigPath.Contains("HLT_DoubleMu3_PFMET50"))passHLT_DoubleMu3_PFMET50 =true;
                if(TrigPath.Contains("HLT_TripleMu_5_3_3"))passHLT_TripleMu_5_3_3 =true;
                if(TrigPath.Contains("HLT_Ele8_CaloIdM_TrackIdM_PFJet30"))passHLT_Ele8_CaloIdM_TrackIdM_PFJet30 =true;
                if(TrigPath.Contains("HLT_Ele12_CaloIdM_TrackIdM_PFJet30_v"))passHLT_Ele12_CaloIdM_TrackIdM_PFJet30 =true;
                if(TrigPath.Contains("HLT_Ele12_CaloIdM_TrackIdM_PFJet30_v")&&!IsMC){
                    Prescale_passHLT_Ele12_CaloIdM_TrackIdM_PFJet30=(triggerPrescales->getPrescaleForIndex(i_Trig))*(L1Prescales->getPrescaleForIndex(i_Trig));
                }
                if(TrigPath.Contains("HLT_Ele17_CaloIdM_TrackIdM_PFJet30_v"))passHLT_Ele17_CaloIdM_TrackIdM_PFJet30 =true;
                if(TrigPath.Contains("HLT_Ele17_CaloIdM_TrackIdM_PFJet30_v")&&!IsMC){
                    Prescale_passHLT_Ele17_CaloIdM_TrackIdM_PFJet30=(triggerPrescales->getPrescaleForIndex(i_Trig))*(L1Prescales->getPrescaleForIndex(i_Trig));
                }
                if(TrigPath.Contains("HLT_Ele23_CaloIdM_TrackIdM_PFJet30"))passHLT_Ele23_CaloIdM_TrackIdM_PFJet30 =true;
                if(TrigPath.Contains("HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30"))passHLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30 =true;
                if(TrigPath.Contains("HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v"))passHLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30 =true;
                if(TrigPath.Contains("HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v")&&!IsMC){
                    Prescale_passHLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30=(triggerPrescales->getPrescaleForIndex(i_Trig))*(L1Prescales->getPrescaleForIndex(i_Trig));
                }
                if(TrigPath.Contains("HLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30_v"))passHLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30 =true;
                if(TrigPath.Contains("HLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30_v")&&!IsMC){
                    Prescale_passHLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30=(triggerPrescales->getPrescaleForIndex(i_Trig))*(L1Prescales->getPrescaleForIndex(i_Trig));
                }
                if(TrigPath.Contains("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30"))passHLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30 =true;
                if(TrigPath.Contains("HLT_Ele10_CaloIdM_TrackIdM_CentralPFJet30_BTagCSV_p13"))passHLT_Ele10_CaloIdM_TrackIdM_CentralPFJet30_BTagCSV_p13 =true;
                if(TrigPath.Contains("HLT_Mu10_CentralPFJet30_BTagCSV_p13"))passHLT_Mu10_CentralPFJet30_BTagCSV_p13 =true;
                if(TrigPath.Contains("HLT_Mu8_v")) passHLT_Mu8 =true;
                if(TrigPath.Contains("HLT_Mu8_v")&&!IsMC){
                    Prescale_passHLT_Mu8=(triggerPrescales->getPrescaleForIndex(i_Trig))*(L1Prescales->getPrescaleForIndex(i_Trig));
                }
                if(TrigPath.Contains("HLT_Mu17_v"))passHLT_Mu17 =true;
                if(TrigPath.Contains("HLT_Mu17_v")&&!IsMC){
                    Prescale_passHLT_Mu17=(triggerPrescales->getPrescaleForIndex(i_Trig))*(L1Prescales->getPrescaleForIndex(i_Trig));
                }
                if(TrigPath.Contains("HLT_Mu20"))passHLT_Mu20 =true;
                if(TrigPath.Contains("HLT_Mu27"))passHLT_Mu27 =true;
                if(TrigPath.Contains("HLT_Mu8_TrkIsoVVL_v"))passHLT_Mu8_TrkIsoVVL =true;
                if(TrigPath.Contains("HLT_Mu8_TrkIsoVVL_v")&&!IsMC){
                    Prescale_passHLT_Mu8_TrkIsoVVL=(triggerPrescales->getPrescaleForIndex(i_Trig))*(L1Prescales->getPrescaleForIndex(i_Trig));
                }
                if(TrigPath.Contains("HLT_Mu17_TrkIsoVVL_v"))passHLT_Mu17_TrkIsoVVL =true;
                if(TrigPath.Contains("HLT_Mu17_TrkIsoVVL_v")&&!IsMC){
                    Prescale_passHLT_Mu17_TrkIsoVVL=(triggerPrescales->getPrescaleForIndex(i_Trig))*(L1Prescales->getPrescaleForIndex(i_Trig));
                }
                
                
                
                
                if(TrigPath.Contains("HLT_Ele15_IsoVVVL_PFHT350"))passHLT_Ele15_IsoVVVL_PFHT350=true;
                if(TrigPath.Contains("HLT_Mu15_IsoVVVL_PFHT350"))passHLT_Mu15_IsoVVVL_PFHT350=true;
                
                
                if(TrigPath.Contains("HLT_DoubleMu8_Mass8_PFHT300")) passHLT_DoubleMu8_Mass8_PFHT300=true;
                if(TrigPath.Contains("HLT_DoubleMu8_Mass8_PFHT250")) passHLT_DoubleMu8_Mass8_PFHT250=true;
                if(TrigPath.Contains("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL")) passHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL=true;
                if(TrigPath.Contains("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL")) passHLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL=true;
                if(TrigPath.Contains("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ")) passHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ=true;
                if(TrigPath.Contains("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ")) passHLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ=true;
                if(TrigPath.Contains("HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT300")) passHLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT300=true;
                if(TrigPath.Contains("HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT250")) passHLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT250=true;
                if(TrigPath.Contains("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL")) passHLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL=true;
                if(TrigPath.Contains("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ")) passHLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ=true;
                if(TrigPath.Contains("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL")) passHLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL=true;
                if(TrigPath.Contains("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ")) passHLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ=true;
                
                if(TrigPath.Contains("HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT300"))passHLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT300 =true;
                if(TrigPath.Contains("HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT250"))passHLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT250 =true;
                if(TrigPath.Contains("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL"))passHLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL=true;
                if(TrigPath.Contains("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL"))passHLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL=true;
                if(TrigPath.Contains("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL")) passHLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL=true;
                if(TrigPath.Contains("HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL"))passHLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL=true;
                if(TrigPath.Contains("HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL")) passHLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL=true;
                
                if(TrigPath.Contains("HLT_IsoMu24_eta2p1")) passHLT_IsoMu24_eta2p1=true;
                if(TrigPath.Contains("HLT_IsoMu20")) passHLT_IsoMu20=true;
                if(TrigPath.Contains("HLT_IsoTkMu20")) passHLT_IsoTkMu20=true;
                if(TrigPath.Contains("HLT_IsoMu18")) passHLT_IsoMu18=true;
                if(TrigPath.Contains("HLT_IsoTkMu18")) passHLT_IsoTkMu18=true;
                if(TrigPath.Contains("HLT_IsoMu22")) passHLT_IsoMu22=true;
                if(TrigPath.Contains("HLT_IsoTkMu22")) passHLT_IsoTkMu22=true;
                
                
                if(TrigPath.Contains("HLT_Ele23_WPLoose_Gsf")) passHLT_Ele23_WPLoose_Gsf=true;
                if(TrigPath.Contains("HLT_Ele27_WPLoose_Gsf")) passHLT_Ele27_WPLoose_Gsf=true;
                if(TrigPath.Contains("HLT_Ele27_eta2p1_WPLoose_Gsf")||TrigPath.Contains("HLT_Ele27_eta2p1_WP75_Gsf")) passHLT_Ele27_eta2p1_WPLoose_Gsf=true;
                if(TrigPath.Contains("HLT_Ele32_eta2p1_WPLoose_Gsf")||TrigPath.Contains("HLT_Ele32_eta2p1_WP75_Gsf")) passHLT_Ele32_eta2p1_WPLoose_Gsf=true;
                if(TrigPath.Contains("HLT_TripleMu_12_10_5"))passHLT_TripleMu_12_10_5=true;
                if(TrigPath.Contains("HLT_DiMu9_Ele9_CaloIdL_TrackIdL"))passHLT_DiMu9_Ele9_CaloIdL_TrackIdL=true;
                if(TrigPath.Contains("HLT_Mu8_DiEle12_CaloIdL_TrackIdL"))passHLT_Mu8_DiEle12_CaloIdL_TrackIdL=true;
                if(TrigPath.Contains("HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL"))passHLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL=true;
                if(TrigPath.Contains("HLT_Ele23_CaloIdL_TrackIdL_IsoVL"))passHLT_Ele23_CaloIdL_TrackIdL_IsoVL=true;
                if(TrigPath.Contains("HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW"))passHLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW=true;
                
                
            }
        }
    }
    
    //std::cout<<"EVENT "<<_eventNb<<std::endl;
    if(debug) cout << "Debug 4.1 " <<endl;
    //============ Beamspot ============
    edm::Handle< reco::BeamSpot > theBeamSpot;
    if(isgetbylabel)    iEvent.getByLabel( IT_beamspot, theBeamSpot );
    else iEvent.getByToken(beamspot_token, theBeamSpot );
    if( ! theBeamSpot.isValid() ) ERR( IT_beamspot ) ;
    BeamSpot::Point  BS= theBeamSpot->position();;
    //==================================
    if(debug) cout << "Debug 4.2 " <<endl;
    //============ Primary vertices ============
    //edm::InputTag IT_goodVtx = edm::InputTag("goodOfflinePrimaryVertices");
    //    edm::InputTag IT_goodVtx = edm::InputTag("offlineSlimmedPrimaryVertices");
    edm::Handle<std::vector<Vertex> > theVertices;
    if(isgetbylabel)    iEvent.getByLabel( "goodOfflinePrimaryVertices", theVertices) ;
    else iEvent.getByToken(goodpv_token,theVertices) ;
    if( ! theVertices.isValid() ) ERR(IT_goodVtx ) ;
    int nvertex = theVertices->size();
    _n_PV = nvertex;
    
    //Nvtx->Fill(TMath::Min(nvertex,49));

    Vertex::Point PV(0,0,0); 
    if( nvertex){ PV = theVertices->begin()->position();
    const Vertex* PVtx = &((*theVertices)[0]);
    _PVchi2 = PVtx->chi2();
    _PVerr[0] = PVtx->xError();
    _PVerr[1] = PVtx->yError();
    _PVerr[2] = PVtx->zError();
    }

    //==================================
    if(debug) cout << "Debug 4.5 " <<endl;
    //============ Pat MET ============
    edm::Handle< vector<pat::MET> > ThePFMET;
    if(isgetbylabel) iEvent.getByLabel(IT_pfmet, ThePFMET);
    else iEvent.getByToken(met_token, ThePFMET); 
    if(debug) cout << "Debug 4.6" <<endl;
    if( ! ThePFMET.isValid() ) ERR( IT_pfmet );
    if(debug) cout << "Debug 4.7 " <<endl;
    const vector<pat::MET> *pfmetcol = ThePFMET.product();
    const pat::MET *pfmet;
    pfmet = &(pfmetcol->front());
    _met = pfmet->pt();
    _met_phi = pfmet->phi();

    //============ PF cand ============
    edm::Handle<pat::PackedCandidateCollection> pfcands;
    if(isgetbylabel)    iEvent.getByLabel("packedPFCandidates", pfcands);
    else iEvent.getByToken(pfcands_token ,pfcands);

    //==================================
    
    //============ Pat Muons ============
    edm::Handle< std::vector<pat::Muon> > thePatMuons;
    if(isgetbylabel) iEvent.getByLabel( IT_muon, thePatMuons );
    else iEvent.getByToken(patmuon_token, thePatMuons );
    if( ! thePatMuons.isValid() )  ERR(IT_muon) ;
    //==================================
    
    //============ Pat Electrons ============
    edm::Handle< std::vector<pat::Electron> > thePatElectrons;
    if(isgetbylabel)  iEvent.getByLabel( IT_electron, thePatElectrons );
    else iEvent.getByToken(patelectron_token, thePatElectrons );
    if( ! thePatElectrons.isValid() ) ERR( IT_electron );
    //==================================
    edm::Handle< edm::View<pat::Electron> > thePatElectronsView;
    if(isgetbylabel) iEvent.getByLabel( IT_electron, thePatElectronsView );
    else    iEvent.getByToken(patelectronview_token, thePatElectronsView );
    //if( ! thePatElectrons.isValid() ) ERR( IT_electron );
    //==================================  

    
    //============ Pat Taus ============
    edm::Handle< std::vector<pat::Tau> > thePatTaus;
    if(isgetbylabel) iEvent.getByLabel( IT_tau, thePatTaus );
    else iEvent.getByToken(pattau_token, thePatTaus );
    if( ! thePatTaus.isValid() )  ERR(IT_tau) ;
    //==================================



    if(debug) cout << "Debug 4.8 " <<endl;
    //============ Conversions ============
    edm::Handle< std::vector<reco::Conversion> > theConversions;
    if(isgetbylabel)    iEvent.getByLabel("reducedEgamma","reducedConversions", theConversions);
    else iEvent.getByToken(conv_token,theConversions);

    //==================================
    
    //============ Pat Jets ============
    edm::Handle< std::vector< pat::Jet> > thePatJets;
    if(isgetbylabel)iEvent.getByLabel(IT_jet , thePatJets );
    else iEvent.getByToken(jet_token,thePatJets );
    if( ! thePatJets.isValid() ) ERR(IT_jet);
    //==================================
    
    //============ Rho ============
    edm::Handle<double> rhoJets;
    if(isgetbylabel) iEvent.getByLabel(edm::InputTag("fixedGridRhoFastjetAll","") , rhoJets);//kt6PFJets    
    else iEvent.getByToken(rhoJets_token,rhoJets);
    
    myRhoJets = *rhoJets;
    //==================================
    

    //============ RhoNC ============                             
    edm::Handle<double> rhoJetsNC;
    if(isgetbylabel)     iEvent.getByLabel(edm::InputTag("fixedGridRhoFastjetCentralNeutral","") , rhoJetsNC);//kt6PFJets         
    else     iEvent.getByToken(rhoJetsNC_token,rhoJetsNC);

    myRhoJetsNC = *rhoJetsNC;
    //===============================
    if(debug) cout << "Debug 4.9 " <<endl;

    // Get MVA values and categories (optional)
    //============ MVA ============                                                                               
    edm::Handle<edm::ValueMap<float> > mvaValues;
    edm::Handle<edm::ValueMap<int> > mvaCategories;
    iEvent.getByToken(mvaValuesMapToken_,mvaValues);
    iEvent.getByToken(mvaCategoriesMapToken_,mvaCategories);
    //================================== 

    if(debug) cout << "Debug 4.10 " <<endl;

    //==================================    
    //MET calcuation
    //in miniAOD v2
    double rawmetX = pfmet->uncorPx();
    double rawmetY = pfmet->uncorPy();
    
    double corrMetx = rawmetX;
    double corrMety = rawmetY;
    
    double corrMetx_up = rawmetX;
    double corrMety_up = rawmetY;
    double corrMetx_down = rawmetX;
    double corrMety_down = rawmetY;
    
    //    cout << "run,LS,event nb "<<_runNb <<":"<<_lumiBlock<<":"<<_eventNb<<endl;
    //cout << "Raw MET: " <<  sqrt(rawmetX*rawmetX+rawmetY*rawmetY)<<endl;

    for( std::vector<pat::Jet>::const_iterator jet = (*thePatJets).begin(); jet != (*thePatJets).end(); jet++ ) {

      TLorentzVector pJet;
      pJet.SetPtEtaPhiE(((&*jet)->correctedP4("Uncorrected")).Pt(), ((&*jet)->correctedP4("Uncorrected")).Eta(),((&*jet)->correctedP4("Uncorrected")).Phi(), ((&*jet)->correctedP4("Uncorrected")).E());
      //cout<<"**********"<<endl;
      //cout << "jet pt, eta, phi, neutral EM, charged EM " <<  (&*jet)->pt() << ", " <<  (&*jet)->eta()<<", "<< (&*jet)->phi()<< ", " << (&*jet)->neutralEmEnergy() <<", "<<(&*jet)->chargedEmEnergy()<< endl; 
      //      cout << "jet pt, eta, phi (uncorrected): "<<  pJet.Pt()<<", "<<pJet.Eta()<<", "<< pJet.Phi()<<endl;

      const std::vector<reco::CandidatePtr> & cands = (&*jet)->daughterPtrVector();
      for ( std::vector<reco::CandidatePtr>::const_iterator cand = cands.begin();
	    cand != cands.end(); ++cand ) {
          const reco::PFCandidate *pfcand = dynamic_cast<const reco::PFCandidate *>(cand->get());
          const reco::Candidate *mu = (pfcand != 0 ? ( pfcand->muonRef().isNonnull() ? pfcand->muonRef().get() : 0) : cand->get());
          if ( mu != 0  && ( mu->isGlobalMuon() || mu->isStandAloneMuon() ) ) {
              //	  reco::Candidate::LorentzVector muonP4 = (*cand)->p4();
              TLorentzVector pMu;
              pMu.SetPtEtaPhiE((*cand)->pt(),(*cand)->eta(),(*cand)->phi(),(*cand)->energy());
              //	  rawJetP4 -= muonP4;
              //  cout << "mu found new method: "<< muonP4.Pt()<< ", "<<muonP4.Eta()<< ", "<<muonP4.Phi()<<endl;
              pJet-= pMu;
          }
      }

      //      cout << "jet pt, eta, phi (uncorrected, but muons removed): "<<  pJet.Pt()<<", "<<pJet.Eta()<<", "<< pJet.Phi()<<endl;
      if(debug) cout << "Debug 4.11 " <<endl;
      std::pair <double, double> corr = fMetCorrector->getCorrections(
                                                                      ((&*jet)->correctedP4("Uncorrected")).Pt(),
                                                                      ((&*jet)->correctedP4("Uncorrected")).Eta(),
                                                                      pJet.Pt(),
                                                                      pJet.Phi(),
                                                                      (&*jet)->neutralEmEnergyFraction()+(&*jet)->chargedEmEnergyFraction(),
                                                                      myRhoJets,
                                                                      (&*jet)->jetArea());
        corrMetx += corr.first ;
        corrMety += corr.second;
	//	cout << "Corr x,y: "<< corr.first<<", "<<corr.second<<endl;
	
      std::pair <double, double> corr_up = fMetCorrector->getCorrections(
                                                                      ((&*jet)->correctedP4("Uncorrected")).Pt(),
                                                                      ((&*jet)->correctedP4("Uncorrected")).Eta(),
                                                                      pJet.Pt(),
                                                                      pJet.Phi(),
                                                                      (&*jet)->neutralEmEnergyFraction()+(&*jet)->chargedEmEnergyFraction(),
                                                                      myRhoJets,
                                                                      (&*jet)->jetArea(),"up");
      corrMetx_up += corr_up.first ;
      corrMety_up += corr_up.second;
      
      std::pair <double, double> corr_down = fMetCorrector->getCorrections(
                                                                      ((&*jet)->correctedP4("Uncorrected")).Pt(),
                                                                      ((&*jet)->correctedP4("Uncorrected")).Eta(),
                                                                      pJet.Pt(),
                                                                      pJet.Phi(),
                                                                      (&*jet)->neutralEmEnergyFraction()+(&*jet)->chargedEmEnergyFraction(),
                                                                      myRhoJets,
                                                                      (&*jet)->jetArea(),"down");
      corrMetx_down += corr_down.first ;
      corrMety_down += corr_down.second;



    }
    double newmet    = sqrt(corrMetx*corrMetx + corrMety*corrMety);
    double newmetphi = atan2(corrMety, corrMetx);
    
    double newmet_up    = sqrt(corrMetx_up*corrMetx_up + corrMety_up*corrMety_up);
    double newmet_down    = sqrt(corrMetx_down*corrMetx_down + corrMety_down*corrMety_up);
    double newmetphi_up = atan2(corrMety_up, corrMetx_up);
    double newmetphi_down = atan2(corrMety_down, corrMetx_down);
    if(debug) cout << "Debug 4.12 " <<endl;
    
    _met = newmet;
    _met_phi = newmetphi;
    _met_JECup = newmet_up;
    _met_phi_JECup = newmetphi_up;
    _met_JECdown = newmet_down;
    _met_phi_JECdown = newmetphi_down;
    

 
    if(debug) cout << "Debug 6 " <<endl;
    enum decay {
        W_L,  // 0
        W_T_L, // 1
        W_B_L, // 2
        W_B_D_L, //3
        W_B_D_T_L, // 4
        W_B_T_L, // 5
        W_D_L, // 6
        W_D_T_L, //7
        B_L, // 8
        B_D_L, //9
        B_D_T_L, //10
        B_T_L,  // 11
        D_L, //12
        D_T_L, //13
        B_Baryon, // 14
        C_Baryon, //15
        pi_0, //16
        photon_, //17
        F_L, //18
        N_U_L_L // 19
    };
    

    
    std::vector<const pat::Muon* > sMu = ssbLooseMuonSelector( *thePatMuons, _muonMinPt, PV, _looseD0Mu, true);
    std::vector<const pat::Electron* > sEl = ssbMVAElectronSelector( *thePatElectrons, _eleMinPt, PV, _looseD0E, _chargeConsistency, _useconversions, theConversions, BS, false);
    std::vector<const pat::Tau* > sTau = ssbLooseTauSelector( *thePatTaus, _tauMinPt, PV, 10000, true);
    std::vector<double > sTau_dz;
    for (const pat::Tau &tau : *thePatTaus) {
        if (tau.pt() < _tauMinPt) continue;
        if (fabs(tau.eta()) > _tauEta) continue;
        
        Vertex::Point vertex = theVertices->begin()->position(); //PV
        Vertex::Point vtx = (*(tau.leadChargedHadrCand())).vertex();
        TLorentzVector p4;
        p4.SetPtEtaPhiE(tau.pt(),tau.eta(),tau.phi(),tau.energy());
        double dzTau = (vtx.z()-vertex.z()) - ((vtx.x()-vertex.x())*p4.X()+(vtx.y()-vertex.y())*p4.Y())/ p4.Pt() *  p4.Z()/ p4.Pt();
        //std::cout<<"tau sig "<<tau.dxy_Sig()<<std::endl;
        //std::cout<<"tau dxy "<<tau.dxy()<<std::endl;
        //std::cout<<"tau dz "<<dzTau<<std::endl;;
        sTau_dz.push_back(dzTau);
    }
    
    SelectedJetsAll = JetSelectorAll(*thePatJets, 5., 3.0);
    //    cout << "jetptcut " << _jetPtCut<<endl;
    std::vector<const pat::Jet* > SelectedJets = JetSelector(*thePatJets, _jetPtCut, _jetEtaCut);
    
    //    std::cout<<"tot lept size " << sEl.size() + sMu.size()<<std::endl;
    



    
    /*
     *
     *Mu
     *
     *
     */

    //Begin modif takeimai Trigger
    //std::cout<<sMu.size()<<" "<<sEl.size()<<std::endl;
    int _sameSign[2][2] = {{0, 0}, {0, 0}};
    if(debug) cout << "Debug 7 " <<endl;
    int leptonCounter = 0;
    for(unsigned int i = 0 ; i < sMu.size() ;i++ ){

      
        const pat::Muon *iM = sMu[i];
        
	    //if((iM->innerTrack()->ptError())/(iM->innerTrack()->pt()) > 0.2) continue;
        
        //Reduce Ntuple size
        if (getPFIsolation(pfcands, dynamic_cast<const reco::Candidate *>(iM), 0.05, 0.2, 10., false, false, myRhoJetsNC)>0.4) continue;
        
        _leptonIndex = i;
        
        if (leptonCounter == nLeptonsMax) continue;
	InitLepton();
	


	//PassTriggerLeg(muontriggerlegs,iM,iEvent);
	 
	if(Skim=="skim2LrecoandTriggerTP" || Skim=="skim2GenL"){
	passleg1L[leptonCounter]=PassTriggerLeg("hltL3crIsoL1sMu20Eta2p1L1f0L2f10QL3f24QL3trkIsoFiltered0p09",iM,iEvent); 
	hltL1sL1TripleMu553[leptonCounter]=PassTriggerLeg("hltL1sL1TripleMu553","hltL1sTripleMu0IorTripleMu553",iM,iEvent);

	hltL1TripleMu553L2TriMuFiltered3L3TriMuFiltered12105[leptonCounter]=PassTriggerLeg("hltL1TripleMu553L2TriMuFiltered3L3TriMuFiltered12105",iM,iEvent);
	hltL1TripleMu553L2TriMuFiltered3L3TriMuFiltered10105[leptonCounter]=PassTriggerLeg("hltL1TripleMu553L2TriMuFiltered3L3TriMuFiltered10105",iM,iEvent);
	hltL1TripleMu553L2TriMuFiltered3L3TriMuFiltered5[leptonCounter]=PassTriggerLeg("hltL1TripleMu553L2TriMuFiltered3L3TriMuFiltered5",iM,iEvent);
	hltDiMu9Ele9CaloIdLTrackIdLMuonlegL3Filtered9[leptonCounter]=PassTriggerLeg("hltDiMu9Ele9CaloIdLTrackIdLMuonlegL3Filtered9",iM,iEvent);
	hltMu8DiEle12CaloIdLTrackIdLMuonlegL3Filtered8[leptonCounter]=PassTriggerLeg("hltMu8DiEle12CaloIdLTrackIdLMuonlegL3Filtered8",iM,iEvent);
	hltL1sL1Mu6HTT150ORL1Mu8HTT125ORL1HTT175[leptonCounter]=PassTriggerLeg("hltL1sL1Mu8HTT","hltL1sL1Mu6HTT",iM,iEvent) || PassTriggerLeg("hltL1sMu8HTT","hltL1sMu6HTT",iM,iEvent);
	hltDoubleMu8Mass8L3Filtered[leptonCounter]=PassTriggerLeg("hltDoubleMu8Mass8L3Filtered",iM,iEvent);
	hltMuon8L3Filtered0[leptonCounter]=PassTriggerLeg("hltMuon8L3Filtered0",iM,iEvent);
	hltL1sL1DoubleMu103p5ORDoubleMu125[leptonCounter]=PassTriggerLeg("hltL1sL1DoubleMu103p5ORDoubleMu125","hltL1sDoubleMu114IorDoubleMu125",&*iM,iEvent);

	hltDiMuonGlb17Glb8RelTrkIsoFiltered0p4[leptonCounter]=PassTriggerLeg("hltDiMuonGlb17Glb8RelTrkIsoFiltered0p4",iM,iEvent);
	hltL3fL1sDoubleMu103p5L1f0L2f10OneMuL3Filtered17[leptonCounter]=PassTriggerLeg("hltL3fL1sDoubleMu103p5L1f0L2f10OneMuL3Filtered17","hltL3fL1sDoubleMu114L1f0L2f10OneMuL3Filtered17",iM,iEvent);

	hltL3pfL1sDoubleMu103p5L1f0L2pf0L3PreFiltered8[leptonCounter]=PassTriggerLeg("hltL3pfL1sDoubleMu103p5L1f0L2pf0L3PreFiltered8","hltL3pfL1sDoubleMu114L1f0L2pf0L3PreFiltered8",iM,iEvent);

	hltL3fL1sDoubleMu103p5L1f0L2f10L3Filtered17[leptonCounter]=PassTriggerLeg("hltL3fL1sDoubleMu103p5L1f0L2f10L3Filtered17","hltL3fL1sDoubleMu114L1f0L2f10L3Filtered17",iM,iEvent);
	                                                                                                                         
	hltDiMuonGlbFiltered17TrkFiltered8[leptonCounter]=PassTriggerLeg("hltDiMuonGlbFiltered17TrkFiltered8",iM,iEvent);
	hltDiMuonGlb17Trk8RelTrkIsoFiltered0p4[leptonCounter]=PassTriggerLeg("hltDiMuonGlb17Trk8RelTrkIsoFiltered0p4",iM,iEvent);
	hltDiMuonGlb17Glb8RelTrkIsoFiltered0p4DzFiltered0p2[leptonCounter]=PassTriggerLeg("hltDiMuonGlb17Glb8RelTrkIsoFiltered0p4DzFiltered0p2",iM,iEvent);
	hltDiMuonGlb17Trk8RelTrkIsoFiltered0p4DzFiltered0p2[leptonCounter]=PassTriggerLeg("hltDiMuonGlb17Trk8RelTrkIsoFiltered0p4DzFiltered0p2",iM,iEvent);
	


	hltL1sL1Mu20EG10[leptonCounter]=PassTriggerLeg("hltL1sL1Mu20EG10","hltL1sMu20EG10",iM,iEvent);
	hltMu23TrkIsoVVLEle8CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered23[leptonCounter]=PassTriggerLeg("hltMu23TrkIsoVVLEle8CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered23",iM,iEvent);
	hltL1sL1Mu5EG20[leptonCounter]=PassTriggerLeg("hltL1sL1Mu5EG20","hltL1sMu5EG20IorMu5IsoEG18",iM,iEvent);
	hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered8[leptonCounter]=PassTriggerLeg("hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered8",iM,iEvent);

	hltL1sL1Mu12EG10[leptonCounter]=PassTriggerLeg("hltL1sL1Mu12EG10","hltL1sMu12EG10",iM,iEvent);
	hltMu17TrkIsoVVLEle12CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered17[leptonCounter]=PassTriggerLeg("hltMu17TrkIsoVVLEle12CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered17",iM,iEvent);
	hltL1sL1Mu5EG15[leptonCounter]=PassTriggerLeg("hltL1sL1Mu5EG15","hltL1sMu5EG15",iM,iEvent);
	hltMu8TrkIsoVVLEle17CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered8[leptonCounter]=PassTriggerLeg("hltMu8TrkIsoVVLEle17CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered8",iM,iEvent);


	hltL3crIsoL1sMu20Eta2p1L1f0L2f10QL3f24QL3trkIsoFiltered0p09[leptonCounter]=PassTriggerLeg("hltL3crIsoL1sMu20Eta2p1L1f0L2f10QL3f24QL3trkIsoFiltered0p09",iM,iEvent);
	hltL3crIsoL1sMu16L1f0L2f10QL3f20QL3trkIsoFiltered0p09[leptonCounter]=PassTriggerLeg("hltL3crIsoL1sMu16L1f0L2f10QL3f20QL3trkIsoFiltered0p09","hltL3crIsoL1sMu18L1f0L2f10QL3f20QL3trkIsoFiltered0p09",iM,iEvent);

	hltL3crIsoL1sMu16L1f0L2f10QL3f18QL3trkIsoFiltered0p09[leptonCounter]=PassTriggerLeg("hltL3crIsoL1sMu16L1f0L2f10QL3f18QL3trkIsoFiltered0p09",iM,iEvent);
	hltL3fL1sMu16L1f0Tkf18QL3trkIsoFiltered0p09[leptonCounter]=PassTriggerLeg("hltL3fL1sMu16L1f0Tkf18QL3trkIsoFiltered0p09",iM,iEvent);

	hltL3crIsoL1sMu20L1f0L2f10QL3f22QL3trkIsoFiltered0p09[leptonCounter]=PassTriggerLeg("hltL3crIsoL1sMu20L1f0L2f10QL3f22QL3trkIsoFiltered0p09",iM,iEvent);
	hltL3fL1sMu20L1f0Tkf22QL3trkIsoFiltered0p09[leptonCounter]=PassTriggerLeg("hltL3fL1sMu20L1f0Tkf22QL3trkIsoFiltered0p09",iM,iEvent);


	
	hltL1sL1SingleMu16ORSingleMu25[leptonCounter]=PassTriggerLeg("hltL1sL1SingleMu16ORSingleMu25",&*iM,iEvent);
	hltDiMuonGlb27Trk8DzFiltered0p2[leptonCounter]=PassTriggerLeg("hltDiMuonGlb27Trk8DzFiltered0p2",&*iM,iEvent);
	hltL3fL1sMu16orMu25L1f0L2f25L3Filtered27[leptonCounter]=PassTriggerLeg("hltL3fL1sMu16orMu25L1f0L2f25L3Filtered27",&*iM,iEvent);
	hltL3fL1sMu16orMu25L1f0L2f10QL3Filtered30Q[leptonCounter]=PassTriggerLeg("hltL3fL1sMu16orMu25L1f0L2f10QL3Filtered30Q",&*iM,iEvent) || PassTriggerLeg("hltL3fL1sMu16orMu25L1f0L2f16QL3Filtered30Q",&*iM,iEvent);
	hltL3fL1sMu22orMu25orMu20EG15orMu5EG20L1f0L2f10QL3Filtered30Q[leptonCounter]=PassTriggerLeg("ltL3fL1sMu22orMu25orMu20EG15orMu5EG20L1f0L2f10QL3Filtered30Q",&*iM,iEvent) ;
	
	hltL1sSingleMu16[leptonCounter]=PassTriggerLeg("hltL1sSingleMu16",&*iM,iEvent);
	hltL1sSingleMu18[leptonCounter]=PassTriggerLeg("hltL1sSingleMu18",&*iM,iEvent);
	hltL1sSingleMu20[leptonCounter]=PassTriggerLeg("hltL1sSingleMu20",&*iM,iEvent);
	hltL1sSingleMu20erlorSingleMu22lorSingleMu25[leptonCounter]=PassTriggerLeg("hltL1sSingleMu20erlorSingleMu22lorSingleMu25",&*iM,iEvent);

	hltL3fL1sL1DoubleMu0ETM40lorDoubleMu0ETM55[leptonCounter]=PassTriggerLeg("hltL3fL1sL1DoubleMu0ETM40lorDoubleMu0ETM55",&*iM,iEvent);
	hltL3crIsoL1sSingleMu18erIorSingleMu20erL1f0L2f10QL3f19QL3trkIsoFiltered0p09[leptonCounter]=PassTriggerLeg("hltL3crIsoL1sSingleMu18erIorSingleMu20erL1f0L2f10QL3f19QL3trkIsoFiltered0p09",&*iM,iEvent);
	
	}
	//end modif takeimai


        _flavors[leptonCounter] = 1;
        _charges[leptonCounter] = iM->charge();
        _isolation[leptonCounter] = pfRelIso(iM,myRhoJetsNC);
        
        _lchargePixSc[leptonCounter] =0;
        _lchargeGSF[leptonCounter] = 0;
        _lchargeCTF[leptonCounter] = 0;
        _lnmisshits[leptonCounter] = -10;
        _passconvveto[leptonCounter] = false;
        _lsietaieta[leptonCounter] = -1;
        _lfull5x5sietaieta[leptonCounter] = -1;
        _lhovere[leptonCounter] = -1;
        _ldphiin[leptonCounter] = -1;
        _ldetain[leptonCounter] = -1;
        _l1oemin1op[leptonCounter] = -1;
        
        double chargedHadronIso = iM->pfIsolationR03().sumChargedHadronPt;
        double neutralHadronIso = iM->pfIsolationR03().sumNeutralHadronEt;
        double photonIso = iM->pfIsolationR03().sumPhotonEt;
	
        double beta = iM->pfIsolationR03().sumPUPt;
        double pfRelIsoMu  = ( chargedHadronIso + TMath::Max ( 0.0 ,neutralHadronIso + photonIso - 0.5 * beta ) )/iM->pt() ;
	
	    _muonSegmentComp[leptonCounter] = iM->segmentCompatibility();
        _pfisocharged[leptonCounter]= chargedHadronIso/iM->pt();
        _pfisoneutral[leptonCounter] = neutralHadronIso/iM->pt();
        _pfisophoton[leptonCounter] = photonIso/iM->pt();
        
        _ipPV[leptonCounter] =iM->innerTrack()->dxy(PV);
        _ipPVerr[leptonCounter] = iM->innerTrack()->dxyError();
        
        _ipZPV[leptonCounter] = iM->innerTrack()->dz(PV);
        _ipZPVerr[leptonCounter] = iM->innerTrack()->dzError();
        
        _3dIP[leptonCounter]    = iM->dB(pat::Muon::PV3D);
        _3dIPerr[leptonCounter] = iM->edB(pat::Muon::PV3D);
        _3dIPsig[leptonCounter] = fabs(_3dIP[leptonCounter]/_3dIPerr[leptonCounter]);
        
        bool goodGlb = iM->isGlobalMuon() && iM->globalTrack()->normalizedChi2() < 3
        && iM->combinedQuality().chi2LocalPosition < 12 && iM->combinedQuality().trkKink < 20;
        bool good = iM->innerTrack()->validFraction() >= 0.8 && iM->segmentCompatibility() >= (goodGlb ? 0.303 : 0.451);
        bool isMedium = iM->isLooseMuon() && good;
        
        _isMediumMuon[leptonCounter] = isMedium;
        _istightID[leptonCounter] = good;
        _istrigemulID[leptonCounter] =triggerEmulator(iM);
        _istrigemulISO[leptonCounter] =isoTriggerEmulator(iM);
	
        ((TLorentzVector *)_leptonP4->At(leptonCounter))->SetPtEtaPhiE(iM->pt(), iM->eta(), iM->phi(), iM->energy());
        
        _mt[leptonCounter] = MT_calc(*((TLorentzVector *)_leptonP4->At(leptonCounter)), _met, _met_phi);
        
        _lPt[leptonCounter] = ((TLorentzVector *)_leptonP4->At(leptonCounter))->Pt();
        _lEta[leptonCounter] = ((TLorentzVector *)_leptonP4->At(leptonCounter))->Eta();
        _lPhi[leptonCounter] = ((TLorentzVector *)_leptonP4->At(leptonCounter))->Phi();
        _lE[leptonCounter] = ((TLorentzVector *)_leptonP4->At(leptonCounter))->E();
        
        _miniisolation_0p2[leptonCounter] = getPFIsolation(pfcands, dynamic_cast<const reco::Candidate *>(iM), 0.05, 0.2, 10., false, false, myRhoJetsNC);
	_miniisolationcharged_0p2[leptonCounter] = getPFIsolation(pfcands, dynamic_cast<const reco::Candidate *>(iM), 0.05, 0.2, 10., false, true, myRhoJetsNC);
        _miniisolation_0p3[leptonCounter] = getPFIsolation(pfcands, dynamic_cast<const reco::Candidate *>(iM), 0.05, 0.3, 10., false, false, myRhoJetsNC);
        _miniisoneutral[leptonCounter] = _miniisolation_0p2[leptonCounter]-_miniisolationcharged_0p2[leptonCounter];
        
	
        fillCloseJetVars(leptonCounter,PV);
   
	if( nvertex){
        
        const Vertex PVtx = ((*theVertices)[0]);
        _issoftmuon[leptonCounter] =SoftTightID(iM,PV, PVtx ) ;}
        _issoftmuonwithiso[leptonCounter] =_issoftmuon[leptonCounter] && pfRelIsoMu <0.5 && pfRelIsoMu*iM->pt()<5;
        
        _isloose[leptonCounter] = ( _istightID[leptonCounter] && _miniisolation_0p2[leptonCounter] < 0.4);
        _istightIso[leptonCounter] =TightIso("RA5",iM, _miniisolation_0p2[leptonCounter],_ptratio[leptonCounter],_ptrel2[leptonCounter]);
        _istightID[leptonCounter]= TightID("RA5",iM,PV,_miniisolation_0p2[leptonCounter]);
        
        
        _istightIDWP2016[leptonCounter] =_istightIso[leptonCounter]&&_istightID[leptonCounter];
        _istightIDWP2016_noIP3D[leptonCounter] =_istightIso[leptonCounter]&&TightID("RA5",iM,PV,_miniisolation_0p2[leptonCounter],false);
        _istightIDWP2016_IsoEmul[leptonCounter]=_istightIDWP2016[leptonCounter]  && _istrigemulISO[leptonCounter];
        _istightIDWP2016_RA7[leptonCounter]= TightIso("RA7",iM, _miniisolation_0p2[leptonCounter],_ptratio[leptonCounter],_ptrel2[leptonCounter])&&TightID("RA7",iM,PV,_miniisolation_0p2[leptonCounter]);
        _istightIDWP2016_EWK[leptonCounter]= TightID_EWK("EWK",iM,PV,_lmva[leptonCounter],_miniisolation_0p2[leptonCounter]);
        
        _isVetoIDWP2016_RA7[leptonCounter]= VetoID(iM,PV, _miniisolation_0p2[leptonCounter]);
        _isVetoIDWP2016_EWK[leptonCounter]=VetoID_EWK(iM,PV, _miniisolation_0p2[leptonCounter]);
        
        _isFOIDWP2016[leptonCounter] =FOID("RA5",iM,PV, _miniisolation_0p2[leptonCounter]);
        _isFOIDWP2016_RA7[leptonCounter] =FOID("RA7",iM,PV, _miniisolation_0p2[leptonCounter]);
        _isFOIDWP2016_EWK[leptonCounter]=FOID_EWK("EWK",iM,PV,_lmva[leptonCounter],_ptratio[leptonCounter],_closeJetCSVAll[leptonCounter],_miniisolation_0p2[leptonCounter]);
        _isFOIDWP2016_IsoEmul[leptonCounter] = _isFOIDWP2016[leptonCounter] && _istrigemulISO[leptonCounter];
        _isFOIDWP2016_noIP3D[leptonCounter] =FOID("RA5",iM,PV, _miniisolation_0p2[leptonCounter],false);
        
        _islooseMT2[leptonCounter] =MuonLooseID_MiniIso(iM,PV, _miniisolation_0p2[leptonCounter]);


        if (IsMC) {
            //**************************************************************************************
            // MC
            //**************************************************************************************
            
            MatchToGenParticle(iEvent,leptonCounter,iM->pt(),iM->eta(), iM->phi(), iM->pdgId(),genpartalreadymatched);
            
            fillIsoMCVars(leptonCounter);
            if (_closeIndex[leptonCounter] >=0)
                matchCloseJet(leptonCounter);
            
            if (_regression) {
                std::vector <double> regVars = RegressionVars(SelectedJetsAll[_closeIndex[leptonCounter]], _closeJetPtAllMC[leptonCounter], iM);
                for (int k=0; k!=15; ++k) {
                    _regVars[k] = regVars[k];
                }
                _regVars[11] = fMetCorrector->getJECUncertainty(SelectedJetsAll[_closeIndex[leptonCounter]]->pt(),SelectedJetsAll[_closeIndex[leptonCounter]]->eta());
            }
            
            //**************************************************************************************
            // MC *
            //**************************************************************************************
        }
        
        if (_regression)
            fillRegVars(SelectedJetsAll.at(_closeIndex[leptonCounter]), _closeJetPtAllMC[leptonCounter], iM);
        
        _sameSign[_istight[leptonCounter]][int((_charges[leptonCounter]+1)/2)]++;
        
        _mllZ[leptonCounter] = 9999.;
        _mllG[leptonCounter] = 9999.;
        leptonCounter++;
        
    }

    _nMu = leptonCounter;
    //std::cout<<leptonCounter<<" "<<((TLorentzVector *)_leptonP4->At(leptonCounter))->Pt() <<std::endl;
    


    /*
     *
     *ELE
     *
     *
     */

    //    for(unsigned int i = 0 ; i < sEl.size() ;i++ ){
 
      for (size_t i = 0; i < thePatElectronsView->size(); ++i){
          
          const auto iE = thePatElectronsView->ptrAt(i);
          if (!ssbMVAElectronSelectorPassed(iE, _eleMinPt, PV,_looseD0E, _chargeConsistency, _useconversions, theConversions, BS, false)) continue;
          
          //Reduce Ntuple size
          if (getPFIsolation(pfcands, dynamic_cast<const reco::Candidate *>(&*iE), 0.05, 0.2, 10., false, false, myRhoJetsNC)>0.4) continue;
          
          if(leptonCounter == nLeptonsMax ) cout << "leptonCounter == nLeptonsMax" <<endl;
        
          _leptonIndex = i;
          if (leptonCounter == nLeptonsMax) continue;
	
          InitLepton();
        

	if(Skim=="skim2LrecoandTriggerTP" || Skim=="skim2GenL"){
	passleg1L[leptonCounter]=PassTriggerLeg("hltEle32WP75GsfTrackIsoFilter",&*iE,iEvent) || PassTriggerLeg("hltEle32WPLooseGsfTrackIsoFilter",&*iE,iEvent); 
	hltDiMu9Ele9CaloIdLTrackIdLElectronlegDphiFilter[leptonCounter]=PassTriggerLeg("hltDiMu9Ele9CaloIdLTrackIdLElectronlegDphiFilter",&*iE,iEvent);
	hltMu8DiEle12CaloIdLTrackIdLElectronlegDphiFilter[leptonCounter]=PassTriggerLeg("hltMu8DiEle12CaloIdLTrackIdLElectronlegDphiFilter",&*iE,iEvent);
	
	hltL1sL1TripleEG14108[leptonCounter]=PassTriggerLeg("hltL1sL1TripleEG14108","hltL1sTripleEG14108",&*iE,iEvent);
	hltEle16Ele12Ele8CaloIdLTrackIdLDphiLeg1Filter[leptonCounter]=PassTriggerLeg("hltEle16Ele12Ele8CaloIdLTrackIdLDphiLeg1Filter",&*iE,iEvent);
	hltEle16Ele12Ele8CaloIdLTrackIdLDphiLeg2Filter[leptonCounter]=PassTriggerLeg("hltEle16Ele12Ele8CaloIdLTrackIdLDphiLeg2Filter",&*iE,iEvent);
	hltEle16Ele12Ele8CaloIdLTrackIdLDphiLeg3Filter[leptonCounter]=PassTriggerLeg("hltEle16Ele12Ele8CaloIdLTrackIdLDphiLeg3Filter",&*iE,iEvent);
	
                                                                                                                       
	hltL1sL1DoubleEG6HTT150orL1HTT175[leptonCounter]=PassTriggerLeg("hltL1sL1DoubleEG6HTT","hltL1sDoubleEG6HTT",&*iE,iEvent);
	hltDoubleEle8CaloIdMGsfTrackIdMDphiFilter[leptonCounter]=PassTriggerLeg("hltDoubleEle8CaloIdMGsfTrackIdMDphiFilter",&*iE,iEvent);
	hltMu8Ele8CaloIdMGsfTrackIdMDphiFilter[leptonCounter]=PassTriggerLeg("hltMu8Ele8CaloIdMGsfTrackIdMDphiFilter",&*iE,iEvent);
	
	
	hltL1sL1DoubleEG2210[leptonCounter]=PassTriggerLeg("hltL1sL1DoubleEG2210","hltL1sSingleEG40IorSingleIsoEG22erIorSingleIsoEG24IorDoubleEG1817IorDoubleEG2212",&*iE,iEvent);
	hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter[leptonCounter]=PassTriggerLeg("hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter",&*iE,iEvent);
	hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter[leptonCounter]=PassTriggerLeg("hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter",&*iE,iEvent);
	hltEle23Ele12CaloIdLTrackIdLIsoVLDZFilter[leptonCounter]=PassTriggerLeg("hltEle23Ele12CaloIdLTrackIdLIsoVLDZFilter",&*iE,iEvent);

	hltL1sL1DoubleEG1510[leptonCounter]=PassTriggerLeg("hltL1sL1DoubleEG1510","hltL1sEle23Ele12CaloIdLTrackIdLIsoVL",&*iE,iEvent);
	hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter[leptonCounter]=PassTriggerLeg("hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter",&*iE,iEvent);
	hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter[leptonCounter]=PassTriggerLeg("hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter",&*iE,iEvent);
	hltEle17Ele12CaloIdLTrackIdLIsoVLDZFilter[leptonCounter]=PassTriggerLeg("hltEle17Ele12CaloIdLTrackIdLIsoVLDZFilter",&*iE,iEvent);

	hltL1sL1Mu20EG10[leptonCounter]=PassTriggerLeg("hltL1sL1Mu20EG10","hltL1sMu20EG10",&*iE,iEvent);
	hltMu23TrkIsoVVLEle8CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter[leptonCounter]=PassTriggerLeg("hltMu23TrkIsoVVLEle8CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter",&*iE,iEvent);
	hltL1sL1Mu5EG20[leptonCounter]=PassTriggerLeg("hltL1sL1Mu5EG20","hltL1sMu5EG20IorMu5IsoEG18",&*iE,iEvent);
	hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter[leptonCounter]=PassTriggerLeg("hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter",&*iE,iEvent);


	hltL1sL1Mu12EG10[leptonCounter]=PassTriggerLeg("hltL1sL1Mu12EG10","hltL1sMu12EG10",&*iE,iEvent);
	hltMu17TrkIsoVVLEle12CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter[leptonCounter]=PassTriggerLeg("hltMu17TrkIsoVVLEle12CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter",&*iE,iEvent);
	hltL1sL1Mu5EG15[leptonCounter]=PassTriggerLeg("hltL1sL1Mu5EG15","hltL1sMu5EG15",&*iE,iEvent);
	hltMu8TrkIsoVVLEle17CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter[leptonCounter]=PassTriggerLeg("hltMu8TrkIsoVVLEle17CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter",&*iE,iEvent);


	hltEle27noerWPLooseGsfTrackIsoFilter[leptonCounter]=PassTriggerLeg("hltEle27noerWPLooseGsfTrackIsoFilter",&*iE,iEvent);
	hltEle27WPLooseGsfTrackIsoFilter[leptonCounter]=PassTriggerLeg("hltEle27WPLooseGsfTrackIsoFilter",&*iE,iEvent)||PassTriggerLeg("hltEle27WP75GsfTrackIsoFilter",&*iE,iEvent);
	hltEle32WPLooseGsfTrackIsoFilter[leptonCounter]=PassTriggerLeg("hltEle32WPLooseGsfTrackIsoFilter",&*iE,iEvent)||PassTriggerLeg("hltEle32WP75GsfTrackIsoFilter",&*iE,iEvent);
	
	hltEle23WPLooseGsfTrackIsoFilter[leptonCounter]=PassTriggerLeg("hltEle23WPLooseGsfTrackIsoFilter",&*iE,iEvent) ||PassTriggerLeg("hltEle23CaloIdLTrackIdLIsoVLTrackIsoFilter",&*iE,iEvent);

	hltEle23CaloIdLTrackIdLIsoVLTrackIsoFilter[leptonCounter]=PassTriggerLeg("hltEle23CaloIdLTrackIdLIsoVLTrackIsoFilter",&*iE,iEvent);
	hltDiEle33CaloIdLNewPixelMatchUnseededFilter[leptonCounter]=PassTriggerLeg("hltDiEle33CaloIdLNewPixelMatchUnseededFilter","hltDiEle33CaloIdLGsfTrkIdVLMWPMS2UnseededFilter",&*iE,iEvent);
	                                                                            
	hltEle30CaloIdLGsfTrkIdVLDPhiUnseededFilter[leptonCounter]=PassTriggerLeg("hltEle30CaloIdLGsfTrkIdVLDPhiUnseededFilter",&*iE,iEvent);

	hltDiEle33CaloIdLGsfTrkIdVLDPhiUnseededFilter[leptonCounter]=PassTriggerLeg("hltDiEle33CaloIdLGsfTrkIdVLDPhiUnseededFilter",&*iE,iEvent);

hltL1sSingleEG40IorSingleIsoEG22erIorSingleIsoEG24erIorSingleIsoEG24IorSingleIsoEG26[leptonCounter]=PassTriggerLeg("hltL1sSingleEG40IorSingleIsoEG22erIorSingleIsoEG24erIorSingleIsoEG24IorSingleIsoEG26",&*iE,iEvent);
	hltEle22WPLooseL1SingleIsoEG20erGsfTrackIsoFilter[leptonCounter]=PassTriggerLeg("hltEle22WPLooseL1SingleIsoEG20erGsfTrackIsoFilter",&*iE,iEvent);
	}

	

	//Flavour and charges
	_flavors[leptonCounter] = 0;
        _charges[leptonCounter] = iE->charge();
	_lchargePixSc[leptonCounter] =  iE->scPixCharge();
	_lchargeGSF[leptonCounter] = iE->isGsfScPixChargeConsistent() ? _lchargePixSc[leptonCounter]: - _lchargePixSc[leptonCounter] ;
	_lchargeCTF[leptonCounter] = iE->isGsfCtfChargeConsistent()? _lchargeGSF[leptonCounter] : -_lchargeGSF[leptonCounter];

	const reco::GsfTrackRef gsfTrack = iE->gsfTrack();
	_lnmisshits[leptonCounter] = gsfTrack->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS);
	
	bool vtxFitConversion = ConversionTools::hasMatchedConversion(reco::GsfElectron (*iE), theConversions, BS);
	_passconvveto[leptonCounter] = !vtxFitConversion;

	_lsietaieta[leptonCounter] =  iE->sigmaIetaIeta();
	_lfull5x5sietaieta[leptonCounter] =  iE->full5x5_sigmaIetaIeta();
	_lhovere[leptonCounter] =  iE->hadronicOverEm();
	_ldphiin[leptonCounter] =iE->deltaPhiSuperClusterTrackAtVtx();
	_ldetain[leptonCounter] =iE->deltaEtaSuperClusterTrackAtVtx();
	_l1oemin1op[leptonCounter] = TMath::Abs(1.0/iE->ecalEnergy() - iE->eSuperClusterOverP()/iE->ecalEnergy()) ;

        _isolation[leptonCounter] = pfRelIso(&*iE, myRhoJetsNC);
        _ipPV[leptonCounter] =iE->gsfTrack()->dxy(PV);
	_ipZPV[leptonCounter] = iE->gsfTrack()->dz(PV);
        

	//Iso
        _pfisocharged[leptonCounter]= iE->chargedHadronIso()/iE->pt();
        _pfisoneutral[leptonCounter] = iE->neutralHadronIso()/iE->pt();
        _pfisophoton[leptonCounter]= iE->photonIso()/iE->pt();
        
        _3dIP[leptonCounter]    = iE->dB(pat::Electron::PV3D);
        _3dIPerr[leptonCounter] = iE->edB(pat::Electron::PV3D);
        _3dIPsig[leptonCounter] = fabs(_3dIP[leptonCounter]/_3dIPerr[leptonCounter]);

        _mvaValue[leptonCounter] = (*mvaValues)[iE]; //myMVATrig->mvaValue(*iE,false);
          
          //Loose WP
          bool passed = false;
          if (TMath::Abs(iE->superCluster()->eta()) < 0.8 ) {
              passed = _mvaValue[leptonCounter] > looseMVA[0][0];
          }
          else if (TMath::Abs(iE->superCluster()->eta()) < 1.479 ) {
              passed = _mvaValue[leptonCounter] > looseMVA[1][0];
          }
          else {
              passed = _mvaValue[leptonCounter] > looseMVA[2][0];
          }
          /*
          if ((iE->pt())<10) {
              if (TMath::Abs(iE->superCluster()->eta()) < 0.8 ) {
                  passed = _mvaValue[leptonCounter] > looseMVA[0][0][0];
              }
              else if (TMath::Abs(iE->superCluster()->eta()) < 1.479 ) {
                  passed = _mvaValue[leptonCounter] > looseMVA[1][0][0];
              }
              else {
                  passed = _mvaValue[leptonCounter] > looseMVA[2][0][0];
              }
          }
          else if ((iE->pt())<20&&(iE->pt())>=10) {
              if (TMath::Abs(iE->superCluster()->eta()) < 0.8 ) {
                  passed = _mvaValue[leptonCounter] > looseMVA[0][0][1];
              }
              else if (TMath::Abs(iE->superCluster()->eta()) < 1.479 ) {
                  passed = _mvaValue[leptonCounter] > looseMVA[1][0][1];
              }
              else {
                  passed = _mvaValue[leptonCounter] > looseMVA[2][0][1];
              }
          }
          else if ((iE->pt())>=20) {
              if (TMath::Abs(iE->superCluster()->eta()) < 0.8 ) {
                  passed = _mvaValue[leptonCounter] > looseMVA[0][0][2];
              }
              else if (TMath::Abs(iE->superCluster()->eta()) < 1.479 ) {
                  passed = _mvaValue[leptonCounter] > looseMVA[1][0][2];
              }
              else {
                  passed = _mvaValue[leptonCounter] > looseMVA[2][0][2];
              }
          }*/
          _isloose[leptonCounter] =passed;
          
          passed = false;
          if (TMath::Abs(iE->superCluster()->eta()) < 0.8 ) {
              passed = _mvaValue[leptonCounter] > looseMVA[0][1];
          }
          else if (TMath::Abs(iE->superCluster()->eta()) < 1.479 ) {
              passed = _mvaValue[leptonCounter] > looseMVA[1][1];
          }
          else {
              passed = _mvaValue[leptonCounter] > looseMVA[2][1];
          }
          /*
          if ((iE->pt())<10) {
              if (TMath::Abs(iE->superCluster()->eta()) < 0.8 ) {
                  passed = _mvaValue[leptonCounter] > looseMVA[0][1][0];
              }
              else if (TMath::Abs(iE->superCluster()->eta()) < 1.479 ) {
                  passed = _mvaValue[leptonCounter] > looseMVA[1][1][0];
              }
              else {
                  passed = _mvaValue[leptonCounter] > looseMVA[2][1][0];
              }
          }
          else if ((iE->pt())<20&&(iE->pt())>=10) {
              if (TMath::Abs(iE->superCluster()->eta()) < 0.8 ) {
                  passed = _mvaValue[leptonCounter] > looseMVA[0][1][1];
              }
              else if (TMath::Abs(iE->superCluster()->eta()) < 1.479 ) {
                  passed = _mvaValue[leptonCounter] > looseMVA[1][1][1];
              }
              else {
                  passed = _mvaValue[leptonCounter] > looseMVA[2][1][1];
              }
          }
          else if ((iE->pt())>=20) {
              if (TMath::Abs(iE->superCluster()->eta()) < 0.8 ) {
                  passed = _mvaValue[leptonCounter] > looseMVA[0][1][2];
              }
              else if (TMath::Abs(iE->superCluster()->eta()) < 1.479 ) {
                  passed = _mvaValue[leptonCounter] > looseMVA[1][1][2];
              }
              else {
                  passed = _mvaValue[leptonCounter] > looseMVA[2][1][2];
              }
          }
          */
          _istightID[leptonCounter] = passed;
          
        
        ((TLorentzVector *)_leptonP4->At(leptonCounter))->SetPtEtaPhiE(iE->pt(), iE->eta(), iE->phi(), iE->energy());
        
        _mt[leptonCounter] = MT_calc(*((TLorentzVector *)_leptonP4->At(leptonCounter)), _met, _met_phi);
        
        _lPt[leptonCounter] = ((TLorentzVector *)_leptonP4->At(leptonCounter))->Pt();
        _lEta[leptonCounter] = ((TLorentzVector *)_leptonP4->At(leptonCounter))->Eta();
        _lPhi[leptonCounter] = ((TLorentzVector *)_leptonP4->At(leptonCounter))->Phi();
        _lE[leptonCounter] = ((TLorentzVector *)_leptonP4->At(leptonCounter))->E();
        
          
        _miniisolation_0p2[leptonCounter] = getPFIsolation(pfcands, dynamic_cast<const reco::Candidate *>(&*iE), 0.05, 0.2, 10., false, false, myRhoJetsNC);
        _miniisolation_0p3[leptonCounter] = getPFIsolation(pfcands, dynamic_cast<const reco::Candidate *>(&*iE), 0.05, 0.3, 10., false, false, myRhoJetsNC);
        _miniisolationcharged_0p2[leptonCounter] = getPFIsolation(pfcands, dynamic_cast<const reco::Candidate *>(&*iE), 0.05, 0.2, 10., false, true, myRhoJetsNC);
        _miniisoneutral[leptonCounter] = _miniisolation_0p2[leptonCounter]-_miniisolationcharged_0p2[leptonCounter];
        
        fillCloseJetVars(leptonCounter,PV);
        
        _isloose[leptonCounter] = ( _isloose[leptonCounter]&& _miniisolation_0p2[leptonCounter] < 0.4);

	_istightIso[leptonCounter] =TightIso("RA5",(&*iE),_mvaValue[leptonCounter],_miniisolation_0p2[leptonCounter],_ptratio[leptonCounter],_ptrel2[leptonCounter]);
	_istightID[leptonCounter] =TightID("RA5",&*iE,theConversions,BS,PV,_mvaValue[leptonCounter], _miniisolation_0p2[leptonCounter]);
    _istrigemulID[leptonCounter] =triggerEmulatorReturned(&*iE);
	_istrigemulISO[leptonCounter] =isoTriggerEmulator(&*iE,_mvaValue[leptonCounter]);

	_istightIDWP2016[leptonCounter] = _istightIso[leptonCounter]&&_istightID[leptonCounter];
	_istightIDWP2016_noIP3D[leptonCounter] =_istightIso[leptonCounter]&&TightID("RA5",&*iE,theConversions,BS,PV,_mvaValue[leptonCounter], _miniisolation_0p2[leptonCounter], false);
	_istightIDWP2016_IsoEmul[leptonCounter]=_istightIDWP2016[leptonCounter]  && _istrigemulISO[leptonCounter];
	_istightIDWP2016_RA7[leptonCounter] =TightIso("RA7",(&*iE),_mvaValue[leptonCounter],_miniisolation_0p2[leptonCounter],_ptratio[leptonCounter],_ptrel2[leptonCounter])&&TightID("RA7",&*iE,theConversions,BS,PV,_mvaValue[leptonCounter], _miniisolation_0p2[leptonCounter]);
    _istightIDWP2016_EWK[leptonCounter] =TightID_EWK("EWK",&*iE,theConversions,BS,PV,_mvaValue[leptonCounter],_lmva[leptonCounter],_miniisolation_0p2[leptonCounter]);
	
	_isFOIDWP2016[leptonCounter] =FOID("RA5",&*iE,theConversions,BS,PV,_mvaValue[leptonCounter],_miniisolation_0p2[leptonCounter]);
	_isFOIDWP2016_IsoEmul[leptonCounter] = _isFOIDWP2016[leptonCounter] && _istrigemulISO[leptonCounter];
	_isFOIDWP2016_noIP3D[leptonCounter] =FOID("RA5",&*iE,theConversions,BS,PV,_mvaValue[leptonCounter],_miniisolation_0p2[leptonCounter],false);
	_isFOIDWP2016_RA7[leptonCounter]=FOID("RA7",&*iE,theConversions,BS,PV,_mvaValue[leptonCounter],_miniisolation_0p2[leptonCounter]);
    _isFOIDWP2016_EWK[leptonCounter]=FOID_EWK("EWK",&*iE,theConversions,BS,PV,_mvaValue[leptonCounter],_lmva[leptonCounter],_ptratio[leptonCounter],_closeJetCSVAll[leptonCounter],_closeJetPtAll[leptonCounter],_miniisolation_0p2[leptonCounter]);
          
	_isVetoIDWP2016_RA7[leptonCounter]=VetoID(&*iE,theConversions, BS,PV, _mvaValue[leptonCounter],_miniisolation_0p2[leptonCounter]);
    _isVetoIDWP2016_EWK[leptonCounter]=VetoID_EWK(&*iE,theConversions, BS,PV, _mvaValue[leptonCounter],_miniisolation_0p2[leptonCounter]) ;
	_islooseMT2[leptonCounter]=EleCutBasedLooseID_MiniIso(&*iE,theConversions, BS,PV,_miniisolation_0p2[leptonCounter]) ;
	_iscutbasedWPtight[leptonCounter]=EleCutBasedTightID(&*iE,theConversions, BS,PV) ;
        if (IsMC) {
            //**************************************************************************************
            // MC
            //**************************************************************************************
	  MatchToGenParticle(iEvent,leptonCounter,iE->pt(),iE->eta(), iE->phi(), iE->pdgId(),genpartalreadymatched);
            
            fillIsoMCVars(leptonCounter);
            if (_closeIndex[leptonCounter] >=0)
                matchCloseJet(leptonCounter);
            
            //**************************************************************************************
            // MC *
            //**************************************************************************************
        }
        
        if (_regression)
            fillRegVars(SelectedJetsAll.at(_closeIndex[leptonCounter]), _closeJetPtAllMC[leptonCounter], &*iE);
        
        _sameSign[_istight[leptonCounter]][int((_charges[leptonCounter]+1)/2)]++;
        _mllZ[leptonCounter] = 9999.;
        _mllG[leptonCounter] = 9999.;

        leptonCounter++;
        
    }
    _nEle = leptonCounter - _nMu;



    /*
     *
     *Tau
     *
     *
     */

    for(unsigned int i = 0 ; i < sTau.size() ;i++ ){
        const pat::Tau *iT = sTau[i];
        _leptonIndex = i;
        
        if (leptonCounter == nLeptonsMax) continue;
	    InitLepton();

	    //Flavour and charges
	    _flavors[leptonCounter] = 2;
        _charges[leptonCounter] = iT->charge();
        _isolation[leptonCounter] = 0;
        _ipPV[leptonCounter] = iT->dxy();
        
        _tau_dz[leptonCounter]=sTau_dz[i];
        _3dIP[leptonCounter]    = sqrt(iT->dxy()*iT->dxy() + sTau_dz[i]*sTau_dz[i]);
        _3dIPerr[leptonCounter] = sqrt(iT->dxy_error()*iT->dxy_error() + iT->dzError()*iT->dzError());
        _3dIPsig[leptonCounter] = _3dIP[leptonCounter]/_3dIPerr[leptonCounter];
        
        ((TLorentzVector *)_leptonP4->At(leptonCounter))->SetPtEtaPhiE(iT->pt(), iT->eta(), iT->phi(), iT->energy());
    
        _mt[leptonCounter] = MT_calc(*((TLorentzVector *)_leptonP4->At(leptonCounter)), _met, _met_phi);
        _lPt[leptonCounter] = ((TLorentzVector *)_leptonP4->At(leptonCounter))->Pt();
        _lEta[leptonCounter] = ((TLorentzVector *)_leptonP4->At(leptonCounter))->Eta();
        _lPhi[leptonCounter] = ((TLorentzVector *)_leptonP4->At(leptonCounter))->Phi();
        _lE[leptonCounter] = ((TLorentzVector *)_leptonP4->At(leptonCounter))->E();
        

	_istightIso[leptonCounter] =TightIso(iT);
	_istightID[leptonCounter] =TightID("",iT,PV,false);
	_istightIDWP2016[leptonCounter] = _istightIso[leptonCounter]&&_istightID[leptonCounter];
	_istightIDWP2016_RA7[leptonCounter] =  _istightIso[leptonCounter]&&_istightID[leptonCounter];
    _istightIDWP2016_EWK[leptonCounter] =  _istightIso[leptonCounter]&&_istightID[leptonCounter];
	_isFOIDWP2016[leptonCounter] =FOID("",iT,PV,false);
	_isFOIDWP2016_RA7[leptonCounter]=_isFOIDWP2016[leptonCounter] ;
    _isFOIDWP2016_EWK[leptonCounter]=_isFOIDWP2016[leptonCounter] ;

	
        //_TauIs_againstElectronLooseMVA5[leptonCounter]=iT->tauID("againstElectronLooseMVA5");
        _TauIs_againstMuonLoose3[leptonCounter]=iT->tauID("againstMuonLoose3");
        _TauIs_byLooseCombinedIsolationDeltaBetaCorr3Hits[leptonCounter]=iT->tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits");
        _TauIs_decayModeFindingNewDMs[leptonCounter]=iT->tauID("decayModeFindingNewDMs");
        _TauIs_decayModeFinding[leptonCounter]=iT->tauID("decayModeFinding");
        _TauIs_byLooseIsolationMVArun2v1DBdR03oldDMwLT[leptonCounter]=iT->tauID("byLooseIsolationMVArun2v1DBdR03oldDMwLT");
        _TauIs_byMediumIsolationMVArun2v1DBdR03oldDMwLT[leptonCounter]=iT->tauID("byMediumIsolationMVArun2v1DBdR03oldDMwLT");
        _TauIs_byVLooseIsolationMVArun2v1DBoldDMwLT[leptonCounter]=iT->tauID("byVLooseIsolationMVArun2v1DBoldDMwLT");
        _TauIs_byTightIsolationMVArun2v1DBoldDMwLT[leptonCounter]=iT->tauID("byTightIsolationMVArun2v1DBoldDMwLT");
        
        if (IsMC) {
            //**************************************************************************************
            // MC
            //**************************************************************************************
	  MatchToGenParticle(iEvent,leptonCounter, iT->pt(),iT->eta(), iT->phi(), iT->pdgId(),genpartalreadymatched);
	  
	  fillIsoMCVars(leptonCounter);
	  if (_closeIndex[leptonCounter] >=0)
	    matchCloseJet(leptonCounter);
            
            //**************************************************************************************
            // MC *
            //**************************************************************************************
        }
        
	
        _sameSign[_istight[leptonCounter]][int((_charges[leptonCounter]+1)/2)]++;
        _mllZ[leptonCounter] = 9999.;
        _mllG[leptonCounter] = 9999.;

        leptonCounter++;
        


    }

    
    _nTau=leptonCounter-_nMu-_nEle;
    
    _nLeptons = leptonCounter;


    _n_Jets = 0; 
    _n_bJets = 0;
    _n_Jets40 = 0;
    HT = 0;
    TLorentzVector jt;
    for(unsigned int i = 0 ; i < SelectedJets.size() ;i++ ){
        if (fabs(SelectedJets[i]->eta()) > 2.4) continue;
	double uncPt = (SelectedJets[i]->correctedP4("Uncorrected")).Pt();
    double uncEta = (SelectedJets[i]->correctedP4("Uncorrected")).Eta();
	//double uncPhi = (SelectedJets[i]->correctedP4("Uncorrected")).Phi();

	double corr = fMetCorrector->getJetCorrectionRawPt(uncPt, uncEta, myRhoJets, SelectedJets[i]->jetArea(),_corrLevel);
	if (uncPt*corr < jetptcut) continue;

        if(_n_Jets>29)cout << "njets is " << _n_Jets<<endl;
	_jetJECuncty.push_back(fMetCorrector->getJECUncertainty(uncPt, uncEta) );
	_jetEta.push_back(SelectedJets[i]->eta());
	_jetPhi.push_back(SelectedJets[i]->phi());
	_jetPt.push_back(uncPt*corr);
	_jetM.push_back(SelectedJets[i]->mass());
	_jetpFlav.push_back(SelectedJets[i]->partonFlavour());
	_jethFlav.push_back(SelectedJets[i]->hadronFlavour());

	//See here:https://twiki.cern.ch/twiki/bin/view/CMS/BTagSFMethods 
	//I use method 1.a
	//For max flexibility I save the SF and eff per jet. 
	//This allows me to recompute the event weight from the ntuples
	//	cout << "jet pt eta " << _jetEta[i]<<","<< _jetPt[i]<<","<<_jethFlav[i]<<endl;
	//cout << reader->eval(BTagEntry::FLAV_B, _jetEta[i], _jetPt[i])<<endl; 
	//cout << btagEffs2D[0]->GetBinContent(btagEffs2D[0]->FindBin(_jetPt[i],fabs(_jetEta[i])))<< endl;
	double ptoverflow = _jetPt[i]< btagEffs2D[0]->GetXaxis()->GetBinLowEdge(btagEffs2D[0]->GetNbinsX()+1)? _jetPt[i]: btagEffs2D[0]->GetXaxis()->GetBinLowEdge(btagEffs2D[0]->GetNbinsX()+1) -0.1;
	//cout << "ptoverflow: " << ptoverflow << endl;
	//cout << btagEffs2D[2]->GetEntries()<<endl;

	if ( _jethFlav[i] ==5) {
	  _jetbtagSF.push_back(reader->eval_auto_bounds("central",BTagEntry::FLAV_B, _jetEta[i], ptoverflow));
	  _jetbtagSF.push_back(reader_up->eval_auto_bounds("up", BTagEntry::FLAV_B, _jetEta[i], ptoverflow));
	  _jetbtagSF_down.push_back(reader_do->eval_auto_bounds("down", BTagEntry::FLAV_B, _jetEta[i], ptoverflow));
	  _jetbtagEff.push_back(btagEffs2D[0]->GetBinContent(btagEffs2D[0]->FindBin(ptoverflow,fabs(_jetEta[i]))));
	}
	else if ( _jethFlav[i] ==4) {
	  _jetbtagSF.push_back(reader->eval_auto_bounds("central",BTagEntry::FLAV_C, _jetEta[i], ptoverflow));
	  _jetbtagSF_up.push_back(reader_up->eval_auto_bounds("up", BTagEntry::FLAV_C, _jetEta[i], ptoverflow));
	  _jetbtagSF_down.push_back(reader_do->eval_auto_bounds("down",BTagEntry::FLAV_C, _jetEta[i], ptoverflow));
	  _jetbtagEff.push_back(btagEffs2D[1]->GetBinContent(btagEffs2D[1]->FindBin(ptoverflow,fabs(_jetEta[i]))));
	}
	else {
	  _jetbtagSF.push_back(1);//(reader->eval_auto_bounds(BTagEntry::FLAV_UDSG, _jetEta[i], ptoverflow));
	  _jetbtagSF_up.push_back(1);//(reader_up->eval_auto_bounds(BTagEntry::FLAV_UDSG, _jetEta[i], ptoverflow));
	  _jetbtagSF_down.push_back(1);//(reader_do->eval_auto_bounds(BTagEntry::FLAV_UDSG, _jetEta[i], ptoverflow));
	  _jetbtagEff.push_back(btagEffs2D[2]->GetBinContent(btagEffs2D[2]->FindBin(ptoverflow,fabs(_jetEta[i]))));
	}

        
        jt.SetPtEtaPhiM(_jetPt[_n_Jets],_jetEta[_n_Jets],_jetPhi[_n_Jets],SelectedJets[i]->mass());
        
	//cout << "debug before btag " << endl ;
	_csv.push_back(SelectedJets[i]->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
	                                                
        //To be replaced by:                           "pfCombinedInclusiveSecondaryVertexV2BJetTags"
	//cout <<"debug after btag " <<endl ;

        if(_csv[_n_Jets] > btag_cut) {
	  //            _bTagged[_n_Jets] = true;
	  _bTagged.push_back(true); //[_n_Jets]
            _n_bJets++;
        } else _bTagged.push_back(false);
	//else _bTagged[_n_Jets] = false;
        
        if (_jetPt[_n_Jets] > 40) {
            HT+= _jetPt[_n_Jets];
            _n_Jets40++;
        }
	//        _jetDeltaRloose[_n_Jets] = 99999;
	_jetDeltaRloose.push_back(99999);
        _n_Jets++;
    }
    
    
    //_jetDeltaRloose and Mll
    //Muons
    for( std::vector<pat::Muon>::const_iterator mu = thePatMuons->begin() ; mu != thePatMuons->end() ; mu++ ) {
        
      if ( mu->pt()  < _muonMinPt) continue;
        if ( TMath::Abs( mu->eta() ) > 2.4 ) continue;
        
        if ( !(mu->isGlobalMuon() || mu->isTrackerMuon() )) continue;
        if ( !(mu->isPFMuon()) ) continue;
        
        const reco::TrackRef innerTrack = mu->innerTrack();
        if( innerTrack.isNull() ) continue;
        
        if(TMath::Abs(innerTrack->dxy(PV)) > _looseD0Mu  ) continue;
        if(TMath::Abs(innerTrack->dz(PV))   > 0.1  ) continue;
        
        double iso = getPFIsolation(pfcands, dynamic_cast<const reco::Candidate *>(&*mu), 0.05, 0.2, 10., false, false, myRhoJetsNC);
        if (iso < 0.4) {
            TLorentzVector mu1; mu1.SetPtEtaPhiE(mu->pt(),mu->eta(),mu->phi(),mu->energy());
            
            if (mu->pt() > _muonMinPt) {
            //if (mu->pt() > _muonMinPt) {
                for(int i = 0 ; i != _n_Jets ;i++ ){
                    jt.SetPtEtaPhiM(_jetPt[i],_jetEta[i],_jetPhi[i],_jetM[i]);
                    double dR = mu1.DeltaR( jt );
                    //std::cout<<dR<<std::endl;
                    if( dR < _jetDeltaRloose[i] ) {
                        _jetDeltaRloose[i] = dR;
                    }
                }
            }
            for(int i = 0 ; i != _nMu ;i++ ) {
                if (_charges[i] != mu->charge()) {
                    double mdiL = Mll_calc( *((TLorentzVector*)_leptonP4->At(i)), mu1);
                    if (mdiL <  _mllG[i]) _mllG[i] = mdiL;
                    if (fabs(mdiL - 91) < fabs(_mllZ[i] - 91)) _mllZ[i] = mdiL;
                }
            }
        }
    }
    
    //Ele
    for (size_t i = 0; i < thePatElectronsView->size(); ++i){
        const auto el = thePatElectronsView->ptrAt(i);
        
        if ( el->pt()  <  7.) continue;
        if ( TMath::Abs( el->superCluster()->eta() ) > 2.5 ) continue;
        
        const reco::GsfTrackRef gsfTrack = el->gsfTrack();
        
        if (!gsfTrack.isNonnull()) {
            continue;
        }
        if( TMath::Abs(gsfTrack->dxy(PV)) > 0.05/*_looseD0E*/  )  continue;
        if( TMath::Abs(gsfTrack->dz(PV)) > 0.1  ) continue;
        
        bool vtxFitConversion = ConversionTools::hasMatchedConversion(reco::GsfElectron (*el), theConversions, BS);
        if( vtxFitConversion )  continue;
        if( gsfTrack->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) > 1 ) continue;
        
        double mvaValueE = (*mvaValues)[el];
        //std::cout<<"MVA value "<<mvaValueE<<"; eta "<< iE->eta()<<std::endl;
        
        bool passed = false;
        if (TMath::Abs(el->superCluster()->eta()) < 0.8 ) {
            passed = mvaValueE > looseMVA[0][0];
        }
        else if (TMath::Abs(el->superCluster()->eta()) < 1.479 ) {
            passed = mvaValueE > looseMVA[1][0];
        }
        else {
            passed = mvaValueE > looseMVA[2][0];
        }
        /*
        if ((el->pt())<10) {
            if (TMath::Abs(el->superCluster()->eta()) < 0.8 ) {
                passed = mvaValueE > looseMVA[0][0][0];
            }
            else if (TMath::Abs(el->superCluster()->eta()) < 1.479 ) {
                passed = mvaValueE > looseMVA[1][0][0];
            }
            else {
                passed = mvaValueE > looseMVA[2][0][0];
            }
        }
        else if ((el->pt())<20&&(el->pt())>=10) {
            if (TMath::Abs(el->superCluster()->eta()) < 0.8 ) {
                passed = mvaValueE > looseMVA[0][0][1];
            }
            else if (TMath::Abs(el->superCluster()->eta()) < 1.479 ) {
                passed = mvaValueE > looseMVA[1][0][1];
            }
            else {
                passed = mvaValueE > looseMVA[2][0][1];
            }
        }
        else if ((el->pt())>=20) {
            if (TMath::Abs(el->superCluster()->eta()) < 0.8 ) {
                passed = mvaValueE > looseMVA[0][0][2];
            }
            else if (TMath::Abs(el->superCluster()->eta()) < 1.479 ) {
                passed = mvaValueE > looseMVA[1][0][2];
            }
            else {
                passed = mvaValueE > looseMVA[2][0][2];
            }
        }
        */
        if (!passed) continue;
        
        double iso = getPFIsolation(pfcands, dynamic_cast<const reco::Candidate *>(&*el), 0.05, 0.2, 10., false, false, myRhoJetsNC);
        if (iso < 0.4) {
            TLorentzVector mu1; mu1.SetPtEtaPhiE(el->pt(),el->eta(),el->phi(),el->energy());
            if (el->pt() > _muonMinPt) {
                //if (el->pt() > _muonMinPt) {
                for(int i = 0 ; i != _n_Jets ;i++ ){
                    jt.SetPtEtaPhiM(_jetPt[i],_jetEta[i],_jetPhi[i],_jetM[i]);
                    double dR = mu1.DeltaR( jt );
                    //std::cout<<dR<<std::endl;
                    if( dR < _jetDeltaRloose[i] ) {
                        _jetDeltaRloose[i] = dR;
                    }
                }
            }
            for(int i = _nMu ; i != _nLeptons ;i++ ) {
                if (_charges[i] != el->charge()) {
                    double mdiL = Mll_calc( *((TLorentzVector*)_leptonP4->At(i)), mu1);
                    if (mdiL <  _mllG[i]) _mllG[i] = mdiL;
                    if (fabs(mdiL - 91) < fabs(_mllZ[i] - 91)) _mllZ[i] = mdiL;
                }
            }
        }
        
    }//Ele
    
    for (int i=0; i!=_nLeptons; ++i) {
      vector < double> mymasses; 
      for (int j=0; j!=_nLeptons; ++j) {
        double themass =  Mll_calc( *((TLorentzVector*)_leptonP4->At(i)), *((TLorentzVector*)_leptonP4->At(j)));
	mymasses.push_back(themass); 
	  }
      _mllv.push_back(mymasses);
      mymasses.clear();
    }
    if(debug) cout << "Debug 10 " << _nLeptons<< endl;
    
  
    
    
    /*
     *
     *Apply skims
     *
     */
    
    if(Skim == "Zprime"){
      TLorentzVector lep1,lep2; 
      double highestmass(0);
      for( unsigned int i  =0; i< _lPt.size();i++){
	if(_lPt[i]<20)continue;
	for( unsigned int j  =i+1; j< _lPt.size();j++){
	  if(_lPt[j]<20)continue;
	  lep1.SetPtEtaPhiM(_lPt[i],_lEta[i],_lPhi[i],0.000511);
	  lep2.SetPtEtaPhiM(_lPt[j],_lEta[j],_lPhi[j],0.000511);
	  if((lep1+lep2).Mag()>highestmass)highestmass= (lep1+lep2).Mag();
	}
      }
      if(highestmass<1000) return;
      //cout << "High Mass: " << highestmass << endl;
    }


    if(Skim == "SkimRA7"){
      int nfoid_ra7(0);
      if(_met<30) return;
      for( unsigned int i  =0; i< _lPt.size();i++){
        if(_isFOIDWP2016_RA7[i] && _lPt[i]>10) nfoid_ra7++;
      }
      if(nfoid_ra7<3) return;
    }

    if(Skim == "SignalSkimRA5"){
      if(_met<30) return;
      int ntightl (0), njpt35(0);
      for( unsigned int i  =0; i< _lPt.size();i++){
	if(_istightIDWP2016[i] && _lPt[i]>10) ntightl++; 
      }
      if(ntightl<2) return;

      for( unsigned int i  =0; i< _jetPt.size();i++){
	bool isjetclean = true;
	if(_jetPt[i] <35|| abs(_jetEta[i])>2.4) continue;
	for( unsigned int j=0; j< _lPt.size();j++){
	  if(_lPt[j]<10 || !_isFOIDWP2016[j]) continue;
	  double deta = _jetEta[i]-_lEta[j];
	  double dphi = fabs(acos(cos(_jetPhi[i]-_lPhi[j])));
	  if(sqrt(deta*deta+dphi*dphi)<0.1) isjetclean =false;
	}
	
	if(isjetclean) njpt35++;
      }
      if(njpt35<2)return;
    }
    

    if(Skim == "FakeRateAppSkimRA5"){
        if(_met<30) return;
        int nFOl (0), njpt35(0);
        int ntightl(0);
        for( unsigned int i  =0; i< _lPt.size();i++){
            if(_isFOIDWP2016[i] && _lPt[i] >10 ) nFOl++;
            if(_istightIDWP2016[i] && _lPt[i]>10) ntightl++;
        }
        if(nFOl<2|| ntightl == nFOl) return;
        
        for( unsigned int i  =0; i< _jetPt.size();i++){
            bool isjetclean = true;
            if(_jetPt[i] <35|| abs(_jetEta[i])>2.4) continue;
            for( unsigned int j=0; j< _lPt.size();j++){
                if(_lPt[j]<10 || !_isFOIDWP2016[j]) continue;
                double deta = _jetEta[i]-_lEta[j];
                double dphi = fabs(acos(cos(_jetPhi[i]-_lPhi[j])));
                if(sqrt(deta*deta+dphi*dphi)<0.1) isjetclean =false;
            }
            
            if(isjetclean) njpt35++;
        }
        
        if(njpt35<2)return;
    }

    if(Skim == "InSituSkimRA5"){
      //      if(_met<30) return;
      int nFOl (0), njpt35(0);
      int ntightl(0);
      for( unsigned int i  =0; i< _lPt.size();i++){
	if(_isFOIDWP2016_noIP3D[i] && _lPt[i] >10 ) nFOl++; 
	if(_istightIDWP2016[i] && _lPt[i]>10) ntightl++;
      }
      if(nFOl<2|| ntightl == 0) return;

      for( unsigned int i  =0; i< _jetPt.size();i++){
	bool isjetclean = true;
	if(_jetPt[i] <35|| abs(_jetEta[i])>2.4) continue;
	for( unsigned int j=0; j< _lPt.size();j++){
	  if(_lPt[j]<10 || !_isFOIDWP2016[j]) continue;
	  double deta = _jetEta[i]-_lEta[j];
	  double dphi = fabs(acos(cos(_jetPhi[i]-_lPhi[j])));
	  if(sqrt(deta*deta+dphi*dphi)<0.1) isjetclean =false;
	}
	
	if(isjetclean) njpt35++;
      }
      // if(njpt35<2)return;
    }


    if(Skim == "skim2LrecoandTriggerTP"){

      bool highptl= false;
      int nbofgoodl(0);
      for(unsigned int i  =0; i< _miniisolation_0p2.size();i++){
	if( (_isFOIDWP2016_RA7[i] && _lPt[i]>5) || (_flavors[i] ==1 && _issoftmuonwithiso[i] && _lPt[i]<30) )
	nbofgoodl ++;
	if(_istightIDWP2016[i] && _lPt[i]>10) highptl=true ;
      }
      if (nbofgoodl <2 || !highptl) return;
    }
    if(Skim == "2TriggEles" ){
      int nbofgoodeles(0);
      for(unsigned int i  =0; i< _miniisolation_0p2.size();i++){
	if(_miniisolation_0p2[i]<0.2&&_flavors[i]==0&& _istrigemulID[i] ) nbofgoodeles++;
	
      }
      if(nbofgoodeles<2)return;
    }
    
 

    outputTree->Fill();   
}

void SSb13_takeimai::fillMCVars(const GenParticle* mc, const int leptonCounter) {
    
    _lPtmc[leptonCounter] = mc->pt();
    _lEmc[leptonCounter] = mc->energy();
    _lPhimc[leptonCounter] = mc->phi();
    _lEtamc[leptonCounter] = mc->eta();
    
    _chargesMC[leptonCounter] = mc->charge();
    
    _origin[leptonCounter] = GPM.origin(mc);
    _originReduced[leptonCounter] = GPM.originReduced(_origin[leptonCounter]);
    const GenParticle* mcMom = GPM.getMotherParton(mc);
    if (mcMom!=NULL) {
        
        //std::cout<<"Mother: "<<std::endl;
        //std::cout<<mcMom->pdgId()<<" "<<mcMom->pt()<<std::endl;
        //std::cout<<mc->pt()<<std::endl;
        
        _mompt[leptonCounter] = mcMom->pt();
        _momphi[leptonCounter] = mcMom->phi();
        _mometa[leptonCounter] = mcMom->eta();
        _mompdg[leptonCounter] = mcMom->pdgId();
        
        PVmc = mcMom->vertex();
        
        //std::cout<<"d0: "<<_ipPVmc<<" "<<_ipPV<<std::endl;
        
        TLorentzVector Gen0;
        Gen0.SetPtEtaPhiE( 0, 0, 0, 0);
        
        
        _isolationMCraw[leptonCounter] = 0;
        _isolationMCnonu[leptonCounter] = 0;
        std::vector<GenParticleRef> dauts;
        for (unsigned int dau=0; dau!=mcMom->numberOfDaughters(); ++dau) {
            GenParticleRef daut = mcMom->daughterRef(dau);
            dauts.push_back(daut);
        }
        unsigned int counterD = 0;
        //if (_isolation == 0)
        //    std::cout<<"==== "<<chargedHadronIso*iM->pt()<<" "<<neutralHadronIso*iM->pt()<<" "<<photonIso*iM->pt()<<" "<<beta*iM->pt()<<std::endl;
        while (counterD < dauts.size()) {
            if (dauts.at(counterD)->status() == 1) {
                _isolationMCraw[leptonCounter]+=dauts.at(counterD)->pt();
                if ((fabs(dauts.at(counterD)->pdgId())!= 12) && (fabs(dauts.at(counterD)->pdgId())!= 14) && (fabs(dauts.at(counterD)->pdgId())!= 16)) {
                    _isolationMCnonu[leptonCounter]+=dauts.at(counterD)->pt();
                    /*if (_isolation == 0) {
                     TLorentzVector dauV;
                     dauV.SetPtEtaPhiE(dauts.at(counterD)->pt(), dauts.at(counterD)->eta(), dauts.at(counterD)->phi(), dauts.at(counterD)->energy());
                     double deltaR1 = ((TLorentzVector *)_leptonP4->At(leptonCounter))->DeltaR(dauV);
                     std::cout<<counterD<<" "<<dauts.at(counterD)->pdgId()<<" "<<dauts.at(counterD)->pt()<<
                     " "<<dauts.at(counterD)->eta()<<" "<<dauts.at(counterD)->phi()<<" "<<dauts.at(counterD)->energy()<<" in "<<deltaR1<<std::endl;
                     }*/
                } else {
                    TLorentzVector Gen;
                    Gen.SetPtEtaPhiE( dauts.at(counterD)->pt(), dauts.at(counterD)->eta(), dauts.at(counterD)->phi(), dauts.at(counterD)->energy() );
                    Gen0 += Gen;
                }
            } else {
                for (unsigned int dau=0; dau!=dauts.at(counterD)->numberOfDaughters(); ++dau) {
                    GenParticleRef daut = dauts.at(counterD)->daughterRef(dau);
                    dauts.push_back(daut);
                }
            }
            counterD++;
        }
        
        //std::cout<<"PTs: "<<mc->pt()<<" "<<iM->pt()<<std::endl;
        _isolationMCraw[leptonCounter]/=mc->pt();
        _isolationMCraw[leptonCounter]-=1;
        _isolationMCnonu[leptonCounter]/=mc->pt();
        _isolationMCnonu[leptonCounter]-=1;
        
        _nuPtmc[leptonCounter] = Gen0.Pt();
        _nuPhimc[leptonCounter] = Gen0.Phi();
        _nuEtamc[leptonCounter] = Gen0.Eta();
        _nuEmc[leptonCounter] = Gen0.E();
        
        TLorentzVector lmc;
        lmc.SetPtEtaPhiE(_lPtmc[leptonCounter],_lEtamc[leptonCounter],_lPhimc[leptonCounter],_lEmc[leptonCounter]);
        _mtmc[leptonCounter] = MT_calc(lmc, _nuPtmc[leptonCounter], _nuPhimc[leptonCounter]);
    } else {
        _mompt[leptonCounter] = 0;
        _momphi[leptonCounter] = 0;
        _mometa[leptonCounter] = 0;
        _mompdg[leptonCounter] = 0;
    }    
}
void SSb13_takeimai::fillCloseJetVars(const int leptonCounter,Vertex::Point PV) {
    _closeJetPtAll[leptonCounter] = ((TLorentzVector *)_leptonP4->At(leptonCounter))->Pt();
    _closeJetAngAll[leptonCounter] = 10000;
    _ptRelAll[leptonCounter] = 0;
    _ptrel[leptonCounter] = 0;
    _ptrel2[leptonCounter] = 0;
    _ptratio[leptonCounter] = 1.;
    _ljetmultiplicity[leptonCounter] = 0;
    _closeIndex[leptonCounter] = 0;
    TLorentzVector pJet;
    //   cout <<  endl;
    //cout << "Lepton pt, eta, phi: " << ((TLorentzVector *)_leptonP4->At(leptonCounter))->Pt()<<", "<<  ((TLorentzVector *)_leptonP4->At(leptonCounter))->Eta()<<", "<< ((TLorentzVector *)_leptonP4->At(leptonCounter))->Phi()<<endl;
    for(unsigned int k = 0 ; k < SelectedJetsAll.size() ;k++ ){
      

      double uncorrPt = (SelectedJetsAll[k]->correctedP4("Uncorrected")).Pt();
      double corr = fMetCorrector->getJetCorrectionRawPt( uncorrPt, (SelectedJetsAll[k]->correctedP4("Uncorrected")).Eta(), myRhoJets, SelectedJetsAll[k]->jetArea(),"L1FastJet");
      double corr2 = fMetCorrector->getJetCorrectionRawPt( uncorrPt, (SelectedJetsAll[k]->correctedP4("Uncorrected")).Eta(), myRhoJets, SelectedJetsAll[k]->jetArea(),_corrLevel);
      pJet.SetPtEtaPhiE( corr*uncorrPt, SelectedJetsAll[k]->eta(), SelectedJetsAll[k]->phi(), corr*(SelectedJetsAll[k]->correctedP4("Uncorrected")).E());
      pJet-=*((TLorentzVector *)_leptonP4->At(leptonCounter));
      pJet*=corr2/corr;
      pJet+=*((TLorentzVector *)_leptonP4->At(leptonCounter));
      double ang = ((TLorentzVector *)_leptonP4->At(leptonCounter))->DeltaR( pJet );

      TLorentzVector lp4 = *((TLorentzVector *)_leptonP4->At(leptonCounter));


      //Trk multiplicity
      int trkmultiplicity = 0;
      const std::vector<reco::CandidatePtr> & cands = (&*SelectedJetsAll[k])->daughterPtrVector();
      for ( std::vector<reco::CandidatePtr>::const_iterator cand = cands.begin();   cand != cands.end(); ++cand ) {

	const pat::PackedCandidate *pfcand = dynamic_cast<const pat::PackedCandidate *>(cand->get());
	if(pfcand->charge() ==0)continue;
	if(pfcand->fromPV() <= 1)continue;
	if(pfcand->pseudoTrack().hitPattern().numberOfValidHits()<8)continue;
	if(pfcand->pseudoTrack().hitPattern().numberOfValidPixelHits()<2)continue;
	if(pfcand->pseudoTrack().normalizedChi2()>=5)continue;
	if(std::fabs(pfcand->pseudoTrack().dxy(PV))>=0.2)continue;
	if(std::fabs(pfcand->pseudoTrack().dz(PV))>=17)continue;
	TLorentzVector candp4; 

	candp4.SetPtEtaPhiE(pfcand->pt(),pfcand->eta(),pfcand->phi(),pfcand->energy());
	double dRjet_cand = candp4.DeltaR( pJet );
	if(dRjet_cand>0.4)continue;
	trkmultiplicity++;
      }
      
      
      //      cout << "Jet pt, eta, phi (uncorr): " <<  (SelectedJetsAll[k]->correctedP4("Uncorrected")).Pt() <<", "<< (SelectedJetsAll[k]->correctedP4("Uncorrected")).Eta()<<", "
      //   <<(SelectedJetsAll[k]->correctedP4("Uncorrected")).Phi()<<endl;
      //double rhoforjetl = (pJet-(TLorentzVector *)_leptonP4->At(leptonCounter)).Rho();
      // cout << "(jet(l1l2l3 corr)-lepton).Rho() " << (pJet-lp4).Rho()<<endl;
	

      if (ang < _closeJetAngAll[leptonCounter]) {
	_closeJetAngAll[leptonCounter] = ang;
	_closeJetPtAll[leptonCounter] = pJet.Pt();
	_closeJetEtaAll[leptonCounter] = SelectedJetsAll[k]->eta();
	_closeJetPhiAll[leptonCounter] = SelectedJetsAll[k]->phi();
	_closeJetEAll[leptonCounter] = pJet.E();
	//	_closeJetMAll[leptonCounter] = pJet.M();
	_closeJetNconstAll[leptonCounter] = SelectedJetsAll[k]->numberOfDaughters();
	_closeJetCSVAll[leptonCounter] = SelectedJetsAll[k]->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
          
          //Calculate _trackSelectionMultiplicity
          int trackSelectionMult = 0;
          
          for(int unsigned it=0; it!=SelectedJetsAll[k]->numberOfDaughters(); ++it) {
              const pat::PackedCandidate * icand = dynamic_cast<const pat::PackedCandidate *> (SelectedJetsAll[k]->daughter(it));
              const reco::Track& trk = icand->pseudoTrack();
              Double_t deta = trk.eta() - pJet.Eta();
              Double_t dphi = TVector2::Phi_mpi_pi(trk.phi() - pJet.Phi());
              bool trackSelection =  trk.pt()>1 && trk.charge() != 0 && (TMath::Sqrt( deta*deta+dphi*dphi ) < 0.4) && (icand->fromPV() > 1) && (trk.hitPattern().numberOfValidHits()>=8) && (trk.hitPattern().numberOfValidPixelHits()>=2) && (trk.normalizedChi2()<5) && std::fabs(trk.dxy(PV))<0.2 && std::fabs(trk.dz(PV))<17;
              if(trackSelection) trackSelectionMult++;
          }
          
          _trackSelectionMultiplicity[leptonCounter] = trackSelectionMult;
          

	_ptRelAll[leptonCounter] = ((TLorentzVector *)_leptonP4->At(leptonCounter))->Perp(pJet.Vect());
	_ptrel[leptonCounter] = ((TLorentzVector *)_leptonP4->At(leptonCounter))->Perp(pJet.Vect() - ((TLorentzVector *)_leptonP4->At(leptonCounter))->Vect());
	pJet-=*((TLorentzVector *)_leptonP4->At(leptonCounter));
	_ptrel2[leptonCounter] = ((TLorentzVector *)_leptonP4->At(leptonCounter))->Perp(pJet.Vect());
	
	_ptratio[leptonCounter] = ((TLorentzVector *)_leptonP4->At(leptonCounter))->Pt()/(_closeJetPtAll[leptonCounter]);
	//std::cout<<_ptrel[leptonCounter]<<std::endl;                                                                                        
	_closeIndex[leptonCounter] = k;
	_ljetmultiplicity[leptonCounter]=trkmultiplicity;

        }
    }
    
    if (_closeJetAngAll[leptonCounter] > 0.4) {
        _closeJetPtAll[leptonCounter] = ((TLorentzVector *)_leptonP4->At(leptonCounter))->Pt();
        _closeJetEAll[leptonCounter] = ((TLorentzVector *)_leptonP4->At(leptonCounter))->E();
        _closeJetPhiAll[leptonCounter] = ((TLorentzVector *)_leptonP4->At(leptonCounter))->Phi();
        _closeJetEtaAll[leptonCounter] = ((TLorentzVector *)_leptonP4->At(leptonCounter))->Eta();
        _closeJetNconstAll[leptonCounter] = 1;
        _closeJetCSVAll[leptonCounter] = 0;
        _closeJetAngAll[leptonCounter] = 0;
        _ptRelAll[leptonCounter] = 0;
        _ptrel[leptonCounter] = 0;
	_ptrel2[leptonCounter] = 0;

        _ptratio[leptonCounter] = 1.;
        _closeIndex[leptonCounter] = -1;
	_ljetmultiplicity[leptonCounter] =0;
    }
    
    //std::cout<<_closeIndex[leptonCounter]<<" "<<_closeJetPtAll[leptonCounter]<<" "<<_closeJetAngAll[leptonCounter]<<" "<<
    //_ptrel[leptonCounter]<<" "<<_ptratio[leptonCounter]<<std::endl;

    if (_miniisolation_0p2[leptonCounter] < _multiConst[0][0] && (_ptratio[leptonCounter] > _multiConst[0][1] || _ptrel2[leptonCounter] > _multiConst[0][2]) ) _multiisolation_L[leptonCounter] = true;     
    if (_miniisolation_0p2[leptonCounter] < _multiConst[1][0] && (_ptratio[leptonCounter] > _multiConst[1][1] || _ptrel2[leptonCounter] > _multiConst[1][2]) ) _multiisolation_M[leptonCounter] = true;     
    if (_miniisolation_0p2[leptonCounter] < _multiConst[2][0] && (_ptratio[leptonCounter] > _multiConst[2][1] || _ptrel2[leptonCounter] > _multiConst[2][2]) ) _multiisolation_T[leptonCounter] = true;
    
    /*
     *LeptonMVA with ISO
     */
    LepGood_pt_forMVA=_lPt[leptonCounter];
    LepGood_eta_forMVA=_lEta[leptonCounter];
    LepGood_jetNDauChargedMVASel_forMVA=_ljetmultiplicity[leptonCounter];
    LepGood_miniRelIsoCharged_forMVA= (_miniisolationcharged_0p2[leptonCounter] <0 )?0: _miniisolationcharged_0p2[leptonCounter];
    LepGood_miniRelIsoNeutral_forMVA= (_miniisolation_0p2[leptonCounter]-LepGood_miniRelIsoCharged_forMVA <0 )?0: _miniisolation_0p2[leptonCounter]-LepGood_miniRelIsoCharged_forMVA;
    LepGood_jetPtRelv2_forMVA=_ptrel2[leptonCounter];
    LepGood_jetPtRatio_forMVA=TMath::Min(_ptratio[leptonCounter],1.5);
    LepGood_jetBTagCSV_forMVA=TMath::Max(_closeJetCSVAll[leptonCounter],0.);
    LepGood_sip3d_forMVA=_3dIPsig[leptonCounter];
    LepGood_dxy_forMVA=TMath::Log(fabs(_ipPV[leptonCounter]));
    LepGood_dz_forMVA=TMath::Log(fabs(_ipZPV[leptonCounter]));
    LepGood_segmentCompatibility_forMVA= _muonSegmentComp[leptonCounter];
    LepGood_mvaIdSpring16GP_forMVA=_mvaValue[leptonCounter];

    float BDToutput= -99.;
    if(_flavors[leptonCounter] ==0&&usemva )BDToutput= (float)tmvareaderele->EvaluateMVA( "BDTG method"           );
    if(_flavors[leptonCounter] ==1&&usemva )BDToutput= (float)tmvareadermu->EvaluateMVA( "BDTG method"           );

    _lmvawithiso[leptonCounter] =BDToutput;
    
    /*
     *LeptonMVA
     */
    LepGood_pt_forMVA=_lPt[leptonCounter];
    LepGood_eta_forMVA=_lEta[leptonCounter];
    LepGood_jetNDauChargedMVASel_forMVA=_trackSelectionMultiplicity[leptonCounter];
    LepGood_miniRelIsoCharged_forMVA= (_miniisolationcharged_0p2[leptonCounter] <0 )?0: _miniisolationcharged_0p2[leptonCounter];
    LepGood_miniRelIsoNeutral_forMVA= (_miniisolation_0p2[leptonCounter]-LepGood_miniRelIsoCharged_forMVA <0 )?0: _miniisolation_0p2[leptonCounter]-LepGood_miniRelIsoCharged_forMVA;
    LepGood_jetPtRelv2_forMVA=_ptrel2[leptonCounter];
    LepGood_jetPtRatio_forMVA=TMath::Min(_ptratio[leptonCounter],1.5);
    LepGood_jetBTagCSV_forMVA=TMath::Max(_closeJetCSVAll[leptonCounter],0.);
    LepGood_sip3d_forMVA=_3dIPsig[leptonCounter];
    LepGood_dxy_forMVA=TMath::Log(fabs(_ipPV[leptonCounter]));
    LepGood_dz_forMVA=TMath::Log(fabs(_ipZPV[leptonCounter]));
    LepGood_segmentCompatibility_forMVA= _muonSegmentComp[leptonCounter];
    LepGood_mvaIdSpring16GP_forMVA=_mvaValue[leptonCounter];
    
    if(_flavors[leptonCounter] ==0&&usemva )BDToutput= (float)tmvareaderele->EvaluateMVA( "BDTG method"           );
    if(_flavors[leptonCounter] ==1&&usemva )BDToutput= (float)tmvareadermu->EvaluateMVA( "BDTG method"           );
    
    _lmva[leptonCounter] =BDToutput;

    
}

void SSb13_takeimai::matchCloseJet(const int leptonCounter) {
    double minDeltaR3 = 9999;
    double minDeltaR2 = 9999;
    TLorentzVector Gen1, Gen2;
    const GenParticle* mc3a = 0;
    const GenParticle* mc2a = 0;
    Gen1.SetPtEtaPhiE(SelectedJetsAll.at(_closeIndex[leptonCounter])->pt(),SelectedJetsAll.at(_closeIndex[leptonCounter])->eta(),SelectedJetsAll.at(_closeIndex[leptonCounter])->phi(),SelectedJetsAll.at(_closeIndex[leptonCounter])->energy());
    
    for(GenParticleCollection::const_reverse_iterator p = TheGenParticles->rbegin() ; p != TheGenParticles->rend() ; p++ ) {
        int id = TMath::Abs(p->pdgId());
        
        if ((id > 0 && id < 6) || (id == 21) || (id == 22)) {
            if (p->status() != 2) {
                if (fabs(p->eta()) <= 10 )
                    Gen2.SetPtEtaPhiE(p->pt(),p->eta(),p->phi(),p->energy());
                else continue;
                //Gen2.SetPtEtaPhiM(p->pt(),0.00001,p->phi(),0);
                double deltaRcur = Gen1.DeltaR(Gen2);
                if (deltaRcur < minDeltaR3) {
                    mc3a = &*p;
                    minDeltaR3 = deltaRcur;
                }
                
            } else if (p->status() == 2) {
                if (fabs(p->eta()) <= 10 )
                    Gen2.SetPtEtaPhiE(p->pt(),p->eta(),p->phi(),p->energy());
                else continue;
                //Gen2.SetPtEtaPhiM(p->pt(),0.00001,p->phi(),0);
                double deltaRcur = Gen1.DeltaR(Gen2);
                if (deltaRcur < minDeltaR2) {
                    mc2a = &*p;
                    minDeltaR2 = deltaRcur;
                }
            }
        }
    }
    
    if ((minDeltaR3 < 0.5 || minDeltaR2 == 9999) && mc3a!=0) {
        _closeJetPtAllMC[leptonCounter] = mc3a->pt();
        _closeJetPtAllstatus[leptonCounter] = 3;
        _partonIdMatched[leptonCounter] = mc3a->pdgId();
    } else if (mc2a!=0) {
        _closeJetPtAllMC[leptonCounter] = mc2a->pt();
        _closeJetPtAllstatus[leptonCounter] = 2;
        _partonIdMatched[leptonCounter] = mc2a->pdgId();
    } else {
        _closeJetPtAllMC[leptonCounter] = 0;
        _closeJetPtAllstatus[leptonCounter] = 0;
        _partonIdMatched[leptonCounter] = 0;
    }
}
void SSb13_takeimai::fillIsoMCVars(const int leptonCounter) {
    if( TheGenParticles.isValid() )
    {
        //if (_isolation[0] == 0) {
        //    std::cout<<"==="<<std::endl;
        //}
        _isolationMCdr03[leptonCounter] = 0;
        _isolationMCdr03nonu[leptonCounter] = 0;
        for(GenParticleCollection::const_reverse_iterator p = TheGenParticles->rbegin() ; p != TheGenParticles->rend() ; p++ ) {
            //int id = TMath::Abs(p->pdgId());
            if (p->status() == 1) {
                TLorentzVector pmc; pmc.SetPtEtaPhiM( p->pt(), p->eta(), p->phi(), p->mass() );
                double ang = ((TLorentzVector *)_leptonP4->At(leptonCounter))->DeltaR( pmc );
                if (ang < 0.3) {
                    _isolationMCdr03[leptonCounter]+=p->pt();
                    if (fabs(p->pdgId())!= 12 && fabs(p->pdgId())!= 14 && fabs(p->pdgId())!= 16) {
                        _isolationMCdr03nonu[leptonCounter]+=p->pt();
                        //if (_isolation[0] == 0) {
                        //    std::cout<<" "<<p->pdgId()<<" "<<p->pt()<<
                        //    " "<<p->eta()<<" "<<p->phi()<<" "<<p->energy()<<" in "<<ang<<std::endl;
                        //}
                    }
                }
            }
        }
        double leptpt = ((TLorentzVector *)_leptonP4->At(leptonCounter))->Pt();
        _isolationMCdr03[leptonCounter]/=leptpt;
        _isolationMCdr03nonu[leptonCounter]/=leptpt;
        _isolationMCdr03[leptonCounter]-=1;
        _isolationMCdr03nonu[leptonCounter]-=1;
        //std::cout<<"*************"<<std::endl;
        //for (int i=0; i!=4; ++i)
        //    std::cout<<_isolationMC[i]<<" ";
        //std::cout<<_isolation[0]<<" "<<_closeJetPtAllMC[0]/leptpt-1<<" "<<_closeJetPtAll[0]/leptpt-1<<std::endl;
    }
}


void SSb13_takeimai::fillRegVars(const pat::Jet *jet, double genpt, const pat::Muon* mu) {
    
    hJet_ptRaw = (jet->correctedP4("Uncorrected")).Pt();
    hJet_genPt = genpt;
    hJet_pt = jet->pt();
    hJet_phi = jet->phi();
    hJet_eta = jet->eta();
    hJet_e = jet->energy();
    
    hJet_ptLeadTrack = 0;
    
    const reco::TrackRefVector &tracks =  jet->associatedTracks();
    for (unsigned int k=0; k!= tracks.size(); ++k) {
        if(tracks[k]->pt() > hJet_ptLeadTrack)
            hJet_ptLeadTrack = tracks[k]->pt();
    }
    
    hJet_vtx3dL = 0;
    hJet_vtx3deL = 0;
    hJet_vtxMass = 0;
    hJet_vtxPt = 0;
    
    const reco::SecondaryVertexTagInfo* scdVtx = jet->tagInfoSecondaryVertex("secondaryVertex");
    
    if (scdVtx) {
        //std::cout<<"Vertetx info: "<<scdVtx->nVertices()<<std::endl;
        if (scdVtx->nVertices()) {
            const reco::Vertex &sv1 = scdVtx->secondaryVertex(0);
            if (!sv1.isFake()) {
                Measurement1D distance1 = scdVtx->flightDistance(0, true);
                hJet_vtx3dL = distance1.value();
                hJet_vtx3deL = distance1.error();
                
                math::XYZTLorentzVectorD p4vtx = sv1.p4();
                hJet_vtxMass = p4vtx.M();
                hJet_vtxPt = p4vtx.Pt();
            }
        }
    }
    
    
    hJet_cef = jet->chargedHadronEnergyFraction() + jet->chargedEmEnergyFraction();
    hJet_nconstituents = jet->getPFConstituents().size();
    hJet_JECUnc = fMetCorrector->getJECUncertainty(jet->pt(),jet->eta());
    
    
    TLorentzVector pJet; pJet.SetPtEtaPhiE( jet->pt(), jet->eta(), jet->phi(), jet->energy() );
    TLorentzVector pLep; pLep.SetPtEtaPhiE( mu->pt(), mu->eta(), mu->phi(), mu->energy() );
    
    hJet_SoftLeptptRel = pLep.Perp(pJet.Vect());
    hJet_SoftLeptPt = mu->pt();
    hJet_SoftLeptdR = pLep.DeltaR(pJet);
    
    hJet_SoftLeptIdlooseMu = 1;
    hJet_SoftLeptId95 = 1;
}

void SSb13_takeimai::fillRegVars(const pat::Jet *jet, double genpt, const pat::Electron* mu) {
    
    hJet_ptRaw = (jet->correctedP4("Uncorrected")).Pt();
    hJet_genPt = genpt;
    hJet_pt = jet->pt();
    hJet_phi = jet->phi();
    hJet_eta = jet->eta();
    hJet_e = jet->energy();
    
    hJet_ptLeadTrack = 0;
    
    const reco::TrackRefVector &tracks =  jet->associatedTracks();
    for (unsigned int k=0; k!= tracks.size(); ++k) {
        if(tracks[k]->pt() > hJet_ptLeadTrack)
            hJet_ptLeadTrack = tracks[k]->pt();
    }
    
    hJet_vtx3dL = 0;
    hJet_vtx3deL = 0;
    hJet_vtxMass = 0;
    hJet_vtxPt = 0;
    
    const reco::SecondaryVertexTagInfo* scdVtx = jet->tagInfoSecondaryVertex("secondaryVertex");
    
    if (scdVtx) {
        //std::cout<<"Vertetx info: "<<scdVtx->nVertices()<<std::endl;
        if (scdVtx->nVertices()) {
            const reco::Vertex &sv1 = scdVtx->secondaryVertex(0);
            if (!sv1.isFake()) {
                Measurement1D distance1 = scdVtx->flightDistance(0, true);
                hJet_vtx3dL = distance1.value();
                hJet_vtx3deL = distance1.error();
                
                math::XYZTLorentzVectorD p4vtx = sv1.p4();
                hJet_vtxMass = p4vtx.M();
                hJet_vtxPt = p4vtx.Pt();
            }
        }
    }
    
    
    hJet_cef = jet->chargedHadronEnergyFraction() + jet->chargedEmEnergyFraction();
    hJet_nconstituents = jet->getPFConstituents().size();
    hJet_JECUnc = fMetCorrector->getJECUncertainty(jet->pt(),jet->eta());
    
    
    TLorentzVector pJet; pJet.SetPtEtaPhiE( jet->pt(), jet->eta(), jet->phi(), jet->energy() );
    TLorentzVector pLep; pLep.SetPtEtaPhiE( mu->pt(), mu->eta(), mu->phi(), mu->energy() );
    
    hJet_SoftLeptptRel = pLep.Perp(pJet.Vect());
    hJet_SoftLeptPt = mu->pt();
    hJet_SoftLeptdR = pLep.DeltaR(pJet);
    
    hJet_SoftLeptIdlooseMu = 1;
    hJet_SoftLeptId95 = 1;
}

void SSb13_takeimai::bookTree() {
  bool istriggerstudy =   (Skim=="skim2LrecoandTriggerTP"); 

  bool frmeasurement =  (Skim=="skim1Lreco");

    //outputTree->Branch("_leptonP4", "TClonesArray", &_leptonP4, 32000, 0);
    //outputTree->Branch("_jetP4", "TClonesArray", &_jetP4, 32000, 0);
    outputTree->Branch("_eventNb",   &_eventNb,   "_eventNb/l");
    outputTree->Branch("_runNb",     &_runNb,     "_runNb/l");
    outputTree->Branch("_lumiBlock", &_lumiBlock, "_lumiBlock/l");
    
    if(!istriggerstudy&&!frmeasurement) outputTree->Branch("passmetfilters",&passmetfilters,"passmetfilters/O");
    /*
    if(!istriggerstudy&&!frmeasurement) outputTree->Branch("passcsctighthalo2015", &passcsctighthalo2015, "passcsctighthalo2015/O");
    if(!istriggerstudy&&!frmeasurement) outputTree->Branch("passglobaltighthalo2016", &passglobaltighthalo2016, "passglobaltighthalo2016/O");
    if(!istriggerstudy&&!frmeasurement) outputTree->Branch("passglobalsupertighthalo2016", &passglobalsupertighthalo2016, "passglobalsupertighthalo2016/O");
     */

    /*
    if(!istriggerstudy&&!frmeasurement) outputTree->Branch("passHLT_PFMET170_NotCleaned",&passHLT_PFMET170_NotCleaned,"passHLT_PFMET170_NotCleaned/O");
    if(!istriggerstudy&&!frmeasurement) outputTree->Branch("passHLT_PFMET170_HBHECleaned",&passHLT_PFMET170_HBHECleaned,"passHLT_PFMET170_HBHECleaned/O");
    if(!istriggerstudy&&!frmeasurement) outputTree->Branch("passHLT_PFMET170_BeamHaloCleaned",&passHLT_PFMET170_BeamHaloCleaned,"passHLT_PFMET170_BeamHaloCleaned/O");
    if(!frmeasurement)outputTree->Branch("passHLT_PFMET120_Mu5",&passHLT_PFMET120_Mu5,"passHLT_PFMET120_Mu5/O");
     */
    //if(!istriggerstudy&&!frmeasurement) outputTree->Branch("passHLT_Ele15_IsoVVVL_PFHT350_PFMET50",&passHLT_Ele15_IsoVVVL_PFHT350_PFMET50,"passHLT_Ele15_IsoVVVL_PFHT350_PFMET50/O");
    //if(!istriggerstudy&&!frmeasurement) outputTree->Branch("passHLT_IsoMu16_eta2p1_MET30",&passHLT_IsoMu16_eta2p1_MET30,"passHLT_IsoMu16_eta2p1_MET30/O");
    if(!istriggerstudy&&!frmeasurement) outputTree->Branch("passHLT_IsoMu19_eta2p1_LooseIsoPFTau20_v",&passHLT_IsoMu19_eta2p1_LooseIsoPFTau20_v,"passHLT_IsoMu19_eta2p1_LooseIsoPFTau20_v/O");
    if(!istriggerstudy&&!frmeasurement) outputTree->Branch("passHLT_Ele22_eta2p1_WPLoose_GSF_LooseIsoPFtau20_SingleL1_v",&passHLT_Ele22_eta2p1_WPLoose_GSF_LooseIsoPFtau20_SingleL1_v,"passHLT_Ele22_eta2p1_WPLoose_GSF_LooseIsoPFtau20_SingleL1_v/O");
    /*
    if(!istriggerstudy&&!frmeasurement) outputTree->Branch("passHLT_Mu17_Mu8_SameSign_DZ",&passHLT_Mu17_Mu8_SameSign_DZ,"passHLT_Mu17_Mu8_SameSign_DZ/O");
    if(!frmeasurement)outputTree->Branch("passHLT_DoubleMu3_PFMET50",&passHLT_DoubleMu3_PFMET50,"passHLT_DoubleMu3_PFMET50/O");
    if(!frmeasurement)outputTree->Branch("passHLT_TripleMu_5_3_3",&passHLT_TripleMu_5_3_3,"passHLT_TripleMu_5_3_3/O");
    if(!istriggerstudy) outputTree->Branch("passHLT_Ele8_CaloIdM_TrackIdM_PFJet30",&passHLT_Ele8_CaloIdM_TrackIdM_PFJet30,"passHLT_Ele8_CaloIdM_TrackIdM_PFJet30/O");
     */
    if(!istriggerstudy) outputTree->Branch("passHLT_Ele12_CaloIdM_TrackIdM_PFJet30",&passHLT_Ele12_CaloIdM_TrackIdM_PFJet30,"passHLT_Ele12_CaloIdM_TrackIdM_PFJet30/O");
    if(!istriggerstudy) outputTree->Branch("passHLT_Ele17_CaloIdM_TrackIdM_PFJet30",&passHLT_Ele17_CaloIdM_TrackIdM_PFJet30,"passHLT_Ele17_CaloIdM_TrackIdM_PFJet30/O");
    if(!istriggerstudy) outputTree->Branch("passHLT_Ele23_CaloIdM_TrackIdM_PFJet30",&passHLT_Ele23_CaloIdM_TrackIdM_PFJet30,"passHLT_Ele23_CaloIdM_TrackIdM_PFJet30/O");
    if(!istriggerstudy) outputTree->Branch("passHLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30",&passHLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30,"passHLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30/O");
    if(!istriggerstudy) outputTree->Branch("passHLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30",&passHLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30,"passHLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30/O");
    if(!istriggerstudy) outputTree->Branch("passHLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30",&passHLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30,"passHLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30/O");
    if(!istriggerstudy) outputTree->Branch("passHLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30",&passHLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30,"passHLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30/O");
    /*
    if(!istriggerstudy) outputTree->Branch("passHLT_Ele15_IsoVVVL_PFHT350",&passHLT_Ele15_IsoVVVL_PFHT350,"passHLT_Ele15_IsoVVVL_PFHT350/O");
    if(!istriggerstudy) outputTree->Branch("passHLT_Ele10_CaloIdM_TrackIdM_CentralPFJet30_BTagCSV_p13",&passHLT_Ele10_CaloIdM_TrackIdM_CentralPFJet30_BTagCSV_p13,"passHLT_Ele10_CaloIdM_TrackIdM_CentralPFJet30_BTagCSV_p13/O");
    */
    /*
    if(!istriggerstudy) outputTree->Branch("passHLT_Mu10_CentralPFJet30_BTagCSV_p13",&passHLT_Mu10_CentralPFJet30_BTagCSV_p13,"passHLT_Mu10_CentralPFJet30_BTagCSV_p13/O");
     */
    if(!istriggerstudy) outputTree->Branch("passHLT_Mu8",&passHLT_Mu8,"passHLT_Mu8/O");
    if(!istriggerstudy) outputTree->Branch("passHLT_Mu17",&passHLT_Mu17,"passHLT_Mu17/O");
    if(!istriggerstudy) outputTree->Branch("passHLT_Mu20",&passHLT_Mu20,"passHLT_Mu20/O");
    if(!istriggerstudy) outputTree->Branch("passHLT_Mu27",&passHLT_Mu27,"passHLT_Mu27/O");
    if(!istriggerstudy) outputTree->Branch("passHLT_Mu8_TrkIsoVVL",&passHLT_Mu8_TrkIsoVVL,"passHLT_Mu8_TrkIsoVVL/O");
    if(!istriggerstudy) outputTree->Branch("passHLT_Mu17_TrkIsoVVL",&passHLT_Mu17_TrkIsoVVL,"passHLT_Mu17_TrkIsoVVL/O");
    //if(!istriggerstudy) outputTree->Branch("passHLT_Mu15_IsoVVVL_PFHT350",&passHLT_Mu15_IsoVVVL_PFHT350,"passHLT_Mu15_IsoVVVL_PFHT350/O");


    
    //if(!frmeasurement)outputTree->Branch("passHLT_DoubleMu8_Mass8_PFHT300",&passHLT_DoubleMu8_Mass8_PFHT300  ,"passHLT_DoubleMu8_Mass8_PFHT300/O");
    //if(!frmeasurement)outputTree->Branch("passHLT_DoubleMu8_Mass8_PFHT250",&passHLT_DoubleMu8_Mass8_PFHT250  ,"passHLT_DoubleMu8_Mass8_PFHT250/O");
    if(!frmeasurement)outputTree->Branch("passHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL",&passHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL ,"passHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL/O");
    if(!frmeasurement)outputTree->Branch("passHLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL",&passHLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL ,"passHLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL/O");
    if(!frmeasurement)outputTree->Branch("passHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ",&passHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ ,"passHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ/O");
    if(!frmeasurement)outputTree->Branch("passHLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ",&passHLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ ,"passHLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ/O");
    //if(!frmeasurement)outputTree->Branch("passHLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT300",&passHLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT300 ,"passHLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT300/O");
    //if(!frmeasurement)outputTree->Branch("passHLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT250",&passHLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT250 ,"passHLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT250/O");
    if(!frmeasurement)outputTree->Branch("passHLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL",&passHLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL ,"passHLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL/O");
    if(!frmeasurement)outputTree->Branch("passHLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",&passHLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ ,"passHLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ/O");
    //if(!frmeasurement)outputTree->Branch("passHLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL",&passHLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL ,"passHLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL/O");
    //if(!frmeasurement)outputTree->Branch("passHLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",&passHLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ ,"passHLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ/O");

    //if(!frmeasurement)outputTree->Branch("passHLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT300",&passHLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT300 ,"passHLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT300/O");
    //if(!frmeasurement)outputTree->Branch("passHLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT250",&passHLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT250 ,"passHLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT250/O");
    if(!frmeasurement)outputTree->Branch("passHLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL",&passHLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL ,"passHLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL/O");
    if(!frmeasurement)outputTree->Branch("passHLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL",&passHLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL ,"passHLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL/O");
    if(!frmeasurement)outputTree->Branch("passHLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL",&passHLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL ,"passHLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL/O");
    //if(!frmeasurement)outputTree->Branch("passHLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL",&passHLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL ,"passHLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL/O");
    //if(!frmeasurement)outputTree->Branch("passHLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL",&passHLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL ,"passHLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL/O");


    /*
    if(!frmeasurement)outputTree->Branch("passHLT_IsoMu24_eta2p1",&passHLT_IsoMu24_eta2p1 ,"passHLT_IsoMu24_eta2p1/O");
    if(!frmeasurement)outputTree->Branch("passHLT_IsoMu20",&passHLT_IsoMu20 ,"passHLT_IsoMu20/O");
    if(!frmeasurement)outputTree->Branch("passHLT_IsoTkMu20",&passHLT_IsoTkMu20 ,"passHLT_IsoTkMu20/O");
    if(!frmeasurement)outputTree->Branch("passHLT_IsoMu18",&passHLT_IsoMu18 ,"passHLT_IsoMu18/O");
    if(!frmeasurement)outputTree->Branch("passHLT_IsoTkMu18",&passHLT_IsoTkMu18 ,"passHLT_IsoTkMu18/O");
     */
    if(!frmeasurement)outputTree->Branch("passHLT_IsoMu22",&passHLT_IsoMu22 ,"passHLT_IsoMu22/O");
    if(!frmeasurement)outputTree->Branch("passHLT_IsoTkMu22",&passHLT_IsoTkMu22 ,"passHLT_IsoTkMu22/O");


    //if(!frmeasurement)outputTree->Branch("passHLT_Ele23_WPLoose_Gsf",&passHLT_Ele23_WPLoose_Gsf ,"passHLT_Ele23_WPLoose_Gsf/O");
    //if(!frmeasurement)outputTree->Branch("passHLT_Ele27_WPLoose_Gsf",&passHLT_Ele27_WPLoose_Gsf ,"passHLT_Ele27_WPLoose_Gsf/O");
    if(!frmeasurement)outputTree->Branch("passHLT_Ele27_eta2p1_WPLoose_Gsf",&passHLT_Ele27_eta2p1_WPLoose_Gsf ,"passHLT_Ele27_eta2p1_WPLoose_Gsf/O");
     /*
    if(!frmeasurement)outputTree->Branch("passHLT_Ele32_eta2p1_WPLoose_Gsf",&passHLT_Ele32_eta2p1_WPLoose_Gsf ,"passHLT_Ele32_eta2p1_WPLoose_Gsf/O");
    if(!frmeasurement)outputTree->Branch("passHLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW",&passHLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW ,"passHLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW/O");
    if(!frmeasurement)outputTree->Branch("passHLT_Ele23_CaloIdL_TrackIdL_IsoVL",&passHLT_Ele23_CaloIdL_TrackIdL_IsoVL,"passHLT_Ele23_CaloIdL_TrackIdL_IsoVL/O");


    if(!frmeasurement)outputTree->Branch("passHLT_TripleMu_12_10_5",&passHLT_TripleMu_12_10_5,"passHLT_TripleMu_12_10_5/O");
    if(!frmeasurement)outputTree->Branch("passHLT_DiMu9_Ele9_CaloIdL_TrackIdL",&passHLT_DiMu9_Ele9_CaloIdL_TrackIdL,"passHLT_DiMu9_Ele9_CaloIdL_TrackIdL/O");
    if(!frmeasurement)outputTree->Branch("passHLT_Mu8_DiEle12_CaloIdL_TrackIdL",&passHLT_Mu8_DiEle12_CaloIdL_TrackIdL,"passHLT_Mu8_DiEle12_CaloIdL_TrackIdL/O");
    if(!frmeasurement)outputTree->Branch("passHLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL",&passHLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL,"passHLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL/O");
     */
    
    //Prescale factors
    outputTree->Branch("Prescale_passHLT_Mu8", &Prescale_passHLT_Mu8, "Prescale_passHLT_Mu8/D");
    outputTree->Branch("Prescale_passHLT_Mu17", &Prescale_passHLT_Mu17, "Prescale_passHLT_Mu17/D");
    outputTree->Branch("Prescale_passHLT_Mu8_TrkIsoVVL", &Prescale_passHLT_Mu8_TrkIsoVVL, "Prescale_passHLT_Mu8_TrkIsoVVL/D");
    outputTree->Branch("Prescale_passHLT_Mu17_TrkIsoVVL", &Prescale_passHLT_Mu17_TrkIsoVVL, "Prescale_passHLT_Mu17_TrkIsoVVL/D");
    outputTree->Branch("Prescale_passHLT_Ele12_CaloIdM_TrackIdM_PFJet30", &Prescale_passHLT_Ele12_CaloIdM_TrackIdM_PFJet30, "Prescale_passHLT_Ele12_CaloIdM_TrackIdM_PFJet30/D");
    outputTree->Branch("Prescale_passHLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30", &Prescale_passHLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30, "Prescale_passHLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30/D");
    outputTree->Branch("Prescale_passHLT_Ele17_CaloIdM_TrackIdM_PFJet30", &Prescale_passHLT_Ele17_CaloIdM_TrackIdM_PFJet30, "Prescale_passHLT_Ele17_CaloIdM_TrackIdM_PFJet30/D");
    outputTree->Branch("Prescale_passHLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30", &Prescale_passHLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30, "Prescale_passHLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30/D");


    outputTree->Branch("_weight", &_weight, "_weight/D");
    if(!istriggerstudy&&!frmeasurement) outputTree->Branch("_genHT", &_genHT, "_genHT/D");
    
    outputTree->Branch("_nLeptons", &_nLeptons, "_nLeptons/I");
    //if(!frmeasurement)outputTree->Branch("nbofL1EG", &nbofL1EG, "nbofL1EG/I");
    
    
    /*
    if(!istriggerstudy&&!frmeasurement) outputTree->Branch("_index1", &_index1, "_index1/I");
    if(!istriggerstudy&&!frmeasurement) outputTree->Branch("_index2", &_index2, "_index2/I");

    if(!istriggerstudy&&!frmeasurement) outputTree->Branch("_sb", &_sb, "_sb/O");
    if(!istriggerstudy&&!frmeasurement) outputTree->Branch("_doubleF", &_doubleF, "_doubleF/O");
     */
    
    //outputTree->Branch("_PVchi2", &_PVchi2, "_PVchi2/D");
    //outputTree->Branch("_PVerr", &_PVerr, "_PVerr[3]/D");
    
    if (_regression) {
        outputTree->Branch("_regVars", &_regVars, "_regVars[15]/D");
        
        outputTree->Branch("hJet_ptRaw", &hJet_ptRaw, "hJet_ptRaw/D");
        outputTree->Branch("hJet_genPt", &hJet_genPt, "hJet_genPt/D");
        outputTree->Branch("hJet_pt", &hJet_pt, "hJet_pt/D");
        outputTree->Branch("hJet_phi", &hJet_phi, "hJet_phi/D");
        outputTree->Branch("hJet_eta", &hJet_eta, "hJet_eta/D");
        outputTree->Branch("hJet_e", &hJet_e, "hJet_e/D");
        
        outputTree->Branch("hJet_ptLeadTrack", &hJet_ptLeadTrack, "hJet_ptLeadTrack/D");
        
        outputTree->Branch("hJet_vtx3dL", &hJet_vtx3dL, "hJet_vtx3dL/D");
        outputTree->Branch("hJet_vtx3deL", &hJet_vtx3deL, "hJet_vtx3deL/D");
        outputTree->Branch("hJet_vtxMass", &hJet_vtxMass, "hJet_vtxMass/D");
        outputTree->Branch("hJet_vtxPt", &hJet_vtxPt, "hJet_vtxPt/D");
        
        outputTree->Branch("hJet_cef", &hJet_cef, "hJet_cef/D");
        
        outputTree->Branch("hJet_nconstituents", &hJet_nconstituents, "hJet_nconstituents/D");
        outputTree->Branch("hJet_JECUnc", &hJet_JECUnc, "hJet_JECUnc/D");
        
        outputTree->Branch("hJet_SoftLeptptRel", &hJet_SoftLeptptRel, "hJet_SoftLeptptRel/D");
        outputTree->Branch("hJet_SoftLeptPt", &hJet_SoftLeptPt, "hJet_SoftLeptPt/D");
        outputTree->Branch("hJet_SoftLeptdR", &hJet_SoftLeptdR, "hJet_SoftLeptdR/D");
        
        outputTree->Branch("hJet_SoftLeptIdlooseMu", &hJet_SoftLeptIdlooseMu, "hJet_SoftLeptIdlooseMu/D");
        outputTree->Branch("hJet_SoftLeptId95", &hJet_SoftLeptId95, "hJet_SoftLeptId95/D");
    }
    
    outputTree->Branch("_n_PV", &_n_PV, "_n_PV/I");
    outputTree->Branch("trueNVtx", &trueNVtx,"trueNVtx/I");
    //if(!istriggerstudy&&!frmeasurement) outputTree->Branch("nVtxNow", &nVtxNow,"nVtxNow/I");
    //if(!istriggerstudy&&!frmeasurement) outputTree->Branch("nVtxBefore", &nVtxBefore,"nVtxBefore/I");
    //if(!istriggerstudy&&!frmeasurement) outputTree->Branch("nVtxAfter", &nVtxAfter,"nVtxAfter/I");
    outputTree->Branch("_met", &_met, "_met/D");
    outputTree->Branch("_met_phi", &_met_phi, "_met_phi/D");
    /*
    if(!istriggerstudy&&!frmeasurement) outputTree->Branch("_met_JECup", &_met_JECup, "_met_JECup/D");
    if(!istriggerstudy&&!frmeasurement) outputTree->Branch("_met_JECdown", &_met_JECdown, "_met_JECdown/D");
   
    if(!istriggerstudy&&!frmeasurement) outputTree->Branch("_met_phi_JECup", &_met_phi_JECup, "_met_phi_JECup/D");
    if(!istriggerstudy&&!frmeasurement) outputTree->Branch("_met_phi_JECdown", &_met_phi_JECdown, "_met_phi_JECdown/D");
     */

    if(!istriggerstudy&&!frmeasurement) outputTree->Branch("HT", &HT, "HT/D");
    
    /*
    if(!istriggerstudy&&!frmeasurement) outputTree->Branch("_genmet", &_genmet, "_genmet/D");
    if(!istriggerstudy&&!frmeasurement) outputTree->Branch("_genmet_phi", &_genmet_phi, "_genmet_phi/D");
    if(!istriggerstudy&&!frmeasurement) outputTree->Branch("_mgluino", &_mgluino, "_mgluino/D");
    if(!istriggerstudy&&!frmeasurement) outputTree->Branch("_mLSP", &_mLSP, "_mLSP/D");
    if(!istriggerstudy&&!frmeasurement) outputTree->Branch("_mchargino", &_mchargino, "_mchargino/D");
    if(!istriggerstudy&&!frmeasurement) outputTree->Branch("_mstop", &_mstop, "_mstop/D");
     */

    outputTree->Branch("_n_bJets", &_n_bJets, "_n_bJets/I");
    
    outputTree->Branch("_n_Jets", &_n_Jets, "_n_Jets/I");
    outputTree->Branch("_n_Jets40", &_n_Jets40, "_n_Jets40/I");
    
    //if(!istriggerstudy&&!frmeasurement) outputTree->Branch("_closeIndex",&_closeIndex);
    //if(!istriggerstudy&&!frmeasurement) outputTree->Branch("_indeces",&_indeces);
    outputTree->Branch("_flavors",&_flavors);
    outputTree->Branch("_charges",&_charges);
    //if(!istriggerstudy&&!frmeasurement) outputTree->Branch("_chargesMC",&_chargesMC);
    //if(!istriggerstudy)outputTree->Branch("_isolation",&_isolation);
    outputTree->Branch("_miniisolation_0p2",&_miniisolation_0p2);
    if(!frmeasurement)outputTree->Branch("_miniisolationcharged_0p2",&_miniisolationcharged_0p2);
    if(!frmeasurement)outputTree->Branch("_miniisoneutral",&_miniisoneutral);
    if(!frmeasurement)outputTree->Branch("_miniisolation_0p3",&_miniisolation_0p3);
    /*
    if(!istriggerstudy&&!frmeasurement) outputTree->Branch("_multiisolation_T",&_multiisolation_T);
    if(!istriggerstudy&&!frmeasurement) outputTree->Branch("_multiisolation_M",&_multiisolation_M);
    if(!istriggerstudy&&!frmeasurement) outputTree->Branch("_multiisolation_L",&_multiisolation_L);
     */
    if(!frmeasurement)outputTree->Branch("_muonSegmentComp",&_muonSegmentComp);
    if(!frmeasurement)outputTree->Branch("_trackSelectionMultiplicity",&_trackSelectionMultiplicity);
    if(!frmeasurement)outputTree->Branch("_lmvawithiso",&_lmvawithiso);
    if(!frmeasurement)outputTree->Branch("_lmva",&_lmva);
    /*
    if(!frmeasurement)outputTree->Branch("_pfisocharged",&_pfisocharged);
    if(!frmeasurement)outputTree->Branch("_pfisoneutral",&_pfisoneutral);
    if(!frmeasurement)outputTree->Branch("_pfisophoton",&_pfisophoton);
    if(!istriggerstudy&&!frmeasurement) outputTree->Branch("_isolationMCraw",&_isolationMCraw);
    if(!istriggerstudy&&!frmeasurement) outputTree->Branch("_isolationMCnonu",&_isolationMCnonu);
    if(!istriggerstudy&&!frmeasurement) outputTree->Branch("_isolationMCdr03",&_isolationMCdr03);
    if(!istriggerstudy&&!frmeasurement) outputTree->Branch("_isolationMCdr03nonu",&_isolationMCdr03nonu);
     */
    outputTree->Branch("_ptrel",&_ptrel);
    if(!istriggerstudy)outputTree->Branch("_ptrel2",&_ptrel2);
    outputTree->Branch("_ptratio",&_ptratio);
    //if(!frmeasurement)outputTree->Branch("_ljetmultiplicity",&_ljetmultiplicity);
    outputTree->Branch("_mvaValue",&_mvaValue);
    /*
    if(!istriggerstudy&&!frmeasurement) outputTree->Branch("_mt",&_mt);
    if(!istriggerstudy&&!frmeasurement) outputTree->Branch("_mllZ",&_mllZ);
    if(!istriggerstudy&&!frmeasurement) outputTree->Branch("_mllG",&_mllG);
    if(!istriggerstudy&&!frmeasurement) outputTree->Branch("_mll",&_mll);
    if(!istriggerstudy&&!frmeasurement) outputTree->Branch("_mllv",&_mllv);
    if(!istriggerstudy) outputTree->Branch("_origin",&_origin);
    if(!istriggerstudy&&!frmeasurement) outputTree->Branch("_originReduced",&_originReduced);
     */
    outputTree->Branch("_ipPV",&_ipPV);
    //outputTree->Branch("_ipPVerr",&_ipPVerr);
    //if(!istriggerstudy)outputTree->Branch("_ipPVmc",&_ipPVmc);
    outputTree->Branch("_ipZPV",&_ipZPV);
    //outputTree->Branch("_ipZPVerr",&_ipZPVerr);
    outputTree->Branch("_3dIP",&_3dIP);
    //outputTree->Branch("_3dIPerr",&_3dIPerr);
    outputTree->Branch("_3dIPsig",&_3dIPsig);
    if(!istriggerstudy) outputTree->Branch("_closeJetPtAll",&_closeJetPtAll);
    /*
    if(!istriggerstudy) outputTree->Branch("_closeJetEtaAll",&_closeJetEtaAll);
    if(!istriggerstudy) outputTree->Branch("_closeJetPhiAll",&_closeJetPhiAll);
    if(!istriggerstudy&&!frmeasurement) outputTree->Branch("_closeJetEAll",&_closeJetEAll);
    */
    if(!istriggerstudy&&!frmeasurement) outputTree->Branch("_closeJetCSVAll",&_closeJetCSVAll);
    /*
    if(!istriggerstudy&&!frmeasurement) outputTree->Branch("_closeJetNconstAll",&_closeJetNconstAll);
     */
    if(!istriggerstudy&&!frmeasurement) outputTree->Branch("_closeJetAngAll",&_closeJetAngAll);
    if(!istriggerstudy&&!frmeasurement) outputTree->Branch("_ptRelAll",&_ptRelAll);
    /*
    if(!istriggerstudy&&!frmeasurement) outputTree->Branch("_closeJetPtAllMC",&_closeJetPtAllMC);
    if(!istriggerstudy&&!frmeasurement) outputTree->Branch("_closeJetPtAllstatus",&_closeJetPtAllstatus);
    if(!istriggerstudy&&!frmeasurement) outputTree->Branch("_partonIdMatched",&_partonIdMatched);
    if(!istriggerstudy&&!frmeasurement) outputTree->Branch("_sameParton",&_sameParton);
    if(!frmeasurement)outputTree->Branch("_isloose",&_isloose);
    if(!frmeasurement)outputTree->Branch("_issoftmuon",&_issoftmuon);
    if(!frmeasurement)outputTree->Branch("_issoftmuonwithiso",&_issoftmuonwithiso);
     */

    if(!frmeasurement)outputTree->Branch("_istight",&_istight);
    if(!frmeasurement)outputTree->Branch("_isMediumMuon",&_isMediumMuon);
    if(!frmeasurement)outputTree->Branch("_istightIso",&_istightIso);
    if(!frmeasurement)outputTree->Branch("_istightID",&_istightID);
    //outputTree->Branch("_isFOIDWP2016_noIP3D", &_isFOIDWP2016_noIP3D);
    //outputTree->Branch("_istightIDWP2016_noIP3D", &_istightIDWP2016_noIP3D);
    outputTree->Branch("_istightIDWP2016",&_istightIDWP2016);
    //outputTree->Branch("_istightIDWP2016_IsoEmul",&_istightIDWP2016_IsoEmul);
    outputTree->Branch("_isFOIDWP2016",&_isFOIDWP2016);
    //outputTree->Branch("_isFOIDWP2016_IsoEmul",&_isFOIDWP2016_IsoEmul);
    outputTree->Branch("_isVetoIDWP2016_RA7",&_isVetoIDWP2016_RA7);
    outputTree->Branch("_istightIDWP2016_RA7",&_istightIDWP2016_RA7);
    outputTree->Branch("_isFOIDWP2016_RA7",&_isFOIDWP2016_RA7);
    outputTree->Branch("_isVetoIDWP2016_EWK",&_isVetoIDWP2016_EWK);
    outputTree->Branch("_isFOIDWP2016_EWK",&_isFOIDWP2016_EWK);
    outputTree->Branch("_istightIDWP2016_EWK",&_istightIDWP2016_EWK);
    //if(!frmeasurement)outputTree->Branch("_islooseMT2",&_islooseMT2);
    //if(!frmeasurement)outputTree->Branch("_iscutbasedWPtight",&_iscutbasedWPtight);
    outputTree->Branch("_istrigemulID",&_istrigemulID);
    outputTree->Branch("_istrigemulISO",&_istrigemulISO);
    /*
    if(!istriggerstudy&&!frmeasurement) outputTree->Branch("_mompt",&_mompt);
    if(!istriggerstudy&&!frmeasurement) outputTree->Branch("_momphi",&_momphi);
    if(!istriggerstudy&&!frmeasurement) outputTree->Branch("_mometa",&_mometa);
    if(!istriggerstudy&&!frmeasurement) outputTree->Branch("_mompdg",&_mompdg);
     */
    outputTree->Branch("_lPt",&_lPt);
    outputTree->Branch("_lEta",&_lEta);
    outputTree->Branch("_lPhi",&_lPhi);
    outputTree->Branch("_lE",&_lE);
    /*
    if(!istriggerstudy&&!frmeasurement) outputTree->Branch("_lPtmc",&_lPtmc);
    if(!istriggerstudy&&!frmeasurement) outputTree->Branch("_lEtamc",&_lEtamc);
    if(!istriggerstudy&&!frmeasurement) outputTree->Branch("_lPhimc",&_lPhimc);
    if(!istriggerstudy&&!frmeasurement) outputTree->Branch("_lEmc",&_lEmc);
    
    if(!istriggerstudy&&!frmeasurement) outputTree->Branch("_lgenPt",&_lgenPt);
    if(!istriggerstudy&&!frmeasurement) outputTree->Branch("_lgenEta",&_lgenEta);
    if(!istriggerstudy&&!frmeasurement) outputTree->Branch("_lgenPhi",&_lgenPhi);
    if(!istriggerstudy&&!frmeasurement) outputTree->Branch("_lgenpdgId",&_lgenpdgId);
     */

    outputTree->Branch("_genWPhi",&_genWPhi);
    outputTree->Branch("_genZPhi",&_genZPhi);
    outputTree->Branch("_genWEta",&_genWEta);
    outputTree->Branch("_genZEta",&_genZEta);
    outputTree->Branch("_genWPt",&_genWPt);
    outputTree->Branch("_genZPt",&_genZPt);
    outputTree->Branch("_genWE",&_genWE);
    outputTree->Branch("_genZE",&_genZE);
    
    /*
    if(!istriggerstudy&&!frmeasurement) outputTree->Branch("_nuPtmc",&_nuPtmc);
    if(!istriggerstudy&&!frmeasurement) outputTree->Branch("_nuEtamc",&_nuEtamc);
    if(!istriggerstudy&&!frmeasurement) outputTree->Branch("_nuPhimc",&_nuPhimc);
    if(!istriggerstudy&&!frmeasurement) outputTree->Branch("_nuEmc",&_nuEmc);
    if(!istriggerstudy&&!frmeasurement) outputTree->Branch("_mtmc",&_mtmc);
    if(!frmeasurement)outputTree->Branch("_lnmisshits",&_lnmisshits);
    if(!frmeasurement)outputTree->Branch("_passconvveto",&_passconvveto);
    if(!istriggerstudy&&!frmeasurement) outputTree->Branch("_lsietaieta",&_lsietaieta);
    if(!istriggerstudy&&!frmeasurement) outputTree->Branch("_lfull5x5sietaieta",&_lfull5x5sietaieta);
    if(!istriggerstudy&&!frmeasurement) outputTree->Branch("_lhovere",&_lhovere);
    if(!istriggerstudy&&!frmeasurement) outputTree->Branch("_ldphiin",&_ldphiin);    
    if(!istriggerstudy&&!frmeasurement) outputTree->Branch("_ldetain",&_ldetain);
    if(!istriggerstudy&&!frmeasurement) outputTree->Branch("_l1oemin1op",&_l1oemin1op);
    */
    
    //Tau leptons
    if(!istriggerstudy&&!frmeasurement)outputTree->Branch("_tau_dz",&_tau_dz);
    if(!istriggerstudy&&!frmeasurement)outputTree->Branch("_TauIs_againstMuonLoose3",&_TauIs_againstMuonLoose3);
    if(!istriggerstudy&&!frmeasurement)outputTree->Branch("_TauIs_byLooseCombinedIsolationDeltaBetaCorr3Hits",&_TauIs_byLooseCombinedIsolationDeltaBetaCorr3Hits);
    if(!istriggerstudy&&!frmeasurement)outputTree->Branch("_TauIs_decayModeFindingNewDMs",&_TauIs_decayModeFindingNewDMs);
    if(!istriggerstudy&&!frmeasurement)outputTree->Branch("_TauIs_decayModeFinding",&_TauIs_decayModeFinding);
    if(!istriggerstudy&&!frmeasurement)outputTree->Branch("_TauIs_byLooseIsolationMVArun2v1DBdR03oldDMwLT",&_TauIs_byLooseIsolationMVArun2v1DBdR03oldDMwLT);
    if(!istriggerstudy&&!frmeasurement)outputTree->Branch("_TauIs_byMediumIsolationMVArun2v1DBdR03oldDMwLT",&_TauIs_byMediumIsolationMVArun2v1DBdR03oldDMwLT);
    if(!istriggerstudy&&!frmeasurement)outputTree->Branch("_TauIs_byVLooseIsolationMVArun2v1DBoldDMwLT",&_TauIs_byVLooseIsolationMVArun2v1DBoldDMwLT);
    if(!istriggerstudy&&!frmeasurement)outputTree->Branch("_TauIs_byTightIsolationMVArun2v1DBoldDMwLT",&_TauIs_byTightIsolationMVArun2v1DBoldDMwLT);


    /*
    if(!frmeasurement)outputTree->Branch("_lchargeGSF",&_lchargeGSF);
    if(!frmeasurement)outputTree->Branch("_lchargeCTF",&_lchargeCTF);
    if(!frmeasurement)outputTree->Branch("_lchargePixSc",&_lchargePixSc);
     */
    


    //if(!istriggerstudy&&!frmeasurement) outputTree->Branch("_jetJECuncty",&_jetJECuncty);
    outputTree->Branch("_bTagged",&_bTagged);
    outputTree->Branch("_jetEta",&_jetEta);
    outputTree->Branch("_jetPhi",&_jetPhi);
    
    
    outputTree->Branch("_jetPt",&_jetPt);
    outputTree->Branch("_jetM",&_jetM);
    outputTree->Branch("_csv",&_csv);
    /*
    if(!istriggerstudy&&!frmeasurement) outputTree->Branch("_jetpFlav",&_jetpFlav);
    if(!istriggerstudy&&!frmeasurement) outputTree->Branch("_jethFlav",&_jethFlav);
    if(!istriggerstudy&&!frmeasurement) outputTree->Branch("_jetbtagSF",&_jetbtagSF);
    if(!istriggerstudy&&!frmeasurement) outputTree->Branch("_jetbtagSF_up",&_jetbtagSF_up);
    if(!istriggerstudy&&!frmeasurement) outputTree->Branch("_jetbtagSF_down",&_jetbtagSF_down);
    if(!istriggerstudy&&!frmeasurement) outputTree->Branch("_jetbtagEff",&_jetbtagEff);
     */
    
    //if(!istriggerstudy&&!frmeasurement) outputTree->Branch("_jetDeltaRloose",&_jetDeltaRloose);
    
    if(Skim=="skim2LrecoandTriggerTP" || Skim=="skim2GenL"){
    outputTree->Branch("passleg1L",&passleg1L);


    outputTree->Branch("hltL1sL1TripleMu553",&hltL1sL1TripleMu553);
    outputTree->Branch("hltL1TripleMu553L2TriMuFiltered3L3TriMuFiltered12105",&hltL1TripleMu553L2TriMuFiltered3L3TriMuFiltered12105);
    outputTree->Branch("hltL1TripleMu553L2TriMuFiltered3L3TriMuFiltered10105",&hltL1TripleMu553L2TriMuFiltered3L3TriMuFiltered10105);
    outputTree->Branch("hltL1TripleMu553L2TriMuFiltered3L3TriMuFiltered5",&hltL1TripleMu553L2TriMuFiltered3L3TriMuFiltered5);
    outputTree->Branch("hltDiMu9Ele9CaloIdLTrackIdLMuonlegL3Filtered9",&hltDiMu9Ele9CaloIdLTrackIdLMuonlegL3Filtered9);
    outputTree->Branch("hltDiMu9Ele9CaloIdLTrackIdLElectronlegDphiFilter",&hltDiMu9Ele9CaloIdLTrackIdLElectronlegDphiFilter);
    outputTree->Branch("hltMu8DiEle12CaloIdLTrackIdLMuonlegL3Filtered8",&hltMu8DiEle12CaloIdLTrackIdLMuonlegL3Filtered8);
    outputTree->Branch("hltMu8DiEle12CaloIdLTrackIdLElectronlegDphiFilter",&hltMu8DiEle12CaloIdLTrackIdLElectronlegDphiFilter);
    outputTree->Branch("hltL1sL1TripleEG14108",&hltL1sL1TripleEG14108);
    outputTree->Branch("hltEle16Ele12Ele8CaloIdLTrackIdLDphiLeg1Filter",&hltEle16Ele12Ele8CaloIdLTrackIdLDphiLeg1Filter);
    outputTree->Branch("hltEle16Ele12Ele8CaloIdLTrackIdLDphiLeg2Filter",&hltEle16Ele12Ele8CaloIdLTrackIdLDphiLeg2Filter);
    outputTree->Branch("hltEle16Ele12Ele8CaloIdLTrackIdLDphiLeg3Filter",&hltEle16Ele12Ele8CaloIdLTrackIdLDphiLeg3Filter);
    outputTree->Branch("hltL1sL1Mu6HTT150ORL1Mu8HTT125ORL1HTT175",&hltL1sL1Mu6HTT150ORL1Mu8HTT125ORL1HTT175);
    outputTree->Branch("hltDoubleMu8Mass8L3Filtered",&hltDoubleMu8Mass8L3Filtered);
    outputTree->Branch("hltL1sL1DoubleEG6HTT150orL1HTT175",&hltL1sL1DoubleEG6HTT150orL1HTT175);
    outputTree->Branch("hltDoubleEle8CaloIdMGsfTrackIdMDphiFilter",&hltDoubleEle8CaloIdMGsfTrackIdMDphiFilter);
    outputTree->Branch("hltMu8Ele8CaloIdMGsfTrackIdMDphiFilter",&hltMu8Ele8CaloIdMGsfTrackIdMDphiFilter);
    outputTree->Branch("hltMuon8L3Filtered0",&hltMuon8L3Filtered0);
    outputTree->Branch("hltL1sL1DoubleMu103p5ORDoubleMu125",&hltL1sL1DoubleMu103p5ORDoubleMu125);
    outputTree->Branch("hltDiMuonGlb17Glb8RelTrkIsoFiltered0p4",&hltDiMuonGlb17Glb8RelTrkIsoFiltered0p4);
    outputTree->Branch("hltL3fL1sDoubleMu103p5L1f0L2f10OneMuL3Filtered17",&hltL3fL1sDoubleMu103p5L1f0L2f10OneMuL3Filtered17);
    outputTree->Branch("hltL3pfL1sDoubleMu103p5L1f0L2pf0L3PreFiltered8",&hltL3pfL1sDoubleMu103p5L1f0L2pf0L3PreFiltered8);
    outputTree->Branch("hltL3fL1sDoubleMu103p5L1f0L2f10L3Filtered17",&hltL3fL1sDoubleMu103p5L1f0L2f10L3Filtered17);
    outputTree->Branch("hltDiMuonGlbFiltered17TrkFiltered8",&hltDiMuonGlbFiltered17TrkFiltered8);
    outputTree->Branch("hltDiMuonGlb17Trk8RelTrkIsoFiltered0p4",&hltDiMuonGlb17Trk8RelTrkIsoFiltered0p4);
    outputTree->Branch("hltDiMuonGlb17Glb8RelTrkIsoFiltered0p4DzFiltered0p2",&hltDiMuonGlb17Glb8RelTrkIsoFiltered0p4DzFiltered0p2);
    outputTree->Branch("hltDiMuonGlb17Trk8RelTrkIsoFiltered0p4DzFiltered0p2",&hltDiMuonGlb17Trk8RelTrkIsoFiltered0p4DzFiltered0p2);
    outputTree->Branch("hltL1sL1DoubleEG2210",&hltL1sL1DoubleEG2210);
    outputTree->Branch("hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter",&hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter);
    outputTree->Branch("hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter",&hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter);
    outputTree->Branch("hltEle23Ele12CaloIdLTrackIdLIsoVLDZFilter",&hltEle23Ele12CaloIdLTrackIdLIsoVLDZFilter);
    outputTree->Branch("hltL1sL1DoubleEG1510",&hltL1sL1DoubleEG1510);
    outputTree->Branch("hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter",&hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter);
    outputTree->Branch("hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter",&hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter);
    outputTree->Branch("hltEle17Ele12CaloIdLTrackIdLIsoVLDZFilter",&hltEle17Ele12CaloIdLTrackIdLIsoVLDZFilter);

    outputTree->Branch("hltL1sL1Mu20EG10",&hltL1sL1Mu20EG10);
    outputTree->Branch("hltMu23TrkIsoVVLEle8CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered23",&hltMu23TrkIsoVVLEle8CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered23);
    outputTree->Branch("hltMu23TrkIsoVVLEle8CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter",&hltMu23TrkIsoVVLEle8CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter);

    outputTree->Branch("hltL1sL1Mu12EG10",&hltL1sL1Mu12EG10);
    outputTree->Branch("hltMu17TrkIsoVVLEle12CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered17",&hltMu17TrkIsoVVLEle12CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered17);
    outputTree->Branch("hltMu17TrkIsoVVLEle12CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter",&hltMu17TrkIsoVVLEle12CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter);


    outputTree->Branch("hltL1sL1Mu5EG20",&hltL1sL1Mu5EG20);
    outputTree->Branch("hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered8",&hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered8);
    outputTree->Branch("hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter",&hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter);

    outputTree->Branch("hltL1sL1Mu5EG15",&hltL1sL1Mu5EG15);
    outputTree->Branch("hltMu8TrkIsoVVLEle17CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered8",&hltMu8TrkIsoVVLEle17CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered8);
    outputTree->Branch("hltMu8TrkIsoVVLEle17CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter",&hltMu8TrkIsoVVLEle17CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter);




    outputTree->Branch("hltL3crIsoL1sMu20Eta2p1L1f0L2f10QL3f24QL3trkIsoFiltered0p09",&hltL3crIsoL1sMu20Eta2p1L1f0L2f10QL3f24QL3trkIsoFiltered0p09);

    outputTree->Branch("hltL3crIsoL1sMu16L1f0L2f10QL3f18QL3trkIsoFiltered0p09",&hltL3crIsoL1sMu16L1f0L2f10QL3f18QL3trkIsoFiltered0p09);
    outputTree->Branch("hltL3fL1sMu16L1f0Tkf18QL3trkIsoFiltered0p09",&hltL3fL1sMu16L1f0Tkf18QL3trkIsoFiltered0p09);
    outputTree->Branch("hltL3crIsoL1sMu20L1f0L2f10QL3f22QL3trkIsoFiltered0p09",&hltL3crIsoL1sMu20L1f0L2f10QL3f22QL3trkIsoFiltered0p09);
    outputTree->Branch("hltL3fL1sMu20L1f0Tkf22QL3trkIsoFiltered0p09",&hltL3fL1sMu20L1f0Tkf22QL3trkIsoFiltered0p09);

    outputTree->Branch("hltL3crIsoL1sMu16L1f0L2f10QL3f20QL3trkIsoFiltered0p09",&hltL3crIsoL1sMu16L1f0L2f10QL3f20QL3trkIsoFiltered0p09);
    outputTree->Branch("hltEle23WPLooseGsfTrackIsoFilter",&hltEle23WPLooseGsfTrackIsoFilter);
    outputTree->Branch("hltEle27WPLooseGsfTrackIsoFilter",&hltEle27WPLooseGsfTrackIsoFilter);
    outputTree->Branch("hltEle27noerWPLooseGsfTrackIsoFilter",&hltEle27noerWPLooseGsfTrackIsoFilter);
    outputTree->Branch("hltEle32WPLooseGsfTrackIsoFilter",&hltEle32WPLooseGsfTrackIsoFilter);
    outputTree->Branch("hltEle23CaloIdLTrackIdLIsoVLTrackIsoFilter",&hltEle23CaloIdLTrackIdLIsoVLTrackIsoFilter);
    outputTree->Branch("hltDiEle33CaloIdLNewPixelMatchUnseededFilter",&hltDiEle33CaloIdLNewPixelMatchUnseededFilter);
    
    outputTree->Branch("hltL1sL1SingleMu16ORSingleMu25",&hltL1sL1SingleMu16ORSingleMu25);
    outputTree->Branch("hltDiMuonGlb27Trk8DzFiltered0p2",&hltDiMuonGlb27Trk8DzFiltered0p2);
    outputTree->Branch("hltL3fL1sMu16orMu25L1f0L2f25L3Filtered27",&hltL3fL1sMu16orMu25L1f0L2f25L3Filtered27);
    outputTree->Branch("hltL3fL1sMu16orMu25L1f0L2f10QL3Filtered30Q",&hltL3fL1sMu16orMu25L1f0L2f10QL3Filtered30Q);
    outputTree->Branch("hltEle30CaloIdLGsfTrkIdVLDPhiUnseededFilter",&hltEle30CaloIdLGsfTrkIdVLDPhiUnseededFilter);
    outputTree->Branch("hltDiEle33CaloIdLGsfTrkIdVLDPhiUnseededFilter",&hltDiEle33CaloIdLGsfTrkIdVLDPhiUnseededFilter);
    outputTree->Branch("hltL3fL1sMu22orMu25orMu20EG15orMu5EG20L1f0L2f10QL3Filtered30Q",&hltL3fL1sMu22orMu25orMu20EG15orMu5EG20L1f0L2f10QL3Filtered30Q);

    outputTree->Branch("hltL3fL1sL1DoubleMu0ETM40lorDoubleMu0ETM55", &hltL3fL1sL1DoubleMu0ETM40lorDoubleMu0ETM55) ;
    outputTree->Branch("hltEle22WPLooseL1SingleIsoEG20erGsfTrackIsoFilter", &hltEle22WPLooseL1SingleIsoEG20erGsfTrackIsoFilter);
    outputTree->Branch("hltL3crIsoL1sSingleMu18erIorSingleMu20erL1f0L2f10QL3f19QL3trkIsoFiltered0p09", &hltL3crIsoL1sSingleMu18erIorSingleMu20erL1f0L2f10QL3f19QL3trkIsoFiltered0p09);
    outputTree->Branch("hltL1sSingleMu16",&hltL1sSingleMu16);
    outputTree->Branch("hltL1sSingleMu18",&hltL1sSingleMu18);
    outputTree->Branch("hltL1sSingleMu20",&hltL1sSingleMu20);
    outputTree->Branch("hltL1sSingleMu20erlorSingleMu22lorSingleMu25",&hltL1sSingleMu20erlorSingleMu22lorSingleMu25);

    outputTree->Branch("hltL1sSingleEG40IorSingleIsoEG22erIorSingleIsoEG24erIorSingleIsoEG24IorSingleIsoEG26",&hltL1sSingleEG40IorSingleIsoEG22erIorSingleIsoEG24erIorSingleIsoEG24IorSingleIsoEG26);

    }
}




int SSb13_takeimai::GenParticleInfo(const edm::Event& iEvent){

  _mgluino=0; 
  _mLSP =0;
  _mchargino=0; 
  _mstop=0;


  if(!IsMC) return 0;

  
  edm::InputTag pileupSrc_(IT_pileup);

  Handle<std::vector<PileupSummaryInfo> > puInfo;
  if(isgetbylabel)   iEvent.getByLabel(pileupSrc_, puInfo);
  else iEvent.getByToken(puInfo_token, puInfo); 
  
  vector<PileupSummaryInfo>::const_iterator pvi;
  //  int trueNVtx(0),nVtxNow(0),nVtxBefore(0),nVtxAfter(0);
  for (pvi = puInfo->begin(); pvi != puInfo->end(); ++pvi) {
    if (pvi->getBunchCrossing() == -1) nVtxBefore = pvi->getPU_NumInteractions();
    else if (pvi->getBunchCrossing() == 0) {
      trueNVtx = pvi->getTrueNumInteractions(); // from Fall11 onwards
      nVtxNow = pvi->getPU_NumInteractions();
      Nvtx->Fill(TMath::Min(trueNVtx,49));
    }
    else if (pvi->getBunchCrossing() == 1) nVtxAfter = pvi->getPU_NumInteractions();
    //cout << " Pileup Information: bunchXing, nvtx: " << pvi->getBunchCrossing() << " " << pvi->getTrueNumInteractions() << endl;
  }

  if(isgetbylabel) iEvent.getByLabel("prunedGenParticles", TheGenParticles);
  else iEvent.getByToken(genpart_token, TheGenParticles);

  for(unsigned int i = 0; i < TheGenParticles->size();i++){
    const reco::GenParticle &genp = (*TheGenParticles)[i];

    if(fabs(genp.pdgId())==1000021 && genp.mass()> _mgluino) _mgluino = genp.mass();
    if(fabs(genp.pdgId())==1000022 && genp.mass()>_mLSP ) _mLSP = genp.mass();
    if(fabs(genp.pdgId())==1000024 && genp.mass()>_mchargino) _mchargino = genp.mass();
    if(fabs(genp.pdgId())==1000006 && genp.mass()>_mstop) _mstop = genp.mass();
    /*    if(fabs(genp.pdgId())>1000000){
      cout << "GenParticle: PdgId, pt, eta, phi, M: "
	   << genp.pdgId() << ", "<< genp.pt() << ", "<<genp.eta() << ", "<<genp.phi() << ", "<<genp.mass()<<endl;
	   }*/
    

  }

  std::vector<const GenParticle*> vGenElectrons, vGenMuons, vGenNPElectrons, vGenNPMuons, vGenW;
  if( TheGenParticles.isValid() ) {
    GPM.SetCollection(TheGenParticles);
    GPM.Classify();
    vGenMuons = GPM.filterByStatus(GPM.getPromptMuons(),statusgenparticles);
    vGenElectrons = GPM.filterByStatus(GPM.getPromptElectrons(),statusgenparticles);
  }
  for(unsigned int i = 0; i < vGenMuons.size();i++){
  const reco::GenParticle *iM = vGenMuons[i];
    _lgenPt.push_back(iM->pt());
    _lgenEta.push_back(iM->eta());
    _lgenPhi.push_back(iM->phi());
    _lgenpdgId.push_back(iM->pdgId());

    if(vGenMuons.size()+vGenElectrons.size()>20) {cout << "Gen muon. Pt: " << iM->pt()<< 
      " , Eta: " <<iM->eta()<< 
      " , Phi: " <<iM->phi()<<  
      " , PdgId: " <<iM->pdgId()<<  
      " , Status: " <<iM->status()<<  
      endl;
      TLorentzVector genmomv; 
      genmomv.SetPtEtaPhiE(iM->mother()->pt(),iM->mother()->eta(),iM->mother()->phi(),iM->mother()->energy());
      cout <<"Mother, pdgid, status: " << iM->mother()->pdgId() <<","<<iM->mother()->status()
	   << endl << "pt, eta, phi, M: " << iM->mother()->pt()<<","<<iM->mother()->eta()<<","<<iM->mother()->phi()<<","<<genmomv.M()
	   << endl;}
  }

  for(unsigned int i = 0; i <vGenElectrons.size();i++){
   const reco::GenParticle *iE = vGenElectrons[i];
   _lgenPt.push_back(iE->pt());
   _lgenEta.push_back(iE->eta());
   _lgenPhi.push_back(iE->phi());
   _lgenpdgId.push_back(iE->pdgId());
   if(vGenMuons.size()+vGenElectrons.size()>20  ) {cout << "Gen ele. Pt: " << iE->pt()<< 
      " , Eta: " <<iE->eta()<< 
      " , Phi: " <<iE->phi()<<  
      " , PdgId: " <<iE->pdgId()<<  
      " , Status: " <<iE->status()<<  
      endl;
     TLorentzVector genmomv;
     genmomv.SetPtEtaPhiE(iE->mother()->pt(),iE->mother()->eta(),iE->mother()->phi(),iE->mother()->energy());
     cout <<"Mother, pdgid, status: " << iE->mother()->pdgId() <<","<<iE->mother()->status() 
	  << endl << "pt, eta, phi, M: " << iE->mother()->pt()<<","<<iE->mother()->eta()<<","<<iE->mother()->phi()<<","<< genmomv.M()
	  << endl;}
  }

  return TheGenParticles->size();
}


void SSb13_takeimai::MatchToGenParticle( const edm::Event& iEvent,const int &leptonctr, const double &recolpt, const double &recoleta, const double &recolphi, const int &recolpdgid, vector <bool> & checkusedpart){
  int genpart_i = 0;

  double dR_closegenl = 0.2;
  double dR_closegenp = 0.2;
  bool quarkcand[3] ={false,false,false};
  bool nonpromptlepton =false;
  enum quarktype{udsg=0,c, b};

  if(isgetbylabel) iEvent.getByLabel("prunedGenParticles", TheGenParticles);
  else iEvent.getByToken(genpart_token, TheGenParticles);

  //  cout <<endl<<endl; 
  
  //  cout << "Lepton Pt, Eta, Phi, PdgId " << recolpt<<", "<<recoleta<<", "<< recolphi<<", "<<recolpdgid<<endl<<endl;
  for(unsigned int i = 0; i < TheGenParticles->size();i++){
    const reco::GenParticle &genp = (*TheGenParticles)[i];
    int provcode = _origin[leptonctr];

    int thepdgid = genp.pdgId();
    if(genp.pt()<1) continue;
    //Match in dr<0.2
    double dR = deltaR( genp.eta(),genp.phi(), recoleta,recolphi);
    if(dR>0.2)continue; 

    //cout << "Genparticle Pt, dR, PdgId, Status, IsPrompt " 
    //	 <<genp.pt()<<", "<<dR/*<<genp.eta()<<", "<< genp.phi()*/<<", "<< genp.pdgId()<<", "<<genp.status()<<", " << (genp.isPromptFinalState() || genp.isDirectPromptTauDecayProductFinalState() || genp.isHardProcess()) <<endl;
    //
    if(genp.status()!=1 &&genp.status()!=2&&genp.status()!=71) continue; 
    
    if(abs(thepdgid)<=3||abs(thepdgid)==21)quarkcand[udsg]=true;
    if(abs(thepdgid)==4)quarkcand[c]=true;
    if(abs(thepdgid)==5)quarkcand[b]=true;
    
    
    if(dR> dR_closegenp && dR> dR_closegenl) continue; 
    if(abs(thepdgid) ==12 || abs(thepdgid) ==14 || abs(thepdgid) ==16) continue;
    //If the genparticle has already been used for another matching, discard it. 
    if(checkusedpart[i] ) continue;
    


    bool isgenl = false;

    if(genp.status()==1 && abs(thepdgid) == abs(recolpdgid) ) isgenl = true;
    if(genp.status()==1 && abs(thepdgid) == 15) isgenl = true;
    if(genp.status()==2 && abs(thepdgid) == 15) isgenl = true;
    if(genp.status()==2 && abs(thepdgid) != 15) continue;
    //If a genl is found, don't check non gen leptons
    
    if(!isgenl && abs(provcode)<=3 ) continue; 
    
    if(isgenl && dR>dR_closegenl) continue;
    if(!isgenl && dR>dR_closegenp&&!nonpromptlepton) continue;
    

    bool isprompt =  genp.isPromptFinalState() || genp.isDirectPromptTauDecayProductFinalState() || genp.isHardProcess() ;
    
    nonpromptlepton = !isprompt && abs(thepdgid) == abs(recolpdgid) ;    

    //Find mother 
    reco::GenParticle genmom ;
   
    for(unsigned int j = 0; j < genp.numberOfMothers();j++){
      const reco::GenParticle genmothercand = *(genp.motherRef(j));
      if(genp. numberOfMothers() ==1 || genmothercand.pdgId() == thepdgid){
	genmom =genmothercand;
      	break;
      }  
    }
    //Repeat previous step until one doesn't find intermediate particles    

    bool keeplooping = true;
    while(genmom.pdgId() == thepdgid&&keeplooping){
      if(genmom.numberOfMothers()==0) break;
      // cout << "****"<<endl;
      keeplooping = false;
      for(unsigned int j = 0; j < genmom.numberOfMothers();j++){
	
	const reco::GenParticle genmothercand = *(genmom.motherRef(j));

	//cout << "deb 1 "<<thepdgid<<","<< genmothercand.pdgId()<<endl;
	if(genmom.numberOfMothers() ==1||genmothercand.pdgId() == thepdgid){
	  genmom =genmothercand;
	  keeplooping = true;
	  break;
	}  

      }  
    }


    //Find grandmother 
    
    reco::GenParticle gengrandmom ;
   
    for(unsigned int j = 0; j < genmom.numberOfMothers();j++){
      const reco::GenParticle gengrandmothercand = *(genmom.motherRef(j));

      if(genmom. numberOfMothers() ==1 || gengrandmothercand.pdgId() == genmom.pdgId()){
	gengrandmom = gengrandmothercand;
	break; 
      }  
    }
    //Repeat previous step until one doesn't find intermediate particles    
    keeplooping = true;
    while(gengrandmom.pdgId() == genmom.pdgId()&&keeplooping){
      if(gengrandmom.numberOfMothers()==0) break;
      //      	cout << "****"<<endl;
	keeplooping = false;

      for(unsigned int j = 0; j < gengrandmom.numberOfMothers();j++){
	const reco::GenParticle gengrandmothercand = *(gengrandmom.motherRef(j));
	//cout << "deb 2 "<<thepdgid<<","<< gengrandmothercand.pdgId()<<endl;

	if(gengrandmom.numberOfMothers() ==1 || gengrandmothercand.pdgId() == genmom.pdgId()){
	  gengrandmom = gengrandmothercand;
	  keeplooping =  true;
	  break;
	}  
      }  
    }


    int pdgidmom = genmom.pdgId();
    int pdgidgrandmom = gengrandmom.pdgId();
      
    bool bestparticle = true;

    bool isgenemu = (abs(thepdgid)==11 ||  abs(thepdgid)==13)&&abs(pdgidmom) !=15 &&abs(pdgidgrandmom) !=15;// || abs(thepdgid)==15 ;
    bool isemufromtau = (abs(thepdgid)==11 ||  abs(thepdgid)==13)&& (abs(pdgidmom)==15|| abs(pdgidgrandmom)==15);
    //e/mu cases:
    if(isprompt &&isgenemu &&thepdgid==recolpdgid ) provcode =0;   
    //Charge flip:
    else if(isprompt &&isgenemu && thepdgid==-recolpdgid) provcode =1;  
    else if(isgenemu &&thepdgid==recolpdgid &&(pdgidmom == 23 || pdgidmom == 24 || pdgidmom ==-24  || pdgidmom == 1000024  || pdgidmom ==  -1000024  )&&thepdgid*recolpdgid>0) provcode =0;
    else if(isgenemu &&thepdgid==recolpdgid &&(pdgidmom == 23 || pdgidmom == 24 || pdgidmom ==-24  || pdgidmom == 1000024  || pdgidmom ==  -1000024  )&&thepdgid*recolpdgid<0) provcode =1;


    //Leptonic Taus:
    else if(isprompt &&isemufromtau&&thepdgid*recolpdgid >0) provcode =2;
    else if(isprompt &&isemufromtau&&thepdgid*recolpdgid<0 ) provcode =3;
    //Hadronic Taus:
    else if( abs(thepdgid)==15&& thepdgid==recolpdgid )  provcode =2;
    else if( abs(thepdgid)==15&& thepdgid==-recolpdgid)  provcode =3;

    // (here we verify the lepton comes from a W/Z, I don't know why...): 
    else if(( abs(thepdgid)==15 /*|| abs(thepdgid) == abs(recolpdgid)*/ )&& 
	(pdgidmom == 23 || pdgidmom == 24 || pdgidmom ==-24  || pdgidmom == 1000024  || pdgidmom ==  -1000024  ) 
	&&thepdgid*recolpdgid<0)  provcode =3; 
    else if( abs(thepdgid)==15&& 
	(pdgidmom == 23 || pdgidmom == 24 || pdgidmom ==-24  || pdgidmom == 1000024  || pdgidmom ==  -1000024  )
        &&thepdgid*recolpdgid>0)  provcode =2; 

    //Sometimes for taus, the mother is also a tau and one needs to check the grandmother
    else if( abs(thepdgid)==15 && abs(thepdgid) == abs(recolpdgid) &&  abs(pdgidmom) ==15 &&
	(pdgidgrandmom == 23 || pdgidgrandmom == 24 || pdgidgrandmom ==-24  || pdgidgrandmom == 1000024  || pdgidgrandmom ==  -1000024  )
	&&thepdgid*recolpdgid<0)  provcode =3; 
    
    else if( abs(thepdgid)==15&& abs(thepdgid) == abs(recolpdgid) &&  abs(pdgidmom) ==15 &&
	(pdgidgrandmom == 23 || pdgidgrandmom == 24 || pdgidgrandmom ==-24  || pdgidgrandmom == 1000024  || pdgidgrandmom ==  -1000024  )
	&&thepdgid*recolpdgid>0)  provcode =2; 

    //Photon conversion in tracker
    else if(isprompt && thepdgid==22) provcode = 5;
    //Internal photon conversion
    else if(isprompt && abs(thepdgid)==abs(recolpdgid) &&pdgidmom==22) provcode = 4;

    else if(quarkcand[b]&& nonpromptlepton  &&  abs(provcode)>=4)  provcode =51;
    else if(quarkcand[c] && nonpromptlepton &&  abs(provcode)>=4  ) provcode =41;
    else if(quarkcand[udsg] && nonpromptlepton  &&  abs(provcode)>=4) provcode =11;
    else if(abs(provcode)>=4 && _origin[leptonctr]!=51&& _origin[leptonctr]!=41&& _origin[leptonctr]!=11){
      //Check if particle is a b quark/hadron
      if(abs(thepdgid)==5||abs(thepdgid)==511||abs(thepdgid)==521||abs(thepdgid)==531||abs(thepdgid)==533||abs(thepdgid)==535||abs(thepdgid)==551||abs(thepdgid)==553) provcode =50; 
      //Check if particle is a c quark/hadron
      else if(_origin[leptonctr]!=50&&(abs(thepdgid)==4||abs(thepdgid)==411||abs(thepdgid)==421||abs(thepdgid)==431||abs(thepdgid)==433||abs(thepdgid)==441||abs(thepdgid)==443)) provcode =40;
      //Check if particle is a light quark/hadron
      else if(_origin[leptonctr]!=50&&_origin[leptonctr]!=40&&(abs(thepdgid)==1||abs(thepdgid)==2||abs(thepdgid)==3||abs(thepdgid)==111||abs(thepdgid)==211||abs(thepdgid)==130||abs(thepdgid)==210||abs(thepdgid)==321||abs(thepdgid)==310||abs(thepdgid)==21)) provcode =10;
      else bestparticle=false;
    }
    else bestparticle=false;

    //    if(abs(thepdgid)==22) cout <<"pdgid 22, " <<bestparticle<<", " << provcode<<endl;
    if(!bestparticle)continue;
    if(isgenl)dR_closegenl =dR; 
    else dR_closegenp =dR;
    /*
    if(provcode==1&&(abs(recolpdgid)==13 ||abs(recolpdgid)==11))cout <<"flip" <<endl;
    if(provcode==1&&(abs(recolpdgid)==13 ||abs(recolpdgid)==11)) cout <<"gen "<<  genp.pt()<<", "<<genp.eta()<<", "<<genp.phi()<<", "<<genp.pdgId()<<", "<< isprompt<<", " <<genp.status()<<endl;
    if(provcode==1&&(abs(recolpdgid)==13 ||abs(recolpdgid)==11)) cout << "reco " << recoleta<<", "<< recolphi<<", "<<recolpdgid<<endl;*/
    genpart_i =(int) i; 



    _lPtmc[leptonctr]=genp.pt();
    _lEtamc[leptonctr]=genp.eta();
    _lPhimc[leptonctr]=genp.phi();
    _lEmc[leptonctr]=genp.energy();
    _chargesMC[leptonctr]=genp.charge();
    _origin[leptonctr]=provcode;

    //    cout << "genmom/gmom " << genmom.pdgId()<< ", " <<gengrandmom.pdgId()<<endl;



  }
  checkusedpart[genpart_i]=true;
    //  GenParticle genmom ;


    /*    for(unsigned int j = 0; j < genp.numberOfMothers();i++){
      const reco::GenParticle &genmother = (*TheGenParticles)[i];
      ->motherRef(i)
    for (unsigned int i =0; i < genp. numberOfMothers();i++){
      if(genp. numberOfMothers() ==1) genmom = genmom[i];
      if(genmom[i].pdgId() == genp.pdgId){
	genmom = genmom[i];
	break;
      }
    }
    if(genmom !=0 && genmom.pdgId ==genp.pdgId) {
      while(genmom.pdgId ==genp.pdgId ){
	for (unsigned int i =0; i < genmomp. numberOfMothers();i++){
	  if()
	
	}
      }
    
      }*/
    
    
    //    if(genp.numberOfMothers()==1)

    //genmom 
    /*    while(genmom.pdgId() == genp. ){
      
    }
    for (unsigned int i =0; i < genp. numberOfMothers()){
      mom = gp.mother(i)
	    if mom.pdgId() == gp.pdgId():
    ret += realGenMothers(mom)
	      else:
		ret.append(mom)
    return ret

    int motherpdgid = gen->mother()->pdgId();
    while(gen->mother() )



    ->mother()

    }*/
  /*  cout<< _origin[leptonctr]<<endl;
  cout<< "Ptgen: "<< _lPtmc[leptonctr]<<endl;
  cout<<"**************"<<endl;
  std::string mylol;*/
   //   if(_origin[leptonctr]==1) cin >> mylol;
}

bool SSb13_takeimai::ApplySkim(const edm::Event& iEvent,std::string skim){
  if(skim == ""  ) return true;
  
  else if(skim == "MET100"  ) {
  //============ Pat MET ============                                                                                                                                                                    
  edm::Handle< vector<pat::MET> > ThePFMET;
  if(isgetbylabel)     iEvent.getByLabel(IT_pfmet, ThePFMET);
  else iEvent.getByToken(met_token,ThePFMET);
  if( ! ThePFMET.isValid() ) ERR( IT_pfmet );
  const vector<pat::MET> *pfmetcol = ThePFMET.product();
  const pat::MET *pfmet;
  pfmet = &(pfmetcol->front());
  if(pfmet->pt()<100 ) return false;
  return true;}


  else if(skim=="SignalSkimRA5"||skim=="FakeRateAppSkimRA5"||skim=="InSituSkimRA5"||skim=="SkimRA7"){
    //============ Pat MET ============                                                                                                          
    edm::Handle< vector<pat::MET> > ThePFMET;
    if(isgetbylabel)     iEvent.getByLabel(IT_pfmet, ThePFMET);
    else     iEvent.getByToken(met_token, ThePFMET);
    if( ! ThePFMET.isValid() ) ERR( IT_pfmet );
    const vector<pat::MET> *pfmetcol = ThePFMET.product();
    const pat::MET *pfmet;
    pfmet = &(pfmetcol->front());
    if (pfmet->pt()<30) return false;
    
    edm::Handle<TriggerResults> trigResults;
    if(isgetbylabel)    iEvent.getByLabel(IT_hltresults, trigResults);
    else     iEvent.getByToken(trgresults_token, trigResults);
    bool passtriggers = IsMC? true: false;
    if( !trigResults.failedToGet() ) {
      int N_Triggers = trigResults->size();

      const edm::TriggerNames & trigName = iEvent.triggerNames(*trigResults);
      for( int i_Trig = 0; i_Trig < N_Triggers; ++i_Trig ) {
	if (trigResults.product()->accept(i_Trig)) {
	  TString TrigPath =trigName.triggerName(i_Trig);
	  //      if(TrigPath.Contains("HLT_IsoMu24_eta2p1")) return true;                                                                      
	  if( TrigPath.Contains("HLT_DoubleMu8_Mass8_PFHT300")) passtriggers = true;
	  //if(TrigPath.Contains("HLT_DoubleMu8_Mass8_PFHT250")) passtriggers = true;
	  if(TrigPath.Contains("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL")) passtriggers = true;
	  if(TrigPath.Contains("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL")) passtriggers = true;
	  if(TrigPath.Contains("HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT300")) passtriggers = true;
	  //	  if(TrigPath.Contains("HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT250")) passtriggers = true;
	  if(TrigPath.Contains("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ")) passtriggers = true;
	  if(TrigPath.Contains("HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT300")) passtriggers = true;
	  //if(TrigPath.Contains("HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT250")) passtriggers = true;
	  if(TrigPath.Contains("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL")) passtriggers = true;
	  if(TrigPath.Contains("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL")) passtriggers = true;
	
	}
      }
    }
    edm::Handle< std::vector<pat::Muon> > thePatMuons;
    if(isgetbylabel)     iEvent.getByLabel( IT_muon, thePatMuons );
    else     iEvent.getByToken(patmuon_token,thePatMuons);
    if( ! thePatMuons.isValid() )  ERR(IT_muon) ;
    edm::Handle< std::vector<pat::Electron> > thePatElectrons;
    if(isgetbylabel)     iEvent.getByLabel( IT_electron, thePatElectrons );
    else iEvent.getByToken(patelectron_token,thePatElectrons);
    if( ! thePatElectrons.isValid() ) ERR( IT_electron );
    int leptonsize = thePatMuons->size()+thePatElectrons->size();
    if((!passtriggers&&!IsFS)||leptonsize<2) return false;
    else return true;
  }
  else  if(skim == "skim2GenL" || skim == "skim3GenL" ) {
    if(!IsMC) return true;
    if(isgetbylabel)     iEvent.getByLabel("prunedGenParticles", TheGenParticles);
    else     iEvent.getByToken(genpart_token, TheGenParticles);
    std::vector<const GenParticle*> vGenElectrons, vGenMuons, vGenNPElectrons, vGenNPMuons, vGenW;
    if( TheGenParticles.isValid() ) {
      GPM.SetCollection(TheGenParticles);
      GPM.Classify();
      vGenMuons = GPM.filterByStatus(GPM.getPromptMuons(),statusgenparticles);
      vGenElectrons = GPM.filterByStatus(GPM.getPromptElectrons(),statusgenparticles);
    }
    if(vGenMuons.size() + vGenElectrons.size()<2) return false;    
    if(vGenMuons.size() + vGenElectrons.size()<3 &&  skim == "skim3GenL" ) return false;    
    return true;
  }
  
  else  if(skim == "skim1Lreco" ||skim == "skim2Lreco" || skim == "skim3Lreco" || skim == "skim2LrecoPFMet30" ||  skim == "skim2LrecoandTriggerTP" ||  skim == "skim3LrecoandTriggerTP" ) {
    //============ Beamspot ============                                                                                                         
    edm::Handle< reco::BeamSpot > theBeamSpot;
    if(isgetbylabel)     iEvent.getByLabel( IT_beamspot, theBeamSpot );
    else iEvent.getByToken(beamspot_token, theBeamSpot );
    if( ! theBeamSpot.isValid() ) ERR( IT_beamspot ) ;
    BeamSpot::Point  BS= theBeamSpot->position();
    //============ Primary vertices ============  
    //edm::InputTag IT_goodVtx = edm::InputTag("offlineSlimmedPrimaryVertices");
    edm::Handle<std::vector<Vertex> > theVertices;
    if(isgetbylabel) iEvent.getByLabel( "goodOfflinePrimaryVertices", theVertices) ;
    else iEvent.getByToken(goodpv_token,theVertices) ;
    

    if( ! theVertices.isValid() ) ERR(IT_goodVtx ) ;
    int nvertex = theVertices->size();
    Vertex::Point PV(0,0,0);
    if( nvertex) PV = theVertices->begin()->position();
    else{
      cout << "[WARNING]: No candidate primary vertices passed the quality cuts, so skipping event" << endl;
      return false;
    } 
    //============ Pat Muons ============                                                                                                        
    edm::Handle< std::vector<pat::Muon> > thePatMuons;
    if(isgetbylabel)     iEvent.getByLabel( IT_muon, thePatMuons );
    else iEvent.getByToken(patmuon_token, thePatMuons); 
    if( ! thePatMuons.isValid() )  ERR(IT_muon) ;
    //============ Pat Electrons ============                                                                                                    
    edm::Handle< std::vector<pat::Electron> > thePatElectrons;
    if(isgetbylabel)     iEvent.getByLabel( IT_electron, thePatElectrons );
    else iEvent.getByToken(patelectron_token, thePatElectrons); 
    if( ! thePatElectrons.isValid() ) ERR( IT_electron );
      //============ Pat Taus ============
      edm::Handle< std::vector<pat::Tau> > thePatTaus;
      if(isgetbylabel)  iEvent.getByLabel( IT_tau, thePatTaus );
      else iEvent.getByToken(pattau_token, thePatTaus );
      if( ! thePatTaus.isValid() ) ERR( IT_tau );
    //============ Conversions ============                                                                                                      
    edm::Handle< std::vector<reco::Conversion> > theConversions;
    if(isgetbylabel)     iEvent.getByLabel("reducedEgamma","reducedConversions", theConversions);
    else iEvent.getByToken(conv_token,theConversions);


    if(skim == "skim1Lreco"){
      _eleMinPt = 10;
      _muonMinPt = 10;
      _jetPtCut =25;
    }

    std::vector<const pat::Muon* > sMu = ssbLooseMuonSelector( *thePatMuons, _muonMinPt, PV, _looseD0Mu, true);
    std::vector<const pat::Electron* > sEl = ssbMVAElectronSelector( *thePatElectrons, _eleMinPt, PV, _looseD0E, _chargeConsistency, _useconversions, theConversions , BS, false);
    std::vector<const pat::Tau* > sTau = ssbLooseTauSelector( *thePatTaus, _tauMinPt, PV, 10000, true);
    
    if (sEl.size() + sMu.size() <1)return false;
    if(skim == "skim1Lreco") return true;
    
    if (sEl.size() + sMu.size() + sTau.size()< 2) return false;
    if(sEl.size() + sMu.size() + sTau.size()< 3 && skim == "skim3Lreco") return false;
    
    //    if(sEl.size() + sMu.size() >=3  && skim == "skim3Lreco") return false;   
    
    //============ Pat MET ============                                                                                                          
    edm::Handle< vector<pat::MET> > ThePFMET;
    if(isgetbylabel)     iEvent.getByLabel(IT_pfmet, ThePFMET);
    else iEvent.getByToken(met_token,ThePFMET);
    if( ! ThePFMET.isValid() ) ERR( IT_pfmet );
    const vector<pat::MET> *pfmetcol = ThePFMET.product();
    const pat::MET *pfmet;
    pfmet = &(pfmetcol->front());
    if(pfmet->pt()<30 &&  skim == "skim2LrecoPFMet30" ) return false;

    //Trigger results 
    if(skim == "skim2LrecoandTriggerTP" ) {

      _jetPtCut = 30;
      edm::Handle<TriggerResults> trigResults;
      if(isgetbylabel)       iEvent.getByLabel(IT_hltresults, trigResults);
      else iEvent.getByToken(trgresults_token, trigResults);

      if( !trigResults.failedToGet() ) {
	int N_Triggers = trigResults->size();
	
	const edm::TriggerNames & trigName = iEvent.triggerNames(*trigResults);
	for( int i_Trig = 0; i_Trig < N_Triggers; ++i_Trig ) {
	  if (trigResults.product()->accept(i_Trig)) {
	    TString TrigPath =trigName.triggerName(i_Trig);
	    //	    if(TrigPath.Contains("HLT_IsoMu24_eta2p1")) return true;
	    if(TrigPath.Contains("HLT_IsoMu18")) return true;
	    if(TrigPath.Contains("HLT_IsoMu20")) return true;
	    if(TrigPath.Contains("HLT_IsoMu22")) return true;
	    if(TrigPath.Contains("HLT_IsoTkMu18")) return true;
	    if(TrigPath.Contains("HLT_IsoTkMu20")) return true;
	    if(TrigPath.Contains("HLT_IsoTkMu22")) return true;

	    if(TrigPath.Contains("HLT_Ele23_WPLoose_Gsf"))return true;
	    if(TrigPath.Contains("HLT_Ele27_WPLoose_Gsf"))return true;
	    if(TrigPath.Contains("HLT_Ele23_CaloIdL_TrackIdL_IsoVL"))return true;
	    //if(TrigPath.Contains("HLT_Ele27_eta2p1_WPLoose_Gsf"))return true;
	    //if(TrigPath.Contains("HLT_Ele32_eta2p1_WPLoose_Gsf"))return true;
	    //if(TrigPath.Contains("HLT_Ele27_eta2p1_WP75_Gsf"))return true;
	    //if(TrigPath.Contains("HLT_Ele32_eta2p1_WP75_Gsf"))return true;
	    if(TrigPath.Contains("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL"))return true;
            if(TrigPath.Contains("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ"))return true;
	    if(TrigPath.Contains("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL"))return true;
            if(TrigPath.Contains("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ"))return true;

            //if(TrigPath.Contains("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL"))return true;
            //if(TrigPath.Contains("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ"))return true;
            if(TrigPath.Contains("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL"))return true;
            if(TrigPath.Contains("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ"))return true;
            if(TrigPath.Contains("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL"))return true;
            if(TrigPath.Contains("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ"))return true;


            if(TrigPath.Contains("HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL")) return true;
	    if(TrigPath.Contains("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL")) return true;
	    if(TrigPath.Contains("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL"))  return true;

            
	    if(TrigPath.Contains("HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL"))return true;
	    if(TrigPath.Contains("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL"))return true;
	  }
	}
      }
      return false;
    }
    
    if(skim == "skim3LrecoandTriggerTP" ) {
      if(sEl.size() + sMu.size() < 3)return false;
      edm::Handle<TriggerResults> trigResults;
      if(isgetbylabel)       iEvent.getByLabel(IT_hltresults, trigResults);
      else iEvent.getByToken(trgresults_token, trigResults);
      if( !trigResults.failedToGet() ) {
        int N_Triggers = trigResults->size();
	
        const edm::TriggerNames & trigName = iEvent.triggerNames(*trigResults);
        for( int i_Trig = 0; i_Trig < N_Triggers; ++i_Trig ) {
          if (trigResults.product()->accept(i_Trig)) {
            TString TrigPath =trigName.triggerName(i_Trig);

	    if(TrigPath.Contains("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL"))return true;
            if(TrigPath.Contains("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ"))return true;
            if(TrigPath.Contains("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL"))return true;
            if(TrigPath.Contains("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ"))return true;

            if(TrigPath.Contains("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL"))return true;
            if(TrigPath.Contains("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ"))return true;
            if(TrigPath.Contains("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL"))return true;
            if(TrigPath.Contains("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ"))return true;


            if(TrigPath.Contains("HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL")) return true;
            if(TrigPath.Contains("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL")) return true;
            if(TrigPath.Contains("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL"))  return true;


            if(TrigPath.Contains("HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL"))return true;
            if(TrigPath.Contains("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL"))return true;


          }
        }
      }
      return false;
    }
    
    
    
    return true;

  } 
  
  else if(skim == "triggerTP" ) {
    edm::Handle<TriggerResults> trigResults;
    if(isgetbylabel) iEvent.getByLabel(IT_hltresults, trigResults);
    else    iEvent.getByToken(trgresults_token, trigResults);
    
    if( !trigResults.failedToGet() ) {
      int N_Triggers = trigResults->size();
      
      const edm::TriggerNames & trigName = iEvent.triggerNames(*trigResults);
      for( int i_Trig = 0; i_Trig < N_Triggers; ++i_Trig ) {
	if (trigResults.product()->accept(i_Trig)) {
	  TString TrigPath =trigName.triggerName(i_Trig);
	    if(TrigPath.Contains("HLT_IsoMu20")) return true;
	    if(TrigPath.Contains("HLT_Ele23_WPLoose_Gsf"))return true;
	    if(TrigPath.Contains("HLT_Ele23_CaloIdL_TrackIdL_IsoVL"))return true;
	    //if(TrigPath.Contains("HLT_Ele27_eta2p1_WPLoose_Gsf"))return true;
	    //if(TrigPath.Contains("HLT_Ele32_eta2p1_WPLoose_Gsf"))return true;
	    //if(TrigPath.Contains("HLT_Ele27_eta2p1_WP75_Gsf"))return true;
	    //if(TrigPath.Contains("HLT_Ele32_eta2p1_WP75_Gsf"))return true;
	    if(TrigPath.Contains("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL"))return true;
            if(TrigPath.Contains("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ"))return true;
            //if(TrigPath.Contains("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL"))return true;
            //if(TrigPath.Contains("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ"))return true;
            if(TrigPath.Contains("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL"))return true;
            if(TrigPath.Contains("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ"))return true;
            if(TrigPath.Contains("HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL")) return true;
            if(TrigPath.Contains("HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL"))return true;


	}
      }
    }
    return false;
  }
  
  

  return true;
} 



bool SSb13_takeimai::PassTriggerLeg(std::string triggerlegstring, std::string triggerlegstringalt,  const pat::Muon *muonit, const edm::Event& iEvent ){

  edm::Handle<edm::TriggerResults> triggerBits;
  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  edm::InputTag triggerBits_("TriggerResults","","HLT");
  edm::InputTag  triggerObjects_("selectedPatTrigger");
  if(isgetbylabel) iEvent.getByLabel(triggerObjects_, triggerObjects);
  else iEvent.getByToken(trigobject_token, triggerObjects);

  if(isgetbylabel)  iEvent.getByLabel(triggerBits_, triggerBits);
  else iEvent.getByToken(trgresults_token, triggerBits);  
  
  const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
  
  for (pat::TriggerObjectStandAlone obj : *triggerObjects) { // note: not "const &" since we want to call unpackPathNames
    obj.unpackPathNames(names);
    for (unsigned h = 0; h < obj.filterLabels().size(); ++h){

      string myfillabl=obj.filterLabels()[h];
      if( (myfillabl.find(triggerlegstring)!=string::npos  || myfillabl.find(triggerlegstringalt)!=string::npos)&& deltaR(muonit->eta(),muonit->phi(), obj.eta(),obj.phi())<0.4 ) return true;
      if( (myfillabl.find(triggerlegstring)!=string::npos  || myfillabl.find(triggerlegstringalt)!=string::npos)&& myfillabl.find("hltL1s")!=string::npos &&deltaR(muonit->eta(),muonit->phi(), obj.eta(),obj.phi())<0.5 ) return true;// For L1, looser dr matching


    }
  }
  return false;
}

bool SSb13_takeimai::PassTriggerLeg(std::string triggerlegstring, std::string triggerlegstringalt,const pat::Electron *eleit, const edm::Event& iEvent ){

  edm::Handle<edm::TriggerResults> triggerBits;
  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  edm::InputTag triggerBits_("TriggerResults","","HLT");
  edm::InputTag  triggerObjects_("selectedPatTrigger");
  if(isgetbylabel) iEvent.getByLabel(triggerObjects_, triggerObjects);
  else iEvent.getByToken(trigobject_token, triggerObjects);

  if(isgetbylabel)  iEvent.getByLabel(triggerBits_, triggerBits);
  else iEvent.getByToken(trgresults_token, triggerBits);

  const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
  
  //   cout << "Filters: "<< endl;
  for (pat::TriggerObjectStandAlone obj : *triggerObjects) { // note: not "const &" since we want to call unpackPathNames
    obj.unpackPathNames(names);
    for (unsigned h = 0; h < obj.filterLabels().size(); ++h){

      string myfillabl=obj.filterLabels()[h];
      //        cout << myfillabl<<endl;
      
      if( (myfillabl.find(triggerlegstring)!=string::npos  || myfillabl.find(triggerlegstringalt)!=string::npos)&& deltaR(eleit->eta(),eleit->phi(), obj.eta(),obj.phi())<0.4 ) return true;
      if( (myfillabl.find(triggerlegstring)!=string::npos  || myfillabl.find(triggerlegstringalt)!=string::npos)&& myfillabl.find("L1sL1")!=string::npos &&deltaR(eleit->eta(),eleit->phi(), obj.eta(),obj.phi())<0.5 ) return true;// For L1, looser dr matching
      
    }
  }
  return false;
}


void SSb13_takeimai::InitLepton()
{

  _closeIndex.push_back(-99);
  _indeces.push_back(-99);
  _flavors.push_back(-99);
  _charges.push_back(-99);
  _chargesMC.push_back(-99);
  _isolation.push_back(-99.);
  _miniisolation_0p2.push_back(-99.);
  _miniisolationcharged_0p2.push_back(-99.);
  _miniisoneutral.push_back(-99);
  _miniisolation_0p3.push_back(-99.);
  _multiisolation_T.push_back(false);
  _multiisolation_M.push_back(false);
  _multiisolation_L.push_back(false);
  _trackSelectionMultiplicity.push_back(-99);
  _muonSegmentComp.push_back(-99.);
  _lmvawithiso.push_back(-99.);
  _lmva.push_back(-99.);
  _pfisocharged.push_back(-99.);
  _pfisoneutral.push_back(-99.);
  _pfisophoton.push_back(-99.);
  _isolationMCraw.push_back(-99.);
  _isolationMCnonu.push_back(-99.);
  _isolationMCdr03.push_back(-99.);
  _isolationMCdr03nonu.push_back(-99.);
  _ptrel.push_back(-99.);
  _ptrel2.push_back(-99.);
  _ptratio.push_back(-99.);
  _ljetmultiplicity.push_back(-99.);
  _mvaValue.push_back(-99.);
  _mt.push_back(-99.);
  _mllZ.push_back(-99.);
  _mllG.push_back(-99.);
  _mll.push_back(-99.);
  _origin.push_back(-99);
  _originReduced.push_back(-99);
  _ipPV.push_back(-99.);
  _ipPVerr.push_back(-99.);
  _ipPVmc.push_back(-99.);
  _ipZPV.push_back(-99.);
  _ipZPVerr.push_back(-99.);
  _3dIP.push_back(-99.);
  _3dIPerr.push_back(-99.);
  _3dIPsig.push_back(-99.);
  _closeJetPtAll.push_back(-99.);
  _closeJetEtaAll.push_back(-99.);
  _closeJetPhiAll.push_back(-99.);
  _closeJetEAll.push_back(-99.);
  _closeJetCSVAll.push_back(-99.);
  _closeJetNconstAll.push_back(-99);
  _closeJetAngAll.push_back(-99.);
  _ptRelAll.push_back(-99.);
  _closeJetPtAllMC.push_back(-99.);
  _closeJetPtAllstatus.push_back(-99.);
  _partonIdMatched.push_back(-99);
  _sameParton.push_back(false);
  _isloose.push_back(false);
  _issoftmuon.push_back(false);
  _issoftmuonwithiso.push_back(false);
  _istight.push_back(false);
  _isMediumMuon.push_back(false);
  _istightIso.push_back(false);
  _istightID.push_back(false);
  _istightIDWP2016.push_back(false);
  _istightIDWP2016_IsoEmul.push_back(false);
  _isFOIDWP2016.push_back(false);
  _isFOIDWP2016_IsoEmul.push_back(false);
  _isVetoIDWP2016_RA7.push_back(false);
  _istightIDWP2016_noIP3D.push_back(false);
  _isFOIDWP2016_noIP3D.push_back(false);
  _istightIDWP2016_RA7.push_back(false);
  _isFOIDWP2016_RA7.push_back(false);
  _isVetoIDWP2016_EWK.push_back(false);
  _isFOIDWP2016_EWK.push_back(false);
  _istightIDWP2016_EWK.push_back(false);
  _islooseMT2.push_back(false);
  _iscutbasedWPtight.push_back(false);
  _istrigemulID.push_back(false);
  _istrigemulISO.push_back(false);
  _mompt.push_back(-99.);
  _momphi.push_back(-99.);
  _mometa.push_back(-99.);
  _mompdg.push_back(-99);
  _lPt.push_back(-99.);
  _lEta.push_back(-99.);
  _lPhi.push_back(-99.);
  _lE.push_back(-99.);


  _lPtmc.push_back(-99.);
  _lEtamc.push_back(-99.);
  _lPhimc.push_back(-99.);
  _lEmc.push_back(-99.);
  _nuPtmc.push_back(-99.);
  _nuEtamc.push_back(-99.);
  _nuPhimc.push_back(-99.);
  _nuEmc.push_back(-99.);
  _mtmc.push_back(-99.);
  _lnmisshits.push_back(-99);
  _passconvveto.push_back(false);
  _lsietaieta.push_back(-99);
  _lfull5x5sietaieta.push_back(-99);
  _lhovere.push_back(-99);

  _ldphiin.push_back(-99);
  _ldetain.push_back(-99);
  _l1oemin1op.push_back(-99);

    _tau_dz.push_back(-99);
    _TauIs_againstMuonLoose3.push_back(false);
    _TauIs_byLooseCombinedIsolationDeltaBetaCorr3Hits.push_back(false);
    _TauIs_decayModeFindingNewDMs.push_back(false);
    _TauIs_decayModeFinding.push_back(false);
    _TauIs_byLooseIsolationMVArun2v1DBdR03oldDMwLT.push_back(false);
    _TauIs_byMediumIsolationMVArun2v1DBdR03oldDMwLT.push_back(false);
    _TauIs_byVLooseIsolationMVArun2v1DBoldDMwLT.push_back(false);
    _TauIs_byTightIsolationMVArun2v1DBoldDMwLT.push_back(false);

  _lchargeGSF.push_back(-99);
  _lchargeCTF.push_back(-99);
  _lchargePixSc.push_back(-99);
  passleg1L.push_back(false);
  hltL1sL1TripleMu553.push_back(false);
  hltL1TripleMu553L2TriMuFiltered3L3TriMuFiltered12105.push_back(false);
  hltL1TripleMu553L2TriMuFiltered3L3TriMuFiltered10105.push_back(false);
  hltL1TripleMu553L2TriMuFiltered3L3TriMuFiltered5.push_back(false);
  hltDiMu9Ele9CaloIdLTrackIdLMuonlegL3Filtered9.push_back(false);
  hltDiMu9Ele9CaloIdLTrackIdLElectronlegDphiFilter.push_back(false);
  hltMu8DiEle12CaloIdLTrackIdLMuonlegL3Filtered8.push_back(false);
  hltMu8DiEle12CaloIdLTrackIdLElectronlegDphiFilter.push_back(false);
  hltL1sL1TripleEG14108.push_back(false);
  hltEle16Ele12Ele8CaloIdLTrackIdLDphiLeg1Filter.push_back(false);
  hltEle16Ele12Ele8CaloIdLTrackIdLDphiLeg2Filter.push_back(false);
  hltEle16Ele12Ele8CaloIdLTrackIdLDphiLeg3Filter.push_back(false);
  hltL1sL1Mu6HTT150ORL1Mu8HTT125ORL1HTT175.push_back(false);
  hltDoubleMu8Mass8L3Filtered.push_back(false);
  hltL1sL1DoubleEG6HTT150orL1HTT175.push_back(false);
  hltDoubleEle8CaloIdMGsfTrackIdMDphiFilter.push_back(false);
  hltMu8Ele8CaloIdMGsfTrackIdMDphiFilter.push_back(false);
  hltMuon8L3Filtered0.push_back(false);
  hltL1sL1DoubleMu103p5ORDoubleMu125.push_back(false);
  hltDiMuonGlb17Glb8RelTrkIsoFiltered0p4.push_back(false);
  hltL3fL1sDoubleMu103p5L1f0L2f10OneMuL3Filtered17.push_back(false);
  hltL3pfL1sDoubleMu103p5L1f0L2pf0L3PreFiltered8.push_back(false);
  hltL3fL1sDoubleMu103p5L1f0L2f10L3Filtered17.push_back(false);
  hltDiMuonGlbFiltered17TrkFiltered8.push_back(false);
  hltDiMuonGlb17Trk8RelTrkIsoFiltered0p4.push_back(false);
  hltDiMuonGlb17Glb8RelTrkIsoFiltered0p4DzFiltered0p2.push_back(false);
  hltDiMuonGlb17Trk8RelTrkIsoFiltered0p4DzFiltered0p2.push_back(false);
  hltL1sL1DoubleEG2210.push_back(false);
  hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter.push_back(false);
  hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter.push_back(false);
  hltEle23Ele12CaloIdLTrackIdLIsoVLDZFilter.push_back(false);
  hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter.push_back(false);
  hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter.push_back(false);
  hltEle17Ele12CaloIdLTrackIdLIsoVLDZFilter.push_back(false);
  hltL1sL1DoubleEG1510.push_back(false);
  hltL1sL1Mu20EG10.push_back(false);
  hltMu23TrkIsoVVLEle8CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered23.push_back(false);
  hltMu23TrkIsoVVLEle8CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter.push_back(false);
  hltL1sL1Mu12EG10.push_back(false);
  hltMu17TrkIsoVVLEle12CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered17.push_back(false);
  hltMu17TrkIsoVVLEle12CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter.push_back(false);

  hltL1sL1Mu5EG20.push_back(false);
  hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered8.push_back(false);
  hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter.push_back(false);

  hltL1sL1Mu5EG15.push_back(false);
  hltMu8TrkIsoVVLEle17CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered8.push_back(false);
  hltMu8TrkIsoVVLEle17CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter.push_back(false);

  hltL3crIsoL1sMu20Eta2p1L1f0L2f10QL3f24QL3trkIsoFiltered0p09.push_back(false);
  hltL3crIsoL1sMu16L1f0L2f10QL3f20QL3trkIsoFiltered0p09.push_back(false);
  hltL3crIsoL1sMu16L1f0L2f10QL3f18QL3trkIsoFiltered0p09.push_back(false);
  hltL3fL1sMu16L1f0Tkf18QL3trkIsoFiltered0p09.push_back(false);
  hltL3crIsoL1sMu20L1f0L2f10QL3f22QL3trkIsoFiltered0p09.push_back(false);
  hltL3fL1sMu20L1f0Tkf22QL3trkIsoFiltered0p09.push_back(false);

  hltEle23WPLooseGsfTrackIsoFilter.push_back(false);
  hltEle27WPLooseGsfTrackIsoFilter.push_back(false);
  hltEle27noerWPLooseGsfTrackIsoFilter.push_back(false);
  hltEle32WPLooseGsfTrackIsoFilter.push_back(false);
  hltEle23CaloIdLTrackIdLIsoVLTrackIsoFilter.push_back(false);
  hltDiEle33CaloIdLNewPixelMatchUnseededFilter.push_back(false);

  hltL1sL1SingleMu16ORSingleMu25.push_back(false);
  hltDiMuonGlb27Trk8DzFiltered0p2.push_back(false);
  hltL3fL1sMu16orMu25L1f0L2f25L3Filtered27.push_back(false);
  hltL3fL1sMu16orMu25L1f0L2f10QL3Filtered30Q.push_back(false);
  hltEle30CaloIdLGsfTrkIdVLDPhiUnseededFilter.push_back(false);
  hltDiEle33CaloIdLGsfTrkIdVLDPhiUnseededFilter.push_back(false);
  hltL3fL1sL1DoubleMu0ETM40lorDoubleMu0ETM55.push_back(false);
  hltEle22WPLooseL1SingleIsoEG20erGsfTrackIsoFilter.push_back(false);
  hltL3crIsoL1sSingleMu18erIorSingleMu20erL1f0L2f10QL3f19QL3trkIsoFiltered0p09.push_back(false);

  hltL3fL1sMu22orMu25orMu20EG15orMu5EG20L1f0L2f10QL3Filtered30Q.push_back(false);
  
  hltL1sSingleMu16.push_back(false);
  hltL1sSingleMu18.push_back(false);
  hltL1sSingleMu20.push_back(false);
  hltL1sSingleMu20erlorSingleMu22lorSingleMu25.push_back(false);
  hltL1sSingleEG40IorSingleIsoEG22erIorSingleIsoEG24erIorSingleIsoEG24IorSingleIsoEG26.push_back(false);





}

void SSb13_takeimai::ClearStuff(){


  _closeIndex.clear();
  _indeces.clear();
  _flavors.clear();
  _charges.clear();
  _chargesMC.clear();
  _isolation.clear();
  _miniisolation_0p2.clear();
  _miniisolationcharged_0p2.clear();
  _miniisoneutral.clear();
  _miniisolation_0p3.clear();
  _multiisolation_T.clear();
  _multiisolation_M.clear();
  _multiisolation_L.clear();
  _trackSelectionMultiplicity.clear();
  _muonSegmentComp.clear();
  _lmvawithiso.clear();
  _lmva.clear();
  _pfisocharged.clear();
  _pfisoneutral.clear();
  _pfisophoton.clear();
  _isolationMCraw.clear();
  _isolationMCnonu.clear();
  _isolationMCdr03.clear();
  _isolationMCdr03nonu.clear();
  _ptrel.clear();
  _ptrel2.clear();
  _ptratio.clear();
  _ljetmultiplicity.clear();
  _mvaValue.clear();
  _mt.clear();
  _mllZ.clear();
  _mllG.clear();
  _mll.clear();
  _origin.clear();
  _originReduced.clear();
  _ipPV.clear();
  _ipPVerr.clear();
  _ipPVmc.clear();
  _ipZPV.clear();
  _ipZPVerr.clear();
  _3dIP.clear();
  _3dIPerr.clear();
  _3dIPsig.clear();
  _closeJetPtAll.clear();
  _closeJetEtaAll.clear();
  _closeJetPhiAll.clear();
  _closeJetEAll.clear();
  _closeJetCSVAll.clear();
  _closeJetNconstAll.clear();
  _closeJetAngAll.clear();
  _ptRelAll.clear();
  _closeJetPtAllMC.clear();
  _closeJetPtAllstatus.clear();
  _partonIdMatched.clear();
  _sameParton.clear();
  _isloose.clear();
  _issoftmuon.clear();
  _issoftmuonwithiso.clear();

  _istight.clear();
  _isMediumMuon.clear();
  _istightIso.clear();
  _istightID.clear();

  _istightIDWP2016.clear();
  _istightIDWP2016_IsoEmul.clear();
  _isFOIDWP2016.clear();
  _isFOIDWP2016_IsoEmul.clear();
  _isVetoIDWP2016_RA7.clear();
  _isFOIDWP2016_noIP3D.clear();
  _istightIDWP2016_noIP3D.clear();
  _istightIDWP2016_RA7.clear();
  _isFOIDWP2016_RA7.clear();
  _isVetoIDWP2016_EWK.clear();
  _isFOIDWP2016_EWK.clear();
  _istightIDWP2016_EWK.clear();
  _islooseMT2.clear();
  _iscutbasedWPtight.clear();

  _istrigemulID.clear();
  _istrigemulISO.clear();

  _mompt.clear();
  _momphi.clear();
  _mometa.clear();
  _mompdg.clear();
  _lPt.clear();
  _lEta.clear();
  _lPhi.clear();
  _lE.clear();
  _lPtmc.clear();
  _lEtamc.clear();
  _lPhimc.clear();
  _lEmc.clear();

  _lgenPt.clear();
  _lgenEta.clear();
  _lgenPhi.clear();
  _lgenpdgId.clear();

    _genWPhi.clear();
    _genWEta.clear();
    _genWPt.clear();
    _genWE.clear();
    _genZPhi.clear();
    _genZEta.clear();
    _genZPt.clear();
    _genZE.clear();
    
  _nuPtmc.clear();
  _nuEtamc.clear();
  _nuPhimc.clear();
  _nuEmc.clear();
  _mtmc.clear();
  _lnmisshits.clear();
  _passconvveto.clear();
  _lsietaieta.clear();
  _lfull5x5sietaieta.clear();
  _lhovere.clear();
  _ldphiin.clear();
  _ldetain.clear();
  _l1oemin1op.clear();
  _lchargeGSF.clear();
  _lchargeCTF.clear();
  _lchargePixSc.clear();
  _jetJECuncty.clear();
  _bTagged.clear();
  _jetEta.clear();
  _jetPhi.clear();
  _jetPt.clear();
  _jetM.clear();
  _csv.clear();
  _jetpFlav.clear();
  _jethFlav.clear();
  _jetbtagSF.clear();
  _jetbtagSF_up.clear();
  _jetbtagSF_down.clear();
  _jetbtagEff.clear();
  _jetDeltaRloose.clear();
  _mllv.clear();
  passleg1L.clear();
    
    _tau_dz.clear();
    _TauIs_againstMuonLoose3.clear();
    _TauIs_byLooseCombinedIsolationDeltaBetaCorr3Hits.clear();
    _TauIs_decayModeFindingNewDMs.clear();
    _TauIs_decayModeFinding.clear();
    _TauIs_byLooseIsolationMVArun2v1DBdR03oldDMwLT.clear();
    _TauIs_byMediumIsolationMVArun2v1DBdR03oldDMwLT.clear();
    _TauIs_byVLooseIsolationMVArun2v1DBoldDMwLT.clear();
    _TauIs_byTightIsolationMVArun2v1DBoldDMwLT.clear();
    
    
  hltL1sL1TripleMu553.clear();
  hltL1TripleMu553L2TriMuFiltered3L3TriMuFiltered12105.clear();
  hltL1TripleMu553L2TriMuFiltered3L3TriMuFiltered10105.clear();
  hltL1TripleMu553L2TriMuFiltered3L3TriMuFiltered5.clear();
  hltDiMu9Ele9CaloIdLTrackIdLMuonlegL3Filtered9.clear();
  hltDiMu9Ele9CaloIdLTrackIdLElectronlegDphiFilter.clear();
  hltMu8DiEle12CaloIdLTrackIdLMuonlegL3Filtered8.clear();
  hltMu8DiEle12CaloIdLTrackIdLElectronlegDphiFilter.clear();
  hltL1sL1TripleEG14108.clear();
  hltEle16Ele12Ele8CaloIdLTrackIdLDphiLeg1Filter.clear();
  hltEle16Ele12Ele8CaloIdLTrackIdLDphiLeg2Filter.clear();
  hltEle16Ele12Ele8CaloIdLTrackIdLDphiLeg3Filter.clear();
  hltL1sL1Mu6HTT150ORL1Mu8HTT125ORL1HTT175.clear();
  hltDoubleMu8Mass8L3Filtered.clear();
  hltL1sL1DoubleEG6HTT150orL1HTT175.clear();
  hltDoubleEle8CaloIdMGsfTrackIdMDphiFilter.clear();
  hltMu8Ele8CaloIdMGsfTrackIdMDphiFilter.clear();
  hltMuon8L3Filtered0.clear();
  hltL1sL1DoubleMu103p5ORDoubleMu125.clear();
  hltDiMuonGlb17Glb8RelTrkIsoFiltered0p4.clear();
  hltL3fL1sDoubleMu103p5L1f0L2f10OneMuL3Filtered17.clear();
  hltL3pfL1sDoubleMu103p5L1f0L2pf0L3PreFiltered8.clear();
  hltL3fL1sDoubleMu103p5L1f0L2f10L3Filtered17.clear();
  hltDiMuonGlbFiltered17TrkFiltered8.clear();
  hltDiMuonGlb17Trk8RelTrkIsoFiltered0p4.clear();
  hltDiMuonGlb17Glb8RelTrkIsoFiltered0p4DzFiltered0p2.clear();
  hltDiMuonGlb17Trk8RelTrkIsoFiltered0p4DzFiltered0p2.clear();
  hltL1sL1DoubleEG2210.clear();
  hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter.clear();
  hltEle23Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter.clear();
  hltEle23Ele12CaloIdLTrackIdLIsoVLDZFilter.clear();
  hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg1Filter.clear();
  hltEle17Ele12CaloIdLTrackIdLIsoVLTrackIsoLeg2Filter.clear();
  hltEle17Ele12CaloIdLTrackIdLIsoVLDZFilter.clear();
  hltL1sL1DoubleEG1510.clear();
  hltL1sL1Mu20EG10.clear();
  hltMu23TrkIsoVVLEle8CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered23.clear();
  hltMu23TrkIsoVVLEle8CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter.clear();
  hltL1sL1Mu12EG10.clear();
  hltMu17TrkIsoVVLEle12CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered17.clear();
  hltMu17TrkIsoVVLEle12CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter.clear();
  hltL1sL1Mu5EG20.clear();
  hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered8.clear();
  hltMu8TrkIsoVVLEle23CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter.clear();
  hltL1sL1Mu5EG15.clear();
  hltMu8TrkIsoVVLEle17CaloIdLTrackIdLIsoVLMuonlegL3IsoFiltered8.clear();
  hltMu8TrkIsoVVLEle17CaloIdLTrackIdLIsoVLElectronlegTrackIsoFilter.clear();

  hltL3crIsoL1sMu20Eta2p1L1f0L2f10QL3f24QL3trkIsoFiltered0p09.clear();
  hltL3crIsoL1sMu16L1f0L2f10QL3f18QL3trkIsoFiltered0p09.clear();
  hltL3fL1sMu16L1f0Tkf18QL3trkIsoFiltered0p09.clear();
  hltL3crIsoL1sMu20L1f0L2f10QL3f22QL3trkIsoFiltered0p09.clear();
  hltL3fL1sMu20L1f0Tkf22QL3trkIsoFiltered0p09.clear();


  hltEle27WPLooseGsfTrackIsoFilter.clear();
  hltEle27noerWPLooseGsfTrackIsoFilter.clear();
  hltEle32WPLooseGsfTrackIsoFilter.clear();
  hltL3crIsoL1sMu16L1f0L2f10QL3f20QL3trkIsoFiltered0p09.clear();
  hltEle23WPLooseGsfTrackIsoFilter.clear();
  hltEle23CaloIdLTrackIdLIsoVLTrackIsoFilter.clear();
  hltDiEle33CaloIdLNewPixelMatchUnseededFilter.clear();

  hltL1sL1SingleMu16ORSingleMu25.clear();
  hltDiMuonGlb27Trk8DzFiltered0p2.clear();
  hltL3fL1sMu16orMu25L1f0L2f25L3Filtered27.clear();
  hltL3fL1sMu16orMu25L1f0L2f10QL3Filtered30Q.clear();
  hltEle30CaloIdLGsfTrkIdVLDPhiUnseededFilter.clear();
  hltDiEle33CaloIdLGsfTrkIdVLDPhiUnseededFilter.clear();
  hltL3fL1sL1DoubleMu0ETM40lorDoubleMu0ETM55.clear();
  hltEle22WPLooseL1SingleIsoEG20erGsfTrackIsoFilter.clear();
  hltL3crIsoL1sSingleMu18erIorSingleMu20erL1f0L2f10QL3f19QL3trkIsoFiltered0p09.clear();

  
  hltL3fL1sMu22orMu25orMu20EG15orMu5EG20L1f0L2f10QL3Filtered30Q.clear();


  hltL1sSingleMu16.clear();
  hltL1sSingleMu18.clear();
  hltL1sSingleMu20.clear();
  hltL1sSingleMu20erlorSingleMu22lorSingleMu25.clear();
  hltL1sSingleEG40IorSingleIsoEG22erIorSingleIsoEG24erIorSingleIsoEG24IorSingleIsoEG26.clear();


}

DEFINE_FWK_MODULE(SSb13_takeimai);
