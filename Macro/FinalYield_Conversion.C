#define FinalYield_cxx
#include "FinalYield.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <stdio.h>
#include "TLorentzVector.h"
#include "TF1.h"
#include "TMath.h"

#include <fstream>

/*
 /Dissertation Ch 13 plots
 /Final results
 */
/*
 *First .L mt2_bisect.cpp
 *Then .L .C
 */
/*
 Uncertainty up/down hists:
 1.pile-up
 2.b-tagging
 3.lepton SF
 */
/*
 Added on uncertainties:
 1. Luminosity
 2. Trigger
 3. Processes
 */
 
using namespace std;

//MC/Data Trees
//Signal regions
/*
 MC:
 WZ
 Nonprompt: TTDilept, TTSingleLeptFromTbar, TTSingleLeptFromT
 ZZ/H: ZZTo4L
 Conversion: ZGamma
 VVV: ZZZ, WZZ, WWZ, WWW
 TTbar+X: TTW(2), TTZToLLNuNu, TTZToLL, ttHtoNonbb
 Total: 15
 */
/*
 const double xsection[nbMCsamples]={
 4.4297,//WZ
 87.315,//TT_Dilept
 182.175,//TT_SingleLeptFromTBar
 182.175,//TT_SingleLeptFromT
 1.256,//ZZ/H
 131.3,//ZGamma
 0.01398,//ZZZ
 0.05565,//WZZ
 0.1651,//WWZ
 0.2086,//WWW
 0.2043,//TTW(1)
 0.2043,//TTW(2)
 0.2529,//TTZToLLNuNu
 0.0493,//TTZToLL
 0.2151//ttHToNonbb
 }
 
 TString samplefilesMC[nbMCsamples] = {
 "WZTo3LNu_skim3Lreco.root",//WZ
 "TTJets_DiLept_skim3Lreco.root",//TT_Dilept
 "TTJets_SingleLeptFromTbar_skim3Lreco.root",//TT_SingleLeptFromTbar
 "TTJets_SingleLeptFromT_skim3Lreco.root",//TT_SingleLeptFromT
 "ZZTo4L_skim3Lreco.root",//ZZ/H
 "ZGTo2LG_skim3Lreco.root",//ZGamma
 "ZZZ_skim3Lreco.root",//ZZZ
 "WZZ_skim3Lreco.root",//WZZ
 "WWZ_skim3Lreco.root",//WWZ
 "WWW_skim3Lreco.root",//WWW
 "TTWJetsToLNu_1_skim3Lreco.root",//TTW(1)
 "TTWJetsToLNu_2_skim3Lreco.root",//TTW(2)
 "TTZToLLNuNu_skim3Lreco.root",//TTZToLLNuNu
 "TTZToLL_skim3Lreco.root",//TTZToLL
 "ttHToNonbb_skim3Lreco.root"//ttHToNonbb
 };
 */

const int nbMCsamples = 2;
const int nSRCat=6;

const int Ptbin[nSRCat] = {20,20,20,20,30,20};
const int METbin[nSRCat] = {10,10,10,10,10,10};
const int MTbin[nSRCat] = {10,10,30,30,30,30};
const int Mllbin[nSRCat] = {28,28,28,28,28,28};
const int M3lbin[nSRCat] = {40,40,40,40,40,40};
const int SumQbin[nSRCat] = {7,7,7,7,7,7};

const double Ptlowerbound = 0;
const double Ptupperbound[nSRCat] = {200,200,200,200,300,200};
const double METlowerbound = 50;
const double METupperbound[nSRCat] = {300,300,300,300,300,300};
const double MTlowerbound = 0;
const double MTupperbound[nSRCat] = {200,200,300,300,300,300};
const double Mlllowerbound = 0;
const double Mllupperbound[nSRCat] = {220,220,220,220,220,220};
const double M3llowerbound = 0;
const double M3lupperbound[nSRCat] = {200,200,200,200,200,200};
const double SumQlowerbound = -3;
const double SumQupperbound[nSRCat] = {4,4,4,4,4,4};

/*
 Const regions
 */
const int nSR[nSRCat]={44,6,18,16,12,12};

const double xsection[nbMCsamples]={
    131.3,//ZGamma
    585.8//WGamma
};

const double intlumi=35900;


TString samplefilesMC[nbMCsamples] = {
    "ZGTo2LG_skim3Lreco.root",//ZGamma
    "WGToLNuG_skim3Lreco.root",//WGamma
};

void FinalYield::Loop()
{
    
    //////////////////////////////
    ////Set up pileup correction//
    //////////////////////////////
    //Read the pile up file
    //Central/Up/Down
    const int nbVtx=50;
    double DataVtxDist[nbVtx]={0};
    double DataVtxDist_Up[nbVtx]={0};
    double DataVtxDist_Down[nbVtx]={0};
    TH1D* pileup = new TH1D("pileup","",nbVtx,0,nbVtx);
    TH1D* pileup_up = new TH1D("pileup_up","",nbVtx,0,nbVtx);
    TH1D* pileup_down = new TH1D("pileup_down","",nbVtx,0,nbVtx);
    
    TFile* DataVtx;
    
    DataVtx=new TFile("DataPileup_20161206.root","open");
    pileup->SetDirectory(0);
    DataVtx->cd();
    pileup->Read("pileup");
    double totalVtxData=pileup->Integral();
    for(int i=0;i<nbVtx;i++){
        DataVtxDist[i]=(pileup->GetBinContent(i+1))/totalVtxData;
    }
    DataVtx->Close();
    
    DataVtx=new TFile("DataPileup_20161206_Up.root","open");
    pileup_up->SetDirectory(0);
    DataVtx->cd();
    pileup_up->Read("pileup");
    double totalVtxData_Up=pileup_up->Integral();
    for(int i=0;i<nbVtx;i++){
        DataVtxDist_Up[i]=(pileup_up->GetBinContent(i+1))/totalVtxData_Up;
    }
    DataVtx->Close();
    
    DataVtx=new TFile("DataPileup_20161206_Down.root","open");
    pileup_down->SetDirectory(0);
    DataVtx->cd();
    pileup_down->Read("pileup");
    double totalVtxData_Down=pileup_down->Integral();
    for(int i=0;i<nbVtx;i++){
        DataVtxDist_Down[i]=(pileup_down->GetBinContent(i+1))/totalVtxData_Down;
    }
    DataVtx->Close();
    
    for(int i=0;i<nbVtx;i++){
        cout<<"Data PU distribution (C/U/D): "<<i<<" "<<DataVtxDist[i]<<" "<<DataVtxDist_Up[i]<<" "<<DataVtxDist_Down[i]<<endl;
    }
    
    //Now read the MC Nvtx distribution histogram
    TFile * MCVtx[nbMCsamples];
    TH1F* MCVtxhist[nbMCsamples];
    double MCVtxDist[nbMCsamples][nbVtx]={0};
    double reweightSF[nbMCsamples][nbVtx]={0};
    double reweightSF_Up[nbMCsamples][nbVtx]={0};
    double reweightSF_Down[nbMCsamples][nbVtx]={0};
    
    double totalVtxMC[nbMCsamples]={0};
    for (int iSam=0; iSam<nbMCsamples; iSam++) {
        MCVtxhist[iSam]=new TH1F("","",nbVtx,0,nbVtx);
        MCVtx[iSam]=new TFile("/cms/data/store/user/t2/users/takeimai/output/SRYield/2016MC/"+samplefilesMC[iSam],"open");
        MCVtx[iSam]->cd("FakeElectrons");
        MCVtxhist[iSam]->SetDirectory(0);
        MCVtxhist[iSam]->Read("N_{vtx}");
        totalVtxMC[iSam]=MCVtxhist[iSam]->Integral();
        for(int i=0;i<nbVtx;i++){
            MCVtxDist[iSam][i]=(MCVtxhist[iSam]->GetBinContent(i+1))/totalVtxMC[iSam];
            if(MCVtxDist[iSam][i]!=0) {
                reweightSF[iSam][i]=DataVtxDist[i]/MCVtxDist[iSam][i];
                reweightSF_Up[iSam][i]=DataVtxDist_Up[i]/MCVtxDist[iSam][i];
                reweightSF_Down[iSam][i]=DataVtxDist_Down[i]/MCVtxDist[iSam][i];
            }
            
            cout<<"PU Reweighting Factors (C/U/D): "<<samplefilesMC[iSam]<<" "<<i<<" "<<reweightSF[iSam][i]<<" "<<reweightSF_Up[iSam][i]<<" "<<reweightSF_Down[iSam][i]<<endl;
        }
    }
    //////////////////////////////
    ////Set up pileup correction//
    //////////////////////////////
    //////////////////////////////
    //////////////////////////////
    
    
    //////////////////////////////
    ////Set up leoton SF//////////
    //////////////////////////////
    //Only Applied to MC
    TFile* SFmap;
    //Electron
    SFmap=new TFile("/cms/data/store/user/t2/users/takeimai/CMSSW_8_0_13/src/SUSYAnalyzer/PatAnalyzer/test/Macro/ScaleFactor/el_SF_FullSim.root","open");
    TH2D* el_SF_FullSim = (TH2D*)SFmap->Get("GsfElectronToLeptonMvaMIDEmuTightIP2DSIP3D8mini04");
    SFmap=new TFile("/cms/data/store/user/t2/users/takeimai/CMSSW_8_0_13/src/SUSYAnalyzer/PatAnalyzer/test/Macro/ScaleFactor/el_RecoSF_FullSim.root","open");
    TH2D* el_RecoSF_FullSim = (TH2D*)SFmap->Get("EGamma_SF2D");
    //Muon
    SFmap=new TFile("/cms/data/store/user/t2/users/takeimai/CMSSW_8_0_13/src/SUSYAnalyzer/PatAnalyzer/test/Macro/ScaleFactor/mu_dxydz_SF_FullSim.root","open");
    TH2D* mu_dxydz_SF_FullSim = (TH2D*)SFmap->Get("SF");
    SFmap=new TFile("/cms/data/store/user/t2/users/takeimai/CMSSW_8_0_13/src/SUSYAnalyzer/PatAnalyzer/test/Macro/ScaleFactor/mu_MediumID_SF_FullSim.root","open");
    TH2D* mu_MediumID_SF_FullSim = (TH2D*)SFmap->Get("SF");
    SFmap=new TFile("/cms/data/store/user/t2/users/takeimai/CMSSW_8_0_13/src/SUSYAnalyzer/PatAnalyzer/test/Macro/ScaleFactor/mu_miniIso_SF_FullSim.root","open");
    TH2D* mu_miniIso_SF_FullSim = (TH2D*)SFmap->Get("SF");
    SFmap=new TFile("/cms/data/store/user/t2/users/takeimai/CMSSW_8_0_13/src/SUSYAnalyzer/PatAnalyzer/test/Macro/ScaleFactor/mu_LepMVAVM_SF_FullSim.root","open");
    TH2D* mu_LepMVAVM_SF_FullSim = (TH2D*)SFmap->Get("SF");
    SFmap=new TFile("/cms/data/store/user/t2/users/takeimai/CMSSW_8_0_13/src/SUSYAnalyzer/PatAnalyzer/test/Macro/ScaleFactor/mu_RecoSF_FullSim.root","open");
    TGraphAsymmErrors* mu_RecoSF_FullSim_Highpt = (TGraphAsymmErrors*)SFmap->Get("ratio_eff_eta3_dr030e030_corr");//Pt>10
    TGraphAsymmErrors* mu_RecoSF_FullSim_Lowpt = (TGraphAsymmErrors*)SFmap->Get("ratio_eff_eta3_tk0_dr030e030_corr");//Pt<10
    //Tau
    //No input for FullSim tau Lep SF
    //////////////////////////////
    ////Set up leoton SF//////////
    //////////////////////////////
    
    //Set up histograms
    TString name;
    TString SRCat[nSRCat]={
        "SR_A",
        "SR_B",
        "SR_C",
        "SR_D",
        "SR_E",
        "SR_F"
    };
    
    TH1F* MCYield_Conversion[nSRCat];
    TH1F* MCYield_Conversion_JECUp[nSRCat];
    TH1F* MCYield_Conversion_JECDown[nSRCat];
    TH1F* MCYield_Conversion_BTagUp[nSRCat];
    TH1F* MCYield_Conversion_BTagDown[nSRCat];
    TH1F* MCYield_Conversion_LepSFUp[nSRCat];
    TH1F* MCYield_Conversion_LepSFDown[nSRCat];
    TH1F* MCYield_Conversion_PUUp[nSRCat];
    TH1F* MCYield_Conversion_PUDown[nSRCat];
    
    TH1F* MCYield_Conversion_Pt1[nSRCat];
    TH1F* MCYield_Conversion_Pt2[nSRCat];
    TH1F* MCYield_Conversion_Pt3[nSRCat];
    TH1F* MCYield_Conversion_Mll[nSRCat];
    TH1F* MCYield_Conversion_MT[nSRCat];
    TH1F* MCYield_Conversion_MET[nSRCat];
    TH1F* MCYield_Conversion_M3l[nSRCat];
    TH1F* MCYield_Conversion_SumQ[nSRCat];
    
    for (int i=0; i<nSRCat; i++) {
        name="MCBG_Conversion_"+SRCat[i];
        MCYield_Conversion[i]=new TH1F(name,"",nSR[i],0.5,nSR[i]+0.5);
        MCYield_Conversion[i]->Sumw2();
        name="MCBG_Conversion_JECUp_"+SRCat[i];
        MCYield_Conversion_JECUp[i]=new TH1F(name,"",nSR[i],0.5,nSR[i]+0.5);
        MCYield_Conversion_JECUp[i]->Sumw2();
        name="MCBG_Conversion_JECDown_"+SRCat[i];
        MCYield_Conversion_JECDown[i]=new TH1F(name,"",nSR[i],0.5,nSR[i]+0.5);
        MCYield_Conversion_JECDown[i]->Sumw2();
        name="MCBG_Conversion_BTagUp_"+SRCat[i];
        MCYield_Conversion_BTagUp[i]=new TH1F(name,"",nSR[i],0.5,nSR[i]+0.5);
        MCYield_Conversion_BTagUp[i]->Sumw2();
        name="MCBG_Conversion_BTagDown_"+SRCat[i];
        MCYield_Conversion_BTagDown[i]=new TH1F(name,"",nSR[i],0.5,nSR[i]+0.5);
        MCYield_Conversion_BTagDown[i]->Sumw2();
        name="MCBG_Conversion_LepSFUp_"+SRCat[i];
        MCYield_Conversion_LepSFUp[i]=new TH1F(name,"",nSR[i],0.5,nSR[i]+0.5);
        MCYield_Conversion_LepSFUp[i]->Sumw2();
        name="MCBG_Conversion_LepSFDown_"+SRCat[i];
        MCYield_Conversion_LepSFDown[i]=new TH1F(name,"",nSR[i],0.5,nSR[i]+0.5);
        MCYield_Conversion_LepSFDown[i]->Sumw2();
        name="MCBG_Conversion_PUUp_"+SRCat[i];
        MCYield_Conversion_PUUp[i]=new TH1F(name,"",nSR[i],0.5,nSR[i]+0.5);
        MCYield_Conversion_PUUp[i]->Sumw2();
        name="MCBG_Conversion_PUDown_"+SRCat[i];
        MCYield_Conversion_PUDown[i]=new TH1F(name,"",nSR[i],0.5,nSR[i]+0.5);
        MCYield_Conversion_PUDown[i]->Sumw2();
        
        name="MCBG_Conversion_Pt1_"+SRCat[i];
        MCYield_Conversion_Pt1[i]=new TH1F(name,"",Ptbin[i],Ptlowerbound,Ptupperbound[i]);
        MCYield_Conversion_Pt1[i]->Sumw2();
        name="MCBG_Conversion_Pt2_"+SRCat[i];
        MCYield_Conversion_Pt2[i]=new TH1F(name,"",Ptbin[i],Ptlowerbound,Ptupperbound[i]);
        MCYield_Conversion_Pt2[i]->Sumw2();
        name="MCBG_Conversion_Pt3_"+SRCat[i];
        MCYield_Conversion_Pt3[i]=new TH1F(name,"",Ptbin[i],Ptlowerbound,Ptupperbound[i]);
        MCYield_Conversion_Pt3[i]->Sumw2();
        name="MCBG_Conversion_Mll_"+SRCat[i];
        MCYield_Conversion_Mll[i]=new TH1F(name,"",Mllbin[i],Mlllowerbound,Mllupperbound[i]);
        MCYield_Conversion_Mll[i]->Sumw2();
        name="MCBG_Conversion_MT_"+SRCat[i];
        MCYield_Conversion_MT[i]=new TH1F(name,"",MTbin[i],MTlowerbound,MTupperbound[i]);
        MCYield_Conversion_MT[i]->Sumw2();
        name="MCBG_Conversion_MET_"+SRCat[i];
        MCYield_Conversion_MET[i]=new TH1F(name,"",METbin[i],METlowerbound,METupperbound[i]);
        MCYield_Conversion_MET[i]->Sumw2();
        name="MCBG_Conversion_M3l_"+SRCat[i];
        MCYield_Conversion_M3l[i]=new TH1F(name,"",M3lbin[i],M3llowerbound,M3lupperbound[i]);
        MCYield_Conversion_M3l[i]->Sumw2();
        name="MCBG_Conversion_SumQ_"+SRCat[i];
        MCYield_Conversion_SumQ[i]=new TH1F(name,"",SumQbin[i],SumQlowerbound,SumQupperbound[i]);
        MCYield_Conversion_SumQ[i]->Sumw2();
    }
    
    ////////////////////////
    //MVA Setup/////////////
    ////////////////////////
    //NO NEED ANYMORE
    //_lmva matches the calculation here
    /*
     TMVA::Reader *readerEle = new TMVA::Reader( "!Color:!Silent" );
     TMVA::Reader *readerMu = new TMVA::Reader( "!Color:!Silent" );
     
     Float_t LepGood_pt, LepGood_eta, LepGood_jetNDauChargedMVASel,
     LepGood_miniRelIsoCharged, LepGood_miniRelIsoNeutral,
     LepGood_jetPtRelv2, LepGood_jetPtRatio,
     LepGood_jetBTagCSV,
     LepGood_sip3d, LepGood_dxy, LepGood_dz,
     LepGood_segmentCompatibility,
     LepGood_mvaIdSpring16GP;
     
     readerEle->AddVariable( "LepGood_pt", &LepGood_pt );
     readerEle->AddVariable( "LepGood_eta", &LepGood_eta );
     readerEle->AddVariable( "LepGood_jetNDauChargedMVASel", &LepGood_jetNDauChargedMVASel );
     readerEle->AddVariable( "LepGood_miniRelIsoCharged", &LepGood_miniRelIsoCharged );
     readerEle->AddVariable( "LepGood_miniRelIsoNeutral", &LepGood_miniRelIsoNeutral );
     readerEle->AddVariable( "LepGood_jetPtRelv2", &LepGood_jetPtRelv2 );
     readerEle->AddVariable( "min(LepGood_jetPtRatiov2,1.5)", &LepGood_jetPtRatio );
     readerEle->AddVariable( "max(LepGood_jetBTagCSV,0)", &LepGood_jetBTagCSV );
     readerEle->AddVariable( "LepGood_sip3d", &LepGood_sip3d );
     readerEle->AddVariable( "log(abs(LepGood_dxy))", &LepGood_dxy );
     readerEle->AddVariable( "log(abs(LepGood_dz))", &LepGood_dz );
     readerEle->AddVariable( "LepGood_mvaIdSpring16GP", &LepGood_mvaIdSpring16GP );
     
     readerEle->BookMVA( "BDTG method", "el_BDTG.weights.xml" );
     
     readerMu->AddVariable( "LepGood_pt", &LepGood_pt );
     readerMu->AddVariable( "LepGood_eta", &LepGood_eta );
     readerMu->AddVariable( "LepGood_jetNDauChargedMVASel", &LepGood_jetNDauChargedMVASel );
     readerMu->AddVariable( "LepGood_miniRelIsoCharged", &LepGood_miniRelIsoCharged );
     readerMu->AddVariable( "LepGood_miniRelIsoNeutral", &LepGood_miniRelIsoNeutral );
     readerMu->AddVariable( "LepGood_jetPtRelv2", &LepGood_jetPtRelv2 );
     readerMu->AddVariable( "min(LepGood_jetPtRatiov2,1.5)", &LepGood_jetPtRatio );
     readerMu->AddVariable( "max(LepGood_jetBTagCSV,0)", &LepGood_jetBTagCSV );
     readerMu->AddVariable( "LepGood_sip3d", &LepGood_sip3d );
     readerMu->AddVariable( "log(abs(LepGood_dxy))", &LepGood_dxy );
     readerMu->AddVariable( "log(abs(LepGood_dz))", &LepGood_dz );
     readerMu->AddVariable( "LepGood_segmentCompatibility", &LepGood_segmentCompatibility );
     
     readerMu->BookMVA( "BDTG method", "mu_BDTG.weights.xml" );
     */
    ////////////////////////
    //MVA Setup/////////////
    ////////////////////////
    
    
    ///////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////Run on the samples /////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////
    //Now loop over MC samples
    TFile * inputMC[nbMCsamples];
    
    for(int iSam = 0;  iSam <nbMCsamples; iSam++){//loop over all MC samples
        
        inputMC[iSam] = new TFile("/cms/data/store/user/t2/users/takeimai/output/SRYield/2016MC/"+samplefilesMC[iSam],"open");
        inputMC[iSam]->cd("FakeElectrons");
        TH1F* eventcounter=new TH1F("hCounter","",5,0,5);
        eventcounter->SetDirectory(0);
        eventcounter->Read("hCounter");
        int nbevents=eventcounter->GetEntries();//total nb of events
        inputMC[iSam]->Close();
        
        inputMC[iSam] = new TFile("/cms/data/store/user/t2/users/takeimai/output/SRYield/2016MC/"+samplefilesMC[iSam],"open");
        TTree *thetree_mc = (TTree*)(inputMC[iSam])->Get("FakeElectrons/fakeTree");
        Init(thetree_mc);
        Long64_t nentries_mc = (*thetree_mc).GetEntries();//Get the number of generated events
        
        cout<<"Processing MC sample "<<samplefilesMC[iSam]<<endl;
        cout<<"Total number of generated events: "<<nbevents<<endl;
        
        for (Long64_t jentry=0; jentry<nentries_mc;jentry++) {//loop over all the events//First loop
            
            LoadTree(jentry);
            thetree_mc->GetEntry(jentry);
            cout<<jentry<<endl;
            
            ////////////////////
            //SELECTION BEGINS//
            ////////////////////
            //Set up B-Tag and MET flag
            bool bJetandMETflag=false;
            bool bJetandMETflag_JECUp=false;
            bool bJetandMETflag_JECDown=false;
            
            ////////////////////////////////
            //Calculate the weight to fill//
            ////////////////////////////////
            //The variable _weight need to be normalized for some MC samples
            if(_weight>0) {
                _weight=1;
            }
            else{
                _weight=-1;
            }
            double W = _weight*intlumi*xsection[iSam]/nbevents*reweightSF[iSam][trueNVtx];
            double W_PUUp = _weight*intlumi*xsection[iSam]/nbevents*reweightSF_Up[iSam][trueNVtx];
            double W_PUDown = _weight*intlumi*xsection[iSam]/nbevents*reweightSF_Down[iSam][trueNVtx];
            
            /////////////////////////////////////////////////////////
            //Apply trigger selection////////////////////////////////
            /////////////////////////////////////////////////////////
            if(!TriggerSelectionMC(samplefilesMC[iSam])) continue;
            
            
            /////////////////////////////////////////////////////////
            //Apply immediate vetos//////////////////////////////////
            /////////////////////////////////////////////////////////
            //Met filter
            if (!passmetfilters) continue;
            //b-Jet veto
            //MET CUT
            bJetandMETflag=bJetandMETcut(_met,0);
            bJetandMETflag_JECUp=bJetandMETcut(_met_JECup,1);
            bJetandMETflag_JECDown=bJetandMETcut(_met_JECdown,-1);
            if (bJetandMETflag&&bJetandMETflag_JECUp&&bJetandMETflag_JECDown) continue;
            /////////////////////////////////////////////////////////
            //Apply immediate vetos//////////////////////////////////
            /////////////////////////////////////////////////////////
            
            
            /////////////////////////////////////////////////////////
            ///////////////SELECT 3 LEPTONS//////////////////////////
            /////////////////////////////////////////////////////////
            /*
             6 SR TYPES
             A: 3l OSSF
             B: 3l NoOSSF
             C: 2l OSSF + 1tau
             D: 2l OS-e/mu pair + 1tau
             E: 2l SSSF + 1tau
             F: 1l + 2tau
             */
            
            //>=3 FO leptons, only the hardest three pass the Tight selection
            int chosenlept[3] = { 9999,9999,9999 };
            double chosenPt[3] = { 0 };
            
            int QualifiedLepton=0;
            vector<int> LeptonID;
            vector<double> LeptonPt;
            LeptonID.clear();
            LeptonPt.clear();
            
            for (int lept = 0; lept <(*_lPt).size(); lept++){//loops over all the leptons//First count light leptons
                if (!FOSelection(lept)) continue;//FOID for both light leptons and taus
                int lptsize=(*_lPt).size();
                if (!DeltaRSelection(lept,lptsize)) continue;
                LeptonID.push_back(lept);
                LeptonPt.push_back((*_lPt)[lept]);
                QualifiedLepton++;
            }
            
            //Skip events with <3 FO leptons
            if (QualifiedLepton<3) continue;
            
            //Choose three leptons with HIGHEST Pt
            for (int lept1 = 0; lept1 <LeptonID.size(); lept1++){//loops over all the selected leptons
                if (((*_lPt)[LeptonID[lept1]]>chosenPt[0]) && (LeptonID[lept1] != chosenlept[1]) && (LeptonID[lept1] != chosenlept[2])) {
                    chosenlept[0] = LeptonID[lept1];
                    chosenPt[0] = (*_lPt)[LeptonID[lept1]];
                }
                
                for (int lept2 = 0; lept2 <LeptonID.size(); lept2++){
                    if (LeptonID[lept2] == LeptonID[lept1]) continue;
                    if (((*_lPt)[LeptonID[lept2]]>chosenPt[1]) && (LeptonID[lept2] != chosenlept[0]) && (LeptonID[lept2] != chosenlept[2])) {
                        chosenlept[1] = LeptonID[lept2];
                        chosenPt[1] = (*_lPt)[LeptonID[lept2]];
                    }
                    
                    for (int lept3 = 0; lept3 <LeptonID.size(); lept3++){
                        if (LeptonID[lept3] == LeptonID[lept2] || lept3 == lept1) continue;
                        if (((*_lPt)[LeptonID[lept3]]>chosenPt[2]) && (LeptonID[lept3] != chosenlept[0]) && (LeptonID[lept3] != chosenlept[1])) {
                            chosenlept[2] = LeptonID[lept3];
                            chosenPt[2] = (*_lPt)[LeptonID[lept3]];
                        }
                    }
                }
            }
            
            //////////////////////////////////////
            //The hardest three FO leptons should pass the Tight selection
            //None of the other FO leptons should pass the Tight selection
            bool eventselection=true;
            for (int lept=0; lept<LeptonID.size(); lept++) {
                if ((LeptonID[lept]==chosenlept[0]||LeptonID[lept]==chosenlept[1]||LeptonID[lept]==chosenlept[2])&&!ObjectSelection(LeptonID[lept])) {
                    eventselection=false;
                    continue;
                }
                if ((LeptonID[lept]!=chosenlept[0]&&LeptonID[lept]!=chosenlept[1]&&LeptonID[lept]!=chosenlept[2])&&ObjectSelection(LeptonID[lept])) {
                    eventselection=false;
                    continue;
                }
            }
            if (!eventselection) continue;
            //////////////////////////////////////
            
            /*
             //Identify tau leptons
             //Calculate SumQ
             */
            int nTau=0;
            double SumQ=0;
            for (int i=0; i<3; i++) {
                SumQ+=(*_charges)[chosenlept[i]];
                if((*_flavors)[chosenlept[i]]==2) nTau++;
            }
            if (nTau>2) continue;
            
            /*
             *Lepton/B-Tagging SF factor
             */
            ////////////////////////////////////////////////////////////////////////////////////
            //Incorporate Lepton SF
            double TotalLeptonSF=1;
            double TotalLeptonSF_LepSFUp=1;
            double TotalLeptonSF_LepSFDown=1;
            for (int i=0; i<3; i++) {
                double LeptonSF=1;
                double LeptonSF_LepSFUp=1;
                double LeptonSF_LepSFDown=1;
                //Ele
                if((*_flavors)[chosenlept[i]]==0) {
                    LeptonSF=CalcLepSF_ele_FullSim(el_SF_FullSim,el_RecoSF_FullSim,chosenlept[i],0);
                    LeptonSF_LepSFUp=CalcLepSF_ele_FullSim(el_SF_FullSim,el_RecoSF_FullSim,chosenlept[i],1);
                    LeptonSF_LepSFDown=CalcLepSF_ele_FullSim(el_SF_FullSim,el_RecoSF_FullSim,chosenlept[i],-1);
                }
                //Mu
                else if((*_flavors)[chosenlept[i]]==1) {
                    if ((*_lPt)[chosenlept[i]]<10) {//Low Pt Muon
                        LeptonSF=CalcLepSF_mu_FullSim(mu_dxydz_SF_FullSim,mu_MediumID_SF_FullSim,mu_miniIso_SF_FullSim,mu_LepMVAVM_SF_FullSim,mu_RecoSF_FullSim_Lowpt,chosenlept[i],0);
                        LeptonSF_LepSFUp=CalcLepSF_mu_FullSim(mu_dxydz_SF_FullSim,mu_MediumID_SF_FullSim,mu_miniIso_SF_FullSim,mu_LepMVAVM_SF_FullSim,mu_RecoSF_FullSim_Lowpt,chosenlept[i],1);
                        LeptonSF_LepSFDown=CalcLepSF_mu_FullSim(mu_dxydz_SF_FullSim,mu_MediumID_SF_FullSim,mu_miniIso_SF_FullSim,mu_LepMVAVM_SF_FullSim,mu_RecoSF_FullSim_Lowpt,chosenlept[i],-1);
                    }
                    else {//High Pt Muon
                        LeptonSF=CalcLepSF_mu_FullSim(mu_dxydz_SF_FullSim,mu_MediumID_SF_FullSim,mu_miniIso_SF_FullSim,mu_LepMVAVM_SF_FullSim,mu_RecoSF_FullSim_Highpt,chosenlept[i],0);
                        LeptonSF_LepSFUp=CalcLepSF_mu_FullSim(mu_dxydz_SF_FullSim,mu_MediumID_SF_FullSim,mu_miniIso_SF_FullSim,mu_LepMVAVM_SF_FullSim,mu_RecoSF_FullSim_Highpt,chosenlept[i],1);
                        LeptonSF_LepSFDown=CalcLepSF_mu_FullSim(mu_dxydz_SF_FullSim,mu_MediumID_SF_FullSim,mu_miniIso_SF_FullSim,mu_LepMVAVM_SF_FullSim,mu_RecoSF_FullSim_Highpt,chosenlept[i],-1);
                    }
                }
                //Tau
                else if((*_flavors)[chosenlept[i]]==2) {
                    LeptonSF=CalcLepSF_tau_FullSim(0);
                    LeptonSF_LepSFUp=CalcLepSF_tau_FullSim(1);
                    LeptonSF_LepSFDown=CalcLepSF_tau_FullSim(-1);
                }
                
                if (LeptonSF==0||TMath::IsNaN(LeptonSF)) LeptonSF=1;
                if (LeptonSF_LepSFUp==0||TMath::IsNaN(LeptonSF_LepSFUp)) LeptonSF_LepSFUp=1;
                if (LeptonSF_LepSFDown==0||TMath::IsNaN(LeptonSF_LepSFDown)) LeptonSF_LepSFDown=1;
                TotalLeptonSF=TotalLeptonSF*LeptonSF;
                TotalLeptonSF_LepSFUp=TotalLeptonSF_LepSFUp*LeptonSF_LepSFUp;
                TotalLeptonSF_LepSFDown=TotalLeptonSF_LepSFDown*LeptonSF_LepSFDown;
            }
            ////////////////////////////////////////////////////////////////////////////////////
            //Incorporate B-Tagging SF
            double TotalBTagSF=1;
            double TotalBTagSF_BTagUp=1;
            double TotalBTagSF_BTagDown=1;
            TotalBTagSF=CalcBTagSF(0);
            TotalBTagSF_BTagUp=CalcBTagSF(1);
            TotalBTagSF_BTagDown=CalcBTagSF(-1);
            if (TotalBTagSF==0||TMath::IsNaN(TotalBTagSF)) TotalBTagSF=1;
            if (TotalBTagSF_BTagUp==0||TMath::IsNaN(TotalBTagSF_BTagUp)) TotalBTagSF_BTagUp=1;
            if (TotalBTagSF_BTagDown==0||TMath::IsNaN(TotalBTagSF_BTagDown)) TotalBTagSF_BTagDown=1;
            
            double W_BTagUp=W*TotalLeptonSF*TotalBTagSF_BTagUp;
            double W_BTagDown=W*TotalLeptonSF*TotalBTagSF_BTagDown;
            double W_LepSFUp=W*TotalLeptonSF_LepSFUp*TotalBTagSF;
            double W_LepSFDown=W*TotalLeptonSF_LepSFDown*TotalBTagSF;
            W=W*TotalLeptonSF*TotalBTagSF;
            W_PUUp=W_PUUp*TotalLeptonSF*TotalBTagSF;
            W_PUDown=W_PUDown*TotalLeptonSF*TotalBTagSF;
            
            
            
            /*
             *Categorize events
             */
            ////////////////////////////////////////////////////////////////////////////////////
            /*
             *0Tau
             */
            if (nTau==0) {
                //Cat A(OSSF) and B(NoOSSF)
                //Event Requirements
                /*
                 1. 3 Tight leptons
                 2. Pt > 20(25), 10(15), 10 for mu(ele)
                 3. If in an event with a leading muon other leptons are electrons or taus, then the leading muon is required to have Pt > 25 GeV
                 4. OSSF: |TriInvMass-Z|<15
                 5. MET>50
                 */
                
                ////Pt cut on 3 leptons////
                double leadingPt=0;
                double subleadingPt=0;
                double trailingPt=0;
                int leadinglep=99;
                int subleadinglep=99;
                int trailinglep=99;
                
                for (int i=0; i<3; i++) {
                    if(chosenPt[i]>leadingPt){
                        trailingPt=subleadingPt;
                        trailinglep=subleadinglep;
                        subleadingPt=leadingPt;
                        subleadinglep=leadinglep;
                        leadingPt=chosenPt[i];
                        leadinglep=chosenlept[i];
                    }
                    else if(chosenPt[i]>subleadingPt && chosenPt[i]<=leadingPt){
                        trailingPt=subleadingPt;
                        trailinglep=subleadinglep;
                        subleadingPt=chosenPt[i];
                        subleadinglep=chosenlept[i];
                    }
                    else if(chosenPt[i]>trailingPt && chosenPt[i]<=subleadingPt){
                        trailingPt=chosenPt[i];
                        trailinglep=chosenlept[i];
                    }
                }
                
                if ((*_flavors)[leadinglep]==1&&leadingPt<20) continue;
                if ((*_flavors)[leadinglep]==0&&leadingPt<25) continue;
                if ((*_flavors)[subleadinglep]==1&&subleadingPt<10) continue;
                if ((*_flavors)[subleadinglep]==0&&subleadingPt<15) continue;
                if (trailingPt<10) continue;
                if ((*_flavors)[leadinglep]==1&&(*_flavors)[subleadinglep]!=1&&(*_flavors)[trailinglep]!=1&&leadingPt<25) continue;
                
                
                
                ///OSSF/NoOSSF///
                bool CatA = false, CatB = false;
                bool OSpairflag = false;
                double InvMassPair=0;
                double DistPeak=9999;
                double Mll=0;
                double Mt=0;
                double Mt_JECUp=0;
                double Mt_JECDown=0;
                int SRIndex=9999;
                
                if (((*_charges)[chosenlept[0]] == (*_charges)[chosenlept[1]]) && ((*_charges)[chosenlept[0]] == (*_charges)[chosenlept[2]])){
                    //For 3 SS leptons, Mt takes the smallest value and Mll=0, event goes to CatB
                    Mt=9999;
                    Mt_JECUp=9999;
                    Mt_JECDown=9999;
                    Mll=0;
                    for (int i=0; i<3; i++) {
                        double Mt_temp=CalcTransMass(chosenlept[i],_met,_met_phi);
                        if (Mt_temp<Mt) Mt=Mt_temp;
                        double Mt_temp_JECUp=CalcTransMass(chosenlept[i],_met_JECup,_met_phi_JECup);
                        if (Mt_temp_JECUp<Mt_JECUp) Mt_JECUp=Mt_temp_JECUp;
                        double Mt_temp_JECDown=CalcTransMass(chosenlept[i],_met_JECdown,_met_phi_JECdown);
                        if (Mt_temp_JECDown<Mt_JECDown) Mt_JECDown=Mt_temp_JECDown;
                    }
                    CatB=true;
                }
                else{
                    //Three leptons with the same flavor
                    if(((*_flavors)[chosenlept[0]]==(*_flavors)[chosenlept[1]])&&((*_flavors)[chosenlept[0]]==(*_flavors)[chosenlept[2]])) {
                        for (int i=0; i<3; i++) {
                            for (int j=0; j<3; j++) {
                                if(i==j) continue;
                                if(i>j) continue;
                                if((*_charges)[chosenlept[i]]==(*_charges)[chosenlept[j]]) continue;
                                InvMassPair=CalcInvariantMass(chosenlept[i],chosenlept[j]);
                                if (abs(InvMassPair-91)<DistPeak) {
                                    //Choose the lepton with inv. mass closest to Z mass as the Mll variable
                                    DistPeak=abs(InvMassPair-91);
                                    Mll=InvMassPair;
                                    for (int k=0; k<3; k++) {
                                        if((k==i)||(k==j)) continue;
                                        Mt=CalcTransMass(chosenlept[k],_met,_met_phi);
                                        Mt_JECUp=CalcTransMass(chosenlept[k],_met_JECup,_met_phi_JECup);
                                        Mt_JECDown=CalcTransMass(chosenlept[k],_met_JECdown,_met_phi_JECdown);
                                    }
                                }
                            }
                        }
                        CatA=true;
                        OSpairflag=true;
                    }
                    //Three leptons with different flavors (uue/uee)
                    else {
                        //aab type
                        if ((*_flavors)[chosenlept[0]]==(*_flavors)[chosenlept[1]]) {
                            if ((*_charges)[chosenlept[0]]==(*_charges)[chosenlept[1]]) {//aa have the same sign
                                for (int i=0; i<3; i++) {
                                    if(i==2) continue;
                                    InvMassPair=CalcInvariantMass(chosenlept[i],chosenlept[2]);
                                    if (abs(InvMassPair-50)<DistPeak) {
                                        //Choose the lepton with inv. mass closest to Z mass as the Mll variable
                                        DistPeak=abs(InvMassPair-50);
                                        Mll=InvMassPair;
                                        for (int k=0; k<3; k++) {
                                            if((k==i)||(k==2)) continue;
                                            Mt=CalcTransMass(chosenlept[k],_met,_met_phi);
                                            Mt_JECUp=CalcTransMass(chosenlept[k],_met_JECup,_met_phi_JECup);
                                            Mt_JECDown=CalcTransMass(chosenlept[k],_met_JECdown,_met_phi_JECdown);
                                        }
                                    }
                                }
                                CatB=true;
                                OSpairflag=true;
                            }
                            if ((*_charges)[chosenlept[0]]!=(*_charges)[chosenlept[1]]) {//aa have opposite sign
                                Mll=CalcInvariantMass(chosenlept[0],chosenlept[1]);
                                Mt=CalcTransMass(chosenlept[2],_met,_met_phi);
                                Mt_JECUp=CalcTransMass(chosenlept[2],_met_JECup,_met_phi_JECup);
                                Mt_JECDown=CalcTransMass(chosenlept[2],_met_JECdown,_met_phi_JECdown);
                                CatA=true;
                                OSpairflag=true;
                            }
                        }
                        //aba type
                        else if((*_flavors)[chosenlept[0]]==(*_flavors)[chosenlept[2]]){
                            if ((*_charges)[chosenlept[0]]==(*_charges)[chosenlept[2]]) {//aa have the same sign
                                for (int i=0; i<3; i++) {
                                    if(i==1) continue;
                                    InvMassPair=CalcInvariantMass(chosenlept[i],chosenlept[1]);
                                    if (abs(InvMassPair-50)<DistPeak) {
                                        //Choose the lepton with inv. mass closest to Z mass as the Mll variable
                                        DistPeak=abs(InvMassPair-50);
                                        Mll=InvMassPair;
                                        for (int k=0; k<3; k++) {
                                            if((k==i)||(k==1)) continue;
                                            Mt=CalcTransMass(chosenlept[k],_met,_met_phi);
                                            Mt_JECUp=CalcTransMass(chosenlept[k],_met_JECup,_met_phi_JECup);
                                            Mt_JECDown=CalcTransMass(chosenlept[k],_met_JECdown,_met_phi_JECdown);
                                        }
                                    }
                                }
                                CatB=true;
                                OSpairflag=true;
                            }
                            if ((*_charges)[chosenlept[0]]!=(*_charges)[chosenlept[2]]) {//aa have opposite sign
                                Mll=CalcInvariantMass(chosenlept[0],chosenlept[2]);
                                Mt=CalcTransMass(chosenlept[1],_met,_met_phi);
                                Mt_JECUp=CalcTransMass(chosenlept[1],_met_JECup,_met_phi_JECup);
                                Mt_JECDown=CalcTransMass(chosenlept[1],_met_JECdown,_met_phi_JECdown);
                                CatA=true;
                                OSpairflag=true;
                            }
                        }
                        //baa type
                        else if((*_flavors)[chosenlept[1]]==(*_flavors)[chosenlept[2]]){
                            if ((*_charges)[chosenlept[1]]==(*_charges)[chosenlept[2]]) {//aa have the same sign
                                for (int i=0; i<3; i++) {
                                    if(i==0) continue;
                                    InvMassPair=CalcInvariantMass(chosenlept[i],chosenlept[0]);
                                    if (abs(InvMassPair-50)<DistPeak) {
                                        //Choose the lepton with inv. mass closest to Z mass as the Mll variable
                                        DistPeak=abs(InvMassPair-50);
                                        Mll=InvMassPair;
                                        for (int k=0; k<3; k++) {
                                            if((k==i)||(k==0)) continue;
                                            Mt=CalcTransMass(chosenlept[k],_met,_met_phi);
                                            Mt_JECUp=CalcTransMass(chosenlept[k],_met_JECup,_met_phi_JECup);
                                            Mt_JECDown=CalcTransMass(chosenlept[k],_met_JECdown,_met_phi_JECdown);
                                        }
                                    }
                                }
                                CatB=true;
                                OSpairflag=true;
                            }
                            if ((*_charges)[chosenlept[1]]!=(*_charges)[chosenlept[2]]) {//aa have the opposite sign
                                Mll=CalcInvariantMass(chosenlept[1],chosenlept[2]);
                                Mt=CalcTransMass(chosenlept[0],_met,_met_phi);
                                Mt_JECUp=CalcTransMass(chosenlept[0],_met_JECup,_met_phi_JECup);
                                Mt_JECDown=CalcTransMass(chosenlept[0],_met_JECdown,_met_phi_JECdown);
                                CatA=true;
                                OSpairflag=true;
                            }
                        }
                    }
                }//Categorize the event by flavor/charge combination
                
                /*
                 *
                 */
                ////////////////////////////////////////////////////////////////
                ///////Fill The Hist////////////////////////////////////////////
                ////////////////////////////////////////////////////////////////
                /*
                 *
                 */
                ///////Exclude Events With Low Mll////////
                if((!CatA)&&(!CatB)) continue;
                if (Mll<12&&OSpairflag) continue;
                
                ///////OSSF: |TriInvMass-Z|<15/////////
                double M3l=CalcTrilepInvariantMass(chosenlept[0],chosenlept[1],chosenlept[2]);
                if(CatA){
                    if(abs(M3l-91)<15) continue;
                }
                
                
                //JECUp
                if (!bJetandMETflag_JECUp) {
                    if (CatA) {
                        SRIndex=CalcSRIndex_A(_met_JECup,Mll,Mt_JECUp);
                        MCYield_Conversion_JECUp[0]->Fill(SRIndex,W);
                    }
                    if (CatB) {
                        SRIndex=CalcSRIndex_B(_met_JECup,Mll,Mt_JECUp);
                        MCYield_Conversion_JECUp[1]->Fill(SRIndex,W);
                    }
                }
                //JECDown
                if (!bJetandMETflag_JECDown) {
                    if (CatA) {
                        SRIndex=CalcSRIndex_A(_met_JECdown,Mll,Mt_JECDown);
                        MCYield_Conversion_JECDown[0]->Fill(SRIndex,W);
                    }
                    if (CatB) {
                        SRIndex=CalcSRIndex_B(_met_JECdown,Mll,Mt_JECDown);
                        MCYield_Conversion_JECDown[1]->Fill(SRIndex,W);
                    }
                }
                //NO JEC
                if (!bJetandMETflag) {
                    ///////////////////////////////////////////////
                    //Fill SR Hist Before Overflow Bin Adjustment//
                    ///////////////////////////////////////////////
                    if (CatA) {
                        SRIndex=CalcSRIndex_A(_met,Mll,Mt);
                        MCYield_Conversion[0]->Fill(SRIndex,W);
                        MCYield_Conversion_BTagUp[0]->Fill(SRIndex,W_BTagUp);
                        MCYield_Conversion_BTagDown[0]->Fill(SRIndex,W_BTagDown);
                        MCYield_Conversion_LepSFUp[0]->Fill(SRIndex,W_LepSFUp);
                        MCYield_Conversion_LepSFDown[0]->Fill(SRIndex,W_LepSFDown);
                        MCYield_Conversion_PUUp[0]->Fill(SRIndex,W_PUUp);
                        MCYield_Conversion_PUDown[0]->Fill(SRIndex,W_PUDown);
                    }
                    if (CatB) {
                        SRIndex=CalcSRIndex_B(_met,Mll,Mt);
                        MCYield_Conversion[1]->Fill(SRIndex,W);
                        MCYield_Conversion_BTagUp[1]->Fill(SRIndex,W_BTagUp);
                        MCYield_Conversion_BTagDown[1]->Fill(SRIndex,W_BTagDown);
                        MCYield_Conversion_LepSFUp[1]->Fill(SRIndex,W_LepSFUp);
                        MCYield_Conversion_LepSFDown[1]->Fill(SRIndex,W_LepSFDown);
                        MCYield_Conversion_PUUp[1]->Fill(SRIndex,W_PUUp);
                        MCYield_Conversion_PUDown[1]->Fill(SRIndex,W_PUDown);
                    }
                    ///////////////////////////////////////////////
                    //Fill SR Hist Before Overflow Bin Adjustment//
                    ///////////////////////////////////////////////
                    ////////////////////
                    ///Overflow Bin/////
                    ////////////////////
                    int CatNum = 9999;
                    if (CatA) {
                        CatNum=0;
                    }
                    if (CatB) {
                        CatNum=1;
                    }
                    
                    double leadingPt_FL=leadingPt;
                    double subleadingPt_FL=subleadingPt;
                    double trailingPt_FL=trailingPt;
                    double _met_FL=_met;
                    double Mt_FL=Mt;
                    double Mll_FL=Mll;
                    double M3l_FL=M3l;
                    if (leadingPt_FL>Ptupperbound[CatNum]) leadingPt_FL=Ptupperbound[CatNum]-1;
                    if (subleadingPt_FL>Ptupperbound[CatNum]) subleadingPt_FL=Ptupperbound[CatNum]-1;
                    if (trailingPt_FL>Ptupperbound[CatNum]) trailingPt_FL=Ptupperbound[CatNum]-1;
                    if (_met_FL>METupperbound[CatNum]) _met_FL=METupperbound[CatNum]-1;
                    if (Mt_FL>MTupperbound[CatNum]) Mt_FL=MTupperbound[CatNum]-1;
                    if (Mll_FL>Mllupperbound[CatNum]) Mll_FL=Mllupperbound[CatNum]-1;
                    if (M3l_FL>M3lupperbound[CatNum]) M3l_FL=M3lupperbound[CatNum]-1;
                    
                    //////FILL HIST////////
                    if (CatA) {
                        MCYield_Conversion_Pt1[0]->Fill(leadingPt_FL,W);
                        MCYield_Conversion_Pt2[0]->Fill(subleadingPt_FL,W);
                        MCYield_Conversion_Pt3[0]->Fill(trailingPt_FL,W);
                        MCYield_Conversion_Mll[0]->Fill(Mll_FL,W);
                        MCYield_Conversion_MT[0]->Fill(Mt_FL,W);
                        MCYield_Conversion_MET[0]->Fill(_met_FL,W);
                        MCYield_Conversion_M3l[0]->Fill(M3l_FL,W);
                        MCYield_Conversion_SumQ[0]->Fill(SumQ,W);
                    }
                    if (CatB) {
                        MCYield_Conversion_Pt1[1]->Fill(leadingPt_FL,W);
                        MCYield_Conversion_Pt2[1]->Fill(subleadingPt_FL,W);
                        MCYield_Conversion_Pt3[1]->Fill(trailingPt_FL,W);
                        MCYield_Conversion_Mll[1]->Fill(Mll_FL,W);
                        MCYield_Conversion_MT[1]->Fill(Mt_FL,W);
                        MCYield_Conversion_MET[1]->Fill(_met_FL,W);
                        MCYield_Conversion_M3l[1]->Fill(M3l_FL,W);
                        MCYield_Conversion_SumQ[1]->Fill(SumQ,W);
                    }
                    
                }//NO JEC
                
                
            }//if(nTau==0)
            /*
             *0Tau
             */
            ////////////////////////////////////////////////////////////////////////////////////
            
            ////////////////////////////////////////////////////////////////////////////////////
            /*
             *1Tau
             */
            if (nTau==1) {
                //Cat C-E
                //Event Requirements
                /*
                 1. 3 Tight leptons
                 2. Pt > 20(25), 10(15), 10 for mu(ele), Pt >20 for tau
                 3. If in an event with a leading muon other leptons are electrons or taus, then the leading muon is required to have Pt > 25 GeV
                 4. OSSF: |TriInvMass-Z|<15
                 5. MET>50
                 */
                
                ////Pt cut on 3 leptons////
                double leadingPt=0;
                double subleadingPt=0;
                double trailingPt=0;
                int leadinglep=99;
                int subleadinglep=99;
                int trailinglep=99;
                
                for (int i=0; i<3; i++) {
                    if(chosenPt[i]>leadingPt){
                        trailingPt=subleadingPt;
                        trailinglep=subleadinglep;
                        subleadingPt=leadingPt;
                        subleadinglep=leadinglep;
                        leadingPt=chosenPt[i];
                        leadinglep=chosenlept[i];
                    }
                    else if(chosenPt[i]>subleadingPt && chosenPt[i]<=leadingPt){
                        trailingPt=subleadingPt;
                        trailinglep=subleadinglep;
                        subleadingPt=chosenPt[i];
                        subleadinglep=chosenlept[i];
                    }
                    else if(chosenPt[i]>trailingPt && chosenPt[i]<=subleadingPt){
                        trailingPt=chosenPt[i];
                        trailinglep=chosenlept[i];
                    }
                }
                
                if ((*_flavors)[leadinglep]==1&&leadingPt<20) continue;
                if ((*_flavors)[leadinglep]==0&&leadingPt<25) continue;
                if ((*_flavors)[subleadinglep]==1&&subleadingPt<10) continue;
                if ((*_flavors)[subleadinglep]==0&&subleadingPt<15) continue;
                if (trailingPt<10) continue;
                if ((*_flavors)[leadinglep]==1&&(*_flavors)[subleadinglep]!=1&&(*_flavors)[trailinglep]!=1&&leadingPt<25) continue;
                
                
                /*
                 C: 1tau + OSSF
                 D: 1tau + OS e/mu pair
                 E: 1tau + SS pair
                 
                 Mt2 is used instead of Mt: mt2_bisect.cpp
                 C: OSSF+1tau: both light leptons
                 D: OSOF+1tau: the eμ pair
                 E: SS+tau: the leading light lepton and the tau (could be of same sign)
                 */
                bool CatC=false, CatD=false, CatE=false;
                bool OSpairflag=false;
                double InvMassPair=0;
                double DistPeak=9999;
                double Mll=0;
                double Mt2=0;
                double Mt2_JECUp=0;
                double Mt2_JECDown=0;
                int SRIndex=9999;
                
                double ParticleMass[3]={0.00051, 0.105, 1.776};//Mass of Ele/Mu/Tau
                TLorentzVector l1;
                TLorentzVector l2;
                TLorentzVector Emiss;
                TLorentzVector Emiss_JECUp;
                TLorentzVector Emiss_JECDown;
                
                //tau+l+l////////////
                if((*_flavors)[chosenlept[0]]==2){//tau+l+l
                    if ((*_charges)[chosenlept[1]]==(*_charges)[chosenlept[2]]) {//E::2l SS
                        CatE=true;
                        //Calculate Mll
                        if ((*_charges)[chosenlept[0]]==(*_charges)[chosenlept[1]]) {//3 SS Leptons
                            Mll=0;
                        }
                        else {//Tau-Lep pair closest to 60
                            for (int i=0; i<3; i++) {
                                if(i==0) continue;
                                InvMassPair=CalcInvariantMass(chosenlept[i],chosenlept[0]);
                                if (abs(InvMassPair-60)<DistPeak) {
                                    //Choose the lepton with inv. mass closest to Z mass as the Mll variable
                                    DistPeak=abs(InvMassPair-60);
                                    Mll=InvMassPair;
                                }
                            }
                            OSpairflag=true;
                        }
                        //Calculate MT2
                        if (chosenPt[1]>=chosenPt[2]) {//Find Leading Light Lepton
                            l1.SetPtEtaPhiE((*_lPt)[chosenlept[1]],0, (*_lPhi)[chosenlept[1]],(*_lPt)[chosenlept[1]]);
                            l2.SetPtEtaPhiE((*_lPt)[chosenlept[0]],0, (*_lPhi)[chosenlept[0]],(*_lPt)[chosenlept[0]]);
                            Emiss.SetPtEtaPhiE(_met,0, _met_phi,_met);
                            Emiss_JECUp.SetPtEtaPhiE(_met_JECup,0, _met_phi_JECup,_met_JECup);
                            Emiss_JECDown.SetPtEtaPhiE(_met_JECdown,0, _met_phi_JECdown,_met_JECdown);
                            double pa[3]={ParticleMass[(*_flavors)[chosenlept[1]]],l1.Px(),l1.Py()};
                            double pb[3]={ParticleMass[(*_flavors)[chosenlept[0]]],l2.Px(),l2.Py()};
                            double pmiss[3]={0,Emiss.Px(),Emiss.Py()};
                            double pmiss_JECUp[3]={0,Emiss_JECUp.Px(),Emiss_JECUp.Py()};
                            double pmiss_JECDown[3]={0,Emiss_JECDown.Px(),Emiss_JECDown.Py()};
                            Mt2=CalcMt2(pa,pb,pmiss);
                            Mt2_JECUp=CalcMt2(pa,pb,pmiss_JECUp);
                            Mt2_JECDown=CalcMt2(pa,pb,pmiss_JECDown);
                        }
                        else {//Find Leading Light Lepton
                            l1.SetPtEtaPhiE((*_lPt)[chosenlept[2]],0, (*_lPhi)[chosenlept[2]],(*_lPt)[chosenlept[2]]);
                            l2.SetPtEtaPhiE((*_lPt)[chosenlept[0]],0, (*_lPhi)[chosenlept[0]],(*_lPt)[chosenlept[0]]);
                            Emiss.SetPtEtaPhiE(_met,0, _met_phi,_met);
                            Emiss_JECUp.SetPtEtaPhiE(_met_JECup,0, _met_phi_JECup,_met_JECup);
                            Emiss_JECDown.SetPtEtaPhiE(_met_JECdown,0, _met_phi_JECdown,_met_JECdown);
                            double pa[3]={ParticleMass[(*_flavors)[chosenlept[2]]],l1.Px(),l1.Py()};
                            double pb[3]={ParticleMass[(*_flavors)[chosenlept[0]]],l2.Px(),l2.Py()};
                            double pmiss[3]={0,Emiss.Px(),Emiss.Py()};
                            double pmiss_JECUp[3]={0,Emiss_JECUp.Px(),Emiss_JECUp.Py()};
                            double pmiss_JECDown[3]={0,Emiss_JECDown.Px(),Emiss_JECDown.Py()};
                            Mt2=CalcMt2(pa,pb,pmiss);
                            Mt2_JECUp=CalcMt2(pa,pb,pmiss_JECUp);
                            Mt2_JECDown=CalcMt2(pa,pb,pmiss_JECDown);
                        }
                    }
                    else {//C/D::2l OS
                        if ((*_flavors)[chosenlept[1]]==(*_flavors)[chosenlept[2]]) {//C
                            CatC=true;
                            OSpairflag=true;
                        }
                        if ((*_flavors)[chosenlept[1]]!=(*_flavors)[chosenlept[2]]) {//D
                            CatD=true;
                            OSpairflag=true;
                        }
                        //C/D: Same Way to Calculate Mt2/Mll
                        //Calculate Mll
                        Mll=CalcInvariantMass(chosenlept[1],chosenlept[2]);
                        //Calculate Mt2
                        l1.SetPtEtaPhiE((*_lPt)[chosenlept[1]],0, (*_lPhi)[chosenlept[1]],(*_lPt)[chosenlept[1]]);
                        l2.SetPtEtaPhiE((*_lPt)[chosenlept[2]],0, (*_lPhi)[chosenlept[2]],(*_lPt)[chosenlept[2]]);
                        Emiss.SetPtEtaPhiE(_met,0, _met_phi,_met);
                        Emiss_JECUp.SetPtEtaPhiE(_met_JECup,0, _met_phi_JECup,_met_JECup);
                        Emiss_JECDown.SetPtEtaPhiE(_met_JECdown,0, _met_phi_JECdown,_met_JECdown);
                        double pa[3]={ParticleMass[(*_flavors)[chosenlept[1]]],l1.Px(),l1.Py()};
                        double pb[3]={ParticleMass[(*_flavors)[chosenlept[2]]],l2.Px(),l2.Py()};
                        double pmiss[3]={0,Emiss.Px(),Emiss.Py()};
                        double pmiss_JECUp[3]={0,Emiss_JECUp.Px(),Emiss_JECUp.Py()};
                        double pmiss_JECDown[3]={0,Emiss_JECDown.Px(),Emiss_JECDown.Py()};
                        Mt2=CalcMt2(pa,pb,pmiss);
                        Mt2_JECUp=CalcMt2(pa,pb,pmiss_JECUp);
                        Mt2_JECDown=CalcMt2(pa,pb,pmiss_JECDown);
                    }
                }
                //tau+l+l////////////
                
                //l+tau+l////////////
                if((*_flavors)[chosenlept[1]]==2){//l+tau+l
                    if ((*_charges)[chosenlept[0]]==(*_charges)[chosenlept[2]]) {//E::2l SS
                        CatE=true;
                        //Calculate Mll
                        if ((*_charges)[chosenlept[1]]==(*_charges)[chosenlept[0]]) {//3 SS Leptons
                            Mll=0;
                        }
                        else {//Tau-Lep pair closest to 60
                            for (int i=0; i<3; i++) {
                                if(i==1) continue;
                                InvMassPair=CalcInvariantMass(chosenlept[i],chosenlept[1]);
                                if (abs(InvMassPair-60)<DistPeak) {
                                    //Choose the lepton with inv. mass closest to Z mass as the Mll variable
                                    DistPeak=abs(InvMassPair-60);
                                    Mll=InvMassPair;
                                }
                            }
                            OSpairflag=true;
                        }
                        //Calculate MT2
                        if (chosenPt[0]>=chosenPt[2]) {//Find Leading Light Lepton
                            l1.SetPtEtaPhiE((*_lPt)[chosenlept[0]],0, (*_lPhi)[chosenlept[0]],(*_lPt)[chosenlept[0]]);
                            l2.SetPtEtaPhiE((*_lPt)[chosenlept[1]],0, (*_lPhi)[chosenlept[1]],(*_lPt)[chosenlept[1]]);
                            Emiss.SetPtEtaPhiE(_met,0, _met_phi,_met);
                            Emiss_JECUp.SetPtEtaPhiE(_met_JECup,0, _met_phi_JECup,_met_JECup);
                            Emiss_JECDown.SetPtEtaPhiE(_met_JECdown,0, _met_phi_JECdown,_met_JECdown);
                            double pa[3]={ParticleMass[(*_flavors)[chosenlept[0]]],l1.Px(),l1.Py()};
                            double pb[3]={ParticleMass[(*_flavors)[chosenlept[1]]],l2.Px(),l2.Py()};
                            double pmiss[3]={0,Emiss.Px(),Emiss.Py()};
                            double pmiss_JECUp[3]={0,Emiss_JECUp.Px(),Emiss_JECUp.Py()};
                            double pmiss_JECDown[3]={0,Emiss_JECDown.Px(),Emiss_JECDown.Py()};
                            Mt2=CalcMt2(pa,pb,pmiss);
                            Mt2_JECUp=CalcMt2(pa,pb,pmiss_JECUp);
                            Mt2_JECDown=CalcMt2(pa,pb,pmiss_JECDown);
                        }
                        else {//Find Leading Light Lepton
                            l1.SetPtEtaPhiE((*_lPt)[chosenlept[2]],0, (*_lPhi)[chosenlept[2]],(*_lPt)[chosenlept[2]]);
                            l2.SetPtEtaPhiE((*_lPt)[chosenlept[1]],0, (*_lPhi)[chosenlept[1]],(*_lPt)[chosenlept[1]]);
                            Emiss.SetPtEtaPhiE(_met,0, _met_phi,_met);
                            Emiss_JECUp.SetPtEtaPhiE(_met_JECup,0, _met_phi_JECup,_met_JECup);
                            Emiss_JECDown.SetPtEtaPhiE(_met_JECdown,0, _met_phi_JECdown,_met_JECdown);
                            double pa[3]={ParticleMass[(*_flavors)[chosenlept[2]]],l1.Px(),l1.Py()};
                            double pb[3]={ParticleMass[(*_flavors)[chosenlept[1]]],l2.Px(),l2.Py()};
                            double pmiss[3]={0,Emiss.Px(),Emiss.Py()};
                            double pmiss_JECUp[3]={0,Emiss_JECUp.Px(),Emiss_JECUp.Py()};
                            double pmiss_JECDown[3]={0,Emiss_JECDown.Px(),Emiss_JECDown.Py()};
                            Mt2=CalcMt2(pa,pb,pmiss);
                            Mt2_JECUp=CalcMt2(pa,pb,pmiss_JECUp);
                            Mt2_JECDown=CalcMt2(pa,pb,pmiss_JECDown);
                        }
                    }
                    else {//C/D::2l OS
                        if ((*_flavors)[chosenlept[0]]==(*_flavors)[chosenlept[2]]) {//C
                            CatC=true;
                            OSpairflag=true;
                        }
                        if ((*_flavors)[chosenlept[0]]!=(*_flavors)[chosenlept[2]]) {//D
                            CatD=true;
                            OSpairflag=true;
                        }
                        //C/D: Same Way to Calculate Mt2/Mll
                        //Calculate Mll
                        Mll=CalcInvariantMass(chosenlept[0],chosenlept[2]);
                        //Calculate Mt2
                        l1.SetPtEtaPhiE((*_lPt)[chosenlept[0]],0, (*_lPhi)[chosenlept[0]],(*_lPt)[chosenlept[0]]);
                        l2.SetPtEtaPhiE((*_lPt)[chosenlept[2]],0, (*_lPhi)[chosenlept[2]],(*_lPt)[chosenlept[2]]);
                        Emiss.SetPtEtaPhiE(_met,0, _met_phi,_met);
                        Emiss_JECUp.SetPtEtaPhiE(_met_JECup,0, _met_phi_JECup,_met_JECup);
                        Emiss_JECDown.SetPtEtaPhiE(_met_JECdown,0, _met_phi_JECdown,_met_JECdown);
                        double pa[3]={ParticleMass[(*_flavors)[chosenlept[0]]],l1.Px(),l1.Py()};
                        double pb[3]={ParticleMass[(*_flavors)[chosenlept[2]]],l2.Px(),l2.Py()};
                        double pmiss[3]={0,Emiss.Px(),Emiss.Py()};
                        double pmiss_JECUp[3]={0,Emiss_JECUp.Px(),Emiss_JECUp.Py()};
                        double pmiss_JECDown[3]={0,Emiss_JECDown.Px(),Emiss_JECDown.Py()};
                        Mt2=CalcMt2(pa,pb,pmiss);
                        Mt2_JECUp=CalcMt2(pa,pb,pmiss_JECUp);
                        Mt2_JECDown=CalcMt2(pa,pb,pmiss_JECDown);
                    }
                }
                //l+tau+l////////////
                
                //l+l+tau////////////
                if((*_flavors)[chosenlept[2]]==2){//l+l+tau
                    if ((*_charges)[chosenlept[0]]==(*_charges)[chosenlept[1]]) {//E::2l SS
                        CatE=true;
                        //Calculate Mll
                        if ((*_charges)[chosenlept[2]]==(*_charges)[chosenlept[0]]) {//3 SS Leptons
                            Mll=0;
                        }
                        else {//Tau-Lep pair closest to 60
                            for (int i=0; i<3; i++) {
                                if(i==2) continue;
                                InvMassPair=CalcInvariantMass(chosenlept[i],chosenlept[2]);
                                if (abs(InvMassPair-60)<DistPeak) {
                                    //Choose the lepton with inv. mass closest to Z mass as the Mll variable
                                    DistPeak=abs(InvMassPair-60);
                                    Mll=InvMassPair;
                                }
                            }
                            OSpairflag=true;
                        }
                        //Calculate MT2
                        if (chosenPt[0]>=chosenPt[1]) {//Find Leading Light Lepton
                            l1.SetPtEtaPhiE((*_lPt)[chosenlept[0]],0, (*_lPhi)[chosenlept[0]],(*_lPt)[chosenlept[0]]);
                            l2.SetPtEtaPhiE((*_lPt)[chosenlept[2]],0, (*_lPhi)[chosenlept[2]],(*_lPt)[chosenlept[2]]);
                            Emiss.SetPtEtaPhiE(_met,0, _met_phi,_met);
                            Emiss_JECUp.SetPtEtaPhiE(_met_JECup,0, _met_phi_JECup,_met_JECup);
                            Emiss_JECDown.SetPtEtaPhiE(_met_JECdown,0, _met_phi_JECdown,_met_JECdown);
                            double pa[3]={ParticleMass[(*_flavors)[chosenlept[0]]],l1.Px(),l1.Py()};
                            double pb[3]={ParticleMass[(*_flavors)[chosenlept[2]]],l2.Px(),l2.Py()};
                            double pmiss[3]={0,Emiss.Px(),Emiss.Py()};
                            double pmiss_JECUp[3]={0,Emiss_JECUp.Px(),Emiss_JECUp.Py()};
                            double pmiss_JECDown[3]={0,Emiss_JECDown.Px(),Emiss_JECDown.Py()};
                            Mt2=CalcMt2(pa,pb,pmiss);
                            Mt2_JECUp=CalcMt2(pa,pb,pmiss_JECUp);
                            Mt2_JECDown=CalcMt2(pa,pb,pmiss_JECDown);
                        }
                        else {//Find Leading Light Lepton
                            l1.SetPtEtaPhiE((*_lPt)[chosenlept[1]],0, (*_lPhi)[chosenlept[1]],(*_lPt)[chosenlept[1]]);
                            l2.SetPtEtaPhiE((*_lPt)[chosenlept[2]],0, (*_lPhi)[chosenlept[2]],(*_lPt)[chosenlept[2]]);
                            Emiss.SetPtEtaPhiE(_met,0, _met_phi,_met);
                            Emiss_JECUp.SetPtEtaPhiE(_met_JECup,0, _met_phi_JECup,_met_JECup);
                            Emiss_JECDown.SetPtEtaPhiE(_met_JECdown,0, _met_phi_JECdown,_met_JECdown);
                            double pa[3]={ParticleMass[(*_flavors)[chosenlept[1]]],l1.Px(),l1.Py()};
                            double pb[3]={ParticleMass[(*_flavors)[chosenlept[2]]],l2.Px(),l2.Py()};
                            double pmiss[3]={0,Emiss.Px(),Emiss.Py()};
                            double pmiss_JECUp[3]={0,Emiss_JECUp.Px(),Emiss_JECUp.Py()};
                            double pmiss_JECDown[3]={0,Emiss_JECDown.Px(),Emiss_JECDown.Py()};
                            Mt2=CalcMt2(pa,pb,pmiss);
                            Mt2_JECUp=CalcMt2(pa,pb,pmiss_JECUp);
                            Mt2_JECDown=CalcMt2(pa,pb,pmiss_JECDown);
                        }
                    }
                    else {//C/D::2l OS
                        if ((*_flavors)[chosenlept[0]]==(*_flavors)[chosenlept[1]]) {//C
                            CatC=true;
                            OSpairflag=true;
                        }
                        if ((*_flavors)[chosenlept[0]]!=(*_flavors)[chosenlept[1]]) {//D
                            CatD=true;
                            OSpairflag=true;
                        }
                        //C/D: Same Way to Calculate Mt2/Mll
                        //Calculate Mll
                        Mll=CalcInvariantMass(chosenlept[0],chosenlept[1]);
                        //Calculate Mt2
                        l1.SetPtEtaPhiE((*_lPt)[chosenlept[0]],0, (*_lPhi)[chosenlept[0]],(*_lPt)[chosenlept[0]]);
                        l2.SetPtEtaPhiE((*_lPt)[chosenlept[1]],0, (*_lPhi)[chosenlept[1]],(*_lPt)[chosenlept[1]]);
                        Emiss.SetPtEtaPhiE(_met,0, _met_phi,_met);
                        Emiss_JECUp.SetPtEtaPhiE(_met_JECup,0, _met_phi_JECup,_met_JECup);
                        Emiss_JECDown.SetPtEtaPhiE(_met_JECdown,0, _met_phi_JECdown,_met_JECdown);
                        double pa[3]={ParticleMass[(*_flavors)[chosenlept[0]]],l1.Px(),l1.Py()};
                        double pb[3]={ParticleMass[(*_flavors)[chosenlept[1]]],l2.Px(),l2.Py()};
                        double pmiss[3]={0,Emiss.Px(),Emiss.Py()};
                        double pmiss_JECUp[3]={0,Emiss_JECUp.Px(),Emiss_JECUp.Py()};
                        double pmiss_JECDown[3]={0,Emiss_JECDown.Px(),Emiss_JECDown.Py()};
                        Mt2=CalcMt2(pa,pb,pmiss);
                        Mt2_JECUp=CalcMt2(pa,pb,pmiss_JECUp);
                        Mt2_JECDown=CalcMt2(pa,pb,pmiss_JECDown);
                    }
                }
                //l+l+tau////////////
                
                /*
                 *
                 */
                ////////////////////////////////////////////////////////////////
                ///////Fill The Hist////////////////////////////////////////////
                ////////////////////////////////////////////////////////////////
                /*
                 *
                 */
                
                //////Exclude Low Mll//////////////////
                if ((!CatC)&&(!CatD)&&(!CatE)) continue;
                if (Mll<12&&OSpairflag) continue;
                
                ///////OSSF: |TriInvMass-Z|<15/////////
                double M3l=CalcTrilepInvariantMass(chosenlept[0],chosenlept[1],chosenlept[2]);
                if(CatC){
                    if(abs(M3l-91)<15) continue;
                }
                
               
                //JECUp
                if (!bJetandMETflag_JECUp) {
                    if (CatC) {
                        SRIndex=CalcSRIndex_C(_met_JECup,Mll,Mt2_JECUp);
                        MCYield_Conversion_JECUp[2]->Fill(SRIndex,W);
                    }
                    if (CatD) {
                        SRIndex=CalcSRIndex_D(_met_JECup,Mll,Mt2_JECUp);
                        MCYield_Conversion_JECUp[3]->Fill(SRIndex,W);
                    }
                    if (CatE) {
                        SRIndex=CalcSRIndex_E(_met_JECup,Mll,Mt2_JECUp);
                        MCYield_Conversion_JECUp[4]->Fill(SRIndex,W);
                    }
                }
                //JECDown
                if (!bJetandMETflag_JECDown) {
                    if (CatC) {
                        SRIndex=CalcSRIndex_C(_met_JECdown,Mll,Mt2_JECDown);
                        MCYield_Conversion_JECDown[2]->Fill(SRIndex,W);
                    }
                    if (CatD) {
                        SRIndex=CalcSRIndex_D(_met_JECdown,Mll,Mt2_JECDown);
                        MCYield_Conversion_JECDown[3]->Fill(SRIndex,W);
                    }
                    if (CatE) {
                        SRIndex=CalcSRIndex_E(_met_JECdown,Mll,Mt2_JECDown);
                        MCYield_Conversion_JECDown[4]->Fill(SRIndex,W);
                    }
                }
                //NO JEC
                if (!bJetandMETflag) {
                    ///////////////////////////////////////////////
                    //Fill SR Hist Before Overflow Bin Adjustment//
                    ///////////////////////////////////////////////
                    if (CatC) {
                        SRIndex=CalcSRIndex_C(_met,Mll,Mt2);
                        MCYield_Conversion[2]->Fill(SRIndex,W);
                        MCYield_Conversion_BTagUp[2]->Fill(SRIndex,W_BTagUp);
                        MCYield_Conversion_BTagDown[2]->Fill(SRIndex,W_BTagDown);
                        MCYield_Conversion_LepSFUp[2]->Fill(SRIndex,W_LepSFUp);
                        MCYield_Conversion_LepSFDown[2]->Fill(SRIndex,W_LepSFDown);
                        MCYield_Conversion_PUUp[2]->Fill(SRIndex,W_PUUp);
                        MCYield_Conversion_PUDown[2]->Fill(SRIndex,W_PUDown);
                    }
                    if (CatD) {
                        SRIndex=CalcSRIndex_D(_met,Mll,Mt2);
                        MCYield_Conversion[3]->Fill(SRIndex,W);
                        MCYield_Conversion_BTagUp[3]->Fill(SRIndex,W_BTagUp);
                        MCYield_Conversion_BTagDown[3]->Fill(SRIndex,W_BTagDown);
                        MCYield_Conversion_LepSFUp[3]->Fill(SRIndex,W_LepSFUp);
                        MCYield_Conversion_LepSFDown[3]->Fill(SRIndex,W_LepSFDown);
                        MCYield_Conversion_PUUp[3]->Fill(SRIndex,W_PUUp);
                        MCYield_Conversion_PUDown[3]->Fill(SRIndex,W_PUDown);
                    }
                    if (CatE) {
                        SRIndex=CalcSRIndex_E(_met,Mll,Mt2);
                        MCYield_Conversion[4]->Fill(SRIndex,W);
                        MCYield_Conversion_BTagUp[4]->Fill(SRIndex,W_BTagUp);
                        MCYield_Conversion_BTagDown[4]->Fill(SRIndex,W_BTagDown);
                        MCYield_Conversion_LepSFUp[4]->Fill(SRIndex,W_LepSFUp);
                        MCYield_Conversion_LepSFDown[4]->Fill(SRIndex,W_LepSFDown);
                        MCYield_Conversion_PUUp[4]->Fill(SRIndex,W_PUUp);
                        MCYield_Conversion_PUDown[4]->Fill(SRIndex,W_PUDown);
                    }
                    ///////////////////////////////////////////////
                    //Fill SR Hist Before Overflow Bin Adjustment//
                    ///////////////////////////////////////////////
                    ////////////////////
                    ///Overflow Bin/////
                    ////////////////////
                    int CatNum = 9999;
                    if (CatC) {
                        CatNum=2;
                    }
                    if (CatD) {
                        CatNum=3;
                    }
                    if (CatE) {
                        CatNum=4;
                    }
                    
                    double leadingPt_FL=leadingPt;
                    double subleadingPt_FL=subleadingPt;
                    double trailingPt_FL=trailingPt;
                    double _met_FL=_met;
                    double Mt2_FL=Mt2;
                    double Mll_FL=Mll;
                    double M3l_FL=M3l;
                    if (leadingPt_FL>Ptupperbound[CatNum]) leadingPt_FL=Ptupperbound[CatNum]-1;
                    if (subleadingPt_FL>Ptupperbound[CatNum]) subleadingPt_FL=Ptupperbound[CatNum]-1;
                    if (trailingPt_FL>Ptupperbound[CatNum]) trailingPt_FL=Ptupperbound[CatNum]-1;
                    if (_met_FL>METupperbound[CatNum]) _met_FL=METupperbound[CatNum]-1;
                    if (Mt2_FL>MTupperbound[CatNum]) Mt2_FL=MTupperbound[CatNum]-1;
                    if (Mll_FL>Mllupperbound[CatNum]) Mll_FL=Mllupperbound[CatNum]-1;
                    if (M3l_FL>M3lupperbound[CatNum]) M3l_FL=M3lupperbound[CatNum]-1;
                    
                    if (CatC) {
                        MCYield_Conversion_Pt1[2]->Fill(leadingPt_FL,W);
                        MCYield_Conversion_Pt2[2]->Fill(subleadingPt_FL,W);
                        MCYield_Conversion_Pt3[2]->Fill(trailingPt_FL,W);
                        MCYield_Conversion_Mll[2]->Fill(Mll_FL,W);
                        MCYield_Conversion_MT[2]->Fill(Mt2_FL,W);
                        MCYield_Conversion_MET[2]->Fill(_met_FL,W);
                        MCYield_Conversion_M3l[2]->Fill(M3l_FL,W);
                        MCYield_Conversion_SumQ[2]->Fill(SumQ,W);
                    }
                    if (CatD) {
                        MCYield_Conversion_Pt1[3]->Fill(leadingPt_FL,W);
                        MCYield_Conversion_Pt2[3]->Fill(subleadingPt_FL,W);
                        MCYield_Conversion_Pt3[3]->Fill(trailingPt_FL,W);
                        MCYield_Conversion_Mll[3]->Fill(Mll_FL,W);
                        MCYield_Conversion_MT[3]->Fill(Mt2_FL,W);
                        MCYield_Conversion_MET[3]->Fill(_met_FL,W);
                        MCYield_Conversion_M3l[3]->Fill(M3l_FL,W);
                        MCYield_Conversion_SumQ[3]->Fill(SumQ,W);
                    }
                    if (CatE) {
                        MCYield_Conversion_Pt1[4]->Fill(leadingPt_FL,W);
                        MCYield_Conversion_Pt2[4]->Fill(subleadingPt_FL,W);
                        MCYield_Conversion_Pt3[4]->Fill(trailingPt_FL,W);
                        MCYield_Conversion_Mll[4]->Fill(Mll_FL,W);
                        MCYield_Conversion_MT[4]->Fill(Mt2_FL,W);
                        MCYield_Conversion_MET[4]->Fill(_met_FL,W);
                        MCYield_Conversion_M3l[4]->Fill(M3l_FL,W);
                        MCYield_Conversion_SumQ[4]->Fill(SumQ,W);
                    }
                }//NO JEC
                
            }//if(nTau==1)
            /*
             *1Tau
             */
            ////////////////////////////////////////////////////////////////////////////////////
            
            ////////////////////////////////////////////////////////////////////////////////////
            /*
             *2Tau
             */
            if (nTau==2) {
                //Cat F
                //Event Requirements
                /*
                 1. 3 Tight leptons
                 2. Pt > 20 for 2 taus and Pt > 25(30) for the mu(ele)
                 3. All Lepton eta<2.1
                 4. MET>50
                 */
                
                ////Pt cut on 3 leptons////
                //Calculate Leading/Subleading/Trailing Pt
                 double leadingPt=0;
                 double subleadingPt=0;
                 double trailingPt=0;
                 int leadinglep=99;
                 int subleadinglep=99;
                 int trailinglep=99;
                 
                 for (int i=0; i<3; i++) {
                     if(chosenPt[i]>leadingPt){
                         trailingPt=subleadingPt;
                         trailinglep=subleadinglep;
                         subleadingPt=leadingPt;
                         subleadinglep=leadinglep;
                         leadingPt=chosenPt[i];
                         leadinglep=chosenlept[i];
                     }
                     else if(chosenPt[i]>subleadingPt && chosenPt[i]<=leadingPt){
                         trailingPt=subleadingPt;
                         trailinglep=subleadinglep;
                         subleadingPt=chosenPt[i];
                         subleadinglep=chosenlept[i];
                     }
                     else if(chosenPt[i]>trailingPt && chosenPt[i]<=subleadingPt){
                         trailingPt=chosenPt[i];
                         trailinglep=chosenlept[i];
                     }
                 }
                
                bool TakeEvent=true;
                for (int i=0; i<3; i++) {
                    if ((*_flavors)[chosenlept[i]]==2&&(abs((*_lEta)[chosenlept[i]])>=2.1||chosenPt[i]<20)) {
                        TakeEvent=false;
                        break;
                    }
                    if ((*_flavors)[chosenlept[i]]==1&&(abs((*_lEta)[chosenlept[i]])>=2.1||chosenPt[i]<25)) {
                        TakeEvent=false;
                        break;
                    }
                    if ((*_flavors)[chosenlept[i]]==0&&(abs((*_lEta)[chosenlept[i]])>=2.1||chosenPt[i]<30)) {
                        TakeEvent=false;
                        break;
                    }
                }
                if (!TakeEvent) continue;
                
                
                /*
                 *CAT F
                 Mt2: F: 1l+2tau: the light lepton and the leading tau (could be of same sign)
                 */
                bool CatF=false;
                bool OSpairflag=false;
                double InvMassPair=0;
                double DistPeak=9999;
                double Mll=0;
                double Mt2=0;
                double Mt2_JECUp=0;
                double Mt2_JECDown=0;
                int SRIndex=9999;
                
                double ParticleMass[3]={0.00051, 0.105, 1.776};//Mass of Ele/Mu/Tau
                TLorentzVector l1;
                TLorentzVector l2;
                TLorentzVector Emiss;
                TLorentzVector Emiss_JECUp;
                TLorentzVector Emiss_JECDown;
                
                //l+tau+tau////////////
                if((*_flavors)[chosenlept[0]]!=2){
                    CatF=true;
                    //Calculate Mt2
                    if (chosenPt[1]>=chosenPt[2]) {//Find Leading Tau Lepton
                        l1.SetPtEtaPhiE((*_lPt)[chosenlept[1]],0, (*_lPhi)[chosenlept[1]],(*_lPt)[chosenlept[1]]);
                        l2.SetPtEtaPhiE((*_lPt)[chosenlept[0]],0, (*_lPhi)[chosenlept[0]],(*_lPt)[chosenlept[0]]);
                        Emiss.SetPtEtaPhiE(_met,0, _met_phi,_met);
                        Emiss_JECUp.SetPtEtaPhiE(_met_JECup,0, _met_phi_JECup,_met_JECup);
                        Emiss_JECDown.SetPtEtaPhiE(_met_JECdown,0, _met_phi_JECdown,_met_JECdown);
                        double pa[3]={ParticleMass[(*_flavors)[chosenlept[1]]],l1.Px(),l1.Py()};
                        double pb[3]={ParticleMass[(*_flavors)[chosenlept[0]]],l2.Px(),l2.Py()};
                        double pmiss[3]={0,Emiss.Px(),Emiss.Py()};
                        double pmiss_JECUp[3]={0,Emiss_JECUp.Px(),Emiss_JECUp.Py()};
                        double pmiss_JECDown[3]={0,Emiss_JECDown.Px(),Emiss_JECDown.Py()};
                        Mt2=CalcMt2(pa,pb,pmiss);
                        Mt2_JECUp=CalcMt2(pa,pb,pmiss_JECUp);
                        Mt2_JECDown=CalcMt2(pa,pb,pmiss_JECDown);
                    }
                    else {//Find Leading Tau Lepton
                        l1.SetPtEtaPhiE((*_lPt)[chosenlept[2]],0, (*_lPhi)[chosenlept[2]],(*_lPt)[chosenlept[2]]);
                        l2.SetPtEtaPhiE((*_lPt)[chosenlept[0]],0, (*_lPhi)[chosenlept[0]],(*_lPt)[chosenlept[0]]);
                        Emiss.SetPtEtaPhiE(_met,0, _met_phi,_met);
                        Emiss_JECUp.SetPtEtaPhiE(_met_JECup,0, _met_phi_JECup,_met_JECup);
                        Emiss_JECDown.SetPtEtaPhiE(_met_JECdown,0, _met_phi_JECdown,_met_JECdown);
                        double pa[3]={ParticleMass[(*_flavors)[chosenlept[2]]],l1.Px(),l1.Py()};
                        double pb[3]={ParticleMass[(*_flavors)[chosenlept[0]]],l2.Px(),l2.Py()};
                        double pmiss[3]={0,Emiss.Px(),Emiss.Py()};
                        double pmiss_JECUp[3]={0,Emiss_JECUp.Px(),Emiss_JECUp.Py()};
                        double pmiss_JECDown[3]={0,Emiss_JECDown.Px(),Emiss_JECDown.Py()};
                        Mt2=CalcMt2(pa,pb,pmiss);
                        Mt2_JECUp=CalcMt2(pa,pb,pmiss_JECUp);
                        Mt2_JECDown=CalcMt2(pa,pb,pmiss_JECDown);
                    }
                    //Calculate Mll
                    if ((*_charges)[chosenlept[0]]==(*_charges)[chosenlept[1]]&&(*_charges)[chosenlept[0]]==(*_charges)[chosenlept[2]]) {
                        Mll=0;
                    }
                    else {
                        OSpairflag=true;
                        if ((*_charges)[chosenlept[1]]==(*_charges)[chosenlept[2]]) {//tau1=tau2
                            for (int i=0; i<3; i++) {
                                if(i==0) continue;
                                InvMassPair=CalcInvariantMass(chosenlept[i],chosenlept[0]);
                                if (abs(InvMassPair-60)<DistPeak) {
                                    //Choose the lepton with inv. mass closest to Z mass as the Mll variable
                                    DistPeak=abs(InvMassPair-60);
                                    Mll=InvMassPair;
                                }
                            }
                        }
                        else{//tau1!=tau2 //Calculate Mll using OSSF tau
                            Mll=CalcInvariantMass(chosenlept[1],chosenlept[2]);
                        }
                    }
                }
                //l+tau+tau////////////
                
                //tau+l+tau////////////
                if((*_flavors)[chosenlept[1]]!=2){
                    CatF=true;
                    //Calculate Mt2
                    if (chosenPt[0]>=chosenPt[2]) {//Find Leading Tau Lepton
                        l1.SetPtEtaPhiE((*_lPt)[chosenlept[0]],0, (*_lPhi)[chosenlept[0]],(*_lPt)[chosenlept[0]]);
                        l2.SetPtEtaPhiE((*_lPt)[chosenlept[1]],0, (*_lPhi)[chosenlept[1]],(*_lPt)[chosenlept[1]]);
                        Emiss.SetPtEtaPhiE(_met,0, _met_phi,_met);
                        Emiss_JECUp.SetPtEtaPhiE(_met_JECup,0, _met_phi_JECup,_met_JECup);
                        Emiss_JECDown.SetPtEtaPhiE(_met_JECdown,0, _met_phi_JECdown,_met_JECdown);
                        double pa[3]={ParticleMass[(*_flavors)[chosenlept[0]]],l1.Px(),l1.Py()};
                        double pb[3]={ParticleMass[(*_flavors)[chosenlept[1]]],l2.Px(),l2.Py()};
                        double pmiss[3]={0,Emiss.Px(),Emiss.Py()};
                        double pmiss_JECUp[3]={0,Emiss_JECUp.Px(),Emiss_JECUp.Py()};
                        double pmiss_JECDown[3]={0,Emiss_JECDown.Px(),Emiss_JECDown.Py()};
                        Mt2=CalcMt2(pa,pb,pmiss);
                        Mt2_JECUp=CalcMt2(pa,pb,pmiss_JECUp);
                        Mt2_JECDown=CalcMt2(pa,pb,pmiss_JECDown);
                    }
                    else {//Find Leading Tau Lepton
                        l1.SetPtEtaPhiE((*_lPt)[chosenlept[2]],0, (*_lPhi)[chosenlept[2]],(*_lPt)[chosenlept[2]]);
                        l2.SetPtEtaPhiE((*_lPt)[chosenlept[1]],0, (*_lPhi)[chosenlept[1]],(*_lPt)[chosenlept[1]]);
                        Emiss.SetPtEtaPhiE(_met,0, _met_phi,_met);
                        Emiss_JECUp.SetPtEtaPhiE(_met_JECup,0, _met_phi_JECup,_met_JECup);
                        Emiss_JECDown.SetPtEtaPhiE(_met_JECdown,0, _met_phi_JECdown,_met_JECdown);
                        double pa[3]={ParticleMass[(*_flavors)[chosenlept[2]]],l1.Px(),l1.Py()};
                        double pb[3]={ParticleMass[(*_flavors)[chosenlept[1]]],l2.Px(),l2.Py()};
                        double pmiss[3]={0,Emiss.Px(),Emiss.Py()};
                        double pmiss_JECUp[3]={0,Emiss_JECUp.Px(),Emiss_JECUp.Py()};
                        double pmiss_JECDown[3]={0,Emiss_JECDown.Px(),Emiss_JECDown.Py()};
                        Mt2=CalcMt2(pa,pb,pmiss);
                        Mt2_JECUp=CalcMt2(pa,pb,pmiss_JECUp);
                        Mt2_JECDown=CalcMt2(pa,pb,pmiss_JECDown);
                    }
                    //Calculate Mll
                    if ((*_charges)[chosenlept[0]]==(*_charges)[chosenlept[1]]&&(*_charges)[chosenlept[0]]==(*_charges)[chosenlept[2]]) {
                        Mll=0;
                    }
                    else {
                        OSpairflag=true;
                        if ((*_charges)[chosenlept[0]]==(*_charges)[chosenlept[2]]) {//tau1=tau2
                            for (int i=0; i<3; i++) {
                                if(i==1) continue;
                                InvMassPair=CalcInvariantMass(chosenlept[i],chosenlept[1]);
                                if (abs(InvMassPair-60)<DistPeak) {
                                    //Choose the lepton with inv. mass closest to Z mass as the Mll variable
                                    DistPeak=abs(InvMassPair-60);
                                    Mll=InvMassPair;
                                }
                            }
                        }
                        else{//tau1!=tau2 //Calculate Mll using OSSF tau
                            Mll=CalcInvariantMass(chosenlept[0],chosenlept[2]);
                        }
                    }
                }
                //tau+l+tau////////////
                
                //tau+tau+l////////////
                if((*_flavors)[chosenlept[2]]!=2){
                    CatF=true;
                    //Calculate Mt2
                    if (chosenPt[0]>=chosenPt[1]) {//Find Leading Tau Lepton
                        l1.SetPtEtaPhiE((*_lPt)[chosenlept[0]],0, (*_lPhi)[chosenlept[0]],(*_lPt)[chosenlept[0]]);
                        l2.SetPtEtaPhiE((*_lPt)[chosenlept[2]],0, (*_lPhi)[chosenlept[2]],(*_lPt)[chosenlept[2]]);
                        Emiss.SetPtEtaPhiE(_met,0, _met_phi,_met);
                        Emiss_JECUp.SetPtEtaPhiE(_met_JECup,0, _met_phi_JECup,_met_JECup);
                        Emiss_JECDown.SetPtEtaPhiE(_met_JECdown,0, _met_phi_JECdown,_met_JECdown);
                        double pa[3]={ParticleMass[(*_flavors)[chosenlept[0]]],l1.Px(),l1.Py()};
                        double pb[3]={ParticleMass[(*_flavors)[chosenlept[2]]],l2.Px(),l2.Py()};
                        double pmiss[3]={0,Emiss.Px(),Emiss.Py()};
                        double pmiss_JECUp[3]={0,Emiss_JECUp.Px(),Emiss_JECUp.Py()};
                        double pmiss_JECDown[3]={0,Emiss_JECDown.Px(),Emiss_JECDown.Py()};
                        Mt2=CalcMt2(pa,pb,pmiss);
                        Mt2_JECUp=CalcMt2(pa,pb,pmiss_JECUp);
                        Mt2_JECDown=CalcMt2(pa,pb,pmiss_JECDown);
                    }
                    else {//Find Leading Tau Lepton
                        l1.SetPtEtaPhiE((*_lPt)[chosenlept[1]],0, (*_lPhi)[chosenlept[1]],(*_lPt)[chosenlept[1]]);
                        l2.SetPtEtaPhiE((*_lPt)[chosenlept[2]],0, (*_lPhi)[chosenlept[2]],(*_lPt)[chosenlept[2]]);
                        Emiss.SetPtEtaPhiE(_met,0, _met_phi,_met);
                        Emiss_JECUp.SetPtEtaPhiE(_met_JECup,0, _met_phi_JECup,_met_JECup);
                        Emiss_JECDown.SetPtEtaPhiE(_met_JECdown,0, _met_phi_JECdown,_met_JECdown);
                        double pa[3]={ParticleMass[(*_flavors)[chosenlept[1]]],l1.Px(),l1.Py()};
                        double pb[3]={ParticleMass[(*_flavors)[chosenlept[2]]],l2.Px(),l2.Py()};
                        double pmiss[3]={0,Emiss.Px(),Emiss.Py()};
                        double pmiss_JECUp[3]={0,Emiss_JECUp.Px(),Emiss_JECUp.Py()};
                        double pmiss_JECDown[3]={0,Emiss_JECDown.Px(),Emiss_JECDown.Py()};
                        Mt2=CalcMt2(pa,pb,pmiss);
                        Mt2_JECUp=CalcMt2(pa,pb,pmiss_JECUp);
                        Mt2_JECDown=CalcMt2(pa,pb,pmiss_JECDown);
                    }
                    //Calculate Mll
                    if ((*_charges)[chosenlept[0]]==(*_charges)[chosenlept[1]]&&(*_charges)[chosenlept[0]]==(*_charges)[chosenlept[2]]) {
                        Mll=0;
                    }
                    else {
                        OSpairflag=true;
                        if ((*_charges)[chosenlept[0]]==(*_charges)[chosenlept[1]]) {//tau1=tau2
                            for (int i=0; i<3; i++) {
                                if(i==2) continue;
                                InvMassPair=CalcInvariantMass(chosenlept[i],chosenlept[2]);
                                if (abs(InvMassPair-60)<DistPeak) {
                                    //Choose the lepton with inv. mass closest to Z mass as the Mll variable
                                    DistPeak=abs(InvMassPair-60);
                                    Mll=InvMassPair;
                                }
                            }
                        }
                        else{//tau1!=tau2 //Calculate Mll using OSSF tau
                            Mll=CalcInvariantMass(chosenlept[0],chosenlept[1]);
                        }
                    }
                }
                //tau+tau+l////////////
                
                /*
                 *
                 */
                ////////////////////////////////////////////////////////////////
                ///////Fill The Hist////////////////////////////////////////////
                ////////////////////////////////////////////////////////////////
                /*
                 *
                 */
                //////Exclude Low Mll//////////////////
                if (!CatF) continue;
                if (Mll<12&&OSpairflag) continue;
                
                double M3l=CalcTrilepInvariantMass(chosenlept[0],chosenlept[1],chosenlept[2]);
                
                
                //JECUp
                if (!bJetandMETflag_JECUp) {
                    if (CatF) {
                        SRIndex=CalcSRIndex_F(_met_JECup,Mll,Mt2_JECUp);
                        MCYield_Conversion_JECUp[5]->Fill(SRIndex,W);
                    }
                }
                //JECDown
                if (!bJetandMETflag_JECDown) {
                    if (CatF) {
                        SRIndex=CalcSRIndex_F(_met_JECdown,Mll,Mt2_JECDown);
                        MCYield_Conversion_JECDown[5]->Fill(SRIndex,W);
                    }
                }
                //NO JEC
                if (!bJetandMETflag) {
                    ///////////////////////////////////////////////
                    //Fill SR Hist Before Overflow Bin Adjustment//
                    ///////////////////////////////////////////////
                    if (CatF) {
                        SRIndex=CalcSRIndex_F(_met,Mll,Mt2);
                        MCYield_Conversion[5]->Fill(SRIndex,W);
                        MCYield_Conversion_BTagUp[5]->Fill(SRIndex,W_BTagUp);
                        MCYield_Conversion_BTagDown[5]->Fill(SRIndex,W_BTagDown);
                        MCYield_Conversion_LepSFUp[5]->Fill(SRIndex,W_LepSFUp);
                        MCYield_Conversion_LepSFDown[5]->Fill(SRIndex,W_LepSFDown);
                        MCYield_Conversion_PUUp[5]->Fill(SRIndex,W_PUUp);
                        MCYield_Conversion_PUDown[5]->Fill(SRIndex,W_PUDown);
                    }
                    ///////////////////////////////////////////////
                    //Fill SR Hist Before Overflow Bin Adjustment//
                    ///////////////////////////////////////////////
                    ////////////////////
                    ///Overflow Bin/////
                    ////////////////////
                    int CatNum = 9999;
                    if (CatF) {
                        CatNum=5;
                    }
                    
                    double leadingPt_FL=leadingPt;
                    double subleadingPt_FL=subleadingPt;
                    double trailingPt_FL=trailingPt;
                    double _met_FL=_met;
                    double Mt2_FL=Mt2;
                    double Mll_FL=Mll;
                    double M3l_FL=M3l;
                    if (leadingPt_FL>Ptupperbound[CatNum]) leadingPt_FL=Ptupperbound[CatNum]-1;
                    if (subleadingPt_FL>Ptupperbound[CatNum]) subleadingPt_FL=Ptupperbound[CatNum]-1;
                    if (trailingPt_FL>Ptupperbound[CatNum]) trailingPt_FL=Ptupperbound[CatNum]-1;
                    if (_met_FL>METupperbound[CatNum]) _met_FL=METupperbound[CatNum]-1;
                    if (Mt2_FL>MTupperbound[CatNum]) Mt2_FL=MTupperbound[CatNum]-1;
                    if (Mll_FL>Mllupperbound[CatNum]) Mll_FL=Mllupperbound[CatNum]-1;
                    if (M3l_FL>M3lupperbound[CatNum]) M3l_FL=M3lupperbound[CatNum]-1;
                    
                    if (CatF) {
                        MCYield_Conversion_Pt1[5]->Fill(leadingPt_FL,W);
                        MCYield_Conversion_Pt2[5]->Fill(subleadingPt_FL,W);
                        MCYield_Conversion_Pt3[5]->Fill(trailingPt_FL,W);
                        MCYield_Conversion_Mll[5]->Fill(Mll_FL,W);
                        MCYield_Conversion_MT[5]->Fill(Mt2_FL,W);
                        MCYield_Conversion_MET[5]->Fill(_met_FL,W);
                        MCYield_Conversion_M3l[5]->Fill(M3l_FL,W);
                        MCYield_Conversion_SumQ[5]->Fill(SumQ,W);
                    }
                }//NO JEC
                
            }//if(nTau==2)
            /*
             *2Tau
             */
            ////////////////////////////////////////////////////////////////////////////////////

            
        }//loop over events/entries
        
        
        cout<<endl;
        
    }//loop over samples

    ////////////////////////////////////////////////////////////////////////////////////
    ////Assign uncertainties and prepare hists for limit setting and plots//////////////
    ////////////////////////////////////////////////////////////////////////////////////
    
    //Hists for limit setting
    //None for conversion
    //Hists for limit setting
    
    //Hists for plots
    TH1F* Plot_MCYield_Conversion[nSRCat];
    TH1F* Plot_MCYield_Conversion_Pt1[nSRCat];
    TH1F* Plot_MCYield_Conversion_Pt2[nSRCat];
    TH1F* Plot_MCYield_Conversion_Pt3[nSRCat];
    TH1F* Plot_MCYield_Conversion_Mll[nSRCat];
    TH1F* Plot_MCYield_Conversion_MT[nSRCat];
    TH1F* Plot_MCYield_Conversion_MET[nSRCat];
    TH1F* Plot_MCYield_Conversion_M3l[nSRCat];
    TH1F* Plot_MCYield_Conversion_SumQ[nSRCat];
    
    TH1F* Plot_StatUnc_MCYield_Conversion[nSRCat];
    TH1F* Plot_StatUnc_MCYield_Conversion_Pt1[nSRCat];
    TH1F* Plot_StatUnc_MCYield_Conversion_Pt2[nSRCat];
    TH1F* Plot_StatUnc_MCYield_Conversion_Pt3[nSRCat];
    TH1F* Plot_StatUnc_MCYield_Conversion_Mll[nSRCat];
    TH1F* Plot_StatUnc_MCYield_Conversion_MT[nSRCat];
    TH1F* Plot_StatUnc_MCYield_Conversion_MET[nSRCat];
    TH1F* Plot_StatUnc_MCYield_Conversion_M3l[nSRCat];
    TH1F* Plot_StatUnc_MCYield_Conversion_SumQ[nSRCat];
    
    TH1F* Plot_TotUnc_MCYield_Conversion[nSRCat];
    TH1F* Plot_TotUnc_MCYield_Conversion_Pt1[nSRCat];
    TH1F* Plot_TotUnc_MCYield_Conversion_Pt2[nSRCat];
    TH1F* Plot_TotUnc_MCYield_Conversion_Pt3[nSRCat];
    TH1F* Plot_TotUnc_MCYield_Conversion_Mll[nSRCat];
    TH1F* Plot_TotUnc_MCYield_Conversion_MT[nSRCat];
    TH1F* Plot_TotUnc_MCYield_Conversion_MET[nSRCat];
    TH1F* Plot_TotUnc_MCYield_Conversion_M3l[nSRCat];
    TH1F* Plot_TotUnc_MCYield_Conversion_SumQ[nSRCat];
    //Hists for plots
    
    for (int i=0; i<nSRCat; i++) {
        
        Plot_MCYield_Conversion[i]=(TH1F*)MCYield_Conversion[i]->Clone();
        Plot_MCYield_Conversion_Pt1[i]=(TH1F*)MCYield_Conversion_Pt1[i]->Clone();
        Plot_MCYield_Conversion_Pt2[i]=(TH1F*)MCYield_Conversion_Pt2[i]->Clone();
        Plot_MCYield_Conversion_Pt3[i]=(TH1F*)MCYield_Conversion_Pt3[i]->Clone();
        Plot_MCYield_Conversion_Mll[i]=(TH1F*)MCYield_Conversion_Mll[i]->Clone();
        Plot_MCYield_Conversion_MT[i]=(TH1F*)MCYield_Conversion_MT[i]->Clone();
        Plot_MCYield_Conversion_MET[i]=(TH1F*)MCYield_Conversion_MET[i]->Clone();
        Plot_MCYield_Conversion_M3l[i]=(TH1F*)MCYield_Conversion_M3l[i]->Clone();
        Plot_MCYield_Conversion_SumQ[i]=(TH1F*)MCYield_Conversion_SumQ[i]->Clone();
        
        name="Plot_StatUnc_MCBG_Conversion_"+SRCat[i];
        Plot_StatUnc_MCYield_Conversion[i]=new TH1F(name,"",nSR[i],0.5,nSR[i]+0.5);
        Plot_StatUnc_MCYield_Conversion[i]->Sumw2();
        name="Plot_StatUnc_MCBG_Conversion_Pt1_"+SRCat[i];
        Plot_StatUnc_MCYield_Conversion_Pt1[i]=new TH1F(name,"",Ptbin[i],Ptlowerbound,Ptupperbound[i]);
        Plot_StatUnc_MCYield_Conversion_Pt1[i]->Sumw2();
        name="Plot_StatUnc_MCBG_Conversion_Pt2_"+SRCat[i];
        Plot_StatUnc_MCYield_Conversion_Pt2[i]=new TH1F(name,"",Ptbin[i],Ptlowerbound,Ptupperbound[i]);
        Plot_StatUnc_MCYield_Conversion_Pt2[i]->Sumw2();
        name="Plot_StatUnc_MCBG_Conversion_Pt3_"+SRCat[i];
        Plot_StatUnc_MCYield_Conversion_Pt3[i]=new TH1F(name,"",Ptbin[i],Ptlowerbound,Ptupperbound[i]);
        Plot_StatUnc_MCYield_Conversion_Pt3[i]->Sumw2();
        name="Plot_StatUnc_MCBG_Conversion_Mll_"+SRCat[i];
        Plot_StatUnc_MCYield_Conversion_Mll[i]=new TH1F(name,"",Mllbin[i],Mlllowerbound,Mllupperbound[i]);
        Plot_StatUnc_MCYield_Conversion_Mll[i]->Sumw2();
        name="Plot_StatUnc_MCBG_Conversion_MT_"+SRCat[i];
        Plot_StatUnc_MCYield_Conversion_MT[i]=new TH1F(name,"",MTbin[i],MTlowerbound,MTupperbound[i]);
        Plot_StatUnc_MCYield_Conversion_MT[i]->Sumw2();
        name="Plot_StatUnc_MCBG_Conversion_MET_"+SRCat[i];
        Plot_StatUnc_MCYield_Conversion_MET[i]=new TH1F(name,"",METbin[i],METlowerbound,METupperbound[i]);
        Plot_StatUnc_MCYield_Conversion_MET[i]->Sumw2();
        name="Plot_StatUnc_MCBG_Conversion_M3l_"+SRCat[i];
        Plot_StatUnc_MCYield_Conversion_M3l[i]=new TH1F(name,"",M3lbin[i],M3llowerbound,M3lupperbound[i]);
        Plot_StatUnc_MCYield_Conversion_M3l[i]->Sumw2();
        name="Plot_StatUnc_MCBG_Conversion_SumQ_"+SRCat[i];
        Plot_StatUnc_MCYield_Conversion_SumQ[i]=new TH1F(name,"",SumQbin[i],SumQlowerbound,SumQupperbound[i]);
        Plot_StatUnc_MCYield_Conversion_SumQ[i]->Sumw2();
        
        name="Plot_TotUnc_MCBG_Conversion_"+SRCat[i];
        Plot_TotUnc_MCYield_Conversion[i]=new TH1F(name,"",nSR[i],0.5,nSR[i]+0.5);
        Plot_TotUnc_MCYield_Conversion[i]->Sumw2();
        name="Plot_TotUnc_MCBG_Conversion_Pt1_"+SRCat[i];
        Plot_TotUnc_MCYield_Conversion_Pt1[i]=new TH1F(name,"",Ptbin[i],Ptlowerbound,Ptupperbound[i]);
        Plot_TotUnc_MCYield_Conversion_Pt1[i]->Sumw2();
        name="Plot_TotUnc_MCBG_Conversion_Pt2_"+SRCat[i];
        Plot_TotUnc_MCYield_Conversion_Pt2[i]=new TH1F(name,"",Ptbin[i],Ptlowerbound,Ptupperbound[i]);
        Plot_TotUnc_MCYield_Conversion_Pt2[i]->Sumw2();
        name="Plot_TotUnc_MCBG_Conversion_Pt3_"+SRCat[i];
        Plot_TotUnc_MCYield_Conversion_Pt3[i]=new TH1F(name,"",Ptbin[i],Ptlowerbound,Ptupperbound[i]);
        Plot_TotUnc_MCYield_Conversion_Pt3[i]->Sumw2();
        name="Plot_TotUnc_MCBG_Conversion_Mll_"+SRCat[i];
        Plot_TotUnc_MCYield_Conversion_Mll[i]=new TH1F(name,"",Mllbin[i],Mlllowerbound,Mllupperbound[i]);
        Plot_TotUnc_MCYield_Conversion_Mll[i]->Sumw2();
        name="Plot_TotUnc_MCBG_Conversion_MT_"+SRCat[i];
        Plot_TotUnc_MCYield_Conversion_MT[i]=new TH1F(name,"",MTbin[i],MTlowerbound,MTupperbound[i]);
        Plot_TotUnc_MCYield_Conversion_MT[i]->Sumw2();
        name="Plot_TotUnc_MCBG_Conversion_MET_"+SRCat[i];
        Plot_TotUnc_MCYield_Conversion_MET[i]=new TH1F(name,"",METbin[i],METlowerbound,METupperbound[i]);
        Plot_TotUnc_MCYield_Conversion_MET[i]->Sumw2();
        name="Plot_TotUnc_MCBG_Conversion_M3l_"+SRCat[i];
        Plot_TotUnc_MCYield_Conversion_M3l[i]=new TH1F(name,"",M3lbin[i],M3llowerbound,M3lupperbound[i]);
        Plot_TotUnc_MCYield_Conversion_M3l[i]->Sumw2();
        name="Plot_TotUnc_MCBG_Conversion_SumQ_"+SRCat[i];
        Plot_TotUnc_MCYield_Conversion_SumQ[i]=new TH1F(name,"",SumQbin[i],SumQlowerbound,SumQupperbound[i]);
        Plot_TotUnc_MCYield_Conversion_SumQ[i]->Sumw2();
    }
    
    /////////////////////////////////////////////////////////////////////////
    //Assign uncertainties///////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////
    for (int i=0; i<nSRCat; i++) {
        
        
        /*
         ttX: 15%
         ZZ/H: 25%
         VVV: 50%
         Non-prompt: 30%
         WZ: Shape+9%(SR)/5%+9%(Var)
         Conversion: 15%
         */
        int nbinSR=nSR[i];
        int nbinPt=Ptbin[i];
        int nbinMll=Mllbin[i];
        int nbinMT=MTbin[i];
        int nbinMET=METbin[i];
        int nbinM3l=M3lbin[i];
        int nbinSumQ=SumQbin[i];
        double uncertainty=0.15;
        //Hists for limit setting
        /*
         *None for conversion
         */
        
        //Hists for plots
        /*
         *All uncertainties
         */
        //SR
        AssignUncForPlot_SRBG_NonWZ(Plot_MCYield_Conversion[i], uncertainty, nbinSR, MCYield_Conversion_JECUp[i], MCYield_Conversion_JECDown[i], MCYield_Conversion_BTagUp[i],MCYield_Conversion_BTagDown[i],MCYield_Conversion_LepSFUp[i],MCYield_Conversion_LepSFDown[i],MCYield_Conversion_PUUp[i],MCYield_Conversion_PUDown[i],Plot_StatUnc_MCYield_Conversion[i],Plot_TotUnc_MCYield_Conversion[i]);
        //Var
        AssignUncForPlot_VarBG_NonWZ(Plot_MCYield_Conversion_Pt1[i], uncertainty, nbinPt, Plot_StatUnc_MCYield_Conversion_Pt1[i],Plot_TotUnc_MCYield_Conversion_Pt1[i]);
        AssignUncForPlot_VarBG_NonWZ(Plot_MCYield_Conversion_Pt2[i], uncertainty, nbinPt, Plot_StatUnc_MCYield_Conversion_Pt2[i],Plot_TotUnc_MCYield_Conversion_Pt2[i]);
        AssignUncForPlot_VarBG_NonWZ(Plot_MCYield_Conversion_Pt3[i], uncertainty, nbinPt, Plot_StatUnc_MCYield_Conversion_Pt3[i],Plot_TotUnc_MCYield_Conversion_Pt3[i]);
        AssignUncForPlot_VarBG_NonWZ(Plot_MCYield_Conversion_Mll[i], uncertainty, nbinMll, Plot_StatUnc_MCYield_Conversion_Mll[i],Plot_TotUnc_MCYield_Conversion_Mll[i]);
        AssignUncForPlot_VarBG_NonWZ(Plot_MCYield_Conversion_MT[i], uncertainty, nbinMT, Plot_StatUnc_MCYield_Conversion_MT[i],Plot_TotUnc_MCYield_Conversion_MT[i]);
        AssignUncForPlot_VarBG_NonWZ(Plot_MCYield_Conversion_MET[i], uncertainty, nbinMET, Plot_StatUnc_MCYield_Conversion_MET[i],Plot_TotUnc_MCYield_Conversion_MET[i]);
        AssignUncForPlot_VarBG_NonWZ(Plot_MCYield_Conversion_M3l[i], uncertainty, nbinM3l, Plot_StatUnc_MCYield_Conversion_M3l[i],Plot_TotUnc_MCYield_Conversion_M3l[i]);
        AssignUncForPlot_VarBG_NonWZ(Plot_MCYield_Conversion_SumQ[i], uncertainty, nbinSumQ, Plot_StatUnc_MCYield_Conversion_SumQ[i],Plot_TotUnc_MCYield_Conversion_SumQ[i]);
        
    }
    
    
    
    /////////////////////////////////////////////////////////////////////////
    //Write to output files//////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////
    TFile *outputRootFile;
    //Original Hists
    outputRootFile=new TFile("FinalYield_Conversion.root","CREATE");
    for (int i=0; i<nSRCat; i++) {
        MCYield_Conversion[i]->Write();
        MCYield_Conversion_JECUp[i]->Write();
        MCYield_Conversion_JECDown[i]->Write();
        MCYield_Conversion_BTagUp[i]->Write();
        MCYield_Conversion_BTagDown[i]->Write();
        MCYield_Conversion_LepSFUp[i]->Write();
        MCYield_Conversion_LepSFDown[i]->Write();
        MCYield_Conversion_PUUp[i]->Write();
        MCYield_Conversion_PUDown[i]->Write();
        
        MCYield_Conversion_Pt1[i]->Write();
        MCYield_Conversion_Pt2[i]->Write();
        MCYield_Conversion_Pt3[i]->Write();
        MCYield_Conversion_Mll[i]->Write();
        MCYield_Conversion_MT[i]->Write();
        MCYield_Conversion_MET[i]->Write();
        MCYield_Conversion_M3l[i]->Write();
        MCYield_Conversion_SumQ[i]->Write();
        
    }
    outputRootFile->Close();
    //For limit setting
    outputRootFile=new TFile("FinalYield_Conversion_Limit.root","CREATE");
    for (int i=0; i<nSRCat; i++) {
        MCYield_Conversion[i]->Write();
        MCYield_Conversion_JECUp[i]->Write();
        MCYield_Conversion_JECDown[i]->Write();
        MCYield_Conversion_BTagUp[i]->Write();
        MCYield_Conversion_BTagDown[i]->Write();
        MCYield_Conversion_LepSFUp[i]->Write();
        MCYield_Conversion_LepSFDown[i]->Write();
        MCYield_Conversion_PUUp[i]->Write();
        MCYield_Conversion_PUDown[i]->Write();
    }
    outputRootFile->Close();
    //For plots
    outputRootFile=new TFile("FinalYield_Conversion_Plot.root","CREATE");
    for (int i=0; i<nSRCat; i++) {
        Plot_MCYield_Conversion[i]->Write();
        Plot_MCYield_Conversion_Pt1[i]->Write();
        Plot_MCYield_Conversion_Pt2[i]->Write();
        Plot_MCYield_Conversion_Pt3[i]->Write();
        Plot_MCYield_Conversion_Mll[i]->Write();
        Plot_MCYield_Conversion_MT[i]->Write();
        Plot_MCYield_Conversion_MET[i]->Write();
        Plot_MCYield_Conversion_M3l[i]->Write();
        Plot_MCYield_Conversion_SumQ[i]->Write();
        
        Plot_StatUnc_MCYield_Conversion[i]->Write();
        Plot_StatUnc_MCYield_Conversion_Pt1[i]->Write();
        Plot_StatUnc_MCYield_Conversion_Pt2[i]->Write();
        Plot_StatUnc_MCYield_Conversion_Pt3[i]->Write();
        Plot_StatUnc_MCYield_Conversion_Mll[i]->Write();
        Plot_StatUnc_MCYield_Conversion_MT[i]->Write();
        Plot_StatUnc_MCYield_Conversion_MET[i]->Write();
        Plot_StatUnc_MCYield_Conversion_M3l[i]->Write();
        Plot_StatUnc_MCYield_Conversion_SumQ[i]->Write();
        
        Plot_TotUnc_MCYield_Conversion[i]->Write();
        Plot_TotUnc_MCYield_Conversion_Pt1[i]->Write();
        Plot_TotUnc_MCYield_Conversion_Pt2[i]->Write();
        Plot_TotUnc_MCYield_Conversion_Pt3[i]->Write();
        Plot_TotUnc_MCYield_Conversion_Mll[i]->Write();
        Plot_TotUnc_MCYield_Conversion_MT[i]->Write();
        Plot_TotUnc_MCYield_Conversion_MET[i]->Write();
        Plot_TotUnc_MCYield_Conversion_M3l[i]->Write();
        Plot_TotUnc_MCYield_Conversion_SumQ[i]->Write();
    }
    outputRootFile->Close();
    

}//void FinalYield::loop


///////////////////////////////////////////////////////////////////////////////////////
//////////////////////////Appendix- Defined Functions//////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
//Trigger Selection
bool FinalYield::TriggerSelectionMC(TString& MCsample){
    /*
     TRIGGERS USED:
     
     @SingleElectron:
     passHLT_Ele27_WPTight_Gsf
     
     
     @SingleMuon:
     passHLT_IsoMu24
     passHLT_IsoTkMu24
     
     @DoubleEG:
     passHLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ
     
     @DoubleMuon:
     passHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL
     passHLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL
     passHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ (Only for RunH)
     passHLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ (Only for RunH: RunNo>280919)
     
     @MuonEG:
     passHLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL
     passHLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ (Only for RunH)
     passHLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL
     passHLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ (Only for RunH)
     */
    
    bool passtriggerselection=false;
    
    //NO REPEATED TRIGGER DROP
    if (passHLT_Ele27_WPTight_Gsf||
        passHLT_IsoMu24||passHLT_IsoTkMu24||
        passHLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ||
        passHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL||passHLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL||
        passHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ||passHLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ||
        passHLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL||passHLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL||
        passHLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ||passHLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ) {
        passtriggerselection=true;
    }

    return passtriggerselection;
}

//b-Jet and MET cut
bool FinalYield::bJetandMETcut(double& MET, const double& VarCoef){
    bool flag=false;
    if (MET<50) {
        flag=true;
    }
    for (int jetnum = 0; jetnum <(*_jetPt).size(); jetnum++){
        if ((*_bTagged)[jetnum]==1&&
            ((*_jetPt)[jetnum]+VarCoef*(*_jetJECuncty)[jetnum])>25) {
            flag=true;
            continue;
        }
    }
    return flag;
}

//Object Selection
bool FinalYield::ObjectSelection(int& n){
    bool passobjselection=false;
    if((*_flavors)[n]!=2){
        if ((*_istightIDWP2016_EWK)[n]&&(*_lPt)[n]>10) passobjselection=true;
    }
    else if((*_flavors)[n]==2){
        //if (abs((*_tau_dz)[n])<0.2&&abs((*_ipPV)[n])<1000&&(*_lPt)[n]>20&&abs((*_lEta)[n])<2.3&&(*_TauIs_decayModeFinding)[n]&&(*_TauIs_byLooseIsolationMVArun2v1DBdR03oldDMwLT)[n])
        //if ((*_lPt)[n]>20&&abs((*_lEta)[n])<2.3&&(*_TauIs_byTightIsolationMVArun2v1DBoldDMwLT)[n])passobjselection=true;
        if ((*_lPt)[n]>20&&abs((*_lEta)[n])<2.3&&(*_istightIDWP2016_EWK)[n])passobjselection=true;//Tight tau selected by _istightID variable
    }
    return passobjselection;
}

//FO Selection
bool FinalYield::FOSelection(int& n){
    bool passFOselection=false;
    if((*_flavors)[n]!=2){
        if ((*_isFOIDWP2016_EWK)[n]&&(*_lPt)[n]>10) passFOselection=true;
    }
    else if((*_flavors)[n]==2){
        if ((*_lPt)[n]>20&&abs((*_lEta)[n])<2.3&&(*_isFOIDWP2016_EWK)[n]) passFOselection=true;
    }
    return passFOselection;
}

//DeltaRSelection
bool FinalYield::DeltaRSelection(int& n,int& l){
    bool passdRselection=true;
    //electrons should not be overlapped (within dR=0.05) with a LOOSE muon (Pt>5)
    if((*_flavors)[n]==0){
        for (int i=0; i<l; i++){
            if((*_isVetoIDWP2016_EWK)[n]&&(*_lPt)[n]>5&&(*_flavors)[i]==1){
                double pi=3.1415926;
                double dphi=(abs((*_lPhi)[n]-(*_lPhi)[i])>pi)?(2*pi-abs((*_lPhi)[n]-(*_lPhi)[i])):abs((*_lPhi)[n]-(*_lPhi)[i]);
                double deta=abs((*_lEta)[n]-(*_lEta)[i]);
                double dr=sqrt(dphi*dphi+deta*deta);
                if(dr<0.05) {
                    passdRselection=false;
                    break;
                }
            }
        }
    }
    //remove taus within DeltaR < 0.4 of LOOSE (AN-16-399 P27) muons and electrons (Pt>10)
    if((*_flavors)[n]==2){
        for (int i=0; i<l; i++){
            if((*_isVetoIDWP2016_EWK)[i]&&(*_flavors)[i]!=2){//Only check for mu/ele
                double pi=3.1415926;
                double dphi=(abs((*_lPhi)[n]-(*_lPhi)[i])>pi)?(2*pi-abs((*_lPhi)[n]-(*_lPhi)[i])):abs((*_lPhi)[n]-(*_lPhi)[i]);
                double deta=abs((*_lEta)[n]-(*_lEta)[i]);
                double dr=sqrt(dphi*dphi+deta*deta);
                if(dr<0.4) {
                    passdRselection=false;
                    break;
                }
            }
        }
    }
    return passdRselection;
}

//Define the function calculating the invariant mass of di-lepton pair
double FinalYield::CalcInvariantMass(int& n,int& l){
    TLorentzVector lep1;
    TLorentzVector lep2;
    
    Float_t et1 = (*_lPt)[n];
    Float_t et2 = (*_lPt)[l];
    lep1.SetPtEtaPhiE(et1, (*_lEta)[n], (*_lPhi)[n], (et1 * cosh((*_lEta)[n])));
    lep2.SetPtEtaPhiE(et2, (*_lEta)[l], (*_lPhi)[l], (et2 * cosh((*_lEta)[l])));
    
    return (lep1+lep2).Mag();
}

double FinalYield::CalcTrilepInvariantMass(int& n,int& l, int& m){
    TLorentzVector lep1;
    TLorentzVector lep2;
    TLorentzVector lep3;
    
    Float_t et1 = (*_lPt)[n];
    Float_t et2 = (*_lPt)[l];
    Float_t et3 = (*_lPt)[m];
    lep1.SetPtEtaPhiE(et1, (*_lEta)[n], (*_lPhi)[n], (et1 * cosh((*_lEta)[n])));
    lep2.SetPtEtaPhiE(et2, (*_lEta)[l], (*_lPhi)[l], (et2 * cosh((*_lEta)[l])));
    lep3.SetPtEtaPhiE(et3, (*_lEta)[m], (*_lPhi)[m], (et3 * cosh((*_lEta)[m])));
    
    return (lep1+lep2+lep3).Mag();
}

//Define the function calculating the transverse mass
double FinalYield::CalcTransMass(int& n, double& MET, double& METPhi){
    double TransMass=sqrt(2*(*_lPt)[n]*MET*(1-cos((*_lPhi)[n]-METPhi)));
    return TransMass;
}

//Define the function calculating MT2
double FinalYield::CalcMt2(double* pa, double* pb, double* pmiss){
    double Mt2;
    mt2_calculator.set_momenta(pa,pb,pmiss);
    mt2_calculator.set_mn(0);
    Mt2=mt2_calculator.get_mt2();
    return Mt2;
}

//Define the function calculating dR(jet,FO)
double FinalYield::CalcdR(int& n, int& l){
    double pi=3.1415926;
    double dphi=(abs((*_jetPhi)[n]-(*_lPhi)[l])>pi)?(2*pi-abs((*_jetPhi)[n]-(*_lPhi)[l])):abs((*_jetPhi)[n]-(*_lPhi)[l]);
    double deta=abs((*_jetEta)[n]-(*_lEta)[l]);
    double dr=sqrt(dphi*dphi+deta*deta);
    return dr;
}

//Define the function calculating recoil vector
/*
double FinalYield::Calcu1(double& uphi, double& bosonphi, double& umag){
    double pi=3.1415926;
    double dphi=(abs(uphi-bosonphi)>pi)?(2*pi-abs(uphi-bosonphi)):abs(uphi-bosonphi);
    double u1=umag*cos(dphi);
    return u1;
}
double FinalYield::Calcu2(double& uphi, double& bosonphi, double& umag){
    double pi=3.1415926;
    double dphi=(abs(uphi-bosonphi)>pi)?(2*pi-abs(uphi-bosonphi)):abs(uphi-bosonphi);
    double u2=umag*sin(dphi);
    if(uphi>=0){
        if(bosonphi>=0){
            if(bosonphi<uphi)u2=-u2;
        }
        if(bosonphi<0){
            if(bosonphi>(uphi-pi))u2=-u2;
        }
    }
    if(uphi<0){
        if(bosonphi>=0){
            if(bosonphi>(uphi+pi))u2=-u2;
        }
        if(bosonphi<0){
            if(bosonphi<uphi)u2=-u2;
        }
    }
    return u2;
}
*/

//Function Calculating B-Tagging SF
double FinalYield::CalcBTagSF(const double& VarCoef){
    double PMC=1;
    double PData=1;
    
    //Central
    if (VarCoef==0) {
        for (int jetnum = 0; jetnum <(*_jetPt).size(); jetnum++){
            //A B-Tag Jet
            if ((*_bTagged)[jetnum]==1&&
                (*_jetPt)[jetnum]>25) {
                PMC=PMC*(*_jetbtagEff)[jetnum];
                PData=PData*(*_jetbtagSF)[jetnum]*(*_jetbtagEff)[jetnum];
            }
            //Not A B-Tag Jet
            else {
                PMC=PMC*(1-(*_jetbtagEff)[jetnum]);
                PData=PData*(1-(*_jetbtagSF)[jetnum]*(*_jetbtagEff)[jetnum]);
            }
        }
    }
    //Up
    if (VarCoef==1) {
        for (int jetnum = 0; jetnum <(*_jetPt).size(); jetnum++){
            //A B-Tag Jet
            if ((*_bTagged)[jetnum]==1&&
                (*_jetPt)[jetnum]>25) {
                PMC=PMC*(*_jetbtagEff)[jetnum];
                PData=PData*(*_jetbtagSF_up)[jetnum]*(*_jetbtagEff)[jetnum];
            }
            //Not A B-Tag Jet
            else {
                PMC=PMC*(1-(*_jetbtagEff)[jetnum]);
                PData=PData*(1-(*_jetbtagSF_up)[jetnum]*(*_jetbtagEff)[jetnum]);
            }
        }
    }
    //Down
    if (VarCoef==-1) {
        for (int jetnum = 0; jetnum <(*_jetPt).size(); jetnum++){
            //A B-Tag Jet
            if ((*_bTagged)[jetnum]==1&&
                (*_jetPt)[jetnum]>25) {
                PMC=PMC*(*_jetbtagEff)[jetnum];
                PData=PData*(*_jetbtagSF_down)[jetnum]*(*_jetbtagEff)[jetnum];
            }
            //Not A B-Tag Jet
            else {
                PMC=PMC*(1-(*_jetbtagEff)[jetnum]);
                PData=PData*(1-(*_jetbtagSF_down)[jetnum]*(*_jetbtagEff)[jetnum]);
            }
        }
    }
    
    if (PMC==0||TMath::IsNaN(PMC)) PMC=1;
    if (PData==0||TMath::IsNaN(PData)) PData=1;
    
    return (PData/PMC);
}

//Function Calculating ele LepSF (Data-FullSim)
double FinalYield::CalcLepSF_ele_FullSim(TH2D* SF, TH2D* RecoSF, int& n, const double& VarCoef){
    
    int binx,biny;
    
    //LepSF
    //For current SF xbin=6(1-6), ybin=5(6-10)
    binx=SF->GetXaxis()->FindBin((*_lPt)[n]);
    if (binx<1) binx=1;
    if (binx>6) binx=6;//Control Overflow
    biny=SF->GetYaxis()->FindBin(fabs((*_lEta)[n]));
    if (biny<6) biny=6;
    if (biny>10) biny=10;//Control Overflow
    double eleSF=SF->GetBinContent(binx,biny)+VarCoef*(SF->GetBinError(binx,biny));
    
    //RecoSF
    //For current RecoSF xbin=30(1-30), ybin=1
    binx=RecoSF->GetXaxis()->FindBin((*_lEta)[n]);
    if (binx<1) binx=1;
    if (binx>30) binx=30;
    double eleRecoSF=RecoSF->GetBinContent(binx,1)+VarCoef*(RecoSF->GetBinError(binx,1));
    //An additional 1% systematic uncertainty should be applied to electrons with pT < 20 or > 80 GeV
    if ((*_lPt)[n]<20||(*_lPt)[n]>80) eleRecoSF+=0.01*VarCoef*(RecoSF->GetBinContent(binx,1));
    
    return eleSF*eleRecoSF;
    
}

//Function Calculating ele LepSF (FullSim-FastSim)
double FinalYield::CalcLepSF_ele_FastSim(TH2D* SF, int& n, const double& VarCoef){
    
    //Statistical uncertainty in hist2D neglected
    //2% uncertainty each leg
    int binx,biny;
    
    //LepSF
    //For current SF xbin=5(1-5), ybin=3(1-3)
    binx=SF->GetXaxis()->FindBin((*_lPt)[n]);
    if (binx<1) binx=1;
    if (binx>5) binx=5;//Control Overflow
    biny=SF->GetYaxis()->FindBin(fabs((*_lEta)[n]));
    if (biny<1) biny=1;
    if (biny>3) biny=3;//Control Overflow
    double eleSF=SF->GetBinContent(binx,biny)+VarCoef*0.02*(SF->GetBinContent(binx,biny));
    
    return eleSF;
    
}

//Function Calculating mu LepSF (Data-FullSim)
double FinalYield::CalcLepSF_mu_FullSim(TH2D* SF1, TH2D* SF2, TH2D* SF3, TH2D* SF4, TGraphAsymmErrors* RecoSF, int& n, const double& VarCoef){
    
    /*
     *Input map order: mu_dxydz_SF_FullSim, mu_MediumID_SF_FullSim, mu_miniIso_SF_FullSim, mu_LepMVAVM_SF_FullSim, RecoSF_FullSim
     */
    int binx,biny;
    
    //SF1: mu_dxydz_SF_FullSim
    //xbin:7(1-7)
    //ybin:4(1-4)
    binx=SF1->GetXaxis()->FindBin((*_lPt)[n]);
    if (binx<1) binx=1;
    if (binx>7) binx=7;//Control Overflow
    biny=SF1->GetYaxis()->FindBin(fabs((*_lEta)[n]));
    if (biny<1) biny=1;
    if (biny>4) biny=4;//Control Overflow
    double muSF1=SF1->GetBinContent(binx,biny)+VarCoef*(SF1->GetBinError(binx,biny));
    
    //SF2: mu_MediumID_SF_FullSim
    //xbin:7(1-7)
    //ybin:4(1-4)
    binx=SF2->GetXaxis()->FindBin((*_lPt)[n]);
    if (binx<1) binx=1;
    if (binx>7) binx=7;//Control Overflow
    biny=SF2->GetYaxis()->FindBin(fabs((*_lEta)[n]));
    if (biny<1) biny=1;
    if (biny>4) biny=4;//Control Overflow
    double muSF2=SF2->GetBinContent(binx,biny)+VarCoef*(SF2->GetBinError(binx,biny));
    
    //SF3: mu_miniIso_SF_FullSim
    //xbin:7(1-7)
    //ybin:4(1-4)
    binx=SF3->GetXaxis()->FindBin((*_lPt)[n]);
    if (binx<1) binx=1;
    if (binx>7) binx=7;//Control Overflow
    biny=SF3->GetYaxis()->FindBin(fabs((*_lEta)[n]));
    if (biny<1) biny=1;
    if (biny>4) biny=4;//Control Overflow
    double muSF3=SF3->GetBinContent(binx,biny)+VarCoef*(SF3->GetBinError(binx,biny));
    
    //SF4: mu_LepMVAVM_SF_FullSim
    //xbin:7(1-7)
    //ybin:4(1-4)
    binx=SF4->GetXaxis()->FindBin((*_lPt)[n]);
    if (binx<1) binx=1;
    if (binx>7) binx=7;//Control Overflow
    biny=SF4->GetYaxis()->FindBin(fabs((*_lEta)[n]));
    if (biny<1) biny=1;
    if (biny>4) biny=4;//Control Overflow
    double muSF4=SF4->GetBinContent(binx,biny)+VarCoef*(SF4->GetBinError(binx,biny));
    
    //RecoSF_FullSim: ratio_eff_eta3_dr030e030_corr (Pt>10)
    //RecoSF_FullSim: ratio_eff_eta3_tk0_dr030e030_corr (Pt<10)
    //Convert Graph Into x Array (15 data points)
    const int np=15;
    double x,y;
    double eta[np],etaerror[np],sf[np],sferror[np];
    for (int i=0; i<np; i++) {
        RecoSF->GetPoint(i,x,y);
        eta[i]=x;
        sf[i]=y;
        etaerror[i]=RecoSF->GetErrorX(i);
        sferror[i]=RecoSF->GetErrorY(i);
    }
    //Now find where the lepton eta falls
    int iup,index;
    for (int i=0; i<np; i++) {
        if ((*_lEta)[n]>eta[np-1]) {
            iup=np;
            break;
        }
        else if((*_lEta)[n]<=eta[i]){
            iup=i;
            break;
        }
    }
    //Falls between eta[i-1],eta[i]
    if (iup==0) {
        index=iup;
    }
    else if (iup==np){
        index=iup-1;
    }
    else {
        if ((*_lEta)[n]<=(eta[iup-1]+etaerror[iup-1])) index=iup-1;
        else index=iup;
    }
    //Calculate RecoSF
    double muRecoSF=sf[index]+VarCoef*sferror[index];
    
    return muSF1*muSF2*muSF3*muSF4*muRecoSF;
}

//Function Calculating mu LepSF (FullSim-FastSim)
double FinalYield::CalcLepSF_mu_FastSim(TH2D* SF1, TH2D* SF2, TH2D* SF3, TH2D* SF4, int& n, const double& VarCoef){
    
    /*
     *Input map order: mu_dxydz_SF_FastSim, mu_MediumID_SF_FastSim, mu_miniIso_SF_FastSim, mu_LepMVAVM_SF_FastSim
     */
    int binx,biny;
    
    //Statistical uncertainty in hist2D neglected
    //2% uncertainty each leg
    
    //SF1: mu_dxydz_SF_FastSim
    //xbin:5(1-5)
    //ybin:4(1-4)
    binx=SF1->GetXaxis()->FindBin((*_lPt)[n]);
    if (binx<1) binx=1;
    if (binx>5) binx=5;//Control Overflow
    biny=SF1->GetYaxis()->FindBin(fabs((*_lEta)[n]));
    if (biny<1) biny=1;
    if (biny>4) biny=4;//Control Overflow
    double muSF1=SF1->GetBinContent(binx,biny)+VarCoef*0.02*(SF1->GetBinContent(binx,biny));
    
    //SF2: mu_MediumID_SF_FastSim
    //xbin:5(1-5)
    //ybin:4(1-4)
    binx=SF2->GetXaxis()->FindBin((*_lPt)[n]);
    if (binx<1) binx=1;
    if (binx>5) binx=5;//Control Overflow
    biny=SF2->GetYaxis()->FindBin(fabs((*_lEta)[n]));
    if (biny<1) biny=1;
    if (biny>4) biny=4;//Control Overflow
    double muSF2=SF2->GetBinContent(binx,biny)+VarCoef*0.02*(SF2->GetBinContent(binx,biny));
    
    //SF3: mu_miniIso_SF_FastSim
    //xbin:5(1-5)
    //ybin:4(1-4)
    binx=SF3->GetXaxis()->FindBin((*_lPt)[n]);
    if (binx<1) binx=1;
    if (binx>5) binx=5;//Control Overflow
    biny=SF3->GetYaxis()->FindBin(fabs((*_lEta)[n]));
    if (biny<1) biny=1;
    if (biny>4) biny=4;//Control Overflow
    double muSF3=SF3->GetBinContent(binx,biny)+VarCoef*0.02*(SF3->GetBinContent(binx,biny));
    
    //SF4: mu_LepMVAVM_SF_FastSim
    //xbin:5(1-5)
    //ybin:4(1-4)
    binx=SF4->GetXaxis()->FindBin((*_lPt)[n]);
    if (binx<1) binx=1;
    if (binx>5) binx=5;//Control Overflow
    biny=SF4->GetYaxis()->FindBin(fabs((*_lEta)[n]));
    if (biny<1) biny=1;
    if (biny>4) biny=4;//Control Overflow
    double muSF4=SF4->GetBinContent(binx,biny)+VarCoef*0.02*(SF4->GetBinContent(binx,biny));
    
    return muSF1*muSF2*muSF3*muSF4;
}

//Function Calculating tau LepSF (Data-FullSim)
double FinalYield::CalcLepSF_tau_FullSim(const double& VarCoef){
    
    //0.95 with 5% uncertainty for 2016 Reco MiniAOD data
    double tauSF=0.95+VarCoef*0.05;
    return tauSF;
    
}

//Function Calculating tau LepSF (FullSim-FastSim)
double FinalYield::CalcLepSF_tau_FastSim(TH2D* SF, int& n, const double& VarCoef){
    
    int binx,biny;
    //tau_SF
    //xbin:5(1-5)
    //ybin:4(1-4)
    binx=SF->GetXaxis()->FindBin((*_lPt)[n]);
    if (binx<1) binx=1;
    if (binx>5) binx=5;//Control Overflow
    biny=SF->GetYaxis()->FindBin(fabs((*_lEta)[n]));
    if (biny<1) biny=1;
    if (biny>4) biny=4;//Control Overflow
    double tauSF=SF->GetBinContent(binx,biny)+VarCoef*(SF->GetBinError(binx,biny));
    
    return tauSF;
    
}

//Function Calculating SR Index
int FinalYield::CalcSRIndex_A(double& MET, double& Mll, double& MT){
    int index=9999;
    
    if ( 50<=MET && MET<100 && Mll<75 && MT>=0 && MT<100) index=1;
    if (100<=MET && MET<150 && Mll<75 && MT>=0 && MT<100) index=2;
    if (150<=MET && MET<200 && Mll<75 && MT>=0 && MT<100) index=3;
    if (200<=MET && MET<250 && Mll<75 && MT>=0 && MT<100) index=4;
    if (250<=MET &&            Mll<75 && MT>=0 && MT<100) index=5;
    
    if ( 50<=MET && MET<100 && Mll<75 && MT>=100 && MT<160) index=6;
    if (100<=MET && MET<150 && Mll<75 && MT>=100 && MT<160) index=7;
    if (150<=MET && MET<200 && Mll<75 && MT>=100 && MT<160) index=8;
    if (200<=MET &&            Mll<75 && MT>=100 && MT<160) index=9;
    
    if ( 50<=MET && MET<100 && Mll<75 && 160<=MT) index=10;
    if (100<=MET && MET<150 && Mll<75 && 160<=MT) index=11;
    if (150<=MET && MET<200 && Mll<75 && 160<=MT) index=12;
    if (200<=MET && MET<250 && Mll<75 && 160<=MT) index=13;
    if (250<=MET &&            Mll<75 && 160<=MT) index=14;
    
    if ( 50<=MET && MET<100 && Mll>=75 && Mll<105 && MT>=0 && MT<100) index=15;
    if (100<=MET && MET<150 && Mll>=75 && Mll<105 && MT>=0 && MT<100) index=16;
    if (150<=MET && MET<200 && Mll>=75 && Mll<105 && MT>=0 && MT<100) index=17;
    if (200<=MET && MET<250 && Mll>=75 && Mll<105 && MT>=0 && MT<100) index=18;
    if (250<=MET && MET<400 && Mll>=75 && Mll<105 && MT>=0 && MT<100) index=19;
    if (400<=MET && MET<550 && Mll>=75 && Mll<105 && MT>=0 && MT<100) index=20;
    if (550<=MET &&            Mll>=75 && Mll<105 && MT>=0 && MT<100) index=21;
    
    if ( 50<=MET && MET<100 && Mll>=75 && Mll<105 && MT>=100 && MT<160) index=22;
    if (100<=MET && MET<150 && Mll>=75 && Mll<105 && MT>=100 && MT<160) index=23;
    if (150<=MET && MET<200 && Mll>=75 && Mll<105 && MT>=100 && MT<160) index=24;
    if (200<=MET &&            Mll>=75 && Mll<105 && MT>=100 && MT<160) index=25;
    
    if ( 50<=MET && MET<100 && Mll>=75 && Mll<105 && 160<=MT) index=26;
    if (100<=MET && MET<150 && Mll>=75 && Mll<105 && 160<=MT) index=27;
    if (150<=MET && MET<200 && Mll>=75 && Mll<105 && 160<=MT) index=28;
    if (200<=MET && MET<250 && Mll>=75 && Mll<105 && 160<=MT) index=29;
    if (250<=MET && MET<400 && Mll>=75 && Mll<105 && 160<=MT) index=30;
    if (400<=MET &&            Mll>=75 && Mll<105 && 160<=MT) index=31;
    
    if ( 50<=MET && MET<100 && Mll>=105 && MT>=0 && MT<100) index=32;
    if (100<=MET && MET<150 && Mll>=105 && MT>=0 && MT<100) index=33;
    if (150<=MET && MET<200 && Mll>=105 && MT>=0 && MT<100) index=34;
    if (200<=MET && MET<250 && Mll>=105 && MT>=0 && MT<100) index=35;
    if (250<=MET &&            Mll>=105 && MT>=0 && MT<100) index=36;

    if ( 50<=MET && MET<100 && Mll>=105 && MT>=100 && MT<160) index=37;
    if (100<=MET && MET<150 && Mll>=105 && MT>=100 && MT<160) index=38;
    if (150<=MET && MET<200 && Mll>=105 && MT>=100 && MT<160) index=39;
    if (200<=MET &&            Mll>=105 && MT>=100 && MT<160) index=40;
    
    if ( 50<=MET && MET<100 && Mll>=105 && 160<=MT) index=41;
    if (100<=MET && MET<150 && Mll>=105 && 160<=MT) index=42;
    if (150<=MET && MET<200 && Mll>=105 && 160<=MT) index=43;
    if (200<=MET &&            Mll>=105 && 160<=MT) index=44;
    
    return index;
}

int FinalYield::CalcSRIndex_B(double& MET, double& Mll, double& MT){
    int index=9999;
    
    if (MT>=0 && MT <120) {
        if ( 50<=MET && MET<100 && Mll<100) index=1;
        if ( MET>=100           && Mll<100) index=2;
        if ( 50<=MET && MET<100 && Mll>=100) index=4;
        if ( MET>=100           && Mll>=100) index=5;
    }
    
    if (MT>=120) {
        if ( 50<=MET && Mll<100) index=3;
        if ( 50<=MET && Mll>=100) index=6;
    }
    
    return index;
}

int FinalYield::CalcSRIndex_C(double& MET, double& Mll, double& MT){
    int index=9999;
    
    if ( 50<=MET && MET<100 && Mll<75 && MT>=0 && MT<100) index=1;
    if (100<=MET && MET<150 && Mll<75 && MT>=0 && MT<100) index=2;
    if (150<=MET && MET<200 && Mll<75 && MT>=0 && MT<100) index=3;
    if (200<=MET && MET<250 && Mll<75 && MT>=0 && MT<100) index=4;
    if (250<=MET &&            Mll<75 && MT>=0 && MT<100) index=5;
    
    if ( 50<=MET && MET<200 && Mll<75 && MT>=100) index=17;
    if (200<=MET &&            Mll<75 && MT>=100) index=18;
    
    if ( 50<=MET && MET<100 && Mll>=75 && Mll<105 && MT>=0 && MT<100) index=6;
    if (100<=MET && MET<150 && Mll>=75 && Mll<105 && MT>=0 && MT<100) index=7;
    if (150<=MET && MET<200 && Mll>=75 && Mll<105 && MT>=0 && MT<100) index=8;
    if (200<=MET && MET<300 && Mll>=75 && Mll<105 && MT>=0 && MT<100) index=9;
    if (300<=MET && MET<400 && Mll>=75 && Mll<105 && MT>=0 && MT<100) index=10;
    if (400<=MET &&            Mll>=75 && Mll<105 && MT>=0 && MT<100) index=11;
    
    if ( 50<=MET && MET<100 && Mll>=75 && Mll<105 && MT>=100) index=6;
    if (100<=MET && MET<150 && Mll>=75 && Mll<105 && MT>=100) index=7;
    if (150<=MET && MET<200 && Mll>=75 && Mll<105 && MT>=100) index=8;
    if (200<=MET && MET<300 && Mll>=75 && Mll<105 && MT>=100) index=9;
    if (300<=MET && MET<400 && Mll>=75 && Mll<105 && MT>=100) index=10;
    if (400<=MET &&            Mll>=75 && Mll<105 && MT>=100) index=11;
    
    if ( 50<=MET && MET<100 && Mll>=105 && MT>=0 && MT<100) index=12;
    if (100<=MET && MET<150 && Mll>=105 && MT>=0 && MT<100) index=13;
    if (150<=MET && MET<200 && Mll>=105 && MT>=0 && MT<100) index=14;
    if (200<=MET && MET<250 && Mll>=105 && MT>=0 && MT<100) index=15;
    if (250<=MET &&            Mll>=105 && MT>=0 && MT<100) index=16;
    
    if ( 50<=MET && MET<200 && Mll>=105 && MT>=100) index=17;
    if (200<=MET &&            Mll>=105 && MT>=100) index=18;
    
    return index;
}

int FinalYield::CalcSRIndex_D(double& MET, double& Mll, double& MT){
    int index=9999;
    
    if ( 50<=MET && MET<100 && Mll<60 && MT>=0 && MT<100) index=1;
    if (100<=MET && MET<150 && Mll<60 && MT>=0 && MT<100) index=2;
    if (150<=MET && MET<200 && Mll<60 && MT>=0 && MT<100) index=3;
    if (200<=MET && MET<250 && Mll<60 && MT>=0 && MT<100) index=4;
    if (250<=MET &&            Mll<60 && MT>=0 && MT<100) index=5;
    
    if ( 50<=MET && MET<100 && Mll>=60 && Mll<100 && MT>=0 && MT<100) index=6;
    if (100<=MET && MET<150 && Mll>=60 && Mll<100 && MT>=0 && MT<100) index=7;
    if (150<=MET && MET<200 && Mll>=60 && Mll<100 && MT>=0 && MT<100) index=8;
    if (200<=MET && MET<250 && Mll>=60 && Mll<100 && MT>=0 && MT<100) index=9;
    if (250<=MET &&            Mll>=60 && Mll<100 && MT>=0 && MT<100) index=10;
    
    if ( 50<=MET && MET<100 && Mll>=100 && MT>=0 && MT<100) index=11;
    if (100<=MET && MET<150 && Mll>=100 && MT>=0 && MT<100) index=12;
    if (150<=MET && MET<200 && Mll>=100 && MT>=0 && MT<100) index=13;
    if (200<=MET &&            Mll>=100 && MT>=0 && MT<100) index=14;
    
    if ( 50<=MET && MET<200  && MT>=100) index=15;
    if (200<=MET &&             MT>=100) index=16;
    
    return index;
}

int FinalYield::CalcSRIndex_E(double& MET, double& Mll, double& MT){
    int index=9999;
    
    if ( 50<=MET && MET<100 && Mll<60 && MT>=0 && MT<100) index=1;
    if (100<=MET && MET<150 && Mll<60 && MT>=0 && MT<100) index=2;
    if (150<=MET && MET<200 && Mll<60 && MT>=0 && MT<100) index=3;
    if (200<=MET && MET<250 && Mll<60 && MT>=0 && MT<100) index=4;
    if (250<=MET &&            Mll<60 && MT>=0 && MT<100) index=5;
    
    if ( 50<=MET && MET<100 && Mll>=60 && Mll<100 && MT>=0 && MT<100) index=6;
    if (100<=MET && MET<150 && Mll>=60 && Mll<100 && MT>=0 && MT<100) index=7;
    if (150<=MET && MET<200 && Mll>=60 && Mll<100 && MT>=0 && MT<100) index=8;
    if (200<=MET && MET<250 && Mll>=60 && Mll<100 && MT>=0 && MT<100) index=9;
    if (250<=MET &&            Mll>=60 && Mll<100 && MT>=0 && MT<100) index=10;
    
    if (Mll>=100 && MT>=0 && MT<100) index=11;
    
    if (50<=MET && MT>=100) index=12;
    
    return index;
}

int FinalYield::CalcSRIndex_F(double& MET, double& Mll, double& MT){
    int index=9999;
    
    if ( 50<=MET && MET<100 && Mll<100 && MT>=0 && MT<100) index=1;
    if (100<=MET && MET<150 && Mll<100 && MT>=0 && MT<100) index=2;
    if (150<=MET && MET<200 && Mll<100 && MT>=0 && MT<100) index=3;
    if (200<=MET && MET<250 && Mll<100 && MT>=0 && MT<100) index=4;
    if (250<=MET && MET<300 && Mll<100 && MT>=0 && MT<100) index=5;
    if (300<=MET &&            Mll<100 && MT>=0 && MT<100) index=6;
    
    if ( 50<=MET && MET<100 && Mll>=100 && MT>=0 && MT<100) index=7;
    if (100<=MET && MET<150 && Mll>=100 && MT>=0 && MT<100) index=8;
    if (150<=MET && MET<200 && Mll>=100 && MT>=0 && MT<100) index=9;
    if (200<=MET &&            Mll>=100 && MT>=0 && MT<100) index=10;
    
    if ( 50<=MET && MET<200  && MT>=100) index=11;
    if (200<=MET &&             MT>=100) index=12;
    
    return index;
}

//Assign Uncertainty
void FinalYield::AssignUncForPlot_SRBG_WZ(TH1F* Hist, int& nbin, TH1F* HistUp1, TH1F* HistDown1, TH1F* HistUp2, TH1F* HistDown2, TH1F* HistUp3, TH1F* HistDown3, TH1F* HistUp4, TH1F* HistDown4, TH1F* HistStatUnc, TH1F* HistTotalUnc){
    
    //SR A
    if (nbin==44) {
        //Statistical uncertainty + up/down uncertainty + shape uncertainty + 9% normalization uncertainty
        double ShapeUnc[44]={
            0.05,0.10,0.20,0.20,0.20,0.05,0.25,0.45,0.45,0.30,0.50,0.40,0.45,0.45,
            0.05,0.10,0.20,0.20,0.20,0.20,0.20,0.05,0.25,0.45,0.45,0.30,0.50,0.40,0.45,0.45,0.45,
            0.05,0.10,0.20,0.20,0.20,0.05,0.25,0.45,0.45,0.30,0.50,0.40,0.45
        };
        for (int i=1; i<=nbin; i++) {
            
            double StatUnc=Hist->GetBinError(i);
            double BinVal=Hist->GetBinContent(i);
            
            double CumUncSq=0;
            //Stat Unc
            CumUncSq=CumUncSq+StatUnc*StatUnc;
            //IntLumi Unc: 2.5%
            CumUncSq=CumUncSq+(0.025*BinVal)*(0.025*BinVal);
            //Trigger Eff Unc: 3%
            CumUncSq=CumUncSq+(0.03*BinVal)*(0.03*BinVal);
            //Up/Down Unc:PU/JEC/SF/b-Tagging
            double AveUnc1=0.5*(abs(HistUp1->GetBinContent(i) - BinVal)+abs(HistDown1->GetBinContent(i) - BinVal));
            double AveUnc2=0.5*(abs(HistUp2->GetBinContent(i) - BinVal)+abs(HistDown2->GetBinContent(i) - BinVal));
            double AveUnc3=0.5*(abs(HistUp3->GetBinContent(i) - BinVal)+abs(HistDown3->GetBinContent(i) - BinVal));
            double AveUnc4=0.5*(abs(HistUp4->GetBinContent(i) - BinVal)+abs(HistDown4->GetBinContent(i) - BinVal));
            CumUncSq=CumUncSq+AveUnc1*AveUnc1+AveUnc2*AveUnc2+AveUnc3*AveUnc3+AveUnc4*AveUnc4;
            //Shape and normalization
            double TotalUnc=sqrt(CumUncSq+(0.09*BinVal)*(0.09*BinVal)+(ShapeUnc[i-1]*BinVal)*(ShapeUnc[i-1]*BinVal));
            //Assign uncertainty
            HistStatUnc->SetBinError(i,StatUnc);
            HistTotalUnc->SetBinError(i,TotalUnc);
            Hist->SetBinError(i,TotalUnc);
        }
    }
    //SR B-F
    else {
        //Statistical uncertainty + up/down uncertainty + 9% normalization uncertainty
        for (int i=1; i<=nbin; i++) {
            
            double StatUnc=Hist->GetBinError(i);
            double BinVal=Hist->GetBinContent(i);
            
            double CumUncSq=0;
            //Stat Unc
            CumUncSq=CumUncSq+StatUnc*StatUnc;
            //IntLumi Unc: 2.5%
            CumUncSq=CumUncSq+(0.025*BinVal)*(0.025*BinVal);
            //Trigger Eff Unc: 3%
            CumUncSq=CumUncSq+(0.03*BinVal)*(0.03*BinVal);
            //Up/Down Unc:PU/JEC/SF/b-Tagging
            double AveUnc1=0.5*(abs(HistUp1->GetBinContent(i) - BinVal)+abs(HistDown1->GetBinContent(i) - BinVal));
            double AveUnc2=0.5*(abs(HistUp2->GetBinContent(i) - BinVal)+abs(HistDown2->GetBinContent(i) - BinVal));
            double AveUnc3=0.5*(abs(HistUp3->GetBinContent(i) - BinVal)+abs(HistDown3->GetBinContent(i) - BinVal));
            double AveUnc4=0.5*(abs(HistUp4->GetBinContent(i) - BinVal)+abs(HistDown4->GetBinContent(i) - BinVal));
            CumUncSq=CumUncSq+AveUnc1*AveUnc1+AveUnc2*AveUnc2+AveUnc3*AveUnc3+AveUnc4*AveUnc4;
            //Shape and normalization
            double TotalUnc=sqrt(CumUncSq+(0.09*BinVal)*(0.09*BinVal));
            //Assign uncertainty
            HistStatUnc->SetBinError(i,StatUnc);
            HistTotalUnc->SetBinError(i,TotalUnc);
            Hist->SetBinError(i,TotalUnc);
        }
    }
    return;
}

void FinalYield::AssignUncForPlot_SRBG_NonWZ(TH1F* Hist, double uncertainty, int& nbin, TH1F* HistUp1, TH1F* HistDown1, TH1F* HistUp2, TH1F* HistDown2, TH1F* HistUp3, TH1F* HistDown3, TH1F* HistUp4, TH1F* HistDown4, TH1F* HistStatUnc, TH1F* HistTotalUnc){
    
    //Statistical uncertainty + up/down uncertainty + process-specific uncertainty
    for (int i=1; i<=nbin; i++) {
        
        double StatUnc=Hist->GetBinError(i);
        double BinVal=Hist->GetBinContent(i);
        
        double CumUncSq=0;
        //Stat Unc
        CumUncSq=CumUncSq+StatUnc*StatUnc;
        //IntLumi Unc: 2.5%
        CumUncSq=CumUncSq+(0.025*BinVal)*(0.025*BinVal);
        //Trigger Eff Unc: 3%
        CumUncSq=CumUncSq+(0.03*BinVal)*(0.03*BinVal);
        //Up/Down Unc:PU/JEC/SF/b-Tagging
        double AveUnc1=0.5*(abs(HistUp1->GetBinContent(i) - BinVal)+abs(HistDown1->GetBinContent(i) - BinVal));
        double AveUnc2=0.5*(abs(HistUp2->GetBinContent(i) - BinVal)+abs(HistDown2->GetBinContent(i) - BinVal));
        double AveUnc3=0.5*(abs(HistUp3->GetBinContent(i) - BinVal)+abs(HistDown3->GetBinContent(i) - BinVal));
        double AveUnc4=0.5*(abs(HistUp4->GetBinContent(i) - BinVal)+abs(HistDown4->GetBinContent(i) - BinVal));
        CumUncSq=CumUncSq+AveUnc1*AveUnc1+AveUnc2*AveUnc2+AveUnc3*AveUnc3+AveUnc4*AveUnc4;
        //Shape and normalization
        double TotalUnc=sqrt(CumUncSq+(uncertainty*BinVal)*(uncertainty*BinVal));
        //Assign uncertainty
        HistStatUnc->SetBinError(i,StatUnc);
        HistTotalUnc->SetBinError(i,TotalUnc);
        Hist->SetBinError(i,TotalUnc);
    }
    return;
}

void FinalYield::AssignUncForPlot_VarBG_WZ(TH1F* Hist, double uncertainty, int& nbin, TH1F* HistStatUnc, TH1F* HistTotalUnc){
    
    //Statistical uncertainty + up/down uncertainty + process-specific uncertainty
    for (int i=1; i<=nbin; i++) {
        
        double StatUnc=Hist->GetBinError(i);
        double BinVal=Hist->GetBinContent(i);
        
        double CumUncSq=0;
        //Stat Unc
        CumUncSq=CumUncSq+StatUnc*StatUnc;
        //IntLumi Unc: 2.5%
        CumUncSq=CumUncSq+(0.025*BinVal)*(0.025*BinVal);
        //Trigger Eff Unc: 3%
        CumUncSq=CumUncSq+(0.03*BinVal)*(0.03*BinVal);
        //Up/Down Unc
        //PU Unc: 3%
        CumUncSq=CumUncSq+(0.03*BinVal)*(0.03*BinVal);
        //JEC Unc: 5%
        CumUncSq=CumUncSq+(0.05*BinVal)*(0.05*BinVal);
        //SF Unc: 4%
        CumUncSq=CumUncSq+(0.04*BinVal)*(0.04*BinVal);
        //b-Tagging Unc: 5%
        CumUncSq=CumUncSq+(0.05*BinVal)*(0.05*BinVal);
        
        //Shape and normalization
        double TotalUnc=sqrt(CumUncSq+(uncertainty*BinVal)*(uncertainty*BinVal)+(0.09*BinVal)*(0.09*BinVal));
        //Assign uncertainty
        HistStatUnc->SetBinError(i,StatUnc);
        HistTotalUnc->SetBinError(i,TotalUnc);
        Hist->SetBinError(i,TotalUnc);
    }
    return;
}

void FinalYield::AssignUncForPlot_VarBG_NonWZ(TH1F* Hist, double uncertainty, int& nbin, TH1F* HistStatUnc, TH1F* HistTotalUnc){
    
    //Statistical uncertainty + up/down uncertainty + process-specific uncertainty + 9% normalization uncertainty
    for (int i=1; i<=nbin; i++) {
        
        double StatUnc=Hist->GetBinError(i);
        double BinVal=Hist->GetBinContent(i);
        
        double CumUncSq=0;
        //Stat Unc
        CumUncSq=CumUncSq+StatUnc*StatUnc;
        //IntLumi Unc: 2.5%
        CumUncSq=CumUncSq+(0.025*BinVal)*(0.025*BinVal);
        //Trigger Eff Unc: 3%
        CumUncSq=CumUncSq+(0.03*BinVal)*(0.03*BinVal);
        //Up/Down Unc
        //PU Unc: 3%
        CumUncSq=CumUncSq+(0.03*BinVal)*(0.03*BinVal);
        //JEC Unc: 5%
        CumUncSq=CumUncSq+(0.05*BinVal)*(0.05*BinVal);
        //SF Unc: 4%
        CumUncSq=CumUncSq+(0.04*BinVal)*(0.04*BinVal);
        //b-Tagging Unc: 5%
        CumUncSq=CumUncSq+(0.05*BinVal)*(0.05*BinVal);
        
        //Shape and normalization
        double TotalUnc=sqrt(CumUncSq+(uncertainty*BinVal)*(uncertainty*BinVal));
        //Assign uncertainty
        HistStatUnc->SetBinError(i,StatUnc);
        HistTotalUnc->SetBinError(i,TotalUnc);
        Hist->SetBinError(i,TotalUnc);
    }
    return;
}
