#define FinalYield_Signal_cxx
#include "FinalYield_Signal.h"
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


//Set up pars
const int nSRCat=6;
const int nbCNMassPoint=77;//100-2000:25
const int nbLSPMassPoint=81;//0-2000:25

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

const double intlumi=35900;

/*
 *Full Signal Sets
 */
const int nbMCsamples = 9;//9 Signal samples
TString samplefilesMC[nbMCsamples] = {
    "TChiSlepSnu_tauenriched_x0p05-skim3Lreco.root",
    "TChiSlepSnu_tauenriched_x0p5-skim3Lreco.root",
    "TChiSlepSnu_tauenriched_x0p95-skim3Lreco.root",
    "TChiSlepSnu_x0p05-skim3Lreco.root",
    "TChiSlepSnu_x0p5-skim3Lreco.root",
    "TChiSlepSnu_x0p95-skim3Lreco.root",
    "TChiStauStau_x0p5-skim3Lreco.root",
    "TChiWH_WToLNu_HToVVTauTau-skim3Lreco.root",
    "TChiWZ_ZToLL_mZMin-0p1-skim3Lreco.root"
};

TString SignalProcess[nbMCsamples]={
    "TChiSlepSnu_tauenriched_x0p05",
    "TChiSlepSnu_tauenriched_x0p5",
    "TChiSlepSnu_tauenriched_x0p95",
    "TChiSlepSnu_x0p05",
    "TChiSlepSnu_x0p5",
    "TChiSlepSnu_x0p95",
    "TChiStauStau_x0p5",
    "TChiWH_WToLNu_HToVVTauTau",
    "TChiWZ_ZToLL_mZMin-0p1"
};

//Branching Fraction
/*
 *Flavor-democratic: 50%
 *Tau-enriched: 100%
 *Tau-dominated: 100%
 *WH->3l: 2.9%
 *WZ->3l: 3.3%
 */
//const double BR[nbMCsamples]={1,1,1,0.5,0.5,0.5,1,0.029,0.033};
const double BR[nbMCsamples]={1,1,1,0.5,0.5,0.5,0.5,0.029,0.100};//For limit setting
/*
 *Full Signal Sets
 */

void FinalYield_Signal::Loop()
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
        MCVtx[iSam]=new TFile("/cms/data/store/user/t2/users/takeimai/output/SRYield/2016Signal/"+samplefilesMC[iSam],"open");
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
    SFmap=new TFile("/cms/data/store/user/t2/users/takeimai/CMSSW_8_0_13/src/SUSYAnalyzer/PatAnalyzer/test/Macro/ScaleFactor/el_SF_FastSim.root","open");
    TH2D* el_SF_FastSim = (TH2D*)SFmap->Get("histo2D");
    SFmap=new TFile("/cms/data/store/user/t2/users/takeimai/CMSSW_8_0_13/src/SUSYAnalyzer/PatAnalyzer/test/Macro/ScaleFactor/el_RecoSF_FullSim.root","open");
    TH2D* el_RecoSF_FullSim = (TH2D*)SFmap->Get("EGamma_SF2D");
    //Muon
    SFmap=new TFile("/cms/data/store/user/t2/users/takeimai/CMSSW_8_0_13/src/SUSYAnalyzer/PatAnalyzer/test/Macro/ScaleFactor/mu_dxydz_SF_FullSim.root","open");
    TH2D* mu_dxydz_SF_FullSim = (TH2D*)SFmap->Get("SF");
    SFmap=new TFile("/cms/data/store/user/t2/users/takeimai/CMSSW_8_0_13/src/SUSYAnalyzer/PatAnalyzer/test/Macro/ScaleFactor/mu_dxydz_SF_FastSim.root","open");
    TH2D* mu_dxydz_SF_FastSim = (TH2D*)SFmap->Get("histo2D");
    
    SFmap=new TFile("/cms/data/store/user/t2/users/takeimai/CMSSW_8_0_13/src/SUSYAnalyzer/PatAnalyzer/test/Macro/ScaleFactor/mu_MediumID_SF_FullSim.root","open");
    TH2D* mu_MediumID_SF_FullSim = (TH2D*)SFmap->Get("SF");
    SFmap=new TFile("/cms/data/store/user/t2/users/takeimai/CMSSW_8_0_13/src/SUSYAnalyzer/PatAnalyzer/test/Macro/ScaleFactor/mu_MediumID_SF_FastSim.root","open");
    TH2D* mu_MediumID_SF_FastSim = (TH2D*)SFmap->Get("histo2D");
    
    SFmap=new TFile("/cms/data/store/user/t2/users/takeimai/CMSSW_8_0_13/src/SUSYAnalyzer/PatAnalyzer/test/Macro/ScaleFactor/mu_miniIso_SF_FullSim.root","open");
    TH2D* mu_miniIso_SF_FullSim = (TH2D*)SFmap->Get("SF");
    SFmap=new TFile("/cms/data/store/user/t2/users/takeimai/CMSSW_8_0_13/src/SUSYAnalyzer/PatAnalyzer/test/Macro/ScaleFactor/mu_miniIso_SF_FastSim.root","open");
    TH2D* mu_miniIso_SF_FastSim = (TH2D*)SFmap->Get("histo2D");
    
    
    SFmap=new TFile("/cms/data/store/user/t2/users/takeimai/CMSSW_8_0_13/src/SUSYAnalyzer/PatAnalyzer/test/Macro/ScaleFactor/mu_LepMVAVM_SF_FullSim.root","open");
    TH2D* mu_LepMVAVM_SF_FullSim = (TH2D*)SFmap->Get("SF");
    SFmap=new TFile("/cms/data/store/user/t2/users/takeimai/CMSSW_8_0_13/src/SUSYAnalyzer/PatAnalyzer/test/Macro/ScaleFactor/mu_LepMVAVM_SF_FastSim.root","open");
    TH2D* mu_LepMVAVM_SF_FastSim = (TH2D*)SFmap->Get("histo2D");
    
    SFmap=new TFile("/cms/data/store/user/t2/users/takeimai/CMSSW_8_0_13/src/SUSYAnalyzer/PatAnalyzer/test/Macro/ScaleFactor/mu_RecoSF_FullSim.root","open");
    TGraphAsymmErrors* mu_RecoSF_FullSim_Highpt = (TGraphAsymmErrors*)SFmap->Get("ratio_eff_eta3_dr030e030_corr");//Pt>10
    TGraphAsymmErrors* mu_RecoSF_FullSim_Lowpt = (TGraphAsymmErrors*)SFmap->Get("ratio_eff_eta3_tk0_dr030e030_corr");//Pt<10
    //Tau
    //No input for FullSim tau Lep SF
    SFmap=new TFile("/cms/data/store/user/t2/users/takeimai/CMSSW_8_0_13/src/SUSYAnalyzer/PatAnalyzer/test/Macro/ScaleFactor/tau_SF_FastSim.root","open");
    TH2D* tau_SF_FastSim = (TH2D*)SFmap->Get("histo2D");
    //////////////////////////////
    ////Set up leoton SF//////////
    //////////////////////////////
    
    //Set up histograms par
    TString name;
    TString SRCat[nSRCat]={
        "SR_A",
        "SR_B",
        "SR_C",
        "SR_D",
        "SR_E",
        "SR_F"
    };
    
    //Histogram per model per mass combination in each SR
    //One root output per sample
    TH1F* MCYield_Signal[nSRCat][nbCNMassPoint*nbLSPMassPoint];
    TH1F* MCYield_Signal_JECUp[nSRCat][nbCNMassPoint*nbLSPMassPoint];
    TH1F* MCYield_Signal_JECDown[nSRCat][nbCNMassPoint*nbLSPMassPoint];
    TH1F* MCYield_Signal_BTagUp[nSRCat][nbCNMassPoint*nbLSPMassPoint];
    TH1F* MCYield_Signal_BTagDown[nSRCat][nbCNMassPoint*nbLSPMassPoint];
    TH1F* MCYield_Signal_LepSFUp[nSRCat][nbCNMassPoint*nbLSPMassPoint];
    TH1F* MCYield_Signal_LepSFDown[nSRCat][nbCNMassPoint*nbLSPMassPoint];
    TH1F* MCYield_Signal_PUUp[nSRCat][nbCNMassPoint*nbLSPMassPoint];
    TH1F* MCYield_Signal_PUDown[nSRCat][nbCNMassPoint*nbLSPMassPoint];
    
    TH1F* MCYield_Signal_Pt1[nSRCat][nbCNMassPoint*nbLSPMassPoint];
    TH1F* MCYield_Signal_Pt2[nSRCat][nbCNMassPoint*nbLSPMassPoint];
    TH1F* MCYield_Signal_Pt3[nSRCat][nbCNMassPoint*nbLSPMassPoint];
    TH1F* MCYield_Signal_Mll[nSRCat][nbCNMassPoint*nbLSPMassPoint];
    TH1F* MCYield_Signal_MT[nSRCat][nbCNMassPoint*nbLSPMassPoint];
    TH1F* MCYield_Signal_MET[nSRCat][nbCNMassPoint*nbLSPMassPoint];
    TH1F* MCYield_Signal_M3l[nSRCat][nbCNMassPoint*nbLSPMassPoint];
    TH1F* MCYield_Signal_SumQ[nSRCat][nbCNMassPoint*nbLSPMassPoint];
    
    //Set CSV files for event counter
    //ofstream EventStat;
    ifstream EventStat;
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
        
        
        //inputMC[iSam] = new TFile("/cms/data/store/user/t2/users/takeimai/output/SRYield/2016Signal/"+samplefilesMC[iSam],"open");
        //inputMC[iSam]->cd("FakeElectrons");
        /*
        TH1F* eventcounter=new TH1F("hCounter","",5,0,5);
        eventcounter->SetDirectory(0);
        eventcounter->Read("hCounter");
        int nbevents=eventcounter->GetEntries();//total nb of events
         */
        //Read model event counter
        /*
        int EventPerModel[nbCNMassPoint][nbLSPMassPoint]={0};
        TH1F* ModelCounter[nbCNMassPoint][nbLSPMassPoint];
        EventStat.open(Form("TotalEventCount_Process_%d.csv",iSam));
        for (int i=0; i<nbCNMassPoint; i++) for (int j=0; j<nbLSPMassPoint; j++){
            ModelCounter[i][j]=new TH1F(Form("ModelCounter_%d_%d",i,j),"",5,0,5);
            ModelCounter[i][j]->SetDirectory(0);
            ModelCounter[i][j]->Read(Form("ModelCounter_%d_%d",i,j));
            EventPerModel[i][j]=ModelCounter[i][j]->GetEntries();
            EventStat<<i<<","<<j<<","<<EventPerModel[i][j]<<endl;
            //ModelCounter[i][j]->Delete();
        }
        EventStat.close();
        */
        //Reading Ends
        //inputMC[iSam]->Close();
        
        //Read model event counter from csv
        int EventPerModel[nbCNMassPoint][nbLSPMassPoint]={0};
        EventStat.open("SigModelCount/"+SignalProcess[iSam]+".csv");
        string ReadLine,ReadColumn;
        if (EventStat.is_open()){
            while (getline(EventStat,ReadLine)){
                stringstream LineStream(ReadLine);
                int count=0;
                int ReadOut[3]={0};
                while(getline(LineStream,ReadColumn,',')){
                    ReadOut[count]=stoi(ReadColumn);
                    count++;
                }
                EventPerModel[(ReadOut[0]-100)/25][ReadOut[1]/25]=ReadOut[2];
            }
        }
        EventStat.close();
        //Reading Ends
        
        /*
         *Set Up Histograms
         */
        //Histogram per model per mass combination in each SR
        for (int i=0; i<nSRCat; i++) for (int k=0; k<nbCNMassPoint; k++) for (int l=0; l<nbLSPMassPoint; l++){
            
            if (EventPerModel[k][l]==0) continue;
            
            name="Signal_"+SignalProcess[iSam]+"_"+SRCat[i]+ Form("_CN%d_LSP%d",100+k*25,l*25);
            MCYield_Signal[i][k*nbLSPMassPoint+l]=new TH1F(name,"",nSR[i],0.5,nSR[i]+0.5);
            MCYield_Signal[i][k*nbLSPMassPoint+l]->Sumw2();
            name="Signal_JECUp_"+SignalProcess[iSam]+"_"+SRCat[i]+ Form("_CN%d_LSP%d",100+k*25,l*25);
            MCYield_Signal_JECUp[i][k*nbLSPMassPoint+l]=new TH1F(name,"",nSR[i],0.5,nSR[i]+0.5);
            MCYield_Signal_JECUp[i][k*nbLSPMassPoint+l]->Sumw2();
            name="Signal_JECDown_"+SignalProcess[iSam]+"_"+SRCat[i]+ Form("_CN%d_LSP%d",100+k*25,l*25);
            MCYield_Signal_JECDown[i][k*nbLSPMassPoint+l]=new TH1F(name,"",nSR[i],0.5,nSR[i]+0.5);
            MCYield_Signal_JECDown[i][k*nbLSPMassPoint+l]->Sumw2();
            name="Signal_BTagUp_"+SignalProcess[iSam]+"_"+SRCat[i]+ Form("_CN%d_LSP%d",100+k*25,l*25);
            MCYield_Signal_BTagUp[i][k*nbLSPMassPoint+l]=new TH1F(name,"",nSR[i],0.5,nSR[i]+0.5);
            MCYield_Signal_BTagUp[i][k*nbLSPMassPoint+l]->Sumw2();
            name="Signal_BTagDown_"+SignalProcess[iSam]+"_"+SRCat[i]+ Form("_CN%d_LSP%d",100+k*25,l*25);
            MCYield_Signal_BTagDown[i][k*nbLSPMassPoint+l]=new TH1F(name,"",nSR[i],0.5,nSR[i]+0.5);
            MCYield_Signal_BTagDown[i][k*nbLSPMassPoint+l]->Sumw2();
            name="Signal_LepSFUp_"+SignalProcess[iSam]+"_"+SRCat[i]+ Form("_CN%d_LSP%d",100+k*25,l*25);
            MCYield_Signal_LepSFUp[i][k*nbLSPMassPoint+l]=new TH1F(name,"",nSR[i],0.5,nSR[i]+0.5);
            MCYield_Signal_LepSFUp[i][k*nbLSPMassPoint+l]->Sumw2();
            name="Signal_LepSFDown_"+SignalProcess[iSam]+"_"+SRCat[i]+ Form("_CN%d_LSP%d",100+k*25,l*25);
            MCYield_Signal_LepSFDown[i][k*nbLSPMassPoint+l]=new TH1F(name,"",nSR[i],0.5,nSR[i]+0.5);
            MCYield_Signal_LepSFDown[i][k*nbLSPMassPoint+l]->Sumw2();
            name="Signal_PUUp_"+SignalProcess[iSam]+"_"+SRCat[i]+ Form("_CN%d_LSP%d",100+k*25,l*25);
            MCYield_Signal_PUUp[i][k*nbLSPMassPoint+l]=new TH1F(name,"",nSR[i],0.5,nSR[i]+0.5);
            MCYield_Signal_PUUp[i][k*nbLSPMassPoint+l]->Sumw2();
            name="Signal_PUDown_"+SignalProcess[iSam]+"_"+SRCat[i]+ Form("_CN%d_LSP%d",100+k*25,l*25);
            MCYield_Signal_PUDown[i][k*nbLSPMassPoint+l]=new TH1F(name,"",nSR[i],0.5,nSR[i]+0.5);
            MCYield_Signal_PUDown[i][k*nbLSPMassPoint+l]->Sumw2();
            
            name="Signal_Pt1_"+SignalProcess[iSam]+"_"+SRCat[i]+ Form("_CN%d_LSP%d",100+k*25,l*25);
            MCYield_Signal_Pt1[i][k*nbLSPMassPoint+l]=new TH1F(name,"",Ptbin[i],Ptlowerbound,Ptupperbound[i]);
            MCYield_Signal_Pt1[i][k*nbLSPMassPoint+l]->Sumw2();
            name="Signal_Pt2_"+SignalProcess[iSam]+"_"+SRCat[i]+ Form("_CN%d_LSP%d",100+k*25,l*25);
            MCYield_Signal_Pt2[i][k*nbLSPMassPoint+l]=new TH1F(name,"",Ptbin[i],Ptlowerbound,Ptupperbound[i]);
            MCYield_Signal_Pt2[i][k*nbLSPMassPoint+l]->Sumw2();
            name="Signal_Pt3_"+SignalProcess[iSam]+"_"+SRCat[i]+ Form("_CN%d_LSP%d",100+k*25,l*25);
            MCYield_Signal_Pt3[i][k*nbLSPMassPoint+l]=new TH1F(name,"",Ptbin[i],Ptlowerbound,Ptupperbound[i]);
            MCYield_Signal_Pt3[i][k*nbLSPMassPoint+l]->Sumw2();
            name="Signal_Mll_"+SignalProcess[iSam]+"_"+SRCat[i]+ Form("_CN%d_LSP%d",100+k*25,l*25);
            MCYield_Signal_Mll[i][k*nbLSPMassPoint+l]=new TH1F(name,"",Mllbin[i],Mlllowerbound,Mllupperbound[i]);
            MCYield_Signal_Mll[i][k*nbLSPMassPoint+l]->Sumw2();
            name="Signal_MT_"+SignalProcess[iSam]+"_"+SRCat[i]+ Form("_CN%d_LSP%d",100+k*25,l*25);
            MCYield_Signal_MT[i][k*nbLSPMassPoint+l]=new TH1F(name,"",MTbin[i],MTlowerbound,MTupperbound[i]);
            MCYield_Signal_MT[i][k*nbLSPMassPoint+l]->Sumw2();
            name="Signal_MET_"+SignalProcess[iSam]+"_"+SRCat[i]+ Form("_CN%d_LSP%d",100+k*25,l*25);
            MCYield_Signal_MET[i][k*nbLSPMassPoint+l]=new TH1F(name,"",METbin[i],METlowerbound,METupperbound[i]);
            MCYield_Signal_MET[i][k*nbLSPMassPoint+l]->Sumw2();
            name="Signal_M3l_"+SignalProcess[iSam]+"_"+SRCat[i]+ Form("_CN%d_LSP%d",100+k*25,l*25);
            MCYield_Signal_M3l[i][k*nbLSPMassPoint+l]=new TH1F(name,"",M3lbin[i],M3llowerbound,M3lupperbound[i]);
            MCYield_Signal_M3l[i][k*nbLSPMassPoint+l]->Sumw2();
            name="Signal_SumQ_"+SignalProcess[iSam]+"_"+SRCat[i]+ Form("_CN%d_LSP%d",100+k*25,l*25);
            MCYield_Signal_SumQ[i][k*nbLSPMassPoint+l]=new TH1F(name,"",SumQbin[i],SumQlowerbound,SumQupperbound[i]);
            MCYield_Signal_SumQ[i][k*nbLSPMassPoint+l]->Sumw2();
        }
        /*
         *
         */
        
        inputMC[iSam] = new TFile("/cms/data/store/user/t2/users/takeimai/output/SRYield/2016Signal/"+samplefilesMC[iSam],"open");
        TTree *thetree_mc = (TTree*)(inputMC[iSam])->Get("FakeElectrons/fakeTree");
        Init(thetree_mc);
        Long64_t nentries_mc = (*thetree_mc).GetEntries();//Get the number of generated events
        
        cout<<"Processing Signal sample "<<samplefilesMC[iSam]<<endl;
        //cout<<"Total number of generated events: "<<nbevents<<endl;
        
        for (Long64_t jentry=0; jentry<nentries_mc;jentry++) {//loop over all the events//First loop
            
            LoadTree(jentry);
            thetree_mc->GetEntry(jentry);
            
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
            //X-section of signal samples: from mchargino/mneutralino
            //Normalization on different mass point
            //Count total events per model per chargino/neutralino2&LSP mass combination
            int MassIndexCN=FindCNIndex(_mneutralino2);
            int MassIndexLSP=FindLSPIndex(_mLSP);
            int nbevents=EventPerModel[MassIndexCN][MassIndexLSP];
            double xsection=FindSignalXSection(_mchargino,_mneutralino2)*BR[iSam];
            cout<<"Event: "<<jentry<<" "<<"mCN: "<<_mneutralino2<<" "<<MassIndexCN<<" "<<"mLSP: "<<_mLSP<<" "<<MassIndexLSP<<" "<<EventPerModel[MassIndexCN][MassIndexLSP]<<endl;
            if (nbevents==0) continue;
            double W = _weight*intlumi*xsection/nbevents*reweightSF[iSam][trueNVtx];
            double W_PUUp = _weight*intlumi*xsection/nbevents*reweightSF_Up[iSam][trueNVtx];
            double W_PUDown = _weight*intlumi*xsection/nbevents*reweightSF_Down[iSam][trueNVtx];
            
            
            /////////////////////////////////////////////////////////
            //Apply trigger selection////////////////////////////////
            /////////////////////////////////////////////////////////
            //It seems FASTSIM samples are not triggered
            //if(!TriggerSelectionMC(samplefilesMC[iSam])) continue;
            
            
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
                    LeptonSF=CalcLepSF_ele_FullSim(el_SF_FullSim,el_RecoSF_FullSim,chosenlept[i],0)*CalcLepSF_ele_FastSim(el_SF_FastSim,chosenlept[i],0);
                    LeptonSF_LepSFUp=CalcLepSF_ele_FullSim(el_SF_FullSim,el_RecoSF_FullSim,chosenlept[i],1)*CalcLepSF_ele_FastSim(el_SF_FastSim,chosenlept[i],1);
                    LeptonSF_LepSFDown=CalcLepSF_ele_FullSim(el_SF_FullSim,el_RecoSF_FullSim,chosenlept[i],-1)*CalcLepSF_ele_FastSim(el_SF_FastSim,chosenlept[i],-1);
                }
                //Mu
                else if((*_flavors)[chosenlept[i]]==1) {
                    if ((*_lPt)[chosenlept[i]]<10) {//Low Pt Muon
                        LeptonSF=CalcLepSF_mu_FullSim(mu_dxydz_SF_FullSim,mu_MediumID_SF_FullSim,mu_miniIso_SF_FullSim,mu_LepMVAVM_SF_FullSim,mu_RecoSF_FullSim_Lowpt,chosenlept[i],0)*CalcLepSF_mu_FastSim(mu_dxydz_SF_FastSim,mu_MediumID_SF_FastSim,mu_miniIso_SF_FastSim,mu_LepMVAVM_SF_FastSim,chosenlept[i],0);
                        LeptonSF_LepSFUp=CalcLepSF_mu_FullSim(mu_dxydz_SF_FullSim,mu_MediumID_SF_FullSim,mu_miniIso_SF_FullSim,mu_LepMVAVM_SF_FullSim,mu_RecoSF_FullSim_Lowpt,chosenlept[i],1)*CalcLepSF_mu_FastSim(mu_dxydz_SF_FastSim,mu_MediumID_SF_FastSim,mu_miniIso_SF_FastSim,mu_LepMVAVM_SF_FastSim,chosenlept[i],1);
                        LeptonSF_LepSFDown=CalcLepSF_mu_FullSim(mu_dxydz_SF_FullSim,mu_MediumID_SF_FullSim,mu_miniIso_SF_FullSim,mu_LepMVAVM_SF_FullSim,mu_RecoSF_FullSim_Lowpt,chosenlept[i],-1)*CalcLepSF_mu_FastSim(mu_dxydz_SF_FastSim,mu_MediumID_SF_FastSim,mu_miniIso_SF_FastSim,mu_LepMVAVM_SF_FastSim,chosenlept[i],-1);
                    }
                    else {//High Pt Muon
                        LeptonSF=CalcLepSF_mu_FullSim(mu_dxydz_SF_FullSim,mu_MediumID_SF_FullSim,mu_miniIso_SF_FullSim,mu_LepMVAVM_SF_FullSim,mu_RecoSF_FullSim_Highpt,chosenlept[i],0)*CalcLepSF_mu_FastSim(mu_dxydz_SF_FastSim,mu_MediumID_SF_FastSim,mu_miniIso_SF_FastSim,mu_LepMVAVM_SF_FastSim,chosenlept[i],0);
                        LeptonSF_LepSFUp=CalcLepSF_mu_FullSim(mu_dxydz_SF_FullSim,mu_MediumID_SF_FullSim,mu_miniIso_SF_FullSim,mu_LepMVAVM_SF_FullSim,mu_RecoSF_FullSim_Highpt,chosenlept[i],1)*CalcLepSF_mu_FastSim(mu_dxydz_SF_FastSim,mu_MediumID_SF_FastSim,mu_miniIso_SF_FastSim,mu_LepMVAVM_SF_FastSim,chosenlept[i],1);
                        LeptonSF_LepSFDown=CalcLepSF_mu_FullSim(mu_dxydz_SF_FullSim,mu_MediumID_SF_FullSim,mu_miniIso_SF_FullSim,mu_LepMVAVM_SF_FullSim,mu_RecoSF_FullSim_Highpt,chosenlept[i],-1)*CalcLepSF_mu_FastSim(mu_dxydz_SF_FastSim,mu_MediumID_SF_FastSim,mu_miniIso_SF_FastSim,mu_LepMVAVM_SF_FastSim,chosenlept[i],-1);
                    }
                }
                //Tau
                else if((*_flavors)[chosenlept[i]]==2) {
                    LeptonSF=CalcLepSF_tau_FullSim(0)*CalcLepSF_tau_FastSim(tau_SF_FastSim,chosenlept[i],0);
                    LeptonSF_LepSFUp=CalcLepSF_tau_FullSim(1)*CalcLepSF_tau_FastSim(tau_SF_FastSim,chosenlept[i],1);
                    LeptonSF_LepSFDown=CalcLepSF_tau_FullSim(-1)*CalcLepSF_tau_FastSim(tau_SF_FastSim,chosenlept[i],-1);
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
                
                
                ////////////////////////////
                //Trigger Efficiencies//////
                //FastSim Samples Only//////
                ////////////////////////////
                //3l events
                int nullpar=99;
                W=W*CalcTriggerEff(nTau, nullpar, subleadinglep, trailinglep);
                W_PUUp=W_PUUp*CalcTriggerEff(nTau, nullpar, subleadinglep, trailinglep);
                W_PUDown=W_PUDown*CalcTriggerEff(nTau, nullpar, subleadinglep, trailinglep);
                
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
                        MCYield_Signal_JECUp[0][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(SRIndex,W);
                    }
                    if (CatB) {
                        SRIndex=CalcSRIndex_B(_met_JECup,Mll,Mt_JECUp);
                        MCYield_Signal_JECUp[1][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(SRIndex,W);
                    }
                }
                //JECDown
                if (!bJetandMETflag_JECDown) {
                    if (CatA) {
                        SRIndex=CalcSRIndex_A(_met_JECdown,Mll,Mt_JECDown);
                        MCYield_Signal_JECDown[0][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(SRIndex,W);
                    }
                    if (CatB) {
                        SRIndex=CalcSRIndex_B(_met_JECdown,Mll,Mt_JECDown);
                        MCYield_Signal_JECDown[1][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(SRIndex,W);
                    }
                }
                //NO JEC
                if (!bJetandMETflag) {
                    ///////////////////////////////////////////////
                    //Fill SR Hist Before Overflow Bin Adjustment//
                    ///////////////////////////////////////////////
                    if (CatA) {
                        SRIndex=CalcSRIndex_A(_met,Mll,Mt);
                        MCYield_Signal[0][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(SRIndex,W);
                        MCYield_Signal_BTagUp[0][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(SRIndex,W_BTagUp);
                        MCYield_Signal_BTagDown[0][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(SRIndex,W_BTagDown);
                        MCYield_Signal_LepSFUp[0][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(SRIndex,W_LepSFUp);
                        MCYield_Signal_LepSFDown[0][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(SRIndex,W_LepSFDown);
                        MCYield_Signal_PUUp[0][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(SRIndex,W_PUUp);
                        MCYield_Signal_PUDown[0][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(SRIndex,W_PUDown);
                    }
                    if (CatB) {
                        SRIndex=CalcSRIndex_B(_met,Mll,Mt);
                        MCYield_Signal[1][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(SRIndex,W);
                        MCYield_Signal_BTagUp[1][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(SRIndex,W_BTagUp);
                        MCYield_Signal_BTagDown[1][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(SRIndex,W_BTagDown);
                        MCYield_Signal_LepSFUp[1][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(SRIndex,W_LepSFUp);
                        MCYield_Signal_LepSFDown[1][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(SRIndex,W_LepSFDown);
                        MCYield_Signal_PUUp[1][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(SRIndex,W_PUUp);
                        MCYield_Signal_PUDown[1][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(SRIndex,W_PUDown);
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
                        MCYield_Signal_Pt1[0][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(leadingPt_FL,W);
                        MCYield_Signal_Pt2[0][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(subleadingPt_FL,W);
                        MCYield_Signal_Pt3[0][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(trailingPt_FL,W);
                        MCYield_Signal_Mll[0][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(Mll_FL,W);
                        MCYield_Signal_MT[0][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(Mt_FL,W);
                        MCYield_Signal_MET[0][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(_met_FL,W);
                        MCYield_Signal_M3l[0][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(M3l_FL,W);
                        MCYield_Signal_SumQ[0][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(SumQ,W);
                    }
                    if (CatB) {
                        MCYield_Signal_Pt1[1][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(leadingPt_FL,W);
                        MCYield_Signal_Pt2[1][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(subleadingPt_FL,W);
                        MCYield_Signal_Pt3[1][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(trailingPt_FL,W);
                        MCYield_Signal_Mll[1][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(Mll_FL,W);
                        MCYield_Signal_MT[1][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(Mt_FL,W);
                        MCYield_Signal_MET[1][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(_met_FL,W);
                        MCYield_Signal_M3l[1][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(M3l_FL,W);
                        MCYield_Signal_SumQ[1][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(SumQ,W);
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
                
                
                ////////////////////////////
                //Trigger Efficiencies//////
                //FastSim Samples Only//////
                ////////////////////////////
                //2l+tau events
                int nullpar=99;
                if((*_flavors)[leadinglep]==2){
                    W=W*CalcTriggerEff(nTau, nullpar, subleadinglep, trailinglep);
                    W_PUUp=W_PUUp*CalcTriggerEff(nTau, nullpar, subleadinglep, trailinglep);
                    W_PUDown=W_PUDown*CalcTriggerEff(nTau, nullpar, subleadinglep, trailinglep);
                }
                if((*_flavors)[subleadinglep]==2){
                    W=W*CalcTriggerEff(nTau, nullpar, leadinglep, trailinglep);
                    W_PUUp=W_PUUp*CalcTriggerEff(nTau, nullpar, leadinglep, trailinglep);
                    W_PUDown=W_PUDown*CalcTriggerEff(nTau, nullpar, leadinglep, trailinglep);
                }
                if((*_flavors)[trailinglep]==2){
                    W=W*CalcTriggerEff(nTau, nullpar, leadinglep, subleadinglep);
                    W_PUUp=W_PUUp*CalcTriggerEff(nTau, nullpar, leadinglep, subleadinglep);
                    W_PUDown=W_PUDown*CalcTriggerEff(nTau, nullpar, leadinglep, subleadinglep);
                }
                
                
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
                        MCYield_Signal_JECUp[2][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(SRIndex,W);
                    }
                    if (CatD) {
                        SRIndex=CalcSRIndex_D(_met_JECup,Mll,Mt2_JECUp);
                        MCYield_Signal_JECUp[3][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(SRIndex,W);
                    }
                    if (CatE) {
                        SRIndex=CalcSRIndex_E(_met_JECup,Mll,Mt2_JECUp);
                        MCYield_Signal_JECUp[4][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(SRIndex,W);
                    }
                }
                //JECDown
                if (!bJetandMETflag_JECDown) {
                    if (CatC) {
                        SRIndex=CalcSRIndex_C(_met_JECdown,Mll,Mt2_JECDown);
                        MCYield_Signal_JECDown[2][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(SRIndex,W);
                    }
                    if (CatD) {
                        SRIndex=CalcSRIndex_D(_met_JECdown,Mll,Mt2_JECDown);
                        MCYield_Signal_JECDown[3][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(SRIndex,W);
                    }
                    if (CatE) {
                        SRIndex=CalcSRIndex_E(_met_JECdown,Mll,Mt2_JECDown);
                        MCYield_Signal_JECDown[4][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(SRIndex,W);
                    }
                }
                //NO JEC
                if (!bJetandMETflag) {
                    ///////////////////////////////////////////////
                    //Fill SR Hist Before Overflow Bin Adjustment//
                    ///////////////////////////////////////////////
                    if (CatC) {
                        SRIndex=CalcSRIndex_C(_met,Mll,Mt2);
                        MCYield_Signal[2][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(SRIndex,W);
                        MCYield_Signal_BTagUp[2][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(SRIndex,W_BTagUp);
                        MCYield_Signal_BTagDown[2][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(SRIndex,W_BTagDown);
                        MCYield_Signal_LepSFUp[2][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(SRIndex,W_LepSFUp);
                        MCYield_Signal_LepSFDown[2][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(SRIndex,W_LepSFDown);
                        MCYield_Signal_PUUp[2][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(SRIndex,W_PUUp);
                        MCYield_Signal_PUDown[2][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(SRIndex,W_PUDown);
                    }
                    if (CatD) {
                        SRIndex=CalcSRIndex_D(_met,Mll,Mt2);
                        MCYield_Signal[3][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(SRIndex,W);
                        MCYield_Signal_BTagUp[3][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(SRIndex,W_BTagUp);
                        MCYield_Signal_BTagDown[3][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(SRIndex,W_BTagDown);
                        MCYield_Signal_LepSFUp[3][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(SRIndex,W_LepSFUp);
                        MCYield_Signal_LepSFDown[3][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(SRIndex,W_LepSFDown);
                        MCYield_Signal_PUUp[3][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(SRIndex,W_PUUp);
                        MCYield_Signal_PUDown[3][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(SRIndex,W_PUDown);
                    }
                    if (CatE) {
                        SRIndex=CalcSRIndex_E(_met,Mll,Mt2);
                        MCYield_Signal[4][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(SRIndex,W);
                        MCYield_Signal_BTagUp[4][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(SRIndex,W_BTagUp);
                        MCYield_Signal_BTagDown[4][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(SRIndex,W_BTagDown);
                        MCYield_Signal_LepSFUp[4][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(SRIndex,W_LepSFUp);
                        MCYield_Signal_LepSFDown[4][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(SRIndex,W_LepSFDown);
                        MCYield_Signal_PUUp[4][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(SRIndex,W_PUUp);
                        MCYield_Signal_PUDown[4][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(SRIndex,W_PUDown);
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
                        MCYield_Signal_Pt1[2][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(leadingPt_FL,W);
                        MCYield_Signal_Pt2[2][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(subleadingPt_FL,W);
                        MCYield_Signal_Pt3[2][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(trailingPt_FL,W);
                        MCYield_Signal_Mll[2][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(Mll_FL,W);
                        MCYield_Signal_MT[2][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(Mt2_FL,W);
                        MCYield_Signal_MET[2][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(_met_FL,W);
                        MCYield_Signal_M3l[2][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(M3l_FL,W);
                        MCYield_Signal_SumQ[2][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(SumQ,W);
                    }
                    if (CatD) {
                        MCYield_Signal_Pt1[3][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(leadingPt_FL,W);
                        MCYield_Signal_Pt2[3][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(subleadingPt_FL,W);
                        MCYield_Signal_Pt3[3][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(trailingPt_FL,W);
                        MCYield_Signal_Mll[3][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(Mll_FL,W);
                        MCYield_Signal_MT[3][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(Mt2_FL,W);
                        MCYield_Signal_MET[3][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(_met_FL,W);
                        MCYield_Signal_M3l[3][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(M3l_FL,W);
                        MCYield_Signal_SumQ[3][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(SumQ,W);
                    }
                    if (CatE) {
                        MCYield_Signal_Pt1[4][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(leadingPt_FL,W);
                        MCYield_Signal_Pt2[4][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(subleadingPt_FL,W);
                        MCYield_Signal_Pt3[4][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(trailingPt_FL,W);
                        MCYield_Signal_Mll[4][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(Mll_FL,W);
                        MCYield_Signal_MT[4][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(Mt2_FL,W);
                        MCYield_Signal_MET[4][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(_met_FL,W);
                        MCYield_Signal_M3l[4][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(M3l_FL,W);
                        MCYield_Signal_SumQ[4][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(SumQ,W);
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
            
                ////////////////////////////
                //Trigger Efficiencies//////
                //FastSim Samples Only//////
                ////////////////////////////
                //l+2tau events
                int nullpar=99;
                for (int i=0; i<3; i++) {
                    if ((*_flavors)[chosenlept[i]]==2) continue;
                    W=W*CalcTriggerEff(nTau, chosenlept[i], nullpar, nullpar);
                    W_PUUp=W_PUUp*CalcTriggerEff(nTau, chosenlept[i], nullpar, nullpar);
                    W_PUDown=W_PUDown*CalcTriggerEff(nTau, chosenlept[i], nullpar, nullpar);
                }
                
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
                        MCYield_Signal_JECUp[5][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(SRIndex,W);
                    }
                }
                //JECDown
                if (!bJetandMETflag_JECDown) {
                    if (CatF) {
                        SRIndex=CalcSRIndex_F(_met_JECdown,Mll,Mt2_JECDown);
                        MCYield_Signal_JECDown[5][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(SRIndex,W);
                    }
                }
                //NO JEC
                if (!bJetandMETflag) {
                    ///////////////////////////////////////////////
                    //Fill SR Hist Before Overflow Bin Adjustment//
                    ///////////////////////////////////////////////
                    if (CatF) {
                        SRIndex=CalcSRIndex_F(_met,Mll,Mt2);
                        MCYield_Signal[5][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(SRIndex,W);
                        MCYield_Signal_BTagUp[5][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(SRIndex,W_BTagUp);
                        MCYield_Signal_BTagDown[5][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(SRIndex,W_BTagDown);
                        MCYield_Signal_LepSFUp[5][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(SRIndex,W_LepSFUp);
                        MCYield_Signal_LepSFDown[5][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(SRIndex,W_LepSFDown);
                        MCYield_Signal_PUUp[5][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(SRIndex,W_PUUp);
                        MCYield_Signal_PUDown[5][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(SRIndex,W_PUDown);
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
                        MCYield_Signal_Pt1[5][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(leadingPt_FL,W);
                        MCYield_Signal_Pt2[5][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(subleadingPt_FL,W);
                        MCYield_Signal_Pt3[5][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(trailingPt_FL,W);
                        MCYield_Signal_Mll[5][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(Mll_FL,W);
                        MCYield_Signal_MT[5][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(Mt2_FL,W);
                        MCYield_Signal_MET[5][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(_met_FL,W);
                        MCYield_Signal_M3l[5][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(M3l_FL,W);
                        MCYield_Signal_SumQ[5][MassIndexCN*nbLSPMassPoint+MassIndexLSP]->Fill(SumQ,W);
                    }
                }//NO JEC

                
            }//if(nTau==2)
            /*
             *2Tau
             */
            ////////////////////////////////////////////////////////////////////////////////////

            
        }//loop over events/entries
        
        /*
         *Write the histograms
         */
        TFile *outputRootFile=new TFile("FinalYield_Signal_"+SignalProcess[iSam]+".root","CREATE");
        for (int i=0; i<nSRCat; i++) for (int k=0; k<nbCNMassPoint; k++) for (int l=0; l<nbLSPMassPoint; l++){
            
            if (EventPerModel[k][l]==0) continue;
            
            if (MCYield_Signal[i][k*nbLSPMassPoint+l]->GetEntries()==0) {
                MCYield_Signal[i][k*nbLSPMassPoint+l]->Delete();
            }
            else{
                MCYield_Signal[i][k*nbLSPMassPoint+l]->Write();
                MCYield_Signal[i][k*nbLSPMassPoint+l]->Delete();
            }
            
            if (MCYield_Signal_JECUp[i][k*nbLSPMassPoint+l]->GetEntries()==0) {
                MCYield_Signal_JECUp[i][k*nbLSPMassPoint+l]->Delete();
            }
            else{
                MCYield_Signal_JECUp[i][k*nbLSPMassPoint+l]->Write();
                MCYield_Signal_JECUp[i][k*nbLSPMassPoint+l]->Delete();
            }
            
            if (MCYield_Signal_JECDown[i][k*nbLSPMassPoint+l]->GetEntries()==0) {
                MCYield_Signal_JECDown[i][k*nbLSPMassPoint+l]->Delete();
            }
            else{
                MCYield_Signal_JECDown[i][k*nbLSPMassPoint+l]->Write();
                MCYield_Signal_JECDown[i][k*nbLSPMassPoint+l]->Delete();
            }

            if (MCYield_Signal_BTagUp[i][k*nbLSPMassPoint+l]->GetEntries()==0) {
                MCYield_Signal_BTagUp[i][k*nbLSPMassPoint+l]->Delete();
            }
            else{
                MCYield_Signal_BTagUp[i][k*nbLSPMassPoint+l]->Write();
                MCYield_Signal_BTagUp[i][k*nbLSPMassPoint+l]->Delete();
            }
            
            if (MCYield_Signal_BTagDown[i][k*nbLSPMassPoint+l]->GetEntries()==0) {
                MCYield_Signal_BTagDown[i][k*nbLSPMassPoint+l]->Delete();
            }
            else{
                MCYield_Signal_BTagDown[i][k*nbLSPMassPoint+l]->Write();
                MCYield_Signal_BTagDown[i][k*nbLSPMassPoint+l]->Delete();
            }
            
            if (MCYield_Signal_LepSFUp[i][k*nbLSPMassPoint+l]->GetEntries()==0) {
                MCYield_Signal_LepSFUp[i][k*nbLSPMassPoint+l]->Delete();
            }
            else{
                MCYield_Signal_LepSFUp[i][k*nbLSPMassPoint+l]->Write();
                MCYield_Signal_LepSFUp[i][k*nbLSPMassPoint+l]->Delete();
            }
            
            if (MCYield_Signal_LepSFDown[i][k*nbLSPMassPoint+l]->GetEntries()==0) {
                MCYield_Signal_LepSFDown[i][k*nbLSPMassPoint+l]->Delete();
            }
            else{
                MCYield_Signal_LepSFDown[i][k*nbLSPMassPoint+l]->Write();
                MCYield_Signal_LepSFDown[i][k*nbLSPMassPoint+l]->Delete();
            }
            
            if (MCYield_Signal_PUUp[i][k*nbLSPMassPoint+l]->GetEntries()==0) {
                MCYield_Signal_PUUp[i][k*nbLSPMassPoint+l]->Delete();
            }
            else{
                MCYield_Signal_PUUp[i][k*nbLSPMassPoint+l]->Write();
                MCYield_Signal_PUUp[i][k*nbLSPMassPoint+l]->Delete();
            }
            
            if (MCYield_Signal_PUDown[i][k*nbLSPMassPoint+l]->GetEntries()==0) {
                MCYield_Signal_PUDown[i][k*nbLSPMassPoint+l]->Delete();
            }
            else{
                MCYield_Signal_PUDown[i][k*nbLSPMassPoint+l]->Write();
                MCYield_Signal_PUDown[i][k*nbLSPMassPoint+l]->Delete();
            }
            
            if (MCYield_Signal_Pt1[i][k*nbLSPMassPoint+l]->GetEntries()==0) {
                MCYield_Signal_Pt1[i][k*nbLSPMassPoint+l]->Delete();
            }
            else{
                MCYield_Signal_Pt1[i][k*nbLSPMassPoint+l]->Write();
                MCYield_Signal_Pt1[i][k*nbLSPMassPoint+l]->Delete();
            }
            
            if (MCYield_Signal_Pt2[i][k*nbLSPMassPoint+l]->GetEntries()==0) {
                MCYield_Signal_Pt2[i][k*nbLSPMassPoint+l]->Delete();
            }
            else{
                MCYield_Signal_Pt2[i][k*nbLSPMassPoint+l]->Write();
                MCYield_Signal_Pt2[i][k*nbLSPMassPoint+l]->Delete();
            }
            
            if (MCYield_Signal_Pt3[i][k*nbLSPMassPoint+l]->GetEntries()==0) {
                MCYield_Signal_Pt3[i][k*nbLSPMassPoint+l]->Delete();
            }
            else{
                MCYield_Signal_Pt3[i][k*nbLSPMassPoint+l]->Write();
                MCYield_Signal_Pt3[i][k*nbLSPMassPoint+l]->Delete();
            }
            
            if (MCYield_Signal_Mll[i][k*nbLSPMassPoint+l]->GetEntries()==0) {
                MCYield_Signal_Mll[i][k*nbLSPMassPoint+l]->Delete();
            }
            else{
                MCYield_Signal_Mll[i][k*nbLSPMassPoint+l]->Write();
                MCYield_Signal_Mll[i][k*nbLSPMassPoint+l]->Delete();
            }
            
            if (MCYield_Signal_MT[i][k*nbLSPMassPoint+l]->GetEntries()==0) {
                MCYield_Signal_MT[i][k*nbLSPMassPoint+l]->Delete();
            }
            else{
                MCYield_Signal_MT[i][k*nbLSPMassPoint+l]->Write();
                MCYield_Signal_MT[i][k*nbLSPMassPoint+l]->Delete();
            }
            
            if (MCYield_Signal_MET[i][k*nbLSPMassPoint+l]->GetEntries()==0) {
                MCYield_Signal_MET[i][k*nbLSPMassPoint+l]->Delete();
            }
            else{
                MCYield_Signal_MET[i][k*nbLSPMassPoint+l]->Write();
                MCYield_Signal_MET[i][k*nbLSPMassPoint+l]->Delete();
            }
            
            if (MCYield_Signal_M3l[i][k*nbLSPMassPoint+l]->GetEntries()==0) {
                MCYield_Signal_M3l[i][k*nbLSPMassPoint+l]->Delete();
            }
            else{
                MCYield_Signal_M3l[i][k*nbLSPMassPoint+l]->Write();
                MCYield_Signal_M3l[i][k*nbLSPMassPoint+l]->Delete();
            }
            
            if (MCYield_Signal_SumQ[i][k*nbLSPMassPoint+l]->GetEntries()==0) {
                MCYield_Signal_SumQ[i][k*nbLSPMassPoint+l]->Delete();
            }
            else{
                MCYield_Signal_SumQ[i][k*nbLSPMassPoint+l]->Write();
                MCYield_Signal_SumQ[i][k*nbLSPMassPoint+l]->Delete();
            }
            
        }
        outputRootFile->Close();
        /*
         *
         */
        
        cout<<endl;
        inputMC[iSam]->Close();
        
    }//loop over samples
    

}//void FinalYield_Signal::loop


///////////////////////////////////////////////////////////////////////////////////////
//////////////////////////Appendix- Defined Functions//////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
//Trigger Selection
bool FinalYield_Signal::TriggerSelectionMC(TString& MCsample){
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
bool FinalYield_Signal::bJetandMETcut(double& MET, const double& VarCoef){
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
bool FinalYield_Signal::ObjectSelection(int& n){
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
bool FinalYield_Signal::FOSelection(int& n){
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
bool FinalYield_Signal::DeltaRSelection(int& n,int& l){
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
double FinalYield_Signal::CalcInvariantMass(int& n,int& l){
    TLorentzVector lep1;
    TLorentzVector lep2;
    
    Float_t et1 = (*_lPt)[n];
    Float_t et2 = (*_lPt)[l];
    lep1.SetPtEtaPhiE(et1, (*_lEta)[n], (*_lPhi)[n], (et1 * cosh((*_lEta)[n])));
    lep2.SetPtEtaPhiE(et2, (*_lEta)[l], (*_lPhi)[l], (et2 * cosh((*_lEta)[l])));
    
    return (lep1+lep2).Mag();
}

double FinalYield_Signal::CalcTrilepInvariantMass(int& n,int& l, int& m){
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
double FinalYield_Signal::CalcTransMass(int& n, double& MET, double& METPhi){
    double TransMass=sqrt(2*(*_lPt)[n]*MET*(1-cos((*_lPhi)[n]-METPhi)));
    return TransMass;
}

//Define the function calculating MT2
double FinalYield_Signal::CalcMt2(double* pa, double* pb, double* pmiss){
    double Mt2;
    mt2_calculator.set_momenta(pa,pb,pmiss);
    mt2_calculator.set_mn(0);
    Mt2=mt2_calculator.get_mt2();
    return Mt2;
}

//Define the function calculating dR(jet,FO)
double FinalYield_Signal::CalcdR(int& n, int& l){
    double pi=3.1415926;
    double dphi=(abs((*_jetPhi)[n]-(*_lPhi)[l])>pi)?(2*pi-abs((*_jetPhi)[n]-(*_lPhi)[l])):abs((*_jetPhi)[n]-(*_lPhi)[l]);
    double deta=abs((*_jetEta)[n]-(*_lEta)[l]);
    double dr=sqrt(dphi*dphi+deta*deta);
    return dr;
}

//Define the function calculating recoil vector
/*
double FinalYield_Signal::Calcu1(double& uphi, double& bosonphi, double& umag){
    double pi=3.1415926;
    double dphi=(abs(uphi-bosonphi)>pi)?(2*pi-abs(uphi-bosonphi)):abs(uphi-bosonphi);
    double u1=umag*cos(dphi);
    return u1;
}
double FinalYield_Signal::Calcu2(double& uphi, double& bosonphi, double& umag){
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
double FinalYield_Signal::CalcBTagSF(const double& VarCoef){
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
double FinalYield_Signal::CalcLepSF_ele_FullSim(TH2D* SF, TH2D* RecoSF, int& n, const double& VarCoef){
    
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
double FinalYield_Signal::CalcLepSF_ele_FastSim(TH2D* SF, int& n, const double& VarCoef){
    
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
double FinalYield_Signal::CalcLepSF_mu_FullSim(TH2D* SF1, TH2D* SF2, TH2D* SF3, TH2D* SF4, TGraphAsymmErrors* RecoSF, int& n, const double& VarCoef){
    
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
double FinalYield_Signal::CalcLepSF_mu_FastSim(TH2D* SF1, TH2D* SF2, TH2D* SF3, TH2D* SF4, int& n, const double& VarCoef){
    
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
double FinalYield_Signal::CalcLepSF_tau_FullSim(const double& VarCoef){
    
    //0.95 with 5% uncertainty for 2016 Reco MiniAOD data
    double tauSF=0.95+VarCoef*0.05;
    return tauSF;
    
}

//Function Calculating tau LepSF (FullSim-FastSim)
double FinalYield_Signal::CalcLepSF_tau_FastSim(TH2D* SF, int& n, const double& VarCoef){
    
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
int FinalYield_Signal::CalcSRIndex_A(double& MET, double& Mll, double& MT){
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

int FinalYield_Signal::CalcSRIndex_B(double& MET, double& Mll, double& MT){
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

int FinalYield_Signal::CalcSRIndex_C(double& MET, double& Mll, double& MT){
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

int FinalYield_Signal::CalcSRIndex_D(double& MET, double& Mll, double& MT){
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

int FinalYield_Signal::CalcSRIndex_E(double& MET, double& Mll, double& MT){
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

int FinalYield_Signal::CalcSRIndex_F(double& MET, double& Mll, double& MT){
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

//Function calculating X-section of signal samples
double FinalYield_Signal::FindSignalXSection(double& mchargino, double& mneutralino2){
    
    //Set up neutralino2 mass point
    int CNindex=0;
    double DCNMass=9999;
    double CNMassPoint[nbCNMassPoint]={0};
    for (int i=0; i<nbCNMassPoint; i++) {
        CNMassPoint[i]=100+i*25;
        if (abs(mneutralino2-CNMassPoint[i])<DCNMass) {
            CNindex=i;
            DCNMass=abs(mneutralino2-CNMassPoint[i]);
        }
    }
    
    //Set up x-section vector
    double xsection[nbCNMassPoint]={
        //fb
        22504,
        9936.7,
        5118.9,
        2912.9,
        1779.1,
        1145.7,
        768.61,
        532.81,
        379.23,
        276.17,
        204.91,
        154.54,
        118.22,
        91.52,
        71.65,
        56.6,
        45.12,
        36.27,
        29.32,
        23.87,
        19.55,
        16.1,
        13.31,
        11.06,
        9.21,
        7.71,
        6.47,
        5.46,
        4.61,
        3.91,
        3.32,
        2.82,
        2.41,
        2.06,
        1.76,
        1.52,
        1.3,
        1.12,
        0.96,
        0.84,
        0.72,
        0.63,
        0.54,
        0.47,
        0.41,
       	0.35,
        0.31,
        0.27,
        0.24,
       	0.21,
       	0.18,
        0.16,
        0.13,
        0.12,
        0.11,
        0.09,
        0.078948,
        0.069138,
        0.06046,
        0.052947,
        0.046398,
        0.040543,
        0.035545,
        0.031169,
        0.027339,
        0.023989,
        0.021048,
        0.018475,
        0.016221,
        0.014244,
        0.012506,
        0.010983,
        0.00968,
        0.0085234,
        0.0074991,
        0.0066073,
        0.0058009
    };
    
    //Find the x-section
    return xsection[CNindex]*0.001;
    
}

int FinalYield_Signal::FindCNIndex(double& mneutralino2){
    int CNindex=0;
    double DCNMass=9999;
    double CNMassPoint[nbCNMassPoint]={0};
    for (int i=0; i<nbCNMassPoint; i++) {
        CNMassPoint[i]=100+i*25;
        if (abs(_mneutralino2-CNMassPoint[i])<DCNMass) {
            CNindex=i;
            DCNMass=abs(_mneutralino2-CNMassPoint[i]);
        }
    }
    return CNindex;
}

int FinalYield_Signal::FindLSPIndex(double& mLSP){
    int LSPindex=0;
    double DLSPMass=9999;
    double LSPMassPoint[nbLSPMassPoint]={0};
    for (int i=0; i<nbLSPMassPoint; i++) {
        LSPMassPoint[i]=i*25;
        if (abs(_mLSP-LSPMassPoint[i])<DLSPMass) {
            LSPindex=i;
            DLSPMass=abs(_mLSP-LSPMassPoint[i]);
        }
    }
    return LSPindex;
}

double FinalYield_Signal::CalcTriggerEff(int& MapType, int& l, int& l1, int& l2){
    
    const int nbPtbin_2tau=12;
    const int nbEtabinele_2tau=4;
    const int nbEtabinmuon_2tau=5;
    const int nbPtbin1_1tau=3;
    const int nbPtbin2_1tau=4;
    const int nbPtbin1_0tau=8;
    const int nbPtbin2_0tau=8;
    
    double Etabinele_2tau[nbEtabinele_2tau+1]={0.0,0.8,1.479,2.0,2.5};//electron
    double Etabinmuon_2tau[nbEtabinmuon_2tau+1]={0.0,0.8,1.25,1.6,2.1,2.4};//muon
    double Ptbin_2tau[nbPtbin_2tau+1]={10,15,20,25,30,35,40,50,75,100,200,300,450};//electron/muon
    //double Ptbin1_1tau[nbPtbin1_1tau+1]={20,30,40,50};
    double Ptbin1_1tau[nbPtbin1_1tau+1]={10,30,40,50};
    double Ptbin2_1tau[nbPtbin2_1tau+1]={10,20,30,40,50};
    double Ptbin1_0tau[nbPtbin1_0tau+1]={10,15,20,25,30,35,40,45,50};
    double Ptbin2_0tau[nbPtbin2_0tau+1]={10,15,20,25,30,35,40,45,50};
    
    double Eff_2tau_SE[nbEtabinele_2tau][nbPtbin_2tau]={
        0.006,0.002,0.005,0.376,0.702,0.774,0.841,0.879,0.914,0.932,0.940,0.945,
        0.004,0.003,0.006,0.355,0.710,0.789,0.855,0.890,0.919,0.931,0.944,0.950,
        0.002,0.003,0.010,0.314,0.602,0.688,0.757,0.787,0.818,0.833,0.852,0.875,
        0.001,0.002,0.005,0.230,0.565,0.644,0.705,0.745,0.810,0.834,0.878,0.886,
    };
    double Eff_2tau_SE_Error[nbEtabinele_2tau][nbPtbin_2tau]={
        0.000,0.000,0.000,0.001,0.001,0.000,0.000,0.000,0.001,0.001,0.003,0.006,
        0.000,0.000,0.000,0.001,0.001,0.001,0.000,0.001,0.001,0.001,0.004,0.008,
        0.000,0.000,0.000,0.002,0.001,0.001,0.001,0.001,0.003,0.003,0.010,0.020,
        0.000,0.000,0.000,0.002,0.002,0.001,0.001,0.002,0.004,0.004,0.015,0.036,
    };
    double Eff_2tau_SM[nbEtabinmuon_2tau][nbPtbin_2tau]={
        0.004,0.014,0.226,0.865,0.891,0.905,0.918,0.924,0.930,0.919,0.909,0.901,
        0.003,0.011,0.211,0.852,0.883,0.899,0.912,0.917,0.920,0.904,0.879,0.851,
        0.0001,0.001,0.210,0.844,0.876,0.895,0.911,0.917,0.921,0.915,0.899,0.893,
        0.0001,0.001,0.201,0.793,0.818,0.830,0.841,0.849,0.861,0.852,0.826,0.835,
        0.0001,0.002,0.194,0.735,0.781,0.803,0.821,0.833,0.850,0.837,0.796,0.822,
    };
    double Eff_2tau_SM_Error[nbEtabinmuon_2tau][nbPtbin_2tau]={
        0.000,0.000,0.000,0.000,0.009,0.000,0.000,0.000,0.000,0.001,0.003,0.006,
        0.000,0.000,0.001,0.000,0.000,0.000,0.000,0.000,0.001,0.001,0.005,0.012,
        0.000,0.000,0.001,0.001,0.000,0.000,0.000,0.000,0.001,0.002,0.008,0.018,
        0.000,0.000,0.001,0.001,0.000,0.000,0.000,0.000,0.001,0.002,0.010,0.026,
        0.000,0.000,0.001,0.001,0.001,0.001,0.000,0.001,0.002,0.004,0.020,0.045,
    };
    double Eff_1tau[nbPtbin2_1tau][nbPtbin1_1tau]={
        0.916,0.846,0.926,
        0.857,0.960,0.955,
        -1.00,0.989,0.982,
        -1.00,-1.00,0.982,
    };
    double Eff_1tau_Error[nbPtbin2_1tau][nbPtbin1_1tau]={
        0.022,0.027,0.025,
        0.036,0.016,0.018,
        -1.00,0.016,0.013,
        -1.00,-1.00,0.015,
    };
    double Eff_0tau[nbPtbin2_0tau][nbPtbin1_0tau]={
        0.921,0.950,0.963,0.975,0.992,0.982,0.984,0.993,
        -1.00,0.943,0.961,0.983,0.992,0.988,0.994,0.997,
        -1.00,-1.00,0.980,0.988,0.985,0.988,0.993,0.989,
        -1.00,-1.00,-1.00,0.991,0.994,0.994,0.994,0.993,
        -1.00,-1.00,-1.00,-1.00,0.995,0.997,0.998,0.996,
        -1.00,-1.00,-1.00,-1.00,-1.00,1.000,0.999,0.994,
        -1.00,-1.00,-1.00,-1.00,-1.00,-1.00,1.000,1.000,
        -1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,1.000,
    };
    double Eff_0tau_Error[nbPtbin2_0tau][nbPtbin1_0tau]={
        0.012,0.008,0.007,0.007,0.005,0.007,0.008,0.008,
        -1.00,0.011,0.007,0.005,0.005,0.006,0.005,0.006,
        -1.00,-1.00,0.008,0.005,0.005,0.005,0.005,0.007,
        -1.00,-1.00,-1.00,0.006,0.004,0.004,0.005,0.005,
        -1.00,-1.00,-1.00,-1.00,0.005,0.004,0.004,0.005,
        -1.00,-1.00,-1.00,-1.00,-1.00,0.004,0.004,0.005,
        -1.00,-1.00,-1.00,-1.00,-1.00,-1.00,0.005,0.004,
        -1.00,-1.00,-1.00,-1.00,-1.00,-1.00,-1.00,0.007,
    };
    
    double TriggerEff=1;
    int MapBinx=9999;
    int MapBiny=9999;
    
    if (MapType==2) {
        //2tau+l events
        for (int i=0; i<nbPtbin_2tau; i++) {
            if ((i!=nbPtbin_2tau-1&&Ptbin_2tau[i]<(*_lPt)[l]&&Ptbin_2tau[i+1]>(*_lPt)[l])||(i==nbPtbin_2tau-1&&Ptbin_2tau[i]<(*_lPt)[l])) MapBinx=i;
        }
        if ((*_flavors)[l]==0) {//2tau+e
            for (int i=0; i<nbEtabinele_2tau; i++) {
                if ((i!=nbEtabinele_2tau-1&&Etabinele_2tau[i]<abs((*_lEta)[l])&&Etabinele_2tau[i+1]>abs((*_lEta)[l]))||(i==nbEtabinele_2tau-1&&Etabinele_2tau[i]<abs((*_lEta)[l]))) MapBiny=i;
            }
            TriggerEff=Eff_2tau_SE[MapBiny][MapBinx];
        }
        if ((*_flavors)[l]==1) {//2tau+mu
            for (int i=0; i<nbEtabinmuon_2tau; i++) {
                if ((i!=nbEtabinmuon_2tau-1&&Etabinmuon_2tau[i]<abs((*_lEta)[l])&&Etabinmuon_2tau[i+1]>abs((*_lEta)[l]))||(i==nbEtabinmuon_2tau-1&&Etabinmuon_2tau[i]<abs((*_lEta)[l]))) MapBiny=i;
            }
            TriggerEff=Eff_2tau_SM[MapBiny][MapBinx];
        }
    }
    
    if (MapType==1) {
        //1tau+2l events
        for (int i=0; i<nbPtbin1_1tau; i++) {
            if ((i!=nbPtbin1_1tau-1&&Ptbin1_1tau[i]<(*_lPt)[l1]&&Ptbin1_1tau[i+1]>(*_lPt)[l1])||(i==nbPtbin1_1tau-1&&Ptbin1_1tau[i]<(*_lPt)[l1])) MapBinx=i;
        }
        for (int i=0; i<nbPtbin2_1tau; i++) {
            if ((i!=nbPtbin2_1tau-1&&Ptbin2_1tau[i]<(*_lPt)[l2]&&Ptbin2_1tau[i+1]>(*_lPt)[l2])||(i==nbPtbin2_1tau-1&&Ptbin2_1tau[i]<(*_lPt)[l2])) MapBiny=i;
        }
        TriggerEff=Eff_1tau[MapBiny][MapBinx];
        
    }
    
    if (MapType==0) {
        //0tau+3l events
        for (int i=0; i<nbPtbin1_0tau; i++) {
            if ((i!=nbPtbin1_0tau-1&&Ptbin1_0tau[i]<(*_lPt)[l1]&&Ptbin1_0tau[i+1]>(*_lPt)[l1])||(i==nbPtbin1_0tau-1&&Ptbin1_0tau[i]<(*_lPt)[l1])) MapBinx=i;
        }
        for (int i=0; i<nbPtbin2_0tau; i++) {
            if ((i!=nbPtbin2_0tau-1&&Ptbin2_0tau[i]<(*_lPt)[l2]&&Ptbin2_0tau[i+1]>(*_lPt)[l2])||(i==nbPtbin2_0tau-1&&Ptbin2_0tau[i]<(*_lPt)[l2])) MapBiny=i;
        }
        TriggerEff=Eff_0tau[MapBiny][MapBinx];
    }
    
    if (TriggerEff>1||TriggerEff<0) TriggerEff=1;
    return TriggerEff;
    
}