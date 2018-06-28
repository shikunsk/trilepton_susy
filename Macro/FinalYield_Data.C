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

//Data Trees
/*
 Signal regions
 */
const int nbDatasets = 46; //Single EG/Muon Double EG/Muon MuonEG
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


TString datafiles[nbDatasets] = {
    
    //SingleElectron
    "SingleElectron-Run2016B-03Feb2017-v1-skim3Lreco.root",
    "SingleElectron-Run2016B-03Feb2017-v2-skim3Lreco.root",
    "SingleElectron-Run2016C-03Feb2017-v1-skim3Lreco.root",
    "SingleElectron-Run2016D-03Feb2017-v1-skim3Lreco.root",
    "SingleElectron-Run2016E-03Feb2017-v1-skim3Lreco.root",
    "SingleElectron-Run2016F-03Feb2017-v1-skim3Lreco.root",
    "SingleElectron-Run2016G-03Feb2017-v1-skim3Lreco.root",
    "SingleElectron-Run2016H-03Feb2017-v2-skim3Lreco.root",
    "SingleElectron-Run2016H-03Feb2017-v3-skim3Lreco.root",
    
    //SingleMuon
    "SingleMuon-Run2016B-03Feb2017-v1-skim3Lreco.root",
    "SingleMuon-Run2016B-03Feb2017-v2-skim3Lreco.root",
    "SingleMuon-Run2016C-03Feb2017-v1-skim3Lreco.root",
    "SingleMuon-Run2016D-03Feb2017-v1-skim3Lreco.root",
    "SingleMuon-Run2016E-03Feb2017-v1-skim3Lreco.root",
    "SingleMuon-Run2016F-03Feb2017-v1-skim3Lreco.root",
    "SingleMuon-Run2016G-03Feb2017-v1-skim3Lreco.root",
    "SingleMuon-Run2016H-03Feb2017-v2-skim3Lreco.root",
    "SingleMuon-Run2016H-03Feb2017-v3-skim3Lreco.root",
    
    //DoubleEG
    "DoubleEG-Run2016B-03Feb2017-v1-skim3Lreco.root",
    "DoubleEG-Run2016B-03Feb2017-v2-skim3Lreco.root",
    "DoubleEG-Run2016C-03Feb2017-v1-skim3Lreco.root",
    "DoubleEG-Run2016D-03Feb2017-v1-skim3Lreco.root",
    "DoubleEG-Run2016E-03Feb2017-v1-skim3Lreco.root",
    "DoubleEG-Run2016F-03Feb2017-v1-skim3Lreco.root",
    "DoubleEG-Run2016G-03Feb2017-v1-skim3Lreco.root",
    "DoubleEG-Run2016H-03Feb2017-v2-skim3Lreco.root",
    "DoubleEG-Run2016H-03Feb2017-v3-skim3Lreco.root",
    
    //DoubleMuon
    "DoubleMuon-Run2016B-03Feb2017-v1-skim3Lreco.root",
    "DoubleMuon-Run2016B-03Feb2017-v2-skim3Lreco.root",
    "DoubleMuon-Run2016C-03Feb2017-v1-skim3Lreco.root",
    "DoubleMuon-Run2016C-03Feb2017-v2-skim3Lreco.root",
    "DoubleMuon-Run2016D-03Feb2017-v1-skim3Lreco.root",
    "DoubleMuon-Run2016E-03Feb2017-v1-skim3Lreco.root",
    "DoubleMuon-Run2016F-03Feb2017-v1-skim3Lreco.root",
    "DoubleMuon-Run2016G-03Feb2017-v1-skim3Lreco.root",
    "DoubleMuon-Run2016H-03Feb2017-v2-skim3Lreco.root",
    "DoubleMuon-Run2016H-03Feb2017-v3-skim3Lreco.root",
    
    //MuonEG
    "MuonEG-Run2016B-03Feb2017-v1-skim3Lreco.root",
    "MuonEG-Run2016B-03Feb2017-v2-skim3Lreco.root",
    "MuonEG-Run2016C-03Feb2017-v1-skim3Lreco.root",
    "MuonEG-Run2016D-03Feb2017-v1-skim3Lreco.root",
    "MuonEG-Run2016E-03Feb2017-v1-skim3Lreco.root",
    "MuonEG-Run2016F-03Feb2017-v1-skim3Lreco.root",
    "MuonEG-Run2016G-03Feb2017-v1-skim3Lreco.root",
    "MuonEG-Run2016H-03Feb2017-v2-skim3Lreco.root",
    "MuonEG-Run2016H-03Feb2017-v3-skim3Lreco.root"
    
    
};

void FinalYield::Loop()
{
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
    
    TH1F* DataYield[nSRCat];
    
    TH1F* DataYield_Pt1[nSRCat];
    TH1F* DataYield_Pt2[nSRCat];
    TH1F* DataYield_Pt3[nSRCat];
    TH1F* DataYield_Mll[nSRCat];
    TH1F* DataYield_MT[nSRCat];
    TH1F* DataYield_MET[nSRCat];
    TH1F* DataYield_M3l[nSRCat];
    TH1F* DataYield_SumQ[nSRCat];
    
    for (int i=0; i<nSRCat; i++) {
        name="Data_"+SRCat[i];
        DataYield[i]=new TH1F(name,"",nSR[i],0.5,nSR[i]+0.5);
        DataYield[i]->Sumw2();
        
        name="Data_Pt1_"+SRCat[i];
        DataYield_Pt1[i]=new TH1F(name,"",Ptbin[i],Ptlowerbound,Ptupperbound[i]);
        DataYield_Pt1[i]->Sumw2();
        name="Data_Pt2_"+SRCat[i];
        DataYield_Pt2[i]=new TH1F(name,"",Ptbin[i],Ptlowerbound,Ptupperbound[i]);
        DataYield_Pt2[i]->Sumw2();
        name="Data_Pt3_"+SRCat[i];
        DataYield_Pt3[i]=new TH1F(name,"",Ptbin[i],Ptlowerbound,Ptupperbound[i]);
        DataYield_Pt3[i]->Sumw2();
        name="Data_Mll_"+SRCat[i];
        DataYield_Mll[i]=new TH1F(name,"",Mllbin[i],Mlllowerbound,Mllupperbound[i]);
        DataYield_Mll[i]->Sumw2();
        name="Data_MT_"+SRCat[i];
        DataYield_MT[i]=new TH1F(name,"",MTbin[i],MTlowerbound,MTupperbound[i]);
        DataYield_MT[i]->Sumw2();
        name="Data_MET_"+SRCat[i];
        DataYield_MET[i]=new TH1F(name,"",METbin[i],METlowerbound,METupperbound[i]);
        DataYield_MET[i]->Sumw2();
        name="Data_M3l_"+SRCat[i];
        DataYield_M3l[i]=new TH1F(name,"",M3lbin[i],M3llowerbound,M3lupperbound[i]);
        DataYield_M3l[i]->Sumw2();
        name="Data_SumQ_"+SRCat[i];
        DataYield_SumQ[i]=new TH1F(name,"",SumQbin[i],SumQlowerbound,SumQupperbound[i]);
        DataYield_SumQ[i]->Sumw2();
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
    
    
    ////////////////////////
    //RUN DATASETS//////////
    ////////////////////////
    TFile * inputData[nbDatasets];
    
    for(int iSam = 0;  iSam <nbDatasets; iSam++){
        
        inputData[iSam] = new TFile("/cms/data/store/user/t2/users/takeimai/output/SRYield/2016Data/"+datafiles[iSam],"open");
        TTree *thetree_data = (TTree*)(inputData[iSam])->Get("FakeElectrons/fakeTree");
        Init(thetree_data);
        Long64_t nentries_data = (*thetree_data).GetEntries();//Get the number of generated events
        
        cout<<"Processing Dataset: "<<datafiles[iSam]<<endl;
        
        for (Long64_t jentry=0; jentry<nentries_data;jentry++) {//loop over all the events//First loop
            
            LoadTree(jentry);
            thetree_data->GetEntry(jentry);
            cout<<jentry<<endl;
            
            ////////////////////
            //SELECTION BEGINS//
            ////////////////////
            
            
            /////////////////////////////////////////////////////////
            //Apply trigger selection////////////////////////////////
            /////////////////////////////////////////////////////////
            if(!TriggerSelectionData(datafiles[iSam])) continue;
            
            
            /////////////////////////////////////////////////////////
            //Apply immediate vetos//////////////////////////////////
            /////////////////////////////////////////////////////////
            //Met filter
            if (!passmetfilters) continue;
            //b-Jet veto
            //MET CUT
            if (bJetandMETcut(_met,0)) continue;
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
             Priority of Lepton Pair in Mll Calculation: OSSF>NoOSSF; within NoOSSF Light Lep>Tau
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
                /*
                 *As long as there is an OSSF pair, categorize the event as CatA event
                 */
                bool CatA = false, CatB = false;
                bool OSpairflag = false;
                double InvMassPair=0;
                double DistPeak=9999;
                double Mll=0;
                double Mt=0;
                int SRIndex=9999;
                
                if (((*_charges)[chosenlept[0]] == (*_charges)[chosenlept[1]]) && ((*_charges)[chosenlept[0]] == (*_charges)[chosenlept[2]])){
                    //For 3 SS leptons, Mt takes the smallest value and Mll=0, event goes to CatB
                    Mt=9999;
                    Mll=0;
                    for (int i=0; i<3; i++) {
                        double Mt_temp=CalcTransMass(chosenlept[i],_met,_met_phi);
                        if (Mt_temp<Mt) Mt=Mt_temp;
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
                                        }
                                    }
                                }
                                CatB=true;
                                OSpairflag=true;
                            }
                            if ((*_charges)[chosenlept[0]]!=(*_charges)[chosenlept[1]]) {//aa have opposite sign
                                Mll=CalcInvariantMass(chosenlept[0],chosenlept[1]);
                                Mt=CalcTransMass(chosenlept[2],_met,_met_phi);
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
                                        }
                                    }
                                }
                                CatB=true;
                                OSpairflag=true;
                            }
                            if ((*_charges)[chosenlept[0]]!=(*_charges)[chosenlept[2]]) {//aa have opposite sign
                                Mll=CalcInvariantMass(chosenlept[0],chosenlept[2]);
                                Mt=CalcTransMass(chosenlept[1],_met,_met_phi);
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
                                        }
                                    }
                                }
                                CatB=true;
                                OSpairflag=true;
                            }
                            if ((*_charges)[chosenlept[1]]!=(*_charges)[chosenlept[2]]) {//aa have the opposite sign
                                Mll=CalcInvariantMass(chosenlept[1],chosenlept[2]);
                                Mt=CalcTransMass(chosenlept[0],_met,_met_phi);
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
                
                ///////////////////////////////////////////////
                //Fill SR Hist Before Overflow Bin Adjustment//
                ///////////////////////////////////////////////
                if (CatA) {
                    SRIndex=CalcSRIndex_A(_met,Mll,Mt);
                    DataYield[0]->Fill(SRIndex);
                }
                if (CatB) {
                    SRIndex=CalcSRIndex_B(_met,Mll,Mt);
                    DataYield[1]->Fill(SRIndex);
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
                    DataYield_Pt1[0]->Fill(leadingPt_FL);
                    DataYield_Pt2[0]->Fill(subleadingPt_FL);
                    DataYield_Pt3[0]->Fill(trailingPt_FL);
                    DataYield_Mll[0]->Fill(Mll_FL);
                    DataYield_MT[0]->Fill(Mt_FL);
                    DataYield_MET[0]->Fill(_met_FL);
                    DataYield_M3l[0]->Fill(M3l_FL);
                    DataYield_SumQ[0]->Fill(SumQ);
                }
                if (CatB) {
                    DataYield_Pt1[1]->Fill(leadingPt_FL);
                    DataYield_Pt2[1]->Fill(subleadingPt_FL);
                    DataYield_Pt3[1]->Fill(trailingPt_FL);
                    DataYield_Mll[1]->Fill(Mll_FL);
                    DataYield_MT[1]->Fill(Mt_FL);
                    DataYield_MET[1]->Fill(_met_FL);
                    DataYield_M3l[1]->Fill(M3l_FL);
                    DataYield_SumQ[1]->Fill(SumQ);
                }
                
                
                
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
                int SRIndex=9999;
                
                double ParticleMass[3]={0.00051, 0.105, 1.776};//Mass of Ele/Mu/Tau
                TLorentzVector l1;
                TLorentzVector l2;
                TLorentzVector Emiss;
                
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
                            double pa[3]={ParticleMass[(*_flavors)[chosenlept[1]]],l1.Px(),l1.Py()};
                            double pb[3]={ParticleMass[(*_flavors)[chosenlept[0]]],l2.Px(),l2.Py()};
                            double pmiss[3]={0,Emiss.Px(),Emiss.Py()};
                            Mt2=CalcMt2(pa,pb,pmiss);
                        }
                        else {//Find Leading Light Lepton
                            l1.SetPtEtaPhiE((*_lPt)[chosenlept[2]],0, (*_lPhi)[chosenlept[2]],(*_lPt)[chosenlept[2]]);
                            l2.SetPtEtaPhiE((*_lPt)[chosenlept[0]],0, (*_lPhi)[chosenlept[0]],(*_lPt)[chosenlept[0]]);
                            Emiss.SetPtEtaPhiE(_met,0, _met_phi,_met);
                            double pa[3]={ParticleMass[(*_flavors)[chosenlept[2]]],l1.Px(),l1.Py()};
                            double pb[3]={ParticleMass[(*_flavors)[chosenlept[0]]],l2.Px(),l2.Py()};
                            double pmiss[3]={0,Emiss.Px(),Emiss.Py()};
                            Mt2=CalcMt2(pa,pb,pmiss);
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
                        double pa[3]={ParticleMass[(*_flavors)[chosenlept[1]]],l1.Px(),l1.Py()};
                        double pb[3]={ParticleMass[(*_flavors)[chosenlept[2]]],l2.Px(),l2.Py()};
                        double pmiss[3]={0,Emiss.Px(),Emiss.Py()};
                        Mt2=CalcMt2(pa,pb,pmiss);
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
                            double pa[3]={ParticleMass[(*_flavors)[chosenlept[0]]],l1.Px(),l1.Py()};
                            double pb[3]={ParticleMass[(*_flavors)[chosenlept[1]]],l2.Px(),l2.Py()};
                            double pmiss[3]={0,Emiss.Px(),Emiss.Py()};
                            Mt2=CalcMt2(pa,pb,pmiss);
                        }
                        else {//Find Leading Light Lepton
                            l1.SetPtEtaPhiE((*_lPt)[chosenlept[2]],0, (*_lPhi)[chosenlept[2]],(*_lPt)[chosenlept[2]]);
                            l2.SetPtEtaPhiE((*_lPt)[chosenlept[1]],0, (*_lPhi)[chosenlept[1]],(*_lPt)[chosenlept[1]]);
                            Emiss.SetPtEtaPhiE(_met,0, _met_phi,_met);
                            double pa[3]={ParticleMass[(*_flavors)[chosenlept[2]]],l1.Px(),l1.Py()};
                            double pb[3]={ParticleMass[(*_flavors)[chosenlept[1]]],l2.Px(),l2.Py()};
                            double pmiss[3]={0,Emiss.Px(),Emiss.Py()};
                            Mt2=CalcMt2(pa,pb,pmiss);
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
                        double pa[3]={ParticleMass[(*_flavors)[chosenlept[0]]],l1.Px(),l1.Py()};
                        double pb[3]={ParticleMass[(*_flavors)[chosenlept[2]]],l2.Px(),l2.Py()};
                        double pmiss[3]={0,Emiss.Px(),Emiss.Py()};
                        Mt2=CalcMt2(pa,pb,pmiss);
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
                            double pa[3]={ParticleMass[(*_flavors)[chosenlept[0]]],l1.Px(),l1.Py()};
                            double pb[3]={ParticleMass[(*_flavors)[chosenlept[2]]],l2.Px(),l2.Py()};
                            double pmiss[3]={0,Emiss.Px(),Emiss.Py()};
                            Mt2=CalcMt2(pa,pb,pmiss);
                        }
                        else {//Find Leading Light Lepton
                            l1.SetPtEtaPhiE((*_lPt)[chosenlept[1]],0, (*_lPhi)[chosenlept[1]],(*_lPt)[chosenlept[1]]);
                            l2.SetPtEtaPhiE((*_lPt)[chosenlept[2]],0, (*_lPhi)[chosenlept[2]],(*_lPt)[chosenlept[2]]);
                            Emiss.SetPtEtaPhiE(_met,0, _met_phi,_met);
                            double pa[3]={ParticleMass[(*_flavors)[chosenlept[1]]],l1.Px(),l1.Py()};
                            double pb[3]={ParticleMass[(*_flavors)[chosenlept[2]]],l2.Px(),l2.Py()};
                            double pmiss[3]={0,Emiss.Px(),Emiss.Py()};
                            Mt2=CalcMt2(pa,pb,pmiss);
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
                        double pa[3]={ParticleMass[(*_flavors)[chosenlept[0]]],l1.Px(),l1.Py()};
                        double pb[3]={ParticleMass[(*_flavors)[chosenlept[1]]],l2.Px(),l2.Py()};
                        double pmiss[3]={0,Emiss.Px(),Emiss.Py()};
                        Mt2=CalcMt2(pa,pb,pmiss);
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
                
                ///////////////////////////////////////////////
                //Fill SR Hist Before Overflow Bin Adjustment//
                ///////////////////////////////////////////////
                if (CatC) {
                    SRIndex=CalcSRIndex_C(_met,Mll,Mt2);
                    DataYield[2]->Fill(SRIndex);
                }
                if (CatD) {
                    SRIndex=CalcSRIndex_D(_met,Mll,Mt2);
                    DataYield[3]->Fill(SRIndex);
                }
                if (CatE) {
                    SRIndex=CalcSRIndex_E(_met,Mll,Mt2);
                    DataYield[4]->Fill(SRIndex);
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
                    DataYield_Pt1[2]->Fill(leadingPt_FL);
                    DataYield_Pt2[2]->Fill(subleadingPt_FL);
                    DataYield_Pt3[2]->Fill(trailingPt_FL);
                    DataYield_Mll[2]->Fill(Mll_FL);
                    DataYield_MT[2]->Fill(Mt2_FL);
                    DataYield_MET[2]->Fill(_met_FL);
                    DataYield_M3l[2]->Fill(M3l_FL);
                    DataYield_SumQ[2]->Fill(SumQ);
                }
                if (CatD) {
                    DataYield_Pt1[3]->Fill(leadingPt_FL);
                    DataYield_Pt2[3]->Fill(subleadingPt_FL);
                    DataYield_Pt3[3]->Fill(trailingPt_FL);
                    DataYield_Mll[3]->Fill(Mll_FL);
                    DataYield_MT[3]->Fill(Mt2_FL);
                    DataYield_MET[3]->Fill(_met_FL);
                    DataYield_M3l[3]->Fill(M3l_FL);
                    DataYield_SumQ[3]->Fill(SumQ);
                }
                if (CatE) {
                    DataYield_Pt1[4]->Fill(leadingPt_FL);
                    DataYield_Pt2[4]->Fill(subleadingPt_FL);
                    DataYield_Pt3[4]->Fill(trailingPt_FL);
                    DataYield_Mll[4]->Fill(Mll_FL);
                    DataYield_MT[4]->Fill(Mt2_FL);
                    DataYield_MET[4]->Fill(_met_FL);
                    DataYield_M3l[4]->Fill(M3l_FL);
                    DataYield_SumQ[4]->Fill(SumQ);
                }
                
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
                int SRIndex=9999;
                
                double ParticleMass[3]={0.00051, 0.105, 1.776};//Mass of Ele/Mu/Tau
                TLorentzVector l1;
                TLorentzVector l2;
                TLorentzVector Emiss;
                
                //l+tau+tau////////////
                if((*_flavors)[chosenlept[0]]!=2){
                    CatF=true;
                    //Calculate Mt2
                    if (chosenPt[1]>=chosenPt[2]) {//Find Leading Tau Lepton
                        l1.SetPtEtaPhiE((*_lPt)[chosenlept[1]],0, (*_lPhi)[chosenlept[1]],(*_lPt)[chosenlept[1]]);
                        l2.SetPtEtaPhiE((*_lPt)[chosenlept[0]],0, (*_lPhi)[chosenlept[0]],(*_lPt)[chosenlept[0]]);
                        Emiss.SetPtEtaPhiE(_met,0, _met_phi,_met);
                        double pa[3]={ParticleMass[(*_flavors)[chosenlept[1]]],l1.Px(),l1.Py()};
                        double pb[3]={ParticleMass[(*_flavors)[chosenlept[0]]],l2.Px(),l2.Py()};
                        double pmiss[3]={0,Emiss.Px(),Emiss.Py()};
                        Mt2=CalcMt2(pa,pb,pmiss);
                    }
                    else {//Find Leading Tau Lepton
                        l1.SetPtEtaPhiE((*_lPt)[chosenlept[2]],0, (*_lPhi)[chosenlept[2]],(*_lPt)[chosenlept[2]]);
                        l2.SetPtEtaPhiE((*_lPt)[chosenlept[0]],0, (*_lPhi)[chosenlept[0]],(*_lPt)[chosenlept[0]]);
                        Emiss.SetPtEtaPhiE(_met,0, _met_phi,_met);
                        double pa[3]={ParticleMass[(*_flavors)[chosenlept[2]]],l1.Px(),l1.Py()};
                        double pb[3]={ParticleMass[(*_flavors)[chosenlept[0]]],l2.Px(),l2.Py()};
                        double pmiss[3]={0,Emiss.Px(),Emiss.Py()};
                        Mt2=CalcMt2(pa,pb,pmiss);
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
                        double pa[3]={ParticleMass[(*_flavors)[chosenlept[0]]],l1.Px(),l1.Py()};
                        double pb[3]={ParticleMass[(*_flavors)[chosenlept[1]]],l2.Px(),l2.Py()};
                        double pmiss[3]={0,Emiss.Px(),Emiss.Py()};
                        Mt2=CalcMt2(pa,pb,pmiss);
                    }
                    else {//Find Leading Tau Lepton
                        l1.SetPtEtaPhiE((*_lPt)[chosenlept[2]],0, (*_lPhi)[chosenlept[2]],(*_lPt)[chosenlept[2]]);
                        l2.SetPtEtaPhiE((*_lPt)[chosenlept[1]],0, (*_lPhi)[chosenlept[1]],(*_lPt)[chosenlept[1]]);
                        Emiss.SetPtEtaPhiE(_met,0, _met_phi,_met);
                        double pa[3]={ParticleMass[(*_flavors)[chosenlept[2]]],l1.Px(),l1.Py()};
                        double pb[3]={ParticleMass[(*_flavors)[chosenlept[1]]],l2.Px(),l2.Py()};
                        double pmiss[3]={0,Emiss.Px(),Emiss.Py()};
                        Mt2=CalcMt2(pa,pb,pmiss);
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
                        double pa[3]={ParticleMass[(*_flavors)[chosenlept[0]]],l1.Px(),l1.Py()};
                        double pb[3]={ParticleMass[(*_flavors)[chosenlept[2]]],l2.Px(),l2.Py()};
                        double pmiss[3]={0,Emiss.Px(),Emiss.Py()};
                        Mt2=CalcMt2(pa,pb,pmiss);
                    }
                    else {//Find Leading Tau Lepton
                        l1.SetPtEtaPhiE((*_lPt)[chosenlept[1]],0, (*_lPhi)[chosenlept[1]],(*_lPt)[chosenlept[1]]);
                        l2.SetPtEtaPhiE((*_lPt)[chosenlept[2]],0, (*_lPhi)[chosenlept[2]],(*_lPt)[chosenlept[2]]);
                        Emiss.SetPtEtaPhiE(_met,0, _met_phi,_met);
                        double pa[3]={ParticleMass[(*_flavors)[chosenlept[1]]],l1.Px(),l1.Py()};
                        double pb[3]={ParticleMass[(*_flavors)[chosenlept[2]]],l2.Px(),l2.Py()};
                        double pmiss[3]={0,Emiss.Px(),Emiss.Py()};
                        Mt2=CalcMt2(pa,pb,pmiss);
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
                
                ///////////////////////////////////////////////
                //Fill SR Hist Before Overflow Bin Adjustment//
                ///////////////////////////////////////////////
                if (CatF) {
                    SRIndex=CalcSRIndex_F(_met,Mll,Mt2);
                    DataYield[5]->Fill(SRIndex);
                }
                ///////////////////////////////////////////////
                //Fill SR Hist Before Overflow Bin Adjustment//
                ///////////////////////////////////////////////
                ////////////////////
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
                    DataYield_Pt1[5]->Fill(leadingPt_FL);
                    DataYield_Pt2[5]->Fill(subleadingPt_FL);
                    DataYield_Pt3[5]->Fill(trailingPt_FL);
                    DataYield_Mll[5]->Fill(Mll_FL);
                    DataYield_MT[5]->Fill(Mt2_FL);
                    DataYield_MET[5]->Fill(_met_FL);
                    DataYield_M3l[5]->Fill(M3l_FL);
                    DataYield_SumQ[5]->Fill(SumQ);
                }

                
            }//if(nTau==2)
            /*
             *2Tau
             */
            ////////////////////////////////////////////////////////////////////////////////////

            
        }//loop over events/entries
        
        
        cout<<endl;
        
    }//loop over samples

    
    TFile *outputRootFile=new TFile("FinalYield_Data.root","CREATE");
    
    for (int i=0; i<nSRCat; i++) {
        DataYield[i]->Write();
        DataYield_Pt1[i]->Write();
        DataYield_Pt2[i]->Write();
        DataYield_Pt3[i]->Write();
        DataYield_Mll[i]->Write();
        DataYield_MT[i]->Write();
        DataYield_MET[i]->Write();
        DataYield_M3l[i]->Write();
        DataYield_SumQ[i]->Write();
    }
    
    outputRootFile->Close();
    

}//void FinalYield::loop


///////////////////////////////////////////////////////////////////////////////////////
//////////////////////////Appendix- Defined Functions//////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////

//Trigger Selection
bool FinalYield::TriggerSelectionData(TString& Dataset){
    /*
     *EVENTS MAY BELONG TO DIFFERENT DATASETS SO DROP THE TRIGGERS THAT HAVE BEEN PASSED
     *DO NOT USE TAU TRIGGERS FOR NOW
     
     All Categories:
     
     passHLT_Ele27_WPTight_Gsf --SingleEle
     passHLT_IsoMu24 --SingleMu
     passHLT_IsoTkMu24 --SingleMu
     
     eee,eetau:
     passHLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ --DoubleEG
     
     μμμ,μμtau:
     passHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL --DoubleMu
     passHLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL --DoubleMu
     passHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ (Only for RunH) --DoubleMu
     passHLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ (Only for RunH: RunNo>280919) --DoubleMu
     
     eeμ:
     passHLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ --DoubleEG
     passHLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL --MuEG
     passHLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ (Only for RunH) --MuEG
     passHLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL --MuEG
     passHLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ (Only for RunH) --MuEG
     
     
     eμμ
     passHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL  --DoubleMu
     passHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ (Only for RunH)  --DoubleMu
     passHLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL  --DoubleMu
     passHLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ (Only for RunH)  --DoubleMu
     passHLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL  --MuEG
     passHLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ (Only for RunH)  --MuEG
     passHLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL  --MuEG
     passHLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ (Only for RunH)  --MuEG
     
     eμtau
     passHLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL  --MuEG
     passHLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ (Only for RunH)  --MuEG
     passHLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL  --MuEG
     passHLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ (Only for RunH)  --MuEG
     
     /////////////////////////////
     //DROP OTHER TAU LEPTON TRIGGERS///
     /////////////////////////////
     etautau
     passHLT_Ele22_eta2p1_WPLoose_GSF_LooseIsoPFtau20_SingleL1 (Runs_B-E)
     passHLT_Ele24_eta2p1_WPLoose_GSF_LooseIsoPFtau30 (Runs_F-H)
     passHLT_Ele36_eta2p1_WPLoose_GSF_LooseIsoPFtau20_SingleL1 (Runs_F-H)
     
     μtautau
     passHLT_IsoMu19_eta2p1_LooseIsoPFTau20_SingleL1_v
     /////////////////////////////
     //DROP Tau LEPTON TRIGGERS///
     /////////////////////////////
     */
    
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
    
    
    //NOTICE: DATASETS HAVE TO BE LOOPED IN !ORDER! FOR THE FUNCTION TO BE PROCESSED CORRECTLY
    //SingleElectron
    if (Dataset.Contains("SingleElectron")) {
        if (passHLT_Ele27_WPTight_Gsf) passtriggerselection=true;
    }
    
    //SingleMuon
    if (Dataset.Contains("SingleMuon")) {
        if (passHLT_IsoMu24||passHLT_IsoTkMu24) passtriggerselection=true;
        if (passHLT_Ele27_WPTight_Gsf) passtriggerselection=false;//Drop repeated events
    }
    
    //DoubleElectron
    if (Dataset.Contains("DoubleEG")) {
        if (passHLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ) passtriggerselection=true;
        if (passHLT_Ele27_WPTight_Gsf||passHLT_IsoMu24||passHLT_IsoTkMu24) passtriggerselection=false;//Drop repeated events
    }
    
    //DoubleMuon
    if (Dataset.Contains("DoubleMuon")) {
        if (passHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL||passHLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL) passtriggerselection=true;
        if (Dataset.Contains("RunH")){
            if (passHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ||passHLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ) passtriggerselection=true;
        }
        if (passHLT_Ele27_WPTight_Gsf||passHLT_IsoMu24||passHLT_IsoTkMu24||passHLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ) passtriggerselection=false;//Drop repeated events
    }
    
    //MuonEG
    if (Dataset.Contains("MuonEG")) {
        if (passHLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL||passHLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL) passtriggerselection=true;
        if (Dataset.Contains("RunH")){
            if (passHLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ||passHLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ) passtriggerselection=true;
        }
        //Drop repeated events
        if (passHLT_Ele27_WPTight_Gsf||passHLT_IsoMu24||passHLT_IsoTkMu24||passHLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ||passHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL||passHLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL) passtriggerselection=false;//Drop repeated events
        if (Dataset.Contains("RunH")){
            if (passHLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ||passHLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ) passtriggerselection=false;
        }
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
        //if (abs((*_tau_dz)[n])<0.2&&abs((*_ipPV)[n])<1000&&(*_lPt)[n]>20&&abs((*_lEta)[n])<2.3&&(*_TauIs_decayModeFinding)[n]&&(*_TauIs_byLooseIsolationMVArun2v1DBdR03oldDMwLT)[n])passobjselection=true;
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
