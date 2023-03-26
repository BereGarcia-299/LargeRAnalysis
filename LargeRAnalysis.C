#include <assert.h>
#include <cmath>
#include <fstream>
#include <stdio.h>
#include <array>
#include <string>
#include <tgmath.h>
#include <iostream>
#include <map>
#include <TChain.h>
#include <TLegend.h>
#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TMath.h"
#include "TError.h"
#include "TBox.h"
#include "TRandom3.h"
#include "getLogBins.h"
#include <TAttMarker.h>
#include <TCanvas.h>
#include <TPaveText.h>
#include "TLorentzVector.h"
#include "bNec.h"
#include "bFunctions.h"
//#include "largeRAnalysisBinning.h"
#include <TSystem.h>
//Centrality Weights
#include "centWeights.h"
#include "jZWeights.h"
#include "logBinsChris.h"
#include "ghostUtil.h"
#include "centralityFromInput.C"
#include "checkMakeDir.C"
using namespace std;

int LargeRAnalysis(string collisionType = "pp", int jet_Rad = R4,int data_or_mc = mcOrdata::mc,  string extraTag="MC", bool debug =false, int sys_uncrt =-1, string jer_or_jes_Sys = "",bool addpTShapeWeight = true,bool jetRate2015Bins=false, bool rAA2015Binning = false, bool dijet2018bins = true, bool officialLargeRpTBins = false){

  cout << "dijet2018bins: " << dijet2018bins << endl;
  cout << "rAA2015Binning: " << rAA2015Binning << endl;
  cout << "jetRate2015Bins: " << jetRate2015Bins << endl;

  TH1::SetDefaultSumw2();
  bool dgb  = false;
  //If you set sumpTCut_MtchTrthJets to truth then you will collect
  //truth jets matched to reco jets that have passed the sumpT cut. 
  bool sumpTCut_MtchTrthJets = true;
  //ptshape weights info
  double lowEndCutOff = 158; 
  double highEndCutOff =600;
  
  


  cout << "This is the jet radius you will be looking at: " << JetRadius[jet_Rad] << endl;
  
  if((jer_or_jes_Sys == "JER" || jer_or_jes_Sys == "JES") && sys_uncrt != -1) extraTag =Form("%s_%d_",jer_or_jes_Sys.c_str(),sys_uncrt) +extraTag;
  extraTag = Form("eta%d_",(int)(etacut*10)) + extraTag;

  
  if(jetRate2015Bins){
    extraTag = "2015JetRateBins_" + extraTag;
  }else if(rAA2015Binning){
    extraTag = "2015RAARateBins_" + extraTag;
  }else if(dijet2018bins){
    extraTag = "2018DiJetBins_" + extraTag;
  }else if(officialLargeRpTBins){
    extraTag = "2018LargerAnalysis_" + extraTag;
  }
  bool matchReq = true;
  if(data_or_mc == mcOrdata::data)matchReq =false;

  
  //Choosing Number of centralities
  std::vector<int> centBins{0,10,20,30,40,50,60,70,80};
  int number_CentBins = centBins.size()-1; 
  if(dijet2018bins)number_CentBins=4;



  int dataSet = 1; //pp
  bool addCentWeights = false; //No Weights for collission type: pp
  bool addpTShapeWght = false;
  bool sumpTReq = true;
  bool clean_ppJets = true; //Apply Jet Cleaning tool

  if(data_or_mc == mcOrdata::mc)addpTShapeWght=addpTShapeWeight;
  
  if(collisionType == "pp" && jet_Rad==R10 && data_or_mc == mcOrdata::mc && addpTShapeWeight){
    highEndCutOff = 700;
  }
		       
  if(collisionType == "PbPb"){
    dataSet = 0;// PbPb
    sumpTReq = false;
    clean_ppJets = false;
    if(data_or_mc == mc){
      addCentWeights = true;
      addpTShapeWght = addpTShapeWeight;
      highEndCutOff = 700;
    }
  }

  if(!addpTShapeWght){
    extraTag = "NopTShpWghts_" + extraTag ;
   }
  


  //FULL/Half Closure Tests
  bool full_and_ClsrTst = true; //Collect TH1Ds + TH2D for Tests
  
  //------------------
  //Systematics 
  //-------------------
  if(sys_uncrt!=-1){
    cout <<"You go through systematics for JER/JES!" << endl;
    
  }

  //-------------------------------------
  //sumTrackpT Distributions (Only for pp)
  //--------------------------------------
  bool sumpTTrackDis = false;
  float pTRange_SumpTDis[] = {600,800,1000};
  TH1D *sumpTTrack_Dis[totRadii][2];//[totRadii][pTRanges]
   

  //This is for efficiency studies.
  bool effiBins = false;
  if(effiBins){
    //These bins are applies to the
    //Efficiency plots (i.e. sumpT Cuts Eff.)
    const int totBins = 190;
    int total_bins = totBins-1;
    float finerBins_Eff[totRadii][totBins]={};
 
    // float up
    int uperedge[] = {260,140};//bins size of 2 up to this GeV value
    float binfivewidths[] = {350,200};
    for(int ijetr =0; ijetr < totRadii; ijetr++){
      int x_lowedge = 0;
      float size_of_bin = 2;
      for(int iBin = 1; iBin < totBins; iBin++){
	finerBins_Eff[ijetr][iBin] = x_lowedge + size_of_bin;
	x_lowedge = finerBins_Eff[ijetr][iBin];
	cout << "This was savved in your pT bin: " << finerBins_Eff[ijetr][iBin] << endl;
	if(finerBins_Eff[ijetr][iBin] == uperedge[ijetr]){
	  size_of_bin = 5;
	}else if(finerBins_Eff[ijetr][iBin] >binfivewidths[ijetr])
	  size_of_bin = 80;
      }//bin loop
    }//jet radius loop
  }

  //If you want to look at jets that were rejected by sumpT cut
  //We want to see how many fake jets are rejected when we apply the sumpT cut 
  //but also how many we missed
  bool fakeJets_Study=false;


  


  //Binning for this analysis
  //Reco Binning
  std::vector<std::vector<std::vector<double>>> binsGenR(totRadii, std::vector<std::vector<double>>(8));
  //Truth Binning
  std::vector<std::vector<std::vector<double>>> binsGenRTruth(totRadii, std::vector<std::vector<double>>(8));
//std::vector<std::vector<std::vector<double>>> binsGenR;
  // binsGenR.reserve(totRadii);
  // cout << __LINE__ << endl;
  // cout << "totRadii" << totRadii << endl;
  // for(int iRad =0; iRad < totRadii; iRad++){
  //   binsGenR[iRad].reserve(8);
  //   cout << __LINE__ << endl;
  // }

  // cout << binsGenR.size() << " , " << binsGenR[0].size() << " , " << binsGenR[1].size() << endl;
  // cout << "R4/R10 = " << R4 << "/" << R10 << endl;
  
  int pTBinsTotRTruth[2][8] = {{13,14,14,13,13,14,13,13},{22,22,22,19,19,20,19,19}};
  int pTBinsTotR[2][8] = {{11,12,12,11,11,12,11,11},{22,22,22,19,19,20,19,19}};
  //R=0.4
  //Truth Binning
  std::vector<std::vector<double>> binsR4_Truth;

  

  //R=0.4
  //Truth Binning
  binsR4_Truth.reserve(8);
  binsR4_Truth.emplace_back(std::vector<double>{66,83,100, 117, 138, 158, 178, 199, 220, 241, 262, 294, 336, 384, 439, 502, 573, 656, 750, 857, 980, 1120,1300}); //0-10%
  binsR4_Truth.emplace_back(std::vector<double>{66,83,100, 117, 138, 158, 178, 199, 220, 241, 262, 294, 336, 384, 439, 502, 573, 656, 750, 857, 980, 1120, 1300}); //10-20%
  binsR4_Truth.emplace_back(std::vector<double>{75,92,109, 126, 143, 160, 177, 194, 211, 228, 245, 262, 294, 336, 384, 439, 502, 573, 656, 750, 857, 980, 1120}); //20-30%
  binsR4_Truth.emplace_back(std::vector<double>{66,83,100, 117, 138, 158, 177, 198, 230, 262, 294, 336, 384, 439, 502, 573, 656, 750, 857, 1000});//30-40%
  binsR4_Truth.emplace_back(std::vector<double>{66,83,100, 117, 138, 158, 177, 198, 230, 262, 294, 336, 384, 439, 502, 573, 656, 750, 857, 1000}); //40-50%
  binsR4_Truth.emplace_back(std::vector<double>{66,83,100, 117, 138, 158, 177, 198, 220, 241, 262, 294, 336, 384, 439, 502, 573, 656, 750, 857, 1000}); //50-60%
  binsR4_Truth.emplace_back(std::vector<double>{66,83,100, 117, 138, 158, 177, 197, 220, 241, 262, 294, 336, 384, 439, 502, 573, 656, 750, 1000}); //60-70%
  binsR4_Truth.emplace_back(std::vector<double>{66,83,100, 117, 138, 158, 177, 197, 220, 241, 262, 294, 336, 384, 439, 502, 573, 656, 750, 1000}); //70-80%
  //Reco Binning
  std::vector<std::vector<double>> binsR4;
  binsR4.reserve(8);
  binsR4.emplace_back(std::vector<double>{100, 117, 138, 158, 178, 199, 220, 241, 262, 294, 336, 384, 439, 502, 573, 656, 750, 857, 980, 1120,1300}); //0-10%
  binsR4.emplace_back(std::vector<double>{100, 117, 138, 158, 178, 199, 220, 241, 262, 294, 336, 384, 439, 502, 573, 656, 750, 857, 980, 1120, 1300}); //10-20%
  binsR4.emplace_back(std::vector<double>{109, 126, 143, 160, 177, 194, 211, 228, 245, 262, 294, 336, 384, 439, 502, 573, 656, 750, 857, 980, 1120}); //20-30%
  binsR4.emplace_back(std::vector<double>{100, 117, 138, 158, 177, 198, 230, 262, 294, 336, 384, 439, 502, 573, 656, 750, 857, 1000}); //30-40%
  binsR4.emplace_back(std::vector<double>{100, 117, 138, 158, 177, 198, 230, 262, 294, 336, 384, 439, 502, 573, 656, 750, 857, 1000}); //40-50%
  binsR4.emplace_back(std::vector<double>{100, 117, 138, 158, 177, 198, 220, 241, 262, 294, 336, 384, 439, 502, 573, 656, 750, 857, 1000}); //50-60%
  binsR4.emplace_back(std::vector<double>{100, 117, 138, 158, 177, 197, 220, 241, 262, 294, 336, 384, 439, 502, 573, 656, 750, 1000}); //60-70%
  binsR4.emplace_back(std::vector<double>{100, 117, 138, 158, 177, 197, 220, 241, 262, 294, 336, 384, 439, 502, 573, 656, 750, 1000}); //70-80~
        

  //Reco Binning
  //Reco Binning
  std::vector<std::vector<double>> binsR10;
  binsR10.reserve(8);
  binsR10.emplace_back(std::vector<double>{294, 336, 384, 439, 502, 573, 656, 750, 857, 980, 1120,1300}); //0-10%
  binsR10.emplace_back(std::vector<double>{262, 294, 336, 384, 439, 502, 573, 656, 750, 857, 980, 1120, 1300}); //10-20%
  binsR10.emplace_back(std::vector<double>{245, 262, 294, 336, 384, 439, 502, 573, 656, 750, 857, 980, 1120}); //20-30%
  binsR10.emplace_back(std::vector<double>{230, 262, 294, 336, 384, 439, 502, 573, 656, 750, 857, 1000}); //30-40%
  binsR10.emplace_back(std::vector<double>{230, 262, 294, 336, 384, 439, 502, 573, 656, 750, 857, 1000}); //40-50%
  binsR10.emplace_back(std::vector<double>{220, 241, 262, 294, 336, 384, 439, 502, 573, 656, 750, 857, 1000}); //50-60%
  binsR10.emplace_back(std::vector<double>{220, 241, 262, 294, 336, 384, 439, 502, 573, 656, 750, 1000}); //60-70%
  binsR10.emplace_back(std::vector<double>{220, 241, 262, 294, 336, 384, 439, 502, 573, 656, 750, 1000}); //70-80%

  std::vector<std::vector<double>> binsR10_Truth;
  binsR10_Truth.reserve(8);
  binsR10_Truth.emplace_back(std::vector<double>{241, 262, 294, 336, 384, 439, 502, 573, 656, 750, 857, 980, 1120,1300}); //0-10%
  binsR10_Truth.emplace_back(std::vector<double>{220 ,241, 262, 294, 336, 384, 439, 502, 573, 656, 750, 857, 980, 1120, 1300}); //10-20%
  binsR10_Truth.emplace_back(std::vector<double>{211, 228, 245, 262, 294, 336, 384, 439, 502, 573, 656, 750, 857, 980, 1120}); //20-30%
  binsR10_Truth.emplace_back(std::vector<double>{177, 198, 230, 262, 294, 336, 384, 439, 502, 573, 656, 750, 857, 1000}); //30-40%
 binsR10_Truth.emplace_back(std::vector<double>{177,198, 230, 262, 294, 336, 384, 439, 502, 573, 656, 750, 857, 1000}); //40-50%
 binsR10_Truth.emplace_back(std::vector<double>{177, 198, 220, 241, 262, 294, 336, 384, 439, 502, 573, 656, 750, 857, 1000}); //50-60%
  binsR10_Truth.emplace_back(std::vector<double>{177, 197, 220, 241, 262, 294, 336, 384, 439, 502, 573, 656, 750, 1000}); //60-70%
 binsR10_Truth.emplace_back(std::vector<double>{177, 197, 220, 241, 262, 294, 336, 384, 439, 502, 573, 656, 750, 1000}); //70-80%

  for(int iRad =0; iRad < totRadii; iRad++){
    for(int icent =0; icent < 8; icent++){
      if(iRad==R4){
        int binSize =binsR4[icent].size();
        for(int ipt=0; ipt < binSize; ipt++){
          binsGenR[R4][icent].push_back(binsR4[icent][ipt]);
	  binsGenRTruth[R4][icent].push_back(binsR4_Truth[icent][ipt]);
        }
      } else if(iRad==R10){
        int binSize =binsR10[icent].size();
        for(int ipt=0; ipt < binSize; ipt++){
          binsGenR[R10][icent].push_back(binsR10[icent][ipt]);
	  binsGenRTruth[R10][icent].push_back(binsR10_Truth[icent][ipt]);
        }
      } 
    }
  }


  if(debug)cout << __LINE__ << endl;
  //saving pT cut
  double pTCutsTruth[totRadii][8] = {};
  //saving reco cut
  double pTCutsReco[totRadii][8] = {};

  double minpTCuts[totRadii] = {10000,10000};

  cout << binsGenR.size() << " , " << binsGenR[0].size() << endl;
  
  
  
  
   for(int iCentBin = 0; iCentBin < number_CentBins; iCentBin++){
     //cout << "This is the pT cut for R=" << JetRadius[jet_Rad] << " jets Truth/Reco: " << binsGenRTruth[jet_Rad].at(iCentBin)[0] << "/" << binsGenR[jet_Rad].at(iCentBin)[0] << endl;
     if(!dijet2018bins&&(binsGenRTruth[jet_Rad].at(iCentBin)[0] < minpTCuts[jet_Rad]))minpTCuts[jet_Rad] = binsGenRTruth[jet_Rad].at(iCentBin)[0];
     if(dijet2018bins){
       minpTCuts[jet_Rad] = dijetBins[iCentBin][0];
     } 
      if(data_or_mc == mc){
	pTCutsTruth[jet_Rad][iCentBin] = binsGenRTruth[jet_Rad].at(iCentBin)[0];
	cout << "Truth pT cut: " << pTCutsTruth[jet_Rad][iCentBin] << " for cent Bin: " << iCentBin << " and radius: " << JetRadius[jet_Rad] <<endl;
	if(dijet2018bins){
	   pTCutsTruth[jet_Rad][iCentBin] = dijetBins[iCentBin][0];
	}

      }

      
      pTCutsReco[jet_Rad][iCentBin] = binsGenR[jet_Rad].at(iCentBin)[0];

      if(dijet2018bins)pTCutsReco[jet_Rad][iCentBin]=dijetBins[iCentBin][6];
    }//Cent Loop
  

  cout << "These are the minimum pT cuts for reconstructed jets: " << minpTCuts[R10] << "/" << minpTCuts[R4] << endl;

  for(int iCentBin =0; iCentBin < number_CentBins; iCentBin++){
    cout << "This is the pT CUt truth jets for R=0.4: " << pTCutsTruth[1][iCentBin] << endl;
    cout << "This is the pT cut reco jets for R=0.4: " << pTCutsReco[1][iCentBin] << endl;
  }


  //The collisionType variable will be 0 when runing over Pb+Pb 
  //THe collisionType variable will be 1 when running over pp

   if(debug)cout << __LINE__ << endl;

   //Debugging
   TH1D *R4_cent0_10TruthMathedJ = new TH1D("R4_cent0_10TruthMathedJ","R4_cent0_10TruthMathedJ",50,0,10000);
   TH1D *R4_cent0_10RecoMathedJ = new TH1D("R4_cent0_10RecoMathedJ","R4_cent0_10RecoMathedJ",50,0,10000);
   TH2D *R4_cent0_10_pTShpWght = new TH2D("R4_cent0_10_pTShpWght","R4_cent0_10_pTShpWght",50,0,1000,50,0,3.5);




  //int number_CentBins = binsGenR[0].size();
  cout << "Number of cent bins: " << number_CentBins << endl;
  centralityFromInput centTable("/usatlas/u/bereniceg299/data/LargeRJet_Study/NewSourceCode/centrality_cuts_Gv32_proposed_RCMOD2.txt");
  if(debug) centTable.PrintTableTex();
  TH2D *responseMatrix[totRadii][number_CentBins]; //Full Closure
  TH2D *halfResponseMtx[totRadii][number_CentBins];//Half Closure
  TH1D *pTDisHalvesTruth[twoHalves][totRadii][number_CentBins];
  TH1D *pTDisHalvesReco[twoHalves][totRadii][number_CentBins];

  //pT Distributions for each difference the various JZ samples
  bool jz_samples_pTDis = false;
  TH1D *pTDis_JZSamples[totRadii][total_samples][2]; //0 is for truth and 1 is reco
  TH1D* pTDis[totRadii][number_CentBins];
  TH1D *pTDisTruth[totRadii][number_CentBins];
  TH1D *pTDisTruth_All[totRadii][number_CentBins];


  //SumET Distribution
  TH1D* sumETDisData[totRadii];
  TH1D* pTLeadJet_R10[totRadii][number_CentBins];
  //This is just to check the number of counts in each pT Bin
  //This is to help us decide how wide our pT bins should be 
  //To see how much we are limited by statistis 

  //if(fakeJets_Study) number_CentBins=1;
  //This is to study the fake jets that were rejected in pp Data
  TH1D* pTDis_RJ_sumpTCut[totRadii]; 
  TH1D* pTDis_RJ_cleanTool[totRadii];

  TH1D *exp_pTBins[totRadii][number_CentBins];
  
  //pTBin (Reco/Truth) Just for 0-10% Bin for 200-250 Truth pT bin
  TH1D* pTRatio[totRadii];

  //This is the eta-phi plot for 'noisy jets' 
  bool etaPhiPlots = false;
  TH2D *etaphiNoisyJets[totRadii];
  TH2D *etaphiSumpTCut[totRadii];
  
  //Eta/Phi Distributions
  bool etaPhiplotsOnly = false;
  const int paramDis =2;
  const int numCentbins = 8;
  const int truthAndReco = 2; //if running over data you will only have one instead of two
  int totType =2; 
  if(data_or_mc == mcOrdata::data)totType=1;
  TH1D *eta_phi_dis[totRadii][paramDis][numCentbins][truthAndReco];
  TH2D *etaphiPlots[totRadii][numCentbins];//Only collecting it for reco jets
   if(debug)cout << __LINE__ << endl;
  //This is to see how much statistics we have in multi-jet system
  //Check how many jets of R=1.0 jets do we have in pp that pT > 400 GeV
  //Creating This Root file to add in all the jets are passing the eta, cleanning, sumpT, etc cuts (pp data)
  bool collectBadJetIndicies = false;
  TFile *BadJetsFile;
  if(collectBadJetIndicies){
    if(data_or_mc == mcOrdata::data)BadJetsFile = new TFile("multiJets/BadJetsIndexR10.root","RECREATE");
    if(data_or_mc == mcOrdata::mc)BadJetsFile= new TFile("multiJets/BadJetsIndexR10MC.root","RECREATE");
  }
  std::vector<int> *BadJets = NULL; //do not pass the sumpT cut and cleaning cut tool 
  int BadJets_n =0;
  if(debug)cout << __LINE__ << endl;
  TTree* BadJetsTree;
  if(collectBadJetIndicies){
    BadJetsTree= new TTree("BadJetsTree","BadJetsTree");
    BadJetsTree->SetDirectory(BadJetsFile);
  }
 if(debug)cout << __LINE__ << endl;
  //Reading File that contains indicies of good jets (read above)
 
  TChain* GrabBadJetsTree = new TChain("BadJetsTree");
  GrabBadJetsTree->Add("multiJets/BadJetsIndexR10.root");
  
  bool checkMultiJetStats = false;
  vector<int> subLeadingJets;//saving indicies of the two subleading jets
  TH1D *pT_R10_400GeV = new TH1D("pT_R10_400GeV","pT_R10_400GeV",25,400,1000);
  TH1D *pT_R10_400GeV_2pT200 = new TH1D("pT_R10_400GeV_2pT200","pT_R10_400GeV_2pT200",25,400,1000);
  TH1D *pT_R10_400GeV_Only2_or3jets = new TH1D("pT_R10_400GeV_Only2_or3jets","pT_R10_400GeV_Only2_or3jets",25,400,1000);
  const int totLeadingpTRanges = 4;
  double leadingpTRanges[5] = {400,500,600,700,900}; 
  TH1D *GaussianDis[totLeadingpTRanges];
  TH1D *AJObservable[totLeadingpTRanges];

  //This Containe sumpT Values
  TChain *sumpT_tree = new TChain("grapes/sumpTValues_ForJetRadii_Data.root");
  

  if(trigEff && (collisionType == "PbPb")){
    for(int iTrig = 0; iTrig < number_CentBins; iTrig++){
      for(int iCentBin =0; iCentBin < number_CentBins; iCentBin++){
        pTLeadJet_R10[iTrig][iCentBin] = new TH1D(Form("pTLeadDis_%s_Cent_%s",trigNames[iTrig].c_str(),centBins_TrigEff[iCentBin].c_str()),Form("pTLeadDis_%s_Cent_%s",trigNames[iTrig].c_str(),centBins_TrigEff[iCentBin].c_str()), 140, 0, 1200);
      }//
    }//Trigger Loop
  }//Trigger Efficiency Plots for HLT_J85


  
    if(sumpTTrackDis){
      double arraypTBins[] = {0,1,2,3,4,5,10,20,30,40,50};
      for(int ipTRange =0; ipTRange < 2; ipTRange++){
	sumpTTrack_Dis[jet_Rad][ipTRange] = new TH1D(Form("sumpT_Dis_R%d_pp_TightPrimary_%d_%d",JetRadius[jet_Rad],(int)pTRange_SumpTDis[ipTRange],(int)pTRange_SumpTDis[ipTRange+1]),Form("sumpT_Dis_R%d_pp_TightPrimary_%d_%d",JetRadius[jet_Rad],(int)pTRange_SumpTDis[ipTRange],(int)pTRange_SumpTDis[ipTRange+1]),10,arraypTBins);
      }
    }

  if(jz_samples_pTDis){
      for(int iJZSamp =0; iJZSamp < total_samples; iJZSamp++){
	pTDis_JZSamples[jet_Rad][iJZSamp][0] = new TH1D(Form("JZ%d_MatchedTruth_R%d",iJZSamp+2,JetRadius[jet_Rad]),Form("JZ%d_MatchedTruth_R%d",iJZSamp+2,JetRadius[jet_Rad]),30,20,1300);
	pTDis_JZSamples[jet_Rad][iJZSamp][1] = new TH1D(Form("JZ%d_MatchedReco_R%d",iJZSamp+2,JetRadius[jet_Rad]),Form("JZ%d_MatchedReco_R%d",iJZSamp+2,JetRadius[jet_Rad]),30,20,1300);
	
      }
    }
    if(checkMultiJetStats && jet_Rad==R10){
      for(int ipTRange =0; ipTRange < totLeadingpTRanges; ipTRange++){
	GaussianDis[ipTRange] = new TH1D(Form("sumSubLeadingpVector_LeadingJetpT_%d_%d_R10",(int)leadingpTRanges[ipTRange],(int)leadingpTRanges[ipTRange+1]),Form("sumSubLeadingpVector_LeadingJetpT_%d_%d_R10",(int)leadingpTRanges[ipTRange],(int)leadingpTRanges[ipTRange+1]),45,0,2);
	AJObservable[ipTRange] = new TH1D(Form("AJ_LeadingJetpT_%d_%d_R10",(int)leadingpTRanges[ipTRange],(int)leadingpTRanges[ipTRange+1]),Form("AJ_LeadingJetpT_%d_%d_R10",(int)leadingpTRanges[ipTRange],(int)leadingpTRanges[ipTRange+1]),20,-1,1);
      }

    }
    pTRatio[jet_Rad] = new TH1D(Form("200_250GeV_JER_JES_R%d",JetRadius[jet_Rad]),Form("200_250GeV_JER_JES_R%d",JetRadius[jet_Rad]),40,0,3);

    if(sumETHist && collisionType == "PbPb")sumETDisData[jet_Rad]= new TH1D(Form("sumETDisData_R%d",JetRadius[jet_Rad]),Form("sumETDisData_R%d",JetRadius[jet_Rad]),200,0,maxSumETValue);

    if(etaPhiPlots){
      etaphiNoisyJets[jet_Rad] = new TH2D(Form("R%d_EtaPhi_NoisyJets",JetRadius[jet_Rad]),Form("R%d_EtaPhi_NoisyJets",JetRadius[jet_Rad]),55,-3.14,3.14,55,-1.5,1.5);
      etaphiSumpTCut[jet_Rad] = new TH2D(Form("R%d_EtaPhi_sumpTCutNoisyJets",JetRadius[jet_Rad]),Form("R%d_EtaPhi_sumpTCutNoisyJets",JetRadius[jet_Rad]),55,-3.14,3.14,55,-1.5,1.5);
    }
    for(int iCent = 0; iCent <  number_CentBins; iCent++){
      double *array = NULL;
      double *arrayTruth = NULL;
      int numpTBins =0;
      int numpTBinsTruth = 0;

      
      
       string centBin_strg = "NAN";
       if(dijet2018bins){
	   centBin_strg =  centBins_DJ[iCent];
       }else{
	    centBin_strg = centBins_2015Meas[iCent];
       }


      if(jetRate2015Bins){
	array = jetRateMC;arrayTruth = jetRateMC;
	numpTBins =15;
        numpTBinsTruth = 15;
      }else if(rAA2015Binning){
	array = rAA_2015Bins[iCent];arrayTruth =rAA_2015Bins[iCent];
	numpTBins =20;
        numpTBinsTruth = 20;
      }else if(dijet2018bins){
	array = dijetBins[iCent]; arrayTruth = dijetBins[iCent];
	numpTBins =dijetBinsTot;
        numpTBinsTruth = dijetBinsTot;
      }else if(!jetRate2015Bins && !rAA2015Binning && !dijet2018bins){
	if(jet_Rad == R4){	
	  array = binsR4.at(iCent).data();
	  arrayTruth = binsR4_Truth.at(iCent).data();
	  numpTBins = binsR4.at(iCent).size()-1;
	  cout << "TOTAL NUMBER OF BINS FOR RECO: " << numpTBins << endl;
	  numpTBinsTruth = binsR4_Truth.at(iCent).size()-1;
	}

	if(jet_Rad == R10){
	  array = binsR10.at(iCent).data();
	  arrayTruth = binsR10_Truth.at(iCent).data();
	  numpTBins = binsR10.at(iCent).size() -1;
	  numpTBinsTruth = binsR10_Truth.at(iCent).size()-1;
	}

      }
       
      if(fakeJets_Study){
	pTDis_RJ_sumpTCut[jet_Rad] = new TH1D(Form("R%d_RejectedJets_sumpTCut",JetRadius[jet_Rad]),Form("R%d_RejectedJetssumpTCut",JetRadius[jet_Rad]),25,0,1800);
	pTDis_RJ_cleanTool[jet_Rad] = new TH1D(Form("R%d_RejectedJets_CleanTool",JetRadius[jet_Rad]),Form("R%d_RejectedJets_CleanTool",JetRadius[jet_Rad]),25,0,1800);
        pTDis[jet_Rad][iCent] = new TH1D(Form("R%d_CleanJets",JetRadius[jet_Rad],centBins_2015Meas[iCent].c_str()),Form("R%d_CleanJets",JetRadius[jet_Rad]),25,0,1800);
      }else{



	pTDis[jet_Rad][iCent] = new TH1D(Form("R%d_Cent_%s",JetRadius[jet_Rad],centBin_strg.c_str()),Form("R%d_Cent_%s",JetRadius[jet_Rad],centBin_strg.c_str()),numpTBins,array);
       

      }


      cout << "WE SHOULD GO THROUGH THE ETA PHI INIT PART!!!!!***" << endl;


      if(etaPhiplotsOnly){
	//etaphiPlots[jet_Rad][iCent] = new TH2D(Form("etaPhiplot_R%d_Cent_%s",JetRadius[jet_Rad],centBins_2015Meas[iCent].c_str()),Form("etaPhiplot_R%d_Cent_%s",JetRadius[jet_Rad],centBins_2015Meas[iCent].c_str()),45,

	for(int iparam = 0; iparam < totJetParam-1; iparam++){
	  for(int iTypeJet =0; iTypeJet < totType; iTypeJet++){
	    
	    if(debug)cout << "AM I GOING THROUGH HERE??? [jet_Rad][iparam][iCent][iTypeJet]: " << jet_Rad << "/" << iparam << "/" << iCent << "/" << iTypeJet << endl;
	    eta_phi_dis[jet_Rad][iparam][iCent][iTypeJet] = new TH1D(Form("R%d_Cent_%s_%sDis_%s",JetRadius[jet_Rad],centBins_2015Meas[iCent].c_str(),jetParam[iparam+1].c_str(),jetType[iTypeJet].c_str()),Form("R%d_Cent_%s_%sDis_%s",JetRadius[jet_Rad],centBins_2015Meas[iCent].c_str(),jetParam[iparam+1].c_str(),jetType[iTypeJet].c_str()),25,eta_phi_Edges[iparam][startEdge],eta_phi_Edges[iparam][endEdge]);
	    }//parameter loop (only 2 param)
	  }
      }				   

      if(data_or_mc == mc){
	if(full_and_ClsrTst){
	  if(debug)cout << __LINE__ << endl;


	  responseMatrix[jet_Rad][iCent]= new TH2D(Form("FullClsr_RespMatrix_R%d_Cent%s",JetRadius[jet_Rad],centBin_strg.c_str()),Form("FullClsr_RespMatrix_R%d_Cent%s",JetRadius[jet_Rad],centBin_strg.c_str()),numpTBins,array,numpTBinsTruth,arrayTruth);
	  if(debug)cout << __LINE__ << endl;
	  //this goes into first half 
	  halfResponseMtx[jet_Rad][iCent]= new TH2D(Form("HalfClsr_RespMatrix_R%d%s",JetRadius[jet_Rad],centBin_strg.c_str()),Form("HalfClsr_RespMatrix_R%d%s",JetRadius[jet_Rad],centBin_strg.c_str()),numpTBins,array,numpTBinsTruth,arrayTruth);
	  if(debug)cout << __LINE__ << endl;
	 for(int iHalf =0; iHalf <  twoHalves; iHalf++){
	   string tagHalf = "First";
	   if(iHalf>0)tagHalf = "Secnd";
	   if(debug)cout << __LINE__ << endl;
	   pTDisHalvesTruth[iHalf][jet_Rad][iCent]=new TH1D(Form("pTDis%sHalfMatch_TruthJets_R%d%s",tagHalf.c_str(),JetRadius[jet_Rad],centBin_strg.c_str()),Form("pTDis%sHalfMatch_TruthJets_R%d%s",tagHalf.c_str(),JetRadius[jet_Rad],centBin_strg.c_str()),numpTBinsTruth,arrayTruth);
	   if(debug)cout << __LINE__ << endl;
	   pTDisHalvesReco[iHalf][jet_Rad][iCent]=new TH1D(Form("pTDis%sHalfMatch_RecoJets_R%d%s",tagHalf.c_str(),JetRadius[jet_Rad],centBin_strg.c_str()),Form("pTDis%sHalfMatch_RecoJets_R%d%s",tagHalf.c_str(),JetRadius[jet_Rad],centBin_strg.c_str()),numpTBins,array);
	   if(debug)cout << __LINE__ << endl;
	 }//Looping over halves
	  	
	}
        if(matchReq){
	  if(debug)cout << __LINE__ << endl;
	  pTDisTruth[jet_Rad][iCent] = new TH1D(Form("R%d_Cent_%s_TruthJetsMatched",JetRadius[jet_Rad],centBin_strg.c_str()),Form("R%d_Cent_%s_TruthJetsMatched",JetRadius[jet_Rad],centBin_strg.c_str()),numpTBinsTruth,arrayTruth);
	  pTDisTruth_All[jet_Rad][iCent] = new TH1D(Form("R%d_Cent_%s_AllTruthJets",JetRadius[jet_Rad],centBin_strg.c_str()),Form("R%d_Cent_%s_TruthJetsMatched",JetRadius[jet_Rad],centBin_strg.c_str()),numpTBinsTruth,arrayTruth);

	  if(debug)cout << __LINE__ << endl;
	}else if(!matchReq){
	  pTDisTruth[jet_Rad][iCent] = new TH1D(Form("R%d_Cent_%s_ALLTruthJets",JetRadius[jet_Rad],centBin_strg.c_str()),Form("R%d_Cent_%s_ALLTruthJets",JetRadius[jet_Rad],centBin_strg.c_str()),numpTBinsTruth,arrayTruth);
	  
	}
      }
    }//Centrality Loop
  

  ///
  if(data_or_mc == mcOrdata::mc){
    for(int iJZsamp = 0; iJZsamp < total_samples; iJZsamp++){
      if(collisionType=="PbPb")cout << "these are the weights: " << jZWeights[iJZsamp] << endl;
      if(collisionType=="pp")cout << "these are the weights: " << weights_ppMC[iJZsamp] << endl;
    }
  }


  //------Grabbing All Root Files
  if(data_or_mc == mcOrdata::data){
    //Pb+Pb Data Set
    cout <<__LINE__ << endl;
    if(collisionType == "PbPb"){
      cout <<__LINE__ << endl;
      //tree->Add("/gpfs/mnt/atlasgpfs01/usatlas/data/bereniceg299/user.berenice.05172022.PbPb2018Data_GRL_R10_R4.0000000000_myOutput.root/*.root"); //R=0.4 and R=1.0 jets are in TTrees
      tree->Add("/gpfs/mnt/atlasgpfs01/usatlas/data/bereniceg299/Clown/user.berenice.03132023.LargeRJet_PbPbData_R10_R04_R02_.00000000002_myOutput.root/*.root");
      
      cout <<__LINE__ << endl;
      cout << "This is how many entries you have: " << tree->GetEntries() << endl;
    }
    cout <<__LINE__ << endl;
    //if(collisionType == "pp")tree->Add("/gpfs/mnt/atlasgpfs01/usatlas/data/bereniceg299/user.berenice.01212023.LargeRJet_ppData_R10_R04_R02_InsituCalibrationsForALL.00000000001_myOutput.root/*.root");
    if(collisionType == "pp")tree->Add("/gpfs/mnt/atlasgpfs01/usatlas/data/bereniceg299/user.berenice.03132023.LargeRJet_ppData_R10_R04_R02_.00000000000_myOutput.root/*.root");

    cout << "This is how many entries you have: " << tree->GetEntries() << endl;
  }else if(data_or_mc == mc){

    cout <<__LINE__ << endl;
    //if(collisionType=="PbPb")directory = "PbPbMC_Uncert";
    if(collisionType=="PbPb")directory = "PbPbMC_Uncert_1";
    //if(collisionType=="pp")directory = "ppMC";
    
    if(collisionType=="pp")directory = "ppMC_1";
    cout <<__LINE__ << endl;
    int PbPbMC_JZ_range_entries[6] = {};
    cout << __LINE__ << endl;
    for(int iJZ_Samples = 0; iJZ_Samples < total_samples; iJZ_Samples++){
      string PbPbtag = "";
      cout <<"current PbPb tag: "<< PbPbtag << endl;
      //if(collisionType == "PbPb")PbPbtag = "_new";
      
      cout <<"current PbPb tag: "<< PbPbtag << endl;
      cout << Form("~/usatlasdata/%s/JZ%d%s/*.root",directory.c_str(),iJZ_Samples+2,PbPbtag.c_str()) << endl;
      tree->Add(Form("~/usatlasdata/%s/JZ%d%s/*.root",directory.c_str(),iJZ_Samples+2,PbPbtag.c_str()));
      cout <<__LINE__ << endl;
      num_JZ_range_entries[iJZ_Samples+1]=tree->GetEntries();

      cout << "inside jz array: " << tree->GetEntries() << endl;
    }//---Picking up all JZ samples

  }

  //We are uploading root filse that contain weights
  //MinBis -> Bias weights
  TFile *centWeightFiles[totRadii];
  TH1D *centWeightsFunc[totRadii];
  cout <<__LINE__ << endl;
  if(addCentWeights){
    
      //Since we are only using one trigger (HLT_j85) we will use one file and that's the R4 file.
      centWeightFiles[jet_Rad]= new TFile("centWeights_PbPbMC_R4.root","READ");
      centWeightsFunc[jet_Rad]=(TH1D*) centWeightFiles[jet_Rad]->Get("ratio_R4");
    
  }
  //pT shape weights
  TFile* ptshape_wghtsFile[totRadii];
  TF1 *ptshapeWghts[totRadii][number_CentBins];
  string ptshapeweights_tag = "erf";
  if(dijet2018bins)ptshapeweights_tag="substructComp";
if(addpTShapeWght){
    if(collisionType=="PbPb"){
      
      ptshape_wghtsFile[jet_Rad] = new TFile(Form("/usatlas/u/bereniceg299/data/LargeRJet_Study/NewSourceCode/LargeRAnalysis/pt_shp_wghts/%s_ptshape_weights_pbpb_R%d.root",ptshapeweights_tag.c_str(),JetRadius[jet_Rad]),"READ");
      
    }else if(collisionType=="pp"){
      ptshape_wghtsFile[jet_Rad] = new TFile(Form("/usatlas/u/bereniceg299/data/LargeRJet_Study/NewSourceCode/LargeRAnalysis/pt_shp_wghts/%s_ptshape_weights_pp_R%d.root",ptshapeweights_tag.c_str(),JetRadius[jet_Rad]),"READ"); //this is not the correcfile.
      
    }

      
     for(int iCentBin =0; iCentBin < number_CentBins;iCentBin++){

       string centBin_strg = "NAN";
       if(dijet2018bins){
	 centBin_strg =  centBins_DJ[iCentBin];
       }else{
	 centBin_strg = centBins_2015Meas[iCentBin];
       }

       ptshapeWghts[jet_Rad][iCentBin] = (TF1*) ptshape_wghtsFile[jet_Rad]->Get(Form("fitfun_ptshape_weights_Cent_%s",centBin_strg.c_str())); 
     }//pt shape weights
    
    
  }
  //This is for MultiJets
  if(collectBadJetIndicies){
    BadJetsTree->Branch("BadJets",&BadJets);
    BadJetsTree->Branch("BadJets_n", &BadJets_n, "BadJets_n/I");
  }

  //MultiJetStudy
  if(checkMultiJetStats){
    GrabBadJetsTree->SetBranchAddress("BadJets",&BadJets);
    GrabBadJetsTree->SetBranchAddress("BadJets_n", &BadJets_n);
  }
  
  //Turning Off all Branches (So it will not load all)
  tree->SetBranchStatus("*",0);


  //Turning on Event Number Branch
  tree->SetBranchStatus("EventNumber",1);
  tree->SetBranchAddress("EventNumber", &eventNumber);

  //Turning on Run Number Branch
  tree->SetBranchStatus("RunNumber",1);
  tree->SetBranchAddress("RunNumber", &runNumber);


  //Turning On Braches that you will use
  tree->SetBranchStatus("nvert",1);

  //---Verticies
  tree->SetBranchAddress("nvert", &nvert);


  if(collisionType == "PbPb"){
    if(data_or_mc == mcOrdata::data){
      tree->SetBranchStatus("HLT_j85_ion_L1J30",1);
      tree->SetBranchAddress("HLT_j85_ion_L1J30",&HLT_j85_ion_L1J30);
    }

    
    tree->SetBranchStatus("fcalA_et",1);
    tree->SetBranchStatus("fcalC_et",1);
    //FCal
    tree->SetBranchAddress("fcalA_et",&fcalEnergyA);
    tree->SetBranchAddress("fcalC_et",&fcalEnergyC);
  }else{
    tree->SetBranchStatus("vert_type",1);
    tree->SetBranchStatus("vert_sumpt",1);
    tree->SetBranchStatus("vert_ntrk", 1);
    //---Tracks
    tree->SetBranchStatus("ntrk",1);
    tree->SetBranchStatus("trk_eta", 1);
    tree->SetBranchStatus("trk_phi", 1);
    tree->SetBranchStatus("trk_pt", 1);
    if(data_or_mc == mcOrdata::data){
      //Triggers
      tree->SetBranchStatus("HLT_j85", 1);
      tree->SetBranchStatus("HLT_j110_a10_lcw_subjes_L1J30", 1);

      tree->SetBranchAddress("HLT_j85", &HLT_j85);
      tree->SetBranchAddress("HLT_j110_a10_lcw_subjes_L1J30", &HLT_j110_a10_lcw_subjes_L1J30);
      
      //Jet Cleaning Tool
      for(int iRad = 0; iRad < totRadii; iRad++){
	tree->SetBranchStatus(Form("LooseBad_%d",JetRadius[iRad]),1);
      	tree->SetBranchAddress(Form("LooseBad_%d",JetRadius[iRad]),&looseBadR[iRad]);
      }//jet radii loop      

    }
    //---Verticies 
    tree->SetBranchAddress("vert_type", vert_type);
    tree->SetBranchAddress("vert_sumpt", vert_sumpt);
    tree->SetBranchAddress("vert_ntrk", vert_ntrk);
    //---Tracks
    tree->SetBranchAddress("ntrk",&ntrk);
    tree->SetBranchAddress("trk_eta", &trk_eta);
    tree->SetBranchAddress("trk_phi", &trk_phi);
    tree->SetBranchAddress("trk_pt", &trk_pt);


  }

  string mc_truthTag = "";
  string mc_recoTag = "";
  string mc_recoTag_prfx = "";
  string data_tag = "";
  string jet_calib = "";

  if(data_or_mc == mcOrdata::data){
    jet_calib = "Insitu";
    data_tag = "hi";
  }else{
    mc_truthTag = "_truth";
    //mc_recoTag = "_reco";
    mc_recoTag = "";
    //mc_recoTag_prfx = "m_b_";
    mc_recoTag_prfx = "";
    data_tag = "hi";
    //jet_calib = "etaJES";
    jet_calib = "etajes";
  }

  
     
    if(debug)cout << "This is the truth jet number variable you grabbed for reco: " << Form("%sakt%d%s%s_jet_n",mc_recoTag_prfx.c_str(),JetRadius[jet_Rad],data_tag.c_str(),mc_recoTag.c_str()) << endl;
    tree->SetBranchStatus(Form("%sakt%d%s%s_jet_n",mc_recoTag_prfx.c_str(),JetRadius[jet_Rad],data_tag.c_str(),mc_recoTag.c_str()),1); //Turning on branch
    cout << "This is the jet radius number of jets branh: " << Form("%sakt%d%s%s_jet_n",mc_recoTag_prfx.c_str(),JetRadius[jet_Rad],data_tag.c_str(),mc_recoTag.c_str()) << endl;

    tree->SetBranchAddress(Form("%sakt%d%s%s_jet_n",mc_recoTag_prfx.c_str(),JetRadius[jet_Rad],data_tag.c_str(),mc_recoTag.c_str()),&akt_jet_n[jet_Rad]);

    if(checkMultiJetStats && jet_Rad==R10){
      tree->SetBranchStatus("akt10hi_jet_Insitu_calib_e",1);
      tree->SetBranchAddress("akt10hi_jet_Insitu_calib_e",&jetEnergy);
    }
    if(data_or_mc == mc){
      cout << "This is the reco jet number: " << Form("akt%d_truth_jet_n",JetRadius[jet_Rad]) << endl;

      tree->SetBranchStatus(Form("akt%d_truth_jet_n",JetRadius[jet_Rad]),1);
      tree->SetBranchAddress(Form("akt%d_truth_jet_n",JetRadius[jet_Rad]),&akt_truth_jet_n[jet_Rad]);
    }


    if(collisionType == "pp"){

      tree->SetBranchStatus(Form("akt%dhi_jet_sumTightPrimaryTrackpT",JetRadius[jet_Rad]),1);
      tree->SetBranchAddress(Form("akt%dhi_jet_sumTightPrimaryTrackpT",JetRadius[jet_Rad]),&akt_jet_sumTightPrimaryTrackpT[jet_Rad]);



    }
    

    for(int iParam =0; iParam < totJetParam; iParam++){
      
      string jetRecoBranch = Form("akt%d%s%s_jet_%s_calib_%s",JetRadius[jet_Rad],data_tag.c_str(),mc_recoTag.c_str(),jet_calib.c_str(),jetParam[iParam].c_str()); //branch name in data
 
      if(data_or_mc == mcOrdata::mc)jetRecoBranch = Form("akt%dhi_etajes_jet_%s",JetRadius[jet_Rad],jetParam[iParam].c_str()); //branch name in mc
      
      if(iParam == pTPar && sys_uncrt != -1){
	
	jetRecoBranch = Form("akt%dhi_etajes_jet_pt_sys_%s_%d",JetRadius[jet_Rad],jer_or_jes_Sys.c_str(),sys_uncrt);  
	  cout << "This is the branch you are looking for the systematical uncertainties: " << jetRecoBranch << endl;
      }

      if(debug)cout << "This is the jet branch used: " << jetRecoBranch.c_str() << endl;

      tree->SetBranchStatus(jetRecoBranch.c_str(),1);//Turning On Branch
      tree->SetBranchAddress(jetRecoBranch.c_str(), &akt_jet_kin[jet_Rad][iParam]);
      
      if(data_or_mc == mc){
        tree->SetBranchStatus(Form("akt%d_truth_jet_%s",JetRadius[jet_Rad],jetParam[iParam].c_str()),1);
        tree->SetBranchAddress(Form("akt%d_truth_jet_%s",JetRadius[jet_Rad],jetParam[iParam].c_str()),&akt_truth_jet_kin[jet_Rad][iParam]);
      }


    }//Parameter Loop
  
  if(debug)cout << __LINE__ << endl;
  TRandom3* rand_p;
  if(debug)cout << __LINE__ << endl;
  if(data_or_mc == mcOrdata::mc){
    if(debug)cout << __LINE__ << endl;
    if(collisionType=="PbPb"){
      rand_p = new TRandom3(12201);
    }else if(collisionType=="pp"){
      rand_p = new TRandom3(20211);
    }
  }

  if(debug)cout << __LINE__ << endl;
  if(debug)cout << "THis is the number of entries(): " << tree->GetEntries() << endl;
  int events_tot= tree->GetEntries();
  if(debug)cout << __LINE__ << endl;
  //events_tot =3249974;
  if(debug)cout << __LINE__ << endl;
  if(debug)cout << "This is the total amount of entries: " << events_tot << endl;
  if(debug)cout << __LINE__ << endl;
  //events_tot = 5000;
  for(int iEvent = 0; iEvent < events_tot; iEvent++){    
    
    float randomValue = 0;
    int corrHalf = 1;
    BadJets_n =0;
    
    tree->GetEntry(iEvent);

    if(data_or_mc == mcOrdata::mc){
      randomValue = rand_p->Uniform(0,1);
      
      if(randomValue>0.5)corrHalf =0;
    }


    
    if (iEvent % 1000000 == 0) printf("entry: %lu\n", iEvent);

    
    if(checkMultiJetStats)GrabBadJetsTree->GetEntry(iEvent);
   
    
    //cout << "This is the eventnumber: " << eventNumber << endl;
    
    // if(1262796547 != eventNumber){
    //     continue;
    //  }
    int numJetsTot= akt_jet_kin[jet_Rad][pTPar]->size();

    // cout << "This is the run number: " << runNumber << endl; 
    // cout << "This is the event number: " << eventNumber << endl;
    // cout << "This is the jet radius: " << akt_jet_n[jet_Rad] << endl;
    // cout << "This is the number of recojets in this event: " << akt_jet_n[jet_Rad] << endl;
    // cout << "This is the nvert: " << nvert << endl;
    

    if(nvert < 2)continue;
    
    // //FCal
    float sumET = 0;
    double centWeight = 1;

    if(collisionType == "PbPb"){

       sumET = fcalEnergyA + fcalEnergyC;
       if(debug)cout << "this is the value of your sumet:" << endl << endl;
       if(sumET > 5020)continue; //this is to get rid of pile up
       if(debug)cout << "NOT PILE UP!" << endl<<endl;

    }








    if(debug)cout << __LINE__ << endl;
    //Checking What JZ Sample This Event is in
    //And adding in the JZ weight
    int iJZBin = -1;
    double jZWeight = 1;

    if(data_or_mc == mc){

      for(int iJZ = 0; iJZ < total_samples; iJZ++){
        if(debug)cout << "This is the JZ sample for loop!" << endl;
        if( (iEvent < num_JZ_range_entries[iJZ]) || (iEvent > num_JZ_range_entries[iJZ+1])) continue;
        iJZBin = iJZ;
        break;
      }

      if(iJZBin == -1)continue;

    }
    if(debug)cout << __LINE__ << endl; 

    
    //JZ weight
    double jZweight = 1;
    if(data_or_mc == mc){
        if(debug)cout << "Grabbing the JZ weight! " << endl;

        if(collisionType == "PbPb"){
          jZweight = jZWeights[iJZBin];
        }else if(collisionType == "pp"){
          jZweight = weights_ppMC[iJZBin]; 
        }
	if(debug)cout << "This is the jz weight: " << jZweight << endl;
        if(debug)cout << "This is the JZ bin: " << iJZBin << endl;


      }
      


     
    // centrality
    Int_t centBin = -1;
    Double_t cent = -1;
    if(collisionType=="PbPb"){
      cent = centTable.GetCent(sumET);
      centBin = ghostPos(centBins, cent, true, false);
      if(debug){
	cout << "This is the cent bin: " << centBin << endl;
        cout << "centrality: " << cent << endl;
      }
    }
    else centBin = 0;
    if(centBin < 0 && collisionType=="PbPb") continue;

     if(data_or_mc == mc && collisionType == "PbPb" && addCentWeights){
       if(sumET < 4700)centWeight = centWeightsFunc[jet_Rad]->GetBinContent(centWeightsFunc[jet_Rad]->FindBin(sumET));
       if(sumET > 4700)centWeight = centWeightsFunc[jet_Rad]->GetBinContent(centWeightsFunc[jet_Rad]->FindBin(4700)); 
     }

     if(dijet2018bins){
     
       if(centBin == centIndex0_10)centBin=0;
       if(centBin > centIndex0_10 && centBin <= centIndex20_30)centBin=1;
       if(centBin >= centIndex30_40 && centBin <= centIndex40_50)centBin=2;
       if(centBin >= centIndex50_60 && centBin <= centIndex70_80)centBin=3;

     }
     
  //centIndex30_40,centIndex40_50,centIndex50_60




    if(debug)cout<< __LINE__ << endl;

    
      if(debug)cout << "This is the jet radius we are looking at: " << jet_Rad << endl;
      
      if((collisionType == "PbPb") && (data_or_mc == mcOrdata::data)){
	if(!HLT_j85_ion_L1J30){continue;} //We decided to use this HI trigger for both R=1.0 and 0.4
      }else if(collisionType =="pp" && data_or_mc != mc){
        if(jet_Rad==R10 && !HLT_j110_a10_lcw_subjes_L1J30)continue; //Only looking R=1.0 jets when this trigger is fired
        if(jet_Rad==R4 && !HLT_j85)continue; //Only looking at R=0.4 jets when this trigger is fired
      }


      



      
      double leading_pT_HLTa10 = -1;
      double leading_pT_HLTj85 = -1;




      if(jet_Rad==R4 && HLT_j85_ion_L1J30 && HLT_j150_a10_ion_L1J50 && trigEff){


        if(debug)cout << "should not go through here!!" << endl;
        if(debug)cout << __LINE__ << endl;
        for(int iJetR10 =0; iJetR10 < akt_jet_n[R10]; iJetR10++){

          if(akt_jet_kin[R10][etaPar]->at(iJetR10) > etacut || akt_jet_kin[R10][etaPar]->at(iJetR10) < -etacut){continue;}


          if(leading_pT_HLTj85 < akt_jet_kin[R10][pTPar]->at(iJetR10)){
            leading_pT_HLTj85 = akt_jet_kin[R10][pTPar]->at(iJetR10);
          }//HLT_j85 trigger fired

        }//R10 Loop

        if(leading_pT_HLTj85 != -1){
          //Saving the leading jet
          if(centBin == 3 || centBin ==4 ){
            pTLeadJet_R10[1][cent_bin_30_50]->Fill(leading_pT_HLTj85);
          }else if(centBin == 5 || centBin == 6 || centBin == 7){
            pTLeadJet_R10[1][cent_bin_50_80]->Fill(leading_pT_HLTj85);
          }else{
            pTLeadJet_R10[1][centBin]->Fill(leading_pT_HLTj85);
          }
        }
      }//Looking at R10 jets that fired the HLT_j85 trigger
      
      //indicies that will be saved if truth found a match
      map<int,int>TruthMatchIndicies;
            
      //give index of reco and true if it was match to give you the truth matched jet.

      //this map will save the indcies of reco jets
      //and return a bool (true if it was matched and false if it did not)
      
      
      if(data_or_mc == mc){

	if(debug) cout << " WE CANNOT BE IN HERE!!!" <<endl; 
	
	bool UniqueMatches[200] = {};
	if(debug)cout << __LINE__ << endl;
	if(debug)cout << "This is how many truth jets you have in this event:  " << akt_truth_jet_n[jet_Rad] << endl;
	if(debug)cout << "This is the jet radius we are looking at: " << jet_Rad << endl;
        if(debug){
	  cout << "Printing all truth jets in this event." << endl;
	  for(int itru =0; itru < akt_truth_jet_n[jet_Rad]; itru++){
	    cout << "TRUTH JET Index/pt/eta/phi: " << itru << "/"<< akt_truth_jet_kin[jet_Rad][pTPar]->at(itru) << "/" << akt_truth_jet_kin[jet_Rad][etaPar]->at(itru)<< "/" << akt_truth_jet_kin[jet_Rad][phiPar]->at(itru) << endl;
	  }
	}
	for (int iJetTruth = 0; iJetTruth < akt_truth_jet_n[jet_Rad]; iJetTruth++){
	  if(debug)cout << "This is the Truth Index/pT/eta/phi: " << iJetTruth << "/" <<akt_truth_jet_kin[jet_Rad][pTPar]->at(iJetTruth)  << " / " << akt_truth_jet_kin[jet_Rad][etaPar]->at(iJetTruth) << " / " << akt_truth_jet_kin[jet_Rad][phiPar]->at(iJetTruth) << endl; 
	  
	  if( (akt_truth_jet_kin[jet_Rad][etaPar]->at(iJetTruth) > etacut) || (akt_truth_jet_kin[jet_Rad][etaPar]->at(iJetTruth)< -etacut) ) continue;//eta cut
	  if(debug)cout << "Truth jet passed the eta cut! " << endl;
	  
	  
	  int pTTruthCut = 0;
	  int pTRecoCut = 0;
	  
	  if(collisionType == "PbPb"){
	     pTTruthCut = pTCutsTruth[jet_Rad][centBin];
	     pTRecoCut = pTCutsReco[jet_Rad][centBin];
	     
	     
	     if(debug)cout << "listing our truth/reco pT cuts: " << pTTruthCut << "/"<< pTRecoCut <<  endl; 
	    }
	    if(debug)cout << __LINE__ << endl;
	    
	  
	  if(debug)cout << "will this pT truth jet pass pT cut? (Index/pT)" << iJetTruth << "/" << akt_truth_jet_kin[jet_Rad][pTPar]->at(iJetTruth) << endl;
	  if(pTTruthCut > akt_truth_jet_kin[jet_Rad][pTPar]->at(iJetTruth) && collisionType=="PbPb")continue;
	  double ptShpWeight =1.0;
	  if(debug)cout << __LINE__ << endl;
	  if(collisionType == "PbPb"){
	    if(addpTShapeWeight){
	      if(debug)cout << __LINE__ << endl;
	      if(akt_truth_jet_kin[jet_Rad][pTPar]->at(iJetTruth)< lowEndCutOff){
		if(debug)cout << __LINE__ << endl;
		ptShpWeight=ptshapeWghts[jet_Rad][centBin]->Eval(lowEndCutOff+1);
		if(debug)cout << __LINE__ << endl;
	      }else if(akt_truth_jet_kin[jet_Rad][pTPar]->at(iJetTruth) > highEndCutOff){
		ptShpWeight=ptshapeWghts[jet_Rad][centBin]->Eval(lowEndCutOff-1);
	      }else{
		if(debug)cout << __LINE__ << endl;
		ptShpWeight=ptshapeWghts[jet_Rad][centBin]->Eval(akt_truth_jet_kin[jet_Rad][pTPar]->at(iJetTruth)); //instead of inputing reco we input truth

	      }
	    }
	    if(debug)cout << __LINE__ << endl;
	    pTDisTruth_All[jet_Rad][centBin]->Fill(akt_truth_jet_kin[jet_Rad][pTPar]->at(iJetTruth),centWeight*jZweight*ptShpWeight);
	  }else if(collisionType == "pp"){
	    for(int iCentBin=0; iCentBin <  number_CentBins; iCentBin++){
	      if(addpTShapeWeight){
		if(akt_truth_jet_kin[jet_Rad][pTPar]->at(iJetTruth)< lowEndCutOff){
		  ptShpWeight=ptshapeWghts[jet_Rad][iCentBin]->Eval(lowEndCutOff+1);
		}else if(akt_truth_jet_kin[jet_Rad][pTPar]->at(iJetTruth) > highEndCutOff){
		  ptShpWeight=ptshapeWghts[jet_Rad][iCentBin]->Eval(lowEndCutOff-1);
		}else{
		  ptShpWeight=ptshapeWghts[jet_Rad][iCentBin]->Eval(akt_truth_jet_kin[jet_Rad][pTPar]->at(iJetTruth)); //instead of inputing reco we input truth

		}
	      }
	      pTDisTruth_All[jet_Rad][iCentBin]->Fill(akt_truth_jet_kin[jet_Rad][pTPar]->at(iJetTruth),centWeight*jZweight*ptShpWeight);
	    }
	  }
	  
	  if(debug)cout << "YES!" << endl;
	  if(!matchReq){
	    //pass if MATCHING is NOT required
	    if(collisionType=="PbPb"){
	      
	      pTDisTruth[jet_Rad][centBin]->Fill(akt_truth_jet_kin[jet_Rad][pTPar]->at(iJetTruth),centWeight*jZweight);  
	    }else if(collisionType == "pp"){
	      if(debug)cout << "This truth jet will be used to fill the histogram that looks at all truth jets! " << endl;
	      if(debug)cout << "Index/pT/eta/phi: " << iJetTruth << "/" <<akt_truth_jet_kin[jet_Rad][pTPar]->at(iJetTruth) << "/" << akt_truth_jet_kin[jet_Rad][etaPar]->at(iJetTruth) << "/" << akt_truth_jet_kin[jet_Rad][etaPar]->at(iJetTruth) << "/" << akt_truth_jet_kin[jet_Rad][phiPar]->at(iJetTruth) << endl;
	      
 	      for(int iCentBin =0; iCentBin < number_CentBins; iCentBin++){
		if(pTCutsTruth[jet_Rad][iCentBin] > akt_truth_jet_kin[jet_Rad][pTPar]->at(iJetTruth))continue;
		
		pTDisTruth[jet_Rad][iCentBin]->Fill(akt_truth_jet_kin[jet_Rad][pTPar]->at(iJetTruth),jZweight);
	      }//centBin loop
	    }//pp collisions type
	  
	  }else if(matchReq){
	    if(debug)cout << "MATCH IS REQUIRED!" << endl;
	    //pass if MATCHING is REQUIRED required
	    int ind_recoJet = -1;
	    float delR_Smallest = 100000;
	    if(debug)cout << __LINE__ << endl;

	    if(debug)cout << "This is the number of reco jets that we have for this event: " << numJetsTot << endl;
	    
	    for(int iRecoJet = 0; iRecoJet < numJetsTot; iRecoJet++){
	      if(UniqueMatches[iRecoJet]){
		if(debug)cout << "Looks like this reco alsread found a truth match! " << endl;
		std::map<int,int>::iterator iter = TruthMatchIndicies.find(iRecoJet);
		int trth = iter->second;
		if(debug)cout << "This is the truth jet it was matched to (index/pt/eta/phi): " << trth << "/" << akt_truth_jet_kin[jet_Rad][pTPar]->at(trth) << "/" << akt_truth_jet_kin[jet_Rad][etaPar]->at(trth) << "/" << akt_truth_jet_kin[jet_Rad][phiPar]->at(trth) << endl;
		

		continue; //This reco jet already found a date.
	      
	      }
	      if(debug)cout << "This is the Reco Info! (Index/pT/eta/phi) " << iRecoJet << "/" << akt_jet_kin[jet_Rad][pTPar]->at(iRecoJet) << "/" << akt_jet_kin[jet_Rad][etaPar]->at(iRecoJet) << "/" << akt_jet_kin[jet_Rad][phiPar]->at(iRecoJet) << endl;
	      if((akt_jet_kin[jet_Rad][etaPar]->at(iRecoJet) > etacut) || (akt_jet_kin[jet_Rad][etaPar]->at(iRecoJet)< -etacut)) continue;
	      if(debug)cout << "This reco jet passed out eta cut!" <<endl;
	      if(pTRecoCut > akt_jet_kin[jet_Rad][pTPar]->at(iRecoJet) && collisionType=="PbPb")continue;	    
	      if(minpTCuts[jet_Rad] >akt_jet_kin[jet_Rad][pTPar]->at(iRecoJet) && collisionType == "pp")continue;

	      //Dhanush had an upper cut at 1TeV
	      if(dijet2018bins && (akt_jet_kin[jet_Rad][pTPar]->at(iRecoJet) > 1000))continue;

	      if(debug){
		cout <<"This reco jet PASEED PT CUT! " << endl;
		cout << "This is the truth eta: " << akt_truth_jet_kin[jet_Rad][etaPar]->at(iJetTruth) << endl;
		cout << "This is the reco eta: " <<  akt_jet_kin[jet_Rad][etaPar]->at(iRecoJet) << endl;
		cout << "This is the reco phi: " <<  akt_jet_kin[jet_Rad][phiPar]->at(iRecoJet) << endl;
		cout << "This is the truth phi: " <<  akt_truth_jet_kin[jet_Rad][phiPar]->at(iJetTruth) << endl;
	      }
	      float delta_R = deltaR_calc(akt_truth_jet_kin[jet_Rad][etaPar]->at(iJetTruth),
					  akt_jet_kin[jet_Rad][etaPar]->at(iRecoJet),
					  akt_truth_jet_kin[jet_Rad][phiPar]->at(iJetTruth),
					  akt_jet_kin[jet_Rad][phiPar]->at(iRecoJet));
	      
	      if(debug)cout << "this is the deltaR number: " << delta_R << endl;
	      
              if(delta_R > deltaR_MatchCut[jet_Rad])continue;
	      if(delta_R > delR_Smallest) continue;
	     
	      if(debug)cout << "saving indicies for now! but we think we found a match! " << endl;
              ind_recoJet = iRecoJet;
              delR_Smallest = delta_R;
	      if(debug)cout << __LINE__ << endl;
	      //after looping over ALL reco jets lets sacve indicies in map matched and not matched indicies!
	      
	    }//Reco Jet
	    
	    
	    if(debug)cout << __LINE__ << endl;

	    if(ind_recoJet==-1){
	      
	       continue;
	    }//no match found
	    
	    TruthMatchIndicies.insert({ind_recoJet,iJetTruth});

	    if(debug)cout << "saved these indicies! " << endl;
	    if(debug)cout << "truth index: " << iJetTruth << endl;
	    if(centBin==0 && jet_Rad == 1 && false){
	      cout << "R=0.4 truth MATCHED jet. This is the pT: "<< akt_truth_jet_kin[jet_Rad][pTPar]->at(iJetTruth) << endl;
	      cout << "R=0.4 reco MATCHED jet. This is the index/pT/eta/phi: " << ind_recoJet << "/" << akt_jet_kin[jet_Rad][pTPar]->at(ind_recoJet) << "/" << akt_jet_kin[jet_Rad][etaPar]->at(ind_recoJet) << "/" << akt_jet_kin[jet_Rad][phiPar]->at(ind_recoJet) << endl;
	      cout << endl;
	     
	      cout << "Printing all of the truth jets... which ones did not make it to the cut? " << endl;
	      cout << "THis is the total number of truth jets: " << akt_truth_jet_n[jet_Rad] << endl;
	     
	      for(int itruthJet =0; itruthJet < akt_truth_jet_n[jet_Rad]; itruthJet++){
		cout << "TRUTH index/pT/eta/phi: " << itruthJet << "/" << akt_truth_jet_kin[jet_Rad][pTPar]->at(itruthJet) << "/" << akt_truth_jet_kin[jet_Rad][etaPar]->at(itruthJet) << "/" << akt_truth_jet_kin[jet_Rad][phiPar]->at(itruthJet) <<  endl;
	      }
	      cout << endl;
	      cout << "Printing reco jet information... " << endl;
	      for(int irecoJet=0; irecoJet< numJetsTot; irecoJet++){
		cout << "RECO JET index/pT/eta/phi: " << irecoJet << "/" << akt_jet_kin[jet_Rad][pTPar]->at(irecoJet) << "/" << akt_jet_kin[jet_Rad][etaPar]->at(irecoJet) << "/" << akt_jet_kin[jet_Rad][phiPar]->at(irecoJet) << endl;
		for(int itruthJet=0; itruthJet < akt_truth_jet_n[jet_Rad]; itruthJet++){
		  cout << "this is the delta r between truth index: " << itruthJet << " and reco index " << irecoJet << endl;
		  float delta_R = deltaR_calc(akt_truth_jet_kin[jet_Rad][etaPar]->at(itruthJet),
					      akt_jet_kin[jet_Rad][etaPar]->at(irecoJet),
					      akt_truth_jet_kin[jet_Rad][phiPar]->at(itruthJet),
					      akt_jet_kin[jet_Rad][phiPar]->at(irecoJet));
		  cout << "delta R = " << delta_R << endl;
		}
		
	      }

	    }
	    if(debug)cout << "reco index: " << ind_recoJet << endl;
	    
	    if(debug)cout << __LINE__ << endl; 
	  
	    if(debug)cout << __LINE__ << endl;

	    
	  
	  // if(collisionType=="PbPb"){
	  //     if(debug)cout << "listing reco jets that were matched" << endl;
	      
	  //     if(debug)cout << __LINE__ << endl;
	  //     UniqueMatches[ind_recoJet] = true; //This reco jet cannot be matched twice.
	  //     //pTDisTruth[jet_Rad][centBin]->Fill(akt_truth_jet_kin[jet_Rad][pTPar]->at(iJetTruth),centWeight*jZweight);
	  //     if(debug)cout << __LINE__ << endl;
	  //     pTDis_RecoJetsMatched_PbPbMC[jet_Rad][centBin]->Fill(akt_jet_kin[jet_Rad][pTPar]->at(ind_recoJet),centWeight*jZweight);
	  // }else if(collisionType == "pp"){
	  //   for(int iCentBin =0; iCentBin < number_CentBins; iCentBin++){
	  //     pTDisTruth[jet_Rad][iCentBin]->Fill(akt_truth_jet_kin[jet_Rad][pTPar]->at(iJetTruth),jZweight);
	  //   }//centBin loop
	  // }//pp collissions
	    
	  if(debug)cout << __LINE__ << endl;
	    
	  
	  if(jz_samples_pTDis){
	    pTDis_JZSamples[jet_Rad][iJZBin][1]->Fill(akt_jet_kin[jet_Rad][pTPar]->at(ind_recoJet),jZweight);
	    pTDis_JZSamples[jet_Rad][iJZBin][0]->Fill(akt_truth_jet_kin[jet_Rad][pTPar]->at(iJetTruth),jZweight);
	   }

	  }//match requirement
	}//Truth jet loop.
	if(debug)cout << __LINE__ << endl;
      }//ONLY MC

      if(debug)cout << __LINE__ << endl;

      bool passOnlyOnce = false;
      
      
        
      
      if(debug)cout << __LINE__ << endl;
      bool UniqueMatchToTrthJet[200]={};
      
      //DEBUG
      if(debug & false){
	cout << "JET RADIUS YOU ARE LOOKING AT:" << jet_Rad << endl;
	cout << "This is how many reco jets you have: " << numJetsTot << endl;
	cout << "Printing all reco jets " << endl;
	
	
	for(int irec =0; irec < numJetsTot; irec++){
	  cout << "irec : " << irec << endl;

	  cout << "This is the pT of jet: "  << akt_jet_kin[jet_Rad][pTPar]->at(irec) << endl;
	  cout << "This is the eta of your jet: " << akt_jet_kin[jet_Rad][etaPar]->at(irec) << endl;
	  cout << "This is the phi of your jet: " << akt_jet_kin[jet_Rad][phiPar]->at(irec) << endl;

	  cout << "RECO (index/pT/eta/phi): " << irec << "/" << akt_jet_kin[jet_Rad][pTPar]->at(irec) << "/" << akt_jet_kin[jet_Rad][etaPar]->at(irec) <<  "/" << akt_jet_kin[jet_Rad][phiPar]->at(irec) <<endl;
	}
	cout << endl;
	cout << "This is how mamny truth jets you have: " << akt_truth_jet_n[jet_Rad] << endl;
	
      
	cout << "Printing all truth jets " << endl;
	
	for(int itru = 0; itru < akt_truth_jet_n[jet_Rad]; itru++){
	  cout << "TRUTH (Index/pT/eta/phi): " << akt_truth_jet_kin[jet_Rad][pTPar]->at(itru) << "/" << akt_truth_jet_kin[jet_Rad][etaPar]->at(itru) << "/" << akt_truth_jet_kin[jet_Rad][phiPar]->at(itru) << endl;
	}

	cout << "This is how many matched you have: " << TruthMatchIndicies.size() << endl;
	cout << "Printing reco jets that were matched..... " << endl;
	for(auto it = TruthMatchIndicies.cbegin(); it != TruthMatchIndicies.cend(); ++it)
	  {
	    cout << "RECO JET: INDEX/pT/eta/phi    -----      TRUTH JET INDEX: INDEX/pT/eta/phi    " << endl;
	    int truthin = it->second;
	    int recoind = it->first;
	  
	    std::cout << it->first << " / " << akt_jet_kin[jet_Rad][pTPar]->at(recoind) << " / " << akt_jet_kin[jet_Rad][etaPar]->at(recoind) << " / " << akt_jet_kin[jet_Rad][phiPar]->at(recoind)  << "   -------   " << it->second << " / " << akt_truth_jet_kin[jet_Rad][pTPar]->at(truthin) << " / " << akt_truth_jet_kin[jet_Rad][etaPar]->at(truthin) << " / " << akt_truth_jet_kin[jet_Rad][phiPar]->at(truthin) << "\n";
	  
	  }
      }//Debug

      // if(debug){
      // 	cout << "printing out matching pairs." << endl;
      // 	cout << "RECO INDEX, TRUTH INDEX" << endl;
      // 	for (auto const &pair: TruthMatchIndicies) {
      // 	  std::cout << "{" << pair.first << ": " << pair.second << "}\n";
      // 	}
      // }


      // cout << "This is the jet radius we are looking at: " << jet_Rad << endl;
      // cout << "This is the total number of jets we have: " << numJetsTot << endl;
      for(int iJetReco = 0; iJetReco < numJetsTot; iJetReco++){
	if(debug)cout << "This reco jet (Index/pT/eta/phi): " << iJetReco << "/"<< akt_jet_kin[jet_Rad][pTPar]->at(iJetReco) << " / " << akt_jet_kin[jet_Rad][etaPar]->at(iJetReco) << " / " << akt_jet_kin[jet_Rad][phiPar]->at(iJetReco) << endl;
	
	if(collisionType == "pp" && data_or_mc == mcOrdata::data && clean_ppJets){

	  if (looseBadR[jet_Rad] == nullptr) {
	    if(debug)std::cerr << "Null pointer encountered at index " << jet_Rad << std::endl;
	    continue; // continue to next iteration of loop
	  }
	  bool badjet = looseBadR[jet_Rad]->at(iJetReco);


	    if(!badjet){	      
	      if(fakeJets_Study && minpTCuts[jet_Rad] < akt_jet_kin[jet_Rad][pTPar]->at(iJetReco))pTDis_RJ_cleanTool[jet_Rad]->Fill(akt_jet_kin[jet_Rad][pTPar]->at(iJetReco));
	      continue;
	    }
	  }//Cleaning Tool for pp data ONLY!

	
	

	if(matchReq){
	  int truth_index = -1;
	  // cout << endl;
	  //cout << "This is the event of iEvent: " << iEvent << endl;
	  // cout << "This is the curent value of truth_index: " << truth_index << endl;
	  std::map<int,int>::iterator iter = TruthMatchIndicies.find(iJetReco);
	  
	  

	  if (iter == TruthMatchIndicies.end()) continue;
	  else truth_index = iter->second;
	  if(debug)cout << "we have grabbed the element: " << iter->second << "<--- this is suppose to be the matched pair to this reco jet!" << endl;
	  if(debug)cout << "if its < 0 then we must move on top next jet! " << endl;
	  if (truth_index < 0) continue;
	  double ptShpWeight = 1.0;
	  double ptShpWeight_2D =1.0;
	  if(addpTShapeWght){
	       
	       if(akt_truth_jet_kin[jet_Rad][pTPar]->at(truth_index) < lowEndCutOff){
		  ptShpWeight=ptshapeWghts[jet_Rad][centBin]->Eval(lowEndCutOff); 
	       }else if(akt_truth_jet_kin[jet_Rad][pTPar]->at(truth_index) > highEndCutOff){
		 ptShpWeight=ptshapeWghts[jet_Rad][centBin]->Eval(highEndCutOff);
	       }else{		
		  ptShpWeight=ptshapeWghts[jet_Rad][centBin]->Eval(akt_truth_jet_kin[jet_Rad][pTPar]->at(truth_index)); //instead of inputing reco we input truth
		  //if(centBin==0 && akt_truth_jet_kin[jet_Rad][pTPar]->at(truth_index) < 70)cout << "This is the value of ptShpWeight: " << ptShpWeight << endl;
	       }

	  }
	  
	  if(UniqueMatchToTrthJet[truth_index]==false){
	    //cout << "Guess not! This truth jet can no longer be matched to a reco jet." << endl;
	     

	    if(collisionType=="PbPb"){
	      //only Pb+Pb Collisions
	      
	     	if(centBin==0 && akt_truth_jet_kin[jet_Rad][pTPar]->at(truth_index) > 700  && jet_Rad==1 && false){
		cout << "This is the cent 0-10% bin. This pT truth made it to the fill histo: " << akt_truth_jet_kin[jet_Rad][pTPar]->at(truth_index) << endl;
		cout << "This is the eta and phi: " << akt_truth_jet_kin[jet_Rad][etaPar]->at(truth_index) << "/" << akt_truth_jet_kin[jet_Rad][phiPar]->at(truth_index) << endl;
	      
		cout << "TRUTH JET made it to the fill histo pT: " << akt_truth_jet_kin[jet_Rad][pTPar]->at(truth_index) << endl;
		cout << "Reco Jet made it to the fill histo pT: " << akt_jet_kin[jet_Rad][pTPar]->at(iJetReco) << endl;
		cout << "This is the pT shape weight: " << ptShpWeight << endl;
		}
		
	       
		  pTDisTruth[jet_Rad][centBin]->Fill(akt_truth_jet_kin[jet_Rad][pTPar]->at(truth_index),centWeight*jZweight*ptShpWeight);
		
		

	      if(centBin==0 && jet_Rad==1 && akt_truth_jet_kin[jet_Rad][pTPar]->at(truth_index) < 268){
		//cout << "This is the pT of your matched jet: " << akt_truth_jet_kin[jet_Rad][pTPar]->at(truth_index) << endl;
		// cout << "This is the pT of the reco jet matched: " << akt_jet_kin[jet_Rad][pTPar]->at(iJetReco) << endl;
		// cout << "this is the pt shape weight: " << ptShpWeight << endl;
		R4_cent0_10TruthMathedJ->Fill(akt_truth_jet_kin[jet_Rad][pTPar]->at(truth_index),jZweight);
		R4_cent0_10RecoMathedJ->Fill(akt_jet_kin[jet_Rad][pTPar]->at(iJetReco),jZweight);
		R4_cent0_10_pTShpWght->Fill(akt_jet_kin[jet_Rad][pTPar]->at(iJetReco),ptShpWeight,jZweight);
	      }
	     
	      pTDis[jet_Rad][centBin]->Fill(akt_jet_kin[jet_Rad][pTPar]->at(iJetReco),centWeight*jZweight*ptShpWeight);
	      if(jz_samples_pTDis){
		pTDis_JZSamples[jet_Rad][iJZBin][1]->Fill(akt_jet_kin[jet_Rad][pTPar]->at(iJetReco),jZweight);
		pTDis_JZSamples[jet_Rad][iJZBin][0]->Fill(akt_truth_jet_kin[jet_Rad][pTPar]->at(truth_index),jZweight);
	      }
	      if(etaPhiplotsOnly){
		//etaphiPlots[jet_Rad][centBin]->Fill(akt_jet_kin[jet_Rad][etaPar]->at(iJetReco),akt_jet_kin[jet_Rad][phiPar]->at(iJetReco),centWeight*jZweight);
		if(debug)cout << __LINE__ << endl;
		for(int iParam =1; iParam < totJetParam; iParam++){
		  if(debug)cout <<__LINE__ << endl;
		  cout << "*******[jet_Rad][iParam][centBin][truth]: " << jet_Rad << "/ " << iParam << "/" << centBin << "/" << truth << endl;
		  eta_phi_dis[jet_Rad][iParam][centBin][truth]->Fill(akt_truth_jet_kin[jet_Rad][iParam]->at(truth_index),centWeight*jZweight);
		  
		  eta_phi_dis[jet_Rad][iParam][centBin][recon]->Fill(akt_jet_kin[jet_Rad][iParam]->at(iJetReco),centWeight*jZweight);
		  if(debug)cout <<__LINE__ << endl;
		}
		if(debug)cout <<__LINE__ << endl;
	      }
	      if(debug)cout <<__LINE__ << endl;
	      if(akt_truth_jet_kin[jet_Rad][pTPar]->at(truth_index) > 100 && akt_truth_jet_kin[jet_Rad][pTPar]->at(truth_index) < 200 && centBin == 0){
		if(debug)cout <<__LINE__ << endl;
		pTRatio[jet_Rad]->Fill(akt_jet_kin[jet_Rad][pTPar]->at(iJetReco)/akt_truth_jet_kin[jet_Rad][pTPar]->at(truth_index),centWeight*jZweight);
		if(debug)cout <<__LINE__ << endl;
	      }
	      if(debug)cout <<__LINE__ << endl;
	      if(full_and_ClsrTst){
		if(debug)cout <<__LINE__ << endl;
		responseMatrix[jet_Rad][centBin]->Fill(akt_jet_kin[jet_Rad][pTPar]->at(iJetReco),akt_truth_jet_kin[jet_Rad][pTPar]->at(truth_index),centWeight*jZweight*ptShpWeight);//Resp. Mtx. Full Closure Test
		
		if(corrHalf==0)halfResponseMtx[jet_Rad][centBin]->Fill(akt_jet_kin[jet_Rad][pTPar]->at(iJetReco),akt_truth_jet_kin[jet_Rad][pTPar]->at(truth_index),centWeight*jZweight*ptShpWeight);
		
		pTDisHalvesTruth[corrHalf][jet_Rad][centBin]->Fill(akt_truth_jet_kin[jet_Rad][pTPar]->at(truth_index),centWeight*jZweight*ptShpWeight);
		pTDisHalvesReco[corrHalf][jet_Rad][centBin]->Fill(akt_jet_kin[jet_Rad][pTPar]->at(iJetReco),centWeight*jZweight*ptShpWeight);
	
	      }//Only for Full and Half closure test
	    }else if(collisionType == "pp"){
	      if(debug)cout << "This is the TRUTH JET that was used to fill histogram " << endl;
	      if(debug)cout << "Index/pT/eta/phi: " << truth_index << "/" <<akt_truth_jet_kin[jet_Rad][pTPar]->at(truth_index) << "/" << akt_truth_jet_kin[jet_Rad][etaPar]->at(truth_index) << "/" << akt_truth_jet_kin[jet_Rad][phiPar]->at(truth_index) << endl;
	      
	      // cout << "THIS IS EVENT: " << iEvent << endl;
	      // 	cout << "This truth jet was matched (index/pt/eta/phi): " << truth_index  << "/" << akt_truth_jet_kin[jet_Rad][pTPar]->at(truth_index) << "/" << akt_truth_jet_kin[jet_Rad][etaPar]->at(truth_index) << "/" << akt_truth_jet_kin[jet_Rad][phiPar]->at(truth_index) << endl;
	      // 	cout << "This reco jet was matched (index/pt/eta/phi/sumpttrack): " << iJetReco << "/" << akt_jet_kin[jet_Rad][pTPar]->at(iJetReco) << "/" << akt_jet_kin[jet_Rad][etaPar]->at(iJetReco) << "/" << akt_jet_kin[jet_Rad][phiPar]->at(iJetReco) << "/" << akt_jet_sumTightPrimaryTrackpT[jet_Rad]->at(iJetReco) <<endl;

		// float delta_R = deltaR_calc(akt_truth_jet_kin[jet_Rad][etaPar]->at(truth_index),
		// 			    akt_jet_kin[jet_Rad][etaPar]->at(iJetReco),
		// 			    akt_truth_jet_kin[jet_Rad][phiPar]->at(truth_index),
		// 			    akt_jet_kin[jet_Rad][phiPar]->at(iJetReco));
		// cout << "This is the delta R between the matched pair: "<< delta_R << endl;
		// float deltaRNum = deltaR_calc(-0.729472,
                //                             -0.718895,
                //                             -0.954941,
		// 			      -1.22764);
		// cout << "dR(44Truth, 105Reco): " << deltaRNum << endl;
	      for(int iCentBin =0; iCentBin < number_CentBins; iCentBin++){
		if(pTCutsReco[jet_Rad][iCentBin] > akt_jet_kin[jet_Rad][pTPar]->at(iJetReco))continue;
		if(dijet2018bins && (akt_jet_kin[jet_Rad][pTPar]->at(iJetReco)>1000))continue;

	        if(pTCutsTruth[jet_Rad][iCentBin]> akt_truth_jet_kin[jet_Rad][pTPar]->at(truth_index))continue;
		
	     
		if(!sumpTCut_MtchTrthJets){
		  pTDisTruth[jet_Rad][iCentBin]->Fill(akt_truth_jet_kin[jet_Rad][pTPar]->at(truth_index),jZweight*ptShpWeight);
		}

		

		
	      }//centBin loop
	    }//only pp collisions
	    UniqueMatchToTrthJet[truth_index]=true;
	  }//Unique matches only
	}//Matching Req.


	
	if(debug)cout << __LINE__ << endl;
	
	
        if(akt_jet_kin[jet_Rad][etaPar]->at(iJetReco) > etacut || akt_jet_kin[jet_Rad][etaPar]->at(iJetReco) < -etacut){continue;}
        
        // if(collisionType == "PbPb" && data_or_mc == mc){
        //   if(debug)cout<< __LINE__ << endl;
        //   if((jet_Rad==R4)&&(akt_jet_kin[jet_Rad][pTPar]->at(iJetReco) > 400) && (iJZBin==0))continue;
        //   if((jet_Rad==R4)&&(akt_jet_kin[jet_Rad][pTPar]->at(iJetReco) > 900) && (iJZBin==2))continue;
        //   if(debug)cout << __LINE__ << endl;
        //   if((jet_Rad==R10)&&(akt_jet_kin[jet_Rad][pTPar]->at(iJetReco) > 750) && (iJZBin==1))continue;
        //   if((jet_Rad==R10)&&(akt_jet_kin[jet_Rad][pTPar]->at(iJetReco) > 500) && (iJZBin==0))continue;
        //   if((jet_Rad==R10)&&(akt_jet_kin[jet_Rad][pTPar]->at(iJetReco) > 1200) && (iJZBin==2))continue;
        // }


	if(debug)cout << "This  is the iJetReco value: " << iJetReco << endl;
        if((jet_Rad==R10) && (leading_pT_HLTa10< akt_jet_kin[jet_Rad][pTPar]->at(iJetReco)) && trigEff){
          leading_pT_HLTa10 = akt_jet_kin[jet_Rad][pTPar]->at(iJetReco); //HLT_j150_a10* fired
          if(debug)cout << "should not go through here!!" << endl;
          if(debug)cout << __LINE__ << endl;
        }//Picking Leading jet
	if(debug)cout << __LINE__ << endl;
	if(debug)cout << "THis is the pT of your jet: " << akt_jet_kin[jet_Rad][pTPar]->at(iJetReco) << endl;
	if(debug)cout << "It must pass this pT cut: " << pTCutsReco[jet_Rad][centBin] << endl;
        if(collisionType == "PbPb" && (akt_jet_kin[jet_Rad][pTPar]->at(iJetReco) < pTCutsReco[jet_Rad][centBin])) continue; //pT cut is require for reco jet!w
	
	if(collisionType == "pp" && mcOrdata::data == data_or_mc  && akt_jet_kin[jet_Rad][pTPar]->at(iJetReco) > 800 && etaPhiPlots){
	  etaphiNoisyJets[jet_Rad]->Fill(akt_jet_kin[jet_Rad][phiPar]->at(iJetReco),akt_jet_kin[jet_Rad][etaPar]->at(iJetReco));
	}

        
	  
	  

          
	if(debug)cout << "*******THIS IS THE THE pT of your jet: " << akt_jet_kin[jet_Rad][pTPar]->at(iJetReco) << endl;
	if(debug)cout << "this is the minimum pT cut: " << minpTCuts[jet_Rad] << endl;
	if(minpTCuts[jet_Rad] > akt_jet_kin[jet_Rad][pTPar]->at(iJetReco))continue;
	
	  
	if(debug)cout << __LINE__ << endl;
	if(debug)cout << "This  is the iJetReco value: " << iJetReco << endl;
        if(collisionType == "pp"){
	  //I would only go through here if you wanted to match jets!
          //make sure to require matchreq to be true 
	  if(debug)cout << "It must meet the pT cut: " << sumpTCuts[jet_Rad] << endl;

	  if(sumpTTrackDis){
	    
	    for(int ipTRange=0; ipTRange < 2; ipTRange++){
	      if(akt_jet_kin[jet_Rad][pTPar]->at(iJetReco) > pTRange_SumpTDis[ipTRange] && akt_jet_kin[jet_Rad][pTPar]->at(iJetReco) < pTRange_SumpTDis[ipTRange+1]){
		sumpTTrack_Dis[jet_Rad][ipTRange]->Fill(akt_jet_sumTightPrimaryTrackpT[jet_Rad]->at(iJetReco),jZweight);
	     }
	    }//pTRange Loop
	  }

	  

	  if(akt_jet_sumTightPrimaryTrackpT[jet_Rad]->at(iJetReco) < sumpTCuts[jet_Rad] && sumpTReq){
	    //Collect Bad Jet Indicies
	    if(collectBadJetIndicies && jet_Rad==R10){
	      // cout << "This is entry: " << iEvent << endl;
	      // cout << "Saving this bad jet index: " << iJetReco << endl;
	      BadJets->push_back(iJetReco);
	      BadJets_n++;
	    }
	   if(fakeJets_Study)pTDis_RJ_sumpTCut[jet_Rad]->Fill(akt_jet_kin[jet_Rad][pTPar]->at(iJetReco));
	     continue; //applying sumpTcut
	   }
	
	  if(debug)cout << __LINE__ << endl;
	  

	  if(etaPhiPlots && akt_jet_kin[jet_Rad][pTPar]->at(iJetReco) > 800 && mcOrdata::data == data_or_mc){
	    etaphiSumpTCut[jet_Rad]->Fill(akt_jet_kin[jet_Rad][phiPar]->at(iJetReco),akt_jet_kin[jet_Rad][etaPar]->at(iJetReco));
	  }
	  if(debug)cout << "****STILL HERE******" << endl;

	
	  //This is to see the MultiJet statistics
	  if(debug)cout << __LINE__ << endl;
	
	  if(checkMultiJetStats==true && collisionType == "pp" && akt_jet_kin[jet_Rad][pTPar]->at(iJetReco) > 400 && jet_Rad==R10){


	    // if(akt_jet_kin[jet_Rad][etaPar]->at(iJetReco)<-etacut || etacut < akt_jet_kin[jet_Rad][etaPar]->at(iJetReco))continue;
	    if(debug)cout << __LINE__ << endl;
	    if(debug)cout << "This jet is above 400 gev! and passed eta cut!" << endl;
	    int iCountJets_pT200GeV = 0;
	    subLeadingJets.clear();
	    if(debug)cout << "This is the entry: " << iEvent << endl;
	    int indexLeading_jet = -1;
	    double pt_leading =-1.0;
	    if(debug)cout << "is this the leading jet?" << endl;
	    if(debug)cout << "This is how many jets you have: " << numJetsTot << endl;
	    for(int ijet =0; ijet< numJetsTot; ijet++){
	      if(debug)cout << "index/pt/eta/phi/sumpt/cleanjet: " <<  ijet << "/"  << akt_jet_kin[jet_Rad][pTPar]->at(ijet) << "/" << akt_jet_kin[jet_Rad][etaPar]->at(ijet) << "/" << akt_jet_kin[jet_Rad][phiPar]->at(ijet) << "/" << akt_jet_sumTightPrimaryTrackpT[jet_Rad]->at(ijet) << "/" << looseBadR[jet_Rad]->at(ijet)  << endl;
	      double temp_jet_pt = akt_jet_kin[jet_Rad][pTPar]->at(ijet);
	      if(temp_jet_pt >pt_leading ){
		pt_leading  = temp_jet_pt;
		indexLeading_jet = ijet;
	      } 
	    }
	    
	    if(indexLeading_jet != iJetReco){
	      if(debug)cout << "guess its not!" << endl;
	      continue;
	    }

	    if(debug)cout << "This is how many jets you have in this event: " << numJetsTot << endl;
	    
	    if(numJetsTot==2 || numJetsTot==3){
	      pT_R10_400GeV_Only2_or3jets->Fill(akt_jet_kin[jet_Rad][pTPar]->at(iJetReco));
	    } 


	    if(numJetsTot != 3){
	      continue;
	    } 
	    
	    if(debug)cout << "Continue on! we fill histogram because this is the leading jet!" << endl;
	    if(numJetsTot==3)pT_R10_400GeV->Fill(akt_jet_kin[jet_Rad][pTPar]->at(iJetReco));

	    if(debug)cout << "This even has more than or equal to three jets" << endl;
	    if(debug)cout << "This is the number of R=1.0 jets in this event: " << numJetsTot << endl;
	    for(int ijetR10 =0; ijetR10 < numJetsTot; ijetR10++){
	      if(debug)cout << "This is the index we are currently looking at: " << ijetR10 << endl;
	      
	      
	      if(ijetR10 == iJetReco)continue;
	      if(debug)cout << "This is a different one than the one that made our 400 GeV cut..." << endl;
	      

	      if(akt_jet_kin[jet_Rad][etaPar]->at(ijetR10)<-etacut || etacut < akt_jet_kin[jet_Rad][etaPar]->at(ijetR10))continue;
	      
	      if(debug)cout << "This is the sumpT cut it had to pass: " << sumpTCuts[jet_Rad] << endl;
	      if(akt_jet_sumTightPrimaryTrackpT[jet_Rad]->at(ijetR10) < sumpTCuts[jet_Rad])continue;
	      
	      if(!looseBadR[jet_Rad]->at(ijetR10))continue;
	      if(debug)cout << "This good jet passed our eta cut, sunmpT, and clean cut! " << endl;
	      if(akt_jet_kin[jet_Rad][pTPar]->at(ijetR10) < 200)continue;
	      if(debug)cout << "This good jet passed the pT cut! (above 200 GeV!)" << endl;
	      iCountJets_pT200GeV++;
	      subLeadingJets.push_back(ijetR10);
	      if(debug)cout << "currently we have: " << iCountJets_pT200GeV << " jets that passed our 50 GeV cut" << endl;
		
	      if(debug)cout << "GOOD JET pT/eta/phi: " << akt_jet_kin[jet_Rad][pTPar]->at(ijetR10) << "/" << akt_jet_kin[jet_Rad][etaPar]->at(ijetR10) << "/" << akt_jet_kin[jet_Rad][phiPar]->at(ijetR10) << endl;
	      
	      
	    }//Loop over R=1.0 jets
	
	    if(iCountJets_pT200GeV==2){
	      if(debug)cout << "Saving this pT: " << akt_jet_kin[jet_Rad][pTPar]->at(iJetReco) << " and its index " << iJetReco << endl;
	      pT_R10_400GeV_2pT200->Fill(akt_jet_kin[jet_Rad][pTPar]->at(iJetReco));
	      if(debug)cout << __LINE__ << endl;

	      int subjetIndex_1 = subLeadingJets.at(0);
	      
	      int subjetIndex_2 = subLeadingJets.at(1);
	      
	      TLorentzVector subjet_1; TLorentzVector subjet_2;
	      TLorentzVector vectsum; 
	      
	      subjet_1.SetPtEtaPhiE(akt_jet_kin[jet_Rad][pTPar]->at(subjetIndex_1),akt_jet_kin[jet_Rad][etaPar]->at(subjetIndex_1),akt_jet_kin[jet_Rad][phiPar]->at(subjetIndex_1),jetEnergy->at(subjetIndex_1));
	      
	      subjet_2.SetPtEtaPhiE(akt_jet_kin[jet_Rad][pTPar]->at(subjetIndex_2),akt_jet_kin[jet_Rad][etaPar]->at(subjetIndex_2),akt_jet_kin[jet_Rad][phiPar]->at(subjetIndex_2),jetEnergy->at(subjetIndex_2));
	      
	      vectsum = subjet_1 + subjet_2;
	      
 	      for(int ipTRange =0; ipTRange < totLeadingpTRanges; ipTRange++){
		if(leadingpTRanges[ipTRange] > akt_jet_kin[jet_Rad][pTPar]->at(iJetReco) || akt_jet_kin[jet_Rad][pTPar]->at(iJetReco) > leadingpTRanges[ipTRange+1])continue;
		
		  GaussianDis[ipTRange]->Fill(vectsum.Pt()/akt_jet_kin[jet_Rad][pTPar]->at(iJetReco));
		  AJObservable[ipTRange]->Fill((akt_jet_kin[jet_Rad][pTPar]->at(subjetIndex_1)- akt_jet_kin[jet_Rad][pTPar]->at(subjetIndex_2))/akt_jet_kin[jet_Rad][pTPar]->at(iJetReco));
	      }
	      
	      

	    }
	  }
	     
	  
	  bool paseUno = false;
	  for(int iCentBin =0; iCentBin < number_CentBins; iCentBin++){
            if(debug)cout << __LINE__ << endl;
	    if(debug)cout << "This is the pT of the reco jet: " << akt_jet_kin[jet_Rad][pTPar]->at(iJetReco) << endl;
	    if(akt_jet_kin[jet_Rad][pTPar]->at(iJetReco) < pTCutsReco[jet_Rad][iCentBin])continue;

	    //Dhanush had an upper cut at 1TeV
	    if(dijet2018bins && (akt_jet_kin[jet_Rad][pTPar]->at(iJetReco) > 1000))continue;

	    if(debug)cout << "This is the pT cut that the jet had to pass: " << pTCutsReco[jet_Rad][iCentBin] << endl;
	    if(debug)cout << "This reco jet passed the minumum pT cut!" << endl; 
	    if(debug)cout << __LINE__ << endl;
	     if(data_or_mc == mcOrdata::mc){
	      std::map<int,int>::iterator iter = TruthMatchIndicies.find(iJetReco);
	      if(debug)cout << "Here is the matcheched truth jet index... " << endl;
	      int truthmatch=iter->second; 
	      if(debug)cout << "truth jet index: " << truthmatch << endl;
	       if(debug)cout << __LINE__ << endl;
	       if(debug)cout << "This is the pT of our truth jet: " << akt_truth_jet_kin[jet_Rad][pTPar]->at(truthmatch) << endl;
	       if(debug)cout << "This is the pT cut that the truth jet must pass: " << pTCutsTruth[jet_Rad][iCentBin] << endl;
	      if(akt_truth_jet_kin[jet_Rad][pTPar]->at(truthmatch) < pTCutsTruth[jet_Rad][iCentBin])continue;
	      if(debug)cout << "It passed! " << endl;
	       if(debug)cout << __LINE__ << endl;
	      
	       if(sumpTCut_MtchTrthJets){
		if(debug)cout << __LINE__ << endl;
		if(jet_Rad==R4 && iCentBin == 7 && akt_jet_kin[jet_Rad][pTPar]->at(iJetReco) > 120 && akt_jet_kin[jet_Rad][pTPar]->at(iJetReco) < 130 && iEvent == 181137 && jet_Rad == R4){
		  cout << "EVENT: " << iEvent << endl;
		  cout << "TRUTH index/pT/eta/phi: " << truthmatch << "/" << akt_truth_jet_kin[jet_Rad][pTPar]->at(truthmatch) << "/" << akt_truth_jet_kin[jet_Rad][etaPar]->at(truthmatch) << "/" << akt_truth_jet_kin[jet_Rad][phiPar]->at(truthmatch) << endl; 
		  cout << "RECO index/pT/eta/phi: " << iJetReco << "/" << akt_jet_kin[jet_Rad][pTPar]->at(iJetReco) << "/" << akt_jet_kin[jet_Rad][etaPar]->at(iJetReco) << "/" << akt_jet_kin[jet_Rad][phiPar]->at(iJetReco) << endl;
		  float delta_R_num = deltaR_calc(akt_truth_jet_kin[jet_Rad][etaPar]->at(truthmatch),
					      akt_jet_kin[jet_Rad][etaPar]->at(iJetReco),
					      akt_truth_jet_kin[jet_Rad][phiPar]->at(truthmatch),
					      akt_jet_kin[jet_Rad][phiPar]->at(iJetReco));

		  cout << "DISTANCE BETWEEN MATCHING PAIR: " << delta_R_num << endl;

		  cout << "PRINTING ALL OF RECO JETS... " << endl;

		  for(int irecojet =0 ; irecojet < numJetsTot; irecojet++){
		    cout << "RECO JET index/pT/eta/phi/sumpT: " << irecojet << "/" << akt_jet_kin[jet_Rad][pTPar]->at(irecojet) << "/" << akt_jet_kin[jet_Rad][etaPar]->at(irecojet) << "/" << akt_jet_kin[jet_Rad][phiPar]->at(irecojet) <<  endl;
		    //if(irecojet==0){
		      for(int truthjet =0; truthjet < akt_truth_jet_n[jet_Rad]; truthjet++){
			cout << "THIS IS THE DISTANCE BETWEEN RECO JET:  " << irecojet << " AND TRUTH JET: " << truthjet << endl;
			cout << "DISTANCE DELTAR: ";

			float delta_R_num = deltaR_calc(akt_truth_jet_kin[jet_Rad][etaPar]->at(truthjet),
							akt_jet_kin[jet_Rad][etaPar]->at(irecojet),
							akt_truth_jet_kin[jet_Rad][phiPar]->at(truthjet),
							akt_jet_kin[jet_Rad][phiPar]->at(irecojet));
			cout << delta_R_num << endl;

		      }//truth jet loop
		      //}

		  }//reco jet
		  
		  for(int itruthjet =0; itruthjet <akt_truth_jet_n[jet_Rad]; itruthjet++){
		    cout << "TRUTH JET index/pT/eta/phi: " << itruthjet << akt_truth_jet_kin[jet_Rad][pTPar]->at(itruthjet) << "/" << akt_truth_jet_kin[jet_Rad][etaPar]->at(itruthjet) << "/" << akt_truth_jet_kin[jet_Rad][phiPar]->at(itruthjet) << endl; 
		  }//truth jets


		}
		double ptShpWeight = 1.0;
		double ptShpWeight_2D = 1.0;
		if(addpTShapeWght){
		  if(akt_truth_jet_kin[jet_Rad][pTPar]->at(truthmatch)< lowEndCutOff){
		    ptShpWeight=ptshapeWghts[jet_Rad][centBin]->Eval(lowEndCutOff+1);
		  }else if(akt_truth_jet_kin[jet_Rad][pTPar]->at(truthmatch) > highEndCutOff){ 
		    ptShpWeight=ptshapeWghts[jet_Rad][centBin]->Eval(lowEndCutOff-1);
		  }else{		
		    ptShpWeight=ptshapeWghts[jet_Rad][centBin]->Eval(akt_truth_jet_kin[jet_Rad][pTPar]->at(truthmatch)); //instead of inputing reco we input truth
	
		  }		  
		
		}
	       
		pTDisTruth[jet_Rad][iCentBin]->Fill(akt_truth_jet_kin[jet_Rad][pTPar]->at(truthmatch),jZweight*ptShpWeight);
		
		if(full_and_ClsrTst)responseMatrix[jet_Rad][iCentBin]->Fill(akt_jet_kin[jet_Rad][pTPar]->at(iJetReco),akt_truth_jet_kin[jet_Rad][pTPar]->at(truthmatch),jZweight*ptShpWeight);//Resp. Mtx. Full Closure Test
		
		pTDis[jet_Rad][iCentBin]->Fill(akt_jet_kin[jet_Rad][pTPar]->at(iJetReco),jZweight*ptShpWeight);
		
		

		if(etaPhiplotsOnly){
		  //etaphiPlots[jet_Rad][iCentBin]->Fill(akt_jet_kin[jet_Rad][etaPar]->at(iJetReco),akt_jet_kin[jet_Rad][phiPar]->at(iJetReco),centWeight*jZweight);
		  for(int iParam =1; iParam < totJetParam; iParam++){
		    if(debug)cout << "[jet_Rad][iParam-1][iCentBin][truth]: " << jet_Rad << "/" << iParam-1 << "/" << iCentBin << "/ " << truth << endl;
		    eta_phi_dis[jet_Rad][iParam-1][iCentBin][truth]->Fill(akt_truth_jet_kin[jet_Rad][iParam]->at(truthmatch),centWeight*jZweight);
		    eta_phi_dis[jet_Rad][iParam-1][iCentBin][recon]->Fill(akt_jet_kin[jet_Rad][iParam]->at(iJetReco),centWeight*jZweight);
		  }
		}


		if(debug)cout << "Filling histogram with this jet pt reco: " << akt_jet_kin[jet_Rad][pTPar]->at(iJetReco) << endl;  
		 
		if(corrHalf==0 && full_and_ClsrTst){
		
		  halfResponseMtx[jet_Rad][iCentBin]->Fill(akt_jet_kin[jet_Rad][pTPar]->at(iJetReco),akt_truth_jet_kin[jet_Rad][pTPar]->at(truthmatch),jZweight*ptShpWeight);
		
		}
		  
		if(full_and_ClsrTst){
		  pTDisHalvesTruth[corrHalf][jet_Rad][iCentBin]->Fill(akt_truth_jet_kin[jet_Rad][pTPar]->at(truthmatch),jZweight*ptShpWeight);
		  pTDisHalvesReco[corrHalf][jet_Rad][iCentBin]->Fill(akt_jet_kin[jet_Rad][pTPar]->at(iJetReco),jZweight*ptShpWeight);
		}
	       }
		
	     }else if(data_or_mc==mcOrdata::data){
	       //cout << "This is the jet pT that we will use to fill our data histo: " << akt_jet_kin[jet_Rad][pTPar]->at(iJetReco) << endl;
	       pTDis[jet_Rad][iCentBin]->Fill(akt_jet_kin[jet_Rad][pTPar]->at(iJetReco));
	       if(etaPhiplotsOnly){
		 //etaphiPlots[jet_Rad][iCentBin]->Fill(akt_jet_kin[jet_Rad][etaPar]->at(iJetReco),akt_jet_kin[jet_Rad][phiPar]->at(iJetReco),centWeight*jZweight);
		 for(int iParam =1; iParam < totJetParam; iParam++){
		   eta_phi_dis[jet_Rad][iParam-1][iCentBin][recon]->Fill(akt_jet_kin[jet_Rad][iParam]->at(iJetReco),centWeight*jZweight);
		  
		 }
		}

	       	 
	       
	       if(debug)cout << __LINE__ << endl;
	      
	     }

          }//Cent Bin Loop
	  
	}else{
	  if(debug)cout << __LINE__ << endl;
          if(debug)cout << "This  is the iJetReco value: " << iJetReco << endl;
	  if(akt_jet_kin[jet_Rad][pTPar]->at(iJetReco) < pTCutsReco[jet_Rad][centBin])continue;
	  if(dijet2018bins && (akt_jet_kin[jet_Rad][pTPar]->at(iJetReco)>1000))continue;

	  if(debug)cout << __LINE__ << endl;
	  if(debug)cout << "This is event number" << iEvent << endl;
	  if(debug)cout << "This is the pT reco that was matched(Radius/Index/pT/eta/phi): " << jet_Rad<< "/" <<iJetReco << "/" << akt_jet_kin[jet_Rad][pTPar]->at(iJetReco) << "/" << akt_jet_kin[jet_Rad][etaPar]->at(iJetReco) << "/" << akt_jet_kin[jet_Rad][phiPar]->at(iJetReco) << endl;
	  
	  if(data_or_mc == mcOrdata::data){
	    if(debug)cout << "This is the pT of the jet that made it to our histo: " << akt_jet_kin[jet_Rad][pTPar]->at(iJetReco) << endl; 
	    pTDis[jet_Rad][centBin]->Fill(akt_jet_kin[jet_Rad][pTPar]->at(iJetReco),centWeight*jZweight);
	    if(etaPhiplotsOnly){
	      //etaphiPlots[jet_Rad][centBin]->Fill(akt_jet_kin[jet_Rad][etaPar]->at(iJetReco),akt_jet_kin[jet_Rad][phiPar]->at(iJetReco),centWeight*jZweight);
	      for(int iParam =1; iParam < totJetParam; iParam++){
		eta_phi_dis[jet_Rad][iParam-1][centBin][recon]->Fill(akt_jet_kin[jet_Rad][iParam]->at(iJetReco),centWeight*jZweight);
	      }
	    }

	  }
	  }//For Pb+Pb collisions


      }//Jet Reco Loop

      if(leading_pT_HLTa10 != -1 && jet_Rad==R10 && trigEff){
        if(centBin == 3 || centBin ==4){
          pTLeadJet_R10[0][cent_bin_30_50]->Fill(leading_pT_HLTa10); //HLT_j150* fired for R=1.0. Saving leading jet
        }else if(centBin == 5 || centBin == 6 || centBin == 7){
          pTLeadJet_R10[0][cent_bin_50_80]->Fill(leading_pT_HLTa10); //HLT_j150* fired for R=1.0. Saving leading jet
        }else{
          pTLeadJet_R10[0][centBin]->Fill(leading_pT_HLTa10); //HLT_j150* fired for R=1.0. Saving leading jet
        }
      }
    
    if(collectBadJetIndicies){  
      
      // cout << "We will fill this tree! " << endl;
      // cout << "This is the number of good jets: " << goodJets_n << endl;
      // cout << "Printing indicies" << endl;
      // for(int ijet =0; ijet < goodJets_n; ijet++){
      // 	cout << "this is index: " << goodJets.at(ijet) << endl;
      // }

      BadJetsTree->Fill();
      BadJets->clear();
      BadJets_n = 0;
      
    }

    }//Event Loop
  
  if(collectBadJetIndicies){
    cout << __LINE__ << endl;
    BadJetsTree->Write("",TObject::kOverwrite);
    BadJetsFile->Close();
    cout << __LINE__ <<endl;
  }

  string place_in_this_dir = "";
  if(jetRate2015Bins){
    place_in_this_dir = Form("/usatlas/u/bereniceg299/data/LargeRJet_Study/NewSourceCode//2015CentBinsFiles/%sVar/",collisionType.c_str());
  }else if(rAA2015Binning && jer_or_jes_Sys != ""){
    place_in_this_dir = Form("/usatlas/u/bereniceg299/data/LargeRJet_Study/NewSourceCode//2015CentBinsFiles/rAABins_%s/",collisionType.c_str());
  }else if(rAA2015Binning){
    cout << __LINE__ << endl;
    place_in_this_dir = "/usatlas/u/bereniceg299/data/LargeRJet_Study/NewSourceCode/2015CentBinsFiles/";
  }else if(dijet2018bins){
    place_in_this_dir = Form("/usatlas/u/bereniceg299/data/LargeRJet_Study/NewSourceCode/DijetBinning/%s/",collisionType.c_str());
  }else if(officialLargeRpTBins && etacut == 2.1){
    place_in_this_dir = Form("/usatlas/u/bereniceg299/data/LargeRJet_Study/NewSourceCode/largeR2018/%s/",collisionType.c_str());
  }else if(officialLargeRpTBins && etacut == 1.5){
    place_in_this_dir = Form("/usatlas/u/bereniceg299/data/LargeRJet_Study/NewSourceCode/NominalDir/R_%d/%s/",JetRadius[jet_Rad],collisionType.c_str());
  }

    cout << "Final name of the output file: " << Form("%sRawHistograms_R%d_%s_%s.root",place_in_this_dir.c_str(),JetRadius[jet_Rad],collisionType.c_str(),extraTag.c_str()) << endl;
  
  TFile *files_root = new TFile(Form("%sRawHistograms_R%d_%s_%s.root",place_in_this_dir.c_str(),JetRadius[jet_Rad],collisionType.c_str(),extraTag.c_str()),"RECREATE");
   files_root->cd();     

    if(trigEff){
      for(int iTrig =0; iTrig < tot_triggers; iTrig++){
        for(int iCentBin =0; iCentBin < number_CentBins; iCentBin++){
          pTLeadJet_R10[iTrig][iCentBin]->Write("",TObject::kOverwrite);
        }//Cent Bin Loop
      }//Trigger Loop
    }//Triggger Eff PLots
    
    
      if(sumpTTrackDis){
	for(int ipTRange =0; ipTRange < 2; ipTRange++){
	  sumpTTrack_Dis[jet_Rad][ipTRange]->Write("",TObject::kOverwrite);
	}//pTRange Loop
      }

      if(etaPhiplotsOnly){
	 for(int iCent =0; iCent < number_CentBins; iCent++){
	   //etaphiPlots[jet_Rad][iCent]->Write("",TObject::kOverwrite);
	   for(int iparam=0;iparam < totJetParam-1; iparam++){
	    for(int iType=0; iType < totType; iType++){
	      
	      if(debug)cout << "[jet_Rad][iparam][iCent][iType]: " << jet_Rad << "/" << iparam << "/" << iCent << "/" << iType << endl;
	      eta_phi_dis[jet_Rad][iparam][iCent][iType]->Write("",TObject::kOverwrite);
	      
	    }
	  }//cent loop
	}//parameter loop
      }



      if(jz_samples_pTDis){
	for(int iJZSamp =0; iJZSamp < total_samples; iJZSamp++){
	  pTDis_JZSamples[jet_Rad][iJZSamp][1]->Write("",TObject::kOverwrite);
	  pTDis_JZSamples[jet_Rad][iJZSamp][0]->Write("",TObject::kOverwrite);
	}
      }


      if(collisionType=="pp"){
	if(etaPhiPlots){
	  etaphiNoisyJets[jet_Rad]->Write("",TObject::kOverwrite);
	  etaphiSumpTCut[jet_Rad]->Write("",TObject::kOverwrite);
	}
	if(checkMultiJetStats==true){
	  pT_R10_400GeV_2pT200->Write("",TObject::kOverwrite);
	  pT_R10_400GeV->Write("",TObject::kOverwrite);
	  pT_R10_400GeV_Only2_or3jets->Write("",TObject::kOverwrite);
	  for(int ipTRange =0; ipTRange < totLeadingpTRanges; ipTRange++){
	    GaussianDis[ipTRange]->Write("",TObject::kOverwrite);
	    AJObservable[ipTRange]->Write("",TObject::kOverwrite);
	   }
	}

      }
      


      if(sumETHist && collisionType == "PbPb")sumETDisData[jet_Rad]->Write("",TObject::kOverwrite);
      
      if(fakeJets_Study && collisionType=="pp"){
	 pTDis_RJ_sumpTCut[jet_Rad]->Write("",TObject::kOverwrite);
	 pTDis_RJ_cleanTool[jet_Rad]->Write("",TObject::kOverwrite);
      }
      

      pTRatio[jet_Rad]->Write("",TObject::kOverwrite);
      
      for(int iCentBin = 0; iCentBin < number_CentBins; iCentBin++){
	
	pTDis[jet_Rad][iCentBin]->Write("",TObject::kOverwrite);
        if(data_or_mc == mc){
	  
	  if(jet_Rad==1 && iCentBin==0){
	    R4_cent0_10RecoMathedJ->Write("",TObject::kOverwrite);
	    R4_cent0_10TruthMathedJ->Write("",TObject::kOverwrite);
	    R4_cent0_10_pTShpWght->Write("",TObject::kOverwrite);
	  }
	  pTDisTruth[jet_Rad][iCentBin]->Write("",TObject::kOverwrite);
	  pTDisTruth_All[jet_Rad][iCentBin]->Write("",TObject::kOverwrite);
	  if(full_and_ClsrTst){
	    halfResponseMtx[jet_Rad][iCentBin]->Write("",TObject::kOverwrite);
	    responseMatrix[jet_Rad][iCentBin]->Write("",TObject::kOverwrite);
	    for(int iHalf =0; iHalf < twoHalves; iHalf++){
	      pTDisHalvesTruth[iHalf][jet_Rad][iCentBin]->Write("",TObject::kOverwrite);
	      pTDisHalvesReco[iHalf][jet_Rad][iCentBin]->Write("",TObject::kOverwrite);
	    }//looping over halves
	  }//Only for Full/Half clsr tests
	}//ONLY FOR MC

      }//Cent Bin Loop
    

    files_root->Close();


    return 0;
  }
