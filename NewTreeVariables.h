#include <assert.h>
#include <cmath>
#include <fstream>
#include <stdio.h>
#include <vector>
#include <math.h>
#include <tgmath.h>
#include <iostream>
#include <TChain.h>
#include <TLegend.h>
#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TRandom3.h"
#include <TAttMarker.h>
#include <TCanvas.h>
#include <TPaveText.h>
#include <dirent.h>
#include <sys/types.h>
#include "TLorentzVector.h"



 /**
  * Returns a linearly spaced array. The 0th element is lo, and the num-th element is hi.
  */
double* linspace (double lo, double hi, int num) {
  double* arr = new double[num+1];
  double delta = ((double)(hi)-(double)(lo))/(double)(num);
  for (int i = 0; i <= num; i++) {
    arr[i] = lo + i * delta;
  }
  return arr;
}


/**
 181  * Returns a logarithmically spaced array, where the 0th element is lo and the num-th element is hi.
 182  */
double* logspace (double lo, double hi, int num) {
  double loghi = TMath::Log2(hi);
  if (lo == 0) {
    double* arr = linspace(TMath::Log2(hi/(100*num)), loghi, num);
    for (int i = 0; i <= num; i++) {
      arr[i] = TMath::Power(2, arr[i]);
    }
    return arr;
  } else {
    double loglo = TMath::Log2(lo);
    double* arr = linspace(loglo, loghi, num);
    for (int i = 0; i <= num; i++) {
      arr[i] = TMath::Power(2, arr[i]);
    }
    return arr;
  }
}


//-----What Are Ytrue; ing To Run Over? 
const bool mc_pp =  false; ///true if you are running over pp mc
const bool data_pp =  true; //true if you are running over pp data
const bool data_PbPb = false; //true if you are running over Pb+Pb Data
const bool mc_PbPb = false; //true if you are running over Pb+Pb MC

//-----When debugging--------// 
bool debug_reco = false;
bool debug_trth = false;
bool debug = false;
bool debug_filled_histo = false;
bool debug_write = true;
bool debg_ppData = false;

bool debug_ppMC = false;
bool debug_PbPbMC = false;

//Unfolding Test
bool unfold_PbPb = false;

//----Tree
TChain *tree = new TChain("tree");
TChain *tree_pp_Data = new TChain("tree");
TChain *tree_pp_03Data = new TChain("tree");
TChain *tree_pp_MC = new TChain("tree");
TChain *tree_pythia = new TChain("tree");
TChain *tree_PbPbMC = new TChain("tree");
TChain *tree_PbPb_Data =  new TChain("tree");
//------Vectors for branches

/* //Binning for this analysis */
/* int pTBinsTot[] = {20,20,19,17,17,18,16,16}; */
/* const int largeRBins = 12; */
/* std::vector<std::vector<double>> binsnew; */
/* binsnew.reserve(8); */
/* binsnew.emplace_back(std::vector<double>{100, 117, 138, 158, 177, 197, 220, 241, 262, 294, 336, 384, 439, 502, 573, 656, 750, 857, 980, 1120, 1300}); //0-10% */
/* binsnewemplace_back(binsLargeR{100, 117, 138, 158, 177, 197, 220, 241, 262, 294, 336, 384, 439, 502, 573, 656, 750, 857, 980, 1120, 1300}); //10-20% */
/* binsnew.emplace_back(binsLargeR{100, 117, 138, 158, 177, 197, 220, 245, 262, 294, 336, 384, 439, 502, 573, 656, 750, 857, 980, 1120}); //20-30% */
/* binsnew.emplace_back(binsLargeR{100, 117, 138, 158, 177, 198, 230, 262, 294, 336, 384, 439, 502, 573, 656, 750, 857, 1000}); //30-40% */
/* binsnew.emplace_back(binsLargeR{100, 117, 138, 158, 177, 198, 230, 262, 294, 336, 384, 439, 502, 573, 656, 750, 857, 1000}); //40-50% */
/* binsnew.emplace_back(binsLargeR{100, 117, 138, 158, 177, 198, 220, 241, 262, 294, 336, 384, 439, 502, 573, 656, 750, 857, 1000}); //50-60% */
/* binsnew.emplace_back(binsLargeR{100, 117, 138, 158, 177, 197, 220, 241, 262, 294, 336, 384, 439, 502, 573, 656, 750, 1000}); //60-70% */
/* binsnew.emplace_back(binsLargeR{100, 117, 138, 158, 177, 197, 220, 241, 262, 294, 336, 384, 439, 502, 573, 656, 750, 1000}); //70-80% */

/* double pTBins_LargeR[8][13] = { */
/*                               {230, 262, 294, 336, 384, 439, 502, 573, 656, 750, 857, 980, 1120, 1300}, //0-10% */
/* 			      {200, 230, 262, 294, 336, 384, 439, 502, 573, 656, 750, 857, 980, 1120, 1300},//10-20% */
/* 			      {220, 230, 245, 262, 294, 336, 384, 439, 502, 573, 656, 750, 857, 980, 1120},//20-30% */
/* 			      {170, 200 ,230, 245, 262, 294, 336, 384, 439, 502, 573, 656, 750, 857, 1000},//30-40% */
/* 			      {180, 200, 230, 245, 262, 294, 336, 384, 439, 502, 573, 656, 750, 857, 1000}, //40-50% */
/* 			      {180, 200, 220, 245, 262, 294, 336, 384, 439, 502, 573, 656, 750, 857, 1000},//50-60% */
/* 			      {180, 200, 220, 245, 262, 294, 336, 384, 439, 502, 573, 656, 750, 1000}, //60-70% */
/* 			      {180, 200, 220, 245, 262, 294, 336, 384, 439, 502, 573, 656, 750, 1000} //70-80% */
/* }; */



//------------------------------------------------------------------------//
//-------------------------------pp Data----------------------------------//
//------------------------------------------------------------------------//
//--------Radii 
const int totRadii = 2;

//LooseBad
std::vector<bool> *isLooseBadR4 = NULL;

//---Parameers(pT,eta,phi)
const int totJetParam = 3;
//--------Tracks
const float pTtrackCut = 4; //GeV 

const int v_size = 100;
float vert_sumpt[v_size];
int vert_ntrk[v_size];
int vert_type[v_size];
int ntrk =0;

std::vector<float> *trk_phi = NULL;
std::vector<float> *trk_eta = NULL;
std::vector<float> *trk_pt = NULL;



//------Vertix, Triggers, Run Number 
int nvert = 0;
int nvert3 = 0;
bool HLT_j110_a10_lcw_subjes_L1J30_3 = 0;
bool HLT_j110_a10_lcw_subjes_L1J30 = 0;
bool HLT_j85 = false;
bool HLT_j110 = false;
int RunNumber = 0;

//Name of Triggers
const int totTrgs = 2; 
string triggerNames[]= {"HLT_j110_a10_lcw_subjes_L1J30","HLT_j85"};
bool *trigArray[totTrgs];
enum hlt_Trig {hlt_j110_a10,hlt_j85};

//When Initializing Histograms
float InitHistoParam[] = {0,1800,
                          -2.2,2.2,
			  -3.14,3.14};
//-----Jets 
int tot_ppDataJets[totRadii] ={};
std::vector<float> * ppDataJets_Info[totRadii][totJetParam] = {};
std::vector<float> *akt_jet_InsituCalib_Kin[totRadii][totRadii] = {};

//Lumi Numbers for HLT_j110* (R=1.0) and HLT_j85 (R=0.4)
float lumiValues[] = {256e+3,132.199e+3}; //nb^-1
int totTriggers = 2;
int wAndwoWeights = 2;

//-----Jet eta Cut
const float etaCutppData = 1.5; //Make sure to CHANGE this back to 1.5 for your final results!! 
//This is when we want to compare to the 2015 results
bool bins2015RAA = false;
int num2015Bins = 12;
float mea2015bins[] = {50,100,126,158,200,251,316,398,501,631,800,1000,1300};
//This is to compare 2018 Results from CMS Measurements
int cmsBins = 12;
float cmsLargeRBins[] = {50,70,100,150,200,250,300,350,400,500,650,1000,1300};
int cmsBinRespMtx = 12;
float cmsRespMtx[] = {50,70,100,150,200,250,300,350,400,500,620,1000,1300};

//-----Jet pT Cut
const float pTCut = 160; //GeV 
//-----Jet eta Cut
const float etaCut =1.5; //changed from 2.2 to 1.5 

//Jet pT Distribution Histograms
//THe first bin is 0-200, because I will cut off at 200 GeV. 




float bins_pTDis[26] = {};




float jetRaidus[] = {1.0, 0.4}; //DO NOT TOUCH THIS!!!! 
int jet_Radius[] = {1,4};
//-----Enum 
enum Radii {Rad10,Rad4,Rad3,Rad2};


//----Jets per Lumi Block plot
const int Triggers =2;
string nameTrig[] = {"a10","j85"};
const int totRuns = 13;
const int tot_jetRadii = 2; 
//int jetRadius[tot_jetRadii] ={10,4,3,2};

std::map <int,int> runNumber {{340644,1},{340683,2},{340697,3},
			      {340718,4},{340814,5},{340849,6},
			      {340850,7},{340910,8},{340925,9},
			      {340973,10},{341027,11},{341123,12},
			       {341184,13}};

string runNum[] = {"340644","340683","340697","340718","340814","340849","340850","340910","340925","340973","341027","341123","341184"};

int pTCutArray_Trigers[] = {200,160}; //Based on Triggers
float pTCuts[] = {294,100};
float pTRanges_NLumi[]={160,200,260};

//----JetpT Distributions
int sumpTWO_With =2; 
int sumpTCut = 10; //GeV
float sumpTCuts[] ={10,5}; //GeV 
enum SumpT_TF{WOSumpT, WSumpT}; //---with sumpT Cut and without it 
string SumpT[]={"WOSumpT","WSumpT"};

//----Eta Distributions & Eta-Phi Plots For Jets w/ SumpT > 10 GeV
const int ends_of_pTSlices = 5;
const int pTSlices_EtaPlots = 4;
float eta_pT_slices[ends_of_pTSlices] =  {400,600,800,1000,1200}; //-----4 different pT slices 
int sumpT_Trk_pTSlices[pTSlices_EtaPlots]={600,800,1000,1300};

//---sumpT Track distributions
const int  Tot_pTRanges_NoiseJets = 6;  
float  pTRanges_NoisyJets[] = {600,800,900,1100,1300,1500,1800};

//This is used for finner binnning to get the sumpTCut
float sumpTfinnerBins[]={0,1,2,3,4,5,10,15,20,30,40,50}; //11 bins
 

//------------------------------------------------------------------------//
//---------------------------------pp MC----------------------------------//
//------------------------------------------------------------------------//
//----Name of JZ Files 
const int tot_JZDirs = 8;
string names_JZFiles[tot_JZDirs] = {"user.berenice.02032022.LargeRJet_ppMC16_JZ2_no1_Analysis.0000000001_myOutput.root",
			   "user.berenice.02032022.LargeRJet_ppMC16_JZ2_no2_Analysis.0000000002_myOutput.root",
			   "user.berenice.02032022.LargeRJet_ppMC16_JZ3_no1_Analysis.0000000001_myOutput.root",
			   "user.berenice.02032022.LargeRJet_ppMC16_JZ3_no3_Analysis.0000000003_myOutput.root",
			   "user.berenice.02032022.LargeRJet_ppMC16_JZ4_no1_Analysis.0000000001_myOutput.root",
			   "user.berenice.02032022.LargeRJet_ppMC16_JZ4_no2_Analysis.0000000002_myOutput.root",
			   "user.berenice.02032022.LargeRJet_ppMC16_JZ5_no1_Analysis.0000000001_myOutput.root",
			   "user.berenice.02032022.LargeRJet_ppMC16_JZ5_no2_Analysis.0000000002_myOutput.root"};

//----Jets R's
enum Rad {R10,R4};
const int Tot_Radii = 2; ////Change it from 2 to 1
int JetRadius[] = {10,4};

//-----Vectors
//-----Jets
int aktR_jet_n_MCReco = 0;

int aktR_jet_n_MCTrth_4=0;  ///ADDED THIS NEW LINE OF CODE
int aktR_jet_n_MCTrth_10=0;
int aktR_jet_n_MCReco_4=0;
int aktR_jet_n_MCReco_10=0;

int aktR_jet_n_MCTrth = 0;

//------Jet Cuts
const int pTTruthCut = 40; //GeV
const int pTRecoCut = 20; //GeV

//------deltaR
const float deltaR_IsoCut = 2.1; //-----this is the delta R use to make sure my jet is isolated
std::vector<int> isolatedRecoJetsIndex;
std::vector<int> isolatedTruthJetsIndex;
std::vector<int> isolatedRecoJetsIndex_test;
std::vector<int> isolatedTruthJetsIndex_test;

float deltaR_MatchCut[] = {0.65,0.3}; //-----this is the deltaR used to match my truth jet to a reco jet  

//------Save the number of event for all the JZ samples
const int total_samples = 5; //changed it from 4 to 1 (testing purposes)
const int jz_entry_ranges = total_samples+1;
int num_JZ_range_entries[jz_entry_ranges] = {};

//-----pp MC Weights
float sigma_a_pp = (6.399*pow(10,5));
float eps_a_pp = (4.2785*pow(10,-3));
float num_events_a_pp = 7998000;

float sigma_d_pp = (4.7195*pow(10,3));
float eps_d_pp = (5.2994*pow(10,-3));
float num_events_d_pp = 7999000;

float sigma_e_pp = (2.6602*pow(10,1));
float eps_e_pp = (4.5901*pow(10,-3));
float num_events_e_pp = 8000000;

float sigma_f_pp = (2.2476*pow(10,-1));
float eps_f_pp = (2.1846*pow(10,-3));
float num_events_f_pp = 7999000;

float sigma_g_pp = 1.34*1e-3;
float eps_g_pp = 0.0008433528;
float num_events_g_pp = 19349;

//--Calculating Weights
const int tot_weights_ppMC = 5;
float weight_0_pp = (eps_a_pp*sigma_a_pp)/num_events_a_pp;
float weight_1_pp = (eps_d_pp*sigma_d_pp)/num_events_d_pp;
float weight_2_pp = (eps_e_pp*sigma_e_pp)/num_events_e_pp;
float weight_3_pp = (eps_f_pp*sigma_f_pp)/num_events_f_pp;
float weight_4_pp = (eps_g_pp*sigma_g_pp)/num_events_g_pp;

float weights_ppMC[tot_weights_ppMC] = {weight_0_pp,weight_1_pp,weight_2_pp,weight_3_pp,weight_4_pp};

//------Verticies and Tracks
int nvert_ppMC = 0;
float vert_sumpt_reco[v_size];
int vert_ntrk_reco[v_size];
int vert_type_reco[v_size];
int ntrk_reco = 0;
//----MC Reco
std::vector<float> *trk_phi_reco = NULL;
std::vector<float> *trk_eta_reco = NULL;
std::vector<float> *trk_pt_reco = NULL;

//SumpT Distributions
bool make_sumpTDis = false;


//------R=1.0
//-------Total Number of Reco Jets
int akt10_pp_jet_n_reco = 0;
int akt10_pp_jet_n_trth = 0;

//-----R=0.4
//--------Total Number of Truth Jets
int akt4_pp_jet_n_reco = 0;
int akt4_pp_jet_n_trth = 0;

//----Jet Energy Resolution and Jet Energy Scale 
//-----R=1.0
const int totRanges = 12;
TH1D* pTReco_pTTruth_Dis[totRanges];
const int end_of_Ranges = 13;
float pTRanges[end_of_Ranges] = {100,200,300,400,500,600,700,800,900,1000,1100,1200,1300};


const int dj_totRanges = 20;
float pTRanges_DJ[] = {100,112,126,141,158,178,200,224,251,282,316,355,398,447,501,562,631,708,794,891,1000};



//----Total Number of type of Jets
const int tot_TypeJets = 2; //Reco and truth
//----Total Number of Jets
const int total_Jets = 100000;
//----Total Number of Jet Variables
const int tot_JetParam = 3;
//----Turning off or on the branches
int branchStatus[] = {false,true}; //true means it is ON and false is OFF
//----This is to find the type of Jet (Reco or Truth)
enum TruthOrReco {trthJet,recoJet};
string TypeOfJet[] = {"Truth","Reco"}; 

//---Parameters of Jets
enum JetPara {pTPar,etaPar,phiPar};

string jetParam[]= {"pt","eta","phi"};
//-----This is  the truth jet information arra
std::vector<float> * truthJets_Info[Tot_Radii][tot_JetParam] = {};  
//----Reco Jets Information
std::vector<float> * recoJets_Info[Tot_Radii][tot_JetParam] = {};
int total_TruthJets[Tot_Radii] ={};
int total_RecoJets[Tot_Radii] = {};
//----Hold number of jets in event array
int tot_JetsInEvent[tot_TypeJets][Tot_Radii] = {};

//Unfolding
const int totIterations = 18;
const int nToys = 100;

//These are to create your Response Matrices (for unfolding)
float pTCutUnfol[] = {150,100};


//pp MC
//pT Truth Cuts
double pT_TruthCuts_ppMC[2][8] ={{241,241,245,211,230,230,241,241},{158,158,160,158,158,158,158,158}};
//pT Reco Cuts
double pT_RecoCuts_ppMC[2][8] ={{294,294,294,294,294,294,294,294},{199,199,194,198,198,198,197,197}};


//---------------------------------------------------------------------------//
//---------------------------------Pb+Pb MC----------------------------------//
//---------------------------------------------------------------------------//

//Total Number of Jets
const int totalJets = 2;


//SumET Distributions (For both MC and Data)
int maxSumETValue = 5000; //5 TeV

//----JZ Samples
const int numJZSamples =4; 


//Verticies
int num_vert = 0;

//FCal Information
const int tot_FCals = 2;
float fcalEnergy_A;
float fcalEnergy_C;


//Centrality Bins
const int  totCentralityBins = 4;
float sumET[] = {2995.94, 1378.92, 533.608, 66.402};
string centBins[] = {"0_10","10_30","30_50","50_80"};

const int totCentBins_DJ = 5;

string centTimBins[]={"0_10","10_20","20_40","40_60","60_80"};
float Tim_sumETBins[]={2989.31,2046.51,875.41,289.595,63.719};

enum centBin_DJ {cent_0_10,cent_10_20,cent_20_40,cent_40_60,cent_60_80};
enum centBin {cent0_10,cent10_30,cent30_50, cent50_80};

//Cent Bins To Compare to 2015 rAA Meas.
const int num2015MeasBins = 8;
enum meas2015_SumET {centIndex0_10,centIndex10_20,centIndex20_30,centIndex30_40,centIndex40_50,centIndex50_60,centIndex60_70,centIndex70_80};
 
string centBins_2015Meas[]={"0_10","10_20","20_30","30_40","40_50","50_60","60_70","70_80"};
double sumET_Bins2015Meas[] = {2.99594*1e+3,2.05577*1e+3,1.37892*1e+3,0.885172*1e+3,0.533608*1e+3,0.29617*1e+3,0.148625*1e+3,0.066402*1e+3};


//This Is ONLY for MC 
string centBining[2][8] = {{"0_10","10_20","20_30","30_40","40_50","50_60","60_70","70_80"},{"0_10","10_20","20_40","40_60","60_80","","",""}}; 

double sumETVal_CentBins[2][8] = {{2.99594*1e+3,2.05577*1e+3,1.37892*1e+3,0.885172*1e+3,0.533608*1e+3,0.29617*1e+3,0.148625*1e+3,0.066402*1e+3},
				  {2989.31,2046.51,875.41,289.595,63.719,0,0,0}};


//This is for the trigger efficiency plot for HLT_j85 
string centBins_TrigEff[]={"0_10","10_20","20_30","30_50","50_80"};
enum centbins_TrigEff{cent_bin_0_10,cent_bin_10_20,cent_bin_20_30,cent_bin_30_50,cent_bin_50_80};
int trigEff_CentBins = 5;
//-----This is  the truth jet information arra
std::vector<float> * truthJets_InfoPbPbMC[Tot_Radii][tot_JetParam] = {};
//----Reco Jets Information
std::vector<float> * recoJets_InfoPbPbMC[Tot_Radii][tot_JetParam] = {};
int total_TruthJets_PbPbMC[Tot_Radii] ={};
int total_RecoJets_PbPbMC[Tot_Radii] = {};

//Weights
//----JZ6
double eps_e = 0.0008433528;
double sigma_e = 0.00134; //[nb]
double num_e = 19349;


//-----Values needed for weight second slice (1200-800 GeV)
double eps_d = 0.0021805;
double sigma_d = 0.225;
double num_d = 7734425;

//-----Values needed for JZ4 Sample

double eps_c = 0.0045851;
double sigma_c = 2.6602*pow(10,1);
double num_c = 7534207;
//-----Values needed for weight 0 slice (400-160 GeV)
double eps_b = 0.005288;
double sigma_b = 4.7195*pow(10,3);
double num_b = 7727675;

//------Crossections and efficiencies for the first JZ sample (160-60 GeV)
double eps_a = 0.0042785;
double sigma_a = 6.40*pow(10,5);
double num_a = 7736371;


double weights_PbPb[] = {(eps_a*sigma_a)/(num_a),(eps_b*sigma_b)/(num_b),(eps_c*sigma_c)/(num_c),(eps_d*sigma_d)/(num_d), (eps_e*sigma_e)/(num_e)};

//These bins are applies to the 
//Efficiency plots (i.e. sumpT Cuts Eff.)
const int totBins = 30;
float finerBins_Eff[totBins]={};
int x_lowedge = 0;
int size_of_bin = 20;



float OnePerCentralityBins[] = {0,50.2,100.4,150.6,200.8,
                              251,301.2,351.4,401.6,451.8,
                              502,552.2,602.4,652.6,702.8,
                              753,803.2,853.4,903.6,953.8,
                              1004,1054.2,1104.4,1154.6,
                              1204.8,1255,1305.2,1405.6,1455.8,
                              1506,1556.2,1606.4,1656.6,1706.8,
                              1757,1807.2,1857.4,1907.6,1957.8,
                              2008,2058.2,2108.4, 2158.6,2208.8,
                              2259,2309.2,2359.4,2409.6,2459.8,
                              2510,2560.2,2610.4,2660.6,2710.8,
                              2761,2811.2,2861.4,2911.6,2961.8,
                              3012,3062.2,3112.4,3162.6,3212.8,
                              3263,3313.2,3363.4,3413.6,3463.8,
			      3514,3564.2,3614.4,3664.6,3714.8,
			      3765,3815.2,3865.4,3915.6,3965.8,
			      4016,4066.2,4116.4,4166.6,4216.8,
			      4267,4317.2,4367.4,4417.6,4467.8,
			      4518,4568.2,4618.4,4668.6,4718.8,
			      4769,4819.2,4869.4,4919.6,4969.8,
			      5020}; 


//Centrlity Weights
bool addCentWeights = false;


//Pb+Pb MC
//pT Truth Cuts
double pT_TruthCuts[totalJets][num2015MeasBins] ={{241,220,211,177,177,177,177,177},{66,66,75,66,66,66,66,66}};



std::map<int, std::string> centBinmap_2015Meas {
  {0,"0-10%"},
  {1,"10-20%"},
  {2,"20-30%"},
  {3,"30-40%"},
  {4,"40-50%"},
  {5,"50-60%"},
  {6,"60-70%"},
  {7,"70-80%"}
};


//---------------------------------------------------------------------------//
//---------------------------------Pb+Pb Data--------------------------------//
//---------------------------------------------------------------------------//






//TAA For 2015 Cent Bins
std::map<int, double> tAA_2015map {
  {0,23.35*1e-6},
  {1,14.33*1e-6},
  {2,8.76883*1e-6},
  {3,5.08891*1e-6},
  {4,2.74507*1e-6},
  {5,1.35167*1e-6},
  {6,0.601251*1e-6},
  {7,0.23877*1e-6}

};


//CentBins Using
bool tim_sumETBins = false;


//Triggers
//R=1.0
bool HLT_j200_a10_ion_L1J50 = false;
bool HLT_j180_a10_ion_L1J50 = false;
bool HLT_j150_a10_ion_L1J50 = false;

string trigNames[] = {"HLT_j150_a10_ion_L1J50","HLT_j85_ion_L1J30"}; 

//R=0.4
bool HLT_j85_ion_L1J30 = false;

bool trigEff = false; //Trigger Efficiency Plots for using HLT_j85* for R=1.0 jets 


//-----Jets
int totJetsPbPbData = 2;
int tot_PbPbDataJets[totRadii] ={};
std::vector<float> * PbPbDataJets_Info[totRadii][totJetParam] = {};
float pTCut_PbPbData[2][8] = {{294,294,245,230,230,220,220,220},{100,100,109,100,100,100,100,100}};
std::vector<float> * jetEnergy = {};

//FCal Information
float fcalEnergyA;
float fcalEnergyC;


int vertType_PbPbData[v_size];
int nvert_PbPbData = 0;

//
const int withAndwoWeights = 2; 
string withOrwoWeights[] = {"withWeight", "withoutWeight"};
enum withOrwithoutWeight  {withWeight,withoutWeight};


//EtaCut
double etacut = 1.5;

//Unfolding
//pTDistributions 

//This variable is se to true if you want to use
//the 2015 Meas. rAA binning
bool rAA2015Bins = false;
//Jet Rate Binning
bool jet_RateBins = false;
//Large-R Jet Analysis Binning
bool largeRbins = true;

//RAA BINNING FOR YOUR ANALYSIS
const int binsRAA = 23;
//float rAA_Bins[] = {100,112,125,141,158,177,199,223,251,281,316,354,398,501,630,999,1300};
float rAA_Bins[] = {20,30,40,50,70,90,100,112,125,141,158,177,199,223,251,281,316,354,398,501,630,999,1200,1300};


//RAA BINNING FOR THE 2015 MEASUREMENT
const int bins2015 = 20;
const int cent2015Bins = 8;
double rAA_2015Bins[cent2015Bins][21]= {{60,70,80,90,100,112,125,141,158,177,199,223,251,281,316,354,398,501,630,999,1300},//0-10%
					{50,60,70,80,90,100,112,125,141,158,177,199,223,251,281,316,398,630,800,1000,1300}, //10-20%
					{50,60,70,80,90,100,112,125,141,158,177,199,223,251,281,316,398,630,800,1000,1300}, //20-30%
					{40,50,60,70,80,90,100,112,125,141,158,177,199,223,251,316,398,630,800,1000,1300}, //30-40%
					{20,30,40,50,56,63,70,79,89,100,112,125,141,158,177,199,251,398,545,692,839}, //40-50%
                                        {20,30,40,50,56,63,70,79,89,100,112,125,141,158,177,199,251,398,545,692,839}, //50-60%
                                        {20,30,40,50,56,63,70,79,89,100,112,125,141,158,177,199,251,398,545,692,839}, //60-70%
                                        {20,30,40,50,56,63,70,79,89,100,112,125,141,158,177,199,251,398,545,692,839}, //70-80%
};

float x_lo = 100; //GeV
float x_hi = 1300; //GeV
//Response Matrix 
int numbins = 1200;//assuming your pT Ramnge is 100-1300 GeV (1GeV sized bin)

//Jet Rate Bins for 2015 Measurement 
//Data
int jetRateBins = 10;
double jetRate[] = {100,126,158,200,251,316,398,501,631,800,1000};
//MC
int jetRateMCBins = 15;
double jetRateMC[] = {50,63,79,100,126,158,200,251,316,398,501,631,800,1000,1200,1400};

//RAA Points from CMS
double RAA_CMS_val[] = {0.8111111111111111};
double RAA_CMS_pT[] = {754.1194441104637};


//Bin Width for this point ^^ [500,1000]
//For the y-values
double RAA_CMS_sysErh[]= {0.9031746031746031-0.8111111111111111};
double RAA_CMS_sysErl[]= {0.8111111111111111-0.719047619047619};
double RAA_CMS_statEr[]= {0.8111111111111111-0.873015873015873};
//For the x-values
double ex_CMS[]={500};

//This is just to check the number of counts in each pT Bin
//This is to help us decide how wide our pT bins should be
//To see how much we are limited by statistis
const int set_numBins=100;
double bins[set_numBins] = {};
enum edges{lowEdge,highEdge};

//I am defining the edges
float edgepTBins[totRadii][cent2015Bins][2] = {{{294,1300},{262,1300},{245,1300},{230,1300},{230,1300},{220,999},{220,999},{220,999}},
					       {{100,1300},{100,1300},{100,1300},{100,1300},{100,1300},{100,999},{100,999},{100,999}}};


enum startAndend{start,end};
