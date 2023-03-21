#include <assert.h>
#include <cmath>
#include <fstream>
#include <stdio.h>
#include <vector>
#include <math.h>
#include <tgmath.h>
#include <iostream>
#include <TChain.h>
#include <TProfile.h>
#include <TLine.h>
#include <TLegend.h>
#include <TF1.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TRandom3.h"
#include <TStyle.h>
#include <TAttMarker.h>
#include <TCanvas.h>
#include <TPaveText.h>
#include <TGaxis.h>
#include <vector>
#include "NewTreeVariables.h"
#include "bFunctions.h"
#include <typeinfo>


using namespace std;


int rAANominal(bool rAAPlot = true,bool UnfoldedOrNot = true, bool jetRatePbPb=false, bool jetRatepp = false, bool centBins2015 = true,string centBin= "70_80",int jetR=R10, int numIterations=15, bool comp2015 = false, bool justJetRate = false, bool largeR=true, bool debg = true, bool savepdf = false, bool plot_systematics=true,bool NopTShapeWeight=false,bool NoCentWeight = false,bool useOldFiles = false){

   int central_bin = centBinsMap_2015Meas[centBin];
  

   if(debg)cout << __LINE__ << endl;
   Color_t colorR4[8] = {kOrange+3,kYellow+2,kYellow+3,kOrange+2,kTeal+3,kRed+2,kCyan+1,kAzure+4};

  if(debg)cout << __LINE__ << endl;
  if(debg)cout	<< __LINE__ << endl;
  int numCentBins = 0;
  double etaRange = 0;
  double eta_range = 0;
  //y-Axis of plot
  double yMin = 0;
  double yMax = 0;
  //x-axis of plot
  double xMin = 200;
  double xMax = 900;
  //TLegend corners
  double pt1 =0;
  double pt2 =0;
  double pt3 =0;
  double pt4 =0;
  
  double offSet[] = {1,1e+2,1e+4,1e+6,1e+8,1e+10,1e+12,1e+14};
  string offsetTag[] {"","#times 10^{2}","#times 10^{4}","#times 10^{6}","#times 10^{8}","#times 10^{10}","#times 10^{12}","#times 10^{14}"};

  
  
  //Uploading files
  //---------Pb+Pb
  TFile * data_PbPb = NULL;
  TFile * data_UnfPbPb = NULL;


  //Jet Rate in 2015 Pb+Pb Data
  TGraph *meas2015_graph[num2015MeasBins];
  TH1D *meas2015_hist[num2015MeasBins];
  //Systematica Erros
  TH1D *sys2015Errors_h[num2015MeasBins];
  TH1D *sys2015Errors_l[num2015MeasBins];
  //Statistical Erros
  TH1D *stat2015Errs[num2015MeasBins];

  //RAA 2015 Meas
  TH1D *rAAVal_hist[num2015MeasBins];
  TH1D *rAAVal_graph[num2015MeasBins];

  //Systematic Errors For rAA values
  TH1D *rAAsys_ePlus[num2015MeasBins];
  TH1D *rAAsys_eMinus[num2015MeasBins];

  //Stat Erros
  TH1D *rAA_statEr[num2015MeasBins];

  //Dhanush's 2018 RAA Values
  TH1D *rAAVal_2018[num2015MeasBins];
  TH1D *rAAVal_2018Sys[num2015MeasBins];
  TFile *file_2018RAA = new TFile("Dhanush_Files/2018PbPb_inclusive_jetraa.root","READ");

  string collSyst = "";
  if(jetRatePbPb){
    collSyst = "PbPb";
  }else if(jetRatepp){
    collSyst = "pp";
  }else if(rAAPlot){
    collSyst = "RAA";
  }
  
    string NocentWeight_tag = "";
  if(NoCentWeight){
    NocentWeight_tag= "_NoCentWghts";
  }

  string NopTShpWghts_tag ="";
  if(NopTShapeWeight){
    NopTShpWghts_tag = "_NopTShpWght";
  }

  string usingOldFiles_tag = "";
  if(useOldFiles){
    usingOldFiles_tag = "_UsingSeptember2022RootFiles";
  }
  
  string NOMINAL_TAG ="";
  if(!NoCentWeight && !NopTShapeWeight)NOMINAL_TAG = "Nominal";
  
  if(debg)cout	<< __LINE__ << endl;
  
  if(jetR==R4){
     data_PbPb = new TFile("rAADir/rootFiles/RawHistograms_PbPb_Data.root","READ");      
     data_UnfPbPb = new TFile("/Users/berenicegarcia/Desktop/Berenice/Spring_2022/LargeRJets/LargeRJetNewCode/Locally/DiagPlots/rootFiles/NominalLargeR/hist_26282548_05302022_Unfolded_PbPbData_ATLAS_Official_RAA_Binning_17Iters_10000Toys_2018LargerAnalysis_etaRange15Nominal_R4.root","READ");
     
  }else if(jetR==R10){
    if(!useOldFiles){
      data_PbPb = new TFile("/Users/berenicegarcia/Desktop/Berenice/Spring_2022/LargeRJets/LargeRJetNewCode/Locally/DiagPlots/rootFiles/NominalLargeR/RawHistograms_R10_PbPb_NopTShpWghts_2018LargerAnalysis_eta15_Data.root","READ");
      data_UnfPbPb = new TFile(Form("/Users/berenicegarcia/Desktop/Berenice/Spring_2022/LargeRJets/LargeRJetNewCode/Locally/DiagPlots/rootFiles/NominalLargeR/hist_26282548_05302022_Unfolded_PbPbData_ATLAS_Official_RAA_Binning_17Iters_10000Toys_2018LargerAnalysis_etaRange15%s%s%s_R10.root",NOMINAL_TAG.c_str(),NocentWeight_tag.c_str(),NopTShpWghts_tag.c_str()),"READ");
    }else if(useOldFiles){
      data_UnfPbPb = new TFile("/Users/berenicegarcia/Desktop/Berenice/Spring_2022/LargeRJets/LargeRJetNewCode/Locally/rAADir/rootFiles/hist_26282548_05302022_Unfolded_PbPbData_Eta15_ATLAS_Official_RAA_Binning_3Iters_10000Toys_NewpTBins.root","READ");
      data_PbPb = new TFile("/Users/berenicegarcia/Desktop/Berenice/Spring_2022/LargeRJets/LargeRJetNewCode/Locally/rAADir/rootFiles/RawHistograms_PbPb_Data.root","READ");
    }
  }

  if(debg)cout	<< __LINE__ << endl;
    numCentBins = 8;
    etaRange = 1.5*2;
    eta_range = 1.5;
    
    if(rAAPlot){
      pt1 =0.7913741;pt2 =0.6763848;pt3 =0.9909729;pt4 =0.8381924;
      yMin=1e-17;
      yMax=1e+16;
     
    }
    if(!centBins2015){
      pt1 =0.7514735;pt2 =0.2332016;pt3 =0.9823183;pt4 =0.5085639;
      yMin=1e-8;
      yMax=1e+9;
    }
 
    if(jetRatepp){
      numCentBins = 1;
    }
    if(debg)cout	<< __LINE__ << endl;
    TFile *data_pp;
    if(jetR==R4){
      data_pp = new TFile("/Users/berenicegarcia/Desktop/Berenice/Spring_2022/LargeRJets/LargeRJetNewCode/Locally/DiagPlots/rootFiles/NominalLargeR/hist_26282548_05302022_Unfolded_ppData_ATLAS_Official_RAA_Binning_17Iters_10000Toys_2018LargerAnalysis_etaRange15Nominal_R4.root","READ");
      
    }else if(jetR==R10){
      data_pp = new TFile(Form("/Users/berenicegarcia/Desktop/Berenice/Spring_2022/LargeRJets/LargeRJetNewCode/Locally/DiagPlots/rootFiles/NominalLargeR/hist_26282548_05302022_Unfolded_ppData_ATLAS_Official_RAA_Binning_17Iters_10000Toys_2018LargerAnalysis_etaRange15%s%s%s_R10.root",NOMINAL_TAG.c_str(),NocentWeight_tag.c_str(),NopTShpWghts_tag.c_str()),"READ");//ATLAS 2015  rAA Bins (fixed binning)
    }

    if(useOldFiles){
      data_pp = new TFile("/Users/berenicegarcia/Desktop/Berenice/Spring_2022/LargeRJets/LargeRJetNewCode/Locally/ppCrossSec/rootFiles/hist_26282548_05302022_Unfolded_ppData_Eta15_ATLAS_Official_RAA_Binning_3Iters_10000Toys_NewpTBins_TEST.root","READ");
    }

    
  //Print Out
  cout << "This is the cent bin we will be plotting: " << centBin.c_str() << endl; 
  cout << "For radius: " << jetRadius[jetR] << endl;

  //Histograms
  //MC
  TH1D* pTDis_Truth[totRads][numCentBins];
  TH1D* pTDis_Reco[totRads][numCentBins];
  //Data
  //Not Unfolded
  TH1D *pTDis_PbPbData[totRads][numCentBins];
  TH1D *pTDis_ppData[totRads][numCentBins];
  //Unfolded
  TH1D *pTDis_Unfo_PbPbData[totRads][numCentBins][numIterations];
  TH1D *pTDis_Unfo_ppData[totRads][numCentBins];
  
  string names_PbPbData[totRads][numCentBins];
  string names_Unf_PbPbData[totRads][numCentBins][numIterations];

    //Number of Iterations you want to apply to each of your pTDistributions
  //Before you make the rAA plot
  const int centbins2015 = 8;

  //The reason this array stores 3 ios for the Nominal Iter Num/StartBin/EndBin
  int nominalInfo_PbPb[centbins2015][3] = {};
  int nominalInfo_pp[centbins2015][3] = {};
  enum three_param{nominalNum,startBin,endBin};
  
  //Nominal Number of Iterations
  TFile *pp_nominalInfo;
  if(debg)cout	<< __LINE__ << endl;
  TFile *pbpb_nominalInfo;
  if(jetR==R4){
    pp_nominalInfo =new TFile("/Users/berenicegarcia/Desktop/Berenice/Spring_2022/LargeRJets/LargeRJetNewCode/Locally/DiagPlots/rootFiles/NominalLargeR/nominalIter_pp_2018LargerAnalysis_withpTShapeWeights_etaRange15.root","READ");
    pbpb_nominalInfo  = new TFile("/Users/berenicegarcia/Desktop/Berenice/Spring_2022/LargeRJets/LargeRJetNewCode/Locally/DiagPlots/rootFiles/NominalLargeR/nominalIter_PbPb_2018LargerAnalysis_withpTShapeWeights_etaRange15.root","READ");
  }else if(jetR==R10){
    pp_nominalInfo =new TFile("/Users/berenicegarcia/Desktop/Berenice/Spring_2022/LargeRJets/LargeRJetNewCode/Locally/DiagPlots/rootFiles/NominalLargeR/nominalIter_2018LargerAnalysispp_withpTShapeWeights_etaRange15_R10.root","READ");
    pbpb_nominalInfo  = new TFile("/Users/berenicegarcia/Desktop/Berenice/Spring_2022/LargeRJets/LargeRJetNewCode/Locally/DiagPlots/rootFiles/NominalLargeR/nominalIter_2018LargerAnalysisPbPb_withpTShapeWeights_etaRange15_R10.root","READ");
  }
  if(debg)cout	<< __LINE__ << endl;
  
  for(int iCentBin =0; iCentBin < centbins2015; iCentBin++){
    if(debg)cout << __LINE__ << endl;
    //Grab Nominal Iteration number and store in array
    nominalInfo_pp[iCentBin][nominalNum] = ((TH1D*)pp_nominalInfo->Get(Form("ppdata_Cent_%s_R%d_IterNum",centBins_2015Meas[iCentBin].c_str(),jetRadius[jetR])))->GetBinContent(1);
     if(debg)cout << __LINE__ << endl;
    nominalInfo_PbPb[iCentBin][nominalNum] = ((TH1D*)pbpb_nominalInfo->Get(Form("pbpbdata_Cent_%s_R%d_IterNum",centBins_2015Meas[iCentBin].c_str(),jetRadius[jetR])))->GetBinContent(1);
    if(debg)cout   << __LINE__ << endl;
    cout << "This is iCentBin number: " << iCentBin << endl;
    
    //Grab Start Bin number and store in array
    nominalInfo_pp[iCentBin][startBin] = ((TH1D*)pp_nominalInfo->Get(Form("ppdata_Cent_%s_R%d_StartBinNum",centBins_2015Meas[iCentBin].c_str(),jetRadius[jetR])))->GetBinContent(1);
    nominalInfo_PbPb[iCentBin][startBin] = ((TH1D*)pbpb_nominalInfo->Get(Form("pbpbdata_Cent_%s_R%d_StartBinNum",centBins_2015Meas[iCentBin].c_str(),jetRadius[jetR])))->GetBinContent(1);

    cout << "This is the start bin number for this cent: " << nominalInfo_PbPb[iCentBin][startBin] << endl;
    
    //Grab End Bin number and store in array
    nominalInfo_pp[iCentBin][endBin] = ((TH1D*)pp_nominalInfo->Get(Form("ppdata_Cent_%s_R%d_EndBinNum",centBins_2015Meas[iCentBin].c_str(),jetRadius[jetR])))->GetBinContent(1);
    nominalInfo_PbPb[iCentBin][endBin] = ((TH1D*)pbpb_nominalInfo->Get(Form("pbpbdata_Cent_%s_R%d_EndBinNum",centBins_2015Meas[iCentBin].c_str(),jetRadius[jetR])))->GetBinContent(1);
    
  }



  cout << "This is the centbin number 0 and the start bin for this cent bin is: " <<  nominalInfo_PbPb[0][startBin] << endl;
  
  for(int iJetR = 0; iJetR < totRads; iJetR++){
    if(iJetR!=jetR)continue;
     for(int iCentBin=0;iCentBin < numCentBins; iCentBin++){
       
       cout << "************Nominal iteration number pp: " << nominalInfo_pp[iCentBin][nominalNum] << endl;
       cout << "************Nominal iteration number Pb+Pb: " << nominalInfo_PbPb[iCentBin][nominalNum] << endl; 
       cout << "R/iCentBin: " << iJetR << "/" << iCentBin << endl;
       cout << "Grabbed this pp histo: " << Form("Unfolded_ppData_R%d_%dIter_%s",jetRadius[iJetR],nominalInfo_pp[iCentBin][nominalNum],centBins_2015Meas[iCentBin].c_str()) << endl;
       pTDis_ppData[iJetR][iCentBin] = (TH1D*) data_pp->Get(Form("Unfolded_ppData_R%d_%dIter_%s",jetRadius[iJetR],nominalInfo_pp[iCentBin][nominalNum],centBins_2015Meas[iCentBin].c_str()));
       cout << __LINE__ << endl;

       if(iJetR==0 && iCentBin==0)cout << "PRINT CONTENT OF FIRST BIN!! : " << pTDis_ppData[iJetR][iCentBin]->GetBinContent(10) << endl;
       
     }


     cout << __LINE__ << endl;
     cout << "PRINT CONTENT OF FIRST BIN!! : " << pTDis_ppData[0][0]->GetBinContent(10) << endl;
    
    for(int iCentBin =0; iCentBin < 8; iCentBin++){
       if(centBins2015){
	 
	pTDis_PbPbData[iJetR][iCentBin] = (TH1D*) data_PbPb->Get(Form("R%d_Cent_%s",jetRadius[iJetR],centBins_2015Meas[iCentBin].c_str()));
	 
	 cout <<	__LINE__ << endl;
	for(int iIter =1; iIter <= numIterations; iIter++){
	   cout <<	__LINE__ << endl;
	  pTDis_Unfo_PbPbData[iJetR][iCentBin][iIter] = (TH1D*) data_UnfPbPb->Get(Form("Unfolded_PbPbData_R%d_%dIter_%s",jetRadius[iJetR],iIter,centBins_2015Meas[iCentBin].c_str()));
	  cout << "iIter: " << iIter << endl;
	  if(iCentBin==0 && iIter ==15 && iJetR==1 && false){
	    cout << "Bin content of bin number 1:" << pTDis_Unfo_PbPbData[iJetR][iCentBin][iIter]->GetBinContent(1) << endl;
	    cout << "Bin center: " << pTDis_Unfo_PbPbData[iJetR][iCentBin][iIter]->GetBinCenter(1) << endl;
	  }
	  cout << __LINE__ << endl;
	}//Looping over Iterations 
	cout << __LINE__ << endl;
       }

    }//Cent Bin Loop
    cout << __LINE__ << endl;
  }//Radii Loop

  cout << __LINE__ << endl;
  cout << "PRINT CONTENT OF FIRST BIN!! : " << pTDis_ppData[0][0]->GetBinContent(10) << endl;
  
  //Normalization Factors
  cout << __LINE__ << endl;
  double numEvents_centBin[] = {7.383e+8*LumNumPbPbData[jetR],7.383e+8*LumNumPbPbData[jetR],2*7.383e+8*LumNumPbPbData[jetR],2*7.383e+8*LumNumPbPbData[jetR],2*7.383e+8*LumNumPbPbData[jetR]};
  cout << __LINE__ << endl;
  double numEventsMinBias = 7.383e+8*LumNumPbPbData[jetR];


  double x_val[centbins2015][20]={};
  double y_val[centbins2015][20]={};
   //Systematic Errors
  double sys_ePlus[centbins2015][20] = {};
  double sys_eMinus[centbins2015][20] = {};
  double valx_error[centbins2015][20] = {};
  //Statistical Errors
  double stat_err[centbins2015][20] = {};

  //These arrays belong to Dhanush's Analysis
  double sys_rAAVals[centbins2015][20] = {};
  double rAAVals_2018[centbins2015][20] = {};
    
  cout << __LINE__ << endl;
  TFile *systematicsFile;
  if(rAAPlot && plot_systematics)systematicsFile= new TFile(Form("Systematics/NominalLargeR/totSystematics/2018LargerAnalysis_FirstTotSysUncert_RAA_R%d.root",jetRadius[jetR]),"READ"); 
  if(jetRatePbPb && plot_systematics)systematicsFile= new TFile(Form("Systematics/NominalLargeR/totSystematics/2018LargerAnalysis_FirstTotSysUncert_PbPb_R%d.root",jetRadius[jetR]),"READ");
  if(jetRatepp && plot_systematics)systematicsFile= new TFile(Form("Systematics/NominalLargeR/totSystematics/2018LargerAnalysis_FirstTotSysUncert_pp_R%d.root",jetRadius[jetR]),"READ");

  cout << __LINE__ << endl;
    
  TCanvas * canv = new TCanvas(Form("RAA_R%d%s%s%s%s",jetRadius[jetR],NOMINAL_TAG.c_str(),NocentWeight_tag.c_str(),NopTShpWghts_tag.c_str(),usingOldFiles_tag.c_str()),Form("RAA_R%d%s%s%s%s",jetRadius[jetR],NOMINAL_TAG.c_str(),NocentWeight_tag.c_str(),NopTShpWghts_tag.c_str(),usingOldFiles_tag.c_str()),765,1618,999,711);
    gStyle->SetOptStat(0);
    gPad->SetTicks(1);
    canv->Range(64.57942,-0.2185758,1269.776,1.332886);
    canv->SetFillColor(0);
    canv->SetBorderMode(0);
    canv->SetBorderSize(2);
    if(rAAPlot){
      canv->SetRightMargin(0.224674);
      canv->SetTopMargin(0.1406706);
      canv->SetBottomMargin(0.1406706);
    }else if(jetRatePbPb){
        canv->SetLeftMargin(0.1514544);
	canv->SetRightMargin(0.1725175);
	canv->SetTopMargin(0.1457726);
	canv->SetBottomMargin(0.1355685);
    }

    cout << __LINE__ << endl;
    cout << __LINE__ << endl;
    canv->SetFrameBorderMode(0);
    canv->SetFrameBorderMode(0);
    cout << __LINE__ << endl;
    if(jetRatePbPb || jetRatepp){
      canv->SetLogx();
      canv->SetLogy();
      cout << __LINE__ << endl;
    }
    cout << "This is the centbin number 0 and the start bin for this cent bin is: " <<  nominalInfo_PbPb[0][startBin] << endl;


    cout << "PRINT CONTENT OF FIRST BIN!! : " << pTDis_ppData[0][0]->GetBinContent(10) << endl;
    
    cout << __LINE__ << endl;
    //Legends
    TLegend *leg;
    if(rAAPlot)leg= new TLegend(0.6359077,0.1924198,0.77332,0.3483965,NULL,"brNDC");
    if(jetRatePbPb)leg= new TLegend(0.1805416,0.1588921,0.3781344,0.3571429,NULL,"brNDC");
    if(jetRatepp)leg= new TLegend(0.1273821,0.1195335,0.4593781,0.1865889,NULL,"brNDC");
    
    for(int iCentBin=0; iCentBin < numCentBins; iCentBin++){
      if((iCentBin%2==0) && rAAPlot)continue;
      //if(iCentBin!=0)continue;

      cout << "This is the centbin number 0 and the start bin for this cent bin is: " <<  nominalInfo_PbPb[0][startBin] << endl;
      
      if(!UnfoldedOrNot){
	
	cout << __LINE__ << endl;
        pTDis_PbPbData[jetR][iCentBin]->Scale(1/numEventsMinBias);
	pTDis_PbPbData[jetR][iCentBin]->Scale(1.,"width");
	pTDis_PbPbData[jetR][iCentBin]->Scale(1/etaRange);
	pTDis_PbPbData[jetR][iCentBin]->Scale(1/tAA_2015map[iCentBin]);
      }else{
	cout << "jetR: " << jetR << endl;
	cout << "iCentBin: " << iCentBin << endl;
	cout << "nominalInfo_PbPb[iCentBin][nominalNum]: " << nominalInfo_PbPb[iCentBin][nominalNum] << endl;
	cout <<      __LINE__ << endl;
	
	pTDis_Unfo_PbPbData[jetR][iCentBin][nominalInfo_PbPb[iCentBin][nominalNum]]->Scale(1/numEventsMinBias);

	pTDis_Unfo_PbPbData[jetR][iCentBin][nominalInfo_PbPb[iCentBin][nominalNum]]->Scale(1.,"width");

	pTDis_Unfo_PbPbData[jetR][iCentBin][nominalInfo_PbPb[iCentBin][nominalNum]]->Scale(1/etaRange);

	pTDis_Unfo_PbPbData[jetR][iCentBin][nominalInfo_PbPb[iCentBin][nominalNum]]->Scale(1/tAA_2015map[iCentBin]);
	cout << "For iCentBin: " << iCentBin << endl;
	cout << "the scaling factor is: " << tAA_2015map[iCentBin]*numEventsMinBias << endl;
	cout << "etaRange: " << etaRange << endl;
	
	
	
      }
      cout <<      __LINE__ << endl;
      ///This isto plot the Sys. Uncert.

       cout << "PRINT CONTENT OF FIRST BIN!! : " << pTDis_ppData[0][0]->GetBinContent(10) << endl;
      
      if(iCentBin==0 || iCentBin==1){
	cout << __LINE__ << endl;
	for(int iCentBin=0;iCentBin < numCentBins; iCentBin++){
	  cout << "iCentBin: " << iCentBin << endl;
	  cout << "iJetR: " << jetR << endl;
	   
	  cout << __LINE__ << endl;
 	  pTDis_ppData[jetR][iCentBin]->Scale(1.,"width");
	    cout << __LINE__ << endl;
	  pTDis_ppData[jetR][iCentBin]->Scale(1/etaRange);
	  cout << __LINE__ << endl;
	  pTDis_ppData[jetR][iCentBin]->Scale(1/ppDataLumiVals[jetR]);
	}
	cout << __LINE__ << endl;
      }



      
      cout << "This is the centbin number 0 and the start bin for this cent bin is: " <<  nominalInfo_PbPb[0][startBin] << endl;

      
      //plotting only jet rates 
      if(jetRatePbPb || jetRatepp){
	
	// Create a new histogram with the same binning as rAA, but without its contents
	TH1D *temp_jetRate;
	TH1D *hist;
	int jetRate_bin_start = -1;
	int jetRate_bin_end =-1;
	if(jetRatePbPb){
	  hist = (TH1D*)pTDis_Unfo_PbPbData[jetR][iCentBin][nominalInfo_PbPb[iCentBin][nominalNum]];
	  temp_jetRate= (TH1D*)pTDis_Unfo_PbPbData[jetR][iCentBin][nominalInfo_PbPb[iCentBin][nominalNum]]->Clone("temp_jetrate");
	  jetRate_bin_start = nominalInfo_PbPb[iCentBin][startBin];
	  jetRate_bin_end = nominalInfo_PbPb[iCentBin][endBin];
	}


	if(jetRatepp){
	  hist = (TH1D*)pTDis_ppData[jetR][iCentBin];
	  temp_jetRate= (TH1D*)pTDis_ppData[jetR][iCentBin]->Clone("temp_jetrate");
	  jetRate_bin_start = nominalInfo_pp[iCentBin][startBin];
	  jetRate_bin_end = nominalInfo_pp[iCentBin][endBin];

	}
	temp_jetRate->Reset();
	
	cout << "*******************this is the start bin: " << jetRate_bin_start << endl;
	cout << "*******************this is the endbin: " << jetRate_bin_end << endl;
	 for (int iBin = jetRate_bin_start+1; iBin <= jetRate_bin_end; iBin++) {
	  
          temp_jetRate->SetBinContent(iBin, hist->GetBinContent(iBin));
          temp_jetRate->SetBinError(iBin, hist->GetBinError(iBin));


	  x_val[iCentBin][iBin] = hist->GetBinCenter(iBin);
          y_val[iCentBin][iBin] = hist->GetBinContent(iBin)*offSet[iCentBin];
	  if(plot_systematics){
	    //going to setup systematics
	    TH1D *sys_hist = (TH1D*)systematicsFile->Get(Form("TotSysUncert_R%d_Cent_%s",jetRadius[jetR],centBins_2015Meas[iCentBin].c_str()));
	    sys_ePlus[iCentBin][iBin] = (sys_hist->GetBinContent(iBin)*hist->GetBinContent(iBin)*offSet[iCentBin])/100;
	    sys_eMinus[iCentBin][iBin] = (sys_hist->GetBinContent(iBin)*hist->GetBinContent(iBin)*offSet[iCentBin])/100;
	    valx_error[iCentBin][iBin] = abs(temp_jetRate->GetBinLowEdge(iBin)-temp_jetRate->GetBinCenter(iBin));
	  }


	 }

	 
	 TGraphAsymmErrors *graph_1;
	 if(plot_systematics){
	   graph_1 = new TGraphAsymmErrors(18,x_val[iCentBin],y_val[iCentBin],valx_error[iCentBin],valx_error[iCentBin],sys_eMinus[iCentBin],sys_ePlus[iCentBin]);

	   graph_1->SetMarkerColor(colorR4[iCentBin]);
	   graph_1->SetLineColor(colorR4[iCentBin]);
	   graph_1->SetMarkerStyle(33);
	   graph_1->SetMarkerSize(1.5);
	   graph_1->SetFillColorAlpha(colorR4[iCentBin],0.5);
	 }
	  
	 
       	 
	 temp_jetRate->SetMarkerColor(colorR4[iCentBin]);
	 temp_jetRate->SetLineColor(colorR4[iCentBin]);
	 temp_jetRate->SetMarkerStyle(33);
	 temp_jetRate->SetMarkerSize(1.5);
	 //Offset
	 temp_jetRate->Scale(offSet[iCentBin]);
	if(iCentBin==0){
	  temp_jetRate->SetTitle("");
	  temp_jetRate->GetXaxis()->SetTitle("p_{T} [GeV]");
	  if(jetRatePbPb){
	    temp_jetRate->GetYaxis()->SetTitle("#frac{1}{<T_{AA}>}#frac{1}{N_{evt}}#frac{d^{2}N_{jet}}{dp_{T}d#eta}");
	    temp_jetRate->SetMaximum(1e+17);
	    temp_jetRate->SetMinimum(1e-13);
	  }
	  if(jetRatepp){
	    temp_jetRate->GetYaxis()->SetTitle("#frac{d^{2}#sigma}{dp_{T}d#eta}");
	    temp_jetRate->SetMaximum(1e-1);
	    temp_jetRate->SetMinimum(4e-7);
	  }
	  temp_jetRate->GetXaxis()->SetMoreLogLabels();
	  temp_jetRate->GetXaxis()->SetTitleOffset(1.4);
	  temp_jetRate->Draw();
	  
	  temp_jetRate->GetXaxis()->SetRangeUser(xMin,xMax);
	 
	}else{
	  
	  temp_jetRate->Draw("same");

	}

	if(plot_systematics)graph_1->Draw("same pe2");
	
	if(jetRatePbPb && plot_systematics){
	  leg->AddEntry(graph_1,Form("%s %s",centBinmap_2015Meas[iCentBin].c_str(),offsetTag[iCentBin].c_str()),"pelf");

	}else{
	  leg->AddEntry(temp_jetRate,Form("%s %s",centBinmap_2015Meas[iCentBin].c_str(),offsetTag[iCentBin].c_str()),"pelf");
	}
	 if(jetRatepp && plot_systematics)leg->AddEntry(graph_1,Form("Using p_{T} binning of Centrality Bin: %s %s",centBinmap_2015Meas[iCentBin].c_str(),offsetTag[iCentBin].c_str()),"pelf");
	
	if((jetRatePbPb && iCentBin==7)|| jetRatepp){
	  leg->SetBorderSize(0);
	  leg->Draw("same");
	  if(jetRatePbPb)Text_Info("pbpbdata",-1, jetRadius[jetR],0.6048144,0.6618076,0.8044132,0.8236152,"",0.03,32,eta_range);
	  if(jetRatepp)Text_Info("ppdata",-1, jetRadius[jetR],0.663992,0.7040816,0.8635908,0.8658892,"",0.03,32,eta_range);
	}
      }


      
      

      
      if(rAAPlot){
	TH1D *rAA = NULL;
      
	cout <<	__LINE__ << endl;
	if(!UnfoldedOrNot)rAA = (TH1D*)pTDis_PbPbData[jetR][iCentBin]->Clone(Form("Tep_Copy_%d",iCentBin));
	cout << __LINE__ << endl;
	if(UnfoldedOrNot){
	  cout << __LINE__ << endl;
	  rAA = (TH1D*)pTDis_Unfo_PbPbData[jetR][iCentBin][nominalInfo_PbPb[iCentBin][nominalNum]]->Clone();
	cout << __LINE__ << endl;
	}
         cout << __LINE__ << endl;  
	rAA->GetXaxis()->SetRangeUser(xMin,xMax);
	rAA->Divide(pTDis_ppData[jetR][iCentBin]);
	 cout << __LINE__ << endl;
	// Create a new histogram with the same binning as rAA, but without its contents
	TH1D *temp_rAA = (TH1D*)rAA->Clone("temp_rAA");
	temp_rAA->Reset();
	 cout << __LINE__ << endl;
	//These arrays will store the systematical uncer. for rAA values
	double rAAVal_y[centbins2015][20] = {};
	double rAAVal_x[centbins2015][20] = {};
	//Systematic Errors
	double sys_ePlus[centbins2015][20] = {};
	double sys_eMinus[centbins2015][20] = {};
	double valx_error[centbins2015][20] = {};
	//Statistical Errors
	double stat_err[centbins2015][20] = {};

	//These arrays belong to Dhanush's Analysis
	double sys_rAAVals[centbins2015][20] = {};
	double rAAVals_2018[centbins2015][20] = {};

	cout << __LINE__ << endl;
	// Set the bin contents of temp_rAA based on rAA within the fiducial region
	cout << "****************Start bin: " << nominalInfo_PbPb[iCentBin][startBin] << endl;
	cout << "*****************End bin: " << nominalInfo_PbPb[iCentBin][endBin] << endl;
	for (int iBin = nominalInfo_PbPb[iCentBin][startBin]+1; iBin <= nominalInfo_PbPb[iCentBin][endBin]; iBin++) {
	  cout << __LINE__ << endl;
	  temp_rAA->SetBinContent(iBin, rAA->GetBinContent(iBin));
	  temp_rAA->SetBinError(iBin, rAA->GetBinError(iBin));

       
	  rAAVal_x[iCentBin][iBin] = rAA->GetBinCenter(iBin);
	  rAAVal_y[iCentBin][iBin] = rAA->GetBinContent(iBin);
	  if(plot_systematics){
	    //going to setup systematics
	    TH1D *sys_hist = (TH1D*)systematicsFile->Get(Form("TotSysUncert_R%d_Cent_%s",jetRadius[jetR],centBins_2015Meas[iCentBin].c_str()));
	    sys_ePlus[iCentBin][iBin] = (sys_hist->GetBinContent(iBin)*rAA->GetBinContent(iBin))/100;
	    sys_eMinus[iCentBin][iBin] = (sys_hist->GetBinContent(iBin)*rAA->GetBinContent(iBin))/100;
	    valx_error[iCentBin][iBin] = abs(sys_hist->GetBinLowEdge(iBin)-sys_hist->GetBinCenter(iBin));
	  }
	}
	

	TGraphAsymmErrors *graph;
	if(plot_systematics){
	  graph= new TGraphAsymmErrors(18,rAAVal_x[iCentBin],rAAVal_y[iCentBin],valx_error[iCentBin],valx_error[iCentBin],sys_eMinus[iCentBin],sys_ePlus[iCentBin]);

	  graph->SetMarkerStyle(20);
     
	  graph->SetMarkerColor(colorR4[iCentBin]);
	  graph->SetLineColor(colorR4[iCentBin]);
	  graph->SetFillColorAlpha(colorR4[iCentBin],0.5);
	}
	cout << __LINE__ << endl;
	temp_rAA->GetXaxis()->SetTitle("p_{T} [GeV]");
	temp_rAA->GetYaxis()->SetTitle("R_{AA}");
	temp_rAA->SetTitle("");
	cout << __LINE__ << endl;
	temp_rAA->SetMarkerStyle(20);
	
	//2015 rAA Values
      
      
      
	if(debg)cout << __LINE__ << endl;
      
      
	temp_rAA->SetMarkerColor(colorR4[iCentBin]);
	temp_rAA->SetLineColor(colorR4[iCentBin]);

  
	if(debg)cout << __LINE__ << endl;
	if((iCentBin==0 || iCentBin==1) && rAAPlot){
	  temp_rAA->SetMaximum(1.5);
	  temp_rAA->SetMinimum(0);
      
	  temp_rAA->Draw("pe");
	  temp_rAA->GetXaxis()->SetRangeUser(xMin, xMax);
	}
	if(plot_systematics)graph->Draw("same pe2");
      
	if(iCentBin!=0){
	  temp_rAA->Draw("same pe");
	  
	}
	if(plot_systematics){
	  leg->AddEntry(graph,Form("%s",centBinmap_2015Meas[iCentBin].c_str()),"pefl");
	}else{
	  leg->AddEntry(temp_rAA,Form("%s",centBinmap_2015Meas[iCentBin].c_str()),"pefl");
	}
	rAA = NULL;
      
      

      leg->SetBorderSize(0);
      
    
      leg->Draw("same");
      
      if(rAAPlot){
	TLine *line1 = new TLine(xMin,1,xMax,1);
	line1->SetLineStyle(9);
	line1->Draw("same");
      }
    
      Text_Info("ppPbPbdata",-1, jetRadius[jetR],0.110331,0.6749271,0.3099298,0.8367347,"",0.03,11,eta_range);
      }    
    
  
  }





  //Save PDF
    if(savepdf){
      canv->SaveAs(Form("jetRates_R%d_%s.pdf",jetRadius[jetR],collSyst.c_str()));
    }
  



  
  
  return 0;


    }
