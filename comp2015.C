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


int comp2015(bool rAAPlot = false,string collsnType = "pp",bool jetRatePbPb=true, bool centBins2015 = true ,int jetR=R4, int centBin = -1, bool comp2015 = true,int xMin_value=158,double xMax_value=1000,bool debug = true){
  
  
  int numCentBins = 8;
  double etaRange = 0; //Using it to scale jet rate histo
  double eta_range = 0; //Using it to input in Text function

  string plotting_today = "rAA";

  string collType_Tag = "pp";

  if(jetRatePbPb==true){
    plotting_today = "jetRate";
   
  }

  if(centBin==-1){
    centBin=0;//this means you are plotting jet rates of pp
  }else if(centBin!=-1){
    collType_Tag = "PbPb";
  }
  
 
  
  
  
  if(comp2015){
    plotting_today = plotting_today+ "_comp2015Meas";
    etaRange = 2.8*2;
    eta_range = 2.8;
  }


  //y-Axis of plot
  double yMin = 0;
  double yMax = 0;
  //TLegend corners
  double pt1 =0;
  double pt2 =0;
  double pt3 =0;
  double pt4 =0;
  
  double offSet[] = {1,1e+2,1e+4,1e+6,1e+8,1e+10,1e+12,1e+14};
  string offsetTag[] {"","#times 10^{2}","#times 10^{4}","#times 10^{6}","#times 10^{8}","#times 10^{10}","#times 10^{12}","#times 10^{14}"};

  
  
  //Uploading files
  //---------Pb+Pb
  TFile * data_File = NULL;
  TFile * data_UnfFile = NULL;

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
  if(debug)cout << __LINE__ << endl;

  //Nominal Iter Info and Fudicial Region Info
  TFile *nom_FGInfoFile = NULL;
  
  int jetR_FR[8][2] = {}; //start and end bins
  int jetR_NominalIter[8] = {}; 
  
  string location = "/Users/berenicegarcia/Desktop/Berenice/Spring_2022/LargeRJets/LargeRJetNewCode/Locally/DiagPlots/";
  string cent2015BinTag = "";
  if(comp2015)cent2015BinTag="_JetRate2015Bins";
  cout << "Using this file to only plot fuducial region: " << Form("%snominalIter_%s%s_withpTShapeWeights_etaRange28.root",location.c_str(),collType_Tag.c_str(),cent2015BinTag.c_str()) << endl;
  
  nom_FGInfoFile = new TFile(Form("%s/nominalIter_%s%s_withpTShapeWeights_etaRange28.root",location.c_str(),collType_Tag.c_str(),cent2015BinTag.c_str()),"READ");

  TFile *file_tot = new TFile("Systematics/2015pTBins_FirstTotSysUncert_pbpbdata_R4.root","READ");


  
  
  //if(comp2015)nom_FGInfoFile = new TFile(Form("%s/nominalIter_%s_withpTShapeWeights.root",location.c_str(),collType_Tag.c_str()),"READ");
  if(debug)cout	<< __LINE__ << endl;
  //Dhanush's 2018 RAA Values
  TH1D *rAAVal_2018[num2015MeasBins];
  TH1D *rAAVal_2018Sys[num2015MeasBins];

  TFile *file_2018RAA = new TFile("Dhanush_Files/2018PbPb_inclusive_jetraa.root","READ");
  if(debug)cout << __LINE__ << endl;


  
   string fileTag = "ppdata";
 
  
   if(collType_Tag=="PbPb")fileTag="pbpbdata";
   if(debug)cout << __LINE__ << endl;
   //Grabbing Nominal iteration Value & fudicial region bins
   cout << "We will grab this hist..,: " << Form("%s_Cent_%s_R%d_%sBinNum",fileTag.c_str(),centBins_2015Meas[centBin].c_str(),jetRadius[jetR],startAndendTag[startAndend::start].c_str()) << endl;
   jetR_FR[centBin][startAndend::start] = ((TH1D*)nom_FGInfoFile->Get(Form("%s_Cent_%s_R%d_%sBinNum",fileTag.c_str(),centBins_2015Meas[centBin].c_str(),jetRadius[jetR],startAndendTag[startAndend::start].c_str())))->GetBinContent(1);
   if(debug)cout << __LINE__ << endl;
   jetR_FR[centBin][startAndend::end] = ((TH1D*)nom_FGInfoFile->Get(Form("%s_Cent_%s_R%d_%sBinNum",fileTag.c_str(),centBins_2015Meas[centBin].c_str(),jetRadius[jetR],startAndendTag[startAndend::end].c_str())))->GetBinContent(1);
   if(debug)cout << __LINE__ << endl;
   jetR_NominalIter[centBin] = ((TH1D*)nom_FGInfoFile->Get(Form("%s_Cent_%s_R%d_IterNum",fileTag.c_str(),centBins_2015Meas[centBin].c_str(),jetRadius[jetR])))->GetBinContent(1);
   if(debug)cout << __LINE__ << endl;

   //Grabbing the 2015 Measurement File and histogram
   int val_plot =3;
   if(collsnType=="PbPb")val_plot = centBin;
   if(debug)cout << __LINE__ << endl;
   if(jetRatePbPb){
    //Jet Rate for Pb+Pb 2015
    if(debug)cout << __LINE__ << endl;
    cout << "******YOU GRABBED THIS STUPID FILE TO COMPARE: " << Form("hep_2015Root/HEPData-ins1673184-v1-Table_%d.root",val_plot+1) << endl;
    TFile *file_2015 = new TFile(Form("hep_2015Root/HEPData-ins1673184-v1-Table_%d.root",val_plot+1),"READ");
    TDirectoryFile * dir = (TDirectoryFile*)file_2015->Get(Form("Table %d",val_plot+1));
    meas2015_graph[centBin] = (TGraph*) dir->Get("Graph1D_y1");
    meas2015_hist[centBin] = (TH1D*) dir->Get("Hist1D_y1");
    
   if(debug)cout << __LINE__ << endl;
    //Systematic Erros
    sys2015Errors_l[centBin] = (TH1D*) dir->Get("Hist1D_y1_e1minus");
    sys2015Errors_h[centBin] = (TH1D*) dir->Get("Hist1D_y1_e1plus");
    stat2015Errs[centBin] = (TH1D*) dir->Get("Hist1D_y1_e2");

   }else if(rAAPlot){
      //RAA Values from the 2015 Measurements
      TFile *file_rAA2015 = new TFile(Form("hep_2015Root/rAA2015/HEPData-ins1673184-v1-Table_%d.root",centBin+19),"READ");
      TDirectoryFile * direc = (TDirectoryFile*)file_rAA2015->Get(Form("Table %d",centBin+19));
      rAAVal_hist[centBin] = (TH1D*) direc->Get("Hist1D_y1");
      rAAVal_graph[centBin] = (TH1D*) direc->Get("Graph1D_y1");
      if(debug)cout	<< __LINE__ << endl;
      //RAA Sys Values from 2015
      rAAsys_ePlus[centBin] = (TH1D*) direc->Get("Hist1D_y1_e1plus");
      rAAsys_eMinus[centBin] = (TH1D*) direc->Get("Hist1D_y1_e1minus");

      //RAA Stat Errot
      rAA_statEr[centBin] = (TH1D*) direc->Get("Hist1D_y1_e2");

      //These are Dhanush's RAA Values from his analysis using 2018 Pb+Pb Data for R=0.4 jets
      rAAVal_2018Sys[centBin] = (TH1D*) file_2018RAA->Get(Form("h_jetpt_unfolded_raa_syst_PbPb (%s)",centBinmap_2015Meas[centBin].c_str())); 
      rAAVal_2018[centBin] = (TH1D*) file_2018RAA->Get(Form("h_jetpt_unfolded_raa_PbPb (%s)",centBinmap_2015Meas[centBin].c_str()));
      cout << "Name of histo: " << Form("h_jetpt_unfolded_raa_PbPb (%s)",centBinmap_2015Meas[centBin].c_str()) << endl;
    }
  
  
  
     if(debug)cout << __LINE__ << endl;
    
    //if(!centBins2015)data_UnfFile = new TFile(Form("",location.c_str()),"READ"); 
    //Using this to compare
     
    data_UnfFile = new TFile(Form("DiagPlots/rootFiles/2015CentBins/hist_26282548_05302022_Unfolded_%sData_ATLAS_Official_RAA_Binning_17Iters_10000Toys_JetRate2015Bins_etaRange28Nominal.root",collType_Tag.c_str()),"READ");

     //data_UnfFile = new TFile(Form("/Users/berenicegarcia/Desktop/debug/hist_26282548_05302022_Unfolded_ppData_Eta15_ATLAS_Official_RAA_Binning_3Iters_10000Toys_NewpTBins.root",collType_Tag.c_str()),"READ");


     
     if(debug)cout << __LINE__ << endl;
    yMin=3e-7;
    yMax=1e-1;

    pt1 =0.7514735;pt2 =0.2332016;pt3 =0.9823183;pt4 =0.5085639;




     if(debug)cout << __LINE__ << endl;

  
  
  
  //Print Out
  //cout << "This is the cent bin we will be plotting: " << centBin.c_str() << endl; 
  cout << "For radius: " << jetRadius[jetR] << endl;

  //Histograms
  //MC
  TH1D* pTDis_Truth[numCentBins];
  TH1D* pTDis_Reco[numCentBins];
  //Data
 
  //Unfolded
  TH1D *pTDis_Unfo_Data[numCentBins];
  
  
  

    
  //int nom_iter_num = jetR_NominalIter[centBin];int nom_iter_num = jetR_NominalIter[centBin];
  int nom_iter_num = 2;
     
  cout << "We grabbeds this histo: " << Form("Unfolded_%sData_R%d_%dIter_%s",collType_Tag.c_str(),jetRadius[jetR],nom_iter_num,centBins_2015Meas[centBin].c_str()) << endl;   
  pTDis_Unfo_Data[centBin] = (TH1D*)data_UnfFile->Get(Form("Unfolded_%sData_R%d_%dIter_%s",collType_Tag.c_str(),jetRadius[jetR],nom_iter_num,centBins_2015Meas[centBin].c_str()));
   if(debug)cout << __LINE__ << endl;
    
  TH1D *sysUncertTot[8];
  cout << __LINE__ << endl;
  cout << "GRABBING THIS: " << Form("TotSysUncert_R%d_Cent_%s",JetRadius[jetR],centBins_2015Meas[centBin].c_str()) << endl;
  sysUncertTot[centBin] = (TH1D*)file_tot->Get(Form("TotSysUncert_R%d_Cent_%s",JetRadius[jetR],centBins_2015Meas[centBin].c_str()));
  
  int numbins = 15;
  double totSysLargeR_yl[15] ={};
  double totSysLargeR_yh[15] ={};
  double totSysLargeR_x[15]={};
  double totSysLargeR_ex[15]={};
  double largeR_y_val[15]={};
  cout << __LINE__ << endl;
  cout << "This is how many bins you have: " << 15 << endl;
    
   if(debug)cout << __LINE__ << endl;

    
  //Normalization Factors
  double numEvents_centBin[] = {7.383e+8*LumNumPbPbData[jetR],7.383e+8*LumNumPbPbData[jetR],2*7.383e+8*LumNumPbPbData[jetR],2*7.383e+8*LumNumPbPbData[jetR],2*7.383e+8*LumNumPbPbData[jetR]};
  double numEventsMinBias = 7.383e+8*LumNumPbPbData[jetR];

  
  //Before you make the rAA plot
  const int centbins2015 = 8;


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
  
  if(rAAPlot){
     if(debug)cout << __LINE__ << endl;
    //Each centrality bins have different number of bins.
    //Storing the total number of bins here
    int totBins2015[8] = {16,12,15,13,15,13,15,13};
    double offset = 14;
      
    for(int iCentBin = 0; iCentBin < numCentBins; iCentBin++){
      
      for(int iBin=1; iBin < totBins2015[iCentBin]+1; iBin++){
	
	
	rAAVal_y[iCentBin][iBin] = rAAVal_hist[iCentBin]->GetBinContent(iBin);
	rAAVal_x[iCentBin][iBin] = rAAVal_hist[iCentBin]->GetBinCenter(iBin)+ offset;
      	
	valx_error[iCentBin][iBin] = abs(rAAVal_hist[iCentBin]->GetBinCenter(iBin)-rAAVal_hist[iCentBin]->GetBinLowEdge(iBin)); 
	
	sys_ePlus[iCentBin][iBin] = rAAsys_ePlus[iCentBin]->GetBinContent(iBin);
	sys_eMinus[iCentBin][iBin] = (-1)*rAAsys_eMinus[iCentBin]->GetBinContent(iBin);

	//Statistical errors
	stat_err[iCentBin][iBin] = rAA_statEr[iCentBin]->GetBinContent(iBin);
	//cout << "This is the stat error for your first bin: " << rAA_statEr[iCentBin]->GetBinContent(iBin) << endl;
      }//Looping over tot bins
    }//Looping over tot. centralities
  }
 
    
    
      TLegend* leg_JetRate = NULL;

      
      TCanvas * canv = new TCanvas(Form("%s%s_R%d_%s",fileTag.c_str(),plotting_today.c_str(),jetRadius[jetR],centBins_2015Meas[centBin].c_str()),Form("%s_R%d_%s",plotting_today.c_str(),jetRadius[jetR],centBins_2015Meas[centBin].c_str()),55,125,1030,884);
      
      leg_JetRate = NULL;
      leg_JetRate  = new TLegend(pt1,pt2,pt3,pt4,NULL,"brNDC");
      leg_JetRate->SetBorderSize(0);
      leg_JetRate->SetTextSize(0.03);
      
      gStyle->SetOptStat(0);
      canv->Range(2.071112,-17.56707,3.700387,25.96283);
      canv->SetFillColor(0);
      canv->SetBorderMode(0);
      canv->SetBorderSize(2);
      
      canv->SetLogx();
      canv->SetLogy();
      if(debug)cout << __LINE__ << endl;
      canv->SetLeftMargin(0.1168831);
      canv->SetRightMargin(0.2577423);
      canv->SetTopMargin(0.2154421);
      canv->SetBottomMargin(0.2017435);
      canv->SetFrameBorderMode(0);

      
      if(debug)cout << __LINE__ << endl; 
      //----------------------------------------FIRST PLOT PAD-----------------------------------------------//
      TPad *pad1 = new TPad(Form("pad_%d",centBin), Form("pad_%d",centBin), 0,0.3667055,0.9980545,0.9976717);
      pad1->Draw();
      pad1->cd();
      pad1->Range(2.067458,-6.407058,3.387166,-0.751512);
      pad1->SetFillColor(0);
      pad1->SetBorderMode(0);
      pad1->SetBorderSize(2);
      pad1->SetLogx();
      pad1->SetLogy();
      pad1->SetRightMargin(0.2933723);
      pad1->SetTopMargin(0.009225075);
      pad1->SetBottomMargin(0.08856085);
      pad1->SetFrameBorderMode(0);
      pad1->SetFrameBorderMode(0);
      gPad->SetTicks(1);
      
      cout << "This is the measurement 2018.... " << endl;
      cout << "Value of unfolded hist bin number 6: " << pTDis_Unfo_Data[centBin]->GetBinContent(6) << endl;
      
      
      pTDis_Unfo_Data[centBin]->Scale(1.,"width");
      cout << "Bin content but now we scale it by the bin width..." << endl;
      cout << "Value: " << pTDis_Unfo_Data[centBin]->GetBinContent(6) << endl;
      
      pTDis_Unfo_Data[centBin]->Scale(1/etaRange);
      cout << "We now scaled it by the eta region (2.8*2): " << pTDis_Unfo_Data[centBin]->GetBinContent(6) << endl;
      

  

     if(collsnType=="pp"){
	pTDis_Unfo_Data[centBin]->Scale(1/ppDataLumiVals[jetR]);
	cout << "scale it by the lumi" << endl;
	cout << "this is the value: " << pTDis_Unfo_Data[centBin]->GetBinContent(6) << endl;
      }else if(collsnType == "PbPb"){
	pTDis_Unfo_Data[centBin]->Scale(1/tAA_2015map[centBin]);
        pTDis_Unfo_Data[centBin]->Scale(1/numEventsMinBias);
      }
	
      //This is to compare systematics between 2015 and 2018 meas.
      double sysError2018_yh[10] ={};
      double sysError2018_yl[10] ={};
      double sysError2015_yl[10] = {};
      double sysError2015_yh[10] = {};
        
            
    double pTval[10] = {};
    double ratio[10] = {};
    //Statistical Erros
    double error_yl[10] = {};
    double error_yh[10]={};
    double error_x[10]={};

    
    //Systematic Errors
    double error_sys_yl[10] = {};
    double error_sys_yh[10] = {};
    //experimenting
    double error_Sys_ylnew[10] ={};
    double error_Sys_yhnew[10] ={};
    double yvalnew[10] = {};
    double xvalnew[10] = {};
    double errorxnew[10]={};

    
      
    int startBin_num = 0;
    bool passOnce = false;
      
    for(int iBin =2; iBin < 10; iBin++){
	  
	 while(!passOnce){
	    cout << "iBin: " << iBin << endl;
	    cout << "This it what is saved in  passOnce var. " << passOnce << endl;
	    cout << "This is what is saved in startBin_num." << startBin_num << endl;
	    
	    cout << "This is the value of the bin center in meas. 2015 w/ the bin number saved in iBin: " << meas2015_hist[centBin]->GetBinCenter(iBin+4) << endl;
	    cout << "This is the bin center of bin number number saved in startBin_num (my meas.): " << pTDis_Unfo_Data[centBin]->GetBinCenter(startBin_num) << endl;
	    if(meas2015_hist[centBin]->GetBinCenter(startBin_num) != pTDis_Unfo_Data[centBin]->GetBinCenter(iBin+4)){
	      startBin_num++;
	    }else{
	      passOnce =true;
	    }
	  }
	 cout << endl;
	 cout << "BIN CENTER OF 2015 MEAS:  " << meas2015_hist[centBin]->GetBinCenter(startBin_num) << endl;
	 cout << "BIN EDGE OF 2015 MEAS: " << meas2015_hist[centBin]->GetBinLowEdge(startBin_num) << endl;
	 cout << "BIN HIGH EDGE OF 2015 MEAS: " << (meas2015_hist[centBin]->GetBinCenter(startBin_num) + meas2015_hist[centBin]->GetBinWidth(startBin_num)/2) << endl;
	 cout << endl;
	 cout << "BIN CENTER OF 2018 MEAS: " << pTDis_Unfo_Data[centBin]->GetBinCenter(iBin+4) << endl;
	 cout << "BIN EDGE OF 2018 MEAS: " << pTDis_Unfo_Data[centBin]->GetBinLowEdge(iBin+4) << endl;
	 cout << "BIN HIGH EDGE OF 2018 MEAS: " << (pTDis_Unfo_Data[centBin]->GetBinCenter(iBin+4) + pTDis_Unfo_Data[centBin]->GetBinWidth(iBin+4)/2) << endl;
	 cout << endl;



	 

	 
	 //For the first plot
	totSysLargeR_yh[iBin] =  (sysUncertTot[centBin]->GetBinContent(iBin+4)*pTDis_Unfo_Data[centBin]->GetBinContent(iBin+4))/100;
	totSysLargeR_yl[iBin] =  ((sysUncertTot[centBin]->GetBinContent(iBin+4))*pTDis_Unfo_Data[centBin]->GetBinContent(iBin+4))/100;
	largeR_y_val[iBin] = pTDis_Unfo_Data[centBin]->GetBinContent(iBin+4);
	totSysLargeR_x[iBin] = pTDis_Unfo_Data[centBin]->GetBinCenter(iBin+4);
	totSysLargeR_ex[iBin] = abs(pTDis_Unfo_Data[centBin]->GetBinCenter(iBin+4)-pTDis_Unfo_Data[centBin]->GetBinLowEdge(iBin+4));
	


	//For the second plot 
	 error_Sys_ylnew[iBin] = (-1)*sys2015Errors_l[centBin]->GetBinContent(startBin_num);
	 error_Sys_yhnew[iBin] =sys2015Errors_h[centBin]->GetBinContent(startBin_num);
	 sysError2015_yl[iBin] = (-1)*sys2015Errors_l[centBin]->GetBinContent(startBin_num)/meas2015_hist[centBin]->GetBinContent(startBin_num);
	 sysError2015_yh[iBin] = sys2015Errors_h[centBin]->GetBinContent(startBin_num)/meas2015_hist[centBin]->GetBinContent(startBin_num);
	sysError2018_yh[iBin] = sysUncertTot[centBin]->GetBinContent(iBin+4)/100;
	sysError2018_yl[iBin] = sysUncertTot[centBin]->GetBinContent(iBin+4)/100;
 
	cout << "for bin center: " << pTDis_Unfo_Data[centBin]->GetBinCenter(iBin+4) << endl;
	cout << "This is our error centered at 1: " << sysError2018_yh[iBin] << endl;
	
	 yvalnew[iBin] = meas2015_hist[centBin]->GetBinContent(startBin_num);
	 xvalnew[iBin] = meas2015_hist[centBin]->GetBinCenter(startBin_num);
	 errorxnew[iBin]=abs(meas2015_hist[centBin]->GetBinCenter(startBin_num)-meas2015_hist[centBin]->GetBinLowEdge(startBin_num));
	 

	 
	 //For the third plot
	  pTval[iBin]=pTDis_Unfo_Data[centBin]->GetBinCenter(iBin+4);
	  ratio[iBin]=pTDis_Unfo_Data[centBin]->GetBinContent(iBin+4)/meas2015_hist[centBin]->GetBinContent(startBin_num);
	  cout << "THIS IS THE Y MEAS OF 2018: " << pTDis_Unfo_Data[centBin]->GetBinContent(iBin+4) << endl;
	  cout << "THIS IS THE Y MEAS OF 2015: " << meas2015_hist[centBin]->GetBinContent(startBin_num) << endl;
	  cout << "THIS IS THE BIN CENTER OF 2018: " << pTDis_Unfo_Data[centBin]->GetBinCenter(iBin+4) << endl;
	  cout << "THIS IS THE BIN CENTER OF 2015: " << meas2015_hist[centBin]->GetBinCenter(startBin_num) << endl;
	  cout << "IT SHOULD BE THE SAME" << endl;
	  cout << "THIS IS THE RATIO OF THE BETWEEN THE TWO: " << ratio[iBin] << endl;
	    
	  double up_er=pow(pTDis_Unfo_Data[centBin]->GetBinErrorUp(iBin+4)/pTDis_Unfo_Data[centBin]->GetBinContent(iBin+4),2) + pow(stat2015Errs[centBin]->GetBinContent(startBin_num)/meas2015_hist[centBin]->GetBinContent(startBin_num),2);

	  double low_er=pow(pTDis_Unfo_Data[centBin]->GetBinErrorLow(iBin+4)/pTDis_Unfo_Data[centBin]->GetBinContent(iBin+4),2) + pow(stat2015Errs[centBin]->GetBinContent(startBin_num)/meas2015_hist[centBin]->GetBinContent(startBin_num),2);

	   
	  error_yh[iBin] = ratio[iBin]*sqrt(up_er);
	  error_yl[iBin] = ratio[iBin]*sqrt(low_er);

	  
	  error_sys_yl[iBin] = ratio[iBin]*(-1)*sqrt(pow(sys2015Errors_l[centBin]->GetBinContent(startBin_num)/meas2015_hist[centBin]->GetBinContent(startBin_num),2)+pow(totSysLargeR_yl[iBin]/largeR_y_val[iBin],2));
	  error_sys_yh[iBin] = ratio[iBin]*sqrt(pow(sys2015Errors_h[centBin]->GetBinContent(startBin_num)/meas2015_hist[centBin]->GetBinContent(startBin_num),2)+pow(totSysLargeR_yh[iBin]/largeR_y_val[iBin],2));
	  
	  
	  error_x[iBin] = abs(pTDis_Unfo_Data[centBin]->GetBinCenter(iBin+4) - pTDis_Unfo_Data[centBin]->GetBinLowEdge(iBin+4));
	  startBin_num++;
	  
	}//Bin Loop







      
      
    auto gSys = new TGraphAsymmErrors(10,xvalnew,yvalnew,errorxnew,errorxnew,error_Sys_ylnew,error_Sys_yhnew);
    auto largeRgSys = new TGraphAsymmErrors(15,totSysLargeR_x,largeR_y_val,totSysLargeR_ex,totSysLargeR_ex,totSysLargeR_yl,totSysLargeR_yh);
      if(debug)cout << __LINE__	<< endl;
      gSys->SetFillStyle(3002);largeRgSys->SetFillStyle(3002);
      gSys->SetFillColor(kCyan+1);largeRgSys->SetFillColor(kYellow-2);
      gSys->SetMarkerStyle(20);largeRgSys->SetMarkerStyle(33);
      gSys->SetLineColor(kCyan+2);largeRgSys->SetLineColor(kYellow-2);
      gSys->SetMarkerColor(kCyan+2);largeRgSys->SetMarkerColor(kYellow-2);
      
      pTDis_Unfo_Data[centBin]->SetTitle("");
      pTDis_Unfo_Data[centBin]->GetXaxis()->SetTitle("p_{T} [GeV]");
      if(collsnType=="PbPb"){
	pTDis_Unfo_Data[centBin]->GetYaxis()->SetTitle("#frac{1}{<T_{AA}>}#frac{1}{N_{evt}}#frac{d^{2}N_{jet}}{dp_{T}d#eta}");
      }else if(collsnType=="pp"){
	pTDis_Unfo_Data[centBin]->GetYaxis()->SetTitle("#frac{d^{2}#sigma}{dp_{T}d#eta}");
      }
      
      pTDis_Unfo_Data[centBin]->GetXaxis()->CenterTitle(true);
      pTDis_Unfo_Data[centBin]->GetYaxis()->CenterTitle(true);
      
      
      pTDis_Unfo_Data[centBin]->GetXaxis()->SetMoreLogLabels();
      pTDis_Unfo_Data[centBin]->GetYaxis()->SetTitleSize(0.032);
      pTDis_Unfo_Data[centBin]->GetXaxis()->SetTitleSize(0.085);
      
      pTDis_Unfo_Data[centBin]->GetYaxis()->SetLabelSize(0.03);
      pTDis_Unfo_Data[centBin]->GetXaxis()->SetLabelSize(0.025);
	
      pTDis_Unfo_Data[centBin]->GetXaxis()->SetLabelSize(0);
      pTDis_Unfo_Data[centBin]->Draw();

      
      pTDis_Unfo_Data[centBin]->GetXaxis()->SetRangeUser(xMin_value,xMax_value);

      
      if(debug)cout << __LINE__	<< endl;
      if(comp2015){
	if(debug)cout << __LINE__	<< endl;
	meas2015_hist[centBin]->SetMarkerStyle(20);
	meas2015_hist[centBin]->SetMarkerColor(kCyan+2);
	meas2015_hist[centBin]->Draw("same HIST P");
      	if(debug)cout << __LINE__	<< endl;
	meas2015_graph[centBin]->SetLineColor(kCyan+2);
      	meas2015_graph[centBin]->Draw("P same");
	gSys->Draw("same pe2");largeRgSys->Draw("same pe2");
	if(debug)cout << __LINE__	<< endl;
      }

      if(debug)cout << __LINE__	<< endl;
      
      string tag_collision = "Pb+Pb";
      string tag_centBin = centBinmap_2015Meas[centBin].c_str();
      int yearDataTaken = 2018;
    
      if(collsnType=="pp"){
	tag_collision = "pp";
	tag_centBin = "";
	yearDataTaken = 2017;
      }
      
      
      leg_JetRate->AddEntry(largeRgSys,Form("%s %d Data %s",tag_collision.c_str(),yearDataTaken,tag_centBin.c_str()),"pealf");

      leg_JetRate->AddEntry(gSys,Form("%s 2015 Meas.%s",tag_collision.c_str(),tag_centBin.c_str()),"pealf");

	  
        
	
      pTDis_Unfo_Data[centBin]->SetLineColor(kYellow+3);
      pTDis_Unfo_Data[centBin]->SetMarkerColor(kYellow+3);
      pTDis_Unfo_Data[centBin]->SetMarkerStyle(33);
      pTDis_Unfo_Data[centBin]->SetMarkerSize(1.5);
      pTDis_Unfo_Data[centBin]->Draw("same");
	  
      
      

    leg_JetRate->Draw("same");
    Text_Info("ppdata",-1, 4,0.762279,0.6126482,0.9626719,0.7720685,"",0.038,11,eta_range);


    // Go back to the main canvas before defining pad2
    canv->cd();


    //---------------------------------Second Plot------------------------------------------------//
    
    TPad *pad = new TPad("", "",0.001945525,0.002328289,0.7140078,0.4225844);
   pad->Draw();
   pad->cd();
   pad->Range(2.070031,-1.245294,3.011576,1.515294);
   pad->SetFillColor(0);
   
   pad->SetBorderMode(0);
   pad->SetBorderSize(2);
   pad->SetLogx();
   pad->SetGridx();
   pad->SetLeftMargin(0.136612);
   pad->SetRightMargin(0.01229508);
   pad->SetTopMargin(0.005540166);
   pad->SetBottomMargin(0.5235457);
   pad->SetFrameBorderMode(0);
   pad->SetFrameBorderMode(0); 
   gPad->SetTicks(1);
   pad->SetTitle("");

   TLegend *legend_SysUncert  = new TLegend(0.7470817,0.3166473,0.9581712,0.3573923,NULL,"brNDC");
   legend_SysUncert->SetBorderSize(0);
   legend_SysUncert->SetTextSize(0.021);
   
   double array_ones[10]={1,1,1,1,1,1,1,1,1,1};
   //std::fill_n(array_ones, sizeof(array_ones), 1);
    
     
    numbins = 9;
    
    TLine *line_1 = new TLine(xMin_value,1,xMax_value,1);
    
    auto plot2018 = new TGraphAsymmErrors(10,pTval,array_ones,error_x,error_x,sysError2018_yh,sysError2018_yl);  
    auto plot2015 = new TGraphAsymmErrors(10,xvalnew,array_ones,errorxnew,errorxnew,sysError2015_yh,sysError2015_yl);
    
    
    

    plot2018->SetTitle("");
    plot2018->GetXaxis()->SetMoreLogLabels(); 
    plot2018->SetMaximum(1.2); 
    plot2018->SetMinimum(0.7); 
    plot2018->SetFillColorAlpha(kYellow-2, 0.5);
    
    
    plot2018->GetYaxis()->SetLabelSize(0.038);
   
    plot2018->GetYaxis()->SetTitleSize(0.079);
    plot2018->GetXaxis()->SetTitleSize(0.089);
    plot2018->GetYaxis()->SetTitleOffset(0.3);

    plot2018->Draw("ape2");
    plot2018->GetXaxis()->SetLimits(xMin_value,xMax_value);
    plot2015->SetFillColorAlpha(kCyan+2, 0.5);
   
    plot2015->Draw("same pe2");
    line_1->Draw("same");
    
    canv->cd();
    
    
    legend_SysUncert->AddEntry(plot2018,Form("%d Systematic Uncert.",yearDataTaken),"f");
    legend_SysUncert->AddEntry(plot2015,"2015 Systematic Uncert.","f");
    legend_SysUncert->Draw("same");

    //----------------------------------Third Plot---------------------------------------------//
    
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0.02, 1, 0.22);
    pad2->Draw();
    pad2->cd();
    pad2->Range(2.067336,-0.3211403,3.380548,0.000000061754);
    pad2->SetFillColor(0);
    pad2->SetBorderMode(0);
    pad2->SetBorderSize(2);
    pad2->SetGridx();
    pad2->SetRightMargin(0.2947471);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.2026775);
    pad2->SetFrameBorderMode(0);
    pad2->SetFrameBorderMode(0);
    pad2->SetLogx();
    pad2->SetTitle("");
    gPad->SetTicks(1);
    TLegend *leg_ratio  = NULL;

    if(comp2015)leg_ratio=new TLegend(0.731336,0.5061924,0.9621807,0.7644269,NULL,"brNDC");
    if(!comp2015)leg_ratio=new TLegend(0.7141454,0.4007905,0.9449902,0.8645586,NULL,"brNDC");
    leg_ratio->SetBorderSize(0);
     

     cout << __LINE__ << endl;
    
    for(int iIter = 1; iIter< 5; iIter++){
      TH1D * ratio = new TH1D(Form("unfold_NOI_%d",iIter),Form("unfold_NOI_%d",iIter),jetRateMCBins,jetRateMC);
      ratio = (TH1D*) pTDis_Unfo_Data[centBin]->Clone();
      ratio->Divide(meas2015_hist[centBin]);
      if(iIter==1){
	ratio->GetXaxis()->SetMoreLogLabels();
	ratio->SetMaximum(1.5);
	ratio->SetMinimum(0.2);
	
	ratio->GetYaxis()->SetLabelSize(0.08);
        ratio->GetXaxis()->SetLabelSize(0.089);
	ratio->GetYaxis()->SetTitleSize(0.089);
	ratio->GetXaxis()->SetTitleSize(0.089);
	ratio->GetYaxis()->SetTitleOffset(0.3);
	
	ratio->GetXaxis()->SetTitle("p_{T} [GeV]");
	ratio->GetYaxis()->SetTitle("Ratio");
	if(comp2015){
	  ratio->SetMarkerColor(kWhite);
	  ratio->SetFillColorAlpha(kBlue, 0.35);
	  ratio->Draw();
	  ratio->GetXaxis()->SetRangeUser(xMin_value,xMax_value);
	}
      }
      
    }

     
	auto g = new TGraphAsymmErrors(10,pTval,ratio,error_x,error_x,error_yl,error_yh);
	auto g_Sys = new TGraphAsymmErrors(10,pTval,ratio,error_x,error_x,error_sys_yl,error_sys_yl);
	
	g_Sys->SetFillStyle(3002);
	g_Sys->SetFillColor(kCyan+1);
	g_Sys->SetMarkerStyle(33);
        g_Sys->SetMarkerSize(1.5);
	g_Sys->SetLineColor(kCyan);
	g_Sys->SetMarkerColor(kCyan+2);
        g_Sys->SetLineColor(kCyan+2);

	
	g->SetMarkerStyle(33);
	g->SetMarkerSize(1.5);
	g->SetMarkerColor(kCyan+2);
	g->SetLineColor(kCyan+2);
	if(comp2015){
	  g->Draw("PE");
	  g->Draw("same PE");
	  g_Sys->Draw("same pe2");
	  leg_ratio->AddEntry(g_Sys,"2018 Meas./2015 Meas.","pealf");
	}
	leg_ratio->SetTextSize(0.099);
	leg_ratio->Draw("same");
        
	TLine *line_at1 = new TLine(xMin_value,1,xMax_value,1);
	line_at1->Draw("same");
 
    
    

      









  



  
  
  return 0;


}
