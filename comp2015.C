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


int comp2015(string plot_type = "RAA",double eta_range_val = 2.8, bool jetRate2015Bins = false,bool rAA2015Bins=true,int jetR=R4, int centBin = -1, bool compDhanush = false ,int xMin_value=158,double xMax_value=1000,bool debug = true){

  
  int numCentBins = 8;
  
  string plotting_today = "";

  string collType_Tag = "pp";
  
   string fileTag = "ppdata";
 
  if(jetRate2015Bins==true){
    plotting_today = "jetRate";
    plotting_today = plotting_today+ "_jetRate2015BinsMeas";
   
  }else if(rAA2015Bins==true){
    plotting_today = plotting_today+ "_2015RAARateBins";
  }
    
  if(centBin==-1){
    centBin=0;//this means you are plotting jet rates of pp
  }else if(centBin!=-1){
    collType_Tag = "PbPb";
  }
  

  
  if(collType_Tag=="PbPb" || plot_type == "RAA"){
    fileTag="pbpbdata";
    collType_Tag = "PbPb";
  }
  
  
  //Uploading files 
  TFile * data_File = NULL;
  TFile * data_UnfFile = NULL;

  //When runnning code to plot RAA this will be for the pp files
  //The above will be for Pb+Pb
  TFile * data_UnfFile_pp = NULL;

  //2015 comparison plots
  TFile *file_2015;
  TDirectoryFile * dir;

  //File that contains the total systematics
  cout << "LET US GRAB THIS: " << Form("Systematics/2015pTBins_FirstTotSysUncert_%s_R4.root",plot_type.c_str()) << endl;
  TFile *file_tot = new TFile(Form("Systematics/2015pTBins_FirstTotSysUncert_%s_R4.root",plot_type.c_str()),"READ");

  
  //Jet Rate in 2015 Pb+Pb Data
  TGraph *meas2015_graph[num2015MeasBins];
  TH1D *meas2015_hist[num2015MeasBins];
  //Systematica Erros
  TH1D *sys2015Errors_h[num2015MeasBins];
  TH1D *sys2015Errors_l[num2015MeasBins];
  //Statistical Erros
  TH1D *stat2015Errs[num2015MeasBins];

  //Nominal Iter Info and Fudicial Region Info
  TFile *nom_FGInfoFile = NULL;
  TFile *nom_FGInfoFile_pp = NULL;
  
  int jetR_FR[8][2] = {}; //start and end bins
  int jetR_NominalIter[8] = {}; 
  int jetR_NominalIter_pp[8] = {};
  
  string location = "/Users/berenicegarcia/Desktop/Berenice/Spring_2022/LargeRJets/LargeRJetNewCode/Locally/DiagPlots/";
  //this will be used when plotting RAA
  string location_rAAhists = "/Users/berenicegarcia/Desktop/Berenice/Spring_2022/LargeRJets/LargeRJetNewCode/Locally/DiagPlots/rootFiles/2015CentBins/rAA2015Bins";


  string pT2015BinsTag = "";

  
  if(jetRate2015Bins)pT2015BinsTag="_JetRate2015Bins";
  if(rAA2015Bins){
    pT2015BinsTag="_2015RAARateBins";
    location = location + "rootFiles/2015CentBins/rAA2015Bins/";
  }
  
  nom_FGInfoFile = new TFile(Form("%s/nominalIter_%s%s_withpTShapeWeights_etaRange%d.root",location.c_str(),collType_Tag.c_str(),pT2015BinsTag.c_str(),(int)(eta_range_val*10.0)),"READ");

  if(plot_type == "RAA"){
    nom_FGInfoFile_pp = new TFile(Form("%s/nominalIter_pp%s_withpTShapeWeights_etaRange%d.root",location_rAAhists.c_str(),pT2015BinsTag.c_str(),(int)(eta_range_val*10.0)),"READ");
  }

   
   ///////////////////////////////////////////////////////////
   //Grabbing Nominal iteration Value & fudicial region bins//
   ///////////////////////////////////////////////////////////
   cout << "We will grab this hist..: " << Form("%s_Cent_%s_R%d_%sBinNum",fileTag.c_str(),centBins_2015Meas[centBin].c_str(),jetRadius[jetR],startAndendTag[startAndend::start].c_str()) << endl;
   jetR_FR[centBin][startAndend::start] = ((TH1D*)nom_FGInfoFile->Get(Form("%s_Cent_%s_R%d_%sBinNum",fileTag.c_str(),centBins_2015Meas[centBin].c_str(),jetRadius[jetR],startAndendTag[startAndend::start].c_str())))->GetBinContent(1);
   if(debug)cout << __LINE__ << endl;
   jetR_FR[centBin][startAndend::end] = ((TH1D*)nom_FGInfoFile->Get(Form("%s_Cent_%s_R%d_%sBinNum",fileTag.c_str(),centBins_2015Meas[centBin].c_str(),jetRadius[jetR],startAndendTag[startAndend::end].c_str())))->GetBinContent(1);
    if(debug)cout << __LINE__ <<	endl;
   jetR_NominalIter[centBin] = ((TH1D*)nom_FGInfoFile->Get(Form("%s_Cent_%s_R%d_IterNum",fileTag.c_str(),centBins_2015Meas[centBin].c_str(),jetRadius[jetR])))->GetBinContent(1);
   if(debug)cout << __LINE__ <<	endl;
   if(plot_type == "RAA"){
     cout << "Trying to grab: " << Form("ppdata_Cent_%s_R%d_IterNum",centBins_2015Meas[centBin].c_str(),jetRadius[jetR]) << endl;
     jetR_NominalIter_pp[centBin] = ((TH1D*)nom_FGInfoFile_pp->Get(Form("ppdata_Cent_%s_R%d_IterNum",centBins_2015Meas[centBin].c_str(),jetRadius[jetR])))->GetBinContent(1);
      if(debug)cout << __LINE__ <<	endl;
   }

    if(debug)cout << __LINE__ <<	endl;
   //Grabbing the 2015 Measurement File and histogram
   int val_plot =3;
   if(plot_type=="PbPb")val_plot = centBin;

   
   
   if(plot_type == "PbPb" || plot_type == "pp"){
     
     //Jet Rate for Pb+Pb 2015
     file_2015 = new TFile(Form("hep_2015Root/HEPData-ins1673184-v1-Table_%d.root",val_plot+1),"READ");
     dir = (TDirectoryFile*)file_2015->Get(Form("Table %d",val_plot+1));
   

   }else if(plot_type == "RAA"){
     //RAA Values from the 2015 Measurements  
     file_2015 = new TFile(Form("hep_2015Root/rAA2015/HEPData-ins1673184-v1-Table_%d.root",centBin+19),"READ");
     dir = (TDirectoryFile*)file_2015->Get(Form("Table %d",centBin+19));
      
  
   }


     /* if(compDhanush){ */
	
     /* 	////////////////////////////// */
     /* 	//Dhanush's 2018 RAA Values/// */
     /* 	////////////////////////////// */
  
     /* 	TFile *file_2018RAA = new TFile("Dhanush_Files/2018PbPb_inclusive_jetraa.root","READ"); */
  
     /* 	//These are Dhanush's RAA Values from his analysis using 2018 Pb+Pb Data for R=0.4 jets */
     /* 	meas2015_hist[centBin] = (TH1D*) file_2018RAA->Get(Form("h_jetpt_unfolded_raa_PbPb (%s)",centBinmap_2015Meas[centBin].c_str())); */
     /*  } */

   

   
   
      meas2015_hist[centBin] = (TH1D*) dir->Get("Hist1D_y1");
      meas2015_graph[centBin] = (TGraph*) dir->Get("Graph1D_y1");
      
      //Sys Values from 2015
      sys2015Errors_h[centBin] = (TH1D*) dir->Get("Hist1D_y1_e1plus");
      sys2015Errors_l[centBin] = (TH1D*) dir->Get("Hist1D_y1_e1minus");

      // Stat Errot
      stat2015Errs[centBin] = (TH1D*) dir->Get("Hist1D_y1_e2");






   
   //////////////////////////////////////////////////
   ////Files that containt the nominal histograms////
   /////////////////////////////////////////////////
   
   data_UnfFile = new TFile(Form("%s/hist_26282548_05302022_Unfolded_%sData_ATLAS_Official_RAA_Binning_17Iters_10000Toys%s_etaRange%dNominal.root",location.c_str(),collType_Tag.c_str(),pT2015BinsTag.c_str(),(int)(eta_range_val*10.0)),"READ");
    
    if(plot_type == "RAA"){
      data_UnfFile_pp = new TFile(Form("%s/hist_26282548_05302022_Unfolded_ppData_ATLAS_Official_RAA_Binning_17Iters_10000Toys%s_etaRange%dNominal.root",location_rAAhists.c_str(),pT2015BinsTag.c_str(),(int)(eta_range_val*10.0)),"READ");
    }
    
     
    
  
  //Print Out
  //cout << "This is the cent bin we will be plotting: " << centBin.c_str() << endl; 
  cout << "For radius: " << jetRadius[jetR] << endl;

  //Data 
  //Unfolded
  TH1D *pTDis_Unfo_Data[numCentBins];
  TH1D *pTDis_Unfo_Data_pp[numCentBins];
  
  int nom_iter_num = jetR_NominalIter[centBin];
  
  cout << "We grabbeds this histo: " << Form("Unfolded_%sData_R%d_%dIter_%s",collType_Tag.c_str(),jetRadius[jetR],nom_iter_num,centBins_2015Meas[centBin].c_str()) << endl;   
  pTDis_Unfo_Data[centBin] = (TH1D*)data_UnfFile->Get(Form("Unfolded_%sData_R%d_%dIter_%s",collType_Tag.c_str(),jetRadius[jetR],nom_iter_num,centBins_2015Meas[centBin].c_str()));

  if(plot_type=="RAA"){
    cout << "Grab this histo: " << Form("Unfolded_ppData_R%d_%dIter_%s",jetRadius[jetR],nom_iter_num,centBins_2015Meas[centBin].c_str()) << endl;
    pTDis_Unfo_Data_pp[centBin] = (TH1D*)data_UnfFile_pp->Get(Form("Unfolded_ppData_R%d_%dIter_%s",jetRadius[jetR],nom_iter_num,centBins_2015Meas[centBin].c_str()));
  }
    
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
  
    
    
      TLegend* leg_JetRate = NULL;

       if(debug)cout << __LINE__ << endl;
      TCanvas * canv = new TCanvas(Form("%s%s_R%d_%s",fileTag.c_str(),plotting_today.c_str(),jetRadius[jetR],centBins_2015Meas[centBin].c_str()),Form("%s_R%d_%s",plotting_today.c_str(),jetRadius[jetR],centBins_2015Meas[centBin].c_str()),55,125,1030,884);
       if(debug)cout << __LINE__ << endl;
      leg_JetRate = NULL;
      leg_JetRate  = new TLegend(0.7514735,0.2332016,0.9823183,0.5085639,NULL,"brNDC");
      leg_JetRate->SetBorderSize(0);
      leg_JetRate->SetTextSize(0.03);
       if(debug)cout << __LINE__ << endl;
      gStyle->SetOptStat(0);
      canv->Range(2.071112,-17.56707,3.700387,25.96283);
      canv->SetFillColor(0);
      canv->SetBorderMode(0);
      canv->SetBorderSize(2);
       if(debug)cout << __LINE__ << endl;
      canv->SetLogx();
      canv->SetLogy();
      if(debug)cout << __LINE__ << endl;
      canv->SetLeftMargin(0.1168831);
      canv->SetRightMargin(0.2577423);
      canv->SetTopMargin(0.2154421);
      canv->SetBottomMargin(0.2017435);
      canv->SetFrameBorderMode(0);
       if(debug)cout << __LINE__ << endl;
      
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
      if(plot_type!="RAA")pad1->SetLogy();
      pad1->SetRightMargin(0.2933723);
      pad1->SetTopMargin(0.009225075);
      pad1->SetBottomMargin(0.08856085);
      pad1->SetFrameBorderMode(0);
      pad1->SetFrameBorderMode(0);
      gPad->SetTicks(1);
      

      double scale_factor = -1;
       if(debug)cout << __LINE__ << endl;
      
     if(plot_type=="pp"){
       scale_factor = ppDataLumiVals[jetR];
      }else if(plot_type == "PbPb" || plot_type == "RAA"){
        scale_factor = tAA_2015map[centBin]*numEventsMinBias; 
      }
      if(debug)cout << __LINE__ << endl;
     
     pTDis_Unfo_Data[centBin]->Scale(1.,"width");
     pTDis_Unfo_Data[centBin]->Scale(1/eta_range_val);
     pTDis_Unfo_Data[centBin]->Scale(1/scale_factor);
      if(debug)cout << __LINE__ << endl;
     if(plot_type=="RAA"){
       if(debug)cout << __LINE__ << endl;
       pTDis_Unfo_Data_pp[centBin]->Scale(1/ppDataLumiVals[jetR]);
       pTDis_Unfo_Data_pp[centBin]->Scale(1.,"width");
       pTDis_Unfo_Data_pp[centBin]->Scale(1/eta_range_val);
       pTDis_Unfo_Data[centBin]->Divide(pTDis_Unfo_Data_pp[centBin]);
       if(debug)cout << __LINE__ << endl;
     }
     
     
      //This is to compare systematics between 2015 and 2018 meas.
      double sysError2018_yh[16] ={};
      double sysError2018_yl[16] ={};
      double sysError2015_yl[16] = {};
      double sysError2015_yh[16] = {};
        
            
    double pTval[16] = {};
    double ratio[16] = {};
    //Statistical Erros
    double error_yl[16] = {};
    double error_yh[16]={};
    double error_x[16]={};

    
    //Systematic Errors
    double error_sys_yl[16] = {};
    double error_sys_yh[16] = {};
    //experimenting
    double error_Sys_ylnew[16] ={};
    double error_Sys_yhnew[16] ={};
    double yvalnew[16] = {};
    double xvalnew[16] = {};
    double errorxnew[16]={};

    
      
    int startBin_num = 0;
    bool passOnce = false;
    int numbinstot = 10;
    if(rAA2015Bins) numbinstot=16;
    for(int iBin =2; iBin < numbinstot; iBin++){
	  
	 while(!passOnce){
	    if(meas2015_hist[centBin]->GetBinCenter(startBin_num) != pTDis_Unfo_Data[centBin]->GetBinCenter(iBin+4)){
	      startBin_num++;
	    }else{
	      passOnce =true;
	    }
	  }

	 
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







      
      
    auto gSys = new TGraphAsymmErrors(numbinstot,xvalnew,yvalnew,errorxnew,errorxnew,error_Sys_ylnew,error_Sys_yhnew);
    auto largeRgSys = new TGraphAsymmErrors(15,totSysLargeR_x,largeR_y_val,totSysLargeR_ex,totSysLargeR_ex,totSysLargeR_yl,totSysLargeR_yh);
      if(debug)cout << __LINE__	<< endl;
      gSys->SetFillStyle(3002);largeRgSys->SetFillStyle(3002);
      gSys->SetFillColor(kCyan+1);largeRgSys->SetFillColor(kYellow-2);
      gSys->SetMarkerStyle(20);largeRgSys->SetMarkerStyle(33);
      gSys->SetLineColor(kCyan+2);largeRgSys->SetLineColor(kYellow-2);
      gSys->SetMarkerColor(kCyan+2);largeRgSys->SetMarkerColor(kYellow-2);
      
      pTDis_Unfo_Data[centBin]->SetTitle("");
      pTDis_Unfo_Data[centBin]->GetXaxis()->SetTitle("p_{T} [GeV]");
      if(plot_type=="PbPb"){
	pTDis_Unfo_Data[centBin]->GetYaxis()->SetTitle("#frac{1}{<T_{AA}>}#frac{1}{N_{evt}}#frac{d^{2}N_{jet}}{dp_{T}d#eta}");
      }else if(plot_type=="pp"){
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
      
      if(debug)cout << __LINE__	<< endl;
      meas2015_hist[centBin]->SetMarkerStyle(20);
      meas2015_hist[centBin]->SetMarkerColor(kCyan+2);
      meas2015_hist[centBin]->Draw("same HIST P");
      if(debug)cout << __LINE__	<< endl;
      meas2015_graph[centBin]->SetLineColor(kCyan+2);
      meas2015_graph[centBin]->Draw("P same");
      gSys->Draw("same pe2");largeRgSys->Draw("same pe2");
      if(debug)cout << __LINE__	<< endl;
      

      if(debug)cout << __LINE__	<< endl;
      
      string tag_collision = "Pb+Pb";
      string tag_centBin = centBinmap_2015Meas[centBin].c_str();
      int yearDataTaken = 2018;
    
      if(plot_type=="pp"){
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
    Text_Info(Form("%s",fileTag.c_str()),-1, 4,0.762279,0.6126482,0.9626719,0.7720685,"",0.038,11,eta_range_val);


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
    
    auto plot2018 = new TGraphAsymmErrors(18,pTval,array_ones,error_x,error_x,sysError2018_yh,sysError2018_yl);  
    auto plot2015 = new TGraphAsymmErrors(18,xvalnew,array_ones,errorxnew,errorxnew,sysError2015_yh,sysError2015_yl);
    
    
    

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
    int numbins_val = 5;
     
    if(jetRate2015Bins)leg_ratio=new TLegend(0.731336,0.5061924,0.9621807,0.7644269,NULL,"brNDC");
    if(!jetRate2015Bins){
      leg_ratio=new TLegend(0.7141454,0.4007905,0.9449902,0.8645586,NULL,"brNDC");
      numbins_val=9;
    }
     leg_ratio->SetBorderSize(0);
     

     cout << __LINE__ << endl;

     TH1D * ratio_hist;
     
    for(int iIter = 1; iIter< numbins_val; iIter++){
      
     if(jetRate2015Bins){
       ratio_hist=new TH1D(Form("unfold_NOI_%d",iIter),Form("unfold_NOI_%d",iIter),jetRateMCBins,jetRateMC);
     }
     if(rAA2015Bins){
       ratio_hist=new TH1D(Form("unfold_NOI_%d",iIter),Form("unfold_NOI_%d",iIter),bins2015,rAA_2015Bins[centBin]);
      
     }
      
      ratio_hist = (TH1D*) pTDis_Unfo_Data[centBin]->Clone();
      ratio_hist->Divide(meas2015_hist[centBin]);
      if(iIter==1){
	ratio_hist->GetXaxis()->SetMoreLogLabels();
	ratio_hist->SetMaximum(1.5);
	ratio_hist->SetMinimum(0.2);
	
	ratio_hist->GetYaxis()->SetLabelSize(0.08);
        ratio_hist->GetXaxis()->SetLabelSize(0.089);
	ratio_hist->GetYaxis()->SetTitleSize(0.089);
	ratio_hist->GetXaxis()->SetTitleSize(0.089);
	ratio_hist->GetYaxis()->SetTitleOffset(0.3);
	
	ratio_hist->GetXaxis()->SetTitle("p_{T} [GeV]");
	ratio_hist->GetYaxis()->SetTitle("Ratio_Hist");

	ratio_hist->SetMarkerColor(kWhite);
	ratio_hist->SetFillColorAlpha(kBlue, 0.35);
	ratio_hist->Draw();
	ratio_hist->GetXaxis()->SetRangeUser(xMin_value,xMax_value);
	
      }
      
    }

     
	auto g = new TGraphAsymmErrors(numbinstot,pTval,ratio,error_x,error_x,error_yl,error_yh);
	auto g_Sys = new TGraphAsymmErrors(numbinstot,pTval,ratio,error_x,error_x,error_sys_yl,error_sys_yl);
	
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
	
	g->Draw("PE");
	g->Draw("same PE");
	g_Sys->Draw("same pe2");
	leg_ratio->AddEntry(g_Sys,"2018 Meas./2015 Meas.","pealf");
	
	leg_ratio->SetTextSize(0.099);
	leg_ratio->Draw("same");
        
	TLine *line_at1 = new TLine(xMin_value,1,xMax_value,1);
	line_at1->Draw("same");
 
    
    

      









  



  
  
  return 0;


}
