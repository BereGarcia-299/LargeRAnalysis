#include <assert.h>
#include <cmath>
#include <string>
#include <fstream>
#include <stdio.h>
#include <math.h>
#include <tgmath.h>
#include <iostream>
#include <array>
#include <TChain.h>
#include <TLegend.h>
#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TPad.h"
#include "TMath.h"
#include "TF1.h"
#include "TRandom3.h"
#include <TAttMarker.h>
#include <TCanvas.h>
#include <TPaveText.h>
#include <dirent.h>
#include <sys/types.h>
#include "TLorentzVector.h"
#include "jer_jes_valuesR4.h"
#include "TLatex.h"


using namespace std;


float pi_dos = 360;
float pi_mi = 180;
float original = 3.14159;

float delta_phi(float phi_uno,float phi_dos){

  float final_phi = fmod(abs(phi_uno - phi_dos), 2*original);
  float phi = final_phi > original  ? 2*original - final_phi : final_phi;
  return phi;

}//-------delta phi function




float deltaR_calc(float eta_1, float eta_2,float phi_1, float phi_2 ){
  float deltaPhi = delta_phi(phi_1,phi_2);
  float deltaEta = eta_1 - eta_2;
  return sqrt(pow(deltaPhi,2) + pow(deltaEta,2));
}//------Calculates deltaR


void MakeTH1D_Eta_Array(TH1D *array_plot[], string name_plot, int bins, int pTSLices[],int total_pTslices){

  for(int ipTSlice = 0; ipTSlice < total_pTslices; ipTSlice++){
   array_plot[ipTSlice]=new TH1D(Form("%s_%i_%i",name_plot.c_str(),pTSLices[ipTSlice],pTSLices[ipTSlice+1]),Form("%s_%i_%i",name_plot.c_str(),pTSLices[ipTSlice],pTSLices[ipTSlice+1]),bins,-2.8,2.8);
  }//---Looping over pTSlices 
  
  
}

void  MakeTH2D_EtaPhi_Array(TH2D *array_plot[], string name_plot, int bins, int pTSLices[],int total_pTslices){
 
  for(int ipTSlice = 0; ipTSlice < total_pTslices; ipTSlice++){
    array_plot[ipTSlice]=new TH2D(Form("%s_%i_%i",name_plot.c_str(),pTSLices[ipTSlice],pTSLices[ipTSlice+1]),Form("%s_%i_%i",name_plot.c_str(),pTSLices[ipTSlice],pTSLices[ipTSlice+1]),bins,-2.8,2.8,bins,-3.14,3.14);
  }//---Looping over pTSlices

 
}


void Make_Pre_JER_Or_JES(TH1D* array_plots[],string name_plots, float  pTRanges[], int tot_ptranges,int bins,float x_low=0.5, float x_high=1.5){
  for(int ipTRange = 0; ipTRange < tot_ptranges; ipTRange++){
    array_plots[ipTRange] = new TH1D(Form("%s_%i_%iGeV",name_plots.c_str(),(int)pTRanges[ipTRange],(int)pTRanges[ipTRange+1]),Form("%s_%i_%iGeV",name_plots.c_str(),(int)pTRanges[ipTRange],(int)pTRanges[ipTRange+1]),bins,x_low,x_high);
  }
}


void write_TH1D_1DArray_to_file(TH1D* array[],int size ,string name_plot, int pTSlice[]){
  for(int index = 0; index < size; index++){
    array[index]->Write(Form("%s_%i_%i",name_plot.c_str(),pTSlice[index],pTSlice[index+1]),TObject::kOverwrite);
  }
}//------Writing out histograms to root file (For 1D Arrays)


void Text_Info(string sampleOrdata,int pT, int size_of_rad){
  bool nopT_info = false;
  if(pT == -1){
    nopT_info = true;
  }
    
  TPaveText *text_runs_numberjets = new TPaveText(0.6153846,0.7705263,0.8143813,0.8589474,"brNDC");
  text_runs_numberjets->SetTextSize(0.03);
  text_runs_numberjets->SetFillColor(0);
  text_runs_numberjets->SetBorderSize(0);
  text_runs_numberjets->SetShadowColor(0);
  text_runs_numberjets->SetTextAlign(33);
  text_runs_numberjets->SetTextFont(42);
  
  if(sampleOrdata == "ppmc"){
    text_runs_numberjets->AddText("pp PYTHIA");
  }else if(sampleOrdata == "pp data"){
    text_runs_numberjets->AddText("pp #sqrt{s}=5.02 TeV");
  }

  if(sampleOrdata == "pbpbmc"){
    text_runs_numberjets->AddText("PYTHIA + Pb+Pb Overlay");
  }
  if(sampleOrdata == "pbpbdata"){
    text_runs_numberjets->AddText("Pb+Pb 2018 Data");
  }

  if(!nopT_info){
    text_runs_numberjets->AddText(Form("p_{T} > %d GeV",pT));
  }
  text_runs_numberjets->AddText("|#eta| < 1.5");
  text_runs_numberjets->AddText(Form("R=0.%d",size_of_rad));
  text_runs_numberjets->Draw("same");
}

void Text_Info(string sampleOrdata,int pT, int size_of_rad, float v1, float v2, float v3, float v4, string centBin ="", float textSize = 0.056, int align = 13, float etaCut=1.5, bool transparent = false){
  bool nopT_info = false;
  if(pT == -1){
    nopT_info = true;
  }

  TPaveText *text_runs_numberjets = new TPaveText(v1,v2,v3,v4,"brNDC");
  text_runs_numberjets->SetTextSize(textSize);
  int ci =0;
  TColor *color;
  if(transparent){
    ci = 1181;
    color = new TColor(ci, 1, 1, 1, " ", 0);
  }
  
  text_runs_numberjets->SetFillColor(ci);
  if(transparent){
    ci = TColor::GetColor("#000000");
    text_runs_numberjets->SetLineColor(ci);
  }
  text_runs_numberjets->SetBorderSize(0);
  text_runs_numberjets->SetShadowColor(0);
  text_runs_numberjets->SetTextAlign(align);
  text_runs_numberjets->SetTextFont(42);
 
  text_runs_numberjets->AddText("#bf{#it{ATLAS}} Internal");

  if(sampleOrdata == "ppPbPbdata"){
    text_runs_numberjets->AddText("Pb+Pb 2018 Data");
    text_runs_numberjets->AddText("pp 2017 Data");
  }
  
  if(sampleOrdata == "ppmc"){
    text_runs_numberjets->AddText("pp PYTHIA");
  }else if(sampleOrdata == "ppdata"){
    text_runs_numberjets->AddText("pp #sqrt{s}=5.02 TeV");
  }

  

  if(sampleOrdata == "pbpbmc"){
    text_runs_numberjets->AddText("PYTHIA + Pb+Pb Overlay");
  }else if(sampleOrdata == "pbpbdata"){
    text_runs_numberjets->AddText("Pb+Pb 2018 Data");
  }

  
  if(!nopT_info){
    text_runs_numberjets->AddText(Form("p_{T} > %d GeV",pT));
  }
  if(etaCut!=-1)text_runs_numberjets->AddText(Form("|#eta| < %.1f",etaCut));
  if((sampleOrdata == "pbpbmc" || sampleOrdata == "pbpbdata") && (centBin != "")){
    text_runs_numberjets->AddText(Form("Centrality: %s",centBin.c_str()));
  }
  if(size_of_rad<10 && size_of_rad!=-1){
    text_runs_numberjets->AddText(Form("Anti-k_{T} R=0.%d",size_of_rad));
  }else if(size_of_rad!=-1){
    text_runs_numberjets->AddText("Anti-k_{T} R=1.0");
  }

  
  text_runs_numberjets->Draw("same");
}

void Text_AtlasInt(float v1, float v2, float v3,  float v4, int align=33){

  TPaveText *text_runs_numberjets = new TPaveText(v1,v2,v3,v4,"brNDC");
  text_runs_numberjets->SetTextSize(0.035);
  text_runs_numberjets->SetFillColor(0);
  text_runs_numberjets->SetBorderSize(0);
  text_runs_numberjets->SetTextFont(42);
  text_runs_numberjets->SetShadowColor(0);
  text_runs_numberjets->SetTextAlign(align);
  text_runs_numberjets->AddText("pp #sqrt{s_{NN}}=5.02 TeV");
  text_runs_numberjets->Draw("same");
}

void Text_AtlasIt(float xpos, float ypos,int align=33){
  TLatex *tex = new TLatex(xpos,ypos,"#bf{#it{ATLAS}} Internal");
  tex->SetTextFont(42);
  tex->SetTextSize(.04);
  tex->SetNDC();
  tex->SetTextColor(kBlack);
  tex->SetLineWidth(1);
  tex->SetTextAlign(align);
  tex->Draw("same");

}


void Text_AntikTJetInfo(float v1, float v2, float v3,  float v4, float radius, int align=33){

  TPaveText *text_runs_numberjets = new TPaveText(v1,v2,v3,v4,"brNDC");
  text_runs_numberjets->SetTextSize(0.035);
  text_runs_numberjets->SetFillColor(0);
  text_runs_numberjets->SetBorderSize(0);
  text_runs_numberjets->SetShadowColor(0);
  text_runs_numberjets->SetTextFont(42);
  text_runs_numberjets->SetTextAlign(align);
  text_runs_numberjets->AddText(Form("anti-k_{t} R=%.1f",radius));
  text_runs_numberjets->AddText("|#eta| < 1.5");
  text_runs_numberjets->Draw("same");
}

void drawText(const char *text, float xp, float yp, bool isRightAlign=0, int textColor=kBlack, double textSize=0.04, int textFont = 42, bool isNDC=true){
  TLatex *tex = new TLatex(xp,yp,text);
  tex->SetTextFont(textFont);
  tex->SetTextSize(textSize);
  tex->SetTextColor(textColor);
  tex->SetLineWidth(1);
  if(isNDC) tex->SetNDC();
  if(isRightAlign) tex->SetTextAlign(31);
  tex->Draw();
}



/////--------------------------------------------------------------------------------------------------------------/////
/////--------------------------------------------------------------------------------------------------------------/////
/////------------------------------------Dedicated to making JER and JES plots-------------------------------------/////
/////--------------------------------------------------------------------------------------------------------------/////
/////--------------------------------------------------------------------------------------------------------------/////

void MakeGaussianDis(string name_canvas,int numRows, int numColumns, TH1D* array[], float pTRanges[], TF1* tfArray[], string namefits, int JetRadius, bool ppMC = true, string centBin = ""){

  //Second Fit
  TF1 * fits_2ndFit[dj_totRanges];

  TCanvas *new_canvas= new TCanvas(Form("GausDis_%s_JetR%d",name_canvas.c_str(),JetRadius),Form("%s_JetR%d",name_canvas.c_str(),JetRadius),1500,1000);

  new_canvas->Divide(numColumns,numRows);

  for(int iPlot = 0; iPlot < dj_totRanges; iPlot++){
    new_canvas->cd(iPlot+1);
    
    array[iPlot]->SetMarkerStyle(20);

    if(JetRadius==4)array[iPlot]->SetMarkerColor(kBlue);
    if(JetRadius==10)array[iPlot]->SetMarkerColor(kRed);

    array[iPlot]->GetXaxis()->SetTitle("p^{Reco}_{T}/p^{Truth}_{T}");
    array[iPlot]->SetTitle(" ");
    
    array[iPlot]->Fit(Form("%s_%i_%iGeV",namefits.c_str(),(int)pTRanges[iPlot],(int)pTRanges[iPlot+1]),"M Q N","");
    
    float x_min = -tfArray[iPlot]->GetParameter(sigma_root)*1.5 + tfArray[iPlot]->GetParameter(mean_root);
    float x_max = tfArray[iPlot]->GetParameter(sigma_root)*1.5 + tfArray[iPlot]->GetParameter(mean_root);
    
    
    fits_2ndFit[iPlot] = new TF1(Form("fit_Jet_R%d_2ndFit_%i_%iGeV",JetRadius,(int)pTRanges_DJ[iPlot],(int)pTRanges_DJ[iPlot+1]),"gaus",x_min,x_max);
    
    array[iPlot]->Fit(Form("fit_Jet_R%d_2ndFit_%i_%iGeV",JetRadius,(int)pTRanges_DJ[iPlot],(int)pTRanges_DJ[iPlot+1]),"M Q N","",x_min,x_max);
    
    
    
    fits_2ndFit[iPlot]->SetLineStyle(9);
    fits_2ndFit[iPlot]->SetLineWidth(3);
    tfArray[iPlot]->SetLineStyle(1);
    tfArray[iPlot]->SetLineWidth(3);
   
    tfArray[iPlot]->SetLineColor(kBlack);
    if(JetRadius==4)fits_2ndFit[iPlot]->SetLineColor(kPink-9);
    if(JetRadius==10)fits_2ndFit[iPlot]->SetLineColor(kGreen+1);
    
    array[iPlot]->SetStats(0);
    array[iPlot]->Draw();
    tfArray[iPlot]->Draw("same");
   
    fits_2ndFit[iPlot]->Draw("same");
    if(iPlot == 0 && ppMC)Text_Info("ppmc",-1, JetRadius,0.1955066,0.5860293,0.3950689,0.8621471);
    if(iPlot == 0 && !ppMC)Text_Info("pbpbmc",-1,JetRadius,0.1298769,0.4949833,0.3301439,0.8963211,centBin.c_str());
    TPaveText *pTranges_text= new TPaveText(0.7952084,0.7714604,0.9917668,0.9721293,"brNDC");
    pTranges_text->SetTextSize(0.055);
    pTranges_text->SetTextFont(42);
    pTranges_text->SetFillColor(0);
    pTranges_text->SetBorderSize(0);
    pTranges_text->SetShadowColor(0);
    pTranges_text->SetTextAlign(33);

    if(centBin != "" && iPlot==0){
      pTranges_text->AddText(Form("%s",centBin.c_str()));
     }

    pTranges_text->AddText(Form("#bf{%d <  p^{Truth}_{T} < %d GeV}",(int)pTRanges[iPlot],(int)pTRanges[iPlot+1]));
    pTranges_text->Draw("same");
  }

}//-----Plots the Gaussian Distributions for jets that were matched (Truth and Redco jets) 



  void MakeGaussianDis(string name_canvas,int numRows, int numColumns, TH1D* array[], TH1D* extra_array[], float pTRanges[], TF1* arrayTF[], TF1* extra_TF[], string name_fits_of_array, string name_fits_of_extra_array,string label_for_array, string label_for_extra_array,string color, int size_of_jet){
    cout << __LINE__ << endl;
    
  TCanvas *new_canvas= new TCanvas(name_canvas.c_str(),name_canvas.c_str(),1000,1000);
  new_canvas->SetLogy();
  new_canvas->Divide(numColumns,numRows);
  cout << __LINE__ <<	endl;
  
  for(int iPlot = 0; iPlot < numRows*numColumns; iPlot++){
    new_canvas->cd(iPlot+1);
   cout << __LINE__ <<	endl;
    //-----Setting Marking Style
    array[iPlot]->SetMarkerStyle(21);
    extra_array[iPlot]->SetMarkerStyle(20);
    cout << __LINE__ <<	endl;
    //-----Color of Markers
    array[iPlot]->SetMarkerColor(kBlack);
    if(color == "red"){
     extra_array[iPlot]->SetMarkerColor(kRed);
     extra_array[iPlot]->SetLineColor(kRed);
    
    }else if(color=="blue"){
     extra_array[iPlot]->SetMarkerColor(kBlue); 
     extra_array[iPlot]->SetLineColor(kBlue);
    
    }
    cout << __LINE__ <<	endl;
    //-----Setting Colors
    array[iPlot]->SetLineColor(kBlack);
    cout << __LINE__ <<	endl;
    //-----Setting Axis Names and title
    array[iPlot]->GetXaxis()->SetTitle("p^{Reco}_{T}/p^{Truth}_{T}");
    array[iPlot]->SetTitle(" ");
    cout << __LINE__ <<	endl;
    array[iPlot]->Fit(Form("%s_%i_%iGeV",name_fits_of_array.c_str(),(int)pTRanges[iPlot],(int)pTRanges[iPlot+1]),"M Q N","",xMin_JetIsoCut[iPlot],xMax_JetIsoCut[iPlot]);
    
    
    extra_array[iPlot]->Fit(Form("%s_%i_%iGeV",name_fits_of_extra_array.c_str(),(int)pTRanges[iPlot],(int)pTRanges[iPlot+1]),"M Q N","",xMin_NoJetIsoCut[iPlot],xMax_NoJetIsoCut[iPlot]);
    
    cout << __LINE__ <<	endl;
    //---Setting Fitting Line Style a& Color
    arrayTF[iPlot]->SetLineStyle(2);
    arrayTF[iPlot]->SetLineColor(kBlack);

    extra_TF[iPlot]->SetLineStyle(2);
    if(color == "red"){
      extra_TF[iPlot]->SetLineColor(kRed);
    }else if(color == "blue"){
      extra_TF[iPlot]->SetLineColor(kBlue);
    }
    array[iPlot]->SetStats(0);

    array[iPlot]->Draw();
    cout << __LINE__ <<	endl;
    arrayTF[iPlot]->Draw("same");
    extra_TF[iPlot]->Draw("same");
    extra_array[iPlot]->Draw("same");
    cout << __LINE__ <<	endl;
    TPaveText *pTranges_text= new TPaveText(0.6583596,0.7849544,0.8576102,0.8837386,"brNDC");
    pTranges_text->SetTextSize(0.05);
    pTranges_text->SetFillColor(0);
    pTranges_text->SetBorderSize(0);
    pTranges_text->SetShadowColor(0);
    pTranges_text->SetTextAlign(21);
    pTranges_text->AddText(Form("%d < Truth  p_{T} < %d GeV",(int)pTRanges[iPlot],(int)pTRanges[iPlot+1]));
    pTranges_text->Draw("same");
    if(iPlot == 0){
      TLegend *ptRanges_legends= new TLegend(0.3411161,0.1609043,0.6929791,0.2606383,NULL,"brNDC");
      ptRanges_legends->AddEntry(array[iPlot],label_for_array.c_str(),"pe");
      ptRanges_legends->AddEntry(extra_array[iPlot],label_for_extra_array.c_str(),"pe");
      ptRanges_legends->SetBorderSize(0);
      ptRanges_legends->Draw("same");
      cout << "Only went through here once!" << endl;
      cout << __LINE__ <<	endl;
      Text_Info("ppmc",-1,size_of_jet);
       
    }
  }
}//-----Plots the Gaussian Distributions for jets that were matched (Truth and Redco jets) 




void JER_Or_JES(TH1D* gausHists[],string gaussHistNames[], int totCentBins, int totpTRanges, float pTRangesArray[],int jetRadius, int sigmaOrmean, float y_low=0.5, float y_high=1.5){

  const int totGauss = 120;
  const int totptRanges = 20;

  
  cout << __LINE__ << endl;
  float pt_bins[totptRanges] = {};
  float pT_bin_errors[totptRanges]= {};
  cout << __LINE__ << endl;
  float mid_pt_point =0;
  
  cout << __LINE__ << endl;
  for(int ipTBins=0; ipTBins < 20; ipTBins++){
    pT_bin_errors[ipTBins] = (pTRangesArray[ipTBins+1] - pTRangesArray[ipTBins])/2.0;

    mid_pt_point = ((pTRangesArray[ipTBins+1] - pTRangesArray[ipTBins])/2.0) + pTRangesArray[ipTBins];
    
    pt_bins[ipTBins] = mid_pt_point;
    cout << "pT mid point: " <<  mid_pt_point << endl;
  }//pT Bin Loop
  cout << __LINE__ << endl;
  
  //First Fit
  TF1 * fits_1stFit[totGauss];
  //Second Fit
  TF1 * fits_2ndFit[totGauss];
  
  int cent_Bin = 0;
  string centBin_Tag = centTimBins[cent_Bin];
  
  int ipTRange = 0;
  
  float sigma_or_mean_val[totptRanges]={};
  float sigma_or_mean_errors[totptRanges]={};

  string JER_or_JES = "";
  if(sigmaOrmean == 2)JER_or_JES="JER";
  if(sigmaOrmean == 1)JER_or_JES="JES";
  
  TCanvas *canv = new TCanvas(Form("%s_R%d_CentBins_PbPbMC",JER_or_JES.c_str(),jetRadius),Form("%s_R%d_CentBins_PbPbMC",JER_or_JES.c_str(),jetRadius),800,600);
  canv->SetLeftMargin(0.12);
  gPad->SetTicks(1);
  TLegend * jer_lg = new TLegend(0.3884712,0.1582609,0.8759398,0.293913,NULL,"brNDC");
  jer_lg->SetBorderSize(0);
  
  bool plot = false;

  bool restart = false;
  int counter = 0;
  
  
  //changed it from 120 -> 20
  for(int iGauss =0; iGauss < 120; iGauss++){
    //95
    cout << __LINE__ << endl;
    if(restart){
      
      
      
      //every 20 gausians that you go through
      //you enter a new centrality bin
      cout << "*******SETTING IPTRANGE = 0 && SIGMA ARRAYS TO ZERO" << endl;
     ipTRange = 0;
     float sigma_or_mean_val[totptRanges]={};
     float sigma_or_mean_errors[totptRanges]={};
     restart = false;
     
    }
    cout << __LINE__ << endl;
    cout << "iGauss: " << iGauss << endl;
    cout << "This is the pT Range :" << ipTRange << endl;
   cout << "Centr Bin: " << cent_Bin << endl;
    //Initilizing Fit Lines
    cout << __LINE__ <<	endl;
    
    fits_1stFit[iGauss] = new TF1(Form("fit_1st_Cent_%s_%sGeV",centBin_Tag.c_str(),pTRangeTag[ipTRange].c_str()),"gaus",0,2);
    cout << __LINE__ <<	endl;
    cout << __LINE__ << endl;
    cout << "where are we? " << centBin_Tag.c_str() << "/"<< pTRangeTag[ipTRange].c_str()<< endl;

    cout << __LINE__ <<	endl;
    cout << "iGauss: " << iGauss << endl; 
    gausHists[iGauss]->Fit(Form("fit_1st_Cent_%s_%sGeV",centBin_Tag.c_str(),pTRangeTag[ipTRange].c_str()),"M Q N","");
    cout << __LINE__ <<	endl;
     
    float x_min = -fits_1stFit[iGauss]->GetParameter(sigma_root)*1.5 + fits_1stFit[iGauss]->GetParameter(mean_root);  
    float x_max = fits_1stFit[iGauss]->GetParameter(sigma_root)*1.5 + fits_1stFit[iGauss]->GetParameter(mean_root);
    cout << __LINE__ << endl;
    cout << "sigma/mean of first fit: " << fits_1stFit[iGauss]->GetParameter(sigma_root) << "/" << fits_1stFit[iGauss]->GetParameter(mean_root) << endl;
    
    cout << __LINE__ << endl;
    fits_2ndFit[iGauss] = new TF1(Form("fit_2nd_Cent_%s_%sGeV",centBin_Tag.c_str(),pTRangeTag[ipTRange].c_str()),"gaus",x_min,x_max);

  cout << __LINE__ << endl;
    gausHists[iGauss]->Fit(Form("fit_2nd_Cent_%s_%sGeV",centBin_Tag.c_str(),pTRangeTag[ipTRange].c_str()),"M Q N","",x_min,x_max);
  
    //cout << "sigma/mean of second fit: " << fits_2ndFit[iGauss]->GetParameter(sigma_root) << "/" << fits_2ndFit[iGauss]->GetParameter(mean_root) << endl;
    cout << __LINE__ << endl;
    if(sigmaOrmean==sigma_root){
      sigma_or_mean_errors[ipTRange] = fits_2ndFit[iGauss]->GetParError(sigma_root);
      sigma_or_mean_val[ipTRange] = fits_2ndFit[iGauss]->GetParameter(sigma_root);
    }else{
      sigma_or_mean_errors[ipTRange] = fits_2ndFit[iGauss]->GetParError(mean_root);
      sigma_or_mean_val[ipTRange] = fits_2ndFit[iGauss]->GetParameter(mean_root);

    }
    
    ipTRange++;
    counter++;
  
    if((counter == 20)&&iGauss!=0){
        counter = 0;
      
      	
	
	auto* plot = new TGraphErrors(20,
                                             pt_bins,
                                             sigma_or_mean_val,
                                             pT_bin_errors,
				             sigma_or_mean_errors);

	//if(iGauss==19){
	if(iGauss==19){
	  plot->SetTitle("");
	
	  plot->GetXaxis()->SetTitle("p^{Truth}_{T} [GeV]");
	  if(sigmaOrmean == sigma_root)plot->GetYaxis()->SetTitle("#sigma");
	  if(sigmaOrmean == mean_root){
	    plot->GetYaxis()->SetTitle("#mu");
	   
	  }
	   plot->GetYaxis()->CenterTitle(true);
	  plot->Draw("A P");
	}else{
	  plot->Draw("P SAME");
	}

	plot->SetMarkerStyle(MarkerClosedStyles[cent_Bin]);
	plot->SetMarkerColor(MarkerColors[cent_Bin]);
	plot->SetLineColor(MarkerColors[cent_Bin]);
	//plot->SetMarkerColor(kBlack);
	//plot->SetMarkerStyle(21);
	//plot->SetLineColor(kBlack);
	//plot->SetMarkerSize(1.2);
	if(cent_Bin >2)plot->SetMarkerSize(1.5);
	//plot->SetName(Form("plot_%s",centTimBinsTags[cent_Bin].c_str()));

	if(sigmaOrmean==sigma_root){
	  plot->SetMaximum(y_high);
          plot->SetMinimum(y_low);
	}
	
	if(sigmaOrmean!=sigma_root){
	  plot->SetMaximum(y_high);
	  plot->SetMinimum(y_low);
	}
	if(iGauss < 100){
	  jer_lg->AddEntry(plot,centTimBinsTags[cent_Bin].c_str(),"pe");
	  
	  cent_Bin++;
	  centBin_Tag = centTimBins[cent_Bin];
	  cout << "**************PLOTTED THE FIRST ONE!*************" << endl;
	  restart=true;
	}else{
	  
	  jer_lg->AddEntry(plot,"pp","pe");
	  jer_lg->Draw("SAME");

	  
	}
    }
      
  
    
  }//Gaussian Loop
  /**
    TFile *DiJet_JER = new TFile("Tim_JER_JES_ppPYTHIA/histograms_2018.root");
    string text = "h_jer";
    if(sigmaOrmean==1)text = "h_jes";
    
    
    TH1D* jer_pp = (TH1D*) DiJet_JER->Get(Form("%s_4_5",text.c_str()));
  
    jer_pp->SetMarkerStyle(20);
    jer_pp->SetMarkerColor(kMagenta+1);
    jer_pp->Draw("P same");
    jer_lg->AddEntry(jer_pp,"From DJ Analysis (60-80%)","pe");
    jer_lg->Draw("same");
  */
    if(sigmaOrmean==1){
      TLine *l = new TLine(0,1,1100,1);
      l->SetLineStyle(2);
      l->SetLineWidth(1.3);
      l->Draw("same");
    }
  Text_Info("pbpbmc",-1, jetRadius,0.6904762,0.7269565,0.8909774,0.8869565,"",0.03,33);
  
}





void MAKE_JER_JES_W_JetIsoCut(TH1D* hist_array1[], TH1D* hist_array2[], TF1* fits_Array1[], TF1* fits_Array2[],string name_fits_of_array, string name_fits_of_extra_array,float pTRanges[], int size_of_jet, bool include_jerjes_values = false, bool include_isoJetCutPts = true){

   //Second Fit
  TF1 * fits_2ndFit[dj_totRanges];
  
  const int tot_fit_arrays = 2; 
  const int tot_fits = 20;
  const int tot_pTRanges = 20;
  const int num_of_parameters = 2; //use 0 index to grab sigmas and 1 for mean values
  
  float param_values[num_of_parameters][tot_fit_arrays][tot_fits]= {}; //This holds the sigma and mean for the first fit array
  float param_errors[num_of_parameters][tot_fit_arrays][tot_fits]={}; //This hols the sigma error values
  
  
  //--putting fits into one array
   
  for(int iArray= 0; iArray < tot_fit_arrays; iArray++){
    for(int ipTRange = 0; ipTRange < dj_totRanges;ipTRange++){
      
        if(iArray == without_IsoJetCt){
	  
	  hist_array1[ipTRange]->Fit(Form("%s_%i_%iGeV",name_fits_of_array.c_str(),(int)pTRanges[ipTRange],(int)pTRanges[ipTRange+1]),
				     "M Q N","");
	   float x_min = -fits_Array1[ipTRange]->GetParameter(sigma_root)*1.5 + fits_Array1[ipTRange]->GetParameter(mean_root);
           float x_max = fits_Array1[ipTRange]->GetParameter(sigma_root)*1.5 + fits_Array1[ipTRange]->GetParameter(mean_root);


	   fits_2ndFit[ipTRange] = new TF1(Form("fit_Jet_R4_2ndFit_%i_%iGeV",(int)pTRanges_DJ[ipTRange],(int)pTRanges_DJ[ipTRange+1]),"gaus",x_min,x_max);
	   hist_array1[ipTRange]->Fit(Form("fit_Jet_R4_2ndFit_%i_%iGeV",(int)pTRanges[ipTRange],(int)pTRanges[ipTRange+1]),"M Q N","",x_min,x_max);
	   
	   
	for(int iParam = 0; iParam < num_of_parameters; iParam++){         
	  if(iParam == sigma_index){
	    param_errors[iParam][iArray][ipTRange] = fits_2ndFit[ipTRange]->GetParError(sigma_root);
	    param_values[iParam][iArray][ipTRange] = fits_2ndFit[ipTRange]->GetParameter(sigma_root);
	  ;
	  }else if(iParam == mean_index){
	    
	    param_errors[iParam][iArray][ipTRange] = fits_2ndFit[ipTRange]->GetParError(mean_root);
	    param_values[iParam][iArray][ipTRange] = fits_2ndFit[ipTRange]->GetParameter(mean_root);
	    
	  }
        }//looping over total number of parameters
	}else if(iArray == with_IsoJetCt && include_isoJetCutPts){
	  //0.8,1.2
	  hist_array2[ipTRange]->Fit(Form("%s_%i_%iGeV",name_fits_of_extra_array.c_str(),(int)pTRanges[ipTRange],(int)pTRanges[ipTRange+1]),
				   "M Q N","",xMin_JetIsoCut[ipTRange],xMax_JetIsoCut[ipTRange]);

	  
	  for(int iParam = 0; iParam < num_of_parameters; iParam++){
	    if(iParam == sigma_index){
	      param_errors[iParam][iArray][ipTRange] = fits_Array2[ipTRange]->GetParError(sigma_root);
	      param_values[iParam][iArray][ipTRange] = fits_Array2[ipTRange]->GetParameter(sigma_root);
	        
	    }else if(iParam == mean_index){
	      param_errors[iParam][iArray][ipTRange] = fits_Array2[ipTRange]->GetParError(mean_root);
	      param_values[iParam][iArray][ipTRange] = fits_Array2[ipTRange]->GetParameter(mean_root);
	       
	   }
	  }//Looping over total number of parameters
	}

      
    }//looping over all of the pT bins
  }//looping over with or without Iso cut arrays
  
    
  
  //------Plotting
  const int tot_bins = 20;
 
  float pt_bins[tot_bins] = {};
  float pT_bin_errors[tot_bins]= {}; 
 
  float step = 0;   
  float mid_pt_point = 0;


    

    
  for(int iBins=0; iBins < dj_totRanges; iBins++){
    pT_bin_errors[iBins] = (pTRanges[iBins+1] - pTRanges[iBins])/2.0;
    
    mid_pt_point = ((pTRanges[iBins+1] - pTRanges[iBins])/2.0) + pTRanges[iBins];  
    
    pt_bins[iBins] = mid_pt_point; 
    
  }//looping over pt ranges to create an array that hols the pT mid points for the x-axis


  for(int iParam = 0; iParam < num_of_parameters; iParam++){
    TCanvas * jer_cav = new TCanvas(Form("%i_JetIso_R%i",iParam,size_of_jet),
				    Form("%i_JetIso_R%i",iParam,size_of_jet),600,500);
    TLegend * jer_lg = new TLegend(0.5610368,0.5663158,0.8670569,0.7031579,NULL,"brNDC");

    for(int NoIsoCut_Or_IsoCut=0 ;NoIsoCut_Or_IsoCut < 2; NoIsoCut_Or_IsoCut++){
    
      //-----W/ Iso cut = 0
      //----W/O Iso Cut = 1
    
      if(NoIsoCut_Or_IsoCut == 1){
	cout << "printing out all of the JER values! " << endl;
	for(int ipTRange =0; ipTRange < 20; ipTRange++ ){
	  cout << "sigma: " << param_values[0][1][ipTRange] << endl;
	  //cout << "sigma other:" << param_values[0][0][ipTRange] << endl;
	  //param_values[iParam][iArray][ipTRange] 
	}
      }
      auto* plot = new TGraphErrors(dj_totRanges,
					     pt_bins,
					     param_values[iParam][NoIsoCut_Or_IsoCut],
					     pT_bin_errors,
				             param_errors[iParam][NoIsoCut_Or_IsoCut]);



     

      
      plot->GetXaxis()->SetTitle("");

      if(iParam == sigma_index){
	plot->GetYaxis()->SetTitle("JER");
       if(size_of_jet == 4){
         plot->GetYaxis()->SetRangeUser(0.4,0.12);
        }else if(size_of_jet == 10){
         plot->GetYaxis()->SetRangeUser(0.4,0.12);
        }

      }else{
	plot->GetYaxis()->SetTitle("JES");
        if(size_of_jet == 4){
          plot->GetYaxis()->SetRangeUser(0.99,1.01);
        }else if(size_of_jet == 10){
	  plot->GetYaxis()->SetRangeUser(0.99,1.01);
        }
	
      }
      
       plot->GetXaxis()->SetTitle("p^{Truth}_{T}");
    
      if(NoIsoCut_Or_IsoCut == without_IsoJetCt){

	if(size_of_jet == 4){
         
         plot->SetMarkerColor(MyColorsJER_JES_R4[NoIsoCut_Or_IsoCut]);
         plot->SetLineColor(MyColorsJER_JES_R4[NoIsoCut_Or_IsoCut]);
	}else if(size_of_jet == 10){
         
         plot->SetMarkerColor(MyColorsJER_JES_R10[NoIsoCut_Or_IsoCut]);
         plot->SetLineColor(MyColorsJER_JES_R10[NoIsoCut_Or_IsoCut]);
        }
	plot->SetMarkerStyle(21);
	plot->SetMarkerSize(1.2);

	 
	 plot->SetTitle("");
	 if(include_isoJetCutPts)plot->Draw("P same");
	 if(!include_isoJetCutPts)plot->Draw("A P");
	 jer_lg->AddEntry(plot,"w/o Jet Isolation Cut");


	if(include_jerjes_values && size_of_jet == 4 && iParam == sigma_index){
	  //auto * add_plot = new TGraphErrors(dj_totRanges,
	  //				   JER_pT_bins_R4,
	  //				   JER_values_R4,
	  //				   error_x,
	  //				   error_y);
	  TFile *DiJet_JER = new TFile("Tim_JER_JES_ppPYTHIA/histograms_2018.root");

          TH1D* jer_pp = (TH1D*) DiJet_JER->Get("h_jer_5_5");
	  
	   jer_pp->SetMarkerStyle(20);
	   jer_pp->SetMarkerColor(kMagenta+1);
	   jer_pp->Draw("P same");
	   jer_lg->AddEntry(jer_pp,"From DJ Analysis");
	 }else if(include_jerjes_values && size_of_jet == 4 && iParam == mean_index){
	  //auto * add_plot = new TGraphErrors(dj_totRanges,
	  //                               JES_pT_bins_R4,
	  //                               JES_values_R4,
	  //                               error_x,
	  //                               error_y);


	   //Tim's JER & JES value for pp PYTHIA
          TFile *DiJet_JER_JES = new TFile("Tim_JER_JES_ppPYTHIA/histograms_2018.root");

	  TH1D* jes_pp = (TH1D*) DiJet_JER_JES->Get("h_jes_5_5");
	  
	  jes_pp->SetMarkerStyle(20);
          jes_pp->SetMarkerColor(kMagenta+1);
          jes_pp->Draw("P same");
          jer_lg->AddEntry(jes_pp,"From DJ Analysis");

	 }

	
      }else if(NoIsoCut_Or_IsoCut == with_IsoJetCt && include_isoJetCutPts){
	if(size_of_jet == 4){
         
         plot->SetMarkerColor(MyColorsJER_JES_R4[NoIsoCut_Or_IsoCut]);
         plot->SetLineColor(MyColorsJER_JES_R4[NoIsoCut_Or_IsoCut]);
        }else if(size_of_jet == 10){
	 
         plot->SetMarkerColor(MyColorsJER_JES_R10[NoIsoCut_Or_IsoCut]);
         plot->SetLineColor(MyColorsJER_JES_R10[NoIsoCut_Or_IsoCut]);
       	}
	plot->SetMarkerStyle(21);
	
	plot->Draw("A P");

	jer_lg->AddEntry(plot,"w/ Jet Isolation Cut");
      
      }
    
    }//looping over two plots (JER for w/ and w/o Jet Isolation cut)
  


   jer_lg->SetBorderSize(0);
   jer_lg->Draw("same");


   Text_Info("ppmc",-1, size_of_jet,0.6904762,0.7269565,0.8909774,0.8869565,"",0.03,33,2.1);
   //Text_AtlasInt(0.1438127,0.7747368,0.3428094,0.8631579);
   //Text_Info("ppmc",-1, size_of_jet);
  }//---looping over total number of parameters
}//looping over pT ranges









///--------------------------------------------------------------------------------------------------------////
///--------------------------------------------------------------------------------------------------------////
////----------------------------------END OF JES AND JER FUNCTIONS------------------------------------------////
///--------------------------------------------------------------------------------------------------------////
///--------------------------------------------------------------------------------------------------------////




void MeanSigmaText(double mean, double sigma,Color_t color, const int CutORNoCut){
 
  TPaveText *text_runs_numberjets = new TPaveText(CoordinatesMeanSigmaTxt[CutORNoCut][0],
      CoordinatesMeanSigmaTxt[CutORNoCut][1],
      CoordinatesMeanSigmaTxt[CutORNoCut][2],
      CoordinatesMeanSigmaTxt[CutORNoCut][3],"brNDC");

  text_runs_numberjets->SetTextSize(0.03);
  text_runs_numberjets->SetFillColor(0);
  text_runs_numberjets->SetBorderSize(0);
  text_runs_numberjets->SetShadowColor(0);
  text_runs_numberjets->SetTextAlign(21);
  text_runs_numberjets->SetTextColor(color);
  text_runs_numberjets->AddText(Form("%.2f +/- %.2f",mean,sigma));
  
  text_runs_numberjets->Draw("same");
}

void smallTexts(string StrArray[], int textLines, string labels[]){
  TPaveText *text_runs_numberjets = new TPaveText(0.1605547,0.2697424,0.3606664,0.3683493,"brNDC");
  text_runs_numberjets->SetTextSize(0.03);
  text_runs_numberjets->SetFillColor(0);
  text_runs_numberjets->SetBorderSize(0);
  text_runs_numberjets->SetShadowColor(0);
  text_runs_numberjets->SetTextAlign(21);
  
  for(int itext =0; itext < textLines; itext++){
    text_runs_numberjets->AddText(Form(" %s %s",labels[itext].c_str(),StrArray[itext].c_str()));
  }
  text_runs_numberjets->Draw("same");
}

void MakeTable_Latex(int num_columns, int num_rows, float ColArray[4][2], string ColumnNames[]){
  string column_name = "";
  string column_sect = "";

  cout << "\begin{center}" << endl;

  for(int iCol = 0; iCol< num_columns; iCol++){
    column_sect = " c " + column_sect; 
    column_name  = column_name + " & " + ColumnNames[iCol];
  }//---Column

  cout << "\\begin{tabular}{||" << column_sect << "||}" << endl;
  cout << " \\hline" << endl;
  cout << " " <<  column_name << " \\ [0.5ex]" <<  endl;
  cout << " \\hline\\hline" << endl;


  for(int iRow = 0; iRow < num_rows; iRow++){
    string raw_print = "";
    for(int iCol = 0; iCol < num_columns; iCol++){
      raw_print = (string)Form("%f",ColArray[iCol][iRow]) + " & " + raw_print;
    }//-Columns
    cout << " " + raw_print << " \\ " << endl;
    cout << " \\hline" << endl;
  }//-Rows 

  cout << "\\end{tabular}" << endl;

  cout << "\\end{center}" << endl;
 
}

//-------------------------------------------------------------//

// void Ratio_Plots(TH1D * param_Datadist1[],string nameDataDis[], TH1D* param_MCdis2[],float jetRadius ,string param,float y_low = 0 , float y_high = 2, bool sumpTCutData = true,bool mc_MatchOrNoMatch = false){
void Ratio_Plots(TH1D * param_Datadist1[], TH1D* param_MCdis2[],float jetRadius ,string param,float y_low = 0 , float y_high = 2, bool sumpTCutData = true,bool cleantool = true,bool sumpTANDcleanjets= false,bool mc_MatchOrNoMatch = false){


  int indexPlotNoMatchReq = 0;
  
  if(mc_MatchOrNoMatch == false){
    //This means we want to plot mc pt distributions w/ matching requirement to truth jets
    //Plot all of the Reco jets that passed pT cut and eta cut
    cout << "We will plot all of reco jets! No matching requirement..." << endl; 
    indexPlotNoMatchReq = 1;
  }
  

  int iJetIndex = -1;
  cout << "This is the jet radius: " << jetRadius << endl; 
  
  for(int iJetR =0; iJetR < Tot_Radii; iJetR++){
    cout << "RadiusSize[iJetR]: " << RadiusSize[iJetR] << endl;
    cout << "(float)RadiusSize[iJetR]: " << (float)RadiusSize[iJetR] << endl;
    cout << "jetRadius: " << jetRadius << endl;
    cout << "(float)jetRadius: " << (float)jetRadius << endl;
    cout << __LINE__ << endl;
    cout << "(float)RadiusSize[iJetR]==(float)jetRadius: " << ((float)RadiusSize[iJetR]==(float)jetRadius) << endl;
    if((float)RadiusSize[iJetR]==(float)jetRadius){
	 iJetIndex = iJetR;
	   break;
	   cout << __LINE__ << endl;
       }
   }//Jet Radius Loop
      

  cout << __LINE__ << endl;
  TCanvas *canv = new TCanvas(Form("%sDis_R%d",param.c_str(),(int)jetRadius), Form("%sDis_R%d",param.c_str(),(int)jetRadius), 800, 800);

  TLegend *ptRanges_legends= new TLegend(0.5350877,0.7308756,0.8859649,0.878341,NULL,"brNDC");
  
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
  pad1->SetBottomMargin(0); // Upper and lower plot are joined}
  if("pT" == param){
    pad1->SetLogy(1); 
    cout << __LINE__ << endl;
  }
  pad1->Draw();             // Draw the upper pad: pad
  pad1->cd();               // pad1 becomes the current pad

  cout << __LINE__ << endl;
  
  int index_data = 0;
  int index_mc = 0;
  int data_tot = 2;
  // If data_tot equatl to 2 then that means I want to plot the Data  pT distributions for:
  // (Data) w/o SumpT 
  // (Data) w/ SumpT
  // (MC) w/o SumpT
  
  if(!sumpTCutData && !cleantool && !sumpTANDcleanjets)data_tot = 1;
  //If tot equals to 1 then that means I only want to plot the following distributions:
  //(Data) w/o SumpT
  //(MC) w/o SumpT
  cout << __LINE__ << endl;
    
  param_Datadist1[index_data]->SetTitle("");
  param_Datadist1[index_data]->GetYaxis()->SetTitle("N/Lumi [nb]");
  param_Datadist1[index_data]->GetXaxis()->SetTitle("p_{T} [GeV]");
  param_Datadist1[index_data]->SetStats(0);
  param_Datadist1[index_data]->GetYaxis()->SetRangeUser(1e-9,8);
  cout << __LINE__ << endl;
  if(jetRadius == 1.0){
    //w/o SumpT Cut
    param_Datadist1[WOSumpT]->SetLineColor(kRed+2);
    param_Datadist1[WOSumpT]->SetMarkerColor(kRed+2);
    param_Datadist1[WOSumpT]->SetMarkerStyle(21);
    param_Datadist1[WOSumpT]->Scale(1/ppDataLumiVals[0]);
    cout << __LINE__ << endl;
    
    if(sumpTCutData || cleantool || sumpTANDcleanjets){
      cout << __LINE__ << endl;
      //w/ SumpT Cut
      param_Datadist1[WSumpT]->SetLineColor(kRed);
      param_Datadist1[WSumpT]->SetMarkerColor(kRed);
      param_Datadist1[WSumpT]->SetMarkerStyle(20);
      param_Datadist1[WSumpT]->Scale(1/ppDataLumiVals[0]);
      cout << __LINE__ << endl;
    }
    
  }else if(jetRadius == (float)0.4){
    cout << __LINE__ << endl;
    //w/o SumpT Cut
    param_Datadist1[WOSumpT]->SetMaximum(8);
    param_Datadist1[WOSumpT]->SetMinimum(1e-10);
    param_Datadist1[WOSumpT]->SetLineColor(kAzure-4);
    param_Datadist1[WOSumpT]->SetMarkerColor(kAzure-4);
    param_Datadist1[WOSumpT]->SetMarkerStyle(21);
    param_Datadist1[WOSumpT]->Scale(1/ppData_LumiR4);
    cout << __LINE__ << endl;
    
    if(sumpTCutData || cleantool || sumpTANDcleanjets){
    //w/ SumpT Cut
    cout << __LINE__ << endl;  
    param_Datadist1[WSumpT]->SetLineColor(kBlue);
    param_Datadist1[WSumpT]->SetMarkerColor(kBlue);
    param_Datadist1[WSumpT]->SetMarkerStyle(20);
    param_Datadist1[WSumpT]->Scale(1/ppData_LumiR4);
    cout << __LINE__ << endl;
    }
  }
  param_MCdis2[indexPlotNoMatchReq]->SetLineColor(kBlack);
  param_MCdis2[indexPlotNoMatchReq]->SetMarkerColor(kBlack);
  param_MCdis2[indexPlotNoMatchReq]->SetMarkerStyle(20);
 
  cout << __LINE__ << endl;
  cout << __LINE__ << endl;
  for(int iDataPlots =0; iDataPlots< data_tot; iDataPlots++){
    if(iDataPlots==0){
      param_Datadist1[iDataPlots]->Draw();
      ptRanges_legends->AddEntry( param_Datadist1[iDataPlots],"pp Data","pe");
      param_Datadist1[index_data]->GetXaxis()->SetRangeUser(230,1800);
      param_Datadist1[index_data]->GetYaxis()->SetRangeUser(1e-8,3000);
    }
    cout << "iDataPlots / " << iDataPlots << endl;
    cout << "sumpTCutData: " << sumpTCutData << endl;
    cout << "cleantool:  " << cleantool << endl;
    cout << "sumpTANDcleanjets: " << sumpTANDcleanjets << endl;
    if(iDataPlots!=0 && sumpTCutData){
      param_Datadist1[iDataPlots]->Draw("same");
      cout << "This is the index of your jet: " << iJetIndex << endl;
      cout << "This is the sumpTCut for this jet: " << sumpTCuts[iJetIndex] << endl;
      ptRanges_legends->AddEntry(param_Datadist1[iDataPlots],Form("pp Data w/ #sum p^{trk}_{T} > %d GeV",sumpTCuts[iJetIndex]),"pe");
    }else if(iDataPlots!=0 && cleantool){
      param_Datadist1[iDataPlots]->Draw("same");
      ptRanges_legends->AddEntry(param_Datadist1[iDataPlots],"pp Data w/ Jet Cleaning Tool (LooseBad)","pe");
      cout << __LINE__ << endl;
    }else if(iDataPlots!=0 && sumpTANDcleanjets){
      param_Datadist1[iDataPlots]->Draw("same");
      ptRanges_legends->AddEntry(param_Datadist1[iDataPlots],Form("pp Data w/ Jet Cleaning Tool (LooseBad) + #sum p^{trk}_{T} > %d GeV",sumpTCuts[iJetIndex]),"pe");
    }
  }
  
  cout << "plotting param_MC_dis2 and this as its index: " << indexPlotNoMatchReq << endl;
  param_MCdis2[indexPlotNoMatchReq]->Draw("same");

  if(indexPlotNoMatchReq == 1){
    ptRanges_legends->AddEntry(param_MCdis2[indexPlotNoMatchReq],"pp PYTHIA Reco (No Matching Req.)","pe");
  }else{

    ptRanges_legends->AddEntry(param_MCdis2[indexPlotNoMatchReq],"pp PYTHIA Reco (Matching Req.)","pe");
    
  }

   ptRanges_legends->SetBorderSize(0);
   ptRanges_legends->Draw("same");

   /**
  Text_AtlasInt(0.7005013,0.8119816,0.8984962,0.8543779);
  Text_AntikTJetInfo(0.6967419,0.7400922,0.89599,0.7935484, jetRadius);
  Text_AtlasIt(0.89,0.85);
   */
   Text_Info("pp",-1, 10, 0.1265664,0.7456221,0.3245614,0.8746544, "", 0.031,11);
  
  canv->cd();          // Go back to the main canvas before defining pad2
  
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.2);
  pad2->SetGridx(); // vertical grid
  pad2->Draw();
  pad2->cd();

  TLegend *ratio_legends= new TLegend(0.6629073,0.7341935,0.8408521,0.9406452,NULL,"brNDC");
  ratio_legends->SetTextSize(0.065);
  ratio_legends->SetBorderSize(0);
  TH1D *data_hist[2];
  data_hist[0] = (TH1D*) param_Datadist1[0]->Clone("");
  data_hist[1] = (TH1D*) param_Datadist1[1]->Clone("");  
  
  // data_hist[0] = (TH1D*) param_Datadist1[0]->Clone(nameDataDis[0].c_str());
  // data_hist[1] = (TH1D*) param_Datadist1[1]->Clone(nameDataDis[1].c_str());  
  
  for(int iRatioPlots = 0; iRatioPlots < data_tot; iRatioPlots++){
    
    
     data_hist[iRatioPlots]->Divide(param_MCdis2[indexPlotNoMatchReq]); 
     //R=1.0 Jets Color Scheme & Marker Style
     if((jetRadius == 1.0) && (iRatioPlots==WOSumpT) ){
       //w/o SumpT Cut
       data_hist[iRatioPlots]->SetLineColor(kRed+2); 
       data_hist[iRatioPlots]->SetMarkerColor(kRed+2);
       data_hist[iRatioPlots]->SetMarkerStyle(21);
     }else if((jetRadius == 1.0) && (iRatioPlots==WSumpT)){
       //w/ SumpT Cut
       data_hist[iRatioPlots]->SetLineColor(kRed);
       data_hist[iRatioPlots]->SetMarkerColor(kRed);
       data_hist[iRatioPlots]->SetMarkerStyle(20);
     }
    
     //R=0.4 Color Scheme
     if((jetRadius == (float)0.4) && (iRatioPlots==WOSumpT) ){
       //w/o SumpT Cut
       data_hist[iRatioPlots]->SetLineColor(kAzure-4);
       data_hist[iRatioPlots]->SetMarkerColor(kAzure-4);
        data_hist[iRatioPlots]->SetMarkerStyle(21);
     }else if((jetRadius == (float)0.4) && (iRatioPlots==WSumpT)){
       //w/ SumpT Cut
       data_hist[iRatioPlots]->SetLineColor(kBlue);
       data_hist[iRatioPlots]->SetMarkerColor(kBlue);
       data_hist[iRatioPlots]->SetMarkerStyle(20);
     }       
    
     if(iRatioPlots==0){
       data_hist[iRatioPlots]->Draw("ep");
       ratio_legends->AddEntry(data_hist[iRatioPlots],"pp Data/MC","pe");
       data_hist[iRatioPlots]->SetMinimum(y_low);  // Define Y ..
       data_hist[iRatioPlots]->SetMaximum(y_high); // .. range
       data_hist[iRatioPlots]->Sumw2();
       data_hist[iRatioPlots]->SetStats(0);      // No statistics on lower plot
       data_hist[iRatioPlots]->SetTitle("");

        // Y axis ratio plot settings
       data_hist[iRatioPlots]->GetYaxis()->SetTitle("Data/MC");
       data_hist[iRatioPlots]->GetYaxis()->SetNdivisions(505);
       data_hist[iRatioPlots]->GetYaxis()->SetTitleSize(20);
       data_hist[iRatioPlots]->GetYaxis()->SetTitleFont(43);
       data_hist[iRatioPlots]->GetYaxis()->SetTitleOffset(1.55);
       data_hist[iRatioPlots]->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
       data_hist[iRatioPlots]->GetYaxis()->SetLabelSize(15);
       // X axis ratio plot settings
       data_hist[iRatioPlots]->GetXaxis()->SetTitleSize(20);
       data_hist[iRatioPlots]->GetXaxis()->SetTitleFont(43);
       data_hist[iRatioPlots]->GetXaxis()->SetTitleOffset(3.1);
       data_hist[iRatioPlots]->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
       data_hist[iRatioPlots]->GetXaxis()->SetLabelSize(15);

        
       data_hist[iRatioPlots]->GetXaxis()->SetRangeUser(294,1800);

       
     }
     if(iRatioPlots!=0){
       data_hist[iRatioPlots]->Draw("ep same");
       ratio_legends->AddEntry(data_hist[iRatioPlots],Form("pp Data/MC (w/ #sum p^{trk}_{T} > %d GeV)",sumpTCuts[iJetIndex]),"pe");
       TLine *line_atOne = new TLine(294,1.0,1800.0,1.0);
       line_atOne->Draw("same");
     }

     
  }//Looping over ratio plots
  
  ratio_legends->Draw("same");

   
  return;
}


void Ratio_Plots_ForPbPb(TH1D * param_Datadist, string nameDataDis,TH1D* param_MCdis,int jetRadiusVal,string centBin_strg,string param,float numEvents_CentBin,float y_low = 0 , float y_high = 2,bool mc_MatchOrNoMatch = false, double yLow = 1e-12,double yMax = 1e6){

  int indexPlotNoMatchReq = 0;
  float centBin_TAA = tAAMap[centBin_strg];
  
  if(mc_MatchOrNoMatch == false){
    //This means we want to plot mc pt distributions w/ matching requirement to truth jets
    //Plot all of the Reco jets that passed pT cut and eta cut
    cout << "We will plot all of reco jets! No matching requirement..." << endl; 
    indexPlotNoMatchReq = 1;
  }
  

  int iJetIndex = -1;
      

  int jetRadius = jetRMap[jetRadiusVal];
  
  TCanvas *canv = new TCanvas(Form("%sDis_R%d",param.c_str(),(int)jetRadius), Form("%sDis_R%d",param.c_str(),(int)jetRadius), 800, 800);

  TLegend *ptRanges_legends= new TLegend(0.5263158,0.5124424,0.877193,0.6599078,NULL,"brNDC");
  
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
  pad1->SetBottomMargin(0); // Upper and lower plot are joined}
  if("pT" == param){
    pad1->SetLogy(1); 
    
  }
  pad1->Draw();             // Draw the upper pad: pad
  pad1->cd();               // pad1 becomes the current pad

  int iJetR = R10;
  if(jetRadiusVal != R10)iJetR=R4;
    
  param_Datadist->SetTitle("");
  param_Datadist->GetYaxis()->SetTitle("#frac{1}{<T_{AA}>}#frac{1}{N_{evt}} #frac{d^{2}N_{jet}}{dp_{T}d#eta} [nb]");
  param_Datadist->GetXaxis()->SetTitle("p_{T} [GeV]");
  param_Datadist->SetStats(0);
  param_Datadist->GetYaxis()->SetTitleSize(0.03);
 
  param_Datadist->SetLineColor(colores_pTDis[jetRadius]);
  param_Datadist->SetMarkerColor(colores_pTDis[jetRadius]);
  param_Datadist->SetMarkerStyle(21);
  
  
  //param_Datadist->Scale(1/LumNumPbPbData[iJetR]);
  cout << "Num events: " << numEvents_CentBin << endl;
  cout << "tAA value: " << centBin_TAA << endl;
  cout << "eta Range: " << (1.5*2) << endl;
  
  param_Datadist->Scale(1/numEvents_CentBin);
  param_Datadist->Scale(1/centBin_TAA);
  param_Datadist->Scale(1/(1.5*2));
  param_Datadist->Scale(1.,"width");

  cout << "******************************This might be your stupid overflow bin: " << param_Datadist->GetBinContent(21) << endl;
  
  param_MCdis->SetLineColor(kBlack);
  param_MCdis->SetMarkerColor(kBlack);
  param_MCdis->SetMarkerStyle(20);

  param_MCdis->Scale(1.,"width");
  param_MCdis->Scale(1/(1.5*2));
  
  
  param_Datadist->GetXaxis()->SetRangeUser(100,1800);
  
  param_Datadist->SetMaximum(yMax);
  param_Datadist->SetMinimum(yLow);
  
  
  param_Datadist->Draw();
  ptRanges_legends->AddEntry(param_Datadist,"Pb+Pb Data","pe");
  
  param_MCdis->Draw("same");

  if(indexPlotNoMatchReq == 1){
    ptRanges_legends->AddEntry(param_MCdis,"Pb+Pb MC Reco (No Matching Req.)","pe");
  }else{
    ptRanges_legends->AddEntry(param_MCdis,"Pb+Pb MC  Reco (Matching Req.)","pe");
  }

 ptRanges_legends->SetBorderSize(0);
 ptRanges_legends->Draw("same");

  
 Text_Info("pbpbmc",-1, jetRadiusVal,0.6904762,0.7269565,0.8909774,0.8869565,centBin_strg.c_str(),0.03,33);
  

  
  canv->cd();          // Go back to the main canvas before defining pad2
  
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.2);
  pad2->SetGridx(); // vertical grid
  pad2->Draw();
  pad2->cd();
  gPad->SetTicks(1);
  
  TLegend *ratio_legends= new TLegend(0.6629073,0.7341935,0.8408521,0.9406452,NULL,"brNDC");
  ratio_legends->SetTextSize(0.065);
  ratio_legends->SetBorderSize(0);
  TH1D *data_mc_comp= (TH1D*) param_Datadist->Clone(nameDataDis.c_str());
  
  data_mc_comp->Divide(param_MCdis); 
  //R=1.0 Jets Color Scheme & Marker Style
  data_mc_comp->SetLineColor(colores_pTDis[jetRadius]); 
  data_mc_comp->SetMarkerColor(colores_pTDis[jetRadius]);
  data_mc_comp->SetMarkerStyle(21);

    
     
  data_mc_comp->Draw("ep");
  ratio_legends->AddEntry(data_mc_comp,"Pb+Pb Data/MC","pe");
  data_mc_comp->SetMinimum(y_low);  // Define Y ..
  data_mc_comp->SetMaximum(y_high); // .. range
  data_mc_comp->Sumw2();
  data_mc_comp->SetStats(0);      // No statistics on lower plot
  data_mc_comp->SetTitle("");

   // Y axis ratio plot settings
  data_mc_comp->GetYaxis()->SetTitle("Data/MC");
  data_mc_comp->GetYaxis()->SetNdivisions(505);
  data_mc_comp->GetYaxis()->SetTitleSize(20);
  data_mc_comp->GetYaxis()->SetTitleFont(43);
  data_mc_comp->GetYaxis()->SetTitleOffset(1.55);
  data_mc_comp->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  data_mc_comp->GetYaxis()->SetLabelSize(15);
  // X axis ratio plot settings
  data_mc_comp->GetXaxis()->SetTitleSize(20);
  data_mc_comp->GetXaxis()->SetTitleFont(43);
  data_mc_comp->GetXaxis()->SetTitleOffset(3.1);
  data_mc_comp->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  data_mc_comp->GetXaxis()->SetLabelSize(15);

   
  data_mc_comp->GetXaxis()->SetRangeUser(120,1800);

  


  TLine *line_atOne = new TLine(120,1.0,1600.0,1.0);
  line_atOne->SetLineColor(kBlack);
  line_atOne->Draw("same");

     

 
 ratio_legends->Draw("same");

 
 return;

}





///SumpT Distributions Functions

void sumpTDis(TH1D* sumpTDataDis[], TH1D* sumpTMCDis[], float jetRadius, float pTRanges_Arr[], int tot_pTRanges,int y_low=0,int y_high=2, int x_low_Bin_integt = 50, int x_high_Bin_integt = 900, string tagName = ""){
  int jet_index = -1;
  if(jetRadius == (float)1.0)jet_index=R10;
  if(jetRadius == (float)0.4)jet_index=R4;
 cout << __LINE__ << endl;
 cout << "Total pT ranges: " << tot_pTRanges << endl;
  for(int ipTRange = 0; ipTRange < tot_pTRanges; ipTRange++){
    cout << __LINE__ << endl;
    TCanvas *canv = new TCanvas(Form("MatchedJets_sumpTDis_R%d_pTRange_pTRange_%d_%d_%s",JetRadius[jet_index],(int)pTRanges_Arr[ipTRange],(int)pTRanges_Arr[ipTRange+1],tagName.c_str()), Form("MatchedJets_sumpTDis_R%d_pTRange_%d_%d_%s",JetRadius[jet_index],(int)pTRanges_Arr[ipTRange],(int)pTRanges_Arr[ipTRange+1],tagName.c_str()), 800, 800);
    cout <<	__LINE__ << endl;
    TLegend *pt_lgds= new TLegend(0.6817043,0.5511521,0.89599,0.6506912,NULL,"brNDC");
    cout <<	__LINE__ << endl;
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
    pad1->SetBottomMargin(0); // Upper and lower plot are joined
    pad1->SetLogy(1);
    pad1->Draw();             // Draw the upper pad: pad
    pad1->cd();               // pad1 becomes the current pad

    cout << __LINE__ << endl;
    //Making weight
    sumpTDataDis[ipTRange]->Scale(1/ppDataLumiVals[jet_index]); 
    cout << __LINE__ << endl;
    float sumpTDisData_Intgrl = sumpTDataDis[ipTRange]->Integral(sumpTDataDis[ipTRange]->FindBin(x_low_Bin_integt),sumpTDataDis[ipTRange]->FindBin(x_high_Bin_integt));
    cout << __LINE__ << endl;
    float sumpTDisMC_Intgrl = sumpTMCDis[ipTRange]->Integral(sumpTMCDis[ipTRange]->FindBin(x_low_Bin_integt),sumpTMCDis[ipTRange]->FindBin(x_high_Bin_integt));
    cout << __LINE__ << endl;
    cout << "Integral value for data: " << sumpTDisData_Intgrl << endl;
    cout << "Integral value for MC: " << sumpTDisMC_Intgrl << endl;
    
    float weight =sumpTDisData_Intgrl/sumpTDisMC_Intgrl;
    //Data
    cout << "Weight value: " << weight << endl;
    //sumpTDataDis[ipTRange]->Scale(weight);
    sumpTDataDis[ipTRange]->SetTitle("");
    sumpTDataDis[ipTRange]->GetYaxis()->SetTitle("N/Lumi [nb]");
    sumpTDataDis[ipTRange]->GetXaxis()->SetTitle("#sum p^{trk}_{T}");
    sumpTDataDis[ipTRange]->SetStats(0);
    if(dbg_bFun)cout << __LINE__ << endl;
    sumpTDataDis[ipTRange]->SetLineColor(sumpTDisColors[jet_index]);
    sumpTDataDis[ipTRange]->SetMarkerColor(sumpTDisColors[jet_index]);
    sumpTDataDis[ipTRange]->SetMarkerStyle(21);
    
    if(dbg_bFun)cout << __LINE__ << endl;
    //MC
    
    for(int iBin =sumpTMCDis[ipTRange]->FindBin(x_high_Bin_integt); iBin < sumpTMCDis[ipTRange]->FindBin(x_high_Bin_integt) + 1; iBin++){
      sumpTMCDis[ipTRange]->SetBinContent(iBin,sumpTMCDis[ipTRange]->GetBinContent(iBin)*weight);
    }
    
    sumpTMCDis[ipTRange]->SetLineColor(kBlack);
    sumpTMCDis[ipTRange]->SetMarkerColor(kBlack);
    sumpTMCDis[ipTRange]->SetMarkerStyle(20);
    if(dbg_bFun)cout << __LINE__ << endl;
    sumpTDataDis[ipTRange]->Draw();
    sumpTMCDis[ipTRange]->Draw("same");
    
    pt_lgds->AddEntry(sumpTDataDis[ipTRange],"pp Data","pe");
    pt_lgds->AddEntry(sumpTMCDis[ipTRange],"pp PYTHIA Reco","pe");
    pt_lgds->SetBorderSize(0);
    pt_lgds->Draw("same");
    if(dbg_bFun)cout << __LINE__ << endl;

    Text_Info("pp",-1, JetRadius[jet_index], 0.1265664,0.7456221,0.3245614,0.8746544, "", 0.031,11);
    
    //Text_AtlasInt(0.7005013,0.8119816,0.8984962,0.8543779);
    //Text_AntikTJetInfo(0.6967419,0.7400922,0.89599,0.7935484, jetRadius);

    TPaveText *pt_range_info = new TPaveText(0.6929825,0.6764977,0.8922306,0.7133641,"brNDC");
    pt_range_info->SetTextAlign(12);
    pt_range_info->SetTextSize(0.035);
    pt_range_info->SetFillColor(0);
    pt_range_info->SetBorderSize(0);
    pt_range_info->SetShadowColor(0);
    pt_range_info->SetTextAlign(33);
    pt_range_info->AddText(Form("%d < p_{T} < %d GeV",(int)pTRanges_Arr[ipTRange],(int)pTRanges_Arr[ipTRange+1]));
    pt_range_info->Draw("same");


    //Text_AtlasIt(0.89,0.85);
    
    
    canv->cd();
    if(dbg_bFun)cout << __LINE__ << endl;
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.2);
    pad2->SetGridx(); // vertical grid
    pad2->Draw();
    pad2->cd();

    if(dbg_bFun)cout << __LINE__ << endl;
    TH1D *data_hist = (TH1D*) sumpTDataDis[ipTRange]->Clone();
    data_hist->SetLineColor(sumpTDisColors[jet_index]);
    data_hist->SetMarkerColor(sumpTDisColors[jet_index]);
    data_hist->SetMarkerStyle(20);
    data_hist->Divide(sumpTMCDis[ipTRange]);
    data_hist->Draw("ep");
    if(dbg_bFun)cout << __LINE__ << endl;
    data_hist->SetMinimum(y_low);  // Define Y ..
    data_hist->SetMaximum(y_high); // .. range
    data_hist->Sumw2();
    data_hist->SetStats(0);      // No statistics on lower plot
    data_hist->SetTitle("");

    TLine *line_atOne = new TLine(0.0,1.0,900,1.0);
    line_atOne->Draw("same");
    
    // Y axis ratio plot settings
    data_hist->GetYaxis()->SetTitle("Data/MC");
    data_hist->GetYaxis()->SetNdivisions(505);
    data_hist->GetYaxis()->SetTitleSize(20);
    data_hist->GetYaxis()->SetTitleFont(43);
    data_hist->GetYaxis()->SetTitleOffset(1.55);
    data_hist->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    data_hist->GetYaxis()->SetLabelSize(15);
    if(dbg_bFun)cout << __LINE__ << endl;
    // X axis ratio plot settings
    data_hist->GetXaxis()->SetTitleSize(15);
    data_hist->GetXaxis()->SetTitleFont(43);
    data_hist->GetXaxis()->SetTitleOffset(4.4);
    data_hist->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    data_hist->GetXaxis()->SetLabelSize(15);
   if(dbg_bFun)cout << __LINE__ << endl;
  }//looping over pT ranges



}

void pTFor_DiffJZSamp(string reco_or_truth,TH1D* CombineJZ, TH1D* pTDis_Array[], string pp_or_pbpb){
TCanvas *jets = new TCanvas(Form("%sJets",reco_or_truth.c_str()),Form("%sJets",reco_or_truth.c_str()),800,600);
 gStyle->SetOptStat(0);
 cout << __LINE__ << endl;
 TLegend *jz_lg =  new TLegend(0.6879699,0.4713043,0.8734336,0.6678261,NULL,"brNDC");
 jz_lg->SetBorderSize(0);
 cout << __LINE__ << endl;
 for(int ijZSample =0; ijZSample <  4; ijZSample++){ //numJZSamples

   cout << __LINE__ << endl;
   if(ijZSample==0){
     cout << __LINE__ << endl;
     CombineJZ->SetTitle("");
     cout << __LINE__ << endl;
     CombineJZ->SetMarkerStyle(21);
     cout << __LINE__ << endl;
     CombineJZ->SetMarkerColor(jzColors[ijZSample]);
     CombineJZ->GetXaxis()->SetTitle("p^{Truth}_{T}");
     CombineJZ->GetYaxis()->SetTitle("N/Lumi [nb]");
     CombineJZ->Draw();
     cout << __LINE__ << endl;
     jz_lg->AddEntry(CombineJZ,"All JZ Samples","pe");
     cout << __LINE__ << endl;
     
     pTDis_Array[ijZSample]->SetMarkerStyle(33);
     pTDis_Array[ijZSample]->SetMarkerSize(1.5);
     pTDis_Array[ijZSample]->SetMarkerColor(jzColors[ijZSample+1]);
     pTDis_Array[ijZSample]->Draw("same");
     jz_lg->AddEntry(pTDis_Array[ijZSample],Form("JZ%d Sample",ijZSample+2),"pe");
     cout << __LINE__ << endl;
   }
   cout << __LINE__ << endl;
   if(ijZSample!=0){
     cout << __LINE__ << endl;
     pTDis_Array[ijZSample]->SetMarkerSize(1.5);
     pTDis_Array[ijZSample]->SetMarkerStyle(33);
     pTDis_Array[ijZSample]->SetMarkerColor(jzColors[ijZSample+1]);
     pTDis_Array[ijZSample]->SetLineColor(jzColors[ijZSample+1]);
     pTDis_Array[ijZSample]->Draw("same");
     jz_lg->AddEntry(pTDis_Array[ijZSample],Form("JZ%d Sample",ijZSample+2),"pe");
     cout << __LINE__ << endl;
   }
   cout << __LINE__ << endl;
 }//JZ Sample Loop
 
 Text_Info(pp_or_pbpb.c_str(),-1, 4,0.6904762,0.7269565,0.8909774,0.8869565,"",0.03,33);
 jz_lg->Draw("same");

}



void UnfoldRatioPlots(TH1D * recopTDis_NoUnfold,TH1D* recopTDis_Unfold[], int numIterations, TH1D * trthpTDis,TH2D* respMatrix,int jetRadius,int binsTot, double array[],bool drawRespMatrix = false,string halfOrfull="Half_Closure",float yMin=4e-4,float yMax=4e3,float xMin=0,float xMax = 1300,bool zoomIn = false,string sample = "ppmc",string CentBin = "",bool savepdf = false,string extratag = ""){
  string zoomOrNo = "";
  
  float ymin_Ratio = 0;
  float ymax_Ratio = 2.45;
  
  if(zoomIn){
    ymin_Ratio = 0.94;
    ymax_Ratio = 1.031;
    zoomOrNo = "_ZoomIn";
  }  
  
  TCanvas *canv = new TCanvas(Form("pTDis_R%d_%s_Test_%s%s_%s",jetRadius,halfOrfull.c_str(),CentBin.c_str(),zoomOrNo.c_str(),extratag.c_str()), Form("pTDis_R%d_%s_Test_%s%s_%s",jetRadius,halfOrfull.c_str(),CentBin.c_str(),zoomOrNo.c_str(),extratag.c_str()), 800, 800);
  canv->SetLeftMargin(0.12);

  gStyle->SetOptStat(0);
  
  TLegend *lgnd;
  if(CentBin == ""){
    lgnd  = new TLegend(0.6741855,0.6866359,0.887218,0.883871,NULL,"brNDC");
  }else{
    lgnd  = new TLegend(0.6829574,0.6258065,0.89599,0.8230415,NULL,"brNDC");
  }

  
  lgnd->SetBorderSize(0);
  
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
  pad1->SetBottomMargin(0); // Upper and lower plot are joined}
  pad1->SetLogy();
  pad1->RangeAxis(xMin,yMin,xMax,yMax);
  pad1->Draw();             // Draw the upper pad: pad
  pad1->cd();               // pad1 becomes the current pad

  gPad->SetLogx(1);
  
  
  
  for(int iIteration = 0; iIteration < numIterations; iIteration++){
    cout << "iteration: " << iIteration << endl;
    cout << __LINE__ << endl;
    if(iIteration == 0){
      
      cout << __LINE__ << endl;
      trthpTDis->SetTitle("");
      cout << __LINE__ <<	endl;
      trthpTDis->GetYaxis()->SetTitle("#frac{dN}{dp_{T}Lumi} [nb]");
      trthpTDis->GetXaxis()->SetTitle("p_{T} [GeV]");
      trthpTDis->GetYaxis()->CenterTitle(true);
      trthpTDis->SetStats(0);
      
      trthpTDis->GetXaxis()->SetRangeUser(xMin,xMax);
      
      trthpTDis->SetMarkerStyle(21);
      trthpTDis->SetMarkerSize(1.3);
      trthpTDis->SetMarkerColor(kBlack);
      
      
      recopTDis_NoUnfold->SetMarkerStyle(20);
      recopTDis_NoUnfold->SetMarkerColor(unfoldColors[iIteration]);
      recopTDis_NoUnfold->SetLineWidth(2);
      
      recopTDis_NoUnfold->SetLineColor(unfoldColors[iIteration]);


      
      lgnd->AddEntry(trthpTDis,"Truth","pe");
      lgnd->AddEntry(recopTDis_NoUnfold,"Reco","pe");

      recopTDis_NoUnfold->Scale(1.,"width");
      trthpTDis->Scale(1.,"width");
      
      
      trthpTDis->Draw();
      recopTDis_NoUnfold->Draw("same");

      trthpTDis->SetMaximum(yMax);
      trthpTDis->SetMinimum(yMin);
      //trthpTDis->GetXaxis()->SetRangeUser(xMin,xMax);      
    }else{
      cout << __LINE__ << endl;
      recopTDis_Unfold[iIteration]->SetMarkerColor(unfoldColors[iIteration]);
      recopTDis_Unfold[iIteration]->SetLineColor(unfoldColors[iIteration]);
      recopTDis_Unfold[iIteration]->SetMarkerStyle(jzMarkers[iIteration]);
      cout << __LINE__ << endl;
      if(iIteration>1)recopTDis_Unfold[iIteration]->SetMarkerSize(1.5);
      cout << __LINE__ << endl;
      lgnd->AddEntry(recopTDis_Unfold[iIteration],Form("Unfolded Reco (%d)",iIteration+1),"pe");
      
      recopTDis_Unfold[iIteration]->Scale(1.,"width");
      
      recopTDis_Unfold[iIteration]->Draw("same");
    }
    
  }//Iteration
  
   Text_Info(sample.c_str(),-1, jetRadius, 0.1265664,0.7456221,0.3245614,0.8746544, "", 0.031,11);
  
  if(CentBin.c_str() != ""){
    
    TPaveText *Cent_text= new TPaveText(0.6754386,0.8414747,0.8721805,0.8949309,"brNDC");
    Cent_text->SetTextSize(0.035);
    Cent_text->SetTextFont(42);
    Cent_text->SetFillColor(0);
    Cent_text->SetBorderSize(0);
    Cent_text->SetShadowColor(0);
    Cent_text->SetTextAlign(33);
    Cent_text->AddText(Form("Centrality: %s",centBinmap_2015Meas[centBinsMap_2015Meas[CentBin]].c_str()));
    Cent_text->Draw("same");
  }
  
  lgnd->Draw("same");
  
  canv->cd();          // Go back to the main canvas before defining pad2
  
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.2);
  pad2->SetGridx(); // vertical grid
  pad2->Range(xMin,1e-9,xMax,4e3);
  pad2->Draw();
  pad2->cd();
  
  TLegend *ratio_legends= new TLegend(0.7205514,0.2025806,0.8984962,0.4658065,NULL,"brNDC");
  ratio_legends->SetTextSize(0.065);
  ratio_legends->SetBorderSize(0);
  cout << __LINE__ << endl;
  for(int iIterate =0; iIterate < numIterations; iIterate++){
    
    TH1D *ratio = new TH1D(Form("ratio_%d",iIterate),Form("ratio_%d",iIterate),binsTot,array);
    if(iIterate==0){
      cout << __LINE__ << endl;
      ratio->Divide(recopTDis_NoUnfold,trthpTDis,1,1);
      ratio->SetStats(0);
      ratio->Sumw2();
      ratio->SetTitle("");
      ratio->GetXaxis()->SetTitle("p_{T} [GeV]");
      ratio->GetYaxis()->SetTitle("Reco/Truth");
      ratio->SetMarkerColor(kBlack);
      ratio->SetLineColor(kBlack);
      ratio->SetMarkerStyle(21);
      ratio_legends->AddEntry(ratio,"Reco/Truth","pe");

      
       ratio->SetMinimum(ymin_Ratio);  // Define Y ..
       ratio->SetMaximum(ymax_Ratio); // .. range
       ratio->GetXaxis()->SetMoreLogLabels();
       ratio->Sumw2();
       ratio->SetStats(0);      // No statistics on lower plot
       ratio->SetTitle("");

        // Y axis ratio plot settings
       ratio->GetYaxis()->SetTitle("Reco/Truth");
       ratio->GetYaxis()->SetNdivisions(505);
       ratio->GetYaxis()->SetTitleSize(20);
       ratio->GetYaxis()->SetTitleFont(43);
       ratio->GetYaxis()->SetTitleOffset(1.55);
       ratio->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
       ratio->GetYaxis()->SetLabelSize(15);
       // X axis ratio plot settings
       ratio->GetXaxis()->SetNdivisions(505);
       ratio->GetXaxis()->SetTitleSize(20);
       ratio->GetXaxis()->SetTitleFont(43);
       ratio->GetXaxis()->SetTitleOffset(3.1);
       ratio->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
       ratio->GetXaxis()->SetLabelSize(15);
       cout << __LINE__ << endl;

       
       ratio->Draw("pe");
       ratio->GetXaxis()->SetRangeUser(xMin,xMax);
       
       ratio = nullptr;
      
    }else{
      cout << __LINE__ << endl;
      ratio->Divide(recopTDis_Unfold[iIterate],trthpTDis,1,1);
      ratio->SetMarkerStyle(jzMarkers[iIterate]);
      ratio->SetMarkerColor(unfoldColors[iIterate]);
      ratio->SetLineColor(unfoldColors[iIterate]);
      if(iIterate>1)ratio->SetMarkerSize(1.5);
      ratio->Draw("same");
      ratio_legends->AddEntry(ratio,Form("Reco(NOI %d)/Truth",iIterate+1),"pe");
      ratio = nullptr;
    }
    
  }//Iteration Loop 

  TLine *l = new TLine(xMin,1,xMax,1);
  l->Draw("same");
  ratio_legends->Draw("same");
  cout << __LINE__ << endl;
  gPad->SetLogx(1);

  if(savepdf)canv->SaveAs(Form("Unfolding/HalfClosureTest/%s_pTDis_R%d_%s_Test_%s%s_%s.pdf",sample.c_str(),jetRadius,halfOrfull.c_str(),CentBin.c_str(),zoomOrNo.c_str(),extratag.c_str()));
  
  //Responce Matrix
  if(drawRespMatrix){
    cout << __LINE__ << endl;
    TCanvas *respCanv = new TCanvas(Form("RespMatrix_R%d_%s",jetRadius,CentBin.c_str()),Form("RespMatrix_R%d_%s",jetRadius,CentBin.c_str()),700,700);
    gStyle->SetPaintTextFormat("2.1g");
    cout << __LINE__ << endl;
    gPad->SetLogz();
    respMatrix->SetMaximum(1e3);
    cout << __LINE__ << endl;
    respMatrix->SetMinimum(1e-6);
    respMatrix->GetXaxis()->SetLabelSize(0.025);
    respMatrix->GetYaxis()->SetLabelSize(0.025);
    respMatrix->GetZaxis()->SetLabelSize(0.025);
    respMatrix->GetYaxis()->SetTitleSize(0.025);
    respMatrix->GetXaxis()->SetTitleSize(0.025);
    respMatrix->SetTitle("");
    respMatrix->GetYaxis()->SetTitle("p^{Truth}_{T} [GeV]");
    respMatrix->GetXaxis()->SetTitle("p^{Reco}_{T} [GeV]");
    respMatrix->Draw("colz");
    
    Text_Info(sample.c_str(),-1, jetRadius, 0.1265664,0.7456221,0.3245614,0.8746544, "", 0.031,11,1.5,true);
    
   if(CentBin != ""){
     cout << "THis is it: " << CentBin.c_str() << endl;
    TPaveText *Cent_text= new TPaveText(0.6754386,0.8414747,0.8721805,0.8949309,"brNDC");
    Cent_text->SetTextSize(0.035);
    Cent_text->SetTextFont(42);
    int ci = 1181;
    TColor *color = new TColor(ci, 1, 1, 1, " ", 0);
    Cent_text->SetFillColor(ci);
    ci = TColor::GetColor("#000000");
    Cent_text->SetLineColor(ci);
    Cent_text->SetBorderSize(0);
    Cent_text->SetShadowColor(0);
    Cent_text->SetTextAlign(33);
    Cent_text->AddText(Form("Centrality: %s",centBinmap_2015Meas[centBinsMap_2015Meas[CentBin]].c_str()));
    Cent_text->Draw("same");
  }

   respCanv->SaveAs(Form("Unfolding/ResponseMatrices/RespMatrix_R%d_%s%s.pdf",jetRadius,CentBin.c_str(),sample.c_str()));
    
  }

}


void BackfoldRatioPlots(TH1D * recopTDis_NoBackfold,TH1D* recopTDis_Backfold, TH1D* recopTDis_Unfold,TH1D * trthpTDis, int iterationNumber,int jetRadius,string halfOrfull="Half_Closure",float ymin_Ratio = 0, float ymax_Ratio = 2, string sample = "ppmc",string CentBin = ""){

   TCanvas *canv = new TCanvas(Form("pTDis_R%d_%s_BackfoldTest_%s",jetRadius,halfOrfull.c_str(),CentBin.c_str()), Form("pTDis_R%d_%s_BackfoldTest_%s",jetRadius,halfOrfull.c_str(),CentBin.c_str()), 800, 800);
  canv->SetLeftMargin(0.12);

  gStyle->SetOptStat(0);


  TLegend *lgnd;
  
  if(CentBin != ""){
    lgnd  = new TLegend(0.6741855,0.6110599,0.887218,0.8082949,NULL,"brNDC");
  }else{
    lgnd  = new TLegend(0.6829574,0.6258065,0.89599,0.8230415,NULL,"brNDC");
  }


  lgnd->SetBorderSize(0);

  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
  pad1->SetBottomMargin(0); // Upper and lower plot are joined}
  pad1->SetLogy();
  pad1->Draw();             // Draw the upper pad: pad
  pad1->cd();               // pad1 becomes the current pad

  //All Matched Truth Jets
  trthpTDis->SetTitle("");
  trthpTDis->GetYaxis()->SetTitle("#frac{dN}{dp_{T}Lumi} [nb]");
  trthpTDis->GetXaxis()->SetTitle("p_{T} [GeV]");
  trthpTDis->GetYaxis()->CenterTitle(true);
  trthpTDis->SetStats(0);
  trthpTDis->GetYaxis()->SetRangeUser(1e-9,2e3);
  trthpTDis->GetXaxis()->SetRangeUser(100,1300);
 
  trthpTDis->SetMarkerStyle(21);
  trthpTDis->SetMarkerSize(1.3);
  trthpTDis->SetMarkerColor(kBlack);

  //All Matched Reco Jets
  recopTDis_NoBackfold->SetMarkerStyle(20);
  recopTDis_NoBackfold->SetMarkerColor(kAzure+5);
  recopTDis_NoBackfold->SetLineWidth(2);
  recopTDis_NoBackfold->SetLineColor(kAzure+5);

  
  //Unfolded 
  recopTDis_Unfold->SetMarkerStyle(34);
  recopTDis_Unfold->SetMarkerColor(kViolet+6);
  recopTDis_Unfold->SetLineColor(kViolet+6);
  recopTDis_Unfold->SetMarkerSize(1.5);
  
  //Backfolded
  recopTDis_Backfold->SetMarkerStyle(33);
  recopTDis_Backfold->SetMarkerColor(kAzure);
  recopTDis_Backfold->SetMarkerSize(1.5);
  
  for(int iBin = 1; iBin< binsRAA+1; iBin++){
        
    trthpTDis->SetBinContent(iBin,trthpTDis->GetBinContent(iBin)/trthpTDis->GetBinWidth(iBin));
    recopTDis_NoBackfold->SetBinContent(iBin, recopTDis_NoBackfold->GetBinContent(iBin)/recopTDis_NoBackfold->GetBinWidth(iBin));
    recopTDis_Unfold->SetBinContent(iBin,recopTDis_Unfold->GetBinContent(iBin)/recopTDis_Unfold->GetBinWidth(iBin));
    recopTDis_Backfold->SetBinContent(iBin,recopTDis_Backfold->GetBinContent(iBin)/recopTDis_Backfold->GetBinWidth(iBin));
    
  }//Bin loop
      
   trthpTDis->Draw();
   recopTDis_NoBackfold->Draw("same");
   recopTDis_Unfold->Draw("same");
   recopTDis_Backfold->Draw("same");

     
  lgnd->AddEntry(trthpTDis,"Truth","pe");
  lgnd->AddEntry(recopTDis_NoBackfold,"Reco","pe");
  lgnd->AddEntry(recopTDis_Unfold,"Unfolded Reco","pe");
  lgnd->AddEntry(recopTDis_Backfold,"Backfolded","pe");
  
  lgnd->Draw("same");

   
   Text_Info(sample.c_str(),-1, jetRadius, 0.1265664,0.7456221,0.3245614,0.8746544, "", 0.031,11);

   
   TPaveText *Info_text= new TPaveText(0.6754386,0.8414747,0.8721805,0.8949309,"brNDC");
   Info_text->SetTextSize(0.035);
   Info_text->SetTextFont(42);
   Info_text->SetFillColor(0);
   Info_text->SetBorderSize(0);
   Info_text->SetShadowColor(0);
   Info_text->SetTextAlign(33);

   Info_text->AddText(Form("Bayes method with iteration: %d",iterationNumber));
   
   if(CentBin != "")Info_text->AddText(Form("Centrality: %s",centBinTag[CentBin.c_str()].c_str()));
   Info_text->Draw("same"); 
  
   
  canv->cd();          // Go back to the main canvas before defining pad2
  
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.2);
  pad2->SetGridx(); // vertical grid
  pad2->Draw();
  pad2->cd();
  
  TLegend *ratio_legends= new TLegend(0.7205514,0.2025806,0.8984962,0.4658065,NULL,"brNDC");
  ratio_legends->SetTextSize(0.065);
  ratio_legends->SetBorderSize(0);
  //Unfold Reco VS Truth 
  TH1D *ratio_UnfoldRecoTrth = new TH1D("ratio_UnfoldRecoTrth","ratio_UnfoldRecoTrth",binsRAA,rAA_Bins);
  ratio_UnfoldRecoTrth->Divide(recopTDis_Unfold,trthpTDis,1,1,"B");

  ratio_UnfoldRecoTrth->SetStats(0);
  ratio_UnfoldRecoTrth->Sumw2();

  ratio_UnfoldRecoTrth->SetTitle("");
  ratio_UnfoldRecoTrth->GetXaxis()->SetTitle("p_{T} [GeV]");
  
  ratio_UnfoldRecoTrth->SetMarkerColor(kViolet+6);
  ratio_UnfoldRecoTrth->SetLineColor(kViolet+6);

  ratio_UnfoldRecoTrth->SetMarkerSize(1.5);
  
  ratio_UnfoldRecoTrth->SetMarkerStyle(34);
  ratio_legends->AddEntry(ratio_UnfoldRecoTrth,"Unfolded Reco/Truth","pe");
 
  ratio_UnfoldRecoTrth->SetMinimum(ymin_Ratio);  // Define Y ..
  ratio_UnfoldRecoTrth->SetMaximum(ymax_Ratio); // .. range
  ratio_UnfoldRecoTrth->Sumw2();
  ratio_UnfoldRecoTrth->SetStats(0);      // No statistics on lower plot
  ratio_UnfoldRecoTrth->SetTitle("");

     
        // Y axis ratio_UnfoldRecoTrth plot settings
  ratio_UnfoldRecoTrth->GetYaxis()->SetTitle("Ratio");
  ratio_UnfoldRecoTrth->GetYaxis()->SetNdivisions(505);
  ratio_UnfoldRecoTrth->GetYaxis()->SetTitleSize(20);
  ratio_UnfoldRecoTrth->GetYaxis()->SetTitleFont(43);
  ratio_UnfoldRecoTrth->GetYaxis()->SetTitleOffset(1.55);
  ratio_UnfoldRecoTrth->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  ratio_UnfoldRecoTrth->GetYaxis()->SetLabelSize(15);

       // X axis ratio_UnfoldRecoTrth plot settings
   ratio_UnfoldRecoTrth->GetXaxis()->SetTitleSize(20);
   ratio_UnfoldRecoTrth->GetXaxis()->SetTitleFont(43);
   ratio_UnfoldRecoTrth->GetXaxis()->SetTitleOffset(3.1);
   ratio_UnfoldRecoTrth->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   ratio_UnfoldRecoTrth->GetXaxis()->SetLabelSize(15);


   ratio_UnfoldRecoTrth->GetXaxis()->SetRangeUser(100,1300);
   ratio_UnfoldRecoTrth->Draw("pe");
   
   //Backfold Reco VS Reco (w/o Unfold & Backunfold)
   TH1D *ratio_BackfoldRecoReco = new TH1D("ratio_BackfoldRecoReco","ratio_BackfoldRecoReco",binsRAA,rAA_Bins);
   ratio_BackfoldRecoReco->Divide(recopTDis_Backfold,recopTDis_NoBackfold,1,1,"B");
   ratio_BackfoldRecoReco->SetMarkerStyle(33);
   ratio_BackfoldRecoReco->SetMarkerColor(kAzure);
   ratio_BackfoldRecoReco->SetLineColor(kAzure);
   ratio_BackfoldRecoReco->SetMarkerSize(1.5);
   ratio_BackfoldRecoReco->Draw("same pe");
   
   ratio_legends->AddEntry(ratio_BackfoldRecoReco,"Backfolded Reco/Reco","pe");

   ratio_legends->Draw("same");


   TLine *line1 = new TLine(100,1,1300,1);
   line1->SetLineColor(kBlack);
   line1->SetLineStyle(6);
   line1->SetLineWidth(2);
   line1->Draw("same");
   
}



void UnfoldingStatError(TH1D* unfoldedHist[], TH1D *NotUnfolded,int TotIter = 18, int startBin=3 , int totBins=10,int JETRad=10, int setYMax=6000, string sampleOrdata =  "ppdata", string centBin = "",bool zoom = false ,string location = "",bool debug=false,string fileNameAndLoc = "nominalIters_PbPb.root",double eta_range = 1.5){
  
  //unfoldedHist with have a higher iteration number
  //lowIterHist has a lower iteration number 
  
  //This is for Error Calculation
  double statError[24] = {}; //#sigma^{2}_{stat}
  double iterError[24] = {}; //#sigma^{2}_{iter}
  double sumError[24] = {}; //sqrt(#sigma^{2}_{stat} + #sigma^{2}_{iter})
  
  double iterXval[24] = {};
  
  double errorStatTemp = 0;
  double errorIterTemp =0;

  //Storing all histograms--Not Unfolded pT Dis + Unfolded pT Dis --> into 1 array.
  //To keep life simple
  cout << "This is the number of iterations: " << TotIter << endl;
  TH1D *histos[100];
  cout << __LINE__ << endl;
  for(int ihist =0; ihist < TotIter; ihist++){
    if(ihist==0){
      cout << __LINE__ << endl;
      histos[ihist] = NotUnfolded; 
    }else{
      cout << __LINE__ << endl;
      histos[ihist] = unfoldedHist[ihist];
      cout << __LINE__ << endl;
      if(ihist==1)histos[ihist]->SetMarkerColor(kRed);
      cout << __LINE__ << endl; 
    }
    cout << __LINE__ << endl;
  }

  
  cout << "This is the iteration and statistical error!" << endl;
  int optimalIterNum = 0;
  double smallestQSValue = 1000000;
  for(int iIter = 1; iIter < TotIter; iIter++){
    iterXval[iIter-1] = iIter;//x-values 

    if(iIter==2)cout << "********PAY ATTENTION! ITERATION NUMBER TWO INFORMATION!!!!*********" << endl;
    if(iIter==6)cout << "********PAY ATTENTION! ITERATION NUMBER SIX INFORMATION!!!*********" << endl;
    cout << "Tot bins: " << totBins << endl;
    cout << "start bin: " << startBin << endl;
    for(int iBin = startBin; iBin < totBins+1; iBin++){
      cout << "iter: " << iIter << endl;
      cout << "iBin: " << iBin << endl;
      cout << "THIS IS THE BIN CENTER: " << histos[iIter]->GetBinCenter(iBin) << endl;


      if(debug)cout << "Value of stat error: " << errorStatTemp << endl;
      if(debug)cout<< "Value of ietration error: " << errorIterTemp << endl;
      if(debug)cout<< "This is bin number: " << iBin << endl;
      
      //Stat Error
      if(debug)cout<< "This is for the histo after applying " << iIter << " number of iterations." << endl;
      if(debug)cout<< "This is the bin center: " << histos[iIter]->GetBinCenter(iBin) << endl;
      if(debug)cout<< "This is the bin content: " << histos[iIter]->GetBinContent(iBin) << endl;
      if(debug)cout<< "Getting STAT Error Bin: " <<  histos[iIter]->GetBinErrorUp(iBin) << endl;
      if(debug)cout<< "squaring error: " <<  histos[iIter]->GetBinErrorUp(iBin)* histos[iIter]->GetBinErrorUp(iBin) << endl;
      
      if(debug)cout<< "Before adding the stat error this is the value of the previous stat error: " << errorStatTemp << endl;
      errorStatTemp = (histos[iIter]->GetBinErrorUp(iBin)*histos[iIter]->GetBinErrorUp(iBin)) + errorStatTemp;
      if(debug)cout<< "total: " << errorStatTemp << endl;
      if(debug)cout<< __LINE__<< endl;
      if(debug)cout<< "iIter: "  << iIter << endl;
      if(debug)cout<< "Bin content: " << histos[iIter]->GetBinContent(iBin) << endl;
       
      if(debug)cout<< "other bin content: "<< histos[iIter-1]->GetBinContent(iBin) << endl;
      //Iter. Error
      errorIterTemp  = pow(histos[iIter]->GetBinContent(iBin) -  histos[iIter-1]->GetBinContent(iBin),2) + errorIterTemp; 
      if(debug)cout<< __LINE__<< endl;
      if(debug)cout<< "This is the iteration error... " << endl;
      if(debug)cout<< "histos[iIter]->GetBinContent(iBin): " << histos[iIter]->GetBinContent(iBin) << endl;
      if(debug)cout<< "histos[iIter-1]->GetBinContent(iBin): " << histos[iIter-1]->GetBinContent(iBin) << endl;
      if(debug)cout<< " histos[iIter]->GetBinContent(iBin) -  histos[iIter-1]->GetBinContent(iBin) : " << histos[iIter]->GetBinContent(iBin) -  histos[iIter-1]->GetBinContent(iBin) << endl;

      if(debug)cout<< "(histos[iIter]->GetBinContent(iBin) -  histos[iIter-1]->GetBinContent(iBin))^2 " << pow(histos[iIter]->GetBinContent(iBin) -  histos[iIter-1]->GetBinContent(iBin),2) << endl;

      
      
    }//Bin Loop

   
    
    if(debug)cout<< "This is the statistical error:  " << errorStatTemp << endl;
    if(debug)cout<< "This is the iteration error: " << errorIterTemp << endl;
    if(debug)cout<< "This is for iteration number: " << iIter << endl;
    
    statError[iIter-1] = sqrt(errorStatTemp);
    iterError[iIter-1] = sqrt(errorIterTemp);
    if(debug)cout<< "This is the quad. sum: " << sqrt(errorIterTemp + errorStatTemp) << endl << endl;
    sumError[iIter-1] = sqrt(errorIterTemp + errorStatTemp); 

    if(debug)cout<< "This is the quad. sum error: " << sumError[iIter-1] << endl;
     
    errorIterTemp = 0;
    errorStatTemp = 0;
    if(sumError[iIter-1] < smallestQSValue){
      optimalIterNum = iIter;
      smallestQSValue = sumError[iIter-1];
    
    }
  }//Iteration Loop


  cout << "This is the samllest QS value: " << smallestQSValue << endl;
  cout << "This is the optimal number of iterations: " << optimalIterNum << endl;
  //This file will save the optimal number of iterations
  TFile *file = new TFile(Form("%s",fileNameAndLoc.c_str()),"UPDATE");
  file->cd();
  cout << "writing this hist: " << Form("%s_Cent_%s_R%d",sampleOrdata.c_str(),centBinmap_2015MeasTag[centBin].c_str(),JETRad) << endl;
  cout << "THis is the start bin number: " << startBin << endl;
  cout << "This is the end bin number: " << totBins << endl;
  
  TH1D* h = new TH1D(Form("%s_Cent_%s_R%d_IterNum",sampleOrdata.c_str(),centBinmap_2015MeasTag[centBin].c_str(),JETRad),Form("%s_Cent_%s_R%d_IterNum",sampleOrdata.c_str(),centBinmap_2015MeasTag[centBin].c_str(),JETRad),1,0,1);
  TH1D *h_startBin_Number = new TH1D(Form("%s_Cent_%s_R%d_StartBinNum",sampleOrdata.c_str(),centBinmap_2015MeasTag[centBin].c_str(),JETRad),Form("%s_Cent_%s_R%d_Start",sampleOrdata.c_str(),centBinmap_2015MeasTag[centBin].c_str(),JETRad),1,0,1);
  TH1D *h_endBin_Number = new TH1D(Form("%s_Cent_%s_R%d_EndBinNum",sampleOrdata.c_str(),centBinmap_2015MeasTag[centBin].c_str(),JETRad),Form("%s_Cent_%s_R%d_End",sampleOrdata.c_str(),centBinmap_2015MeasTag[centBin].c_str(),JETRad),1,0,1);

  h_startBin_Number->SetBinContent(1,startBin);
  h_endBin_Number->SetBinContent(1,totBins);
  
  h->SetBinContent(1,optimalIterNum);
  h->Write("",TObject::kOverwrite);
  h_startBin_Number->Write("",TObject::kOverwrite);
  h_endBin_Number->Write("",TObject::kOverwrite);
  

  file->Close();

  
  //Normalizing
  double statError_Norm[24] = {};
  double iterError_Norm[24] = {};
  double sumError_Norm[24] = {};

  double statError_normConst = 0;
  double iterError_normConst = 0;
  double sumError_normConst = 0;


  
  for(int iBin =0; iBin < 24; iBin++){
    statError_normConst = statError[iBin] + statError_normConst;
    iterError_normConst = iterError[iBin] + iterError_normConst;
    sumError_normConst = sumError[iBin] + sumError_normConst;

  }


  for(int iBin=0;iBin <24; iBin++){
    statError_Norm[iBin] = statError[iBin]/statError_normConst;
    iterError_Norm[iBin] = iterError[iBin]/iterError_normConst;
    sumError_Norm[iBin] = sumError[iBin]/sumError_normConst;
    
  }


  
  TGraph *statErrorg = new TGraph(TotIter,iterXval,statError);
  TGraph *iterErrorg = new TGraph(TotIter,iterXval,iterError);
  TGraph *sumErrorsg = new TGraph(TotIter,iterXval,sumError);
  string zoomStr = "";
  if(zoom){
    if(debug)cout<< "We finna zoom " << endl;
    zoomStr="_ZoomInto";
  }
  cout << "cent 1: " << centBin << endl;
  cout << "cent: " << centBinmap_2015MeasTag[centBin] << endl;
  string centTag = centBinmap_2015MeasTag[centBin];
  
  TCanvas *cav = new TCanvas(Form("StatAndIterErrors_R%d_%s%s",JETRad,zoomStr.c_str(),centTag.c_str()),Form("StatAndIterErrors_R%d_%s%s",JETRad,zoomStr.c_str(),centTag.c_str()),800,600);
  TLegend *leg = new TLegend(0.6115288,0.4678261,0.8709273,0.6191304,NULL,"brNDC");

  leg->SetBorderSize(0);

  iterErrorg->SetTitle("");
  iterErrorg->GetXaxis()->SetTitle("Number of Iterations");
  iterErrorg->GetYaxis()->SetTitle("#sqrt{#delta^{2}}");
  iterErrorg->GetYaxis()->SetLabelSize(0.025);
  iterErrorg->GetXaxis()->SetLabelSize(0.025);
  
  iterErrorg->SetMarkerStyle(20);
  iterErrorg->SetMarkerColor(kViolet-8);
  
  statErrorg->SetMarkerStyle(20);
  statErrorg->SetMarkerColor(kAzure-8);

  sumErrorsg->SetMarkerStyle(33);
  sumErrorsg->SetMarkerColor(kPink-8);
  sumErrorsg->SetMarkerSize(1.5);

  iterErrorg->GetXaxis()->CenterTitle(true);
  iterErrorg->GetYaxis()->CenterTitle(true);
  
  iterErrorg->SetMaximum(setYMax);
  iterErrorg->SetMinimum(0);

  
  
  iterErrorg->Draw(" A P");
  statErrorg->Draw("same P");
  sumErrorsg->Draw("same P");
  
  leg->AddEntry(statErrorg,"Statistical Error", "p");
  leg->AddEntry(iterErrorg, "Iteration Error", "p");
  leg->AddEntry(sumErrorsg, "Quadrature Sum", "p");
  
  if(!zoom)leg->Draw("same");

  if(!zoom)Text_Info(sampleOrdata.c_str(),-1, JETRad, 0.660401,0.6226087,0.858396,0.7513043, centBin, 0.031,11,eta_range);
  
  
  TLine *l1 = new TLine(11,0,11,setYMax);
  l1->SetLineStyle(2);
  //l1->Draw("same");
 
  if(debug)cout<< "this is true or false: " << zoom << endl;
 
  cav->SaveAs(Form("%s/StatAndIterErrors_R%d_%s%s_%s.pdf",location.c_str(),JETRad,centBinmap_2015MeasTag[centBin].c_str(),zoomStr.c_str(),sampleOrdata.c_str())); 

  
}

void histo_to_TGraph(TH1D* ConvtToTGraph, int startBin, int endBin, double *xArray, double *xStatErr, double *yArray, double *yStatErr){
  //These are all pointers to arrays.
  //xArray,xStatErr,yArray,yStatErr

  
  //This will keep track of the bins we want to store in array
  int binForArray = 0;

  
  for(int iBin = startBin; iBin < endBin + 1; iBin++){

    //Storing Points
    yArray[binForArray] = ConvtToTGraph->GetBinContent(iBin);
    xArray[binForArray] = ConvtToTGraph->GetBinCenter(iBin);
    
    //Storing Error bars
    yStatErr[binForArray] = ConvtToTGraph->GetBinError(iBin);      
    xStatErr[binForArray] = ConvtToTGraph->GetBinCenter(iBin)-ConvtToTGraph->GetBinLowEdge(iBin);
    
    binForArray++;
  
  }//Looping Over Bins of Interest

  
  
}
