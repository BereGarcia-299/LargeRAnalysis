#include <assert.h>
#include <cmath>
#include <string>
#include <fstream>
#include <stdio.h>
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
#include "TF1.h"
#include "TRandom3.h"
#include <TAttMarker.h>
#include <TCanvas.h>
#include <TPaveText.h>
#include <dirent.h>
#include <sys/types.h>
#include "TLorentzVector.h"
#include <TInterpreter.h>
#include <TSystem.h>



float  pi = TMath::Pi();

float delta_phi(float phi_uno,float phi_dos){

  float final_phi = fmod(abs(phi_uno - phi_dos), 2*pi);
  float phi = final_phi > pi  ? 2*pi - final_phi : final_phi;
  return phi;

}//-------delta phi function




float deltaR_calc(float eta_1, float eta_2,float phi_1, float phi_2 ){
  float deltaPhi = delta_phi(phi_1,phi_2);
  float deltaEta = eta_1 - eta_2;
  return sqrt(pow(deltaPhi,2) + pow(deltaEta,2));
}//------Calculates deltaR

Double_t getRTrkJESSysPt(float cent, Double_t jtPt)
{
  Double_t factor = 1 - 0.000153*(80.0 - cent);
  return jtPt*factor;
}



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


void Make_Pre_JER_Or_JES(TH1D* array_plots[],string name_plots, float pTRanges[], int tot_ptranges,int bins,float x_low=0.5, float x_max = 1.5){
  
   for(int ipTRange = 0; ipTRange < tot_ptranges; ipTRange++){
     array_plots[ipTRange] = new TH1D(Form("%s_%i_%iGeV",name_plots.c_str(),(int)pTRanges[ipTRange],(int)pTRanges[ipTRange+1]),Form("%s_%i_%iGeV",name_plots.c_str(),(int)pTRanges[ipTRange],(int)pTRanges[ipTRange+1]),bins,x_low,x_max);
    }//pT range loop

  
}






void write_TH1D_1DArray_to_file(TH1D* array[],int size ,string name_plot, float pTSlice[]){
  for(int index = 0; index < size; index++){
    array[index]->Write(Form("%s_%i_%iGeV",name_plot.c_str(),(int)pTSlice[index],(int)pTSlice[index+1]),TObject::kOverwrite);
  }
}//------Writing out histograms to root file (For 1D Arrays)



void MakeGaussianDis(string name_canvas,int numRows, int numColumns, TH1D* array[], int pTRanges[], TF1* tfArray[]){
  
  TCanvas *new_canvas= new TCanvas(name_canvas.c_str(),name_canvas.c_str(),1000,1000);
  
  new_canvas->Divide(numColumns,numRows);
  
  for(int iPlot = 0; iPlot < numRows*numColumns; iPlot++){
    new_canvas->cd(iPlot+1);
    cout << __LINE__ << endl;
    array[iPlot]->SetMarkerStyle(20);
    cout << __LINE__ << endl;
    array[iPlot]->SetLineColor(kBlue);
    cout << __LINE__ << endl;
    array[iPlot]->GetXaxis()->SetTitle("p^{Reco}_{T}/p^{Truth}_{T}");
    array[iPlot]->SetTitle(" ");
    cout << __LINE__ << endl;
    array[iPlot]->Fit(Form("fit_%i_%iGeV",pTRanges[iPlot],pTRanges[iPlot+1]),"M Q N","",0.895,1.05);
    tfArray[iPlot]->SetLineStyle(2);
    cout << __LINE__ << endl;
    tfArray[iPlot]->SetLineColor(kBlue);
    
    array[iPlot]->Draw();
    
    tfArray[iPlot]->Draw("same");
    
    
  }
}//-----Plots the Gaussian Distributions for jets that were matched (Truth and Redco jets) 



void MakeGaussianDis(string name_canvas,int numRows, int numColumns, TH1D* array[], TH1D* extra_array[], int pTRanges[], TF1* arrayTF[], TF1* extra_TF[]){
  
  TCanvas *new_canvas= new TCanvas(name_canvas.c_str(),name_canvas.c_str(),1000,1000);
  
  new_canvas->Divide(numColumns,numRows);
  
  for(int iPlot = 0; iPlot < numRows*numColumns; iPlot++){
    new_canvas->cd(iPlot+1);
    
    //-----Setting Marking Stayle
    array[iPlot]->SetMarkerStyle(20);
    extra_array[iPlot]->SetMarkerStyle(20);
    
    //-----Setting Colors
    extra_array[iPlot]->SetLineColor(kBlack);
    extra_array[iPlot]->SetLineColor(kRed);
    
    //-----Setting Axis Names and title
    array[iPlot]->GetXaxis()->SetTitle("p^{Reco}_{T}/p^{Truth}_{T}");
    array[iPlot]->SetTitle(" ");


    array[iPlot]->Fit(Form("fit_%i_%iGeV",pTRanges[iPlot],pTRanges[iPlot+1]),"M Q N","",0.8,1.2);
    extra_array[iPlot]->Fit(Form("fit_%i_%iGeV",pTRanges[iPlot],pTRanges[iPlot+1]),"M Q N","",0.8,1.2);

    //---Setting Fitting Line Style a& Color
    arrayTF[iPlot]->SetLineStyle(2);
    arrayTF[iPlot]->SetLineColor(kRed);

    extra_TF[iPlot]->SetLineStyle(2);
    extra_TF[iPlot]->SetLineColor(kBlack);

    
    array[iPlot]->Draw();
    
    arrayTF[iPlot]->Draw("same");
    extra_TF[iPlot]->Draw("same");
    extra_array[iPlot]->Draw("same");
  }
}//-----Plots the Gaussian Distributions for jets that were matched (Truth and Redco jets) 

/**
void CreateSuperJetArray(int R4_totNum,int R10_totNum,std::vector<float> *pT_values_R4, std::vector<float> *pT_values_R10){
  
  for(int iJetR4 = 0; iJetR4 < R4_totNum; iJetR4++){
    
    allJet_InfoArray[trthJet][R4][iJetR4][pTPar] = pT_values_R4->at(iJetR4);

  }//looping over total number of jets

  for(int iJetR10 = 0; iJetR10 < R10_totNum; iJetR10++){

    allJet_InfoArray[recoJet][R10][iJetR10][pTPar] =pT_values_R10->at(iJetR10);

  }//looping over total number of truth jets R=1.0




  }*/

