#include "/Users/berenicegarcia/Desktop/Berenice/Spring_2022/LargeRJets/LargeRJetNewCode/Locally/NewTreeVariables.h"
#include "/Users/berenicegarcia/Desktop/Berenice/Spring_2022/LargeRJets/LargeRJetNewCode/Locally/bFunctions.h"
#include "uncersysNames.h"
#include <stdio.h>
#include "TH1.h"
#include "TFile.h"
#include <TCanvas.h>
#include <TColor.h>



int systUncertHist(string collisionSys = "RAA",int jetR = R4,bool drawLeg = false, int centBin = Cent010, string jer_or_jes = "JES", double xMin =138, double xMax = 857,bool jetRate2015Bins = false, bool rAA2015pTBins = true, bool substructpT2018Bins =true,double eta_range = 2.1){

  TFile *hist_sys;
  string location_var_rootFile = "";
  string extra_tag = "";
  string path = "";
  if(jetRate2015Bins){
    location_var_rootFile = "/Users/berenicegarcia/Desktop/Berenice/Spring_2022/LargeRJets/LargeRJetNewCode/Locally/Systematics/2015Meas/";
    extra_tag = "_JetRate2015Bins";
  }else if(rAA2015pTBins){
     extra_tag = "_2015RAARateBins";
     location_var_rootFile = "/Users/berenicegarcia/Desktop/Berenice/Spring_2022/LargeRJets/LargeRJetNewCode/Locally/Systematics/2015Meas/";
  }else if(substructpT2018Bins){
     extra_tag = "_2018DiJetBins";
     location_var_rootFile = "/Users/berenicegarcia/Desktop/Berenice/Spring_2022/LargeRJets/LargeRJetNewCode/Locally/Systematics/2018Meas/";
  }

  
  string nominalIter = "nominalIter_PbPb_withpTShapeWeights.root";
  if(jetRate2015Bins || rAA2015pTBins)nominalIter=Form("nominalIter_PbPb%s_withpTShapeWeights_etaRange28.root",extra_tag.c_str());
  if(rAA2015pTBins)path="rootFiles/2015CentBins/rAA2015Bins/";
  TFile * fileNominalIterations = new TFile(Form("/Users/berenicegarcia/Desktop/Berenice/Spring_2022/LargeRJets/LargeRJetNewCode/Locally/DiagPlots/%s%s",path.c_str(),nominalIter.c_str()),"READ");
  cout << "THIS IS THE FILE WE ARE USING TO GRAB Nominal Iteration Num/start bin/end bin: " << Form("/Users/berenicegarcia/Desktop/Berenice/Spring_2022/LargeRJets/LargeRJetNewCode/Locally/DiagPlots/%s%s",path.c_str(),nominalIter.c_str()) << endl;
  
  double yMax[2][8] = {{1,1,1,1,1,1,1,1},{0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2}}; //This is to plot JER Values
  double yMin[2][8]  = {{1,1,1,1,1,1,1,1},{-0.15,-0.15,-0.12,-0.15,-0.15,-0.15,-0.15,-0.15}};

  double yMaxJES[2][8] = {{1,1,1,1,1,1,1,1},{0.31,0.151,0.067,0.067,0.067,0.123,0.05,0.1}}; //This is tto plot JES Values
  double yMinJES[2][8] = {{1,1,1,1,1,1,1,1},{-0.05,-0.044,-0.013,-0.01,-0.01,-0.015,-0.015,-0.015}}; //This is tto plot JES Values
  
  
  if(collisionSys=="pbpbdata")hist_sys=new TFile(Form("%ssystematics_pbpbdata%s_R4_%s.root",location_var_rootFile.c_str(),extra_tag.c_str(),jer_or_jes.c_str()),"READ");
  if(collisionSys=="ppdata")hist_sys=new TFile(Form("%ssystematics_pp%s_R4_%s.root",location_var_rootFile.c_str(),extra_tag.c_str(),jer_or_jes.c_str()),"READ");
  if(collisionSys=="RAA"){
    hist_sys=new TFile(Form("%ssystematics_RAA%s_R4_%s.root",location_var_rootFile.c_str(),extra_tag.c_str(),jer_or_jes.c_str()),"READ");
    cout << "Grabbed this systematic file: " << Form("%ssystematics_RAA%s_R4_%s.root",location_var_rootFile.c_str(),extra_tag.c_str(),jer_or_jes.c_str()) << endl;
  }
  int totCentBins = 8;
  cout << __LINE__ <<	endl;
  Color_t colores[] = {kViolet-8,kViolet,kBlue-7,kPink+2,kMagenta-4,kMagenta-5,kRed-4,kRed-1,kOrange+1,kOrange+7,kYellow-4,kYellow-2,kYellow-5,kSpring+8,kGreen-4,kCyan-4,kCyan-6,kGray+1,kMagenta-2, kCyan-2, kYellow+3};

  int setMarkers[] = {21,20,22,24,25,26,27,28,29,40,33,29,43,47,42,35,49,48,41,40,24};

  
  int index = 0; 
  if(jer_or_jes == "JER")index=1;
  cout << __LINE__ <<	endl;
  TCanvas *canv = new TCanvas(Form("Zoom_%s_SystematicUncertainties_R%d_%s2015JetRatepTBins",jer_or_jes.c_str(),jetRadius[jetR],collisionSys.c_str()),Form("Zoom_%s_SystematicUncertainties_R%d_%s%s",jer_or_jes.c_str(),jetRadius[jetR],collisionSys.c_str(),extra_tag.c_str()),172,144,1971,709);

    
    cout << __LINE__ <<	endl;
     gStyle->SetOptStat(0);
     canv->Divide(5,2,0.00001,0.00001);
     cout << __LINE__ <<endl;

    TLegend *leg;
    if(jer_or_jes=="JER"){
      leg = new TLegend(0.08668902,0.2305622,0.5313995,0.827659,NULL,"brNDC");
     
    }else if(jer_or_jes=="JES"){
      leg = new TLegend(0.006895134,0.007911709,0.4548566,0.9921283,NULL,"brNDC");
    }else if(jer_or_jes=="Unfolding"){
      leg = new TLegend(0.006895134,0.007911709,0.4548566,0.9921283,NULL,"brNDC");
    }
    cout << __LINE__ <<endl;
     leg->SetBorderSize(0);
     leg->SetTextSize(0.05);
    cout << __LINE__ <<	endl;


    
  TH1D* sys;
  int offsetCanv = 1;
  cout << __LINE__ <<	endl;
  for(int iCentBin =0; iCentBin < 8; iCentBin++){
    cout << __LINE__ <<   endl;
    //if(iCentBin!=centBin)continue;
    cout << "Grabbing this hist: " << Form("QS_R%d_%s_Cent_%s",jetRadius[jetR],jer_or_jes.c_str(),centBins_2015Meas[iCentBin].c_str()) << endl;
    TH1D *quadSum = (TH1D*) hist_sys->Get(Form("QS_R%d_%s_Cent_%s",jetRadius[jetR],jer_or_jes.c_str(),centBins_2015Meas[iCentBin].c_str()));
    cout << __LINE__ <<   endl;
    //quadSum->Scale(100);
    quadSum->SetLineColor(kBlack);
    //quadSum->SetMarkerStyle(33);
    cout << __LINE__ <<	endl;
    
    int totSystUncer = 21;
    if(jer_or_jes == "JER"){
      totSystUncer = 9;
    }else if(jer_or_jes=="Unfolding"){
      totSystUncer =2;
    }

    if(jer_or_jes == "JES" && collisionSys == "ppdata")totSystUncer=20;

    
    if(iCentBin==4)offsetCanv++;
    cout << __LINE__ << endl;
    canv->cd(iCentBin+offsetCanv);
    gPad->SetTicks(1); 
     for(int iSysUncrt =0; iSysUncrt < totSystUncer; iSysUncrt++){
      if(collisionSys == "ppdata" && jer_or_jes == "JES" && iSysUncrt==17)continue;
      cout << "iSysUncrt: "<< iSysUncrt << endl;
       cout << "We will grab this histogram: " << Form("R%d_%s_%d_Cent_%s",jetRadius[jetR],jer_or_jes.c_str(),iSysUncrt,centBins_2015Meas[iCentBin].c_str()) << endl;
       sys = (TH1D*) hist_sys->Get(Form("R%d_%s_%d_Cent_%s",jetRadius[jetR],jer_or_jes.c_str(),iSysUncrt,centBins_2015Meas[iCentBin].c_str()));
       cout << __LINE__ << endl;
      sys->SetMarkerStyle(setMarkers[iSysUncrt]);
      sys->SetMarkerSize(1.3);
      sys->SetMarkerColor(colores[iSysUncrt]);
      sys->SetLineColor(colores[iSysUncrt]);
      cout << __LINE__ << endl;
      //sys->SetLineWidth(2);
      //sys->Scale(100);
      if(iSysUncrt==0){
	cout << __LINE__ << endl;
	sys->SetTitle("");
	cout << __LINE__ << endl;
	sys->GetXaxis()->SetTitle("p_{T} [GeV]");
	sys->GetYaxis()->SetTitle("Relative Uncertainty [%]");
	sys->GetYaxis()->CenterTitle(true);
	sys->GetYaxis()->SetTitleOffset(1.34);
	sys->GetYaxis()->SetTitleSize(0.039);
	sys->GetXaxis()->SetTitleSize(0.039);
	cout << __LINE__ << endl;
	cout << __LINE__ << endl;
	TH1D *startbin = (TH1D*)fileNominalIterations->Get(Form("pbpbdata_Cent_%s_R%d_StartBinNum",centBins_2015Meas[iCentBin].c_str(),jetRadius[jetR]));
	cout << __LINE__ << endl;
        TH1D *endbin = (TH1D*)fileNominalIterations->Get(Form("pbpbdata_Cent_%s_R%d_EndBinNum",centBins_2015Meas[iCentBin].c_str(),jetRadius[jetR]));
	cout << __LINE__ << endl;
        double startBin = startbin->GetBinContent(1); double endBin = endbin->GetBinContent(1);
	cout << "THis is the start bin number: " << startBin << endl;
	cout << "This is the bin center of that bin: " << sys->GetBinCenter(startBin) << endl;
	cout << "THis is the starting point of cent 0-10%: " << (sys->GetBinCenter(startBin)- sys->GetBinWidth(startBin)/2) << endl;
	cout << "This is the end bin number: " << endBin << endl;
	cout << "This is the bincenter of that bin:  " << sys->GetBinCenter(endBin) << endl;
	cout << "This is the end point: " << (sys->GetBinCenter(endBin) +sys->GetBinWidth(endBin)/2) << endl;
	cout << __LINE__ << endl;
	
	
       sys->GetXaxis()->SetRangeUser(sys->GetBinCenter(startBin)+ sys->GetBinWidth(startBin)/2,sys->GetBinCenter(endBin) +sys->GetBinWidth(endBin)/2);
	
	
	sys->Draw("ph same");

	TLine *line = new TLine(sys->GetBinCenter(startBin) - sys->GetBinWidth(startBin)/2,0,sys->GetBinCenter(endBin) + sys->GetBinWidth(endBin)/2,0);
	line->SetLineStyle(2);
	line->SetLineColor(kBlack);
	line->Draw("same");
	if(jer_or_jes=="JER"){
	  sys->SetMaximum(20);
	  sys->SetMinimum(-0.1);
	}
	if(jer_or_jes=="JES"){
	  sys->SetMaximum(13);
	  sys->SetMinimum(-1);
	}
	if(jer_or_jes=="Unfolding"){
	  sys->SetMaximum(7);
	  sys->SetMinimum(0);
	  
	}
	
	quadSum->SetLineWidth(2);
	quadSum->Draw("ph same");
      }else{
	cout << __LINE__ << endl;
	sys->Draw("ph same");
	cout << __LINE__ << endl;
      }
      cout << __LINE__ << endl;
    string nameUncer = "Unfolding";
    cout << __LINE__ << endl;
    if(jer_or_jes == "JER")nameUncer = jerNames[iSysUncrt];
    if(jer_or_jes == "JES")nameUncer = jesNames[iSysUncrt];
    if(jer_or_jes=="Unfolding")nameUncer =UnfoldingSysUncert[iSysUncrt];
    cout << __LINE__ << endl;
    if(iCentBin==0)leg->AddEntry(sys,Form("%s",nameUncer.c_str()),"lp"); 
    cout << __LINE__ << endl;  
    }//looping over systematics  
    cout << __LINE__ << endl;
    
     if(iCentBin==0)leg->AddEntry(quadSum,Form("Total %s Uncertainty",jer_or_jes.c_str()),"pl");
    cout << __LINE__ << endl;
     TPaveText *centText;
       if(collisionSys=="pbpbdata" || collisionSys=="RAA"){
	 cout << __LINE__ << endl;
	 centText = new TPaveText(0.3186728,0.7834916,0.5192901,0.8691983,"brNDC");
       }else if(collisionSys=="ppdata"){
	 cout << __LINE__ << endl;
	 centText = new TPaveText(0.2935294,0.8070298,0.6109811,0.8947526,"brNDC");
       }
     cout << __LINE__ << endl;  
     centText->SetTextSize(0.05);
     centText->SetFillColor(0);
     centText->SetBorderSize(0);
     centText->SetShadowColor(0);
     centText->SetTextAlign(33);
     centText->SetTextFont(42);
     cout << __LINE__ << endl;
     if(collisionSys=="pbpbdata" || collisionSys=="RAA"){
       centText->AddText(Form("Centrality: %s",centBinmap_2015Meas[iCentBin].c_str()));
       centText->Draw("same");
     }else if(collisionSys=="ppdata"){
       centText->AddText(Form("Centrality: %s p_{T} Binning",centBinmap_2015Meas[iCentBin].c_str()));
       centText->Draw("same");
     }
  }//cent Loop

    canv->cd(5);
    Text_Info(collisionSys.c_str(),-1, jetRadius[jetR],0.1307594,0.5177954,0.3270731,0.8691293,"", 0.081,11,eta_range);

    canv->cd(10);
    
    leg->Draw();

    sys = NULL;
    return 0;

}
