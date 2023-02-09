#include "/Users/berenicegarcia/Desktop/Berenice/Spring_2022/LargeRJets/LargeRJetNewCode/Locally/NewTreeVariables.h"
#include "/Users/berenicegarcia/Desktop/Berenice/Spring_2022/LargeRJets/LargeRJetNewCode/Locally/bFunctions.h"
#include "uncersysNames.h"
#include <stdio.h>
#include "TH1.h"
#include "TFile.h"
#include <TCanvas.h>
#include <TColor.h>

int totalSysUncert(string collisionSys = "pbpbdata", int jetYields =-1, int crosssec =0, int jetR = R4, bool debug = false, string nameFile = "2018pTBins_FirstTotSysUncert_pbpbdata_R4", bool plot_totalSysUncert = false, bool jetRate2015Bins = true, bool rAA2015pTBins=false, bool substrct2018pTBins =true ,double etaRange_val = 2.1){
  const int totCentBins = 8;
  string tagfile = "pbpbdata";
  if(collisionSys=="pp"){
    tagfile = "pp";
  }else if(collisionSys=="RAA"){
    tagfile = "RAA";
  }
  string locatn = "";
  string meas2015 = "";
  if(jetRate2015Bins || rAA2015pTBins){
    locatn = "/Users/berenicegarcia/Desktop/Berenice/Spring_2022/LargeRJets/LargeRJetNewCode/Locally/Systematics/2015Meas/"; 
  }else if(substrct2018pTBins){
    locatn = "/Users/berenicegarcia/Desktop/Berenice/Spring_2022/LargeRJets/LargeRJetNewCode/Locally/Systematics/2018MeasSys/";
  }

  if(jetRate2015Bins){
    meas2015 = "_2015JetRateBins";
  }else if(rAA2015pTBins){
    meas2015 = "_2015RAARateBins";
  }else if(substrct2018pTBins){
    meas2015 = "_2018DiJetBins";
  }
  
  TFile *JER_SysUncertFile = new TFile(Form("%ssystematics_%s%s_R%d_JER.root",locatn.c_str(),tagfile.c_str(),meas2015.c_str(),jetRadius[jetR]),"READ");
  TFile *JES_SysUncertFile = new TFile(Form("%ssystematics_%s%s_R%d_JES.root",locatn.c_str(),tagfile.c_str(),meas2015.c_str(),jetRadius[jetR]),"READ");
  TFile *Unf_SysUncertFile = new TFile(Form("%ssystematics_%s%s_R%d_Unfolding.root",locatn.c_str(),tagfile.c_str(),meas2015.c_str(),jetRadius[jetR]),"READ");
  

  //Grabbing pTRanges
   TFile * fileNominalIterations = new TFile("/Users/berenicegarcia/Desktop/Berenice/Spring_2022/LargeRJets/LargeRJetNewCode/Locally/DiagPlots/nominalIter_PbPb_withpTShapeWeights.root","READ");
   if(rAA2015pTBins)fileNominalIterations = new TFile(Form("/Users/berenicegarcia/Desktop/Berenice/Spring_2022/LargeRJets/LargeRJetNewCode/Locally/DiagPlots/nominalIter_PbPb%s_withpTShapeWeights_etaRange28.root",meas2015.c_str()),"READ");

   cout << "WE ARE USING THIS FILE TO GRAB nominal iter. num/start bin/ end bin: " << Form("/Users/berenicegarcia/Desktop/Berenice/Spring_2022/LargeRJets/LargeRJetNewCode/Locally/DiagPlots/nominalIter_PbPb%s_withpTShapeWeights_etaRange28.root",meas2015.c_str()) << endl;

   
  int numbins = 15;
  if(collisionSys == "RAA")numbins=20;
  
  int pTBinsTotR[2][8] = {{11,12,12,11,11,12,11,11},{20,20,20,17,17,18,17,17}};

  //if(!jetRate2015Bins)numbins = pTBinsTotR[jetR][iCentBin];
  
  double pTRange[2][8][2] = {};


  double tAAUncert[] = {0.124,0.132,0.132,0.120,0.0946,0.0651,0.0394,0.0199};//for each centrality bin

  
  for(int iCentBin =0; iCentBin <totCentBins; iCentBin++ ){
    cout << "For this cent bin: " << iCentBin << endl;
    
    pTRange[jetR][iCentBin][0] = ((TH1D*)fileNominalIterations->Get(Form("pbpbdata_Cent_%s_R%d_StartBinNum",centBins_2015Meas[iCentBin].c_str(),jetRadius[jetR])))->GetBinContent(1);

    cout << "This is your start bin: " << pTRange[jetR][iCentBin][0] << endl;
    
    pTRange[jetR][iCentBin][1] = ((TH1D*)fileNominalIterations->Get(Form("pbpbdata_Cent_%s_R%d_EndBinNum",centBins_2015Meas[iCentBin].c_str(),jetRadius[jetR])))->GetBinContent(1);
    
    cout << "THis is the endbin: " <<  pTRange[jetR][iCentBin][1] << endl;
  }
  
  
   //Binning for this analysis
  std::vector<std::vector<double>> bins;
  bins.reserve(8);

  if(jetR==R10){
    //R = 1.0
    bins.emplace_back(std::vector<double>{262, 294, 336, 384, 439, 502, 573, 656, 750, 857, 980, 1120, 1300}); //10-20%
    bins.emplace_back(std::vector<double>{245, 262, 294, 336, 384, 439, 502, 573, 656, 750, 857, 980, 1120}); //20-30%
    bins.emplace_back(std::vector<double>{230, 262, 294, 336, 384, 439, 502, 573, 656, 750, 857, 1000}); //30-40%
    bins.emplace_back(std::vector<double>{230, 262, 294, 336, 384, 439, 502, 573, 656, 750, 857, 1000}); //40-50%
    bins.emplace_back(std::vector<double>{220, 241, 262, 294, 336, 384, 439, 502, 573, 656, 750, 857, 1000}); //50-60%
    bins.emplace_back(std::vector<double>{220, 241, 262, 294, 336, 384, 439, 502, 573, 656, 750, 1000}); //60-70%
    bins.emplace_back(std::vector<double>{220, 241, 262, 294, 336, 384, 439, 502, 573, 656, 750, 1000}); //70-80%~
  }else if(jetR==R4){
    //R=0.4
    bins.emplace_back(std::vector<double>{100, 117, 138, 158, 178, 199, 220, 241, 262, 294, 336, 384, 439, 502, 573, 656, 750, 857, 980, 1120,1300}); //0-10%
    bins.emplace_back(std::vector<double>{100, 117, 138, 158, 178, 199, 220, 241, 262, 294, 336, 384, 439, 502, 573, 656, 750, 857, 980, 1120, 1300}); //10-20%
    bins.emplace_back(std::vector<double>{109, 126, 143, 160, 177, 194, 211, 228, 245, 262, 294, 336, 384, 439, 502, 573, 656, 750, 857, 980, 1120}); //20-30%
    bins.emplace_back(std::vector<double>{100, 117, 138, 158, 177, 198, 230, 262, 294, 336, 384, 439, 502, 573, 656, 750, 857, 1000}); //30-40%
    bins.emplace_back(std::vector<double>{100, 117, 138, 158, 177, 198, 230, 262, 294, 336, 384, 439, 502, 573, 656, 750, 857, 1000}); //40-50%
    bins.emplace_back(std::vector<double>{100, 117, 138, 158, 177, 198, 220, 241, 262, 294, 336, 384, 439, 502, 573, 656, 750, 857, 1000}); //50-60%
    bins.emplace_back(std::vector<double>{100, 117, 138, 158, 177, 197, 220, 241, 262, 294, 336, 384, 439, 502, 573, 656, 750, 1000}); //60-70%
    bins.emplace_back(std::vector<double>{100, 117, 138, 158, 177, 197, 220, 241, 262, 294, 336, 384, 439, 502, 573, 656, 750, 1000}); //70-80%~
    
  }




  TFile *TotSysUncer = new TFile(Form("/Users/berenicegarcia/Desktop/Berenice/Spring_2022/LargeRJets/LargeRJetNewCode/Locally/Systematics/2015Meas/totSystematics/%s.root",nameFile.c_str()),"RECREATE");

  TH1D *JES_UnsrtTot[num2015MeasBins];
  TH1D *JER_UnsrtTot[num2015MeasBins];
  TH1D *Unf_UnsrtTot[num2015MeasBins];
  TH1D *TotSysUncert[num2015MeasBins];


  

  
  
  TCanvas *canv;
  TLegend *leg;
  
  if(plot_totalSysUncert){
    canv= new TCanvas(Form("%s_TotSystematicUncertainties_R%d",tagfile.c_str(),jetRadius[jetR]),Form("%s_TotSystematicUncertainties_R%d",tagfile.c_str(),jetRadius[jetR]),172,144,1971,709);
   canv->Divide(5,2,0.00001,0.00001);
   leg = new TLegend(0.08070972,0.3304026,0.5251422,0.6286601,NULL,"brNDC");
   gStyle->SetOptStat(0);
  }
  int offsetCanv = 1;

  Color_t totUncertColors[] = {kCyan-2,kMagenta-2,kGray+1,kCyan+2,kYellow+2};

  
  
  
  for(int iCentBin =0; iCentBin < num2015MeasBins; iCentBin++){
       if(iCentBin==4 && plot_totalSysUncert)offsetCanv++;
       cout << __LINE__ << endl;
       if(plot_totalSysUncert)canv->cd(iCentBin+offsetCanv);

    
     
      double *array = NULL; 
      if(!jetRate2015Bins && !rAA2015pTBins)array = bins.at(iCentBin).data();
      if(jetRate2015Bins){
	array=jetRateMC;
      }else if(rAA2015pTBins){
	array=rAA_2015Bins[iCentBin];
	numbins=20;
      }
    JER_UnsrtTot[iCentBin] = (TH1D*) JER_SysUncertFile->Get(Form("QS_R%d_JER_Cent_%s",jetRadius[jetR],centBins_2015Meas[iCentBin].c_str()));

    if(debug)cout << "This is the histo that we grabbed: " << Form("QS_R%d_JER_Cent_%s",jetRadius[jetR],centBins_2015Meas[iCentBin].c_str()) << endl;
    
    JES_UnsrtTot[iCentBin] = (TH1D*) JES_SysUncertFile->Get(Form("QS_R%d_JES_Cent_%s",jetRadius[jetR],centBins_2015Meas[iCentBin].c_str()));
    Unf_UnsrtTot[iCentBin] = (TH1D*) Unf_SysUncertFile->Get(Form("QS_R%d_Unfolding_Cent_%s",jetRadius[jetR],centBins_2015Meas[iCentBin].c_str()));
    
    

    TotSysUncert[iCentBin] = new TH1D(Form("TotSysUncert_R%d_Cent_%s",jetRadius[jetR],centBins_2015Meas[iCentBin].c_str()),Form("TotSysUncert_R%d_Cent_%s",jetRadius[jetR],centBins_2015Meas[iCentBin].c_str()),numbins,array);

    cout << "PRINTINGOUT ARRAY FOR TOT SYSTEMATICS: " << endl;
    cout << "array[0]: " << array[0] << endl;
    cout << "array[1]: " << array[1] << endl;


    
    if(debug)cout << __LINE__ << endl;

    //if(debug)cout << "This is how many bins we have for this plot: " << pTBinsTotR[jetR][iCentBin] << endl;
    cout << "THIS IS THE NUMBER OF BINS: " << numbins << endl;
    for(int iBin =0; iBin < numbins+4; iBin++){
      if(debug)cout << "This is the bin center for JER " << JER_UnsrtTot[iCentBin]->GetBinCenter(iBin) << " for bin number " << iBin << endl;
      if(debug)cout << "This is the bin center for JES "	<< JES_UnsrtTot[iCentBin]->GetBinCenter(iBin) << " for bin number " << iBin << endl;
      if(debug)cout << "This is the bin center of Total " << TotSysUncert[iCentBin]->GetBinCenter(iBin) << " for bin number " << iBin << endl; 
      if(iCentBin==0){
	cout << "This is the ibin " << iBin << endl;
	cout << "This is the bin center of JER:" << JER_UnsrtTot[iCentBin]->GetBinCenter(iBin) << endl;
	cout << "This is the bin content of JER: " << JER_UnsrtTot[iCentBin]->GetBinContent(iBin) << endl;
	cout << "This is the bin center of JES: " << JES_UnsrtTot[iCentBin]->GetBinCenter(iBin) << endl;
	cout << "This is the bin content of JES: " << JES_UnsrtTot[iCentBin]->GetBinContent(iBin) << endl;
      }

      
      double taaUncert = tAAUncert[iCentBin];
      if(collisionSys=="pp" || collisionSys=="RAA")taaUncert =0;//uncertainty of luminosity of pp
      cout << "THIS IS FOR CENT BIN: " << iCentBin << endl;
      cout << "This is JES UNCERT SQARED: " << pow(JES_UnsrtTot[iCentBin]->GetBinContent(iBin),2) << endl;
      cout << "This is the JER UNCERT SQUARED: " << pow(JER_UnsrtTot[iCentBin]->GetBinContent(iBin),2) << endl;
      cout << "This is the UNFOLDING SQUARED: " << pow(Unf_UnsrtTot[iCentBin]->GetBinContent(iBin),2) << endl;
      double totSysUncert = sqrt(pow(JES_UnsrtTot[iCentBin]->GetBinContent(iBin),2) + pow(JER_UnsrtTot[iCentBin]->GetBinContent(iBin),2) + pow(Unf_UnsrtTot[iCentBin]->GetBinContent(iBin),2) + pow(taaUncert,2));
      TotSysUncert[iCentBin]->SetBinContent(iBin,totSysUncert);
      cout << "THIS IS THE TOTAL SYS: " << totSysUncert << " for bin with center: " << Unf_UnsrtTot[iCentBin]->GetBinCenter(iBin) << endl;
      
      if(iCentBin==0)cout << "This is the total uncertainty: " << totSysUncert << endl;
      
    }


    if(plot_totalSysUncert){
       JER_UnsrtTot[iCentBin]->SetLineColor(totUncertColors[1]);
       JES_UnsrtTot[iCentBin]->SetLineColor(totUncertColors[2]);
       //Unf_UnsrtTot[iCentBin]->SetLineColor(totUncertColors[4]);
       JES_UnsrtTot[iCentBin]->SetLineWidth(2);
       JER_UnsrtTot[iCentBin]->SetLineWidth(2);
       //Unf_UnsrtTot[iCentBin]->SetLineWidth(2);
       TotSysUncert[iCentBin]->SetLineWidth(2);
       TotSysUncert[iCentBin]->SetLineColor(totUncertColors[0]);
       TotSysUncert[iCentBin]->SetTitle("");

       TotSysUncert[iCentBin]->GetXaxis()->SetTitle("p_{T} [GeV]");
       if(collisionSys=="pbpbdata")TotSysUncert[iCentBin]->GetYaxis()->SetTitle("Pb+Pb Yields Uncertainty [%]");
       if(collisionSys=="pp")TotSysUncert[iCentBin]->GetYaxis()->SetTitle("pp Cross-Section Uncertainty [%]");
       TotSysUncert[iCentBin]->SetMinimum(0);
       TotSysUncert[iCentBin]->SetMaximum(15);
       cout << "This is the bin number for the start bin: " << pTRange[jetR][iCentBin][0] << endl;
       cout << "This is the bin center of the startbin: " << Unf_UnsrtTot[iCentBin]->GetBinCenter(pTRange[jetR][iCentBin][0]) << endl;
       cout << "This is adding half of its bin width: " << Unf_UnsrtTot[iCentBin]->GetBinCenter(pTRange[jetR][iCentBin][0]) + (Unf_UnsrtTot[iCentBin]->GetBinWidth(pTRange[jetR][iCentBin][0])/2) << endl;
       cout << "This is the bin number for the end bin: " << pTRange[jetR][iCentBin][1] << endl;
       
       TotSysUncert[iCentBin]->GetXaxis()->SetRangeUser(Unf_UnsrtTot[iCentBin]->GetBinCenter(pTRange[jetR][iCentBin][0])- Unf_UnsrtTot[iCentBin]->GetBinWidth(pTRange[jetR][iCentBin][0])/2,Unf_UnsrtTot[iCentBin]->GetBinCenter(pTRange[jetR][iCentBin][1])+Unf_UnsrtTot[iCentBin]->GetBinWidth(pTRange[jetR][iCentBin][1])/2);
       
       if(iCentBin==0){
      
       
	 TotSysUncert[iCentBin]->Draw();
	 JER_UnsrtTot[iCentBin]->Draw("same");
	 JES_UnsrtTot[iCentBin]->Draw("same");
	 Unf_UnsrtTot[iCentBin]->Draw("same");
	 
         leg->AddEntry(JER_UnsrtTot[iCentBin],"JER Syst. Unc.","l");
         leg->AddEntry(JES_UnsrtTot[iCentBin],"JES Syst. Unc.","l");
	 leg->AddEntry(Unf_UnsrtTot[iCentBin],"Unfolding Syst. Unc.","l");
         leg->AddEntry(TotSysUncert[iCentBin],"Total Syst. Unc.","l");
       
	 
       }else if(iCentBin != 0){
	 TotSysUncert[iCentBin]->Draw("same");
	 JER_UnsrtTot[iCentBin]->Draw("same");
	 JES_UnsrtTot[iCentBin]->Draw("same");
	 Unf_UnsrtTot[iCentBin]->Draw("same");
       }

       TPaveText *centText = new TPaveText(0.3186728,0.7834916,0.5192901,0.8691983,"brNDC");
       centText->SetTextSize(0.05);
       centText->SetFillColor(0);
       centText->SetBorderSize(0);
       centText->SetShadowColor(0);
       centText->SetTextAlign(33);
       centText->SetTextFont(42);
       if(collisionSys =="pbpbdata"){
	 centText->AddText(Form("Centrality: %s",centBinmap_2015Meas[iCentBin].c_str()));
         centText->Draw("same");
       }
       
    }else if(!plot_totalSysUncert){
      TotSysUncer->cd();
      TotSysUncert[iCentBin]->Write("",TObject::kOverwrite);
    } 

    
    
 }//Centrlity bin loop

  if(plot_totalSysUncert){
    canv->cd(5);
    string plotting_data = "";
    if(collisionSys == "RAA"){
      plotting_data = "ppPbPbdata";
    }else if(collisionSys=="ppdata"){
      plotting_data = "ppdata";
    }else if(collisionSys == "PbPb"){
      plotting_data = "pbpbdata";
    }
    
    Text_Info(plotting_data.c_str(),-1, jetRadius[jetR],0.1307594,0.5177954,0.3270731,0.8691293,"", 0.081,11,etaRange_val);

    canv->cd(10);
    leg->SetBorderSize(0);
    leg->Draw();
    
  }else if(!plot_totalSysUncert){
    TotSysUncer->Close();
  }
  
  /* for(int iCentBin=0; iCentBin < num2015MeasBins; iCentBin++){ */
  
     
  /*     TPaveText *centText = new TPaveText(0.3186728,0.7834916,0.5192901,0.8691983,"brNDC"); */
  /*     centText->SetTextSize(0.05); */
  /*     centText->SetFillColor(0); */
  /*     centText->SetBorderSize(0); */
  /*     centText->SetShadowColor(0); */
  /*     centText->SetTextAlign(33); */
  /*     centText->SetTextFont(42); */
  /*     centText->AddText(Form("Centrality: %s",centBinmap_2015Meas[iCentBin].c_str())); */
  /*     centText->Draw("same"); */

     
    
  /* } */



  
  return 0;

}
