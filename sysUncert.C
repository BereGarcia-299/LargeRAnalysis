#include <stdio.h>
#include <vector>
#include "TH1.h"
#include "TFile.h"
#include "/usatlas/u/bereniceg299/data/LargeRJet_Study/NewSourceCode/NewTreeVariables.h"

int sysUncert(int jetR = R4, string collisionType = "RAA" , int jes_uncrt = 0, int jer_uncrt = -1, int unf_uncrt =-1,bool debug =true, bool binsJetRate2015Meas = true,bool binsRAA2015Meas=true,int etacut = 28){
 
  string collSys = "pp";
  if(collisionType=="pbpbdata")collSys = "PbPb";
  
  TFile *nominalDataFile;
  
  string sysUncertType = "JES";
  
  TFile *sysUncrt[40]; //This will be used for either pp or Pb+Pb 
  //If RAA_Uncert true then it will use ot for Pb+Pb and use the bottom array for pp
  TFile *sysUncert_pp[40];

  if(debug)cout << __LINE__ << endl;

  int totUncert = 21; //total number of JES uncertanties
  if(jer_uncrt== 0){
    totUncert=9; //total number of JER uncertanties 
    sysUncertType = "JER";
  }else if(unf_uncrt==0){
    totUncert=2;
    sysUncertType = "Unfolding";
  }
  
  if(collisionType == "pp"){
    if(jes_uncrt!=-1)totUncert=20;
    if(jer_uncrt!=-1)totUncert=9;
  }
  
  const int totCentBins = 8;

  int fiducialRegion_WOptShpWeights[Tot_Radii][num2015MeasBins][2]={};

  //Nominal iterations for unfilding systematic
  TFile *fileNominalIter_wPtShpWeights;
  TFile *fileNominalIter_woPtShpWeights;
  
  int numIter_wShpWghts[Tot_Radii][num2015MeasBins];
  int numIter_woShpWghts[Tot_Radii][num2015MeasBins];

  string tag = "";
 

  string location_nominal = "";
  string ppOrPbPbTag = "pp";
  string extratag = "";
  
  location_nominal = "/usatlas/u/bereniceg299/data/LargeRJet_Study/NewSourceCode/SysUncert/";
  string location_var = "UnfoldedData";

  if(binsJetRate2015Meas == true){
    location_nominal = "/usatlas/u/bereniceg299/data/LargeRJet_Study/NewSourceCode/2015CentBinsFiles";
    location_var = Form("/usatlas/u/bereniceg299/data/LargeRJet_Study/NewSourceCode/2015CentBinsFiles/Unfolded/%sDir",collSys.c_str());
    
  }else if(binsRAA2015Meas==true){
    

  }
  string binsUsed = "2018MeasBins";
  if(binsJetRate2015Meas)binsUsed = "JetRate2015Bins";
  
  fileNominalIter_wPtShpWeights = new TFile(Form("%s/nominal/nominalIter_%s_%s_withpTShapeWeights_etaRange%d.root",location_nominal.c_str(),collSys.c_str(),binsUsed.c_str(),etacut),"READ");
  if(unf_uncrt != -1)fileNominalIter_woPtShpWeights = new TFile(Form("%s/nominal/nominalIter_%s_%s_withoutpTShapeWeights_etaRange%d.root",location_nominal.c_str(),collSys.c_str(),binsUsed.c_str(),etacut));

   



    for(int iCentBin=0; iCentBin < num2015MeasBins; iCentBin++){
     string tag = "";
     if(collisionType=="pp")tag="data"; 
     if(debug)cout << __LINE__ << endl;
     if(unf_uncrt != -1){
       if(debug)cout<< __LINE__ << endl;
       TH1D *hist1 = (TH1D*)fileNominalIter_woPtShpWeights->Get(Form("%s%s_Cent_%s_R%d_IterNum",collisionType.c_str(),tag.c_str(),centBins_2015Meas[iCentBin].c_str(),JetRadius[jetR]));
      
       numIter_woShpWghts[jetR][iCentBin]=hist1->GetBinContent(1);
        cout << "For this cent bin" << iCentBin << "this is the nominal iteration (w/o pT shape weight): " << hist1->GetBinContent(1) << endl;
      }
     if(debug)cout<< __LINE__ << endl;
     cout << "Grabbing this hist: " << Form("%s%s_Cent_%s_R%d_IterNum",collisionType.c_str(),tag.c_str(),centBins_2015Meas[iCentBin].c_str(),JetRadius[jetR]) << endl;
      TH1D *hist2 =  (TH1D*) fileNominalIter_wPtShpWeights->Get(Form("%s%s_Cent_%s_R%d_IterNum",collisionType.c_str(),tag.c_str(),centBins_2015Meas[iCentBin].c_str(),JetRadius[jetR]));
      cout << "For this cent bin" << iCentBin << " this is the nominal iteration (w/ pT shape weight): " << hist2->GetBinContent(1) << endl;
      if(debug)cout<< __LINE__ << endl;
       
      numIter_wShpWghts[jetR][iCentBin] =hist2->GetBinContent(1); 
      	
      //w/ pT shape weights
      TH1D *hist3 = (TH1D*) fileNominalIter_wPtShpWeights->Get(Form("%s%s_Cent_%s_R%d_StartBinNum",collisionType.c_str(),tag.c_str(),centBins_2015Meas[iCentBin].c_str(),JetRadius[jetR]));
      fiducialRegion[jetR][iCentBin][0] = hist3->GetBinContent(1);
      TH1D *hist4  = (TH1D*) fileNominalIter_wPtShpWeights->Get(Form("%s%s_Cent_%s_R%d_EndBinNum",collisionType.c_str(),tag.c_str(),centBins_2015Meas[iCentBin].c_str(),JetRadius[jetR]));
      fiducialRegion[jetR][iCentBin][1] = hist4->GetBinContent(1);
      
      if(unf_uncrt==0){
	//w/o pT Shape weights
	TH1D *hist5 = (TH1D*)fileNominalIter_woPtShpWeights->Get(Form("%s%s_Cent_%s_R%d_StartBinNum",collisionType.c_str(),tag.c_str(),centBins_2015Meas[iCentBin].c_str(),JetRadius[jetR]));
	fiducialRegion_WOptShpWeights[jetR][iCentBin][0] = hist5->GetBinContent(1);
	TH1D *hist6 = (TH1D*)fileNominalIter_woPtShpWeights->Get(Form("%s%s_Cent_%s_R%d_StartBinNum",collisionType.c_str(),tag.c_str(),centBins_2015Meas[iCentBin].c_str(),JetRadius[jetR]));
	fiducialRegion_WOptShpWeights[jetR][iCentBin][1] = hist6->GetBinContent(1);
      }
    }
   
    for(int iUncertVal =0; iUncertVal < totUncert; iUncertVal++){
      if(jer_uncrt==0 || jes_uncrt==0){
        
        sysUncrt[iUncertVal] = new TFile(Form("%s/hist_26282548_05302022_Unfolded_%sData_ATLAS_Official_RAA_Binning_17Iters_10000Toys_%sSys_%d_%s_etaRange%d.root",location_var.c_str(),collSys.c_str(),sysUncertType.c_str(),iUncertVal,binsUsed.c_str(),etacut),"READ");
      }else if(unf_uncrt==0){
        if(iUncertVal==1){
	  extratag="UnfoldingSystematic_FiniteMCStats";
	}else if(iUncertVal==0){
	  extratag="Unfold_pTShpWghtSys";
	}
        sysUncrt[iUncertVal] = new TFile(Form("%s/hist_26282548_05302022_Unfolded_%sData_ATLAS_Official_RAA_Binning_17Iters_10000Toys_%s_%s_etaRange%d.root",location_var.c_str(),collSys.c_str(),extratag.c_str(),binsUsed.c_str(),etacut),"READ");
        cout << __LINE__ << endl;
      }
    }//Uncertainty Loop
 
  
  

  
 

  
  //Number of bins
  //int pTBinsTotR[2][8] = {{13,14,14,13,13,14,13,13},{22,22,22,19,19,20,19,19}};
  int pTBinsTotR[2][8] = {{11,12,12,11,11,12,11,11},{20,20,20,17,17,18,17,17}};
  //Luminosity Numbers
  float LumNumPbPbData[] = {1.72317,1.72317}; //(nb^-1) R=1.0 HLT_j150_a10_ion_L1J50 / R=0.4 HLT_j85_ion_L1J30
  
  float ppDataLumiVals[] = {132.199e+3,132.199e+3};

  //numEvebtsMinBias 
  double numEventsMinBias = 7.383e+8*LumNumPbPbData[jetR];
  

  //eta range
  double etaRange = 1.5*2;
  if(etacut==28)etaRange = 2.8*2;

  //Binning for this analysis
  std::vector<std::vector<double>> bins;
  bins.reserve(8);

  if(jetR==R10 && !binsJetRate2015Meas){
    //R = 1.0

    bins.emplace_back(std::vector<double>{241, 262, 294, 336, 384, 439, 502, 573, 656, 750, 857, 980, 1120,1300}); //0-10%
    bins.emplace_back(std::vector<double>{220,241, 262, 294, 336, 384, 439, 502, 573, 656, 750, 857, 980, 1120, 1300}); //10-20%
    bins.emplace_back(std::vector<double>{211, 228, 245, 262, 294, 336, 384, 439, 502, 573, 656, 750, 857, 980, 1120}); //20-30%
    bins.emplace_back(std::vector<double>{177, 198, 230, 262, 294, 336, 384, 439, 502, 573, 656, 750, 857, 1000}); //30-40%
    bins.emplace_back(std::vector<double>{177,198, 230, 262, 294, 336, 384, 439, 502, 573, 656, 750, 857, 1000}); //40-50%
    bins.emplace_back(std::vector<double>{177, 198, 220, 241, 262, 294, 336, 384, 439, 502, 573, 656, 750, 857, 1000}); //50-60%
    bins.emplace_back(std::vector<double>{177, 197, 220, 241, 262, 294, 336, 384, 439, 502, 573, 656, 750, 1000}); //60-70%
    bins.emplace_back(std::vector<double>{177, 197, 220, 241, 262, 294, 336, 384, 439, 502, 573, 656, 750, 1000}); //70-80%~
  }else if(jetR==R4 && !binsJetRate2015Meas){
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





  nominalDataFile = new TFile(Form("%s/nominal/hist_26282548_05302022_Unfolded_%sData_ATLAS_Official_RAA_Binning_17Iters_10000Toys_%s_etaRange%dNominal.root",location_nominal.c_str(),collSys.c_str(),binsUsed.c_str(),etacut),"READ");

  

  
  if(debug)cout <<__LINE__ << endl;
  //File output
  string directory_saved = "";
  if(binsJetRate2015Meas)directory_saved =location_nominal;
  TFile *sysroot = new TFile(Form("%s/systematics/systematics_%s%s_%s_R%d_%s.root",directory_saved.c_str(),collisionType.c_str(),tag.c_str(),binsUsed.c_str(),JetRadius[jetR],sysUncertType.c_str()),"RECREATE");



  double sqrdSum[totCentBins][40] = {};


  for(int iCentBin =0; iCentBin < totCentBins; iCentBin++){
    double *array = NULL;
    int numBins = pTBinsTotR[jetR][iCentBin];
    if(!binsJetRate2015Meas){
      array = bins.at(iCentBin).data();
    }else if(binsJetRate2015Meas){
      array = jetRateMC;
      numBins = jetRateMCBins;
    }
   

    for(int iSysUncert = 0; iSysUncert < totUncert; iSysUncert++){ 
      
      
     
      int iter_wPtShp =0;
      int iter = 0;
     
      iter_wPtShp = numIter_wShpWghts[jetR][iCentBin];
      if(unf_uncrt==0 && iSysUncert==0){
	  
	  iter = numIter_woShpWghts[jetR][iCentBin]; 
	
      }else{
	iter = numIter_wShpWghts[jetR][iCentBin];
      }

      cout << "This is collisionType: " << collisionType.c_str() << endl;
      //cout << "is this true: " << (collisionType.c_str()=="pbpbdata") << endl;
      if((collisionType=="pbpbdata" || collisionType=="PbPb") && unf_uncrt!=0){
	iter_wPtShp = numIter_wShpWghts[jetR][iCentBin];
	cout << "This is the iteration number: " << iter_wPtShp << endl;
      }
      //systematic
      if(collisionType=="pbpbdata")collisionType = "PbPb";
      
      
      
      TH1D *sys_hist = (TH1D*) sysUncrt[iSysUncert]->Get(Form("Unfolded_%sData_R%d_%dIter_%s",collisionType.c_str(),JetRadius[jetR],iter,centBins_2015Meas[iCentBin].c_str()));
      cout << "Grabbed this histogram: " << Form("Unfolded_%sData_R%d_%dIter_%s",collisionType.c_str(),JetRadius[jetR],iter,centBins_2015Meas[iCentBin].c_str()) << endl;
      //cout << "getting bin content of number of bin 11: " << sys_hist->GetBinContent(11) << endl;
      
      //Nominal measurement
      TH1D *nom = (TH1D*)nominalDataFile->Get(Form("Unfolded_%sData_R%d_%dIter_%s",collisionType.c_str(),JetRadius[jetR],iter_wPtShp,centBins_2015Meas[iCentBin].c_str()));
      cout << "Grabbing the nominal hist: " << Form("Unfolded_%sData_R%d_%dIter_%s",collisionType.c_str(),JetRadius[jetR],iter_wPtShp,centBins_2015Meas[iCentBin].c_str()) << endl;
      
      if(debug)cout << __LINE__ << endl;
      
      

      

      sysroot->cd();
      
      if(collisionType=="pbpbdata"){
        if(debug)cout << __LINE__ << endl;	
	sys_hist->Scale(1/numEventsMinBias);
      }else if(collisionType=="pp"){
	 if(debug)cout << __LINE__ << endl;
	sys_hist->Scale(1/ppDataLumiVals[jetR]);
      }
      if(debug)cout << __LINE__ << endl;
      sys_hist->Scale(1.,"width");
      sys_hist->Scale(1/etaRange);
      if(debug)cout << __LINE__ << endl;
      if(collisionType=="pbpbdata")sys_hist->Scale(1/tAA_2015map[iCentBin]);
      sys_hist->Draw();
      if(debug)cout << __LINE__ << endl;
      if(collisionType=="pbpbdata"){
	if(debug)cout << __LINE__ << endl;
	nom->Scale(1/numEventsMinBias);
	nom->Scale(1/tAA_2015map[iCentBin]);
	if(debug)cout << __LINE__ << endl;
      }

      if(debug)cout << __LINE__ << endl;
      if(collisionType=="pp"){
	nom->Scale(1/ppDataLumiVals[jetR]);
      }
      if(debug)cout <<__LINE__ << endl;

      
      
      nom->Scale(1.,"width");
      nom->Scale(1/etaRange);
      

      nom->Write("",TObject::kOverwrite);
      sys_hist->Write(Form("Unfolded_%sData_R%d_%dIter_%s_%s_%d",collisionType.c_str(),JetRadius[jetR],iter_wPtShp,centBins_2015Meas[iCentBin].c_str(),sysUncertType.c_str(),iSysUncert),TObject::kOverwrite);
     
      TH1D *plot = new TH1D(Form("R%d_%s_%d_Cent_%s",JetRadius[jetR],sysUncertType.c_str(),iSysUncert,centBins_2015Meas[iCentBin].c_str()),Form("R%d_%s_%d_Cent_%s",JetRadius[jetR],sysUncertType.c_str(),iSysUncert,centBins_2015Meas[iCentBin].c_str()),numBins,array);
      if(debug)cout << __LINE__ << endl;
      
      for(int iBin = fiducialRegion[jetR][iCentBin][startAndend::start]; iBin < fiducialRegion[jetR][iCentBin][startAndend::end]+1; iBin++){
	if(debug)cout << __LINE__ << endl;
	cout << "This is the iBin: " << iBin << endl;
	cout << "This is the bin center: " << nom->GetBinCenter(iBin) << endl;
	plot->SetBinContent(iBin,abs(nom->GetBinContent(iBin) - sys_hist->GetBinContent(iBin))*100/nom->GetBinContent(iBin));
	
	cout << "Difference between nominal and sys: " << (nom->GetBinContent(iBin) - sys_hist->GetBinContent(iBin)) << endl;
	cout << "with pTshape value: " << sys_hist->GetBinContent(iBin) << endl;
	cout << "nominal value: " << nom->GetBinContent(iBin) << endl;
	cout << "Systematic uncert: " << (double)abs(nom->GetBinContent(iBin) - sys_hist->GetBinContent(iBin))/nom->GetBinContent(iBin) << endl;
	
	sqrdSum[iCentBin][iBin] = sqrdSum[iCentBin][iBin] + pow((nom->GetBinContent(iBin) - sys_hist->GetBinContent(iBin))*100/nom->GetBinContent(iBin),2);
	
	  
      }//Bin Loop

      if(debug)cout << __LINE__ << endl;

    
	
	plot->Write(Form("R%d_%s_%d_Cent_%s",JetRadius[jetR],sysUncertType.c_str(),iSysUncert,centBins_2015Meas[iCentBin].c_str()),TObject::kOverwrite);
    


    }//Systematic Loop

  }//CentBin Loop


  for(int iCentBin =0; iCentBin < totCentBins; iCentBin++){
    double *array = NULL;
    int numBins = pTBinsTotR[jetR][iCentBin];
    if(!binsJetRate2015Meas){
      array = bins.at(iCentBin).data();
    }else if(binsJetRate2015Meas){
      array = jetRateMC;
      numBins = jetRateMCBins;
    }
   

    TH1D *quadsum = new TH1D(Form("QS_R%d_%s_Cent_%s",JetRadius[jetR],sysUncertType.c_str(),centBins_2015Meas[iCentBin].c_str()),Form("QS_R%d_%s_Cent_%s",JetRadius[jetR],sysUncertType.c_str(),centBins_2015Meas[iCentBin].c_str()),numBins,array);
    
    
    
    for(int iBin = fiducialRegion[jetR][iCentBin][startAndend::start]; iBin < fiducialRegion[jetR][iCentBin][startAndend::end]+1; iBin++){
      
      quadsum->SetBinContent(iBin,sqrt(sqrdSum[iCentBin][iBin]));
    }//bin loop

    quadsum->Write("",TObject::kOverwrite);
    
  }


  sysroot->Close();

  return 0;
}
