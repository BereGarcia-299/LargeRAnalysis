#include <stdio.h>
#include <vector>
#include "TH1.h"
#include "TFile.h"
#include "/usatlas/u/bereniceg299/data/LargeRJet_Study/NewSourceCode/NewTreeVariables.h"

int sysUncert(int jetR = R4, string collisionType = "RAA" , int jes_uncrt = -1, int jer_uncrt = -1, int unf_uncrt =0,bool debug =true, bool binsJetRate2015Meas = false,bool binsRAA2015Meas=true,int etacut = 28){
 
  string collSys = "pp";
  string tag_for_iterFile = "ppdata";
  if(collisionType=="pbpbdata" || collisionType=="RAA"){
    collSys = "PbPb";
    tag_for_iterFile = "pbpbdata";
  }
  TFile *nominalDataFile;
  
  //This is when you are calcualting RAA Sys. Uncert.
  TFile *nominalDataFile_pp;

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

  cout << "TOT UNCERTIES: " << totUncert << endl;
  
  const int totCentBins = 8;
  int fiducialRegion[Tot_Radii][num2015MeasBins][2]={};
  int fiducialRegion_WOptShpWeights[Tot_Radii][num2015MeasBins][2]={};

  //Nominal iterations for unfilding systematic
  TFile *fileNominalIter_wPtShpWeights;
  TFile *fileNominalIter_woPtShpWeights;

  
  int numIter_wShpWghts[Tot_Radii][num2015MeasBins];
  int numIter_woShpWghts[Tot_Radii][num2015MeasBins];

  
  //This is when you want RAA sytematics 
  TFile *fileNominalIter_wPtShpWeights_pp;
  int numIter_wShpWghts_pp[Tot_Radii][num2015MeasBins];
  int numIter_woShpWghts_pp[Tot_Radii][num2015MeasBins];
  TFile *fileNominalIter_woPtShpWeights_pp;

  string tag = "";
 

  string location_nominal = "";
  string ppOrPbPbTag = "pp";
  string extratag = "";
  
  location_nominal = "/usatlas/u/bereniceg299/data/LargeRJet_Study/NewSourceCode/SysUncert/";
  string location_var = "UnfoldedData";



  //This is for when you want to calculate the systematics for RAA
  string location_pp_nominal = "";
  string location_pp_var = "";


  if(binsJetRate2015Meas == true){
    location_nominal = "/usatlas/u/bereniceg299/data/LargeRJet_Study/NewSourceCode/2015CentBinsFiles";
    location_var = Form("/usatlas/u/bereniceg299/data/LargeRJet_Study/NewSourceCode/2015CentBinsFiles/Unfolded/%sDir",collSys.c_str());
    
  }else if(binsRAA2015Meas==true){
    //pp
    location_pp_nominal = "/usatlas/u/bereniceg299/data/LargeRJet_Study/NewSourceCode/2015CentBinsFiles/rAABins_pp";
    location_pp_var = "/usatlas/u/bereniceg299/data/LargeRJet_Study/NewSourceCode/2015CentBinsFiles/rAABins_pp/Unfolded";
    //Pb+Pb
    location_nominal = "/usatlas/u/bereniceg299/data/LargeRJet_Study/NewSourceCode/2015CentBinsFiles/rAABins_PbPb";
    location_var = "/usatlas/u/bereniceg299/data/LargeRJet_Study/NewSourceCode/2015CentBinsFiles/rAABins_PbPb/Unfolded";
  }

  string binsUsed = "2018MeasBins";
  if(binsJetRate2015Meas)binsUsed = "JetRate2015Bins";
  if(binsRAA2015Meas==true)binsUsed = "2015RAARateBins";

  //NOMINAL DATA FILE
  nominalDataFile = new TFile(Form("%s/nominal/hist_26282548_05302022_Unfolded_%sData_ATLAS_Official_RAA_Binning_17Iters_10000Toys_%s_etaRange%dNominal.root",location_nominal.c_str(),collSys.c_str(),binsUsed.c_str(),etacut),"READ");

  //Nominal data file for pp when calculating RAA Sys.Ucer.
  if(binsRAA2015Meas==true){
    cout << "Grabbing this file for PP: " << Form("%s/nominal/hist_26282548_05302022_Unfolded_ppData_ATLAS_Official_RAA_Binning_17Iters_10000Toys_%s_etaRange%dNominal.root",location_pp_nominal.c_str(),binsUsed.c_str(),etacut) << endl;
    nominalDataFile_pp = new TFile(Form("%s/nominal/hist_26282548_05302022_Unfolded_ppData_ATLAS_Official_RAA_Binning_17Iters_10000Toys_%s_etaRange%dNominal.root",location_pp_nominal.c_str(),binsUsed.c_str(),etacut),"READ");
  }


  //Pick up file with nominal iter. num. for nominal pT dis.
  fileNominalIter_wPtShpWeights = new TFile(Form("%s/nominal/nominalIter_%s_%s_withpTShapeWeights_etaRange%d.root",location_nominal.c_str(),collSys.c_str(),binsUsed.c_str(),etacut),"READ");

  if(binsRAA2015Meas==true){
    fileNominalIter_wPtShpWeights_pp = new TFile(Form("%s/nominal/nominalIter_pp_%s_withpTShapeWeights_etaRange%d.root",location_pp_nominal.c_str(),binsUsed.c_str(),etacut),"READ");
  }
  
  //Pick up file with nominal iter. num. when 2D unfolding matricies are not pT shape weighted
  if(unf_uncrt != -1){
    cout << "UNFOLD HISTO FILE WE WANT: " << Form("%s/nominal/nominalIter_%s_%s_withoutpTShapeWeights_etaRange%d.root",location_nominal.c_str(),collSys.c_str(),binsUsed.c_str(),etacut) << endl;

    fileNominalIter_woPtShpWeights = new TFile(Form("%s/nominal/nominalIter_%s_%s_withoutpTShapeWeights_etaRange%d.root",location_nominal.c_str(),collSys.c_str(),binsUsed.c_str(),etacut));
    
   if(binsRAA2015Meas==true){
      fileNominalIter_woPtShpWeights_pp = new TFile(Form("%s/nominal/nominalIter_pp_%s_withoutpTShapeWeights_etaRange%d.root",location_pp_nominal.c_str(),binsUsed.c_str(),etacut));
    }
  }
   



    for(int iCentBin=0; iCentBin < num2015MeasBins; iCentBin++){
     string tag = "";
     if(collisionType=="pp")tag="data"; 
     if(debug)cout << __LINE__ << endl;
     if(unf_uncrt != -1){
       if(debug)cout<< __LINE__ << endl;
       cout << "get unfolding histo: " << Form("%s%s_Cent_%s_R%d_IterNum",tag_for_iterFile.c_str(),tag.c_str(),centBins_2015Meas[iCentBin].c_str(),JetRadius[jetR]) << endl;
       if(debug)cout<< __LINE__ << endl;
       TH1D *hist1 = (TH1D*)fileNominalIter_woPtShpWeights->Get(Form("%s%s_Cent_%s_R%d_IterNum",tag_for_iterFile.c_str(),tag.c_str(),centBins_2015Meas[iCentBin].c_str(),JetRadius[jetR]));
       if(debug)cout<< __LINE__ << endl;
       numIter_woShpWghts[jetR][iCentBin]=hist1->GetBinContent(1);
       if(debug)cout<< __LINE__ << endl;
       if(binsRAA2015Meas==true){
	 cout << "LETS GRAB THIS: " << Form("ppdata_Cent_%s_R%d_IterNum",centBins_2015Meas[iCentBin].c_str(),JetRadius[jetR]) << endl;
	 numIter_woShpWghts_pp[jetR][iCentBin] = ((TH1D*)fileNominalIter_woPtShpWeights_pp->Get(Form("ppdata_Cent_%s_R%d_IterNum",centBins_2015Meas[iCentBin].c_str(),JetRadius[jetR])))->GetBinContent(1);
	 if(debug)cout<< __LINE__ << endl;
       }
     }
     if(debug)cout << __LINE__ << endl;
     if(binsRAA2015Meas==true){
       if(debug)cout << __LINE__ << endl;
       numIter_wShpWghts_pp[jetR][iCentBin] = ((TH1D*)fileNominalIter_wPtShpWeights_pp->Get(Form("ppdata_Cent_%s_R%d_IterNum",centBins_2015Meas[iCentBin].c_str(),JetRadius[jetR])))->GetBinContent(1); //grabbing nominal iter num. for pp
       if(debug)cout << __LINE__ << endl;
     }
     if(debug)cout << __LINE__ << endl;
     cout << "Grabbing this hist: " << Form("%s%s_Cent_%s_R%d_IterNum",tag_for_iterFile.c_str(),tag.c_str(),centBins_2015Meas[iCentBin].c_str(),JetRadius[jetR]) << endl;
      TH1D *hist2 =  (TH1D*) fileNominalIter_wPtShpWeights->Get(Form("%s%s_Cent_%s_R%d_IterNum",tag_for_iterFile.c_str(),tag.c_str(),centBins_2015Meas[iCentBin].c_str(),JetRadius[jetR]));
      cout << "For this cent bin" << iCentBin << " this is the nominal iteration (w/ pT shape weight): " << hist2->GetBinContent(1) << endl;
      if(debug)cout << __LINE__ << endl;
       
      numIter_wShpWghts[jetR][iCentBin] =hist2->GetBinContent(1); 
      	
      //w/ pT shape weights
      TH1D *hist3 = (TH1D*) fileNominalIter_wPtShpWeights->Get(Form("%s%s_Cent_%s_R%d_StartBinNum",tag_for_iterFile.c_str(),tag.c_str(),centBins_2015Meas[iCentBin].c_str(),JetRadius[jetR]));
      fiducialRegion[jetR][iCentBin][0] = hist3->GetBinContent(1);
      TH1D *hist4  = (TH1D*) fileNominalIter_wPtShpWeights->Get(Form("%s%s_Cent_%s_R%d_EndBinNum",tag_for_iterFile.c_str(),tag.c_str(),centBins_2015Meas[iCentBin].c_str(),JetRadius[jetR]));
      fiducialRegion[jetR][iCentBin][1] = hist4->GetBinContent(1);
      
      if(unf_uncrt==0){
	//w/o pT Shape weights
	TH1D *hist5 = (TH1D*)fileNominalIter_woPtShpWeights->Get(Form("%s%s_Cent_%s_R%d_StartBinNum",tag_for_iterFile.c_str(),tag.c_str(),centBins_2015Meas[iCentBin].c_str(),JetRadius[jetR]));
	fiducialRegion_WOptShpWeights[jetR][iCentBin][0] = hist5->GetBinContent(1);
	TH1D *hist6 = (TH1D*)fileNominalIter_woPtShpWeights->Get(Form("%s%s_Cent_%s_R%d_StartBinNum",tag_for_iterFile.c_str(),tag.c_str(),centBins_2015Meas[iCentBin].c_str(),JetRadius[jetR]));
	fiducialRegion_WOptShpWeights[jetR][iCentBin][1] = hist6->GetBinContent(1);
    
      }
    }
   
    for(int iUncertVal =0; iUncertVal < totUncert; iUncertVal++){
      if(jer_uncrt==0 || jes_uncrt==0){
        cout << "This is iUncertVal: " << iUncertVal << endl;
        sysUncrt[iUncertVal] = new TFile(Form("%s/hist_26282548_05302022_Unfolded_%sData_ATLAS_Official_RAA_Binning_17Iters_10000Toys_%sSys_%d_%s_etaRange%d.root",location_var.c_str(),collSys.c_str(),sysUncertType.c_str(),iUncertVal,binsUsed.c_str(),etacut),"READ");
	if(binsRAA2015Meas==true){
	 
	  if(iUncertVal!=17 && iUncertVal!=20){
	      sysUncert_pp[iUncertVal] = new TFile(Form("%s/hist_26282548_05302022_Unfolded_ppData_ATLAS_Official_RAA_Binning_17Iters_10000Toys_%sSys_%d_%s_etaRange%d.root",location_pp_var.c_str(),sysUncertType.c_str(),iUncertVal,binsUsed.c_str(),etacut),"READ");
	  }else{
	      cout << "THIS BE iUncertVal: "  << iUncertVal << endl; 
	      sysUncert_pp[iUncertVal] = nominalDataFile_pp;
	  }
	}

      }else if(unf_uncrt==0){
        if(iUncertVal==1){
	  extratag="UnfoldingSystematic_FiniteMCStats";
	}else if(iUncertVal==0){
	  extratag="Unfold_pTShpWghtSys";
	}
        
	sysUncrt[iUncertVal] = new TFile(Form("%s/hist_26282548_05302022_Unfolded_%sData_ATLAS_Official_RAA_Binning_17Iters_10000Toys_%s_%s_etaRange%d.root",location_var.c_str(),collSys.c_str(),extratag.c_str(),binsUsed.c_str(),etacut),"READ");
        if(binsRAA2015Meas){
	  sysUncert_pp[iUncertVal] = new TFile(Form("%s/hist_26282548_05302022_Unfolded_ppData_ATLAS_Official_RAA_Binning_17Iters_10000Toys_%s_%s_etaRange%d.root",location_pp_var.c_str(),extratag.c_str(),binsUsed.c_str(),etacut),"READ");
	}
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





  

  
  if(debug)cout <<__LINE__ << endl;
  //File output
  string directory_saved = "";
  if(binsJetRate2015Meas){
    directory_saved =location_nominal;
  }else if(binsRAA2015Meas){
    directory_saved = "/usatlas/u/bereniceg299/data/LargeRJet_Study/NewSourceCode/2015CentBinsFiles/rAABins_systematics";
  }
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
     
      int iter_wPtShp_pp =0;
      int iter_pp =0; 
      
      iter_wPtShp = numIter_wShpWghts[jetR][iCentBin];
      if(binsRAA2015Meas==true){
	iter_wPtShp_pp = numIter_wShpWghts_pp[jetR][iCentBin];
	cout << "Saving this value in variable iter_wPtShp_pp: " << iter_wPtShp_pp << endl;
      }
      if(unf_uncrt==0 && iSysUncert==0){
	  
	  iter = numIter_woShpWghts[jetR][iCentBin]; 
	  if(binsRAA2015Meas==true){
	   
	    iter_pp = numIter_woShpWghts_pp[jetR][iCentBin];
	    cout << "for the unfolding systematic (pt shp wght): " <<  iter_pp << endl;
	  }
      }else{
	if(binsRAA2015Meas==true)iter_pp = numIter_wShpWghts[jetR][iCentBin];
	iter = numIter_wShpWghts[jetR][iCentBin];
      }
      
      cout << "This is collisionType: " << collisionType.c_str() << endl;
    
      if((collisionType=="pbpbdata" || collisionType=="PbPb") && unf_uncrt!=0){
	iter_wPtShp = numIter_wShpWghts[jetR][iCentBin];
	cout << "This is the iteration number: " << iter_wPtShp << endl;
      }
      //systematic
      if(collisionType=="pbpbdata")collisionType = "PbPb";
      
      
      
      TH1D *sys_hist = (TH1D*) sysUncrt[iSysUncert]->Get(Form("Unfolded_%sData_R%d_%dIter_%s",collSys.c_str(),JetRadius[jetR],iter,centBins_2015Meas[iCentBin].c_str()));
      cout << "Grabbed this histogram: " << Form("Unfolded_%sData_R%d_%dIter_%s",collSys.c_str(),JetRadius[jetR],iter,centBins_2015Meas[iCentBin].c_str()) << endl;
      //cout << "getting bin content of number of bin 11: " << sys_hist->GetBinContent(11) << endl;
      
      //Nominal measurement
      TH1D *nom = (TH1D*)nominalDataFile->Get(Form("Unfolded_%sData_R%d_%dIter_%s",collSys.c_str(),JetRadius[jetR],iter_wPtShp,centBins_2015Meas[iCentBin].c_str()));
      cout << "Grabbing the nominal hist: " << Form("Unfolded_%sData_R%d_%dIter_%s",collSys.c_str(),JetRadius[jetR],iter_wPtShp,centBins_2015Meas[iCentBin].c_str()) << endl;
      if(debug)cout << __LINE__ << endl;
      TH1D *nom_pp;
      TH1D *sys_hist_pp;
      if(debug)cout << __LINE__<< endl;
      if(binsRAA2015Meas){
	if(debug)cout << __LINE__<< endl;
	cout << "grab iteration number: " << iter_wPtShp_pp << endl;
	cout << "Grabbing this PP histo: " << Form("Unfolded_ppData_R%d_%dIter_%s",JetRadius[jetR],iter_wPtShp_pp,centBins_2015Meas[iCentBin].c_str()) << endl;
	nom_pp = (TH1D*)nominalDataFile_pp->Get(Form("Unfolded_ppData_R%d_%dIter_%s",JetRadius[jetR],iter_wPtShp_pp,centBins_2015Meas[iCentBin].c_str()));
	if(debug)cout << __LINE__<< endl;
	cout << "VARIED PP SPECTRA HIST: " << Form("Unfolded_ppData_R%d_%dIter_%s",JetRadius[jetR],iter_pp,centBins_2015Meas[iCentBin].c_str()) << endl;
	cout << "This is for sys num: " << iSysUncert << endl;
	cout << "Grab this pp varied hist: " << Form("Unfolded_ppData_R%d_%dIter_%s",JetRadius[jetR],iter_pp,centBins_2015Meas[iCentBin].c_str()) << endl;
	sys_hist_pp = (TH1D*)sysUncert_pp[iSysUncert]->Get(Form("Unfolded_ppData_R%d_%dIter_%s",JetRadius[jetR],iter_pp,centBins_2015Meas[iCentBin].c_str()));
	if(debug)cout << __LINE__<< endl;
      }
      if(debug)cout << __LINE__<< endl;
      
      sysroot->cd();
      
      //scaling factor for Pb+Pb jet rates
      double scaling_fac = -1.0;
      //scaling factor for pp jet rates when plotting RAA
      if(collisionType=="RAA"){
	if(debug)cout << __LINE__<< endl;
	//pp Nominal pT spectra
	nom_pp->Scale(1/ppDataLumiVals[jetR]);
	if(debug)cout << __LINE__<< endl;
	nom_pp->Scale(1.,"width");
	if(debug)cout << __LINE__<< endl;
	nom_pp->Scale(1/etaRange);
	if(debug)cout << __LINE__<< endl;
	//pp Variated pT spectra
	sys_hist_pp->Scale(1/ppDataLumiVals[jetR]);
	if(debug)cout << __LINE__<< endl;
	sys_hist_pp->Scale(1.,"width");
	if(debug)cout << __LINE__<< endl;
	sys_hist_pp->Scale(1/etaRange);
	if(debug)cout << __LINE__<< endl;
      }


      if(collisionType=="pbpbdata" || collisionType=="RAA"){
        scaling_fac = tAA_2015map[iCentBin]*numEventsMinBias; //Scaling factor for Pb+Pb
      }else if(collisionType=="pp"){
	scaling_fac = ppDataLumiVals[jetR]; //Scaling factor for pp
      }
      
      //Scaling nominal
      nom->Scale(1/scaling_fac);
      nom->Scale(1.,"width");
      nom->Scale(1/etaRange);

      //Scaling variated
      sys_hist->Scale(1/scaling_fac);
      sys_hist->Scale(1.,"width");
      sys_hist->Scale(1/etaRange);
      
      if(collisionType=="RAA"){
	nom->Divide(nom_pp);
	sys_hist->Divide(sys_hist_pp);
      }
      
      sys_hist->Draw();
      
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
    }else if(binsRAA2015Meas){
      array = rAA_2015Bins[iCentBin];
      numBins = bins2015;
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
