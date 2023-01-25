#include <stdio.h>
#include <vector>
#include "TH1.h"
#include "TFile.h"
#include "/usatlas/u/bereniceg299/data/LargeRJet_Study/NewSourceCode/NewTreeVariables.h"

int rAA_sysUncert(int jetR = R4, string collisionType = "RAA" , int jes_uncrt = -1, int jer_uncrt =-1, int unf_uncrt = 0,bool debug =true){
 
  //Set this variable to true if you want to caluclate the uncertainties for RAA
  //Please set 'collisionType' to 'RAA'
  bool RAA_Uncert = true;
  
  TFile *nominalDataFilePbPb;
  TFile *nominalDataFilepp;

  string sysUncertType = "JES";
  
  TFile *sysUncrt_PbPb[40]; //This will be used for either pp or Pb+Pb 
  //If RAA_Uncert true then it will use ot for Pb+Pb and use the bottom array for pp
  TFile *sysUncrt_pp[40];

  

  int totUncert = 21; //total number of JES uncertanties
  if(jer_uncrt== 0){
    totUncert=9; //total number of JER uncertanties 
    sysUncertType = "JER";
  }else if(unf_uncrt==0){
    totUncert=2;
    sysUncertType = "Unfolding";
  }
  
  if(collisionType == "pp"){
    if(jes_uncrt!=-1)totUncert=19;
    if(jer_uncrt!=-1)totUncert=9;
  }
  
  const int totCentBins = 8;

  //Nuber of iterations
  //Nominal
  int numIter_PbPb[Tot_Radii][num2015MeasBins]= {{2,3,4,2,1,1,1,1},{4,2,5,3,3,2,2,3}}; //Jet Yields in Pb+Pb
  int numIter_pp[Tot_Radii][num2015MeasBins]= {{2,3,4,2,1,1,1,1},{5,5,2,4,4,4,4,4}}; //Cross-sections in pp  
  if(debug)cout << __LINE__ << endl;
    //Fudicial Region
  int fiducialRegion[Tot_Radii][num2015MeasBins][2] = {
    {{5,13},{5,15},{5,14},{5,13},{5,13},{5,14},{5,13},{5,13}},
    {{5,21},{5,21},{5,21},{5,18},{5,18},{5,19},{5,18},{5,18}}
  };

   if(debug)cout << __LINE__ << endl; 
  TFile *nominal_Iter_PbPbWptShp;
  TFile *nominal_Iter_PbPbWOptShp;

  TFile *nominal_Iter_ppWptShp;
  TFile *nominal_Iter_ppWOptShp;
						

  int numIter_PbPb_nominal[Tot_Radii][totCentBins] ={};
  int numIter_PbPb_varied[Tot_Radii][totCentBins] ={};
  
  int numIter_pp_nominal[Tot_Radii][totCentBins] ={};
  int numIter_pp_varied[Tot_Radii][totCentBins] ={};


  string FR[2] = {"Start","End"};

  nominal_Iter_PbPbWptShp =new TFile("nominalIter_PbPb_withpTShapeWeights.root","READ");
  nominal_Iter_ppWptShp = new TFile("nominalIter_pp_withpTShapeWeights.root","READ");
  if(unf_uncrt==0){
     if(debug)cout << __LINE__ << endl;
    
    nominal_Iter_PbPbWOptShp = new TFile("nominalIter_PbPb_withoutpTShapeWeights.root","READ");
     if(debug)cout << __LINE__ << endl;
    
    nominal_Iter_ppWOptShp = new TFile("nominalIter_pp_withoutpTShapeWeights.root","READ");
  }						

  

  for(int iCentBin =0; iCentBin <  totCentBins; iCentBin++){
     if(debug)cout << __LINE__ << endl;
     cout << "Trying to grab this file: " << endl;
     cout << Form("pbpbdata_Cent_%s_R%d_IterNum",centBins_2015Meas[iCentBin].c_str(),jetRadius[jetR]) << endl;
       
    numIter_PbPb_nominal[jetR][iCentBin] = ((TH1D*)nominal_Iter_PbPbWptShp->Get(Form("pbpbdata_Cent_%s_R%d_IterNum",centBins_2015Meas[iCentBin].c_str(),jetRadius[jetR])))->GetBinContent(1);
    
    if(debug)cout << __LINE__ << endl;
    if(unf_uncrt==0){
      numIter_PbPb_varied[jetR][iCentBin] =  ((TH1D*)nominal_Iter_PbPbWOptShp->Get(Form("pbpbdata_Cent_%s_R%d_IterNum",centBins_2015Meas[iCentBin].c_str(),jetRadius[jetR])))->GetBinContent(1);
      numIter_pp_varied[jetR][iCentBin] = ((TH1D*)nominal_Iter_ppWOptShp->Get(Form("ppdata_Cent_%s_R%d_IterNum",centBins_2015Meas[iCentBin].c_str(),jetRadius[jetR])))->GetBinContent(1);
    } 
    if(debug)cout << __LINE__ << endl;
    numIter_pp_nominal[jetR][iCentBin] = ((TH1D*)nominal_Iter_ppWptShp->Get(Form("ppdata_Cent_%s_R%d_IterNum",centBins_2015Meas[iCentBin].c_str(),jetRadius[jetR])))->GetBinContent(1);
    
     if(debug)cout << __LINE__ << endl;
    for(int iFudicialReg =0; iFudicialReg < 2; iFudicialReg++){
       if(debug)cout << __LINE__ << endl;
       cout << "This is the histo that we are trying to grab: " << Form("ppdata_Cent_%s_R%d_%sBinNum",centBins_2015Meas[iCentBin].c_str(),jetRadius[jetR],FR[iFudicialReg].c_str()) << endl;
       fiducialRegion[jetR][iCentBin][iFudicialReg] = ((TH1D*)nominal_Iter_ppWptShp->Get(Form("ppdata_Cent_%s_R%d_%sBinNum",centBins_2015Meas[iCentBin].c_str(),jetRadius[jetR],FR[iFudicialReg].c_str())))->GetBinContent(1); 
      if(debug)cout << __LINE__ << endl;
    }//start and end bin
  }//Cent Bin Loop

  



  //
   

  
  //Number of bins
  int pTBinsTotR[2][8] = {{13,14,14,13,13,14,13,13},{22,22,22,19,19,20,19,19}};
  
  //Luminosity Numbers
  float LumNumPbPbData[] = {1.72317,1.72317}; //(nb^-1) R=1.0 HLT_j150_a10_ion_L1J50 / R=0.4 HLT_j85_ion_L1J30
  
  float ppDataLumiVals[] = {132.199e+3,132.199e+3};

  //numEvebtsMinBias 
  double numEventsMinBias = 7.383e+8*LumNumPbPbData[jetR];
  

  //eta range
  double etaRange = 1.5*2;


  //Binning for this analysis
  std::vector<std::vector<double>> bins;
  bins.reserve(8);

  if(jetR==R10){
    //R = 1.0
     if(debug)cout << __LINE__ << endl;
    bins.emplace_back(std::vector<double>{241, 262, 294, 336, 384, 439, 502, 573, 656, 750, 857, 980, 1120,1300}); //0-10%
    bins.emplace_back(std::vector<double>{220,241, 262, 294, 336, 384, 439, 502, 573, 656, 750, 857, 980, 1120, 1300}); //10-20%
    bins.emplace_back(std::vector<double>{211, 228, 245, 262, 294, 336, 384, 439, 502, 573, 656, 750, 857, 980, 1120}); //20-30%
    bins.emplace_back(std::vector<double>{177, 198, 230, 262, 294, 336, 384, 439, 502, 573, 656, 750, 857, 1000}); //30-40%
    bins.emplace_back(std::vector<double>{177,198, 230, 262, 294, 336, 384, 439, 502, 573, 656, 750, 857, 1000}); //40-50%
    bins.emplace_back(std::vector<double>{177, 198, 220, 241, 262, 294, 336, 384, 439, 502, 573, 656, 750, 857, 1000}); //50-60%
    bins.emplace_back(std::vector<double>{177, 197, 220, 241, 262, 294, 336, 384, 439, 502, 573, 656, 750, 1000}); //60-70%
    bins.emplace_back(std::vector<double>{177, 197, 220, 241, 262, 294, 336, 384, 439, 502, 573, 656, 750, 1000}); //70-80%~
  }else if(jetR==R4){
    //R=0.4
    bins.emplace_back(std::vector<double>{66,83,100, 117, 138, 158, 178, 199, 220, 241, 262, 294, 336, 384, 439, 502, 573, 656, 750, 857, 980, 1120,1300}); //0-10%
    bins.emplace_back(std::vector<double>{66,83,100, 117, 138, 158, 178, 199, 220, 241, 262, 294, 336, 384, 439, 502, 573, 656, 750, 857, 980, 1120, 1300}); //10-20%
    bins.emplace_back(std::vector<double>{75,92,109, 126, 143, 160, 177, 194, 211, 228, 245, 262, 294, 336, 384, 439, 502, 573, 656, 750, 857, 980, 1120}); //20-30%
    bins.emplace_back(std::vector<double>{66,83,100, 117, 138, 158, 177, 198, 230, 262, 294, 336, 384, 439, 502, 573, 656, 750, 857, 1000}); //30-40%
    bins.emplace_back(std::vector<double>{66,83,100, 117, 138, 158, 177, 198, 230, 262, 294, 336, 384, 439, 502, 573, 656, 750, 857, 1000}); //40-50%
    bins.emplace_back(std::vector<double>{66,83,100, 117, 138, 158, 177, 198, 220, 241, 262, 294, 336, 384, 439, 502, 573, 656, 750, 857, 1000}); //50-60%
    bins.emplace_back(std::vector<double>{66,83,100, 117, 138, 158, 177, 197, 220, 241, 262, 294, 336, 384, 439, 502, 573, 656, 750, 1000}); //60-70%
    bins.emplace_back(std::vector<double>{66,83,100, 117, 138, 158, 177, 197, 220, 241, 262, 294, 336, 384, 439, 502, 573, 656, 750, 1000}); //70-80%~
  }








  string location = "";
  string ppOrPbPbTag = "pp";
  location = "/usatlas/u/bereniceg299/data/LargeRJet_Study/NewSourceCode/SysUncert/nominal/";

  //Nominal Values
    nominalDataFilePbPb = new TFile(Form("%shist_26282548_05302022_Unfolded_PbPbData_Eta15_ATLAS_Official_RAA_Binning_17Iters_10000Toys_NewpTBinsptSHapeWeightsAppliedTo2DMatricies_NewpTShpWghts.root",location.c_str()),"READ");
    
    nominalDataFilepp = new TFile(Form("%shist_26282548_05302022_Unfolded_ppData_Eta15_ATLAS_Official_RAA_Binning_17Iters_10000Toys_NewpTBinsptSHapeWeightsAppliedTo2DMatricies_NewpTShpWghts.root",location.c_str()),"READ");
     
    //Unfolded w/ JES/JER Systematic Variations
    for(int iUncertVal =0; iUncertVal < totUncert; iUncertVal++){
      if(sysUncertType != "Unfolding"){
	if(iUncertVal <= totUncert-3){
	
	  sysUncrt_pp[iUncertVal] = new TFile(Form("UnfoldedData/hist_26282548_05302022_Unfolded_ppData_Eta15_ATLAS_Official_RAA_Binning_17Iters_10000Toys_NewpTBins%sSys_%dptSHapeWeightsAppliedTo2DMatricies_NewpTShpWghts.root",sysUncertType.c_str(),iUncertVal),"READ");
	}else{
	  sysUncrt_pp[iUncertVal] = nominalDataFilepp;
	}
      }else if(sysUncertType == "Unfolding"){
	string extratag ="";
	if(iUncertVal==1)extratag="_UnfoldingSystematic_FiniteMCStats";
	sysUncrt_pp[iUncertVal] = new TFile(Form("UnfoldedData/hist_26282548_05302022_Unfolded_ppData_Eta15_ATLAS_Official_RAA_Binning_17Iters_10000Toys_NewpTBins%s_NewpTShpWghts.root",extratag.c_str()),"READ");
      }

    }//pp Uncertainty Loop
    int counter =0;
    for(int iUncertValPbPb =0; iUncertValPbPb < totUncert; iUncertValPbPb++){
      if((iUncertValPbPb!=17 || iUncertValPbPb!=20) && sysUncertType != "Unfolding"){
	sysUncrt_PbPb[counter] = new TFile(Form("UnfoldedData/hist_26282548_05302022_Unfolded_PbPbData_Eta15_ATLAS_Official_RAA_Binning_17Iters_10000Toys_NewpTBins%sSys_%dptSHapeWeightsAppliedTo2DMatricies_NewpTShpWghts.root",sysUncertType.c_str(),counter),"READ");
	
	counter++;
      }else if(sysUncertType == "Unfolding"){
	string extratag1="";
	cout << "This is the number of counter: " << counter << endl;
	if(iUncertValPbPb==1)extratag1="_UnfoldingSystematic_FiniteMCStats";
	cout << "UNFOLDING SYSTEMATICS..." << Form("UnfoldedData/hist_26282548_05302022_Unfolded_PbPbData_Eta15_ATLAS_Official_RAA_Binning_17Iters_10000Toys_NewpTBins%s_NewpTShpWghts.root",extratag1.c_str()) << endl; 
	sysUncrt_PbPb[counter] = new TFile(Form("UnfoldedData/hist_26282548_05302022_Unfolded_PbPbData_Eta15_ATLAS_Official_RAA_Binning_17Iters_10000Toys_NewpTBins%s_NewpTShpWghts.root",extratag1.c_str()),"READ");
	counter++;
      }
    }//PbPb Uncert
    if(debug)cout << __LINE__ << endl;

    //Adding JES 
    if(sysUncertType=="JES"){
      sysUncrt_PbPb[counter+1] = new TFile(Form("UnfoldedData/hist_26282548_05302022_Unfolded_PbPbData_Eta15_ATLAS_Official_RAA_Binning_17Iters_10000Toys_NewpTBins%sSys_%dptSHapeWeightsAppliedTo2DMatricies_NewpTShpWghts.root",sysUncertType.c_str(),17),"READ");
      sysUncrt_PbPb[counter+2] = new TFile(Form("UnfoldedData/hist_26282548_05302022_Unfolded_PbPbData_Eta15_ATLAS_Official_RAA_Binning_17Iters_10000Toys_NewpTBins%sSys_%dptSHapeWeightsAppliedTo2DMatricies_NewpTShpWghts.root",sysUncertType.c_str(),20),"READ");
    }else if(sysUncertType=="Unfolding"){
      sysUncrt_PbPb[counter] = new TFile(Form("UnfoldedData/hist_26282548_05302022_Unfolded_PbPbData_Eta15_ATLAS_Official_RAA_Binning_17Iters_10000Toys_NewpTBins_NewpTShpWghts.root"),"READ");
    }



     if(debug)cout << __LINE__ << endl;
  //File output
  TFile *sysroot = new TFile(Form("systematics/RAAsystematics_R%d_%s.root",jetRadius[jetR],sysUncertType.c_str()),"RECREATE");




  double sqrdSum[totCentBins][40] = {};

  if(debug)cout << "This is the uncertainty type: " << sysUncertType.c_str() << endl;

  for(int iCentBin =0; iCentBin < totCentBins; iCentBin++){
    double *array = NULL;
    array = bins.at(iCentBin).data();
     if(debug)cout << __LINE__ << endl;
   

    for(int iSysUncert = 0; iSysUncert < totUncert; iSysUncert++){ 
      cout << "This is JES uncert number : " <<  iSysUncert << endl;
      
      int iter= numIter_PbPb[jetR][iCentBin];
    
      int iterpp = numIter_pp[jetR][iCentBin];
      if(debug)cout << __LINE__ << endl;
      //systematic
      TH1D *sys_hist = (TH1D*) sysUncrt_PbPb[iSysUncert]->Get(Form("Unfolded_PbPbData_R%d_%dIter_%s",jetRadius[jetR],iter,centBins_2015Meas[iCentBin].c_str()));
      if(debug)cout << __LINE__ << endl;
      //Nominal measurement
      TH1D *nom = (TH1D*)nominalDataFilePbPb->Get(Form("Unfolded_PbPbData_R%d_%dIter_%s",jetRadius[jetR],iter,centBins_2015Meas[iCentBin].c_str()));
      if(debug)cout << __LINE__ << endl;
      TH1D *sys_hist_pp;
      if(iSysUncert !=19 || iSysUncert != 20){
	if(debug)cout << __LINE__ << endl;
	sys_hist_pp = (TH1D*) sysUncrt_pp[iSysUncert]->Get(Form("Unfolded_ppData_R%d_%dIter_%s",jetRadius[jetR],iterpp,centBins_2015Meas[iCentBin].c_str()));
	if(debug)cout << __LINE__ << endl;
      }else{
        sys_hist_pp  = (TH1D*) nominalDataFilepp->Get(Form("Unfolded_ppData_R%d_%dIter_%s",jetRadius[jetR],iterpp,centBins_2015Meas[iCentBin].c_str()));
      }
       if(debug)cout << __LINE__ << endl;
      TH1D *nompp = (TH1D*) nominalDataFilepp->Get(Form("Unfolded_ppData_R%d_%dIter_%s",jetRadius[jetR],iterpp,centBins_2015Meas[iCentBin].c_str()));

       if(debug)cout << __LINE__ << endl;

      

      

      sysroot->cd();
       if(debug)cout << __LINE__ << endl;
    	
      sys_hist->Scale(1/numEventsMinBias);
      sys_hist->Scale(1/tAA_2015map[iCentBin]);	 
      sys_hist->Scale(1.,"width");
      sys_hist->Scale(1/etaRange);
      if(debug)cout << __LINE__ << endl;

      nom->Scale(1/numEventsMinBias);
      nom->Scale(1/tAA_2015map[iCentBin]);
      nom->Scale(1.,"width");
      nom->Scale(1/etaRange);

      sys_hist_pp->Scale(1/ppDataLumiVals[jetR]);
      sys_hist_pp->Scale(1.,"width");
      sys_hist_pp->Scale(1/etaRange);
       if(debug)cout << __LINE__ << endl;
     
      nompp->Scale(1/ppDataLumiVals[jetR]);
      nompp->Scale(1.,"width");
      nompp->Scale(1/etaRange);
     
      
      sys_hist->Divide(sys_hist_pp);
      nom->Divide(nompp);

  	
      if(debug)cout << __LINE__ << endl;
     
      TH1D *plot = new TH1D(Form("RAA_R%d_%s_%d_Cent_%s",jetRadius[jetR],sysUncertType.c_str(),iSysUncert,centBins_2015Meas[iCentBin].c_str()),Form("R%d_%s_%d_Cent_%s",jetRadius[jetR],sysUncertType.c_str(),iSysUncert,centBins_2015Meas[iCentBin].c_str()),pTBinsTotR[jetR][iCentBin],array);
      
      for(int iBin = fiducialRegion[jetR][iCentBin][startAndend::start]; iBin < fiducialRegion[jetR][iCentBin][startAndend::end]; iBin++){
	
	float uncert = (abs(sys_hist->GetBinContent(iBin) - nom->GetBinContent(iBin))*100)/nom->GetBinContent(iBin);

	plot->SetBinContent(iBin,uncert);
	 if(debug)cout << __LINE__ << endl;
	sqrdSum[iCentBin][iBin] = sqrdSum[iCentBin][iBin] + pow(uncert,2);
	
	 if(iSysUncert ==0 && iCentBin ==0){
	   cout << "This is the bin number: " << iBin << endl;
	   cout << "This is bin center: (they should match) " << sys_hist->GetBinCenter(iBin) << "/" << nom->GetBinCenter(iBin) << endl;
	   cout << "This is rAA value of nominal: " << nom->GetBinContent(iBin) << endl;
	   cout << "This is the rAA of the sys: " << sys_hist->GetBinContent(iBin) << endl;
	   cout << "This is the uncertainty: " << uncert << endl;
	   //cout << "This would be the relative uncertainty: " << uncert/nom->GetBinConten << endl;
	 }

	
	  
       if(debug)cout << __LINE__ << endl;
      }//Bin Loop

       if(debug)cout << __LINE__ << endl;
	
      plot->Write(Form("RAA_R%d_%s_%d_Cent_%s",jetRadius[jetR],sysUncertType.c_str(),iSysUncert,centBins_2015Meas[iCentBin].c_str()),TObject::kOverwrite);
    
      
     sys_hist->Write(Form("RAA_R%d_%s_%d_Cent_%s_NotSysUncer",jetRadius[jetR],sysUncertType.c_str(),iSysUncert,centBins_2015Meas[iCentBin].c_str()),TObject::kOverwrite);
     nom->Write(Form("RAA_R%d_%s_%d_Cent_%s_NotSysUncer_Nominal",jetRadius[jetR],sysUncertType.c_str(),iSysUncert,centBins_2015Meas[iCentBin].c_str()),TObject::kOverwrite);


    }//Systematic Loop

  }//CentBin Loop

  if(debug)cout << __LINE__ << endl;
  

  for(int iCentBin =0; iCentBin < totCentBins; iCentBin++){
    double *array = NULL;
    array = bins.at(iCentBin).data();
     if(debug)cout << __LINE__ << endl;
    TH1D *quadsum = new TH1D(Form("RAA_QS_R%d_%s_Cent_%s",jetRadius[jetR],sysUncertType.c_str(),centBins_2015Meas[iCentBin].c_str()),Form("QS_R%d_%s_Cent_%s",jetRadius[jetR],sysUncertType.c_str(),centBins_2015Meas[iCentBin].c_str()),pTBinsTotR[jetR][iCentBin],array);
    
     if(debug)cout << __LINE__ << endl;
    
    for(int iBin = fiducialRegion[jetR][iCentBin][startAndend::start]; iBin < fiducialRegion[jetR][iCentBin][startAndend::end]; iBin++){
      cout << "This is bin number: " << iBin << endl;
      cout << "This is the qdr sum =:  " << sqrdSum[iCentBin][iBin] << endl;
      quadsum->SetBinContent(iBin,sqrt(sqrdSum[iCentBin][iBin]));
    }//bin loop

    quadsum->Write("",TObject::kOverwrite);
     if(debug)cout << __LINE__ << endl;
  }


  sysroot->Close();
   if(debug)cout << __LINE__ << endl;
  return 0;
}
