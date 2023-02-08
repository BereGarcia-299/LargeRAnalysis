#include <assert.h>
#include <cmath>
#include <fstream>
#include <stdio.h>
#include <array>
#include <tgmath.h>
#include <iostream>
#include <TChain.h>
#include <TLegend.h>
#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TRandom3.h"
#include <TAttMarker.h>
#include <TCanvas.h>
#include <TPaveText.h>
#include "TLorentzVector.h"
#include "bFunctions.h"
#include "NewTreeVariables.h"
#include <TSystem.h>
//#include "ncollFunctions_5TeV.h"
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"



using namespace std;

int UnfoldingJetsCopy(bool Unfoldpp = false, bool UnfoldPbPb = true, int jet_radius = R4,bool regCentBins = false,int numIter = 17, int etaRange = 21, int uncrt_sys = -1, string jer_or_jes = "",bool systematicUnfd = false, bool jetRate2015Bins = false,bool rAA2015Binning = false, bool dijet2018bins = true){

  //We set pp_Or_PbPb to true if unfolding Pb+Pb Data
  string pp_or_PbPb = "ppData" ;
  string colsnSys = "pp";
  if(UnfoldPbPb){
    pp_or_PbPb="PbPbData";
    colsnSys = "PbPb";
  }
  int index_CentBinnning = 0;
  int number_CentBins = 0;
  bool usepTShapeWeightsMatricies = false;
  usepTShapeWeightsMatricies = false;

  if(!regCentBins){
    //cout << "We will NOT use tim's cent bins " << endl;
    number_CentBins = num2015MeasBins;
    index_CentBinnning = 0;
  }else{
    number_CentBins = totCentBins_DJ;
    index_CentBinnning = 1;
  }
   
  cout << "Total number of cent bins: " << number_CentBins << endl;
  

  //location of output files
  string location = "/usatlas/u/bereniceg299/data/LargeRJet_Study/NewSourceCode/mc_test_output/";
  if(uncrt_sys!=-1)location="/usatlas/u/bereniceg299/data/LargeRJet_Study/NewSourceCode/SysUncert/UnfoldedData/";
  
 if(jetRate2015Bins){
    location="/usatlas/u/bereniceg299/data/LargeRJet_Study/NewSourceCode/2015CentBinsFiles/";
  }else if(rAA2015Binning){
    location=Form("/usatlas/u/bereniceg299/data/LargeRJet_Study/NewSourceCode/2015CentBinsFiles/rAABins_%s/",colsnSys.c_str());
  }else if(dijet2018bins){
    location = Form("/usatlas/u/bereniceg299/data/LargeRJet_Study/NewSourceCode/DijetBinning/%s/",colsnSys.c_str());    
  }

  if(jer_or_jes==""){
    location = location + "nominal/";
  }else if(jer_or_jes!=""){
    location = location + "Unfolded/";
  }
    //These will contain necessary
  //histograms for the unfolding procedure 
  TFile *forUnfolding = NULL;
  
  //pT Distributions that Need Unfolding
  TFile *pTDis_NotUnfolded = NULL;
  
  //Beginning Tag 
  string mc_tag = "_MC";
  string data_tag = "_Data";
  string location_of_dataFile = "";
  string location_of_mcFile = "";

  if(jetRate2015Bins){
    location_of_dataFile = "/usatlas/u/bereniceg299/data/LargeRJet_Study/NewSourceCode/2015CentBinsFiles";
    location_of_mcFile = "/usatlas/u/bereniceg299/data/LargeRJet_Study/NewSourceCode/2015CentBinsFiles/";

    if(jer_or_jes!="")location_of_mcFile = location_of_mcFile + Form("%sVar",colsnSys.c_str());
  }else if(rAA2015Binning){
    location_of_dataFile = "/usatlas/u/bereniceg299/data/LargeRJet_Study/NewSourceCode/2015CentBinsFiles";
    location_of_mcFile = Form("/usatlas/u/bereniceg299/data/LargeRJet_Study/NewSourceCode/2015CentBinsFiles/rAABins_%s",colsnSys.c_str());
  }else if(dijet2018bins){
    location_of_dataFile = Form("/usatlas/u/bereniceg299/data/LargeRJet_Study/NewSourceCode/DijetBinning/%s",colsnSys.c_str());
    location_of_mcFile = Form("/usatlas/u/bereniceg299/data/LargeRJet_Study/NewSourceCode/DijetBinning/%s",colsnSys.c_str());
  }


   if(uncrt_sys!=-1){
    mc_tag = Form("_%s_%d",jer_or_jes.c_str(),uncrt_sys)+mc_tag;
  }

   mc_tag = Form("_eta%d",etaRange)+mc_tag;
   data_tag = Form("_eta%d",etaRange)+data_tag;
 


  if(jetRate2015Bins==true){
    mc_tag = "_2015JetRateBins" + mc_tag;
    data_tag = "_2015JetRateBins" + data_tag;
  }else if(rAA2015Binning==true){
    mc_tag = "_2015RAARateBins" +mc_tag;
    data_tag = "_2015RAARateBins" + data_tag;
  }else if(dijet2018bins){
    mc_tag = "_2018DiJetBins" + mc_tag;
    data_tag = "_2018DiJetBins" + data_tag;
  }

  data_tag = "_NopTShpWghts"+ data_tag;
  if(jer_or_jes == "Unfolding")mc_tag = "_NopTShpWghts" + mc_tag;
  

  mc_tag = Form("RawHistograms_R%d_%s",JetRadius[jet_radius],colsnSys.c_str()) + mc_tag;
  data_tag = Form("RawHistograms_R%d_%s",JetRadius[jet_radius],colsnSys.c_str()) + data_tag;
  
  std::ofstream file("name_files_used.txt"); // open text file  

  
  if (file.is_open()) {
    file << "These are the files you used for the Unfolding Procedure..." << std::endl;
    file << "Raw Data File: " << Form("%s/%s.root",location_of_dataFile.c_str(),data_tag.c_str()) << std::endl;
    file << "Raw MC File: " << Form("%s/%s.root",location_of_mcFile.c_str(),mc_tag.c_str()) << endl;
    file.close(); // close the file
  } else {
    std::cout << "Unable to open file" << std::endl;
  }

  
  //File with histos that have raw spectras
  pTDis_NotUnfolded = new TFile(Form("%s/%s.root",location_of_dataFile.c_str(),data_tag.c_str()),"READ");
    
  //File w/ necessary histos to Unfold Raw Data
  forUnfolding = new TFile(Form("%s/%s.root",location_of_mcFile.c_str(),mc_tag.c_str()),"READ");

    
    

  
  if(Unfoldpp && false){
   
    if(uncrt_sys==-1){

      if(!usepTShapeWeightsMatricies){
        forUnfolding = new TFile("RawHistograms_pp_NopTShpWeights_010523_MC.root","READ"); //Jet Rate (2015) ATLAS Official Bins
      }else if(usepTShapeWeightsMatricies){
        forUnfolding = new TFile("RawHistograms_pp_ptShpWeights_MC.root","READ");
      }

       
    }
     pTDis_NotUnfolded = new TFile("RawHistograms_pp_Data.root","READ");

     


  }else if(UnfoldPbPb && false){
    if(uncrt_sys==-1){
     if(!usepTShapeWeightsMatricies){
	forUnfolding = new TFile("RawHistograms_PbPb_NopTShpWeights_010523_MC.root","READ");
      }else if(usepTShapeWeightsMatricies){
	forUnfolding = new TFile("RawHistograms_PbPb_centWeights_ptShpWeights_MC.root","READ");
      }
      //forUnfolding = new TFile("RawHistograms_PbPb_AddCentWeights_ptShapeWeights_MC.root","READ");
    }
    pTDis_NotUnfolded = new TFile("RawHistograms_PbPb_Data.root","READ");
  }

  //Reco Truth-Matched Jets 
  TH1D *pTDisRecoMatchedJets[number_CentBins]; //Matched to Reco Jets
  TH1D *pTDisTruthMatchedJets[number_CentBins]; //Matched to Truth Jets
  
  TH1D *pTDisRecoNoMatchJets[number_CentBins]; //No matching req.
  TH2D* RespMatrix[number_CentBins];           //Response Matrix
  TH1D *pTDisAllTruthJets[number_CentBins];    //All Truth Jets
  
    
  for(int iCentBin=0;iCentBin < number_CentBins; iCentBin++){
      
      pTDisRecoMatchedJets[iCentBin] = (TH1D*) forUnfolding->Get(Form("R%d_Cent_%s",JetRadius[jet_radius],centBins_2015Meas[iCentBin].c_str()));
      RespMatrix[iCentBin] = (TH2D*) forUnfolding->Get(Form("FullClsr_RespMatrix_R%d_Cent%s",JetRadius[jet_radius],centBins_2015Meas[iCentBin].c_str()));
      pTDisTruthMatchedJets[iCentBin]= (TH1D*) forUnfolding->Get(Form("R%d_Cent_%s_TruthJetsMatched",JetRadius[jet_radius],centBins_2015Meas[iCentBin].c_str()));

     }//Cent Bins Loop
      
  
  //File that contains Not Unfolded  pT dis.
  
  TH1D * notUnfolded[number_CentBins];
  
  //Grabbing NOT unfolded pT Distributions
  for(int iCentBin=0;iCentBin<number_CentBins; iCentBin++){
       
    notUnfolded[iCentBin] = (TH1D*) pTDis_NotUnfolded->Get(Form("R%d_Cent_%s", JetRadius[jet_radius],centBins_2015Meas[iCentBin].c_str()));
  }//Cent Bin Loop
  
  //This is where we will store
  //Unfolded pT Distributions
  //TFile *pTDisR4UnfoldRootFile = new TFile(Form("%shist_26282548_03012022_UnfoldedppDataR4_Eta15_NewCode_CMS_12Bins_New.root",location.c_str()),"RECREATE"); //CMS Binning

  //Systematics Tag
  string uncertSys = "";
  if(uncrt_sys!=-1)uncertSys=Form("_%sSys_%d",jer_or_jes.c_str(),uncrt_sys);
  string ptShapeAppliedMatriciesTag = "";
  if(usepTShapeWeightsMatricies)ptShapeAppliedMatriciesTag="ptSHapeWeightsAppliedTo2DMatricies";
  string unfoldSystematic = "";
  if(systematicUnfd){
    unfoldSystematic = "_UnfoldingSystematic_FiniteMCStats";
  }else if(systematicUnfd==false && jer_or_jes=="Unfolding"){
    unfoldSystematic = "_Unfold_pTShpWghtSys";
  }
  string jetRate2015BinsTag = "";
  if(jetRate2015Bins){
    jetRate2015BinsTag="_JetRate2015Bins";
  }else if(rAA2015Binning){
    jetRate2015BinsTag = "_2015RAARateBins";
  
  }else if(dijet2018bins){
    jetRate2015BinsTag = "_2018DiJetBins";
  }
  string NOMINAL_TAG =""; 
  if(jer_or_jes=="")NOMINAL_TAG = "Nominal";

  cout << "*******This is the file where we will save everything: " << Form("%shist_26282548_05302022_Unfolded_%s_ATLAS_Official_RAA_Binning_17Iters_10000Toys%s%s%s%s_etaRange%d%s.root***",location.c_str(),pp_or_PbPb.c_str(),uncertSys.c_str(),ptShapeAppliedMatriciesTag.c_str(),unfoldSystematic.c_str(),jetRate2015BinsTag.c_str(),etaRange,NOMINAL_TAG.c_str());
  
  TFile *pTDisUnfold_File = new TFile(Form("%shist_26282548_05302022_Unfolded_%s_ATLAS_Official_RAA_Binning_17Iters_10000Toys%s%s%s%s_etaRange%d%s.root",location.c_str(),pp_or_PbPb.c_str(),uncertSys.c_str(),ptShapeAppliedMatriciesTag.c_str(),unfoldSystematic.c_str(),jetRate2015BinsTag.c_str(),etaRange,NOMINAL_TAG.c_str()),"RECREATE"); //ATLAS RAA Binning
  
  string centBins = "NoCentBins";
  
  
  for(int iCentBin =0; iCentBin <number_CentBins; iCentBin++){
    
      if(Unfoldpp && iCentBin > 0 && regCentBins)continue;
      cout << __LINE__ << endl;
      if(UnfoldPbPb || (Unfoldpp && !regCentBins)){
	
	centBins = centBining[index_CentBinnning][iCentBin];
      }
      cout << __LINE__ << endl;
      RooUnfoldResponse *unf_R = nullptr;
       cout << __LINE__ << endl;
      pTDisRecoMatchedJets[iCentBin]->SetMarkerColor(kRed);
      cout << __LINE__ << endl;
      pTDisTruthMatchedJets[iCentBin]->SetMarkerColor(kBlue);
      
      unf_R = new RooUnfoldResponse(pTDisRecoMatchedJets[iCentBin],pTDisTruthMatchedJets[iCentBin],RespMatrix[iCentBin]);
      
      for(int iIter = 0; iIter < numIter; iIter++){
	RooUnfoldBayes *rooUnfold_Data = nullptr;
	rooUnfold_Data = new RooUnfoldBayes(unf_R,notUnfolded[iCentBin], iIter);
	if(!systematicUnfd)rooUnfold_Data->SetNToys(nToys);
	if(systematicUnfd){
	  rooUnfold_Data->SetNToys(50);
	  rooUnfold_Data->IncludeSystematics(2);
	  
	}
	 TH1D *unfolded_Data = nullptr;
	cout << __LINE__ << endl;
	
	if(!systematicUnfd)unfolded_Data = (TH1D*) rooUnfold_Data->Hreco(RooUnfold::kCovToy);
	if(systematicUnfd)unfolded_Data = (TH1D*) rooUnfold_Data->Hreco(RooUnfold::ErrorTreatment::kCovToy);
	unfolded_Data->Write(Form("Unfolded_%s_R%d_%dIter_%s",pp_or_PbPb.c_str(),JetRadius[jet_radius],iIter,centBins.c_str()),TObject::kOverwrite);
	unfolded_Data = nullptr;
      }//CentBin Loop
      
    }//Jet
    

  pTDis_NotUnfolded->Close();
  pTDisUnfold_File->Close();
  

  
  return 0;

}
