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
#include "bNec.h"
#include <TSystem.h>
//#include "ncollFunctions_5TeV.h"
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"



using namespace std;

int UnfoldingJetsCopy(bool Unfoldpp = true, bool UnfoldPbPb = false,bool Inefficciency_Corr = true, bool Fake_Corr = true, bool testUnfoldingPro = true, int jet_radius = R4,bool regCentBins = false,int numIter = 17, int etaRange = 21, int uncrt_sys = -1, string jer_or_jes = "",bool systematicUnfd = false, bool jetRate2015Bins = false,bool rAA2015Binning = false, bool dijet2018bins = true, bool officialLargeRpTBins = false,bool NopTShapeWeight=false,bool NoCentWeight = false){

  bool matchbins = false;


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
   
  if(dijet2018bins)number_CentBins=4;

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
 }else if(officialLargeRpTBins && etaRange ==21){
   location = Form("/usatlas/u/bereniceg299/data/LargeRJet_Study/NewSourceCode/largeR2018/%s/",colsnSys.c_str());
 }else if(officialLargeRpTBins && etaRange==15){
  location = Form("/usatlas/u/bereniceg299/data/LargeRJet_Study/NewSourceCode/NominalDir/R_%d/%s/",JetRadius[jet_radius],colsnSys.c_str());
 }

  if(jer_or_jes=="" && !NopTShapeWeight && !NoCentWeight){
    location = location + "nominal/";
  }else if(jer_or_jes!="" || NoCentWeight || NopTShapeWeight){
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
  }else if(officialLargeRpTBins && etacut ==2.1){
    location_of_dataFile=Form("/usatlas/u/bereniceg299/data/LargeRJet_Study/NewSourceCode/largeR2018/%s",colsnSys.c_str());
    location_of_mcFile = Form("/usatlas/u/bereniceg299/data/LargeRJet_Study/NewSourceCode/largeR2018/%s",colsnSys.c_str());
  }else if(officialLargeRpTBins && etacut ==1.5){
    location_of_dataFile=Form("/usatlas/u/bereniceg299/data/LargeRJet_Study/NewSourceCode/NominalDir/R_%d/%s/",JetRadius[jet_radius],colsnSys.c_str());
    location_of_mcFile = Form("/usatlas/u/bereniceg299/data/LargeRJet_Study/NewSourceCode/NominalDir/R_%d/%s/",JetRadius[jet_radius],colsnSys.c_str());
    
  }


   if(uncrt_sys!=-1){
    mc_tag = Form("_%s_%d",jer_or_jes.c_str(),uncrt_sys)+mc_tag;
  }

   mc_tag = Form("_eta%d",etaRange)+mc_tag;
   data_tag = Form("_eta%d",etaRange)+data_tag;
 



  if(jetRate2015Bins){
    mc_tag = "_2015JetRateBins" + mc_tag;
    data_tag = "_2015JetRateBins" + data_tag;
  }else if(rAA2015Binning){
    mc_tag = "_2015RAARateBins" +mc_tag;
    data_tag = "_2015RAARateBins" + data_tag;
  }else if(dijet2018bins){
    mc_tag = "_2018DiJetBins" + mc_tag;
    data_tag = "_2018DiJetBins" + data_tag;
  }else if(officialLargeRpTBins){
    mc_tag = "_2018LargerAnalysis" + mc_tag;
    data_tag = "_2018LargerAnalysis" + data_tag;
  }

  data_tag = "_NopTShpWghts"+ data_tag;
  if(jer_or_jes == "Unfolding" || NopTShapeWeight)mc_tag = "_NopTShpWghts" + mc_tag;
  

  mc_tag = Form("RawHistograms_R%d_%s",JetRadius[jet_radius],colsnSys.c_str()) + mc_tag;
  
  //mc_tag="RawHistograms_R4_PbPb_2018DiJetBins_eta21_v0_MC";
  if(!testUnfoldingPro){
    //You will unfold data
    data_tag = Form("RawHistograms_R%d_%s",JetRadius[jet_radius],colsnSys.c_str()) + data_tag;
  }else if(testUnfoldingPro){
    //You wil unfold MC Reco jets
    data_tag = mc_tag;
    if(matchbins)data_tag = "RawHistograms_R4_PbPb_NopTShpWghts_2018LargerAnalysis_eta15_MatchTruthAndRecopTBins_MC";

  }
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
  
    
  //When doing the unfolding test you will fill in the matricies in the bottom 
  //that belong to the half closure necessary items
  TH1D *pTDisRecoMatchedJets_HC[number_CentBins]; //Matched to Reco Jets
  TH1D *pTDisTruthMatchedJets_HC[number_CentBins]; //Matched to Truth Jets
  TH2D* RespMatrix_HC[number_CentBins];           //Response Matrix





  for(int iCentBin=0;iCentBin < number_CentBins; iCentBin++){


    string centBin_strg = "NAN";
    if(dijet2018bins){
      centBin_strg =  centBins_DJ[iCentBin];
    }else{
      centBin_strg = centBins_2015Meas[iCentBin];
    }


      
      
      cout << "THis is for iCentbin:" << iCentBin << endl;
      cout << "Histo: " << Form("R%d_Cent_%s",JetRadius[jet_radius],centBin_strg.c_str()) << endl;
      
      RespMatrix[iCentBin] = (TH2D*) forUnfolding->Get(Form("FullClsr_RespMatrix_R%d_Cent%s",JetRadius[jet_radius],centBin_strg.c_str()));

      if(!Inefficciency_Corr){
	pTDisTruthMatchedJets[iCentBin]= (TH1D*) forUnfolding->Get(Form("R%d_Cent_%s_TruthJetsMatched",JetRadius[jet_radius],centBin_strg.c_str()));
	
      }else{
	pTDisTruthMatchedJets[iCentBin]= (TH1D*) forUnfolding->Get(Form("R%d_Cent_%s_AllTruthJets",JetRadius[jet_radius],centBin_strg.c_str()));
      }

      if(!Fake_Corr){
	pTDisRecoMatchedJets[iCentBin] = (TH1D*) forUnfolding->Get(Form("R%d_Cent_%s",JetRadius[jet_radius],centBin_strg.c_str()));
      }else{
	pTDisRecoMatchedJets[iCentBin] = (TH1D*) forUnfolding->Get(Form("R%d_Cent_%s_AllRecoJets",JetRadius[jet_radius],centBin_strg.c_str()));
      }



      if(testUnfoldingPro){
	cout << "Grabbing this histos: " << Form("pTDisFirstHalfMatch_RecoJets_R%d%s",JetRadius[jet_radius],centBin_strg.c_str()) << endl;
	cout << "Grab this histo: " << Form("HalfClsr_RespMatrix_R%d%s",JetRadius[jet_radius],centBin_strg.c_str()) << endl;
	cout << "grab: " << Form("pTDisFirstHalfMatch_TruthJets_R%d%s",JetRadius[jet_radius],centBin_strg.c_str()) << endl;
	pTDisRecoMatchedJets_HC[iCentBin] = (TH1D*) forUnfolding->Get(Form("pTDisFirstHalfMatch_RecoJets_R%d%s",JetRadius[jet_radius],centBin_strg.c_str()));
	RespMatrix_HC[iCentBin] = (TH2D*) forUnfolding->Get(Form("HalfClsr_RespMatrix_R%d%s",JetRadius[jet_radius],centBin_strg.c_str()));
	pTDisTruthMatchedJets_HC[iCentBin]= (TH1D*) forUnfolding->Get(Form("pTDisFirstHalfMatch_TruthJets_R%d%s",JetRadius[jet_radius],centBin_strg.c_str()));
      }
      

     }//Cent Bins Loop
   //File that contains Not Unfolded  pT dis.
  
  TH1D * notUnfolded[number_CentBins]; //this will be for the FULL closure test when running test
  TH1D * notUnfolded_HC[number_CentBins]; //this will be for the HALF closure test when running test
  
  //Grabbing NOT unfolded pT Distributions
  for(int iCentBin=0;iCentBin<number_CentBins; iCentBin++){
    string centBin_strg = "NAN";
    if(dijet2018bins){
      centBin_strg =  centBins_DJ[iCentBin];
    }else{
      centBin_strg = centBins_2015Meas[iCentBin];
    }


       
    notUnfolded[iCentBin] = (TH1D*) pTDis_NotUnfolded->Get(Form("R%d_Cent_%s", JetRadius[jet_radius],centBin_strg.c_str()));
    if(testUnfoldingPro){
      notUnfolded_HC[iCentBin] = (TH1D*) pTDis_NotUnfolded->Get(Form("pTDisSecndHalfMatch_RecoJets_R%d%s", JetRadius[jet_radius],centBin_strg.c_str()));
    }

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
  }else if(officialLargeRpTBins){
    jetRate2015BinsTag = "_2018LargerAnalysis";
  }

  string NocentWeight_tag = "";
  if(NoCentWeight){
    NocentWeight_tag= "_NoCentWghts";
  }

  string NopTShpWghts_tag ="";
  if(NopTShapeWeight){
    NopTShpWghts_tag = "_NopTShpWght";
  }
    


  string NOMINAL_TAG =""; 
  if(jer_or_jes=="" && !testUnfoldingPro && !NoCentWeight && !NopTShapeWeight)NOMINAL_TAG = "Nominal";
  if(testUnfoldingPro)NOMINAL_TAG = "_UnfoldingProcedureTEST";
  cout << "*******This is the file where we will save everything: " << Form("%shist_26282548_05302022_Unfolded_%s_ATLAS_Official_RAA_Binning_17Iters_10000Toys%s%s%s%s_etaRange%d%s.root***",location.c_str(),pp_or_PbPb.c_str(),uncertSys.c_str(),ptShapeAppliedMatriciesTag.c_str(),unfoldSystematic.c_str(),jetRate2015BinsTag.c_str(),etaRange,NOMINAL_TAG.c_str());
  
  TFile *pTDisUnfold_File = new TFile(Form("%shist_26282548_05302022_Unfolded_%s_ATLAS_Official_RAA_Binning_17Iters_10000Toys%s%s%s%s_etaRange%d%s%s%s_R%d.root",location.c_str(),pp_or_PbPb.c_str(),uncertSys.c_str(),ptShapeAppliedMatriciesTag.c_str(),unfoldSystematic.c_str(),jetRate2015BinsTag.c_str(),etaRange,NOMINAL_TAG.c_str(),NocentWeight_tag.c_str(),NopTShpWghts_tag.c_str(),JetRadius[jet_radius]),"RECREATE"); //ATLAS RAA Binning
  
  string centBins = "NoCentBins";
  
  
  for(int iCentBin =0; iCentBin <number_CentBins; iCentBin++){
    
    string centBin_strg = "NAN";
    double scaling_factor = 1.0;


    if(dijet2018bins){
      centBin_strg =  centBins_DJ[iCentBin];
      if(UnfoldPbPb && !testUnfoldingPro)scaling_factor = tAASubstructure_2018map[iCentBin]*numEvents_centBin[iCentBin];
      if(Unfoldpp && !testUnfoldingPro)scaling_factor = ppDataLumiVals[jet_radius];
    }else{
      centBin_strg = centBins_2015Meas[iCentBin];
    }




      if(Unfoldpp && iCentBin > 0 && regCentBins)continue;
      cout << __LINE__ << endl;
      if(UnfoldPbPb || (Unfoldpp && !regCentBins)){
	
	centBins = centBining[index_CentBinnning][iCentBin];
      }
      cout << __LINE__ << endl;
      RooUnfoldResponse *unf_R = nullptr;
      RooUnfoldResponse *unf_R_HC = nullptr;
      cout << __LINE__ << endl;
      cout << "iCentBin: " << iCentBin << endl;
      pTDisRecoMatchedJets[iCentBin]->SetMarkerColor(kRed);
      cout << __LINE__ << endl;
      pTDisTruthMatchedJets[iCentBin]->SetMarkerColor(kBlue);
      cout << __LINE__ << endl;

      unf_R = new RooUnfoldResponse(pTDisRecoMatchedJets[iCentBin],pTDisTruthMatchedJets[iCentBin],RespMatrix[iCentBin]);
      cout << __LINE__ << endl;


      if(testUnfoldingPro){
      	cout << __LINE__ << endl;
      	pTDisRecoMatchedJets_HC[iCentBin]->SetMarkerColor(kRed);
      	cout << __LINE__ << endl;
      	pTDisTruthMatchedJets_HC[iCentBin]->SetMarkerColor(kRed);
       
      	RespMatrix_HC[iCentBin]->Draw("colz");
      	cout << __LINE__ << endl;
      	unf_R_HC = new RooUnfoldResponse(pTDisRecoMatchedJets_HC[iCentBin],pTDisTruthMatchedJets_HC[iCentBin],RespMatrix_HC[iCentBin]);
      	cout << __LINE__ << endl;
      }
      cout << __LINE__ << endl;

      notUnfolded[iCentBin]->Scale(1/scaling_factor);

      for(int iIter = 0; iIter < numIter; iIter++){
      	cout << __LINE__ << endl;
      	RooUnfoldBayes *rooUnfold_Data = nullptr;
      	cout << __LINE__ << endl;
      	RooUnfoldBayes *rooUnfold_Data_HC = nullptr;

      	rooUnfold_Data = new RooUnfoldBayes(unf_R,notUnfolded[iCentBin], iIter);
      	cout << __LINE__ << endl;
	
      	if(testUnfoldingPro){
      	    rooUnfold_Data_HC = new RooUnfoldBayes(unf_R_HC,notUnfolded_HC[iCentBin], iIter);
      	    rooUnfold_Data_HC->SetNToys(nToys);
      	 }

      	if(!systematicUnfd){
      	  rooUnfold_Data->SetNToys(nToys);
	  
      	}
      	cout << __LINE__ << endl;
      	if(systematicUnfd){
      	  rooUnfold_Data->SetNToys(50);
      	  rooUnfold_Data->IncludeSystematics(2);
	  
      	}
      	cout << __LINE__ << endl;
      	 TH1D *unfolded_Data = nullptr;
      	 TH1D *unfolded_Data_HC = nullptr;
      	 cout << __LINE__ << endl;
      	if(!systematicUnfd){
      	  cout << __LINE__ << endl;
      	  unfolded_Data = (TH1D*) rooUnfold_Data->Hreco(RooUnfold::kCovToy);
	  
      	  if(testUnfoldingPro){
      	    cout << __LINE__ << endl;
      	    unfolded_Data_HC = (TH1D*) rooUnfold_Data_HC->Hreco(RooUnfold::kCovToy);
      	    cout << __LINE__ << endl;
      	    rooUnfold_Data_HC->Hreco(RooUnfold::kCovToy);
      	    cout << __LINE__ << endl;
      	  }
      	  cout << __LINE__ << endl;
      	}
      	if(systematicUnfd)unfolded_Data = (TH1D*) rooUnfold_Data->Hreco(RooUnfold::ErrorTreatment::kCovToy);
      	if(!testUnfoldingPro){
      	  cout << __LINE__ << endl;
      	  unfolded_Data->Write(Form("Unfolded_%s_R%d_%dIter_%s",pp_or_PbPb.c_str(),JetRadius[jet_radius],iIter,centBin_strg.c_str()),TObject::kOverwrite);
      	  unfolded_Data = nullptr;
      	}else if(testUnfoldingPro){
      	  cout << __LINE__ << endl;
      	  unfolded_Data->Write(Form("pTDisUnfold_Full_R%d_%dIter_Cent%s",JetRadius[jet_radius],iIter,centBin_strg.c_str()),TObject::kOverwrite);
          unfolded_Data = nullptr;
      	  cout << __LINE__ << endl;
      	  unfolded_Data_HC->Write(Form("pTDisUnfold_Half_R%d_%dIter_Cent%s",JetRadius[jet_radius],iIter,centBin_strg.c_str()),TObject::kOverwrite);
      	  unfolded_Data_HC = nullptr;
      	  cout << __LINE__ << endl;
	  
      }
      }//CentBin Loop
      
  }//Jet
    

  pTDis_NotUnfolded->Close();
  pTDisUnfold_File->Close();
  

  
  return 0;

}
