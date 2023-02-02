#This bash script is meant to test my analysis and compare the RAA results to the 2015 measurement. 

#############################################################
##-------Arguements that belong to my analysis code---------#
#############################################################


#int LargeRAnalysis(string collisionType = "Data", int jet_Rad = R4,int data_or_mc = mcOrdata::data,  string extraTag="test_Data", bool debug =false, int sys_uncrt =-1, string jer_or_jes_Sys = "",bool addpTShapeWeight = false,bool jetRate2015Bins=false, bool rAA2015Binning = true){

path_of_analysis_code="/usatlas/u/bereniceg299/data/LargeRJet_Study/NewSourceCode"


#------1.Run analysis code over pp and Pb+Pb Data and MC
#Setting Variables for Analysis code
jet_Rad=1                #1 is R=0.4 and 0 is R=1.0
data_or_mc=0             #1 is MC and 0 is Data
extraTag="Data"          #Tag this to root file created
debug=false              #You want to debug? If yes, select true, but if not then select false
sys_uncrt=-1             #Systematic Vartiation you want to look at
jer_or_jes_Sys=""        #Systematic type
addpTShapeWeight=false   #Apply pT shape weights
jetRate2015Bins=false    #2015 Jet rate pT Binning
rAA2015Binning=true      #2015 RAA pT Binning

echo "We will start by running your analysis code over the 2017 pp and 2018 Pb+Pb data"

#Data
root -l -b -q ''$path_of_analysis_code'/LargeRAnalysis.C("PbPb",'$jet_Rad','$data_or_mc',"'$extraTag'",'$debug','$sys_uncrt',"'$jer_or_jes_Sys'",'$addpTShapeWeight','$jetRate2015Bins','$rAA2015Binning')' &
root -l -b -q ''$path_of_analysis_code'/LargeRAnalysis.C("pp",'$jet_Rad','$data_or_mc',"'$extraTag'",'$debug','$sys_uncrt',"'$jer_or_jes_Sys'",'$addpTShapeWeight','$jetRate2015Bins','$rAA2015Binning')' &

#MC
#Modifying the varibles below so analysis code will run over MC
extraTag="MC" 
data_or_mc=1  
addpTShapeWeight=true   #Apply pT shape weights

echo "While those run in parallel we will also want to run the MC samples through your analysis code!" 

root -l -b -q ''$path_of_analysis_code'/LargeRAnalysis.C("PbPb",'$jet_Rad','$data_or_mc',"'$extraTag'",'$debug','$sys_uncrt',"'$jer_or_jes_Sys'",'$addpTShapeWeight','$jetRate2015Bins','$rAA2015Binning')' &
root -l -b -q ''$path_of_analysis_code'/LargeRAnalysis.C("pp",'$jet_Rad','$data_or_mc',"'$extraTag'",'$debug','$sys_uncrt',"'$jer_or_jes_Sys'",'$addpTShapeWeight','$jetRate2015Bins','$rAA2015Binning')' &

# #Running over MC to NOT apply pT shape weight to pT distributions and also the 2D response matricies
# addpTShapeWeight=false   #Apply pT shape weights

# root -l -b 'LargeRAnalysis.C("PbPb",'$jet_Rad','$data_or_mc',"'$extraTag'",'$debug','$sys_uncrt',"'$jer_or_jes_Sys'",'$addpTShapeWeight','$jetRate2015Bins','$rAA2015Binning')'> logMC_PbPb.log &
# root -l -b 'LargeRAnalysis.C("pp",'$jet_Rad','$data_or_mc',"'$extraTag'",'$debug','$sys_uncrt',"'$jer_or_jes_Sys'",'$addpTShapeWeight','$jetRate2015Bins','$rAA2015Binning')' > logMC_PbPb.log &


# #We will apply this code so that it will wait until the above commands are done running 
# pid=$!
# wait $pid

# #Running over MC but producing the JER Variations
# collisionType="PbPb"
# addpTShapeWeight=true   #Apply pT shape weights
# jer_or_jes_Sys="JER"
# for i in {0..20}; do
#     root -l -b 'LargeRAnalysis.C("'$collisionType'",'$jet_Rad','$data_or_mc',"'$extraTag'",'$debug','$i',"'$jer_or_jes_Sys'",'$addpTShapeWeight','$jetRate2015Bins','$rAA2015Binning')' &
# done
# #We will apply this code so that it will wait until the above commands are done running
# pid=$!
# wait $pid

# #Running over MC but producing the JER Variations
# collisionType="pp"
# for i in {0..8}; do
#     root -l -b 'LargeRAnalysis.C("'$collisionType'",'$jet_Rad','$data_or_mc',"'$extraTag'",'$debug','$i',"'$jer_or_jes_Sys'",'$addpTShapeWeight','$jetRate2015Bins','$rAA2015Binning')' &
# done


# #We will apply this code so that it will wait until the above commands are done running
# pid=$!
# wait $pid
# #Running over MC to produce JES systematics
# jer_or_jes_Sys="JES"
# for i in {0..19}; do
#     root -l -b 'LargeRAnalysis.C("'$collisionType'",'$jet_Rad','$data_or_mc',"'$extraTag'",'$debug','$i',"'$jer_or_jes_Sys'",'$addpTShapeWeight','$jetRate2015Bins','$rAA2015Binning')' &
# done

# #We will apply this code so that it will wait until the above commands are done running
# pid=$!
# wait $pid

# #Running over MC but producing the JES Variations
# collisionType="Pb+Pb"

# for i in {0..20}; do
#     root -l -b 'LargeRAnalysis.C("'$collisionType'",'$jet_Rad','$data_or_mc',"'$extraTag'",'$debug','$i',"'$jer_or_jes_Sys'",'$addpTShapeWeight','$jetRate2015Bins','$rAA2015Binning')' &
# done


# #We will apply this code so that it will wait until the above commands are done running
# pid=$!
# wait $pid

# #UNFOLDING STEP
# bash unfolding_step.sh
