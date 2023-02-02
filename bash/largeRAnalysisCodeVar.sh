#int LargeRAnalysis(string collisionType = "PbPb", int jet_Rad = R4,int data_or_mc = mcOrdata::mc,  string extraTag="MC", bool debug =false, int sys_uncrt =-1, string jer_or_jes_Sys = "JES", bool jetRate2015Bins=false, bool rAA2015Binning = true){



#Setting Variables for Analysis code
collisionType="PbPb"
jet_Rad=1 #1 is R=0.4 and 0 is R=1.0 
data_or_mc=1 #1 is MC and 0 is Data 
extraTag="MC"
debug=false #You want to debug? If yes, select true, but if not then select false
jer_or_jes_Sys="JER" #Systematic type
jetRate2015Bins=false #2015 Jet rate pT Binning
rAA2015Binning=true #2015 RAA pT Binning






root -l -b 'LargeRAnalysis.C("'$collisionType'",'$jet_Rad','$data_or_mc',"'$extraTag'",'$debug',0,"'$jer_or_jes_Sys'",'$jetRate2015Bins','$rAA2015Binning')' &
root -l -b 'LargeRAnalysis.C("'$collisionType'",'$jet_Rad','$data_or_mc',"'$extraTag'",'$debug',1,"'$jer_or_jes_Sys'",'$jetRate2015Bins','$rAA2015Binning')' &
root -l -b 'LargeRAnalysis.C("'$collisionType'",'$jet_Rad','$data_or_mc',"'$extraTag'",'$debug',2,"'$jer_or_jes_Sys'",'$jetRate2015Bins','$rAA2015Binning')' &
root -l -b 'LargeRAnalysis.C("'$collisionType'",'$jet_Rad','$data_or_mc',"'$extraTag'",'$debug',3,"'$jer_or_jes_Sys'",'$jetRate2015Bins','$rAA2015Binning')' &
root -l -b 'LargeRAnalysis.C("'$collisionType'",'$jet_Rad','$data_or_mc',"'$extraTag'",'$debug',4,"'$jer_or_jes_Sys'",'$jetRate2015Bins','$rAA2015Binning')' &
root -l -b 'LargeRAnalysis.C("'$collisionType'",'$jet_Rad','$data_or_mc',"'$extraTag'",'$debug',5,"'$jer_or_jes_Sys'",'$jetRate2015Bins','$rAA2015Binning')' &
root -l -b 'LargeRAnalysis.C("'$collisionType'",'$jet_Rad','$data_or_mc',"'$extraTag'",'$debug',6,"'$jer_or_jes_Sys'",'$jetRate2015Bins','$rAA2015Binning')' &
root -l -b 'LargeRAnalysis.C("'$collisionType'",'$jet_Rad','$data_or_mc',"'$extraTag'",'$debug',7,"'$jer_or_jes_Sys'",'$jetRate2015Bins','$rAA2015Binning')' &
root -l -b 'LargeRAnalysis.C("'$collisionType'",'$jet_Rad','$data_or_mc',"'$extraTag'",'$debug',8,"'$jer_or_jes_Sys'",'$jetRate2015Bins','$rAA2015Binning')' &
# root -l -b 'LargeRAnalysis.C("'$collisionType'",'$jet_Rad','$data_or_mc',"'$extraTag'",'$debug',9,"'$jer_or_jes_Sys'",'$jetRate2015Bins','$rAA2015Binning')' &
# root -l -b 'LargeRAnalysis.C("'$collisionType'",'$jet_Rad','$data_or_mc',"'$extraTag'",'$debug',10,"'$jer_or_jes_Sys'",'$jetRate2015Bins','$rAA2015Binning')' &
# root -l -b 'LargeRAnalysis.C("'$collisionType'",'$jet_Rad','$data_or_mc',"'$extraTag'",'$debug',11,"'$jer_or_jes_Sys'",'$jetRate2015Bins','$rAA2015Binning')' &
# root -l -b 'LargeRAnalysis.C("'$collisionType'",'$jet_Rad','$data_or_mc',"'$extraTag'",'$debug',12,"'$jer_or_jes_Sys'",'$jetRate2015Bins','$rAA2015Binning')' &
# root -l -b 'LargeRAnalysis.C("'$collisionType'",'$jet_Rad','$data_or_mc',"'$extraTag'",'$debug',13,"'$jer_or_jes_Sys'",'$jetRate2015Bins','$rAA2015Binning')' &
# root -l -b 'LargeRAnalysis.C("'$collisionType'",'$jet_Rad','$data_or_mc',"'$extraTag'",'$debug',14,"'$jer_or_jes_Sys'",'$jetRate2015Bins','$rAA2015Binning')' &
# root -l -b 'LargeRAnalysis.C("'$collisionType'",'$jet_Rad','$data_or_mc',"'$extraTag'",'$debug',15,"'$jer_or_jes_Sys'",'$jetRate2015Bins','$rAA2015Binning')' &
# root -l -b 'LargeRAnalysis.C("'$collisionType'",'$jet_Rad','$data_or_mc',"'$extraTag'",'$debug',16,"'$jer_or_jes_Sys'",'$jetRate2015Bins','$rAA2015Binning')' &
# root -l -b 'LargeRAnalysis.C("'$collisionType'",'$jet_Rad','$data_or_mc',"'$extraTag'",'$debug',17,"'$jer_or_jes_Sys'",'$jetRate2015Bins','$rAA2015Binning')' &
# root -l -b 'LargeRAnalysis.C("'$collisionType'",'$jet_Rad','$data_or_mc',"'$extraTag'",'$debug',18,"'$jer_or_jes_Sys'",'$jetRate2015Bins','$rAA2015Binning')' &
# root -l -b 'LargeRAnalysis.C("'$collisionType'",'$jet_Rad','$data_or_mc',"'$extraTag'",'$debug',19,"'$jer_or_jes_Sys'",'$jetRate2015Bins','$rAA2015Binning')' &
#root -l -b 'LargeRAnalysis.C("'$collisionType'",'$jet_Rad','$data_or_mc',"'$extraTag'",'$debug',20,"'$jer_or_jes_Sys'",'$jetRate2015Bins','$rAA2015Binning')' &




