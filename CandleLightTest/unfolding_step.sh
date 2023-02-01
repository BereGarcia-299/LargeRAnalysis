#This bash script is dedicated to execute the unfolding step of my analysis

#Here are the arguemenets and the name of the function that takes care of this step

#int UnfoldingJetsCopy(bool Unfoldpp = false, bool UnfoldPbPb = true, int jet_radius = R4,bool regCentBins = false,int numIter = 17, int etaRange = 28, int uncrt_sys = 2, string jer_or_jes = "JES",bool systematicUnfd = false, bool jetRate2015Bins = false,bool rAA2015Binning = true){


#Setting Variables to Unfold the 
unfoldpp=false
unfoldPbPb=true
jet_radius=1 #R=0.4
regCentBins=false
numIter=17
etaRange=28
uncrt_sys=0 #variation
jer_or_jes="JES"
systematicUnfd=false #apply unfolding systematic (finite MC stats)
jetRate2015Bins=false #2015 jet rate bins
rAA2015Binning=true #2015 RAA bins

#############################################
##########      Pb+Pb          ##############
############################################

#This will Unfold the JES root files that I have made with my analysis code
for i in {0..20}; do
  root -l -b 'UnfoldingJetsCopy.C('$unfoldpp','$unfoldPbPb','$jet_radius','$regCentBins','$numIter', '$etaRange', '$i',"'$jer_or_jes'",'$systematicUnfd', '$jetRate2015Bins','$rAA2015Binning')' &

done

#We will apply this code so that it will wait until the above commands are done running
pid=$!
wait $pid


#This will Unfold the JER root files that I have made with my analysis code
for i in {0..8}; do
  root -l -b 'UnfoldingJetsCopy.C('$unfoldpp','$unfoldPbPb','$jet_radius','$regCentBins','$numIter', '$etaRange', '$i',"'$jer_or_jes'",'$systematicUnfd', '$jetRate2015Bins','$rAA2015Binning')' &

done

#We will apply this code so that it will wait until the above commands are done running
pid=$!
wait $pid


#We will unfold to get the NOMINAL pT spectras
jer_or_jes=""

root -l -b 'UnfoldingJetsCopy.C('$unfoldpp','$unfoldPbPb','$jet_radius','$regCentBins','$numIter', '$etaRange', -1,"'$jer_or_jes'",'$systematicUnfd', '$jetRate2015Bins','$rAA2015Binning')' &
