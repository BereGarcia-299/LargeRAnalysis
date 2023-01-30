#int UnfoldingJetsCopy(bool Unfoldpp = false, bool UnfoldPbPb = true, int jet_radius = R4,bool regCentBins = false,int numIter = 17, int etaRange = 28, int uncrt_sys = 2, string jer_or_jes = "JES",bool systematicUnfd = false, bool jetRate2015Bins = false,bool rAA2015Binning = true){


#Setting Variables
unfoldpp=true
unfoldPbPb=false
jet_radius=1 #R=0.4
regCentBins=false 
numIter=17
etaRange=28
uncrt_sys=0 #variation
jer_or_jes="JER"
systematicUnfd=false #apply unfolding systematic (finite MC stats)
jetRate2015Bins=false #2015 jet rate bins
rAA2015Binning=true #2015 RAA bins


root -l -b 'UnfoldingJetsCopy.C('$unfoldpp','$unfoldPbPb','$jet_radius','$regCentBins','$numIter', '$etaRange', 0,"'$jer_or_jes'",'$systematicUnfd', '$jetRate2015Bins','$rAA2015Binning')' &
root -l -b 'UnfoldingJetsCopy.C('$unfoldpp','$unfoldPbPb','$jet_radius','$regCentBins','$numIter', '$etaRange', 1,"'$jer_or_jes'",'$systematicUnfd', '$jetRate2015Bins','$rAA2015Binning')' &
root -l -b 'UnfoldingJetsCopy.C('$unfoldpp','$unfoldPbPb','$jet_radius','$regCentBins','$numIter', '$etaRange', 2,"'$jer_or_jes'",'$systematicUnfd', '$jetRate2015Bins','$rAA2015Binning')' &
root -l -b 'UnfoldingJetsCopy.C('$unfoldpp','$unfoldPbPb','$jet_radius','$regCentBins','$numIter', '$etaRange', 3,"'$jer_or_jes'",'$systematicUnfd', '$jetRate2015Bins','$rAA2015Binning')' &
root -l -b 'UnfoldingJetsCopy.C('$unfoldpp','$unfoldPbPb','$jet_radius','$regCentBins','$numIter', '$etaRange', 4,"'$jer_or_jes'",'$systematicUnfd', '$jetRate2015Bins','$rAA2015Binning')' &
root -l -b 'UnfoldingJetsCopy.C('$unfoldpp','$unfoldPbPb','$jet_radius','$regCentBins','$numIter', '$etaRange', 5,"'$jer_or_jes'",'$systematicUnfd', '$jetRate2015Bins','$rAA2015Binning')' &
root -l -b 'UnfoldingJetsCopy.C('$unfoldpp','$unfoldPbPb','$jet_radius','$regCentBins','$numIter', '$etaRange', 6,"'$jer_or_jes'",'$systematicUnfd', '$jetRate2015Bins','$rAA2015Binning')' &
root -l -b 'UnfoldingJetsCopy.C('$unfoldpp','$unfoldPbPb','$jet_radius','$regCentBins','$numIter', '$etaRange', 7,"'$jer_or_jes'",'$systematicUnfd', '$jetRate2015Bins','$rAA2015Binning')' &
root -l -b 'UnfoldingJetsCopy.C('$unfoldpp','$unfoldPbPb','$jet_radius','$regCentBins','$numIter', '$etaRange', 8,"'$jer_or_jes'",'$systematicUnfd', '$jetRate2015Bins','$rAA2015Binning')' &
# root -l -b 'UnfoldingJetsCopy.C('$unfoldpp','$unfoldPbPb','$jet_radius','$regCentBins','$numIter', '$etaRange', 9,"'$jer_or_jes'",'$systematicUnfd', '$jetRate2015Bins','$rAA2015Binning')' &
# root -l -b 'UnfoldingJetsCopy.C('$unfoldpp','$unfoldPbPb','$jet_radius','$regCentBins','$numIter', '$etaRange', 10,"'$jer_or_jes'",'$systematicUnfd', '$jetRate2015Bins','$rAA2015Binning')' &
# root -l -b 'UnfoldingJetsCopy.C('$unfoldpp','$unfoldPbPb','$jet_radius','$regCentBins','$numIter', '$etaRange', 11,"'$jer_or_jes'",'$systematicUnfd', '$jetRate2015Bins','$rAA2015Binning')' &
# root -l -b 'UnfoldingJetsCopy.C('$unfoldpp','$unfoldPbPb','$jet_radius','$regCentBins','$numIter', '$etaRange', 12,"'$jer_or_jes'",'$systematicUnfd', '$jetRate2015Bins','$rAA2015Binning')' &
# root -l -b 'UnfoldingJetsCopy.C('$unfoldpp','$unfoldPbPb','$jet_radius','$regCentBins','$numIter', '$etaRange', 13,"'$jer_or_jes'",'$systematicUnfd', '$jetRate2015Bins','$rAA2015Binning')' &
# root -l -b 'UnfoldingJetsCopy.C('$unfoldpp','$unfoldPbPb','$jet_radius','$regCentBins','$numIter', '$etaRange', 14,"'$jer_or_jes'",'$ystematicUnfd', '$jetRate2015Bins','$rAA2015Binning')' &
# root -l -b 'UnfoldingJetsCopy.C('$unfoldpp','$unfoldPbPb','$jet_radius','$regCentBins','$numIter', '$etaRange', 15,"'$jer_or_jes'",'$systematicUnfd', '$jetRate2015Bins','$rAA2015Binning')' &
# root -l -b 'UnfoldingJetsCopy.C('$unfoldpp','$unfoldPbPb','$jet_radius','$regCentBins','$numIter', '$etaRange', 16,"'$jer_or_jes'",'$systematicUnfd', '$jetRate2015Bins','$rAA2015Binning')' &
# root -l -b 'UnfoldingJetsCopy.C('$unfoldpp','$unfoldPbPb','$jet_radius','$regCentBins','$numIter', '$etaRange', 17,"'$jer_or_jes'",'$systematicUnfd', '$jetRate2015Bins','$rAA2015Binning')' &
# root -l -b 'UnfoldingJetsCopy.C('$unfoldpp','$unfoldPbPb','$jet_radius','$regCentBins','$numIter', '$etaRange', 18,"'$jer_or_jes'",'$systematicUnfd', '$jetRate2015Bins','$rAA2015Binning')' &
# root -l -b 'UnfoldingJetsCopy.C('$unfoldpp','$unfoldPbPb','$jet_radius','$regCentBins','$numIter', '$etaRange', 19,"'$jer_or_jes'",'$systematicUnfd', '$jetRate2015Bins','$rAA2015Binning')' &
#root -l -b 'UnfoldingJetsCopy.C('$unfoldpp','$unfoldPbPb','$jet_radius','$regCentBins','$numIter', '$etaRange', 20,"'$jer_or_jes'",'$systematicUnfd', '$jetRate2015Bins','$rAA2015Binning')' &

#root -l -b -q 'draw_result_v5.C("'$ver'", "jetEnergyScale19", "bayes")' &
#root -l -b -q 'draw_result_v5.C("'$ver'", "jetEnergyScale20", "bayes")' &
#root -l -b -q 'draw_result_v5.C("'$ver'", "jetEnergyScale21", "bayes")' &
#root -l -b -q 'draw_result_v5.C("'$ver'", "jetEnergyScale22", "bayes")' &
#root -l -b -q 'draw_result_v5.C("'$ver'", "bkgPhoSub3", "bayes")' &
#root -l -b -q 'draw_result_v5.C("'$ver'", "bkgPhoSub4", "bayes")' &
#root -l -b -q 'draw_result_v5.C("'$ver'", "bkgPhoSub5", "bayes")' &
