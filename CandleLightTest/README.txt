#############################################################################################
################################### TLANĒX-TLI Test #########################################
############################################################################################


-----Purpose of test: Make the 2018 & 2015 comparison plots for the R=0.4 jets. If everything goes well, then we should see the same level of agreement that we saw previously before making any changes to the analysis code.


------>First Step: Execute the 'first_step.sh' bash script. This bash script will copy all the necesary files to run this test (i.e. analysis code, header files) on to the local git repository. It will also commit and push all the files on to github.

------>Second Step:There are three componenets to this test for when the following bash script is executed—candle_light_test.sh. 
	      1. It will run over pp and Pb+Pb data to make the raw distributions using the Jet Rate OR RAA 2015 binning. In this bash script there are important two variables 'jetRate2015Bins' and 'rAA2015Binning', which you can set to true or false depending what binning you would like to use. 
	      2. Once it is done making all the necessary raw pT spectras it will begin to UNFOLD the raw distributions using the variated jet pT (this is for the JER and JES systematic uncertainties) and also produce the nominal distributions.
	      3. This last step is to calculate the systematic uncertaintities for JES, JER, and Unfolding (MC finite stats and pT shape weighting).   


     
