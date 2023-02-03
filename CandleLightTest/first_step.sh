#make sure we update the code that I will be using in this test into git repository
cp /usatlas/u/bereniceg299/data/LargeRJet_Study/NewSourceCode/LargeRAnalysis.C LargeRAnalysis/LargeRAnalysis.C #Analysis code
cp /usatlas/u/bereniceg299/data/LargeRJet_Study/NewSourceCode/NewTreeVariables.h LargeRAnalysis/NewTreeVariables.C #Header File
cp /usatlas/u/bereniceg299/data/LargeRJet_Study/NewSourceCode/SysUncert/sysUncert.C LargeRAnalysis/sysUncert.C #Calc Systematcis
cp /usatlas/u/bereniceg299/data/LargeRJet_Study/NewSourceCode/SysUncert/UnfoldingJetsCopy.C LargeRAnalysis/sysUncert.C #Unfolding Code
cp /usatlas/u/bereniceg299/data/LargeRJet_Study/NewSourceCode/SysUncert/UnfoldingJetsCopy.C LargeRAnalysis/bNec.h #Header file
#git
git add *
git commit -m "pushing code to GitHub before running cadle test"
git push
