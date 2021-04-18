use "G:\Github\TimeSeries\Research\MergedTest.dta"
generate date =tq(2011q3) +_n-1
format %tq date
tsset date



matrix A = (1,0,0,0 \ .,1,0,0 \ .,.,1,0 \ .,.,.,1)

matrix B = (.,0,0,0 \ 0,.,0,0 \ 0,0,.,0 \ 0,0,0,.)

matrix C = (.,0,0,0 \ .,.,0,0 \ .,.,.,0 \ .,.,.,.)


svar ADJGoldSupply ADJGDAverage ADJSP500 ADJGPrice , lags(1/4) aeq(A) beq(B)
irf create svar1, step(4) set(myGraph1, replace)

irf cgraph (svar1 ADJGoldSupply ADJGoldSupply sirf) (svar1 ADJGoldSupply ADJGDAveragesirf) (svar1 ADJGoldSupply ADJSP500 sirf) (svar1 ADJGoldSupply ADJGPrice sirf)(svar1ADJGDAverage ADJGoldSupply sirf) (svar1 ADJGDAverage ADJGDAverage sirf) (svar1 ADJGDAverage ADJSP500 sirf) (svar1 ADJGDAverage ADJGPrice sirf)(svar1 ADJSP500 ADJGoldSupply sirf)(svar1 ADJSP500 ADJGDAverage sirf) (svar1 ADJSP500 ADJSP500 sirf) (svar1 ADJSP500 ADJGPrice sirf)(svar1 ADJGPrice ADJGoldSupply sirf) (svar1 ADJGPrice ADJGDAverage sirf) (svar1 ADJGPrice ADJSP500 sirf) (svar1 ADJGPrice ADJGPrice sirf),title ("Irfsof VAR Model 2", size(vsmall))
