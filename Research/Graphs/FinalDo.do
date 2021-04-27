
 ** use "C:\Users\gaeta\Documents\GitHub\TimeSeries\Research\MergedTest.dta"
use "G:\Github\TimeSeries\Research\MergedTest.dta"
 
generate date =tq(2011q3) +_n-1
format %tq date
list date in 1/5
tsset date



matrix A = (1,0,0,0 \ .,1,0,0 \ .,.,1,0 \ .,.,.,1)

matrix B = (.,0,0,0 \ 0,.,0,0 \ 0,0,.,0 \ 0,0,0,.)

matrix C = (.,0,0,0 \ .,.,0,0 \ .,.,.,0 \ .,.,.,.)

 
svar ADJGoldSupply ADJGDAverage ADJSP500 ADJGPrice , lags(1/4) aeq(A) beq(B)
irf create svar1, step(4) set(GraphGoldPriceLast, replace)

irf cgraph (svar1 ADJGoldSupply ADJGoldSupply sirf) (svar1 ADJGoldSupply ADJGDAverage sirf) (svar1 ADJGoldSupply ADJSP500 sirf) (svar1 ADJGoldSupply ADJGPrice sirf)(svar1 ADJGDAverage ADJGoldSupply sirf) (svar1 ADJGDAverage ADJGDAverage sirf) (svar1 ADJGDAverage ADJSP500 sirf) (svar1 ADJGDAverage ADJGPrice sirf)(svar1 ADJSP500 ADJGoldSupply sirf)(svar1 ADJSP500 ADJGDAverage sirf) (svar1 ADJSP500 ADJSP500 sirf) (svar1 ADJSP500 ADJGPrice sirf)(svar1 ADJGPrice ADJGoldSupply sirf) (svar1 ADJGPrice ADJGDAverage sirf) (svar1 ADJGPrice ADJSP500 sirf) (svar1 ADJGPrice ADJGPrice sirf),title ("IRF With Gold Price Shock as Last", size(vsmall))

 *** graph export graph.png, name("IRF With Gold Price Shock as Last") figure this out 
 
 svar ADJGoldSupply ADJGDAverage ADJGPrice ADJSP500, lags(1/4) aeq(A) beq(B)
irf create svar1, step(4) set(GraphSP500Last , replace)

irf cgraph (svar1 ADJGoldSupply ADJGoldSupply sirf) (svar1 ADJGoldSupply ADJGDAverage sirf) (svar1 ADJGoldSupply ADJSP500 sirf) (svar1 ADJGoldSupply ADJGPrice sirf)(svar1 ADJGDAverage ADJGoldSupply sirf) (svar1 ADJGDAverage ADJGDAverage sirf) (svar1 ADJGDAverage ADJSP500 sirf) (svar1 ADJGDAverage ADJGPrice sirf)(svar1 ADJGPrice ADJGoldSupply sirf) (svar1 ADJGPrice ADJGDAverage sirf) (svar1 ADJGPrice ADJSP500 sirf) (svar1 ADJGPrice ADJGPrice sirf)(svar1 ADJSP500 ADJGoldSupply sirf)(svar1 ADJSP500 ADJGDAverage sirf) (svar1 ADJSP500 ADJSP500 sirf) (svar1 ADJSP500 ADJGPrice sirf),title ("IRF With S&P500 Shock as Last", size(vsmall))