
***generate VARs***
clear

*Import data
import excel "/Users/USER/Documents/TimeSeries/Project/project_data_full_csv.csv", sheet("FRED Graph") firstrow
*************** PART A*******************
tsset index
tsline ln_Y
tsline CPI
tsline ffr
***Plotting untransformed series****
generate d_oil=d.oil_supply
generate d_oilp = d.ln_oil

**** Graphing the transformed series******
tsline d_oil
tsline econ_act
tsline d_oilp


********************Part B********************
** Run VAR on transformed series
*Short Run Restrictions

matrix A2 = (1,0,0 \ .,1,0 \ .,.,1)
matrix B2 = (.,0,0 \ 0,.,0 \ 0,0,.)

*I put d_lnY first.  

svar  d_oil econ_act d_oilp, lags(1/24) aeq(A2) beq(B2)

*********************************part C/D************************
irf create var2irf, step(24) set(var2,replace)
irf graph sirf, irf(var2irf) yline(0,lcolor(black)) 

svar d_oil econ_act d_oilp, lags(1/24) aeq(A2) beq(B2)

*Plot irfs
irf create var2irf, step(24) set(var2,replace)
irf graph sirf, irf(var2irf) yline(0,lcolor(black)) 
** Plotting those that don't appear in table** 

irf graph sirf, impulse(d_oil) response(econ_act)
irf graph sirf, impulse(d_oil) response(d_oilp)
irf graph sirf, impulse(econ_act) response(d_oilp)


svar ln_P ln_Y ffr, lags(1/8) aeq(A2) beq(B2)
irf create svar1, set(myGraph1, replace)
irf cgraph (svar1 ln_Y ln_Y sirf) (svar1 ln_Y ln_P sirf) (svar1 ln_Y ffr sirf) ///
           (svar1 ln_P ln_Y sirf) (svar1 ln_P ln_P sirf) (svar1 ln_P ffr sirf) ///
		   (svar1 ffr ln_Y sirf) (svar1 ffr ln_P sirf) (svar1 ffr ffr sirf), ///
		   title ("Irfs of VAR Model 2", size(vsmall))

		   
irf cgraph 
(svar1 ADJGoldSupply ADJGoldSupply sirf) (svar1 ADJGoldSupply ADJGDAverage sirf) (svar1 ADJGoldSupply ADJSP500 sirf) (svar1 ADJGoldSupply ADJGPrice sirf)
(svar1 ADJGDAverage ADJGoldSupply sirf) (svar1 ADJGDAverage ADJGDAverage sirf) (svar1 ADJGDAverage ADJSP500 sirf) (svar1 ADJGDAverage ADJGPrice sirf)
(svar1 ADJSP500 ADJGoldSupply sirf) (svar1 ADJSP500 ADJGDAverage sirf) (svar1 ADJSP500 ADJSP500 sirf) (svar1 ADJSP500 ADJGPrice sirf)
(svar1 ADJGPrice ADJGoldSupply sirf) (svar1 ADJGPrice ADJGDAverage sirf) (svar1 ADJGPrice ADJSP500 sirf) (svar1 ADJGPrice ADJGPrice sirf)
,title ("Irfs of VAR Model 2", size(vsmall))
	   



