library(readxl)
library(forecast)
library(timeSeries)
library(stats)
library(dplyr)
library(tseries)

GoldSupply <- read_excel("GitHub/TimeSeries/Research/GoldSupply.xlsx",range = "A1:B45")

Gold_Demand <- read_excel("GitHub/TimeSeries/Research/GoldDemandClean.xlsx")

SP500 <- read_excel("GitHub/TimeSeries/Research/SP500.xls",range = "A11:B50")

Gold_Price <- read_excel("GitHub/TimeSeries/Research/GOLDAMGBD228NLBM.xls",range = "A11:B56")

# Q1- 2010 to Q4-2020 Gold Supply
# Q1-2005 to Q4-2020 Gold Deman
# Q1- 2011 to Q1-2021 SP500
# Q1 - 2010 to Q1-2021 Gold_Price

# Dropping columns to make the same start and end date for the time series. 
GoldSupply <- GoldSupply[-c(1:6),]
Gold_Demand <- Gold_Demand[-c(1:26),]
SP500 <- SP500[-c(39),]
Gold_Price <- Gold_Price[-c(1:6,45),]

# Transform data into TS
GoldSupply <- ts(GoldSupply$`Mine Production`,start = c(2011,3),frequency = 4)
Gold_DemandAverage <- ts(Gold_Demand$`Average of Indexes`,start = c(2011,3),frequency = 4)
Gold_DemandSum <- ts(Gold_Demand$`Sum of Indexes`,start = c(2011,3),frequency = 4)
SP500 <- ts(SP500$SP500_PCH,start = c(2011,3),frequency = 4)
Gold_Price <- ts(Gold_Price$GOLDAMGBD228NLBM,start = c(2011,3),frequency = 4)

# plotting graphs
plot.ts(GoldSupply) # time trend and seasonality 
plot.ts(Gold_DemandAverage) # time trend
plot.ts(Gold_DemandSum) # time trend
plot.ts(SP500) #time trend possibly seasonality?
plot.ts(Gold_Price) #time trend?

# Let us decompose our time series to find evidence of time trend or seasonality
GoldSupplyDecomp<- decompose(GoldSupply)
Gold_DemandAverageDecomp <- decompose(Gold_DemandAverage)
Gold_DemandSumDecomp <- decompose(Gold_DemandSum)
SP500Decomp <- decompose(SP500)
Gold_PriceDecomp <- decompose(Gold_Price)

# Let us plot the decomposed graphs
plot(GoldSupplyDecomp) # confirm time trend maybe seasonality?
plot(Gold_DemandAverageDecomp) #confirm time trend maybe seasonality?
plot(Gold_DemandSumDecomp) # same as above
plot(SP500Decomp) # no real trend but there is a seasonal aspect
plot(Gold_PriceDecomp) # unkown

# Let us test to see if there our data is stationary
DF <- tseries::adf.test(GoldSupply, k = 0) #no unit root 
ADF2 <- tseries::adf.test(GoldSupply, k = 2) #exists a unit root 
ADF <- tseries::adf.test(GoldSupply) #Unit root

DF2 <- tseries::adf.test(Gold_DemandAverage, k = 0) #Unit root
ADF22 <- tseries::adf.test(Gold_DemandAverage, k = 2) # Unit Root
ADF2 <- tseries::adf.test(Gold_DemandAverage) # Unit Root

DF3 <- tseries::adf.test(Gold_DemandSum, k = 0) #Unit root
ADF23 <- tseries::adf.test(Gold_DemandSum, k = 2) # Unit Root
ADF3 <- tseries::adf.test(Gold_DemandSum) # Unit Root

DF4 <- tseries::adf.test(SP500, k = 0) # No unit root
ADF24 <- tseries::adf.test(SP500, k = 2) # Unit Root
ADF4 <- tseries::adf.test(SP500) # Unit Root

DF5 <- tseries::adf.test(Gold_Price, k = 0) # Unit root
ADF25 <- tseries::adf.test(Gold_Price, k = 2) # Unit Root
ADF5 <- tseries::adf.test(Gold_Price) # Unit Root


# Going to adjust for seasonality and time trend
ADJGoldSupply <- GoldSupply - GoldSupplyDecomp$seasonal - GoldSupplyDecomp$trend
ADJGDAverage <- Gold_DemandAverage- Gold_DemandAverageDecomp$seasonal - Gold_DemandAverageDecomp$trend
ADJGDSum <- Gold_DemandSum - Gold_DemandSumDecomp$seasonal - Gold_DemandSumDecomp$trend
ADJSP500 <- SP500 - SP500Decomp$seasonal - SP500Decomp$trend
ADJGPrice <- Gold_Price - Gold_PriceDecomp$seasonal - Gold_PriceDecomp$trend

# Replotting the Graphs
plot.ts(ADJGoldSupply)
plot.ts(ADJGDAverage)
plot.ts(ADJGDSum)
plot.ts(ADJSP500)
plot.ts(ADJGPrice)

# Everything seems good.Exporting the clean datasets to a new data so that I can perform my analysis. 

#### Warning: You will have to adjust the directory in order to save this data ####
write.csv2(ADJGoldSupply, file = "C:/Users/gaeta/Documents/GitHub/TimeSeries/Research/AdjustedGoldSupply.csv")
write.csv2(ADJGDAverage, file = "C:/Users/gaeta/Documents/GitHub/TimeSeries/Research/AdjustedGoldDemandAverage.csv")
write.csv2(ADJGDSum, file = "C:/Users/gaeta/Documents/GitHub/TimeSeries/Research/AdjustedGoldDemandSum.csv")
write.csv2(ADJSP500, file = "C:/Users/gaeta/Documents/GitHub/TimeSeries/Research/AdjustedS&P500.csv")
write.csv2(ADJGPrice, file = "C:/Users/gaeta/Documents/GitHub/TimeSeries/Research/AdjustedGoldPrice.csv")
