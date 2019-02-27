rhessysver = "/PATH/rhessys7.0"

setwd("/PATH/scripts")

# Load these libraries

library(pse)
library(tidyverse)
library(RHESSysIOinR)
library(RHESSysPreprocessing)
library(lubridate)

# define a function - nse
NSE = function (m, o)
{
  err = m - o
  meanobs = mean(o)
  mse = sum(err*err)
  ovar = sum((o-meanobs)*(o-meanobs))
  nse = ifelse(ovar > 0, 1.0 - mse/ovar,0)
  nse
}


# set period for the calibration simulation
startyr = YEAR 
endyr = YEAR 

# startoutyear - when output should start printing
startoutyr = YEAR 
nyrs = endyr-startoutyr

# generate a tec file
tecfile = "../tecfiles/tec.cal"
strg = sprintf("%d 10 1 1 print_daily_on \n%d 10 1 2 print_daily_growth_on",
               startoutyr-1, startoutyr-1)
write(strg, file=tecfile)

# generate parameter sets

# how many parameter sets
nsets=20

# list parameters
factors = c("m","K","pa","po")
# type of distributions they arise from
q = c("qlnorm", "qunif", "qunif","qunif")
# parameters for those distributions
q.arg = list(list(meanlog=1,sd=2), list(min=0.9, max=1000), list(min=0.5, max=2.0), list(min=0.5, max=2.0))


# you can choose different distributions and different parameter sets
# for example

# parameters for those distributions
# factors = c("m","K")
#q = c("qnorm", "qunif")
#q.arg = list(list(mean=0.5,sd=0.01), list(min=100, max=200))
# or for example
#factors = c("m","K","gw1","gw2")
#q = c("qunif","qunif","qunif","qunif")
#q.arg = list(list(min=0.2, max=10), list(min=1, max=1000), list(min=0,max=0.3), list(min=0.5, max=0.99))

#LHS = latin hypercube sample
sens = LHS(NULL,factors,nsets,q,q.arg, nboot=500)
sens.pars = get.data(sens)

# check to make sure values look reasonable
summary(sens.pars)

# set a unique scenarios number for each parameter set
sens.pars$scen = seq(from=1,to=length(sens.pars$m))

# create some data structures to store results
# first one to store simulation output 
# these will just be daily values that you want to keep for ALL simulations, so be cautious
# ecovars = c("scen","date","streamflow") might be best
ecovars = c("scen", "date","precip", "streamflow","trans","evap","psn","lai")


# create a data structure to store calibration metric results 
# expand this to add additional metrics
calmetrics = as.data.frame(matrix(nrow=nsets,ncol=7))
colnames(calmetrics)=c("scen","nse","nselog","rmse","annualbias","perr","meanerr_minmonth")


# readin some observed data
# this is just an example remember RHESSys streamflow is mm/day
obsflow = read.csv("../../flow.csv")
colnames(obsflow)=c("ID","PARAM","bdate","value","err")
obsflow$date = mdy(obsflow$bdate)
obsflow$year = as.integer(year(obsflow$date))
obsflow$month = as.integer(month(obsflow$date))
obsflow$day = as.integer(day(obsflow$date))
obsflow = mkdate(obsflow)

# convert to mm/day divide by basin area
basin.area = 209 #km2
stoday = 60*60*24
obsflow$mm = obsflow$value * stoday / (basin.area * 1000 * 1000)*1000

ggplot(obsflow, aes(date, mm))+geom_line()

obs_wy = obsflow %>% group_by(wy) %>% summarize(totalflow=sum(mm))
ggplot(obs_wy, aes(wy, totalflow))+geom_col()


# set the worldfile, flowtable, header, make sure paths make sense relative to working directory
worldfile = "../worldfiles/NAME.world"
flowtable = "../flowtables/NAME.flow"
headrh =  "../worldfiles/NAME.hdr"

# run rhessys for all the parameters sets
start=1
for (j in 1:nsets)  {


cmd1 = sprintf("%s -t %s -w %s -r %s", rhessysver, tecfile, worldfile, flowtable)
cmd2 = sprintf("-pre ../out/test   -s %f %f -sv %f %f -gw 0.0 1.0 -svalt %f %f", 
sens.pars$m[j], sens.pars$K[j], sens.pars$m[j], sens.pars$K[j], sens.pars$pa[j], sens.pars$po[j])

# you will need to edit this line if you are changing other parameters for example
#cmd2 = sprintf("-pre ../out/test   -s %f %f -sv %f %f -gw 0.0 1.0 -svalt 5.0 1.0", sens.pars$m[j], sens.pars$K[j], sens.pars$m[j], sens.pars$K[j])
cmd3 = sprintf("-st %d 10 1 1 -ed %d 10 1 1 -b  -whdr %s -b -tchange 0 0 -climrepeat -b -g", startyr, endyr, headrh)

cmdall = paste(cmd1,cmd2,cmd3)

# as an alterantive to system you could write all runs to a script and then source that
#write(cmdall, file="newscrpt", append=FALSE)
system(cmdall)

a = readin_rhessys_output("../out/test", g=1, c=1, wy=0)
compa = inner_join(a$bd, obsflow)

# add additional metrics here if needed
calmetrics[j,"rmse"] = sqrt(sum( (compa$streamflow-compa$mm)*(compa$streamflow-compa$mm)) /
                        length(compa$wy))
calmetrics[j,"scen"]=j
calmetrics[j,"nse"] = NSE(m=compa$streamflow, o=compa$mm)
calmetrics[j,"lognse"]=NSE(m=log(compa$streamflow+0.000001), o=log(compa$mm+0.000001))

calmetrics[j,"pbias"] = sum(compa$streamflow-compa$mm)/sum(compa$mm)*100

calmetrics[j,"annualbias"] = (mean(compa$streamflow)-mean(compa$mm))*365

# min monthly flow cal 
compa_mwy = compa %>% group_by(wy, month) %>% summarize_all(funs(mean))
compa_mth = compa_mwy %>% group_by(month) %>% summarize(flow=mean(mm))
minmonth = compa_mth[which.min(compa_mth$flow),"month"]$month
tmp = subset(compa_mwy, month == minmonth)
calmetrics[j,"meanerr_minmonth"] = mean(tmp$streamflow-tmp$mm)


#  if you want to save some of the results
# if not comment all of this out - note calres can be large for long simulation times and/or large nsets
if (j==1) {
  simlen = length(a$bd$year)
  calres = as.data.frame(matrix(nrow=nsets*simlen, ncol=length(ecovars)))
  colnames(calres)=ecovars
  calres$date = ymd(a$bd$date)

}

endj = start+length(a$bd$wy)-1
calres$scen[start:endj] = j
calres$wy[start:endj] = a$bd$wy
calres$date[start:endj]=ymd(a$bd$date)
calres$streamflow[start:endj] = a$bd$streamflow
calres$lai[start:endj] = a$bd$lai
calres$psn[start:endj] = a$bd$psn
calres$evap[start:endj] = a$bd$evap
calres$trans[start:endj] = a$bd$trans
calres$precip[start:endj] = a$bd$precip

start = endj+1
}

# this is the end of the calibration

# some output examples
ggplot(calres, aes(date, trans, col=as.factor(scen)))+geom_line()+
  labs(y="Daily Transpiration mm/day",x="Date")

calres$year = year(calres$date)
calres$month = as.integer(month(calres$date))
calres$day = day(calres$date)
calres = mkdate(calres)

ggplot(calres, aes(wy, psn, col=as.factor(scen)))+stat_summary(fun.y=sum, geom="line")+
  labs(y="Total Annual Psn (kg/m2/yr)",x="Water Year")

calres$month = months(calres$date)
ggplot(calres, aes(as.factor(month), trans))+geom_boxplot()

calresc = left_join(calres, sens.pars)
calresb = left_join(calresc,obsflow[,c("mm","date")], by="date")
str_wy = calresb %>% group_by(scen,wy) %>% summarize(streamflow=sum(streamflow),
                                                    trans=sum(trans), precip=sum(precip), obsstr=sum(mm))
ggplot(str_wy, aes(as.factor(wy), streamflow-obsstr, fill=as.factor(scen)))+geom_col(position="dodge")+
  labs(y="Annual Flow Err mm/yr")

ggplot(calresc, aes(m, streamflow-mm))+geom_points()+labs(main="Daily Error by m parameter")
ggplot(calresc, aes(as.factor(round(K/10)), streamflow-mm))+geom_boxplot()+labs(main="Daily Error by K parameter")

ggplot(calresc, aes(po, lai) )+stat_summary(fun.y=mean, geom="point", cex=2, col="red")+
  labs(main="Mean Simulation LAI across PO parameter")




save.image("cal.RData")




