# streamflow data from USGS stream gauge 10343500

library(tidyverse)
library(lubridate)
library(chron)

sf = read_csv("streamflow.csv", col_names=FALSE)
colnames(sf) = c("Timestamp","Flow_cfs")

sf$Date = as.Date(sf$Timestamp, format="%m/%d/%Y")
sf_daily = sf %>% group_by(Date) %>% summarize_at(vars(Flow_cfs), funs(sum))
sf_daily$month = as.numeric(format(sf_daily$Date, "%m"))
sf_daily$year = as.numeric(format(sf_daily$Date, "%Y"))
sf_daily$wy = with(sf_daily, ifelse(month>=10, year+1, year))

sf_annual = sf_daily %>% group_by(wy) %>% summarize_at(vars(Flow_cfs), funs(sum))

ggplot(sf_daily, aes(x=Date, y=Flow_cfs))+geom_line()

ggplot(sf_annual, aes(x=wy, y=Flow_cfs))+geom_col()
