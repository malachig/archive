library(rjson)
library(ggplot2)
library(RCurl)
library(lubridate)

#define source URLs for data from covidtracking.org API
url_us_daily = "https://covidtracking.com/api/us/daily"
eurl_us_daily <- URLencode(url_us_daily)

#load data as JSON
data_us_daily <- fromJSON(getURL(eurl_us_daily))

#replace nulls with 0
for(i in 1:length(data_us_daily)){
  if (is.null(data_us_daily[[i]]$death)){
    data_us_daily[[i]]$death = 0
  }
}

#convert lists to dataframes
df <- lapply(data_us_daily, function(day) # Loop through each "day"
{
  # Convert each group to a data frame.
  data.frame(matrix(unlist(day), ncol=8, byrow=T))
})

#now you have a list of data frames, connect them together in one single dataframe
df <- do.call(rbind, df)

#overwrite the old dataframe with the clean one
data_us_daily = df

#add column names
names(data_us_daily) = c("Date","StatesTracked","Positive","Negative","Pos+Neg","Pending","Deaths","TotalTests")

#determine missing per day stats (e.g. the number of new tests performed each day)
data_us_daily[,"NewTests"] = 0
current_tests = 0

data_us_daily[,"NewDeaths"] = 0
current_deaths = 0

data_us_daily[,"NewPositives"] = 0
current_positives = 0

data_us_daily[,"NewPosNeg"] = 0
current_posneg = 0

for(i in 1:length(data_us_daily$Date)){
  new_tests = data_us_daily[i,"TotalTests"] - current_tests
  data_us_daily[i,"NewTests"] = new_tests
  current_tests = data_us_daily[i,"TotalTests"]
  
  new_deaths = data_us_daily[i,"Deaths"] - current_deaths
  data_us_daily[i,"NewDeaths"] = new_deaths
  current_deaths = data_us_daily[i,"Deaths"]
  
  new_positives = data_us_daily[i,"Positive"] - current_positives
  data_us_daily[i,"NewPositives"] = new_positives
  current_positives = data_us_daily[i,"Positive"]
  
  new_posneg = data_us_daily[i,"Pos+Neg"] - current_posneg
  data_us_daily[i,"NewPosNeg"] = new_posneg
  current_posneg = data_us_daily[i,"Pos+Neg"]
}

#determine daily and cumulative rates
data_us_daily[,"DailyPosTestRate"] = data_us_daily[,"NewPositives"] / data_us_daily[,"NewPosNeg"]
data_us_daily[,"CumPosTestRate"] = data_us_daily[,"Positive"] / data_us_daily[,"Pos+Neg"]
data_us_daily[,"CumDeathRate"] = data_us_daily[,"Deaths"] / data_us_daily[,"Positive"]

#plot the number of new tests performed each day
p1 = ggplot(data=data_us_daily, aes(x=ymd(Date), y=NewTests, group=1))
p1 = p1 + geom_line()
p1 = p1 + geom_point()
p1 = p1 + theme_bw()
p1 = p1 + theme(axis.text.x=element_text(angle=45, hjust=1)) + theme(plot.title = element_text(hjust = 0.5))
p1 = p1 + ggtitle("Number of COVID-19 tests performed per day in US")
p1 = p1 + ylab("New Tests Performed") + xlab("Date")
p1 = p1 + scale_x_date(date_breaks = "1 day", date_labels =  "%b %d") 
print(p1)

#plot the positive test rate each day (new positive / new positive+negative)
p2 = ggplot(data=data_us_daily, aes(x=ymd(Date), y=DailyPosTestRate, group=1))
p2 = p2 + geom_line()
p2 = p2 + geom_point()
p2 = p2 + theme_bw()
p2 = p2 + theme(axis.text.x=element_text(angle=45, hjust=1)) + theme(plot.title = element_text(hjust = 0.5))
p2 = p2 + ggtitle("Positive test rate observed per day in US")
p2 = p2 + ylab("Daily positive test rate (new pos / new pos+neg)") + xlab("Date")
p2 = p2 + scale_x_date(date_breaks = "1 day", date_labels =  "%b %d") 
print(p2)

#plot the cumulative positive test rate for each day (positive / positive+negative)
p3 = ggplot(data=data_us_daily, aes(x=ymd(Date), y=CumPosTestRate, group=1))
p3 = p3 + geom_line()
p3 = p3 + geom_point()
p3 = p3 + theme_bw()
p3 = p3 + theme(axis.text.x=element_text(angle=45, hjust=1)) + theme(plot.title = element_text(hjust = 0.5))
p3 = p3 + ggtitle("Cumulative positive test rate observed in US")
p3 = p3 + ylab("Cumulative positive test rate (pos / pos+neg)") + xlab("Date")
p3 = p3 + scale_x_date(date_breaks = "1 day", date_labels =  "%b %d") 
print(p3)

#plot the new deaths each day 
p4 = ggplot(data=data_us_daily, aes(x=ymd(Date), y=NewDeaths, group=1))
p4 = p4 + geom_line()
p4 = p4 + geom_point()
p4 = p4 + theme_bw()
p4 = p4 + theme(axis.text.x=element_text(angle=45, hjust=1)) + theme(plot.title = element_text(hjust = 0.5))
p4 = p4 + ggtitle("New deaths observed in US")
p4 = p4 + ylab("New deaths") + xlab("Date")
p4 = p4 + scale_x_date(date_breaks = "1 day", date_labels =  "%b %d") 
print(p4)

#plot the cumulative death rate (cumulative deaths / cumulative positives)
p5 = ggplot(data=data_us_daily, aes(x=ymd(Date), y=CumDeathRate, group=1))
p5 = p5 + geom_line()
p5 = p5 + geom_point()
p5 = p5 + theme_bw()
p5 = p5 + theme(axis.text.x=element_text(angle=45, hjust=1)) + theme(plot.title = element_text(hjust = 0.5))
p5 = p5 + ggtitle("Cumulative death rate observed in US")
p5 = p5 + ylab("Observed death rate (total deaths / total positives)") + xlab("Date")
p5 = p5 + scale_x_date(date_breaks = "1 day", date_labels =  "%b %d") 
print(p5)
