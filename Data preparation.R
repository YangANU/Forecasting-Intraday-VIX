########################
# Tidy-up VIX tick data
########################
library(forecast)
library(rainbow)
library(ftsa)
library(doParallel)

# Read in data
data_raw = read.table("VIX_tick.txt", header=FALSE, sep=",")

data_date = as.data.frame(matrix(, nrow(data_raw), ncol(data_raw)-1))
data_date[,1] = as.Date(data_raw[,1], format = "%m/%d/%Y")

###############################
# Finding weeks based on years
###############################

# Add day of week 
dayofweek = weekdays(as.Date(data_raw[,1], format = "%m/%d/%Y"))
weeknum = as.numeric(format(as.Date(data_raw[,1], format = "%m/%d/%Y"), "%U"))
data_date[,2] = dayofweek
data_date[,3] = weeknum

# Add timing and price of each tick
data_date[,4:5] = data_raw[,2:3]
colnames(data_date) = c("Date", "Day of week", "Week number", "Time", "Tick price")

###########################################################################################
# Use data from the first trading day of January 2017 to the last trading day of June 2017
###########################################################################################

# index of the first and the last of selected data
day_first = which.max(data_date[,1] == "2017-01-03")
day_last = which.max(data_date[,1] == "2017-07-03") - 1

# extract the selected data from the full dataset
data_reduce = data_date[day_first:day_last,]

##########################
# manually fill in values
##########################

# Set maximum number of tick values per day
ndays = length(unique(data_reduce[,1]))
fts_data = matrix(, 1621, ndays)

# function used to select data according to date
select_date <- function(data_mat, date)
{
  return(data_mat[data_mat[,1] == date,])
}

# function used to format time into "%H:%M:%S"

form_time <- function(time)
{
  return(format(time, format = '%H:%M:%S'))
}

# starting point and ending point of trading times
start_data = vector(, ndays)
end_data = vector(, ndays)
len_select = vector(, ndays)
for(i in 1:ndays)
{
  date = unique(data_reduce[,1])[i]
  
  # take the first trading time that is larger than 09:30:00 as the starting point
  # take the last trading time that is less than 16:15:00 as the ending point
  
  start_data[i] = ifelse(form_time(select_date(data_reduce, date)[1,4])>="09:30:00",
                         form_time(select_date(data_reduce, date)[1,4]),
                         form_time(select_date(data_reduce, date)[,4])[which.max(form_time(select_date(data_reduce, date)[,4])>"09:30:00")])
  
  end_data[i] = ifelse(form_time(select_date(data_reduce,date)[dim(select_date(data_reduce,date))[1],4])<="16:15:00",
                       form_time(select_date(data_reduce,date)[dim(select_date(data_reduce,date))[1],4]),
                       form_time(select_date(data_reduce, date)[,4])[which.max(form_time(select_date(data_reduce, date)[,4])>"16:15:00")-1])
  
  # accounting for possible duplicated starting points and ending points
  
  len_select[i] = which.max(form_time(select_date(data_reduce, date)[,4]) == end_data[i]) - which.max(form_time(select_date(data_reduce, date)[,4]) > start_data[i]) + 2
  
  print(i)
}

  
# check if all starting points and ending points are within the specified interval

sum(start_data < "09:30:00")
sum(end_data > "16:15:00")
sum(len_select > 1621) 

# there are duplications for some days
dup_date_ind <- which(len_select > 1621)

# "09:30:00" converted to seconds equals to 34200; "16:15:00" converted to seconds equals to 58500
ticks_interval = format(as.POSIXct('0001-01-01 00:00:00') + seq(34200, 58500, 15), "%H:%M:%S") 

# fill in available data which do not have duplications

fill_data_nodup = setdiff(1:ndays, dup_date_ind)

for (i in 1:length(fill_data_nodup))
{
  
  # index of date which should be filled
  
  ind = fill_data_nodup[i]
  date = unique(data_reduce[,1])[ind]
  
  # indices of starting and ending points in the original matrix
  
  start_data_ind = which.max(form_time(select_date(data_reduce, date)[,4]) > start_data[ind]) - 1
  end_data_ind = which.max(form_time(select_date(data_reduce, date)[,4]) == end_data[ind])
  
  # indices of starting point in the target matrix
  
  start_fill_ind = which.max(start_data[ind] < ticks_interval)-1
  
  # time indices of selected data
  
  data_select_Tind = as.numeric(select_date(data_reduce, date)[start_data_ind:end_data_ind,4])
  
  # time indices of 15-sec integer multiple in the target matrix
  
  fill_Tind = seq(as.numeric(select_date(data_reduce,date)[start_data_ind,4]), 
      as.numeric(select_date(data_reduce,date)[end_data_ind,4]), 15)
  
  ### fill in the target matrix according to time indices, leaving unfilled position as NA
    
  # Step 1, fill in direct match values
  
  match_Tind = as.numeric(na.omit(match(data_select_Tind, fill_Tind)))
  match_target = match_Tind + start_fill_ind - 1
  match_origin = na.omit(match(fill_Tind[match_Tind], data_select_Tind)) + start_data_ind - 1
  
  fts_data[match_target, ind] = select_date(data_reduce, date)[match_origin, 5]
  
  # Step 2, fill in "nomatch" values observed within 5 seconds of any tick times in the target matrix
  
  if(length(match_origin) < (end_data_ind - start_data_ind + 1))
  {
    # index in the original matrix
    nomatch_origin_Tind = setdiff((start_data_ind:end_data_ind), (start_data_ind:end_data_ind)[na.omit(match(fill_Tind[match_Tind], data_select_Tind))])
    nomatch_target_Tind = as.numeric(attributes(na.omit(match(data_select_Tind, fill_Tind)))$na.action)  + start_fill_ind - 1
    
    # index in the target matrix
    for (j in 1:length(nomatch_origin_Tind))
    {
      # check if within 5 seconds of any later tick times
      
      nomatch_origin = nomatch_origin_Tind[j] 
      nomatch_target = which(abs(as.numeric(select_date(data_reduce, date)[nomatch_origin,4]) - fill_Tind) <= 5) + start_fill_ind - 1
      
      # fill in value if a match is found
      
      if(length(nomatch_target)> 0)
      {
        fts_data[nomatch_target ,ind] = select_date(data_reduce, date)[nomatch_origin, 5]
      }
    }
  }
  print(i)
}


# delete duplicate data and fill in trimmed data (not necessary if using six months of data)

for(i in 1:length(dup_date_ind))
{
  # exact date of duplicate data
  
  dup_date = unique(data_reduce[,1])[dup_date_ind[i]]
  fts_data[, dup_date_ind[i]] = rep(NA, 1621)
  
  # remove duplications
  
  dum = select_date(data_reduce, dup_date)
  dum_index = as.numeric(dum[,4])
  dum_dup_index = diff(dum_index)
  dum_trim = dum[-which(dum_dup_index == 0),]
  dum_trim_len = dim(dum_trim)[1]
  
  # take the first trading time that is larger than 09:30:00 as the starting point
	# take the last trading time that is less than 16:00:00 as the ending point
  
  start_trim = ifelse(form_time(dum_trim[1,4])>="09:30:00", form_time(dum_trim[1,4]),
							 form_time(dum_trim[,4])[which.max(form_time(dum_trim[,4])>"09:30:00")])
  end_trim = ifelse(form_time(dum_trim[dum_trim_len,4])<="16:15:00", form_time(dum_trim[dum_trim_len,4]),
						   form_time(dum_trim[,4])[which.max(form_time(dum_trim[,4])>"16:15:00")-1])
   
  # indices of starting and ending points in the original matrix
  
  start_trim_ind = which.max(form_time(dum_trim[,4]) > start_trim) - 1
	end_trim_ind = which.max(form_time(dum_trim[,4]) == end_trim)
  
  # indices of starting points in the target matrix 
  
  start_trim_fill_ind = which.max(start_trim < ticks_interval)-1
  
  # time indices of selected data
  
  data_trim_Tind = as.numeric(dum_trim[start_trim_ind:end_trim_ind,4])
  
  # time indices of 15-sec integer multiple in the target matrix
  
  fill_trim_Tind = seq(as.numeric(dum_trim[start_trim_ind,4]), as.numeric(dum_trim[end_trim_ind,4]), 15)
  
  ### fill in the target matrix according to time indices, leaving unfilled position as NA
  
  # Step 1, fill in direct match values
  
  match_trim_Tind = as.numeric(na.omit(match(data_trim_Tind, fill_trim_Tind)))
	match_trim_target = match_trim_Tind + start_trim_fill_ind - 1
	match_trim_origin = na.omit(match(fill_trim_Tind[match_trim_Tind], data_trim_Tind)) + start_trim_ind - 1
	 
	fts_data[match_trim_target, dup_date_ind[i]] = dum_trim[match_trim_origin, 5]
	
	# Step 2, fill in "nomatch" values observed within 5 seconds of any tick times in the target matrix
	
	if(length(match_trim_origin) < (end_trim_ind - start_trim_ind + 1))
	{
	  # index in the original matrix
	  nomatch_trim_origin_Tind = setdiff(start_trim_ind:end_trim_ind, na.omit(match(fill_trim_Tind[match_trim_Tind], data_trim_Tind)))
	  
	  # index in the target matrix    
	  nomatch_trim_target_Tind = as.numeric(attributes(na.omit(match(data_trim_Tind, fill_trim_Tind)))$na.action)
	  
	  for (j in 1:length(nomatch_trim_origin_Tind))
	  {
	    # check if within 5 seconds of any later tick times
	    
	    nomatch_trim_origin = nomatch_trim_origin_Tind[j]
	    nomatch_trim_target = which(abs(as.numeric(dum_trim[nomatch_trim_origin, 4]) - fill_trim_Tind[head(nomatch_trim_target_Tind, 1):tail(nomatch_trim_target_Tind, 1)]) <= 5)  + head(nomatch_trim_target_Tind, 1) + start_trim_fill_ind - 2
	    
	    # fill in value if a match is found
	    
	    if(length(nomatch_trim_target)> 0)
	    {
	      fts_data[nomatch_trim_target, dup_date_ind[i]] = dum_trim[nomatch_trim_origin, 5]
	    }
	  }
	}
	print(i)
}

write.csv(fts_data, file = "VIX_organised_with_NA.csv", quote=F, row.names=FALSE)

# count NAs in each column
count_na = vector(, ndays)
for(i in 1:ndays)
{
  count_na[i] = sum(is.na(fts_data[,i]))
}


#########################
# fill in missing values
#########################

# use linear interpolation to fill in missing values
fts_nona = data.frame(matrix(, 1621, ncol(fts_data)))
for(j in 1:ncol(fts_data))
{
  fts_nona[,j] = na.interp(fts_data[,j])
  colnames(fts_nona)[j] = noquote(paste("Day", j, sep = "_"))
}

# check if there is still NAs
sum(is.na(fts_nona))

###############
# write output
###############

write.csv(fts_nona, file = "VIX_organised_by_date.csv", quote=F, row.names=FALSE)

##############################
# Convert the price CIDR data 
##############################

trading_time = seq(9.5, 16.25, by = 15/3600)
fts_return = matrix(, 1621, ncol(fts_data))
for(j in 1:ncol(fts_data))
{
  for(i in 2:length(trading_time))
  {
    fts_return[1,j] = 1
    fts_return[i,j] = (log(fts_nona[i,j]) - log(fts_nona[1,j])) * 100
  }
}

# save cumulative return data

write.csv(fts_return, file = "VIX_returns.csv", quote=F, row.names=FALSE)

# make some plots

pdf("fboxplot_VIX_index.pdf",width=12,height=10,pointsize=10)
fboxplot(fts(trading_time, fts_nona), "functional","hdr",c(0.05,0.5),
         legendpos = "topleft", xlab="Fifteen-second time interval", ylab="Daily VIX index", cex=1.5,
         legend = c("30/Jan/2017", "21/Mar/2017", "27/Mar/2017", "17/May/2017", "19/May/2017", "09/Jun/2017", "29/Jun/2017"))
dev.off()

pdf("fboxplot_Cumulative_returns.pdf",width=12,height=10,pointsize=10)
fboxplot(fts(trading_time, fts_return), "functional","hdr",c(0.05,0.5),
         legendpos = "topleft", xlab="Fifteen-second time interval", ylab="Intraday cumulative return", cex=1.5,
         legend = c("30/Jan/2017", "21/Mar/2017", "27/Mar/2017", "17/May/2017", "19/May/2017", "09/Jun/2017", "29/Jun/2017"))
dev.off()


pdf("plot_VIX_index.pdf",width=12,height=10,pointsize=10)
plot(fts(trading_time, fts_nona,xname="Fifteen-second time interval",yname="Daily VIX index"))
dev.off()

pdf("plot_Cumulative_returns.pdf",width=12,height=10,pointsize=10)
plot(fts(trading_time, fts_return,xname="Fifteen-second time interval",yname="Intraday cumulative return"))
dev.off()

###############################################################
# ftsm forecasts via ARIMA and ETS (one-step-ahead experiment)
###############################################################

fts_object = fts(trading_time, fts_return[1:1621,1:(ncol(fts_return)-1)])
ftsm_object_order = head(which(cumsum(ftsm(fts_object)$varprop)>=0.9),1)
arima_forc = forecast(ftsm(fts_object, order = ftsm_object_order), h = 1, method = "arima")

plot(trading_time, arima_forc$mean$y, type = "l", ylim=c(-10,10))
lines(trading_time, fts_return[1:1621,ncol(fts_return)], col=2)



##########################################################################
# Functions for mmeu/mmeo
##########################################################################

mmeu<-function(obs,predicted)
{
  return((1/(length(obs)))*(sum(abs(obs[which(obs>predicted)]-predicted[which(obs>predicted)]))+sum(sqrt(abs(obs[which(obs<predicted)]-predicted[which(obs<predicted)])))))
}
mmeo<-function(obs,predicted)
{
  return((1/(length(obs)))*(sum(abs(obs[which(obs<predicted)]-predicted[which(obs<predicted)]))+sum(sqrt(abs(obs[which(obs>predicted)]-predicted[which(obs>predicted)])))))
}

#####################################
# ratio of change
#####################################

mcpdc<-function(obs,predicted)
{
  return(sum((obs*predicted)>0)/length(obs))
}

































