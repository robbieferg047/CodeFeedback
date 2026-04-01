from datetime import datetime, date, timedelta, time
import schedule
import csv
import matplotlib.pyplot as plt
import scipy.signal as sig
import numpy as np
import pandas as pd
from datetime import datetime, date, timedelta, time
import matplotlib.dates as mdates
import statistics
import json

#Getting savepath and the day:
savepath = "C:\\Users\\robbi\\Documents\\GitHub\mazerex2\\cohort3\\fem9_c_rex2\\"
today=np.datetime64(datetime.today())
today=np.datetime64(date(2026,3,20)) #Date where we had an error file

#Telling script what to read:
data_lindz = pd.read_csv(savepath + str(today)[0:12] + "_autolindz.csv")
frames=[data_lindz]
data_lindz=pd.concat(frames)
data_lindz['Start_Time'] = pd.to_datetime(data_lindz['Start_Time'], format ='mixed')

#Extracting times
savepat = "C:\\Users\\robbi\\Documents\\GitHub\mazerex2\\cohort2\\fem7_c_VRF_rex4\\"
t=data_lindz['Start_Time'].values

#Defining class to be run:
def check():
    if len(data_lindz) >= 1000:
        times = np.diff(t)
        counts = len(times <= 5000) #Counts how many unrealistic read times there are (This number may need work)

        sortedtimes = sorted(times) #Arranges them from biggest to smallest
        tenpercent = int(0.1*len(sortedtimes)) #Extracts the top 10% (If it's in error mode the top 10% (0.1) of read times should be unrealistic)

        sortedtimes = pd.DataFrame(sortedtimes) #Converting to df

        threshold = sortedtimes.nlargest(tenpercent, 0) #Extracting the top 10% largest numbers
        threshold = np.mean(threshold) #The threshold for an unrealistic read is the average
        
        j = -1 #Index for cutting

        if counts > 1000: 
            for timediffs in times:
                if timediffs > threshold:
                    j = j+1 #j will be used to splice the array
                else:
                    break
                print(j)
                
            newfile = data_lindz[j:len(data_lindz)] #Chops the csv file according to how many unrealistic reads until the end of it
            newfile.to_csv(savepat + "test7" + "_autolindz.csv")    
                
check()

#Scheduling check() to run every few hours:
#schedule.every(1).hours.do(check) #Commented for now so it doesn't do that.