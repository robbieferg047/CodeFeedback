import matplotlib.pyplot as plt
import scipy.signal as sig
import numpy as np
import pandas as pd
from datetime import datetime, date, timedelta, time
import matplotlib.dates as mdates
import statistics
import json
import scipy
from scipy import stats
import statsmodels.api as sm

#For loop to format dates so arrays can be sliced (doesn't work otherwise)
savepath="C:\\Users\\robbi\\Documents\\GitHub\mazerex2\\cohort2\\fem4_e_PFChm4di_rex1\\"
time_line = np.array(pd.read_csv(savepath+"TimeLine.csv", header=None)).astype(np.datetime64).reshape(-1,1)
win = 30

def timeline(X): #For loop for splicing arrays:
    for date in time_line[X]:
        if date != 2026:
            break
    return(date)


#Habits analyser: Do VRF mice eat at more stereotyped times of day than control?
def gettimes(X, Y, Z, d): #X = cohort, Y = rig, Z = Animal, d = day
    savepath="C:\\Users\\robbi\\Documents\\GitHub\mazerex2\\cohort2\\fem4_e_PFChm4di_rex1\\"
    time_line = np.array(pd.read_csv(savepath+"TimeLine.csv", header=None)).astype(np.datetime64).reshape(-1,1)

    #If statements for determining day from class: (bsl = baseline, ind = induction)

    #Define date from first baseline day (plus or minus)
    if d == -1:
        start_date = time_line[1] + np.timedelta64(d, 'D')
    else:
        if d == 'day1ind':
            start_date = time_line[2]
            
        if d == 'day2ind':
            start_date = time_line[2] + np.timedelta64(1, 'D')
            
        if d == 'day3ind':
            start_date = time_line[2] + np.timedelta64(2, 'D')

    #Standard stuff
    savepath="C:\\Users\\robbi\\Documents\\GitHub\mazerex2\\" + X + "\\" + Y + "\\"
    tags = np.array(pd.read_csv(savepath + "AnimalTags.csv", header=None)).ravel().tolist()
    
    known_tags = [tags[Z]]

    #concatenate data across days to a pandas dataframe
    df = pd.read_csv(savepath + str(start_date)[2:12] + "_events.csv")
    #df['Pellets'] = df['Pellets'] != 0 #Can include reads where pellets were 0 or not

    df['Start_Time']=pd.to_datetime(df['Start_Time'])

    sorted_df = df.sort_values(by=['Start_Time'], ascending=True)
    sorted_df.reset_index(drop=True, inplace=True)
    df=sorted_df.drop([0,1])

    list_x=[]
    list_weights=[]
    list_weightmeds=[]
    an=-1 #plotting index
    for rfid in known_tags: #for loop across animals
        an=an+1
        print(rfid)
        animal_weights=df[df['Animal']==rfid]['Weight'].values
        animal_times=df[df['Animal']==rfid]['Start_Time'].values
        print(animal_times)
        init_prctile=np.percentile(animal_weights[0:win-1],80)
        keep_i=[]
        rolling_i=[]
        for i in range(win-1): #for loop to exclude outliers in first window
            weight=animal_weights[i]
            if weight<1.1*init_prctile and weight>0.9*init_prctile:
                keep_i.append(i)
                rolling_i.append(i)
        rolling_median=np.median(animal_weights[rolling_i])
        print(keep_i)
        rolling_medians=[]
        for i in range(len(keep_i)):
            rolling_medians.append(rolling_median) 
        for j in range(len(animal_weights[win:])): #for loop rolling through data
            weight=animal_weights[j+win]
            if weight<1.2*rolling_median and weight>0.8*rolling_median:
                keep_i.append(j+win)
                rolling_i.pop(0)
                rolling_i.append(j+win)
                rolling_median = np.median(animal_weights[rolling_i])
                rolling_medians.append(rolling_median) 
        print(keep_i)
        x = animal_times[keep_i]
        y = animal_weights[keep_i]
        
        list_x.append(x)
        list_weightmeds.append(rolling_medians)
        list_weights.append(y)
        df = df.reset_index().rename(columns={"Index": "Start_Time", 0: "Counts"})

    #Extracting times
    times = df['Start_Time']

    #Rounding time values to nearest hour
    times = pd.Series(df['Start_Time'])
    print(times)
    times = times.dt.round('1h')
    times = times.dt.time

    #Counting unique values
    times = pd.Series.value_counts(times)
    #Making df readable
    times = times.reset_index().rename(columns={"Index": "Start_Time", 0: "Counts"})
    times.columns = ['Start_Time', 'Counts'] 
    times = pd.DataFrame(times)

    #There isn't really any easy code to include hours where reads were 0, so I wrote it myself:
    #The idea here is that there should be 00:00:00 - 23:00:00 in each df. Where it isn't there, means there
    #was 0 counts, which has to be inserted manually. Can be done automatically by comparing to a df that contains
    #all times 00:00:00 - 23:00:00.

    times24h = (pd.DataFrame({'Time':pd.date_range(start='00:00:00', end='23:00:00', freq='1h')}))['Time'].dt.time #All 24hr times
    timesarranged = (times.sort_values(by=['Start_Time'], ascending=True)).reset_index() #Arranges the time values from the df
    timesarranged = timesarranged['Start_Time']

    #Variable to search through array
    j = -1

    for row in times24h: #The below three lines are here to reformat the df after it is appended with a zero by the loop.
        times24h = (pd.DataFrame({'Time':pd.date_range(start='00:00:00', end='23:00:00', freq='1h')}))['Time'].dt.time
        timesarranged = (times.sort_values(by=['Start_Time'], ascending=True)).reset_index()
        timesarranged = timesarranged['Start_Time']

        j = j + 1 
            
        try:
            if timesarranged[j] != times24h[j]: #j increases each loop, so each line is successively checked for the next hour value
                zeroes = {
                "Start_Time":[times24h[j]], 
                "Counts":[0]}
                zeroes = pd.DataFrame(zeroes)
                times = pd.concat([times, zeroes]) #If it isn't there, we concat that specific time, and add a 'Count' of 0.
        except: #This is here because the above won't work if 23:00:00 is missing.
            j == 23
            zeroes = {
            "Start_Time":[times24h[j]], 
            "Counts":[0]}
            zeroes = pd.DataFrame(zeroes)
            times = pd.concat([times, zeroes])
            times = (times.sort_values(by=['Start_Time'], ascending=True)).reset_index()
            times['Start_Time'] = times['Start_Time'].astype(str) #Formatting for graphs
            print(times)                    
            return times
        else:
            if j == 23: #This is here because the loop might hit 23, and thus be out of range for the arranged times. In that case, this script just forces the concat
                print(Y, d, times) #So I can see if the script isn't correcting something it should be
                times = (times.sort_values(by=['Start_Time'], ascending=True)).reset_index()
                times['Start_Time'] = times['Start_Time'].astype(str) #Formatting for graphs
                print(times)
                return times
            else:
                continue

#The below is all just running the commands and plotting.
PFCtimesday1ind = gettimes("cohort2", "fem4_e_PFChm4di_rex1", 1, 'day1ind')

controlcol = "#00000092"
VRFcol = "#FF00008D"
PFCcol = "#0040FF8D"
ICcol = "#00FF488D"

#####################################

fig1, ax1 = plt.subplots(figsize=(16, 16))

x1, y1 = PFCtimesday1ind['Start_Time'], PFCtimesday1ind['Counts']
ax1.plot(x1, y1, marker='o', linestyle = '-', color=(PFCcol), label = 'PFC')
ax1.legend()
ax1.axvline('12:00:00', color='black', linestyle='dashed', label = 'midday')
plt.title("Day 1 Induction, PFC Animal 1")
plt.xticks(rotation=45, ha='right')
plt.xlabel("Time")
plt.ylabel("Counts")
plt.show()
