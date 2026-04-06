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

    if d == 'day-4bsl':
        start_date = time_line[1] - np.timedelta64(4, 'D')

    if d == 'day-3bsl':
        start_date = time_line[1] - np.timedelta64(3, 'D') 

    if d == 'day-2bsl':
        start_date = time_line[1] - np.timedelta64(2, 'D') 
        
    if d == 'day-1bsl':
        start_date = time_line[1] - np.timedelta64(1, 'D')

    if d == 'day1bsl':
        start_date = time_line[1] 
        
    if d == 'day2bsl':
        start_date = time_line[1] + np.timedelta64(1, 'D')

    if d == 'day1ind':
        start_date = time_line[2]

    if d == 'day2ind':
        start_date = time_line[2] + np.timedelta64(1, 'D')

    if d == 'day3ind':
        start_date = time_line[2] + np.timedelta64(2, 'D')


    #Standard stuff
    savepath="C:\\Users\\robbi\\Documents\\GitHub\mazerex2\\" + X + "\\" + Y + "\\"
    tags = np.array(pd.read_csv(savepath + "AnimalTags.csv", header=None)).ravel().tolist()

    for tags_ in tags:
        if Z == "all":
            known_tags = tags
            break
        else:
            known_tags = [tags[Z]]


    #concatenate data across days to a pandas dataframe
    df = pd.read_csv(savepath + str(start_date)[2:12] + "_events.csv")

    df['Start_Time']=pd.to_datetime(df['Start_Time'])

    sorted_df = df.sort_values(by=['Start_Time'], ascending=True)
    sorted_df.reset_index(drop=True, inplace=True)
    df=sorted_df.drop([0,1])
    
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
                print(Y, d, 
                      times) #So I can see if the script isn't correcting something it should be
                times = (times.sort_values(by=['Start_Time'], ascending=True)).reset_index()
                times['Start_Time'] = times['Start_Time'].astype(str) #Formatting for graphs
                print(times)
                return times
            else:
                continue

#The below is all just running the commands and plotting.

controltimesdaym4bsl = gettimes("cohort2", "fem5_c_rex2", "all", "day-4bsl")
VRFtimesdaym4bsl = gettimes("cohort2", "fem7_c_VRF_rex4", "all", "day-4bsl")
PFCtimesdaym4bsl = gettimes("cohort2", "fem4_e_PFChm4di_rex1", "all", "day-4bsl")
ICtimesdaym4bsl = gettimes("cohort2", "fem6_e_IChm4di_rex3", "all", "day-4bsl")

controltimesdaym3bsl = gettimes("cohort2", "fem5_c_rex2", "all", "day-3bsl")
VRFtimesdaym3bsl = gettimes("cohort2", "fem7_c_VRF_rex4", "all", "day-3bsl")
PFCtimesdaym3bsl = gettimes("cohort2", "fem4_e_PFChm4di_rex1", "all", "day-3bsl")
ICtimesdaym3bsl = gettimes("cohort2", "fem6_e_IChm4di_rex3", "all", "day-3bsl")

controltimesdaym2bsl = gettimes("cohort2", "fem5_c_rex2", "all", "day-2bsl")
VRFtimesdaym2bsl = gettimes("cohort2", "fem7_c_VRF_rex4", "all", "day-2bsl")
PFCtimesdaym2bsl = gettimes("cohort2", "fem4_e_PFChm4di_rex1", "all", "day-2bsl")
ICtimesdaym2bsl = gettimes("cohort2", "fem6_e_IChm4di_rex3", "all", "day-2bsl")

controltimesdaym1bsl = gettimes("cohort2", "fem5_c_rex2", "all", "day-1bsl")
VRFtimesdaym1bsl = gettimes("cohort2", "fem7_c_VRF_rex4", "all", "day-1bsl")
PFCtimesdaym1bsl = gettimes("cohort2", "fem4_e_PFChm4di_rex1", "all", "day-1bsl")
ICtimesdaym1bsl = gettimes("cohort2", "fem6_e_IChm4di_rex3", "all", "day-1bsl")

controltimesday1bsl = gettimes("cohort2", "fem5_c_rex2", "all", "day1bsl")
VRFtimesday1bsl = gettimes("cohort2", "fem7_c_VRF_rex4", "all", "day1bsl")
PFCtimesday1bsl = gettimes("cohort2", "fem4_e_PFChm4di_rex1", "all", "day1bsl")
ICtimesday1bsl = gettimes("cohort2", "fem6_e_IChm4di_rex3", "all", "day1bsl")

controltimesday2bsl = gettimes("cohort2", "fem5_c_rex2", "all", "day2bsl")
VRFtimesday2bsl = gettimes("cohort2", "fem7_c_VRF_rex4", "all", "day2bsl")
PFCtimesday2bsl = gettimes("cohort2", "fem4_e_PFChm4di_rex1", "all", "day2bsl")
ICtimesday2bsl = gettimes("cohort2", "fem6_e_IChm4di_rex3", "all", "day2bsl")

controltimesday1ind = gettimes("cohort2", "fem5_c_rex2", "all", "day1ind")
VRFtimesday1ind = gettimes("cohort2", "fem7_c_VRF_rex4", "all", "day1ind")
PFCtimesday1ind = gettimes("cohort2", "fem4_e_PFChm4di_rex1", "all", "day1ind")
ICtimesday1ind = gettimes("cohort2", "fem6_e_IChm4di_rex3", "all", "day1ind")

controltimesday2ind = gettimes("cohort2", "fem5_c_rex2", "all", "day2ind")
VRFtimesday2ind = gettimes("cohort2", "fem7_c_VRF_rex4", "all", "day2ind")
PFCtimesday2ind = gettimes("cohort2", "fem4_e_PFChm4di_rex1", "all", "day2ind")
ICtimesday2ind = gettimes("cohort2", "fem6_e_IChm4di_rex3", "all", "day2ind")

controltimesday3ind = gettimes("cohort2", "fem5_c_rex2", "all", "day3ind")
VRFtimesday3ind = gettimes("cohort2", "fem7_c_VRF_rex4", "all", "day3ind")
PFCtimesday3ind = gettimes("cohort2", "fem4_e_PFChm4di_rex1", "all", "day3ind")
ICtimesday3ind = gettimes("cohort2", "fem6_e_IChm4di_rex3", "all", "day3ind")


controlcol = "#00000092"
VRFcol = "#FF00008D"
PFCcol = "#0040FF8D"
ICcol = "#00FF488D"

#####################################

fig1 = plt.figure(figsize=(16, 16))
ax1 = fig1.add_subplot(221)

x1, y1 = controltimesdaym1bsl['Start_Time'], controltimesdaym1bsl['Counts']
x2, y2 = VRFtimesdaym1bsl['Start_Time'], VRFtimesdaym1bsl['Counts']
#x3, y3 = PFCtimesdaym1bsl['Start_Time'], PFCtimesdaym1bsl['Counts']
#x4, y4 = ICtimesdaym1bsl['Start_Time'], ICtimesdaym1bsl['Counts']

ax1.plot(x1, y1, marker='o', linestyle = '-', color=(controlcol), label = 'Control')
ax1.plot(x2, y2, marker='o', linestyle = '-', color=(VRFcol), label = 'VRF')
#ax1.plot(x3, y3, marker='o', linestyle = '-', color=(PFCcol), label = 'PFC')
#ax1.plot(x4, y4, marker='o', linestyle = '-', color=(ICcol), label = 'IC')
ax1.legend()
ax1.axvline('12:00:00', color='black', linestyle='dashed', label = 'midday')
ax1.axes.get_xaxis().set_visible(False)
plt.title("Day -1 baseline")
plt.ylabel("Counts")

ax2 = fig1.add_subplot(222)

x1, y1 = controltimesdaym2bsl['Start_Time'], controltimesdaym2bsl['Counts']
x2, y2 = VRFtimesdaym2bsl['Start_Time'], VRFtimesdaym2bsl['Counts']
#x3, y3 = PFCtimesdaym2bsl['Start_Time'], PFCtimesdaym2bsl['Counts']
#x4, y4 = ICtimesdaym2bsl['Start_Time'], ICtimesdaym2bsl['Counts']

ax2.plot(x1, y1, marker='o', linestyle = '-', color=(controlcol), label = 'Control')
ax2.plot(x2, y2, marker='o', linestyle = '-', color=(VRFcol), label = 'VRF')
#ax2.plot(x3, y3, marker='o', linestyle = '-', color=(PFCcol), label = 'PFC')
#ax2.plot(x4, y4, marker='o', linestyle = '-', color=(ICcol), label = 'IC')
ax2.legend()
ax2.axvline('12:00:00', color='black', linestyle='dashed', label = 'midday')
ax2.axes.get_xaxis().set_visible(False)
plt.title("Day -2 baseline")
plt.ylabel("Counts")

ax3 = fig1.add_subplot(223)

x1, y1 = controltimesdaym3bsl['Start_Time'], controltimesdaym3bsl['Counts']
x2, y2 = VRFtimesdaym3bsl['Start_Time'], VRFtimesdaym3bsl['Counts']
#x3, y3 = PFCtimesdaym3bsl['Start_Time'], PFCtimesdaym3bsl['Counts']
#x4, y4 = ICtimesdaym3bsl['Start_Time'], ICtimesdaym3bsl['Counts']

ax3.plot(x1, y1, marker='o', linestyle = '-', color=(controlcol), label = 'Control')
ax3.plot(x2, y2, marker='o', linestyle = '-', color=(VRFcol), label = 'VRF')
#ax3.plot(x3, y3, marker='o', linestyle = '-', color=(PFCcol), label = 'PFC')
#ax3.plot(x4, y4, marker='o', linestyle = '-', color=(ICcol), label = 'IC')
ax3.legend()
ax3.axvline('12:00:00', color='black', linestyle='dashed', label = 'midday')
plt.title("Day -3 baseline")
plt.xticks(rotation=45, ha='right')
plt.xlabel("Time")
plt.ylabel("Counts")

ax4 = fig1.add_subplot(224)

x1, y1 = controltimesdaym4bsl['Start_Time'], controltimesdaym4bsl['Counts']
x2, y2 = VRFtimesdaym4bsl['Start_Time'], VRFtimesdaym4bsl['Counts']
#x3, y3 = PFCtimesdaym4bsl['Start_Time'], PFCtimesdaym4bsl['Counts']
#x4, y4 = ICtimesdaym4bsl['Start_Time'], ICtimesdaym4bsl['Counts']

ax4.plot(x1, y1, marker='o', linestyle = '-', color=(controlcol), label = 'Control')
ax4.plot(x2, y2, marker='o', linestyle = '-', color=(VRFcol), label = 'VRF')
#ax4.plot(x3, y3, marker='o', linestyle = '-', color=(PFCcol), label = 'PFC')
#ax4.plot(x4, y4, marker='o', linestyle = '-', color=(ICcol), label = 'IC')
ax4.legend()
ax4.axvline('12:00:00', color='black', linestyle='dashed', label = 'midday')
plt.title("Day -4 baseline")
plt.xticks(rotation=45, ha='right')

plt.xlabel("Time")
plt.ylabel("Counts")
plt.show()

##################################################################

fig2 = plt.figure(figsize=(16, 16))
ax1 = fig2.add_subplot(221)

x1, y1 = controltimesday1bsl['Start_Time'], controltimesday1bsl['Counts']
x2, y2 = VRFtimesday1bsl['Start_Time'], VRFtimesday1bsl['Counts']
#x3, y3 = PFCtimesday1bsl['Start_Time'], PFCtimesday1bsl['Counts']
#x4, y4 = ICtimesday1bsl['Start_Time'], ICtimesday1bsl['Counts']

ax1.plot(x1, y1, marker='o', linestyle = '-', color=(controlcol), label = 'Control')
ax1.plot(x2, y2, marker='o', linestyle = '-', color=(VRFcol), label = 'VRF')
#ax1.plot(x3, y3, marker='o', linestyle = '-', color=(PFCcol), label = 'PFC')
#ax1.plot(x4, y4, marker='o', linestyle = '-', color=(ICcol), label = 'IC')
ax1.legend()
ax1.axvline('12:00:00', color='black', linestyle='dashed', label = 'midday')
ax1.axes.get_xaxis().set_visible(False)
plt.title("Day 1 Baseline")
plt.ylabel("Counts")

ax2 = fig2.add_subplot(222)

x1, y1 = controltimesday2bsl['Start_Time'], controltimesday2bsl['Counts']
x2, y2 = VRFtimesday2bsl['Start_Time'], VRFtimesday2bsl['Counts']
#x3, y3 = PFCtimesday2bsl['Start_Time'], PFCtimesday2bsl['Counts']
#x4, y4 = ICtimesday2bsl['Start_Time'], ICtimesday2bsl['Counts']

ax2.plot(x1, y1, marker='o', linestyle = '-', color=(controlcol), label = 'Control')
ax2.plot(x2, y2, marker='o', linestyle = '-', color=(VRFcol), label = 'VRF')
#ax2.plot(x3, y3, marker='o', linestyle = '-', color=(PFCcol), label = 'PFC')
#ax2.plot(x4, y4, marker='o', linestyle = '-', color=(ICcol), label = 'IC')
ax2.legend()
ax2.axvline('12:00:00', color='black', linestyle='dashed', label = 'midday')
ax2.axes.get_xaxis().set_visible(False)
plt.title("Day 2 Baseline")
plt.ylabel("Counts")

ax3 = fig2.add_subplot(223)

x1, y1 = controltimesday1ind['Start_Time'], controltimesday1ind['Counts']
x2, y2 = VRFtimesday1ind['Start_Time'], VRFtimesday1ind['Counts']
#x3, y3 = PFCtimesday1ind['Start_Time'], PFCtimesday1ind['Counts']
#x4, y4 = ICtimesday1ind['Start_Time'], ICtimesday1ind['Counts']

ax3.plot(x1, y1, marker='o', linestyle = '-', color=(controlcol), label = 'Control')
ax3.plot(x2, y2, marker='o', linestyle = '-', color=(VRFcol), label = 'VRF')
#ax3.plot(x3, y3, marker='o', linestyle = '-', color=(PFCcol), label = 'PFC')
#ax3.plot(x4, y4, marker='o', linestyle = '-', color=(ICcol), label = 'IC')
ax3.legend()
ax3.axvline('12:00:00', color='black', linestyle='dashed', label = 'midday')
plt.title("Day 1 Induction")
plt.xticks(rotation=45, ha='right')
plt.xlabel("Time")
plt.ylabel("Counts")

ax4 = fig2.add_subplot(224)

x1, y1 = controltimesday2ind['Start_Time'], controltimesday2ind['Counts']
x2, y2 = VRFtimesday2ind['Start_Time'], VRFtimesday2ind['Counts']
#x3, y3 = PFCtimesday2ind['Start_Time'], PFCtimesday2ind['Counts']
#x4, y4 = ICtimesday2ind['Start_Time'], ICtimesday2ind['Counts']

ax4.plot(x1, y1, marker='o', linestyle = '-', color=(controlcol), label = 'Control')
ax4.plot(x2, y2, marker='o', linestyle = '-', color=(VRFcol), label = 'VRF')
#ax4.plot(x3, y3, marker='o', linestyle = '-', color=(PFCcol), label = 'PFC')
#ax4.plot(x4, y4, marker='o', linestyle = '-', color=(ICcol), label = 'IC')
ax4.legend()
ax4.axvline('12:00:00', color='black', linestyle='dashed', label = 'midday')
plt.title("Day 2 Induction")
plt.xticks(rotation=45, ha='right')

plt.xlabel("Time")
plt.ylabel("Counts")
plt.show()

###############################################
fig2 = plt.figure(figsize=(16, 16))
ax1 = fig2.add_subplot(221)

x1, y1 = controltimesday3ind['Start_Time'], controltimesday3ind['Counts']
x2, y2 = VRFtimesday3ind['Start_Time'], VRFtimesday3ind['Counts']
#x3, y3 = PFCtimesday3ind['Start_Time'], PFCtimesday3ind['Counts']
#x4, y4 = ICtimesday3ind['Start_Time'], ICtimesday3ind['Counts']

ax1.plot(x1, y1, marker='o', linestyle = '-', color=(controlcol), label = 'Control')
ax1.plot(x2, y2, marker='o', linestyle = '-', color=(VRFcol), label = 'VRF')
#ax1.plot(x3, y3, marker='o', linestyle = '-', color=(PFCcol), label = 'PFC')
#ax1.plot(x4, y4, marker='o', linestyle = '-', color=(ICcol), label = 'IC')
ax1.legend()
ax1.axvline('12:00:00', color='black', linestyle='dashed', label = 'midday')
plt.title("Day 3 Induction")
plt.xticks(rotation=45, ha='right')

plt.xlabel("Time")
plt.ylabel("Counts")
plt.show()
