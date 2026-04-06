import matplotlib.pyplot as plt
import scipy.signal as sig
import numpy as np
import pandas as pd
from datetime import datetime, date, timedelta, time
import matplotlib.dates as mdates
import statistics
import json

#"fem7_c_VRF_rex4"
#"fem4_e_PFChm4di_rex1"
#"fem6_e_IChm4di_rex3"
#"fem5_c_rex2"

#cohort2
#cohort3

savepath="C:\\Users\\robbi\\Documents\\GitHub\mazerex2\\cohort2\\fem4_e_PFChm4di_rex1\\"
win=30 #size of rolling median window
bins=58
time_line = np.array(pd.read_csv(savepath+"TimeLine.csv", header=None)).astype(np.datetime64).reshape(-1,1)
start_date=time_line[0]

marker_times=time_line#add important dates here to add vertical lines on last plot
#last_date=np.datetime64(datetime.today()) #OR TYPE DESIRED DATE ON NEXT LINE AND UNCOMMENT IT

last_date=np.datetime64(date(2026,3,2))
datetag=str(pd.to_datetime(last_date).to_pydatetime().date())

#For loop to format dates so arrays can be sliced (doesn't work otherwise)
def timeline(X):
    for date in time_line[X]:
        if date != 2026:
            break
    return(date)

#Class for getting rolling_medians:
def Rolling_Medians(X, Y, Z): #X Defines the cohort, Y defines the rig, Z defines the animal (I.e., Rolling_Medians("cohort 2", "fem4_e_PFChm4di_rex1", 1)) is cohort 2, PFC, first animal
    savepath="C:\\Users\\robbi\\Documents\\GitHub\mazerex2\\" + X + "\\" + Y + "\\"
    time_line = np.array(pd.read_csv(savepath+"TimeLine.csv", header=None)).astype(np.datetime64).reshape(-1,1)
    known_tags = np.array(pd.read_csv(savepath + "AnimalTags.csv", header=None)).ravel().tolist()[Z]
    

    filtermin=14 #lower limit in g
    filtermax=24 #upper limit in g 
    d=(last_date-start_date)/np.timedelta64(1,'D')

    days_to_plot=round(float(d[0]))

    #concatenate data across days to a long pandas dataframe
    data_coll_weight = pd.read_csv(savepath + str(start_date)[2:12] + "_events.csv")
    for j in range(days_to_plot):
            day=start_date+np.timedelta64(j+1,'D') 
            data = pd.read_csv(savepath + str(day)[2:12] + "_events.csv") 
            frames=[data_coll_weight,data]
            data_coll_weight=pd.concat(frames)
    df=data_coll_weight
    df['Start_Time']=pd.to_datetime(df['Start_Time'])
    sorted_df = df.sort_values(by=['Start_Time'], ascending=True)
    del df

    sorted_df.reset_index(drop=True, inplace=True)
    df=sorted_df.drop([0,1])
   
    #Extracting relevant days
    df = df.loc[df['Animal'] == known_tags]
    print(df)

    df = df.loc[df['Start_Time'] <= timeline(6)]
    df = df.loc[df['Start_Time'] >= timeline(1)]

    # create a rolling median filter, roll through each animal's data, and plot
    list_x=[]
    list_weights=[]
    list_weightmeds=[] 
    #Extracting baseline
    dfbl = df.loc[df['Start_Time'] <= timeline(2)]

    #Extracting average weight across the baseline period
    baselineweight = dfbl['Weight'].mean()

    animal_weights = df['Weight'].values
    animal_times = df['Start_Time'].values

    list_daily_w=[]

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
    
    x = animal_times[keep_i]
    y = animal_weights[keep_i]

    y = (((y) / (baselineweight)))*100

    list_x.append(x)
    list_weightmeds.append(rolling_medians)
    list_weights.append(y)

    data={
         "Date":list_x[0],
         "Weight":list_weights[0]
        }

    filtered_df=pd.DataFrame(data)
    daily_avg_w=filtered_df.groupby(pd.Grouper(key='Date', freq='12h', origin=str(start_date)[2:12])).median().reset_index()

    y = daily_avg_w['Weight']
    x = daily_avg_w['Date']

    return x, y

#Grabbing PFC data for one animal
an1 = Rolling_Medians("cohort2", "fem4_e_PFChm4di_rex1", 0)
#an2 = Rolling_Medians("cohort2", "fem4_e_PFChm4di_rex1", 1)
#an3 = Rolling_Medians("cohort2", "fem4_e_PFChm4di_rex1", 2)
#an4 = Rolling_Medians("cohort2", "fem4_e_PFChm4di_rex1", 3)
#an5 = Rolling_Medians("cohort2", "fem4_e_PFChm4di_rex1", 4)

#Plotting
fig1, ax1 = plt.subplots(figsize=(16, 16))
ax1.plot(an1[0], an1[1], marker='s', linestyle = '-', alpha=0.9, color='lightcoral', linewidth=2)
#ax1.plot(an2[0], an2[1], marker='s', linestyle = '-', alpha=0.9, color='darkred', linewidth=2)
#ax1.plot(an3[0], an3[1], marker='s', linestyle = '-', alpha=0.9, color='salmon', linewidth=2)
#ax1.plot(an4[0], an4[1], marker='s', linestyle = '-', alpha=0.9, color='indianred', linewidth=2)
#ax1.plot(an5[0], an5[1], marker='s', linestyle = '-', alpha=0.9, color='crimson', linewidth=2)

ax1.axvline(timeline(1), color='black', linestyle='dashed', label = "Baseline")
ax1.axvline(timeline(2), color='blue', linestyle='dashed', label = "Induction 1")
ax1.axvline(timeline(3), color='blue', linestyle='dashed')
ax1.axvline(timeline(4), color='red', linestyle='dashed', label = "Induction 2")
ax1.axvline(timeline(5), color='red', linestyle='dashed')
ax1.axvline(timeline(6), color='gray', linestyle='dashed')

ax1.set_title("mPFC DREADDs" )
ax1.legend()
ax1.set_ylabel("% Body Weight")
ax1.set_xlabel("Time")
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%m-%d %H:%M:%S'))
ax1.set_xticks(ax1.get_xticks(), ax1.get_xticklabels(), rotation=45, ha='right')
ax1.grid(True)
plt.show()
