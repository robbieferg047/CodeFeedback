import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from datetime import datetime, date, timedelta, time
import matplotlib.dates as mdates

def GetData(cohort, rig): #I.e., Rolling_Medians("cohort 2", "fem4_e_PFChm4di_rex1", 1 is cohort 2, PFC
    print("Processing Data") #The script may take a bit of time to run, so added print commands throughout to show the status
    win = 30
    savepath="C:\\Users\\robbi\\Documents\\GitHub\mazerex2\\" + cohort + "\\" + rig + "\\" #standard stuff

    if rig == 'None':
        savepath="C:\\Users\\robbi\\Documents\\GitHub\mazerex2\\" + cohort + "\\" #For mal2_RI

    known_tags = np.array(pd.read_csv(savepath + "AnimalTags.csv", header=None)).ravel().tolist() 
    time_line = np.array(pd.read_csv(savepath + "TimeLine.csv", header=None)).astype(np.datetime64).reshape(-1,1) 

    def timeline(X):  #for loop to format dates so arrays can be sliced (doesn't work otherwise)
        for date in time_line[X]:
            if date != 2026:
                break
        return(date)
    
    start_date = (time_line[1] - np.timedelta64(24, 'h')) #Starting for 72h before induction 1
    last_date = (time_line[1] + np.timedelta64(8, 'D')) #Going to end of induction 1

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
    df['Animal']=df['Animal'].astype(int)
    sorted_df = df.sort_values(by=['Start_Time'], ascending=True)
    sorted_df = sorted_df.loc[sorted_df['Start_Time'] >= (timeline(1) - np.timedelta64(24, 'h'))] #Cutting off to 72h before induction 1

    sorted_df['Start_Time'] = ((sorted_df['Start_Time']) - (timeline(1) + np.timedelta64(2, 'D'))) #Making all timepoints relative to induction 1
    sorted_df['Start_Time'] = sorted_df['Start_Time']/np.timedelta64(1, 's') #Converting to seconds

    if rig == "None":
        for animal in known_tags:
            if animal == 0:
                continue
            else:
                sorted_df_an = sorted_df.loc[sorted_df['Animal'] == animal]
                sorted_df_an.to_csv("C:\\Users\\robbi\\Desktop\\Masters\\VRF Project\\OSF_repo\\" + cohort + "\\" + f'{str(animal)}.csv')

    else:
        for animal in known_tags:
            if animal == 0:
                continue
            else:
                sorted_df_an = sorted_df.loc[sorted_df['Animal'] == animal]
                sorted_df_an.to_csv("C:\\Users\\robbi\\Desktop\\Masters\\VRF Project\\OSF_repo\\" + cohort + "\\" + rig + "\\" + f'{str(animal)}.csv')

    return sorted_df

cohorts = ["cohort2", "cohort3",  "mal2_RI"]
rigs = ["fem7_c_VRF_rex4", "fem4_e_PFChm4di_rex1", "fem6_e_IChm4di_rex3", "fem5_c_rex2", 
     "fem8_c_VRF_rex1", "fem9_c_rex2", "fem10_e_PFChm4di_rex3", "fem11_e_IChm4di_rex4", "None"]

for cohort in cohorts:
    for rig in rigs:
        try:
            GetData(cohort, rig)
            print("Running Unique Combination!")
        except Exception:
            print("Null combo! Passing...")
            pass

