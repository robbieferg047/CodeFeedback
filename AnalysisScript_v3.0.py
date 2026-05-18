import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from datetime import datetime, date, timedelta, time
import matplotlib.dates as mdates
import scipy.stats as sp
import statsmodels.api as sm
import statsmodels.formula.api as smf
import statistics as stats
import re
import itertools
import seaborn as sb
import statannotations.Annotator as sn

fig1 = plt.figure(figsize=(16, 8)) #Defining figs
ax1 = fig1.add_subplot(221)
ax2 = fig1.add_subplot(222)
ax3 = fig1.add_subplot(223)
ax4 = fig1.add_subplot(224)
fig2 = plt.figure(figsize=(16, 8))
ax5 = fig2.add_subplot(221)
ax6 = fig2.add_subplot(222)
ax7 = fig2.add_subplot(223)
fig3 = plt.figure(figsize=(16, 8))
ax8 = fig3.add_subplot(221)
ax9 = fig3.add_subplot(222)
ax10 = fig3.add_subplot(223)
ax11 = fig3.add_subplot(224)
fig4 = plt.figure(figsize=(16, 8))
ax12 = fig4.add_subplot(221)
ax13 = fig4.add_subplot(222)
ax14 = fig4.add_subplot(223)
ax15 = fig4.add_subplot(224)
fig5 = plt.figure(figsize=(16, 8))
ax16 = fig5.add_subplot(221)
ax17 = fig5.add_subplot(222)
ax18 = fig5.add_subplot(223)
ax19 = fig5.add_subplot(224)
fig6 = plt.figure(figsize=(16, 8))
ax20 = fig6.add_subplot(221)
ax21 = fig6.add_subplot(222)
ax22 = fig6.add_subplot(223)

fig_axes = [ax8, ax9, ax10, ax11, ax12, ax13, ax14, 
    ax15, ax16, ax17, ax18, ax19, ax20, ax21, ax22] #For habits analysis

cohorts = ["cohort2", "cohort3"]
rigs = ["fem7_c_VRF_rex4", "fem4_e_PFChm4di_rex1", "fem6_e_IChm4di_rex3", "fem5_c_rex2", 
     "fem8_c_VRF_rex1", "fem9_c_rex2", "fem10_e_PFChm4di_rex3", "fem11_e_IChm4di_rex4"] #Here, we can list off our cohorts and rigs for ease

bin_size = 12 #For unilateral control of bins variable (Supports 12 and 24 hours)
days_for_plots = bin_size/24

def GetData(cohort, rig): #I.e., Rolling_Medians("cohort 2", "fem4_e_PFChm4di_rex1", 1 is cohort 2, PFC
    print("Processing Data") #The script may take a bit of time to run, so added print commands throughout to show the status
    win = 30
    savepath="C:\\Users\\robbi\\Documents\\GitHub\mazerex2\\" + cohort + "\\" + rig + "\\" #standard stuff
    known_tags = np.array(pd.read_csv(savepath + "AnimalTags.csv", header=None)).ravel().tolist() 
    time_line = np.array(pd.read_csv(savepath+"TimeLine.csv", header=None)).astype(np.datetime64).reshape(-1,1) 

    def timeline(X):  #for loop to format dates so arrays can be sliced (doesn't work otherwise)
        for date in time_line[X]:
            if date != 2026:
                break
        return(date)
    
    start_date = (time_line[1])
    last_date = (time_line[1] + np.timedelta64(14, 'D'))

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

    sorted_df.reset_index(drop=True, inplace=True)
    df=sorted_df.drop([0,1]) #workaround!

    #create a rolling median filter, roll through each animal's data, and plot
    list_x=[]
    list_weights=[]
    list_weightmeds=[]
    an=-1 #plotting index
    for rfid in known_tags: #for loop across animals
        an=an+1
        animal_weights=df[df['Animal']==rfid]['Weight'].values
        animal_times=df[df['Animal']==rfid]['Start_Time'].values
        init_prctile=np.percentile(animal_weights[0:win-1],80)
        keep_i=[]
        rolling_i=[]
        for i in range(win-1): #for loop to exclude outliers in first window
            weight=animal_weights[i]
            if weight<1.1*init_prctile and weight>0.9*init_prctile:
                keep_i.append(i)
                rolling_i.append(i)
        rolling_median=np.median(animal_weights[rolling_i])
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
        list_x.append(x)
        list_weights.append(y)

    an=-1
    df1=df[df['Pellets'] != 0]
    list_m = []
    list_t = []
    for rfid in known_tags: #for loop across animals
        an = an + 1
        p=df1[df1['Animal']==rfid]['Pellets'].values
        t=df1[df1['Animal']==rfid]['Start_Time'].values
        u=df1[df1['Animal']==rfid]['Unit'].values

        intervals=np.diff(t)/1000000000
        y, x = np.histogram(intervals, bins=np.arange(0,120,2))
        
        meal_threshold=30 #define duration (s) of meal based on histograms 
        rows=len(p)-1
        for i in range(rows):
            row=rows-i
            interval=intervals[row-1]
            if u[row]==u[row-1]: #same unit
                if interval<meal_threshold:
                    p[row-1]=p[row-1]+p[row]
                    p=np.delete(p, row)
                    t=np.delete(t, row)
                    u=np.delete(u, row)
        list_m.append(p)
        list_t.append(t)
                   
    allanimals = []
    an=-1

    for rfid in known_tags: #for loop across animals
        an=an+1
        data={
            "Date":list_x[an],
            "Weight":list_weights[an],
            "Animal":int(known_tags[an])
            }
        data_m = {
            "MealSize":list_m[an],
            "MealDates":list_t[an]
            }   
        
        filtered_df=pd.DataFrame(data)
        filtered_df_m=pd.DataFrame(data_m)

        #Rhythmicity analyser:
        mealtimes = filtered_df_m['MealDates'] #Extract mealdates
        mealtimes = mealtimes.dt.round('1h') #Extract just the hour value

        #Counting unique values
        times = pd.Series.value_counts(mealtimes)
        #Making df readable
        times = times.reset_index().rename(columns={"Index": "MealDates", 0: "Counts"})
        times.columns = ['MealDates', 'Counts'] 
        times = (times.sort_values(by=['MealDates'], ascending=True)).reset_index() 
        times = pd.DataFrame(times)

        dailyp_list = []
        dailyp_date = []
        dailyp_animal = []

        count_avglist = []
        datelist = []
        hourlist = []
        hourdf = []

        for days in pd.date_range(start=timeline(1) - np.timedelta64(12, 'h'), 
                                  end=timeline(1) + np.timedelta64(348, 'h'), freq = '12h'):
                filtered_df_m['MealDates'] = filtered_df_m['MealDates'].dt.round('12h')
                pellets = filtered_df_m.loc[filtered_df_m['MealDates'] == days]
                pellets = np.sum(pellets['MealSize'])
                pellets_an = known_tags[an]
                dailyp_list.append(pellets)
                dailyp_date.append(days)
                dailyp_animal.append(pellets_an)
        
        data_p = {
                "TotalPellets":dailyp_list,
                "PelletDates":dailyp_date
                }   
        
        df_p=pd.DataFrame(data_p)

        df = filtered_df.groupby(pd.Grouper(key='Date', freq = f'{bin_size}' + "h", origin=str(time_line[1])[2:12])).median().reset_index()
        df_m = filtered_df_m.groupby(pd.Grouper(key='MealDates', freq = f'{bin_size}' + "h", origin=str(time_line[1])[2:12])).median().reset_index()

        baselineweight = (df.loc[df['Date'] < (timeline(1) + np.timedelta64(2, "D"))])['Weight'].mean() #Pulls baseline weight (I.e., all weight values made during baseline period)
        baselinemealsize = (df_m.loc[df_m['MealDates'] < (timeline(1) + np.timedelta64(2, "D"))])['MealSize'].mean()
        baselinepellets = (df_p.loc[df_p['PelletDates'] < (timeline(1) + np.timedelta64(2, "D"))])['TotalPellets'].mean()
    
        df.dropna(axis='index', how = 'any', inplace = True) #Resetting axis    
        df['Body Weight %'] = (df['Weight']/baselineweight)*100 #Normalising each value
        df['Meal Size %'] = (df_m['MealSize']/baselinemealsize)*100
        df['Total Pellets %'] = (df_p['TotalPellets']/baselinepellets)*100

        #Getting a column with treatment in it, and unique colours for each treatment (probably a more Pythonic way to do this):
        if "PFC" in rig: 
            T = 'PFC hM4Di'
            col = 'dodgerblue'
        if "IC" in rig:
            T = 'IC hM4Di'
            col = 'royalblue'
        if 'VRF' in rig:
            T = 'VRF'
            col = 'red'
        if '_c_rex' in rig:
            T = 'Control'
            col = 'black'
   
        days = np.arange(0, len(df), days_for_plots).tolist() #Goes up in 0.5 depending on number of entries (as RollingMedians returns 12 hour bins)
        days = pd.DataFrame(days, columns=['Days']) #Making it a df

        df.dropna(axis='index', how = 'any', inplace = True) #Resetting axis    
        df['Date'] = days['Days'] #Formatting days column
        df['Treatment'] = T #Formatting treatment column
        df['Cohort'] = cohort #Nice to retain this information
        df['Cage'] = int((rig[-1:])) #This is important as we need to model cage as a random effect in stats test (later) to account for potential pseudoreplication
        df['Colour'] = col
        
        #ax1.plot(df['Date'], df['Weight'], marker='o', linestyle = '-', color=col, alpha = 0.15)
        #ax2.plot(df['Date'], df['MealSize'], marker='o', linestyle = '-', color=col, alpha = 0.15) #(Un)comment to plot individual animals (Makes graph confusing)
        #ax3.plot(df['Date'], df['TotalPellets'], marker='o', linestyle = '-', color=col, alpha = 0.15) #(Un)comment to plot individual animals (Makes graph confusing)    

        day = -1 #Plotting index for habits analysis
        
        #Rhythmicity:
        for date in pd.date_range(start=timeline(1) - np.timedelta64(12, 'h'), end=timeline(1) + np.timedelta64(348, 'h'), freq = '1D').date: #For each day in the experiemnt
            try:
                day = day + 1

                axis = fig_axes[day]
                times_plot = times.loc[times['MealDates'] >= np.datetime64(date)]
                times_plot = times_plot.loc[times_plot['MealDates'] < (np.datetime64(date) + np.timedelta64(24, 'h'))]

                times_plot['MealDates'] = (times_plot['MealDates'].dt.hour) #Return just the hour value
                times_plot['Date'] = day

                count_avglist.append(times_plot['Counts'])
                datelist.append(times_plot['Date'])
                hourlist.append(times_plot['MealDates'])

                #axis.plot(times_plot['MealDates'], times_plot['Counts'], marker='o', linestyle = '-', color=col, alpha = 0)
                axis.set_title("Day: " + f'{day}')
                axis.set_ylabel("Meal Counts")
                if axis == ax10 or axis == ax11 or axis == ax14 or axis == ax15 or axis == ax19 or axis == ax20 or axis == ax22:
                    axis.set_xlabel("Time of Day (Hour)") #Just makes graphs intelligible
            except:
                break

        hourly_average = pd.DataFrame()
        hourly_average['Count_Meals'] = pd.concat(count_avglist)
        hourly_average['Hour'] = pd.concat(hourlist)
        hourly_average['Day'] = pd.concat(datelist)
        hourly_average['Treatment_Meals'] = T
        hourdf.append(hourly_average)

        allanimals.append(df) #Appending lists of animal data into a new object

    allanimals = pd.concat(allanimals)
    allanimals = pd.DataFrame(allanimals) 

    hourdf = pd.concat(hourdf)
    hourdf = pd.DataFrame(hourdf)
    
    allanimals = pd.concat([allanimals, hourdf], axis = 0)
    
    return allanimals #df can be returned for further analysis

def add_dailyavg(data): #Input PFCcohort2 and PFCcohort3 for example to return daily averages across cohorts
    print("Running daily average")
    avglist = []
    avgmeallist = []
    datelist = []
    errorlistw = []
    errorlistm = []
    errorlistp = []
    avgpelletslist = []
    treatmentlist = []

    Date = 0

    for days in np.unique((data['Date']).dropna()): #For loop across days (I.e., the loop will run for all unique day values)
        avgday = data.loc[data['Date'] == Date] #I.e., day 0, day 0.5...
        avgweight = avgday['Body Weight %'].mean() #Returns the average associated with that day
        avgmeal = avgday['Meal Size %'].mean()
        avgpellets = avgday['Total Pellets %'].mean()
        errorw = sp.tstd(avgday['Body Weight %']) #Finding the standard error mean associated with the weight values that day for error bars
        errorm = sp.tstd(avgday['Meal Size %'])
        errorp = sp.tstd(avgday['Total Pellets %'])
        treatment = (avgday['Treatment'].iloc[0])

        datelist.append(Date) #Appending these values to lists
        avglist.append(avgweight)
        avgmeallist.append(avgmeal)
        errorlistm.append(errorm)
        errorlistw.append(errorw)
        errorlistp.append(errorp)
        avgpelletslist.append(avgpellets)
        treatmentlist.append(treatment)

        if avgday['Treatment'].unique() == 'IC hM4Di': #Adds a colour column
            col = 'dodgerblue'
            n = 10
            alpha = 1
        if avgday['Treatment'].unique() == 'PFC hM4Di':
            col = 'royalblue'
            n = 9
            alpha = 1
        if avgday['Treatment'].unique() == 'VRF':
            col = 'red'
            n = 10
            alpha = 1
        if avgday['Treatment'].unique() == 'Control':
            col = 'black'
            n = 9
            alpha = 1

        Date = Date + days_for_plots #Run the loop again for the next bin

    daily_avg = pd.DataFrame()
    hourly_avg = pd.DataFrame()

    daily_avg['Date'] = datelist #Once all days gone through, make a dataframe consisting of all of these lists
    daily_avg['Body Weight %'] = avglist
    daily_avg['Meal Size %'] = avgmeallist
    daily_avg['Daily_SEM_W'] = errorlistw
    daily_avg['Daily_SEM_M'] = errorlistm
    daily_avg['Daily_SEM_P'] = errorlistp
    daily_avg['Total Pellets %'] = avgpelletslist
    daily_avg['Treatment'] = treatmentlist
    daily_avg['Colour'] = col

    hourly_avg['Date'] = (data['Day']).dropna()
    hourly_avg['Hour'] = (data['Hour']).dropna()
    hourly_avg['Count_Meals'] = (data['Count_Meals']).dropna()
    hourly_avg['Treatment_Meals'] = (data['Treatment_Meals']).dropna()

    countlist = []
    hourlist = []
    datelist = []
    axislist = []

    for date in np.arange(0, 15, 1): #For each day in the experiment
        hourly_avg_plt_d = hourly_avg.loc[hourly_avg['Date'] == date] #Grab a date...
        for hour in np.unique((hourly_avg_plt_d['Hour'])):
                hourly_avg_plt_h = hourly_avg_plt_d.loc[hourly_avg_plt_d['Hour'] == hour] #And extract entries from that day's hour values

                zeroes = n - len(hourly_avg_plt_h) #This is a way to include hours and values where zero meals were had. It's disabled for now, but you can turn it on by saying "if zeroes == 0"
                if zeroes == zeroes: #Disabling this zero'er for now.
                    pass
                else:
                    z_hours_list = []
                    z_date_list = []
                    z_counts_list = []
                    z_treatment_list = []

                    for i in range(zeroes):
                        z_hours_list.append(hour)
                        z_date_list.append(date)
                        z_counts_list.append(0)
                        z_treatment_list.append(hourly_avg_plt_h['Treatment_Meals'].iloc[0]) #Add lines containing day and hour no meals were had
                        
                    add_zeroes_df = pd.DataFrame()
                    add_zeroes_df['Hour'] = z_hours_list
                    add_zeroes_df['Date'] = z_date_list
                    add_zeroes_df['Count_Meals'] = z_counts_list
                    add_zeroes_df['Treatment_Meals'] = z_treatment_list
                    hourly_avg_plt_h = pd.concat([hourly_avg_plt_h, add_zeroes_df])
                                    
                axis = fig_axes[date]

                plot_hourly_val = hourly_avg_plt_h['Count_Meals'].mean()

                countlist.append(plot_hourly_val)
                hourlist.append(hour)
                axislist.append(axis)
                datelist.append(date)
              
    hourlydf = pd.DataFrame()
    hourlydf['Count_Meals'] = countlist
    hourlydf['Axis'] = axis
    hourlydf['Hour'] = hourlist
    hourlydf['Date'] = datelist

    j = -1

    for date in np.unique(hourlydf['Date']):
        j = j + 1
        hourlydf_plot = hourlydf.loc[hourlydf['Date'] == date]
        axis = fig_axes[j]
        error = sp.tstd(hourlydf_plot['Count_Meals']) #Plots standard deviation

        axis.plot(hourlydf_plot['Hour'], hourlydf_plot['Count_Meals'],  marker='o', linestyle = '-', color=col, alpha = alpha, label = avgday['Treatment'].iloc[0])
        axis.errorbar(hourlydf_plot['Hour'], hourlydf_plot['Count_Meals'], yerr=error, xerr = None, color = col, ls = None, alpha = alpha)
        ax8.legend()
    return daily_avg

IC = add_dailyavg(pd.concat([GetData("cohort2", "fem6_e_IChm4di_rex3"), 
                            GetData("cohort3", "fem11_e_IChm4di_rex4")]))
VRF = add_dailyavg(pd.concat([GetData("cohort3", "fem8_c_VRF_rex1"), 
                            GetData("cohort2", "fem7_c_VRF_rex4")]))
Control = add_dailyavg(pd.concat([GetData("cohort2", "fem5_c_rex2"), 
                            GetData("cohort3", "fem9_c_rex2")]))
PFC = add_dailyavg(pd.concat([GetData("cohort2", "fem4_e_PFChm4di_rex1"), 
                            GetData("cohort3", "fem10_e_PFChm4di_rex3")]))

#Plotting the figures
ax1.plot(IC['Date'], IC['Body Weight %'], marker='o', linestyle = '-', color=('royalblue'), label = 'IC, n=10')
ax1.plot(VRF['Date'], VRF['Body Weight %'], marker='o', linestyle = '-', color=('red'), label = 'VRF, n=10')
ax1.plot(Control['Date'], Control['Body Weight %'], marker='o', linestyle = '-', color=('black'), label = 'Control, n=9')
ax1.plot(PFC['Date'], PFC['Body Weight %'], marker='o', linestyle = '-', color=('dodgerblue'), label = 'PFC, n=9')
ax1.axvline(2, color='black', linestyle='dashed', label = "Induction") #These axv line commands plot axv lines at days where induction occurs
ax1.axvline(5, color='black', linestyle='dashed')
ax1.axvline(8, color='black', linestyle='dashed')
ax1.axvline(11, color='black', linestyle='dashed')
ax1.errorbar(IC['Date'], IC['Body Weight %'], yerr=IC['Daily_SEM_W'], xerr = None, color = 'royalblue', ls = None)
ax1.errorbar(VRF['Date'], VRF['Body Weight %'], yerr=VRF['Daily_SEM_W'], xerr = None, color = 'red', ls = None)
ax1.errorbar(Control['Date'], Control['Body Weight %'], yerr=Control['Daily_SEM_W'], xerr = None, color = 'black', ls = None)
ax1.errorbar(PFC['Date'], PFC['Body Weight %'], yerr=PFC['Daily_SEM_W'], xerr = None, color = 'dodgerblue', ls = None)
ax1.set_title("Average Body Weight")
ax1.set_ylabel("Body Weight (% of Baseline)")
ax1.grid(True)  

ax2.plot(IC['Date'], IC['Meal Size %'], marker='o', linestyle = '-', color=('royalblue'), label = 'IC, n=10')
ax2.plot(VRF['Date'], VRF['Meal Size %'], marker='o', linestyle = '-', color=('red'), label = 'VRF, n=10')
ax2.plot(Control['Date'], Control['Meal Size %'], marker='o', linestyle = '-', color=('black'), label = 'Control, n=9')
ax2.plot(PFC['Date'], PFC['Meal Size %'], marker='o', linestyle = '-', color=('dodgerblue'), label = 'PFC, n=9')
ax2.errorbar(IC['Date'], IC['Meal Size %'], yerr=IC['Daily_SEM_M'], xerr = None, color = 'royalblue', ls = None)
ax2.errorbar(VRF['Date'], VRF['Meal Size %'], yerr=VRF['Daily_SEM_M'], xerr = None, color = 'red', ls = None)
ax2.errorbar(Control['Date'], Control['Meal Size %'], yerr=Control['Daily_SEM_M'], xerr = None, color = 'black', ls = None)
ax2.errorbar(PFC['Date'], PFC['Meal Size %'], yerr=PFC['Daily_SEM_M'], xerr = None, color = 'dodgerblue', ls = None)
ax2.set_title("Average Meal Size")
ax2.set_ylabel("Meal Size (% of Baseline)")
ax2.axvline(2, color='black', linestyle='dashed', label = "Induction") #These axv line commands plot axv lines at days where induction occurs
ax2.axvline(5, color='black', linestyle='dashed')
ax2.axvline(8, color='black', linestyle='dashed')
ax2.axvline(11, color='black', linestyle='dashed')
ax2.grid(True)  

ax3.plot(IC['Date'], IC['Total Pellets %'], marker='o', linestyle = '-', color=('royalblue'), label = 'IC, n=10')
ax3.plot(VRF['Date'], VRF['Total Pellets %'], marker='o', linestyle = '-', color=('red'), label = 'VRF, n=10')
ax3.plot(Control['Date'], Control['Total Pellets %'], marker='o', linestyle = '-', color=('black'), label = 'Control, n=9')
ax3.plot(PFC['Date'], PFC['Total Pellets %'], marker='o', linestyle = '-', color=('dodgerblue'), label = 'PFC, n=9')
ax3.errorbar(IC['Date'], IC['Total Pellets %'], yerr=IC['Daily_SEM_P'], xerr = None, color = 'royalblue', ls = None)
ax3.errorbar(VRF['Date'], VRF['Total Pellets %'], yerr=VRF['Daily_SEM_P'], xerr = None, color = 'red', ls = None)
ax3.errorbar(Control['Date'], Control['Total Pellets %'], yerr=Control['Daily_SEM_P'], xerr = None, color = 'black', ls = None)
ax3.errorbar(PFC['Date'], PFC['Total Pellets %'], yerr=PFC['Daily_SEM_P'], xerr = None, color = 'dodgerblue', ls = None)
ax3.set_title("Daily Pellets Eaten")
ax3.set_ylabel("Pellets Eaten (% of Baseline)")
ax3.set_xlabel("Time (Days)")
ax3.axvline(2, color='black', linestyle='dashed', label = "Induction") #These axv line commands plot axv lines at days where induction occurs
ax3.axvline(5, color='black', linestyle='dashed')
ax3.axvline(8, color='black', linestyle='dashed')
ax3.axvline(11, color='black', linestyle='dashed')
ax3.grid(True)  

#Derivative of weight (rate of change) (Just to use up the final axis on the figure. Possibly interesting, possibly not)
ICweightdx = np.gradient(IC['Body Weight %'], IC['Date'], )
PFCweightdx = np.gradient(PFC['Body Weight %'], PFC['Date'])
VRFweightdx = np.gradient(VRF['Body Weight %'], VRF['Date'])
Controlweightdx = np.gradient(Control['Body Weight %'], Control['Date'])
ax4.plot(IC['Date'], ICweightdx, marker='o', linestyle = '-', color=('royalblue'), label = 'IC, n=10')
ax4.plot(VRF['Date'], VRFweightdx, marker='o', linestyle = '-', color=('red'), label = 'VRF, n=10')
ax4.plot(Control['Date'], Controlweightdx, marker='o', linestyle = '-', color=('black'), label = 'Control, n=9')
ax4.plot(PFC['Date'], PFCweightdx, marker='o', linestyle = '-', color=('dodgerblue'), label = 'PFC, n=9')
ax4.axvline(2, color='black', linestyle='dashed', label = "Induction") #These axv line commands plot axv lines at days where induction occurs
ax4.axvline(5, color='black', linestyle='dashed')
ax4.axvline(8, color='black', linestyle='dashed')
ax4.axvline(11, color='black', linestyle='dashed')
ax4.set_title("Rate of Change of Body Weight")
ax4.set_ylabel("Δ%BW/Δt")
ax4.set_xlabel("Time (Days)")
ax4.grid(True)  

#Concatenating data across animals into one variable for statistical test
data = []
for cohort in cohorts:
    for rig in rigs:
        try:
            data.append(GetData(cohort, rig))
            print("Running Unique Combination!")
        except Exception:
            print("Null combo! Passing...")
            pass
print("Concatenating complete!")
data_final = pd.concat(data)

def BarChart(day, parameter, axis): #Specify the day of analysis, the parameter you want (% BW, etc), and the axis to plot on (Each parameter gets its own axis)
    data_chartlist = [] #List for t-test later
    for treatment in data_final['Treatment'].dropna().unique(): #Iterate over unique treatments
        data = data_final.loc[data_final['Treatment'] == treatment] 
        data = data.loc[data['Date'] == day]
       
        #data = data.loc[data['Date'] <= day + 0.5] #Can include the other 12 hour bin if desired by uncommenting
        data_chart = pd.DataFrame(data[[parameter, 'Treatment', 'Animal']])

        order = ['Control', 'VRF', 'PFC hM4Di', 'IC hM4Di']

        if treatment == 'Control':
            data['Colour'].iloc[0] = 'gray'
   
        sb.barplot(data_chart, x = data_chart['Treatment'], order = order, 
                   y = parameter, color = data['Colour'].iloc[0], errorbar = 'sd', ax = axis) #Using seaborn to plot barplot
        
        sb.stripplot(x = data_chart['Treatment'], order = order, 
                     y = parameter, data = pd.DataFrame(data_chart[parameter]), color = data['Colour'].iloc[0], alpha = 1, 
                     edgecolor = 'black', linewidth = 1, ax = axis) #Stripplot for individual data points
        data_chartlist.append(data_chart) #Appending a list for t-test

    data_chart = pd.concat(data_chartlist)
    data_chart = pd.DataFrame(data_chart)
    pairs = list(itertools.combinations(np.unique(data_chart['Treatment']), 2)) #All the t-test pairs
   
    annotate = sn.Annotator(axis, pairs, plot = 'barplot', data=data_chart, x = "Treatment", 
                            y = parameter, loc = 'inside', order = order, fontsize = 'small') #Using statannotations to annotate with t-test
    
    annotate.configure(test='Mann-Whitney', verbose = 3, text_format = 'simple', 
                            alpha = 0.01) #Non-parametric to be sensitive to assumption violations across parameters, plus p-value threshold
    
    print("Day:", day, "Parameter:", parameter) #So you can see what you ran t-tests for in the terminal
    annotate.apply_and_annotate() #It should spit out p-values too

    ax6.set_ylim([80, 160]) #Need to set sensible ylims for BW chart, it looks ridiculous otherwise

    axis.set_title("Day: " + f'{day}' + ", " + " Parameter: " f'{parameter}')

analysis_day = 5

BarData = BarChart(analysis_day, 'Meal Size %', ax5)
BarData = BarChart(analysis_day, 'Body Weight %', ax6)
BarData = BarChart(analysis_day, 'Total Pellets %', ax7)

ax1.axvline(analysis_day, color='blue', linestyle='dashed')
ax2.axvline(analysis_day, color='blue', linestyle='dashed')
ax3.axvline(analysis_day, color='blue', linestyle='dashed')
ax4.axvline(analysis_day, color='blue', linestyle='dashed', label = "t-test") #Highlights the day the t-tests are being run for on the graphs

ax4.legend(fontsize = "9")
plt.subplots_adjust(left=0.1, bottom=0.1, right=0.9, 
                    top=0.9, wspace=0.3, hspace=0.3)
plt.show()
