import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from datetime import datetime, date, timedelta, time
import matplotlib.dates as mdates
import scipy.stats as sp
import statsmodels.api as sm
import statsmodels.formula.api as smf

def Rolling_Medians(X, Y): #X defines the cohort, Y defines the rig (I.e., Rolling_Medians("cohort 2", "fem4_e_PFChm4di_rex1", 1)) is cohort 2, PFC
    win = 30
    savepath="C:\\Users\\robbi\\Documents\\GitHub\mazerex2\\" + X + "\\" + Y + "\\" #standard stuff
    known_tags = np.array(pd.read_csv(savepath + "AnimalTags.csv", header=None)).ravel().tolist() 
    time_line = np.array(pd.read_csv(savepath+"TimeLine.csv", header=None)).astype(np.datetime64).reshape(-1,1) 

    def timeline(X):  #for loop to format dates so arrays can be sliced (doesn't work otherwise)
        for date in time_line[X]:
            if date != 2026:
                break
        return(date)
    
    bin_size = 12 #For unilateral control of bins variable

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
        list_weightmeds.append(rolling_medians)
        list_weights.append(y)

    list_daily_w=[]
    allanimals = []
    an=-1
    for rfid in known_tags: #for loop across animals
        an=an+1
        data={
            "Date":list_x[an],
            "Weight":list_weights[an],
            "Animal":int(known_tags[an]) #Adding animal tag
            }   
        filtered_df=pd.DataFrame(data)
        df = filtered_df.groupby(pd.Grouper(key='Date', freq = f'{bin_size}' + "h", origin=str(time_line[1])[2:12])).median().reset_index()
        baselineweight = (df.loc[df['Date'] < (timeline(1) + np.timedelta64(2, "D"))])['Weight'].mean() #Pulls baseline weight (I.e., all weight values made during baseline period)
        df['Weight'] = (df['Weight']/baselineweight)*100 #Normalising each value

        #Getting a column with treatment in it, and unique colours for each treatment (probably a more Pythonic way to do this):
        if "PFC" in Y: 
            T = 'PFC hM4Di'
            col = 'dodgerblue'
        if "IC" in Y:
            T = 'IC hM4Di'
            col = 'royalblue'
        if 'VRF' in Y:
            T = 'VRF'
            col = 'red'
        if '_c_rex' in Y:
            T = 'Control'
            col = 'black'
   
        days = np.arange(0, len(df), 0.5).tolist() #Goes up in 0.5 depending on number of entries (as RollingMedians returns 12 hour bins)
        days = pd.DataFrame(days, columns=['Days']) #Making it a df

        df.dropna(axis='index', how = 'any', inplace = True) #Resetting axis    
        df['Date'] = days['Days'] #Formatting days column
        df['Treatment'] = T #Formatting treatment column
        df['Cohort'] = X #Nice to retain this information
        df['Cage'] = int((Y[-1:]))
        df['Cage'] = pd.Series.to_frame(df['Cage']) #This is important as we need to model cage as a random effect in stats test (later) to account for potential pseudoreplication
        #ax.plot(df['Date'], df['Weight'], marker='o', linestyle = '-', color=col, alpha = 0.15) #(Un)comment to plot individual animals (Makes graph confusing)
        allanimals.append(df) #Appending lists of animal data into a new object
    allanimals = pd.concat(allanimals) #Concatenate into on df
    allanimals = pd.DataFrame(allanimals)
    print(allanimals)
    return allanimals #df can be returned for further analysis

def add_dailyavg(data): #Input PFCcohort2 and PFCcohort3 for example to return daily averages across cohorts
    data = data
    avglist = []
    datelist = []
    errorlist = []
    Date = 0

    for days in np.unique(data['Date']): #For loop across days (I.e., the loop will run for all unique day values)
        avgday = data.loc[data['Date'] == Date] #I.e., day 0, day 0.5...
        avgweight = avgday['Weight'].mean() #Returns the average associated with that day
        error = sp.sem(avgday['Weight']) #Finding the standard error mean associated with the weigth values that day for error bars

        datelist.append(Date) #Appending these values to lists
        avglist.append(avgweight)
        errorlist.append(error)
        Date = Date + 0.5 #Run the loop again for the next 12 hour bin
    daily_avg = pd.DataFrame() #Once all days gone through, make a dataframe consisting of all of these lists
    daily_avg['Weight'] = avglist
    daily_avg['Date'] = datelist
    daily_avg['Daily_SEM'] = errorlist

    print(daily_avg)

    return daily_avg


#Now run the class for all of our treatments *ACROSS COHORTS!*
IC = add_dailyavg(pd.concat([Rolling_Medians("cohort2", "fem6_e_IChm4di_rex3"), 
                            Rolling_Medians("cohort3", "fem11_e_IChm4di_rex4")]))

VRF = add_dailyavg(pd.concat([Rolling_Medians("cohort3", "fem8_c_VRF_rex1"), 
                            Rolling_Medians("cohort2", "fem7_c_VRF_rex4")]))

Control = add_dailyavg(pd.concat([Rolling_Medians("cohort2", "fem5_c_rex2"), 
                            Rolling_Medians("cohort3", "fem9_c_rex2")]))

PFC = add_dailyavg(pd.concat([Rolling_Medians("cohort2", "fem4_e_PFChm4di_rex1"), 
                            Rolling_Medians("cohort3", "fem10_e_PFChm4di_rex3")]))

print(PFC)


#Plotting the figure
ax.plot(IC['Date'], IC['Weight'], marker='o', linestyle = '-', color=('royalblue'), label = 'IC, n=10')
ax.plot(VRF['Date'], VRF['Weight'], marker='o', linestyle = '-', color=('red'), label = 'VRF, n=10')
ax.plot(Control['Date'], Control['Weight'], marker='o', linestyle = '-', color=('black'), label = 'Control, n=9')
ax.plot(PFC['Date'], PFC['Weight'], marker='o', linestyle = '-', color=('dodgerblue'), label = 'PFC, n=9')
ax.axvline(2, color='black', linestyle='dashed', label = "Induction") #These axv line commands plot axv lines at days where induction occurs
ax.axvline(5, color='black', linestyle='dashed')
ax.axvline(8, color='black', linestyle='dashed')
ax.axvline(11, color='black', linestyle='dashed')
ax.errorbar(IC['Date'], IC['Weight'], yerr=IC['Daily_SEM'], xerr = None, color = 'royalblue', ls = None)
ax.errorbar(VRF['Date'], VRF['Weight'], yerr=VRF['Daily_SEM'], xerr = None, color = 'red', ls = None)
ax.errorbar(Control['Date'], Control['Weight'], yerr=Control['Daily_SEM'], xerr = None, color = 'black', ls = None)
ax.errorbar(PFC['Date'], PFC['Weight'], yerr=PFC['Daily_SEM'], xerr = None, color = 'dodgerblue', ls = None)
plt.title("Average Body Weight")
plt.ylabel("% Body Weight")
plt.legend()
plt.xlabel("Time (Days)")
plt.grid(True)  
plt.show()

#Concatenating the data into one variable for statistical test
data_final = pd.concat([Rolling_Medians("cohort2", "fem6_e_IChm4di_rex3"), Rolling_Medians("cohort3", "fem11_e_IChm4di_rex4"),
                    Rolling_Medians("cohort2", "fem5_c_rex2"), Rolling_Medians("cohort3", "fem9_c_rex2"),
                    Rolling_Medians("cohort3", "fem8_c_VRF_rex1"), Rolling_Medians("cohort2", "fem7_c_VRF_rex4"),
                    Rolling_Medians("cohort2", "fem4_e_PFChm4di_rex1"), Rolling_Medians("cohort3", "fem10_e_PFChm4di_rex3")])

#Linear mixed effects model (Considering our study design)
model = smf.mixedlm("Weight ~ Treatment", data_final, groups=data_final["Animal"], re_formula = '0 + Cage') #Models weight with predictor variables as treatment and date. We need to model animal as a random effect to account for potential pseudoreplication
#(For above) We need to model animal AND cage as random effects to account for 
#pseudoreplication. We might consider also modelling cohort too to account for cohort specific effects
model = model.fit()
residuals = model.resid
fitted = model.fittedvalues
summary = model.summary()

#Testing assumptions:
#Q-Q plot for normality of residuals
sp.probplot(residuals, dist="norm", plot=plt)
plt.title("Q-Q Plot")
plt.show()
 
#Residuals vs fitted values to check equal variance of residuals
plt.scatter(fitted, residuals)
plt.axhline(y=0, color='r', linestyle='--')
plt.xlabel("Fitted Values")
plt.ylabel("Residuals")
plt.title("Residuals Plot")
plt.show()

#Both of these look good to me.
print(summary) #Printing model



















