# -*- coding: utf-8 -*-
"""
Created on Sat Oct  7 22:29:07 2023

@author: johnmallon
"""


import pandas as pd
import seaborn
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score 
#import equationsXY_version_01 as eq
import random
from scipy.signal import savgol_filter
import scipy.stats as st


def groupCellsRandomly(df,groupSize,maxCellsToGraph):
    #grab cells in bunches so we can graph them in groups
    cellID = df['cell ID'].unique()
    count = 0
    alreadyPicked = []
    groupsOfCells = []
    breakOut = False
    while count < maxCellsToGraph:
        if breakOut:
            break
        randomList = []
        i = 0
        while i < groupSize:
            n = random.randint(0,len(cellID)-1)
            if not(n in alreadyPicked):
                alreadyPicked.append(n)
                randomList.append(cellID[n])
                i+=1
            if len(alreadyPicked) == len(cellID):  #for when the number of cells in the dataset is less than the groupSize/amount to graph
                breakOut = True
                break
        groupsOfCells.append(randomList)
        count+=groupSize
    
    return groupsOfCells


def smoothDatapoints(subDF):
    #smooth the growth data with a moving average
    window_size =5   # Window size for smoothing (should be odd)
    order =2  # Polynomial order
    y_smooth = savgol_filter(y, window_size, order)
    y_6th = y_smooth[::1]
    x_6th = x[::1]
    #print(x_6th)
    x_norm = [item - x_6th[0] for item in x_6th]
    
    return newX,newY

    
def exponentialEquation(x,y,k,a,c):
    return a*np.exp(k*x)+0  #do we use 0 or a variable?    
    
    
def fitExponential(df):
    lineages = []
    parents = []
    cellIDs = []
    shapes = []
    r2 = []
    growthRates = []
    cellCycleTime= []
    areasBirth = []
    areasDivision = []
    cellsCycleAreaRatio = []
    divisionAreaRatio = []
    growthDF = pd.DataFrame(columns=['Lineage','cell ID','R2','Growth Rate','Doubling Time','Cell Cycle Time','Area at Birth','Area at Division'])
    
    xPoints = []
    yFitPoints = []
    cellIDsss = []
    fitsDF = pd.DataFrame(columns=['cell ID','X','yFit'])
    
    cells = df.groupby(['Lineage','cell ID'])
    for cellID, cell in cells:
        lineages.append(cellID[0])
        cellIDs.append(cellID[1])
        parents.append(cell['Parent'].iloc[0])
        shapes.append(cell['Shape'].iloc[0]) #shapes recorded as whatever the shape was at the start of growth
        time = cell['Relative Time']
        area = cell['Area']
        areasBirth.append(area.iloc[0])
        areasDivision.append(area.iloc[-1])
        
        if not(pd.isna(cell['Parent'].iloc[0])):
            divisionAreaRatio.append(cell['Parent Area'].iloc[0]/cell['Area'].iloc[0])
        else:
            divisionAreaRatio.append(0) #append zero if NaN because we don't know the parent
                
                
        
        
        cellsCycleAreaRatio.append(area.iloc[-1]/area.iloc[0])
        cellCycleTime.append(time.iloc[-1]-time.iloc[0])
        if ((len(time)==len(area)) & (len(time)>1)):
            #initVals = [0.05,3]
            popt, pcov = curve_fit(exponentialEquation,time,area)
            yFit = exponentialEquation(time,*popt)
            rSquaredValue =  r2_score(area, yFit)
            r2.append(rSquaredValue)
            growthRates.append(popt[1])
            xPoints.extend(time)
            yFitPoints.extend(yFit)
            size = len(time)
            cellIDsss.extend([cellID[1]]*size)
            
        else:
            r2.append(0)
            growthRates.append(0)
    
    growthDF['Lineage'] = lineages
    growthDF['Parent'] = parents
    growthDF['cell ID'] = cellIDs
    growthDF['Shape'] = shapes
    growthDF['R2'] = r2
    growthDF['Growth Rate'] = growthRates
    growthDF['Doubling Time'] = np.log(2)/growthDF['Growth Rate']
    growthDF['Cell Cycle Time'] = cellCycleTime
    growthDF['Area at Birth'] = areasBirth
    growthDF['Area at Division'] = areasDivision
    growthDF['Lifetime Area Ratio'] = cellsCycleAreaRatio
    growthDF['Division Area Ratio'] = divisionAreaRatio
    
    fitsDF['X'] = xPoints
    fitsDF['yFit'] = yFitPoints
    fitsDF['cell ID'] =  cellIDsss
    
    return growthDF, fitsDF


def createGraphs(yDataList,lineages,yTitle,directory):
    fig, ax = plt.subplots(dpi=300) #Growth Rates
    GREY_LIGHT = "#b4aea9"
    GREY50 = "#7F7F7F"
    BLACK = "#282724"
    GREY_DARK = "#747473"
    RED_DARK = "#850e00"
    COLOR_SCALE = ["#1B9E77", "#D95F02", "#7570B3"]
  
    # Horizontal positions for the violins. 
    # They are arbitrary numbers. They could have been [-1, 0, 1] for example.
    POSITIONS = list(range(len(yDataList)))
    

    # Horizontal lines, figure out a way to make this universal to any set of dataranges
    #HLINES = list(range(0,maxValue,stepSize))
    # Horizontal lines that are used as scale reference
    #for h in HLINES:
     #   ax.axhline(h, color=GREY50, ls=(0, (5, 5)), alpha=0.8, zorder=0)
    
    
   
    
    # Add violins ----------------------------------------------------
    # bw_method="silverman" means the bandwidth of the kernel density
    # estimator is computed via Silverman's rule of thumb. 
    # More on this in the bonus track ;)

    # The output is stored in 'violins', used to customize their appearence
    seaborn.violinplot(
        
        )
    violins = ax.violinplot(
        yDataList, 
        positions=POSITIONS,
        widths=0.45,
        bw_method="silverman",
        showmeans=False, 
        showmedians=False,
        showextrema=False
        )

    # Customize violins (remove fill, customize line, etc.)
    for pc in violins["bodies"]:
        pc.set_facecolor("none")
        pc.set_edgecolor(BLACK)
        pc.set_linewidth(1.4)
        pc.set_alpha(1)
        pc.set_zorder(5)
    
    # Add boxplots ---------------------------------------------------
    # Note that properties about the median and the box are passed
    # as dictionaries.

    medianprops = dict(
        linewidth=2, 
        color=GREY_DARK,
        solid_capstyle="butt"
        )
    boxprops = dict(
        linewidth=1, 
        color=GREY_DARK
        )

    ax.boxplot(
        yDataList,
        positions=POSITIONS, 
        showfliers = False, # Do not show the outliers beyond the caps.
        showcaps = False,   # Do not show the caps
        medianprops = medianprops,
        whiskerprops = boxprops,
        boxprops = boxprops,
        widths=0.2
        )
    
    # Create jittered version of "x" (which is only 0, 1, and 2)
    # More about this in the bonus track!
    jitter = 0.075
    x_data = [np.array([i] * len(d)) for i, d in enumerate(yDataList)]
    x_jittered = [x + st.t(df=6, scale=jitter).rvs(len(x)) for x in x_data]
    
    # Add jittered dots ----------------------------------------------
    for x, y, color in zip(x_jittered, yDataList, COLOR_SCALE):
        ax.scatter(x, y, s = 100, color=color, alpha=0.25)
    
    
    
    #labels of means
    # Add mean value labels ------------------------------------------
    
    means = [np.mean(y) for y in yDataList]
    
    for i, mean in enumerate(means):
    # Add dot representing the mean
        ax.scatter(i, mean, s=100, color=RED_DARK, zorder=3)
    
    # Add line conecting mean value and its label
        ax.plot([i, i + 0.25], [mean, mean], ls="dashdot", color="black", zorder=3)
    
    # Add mean value label.
        ax.text(
        i + 0.25,
        mean,
        r"$\hat{\mu}_{\rm{mean}} = $" + str(round(mean, 2)),
        fontsize=4,
        va="center",
        bbox = dict(
            facecolor="white",
            edgecolor="black",
            boxstyle="round",
            pad=0.5
            ),
        zorder=10 # to make sure the line is on top
        )
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines["left"].set_color(GREY_LIGHT)
    ax.spines["left"].set_linewidth(2)
    ax.spines["bottom"].set_color(GREY_LIGHT)
    ax.spines["bottom"].set_linewidth(2)
    ax.tick_params(length=0)
    #ax.set_yticks(HLINES)
    #ax.set_yticklabels(HLINES, size=15)
    ax.set_ylabel(yTitle, size=18, weight="bold")
    
    # xlabels accounts for the sample size for each lineage
    xlabels = [f"{lineage}\n(n={len(yDataList[i])})" for i, lineage in enumerate(lineages)]
    ax.set_xticks(POSITIONS)
    ax.set_xticklabels(xlabels, size=15, ha="center", ma="center")
    ax.set_xlabel("Lineages", size=18, weight="bold")
    
    #save figure
    #fig.savefig(directory+yTitle+'.svg',format='svg')
    
    
def main():
    ###############################################
    ####### USER INPUT INFORMATION START HERE ###########
    
    directory = r"C:\Users\bisso\Desktop\\"
    file = "analysisWithParentArea.csv"
    #set colors for graphs
    qualitative_colors = seaborn.color_palette("Set3")
    groupSize = 10
    maxCellsToGraph = 0
    goodnessOfFit = 0.90
    
    ####### USER INPUT INFORMATION END HERE ###########
    ##############################################
    
    
    df = pd.read_csv(directory+file)
    
    ##############################################
    ####FILTER OUT SHORT TRACKS AND BAD TRACKS####
    
    print('track length filter')
    print(df.shape[0])
    #get an idea of how long cell tracks in time are with this plot
    cells = df.groupby('cell ID')
    #cellTracks = cells.size()
    #figure, sx = plt.subplots(dpi=300)
    #seaborn.swarmplot(ax=sx,
     #               y=cellTracks,
      #              edgecolor="k",linewidth=0.75,alpha=0.5,
       #             zorder=2.0)
    
    df = cells.filter(lambda x: len(x) > 5)
    print(df.shape[0])
    
    ####FILTER OUT SHORT TRACKS AND BAD TRACKS####
    ##############################################
    
    
    #graph growth of single cells that have good exponential fits
    growthDF, fitsDF = fitExponential(df)
    groupsOfCells = groupCellsRandomly(df,groupSize,maxCellsToGraph)
    for cells in groupsOfCells:
        figure, ax = plt.subplots(dpi=300)
        figure, bx = plt.subplots(dpi=300)
        for cell in cells:
            subDF = df.loc[df['cell ID'] == cell]
            print(cell)
            print(fitsDF)
            fit = fitsDF.loc[fitsDF['cell ID']==cell]
            growth = growthDF.loc[growthDF['cell ID']==cell]
            if growth['R2'].iloc[0] >= goodnessOfFit:
                
                seaborn.scatterplot(data=subDF,ax=ax,
                                x='Relative Time',y='Area',
                                edgecolor="k",linewidth=0.75,alpha=0.5,
                                zorder=2.0)
                seaborn.lineplot(data=fit,ax=ax,
                                x='X',y='yFit',
                                color="k",linewidth=0.75,alpha=1.0,
                                zorder=1.0)
                
                seaborn.scatterplot(data=subDF,ax=bx,
                                x='Relative Time',y='Solidity',
                                edgecolor="k",linewidth=0.75,alpha=0.5,
                                zorder=2.0)
                bx.set_ylim(0.8,1,0.1)
    
    #separate rods and disks, and look at their shape measurements
    shapes = df['Shape'].unique()
    for shape in shapes:
        subDF = df.loc[df['Shape']==shape]
        lineages = sorted(subDF["Lineage"].unique())
        AR = [subDF[subDF["Lineage"] == lineage]["Aspect Ratio"].values for lineage in lineages] 
        circ = [subDF[subDF["Lineage"] == lineage]["Circularity"].values for lineage in lineages]
        
        createGraphs(AR,lineages,'Aspect Ratio ('+ shape +')',directory)
        createGraphs(circ,lineages,'Circularity ('+ shape +')',directory)
        
    
    #filter out bad fits for graphing
    print('R2 filter')
    print(growthDF.shape[0])
    growthDF = growthDF.loc[(growthDF['R2']>=goodnessOfFit)]
    print(growthDF.shape[0])
    
    #medianGR = growthDF['Growth Rate'].median()
    #medianDT = growthDF['Doubling Time'].median()
    #medianCT = growthDF['Cell Cycle Time'].median()

    fig2, cx = plt.subplots(dpi=300) #Growth Rate vs. Cycle Time
    seaborn.scatterplot(data=growthDF,ax=cx,
                x='Cell Cycle Time',y='Growth Rate',
                color='k',edgecolor="k",linewidth=0.75,alpha=0.5,
                zorder=1.0)
    #fig2.savefig(directory+'grVscycle.svg',format='svg')
    
    #sort data into lineages and graph
    lineages = sorted(growthDF["Lineage"].unique())
    growthRates = [growthDF[growthDF["Lineage"] == lineage]["Growth Rate"].values for lineage in lineages]
    doublingTimes = [growthDF[growthDF["Lineage"] == lineage]["Doubling Time"].values for lineage in lineages]
    cellCycleTimes = [growthDF[growthDF["Lineage"] == lineage]["Cell Cycle Time"].values for lineage in lineages]
    areaBirth = [growthDF[growthDF["Lineage"] == lineage]["Area at Birth"].values for lineage in lineages]
    areaDivision = [growthDF[growthDF["Lineage"] == lineage]["Area at Division"].values for lineage in lineages]
    areaRatio = [growthDF[growthDF["Lineage"] == lineage]["Lifetime Area Ratio"].values for lineage in lineages]
    divisionRatio = [growthDF[growthDF["Lineage"] == lineage]["Division Area Ratio"].values for lineage in lineages]
    
    
    createGraphs(growthRates,lineages,'Growth Rate',directory)
    createGraphs(doublingTimes,lineages,'Doubling Time',directory)
    createGraphs(cellCycleTimes,lineages,'Cell Cycle Time',directory)
    createGraphs(areaBirth,lineages,'Area at Birth',directory)
    createGraphs(areaDivision,lineages,'Area at Division',directory)
    createGraphs(areaRatio,lineages,'Lifetime Area Ratio',directory)
    createGraphs(divisionRatio,lineages,'Division Area Ratio',directory)
    
    #print(growthDF.columns)
    reorder = ['Lineage', 'Parent','cell ID','Shape','R2', 'Growth Rate', 'Doubling Time',
           'Cell Cycle Time', 'Area at Birth', 'Area at Division',
           'Lifetime Area Ratio', 'Division Area Ratio']

    outputDF = growthDF[reorder]
    print(outputDF)
    #save dataframe as a .csv file
    outputDF.to_csv(directory+'growthData.csv',index=False)
    
main()