import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import ROOT

from matplotlib.colors import LogNorm
from matplotlib.backend_bases import MouseButton

def terminal_input():
    run_list = []; ch_list = []
    run_list_string = input("Input run list (use space between numbers): ")
    ch_list_string = input("Input channel list (use space between numbers): ")
    run_list_split = run_list_string.split(" ")
    ch_list_split = ch_list_string.split(" ")

    for i in range(np.size(run_list_split)):
        run_list.append(int(run_list_split[i])) 

    for i in range(np.size(ch_list_split)):
        ch_list.append(int(ch_list_split[i]))
    
    return run_list, ch_list

def tellme(title):
    plt.title(title, fontsize=12)
    plt.draw()

def ntuple_to_pd(path,runs,channels):
    
    mylist=[]
    for run in runs:
        aux=[]
        df1 = ROOT.RDataFrame("ntuple",path+"run%i_NTuple.root" %run)
        Qdf1 = ROOT.RDataFrame("charge",path+"run%i_NTuple.root" %run)
        diccionario = df1.AsNumpy()
        Qdiccionario=Qdf1.AsNumpy()
        aspd = pd.DataFrame(diccionario,columns=None)
        Qaspd=pd.DataFrame(Qdiccionario,columns=None)
        all_pd=pd.DataFrame.join(aspd,Qaspd)
        
        for ch in channels:
            aux.append(all_pd[all_pd["ch"]==ch])
        mylist.append(aux)
    
    return mylist

def charge_map(dfin,run,ch,charge,interactive):

    f90 = dfin[charge[0]]/dfin[charge[1]]
    np_f90 = f90.to_numpy()
    np_amp = dfin["Amp"].to_numpy()
    limits = [[0,0],[0,0]]

    ymin, ymax, xmin, xmax = np.min(np_f90), np.max(np_f90), np.min(np_amp), np.max(np_amp)  
    
    if ymin < 0:
        ymin = 0
    if xmax > 20000:
        xmax = 20000
    if ymax > 1:
        ymax = 1
    if xmin < 0:
        xmin = 0
    
    hist = plt.hist2d(np_amp,np_f90,200,[[xmin,xmax],[ymin,ymax]],norm=LogNorm())    
    plt.xlabel("Amp (ADC counts)", fontsize=12)
    plt.ylabel(charge[0]+"/"+charge[1], fontsize=12)  
    
    plt.colorbar(hist[3])

    if interactive == True:
        
        while True:
            plt.ion()
            tellme("(RUN %i CH %i) Select 2 corners with your mouse"%(run,ch))
            limits = plt.ginput(n=2, timeout=-1, show_clicks=True, mouse_add=MouseButton.LEFT)

            if limits[0][0] < limits[1][0] and limits[0][1] < limits[1][1]:
                amp_min, amp_max = limits[0][0], limits[1][0]
                f90_min, f90_max = limits[0][1], limits[1][1]
            
            elif limits[0][0] > limits[1][0] and limits[0][1] < limits[1][1]:
                amp_min, amp_max = limits[1][0], limits[0][0]
                f90_min, f90_max = limits[0][1], limits[1][1]

            elif limits[0][0] < limits[1][0] and limits[0][1] > limits[1][1]:
                amp_min, amp_max = limits[0][0], limits[1][0]
                f90_min, f90_max = limits[1][1], limits[0][1]

            elif limits[0][0] > limits[1][0] and limits[0][1] > limits[1][1]:
                amp_min, amp_max = limits[1][0], limits[0][0]       
                f90_min, f90_max = limits[1][1], limits[0][1]
            
            else:
                tellme("Same point selected twice! Click to continue")
                plt.waitforbuttonpress()
                continue
            
            tellme('Happy? Key click for yes, mouse click for no')
            
            plt.plot([amp_min,amp_max],[f90_min,f90_min],color="k")
            plt.plot([amp_min,amp_max],[f90_max,f90_max],color="k")
            plt.plot([amp_min,amp_min],[f90_min,f90_max],color="k")
            plt.plot([amp_max,amp_max],[f90_min,f90_max],color="k")

            if plt.waitforbuttonpress(-1):
                break
            
            plt.clf()

            hist = plt.hist2d(np_amp,np_f90,200,[[xmin,xmax],[ymin,ymax]],norm=LogNorm())    
            plt.title("Charge map for rate calculation and particle ID (RUN %i CH %i)"%(run,ch)) 
            plt.xlabel("Amp (ADC counts)", fontsize=12)
            plt.ylabel(charge[0]+"/"+charge[1], fontsize=12)  
            plt.colorbar(hist[3])

    else:
        #___ASK IN TERMINAL FOR F90 RANGE TO COUNT EVENTS___
        print("\nINPUT RANGE FOR EVENT COUNTING\n")
        
        plt.show()
        amp_min = float(input("Enter amp min (0.0-1e5): "))
        amp_max = float(input("Enter amp max (0.0-1e5): "))
        f90_min = float(input("Enter F90 min (0.0-1.0): "))
        f90_max = float(input("Enter F90 max (0.0-1.0): "))
    
    plt.close("all")  

    #___COUNT F90 VALUES INSIDE GIVEN RANGE___
    dfout = dfin.loc[(dfin["Amp"] < amp_max) & (dfin["Amp"] > amp_min) & (f90 < f90_max) & (f90 > f90_min)]

    #___PRINT NUMBER OF SELECTED EVENTS AND CALCULATE RATE___
    time = dfout["time"].to_numpy()
    if np.size(time) > 0:
        total_time = (time[-1]-time[0])*8e-9
        nevents = dfout["evt"].size
        rate = float(format(dfout["evt"].size/total_time,".4f"))
    
    else:
        print("\nERROR: NO EVENTS SELECTED IN CHARGE MAP")
        rate = 0
        nevents = 0

    #plt.ioff()
    return rate,nevents,dfout

def q_vs_amp(dfin,run,ch,charge):
    
    np_q = dfin[charge].to_numpy()
    np_amp = dfin["Amp"].to_numpy()

    plt.scatter(np_amp, np_q, s = 25)
    plt.title("Q vs. Amp (RUN %i CH %i)"%(run,ch)) 
    plt.xlabel("Amp (ADC counts)", fontsize=16)
    plt.ylabel(charge, fontsize=16)
    plt.show()

def lin_fit(dfin,run,ch,charge):

    y_scatt = dfin[charge].to_numpy()
    x_scatt = dfin["Amp"].to_numpy()

    #plt.scatter(x_scatt,y_scatt)

    linear_model = np.polyfit(x_scatt,y_scatt,1)
    linear_model_fn = np.poly1d(linear_model)

    x = np.linspace(np.min(x_scatt),np.max(x_scatt))
    plt.plot(x,linear_model_fn(x),label = "Linear fit run %i and ch %i"%(run,ch))
    
    return linear_model_fn