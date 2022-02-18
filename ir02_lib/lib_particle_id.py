import numpy as np
import matplotlib.pyplot as plt

from matplotlib.colors import LogNorm
from matplotlib.backend_bases import MouseButton
from ir02_lib.lib_tellme import tellme
    
def particle_id(dfin,run,ch,month,charge,automatic="yes",zoom=False):
    
    yes_string = ["yes", "y", "YES", "Y", "Yes"]
    no_string = ["no", "n", "NO", "N", "No"]

    if automatic in yes_string:
        auto = True
        print("\nAUTOMATIC MODE SELECTED\n")
    elif automatic in no_string:
        auto = False
        print("\nMANUAL MODE SELECTED\n")
    
    f90 = dfin[charge[0]]/dfin[charge[1]]
    np_f90 = f90.to_numpy()
    np_amp = dfin["Amp"].to_numpy()
    limits = [[0,0],[0,0]]

    #CLEAN DATA
    clean_f90 = []
    clean_amp = []
    for i in range (0,np.size(np_f90)):
        if 0 < np_f90[i] < 1 and 0 < np_amp[i]:
            clean_f90.append(np_f90[i])
            clean_amp.append(np_amp[i])
        else:
            continue
    np_f90 = np.array(clean_f90)
    np_amp = np.array(clean_amp)

    #ZOOM TO DATA
    if zoom == True:
        ymin, ymax, xmin, xmax = np.min(np_f90), np.max(np_f90), np.min(np_amp), np.max(np_amp)

        if ymin < 0:
            ymin = 0
        if xmax > 20000:
            xmax = 20000
        if ymax > 1:
            ymax = 1
        if xmin < 0:
            xmin = 0
    else:
        ymin, ymax, xmin, xmax = 0, 1, 0, np.max(np_amp)
    
    counts, xedges, yedges = np.histogram2d(np_amp,np_f90, bins=200)
    x_ind, y_ind = np.unravel_index(np.argmax(counts), counts.shape)
    xpeak = (xedges[x_ind]+xedges[x_ind+1])/2; ypeak = (yedges[y_ind]+yedges[y_ind+1])/2
    print(f'The maximum count is {counts[x_ind][y_ind]:.0f} at index ({x_ind}, {y_ind})')
    print(f'Between x values {xedges[x_ind]} and {xedges[x_ind+1]}')
    print(f'and between y values {yedges[y_ind]} and {yedges[y_ind+1]}')
    
    plt.ioff
    with plt.ion():
        
        plthist = plt.hist2d(np_amp,np_f90,200,[[xmin,xmax],[ymin,ymax]], norm = LogNorm()); plt.colorbar(plthist[3])
        plt.title(month+"_RUN%i_CH%i"%(run,ch))
        plt.scatter(xpeak, ypeak, s=50,color='crimson')
        plt.axhline(ypeak, linestyle = '--',color='crimson', label= "Maximum value found at %.2f"%ypeak)
        plt.xlabel("Amp (ADC counts)", fontsize=12)
        plt.ylabel(charge[0]+"/"+charge[1], fontsize=12)
        plt.legend()
        if auto == True:
            lim_counts = counts[x_ind][y_ind]*5e-4

            ycounts_max = counts[x_ind][y_ind]
            ycounts_min = counts[x_ind][y_ind]
            xcounts_max = counts[x_ind][y_ind]
            xcounts_min = counts[x_ind][y_ind]

            for xind_min in np.linspace(x_ind,0):
                if xcounts_min > lim_counts:
                    xcounts_min = counts[int(xind_min)][y_ind]
                    amp_min = xedges[int(xind_min)]
                else:
                    break
            amp_min = xedges[int(xind_min)]

            for xind_max in range(x_ind, np.size(counts[:][y_ind])):
                if xcounts_max > lim_counts:
                    xcounts_max = counts[xind_max][y_ind]
                    amp_max = xedges[xind_max]
                else:
                    break
            amp_max = xedges[int(xind_max)]

            for yind_max in range(y_ind, np.size(counts[x_ind][:])):
                if ycounts_max > lim_counts:
                    ycounts_max = counts[x_ind][yind_max]
                    f90_max = yedges[yind_max]
                else:
                    break
            f90_max = yedges[int(yind_max)]

            for yind_min in np.linspace(y_ind,0):
                if ycounts_min > lim_counts:
                    ycounts_min = counts[x_ind][int(yind_min)]
                    f90_min = yedges[int(yind_min)]
                else:
                    break
            f90_min = yedges[int(yind_min)]
            
            tellme("Happy? Press any key to continue\nor select events yourself")
            plt.plot([amp_min,amp_max],[f90_min,f90_min],color="k")
            plt.plot([amp_min,amp_max],[f90_max,f90_max],color="k")
            plt.plot([amp_min,amp_min],[f90_min,f90_max],color="k")
            plt.plot([amp_max,amp_max],[f90_min,f90_max],color="k")
            
            test = plt.waitforbuttonpress()
            interactive = False
            if test == False:
                interactive = True
            
            if interactive == True:
                
                while True:
                    plt.clf()
                    plthist = plt.hist2d(np_amp,np_f90,200,[[xmin,xmax],[ymin,ymax]],norm=LogNorm())    
                    plt.xlabel("Amp (ADC counts)", fontsize=12)
                    plt.ylabel(charge[0]+"/"+charge[1], fontsize=12)  
                    plt.colorbar(plthist[3])
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
                    
                    plthist = plt.hist2d(np_amp,np_f90,200,[[xmin,xmax],[ymin,ymax]],norm=LogNorm())    
                    plt.title("Charge map for rate calculation and particle ID (RUN %i CH %i)"%(run,ch)) 
                    plt.xlabel("Amp (ADC counts)", fontsize=12)
                    plt.ylabel(charge[0]+"/"+charge[1], fontsize=12)  
                    plt.colorbar(plthist[3])

        else:
            #___ASK IN TERMINAL FOR F90 RANGE TO COUNT EVENTS___
            print("\nINPUT RANGE FOR EVENT COUNTING\n")
            #plt.axhline(.7)
            default = input("Use default cut (0.7)? (yes / no): ")
            
            if default in yes_string:
                amp_min = 0
                amp_max = np.max(np_amp)
                f90_min = 0.7
                f90_max = 1
            
            elif default in no_string:
                amp_min = float(input("Enter amp min (0.0-1e5): "))
                amp_max = float(input("Enter amp max (0.0-1e5): "))
                f90_min = float(input("Enter F90 min (0.0-1.0): "))
                f90_max = float(input("Enter F90 max (0.0-1.0): "))

            else:
                print("WRONG INPUT!")

    plt.clf()
    #plt.close("all")  

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

    return rate,nevents,dfout
