import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import ROOT
import scipy.fft

from scipy.optimize import curve_fit

from matplotlib.colors import LogNorm
from matplotlib.backend_bases import MouseButton
from ROOT import TF2, TH1D, TF1, TFile
#from matplotlib.patches import Rectangle

class my_wvf:
    smooth=False
    kernel_Done=False
    
    def __init__(self,f_type="",file_path="",item_path="",timebin=4e-9,normalize=False,item=np.zeros(1)):
        if f_type == "hist":
            self.wvf = tfile_hist2array(file_path,item_path)
        if f_type =="vector":
            self.wvf = tfile_vect2array(file_path,item_path)
        if f_type =="python":
            self.wvf = item
        if normalize:
            self.wvf=self.wvf/max(self.wvf)
        self.timebin=timebin
        N=len(self.wvf)
        self.wvf_x = x_time=np.linspace(0,N*timebin,N)
        self.doFFT()
        
    def apply_smooth(self,alfa):
        if self.smooth==False:
            self.wvf=expo_average(self.wvf,alfa)
            self.wvf=unweighted_average(self.wvf)
            self.smooth=True
            self.doFFT()
        else:
            print("Already smoothed")
    
    def doFFT(self):
        FFT=myrfft(self.wvf, self.timebin)
        self.wvf_F=FFT[0]
        self.wvf_F_x=FFT[1]
    
    def compute_kernel(self):
        if (not self.kernel_Done):
            print("code to compute the kernel here")
            
            self.kernel_Done=True
        else:
            print("Kernel already computed")
        
    def deconvolve_signal(self):
        #compute kernel in case is not done by this point.
        if not self.kernel_Done:
            self.compute_kernel()
        #deconvolve the signal (convolve with the kernel)
    
    def Resize(self,N):
        aux=self.wvf.copy()
        aux.resize(N)
        self.wvf=aux
        N=len(self.wvf)
        self.wvf_x = x_time=np.linspace(0,N*self.timebin,N)
        self.doFFT()
    
    def subs_baseline(self):
        self.wvf=self.wvf-1500 # Rodrigo tenía puesto 1500 y yo 375
        self.doFFT()
    
    def deconvolve(self,denominator):
        self.wvf_deco_F=self.wvf_F/denominator
        self.wvf_deco=scipy.fft.irfft(self.wvf_deco_F)

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

def charge_map(dfin,run,ch,charge,automatic=True,zoom=False):

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
    print(f'\nThe maximum count is {counts[x_ind][y_ind]:.0f} at index ({x_ind}, {y_ind})')
    print(f'Between x values {xedges[x_ind]} and {xedges[x_ind+1]}')
    print(f'and between y values {yedges[y_ind]} and {yedges[y_ind+1]}')

    fig1 = plt.subplot()
    plthist = plt.hist2d(np_amp,np_f90,200,[[xmin,xmax],[ymin,ymax]],norm=LogNorm())    
    plt.xlabel("Amp (ADC counts)", fontsize=12)
    plt.ylabel(charge[0]+"/"+charge[1], fontsize=12)  
    plt.colorbar(plthist[3])
    plt.scatter(xpeak, ypeak, s=50,color='crimson')

    if automatic == True:
        plt.ion()
        lim_counts = counts[x_ind][y_ind]*0.01
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
        amp_max = xedges[xind_max]

        for yind_max in range(y_ind, np.size(counts[x_ind][:])):
            if ycounts_max > lim_counts:
                ycounts_max = counts[x_ind][yind_max]
                f90_max = yedges[yind_max]
            else:
                break
        f90_max = yedges[yind_max]

        for yind_min in np.linspace(y_ind,0):
            if ycounts_min > lim_counts:
                ycounts_min = counts[x_ind][int(yind_min)]
                f90_min = yedges[int(yind_min)]
            else:
                break
        f90_min = yedges[int(yind_min)]
        
        tellme("Happy? Press any key to continue or click your own selection")
        plt.plot([amp_min,amp_max],[f90_min,f90_min],color="k")
        plt.plot([amp_min,amp_max],[f90_max,f90_max],color="k")
        plt.plot([amp_min,amp_min],[f90_min,f90_max],color="k")
        plt.plot([amp_max,amp_max],[f90_min,f90_max],color="k")
        
        test = plt.waitforbuttonpress()
        if test == False:
            interactive = True
        """
        # ATTEMPT OF INCREMETAL AUTOMATIC
        counter = 0
        while True:
            plt.ion()
            counter = counter +1
            amp_min, amp_max = xpeak-150*counter, xpeak+150*counter
            f90_min, f90_max = ypeak-0.05*counter, ypeak+0.05*counter
            rect = Rectangle((amp_min, f90_min), amp_max-amp_min, f90_max-f90_min,linewidth=1,edgecolor='crimson',facecolor='none')
            fig1.add_patch(rect)
            tellme("Happy with automatic selection? (press key for yes; mouse click for no)")
            
            if plt.waitforbuttonpress(-1):
                    break
        """
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
                """
                # DOES NOT SEEM TO WORK FOR SUCCESIVE ITERATIONS
                rect = Rectangle((amp_min, f90_min), amp_max-amp_min, f90_max-f90_min,linewidth=1,edgecolor='crimson',facecolor='none')
                fig1.add_patch(rect)
                """
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

def tfile_hist2array(tfile,hist_path):
    file=TFile( tfile, 'READ' )
    h=file.Get(hist_path)
    a=root_numpy.hist2array(h)
    return a

def tfile_vect2array(tfile,hist_path):
    file=TFile( tfile, 'READ' )
    item=file.Get(hist_path)
    Pitem=np.asarray(item)
    return Pitem

def myrfft(array, timebin):
    array_F = scipy.fft.rfft(array,len(array))
    xf_array = scipy.fft.rfftfreq(len(array), timebin)
    return array_F,xf_array

def expo_average(vector,alpha):
    v_averaged=np.zeros(len(vector))
    v_averaged[0]=vector[0]
    for i in range (len(vector)-1):
        v_averaged[i+1]=(1-alpha)*v_averaged[i]+alpha*vector[i+1]
    return v_averaged

def unweighted_average(vector):
    v_averaged=np.zeros(len(vector))
    v_averaged[0]=vector[0]
    v_averaged[-1]=vector[-1]

    for i in range (len(vector)-2):
        v_averaged[i+1]=(vector[i]+vector[i+1]+vector[i+2])/3
    return v_averaged

def plot_F(mps):
    for mp in mps:
        plt.plot(mp[1]/units[0],np.abs(mp[0]))

def func(x, a, c):
    return a*np.exp(-x/c)

def deconvolution_noise(sipm_name, laser_name, noise_name, deconv_time = 2400, initial_time = 450, end_time = 1000, ch_sipm = 0, day = ""):
    base_dir = 'DATA/DIC_21/AnalysisROOT/'
#     Load SiPM name
    inSiPMName = base_dir+sipm_name
    inSiPM = ROOT.TFile.Open(inSiPMName ,"READ")
    listkeys_SiPM = inSiPM.GetListOfKeys()
#     Load laser name
    inLaserName = base_dir+laser_name
    inLaser = ROOT.TFile.Open(inLaserName ,"READ")
    listkeys_Laser = inLaser.GetListOfKeys()
    # Load Noise name
    inNoiseName = base_dir+noise_name
    inNoise = ROOT.TFile.Open(inNoiseName ,"READ")
    listkeys_Noise = inNoise.GetListOfKeys()
    
    # Dump histograms
    if ch_sipm == 0:
        sipm = my_wvf("vector",inSiPMName,listkeys_SiPM[0].GetName())
        laser = my_wvf("vector",inLaserName,listkeys_Laser[0].GetName())
        noise = my_wvf("vector",inNoiseName,listkeys_Noise[0].GetName())
    else:
        sipm = my_wvf("vector",inSiPMName,listkeys_SiPM[1].GetName())
        laser = my_wvf("vector",inLaserName,listkeys_Laser[1].GetName())
        noise = my_wvf("vector",inNoiseName,listkeys_Noise[1].GetName())
    
    # Plot signals and apply smooth
    fig=plt.figure(figsize=(6,4), dpi= 150, facecolor='w', edgecolor='k')
    smooth=0.6;
    plt.plot(laser.wvf/max(laser.wvf))
    laser.apply_smooth(smooth)
    plt.plot(laser.wvf/max(laser.wvf))
    plt.plot(sipm.wvf/max(sipm.wvf))
    sipm.apply_smooth(smooth)
    plt.plot(sipm.wvf/max(sipm.wvf))
    plt.semilogy()
    plt.legend(["Laser Raw","Laser Raw+Smooth", "SiPM Raw", "SiPM Raw+Smooth"]); plt.xlabel("Bin number (1bin = 4ns)")
    plt.show()
    
    # Calculate Wiener filter
    wiener = abs(sipm.wvf_F)**2/(abs(sipm.wvf_F)**2+abs(noise.wvf_F)**2)
    wiener_laser = abs(laser.wvf_F)**2/(abs(laser.wvf_F)**2+abs(noise.wvf_F)**2)
    # Calculate deconvolved
    deconvolved = sipm.wvf_F/laser.wvf_F*wiener
#     deconvolved = sipm.wvf_F/laser.wvf_F*wiener*wiener_laser
    
    # Plot freqs
    fig=plt.figure(figsize=(6,4), dpi= 150, facecolor='w', edgecolor='k')
    plt.plot(sipm.wvf_F_x,abs(deconvolved))
    plt.plot(sipm.wvf_F_x,abs(sipm.wvf_F))
    plt.plot(sipm.wvf_F_x,abs(laser.wvf_F))
    plt.grid(which="both"); plt.semilogy(); plt.semilogx();
    plt.legend(["Deconvolved", "SiPM", "Laser"])
    plt.xlabel("Freq (Hz)")
    plt.show()
    
    # Probably tail is problematic
    print("Check deconvolved time!!")
    deconvolved_short=deconvolved[:deconv_time]
    deconvolved_time=scipy.fft.irfft(deconvolved_short)
    fig=plt.figure(figsize=(6,4), dpi= 90, facecolor='w', edgecolor='k')
    plt.plot(abs(deconvolved))
    plt.plot(abs(deconvolved_short))
    plt.semilogy(); plt.show()
    
    deco1=my_wvf("python",item=deconvolved_time)
    deco1_maxindex = np.where(deco1.wvf==max(deco1.wvf))[0][0] # Bin where deconvolved signal reaches max
    deco1_13us = np.where(deco1.wvf_x>1.3e-6)[0][0] # Bin where 1.3 us
    index_diff = deco1_maxindex-deco1_13us # Bin difference in order to place max in 1.3us
    ls = [] # Empty list for tau values
    
    deco1_mod = np.concatenate((deco1.wvf[index_diff:],deco1.wvf[:index_diff]), axis = None)
    fig=plt.figure(figsize=(6,4), dpi= 150, facecolor='w', edgecolor='k')
    plt.plot(deco1.wvf_x,deco1_mod)
    popt, pcov = curve_fit(func, deco1.wvf_x[initial_time:end_time],deco1_mod[initial_time:end_time], p0=(1, 1e-6)) # Fit between 1.8 and 4 us
    y_before = func(deco1.wvf_x[initial_time:end_time], *popt)
    print("T_slow (pre-smooth) =",popt[1])
    ls.append(popt[1]*1e6)
    
    deco1.apply_smooth(0.6)
    deco1_mod = np.concatenate((deco1.wvf[index_diff:],deco1.wvf[:index_diff]), axis = None)
    popt, pcov = curve_fit(func, deco1.wvf_x[initial_time:end_time],deco1_mod[initial_time:end_time], p0=(1, 1e-6))
    y_after = func(deco1.wvf_x[initial_time:end_time], *popt)
    print("T_slow (post-smooth) =",popt[1])
    ls.append(popt[1]*1e6)
    
    plt.plot(deco1.wvf_x,deco1_mod)
    plt.plot(deco1.wvf_x[initial_time:end_time],y_before)
    plt.plot(deco1.wvf_x[initial_time:end_time],y_after)
    plt.xlim([0,5e-6])
    plt.xlabel("Time (s)"); 
    plt.title("Scintillation Profile "+"("+day+")")
    plt.legend(["Pre-Smooth", "Post-Smooth",r'$\tau_{slow} =$'+str(round(ls[0],2)), r'$\tau_{slow} = $'+str(round(ls[1],2))])
    plt.semilogy()
    plt.show()