import numpy as np
import ROOT
import scipy
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit
from scipy.signal import convolve
from scipy import fft, signal
from ROOT import TF2, TH1D, TF1, TFile

class my_wvf:
    smooth=False
    kernel_Done=False
    
    def __init__(self,f_type="",file_path="",item_path="",timebin=4e-9,wvfalize=True,item=np.zeros(1)):
        if f_type == "hist":
            self.wvf = tfile_hist2array(file_path,item_path)
        if f_type =="vector":
            self.wvf = tfile_vect2array(file_path,item_path)
        if f_type =="python":
            self.wvf = item
        if wvfalize:
            self.wvf=self.wvf/max(self.wvf)
        self.timebin=timebin
        N=len(self.wvf)
        self.wvf_x = np.linspace(0,N*timebin,N)
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
        self.wvf_x = np.linspace(0,N*self.timebin,N)
        self.doFFT()
    
    def subs_baseline(self):
        self.wvf=self.wvf-1500 # Rodrigo tenía puesto 1500 y yo 375
        self.doFFT()
    
    def deconvolve(self,denominator):
        self.wvf_deco_F=self.wvf_F/denominator
        self.wvf_deco=scipy.fft.irfft(self.wvf_deco_F)

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

def expo_average(vector,Raw):
    v_averaged=np.zeros(len(vector))
    v_averaged[0]=vector[0]
    for i in range (len(vector)-1):
        v_averaged[i+1]=(1-Raw)*v_averaged[i]+Raw*vector[i+1]
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

def conv_guess3(wvf,t_fast,t_int,t_slow,amp_fast,amp_int,amp_slow):
    resp = amp_fast*np.exp(-wvf.wvf_x/t_fast)+amp_int*np.exp(-wvf.wvf_x/t_int)+amp_slow*np.exp(-wvf.wvf_x/t_slow)
    conv = convolve(wvf.wvf,resp)
    return conv[:len(wvf.wvf_x)]/np.max(conv[:len(wvf.wvf_x)])

def conv_guess2(wvf,t_fast,t_slow,amp_fast,amp_slow):
    resp = amp_fast*np.exp(-wvf.wvf_x/t_fast)+amp_slow*np.exp(-wvf.wvf_x/t_slow)
    conv = convolve(wvf.wvf,resp)
    return conv[:len(wvf.wvf_x)]/np.max(conv[:len(wvf.wvf_x)])

def deconvolution_noise(path, raw_signal, pure_signal, noise_signal, deconv_time = 2400, initial_time = 450, end_time = 1000, smooth = 0.6):
    base_dir = path; dic = {raw_signal, pure_signal, noise_signal}
    
    #Load Raw name
    inRawName = base_dir+raw_signal
    inRaw = ROOT.TFile.Open(inRawName ,"READ")
    listkeys_Raw = inRaw.GetListOfKeys()
    #Load Pure name
    inPureName = base_dir+pure_signal
    inPure = ROOT.TFile.Open(inPureName ,"READ")
    listkeys_Pure = inPure.GetListOfKeys()
    # Load Noise name
    inNoiseName = base_dir+noise_signal
    inNoise = ROOT.TFile.Open(inNoiseName ,"READ")
    listkeys_Noise = inNoise.GetListOfKeys()
    
    Raw = my_wvf("vector",inRawName,listkeys_Raw[0].GetName())
    Pure = my_wvf("vector",inPureName,listkeys_Pure[0].GetName())
    noise = my_wvf("vector",inNoiseName,listkeys_Noise[0].GetName())
    
    plt.plot(Pure.wvf/max(Pure.wvf))
    plt.plot(Raw.wvf/max(Raw.wvf))
    plt.legend(["Pure Raw", "Signal Raw"])
    plt.title("Pure vs. Raw Average Waveforms")
    plt.xlabel("Bin number [1bin = 4ns]")
    plt.ylabel("wvfalized Amplitude")
    plt.show()

    Pure.apply_smooth(smooth)
    Raw.apply_smooth(smooth)
    plt.plot(Pure.wvf/max(Pure.wvf))
    plt.plot(Raw.wvf/max(Raw.wvf))
    #plt.semilogy()
    plt.legend(["Pure Smooth", "Signal Smooth"]); plt.xlabel("Bin number (1bin = 4ns)")
    plt.show()

    shift_sipm = np.roll(np.array(Raw.wvf),np.argmax(Pure.wvf)-np.argmax(Raw.wvf))
    plt.plot(Pure.wvf/max(Pure.wvf))
    plt.plot(shift_sipm/max(shift_sipm))
    #plt.semilogy()
    plt.title("Pure vs. Raw Average Waveforms (Aligned)")
    plt.legend(["Pure Smooth", "Signal Smooth"]); plt.xlabel("Bin number [1bin = 4ns]")
    plt.ylabel("wvfalized Amplitude")
    plt.show()

    # Calculate Wiener filter
    # wiener = abs(Raw.wvf_F)**2/(abs(Raw.wvf_F)**2+abs(noise.wvf_F)**2)
    # wiener_Pure = abs(Pure.wvf_F)**2/(abs(Pure.wvf_F)**2+abs(noise.wvf_F)**2)
    # Calculate deconvolved

    if len(Raw.wvf) == len(Pure.wvf):
        deconvolved = Raw.wvf_F/Pure.wvf_F
        #deconvolved = Raw.wvf_F/Pure.wvf_F*wiener
        #deconvolved = Raw.wvf_F/Pure.wvf_F*wiener*wiener_Pure
    """
    else:
        trunc_Raw_F = myrfft(np.array(Raw.wvf)[:len(Pure.wvf)],4e-9)
        print(trunc_Raw_F)
        deconvolved = np.array(trunc_Raw_F/Pure.wvf_F)
        print(deconvolved)
    """
    # Plot freqs
    #fig=plt.figure(figsize=(6,4), dpi= 150, facecolor='w', edgecolor='k')
    #plt.plot(sipm.wvf_F_x,abs(wiener),label = "Wiener", lw = 1)
    #plt.plot(sipm.wvf_F_x,abs(wiener_Pure),label = "Wiener Pure", lw = 1)
    plt.plot(abs(deconvolved), label = "Deconvolution", lw = 1)
    plt.plot(abs(Raw.wvf_F), label = "Signal", lw = 1)
    plt.plot(abs(Pure.wvf_F), label = "Pure", lw = 1)
    plt.grid(which="both"); plt.semilogy(); plt.semilogx()
    plt.legend()
    plt.xlabel("Freq (Hz)")
    plt.show()
    
    deconvolved_short=deconvolved[:deconv_time]
    deconvolved_time=scipy.fft.irfft(deconvolved_short)
    print(deconv_time)
    """
    # Probably tail is problematic
    print("Check deconvolved time!!")
    fig=plt.figure(figsize=(6,4), dpi= 90, facecolor='w', edgecolor='k')
    plt.plot(abs(deconvolved))
    plt.plot(abs(deconvolved_short))
    plt.semilogy(); plt.show()
    """

    deco1=my_wvf("python",item=deconvolved_time)
    deco1_maxindex = np.where(deco1.wvf==max(deco1.wvf))[0][0]  # Bin where deconvolved signal reaches max
    deco1_13us = np.where(deco1.wvf_x>1.3e-6)[0][0]             # Bin where 1.3 us
    index_diff = deco1_maxindex-deco1_13us                      # Bin difference in order to place max in 1.3us
    ls = []                                                     # Empty list for tau values
    
    deco1_mod = np.concatenate((deco1.wvf[index_diff:],deco1.wvf[:index_diff]), axis = None)
    plt.plot(deco1.wvf_x,deco1_mod, label = "Pre-Smooth")
    
    popt, pcov = curve_fit(func, deco1.wvf_x[initial_time:end_time],deco1_mod[initial_time:end_time], p0=(0, 4e-4)) # Fit between 1.8 and 4 us
    y_before = func(deco1.wvf_x[initial_time:end_time], *popt)
    print(deco1.wvf_x)
    ls.append(popt[1])
    plt.plot(deco1.wvf_x[initial_time:end_time],y_before, label = r'$\tau_{slow} = $'+str(round(ls[0],4)))
    print("T_slow (pre-smooth) =",popt[1])
    
    deco1.apply_smooth(0.6)
    
    deco1_mod = np.concatenate((deco1.wvf[index_diff:],deco1.wvf[:index_diff]), axis = None)
    popt, pcov = curve_fit(func, deco1.wvf_x[initial_time:end_time],deco1_mod[initial_time:end_time], p0=(0, 4e-4))
    y_after = func(deco1.wvf_x[initial_time:end_time], *popt)
    ls.append(popt[1])
    print(y_after)
    plt.plot(deco1.wvf_x[initial_time:end_time],y_after, label = r'$\tau_{slow} = $'+str(round(ls[1],4)))
    print("T_slow (post-smooth) =",popt[1])
    
    
    plt.plot(deco1.wvf_x,deco1_mod, label = "Post-Smooth")

    plt.ylim([-2,8]); plt.xlim([0,5e-6])
    plt.xlabel("Time [s]");plt.ylabel("Amplitude [a.u.]");
    plt.title("Scintillation Profile")
    plt.legend()
    # plt.semilogy()
    plt.show()

def convolution(path, signals, labels):
    
    wvfs = []
    for s in signals:
        #Load Raw name
        inRawName = path+s
        inRaw = ROOT.TFile.Open(inRawName ,"READ")
        listkeys_Raw = inRaw.GetListOfKeys()
        wvf = my_wvf("vector",inRawName,listkeys_Raw[0].GetName())
        # print(listkeys_Raw[0].GetName())
        wvfs.append(wvf)

    shift_wvfs = [wvfs[0].wvf]
    for n in range(1,len(wvfs)):
        if np.argmax(wvfs[n].wvf)-np.argmax(wvfs[0].wvf) > 0:
            shift = np.roll(np.array(wvfs[n].wvf),np.argmax(wvfs[n].wvf)-np.argmax(wvfs[0].wvf))
        else:
            shift = np.roll(np.array(wvfs[n].wvf),np.argmax(wvfs[0].wvf)-np.argmax(wvfs[n].wvf))
        shift_wvfs.append(shift)

    # print(wvfs)
    # print(shift_wvfs)

    t_fast = 1e-9
    t_int = 1e-7    
    t_slow = 1e-6

    n_label = 0
    for shift in shift_wvfs:
        plt.plot(wvfs[0].wvf_x,shift,label = labels[n_label])
        n_label = n_label+1

    # plt.plot(my_guess/np.max(my_guess), label = "Guess")
    plt.xlabel("Time in [s]"); plt.ylabel("Normalized Amplitude")
    plt.semilogy()
    plt.xlim(0,10e-6); plt.ylim(1e-6, 1.2)
    # plt.show()
    
    fit_initials = (t_fast,t_int,t_slow,0.999,0.001)
    fit_limits = ([1e-10,1e-12,1e-8,1e-4,1e-4],[1e-8,1e-5,1e-5,1,0.5])

    popt, pcov = curve_fit(conv_guess,wvfs[0],shift_wvfs[1], p0 = fit_initials, bounds = fit_limits)
    conv = conv_guess(wvfs[0],*popt)
    print("Tau_fast = %.2e [s] - Amp = %.2e"%(popt[0],popt[3]))
    print("Tau_int = %.2e [s] - Amp = %.2e"%(popt[1],1-popt[3]-popt[4]))
    print("Tau_slow = %.2e [s]- Amp = %.2e"%(popt[2],popt[4]))
    plt.plot(wvfs[0].wvf_x,conv,label = "Fitted Convolution")
    plt.legend()
    plt.show()
    """
    
    # plt.figure(dpi=300)
    plt.plot(scnt_prof)
    plt.semilogy(); #plt.semilogx()
    plt.ylim(1e-10,2)
    plt.show()
    
    plt.plot(x_Pure,Pure.wvf/max(Pure.wvf), label = "SPE") 
    plt.plot(x_Raw,Raw.wvf/max(Raw.wvf), label = "Raw")
    plt.plot(my_guess /max(my_guess ), label = "guess")
    plt.title("Pure vs. Raw Average Waveforms")
    # plt.xlabel("Bin number [1bin = 4ns]")
    plt.ylabel("wvfalized Amplitude")
    # plt.xlim(0,5000)
    plt.legend()
    plt.show()
    

    vline1 = 1900
    vline2 = 2900
    plt.axhline(1e-3,color = "grey",ls= ":")
    # plt.axvline(vline1,color = "grey" ,ls= ":",label = "V_line_1 = %i [ns]"%vline1 )
    # plt.axvline(vline2,color = "grey", ls= ":",label = "V_line_2 = %i [ns]"%vline2 )

    # plt.plot(x_Pure2,Pure2.wvf/max(Pure2.wvf), label = "SPE_Sergio")
    # plt.plot(x_Pure2,shift_Pure2/max(shift_Pure2), label = "Laser")
    # plt.plot(x_Pure2,Pure2.wvf, label = "Laser")
    # plt.plot(x_Pure2,shift_Pure2, label = "SPE_Sergio")

    plt.plot(x_Pure1,Pure1.wvf/max(Pure1.wvf), label = "SPE")
    # plt.plot(x_Pure1,Pure1.wvf, label = "Alpha")

    plt.plot(x_Raw,shift_Raw1/max(shift_Raw1), label = "Alpha")
    # plt.plot(x_Raw,Raw.wvf/max(Raw.wvf), label = "SPE OV = 2V")
    # plt.plot(x_Raw,Raw.wvf, label = "SPE")
    # plt.plot(x_Raw,shift_Raw1, label = "Alpha Signal")

    plt.plot(x_Guess,my_guess /max(my_guess ), label = "Guess")
    plt.plot(test,label = "guess")
    # plt.semilogy()
    plt.title("SPE vs. Alpha Signal Average Waveforms")
    plt.xlabel("Time in [ns]")
    plt.ylabel("Amplitude ADC")
    plt.xlim(2000,5000)
    plt.legend()
    plt.show()
    """

