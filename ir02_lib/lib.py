import numpy as np
import ROOT
import scipy
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit
from scipy import fft, signal
from ROOT import TF2, TH1D, TF1, TFile

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
        self.wvf_x = x_time = np.linspace(0,N*timebin,N)
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

def deconvolution_noise(path, alpha_name, laser_name, noise_name, deconv_time = 2400, initial_time = 450, end_time = 1000, smooth = 0.6):
    base_dir = path; dic = {alpha_name, laser_name, noise_name}
    #Load SiPM name
    inAlphaName = base_dir+alpha_name
    inAlpha = ROOT.TFile.Open(inAlphaName ,"READ")
    listkeys_Alpha = inAlpha.GetListOfKeys()
    #Load laser name
    inLaserName = base_dir+laser_name
    inLaser = ROOT.TFile.Open(inLaserName ,"READ")
    listkeys_Laser = inLaser.GetListOfKeys()
    # Load Noise name
    inNoiseName = base_dir+noise_name
    inNoise = ROOT.TFile.Open(inNoiseName ,"READ")
    listkeys_Noise = inNoise.GetListOfKeys()
    
    alpha = my_wvf("vector",inAlphaName,listkeys_Alpha[0].GetName())
    laser = my_wvf("vector",inLaserName,listkeys_Laser[0].GetName())
    noise = my_wvf("vector",inNoiseName,listkeys_Noise[0].GetName())
    
    plt.plot(laser.wvf/max(laser.wvf))
    plt.plot(alpha.wvf/max(alpha.wvf))
    plt.legend(["Laser Raw", "Signal Raw"])
    plt.title("Laser vs. Alpha Average Waveforms")
    plt.xlabel("Bin number [1bin = 4ns]")
    plt.ylabel("Normalized Amplitude")
    plt.show()

    laser.apply_smooth(smooth)
    alpha.apply_smooth(smooth)
    plt.plot(laser.wvf/max(laser.wvf))
    plt.plot(alpha.wvf/max(alpha.wvf))
    #plt.semilogy()
    plt.legend(["Laser Smooth", "Signal Smooth"]); plt.xlabel("Bin number (1bin = 4ns)")
    plt.show()

    shift_sipm = np.roll(np.array(alpha.wvf),np.argmax(laser.wvf)-np.argmax(alpha.wvf))
    plt.plot(laser.wvf/max(laser.wvf))
    plt.plot(shift_sipm/max(shift_sipm))
    #plt.semilogy()
    plt.title("Laser vs. Alpha Average Waveforms (Aligned)")
    plt.legend(["Laser Smooth", "Signal Smooth"]); plt.xlabel("Bin number [1bin = 4ns]")
    plt.ylabel("Normalized Amplitude")
    plt.show()
    
    tau_slow = 1e6
    x = np.linspace(0, 5000, 5002)
    y = np.exp(-x/tau_slow)*1e-2
    shift_exp = np.roll(y,np.argmax(laser.wvf)-np.argmax(y))
    conv = scipy.signal.fftconvolve(np.array(laser.wvf/max(laser.wvf)), shift_exp, "same")
    plt.plot(laser.wvf/max(laser.wvf))
    plt.plot(x,shift_exp)
    plt.plot(x, conv)
    plt.show()

    # Calculate Wiener filter
    wiener = abs(alpha.wvf_F)**2/(abs(alpha.wvf_F)**2+abs(noise.wvf_F)**2)
    wiener_laser = abs(laser.wvf_F)**2/(abs(laser.wvf_F)**2+abs(noise.wvf_F)**2)
    # Calculate deconvolved
    deconvolved = alpha.wvf_F/laser.wvf_F
    #deconvolved = sipm.wvf_F/laser.wvf_F*wiener
    #deconvolved = sipm.wvf_F/laser.wvf_F*wiener*wiener_laser
    
    # Plot freqs
    #fig=plt.figure(figsize=(6,4), dpi= 150, facecolor='w', edgecolor='k')
    #plt.plot(sipm.wvf_F_x,abs(wiener),label = "Wiener", lw = 1)
    #plt.plot(sipm.wvf_F_x,abs(wiener_laser),label = "Wiener Laser", lw = 1)
    plt.plot(alpha.wvf_F_x,abs(deconvolved), label = "Deconvolution", lw = 1)
    plt.plot(alpha.wvf_F_x,abs(alpha.wvf_F), label = "Signal", lw = 1)
    plt.plot(alpha.wvf_F_x,abs(laser.wvf_F), label = "Laser", lw = 1)
    plt.grid(which="both"); plt.semilogy(); plt.semilogx();
    plt.legend()
    plt.xlabel("Freq (Hz)")
    plt.show()
    
    deconvolved_short=deconvolved[:deconv_time]
    deconvolved_time=scipy.fft.irfft(deconvolved_short)
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
    plt.plot(deco1.wvf_x,deco1_mod)
    
    popt, pcov = curve_fit(func, deco1.wvf_x[initial_time:end_time],deco1_mod[initial_time:end_time], p0=(0, 4e-6)) # Fit between 1.8 and 4 us
    y_before = func(deco1.wvf_x[initial_time:end_time], *popt)
    print("T_slow (pre-smooth) =",popt[1])
    ls.append(popt[1])
    
    deco1.apply_smooth(0.6)
    deco1_mod = np.concatenate((deco1.wvf[index_diff:],deco1.wvf[:index_diff]), axis = None)
    
    popt, pcov = curve_fit(func, deco1.wvf_x[initial_time:end_time],deco1_mod[initial_time:end_time], p0=(0, 4e-6))
    y_after = func(deco1.wvf_x[initial_time:end_time], *popt)
    print("T_slow (post-smooth) =",popt[1])
    ls.append(popt[1])
    
    plt.plot(deco1.wvf_x,deco1_mod)
    plt.plot(deco1.wvf_x[initial_time:end_time],y_before)
    plt.plot(deco1.wvf_x[initial_time:end_time],y_after)
    plt.xlim([0,5e-6])
    plt.xlabel("Time [s]");plt.ylabel("Amplitude [a.u.]");
    plt.title("Scintillation Profile")
    plt.legend(["Pre-Smooth", "Post-Smooth",r'$\tau_{slow} = $'+str(round(ls[0],4)), r'$\tau_{slow} = $'+str(round(ls[1],4))])
    #plt.semilogy()
    plt.show()