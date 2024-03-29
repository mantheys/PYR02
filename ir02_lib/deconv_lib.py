import numpy as np
import ROOT
import matplotlib.pyplot as plt
import scipy.fft

from numpy import convolve
from scipy.optimize import curve_fit
from ROOT import TF2, TH1D, TF1, TFile, TCanvas
from ROOT import gROOT
from scipy.special import erf

class my_wvf:
    smooth=True
    kernel_Done=False
    def __init__(self,f_type="",file_path="",item_path="",timebin=4e-9,normalize=False,trim=False,align=False,start=100,cut_i=0,cut_f=0,invert=False,item=np.zeros(1)):
        if f_type == "hist":
            self.wvf = tfile_hist2array(file_path,item_path)
        if f_type =="vector":
            self.wvf = tfile_vect2array(file_path,item_path)
        if f_type =="python":
            self.wvf = item
        if normalize and invert == False:
            self.wvf=self.wvf/max(self.wvf)
        if normalize and invert == True:
            self.wvf=self.wvf/max(-self.wvf)
        if trim == True:
            for i in range(len(self.wvf[np.argmax(self.wvf):])):
                if self.wvf[np.argmax(self.wvf)+i]<0<self.wvf[np.argmax(self.wvf)+i+1]:
                    self.wvf = self.wvf[:np.argmax(self.wvf)+i+1]
                    break
        if align == True and invert == False:
            self.wvf=np.roll(self.wvf,start-np.argmax(self.wvf))
        if align == True and invert == True:
            self.wvf=np.roll(self.wvf,start-np.argmax(-self.wvf))

        self.timebin = timebin
        self.wvf = self.wvf[cut_i:len(self.wvf)-cut_f]
        N=len(self.wvf)
        self.wvf_x = np.linspace(0,N*timebin,N)
        self.doFFT()
        
    def apply_smooth(self,alfa):
        if self.smooth==True:
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

def nonblank_lines(f):
    for l in f:
        line = l.rstrip()
        if line:
            yield line

def import_scint_prof(path,timebin,normalize,trim,align,start,cut_i,cut_f,invert):
    inSiPMName = path
    inSiPM = ROOT.TFile.Open(inSiPMName ,"READ")
    listkeys_SiPM = inSiPM.GetListOfKeys()
    # print(listkeys_SiPM[0].GetName())
    return my_wvf("vector",inSiPMName,listkeys_SiPM[0].GetName(),timebin,normalize,trim,align,start,cut_i,cut_f,invert)

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

def func(x, a, c):
    return a*np.exp(-x/c)

def pdf(x, m, sd, norm ="standard", n=2):
    A=1;mean=m;std=sd
    if norm=="standard":
        A=1/(std * np.sqrt(2 * np.pi))
    else:
        A=norm
    y_out = A*np.exp( - (x - mean)**n / (2 * std**n))
    
    return y_out

def signal_int(name,data,timebin,detector,int_type,th,i_range=10,f_range=1000,out=False):
    detector_list = ["SiPM","PMT","SC"]
    conv_factor = [250,50,1030]
    for det in range(len(detector_list)):
        if detector == detector_list[det]:
            factor = (16384.0/2.0)*conv_factor[det]

    integral = 0
    start_waveform = 0
    end_waveform = len(data)-1
    max_index = np.argmax(data)
    # th = data[max_index]*th

    if int_type == "ALL":
        start_waveform = 0
        end_waveform = len(data)-1
        integral =  1e12*timebin*np.sum(data)/factor

    if int_type == "BASEL":
        th = data[max_index]*th
        for i in range(len(data[max_index:])):
            if data[i+max_index] <= th:
                end_waveform = i+max_index
                break

        for j in range(len(data[:max_index])):
            if data[max_index-j] <= th:
                start_waveform = max_index-j
                break       
        
        integral =  1e12*timebin*np.sum(data[start_waveform:end_waveform])/factor

    if int_type == "MIXED":
        th = data[max_index]*th
        for j in range(len(data[:max_index])):
            if data[max_index-j] <= th:
                start_waveform = max_index-j
                break
        end_time = 2000
        for j in range(len(data[max_index:])):
            if np.mean(data[max_index+j-10:max_index+j]) <= th:
                end_time = max_index+j
                break
        if end_time > len(data)- max_index:
            end_time = len(data)- max_index                
        end_waveform = max_index + end_time

        integral =  1e12*timebin*np.sum(data[start_waveform:end_waveform])/factor    

    if int_type == "FIXED":
        start_waveform = max_index - 20
        end_waveform = max_index + 2000       
        integral =  1e12*timebin*np.sum(data[start_waveform:end_waveform])/factor

    if int_type == "RANGE":
        start_waveform = max_index - i_range
        end_waveform = max_index + f_range
        integral =  1e12*timebin*np.sum(data[start_waveform:end_waveform])/factor
    
    if int_type == "THRLD":
        for k in range(len(data)):
            if data[k]>=th:
                if integral == 0:
                   start_waveform = k 
                end_waveform = k
                integral = integral + 1e12*timebin*np.sum(data[k])/factor
    
    if int_type == "UNDER":
        for i in range(len(data[max_index:])):
            if data[i+max_index] < 0.:
                start_waveform = i+max_index
                break
        for j in range(len(data[start_waveform:])):
            if data[j+start_waveform] < 0. < data[j+start_waveform+1]:
                end_waveform = j+start_waveform+1
                break
        integral =  1e12*timebin*np.sum(data[start_waveform:end_waveform])/factor

    if out == True:
        print("\n%s integrated charge:\n %.4e pC for %s with type %s"%(name,integral,detector,int_type))

    return integral,end_waveform*timebin,start_waveform*timebin

def simple_2exp(x,t_fast,t_slow,amp_fast,amp_slow):
    return (amp_fast/t_fast)*np.exp(-x/t_fast)+(amp_slow/t_slow)*np.exp(-x/t_slow)

def conv_guess3(wvf,t_fast,t_int,t_slow,amp_fast,amp_int,amp_slow):
    resp = amp_fast*np.exp(-wvf.wvf_x/t_fast)+amp_int*np.exp(-wvf.wvf_x/t_int)+amp_slow*np.exp(-wvf.wvf_x/t_slow)
    conv = convolve(wvf.wvf,resp)
    return conv[:len(wvf.wvf_x)]/np.max(conv[:len(wvf.wvf_x)])

def conv_guess2(wvf,t_fast,t_slow,amp_fast,amp_slow):
    resp = (amp_fast/t_fast)*np.exp(-wvf.wvf_x/t_fast)+(amp_slow/t_slow)*np.exp(-wvf.wvf_x/t_slow)
    conv = convolve(wvf.wvf,resp)
    return conv[:len(wvf.wvf_x)]/np.max(conv[:len(wvf.wvf_x)])

def func(T,P,A,SIGMA,TAU,T0):
    # P = 5e-6
    return P+(2*A/TAU)*np.exp((SIGMA/(np.sqrt(2)*TAU))**2-(np.array(T)-T0)/TAU)*(1-erf((SIGMA**2-TAU*(np.array(T)-T0))/(np.sqrt(2)*SIGMA*TAU)))

def func2(x,TAU1,TAU2,A1,A2,SIGMA):
    return func(x,0,A1,SIGMA,TAU1,1e-7) + func(x,0,A2,SIGMA,TAU2,1e-7)

def conv_func2(wvf,TAU1,TAU2,A1,A2,SIGMA):
    # resp = func(wvf.wvf_x,0,A1,SIGMA,TAU1,1e-7) + func(wvf.wvf_x,0,A2,SIGMA,TAU2,1e-7)
    resp = func(wvf[0],0,A1,SIGMA,TAU1,1e-7) + func(wvf[0],0,A2,SIGMA,TAU2,1e-7)
    
    conv = convolve(wvf[1],resp)
    conv = conv/np.max(conv)
    wvf_max = np.argmax(wvf[1])
    conv_max = np.argmax(conv)
    return conv[conv_max-wvf_max:conv_max+len(wvf[1])-wvf_max]

def logconv_func2(wvf,TAU1,TAU2,A1,A2,SIGMA):
    # resp = func(wvf.wvf_x,0,A1,SIGMA,TAU1,1e-7) + func(wvf.wvf_x,0,A2,SIGMA,TAU2,1e-7)
    resp = np.log(func(wvf[0],1e-6,A1,SIGMA,TAU1,1e-7) + func(wvf[0],1e-6,A2,SIGMA,TAU2,1e-7))
    conv = convolve(wvf[1],resp)
    conv = conv/np.max(conv)
    wvf_max = np.argmax(wvf[1])
    conv_max = np.argmax(conv)
    return conv[conv_max-wvf_max:conv_max+len(wvf[1])-wvf_max]

def str2bool(v):
  return v.lower() in ("True", "true", "t", "1", "yes", "Yes")

def import_deconv_runs(path,debug = False):
    paths = ["","","",""]
    with open(path,'r') as f_in:
        for next_line in nonblank_lines(f_in):
            # while True:
                # next_line = file.readline()
            if next_line.startswith("detector:"):
                detector = next_line.split()[1]
            elif next_line.startswith("label:"):
                label_particle = next_line.split()[1]
            elif next_line.startswith("timebin:"):
                timebin = float(next_line.split()[1])
            elif next_line.startswith("int_st:"):
                int_st = next_line.split()[1]
            elif next_line.startswith("path:"):
                try:
                    base_dir = next_line.split()[1]
                except:
                    base_dir = "/"        
                paths[0]=base_dir
            elif next_line.startswith("path_alp:"):
                path_alp = next_line.split()[1]
                paths[1]=path_alp
                print("Signal file: %s"%(base_dir+path_alp))
            elif next_line.startswith("path_las:"):
                path_las = next_line.split()[1]
                paths[2]=path_las
            elif next_line.startswith("path_spe:"):
                path_spe = next_line.split()[1]
                paths[3]=path_spe
            elif next_line.startswith("filename:"):
                filename = next_line.split()[1]
            elif next_line.startswith("shift:"):
                shift = next_line.split()[1]
                shift = str2bool(shift)
            elif next_line.startswith("f_strength:"):
                f_stregth = int(next_line.split()[1])
            elif next_line.startswith("reverse:"):
                reverse = next_line.split()[1]
            elif not next_line:
                break
            elif next_line.startswith("f_array_cut:"):
                f_cut = int(next_line.split()[1])
            elif next_line.startswith("i_array_cut:"):
                i_cut = int(next_line.split()[1])
            
        if debug == True:
            print("Detector: %s"%detector)
            print("Particle Label: %s"%label_particle)        
            print("Timebin: %.2e"%timebin)
            print("Integration Type: %s"%int_st)
            print("Shiftting: %s"%shift)
            print("Filter Strength: %i"%f_stregth)
            print("Array final cut: %i"%f_cut)
            print("Array initial cut: %i"%i_cut)    
        # print(next_line.strip())
        # file.close()

    return detector,timebin,int_st,paths,filename,shift,f_stregth,reverse,f_cut,i_cut,label_particle

def low_pass_filter(array,timebin,cut_off,grade):
    fc = cut_off  # Cut-off frequency of the filter
    d_freq = 1/(timebin*len(array))
    w = fc/(d_freq/2) # Normalize the frequency
    b, a = array.butter(grade, w, 'low')
    output = array.filtfilt(b,a,array,method='gust')
    return output

def terminal_input():
    run_list = []; ch_list = []
    run_list_string = input("Enter run list (use space between numbers): ")
    ch_list_string = input("Enter channel list (use space between numbers): ")
    month = input("Enter month label: ")

    run_list_split = run_list_string.split(" ")
    ch_list_split = ch_list_string.split(" ")

    for i in range(np.size(run_list_split)):
        run_list.append(int(run_list_split[i])) 

    for i in range(np.size(ch_list_split)):
        ch_list.append(int(ch_list_split[i]))
    
    return run_list, ch_list, month

def tellme(title):
    plt.suptitle(title,x=10,y=10,fontsize=12)
    plt.draw()
