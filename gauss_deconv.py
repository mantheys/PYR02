import ROOT
import scipy
import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit
from ir02_lib.deconv_lib import import_scint_prof,pdf,func,my_wvf,signal_int,import_deconv_runs

########################################################################
#______________________IMPORT_SCINTILLATION_FILES______________________#
########################################################################

# WELCOME USER AND IMPORT CONFIG FILE WITH DATA PATHS AND DECONVOLUTION PARAMETERS
print("\n### WELCOME TO THE DECONVOLUTION STUDIES ###")
detector,timebin,int_st,paths,filename,filter_strenght,reverse = import_deconv_runs("deconvolution_input/FEB_2_SC_DAY1_OV1.txt",debug = True)
check = True; autozoom = False; logy = False; shift = True; fit = False; thrld = 1e-2; pro_abs = False; pro_rodrigo = True

# SELECT THE RIGHT TIME BIN FOR ACCURATE PLOT REPRESENATION
s_start = 500; init=0; fnal=500; trm = False
alp = import_scint_prof(paths[0]+paths[1],timebin,normalize=False,trim=trm,align=shift,start=s_start,cut_i=init,cut_f=fnal)
las = import_scint_prof(paths[0]+paths[2],timebin,normalize=False,trim=trm,align=shift,start=s_start,cut_i=init,cut_f=fnal)
"""
for i in range(len(alp.wvf)):
    if len(las.wvf)>len(alp.wvf):
        las = import_scint_prof(paths[0]+paths[2],timebin,normalize=False,trim=True,align=False,start=100,cut_i=init,cut_f=fnal)
        init = init+1
        print("Array len of alpha (%i) vs laser (%i)"%(len(alp.wvf),len(las.wvf)))
    else:
        break
"""
spe = import_scint_prof(paths[0]+paths[3],timebin,normalize=False,trim=False,align=shift,start=100,cut_i=0,cut_f=0)

# CALCULATE INTEGRAL OF RAW DATA
labels = ["Alpha","Laser","SPE"]

# CHECK INTEGRATION METHODS FOR RAW ALPHA SIGNALq
alp_int_all,f_alp_all,i_alp_all       = signal_int(labels[0],alp.wvf,timebin,detector,"ALL"  ,th = thrld,out = True)
alp_int_basel,f_alp_basel,i_alp_basel = signal_int(labels[0],alp.wvf,timebin,detector,"BASEL",th = thrld,out = True)
alp_int_thrld,f_alp_thrld,i_alp_thrld = signal_int(labels[0],alp.wvf,timebin,detector,"THRLD",th = thrld,out = True)

signal_int(labels[1],las.wvf,timebin,detector,"BASEL",th = thrld,out=True)
try:
    signal_int(labels[2],spe.wvf,timebin,detector,int_st)
except:
    print("\nERROR: %s could not be properly integrated with type: %s"%(labels[2],int_st))

# FIRST PLOT TO VISUALIZE DATA (VERTICAL LINES SHOW INTEGRATION LIMITS OF ALPHA SIGNAL)
if check == True:
    plt.plot(alp.wvf_x,alp.wvf,label="Raw Alpha signal")
    plt.plot(las.wvf_x,las.wvf,label="Raw Laser signal")
    plt.plot(spe.wvf_x,spe.wvf,label="Raw SPE signal")
    plt.xlabel("Time in [s]");plt.ylabel("Amplitude in ADC counts")
    plt.axhline(0, ls="--", color = "k",alpha=0.25)
    plt.axvline(f_alp_all,  color = "k",ls = ":");plt.axvline(i_alp_all  ,color = "k",ls = ":",label = "Alpha Integration (ALL) Limits")
    plt.axvline(f_alp_basel,color = "g",ls = ":");plt.axvline(i_alp_basel,color = "g",ls = ":",label = "Alpha Integration (BASEL) Limits")
    plt.axvline(f_alp_thrld,color = "b",ls = ":");plt.axvline(i_alp_thrld,color = "b",ls = ":",label = "Alpha Integration (THRLD) Limits")
    
    if autozoom == True:
        plt.xlim(timebin*(np.argmax(alp.wvf)-100),timebin*(np.argmax(alp.wvf)+1000));plt.ylim()
    
    plt.legend();plt.show()

########################################################################
#________________________ANALYSE_FILTER_EFFECTS________________________#
########################################################################

# GAUSS FUNCTION IS CREATED IN FOURIER SPACE WITH AMP = 1 TO CONSERVE CHARGE
gauss_f = pdf(np.arange(len(alp.wvf_F_x)), m = 0, sd = filter_strenght, norm = 1, n = 2.)
if pro_rodrigo == True:
    gauss_f[0]=0
gauss = scipy.fft.irfft(gauss_f)
gauss = np.roll(gauss,int(len(gauss)/2))

# CHECK GAUSS BELL PROPERTIES
if check == True:
    gauss_int,f_gauss,i_gauss = signal_int("Gauss",gauss,timebin,detector,"BASEL",th = np.max(gauss)*1e-3,out = True)
    signal_int("Gauss",gauss,timebin,detector,"ALL",th = 1e-3,out = True)
    print("\nGaussian FWHM: %.2e [s]"%(timebin*2*np.sqrt(2*np.log(2))*filter_strenght))
    plt.xlabel("Time in [s]");plt.ylabel("Amplitude in a.u.")
    plt.axvline(f_gauss,color = "k", ls = ":");plt.axvline(i_gauss,color = "k", ls = ":")
    plt.plot(alp.wvf_x,gauss)
    
    if autozoom == True:
        plt.xlim(i_gauss,f_gauss)
    
    plt.show()

pr_signal = alp.wvf_F*(gauss_f)
amp_test=scipy.fft.irfft(pr_signal)
signal_int("Filtered Alpha",amp_test,timebin,detector,"ALL",out = True)

if check == True:
    plt.plot(alp.wvf_x,alp.wvf ,label="Raw Alpha signal")
    plt.plot(alp.wvf_x,amp_test,label="Filtered Alpha Signal")
    plt.xlabel("Time in [s]");plt.ylabel("Amplitude in ADC counts")
    
    if logy == True:
        plt.semilogy()
    
    if autozoom == True:
        plt.xlim(alp.wvf_x[np.argmax(alp.wvf)-200],alp.wvf_x[np.argmax(alp.wvf)+200]);plt.ylim(1e-3,2*np.max(alp.wvf))
    
    plt.legend();plt.show()

########################################################################
#_________EVALUATE DECONVOLUTION AND FILTER IN FOURIER SPACE___________#
########################################################################

if check == True:
    # plt.plot(alp.wvf_F_x,abs(gauss_f), label = "Gauss filter")
    plt.plot(las.wvf_F_x,abs(las.wvf_F), label = "Laser freq.")
    plt.plot(las.wvf_F_x,abs(las.wvf_F/np.max(las.wvf_F)), label = "Norm Laser freq.")
    plt.plot(alp.wvf_F_x,abs(alp.wvf_F), label = "Signal freq.")
    plt.plot(alp.wvf_F_x,abs(pr_signal), label = "Filtered signal freq.")
    plt.plot(alp.wvf_F_x,abs(pr_signal/(las.wvf_F/np.max(las.wvf_F))), label = "Deconvolved signal (norm laser)")
    plt.plot(alp.wvf_F_x,abs(pr_signal/las.wvf_F), label = "Deconvolved signal")
    plt.xlabel("Frequency in [Hz]");plt.ylabel("Amplitude in ADC counts")
    plt.semilogy();plt.semilogx()
    
    if autozoom == True:
        plt.ylim(1e-6,np.max(abs(alp.wvf_F))*10)
    
    plt.legend();plt.show()

# AFTER INVERSE FOURIER TRANSFORM YOU ARE REQUIRED TO SHIFT THE SIGNAL ARRAY TO RECONSTRUCT THE SHAPE OF A PULSE
dec=scipy.fft.irfft(pr_signal/(las.wvf_F/np.max(las.wvf_F)))
# dec=scipy.fft.irfft(pr_signal/las.wvf_F)
if reverse == "True":
    dec = dec[::-1]
dec = np.roll(dec,s_start);ref = np.argmax(dec)
if pro_abs == True:
    dec = np.abs(dec)

if check == True:
    conv_int_all,f_conv_all,i_conv_all = signal_int("Deconvolution",dec,timebin,detector,"ALL",out = True)
    # print("\nFull Deconvolution Integral: %.2e"%conv_int)
   
conv_int_range,f_conv_range,i_conv_range = signal_int("Deconvolution",dec,timebin,detector,"RANGE",th = thrld,i_range = np.argmax(alp.wvf)-int(i_alp_basel/timebin),f_range = int(f_alp_basel/timebin)-np.argmax(alp.wvf))
conv_int_thrld,f_conv_thrld,i_conv_thrld = signal_int("Deconvolution",dec,timebin,detector,"THRLD",th = thrld, out= True)
conv_int_ithrld,f_conv_ithrld,i_conv_ithrld = signal_int("Deconvolution",dec,timebin,detector,"I_THRLD",th = 0, out= True)
conv_int_basel,f_conv_basel,i_conv_basel = signal_int("Deconvolution",dec,timebin,detector,"BASEL",th = thrld,out= True)

try:
    pro_alp_int,f_pro_alp,i_pro_alp = signal_int("Filtered Alpha",amp_test,timebin,detector,"BASEL",th =thrld,out = True)
    print("\nCharge increase after filter application:\n %.2f%%\n"%(100*(pro_alp_int/alp_int_basel-1)))
except:
    print("\nERROR: %s could not be properly integrated with type: %s"%("Filtered Alpha","BASEL"))

########################################################################
#________________________PLOT_FINAL_RESULT_____________________________#
########################################################################

# SELECT THE INTERVAL LENGHT FOR FIT
if fit == True:
    roll = 400; ifit = 100; ffit = 400 #(fit starts at bin rollber and ends at bin rollber+ffit)
    popt, pcov = curve_fit(func, np.arange(0,ffit)*timebin, dec[ref+ifit:ref+ifit+ffit], bounds =([0,1e-7],[1e3,1e-5]))
    perr = np.sqrt(np.diag(pcov))
    plt.plot(np.arange(ref+ifit,ref+ifit+ffit)*timebin, func(np.arange(0,ffit)*timebin, *popt), label='Fit: tau_slow = %.2e s'%popt[1])

# PRINT CHARGE INFORMATION IN RELATION TO RAW DATA
print("\nCharge increase after deconvolution:\n %.2f%% with integration type 'ALL'\n"%(100*(conv_int_all/alp_int_all-1)))
print("\nCharge increase after deconvolution:\n %.2f%% with integration type 'BASEL'\n"%(100*(conv_int_basel/alp_int_basel-1)))
print("\nCharge increase after deconvolution:\n %.2f%% with integration type 'RANGE'\n"%(100*(conv_int_range/alp_int_basel-1)))
print("\nCharge increase after deconvolution:\n %.2f%% with integration type 'THRLD'\n"%(100*(conv_int_thrld/alp_int_thrld-1)))
print("\nCharge increase after deconvolution:\n %.2f%% with integration type 'BASEL' substracting 'I_THRLD'\n"%(100*((conv_int_basel+conv_int_ithrld)/alp_int_basel-1)))

# PLOT DECONVOLUTED WAVEFORM AND RAW DATA
plt.plot(alp.wvf_x,alp.wvf,label = "Raw Alpha Signal")
plt.plot(np.arange(len(dec))*timebin,dec,label = "Deconvolution with Gauss filter")
plt.axvline(f_conv_basel,color = "g", ls = ":");plt.axvline(i_conv_basel,color = "g", ls = ":",label = "(BASEL) Alpha Integration Limits")
plt.axvline(f_conv_range,color = "r", ls = ":");plt.axvline(i_conv_range,color = "r", ls = ":",label = "(RANGE) Alpha Integration Limits")
plt.axvline(f_conv_thrld,color = "b", ls = ":");plt.axvline(i_conv_thrld,color = "b", ls = ":",label = "(THRLD) Alpha Integration Limits")
plt.xlabel("Time in [s]");plt.ylabel("Amplitude in ADC counts")

if logy == True:
    plt.grid(); plt.semilogy()

if autozoom == True:
    plt.xlim([0,10e-6]);plt.ylim(1e-3,2*np.max(dec))

plt.legend(); plt.show()

# EXPORT DECONVOLUTED WAVEFORM TO ROOT FORMAT
f = ROOT.TFile("deconvolution_output/"+filename+".root", "RECREATE")
h = ROOT.TH1D(filename, filename,len(dec),0,len(dec)*4e-9)
for i in range(len(dec)):
    h.SetBinContent(i,dec[i]) # Use the following lines to feed the Fill method of the histograms in order to fill
    h.SetBinError(i,dec[i]*1e-2)
h.Write() # Here write the three histograms on the file and close the file
print(f.GetListOfKeys()[0].GetName())
f.Close()