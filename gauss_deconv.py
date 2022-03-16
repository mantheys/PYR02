import ROOT
import sys
import scipy
import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit
from ir02_lib.deconv_lib import import_scint_prof,pdf,func,my_wvf,signal_int,import_deconv_runs

########################################################################
#______________________IMPORT_SCINTILLATION_FILES______________________#
########################################################################

print("\n### WELCOME TO THE DECONVOLUTION STUDIES ###")

detector,timebin,int_st,paths,filename,filter_strenght,reverse = import_deconv_runs("deconvolution_input/FEB_2_SiPM_DAY2_OV1.txt",debug = True)

# SELECT THE RIGHT TIME BIN FOR ACCURATE PLOT REPRESENATION
check = True; autozoom = True; logy = True; norm = False; shift = True; fit = False
s_start = 500; s_cut = 500; thrld = 1e-6
alp = import_scint_prof(paths[0]+paths[1],timebin,normalize=norm,align=shift,start=s_start,cut=s_cut)
las = import_scint_prof(paths[0]+paths[2],timebin,normalize=norm,align=shift,start=s_start,cut=s_cut)
spe = import_scint_prof(paths[0]+paths[3],timebin,normalize=norm,align=False,start=0,cut=0)

roll = 400; ifit = 100; ffit = 400 #(fit starts at bin rollber and ends at bin rollber+ffit)

# CALCULATE INTEGRAL
labels = ["Alpha","Laser","SPE"]
alp_int,f_alp,i_alp = signal_int(labels[0],alp.wvf,timebin,detector,int_st,th = thrld,out = True)
signal_int(labels[0],alp.wvf,timebin,detector,"ALL",th = thrld,out = True)
las_int,f_las,i_las = signal_int(labels[1],las.wvf,timebin,detector,int_st,th = thrld)
try:
    spe_int,f_spe,i_spe = signal_int(labels[2],spe.wvf,timebin,detector,int_st)
except:
    print("\nERROR: %s could not be properly integrated with type: %s"%(labels[2],int_st))

# print("\nNumber of SPE in Laser pulse:\n %.2f SPE"%(las_int/spe_int))

# FIRST PLOT AND INTEGRAL CALCULATION
if check == True:
    plt.plot(alp.wvf_x,alp.wvf,label="Raw Alpha signal")
    plt.plot(las.wvf_x,las.wvf,label="Raw Laser signal")
    plt.plot(spe.wvf_x,spe.wvf,label="Raw SPE signal")
    if autozoom == True:
        plt.xlim(timebin*(np.argmax(alp.wvf)-100),timebin*(np.argmax(alp.wvf)+1000));plt.ylim()
    plt.xlabel("Time in [s]");plt.ylabel("Amplitude in ADC counts")
    plt.axvline(f_alp,color = "k", ls = ":");plt.axvline(i_alp,color = "k", ls = ":",label = "Alpha Integration Limits")
    plt.legend();plt.show()

########################################################################
#____________________DECONVOLUTION_AND_FIT_PARAMETERS__________________#
########################################################################

# SELECT A FILTER STRENGHT THAT ALLOWS YOU TO SEE THE SLOPE OF THE SIGNAL
# THE CORRECPONDING NORM GAUSS FILTER IS DEFINED AND MULTIPLIED IN FOURIER SPACE

gauss_f = pdf(np.arange(len(alp.wvf_F_x)), m = 0, sd = filter_strenght, norm = 1, n = 2.)#;gauss_f[0]=0
gauss = scipy.fft.irfft(gauss_f)
gauss = np.roll(gauss,int(len(gauss)/2))
pr_signal = alp.wvf_F*(gauss_f)

if check == True:
    gauss_int,f_gauss,i_gauss = signal_int("Gauss",gauss,timebin,detector,"BASEL",th = np.max(gauss)*1e-3)
    signal_int("Gauss",gauss,timebin,detector,"ALL",th = 1e-3,out = True)
    print("\nGaussian Width: %.2e [s]"%(f_gauss-i_gauss))
    plt.xlabel("Time in [s]");plt.ylabel("Amplitude in a.u.")
    plt.axvline(f_gauss,color = "k", ls = ":");plt.axvline(i_gauss,color = "k", ls = ":")
    plt.plot(alp.wvf_x,gauss)
    if autozoom == True:
        plt.xlim(i_gauss,f_gauss)
    plt.show()

if check == True:
    plt.plot(alp.wvf_F_x,abs(gauss_f), label = "Gauss filter")
    plt.plot(las.wvf_F_x,abs(las.wvf_F), label = "Laser freq.")
    plt.plot(alp.wvf_F_x,abs(alp.wvf_F), label = "Signal freq.")
    plt.plot(alp.wvf_F_x,abs(pr_signal), label = "Filtered signal freq.")
    plt.xlabel("Frequency in [Hz]");plt.ylabel("Amplitude in ADC counts")
    plt.semilogy();plt.semilogx()
    if autozoom == True:
        plt.ylim(1e-6,np.max(abs(alp.wvf_F))*10)
    plt.legend();plt.show()

# AFTER INVERSE FOURIER TRANSFORM YOU ARE REQUIRED TO SHIFT THE ARRAY TO CENTER THE SIGNAL

amp_test=scipy.fft.irfft(pr_signal)
dec=scipy.fft.irfft(pr_signal/(las.wvf_F/np.max(las.wvf_F)))
if reverse == "True":
    dec = dec[::-1]
dec = np.roll(dec,s_start)
ref = np.argmax(dec)

if check == True:
    conv_int_all,f_conv_all,i_conv_all = signal_int("Deconvolution",dec,timebin,detector,"ALL",out = True)
    # print("\nFull Deconvolution Integral: %.2e"%conv_int)

if detector == "SC":
    for i in range(len(dec[:ref])):
        if dec[ref-i] < 0:
            i_index = i-1
            break
    for j in range(len(alp.wvf[ref:])):
        if alp.wvf[ref+j-1] < 0 < alp.wvf[ref+j]:
            f_index = j
            break
    try:
        conv_int_range,f_conv_range,i_conv_range = signal_int("Deconvolution",dec,timebin,detector,"RANGE",th = thrld,i_range = i_index, f_range = f_index, out= True)
        conv_int_basel,f_conv_basel,i_conv_basel = signal_int("Deconvolution",dec,timebin,detector,"BASEL",th = thrld,i_range = i_index, f_range = f_index, out= True)
    except:
        # conv_int,f_conv,i_conv = signal_int("Convolution",dec,timebin,detector,int_st,th = thrld)
        print("\nERROR: %s could not be properly integrated with type: %s"%("Deconvolution","RANGE"))
else:
    conv_int_basel,f_conv_basel,i_conv_basel = signal_int("Deconvolution",dec,timebin,detector,int_st,th = thrld)
    conv_int_range,f_conv_range,i_conv_range = signal_int("Deconvolution",dec,timebin,detector,"RANGE",th = thrld,i_range = 0,f_range = len(dec)-1)
try:
    pro_alp_int,f_pro_alp,i_pro_alp = signal_int("Filtered Alpha",amp_test,timebin,detector,"BASEL",th =thrld,out = True)
    print("\nCharge increase after filter application:\n %.2f%%\n"%(100*(pro_alp_int/alp_int-1)))
except:
    print("\nERROR: %s could not be properly integrated with type: %s"%("Filtered Alpha","BASEL"))

signal_int("Filtered Alpha",amp_test,timebin,detector,"ALL",out = True)

# SELECT THE INTERVAL LENGHT FOR FIT

if fit == True:
    popt, pcov = curve_fit(func, np.arange(0,ffit)*timebin, dec[ref+ifit:ref+ifit+ffit], bounds =([0,1e-7],[1e3,1e-5]))
    perr = np.sqrt(np.diag(pcov))

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
#________________________PLOT_FINAL_RESULT_____________________________#
########################################################################

print("\nCharge increase after deconvolution:\n %.2f%% with integration type 'RANGE'\n"%(100*(conv_int_range/alp_int-1)))
print("\nCharge increase after deconvolution:\n %.2f%% with integration type 'BASEL'\n"%(100*(conv_int_basel/alp_int-1)))
print("\nCharge increase after deconvolution:\n %.2f%% with integration type 'ALL'\n"%(100*(conv_int_all/alp_int-1)))

f = ROOT.TFile("deconvolution_output/"+filename+".root", "RECREATE")
h = ROOT.TH1D(filename, filename,len(dec),0,len(dec)*4e-9)

for i in range(len(dec)):
    h.SetBinContent(i,dec[i]) # Use the following lines to feed the Fill method of the histograms in order to fill
    h.SetBinError(i,dec[i]*1e-2)

h.Write() # Here write the three histograms on the file and close the file
print(f.GetListOfKeys()[0].GetName())
f.Close()

plt.plot(alp.wvf_x,alp.wvf,label = "Raw Alpha Signal")
# plt.plot(alp.wvf_x,np.roll(alp.wvf,np.argmax(dec)-np.argmax(alp.wvf)),label = "Raw Alpha Signal")
plt.plot(np.arange(len(dec))*timebin,dec,label = "Deconvolution with Gauss filter")
# plt.plot(np.arange(ref+ifit,ref+ifit+ffit)*timebin, func(np.arange(0,ffit)*timebin, *popt), label='Fit: tau_slow = %.2e s'%popt[1])
plt.axvline(f_conv_basel,color = "k", ls = ":");plt.axvline(i_conv_basel,color = "k", ls = ":")
if logy == True:
    plt.grid(); plt.semilogy()
if autozoom == True:
    plt.xlim([0,10e-6]);plt.ylim(1e-3,2*np.max(dec))
plt.xlabel("Time in [s]");plt.ylabel("Amplitude in ADC counts")
plt.legend(); plt.show()