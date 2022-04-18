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
detector,timebin,int_st,paths,filename,shift,filter_strenght,reverse,fnal,s_start = import_deconv_runs("deconvolution_input/FEB_SC_DAY2_OV1.txt",debug = True)
check = False; autozoom = False; logy = True; norm = True; fit = False; thrld = 1e-3; pro_abs = False; pro_rodrigo = False
term_output = check
# SELECT THE RIGHT TIME BIN FOR ACCURATE PLOT REPRESENATION
init=0; trm = False; #fnal=500; s_start = 600;
alp = import_scint_prof(paths[0]+paths[1],timebin,normalize=norm,trim=trm,align=shift,start=s_start,cut_i=init,cut_f=fnal)
las = import_scint_prof(paths[0]+paths[2],timebin,normalize=norm,trim=trm,align=shift,start=s_start,cut_i=init,cut_f=fnal)
spe = import_scint_prof(paths[0]+paths[3],timebin,normalize=norm,trim=trm,align=shift,start=s_start,cut_i=0,cut_f=0)

# CALCULATE INTEGRAL OF RAW DATA
labels = ["Alpha","Laser","SPE"]

# CHECK INTEGRATION METHODS FOR RAW ALPHA SIGNALq
alp_int_all,f_alp_all,i_alp_all       = signal_int(labels[0],alp.wvf,timebin,detector,"ALL"  ,th = thrld,out = term_output)
alp_int_thrld,f_alp_thrld,i_alp_thrld = signal_int(labels[0],alp.wvf,timebin,detector,"THRLD",th = thrld,out = term_output)

if detector == "SiPM":
    alp_int,f_alp,i_alp = signal_int(labels[0],alp.wvf,timebin,detector,"MIXED"  ,th = thrld,out = True)
    signal_int(labels[1],las.wvf,timebin,detector,"MIXED",th = thrld,out=True)
    signal_int(labels[2],spe.wvf,timebin,detector,"MIXED",th = thrld,out=term_output)
    plt.axvline(f_alp,color = "b",ls = ":");plt.axvline(i_alp,color = "b",ls = ":",label = "Alpha Integration (MIXED) Limits")

if detector == "SC":
    alp_int,f_alp,i_alp = signal_int(labels[0],alp.wvf,timebin,detector,"BASEL",th = thrld,out = True)
    alp_int_under,f_alp_under,i_alp_under = signal_int(labels[0],alp.wvf,timebin,detector,"UNDER",th = thrld,out = True)
    signal_int(labels[1],las.wvf,timebin,detector,"BASEL",th = thrld,out=True)
    signal_int(labels[2],spe.wvf,timebin,detector,"BASEL",th = thrld,out=term_output)
    print("\nCharge content of undershoot:\n %.2f%%\n"%(100*(-alp_int_under/alp_int)))
    plt.axvline(f_alp,color = "g",ls = ":");plt.axvline(i_alp,color = "g",ls = ":",label = "Alpha Integration (BASEL) Limits")
    plt.axvline(f_alp_under,color = "k",ls = ":");plt.axvline(i_alp_under,color = "k",ls = ":",label = "Alpha Integration (UNDER) Limits")

# FIRST PLOT TO VISUALIZE DATA (VERTICAL LINES SHOW INTEGRATION LIMITS OF ALPHA SIGNAL)
plt.plot(spe.wvf_x,spe.wvf,c="g",label="Raw SPE signal")
plt.plot(alp.wvf_x,alp.wvf,c="b",label="Raw Alpha signal")
plt.plot(las.wvf_x,las.wvf,c="r",label="Raw Laser signal")

plt.xlabel("Time in [s]");plt.ylabel("Amplitude in ADC counts")
plt.axhline(0, ls="--", color = "k",alpha=0.25)

if autozoom == True:
    plt.xlim(timebin*(np.argmax(alp.wvf)-100),timebin*(np.argmax(alp.wvf)+1000));plt.ylim()
if logy == True:
        plt.semilogy()
plt.legend();plt.show()

print("\n-------------------------------------------------------------")

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
    gauss_int,f_gauss,i_gauss = signal_int("Gauss",gauss,timebin,detector,"BASEL",th = np.max(gauss)*1e-3,out = term_output)
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
# signal_int("Filtered Alpha",amp_test,timebin,detector,"ALL",out = term_output)

if check == True:
    plt.plot(alp.wvf_x,alp.wvf ,label="Raw Alpha signal")
    plt.plot(alp.wvf_x,amp_test,label="Filtered Alpha Signal",c="b")
    plt.xlabel("Time in [s]");plt.ylabel("Amplitude in ADC counts")
    
    if logy == True:
        plt.semilogy()
    
    if autozoom == True:
        plt.xlim(alp.wvf_x[np.argmax(alp.wvf)-200],alp.wvf_x[np.argmax(alp.wvf)+200]);plt.ylim(1e-3,2*np.max(alp.wvf))
    
    plt.legend();plt.show()
    print("\n-------------------------------------------------------------")

########################################################################
#_________EVALUATE DECONVOLUTION AND FILTER IN FOURIER SPACE___________#
########################################################################

if check == True:
    plt.plot(alp.wvf_F_x,abs(gauss_f), label = "Gauss filter")
    # plt.plot(las.wvf_F_x,abs(las.wvf_F), label = "Laser freq.")
    plt.plot(las.wvf_F_x,abs(las.wvf_F/np.max(las.wvf_F)),c="r",label = "Norm Laser freq.")
    plt.plot(alp.wvf_F_x,abs(alp.wvf_F),c="b",label = "Signal freq.")
    # plt.plot(alp.wvf_F_x,abs(pr_signal), label = "Filtered signal freq.")
    plt.plot(alp.wvf_F_x,abs(pr_signal/(las.wvf_F/np.max(las.wvf_F))), label = "Deconvolved signal (norm laser)")
    # plt.plot(alp.wvf_F_x,abs(pr_signal/las.wvf_F), label = "Deconvolved signal")
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

# if check == True:
conv_int_all,f_conv_all,i_conv_all = signal_int("Deconvolution",dec,timebin,detector,"ALL",out = term_output)
    # print("\nFull Deconvolution Integral: %.2e"%conv_int)
   
conv_int_range,f_conv_range,i_conv_range = signal_int("Deconvolution",dec,timebin,detector,"RANGE",th = thrld,i_range = np.argmax(alp.wvf)-int(i_alp/timebin),f_range = int(f_alp/timebin)-np.argmax(alp.wvf))
conv_int_thrld,f_conv_thrld,i_conv_thrld = signal_int("Deconvolution",dec,timebin,detector,"THRLD",th = thrld,out= term_output)
conv_int_basel,f_conv_basel,i_conv_basel = signal_int("Deconvolution",dec,timebin,detector,"BASEL",th = thrld,out= term_output)
conv_int_fixed,f_conv_fixed,i_conv_fixed = signal_int("Deconvolution",dec,timebin,detector,"FIXED",th = thrld,out= term_output)
conv_int_mixed,f_conv_mixed,i_conv_mixed = signal_int("Deconvolution",dec,timebin,detector,"MIXED",th = thrld,out= term_output)

if detector == "SC":
    try:
        pro_alp_int,f_pro_alp,i_pro_alp = signal_int("Filtered Alpha",amp_test,timebin,detector,"BASEL",th =thrld,out = term_output)
        conv_int_under,f_conv_under,i_conv_under = signal_int("Deconvolution",dec,timebin,detector,"UNDER",th = 0, out= term_output)
        print("\nCharge increase after filter application:\n %.4f%% with integration type 'BASEL'\n"%(100*(pro_alp_int/alp_int-1)))
    except:
        print("\nERROR: %s could not be properly integrated with type: %s"%("Filtered Alpha","BASEL"))

if detector == "SiPM":
    try:
        pro_alp_int,f_pro_alp,i_pro_alp = signal_int("Filtered Alpha",amp_test,timebin,detector,"MIXED",th =thrld,out = term_output)
        print("\nCharge increase after filter application:\n %.4f%% with integration type 'MIXED'\n"%(100*(pro_alp_int/alp_int-1)))
    except:
        print("\nERROR: %s could not be properly integrated with type: %s"%("Filtered Alpha","MIXED"))

print("\n-------------------------------------------------------------")

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
if check == True:
    print("\nCharge increase after deconvolution:\n %.2f%% with integration type 'ALL'  \n"%(100*(conv_int_all/alp_int-1)))
    print("\nCharge increase after deconvolution:\n %.2f%% with integration type 'FIXED'\n"%(100*(conv_int_fixed/alp_int-1)))
    print("\nCharge increase after deconvolution:\n %.2f%% with integration type 'BASEL'\n"%(100*(conv_int_basel/alp_int-1)))
    print("\nCharge increase after deconvolution:\n %.2f%% with integration type 'THRLD'\n"%(100*(conv_int_thrld/alp_int-1)))
    print("\nCharge increase after deconvolution:\n %.2f%% with integration type 'BASEL' substracting 'UNDER'\n"%(100*((conv_int_basel+conv_int_under)/alp_int-1)))
    # plt.axvline(f_conv_basel,color = "r", ls = ":");plt.axvline(i_conv_basel,color = "r", ls = ":",label = "(BASEL) Integration Limits")
    # plt.axvline(f_conv_fixed,color = "k", ls = ":");plt.axvline(i_conv_fixed,color = "k", ls = ":",label = "(FIXED) Integration Limits")

print("\nCharge increase after deconvolution:\n %.2f%% with integration type 'MIXED'\n"%(100*(conv_int_mixed/alp_int-1)))
# print("\nCharge increase after deconvolution:\n %.2f%% with integration type 'MIXED'\n"%(100*(1-alp_int/conv_int_mixed)))
print("\nCharge increase after deconvolution:\n %.2f%% with integration type 'RANGE'\n"%(100*(conv_int_range/alp_int-1)))
# print("\nCharge increase after deconvolution:\n %.2f%% with integration type 'RANGE'\n"%(100*(1-alp_int/conv_int_range)))

# PLOT DECONVOLUTED WAVEFORM AND RAW DATA
plt.plot(alp.wvf_x,alp.wvf,label = "Raw Alpha Signal")
plt.plot(np.arange(len(dec))*timebin,dec,label = "Deconvolution with Gauss filter")
plt.axhline(0, ls="--", color = "k",alpha=0.25)
plt.axvline(f_conv_mixed,color = "b", ls = ":");plt.axvline(i_conv_mixed,color = "b", ls = ":",label = "(MIXED) Integration Limits")
plt.axvline(f_conv_range,color = "g", ls = ":");plt.axvline(i_conv_range,color = "g", ls = ":",label = "(RANGE) Integration Limits")
# plt.axvline(f_conv_thrld,color = "b", ls = ":");plt.axvline(i_conv_thrld,color = "b", ls = ":",label = "(THRLD) Integration Limits")
plt.xlabel("Time in [s]");plt.ylabel("Amplitude in ADC counts")

if logy == True:
    plt.grid(); plt.semilogy()

if autozoom == True:
    plt.xlim([0,10e-6]);plt.ylim(1e-3,2*np.max(dec))

plt.legend(); plt.show()

# EXPORT DECONVOLUTED WAVEFORM TO ROOT FORMAT
f = ROOT.TFile("deconvolution_output/"+filename+"_f%i.root"%filter_strenght, "RECREATE")
h = ROOT.TH1D(filename+"_f%i"%filter_strenght, filename+"_f%i"%filter_strenght,len(dec),0,len(dec)*4e-9)
for i in range(len(dec)):
    h.SetBinContent(i,dec[i]) # Use the following lines to feed the Fill method of the histograms in order to fill
    h.SetBinError(i,dec[i]*1e-2)
h.Write() # Here write the three histograms on the file and close the file
print("Generated output file: %s.root"%f.GetListOfKeys()[0].GetName())
f.Close()