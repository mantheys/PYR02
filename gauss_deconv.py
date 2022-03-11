import ROOT
import sys
import scipy
import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit
from ir02_lib.deconv_lib import import_scint_prof, pdf,func,my_wvf,signal_int

########################################################################
#______________________IMPORT_SCINTILLATION_FILES______________________#
########################################################################

print("\n### WELCOME TO THE DECONVOLUTION STUDIES ###")

base_dir   = '/pnfs/ciemat.es/data/neutrinos/Super-cells_LAr/Jan22/AnalysisROOT/' 

# detector = "SiPM"; timebin = 4e-9; int_st = "BASEL";
# filter_strenght = 400
# path_alp = 'run13_ScintProfFirstSignalBin_SiPM.root'
# path_las = 'run56_ScintProfFirstSignalBin_SiPM.root'
# path_spe = 'run2_ScintProfFirstSignalBin_SiPM.root'
# filename = "JAN_SiPM_DAY1"

# path_alp = 'run36_ScintProfFirstSignalBin_SiPM.root'
# path_las = 'run56_ScintProfFirstSignalBin_SiPM.root'
# path_spe = 'run2_ScintProfFirstSignalBin_SiPM.root'
# filename = "JAN_SiPM_DAY2"

# path_alp = 'run79_ScintProfFirstSignalBin_SiPM.root'
# path_las = 'run56_ScintProfFirstSignalBin_SiPM.root'
# path_spe = 'run2_ScintProfFirstSignalBin_SiPM.root'
# filename = "JAN_SiPM_DAY3"

detector = "SC"; timebin = 4e-9; int_st = "BASEL";
filter_strenght = 50
# path_alp = 'run13_ScintProfFirstSignalBin_SC.root'
# path_las = 'Jan_ScintProf_run55_ch6.root'
# path_spe = 'run2_ScintProfFirstSignalBin_SC.root'
# filename = "JAN_SC_DAY1"

path_alp = 'Jan_ScintProf_run36_ch6.root'
path_las = 'Jan_ScintProf_run56_ch6.root'
path_spe = 'run63_ScintProfFirstSignalBin_SC.root'
filename = "JAN_SC_DAY2"

# path_alp = 'run79_ScintProfFirstSignalBin_SC.root'
# path_las = 'Jan_ScintProf_run55_ch6.root'
# path_spe = 'run68_ScintProfFirstSignalBin_SC.root'
# filename = "JAN_SC_DAY3"

# SELECT THE RIGHT TIME BIN FOR ACCURATE PLOT REPRESENATION
check = True; norm = False; shift = True; th = 1e-5; 

alp = import_scint_prof(base_dir+path_alp,timebin,normalize=norm,align=shift)
las = import_scint_prof(base_dir+path_las,timebin,normalize=norm,align=shift)
spe = import_scint_prof(base_dir+path_spe,timebin,normalize=norm,align=shift)

roll = 400; ifit = 100; ffit = 400 #(fit starts at bin rollber and ends at bin rollber+ffit)

# CALCULATE INTEGRAL
alp_int,f_alp,i_alp = signal_int("Alpha",alp.wvf,timebin,detector,int_st,th=np.max(alp.wvf)*th,out = True)
las_int,f_las,i_las = signal_int("Laser",las.wvf,timebin,detector,int_st,th=np.max(las.wvf)*th)
spe_int,f_spe,i_spe = signal_int("SPE"  ,spe.wvf,timebin,detector,int_st,th=np.max(spe.wvf)*th)
# print("\nNumber of SPE in Laser pulse:\n %.2f SPE"%(las_int/spe_int))

# FIRST PLOT AND INTEGRAL CALCULATION
if check == True:
    plt.plot(alp.wvf_x,alp.wvf,label="Raw Alpha signal")
    plt.plot(las.wvf_x,las.wvf,label="Raw Laser signal")
    plt.plot(spe.wvf_x,spe.wvf,label="Raw SPE signal")

    plt.xlim(0,5e-6);plt.ylim()
    plt.xlabel("Time in [s]");plt.ylabel("Amplitude in ADC counts")
    plt.axvline(f_alp*4e-9,color = "k", ls = ":");plt.axvline(i_alp*4e-9,color = "k", ls = ":")
    plt.legend();plt.show()

########################################################################
#____________________DECONVOLUTION_AND_FIT_PARAMETERS__________________#
########################################################################

# SELECT A FILTER STRENGHT THAT ALLOWS YOU TO SEE THE SLOPE OF THE SIGNAL
# THE CORRECPONDING NORM GAUSS FILTER IS DEFINED AND MULTIPLIED IN FOURIER SPACE

gauss_f = pdf(np.arange(len(alp.wvf_F_x)), m = 0, sd = filter_strenght, norm = 1, n = 2.)
pr_signal = alp.wvf_F*(gauss_f)

if check == True:
    plt.plot(alp.wvf_F_x,abs(gauss_f), label = "Gauss filter")
    plt.plot(las.wvf_F_x,abs(las.wvf_F), label = "Laser freq.")
    plt.plot(alp.wvf_F_x,abs(alp.wvf_F), label = "Signal freq.")
    plt.plot(alp.wvf_F_x,abs(pr_signal), label = "Filtered signal freq.")
    plt.xlabel("Frequency in [Hz]]");plt.ylabel("Amplitude in ADC counts")
    plt.semilogy();plt.semilogx()
    plt.legend();plt.show()

# AFTER INVERSE FOURIER TRANSFORM YOU ARE REQUIRED TO SHIFT THE ARRAY TO CENTER THE SIGNAL

amp_test=scipy.fft.irfft(pr_signal)
gauss = np.roll(scipy.fft.irfft(gauss_f),int(len(gauss_f)/2))
dec=scipy.fft.irfft(pr_signal/(las.wvf_F/np.max(las.wvf_F)))
dec = np.roll(dec,400-np.argmax(dec))
ref = np.argmax(dec)

conv_int,f_conv,i_conv = signal_int("Convolution",dec,timebin,detector,int_st,th = np.max(dec)*1e-4)
pro_alp_int,f_pro_alp,i_pro_alp = signal_int("Filtered Alpha",amp_test,timebin,detector,int_st,th = np.max(amp_test)*th,out = True)

# SELECT THE INTERVAL LENGHT FOR FIT

popt, pcov = curve_fit(func, np.arange(0,ffit)*timebin, dec[ref+ifit:ref+ifit+ffit], bounds =([0,1e-7],[1e3,1e-5]))
perr = np.sqrt(np.diag(pcov))

if check == True:
    plt.plot(alp.wvf_x,alp.wvf ,label="Raw Alpha signal")
    plt.plot(alp.wvf_x,amp_test,label="Filtered Alpha Signal")
    plt.xlabel("Time in [s]");plt.ylabel("Amplitude in ADC counts")
    plt.semilogy();plt.xlim(alp.wvf_x[np.argmax(alp.wvf)-200],alp.wvf_x[np.argmax(alp.wvf)+200]);plt.ylim(1e-3,2*np.max(alp.wvf))
    plt.legend();plt.show()

########################################################################
#________________________PLOT_FINAL_RESULT_____________________________#
########################################################################

# plt.plot(alp.wvf_x,np.roll(alp.wvf,np.argmax(dec)-np.argmax(alp.wvf)),label = "Raw Alpha Signal")
plt.plot(alp.wvf_x,np.roll(alp.wvf,np.argmax(dec)-np.argmax(alp.wvf)),label = "Raw Alpha Signal")
plt.plot(np.arange(len(dec))*timebin,dec,label = "Deconvolution with Gauss filter")
plt.plot(np.arange(ref+ifit,ref+ifit+ffit)*timebin, func(np.arange(0,ffit)*timebin, *popt), label='Fit: tau_slow = %.2e s'%popt[1])

plt.axvline(f_conv*4e-9,color = "k", ls = ":");plt.axvline(i_conv*4e-9,color = "k", ls = ":")
plt.grid(); plt.semilogy();plt.xlim([0,10e-6]);plt.ylim(1e-3,2*np.max(dec))
plt.xlabel("Time in [s]");plt.ylabel("Amplitude in ADC counts")
plt.legend(); plt.show()

print("\nCharge increase after deconvolution:\n %.2f%%\n"%(100*(conv_int/alp_int-1)))

f = ROOT.TFile(filename+".root", "RECREATE")
h = ROOT.TH1D(filename, filename,len(dec),0,len(dec)*4e-9)

for i in range(len(dec)):
    h.SetBinContent(i,dec[i]) # Use the following lines to feed the Fill method of the histograms in order to fill
    h.SetBinError(i,dec[i]*1e-2)

h.Write() # Here write the three histograms on the file and close the file
print(f.GetListOfKeys()[0].GetName())
f.Close()