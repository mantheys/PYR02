from cProfile import label
import scipy
import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit
from ir02_lib.deconv_lib import import_scint_prof, pdf,func,my_wvf,signal_int

########################################################################
#______________________IMPORT_SCINTILLATION_FILES______________________#
########################################################################

base_dir = '/pnfs/ciemat.es/data/neutrinos/Super-cells_LAr/Feb22_2/AnalysisROOT/' 
path_laser = 'ScintProfFirstSignalBin_run03_ch0.root'
path_alp   = 'ScintProf_run27_ch0.root'
path_SPE   = 'ScintProfFirstSignalBin_run03_ch0.root'

# SELECT THE RIGHT TIME BIN FOR ACCURATE PLOT REPRESENATION
timebin = 4e-9; norm = False
alpha = import_scint_prof(base_dir+path_alp,timebin,normalize=norm)
laser = import_scint_prof(base_dir+path_laser,timebin,normalize=norm)
SPE = import_scint_prof(base_dir+path_SPE,timebin,normalize=norm)

# CALCULATE INTEGRAL 
x = [alpha.wvf_F_x,alpha.wvf_F,laser.wvf_F]
alpha_int,end_alpha,start_alpha = signal_int("Alpha",alpha.wvf,timebin,"SiPM",th = 1)
laser_int,end_laser,start_laser = signal_int("Laser",laser.wvf,timebin,"SiPM",th = 1)
SPE_int,end_SPE,start_SPE = signal_int("SPE",SPE.wvf,timebin,"SiPM",th = 1)

# FIRST PLOT AND INTEGRAL CALCULATION
plt.plot(alpha.wvf_x,alpha.wvf,label = "Raw Alpha signal")
plt.plot(laser.wvf_x,laser.wvf, label = "Raw Light signal")
plt.axvline(end_alpha*4e-9,color = "k", ls = ":");plt.axvline(start_alpha*4e-9,color = "k", ls = ":")
plt.xlabel("Time in [s]");plt.ylabel("Amplitude in ADC counts")
plt.legend();plt.show()

########################################################################
#____________________DECONVOLUTION_AND_FIT_PARAMETERS__________________#
########################################################################

# SELECT A FILTER STRENGHT THAT ALLOWS YOU TO SEE THE SLOPE OF THE SIGNAL
# THE CORRECPONDING NORM GAUSS FILTER IS DEFINED AND MULTIPLIED IN FOURIER SPACE
filter_strenght = 500
gauss_f = pdf(np.arange(len(alpha.wvf_F_x)), m = 0, sd = filter_strenght, norm = 1, n = 2.)
processed_signal = alpha.wvf_F*(gauss_f)
plt.plot(laser.wvf_F_x,abs(laser.wvf_F), label = "Laser freq.")
plt.plot(alpha.wvf_F_x,abs(alpha.wvf_F), label = "Signal freq.")
plt.plot(alpha.wvf_F_x,abs(processed_signal), label = "Filtered signal freq.")
plt.semilogy();plt.semilogx()
plt.legend();plt.show()

# AFTER INVERSE FOURIER TRANSFORM YOU ARE REQUIRED TO SHIFT THE ARRAY TO CENTER THE SIGNAL
roll_num = 200
dec2=scipy.fft.irfft(processed_signal/laser.wvf_F)
dec2=np.roll(dec2,roll_num)
conv_int,end_conv,start_conv = signal_int("Convolution",dec2,timebin,"SiPM",th = 1e-3)

# SELECT THE INTERVAL LENGHT FOR FIT
fit_space = 400 #(fit starts at bin roll_number and ends at bin roll_number+fit_space)
popt, pcov = curve_fit(func, np.arange(0,fit_space)*timebin, dec2[roll_num:roll_num+fit_space], bounds =([0,1e-7],[1,1e-5]))
print(popt)

########################################################################
#________________________PLOT_FINAL_RESULT_____________________________#
########################################################################

plt.plot(np.arange(len(dec2))*timebin,dec2,label = "Deconvolution with Gauss filter")
plt.plot(np.arange(roll_num,roll_num+fit_space)*timebin, func(np.arange(0,fit_space)*timebin, *popt), label='fit: tau_slow = %.2e s'%popt[1])
plt.axvline(end_conv*4e-9,color = "k", ls = ":");plt.axvline(start_conv*4e-9,color = "k", ls = ":")
plt.grid(); plt.semilogy()
plt.xlabel("Time in [s]");plt.ylabel("Amplitude in a.u.")
plt.xlim([0,7.5e-6])
plt.legend(); plt.show()