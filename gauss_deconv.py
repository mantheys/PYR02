import scipy
import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit
from deconv_lib import import_scint_prof, pdf,func,my_wvf

########################################################################
#______________________IMPORT_SCINTILLATION_FILES______________________#
########################################################################

base_dir = 'ref_deconv/' 
path_luz='ScintProf_run46_ch0.root'
path_alp='ScintProf_run34_ch0.root'

# SELECT THE RIGHT TIME BIN FOR ACCURATE PLOT REPRESENATION
timebin = 4e-9
alpha = import_scint_prof(base_dir+path_alp,timebin,normalize=True)
laser = import_scint_prof(base_dir+path_luz,timebin,normalize=True)

# FORMAT DATA IN COMPACT LIST
x = [alpha.wvf_F_x,alpha.wvf_F,laser.wvf_F]

########################################################################
#____________________DECONVOLUTION_AND_FIT_PARAMETERS__________________#
########################################################################

# SELECT A FILTER STRENGHT THAT ALLOWS YOU TO SEE THE SLOPE OF THE SIGNAL
# THE CORRECPONDING NORM GAUSS FILTER IS DEFINED AND MULTIPLIED IN FOURIER SPACE
filter_strenght = 300
gauss_f = pdf(np.arange(len(x[0])), m=0, sd=filter_strenght, norm=1, n=2.)
processed_signal=x[1]*(gauss_f)

# AFTER INVERSE FOURIER TRANSFORM YOU ARE REQUIRED TO SHIFT THE ARRAY TO CENTER THE SIGNAL
roll_num = 200
dec2=scipy.fft.irfft(processed_signal/x[2])
dec2=np.roll(dec2,roll_num)

# SELECT THE INTERVAL LENGHT FOR FIT
fit_space = 400 #(fit starts at bin roll_number and ends at bin roll_number+fit_space)
popt, pcov = curve_fit(func, np.arange(0,fit_space)*timebin, dec2[roll_num:roll_num+fit_space],p0=(1e-4,1e-6))
print(popt)

########################################################################
#________________________PLOT_FINAL_RESULT_____________________________#
########################################################################

plt.plot(np.arange(len(dec2))*timebin,dec2,label = "Deconvolution with Gauss filter")
plt.plot(np.arange(roll_num,roll_num+fit_space)*timebin, func(np.arange(0,fit_space)*timebin, *popt), label='fit: tau_slow =%.2e s'%popt[1])
plt.grid(); plt.semilogy()
plt.xlabel("Time in [s]");plt.ylabel("Amplitude in a.u.")
plt.xlim([0,7.5e-6])
plt.legend(); plt.show()