import ROOT

import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit
from ir02_lib.deconv_lib import import_scint_prof,my_wvf,conv_guess2, conv_guess3

########################################################################
#______________________IMPORT_SCINTILLATION_FILES______________________#
########################################################################

base_dir = '/pnfs/ciemat.es/data/neutrinos/Super-cells_LAr/Feb22/AnalysisROOT/'

path_luz = 'run46_ScintProfFirstSignalBin_SC.root'; label_luz = "SC SPE Signal"
path_alp = 'run34_ScintProfFirstSignalBin_ch6.root'; label_alp = "SC Alpha Signal"

timebin = 4e-9
alpha = import_scint_prof(base_dir+path_alp,timebin,normalize=True,trim=False,align=True,start=400,cut_i=0,cut_f=0)
laser = import_scint_prof(base_dir+path_luz,timebin,normalize=True,trim=False,align=True,start=400,cut_i=0,cut_f=0)

if np.argmax(alpha.wvf)-np.argmax(laser.wvf) > 0:
    laser.wvf = np.roll(np.array(laser.wvf),np.argmax(alpha.wvf)-np.argmax(laser[0].wvf))
else:
    alpha.wvf = np.roll(np.array(alpha.wvf),np.argmax(laser.wvf)-np.argmax(alpha.wvf))

plt.plot(laser.wvf_x,laser.wvf,label = label_luz)
plt.plot(alpha.wvf_x,alpha.wvf,label = label_alp)

plt.xlabel("Time in [s]"); plt.ylabel("Normalized Amplitude")

########################################################################
#_____________________CONVOLUTION_AND_FIT_PARAMETERS___________________#
########################################################################

int_exp = True
t_fast = 1e-8; t_int = 1e-7; t_slow = 1e-6
amp_fast = 0.9; amp_int = 0.01; amp_slow = 0.1

if int_exp == False:
    fit_initials = (t_fast,t_slow,amp_fast,amp_slow)
    fit_limits = ([1e-9,1e-7,1e-2,1e-2],[1e-7,1e-5,1,1])
    popt, pcov = curve_fit(conv_guess2,laser,alpha.wvf, p0 = fit_initials, bounds = fit_limits)
    conv = conv_guess2(laser,*popt)
    print("Tau_fast = %.2e [s] - Amp = %.2e"%(popt[0],popt[2]))
    print("Tau_slow = %.2e [s]- Amp = %.2e"%(popt[1],popt[3]))

if int_exp == True:
    fit_initials = (t_fast,t_int,t_slow,amp_fast,amp_int,amp_slow)
    fit_limits = ([1e-9,1e-8,1e-7,1e-2,1e-2,1e-2],[1e-7,1e-6,1e-5,1,1,1])
    popt, pcov = curve_fit(conv_guess3,laser,alpha.wvf, p0 = fit_initials, bounds = fit_limits)
    conv = conv_guess3(laser,*popt)
    print("Tau_fast = %.2e [s] - Amp = %.2e"%(popt[0],popt[3]))
    print("Tau_int = %.2e [s] - Amp = %.2e"%(popt[1],popt[4]))
    print("Tau_slow = %.2e [s]- Amp = %.2e"%(popt[2],popt[5]))

########################################################################
#________________________PLOT_FINAL_RESULT_____________________________#
########################################################################

logy = False
if logy == True:
    plt.semilogy(); plt.xlim(0,20e-6); plt.ylim(1e-6, 1.2)
    
plt.plot(laser.wvf_x,conv,label = "Fitted Convolution")
plt.axhline(0,color = "grey", ls = ":")
plt.legend()
plt.show()