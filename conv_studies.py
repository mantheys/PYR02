import ROOT

import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit
from ir02_lib.deconv_lib import import_scint_prof,my_wvf,conv_guess2, conv_guess3, signal_int, simple_2exp, conv_func2, func2

########################################################################
#______________________IMPORT_SCINTILLATION_FILES______________________#
########################################################################

base_dir = '/pnfs/ciemat.es/data/neutrinos/Super-cells_LAr/Feb22_2/AnalysisROOT/'

path_luz = 'run02_ScintProfFirstSignalBin_SC.root'; label_luz = "SC SPE Signal"
path_alp = 'run55_ScintProfFirstSignalBin_SC.root'; label_alp = "SC MUON Signal"

timebin = 4e-9
norm    = True
detector= "SC"
thrld   = 1e-5 
inv     = False
simple  = False
logy    = False
pro_func= True

alpha = import_scint_prof(base_dir+path_alp,timebin,normalize=norm,trim=False,align=True,start=400,cut_i=0,cut_f=0,invert=inv)
laser = import_scint_prof(base_dir+path_luz,timebin,normalize=norm,trim=False,align=True,start=400,cut_i=0,cut_f=0,invert=inv)

luz_int,f_luz,i_luz = signal_int(label_luz,laser.wvf,timebin,detector,"BASEL",th = thrld,out = True)
alp_int,f_alp,i_alp = signal_int(label_alp,alpha.wvf,timebin,detector,"BASEL",th = thrld,out = True)

if np.argmax(alpha.wvf)-np.argmax(laser.wvf) > 0:
    laser.wvf = np.roll(np.array(laser.wvf),np.argmax(alpha.wvf)-np.argmax(laser[0].wvf))
else:
    alpha.wvf = np.roll(np.array(alpha.wvf),np.argmax(laser.wvf)-np.argmax(alpha.wvf))

########################################################################
#_____________________CONVOLUTION_AND_FIT_PARAMETERS___________________#
########################################################################

t_fast = 6e-9; t_int = 1e-7; t_slow = 1e-6
amp_fast = 1e-8; amp_int = 0.01; amp_slow = 1e-8
sigma = 2e-8

# if simple == True:
#     fit_initials = (t_fast,t_slow,amp_fast,amp_slow)
#     fit_limits = ([1e-9,1e-7,1e-9,1e-9],[9e-9,1e-5,1e-7,1e-7])
#     popt, pcov = curve_fit(conv_guess2,laser,alpha.wvf, p0 = fit_initials, bounds = fit_limits)
#     conv = conv_guess2(laser,*popt)
#     print("\nTau_fast = %.2e [s] - Amp = %.2e"%(popt[0],popt[2]))
#     print("Tau_slow = %.2e [s]- Amp = %.2e"%(popt[1],popt[3]))
#     plt.plot(simple_2exp(alpha.wvf_x,*popt),label="Convolution Func.")
#     plt.legend()
#     plt.semilogy()
#     plt.show()

# if int_exp == True:
#     fit_initials = (t_fast,t_int,t_slow,amp_fast,amp_int,amp_slow)
#     fit_limits = ([1e-9,1e-8,1e-7,1e-2,1e-2,1e-2],[1e-7,1e-6,1e-5,1,1,1])
#     popt, pcov = curve_fit(conv_guess3,laser,alpha.wvf, p0 = fit_initials, bounds = fit_limits)
#     conv = conv_guess3(laser,*popt)
#     print("Tau_fast = %.2e [s] - Amp = %.2e"%(popt[0],popt[3]))
#     print("Tau_int = %.2e [s] - Amp = %.2e"%(popt[1],popt[4]))
#     print("Tau_slow = %.2e [s]- Amp = %.2e"%(popt[2],popt[5]))

if pro_func == True:
    fit_initials = (t_fast,t_slow,amp_fast,amp_slow,sigma)
    fit_limits = ([1e-9,1e-7,1e-9,1e-9,1e-9],[9e-8,1e-5,1e-7,1e-7,1e-7])
    popt, pcov = curve_fit(conv_func2,laser,alpha.wvf, p0 = fit_initials, bounds = fit_limits,method="trf")
    conv = conv_func2(laser,*popt)
    conv_int,f_conv,i_conv = signal_int("CONV FUNC",func2(alpha.wvf_x,*popt),timebin,"SiPM","MIXED",th = thrld,out = True)
    print(100*(conv_int/alp_int-1))
    print("\nTau_fast = %.2e [s] - Amp = %.2e"%(popt[0],popt[2]))
    print("Tau_slow = %.2e [s]- Amp = %.2e"%(popt[1],popt[3]))
    print("Sigma = %.2e [s]"%(popt[4]))
    print("Slow Percent. = %.2f"%(100*popt[3]/(popt[3]+popt[2])))
    
    plt.plot(func2(alpha.wvf_x,*popt),label="Convolution Func.")
    plt.legend()
    plt.semilogy()
    plt.show()

########################################################################
#________________________PLOT_FINAL_RESULT_____________________________#
########################################################################

if logy == True:
    plt.semilogy(); plt.xlim(0,20e-6); plt.ylim(1e-6, 1.2)

plt.plot(laser.wvf_x,laser.wvf,label = label_luz)
plt.plot(alpha.wvf_x,alpha.wvf,label = label_alp)
plt.plot(laser.wvf_x,conv,label = "Fitted Convolution")
plt.xlabel("Time in [s]"); plt.ylabel("Normalized Amplitude")  
plt.axhline(0,color = "grey", ls = ":")
plt.legend()
plt.show()