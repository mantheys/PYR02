import ROOT as root

import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit
from ir02_lib.deconv_lib import import_scint_prof,my_wvf,signal_int,conv_func2,func2,import_deconv_runs,logconv_func2

########################################################################
#______________________IMPORT_SCINTILLATION_FILES______________________#
########################################################################

print("\n### WELCOME TO THE CONVOLUTION STUDIES ###\n")

path = "deconvolution_input/GAr/"
decon_runs = "GAr_SiPM_DAY21_OV2_T-180"
output_file = open("scint_fit/"+decon_runs+"_CONV.txt","w")
detector,timebin,int_st,paths,filename,shift,filter_strenght,reverse,fnal,s_start,particle_label = import_deconv_runs(path+decon_runs+".txt",debug = True)

limit       = 1
check       = False
autozoom    = True
logy        = True
norm        = True
fit         = False
thrld       = 1e-4
pro_abs     = False
pro_rodrigo = False
inv         = False
term_output = check

# SELECT THE RIGHT TIME BIN FOR ACCURATE PLOT REPRESENATION
init=0; trm = False; #fnal=500; s_start = 600;
alpha = import_scint_prof(paths[0]+paths[1],timebin,normalize=False,trim=trm,align=shift,start=s_start,cut_i=init,cut_f=fnal,invert=inv)
laser = import_scint_prof(paths[0]+paths[2],timebin,normalize=False,trim=trm,align=shift,start=s_start,cut_i=init,cut_f=fnal,invert=False)
# laser = import_scint_prof(paths[0]+paths[3],timebin,normalize=norm,trim=trm,align=shift,start=s_start,cut_i=init,cut_f=fnal,invert=False)

# alpha = import_scint_prof(base_dir+path_alp,timebin,normalize=False,trim=False,align=True,start=400,cut_i=0,cut_f=0,invert=inv)
# laser = import_scint_prof(base_dir+path_luz,timebin,normalize=False,trim=False,align=True,start=400,cut_i=0,cut_f=0,invert=inv)

label_luz = "LASER"
label_alp = "ALPHA"

if norm == True:
    max_alpha = np.max(alpha.wvf)
    max_laser = np.max(laser.wvf)

    alpha.wvf = alpha.wvf/max_alpha
    laser.wvf = laser.wvf/max_laser

luz_int_all,f_luz,i_luz = signal_int(label_luz,laser.wvf,timebin,detector,"ALL",th = thrld,out = True)
alp_int_all,f_alp,i_alp = signal_int(label_alp,alpha.wvf,timebin,detector,"ALL",th = thrld,out = True)

luz_int,f_luz,i_luz = signal_int(label_luz,laser.wvf,timebin,detector,"BASEL",th = thrld,out = True)
alp_int,f_alp,i_alp = signal_int(label_alp,alpha.wvf,timebin,detector,"BASEL",th = thrld,out = True)

if np.argmax(alpha.wvf)-np.argmax(laser.wvf) > 0:
    laser.wvf = np.roll(np.array(laser.wvf),np.argmax(alpha.wvf)-np.argmax(laser[0].wvf))
else:
    alpha.wvf = np.roll(np.array(alpha.wvf),np.argmax(laser.wvf)-np.argmax(alpha.wvf))

# for i in range(len(alpha.wvf)):
#     if alpha.wvf[i] <= 0:
#         alpha.wvf[i] = 1e-8
#     if np.isnan(alpha.wvf[i]):
#         alpha.wvf[i] = 1e-8
#     if laser.wvf[i] <= 0:
#         laser.wvf[i] = 1e-8
#     if np.isnan(alpha.wvf[i]):
#         laser.wvf[i] = 1e-8

# alpha.wvf = tfile_vect2array("ITALIAN_WVFs.root","hmu")
# laser.wvf = tfile_vect2array("ITALIAN_WVFs.root","hspe1")

# alpha.wvf_x = np.arange(0, len(alpha.wvf)*4e-9, 4e-9)
# laser.wvf_x = alpha.wvf_x

########################################################################
#_____________________CONVOLUTION_AND_FIT_PARAMETERS___________________#
########################################################################

# MUON SC CONFIG
t_fast   = 2e-8; t_fast_low   = 1e-8; t_fast_high   = 4e-8
t_slow   = 1e-6; t_slow_low   = 6e-7; t_slow_high   = 5e-6
amp_fast = 2e-8; amp_fast_low = 1e-8; amp_fast_high = 3e-8
amp_slow = 5e-8; amp_slow_low = 1e-8; amp_slow_high = 9e-8
sigma    = 2e-8; sigma_low    = 9e-9; sigma_high    = 3e-8

# MUON SiPM CONFIG
t_fast   = 2e-8; t_fast_low   = 1e-8; t_fast_high   = 4e-8
t_slow   = 1.2e-6; t_slow_low   = 1e-6; t_slow_high   = 5e-6
amp_fast = 2e-8; amp_fast_low = 1e-8; amp_fast_high = 3e-8
amp_slow = 2e-8; amp_slow_low = 8e-9; amp_slow_high = 9e-8
sigma    = 2e-8; sigma_low    = 9e-9; sigma_high    = 3e-8

fit_initials = (t_fast,t_slow,amp_fast,amp_slow,sigma)
fit_finals = [t_fast,t_slow,amp_fast,amp_slow,sigma]
limits_low = [t_fast_low,t_slow_low,amp_fast_low,amp_slow_low,sigma_low]
limits_high = [t_fast_high,t_slow_high,amp_fast_high,amp_slow_high,sigma_high]
fit_limits = (limits_low,limits_high)

popt, pcov = curve_fit(conv_func2,[laser.wvf_x[:-limit],laser.wvf[:-limit]],alpha.wvf[:-limit], p0 = fit_initials, bounds = fit_limits,method="trf")
perr = np.sqrt(np.diag(pcov))
conv = conv_func2([laser.wvf_x[:-limit],laser.wvf[:-limit]],*popt)
func = func2(alpha.wvf_x,*popt)
conv_int,f_conv,i_conv = signal_int("CONV FUNC",func2(np.arange(0,alpha.wvf_x[-1],5e-10),*popt),timebin,"SiPM","ALL",th = thrld,out = True)

labels = ["TFAST","TSLOW","AFAST","ASLOW","SIGMA"]
print("\n--- FIT VALUES ---")
for i in range(len(fit_initials)):
    fit_finals[i] = popt[i]
    print("%s: %.2E \u00B1 %.2E"%(labels[i],popt[i],perr[i]))
print("------------------\n")

print("SLOW = %.2f%%"%(100*popt[3]/(popt[3]+popt[2])))

########################################################################
#________________________PLOT_FIRST_RESULT_____________________________#
########################################################################

# fig1, axs = plt.subplots(2, 1, sharex=True)

fig1, axs = plt.subplots(2, 1)
plt.title(decon_runs)
fig1.subplots_adjust(hspace=0.25)
fig1.set_figheight(6)
fig1.set_figwidth(6)

# axs[0] = plt.subplot2grid(shape=(1, 1), loc=(0, 0), colspan=3)
# axs[1] = plt.subplot2grid(shape=(1, 1), loc=(3, 0), colspan=3)

axs[0].plot(laser.wvf_x,laser.wvf,label = label_luz)
axs[0].plot(alpha.wvf_x,alpha.wvf,label = label_alp)
axs[0].plot(laser.wvf_x[:-limit],conv,label = "Fitted Convolution")
axs[0].axvline(laser.wvf_x[-limit],color = "grey", ls = ":")
axs[0].set_ylabel("Normalized Amplitude")  
axs[0].axhline(0,color = "grey", ls = ":")
axs[0].set_ylim(1e-4,np.max(alpha.wvf)*1.5)
axs[0].legend()
if logy == True:
    axs[0].semilogy()
    axs[1].semilogy()

axs[1].plot(alpha.wvf_x[np.argmax(alpha.wvf)-np.argmax(func):],func[:np.argmax(func)-np.argmax(alpha.wvf)],label="Convolution Func.")
axs[1].axhline(0,color = "grey", ls = ":")
axs[1].set_ylim(1e-6,10)
axs[1].set_xlabel("Time in [s]"); axs[1].set_ylabel("Convolution signal") 

plt.show()

# output_file.write("%.2E \t\u00B1\t %.2E\n"%(p,perr2[0]))
# output_file.write("%.2E \t\u00B1\t %.2E\n"%(t0,perr1[4]))
output_file.write("%.2E \t\u00B1\t %.2E\n"%(fit_finals[4],perr[4]))
output_file.write("%.2E \t\u00B1\t %.2E\n"%(fit_finals[2],perr[2]))
output_file.write("%.2E \t\u00B1\t %.2E\n"%(fit_finals[0],perr[0]))
output_file.write("%.2E \t\u00B1\t %.2E\n"%(fit_finals[3],perr[3]))
output_file.write("%.2E \t\u00B1\t %.2E"%(fit_finals[1],perr[1]))
output_file.close()
