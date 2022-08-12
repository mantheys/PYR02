import ROOT
import numpy as np
import scipy as sc
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit
from ir02_lib.deconv_lib import import_scint_prof,signal_int
from ir02_lib.fit_lib import funcP,func2sigma,logfunc2sigma

print("\n### WELCOME TO THE SCINTILLATION FIT ###")

path = "output/"
# path = "/pnfs/ciemat.es/data/neutrinos/Super-cells_LAr/Jan22/AnalysisROOT/"
file = "GAr_SiPM1_DAY21_OV2_T-180_f800"
# file = "run80_ScintProfFirstSignalBin_PMT"
signal = import_scint_prof(path+file+".root",timebin = 4e-9,normalize=1,trim=0,align=False,start=0,cut_i=0,cut_f=0,invert=False)
signal.apply_smooth(0.6)
# path = "deconvolution_input/LAr/JAN/"
# file = "JAN_SiPM_DAY2_OV2_MUON"
# output_file = open("scint_fit/"+file+"_CONV.txt","w")
# detector,timebin,int_st,paths,filename,shift,filter_strenght,reverse,fnal,s_start,particle_label = import_deconv_runs(path+file+".txt",debug = True)
# init=0; trm = False; #fnal=500; s_start = 600;
# signal = import_scint_prof(paths[0]+paths[1],timebin,normalize=True,trim=trm,align=shift,start=s_start,cut_i=init,cut_f=fnal,invert=False)

buffer1     = 25
buffer2     = 5000
check       = True
autozoom    = True
logy        = True
thrld       = 1e-2
timebin     = 4e-9
detector    = "SiPM"
term_output = check
# signal.apply_smooth(0.2)

for i in range(len(signal.wvf)):
    if signal.wvf[i] <= 1e-10:
        signal.wvf[i] = 1e-10
    if np.isnan(signal.wvf[i]):
        signal.wvf[i] = 1e-10

output_file = open("output/"+file+"_SCINT_FIT.txt","w")
max_in = np.argmax(signal.wvf)
guess_t0 = signal.wvf_x[np.argmax(signal.wvf)-10]

if buffer2 < max_in:
    buffer2 = max_in+buffer2

pedestal = np.mean(signal.wvf[:max_in-buffer1])
print("Estimated pedestal %.2E"%pedestal)

p = pedestal; p_low = pedestal*0.01; p_high = pedestal*100
a1 = 2e-8; a1_low = 2e-9; a1_high = 5e-8                    
a3 = 8e-9; a3_low = 1e-9; a3_high = 5e-8
t0 = guess_t0; t0_low = guess_t0*0.1; t0_high = guess_t0*10
tau1 = 8e-9; tau1_low = 6e-9; tau1_high = 2e-8
tau3 = 1e-6; tau3_low = tau1_high; tau3_high = 9e-5
sigma = 2e-8; sigma_low = 8e-9; sigma_high = 6e-8

bounds = ([p_low,a1_low,sigma_low,tau1_low,t0_low],[p_high,a1_high,sigma_high,tau1_high,t0_high])
initial = (p,a1,sigma,tau1,t0)
labels = ["PED","AMP","SIGMA","TAU","TIME"]

# FIT PEAK
popt, pcov = curve_fit(lambda T,P,A1,SIGMA,TAU1,T0: funcP(T,P,A1,SIGMA,TAU1,T0),signal.wvf_x[max_in-buffer1:max_in+int(buffer1/2)],signal.wvf[max_in-buffer1:max_in+int(buffer1/2)],p0 = initial, bounds = bounds, method = "trf")
perr1 = np.sqrt(np.diag(pcov))

# PRINT FIRST FIT VALUE
print("\n--- FISRT FIT VALUES (FAST) ---")
for i in range(len(initial)):
    print("%s: \t%.2E \u00B1 %.2E"%(labels[i],popt[i],perr1[i]))
print("-------------------------------")

# EXPORT FIT PARAMETERS
p = popt[0];a1 = popt[1];sigma = popt[2];tau1 = popt[3];t0 = popt[4]
sigma2 = sigma*0.1; sigma2_low = sigma*0.01; sigma2_high = sigma

# CHECK FIRST FIT
plt.subplot(1,2,1)
plt.plot(signal.wvf_x,signal.wvf)
plt.plot(signal.wvf_x[max_in-buffer1:max_in+int(buffer1/2)],funcP(signal.wvf_x[max_in-buffer1:max_in+int(buffer1/2)],*popt))
plt.axvline(signal.wvf_x[buffer2],ls = ":",c = "k",alpha=0.5)
plt.xlabel("Time in s"); plt.ylabel("Amp in a.u.")
plt.semilogy();plt.ylim(1e-6,np.max(signal.wvf)*1.2)

plt.subplot(1,2,2)
plt.plot(signal.wvf_x,signal.wvf)
plt.plot(signal.wvf_x[max_in-buffer1:max_in+int(buffer1/2)],funcP(signal.wvf_x[max_in-buffer1:max_in+int(buffer1/2)],*popt))
plt.axvline(signal.wvf_x[buffer2],ls = ":",c = "k",alpha=0.5)
plt.xlabel("Time in s")
plt.semilogy(); plt.xlim(signal.wvf_x[max_in-50],signal.wvf_x[max_in+250]);plt.ylim(1e-6,np.max(signal.wvf)*1.2)

plt.show()

# USING VALUES FROM FIRST FIT PERFORM SECONG FIT FOR THE SLOW COMPONENT

bounds2 = ([p_low,sigma2_low,a3_low,tau3_low],[p_high,sigma2_high,a3_high,tau3_high])
initial2 = (p/2,sigma2,a3,tau3)
labels2 = ["PED","SIGMA","AMP","TAU"]
popt, pcov = curve_fit(lambda T,P,SIGMA2,A3,TAU3: logfunc2sigma(T,P,a1,SIGMA2,SIGMA2,tau1,t0,A3,TAU3),signal.wvf_x[max_in-buffer1:buffer2],np.log(signal.wvf[max_in-buffer1:buffer2]),p0 = initial2, bounds = bounds2, method = "trf")
perr2 = np.sqrt(np.diag(pcov))
p = popt[0];sigma2 = popt[1];a3 = popt[2];tau3 = popt[3]
param = [p,a1,sigma2,sigma2,tau1,t0,a3,tau3]

output_file.write("%.2E \t\u00B1\t %.2E\n"%(p,perr2[0]))
output_file.write("%.2E \t\u00B1\t %.2E\n"%(t0,perr1[4]))
output_file.write("%.2E \t\u00B1\t %.2E\n"%(sigma2,perr2[1]))
# output_file.write("%.2E \t\u00B1\t %.2E\n"%(sigma,perr1[2]))
output_file.write("%.2E \t\u00B1\t %.2E\n"%(a1,perr1[1]))
output_file.write("%.2E \t\u00B1\t %.2E\n"%(tau1,perr1[3]))
output_file.write("%.2E \t\u00B1\t %.2E\n"%(a3,perr2[2]))
output_file.write("%.2E \t\u00B1\t %.2E"%(tau3,perr2[3]))
# output_file.write("%.2E \t\u00B1\t %.2E\n"%(sigma2,perr2[1]))

print("\n--- SECOND FIT VALUES (SLOW) ---")
for i in range(len(initial2)):
    print("%s: \t%.2E \u00B1 %.2E"%(labels2[i],popt[i],perr2[i]))
print("--------------------------------\n")

def fit(X_ARRAY):
    return func2sigma(X_ARRAY,*param)

int_sig_mix,f_sig_mix,i_sig_mix = signal_int("SIGNAL",signal.wvf,timebin,detector,"MIXED",th = thrld,out = term_output)
int_fit_mix,f_fit_mix,i_fit_mix = signal_int("FIT",fit(signal.wvf_x),timebin,detector,"MIXED",th = thrld,out = term_output)
print("\nCharge increase after fit:\n %.2f%% (integration type 'MIXED')\n"%(100*(int_fit_mix/int_sig_mix-1)))

fig,ax0 = plt.subplots(ncols=2)
plt.title("Scintillation fit for %s"%file)
plt.subplot(1,2,1)
plt.plot(signal.wvf_x,signal.wvf,zorder=0,c = "tab:blue")
plt.plot(signal.wvf_x,fit(signal.wvf_x),c = "tab:orange")
plt.plot(signal.wvf_x[:buffer2],fit(signal.wvf_x[:buffer2]),c = "tab:green")
plt.xlabel("Time in s"); plt.ylabel("Amp in a.u.")
plt.axvline(signal.wvf_x[buffer2],ls = ":",c = "k",alpha=0.5,label = "FIT endpoint")
plt.axvline(f_sig_mix,ls = "--",c = "k",label = "SIGNAL integration limits")
plt.axvline(i_sig_mix,ls = "--",c = "k")
plt.axvline(f_fit_mix,ls = "--",c = "b",label = "FIT integration limits")
plt.axvline(i_fit_mix,ls = "--",c = "b")
plt.semilogy();plt.ylim(1e-6,np.max(signal.wvf)*1.2)
plt.legend()

plt.subplot(1,2,2)
plt.plot(signal.wvf_x,signal.wvf,zorder=0,c = "tab:blue")
plt.plot(signal.wvf_x,fit(signal.wvf_x),c = "tab:orange")
plt.plot(signal.wvf_x[:buffer2],fit(signal.wvf_x[:buffer2]),c = "tab:green")
plt.xlabel("Time in s")
plt.axvline(signal.wvf_x[buffer2],ls = ":",c = "k",alpha=0.5,label = "FIT endpoint")
plt.axvline(f_sig_mix,ls = "--",c = "k",label = "SIGNAL integration limits")
plt.axvline(i_sig_mix,ls = "--",c = "k")
plt.axvline(f_fit_mix,ls = "--",c = "b",label = "FIT integration limits")
plt.axvline(i_fit_mix,ls = "--",c = "b")
plt.semilogy(); plt.xlim(signal.wvf_x[max_in-50],signal.wvf_x[max_in+250]);plt.ylim(1e-6,np.max(signal.wvf)*1.2)

plt.show()
output_file.close()
print("Generated output file: %s"%(file+"_SCINT_FIT.txt"))
