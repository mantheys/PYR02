import ROOT
import numpy as np
import scipy as sc
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit
from ir02_lib.deconv_lib import import_scint_prof,pdf,func,my_wvf,signal_int,import_deconv_runs
from ir02_lib.fit_lib import func,func2,logfunc2,func2sigma,logfunc2sigma

path = "/pnfs/ciemat.es/data/neutrinos/Super-cells_LAr/Jan22/Deconvolución/"
file = "JAN_SC_DAY1_OV2_MUON_f150"

# path = "/pnfs/ciemat.es/data/neutrinos/Eficiencia/GAr_Julio/AnalysisROOT/"
# file = "ScintProfFirstSignalBin_run02_ch0"

buffer1 = 40; buffer2 = 3900; double_sigma = False
signal = import_scint_prof(path+file+".root",timebin = 4e-9,normalize=1,trim=0,align=False,start=0,cut_i=0,cut_f=0,invert=False)
output_file = open("scint_fit/"+file+"_FIT.txt","w")

for i in range(len(signal.wvf)):
    if signal.wvf[i] <= 1e-10:
        signal.wvf[i] = 1e-10
    if np.isnan(signal.wvf[i]):
        signal.wvf[i] = 1e-10

max_in = np.argmax(signal.wvf)
guess_t0 = signal.wvf_x[np.argmax(signal.wvf)-10]

pedestal = np.mean(signal.wvf[-buffer2:])
# print(pedestal)

p = pedestal; p_low = pedestal*0.1; p_high = pedestal*10

a1 = 1e-8; a1_low = 5e-10; a1_high = 9e-7                    
a3 = 9e-8; a3_low = 5e-10; a3_high = 9e-7

t0 = guess_t0; t0_low = guess_t0*0.1; t0_high = guess_t0*10

tau1 = 6e-9; tau1_low = 1e-9; tau1_high = 5e-8
tau3 = 1e-6; tau3_low = 1e-7; tau3_high = 5e-5

sigma = 2e-8; sigma_low = 1e-9; sigma_high = 9e-8
# sigma2 = 2e-8; sigma2_low = 1e-9; sigma2_high = 9e-8

bounds = ([p_low,a1_low,sigma_low,tau1_low,t0_low],[p_high,a1_high,sigma_high,tau1_high,t0_high])
initial = (p,a1,sigma,tau1,t0)
labels = ["PED","AMP","SIGMA","TAU","TIME"]


# FIT PEAK
popt, pcov = curve_fit(lambda T,P,A1,SIGMA,TAU1,T0: func(T,P,A1,SIGMA,TAU1,T0),signal.wvf_x[max_in-buffer1:max_in+int(buffer1/2)],signal.wvf[max_in-buffer1:max_in+int(buffer1/2)],p0 = initial, bounds = bounds, method = "trf")
perr1 = np.sqrt(np.diag(pcov))

# PRINT FIRST FIT VALUE
print("\n--- FISRT FIT VALUES (FAST) ---")
for i in range(len(initial)):
    print("%s: \t%.2E \u00B1 %.2E"%(labels[i],popt[i],perr1[i]))
print("-------------------------------\n")

# EXPORT FIT PARAMETERS
p = popt[0];a1 = popt[1];sigma = popt[2];tau1 = popt[3];t0 = popt[4]

# CHECK FIRST FIT
plt.plot(signal.wvf_x,signal.wvf)
plt.plot(signal.wvf_x[max_in-buffer1:max_in+int(buffer1/2)],func(signal.wvf_x[max_in-buffer1:max_in+int(buffer1/2)],*popt))
plt.axvline(signal.wvf_x[-buffer2],ls = "--",c = "k")
plt.xlabel("Time in s"); plt.ylabel("Amp in a.u.")
plt.semilogy()
plt.show()

# USING VALUES FROM FIRST FIT PERFORM SECONG FIT FOR THE SLOW COMPONENT
if double_sigma == False:
    bounds2 = ([p_low,a3_low,tau3_low],[p_high,a3_high,tau3_high])
    initial2 = (p/2,a3,tau3)
    labels2 = ["PED","AMP","TAU"]
    popt, pcov = curve_fit(lambda T,P,A3,TAU3: logfunc2(T,P,a1,sigma,tau1,t0,A3,TAU3),signal.wvf_x[max_in-buffer1:-buffer2],np.log(signal.wvf[max_in-buffer1:-buffer2]),p0 = initial2, bounds = bounds2, method = "trf")
    perr2 = np.sqrt(np.diag(pcov))
    p = popt[0];a3 = popt[1];tau3 = popt[2]
    param = [p,a1,sigma,tau1,t0,a3,tau3]
    plt.plot(signal.wvf_x[max_in-buffer1:-buffer2],func2(signal.wvf_x[max_in-buffer1:-buffer2],*param))
    output_file.write("%.2E \t\u00B1\t %.2E\n"%(p,perr2[0]))
    output_file.write("%.2E \t\u00B1\t %.2E\n"%(t0,perr1[4]))
    output_file.write("%.2E \t\u00B1\t %.2E\n"%(sigma,perr1[2]))
    output_file.write("%.2E \t\u00B1\t %.2E\n"%(a1,perr1[1]))
    output_file.write("%.2E \t\u00B1\t %.2E\n"%(tau1,perr1[3]))
    output_file.write("%.2E \t\u00B1\t %.2E\n"%(a3,perr2[1]))
    output_file.write("%.2E \t\u00B1\t %.2E"%(tau3,perr2[2]))

else:
    bounds2 = ([p_low,sigma2_low,a3_low,tau3_low],[p_high,sigma2_high,a3_high,tau3_high])
    initial2 = (p/2,sigma2,a3,tau3)
    labels2 = ["PED","SIGMA","AMP","TAU"]
    popt, pcov = curve_fit(lambda T,P,SIGMA2,A3,TAU3: logfunc2sigma(T,P,a1,sigma,SIGMA2,tau1,t0,A3,TAU3),signal.wvf_x[max_in-buffer1:-buffer2],np.log(signal.wvf[max_in-buffer1:-buffer2]),p0 = initial2, bounds = bounds2, method = "trf")
    perr2 = np.sqrt(np.diag(pcov))
    p = popt[0];sigma2 = popt[1];a3 = popt[2];tau3 = popt[3]
    param = [p,a1,sigma,sigma2,tau1,t0,a3,tau3]
    plt.plot(signal.wvf_x[max_in-buffer1:-buffer2],func2sigma(signal.wvf_x[max_in-buffer1:-buffer2],*param))
    output_file.write("%.2E \t\u00B1\t %.2E\n"%(p,perr2[0]))
    output_file.write("%.2E \t\u00B1\t %.2E\n"%(t0,perr1[4]))
    output_file.write("%.2E \t\u00B1\t %.2E\n"%(sigma,perr1[2]))
    output_file.write("%.2E \t\u00B1\t %.2E\n"%(sigma2,perr2[1]))
    output_file.write("%.2E \t\u00B1\t %.2E\n"%(a1,perr1[1]))
    output_file.write("%.2E \t\u00B1\t %.2E\n"%(tau1,perr1[3]))
    output_file.write("%.2E \t\u00B1\t %.2E\n"%(a3,perr2[2]))
    output_file.write("%.2E \t\u00B1\t %.2E"%(tau3,perr2[3]))

print("\n--- SECOND FIT VALUES (SLOW) ---")
for i in range(len(initial2)):
    print("%s: \t%.2E \u00B1 %.2E"%(labels2[i],popt[i],perr2[i]))
print("--------------------------------\n")

plt.plot(signal.wvf_x,signal.wvf,zorder=0)
plt.xlabel("Time in s"); plt.ylabel("Amp in a.u.")
plt.axvline(signal.wvf_x[-buffer2],ls = "--",c = "k")
plt.semilogy()
plt.show()
output_file.close()