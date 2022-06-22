import ROOT
import numpy as np
import scipy as sc
import matplotlib.pyplot as plt

from scipy.special import erf
from scipy.optimize import curve_fit
from ir02_lib.deconv_lib import import_scint_prof,pdf,func,my_wvf,signal_int,import_deconv_runs

def gauss(T,SIGMA):
    return np.exp(-1/2*(np.array(T)/SIGMA)**2)

def func(T,P,A,SIGMA,TAU,T0):
    # P = 5e-6
    return P+A*np.exp((SIGMA/(np.sqrt(2)*TAU))**2-(np.array(T)-T0)/TAU)*(1-erf((SIGMA**2-TAU*(np.array(T)-T0))/(np.sqrt(2)*SIGMA*TAU)))

def func2(T,P,A1,SIGMA,TAU1,T0,A2,TAU2):
    return func(T,P,A1,SIGMA,TAU1,T0) + func(T,P,A2,SIGMA,TAU2,T0)

def convfunc2(T,A1,SIGMA,TAU1,T0,A2,TAU2):
    return np.convolve(gauss(T,SIGMA),func(T,A1,SIGMA,TAU1,T0) + func(T,A2,SIGMA,TAU2,T0),'same')

def func3(T,A1,SIGMA,TAU1,T0,A2,TAU2,A3,TAU3):
    return func(T,A1,SIGMA,TAU1,T0) + func(T,A2,SIGMA,TAU2,T0) + func(T,A3,SIGMA,TAU3,T0)

def logfunc2(T,P,A1,SIGMA,TAU1,T0,A2,TAU2):
    return np.log(func(T,P,A1,SIGMA,TAU1,T0) + func(T,P,A2,SIGMA,TAU2,T0))

def logfunc3(T,A1,SIGMA,TAU1,T0,A2,TAU2,A3,TAU3):
    # p = 1e-6/3
    return np.log(func(T,A1,SIGMA,TAU1,T0) + func(T,A2,SIGMA,TAU2,T0) + func(T,A3,SIGMA,TAU3,T0))

# path = "/pnfs/ciemat.es/data/neutrinos/Super-cells_LAr/Feb22_2/Deconvolución/"
path = "/pnfs/ciemat.es/data/neutrinos/Super-cells_LAr/Feb22_2/AnalysisROOT/"
# file = "FEB_2_SiPM_DAY1_OV3_f100.root"
file = "run104_ScintProfFirstSignalBin_ch4.root"

signal = import_scint_prof(path+file,timebin = 4e-9,normalize=1,trim=0,align=False,start=0,cut_i=0,cut_f=0,invert=False)

for i in range(len(signal.wvf)):
    if signal.wvf[i] <= 1e-10:
        signal.wvf[i] = 1e-10
    if np.isnan(signal.wvf[i]):
        signal.wvf[i] = 1e-10

max_in = np.argmax(signal.wvf)

buffer1 = 50; buffer2 = 3500
guess_t0 = signal.wvf_x[np.argmax(signal.wvf)-10]

pedestal = np.mean(signal.wvf[-buffer2:])
# print(pedestal)

p = pedestal; p_low = pedestal*0.1; p_high = pedestal*10
a1 = 1e-0; a1_low = 1e-3; a1_high = 10
# a2 = 1e-2; a2_low = 1e-4; a2_high = 1
a3 = 5e-3; a3_low = 1e-4; a3_high = 1e-0
t0 = guess_t0; t0_low = guess_t0*0.1; t0_high = guess_t0*10
tau1 = 1e-8; tau1_low = 5e-10; tau1_high = 4e-8
# tau2 = 5e-8; tau2_low = tau1_high; tau2_high = 9e-7
tau3 = 1e-6; tau3_low = 5e-7; tau3_high = 9e-6
sigma = 1e-8; sigma_low = 9e-10; sigma_high = 1e-7

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

# EXPORT FIT PARAMETERS
p = popt[0];a1 = popt[1];sigma = popt[2];tau1 = popt[3];t0 = popt[4]

# CHECK FIRST FIT
plt.plot(signal.wvf_x,signal.wvf)
plt.plot(signal.wvf_x[max_in-buffer1:max_in+int(buffer1/2)],func(signal.wvf_x[max_in-buffer1:max_in+int(buffer1/2)],*popt))
plt.semilogy()
plt.show()

# USING VALUES FROM FIRST FIT PERFORM SECONG FIT FOR THE SLOW COMPONENT
bounds2 = ([p_low,a3_low,tau3_low],[p_high,a3_high,tau3_high])
initial2 = (p/2,a3,tau3)
labels2 = ["PED","AMP","TAU"]

popt, pcov = curve_fit(lambda T,P,A3,TAU3: logfunc2(T,P,a1,sigma,tau1,t0,A3,TAU3),signal.wvf_x[max_in-buffer1:-buffer2],np.log(signal.wvf[max_in-buffer1:-buffer2]),p0 = initial2, bounds = bounds2, method = "trf")
perr2 = np.sqrt(np.diag(pcov))

print("\n--- SECOND FIT VALUES (SLOW) ---")
for i in range(len(initial2)):
    print("%s: \t%.2E \u00B1 %.2E"%(labels2[i],popt[i],perr2[i]))

# EXPORT FIT PARAMETERS
p = popt[0];a3 = popt[1];tau3 = popt[2]
param = [p,a1,sigma,tau1,t0,a3,tau3]
print("\nAMP RATIO %.2f:100"%(100*a3/a1))

# EXPORT FOR EXCEL FORMAT
print("%.2E \t\u00B1\t %.2E"%(p,perr2[0]))
print("%.2E \t\u00B1\t %.2E"%(t0,perr1[4]))
print("%.2E \t\u00B1\t %.2E"%(sigma,perr1[2]))
print("%.2E \t\u00B1\t %.2E"%(a1,perr1[1]))
print("%.2E \t\u00B1\t %.2E"%(tau1,perr1[3]))
print("%.2E \t\u00B1\t %.2E"%(a3,perr2[1]))
print("%.2E \t\u00B1\t %.2E"%(tau3,perr2[2]))

plt.plot(signal.wvf_x,signal.wvf)
plt.plot(signal.wvf_x[max_in-buffer1:-buffer2],func2(signal.wvf_x[max_in-buffer1:-buffer2],*param))
plt.semilogy()
plt.show()

