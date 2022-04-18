import ROOT
import scipy
import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit
from ir02_lib.deconv_lib import import_scint_prof,pdf,func,my_wvf,signal_int,import_deconv_runs

def spe_rise(x,x_start,a,t):
    y = a*np.exp((x+x_start)/t)
    return y

space = 5
init = 200

timebin = 4e-9; norm = False; trm = False; shift = False; s_start = 0;
path = "/pnfs/ciemat.es/data/neutrinos/Super-cells_LAr/Dic21/AnalysisROOT/run63_ScintProfFirstSignalBin_ch4.root"
spe = import_scint_prof(path,timebin,normalize=norm,trim=trm,align=shift,start=s_start,cut_i=0,cut_f=0)

maxim = np.argmax(spe.wvf)-10

popt, pcov = curve_fit(spe_rise,spe.wvf_x[0:init-space],spe.wvf[maxim-init:maxim-space], p0=[])
fit = spe_rise(spe.wvf_x[0:init-space],*popt)

# plt.plot(spe.wvf_x[maxim-init:],spe_fit(spe.wvf_x[0:-(maxim-init)],init,maxim,1,1e-2,1e-3))
plt.plot(spe.wvf_x,spe.wvf)
plt.plot(spe.wvf_x[maxim-init:maxim-space],fit)
plt.xlabel("Time in [s]"); plt.ylabel("Amplitude in ADC counts")
plt.show()