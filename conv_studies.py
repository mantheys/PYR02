import ROOT

import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit
from ir02_lib.lib import my_wvf,conv_guess

path = '/pnfs/ciemat.es/data/neutrinos/Super-cells_LAr/Feb22/AnalysisROOT/'

signal1 = 'ScintProf_run02_ch6.root'; label1 = "SC SPE Signal"
signal2 = 'ScintProf_run13_ch6.root'; label2 = "SC Alpha Signal"

signals = [signal1,signal2]; labels = [label1,label2]

wvfs = []
for s in signals:
    #Load Raw name
    inRawName = path+s
    inRaw = ROOT.TFile.Open(inRawName ,"READ")
    listkeys_Raw = inRaw.GetListOfKeys()
    wvf = my_wvf("vector",inRawName,listkeys_Raw[0].GetName())
    # print(listkeys_Raw[0].GetName())
    wvfs.append(wvf)

shift_wvfs = [wvfs[0].wvf]
for n in range(1,len(wvfs)):
    if np.argmax(wvfs[n].wvf)-np.argmax(wvfs[0].wvf) > 0:
        shift = np.roll(np.array(wvfs[n].wvf),np.argmax(wvfs[n].wvf)-np.argmax(wvfs[0].wvf))
        shift = np.roll(np.array(wvfs[n].wvf),110)
    else:
        shift = np.roll(np.array(wvfs[n].wvf),np.argmax(wvfs[0].wvf)-np.argmax(wvfs[n].wvf))
        shift = np.roll(np.array(wvfs[n].wvf),110)
    shift_wvfs.append(shift)

t_fast = 1e-9; t_int = 1e-16; t_slow = 1e-6

n_label = 0
for shift in shift_wvfs:
    plt.plot(wvfs[0].wvf_x,shift,label = labels[n_label])
    n_label = n_label+1

plt.xlabel("Time in [s]"); plt.ylabel("Normalized Amplitude")
plt.semilogy()
plt.xlim(0,20e-6); plt.ylim(1e-6, 1.2)

fit_initials = (t_fast,t_int,t_slow,0.999,0.001)
fit_limits = ([1e-10,1e-17,1e-8,1e-4,1e-4],[1e-8,1e-5,1e-5,1,0.5])

popt, pcov = curve_fit(conv_guess,wvfs[0],shift_wvfs[1], p0 = fit_initials, bounds = fit_limits)
conv = conv_guess(wvfs[0],*popt)
print("Tau_fast = %.2e [s] - Amp = %.2e"%(popt[0],popt[3]))
print("Tau_int = %.2e [s] - Amp = %.2e"%(popt[1],1-popt[3]-popt[4]))
print("Tau_slow = %.2e [s]- Amp = %.2e"%(popt[2],popt[4]))
# plt.plot(wvfs[0].wvf_x,conv,label = "Fitted Convolution")
plt.axhline(0,color = "grey", ls= ":")
plt.legend()
plt.show()