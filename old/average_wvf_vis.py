import ROOT

import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit
from ir02_lib.deconv_lib import import_scint_prof,my_wvf,tfile_vect2array

# base_dir = '/pnfs/ciemat.es/data/neutrinos/Massive_SC_Testing/Set_5_rep/AnalysisROOT/'
base_dir = ''
path = []; label = []

file = open("WVF_COMPARISON.txt",'r')
lines = file.readlines()
for line in lines:
    names = line.split()
    if names[0] == "#":
        continue
    else:
        path.append(names[0])
        label.append(names[1])
file.close()

check = False; autozoom = False; logy = True; norm = True; fit = False
term_output = check; timebin = 4e-9
SC = []
for i in range(len(path)):
    SC.append(import_scint_prof(base_dir+path[i],timebin,normalize=norm,trim=False,align="True",start=300,cut_i=0,cut_f=0,invert=False))
    plt.plot(SC[i].wvf_x,SC[i].wvf,label = label[i])

plt.grid()
plt.xlabel("Time in s"); plt.ylabel("Amp. in a.u.")
plt.semilogy()
plt.xlim(1e-6,4e-6)
plt.legend()
plt.show()
