import ROOT

import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit
from ir02_lib.deconv_lib import import_scint_prof,my_wvf,conv_guess2, conv_guess3

base_dir = '/pnfs/ciemat.es/data/neutrinos/Massive_SC_Testing/Set_5_rep/AnalysisROOT/'
path = []; label = []

file = open("SC_SET5_REP.txt",'r')
lines = file.readlines()
for line in lines:
    names = line.split()
    path.append(names[0])
    label.append(names[1])
file.close()

check = False; autozoom = False; logy = True; norm = False; fit = False; thrld = 1e-3; pro_abs = False; pro_rodrigo = False
term_output = check; timebin = 4e-9
SC = []
for i in range(len(path)):
    SC.append(import_scint_prof(base_dir+path[i],timebin,normalize=norm,trim=False,align=False,start=0,cut_i=0,cut_f=0))
    plt.plot(SC[i].wvf_x,SC[i].wvf,label = label[i])

plt.xlabel("Time in s")
plt.ylabel("ADC counts")
plt.legend()
plt.show()
