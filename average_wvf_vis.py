import ROOT

import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit
from ir02_lib.deconv_lib import import_scint_prof,my_wvf,tfile_vect2array

# base_dir = '/pnfs/ciemat.es/data/neutrinos/Massive_SC_Testing/Set_5_rep/AnalysisROOT/'
base_dir = '/pnfs/ciemat.es/data/neutrinos/Super-cells_LAr/Jan22/AnalysisROOT/'
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

check = False; autozoom = False; logy = True; norm = True; fit = False; thrld = 1e-3; pro_abs = False; pro_rodrigo = False
term_output = check; timebin = 4e-9
SC = []
for i in range(len(path)):
    SC.append(import_scint_prof(base_dir+path[i],timebin,normalize=norm,trim=False,align="True",start=400,cut_i=0,cut_f=0,invert=False))
    plt.plot(SC[i].wvf_x,SC[i].wvf,label = label[i])

italian_path = "/afs/ciemat.es/user/m/manthey/PYR02/ITALIAN_WVF/averaged_waveforms.root"
# italian_path = "/pnfs/ciemat.es/data/neutrinos/Super-cells_LAr/ItalianWVFs/sphe_waveforms_Ch1.root"
italian_wvf = ROOT.TFile.Open(italian_path, "READ")
graph = italian_wvf.Get("average_normalized_ch1")
value = graph.GetY()
n = graph.GetN()
dtime = 4
# print(value)
h = ROOT.TH1D("graph","graph",n,0,n*dtime)
normsphe = -1e12;
for i in range(0,n):
    h.SetBinContent(i+1,*(graph.GetY()+i))
    # tmp = *(gspe1.GetY()+i);
    # if tmp >= normsphe:
        # normsphe = tmp
# for key in italian_wvf.GetListOfKeys():
#     print(key)
# item = italian_wvf.GetListOfKeys()
# print(item)

# file = ROOT.TFile( italian_path, 'READ' )
# item = italian_wvf.Get("averaged_normalized_ch1")

# print(itemx)
# for entry in italian_wvf:
    # print(entry)
# Pitem = np.asarray(item)
# italian_wvf = np.array(tfile_vect2array(italian_path,item[1].GetName()))
# print(item.GetHistogram("p"))

# print(item.GetObjectInfo(1,1))

# it_wvf = []
# for i in range(len(italian_wvf.wvf)):
#     it_wvf.append(italian_wvf.wvf[i])

# plt.plot(it_wvf)
plt.grid()
plt.xlabel("Time in s"); plt.ylabel("ADC counts")
plt.xlim(1e-6,4e-6)
plt.legend()
plt.show()
