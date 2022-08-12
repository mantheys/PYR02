import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import csv
import scipy.constants as cte
from scipy.optimize import root_scalar
from scipy import special

def Peot(d,T):
    #TABLE 1#
    R = 208.133

    a = 89.98
    b = 0.5406E-3
    c = 1.0434E-3

    return (R*T*d) * (1 + (c*d)/(1-b*d)) - a*d**2

def Teot(d,P):
    #TABLE 1#
    R = 208.133
    # P = 1.22e5
    a = 89.98
    b = 0.5406E-3
    c = 1.0434E-3

    return (P + a*d**2)/(R*d *(1 + (c*d)/(1-b*d)))

P = 1.22e5
T = 533
temp = [93,169,223,298]
pres = [2.03E5, 4.05E5, 6.08E5, 7.09E5, 8.10E5, 1.01E6]

t2d = []; p2d = []
for i in range(len(temp)):
    density = root_scalar(lambda d, *args: Teot(d, *args)-temp[i],x0=1.784,x1=5,args=P)
    t2d.append(density.root)

for i in range(len(pres)):
    density = root_scalar(lambda d, *args: Peot(d, *args)-pres[i],x0=1.784,x1=5,args=T)
    p2d.append(density.root)


Ws = np.array([50.6, 46.3, 36.0, 32.1, 28.9, 25.3])

fig, ax = plt.subplots(1,1, figsize = (10,10), sharex=True)

plt.rcParams.update({'font.size': 2})
plt.rc('font', size=12)          # controls default text sizes
plt.rc('axes', titlesize=12)     # fontsize of the axes title
plt.rc('axes', labelsize=12)     # fontsize of the x and y labels
plt.rc('xtick', labelsize=12)    # fontsize of the tick labels
plt.rc('ytick', labelsize=12)    # fontsize of the tick labels
plt.rc('legend', fontsize=12)    # legend fontsize
plt.rc('figure', titlesize=12)   # fontsize of the figure title
plt.rc('axes.formatter', useoffset=False)

plt.suptitle('Absolute Number of Scintillation Photons emitted by alpha particules in GAr')
colors = ['C2', 'C3', 'C4', 'C5', 'C6']
# fit2 = np.poly1d(np.polyfit(p2d, 5.49E6/Ws, 2))
fit2 = np.poly1d(np.polyfit(p2d, 5.49E6/Ws, 1))
# fit3 = np.poly1d(np.polyfit(p2d, 5.49E6/Ws, 3))
xarray = np.linspace(0,10)

ax.scatter(p2d, 5.49E6/Ws, label='Literature (density)')
# ax.scatter(1.78,200000, label='LAr')
ax.plot(xarray, fit2(xarray), c='C1', label='Fit2_Literature')
# ax.plot(xarray, fit3(xarray), '--', c='C2', label='Fit3_Literature')

for i in np.arange(len(temp)):
    ax.scatter(t2d[i], fit2(t2d[i]), label='IR02: ' + str(temp[i]) + ' K', marker='s', s=50, c=colors[i])
    # ax.annotate(r'N$_p$ = ' + str("{:.2e}".format(fit2(t2d[i]))), (t2d[i], fit2(t2d[i])-0.25E6), c=colors[i])
ax.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
ax.set_xlabel(r'Density(kg/m$^3$)')
ax.set_ylabel(r'Number of photons (N$_p$)')


plt.legend()
plt.show()
