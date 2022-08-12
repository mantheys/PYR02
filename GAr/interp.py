import scipy 
import scipy.interpolate as inter
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

with open("GAr_DATA/sipm375_data.txt",'r') as file:
    wvl = []; pde = []
    while True:
        next_line = file.readline()
        if not next_line:
            break
        else:
            line = next_line.split("\t")
            try:
                wvl.append(int(line[0]))
                pde.append(float(line[1]))
            except:
                continue

wvl = np.array(wvl)
pde = np.array(pde)
y_int = scipy.interpolate.interp1d(wvl, pde)

def lin_func(x,m,n):
    return m*x+n

def gain_interp(ov,gain,value):
    popt, pcov = curve_fit(lin_func,ov,gain)
    perr = np.sqrt(np.diag(pcov))
    offset = -popt[1]/popt[0]     
    return lin_func(value+offset,*popt)

gain = [4.44E+06,6.33E+06,8.23E+06]
ov = [3,4,5]

print(gain_interp(ov,gain,4))

critical_wvl = [127,270,317,340,365,385,470,570]
critical_pde = np.array(y_int(critical_wvl))

sipm375 =  pd.read_csv('GAr_DATA/sipm375_data.txt', delimiter='\t', header = 0)
# jan1_pde = pd.read_csv('GAr_DATA/jan1_pde.txt', delimiter='\t',header = 0)
# jan2_pde = pd.read_csv('GAr_DATA/jan2_pde.txt', delimiter='\t',header = 0)
# feb1_pde = pd.read_csv('GAr_DATA/feb1_pde.txt', delimiter='\t',header = 0)
# feb2_pde = pd.read_csv('GAr_DATA/feb2_pde.txt', delimiter='\t',header = 0)
# apri_pde = pd.read_csv('GAr_DATA/apri_pde.txt', delimiter='\t',header = 0)
# mean_pde = pd.read_csv('GAr_DATA/mean_pde.txt', delimiter='\t',header = 0)
gar1_pde = pd.read_csv('GAr_DATA/GAr1_pde.txt', delimiter='\t',header = 0)
gar3_pde = pd.read_csv('GAr_DATA/GAr3_pde.txt', delimiter='\t',header = 0)

fig, ax = plt.subplots(1,1, figsize = (16,8))
plt.plot(sipm375["Wavelength"], sipm375["375.00"],'o', color = "dimgray", label = "RT")
for i in [127,270,317,340,365,385,470,570]:
    plt.axvline(i, linestyle = "--", color = "grey", alpha = 0.5)

plt.errorbar(critical_wvl, critical_pde,0,5,'s', markersize = 9, color = "black")
# plt.errorbar(gar1_pde["Wavelength"][gar1_pde["Temperatura"]=="-181"], gar1_pde["PDE_375 (LN2)"][gar1_pde["Temperatura"]=="-181"],
#              gar1_pde["DPDE_375 (LN2)"][gar1_pde["Temperatura"]=="-181"],5,'o', markersize = 8, color= "dodgerblue", label = "LN2 - RT")
plt.errorbar(gar3_pde["Wavelength"][gar3_pde["Temperatura"]=="-180"], gar3_pde["PDE_375 (LN2)"][gar3_pde["Temperatura"]=="-180"],
             gar3_pde["DPDE_375 (LN2)"][gar3_pde["Temperatura"]=="-180"],5,'o', markersize = 8, color= "orange", label = "SiPM1 LN2 - RT")
plt.errorbar(gar3_pde["Wavelength"][gar3_pde["Temperatura"]=="-181"], gar3_pde["PDE_375 (LN2)"][gar3_pde["Temperatura"]=="-181"],
             gar3_pde["DPDE_375 (LN2)"][gar3_pde["Temperatura"]=="-181"],5,'o', markersize = 8, color= "red", label = "SiPM2 LN2 - RT")
# plt.errorbar(gar3_pde["Wavelength"][gar3_pde["Temperatura"]=="-145"], gar3_pde["PDE_375 (LN2)"][gar3_pde["Temperatura"]=="-145"],
#              gar3_pde["DPDE_375 (LN2)"][gar3_pde["Temperatura"]=="-145"],5,'o', markersize = 8, color= "green",alpha =0.5, label = "CT(-145) - RT")
# plt.errorbar(gar1_pde["Wavelength"][gar1_pde["Temperatura"]=="-104"], gar1_pde["PDE_375 (LN2)"][gar1_pde["Temperatura"]=="-104"],
#              gar1_pde["DPDE_375 (LN2)"][gar1_pde["Temperatura"]=="-104"],5,'o', markersize = 8, color= "green", label = "CT(-104) - RT")
# plt.errorbar(gar1_pde["Wavelength"][gar1_pde["Temperatura"]=="-50"], gar1_pde["PDE_375 (LN2)"][gar1_pde["Temperatura"]=="-50"],
#              gar1_pde["DPDE_375 (LN2)"][gar1_pde["Temperatura"]=="-50"],5,'o', markersize = 8, color= "orange", label = "CT(-50) - RT")
ax.legend(fontsize = 16)
plt.grid(axis="y",linewidth = 0.5)
plt.grid(axis="x",linewidth = 0.5)

plt.ylim(0,40); plt.ylabel("PDE (%)", size = 20); plt.yticks(size = 20)
plt.xlabel("Wavelength (nm)", size = 20); plt.xticks(np.arange(100, 950, step=50),size = 20)
fig.savefig("myplot.png")
# plt.show()
