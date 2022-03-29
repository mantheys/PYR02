import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit

def lin_func(x,m,n):
    return m*x+n

OV = [30.3,31.3,33.8]
GAIN_0 = [9.7485E+05,9.6374E+05,9.3083E+05,9.7296E+05,9.4447E+05,9.8003E+05,9.6248E+05,9.6116E+05,9.4640E+05,9.4051E+05,9.3022E+05,9.5563E+05]
GAIN_1 = [1.2305E+06,1.2224E+06,1.1930E+06,1.2423E+06,1.2126E+06,1.2505E+06,1.2169E+06,1.2265E+06,1.2182E+06,1.1895E+06,1.1771E+06,1.2173E+06]
GAIN_2 = [1.9029E+06,1.9188E+06,1.8732E+06,1.9414E+06,1.8575E+06,1.9418E+06,1.9051E+06,1.9003E+06,1.8805E+06,1.8600E+06,1.8402E+06,1.8959E+06]

for i in range(12):
    plt.scatter(OV,[GAIN_0[i],GAIN_1[i],GAIN_2[i]],label = "SC  #%i"%(i+1))
    popt, pcov = curve_fit(lin_func, OV, [GAIN_0[i],GAIN_1[i],GAIN_2[i]]); perr = np.sqrt(np.diag(pcov))
    plt.plot([26,27,33.8],[lin_func(26,*popt),lin_func(27,*popt),lin_func(33.8,*popt)])
    print("Breakdown voltage for SC #%i: %.2f"%(i+1,-lin_func(0,*popt)/popt[0]))
# plt.semilogy()
plt.title("Massive SC testing: Set 1")
plt.axhline(0, c = "k", ls = ":", alpha = 0.5)
plt.xlabel("V_op (V)"); plt.ylabel("Gain")
plt.legend();plt.show()

