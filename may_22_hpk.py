import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from scipy.optimize import curve_fit

def lin_func(x,m,n):
    return np.array(m)*x+n

bursts_ov2  = [5.99E+04,3.28E+04,3.84E+04,3.47E+04,1.38E+04,4.83E+04,1.62E+04,5.85E+04,9.29E+04,2.86E+04,1.81E+04,4.34E+04,7.61E+03,7.24E+04,3.91E+04,5.03E+04,4.96E+04,3.07E+04,2.34E+04,2.83E+04,8.46E+03,4.93E+04,1.84E+04,3.79E+04,5.39E+04,1.65E+04,2.84E+04,4.69E+04,1.52E+04,4.79E+04,4.73E+04,3.49E+04,7.26E+04,3.98E+04,3.58E+04,2.22E+04,3.21E+04,2.45E+04,1.63E+04,4.54E+04,2.24E+04,3.95E+04,3.80E+04,1.69E+04,2.77E+04,1.58E+04,1.19E+04,5.13E+04]
bursts_ov25 = [2.88E+04,9.44E+04,4.15E+04,2.51E+04,6.29E+04,5.70E+04,9.86E+04,8.63E+04,4.51E+04,4.37E+04,1.24E+05,2.66E+04,4.53E+04,3.37E+04,1.35E+05,4.68E+04,4.22E+04,3.29E+04,6.81E+04,2.20E+04,4.63E+04,4.22E+04,2.53E+04,5.47E+04,3.04E+04,3.98E+04,2.55E+04,3.07E+04,7.52E+04,4.02E+04,3.09E+04,7.25E+04,1.13E+05,9.16E+04,5.18E+04,6.40E+04,3.19E+04,2.18E+04,1.79E+04,9.57E+04,5.13E+04,3.61E+04,1.32E+05,4.09E+04,2.78E+04,2.20E+04,9.47E+04,2.92E+04,1.33E+04,5.65E+04]
bursts_ov3  = [5.75E+04,6.21E+04,5.27E+04,8.02E+04,6.51E+04,6.35E+04,1.12E+05,9.34E+04,9.55E+04,2.62E+04,1.24E+05,8.23E+04,5.12E+04,6.54E+04,1.27E+05,1.14E+05,4.32E+04,8.63E+04,8.41E+04,6.54E+04,1.10E+05,5.21E+04,5.86E+04,1.58E+05,4.00E+04,5.43E+04,2.51E+04,6.45E+04,6.53E+04,2.91E+04,8.37E+04,5.12E+04,4.34E+04,6.91E+04,5.24E+04,7.23E+04,3.97E+04,5.07E+04,3.17E+04,1.26E+05,5.55E+04,6.20E+04,1.39E+04,4.38E+04,2.35E+04,2.36E+04,4.11E+04,2.38E+04,2.33E+04,6.80E+04]

rate_ov2    = [69.55,47.67,46.53,54.23,39.60,50.54,32.70,47.13,34.62,36.71,66.24,50.47,43.12,61.83,42.75,89.76,32.90,38.09,42.03,32.05,62.33,48.50,38.22,50.55,36.89,47.02,41.05,44.54,44.57,40.53,47.94,44.57,37.81,34.52,49.34,34.91,55.16,36.43,27.48,36.99,27.69,37.54,21.80,30.36,27.37,38.99,37.19,27.99]
rate_ov25   = [60.30,60.58,63.96,46.97,48.61,72.23,52.82,54.93,42.85,29.24,63.41,56.72,73.93,42.95,66.48,68.39,49.45,48.37,36.02,39.28,73.68,66.39,55.50,49.20,34.14,88.33,57.01,60.61,70.00,52.45,62.89,65.63,57.93,51.54,50.58,67.19,47.89,44.06,42.58,52.54,51.73,40.47,51.44,37.97,42.65,50.32,63.32,39.54,47.59,35.53]
rate_ov3    = [73.32,65.99,58.58,76.66,71.21,86.97,58.64,56.79,67.17,38.71,82.28,69.47,61.85,73.10,59.12,74.94,42.28,60.11,61.08,43.32,92.48,64.29,43.73,38.61,43.77,87.02,66.57,88.22,60.41,54.99,74.66,66.77,60.75,60.32,49.12,75.48,68.94,57.87,42.18,53.03,54.78,45.38,75.85,40.08,58.50,48.44,55.00,48.50,45.76,54.22]

mean_bursts_ov2 = np.mean(bursts_ov2)
mean_rate_ov2 = np.mean(rate_ov2)
mean_bursts_ov25 = np.mean(bursts_ov25)
mean_rate_ov25 = np.mean(rate_ov25)
mean_bursts_ov3 = np.mean(bursts_ov3)
mean_rate_ov3 = np.mean(rate_ov3)

fig = plt.figure()
ax = fig.add_subplot(111)

plt.scatter(bursts_ov2,rate_ov2,c = 'tab:blue',alpha=0.5)
popt, pcov = curve_fit(lin_func, bursts_ov2, rate_ov2); perr = np.sqrt(np.diag(pcov))
plt.plot(bursts_ov2,lin_func(bursts_ov2,*popt),c = 'tab:blue',label='OV2 Y=(%.2E$\pm$%.2E)*X + (%.2f$\pm$%.2f)'%(popt[0],perr[0],popt[1],perr[1]))

plt.scatter(bursts_ov25,rate_ov25,c = 'tab:orange',alpha=0.5)
popt, pcov = curve_fit(lin_func, bursts_ov25, rate_ov25); perr = np.sqrt(np.diag(pcov))
plt.plot(bursts_ov25,lin_func(bursts_ov25,*popt),c = 'tab:orange',label='OV2.5 Y=(%.2E$\pm$%.2E)*X + (%.2f$\pm$%.2f)'%(popt[0],perr[0],popt[1],perr[1]))

plt.scatter(bursts_ov3,rate_ov3,c = 'tab:green',alpha=0.5)
popt, pcov = curve_fit(lin_func, bursts_ov3, rate_ov3); perr = np.sqrt(np.diag(pcov))
plt.plot(bursts_ov3,lin_func(bursts_ov3,*popt),c = 'tab:green',label='OV3 Y=(%.2E$\pm$%.2E)*X + (%.2f$\pm$%.2f)'%(popt[0],perr[0],popt[1],perr[1]))

plt.scatter(mean_bursts_ov2,mean_rate_ov2,c='tab:blue',s=100,marker='o')
plt.scatter(mean_bursts_ov25,mean_rate_ov25,c='tab:orange',s=100,marker='o')
plt.scatter(mean_bursts_ov3,mean_rate_ov3,c='tab:green',s=100,marker='o')

plt.xlabel('Bursts Magnitude');plt.ylabel('DCR+B [mH/mm²]')
ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
plt.legend();plt.show()