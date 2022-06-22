import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit

def lin_func(x,m):
    return m*x

x = [0,800]; y = x

OV3_x = [643.3701814,539.2085491,422.050839,358.9362242,280.5093991,229.2327956,82.54439428,49.86118906,27.18239017,11.58066327,1.544088435,0.2133415522]
OV3 = [615.495,559.812,464.226,387.355,337.853,236.035,82.667,52.679,27.898,10.566,1.659,0.197]
DOV3 = [67.74112912,56.77384031,44.4381807,37.79277593,29.5351323,24.1361643,8.691186866,5.249937507,2.862062705,1.219340324,0.1625787099,0.02246295842]
popt_OV3, pcov_OV3 = curve_fit(lin_func, OV3_x, OV3); perr_OV3 = np.sqrt(np.diag(pcov_OV3))
print("OV 3.0 Fitted slope m = %.3f with ∆m = %.3f"%(popt_OV3,perr_OV3))

OV2_5_x = [587.9166733,492.7329639,385.6733377,327.998712,256.3316695,209.4747107,75.42970918,45.56354218,24.83947945,10.58250012,1.411000016,0.1949531689]
OV2_5 = [550.152,508.169,422.611,352.359,307.511,213.474,74.818,47.624,25.228,9.689,1.508,0.182]
DOV2_5 = [57.11142227,47.865083,37.46509301,31.86246248,24.90058011,20.34879975,7.327395477,4.426135226,2.412957591,1.028005601,0.1370674134,0.01893814762]
popt_OV2_5, pcov_OV2_5 = curve_fit(lin_func, OV2_5_x, OV2_5); perr_OV2_5 = np.sqrt(np.diag(pcov_OV2_5))
print("OV 2.5 Fitted slope m = %.3f with ∆m = %.3f"%(popt_OV2_5,perr_OV2_5))

OV2_x = [523.5457617,438.7837029,343.4460197,292.0861805,228.2659521,186.5393549,67.17092123,40.57479653,22.11980843,9.423823711,1.256509828,0.1736077746]
OV2 = [484.569,449.846,375.044,310.922,272.077,188.775,66.489,42.041,22.197,8.456,1.356,0.160]
DOV2 = [58.49997131,49.02882596,38.37598118,32.63713399,25.50598749,20.84353978,7.505546319,4.533747777,2.471623788,1.052999484,0.1403999311,0.01939859049]
popt_OV2, pcov_OV2 = curve_fit(lin_func, OV2_x, OV2); perr_OV2 = np.sqrt(np.diag(pcov_OV2))
print("OV 2.0 Fitted slope m = %.3f with ∆m = %.3f"%(popt_OV2,perr_OV2))

plt.title("SC Linearity measured with Laser (405nm)")
plt.ylabel("Measured PE");plt.xlabel("Theoretical PE")
plt.plot(x, y, label = "Ref. (Measured PE = Theoretical PE)",alpha = 0.5, color = "k")

plt.plot(OV3_x,lin_func(OV3_x,popt_OV3),c='tab:blue')
plt.plot(OV2_5_x,lin_func(OV2_5_x,popt_OV2_5),c='tab:orange')
plt.plot(OV2_x,lin_func(OV2_x,popt_OV2),c='tab:green')

plt.scatter(OV3_x, OV3, label = "Measured PE (OV = 3.0V) Fitted slope %.3f"%popt_OV3)
plt.errorbar(OV3_x, OV3, xerr=DOV3, fmt="o")

plt.scatter(OV2_5_x, OV2_5, label = "Measured PE (OV = 2.5V) Fitted slope %.3f"%popt_OV2_5)
plt.errorbar(OV2_5_x, OV2_5, xerr=DOV2_5, fmt="o")

plt.scatter(OV2_x, OV2, label = "Measured PE (OV = 2.0V) Fitted slope %.3f"%popt_OV2)
plt.errorbar(OV2_x, OV2, xerr=DOV2, fmt="o")

# plt.scatter(RUN2_OV3[0],RUN2_OV3[3],label = "Recovered Data 1");
# plt.scatter(RUN2_OV3[1],RUN2_OV3[4],label = "Recovered Data 2");
# plt.scatter(RUN2_OV3[2],RUN2_OV3[5],label = "Recovered Data 3"); 
# plt.errorbar(RUN2_OV3[1],RUN2_OV3[4], xerr=47.52087456, fmt="o")
plt.legend();plt.show()

# plt.title("SC Linearity measured with Laser (405nm)")
# plt.ylabel("Measured PE");plt.xlabel("Theoretical PE")
# plt.plot(x, y, label = "Ref. (Measured PE = Theoretical PE)",alpha = 0.5, color = "k")
# plt.plot(OV3_x,lin_func(OV3_x,popt_OV3),c='tab:blue')
# plt.plot(OV3_x_new,lin_func(OV3_x_new,popt_OV3_new),c='b')
# plt.scatter(OV3_x, OV3, label = "Measured PE (OV = 3.0V) %.3f"%popt_OV3,c='tab:blue')
# plt.scatter(OV3_x_new, OV3, label = "New Measured PE (OV = 3.0V) %.3f"%popt_OV3_new,c='b')
# plt.errorbar(OV3_x, OV3, xerr=DOV3, fmt="o",c='tab:blue')
# plt.errorbar(OV3_x_new, OV3, xerr=DOV3_new, fmt="o",c='b')
# plt.legend();plt.show()

# plt.title("SC Linearity measured with Laser (405nm)")
# plt.ylabel("Measured PE");plt.xlabel("Theoretical PE")
# plt.plot(x, y, label = "Ref. (Measured PE = Theoretical PE)",alpha = 0.5, color = "k")
# plt.plot(OV2_5_x,lin_func(OV2_5_x,popt_OV2_5),c='tab:orange')
# plt.plot(OV2_5_x_new,lin_func(OV2_5_x_new,popt_OV2_5_new),c='orange')
# plt.scatter(OV2_5_x, OV2_5, label = "Measured PE (OV = 2.5V) %.3f"%popt_OV2_5,c='tab:orange')
# plt.scatter(OV2_5_x_new, OV2_5, label = "New Measured PE (OV = 2.5V) %.3f"%popt_OV2_5_new,c='orange')
# plt.errorbar(OV2_5_x, OV2_5, xerr=DOV2_5, fmt="o",c='tab:orange')
# plt.errorbar(OV2_5_x_new, OV2_5, xerr=DOV2_5_new, fmt="o",c='orange')
# plt.legend();plt.show()

# plt.title("SC Linearity measured with Laser (405nm)")
# plt.ylabel("Measured PE");plt.xlabel("Theoretical PE")
# plt.plot(x, y, label = "Ref. (Measured PE = Theoretical PE)",alpha = 0.5, color = "k")
# plt.plot(OV2_x,lin_func(OV2_x,popt_OV2),c='tab:green')
# plt.plot(OV2_x_new,lin_func(OV2_x_new,popt_OV2_new),c='g')
# plt.scatter(OV2_x, OV2, label = "Measured PE (OV = 2.0V) %.3f"%popt_OV2,c='tab:green')
# plt.scatter(OV2_x_new, OV2, label = "New Measured PE (OV = 2.0V) %.3f"%popt_OV2_new,c='g')
# plt.errorbar(OV2_x, OV2, xerr=DOV2, fmt="o",c='tab:green')
# plt.errorbar(OV2_x_new, OV2, xerr=DOV2_new, fmt="o",c='g')
# plt.legend();plt.show()