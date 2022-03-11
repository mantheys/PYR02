import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit

def lin_func(x,m):
    return m*x

x = [0,625]; y = x

OV5_x = [607.8985049,509.4797369,398.7814192,339.1465759,265.0437481,216.5942373,77.99337817,47.11213413,25.68371183,10.94217309,1.458956412,0.2015791442]
OV5   = [638.61,586.69,467.09,386.65,327.63,247.72,90.19,56.30,29.09,10.10,1.57,0.19]
DOV5  = [45.45777309,38.09815963,29.82029915,25.36089161,19.81958907,16.19660455,5.832232288,3.522977415,1.920590913,0.8182399156,0.1090986554,0.01507379756]
popt_OV5, pcov_OV5 = curve_fit(lin_func, OV5_x, OV5); perr_OV5 = np.sqrt(np.diag(pcov_OV5))
print("OV 5.0 Fitted slope m = %.3f with ∆m = %.3f"%(popt_OV5,perr_OV5))

OV4_x = [584.2524833,489.6620062,383.269629,325.9544604,254.7340827,208.1691598,74.95959361,45.27956745,24.68466742,10.5165447,1.40220596,0.1937381235]
OV4   = [535.71,491.84,384.75,322.88,281.76,212.47,78.40,48.18,25.09,9.22,1.44,0.19]
DOV4  = [16.2976133,13.6590297,10.69123432,9.092438458,7.105759397,5.806839618,2.090983786,1.26306503,0.6885741618,0.2933570393,0.03911427191,0.005404288569]
popt_OV4, pcov_OV4 = curve_fit(lin_func, OV4_x, OV4); perr_OV4 = np.sqrt(np.diag(pcov_OV4))
print("OV 4.0 Fitted slope m = %.3f with ∆m = %.3f"%(popt_OV4,perr_OV4))

OV3_x = [495.3919766,415.1880156,324.9771366,276.3791837,215.9909018,176.5081613,63.5587906,38.39287819,20.93031101,8.917055579,1.188940744,0.1642719794]
OV3   = [455.22,412.60,328.13,271.26,234.56,175.60,64.54,40.25,21.12,7.07,1.25,0.16]
DOV3  = [23.75797087,19.91155538,15.58522889,13.25457195,10.3584753,8.464965019,3.048147662,1.841242742,1.003774269,0.4276434756,0.05701913008,0.007878143139]
popt_OV3, pcov_OV3 = curve_fit(lin_func, OV3_x, OV3); perr_OV3 = np.sqrt(np.diag(pcov_OV3))
print("OV 3.0 Fitted slope m = %.3f with ∆m = %.3f"%(popt_OV3,perr_OV3))

plt.title("SiPM Linearity measured with Laser (405nm)")
plt.ylabel("Measured PE");plt.xlabel("Theoretical PE")
plt.plot(x, y, label = "Ref. (Measured PE = Theoretical PE)",alpha = 0.5, color = "k")
plt.scatter(OV5_x, OV5, label = "Measured PE (OV = 5.0V)")
plt.errorbar(OV5_x, OV5, xerr=DOV5, fmt="o")

plt.scatter(OV4_x, OV4, label = "Measured PE (OV = 4.0V)")
plt.errorbar(OV4_x, OV4, xerr=DOV4, fmt="o")

plt.scatter(OV3_x, OV3, label = "Measured PE (OV = 3.0V)")
plt.errorbar(OV3_x, OV3, xerr=DOV3, fmt="o")

plt.legend()
plt.show()