import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit

def lin_func(x,m):
    return m*x

x = [0,1750]; y = x

OV3_x = [948.53,794.96,622.24,529.19,413.56,337.96,121.70,73.51,40.08,17.07,2.28,0.31]
OV3_x_new = [1179.050737,1000.126363,791.5764211,656.5634605,553.283815,417.3340776,145.3175057,89.5324455,48.09178049,16.78001496,2.6158322,0.2757988084]
OV3 = [906.55,829.76,638.17,530.35,461.03,344.93,125.23,77.20,40.60,14.99,2.37,0.30]
DOV3 = [37.61,31.52,24.67,20.98,16.40,13.40,4.83,2.92,1.59,0.68,0.09,0.01]
DOV3_new = [112.2699529,95.23266141,75.37440477,62.51838573,52.68403292,39.7388134,13.83722431,8.525335782,4.579329591,1.59780358,0.2490811876,0.02626173603]
RUN2_OV3 = [1358.378242,1542.23]
popt_OV3, pcov_OV3 = curve_fit(lin_func, OV3_x, OV3); perr_OV3 = np.sqrt(np.diag(pcov_OV3))
print("OV 3.0 Fitted slope m = %.3f with ∆m = %.3f"%(popt_OV3,perr_OV3))
popt_OV3_new, pcov_OV3_new = curve_fit(lin_func, OV3_x_new, OV3); perr_OV3_new = np.sqrt(np.diag(pcov_OV3_new))
print("OV 3.0 NEW Fitted slope m = %.3f with ∆m = %.3f"%(popt_OV3_new,perr_OV3_new))
print("OV 3.0 Slope differnece ∆m = %.3f"%(popt_OV3-popt_OV3_new))

OV2_5_x = [840.9603542,704.8088729,551.6699924,469.1717816,366.6587144,299.6341742,107.8952134,65.17442745,35.53057497,15.13728638,2.01830485,0.2788624535]
OV2_5_x_new = [1045.210673,886.5969172,701.7205433,582.0335926,490.4777466,369.9603574,128.8217743,79.36916086,42.63263715,14.87522986,2.318895742,0.2444914786]
OV2_5 = [855.47,778.47,610.02,499.58,420.44,316.30,113.20,70.55,36.95,13.47,2.10,0.27]
DOV2_5 = [34.09410105,28.57426609,22.36573029,19.02109897,14.86502806,12.1477282,4.374273164,2.642292831,1.440475769,0.6136938189,0.08182584252,0.01130560391]
DOV2_5_new = [98.60536651,83.6417157,66.2004447,54.90915583,46.27176053,34.90212794,12.1530698,7.487701186,4.021970804,1.403331912,0.2187650494,0.02306537092]
popt_OV2_5, pcov_OV2_5 = curve_fit(lin_func, OV2_5_x, OV2_5); perr_OV2_5 = np.sqrt(np.diag(pcov_OV2_5))
print("OV 2.5 Fitted slope m = %.3f with ∆m = %.3f"%(popt_OV2_5,perr_OV2_5))
popt_OV2_5_new, pcov_OV2_5_new = curve_fit(lin_func, OV2_5_x_new, OV2_5); perr_OV2_5_new = np.sqrt(np.diag(pcov_OV2_5_new))
print("OV 2.5 NEW Fitted slope m = %.3f with ∆m = %.3f"%(popt_OV2_5_new,perr_OV2_5_new))
print("OV 2.5 Slope differnece ∆m = %.3f"%(popt_OV2_5-popt_OV2_5_new))

OV2_x = [690.9497196,579.08496,453.263016,385.4808485,301.2540777,246.1853851,88.64884902,53.54860327,29.19262565,12.43709495,1.658279327,0.229118927]
OV2_x_new = [862.6737762,731.7605246,579.1711915,480.3865193,404.8201005,305.3500187,106.3241788,65.50803151,35.18722017,12.27740114,1.913920894,0.2017931815]
OV2 = [669.90,603.32,475.15,390.13,336.34,250.98,91.79,57.23,30.38,11.07,1.67,0.23]
DOV2 = [4.739460842,3.972142132,3.109086312,2.644145204,2.066404927,1.688669898,0.608072826,0.3673082152,0.2002422206,0.08531029515,0.01137470602,0.001571605215]
DOV2_new = [110.099961,93.39197206,73.91754257,61.31000905,51.66573797,38.97072807,13.56977372,8.360555185,4.490818748,1.566920688,0.2442668615,0.02575414024]
popt_OV2, pcov_OV2 = curve_fit(lin_func, OV2_x, OV2); perr_OV2 = np.sqrt(np.diag(pcov_OV2))
print("OV 2.0 Fitted slope m = %.3f with ∆m = %.3f"%(popt_OV2,perr_OV2))
popt_OV2_new, pcov_OV2_new = curve_fit(lin_func, OV2_x_new, OV2); perr_OV2_new = np.sqrt(np.diag(pcov_OV2_new))
print("OV 2.0 NEW Fitted slope m = %.3f with ∆m = %.3f"%(popt_OV2_new,perr_OV2_new))
print("OV 2.0 Slope differnece ∆m = %.3f"%(popt_OV2-popt_OV2_new))

plt.title("SC Linearity measured with Laser (405nm)")
plt.ylabel("Measured PE");plt.xlabel("Theoretical PE")
plt.plot(x, y, label = "Ref. (Measured PE = Theoretical PE)",alpha = 0.5, color = "k")
plt.scatter(OV3_x, OV3, label = "Measured PE (OV = 3.0V) Fitted slope %.3f"%popt_OV3)
plt.errorbar(OV3_x, OV3, xerr=DOV3, fmt="o")

plt.scatter(OV2_5_x, OV2_5, label = "Measured PE (OV = 2.5V) Fitted slope %.3f"%popt_OV2_5)
plt.errorbar(OV2_5_x, OV2_5, xerr=DOV2_5, fmt="o")

plt.scatter(OV2_x, OV2, label = "Measured PE (OV = 2.0V) Fitted slope %.3f"%popt_OV2)
plt.errorbar(OV2_x, OV2, xerr=DOV2, fmt="o")

plt.scatter(RUN2_OV3[0],RUN2_OV3[1],label = "Recovered Data"); plt.errorbar(RUN2_OV3[0],RUN2_OV3[1], xerr=47.52087456, fmt="o")
plt.legend();plt.show()

plt.title("SC Linearity measured with Laser (405nm)")
plt.ylabel("Measured PE");plt.xlabel("Theoretical PE")
plt.plot(x, y, label = "Ref. (Measured PE = Theoretical PE)",alpha = 0.5, color = "k")
plt.scatter(OV3_x, OV3, label = "Measured PE (OV = 3.0V) %.3f"%popt_OV3)
plt.scatter(OV3_x_new, OV3, label = "New Measured PE (OV = 3.0V) %.3f"%popt_OV3_new)
plt.errorbar(OV3_x, OV3, xerr=DOV3, fmt="o")
plt.errorbar(OV3_x_new, OV3, xerr=DOV3_new, fmt="o")
plt.legend();plt.show()

plt.title("SC Linearity measured with Laser (405nm)")
plt.ylabel("Measured PE");plt.xlabel("Theoretical PE")
plt.plot(x, y, label = "Ref. (Measured PE = Theoretical PE)",alpha = 0.5, color = "k")
plt.scatter(OV2_5_x, OV2_5, label = "Measured PE (OV = 2.5V) %.3f"%popt_OV2_5)
plt.scatter(OV2_5_x_new, OV2_5, label = "New Measured PE (OV = 2.5V) %.3f"%popt_OV2_5_new)
plt.errorbar(OV2_5_x, OV2_5, xerr=DOV2_5, fmt="o")
plt.errorbar(OV2_5_x_new, OV2_5, xerr=DOV2_5_new, fmt="o")
plt.legend();plt.show()

plt.title("SC Linearity measured with Laser (405nm)")
plt.ylabel("Measured PE");plt.xlabel("Theoretical PE")
plt.plot(x, y, label = "Ref. (Measured PE = Theoretical PE)",alpha = 0.5, color = "k")
plt.scatter(OV2_x, OV2, label = "Measured PE (OV = 2.0V) %.3f"%popt_OV2)
plt.scatter(OV2_x_new, OV2, label = "New Measured PE (OV = 2.0V) %.3f"%popt_OV2_new)
plt.errorbar(OV2_x, OV2, xerr=DOV2, fmt="o")
plt.errorbar(OV2_x_new, OV2, xerr=DOV2_new, fmt="o")
plt.legend();plt.show()