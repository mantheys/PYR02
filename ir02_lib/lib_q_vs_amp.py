import numpy as np
import scipy as sc
import matplotlib.pyplot as plt

from matplotlib.colors import LogNorm
from ir02_lib.lib_tellme import tellme
from scipy.optimize import curve_fit

def lin_func(x,m,n):
    return m*x+n

def lin_fit(dfin,run,ch,charge,n):

    y_scatt = dfin[charge].to_numpy()
    x_scatt = dfin["Amp"].to_numpy()

    linear_model = np.polyfit(x_scatt,y_scatt,n)
    linear_model_fn = np.poly1d(linear_model)
    print("\nMean relative deviation for fit: ")
    print(np.std((linear_model_fn(x_scatt)-y_scatt)/linear_model_fn(x_scatt)))
    
    return linear_model_fn

def q_vs_amp(dfin,run,ch,month,charge):
    
    np_q = dfin[charge].to_numpy()
    np_amp = dfin["Amp"].to_numpy()

    plt.ioff
    test = False
    plthist = plt.hist2d(np_amp,np_q,200,[[np.min(np_amp),np.max(np_amp)],[np.min(np_q),np.max(np_q)]], norm = LogNorm(), label = "Linear fit run %i and ch %i"%(run,ch))
    plt.title(month+"_RUN%i_CH%i"%(run,ch))
    plt.colorbar(plthist[3])
    plt.xlabel("Amp (ADC counts)", fontsize=12)
    plt.ylabel("Total Charge (%s) in pC"%charge, fontsize=12)

    with plt.ion():
        tellme("Press any key to show fit")
        plt.waitforbuttonpress(-1); test = False
        #n = 0
        while test == False:
           
            popt, pcov = curve_fit(lin_func,np_amp,np_q)
            perr = np.sqrt(np.diag(pcov))
            print("\nFit parameters:")
            print("Line slope = %.2f +- %.4f (%.4f%%)"%(popt[0],perr[0],100*perr[0]/popt[0]))
            print("y-Axis intersection = %.2f +- %.4f (%.4f%%)\n"%(popt[1],perr[1],100*perr[1]/popt[1]))

            """
            n = n+1
            fit = lin_fit(dfin,run,ch,charge,n)
            x = np.linspace(np.min(np_amp), np.max(np_amp), 1000)
            plt.plot(x, fit(x), label = fit)
            plt.legend()
            """
            tellme("Click to exit")
            plt.plot(np_amp, lin_func(np_amp, *popt))
            test = plt.waitforbuttonpress(-1)

    plt.clf()
    #plt.close("all") 