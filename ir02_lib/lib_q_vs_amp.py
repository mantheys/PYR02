from cProfile import label
import numpy as np
import matplotlib.pyplot as plt

from matplotlib.colors import LogNorm

def lin_fit(dfin,run,ch,charge):

    y_scatt = dfin[charge].to_numpy()
    x_scatt = dfin["Amp"].to_numpy()

    #plt.scatter(x_scatt,y_scatt)

    linear_model = np.polyfit(x_scatt,y_scatt,1)
    linear_model_fn = np.poly1d(linear_model)

    x = np.linspace(np.min(x_scatt),np.max(x_scatt))
    plt.plot(x,linear_model_fn(x),label = "Linear fit run %i and ch %i"%(run,ch))
    
    return linear_model_fn

def q_vs_amp(dfin,run,ch,charge):
    
    np_q = dfin[charge].to_numpy()
    np_amp = dfin["Amp"].to_numpy()

    plt.ioff
    test = False
    plthist = plt.hist2d(np_amp,np_q,200, norm = LogNorm());plt.colorbar(plthist[3])
    plt.xlabel("Amp (ADC counts)", fontsize=12)
    plt.ylabel("Total Charge (%s) in pC"%charge, fontsize=12)
    while test == False:
        with plt.ion():
            fit = lin_fit(dfin,run,ch,charge)
            x = np.linspace(np.min(np_amp), np.max(np_amp), 1000)
            plt.plot(x, fit(x), label = fit)
            plt.legend()
            test = plt.waitforbuttonpress(-1)

    plt.clf()
    plt.close("all") 