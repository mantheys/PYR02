import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erf
from ir02_lib.deconv_lib import import_scint_prof

def gauss(T,SIGMA):
    return np.exp(-1/2*(np.array(T)/SIGMA)**2)

def func(T,A,SIGMA,TAU,T0):
    return (2*A/TAU)*np.exp((SIGMA/(np.sqrt(2)*TAU))**2-(np.array(T)-T0)/TAU)*(1-erf((SIGMA**2-TAU*(np.array(T)-T0))/(np.sqrt(2)*SIGMA*TAU)))

def funcP(T,P,A,SIGMA,TAU,T0):
    return P + func(T,A,SIGMA,TAU,T0)

def logfunc(T,P,A,SIGMA,TAU,T0):
    return np.log10(P+(2*A/TAU)*np.exp((SIGMA/(np.sqrt(2)*TAU))**2-(np.array(T)-T0)/TAU)*(1-erf((SIGMA**2-TAU*(np.array(T)-T0))/(np.sqrt(2)*SIGMA*TAU))))

def func2(T,P,A1,SIGMA,TAU1,T0,A2,TAU2):
    return P + func(T,A1,SIGMA,TAU1,T0) + func(T,A2,SIGMA,TAU2,T0)

def func2sigma(T,P,A1,SIGMA1,SIGMA2,TAU1,T0,A2,TAU2):
    return P + func(T,A1,SIGMA1,TAU1,T0) + func(T,A2,SIGMA2,TAU2,T0)

def logfunc2(T,P,A1,SIGMA,TAU1,T0,A2,TAU2):
    return np.log(P + func(T,A1,SIGMA,TAU1,T0) + func(T,A2,SIGMA,TAU2,T0))

def logfunc3(T,P,A1,SIGMA,TAU1,T0,A2,TAU2,A3,TAU3):
    return np.log(P + func(T,A1,SIGMA,TAU1,T0) + func(T,A2,SIGMA,TAU2,T0) + func(T,A3,SIGMA,TAU3,T0))

def func3(T,P,A1,SIGMA,TAU1,T0,A2,TAU2,A3,TAU3):
    return P + func(T,A1,SIGMA,TAU1,T0) + func(T,A2,SIGMA,TAU2,T0) + func(T,A3,SIGMA,TAU3,T0)

def logfunc4(T,P,A1,SIGMA,TAU1,T0,A2,TAU2,A3,TAU3,A4,TAU4):
    return np.log(P + func(T,A1,SIGMA,TAU1,T0) + func(T,A2,SIGMA,TAU2,T0) + func(T,A3,SIGMA,TAU3,T0) + func(T,A4,SIGMA,TAU4,T0))

def func4(T,P,A1,SIGMA,TAU1,T0,A2,TAU2,A3,TAU3,A4,TAU4):
    return P + func(T,A1,SIGMA,TAU1,T0) + func(T,A2,SIGMA,TAU2,T0) + func(T,A3,SIGMA,TAU3,T0) + func(T,A4,SIGMA,TAU4,T0)

def logfunc2sigma(T,P,A1,SIGMA1,SIGMA2,TAU1,T0,A2,TAU2):
    return np.log(P + func(T,A1,SIGMA1,TAU1,T0) + func(T,A2,SIGMA2,TAU2,T0))

def wvf_vis(PATH,SHOW,ALIGN,SMOOTH,NORM,LOGY,BS,INIT,X_LIM):
    paths = []; files = []; labels = []

    input = open(PATH+".txt",'r')
    lines = input.readlines()
    for line in lines:
        names = line.split()
        if names[0] == "#" or "":
            continue
        else:
            files.append(names[0].split("/")[-1])
            paths.append(names[0])
            labels.append(names[1])
    input.close()

    # print(files)
    check = False; autozoom = False; fit = False
    term_output = check; timebin = 4e-9
    signals = []
    for i in range(len(paths)):
        signals.append(import_scint_prof(paths[i],timebin,normalize=NORM,trim=False,align=ALIGN,start=INIT,cut_i=0,cut_f=0,invert=False))
        if SHOW == True:
            if SMOOTH != 0:
                signals[i].apply_smooth(SMOOTH)
            pedestal = np.mean(signals[i].wvf[:INIT-10])
            plt.plot(signals[i].wvf_x,np.roll(signals[i].wvf,INIT-np.argmax(signals[i].wvf)),label = labels[i])
            if BS == True:
                plt.axhline(pedestal,c="grey",ls=":")

    if SHOW == True:
        plt.grid()
        plt.xlabel("Time in s"); plt.ylabel("Amp. in a.u.")
        if LOGY == True:
            plt.semilogy()
        plt.xlim(0,X_LIM)
        plt.legend()
        plt.show()
    
    return paths,files,labels,signals
