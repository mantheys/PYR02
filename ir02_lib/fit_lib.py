import numpy as np
from scipy.special import erf

def gauss(T,SIGMA):
    return np.exp(-1/2*(np.array(T)/SIGMA)**2)

def func(T,P,A,SIGMA,TAU,T0):
    # P = 5e-6
    return P+(2*A/TAU)*np.exp((SIGMA/(np.sqrt(2)*TAU))**2-(np.array(T)-T0)/TAU)*(1-erf((SIGMA**2-TAU*(np.array(T)-T0))/(np.sqrt(2)*SIGMA*TAU)))

def testfunc(T,P,A,SIGMA,TAU,T0):
    # P = 5e-6
    return P+(2*A/TAU)*np.exp((SIGMA/(np.sqrt(2)*TAU))**2-(np.array(T)-T0)/TAU)*(1-erf((SIGMA**2-TAU*(np.array(T)-T0))/(np.sqrt(2)*SIGMA*TAU)))

def logfunc(T,P,A,SIGMA,TAU,T0):
    # P = 5e-6
    return np.log10(P+(2*A/TAU)*np.exp((SIGMA/(np.sqrt(2)*TAU))**2-(np.array(T)-T0)/TAU)*(1-erf((SIGMA**2-TAU*(np.array(T)-T0))/(np.sqrt(2)*SIGMA*TAU))))

def func2(T,P,A1,SIGMA,TAU1,T0,A2,TAU2):
    # P = P/2
    return func(T,P,A1,SIGMA,TAU1,T0) + func(T,P,A2,SIGMA,TAU2,T0)

def func2sigma(T,P,A1,SIGMA1,SIGMA2,TAU1,T0,A2,TAU2):
    # P = P/2
    return func(T,P,A1,SIGMA1,TAU1,T0) + func(T,P,A2,SIGMA2,TAU2,T0)

def logfunc2(T,P,A1,SIGMA,TAU1,T0,A2,TAU2):
    # P = P/2
    return np.log(func(T,P,A1,SIGMA,TAU1,T0) + func(T,P,A2,SIGMA,TAU2,T0))

def logfunc3(T,P,A1,SIGMA,TAU1,T0,A2,TAU2,A3,TAU3):
    # P = P/2
    return np.log(func(T,P,A1,SIGMA,TAU1,T0) + func(T,P,A2,SIGMA,TAU2,T0) + func(T,P,A3,SIGMA,TAU3,T0))

def func3(T,P,A1,SIGMA,TAU1,T0,A2,TAU2,A3,TAU3):
    # P = P/2
    return func(T,P,A1,SIGMA,TAU1,T0) + func(T,P,A2,SIGMA,TAU2,T0) + func(T,P,A3,SIGMA,TAU3,T0)

def logfunc4(T,P,A1,SIGMA,TAU1,T0,A2,TAU2,A3,TAU3,A4,TAU4):
    # P = P/2
    return np.log(func(T,P,A1,SIGMA,TAU1,T0) + func(T,P,A2,SIGMA,TAU2,T0) + func(T,P,A3,SIGMA,TAU3,T0) + func(T,P,A4,SIGMA,TAU4,T0))

def func4(T,P,A1,SIGMA,TAU1,T0,A2,TAU2,A3,TAU3,A4,TAU4):
    # P = P/2
    return func(T,P,A1,SIGMA,TAU1,T0) + func(T,P,A2,SIGMA,TAU2,T0) + func(T,P,A3,SIGMA,TAU3,T0) + func(T,P,A4,SIGMA,TAU4,T0)

def logfunc2sigma(T,P,A1,SIGMA1,SIGMA2,TAU1,T0,A2,TAU2):
    # P = P/2
    return np.log(func(T,P,A1,SIGMA1,TAU1,T0) + func(T,P,A2,SIGMA2,TAU2,T0))