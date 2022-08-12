import ROOT
from cycler import cycler
import numpy as np
import scipy as sc
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit
from ir02_lib.fit_lib import funcP,logfunc3,func3,wvf_vis
from ir02_lib.deconv_lib import signal_int

print("\n### WELCOME TO THE SiPM FIT ###")

input_source = "LIST" # SELECT BETWEEN "LIST" OR "CONFIG_FILE"
input_source_path = "input/"
paths,files,labels,signals = wvf_vis(input_source_path+"WVF_COMPARISON", SHOW = True, ALIGN = True, SMOOTH = 0, NORM = False, LOGY = True, BS = False, INIT = 500, X_LIM = 50e-6)

counter = 0
for signal in signals:

    # path = "/pnfs/ciemat.es/data/neutrinos/Sergio/Eficiencia/GAr_Julio/GAr_26_07_2022/AnalysisROOT/"
    # file = "ScintProfFirstSignalBin_SiPM_run01_ch0"
    # signal = import_scint_prof(path+file+".root",timebin = 4e-9,normalize=1,trim=0,align=False,start=0,cut_i=0,cut_f=0,invert=False)
    
    # path = "input/LAr/FEB2/"
    # file = "FEB_2_SiPM_DAY1_OV1"
    # detector,timebin,int_st,paths,filename,shift,filter_strenght,reverse,fnal,s_start,particle_label = import_deconv_runs(path+file+".txt",debug = True)
    
    # init=0; trm = False
    # signal = import_scint_prof(paths[0]+paths[2],timebin,normalize=True,trim=trm,align=shift,start=s_start,cut_i=init,cut_f=fnal,invert=False)
    timebin = 4e-9; detector="SiPM"
    path = paths[counter]; file = files[counter]; label = labels[counter]
    int_sig_all,f_sig_all,i_sig_all= signal_int(label,signal.wvf,timebin,detector,"ALL",th = 1e-3,out = True)
    signal.wvf = signal.wvf/np.max(signal.wvf)
    # signal.apply_smooth(SMOOTH)

    buffer1     = 40
    check       = True
    autozoom    = False
    logy        = True
    thrld       = 1e-6
    term_output = check
    
    buffer2     = 5000
    if buffer2 > len(signal.wvf):
        buffer2 = len(signal.wvf)-1

    for i in range(len(signal.wvf)):
        if signal.wvf[i] <= thrld:
            signal.wvf[i] = thrld
        if np.isnan(signal.wvf[i]):
            signal.wvf[i] = thrld

    output_file = open("output/"+file+"_SiPM_FIT.txt","w")
    max_in = np.argmax(signal.wvf)
    guess_t0 = signal.wvf_x[np.argmax(signal.wvf)-10]

    for i in range(len(signal.wvf[max_in:])):
        if signal.wvf[max_in+i] <= 6e-3*signal.wvf[max_in+i] and i+max_in<buffer2:
            buffer2 = i
            break

    pedestal = np.mean(signal.wvf[:max_in-buffer1])
    print("Pedestas estimate: %.2E"%pedestal)

    p = pedestal; p_low = pedestal*1e-1; p_high = pedestal*10
    t0 = guess_t0; t0_low = guess_t0*0.02; t0_high = guess_t0*50
    sigma = 2e-8; sigma_low = 2e-9; sigma_high = 2e-8
    a1 = 2e-8; a1_low = 1e-9;  a1_high = 9e-7                    
    a2 = 2e-9; a2_low = 1e-12; a2_high = 9e-8                    
    a3 = 2e-9; a3_low = 1e-12; a3_high = 9e-8
    tau1 = 8e-8; tau1_low = 6e-9; tau1_high = 1e-7
    tau2 = 8e-7; tau2_low = tau1_high; tau2_high = 1e-6
    tau3 = 8e-6; tau3_low = tau2_high; tau3_high = 1e-5

    bounds = ([p_low,a1_low,sigma_low,tau1_low,t0_low],[p_high,a1_high,sigma_high,tau1_high,t0_high])
    initial = (p,a1,sigma,tau1,t0)
    variables = ["PEDT","AMP1","SIGM","TAU1","TIME"]

    # FIT PEAK
    popt, pcov = curve_fit(lambda T,P,A1,SIGMA,TAU1,T0: funcP(T,P,A1,SIGMA,TAU1,T0),signal.wvf_x[max_in-buffer1:buffer2],signal.wvf[max_in-buffer1:buffer2],p0 = initial, bounds = bounds, method = "trf")
    perr1 = np.sqrt(np.diag(pcov))

    # PRINT FIRST FIT VALUE
    print("\n--- FISRT FIT VALUES (FAST) ---")
    for i in range(len(initial)):
        print("%s: \t%.2E \u00B1 %.2E"%(variables[i],popt[i],perr1[i]))
    print("-------------------------------")

    # EXPORT FIT PARAMETERS
    p = popt[0];a1 = popt[1];sigma = popt[2];tau1 = popt[3];t0 = popt[4]
    sigma2 = sigma*10; sigma2_low = sigma; sigma2_high = sigma*100

    # CHECK FIRST FIT
    fig,ax0 = plt.subplots(nrows=1)
    custom_cycler = (cycler(color=['#a6e22e', '#a1efe4', '#2ca02c', '#66d9ef', '#ae81ff', '#fd5ff0', '#f92672', '#fd971f', '#e6db74', '#cfcfc2']))
    # ax0.set_prop_cycle(custom_cycler)

    plt.plot(signal.wvf_x,signal.wvf,label = labels[counter])
    plt.axhline(np.mean(signal.wvf[:max_in-buffer1]),c="grey",ls=":")
    plt.plot(signal.wvf_x[max_in-buffer1:max_in+int(buffer1/2)],funcP(signal.wvf_x[max_in-buffer1:max_in+int(buffer1/2)],*popt))
    plt.axvline(signal.wvf_x[buffer2],ls = "--",c = "k")
    plt.xlabel("Time in s"); plt.ylabel("Amp in a.u.")
    plt.semilogy()
    plt.legend();plt.show()

    # USING VALUES FROM FIRST FIT PERFORM SECONG FIT FOR THE SLOW COMPONENT

    bounds2 = ([pedestal,sigma2_low,a2_low,tau2_low,a3_low,tau3_low],[pedestal*100,sigma2_high,a2_high,tau2_high,a3_high,tau3_high])
    initial2 = (pedestal*10,sigma2,a2,tau2,a3,tau3)

    variables2 = ["PEDT","SIGM","AMP2","TAU2","AMP3","TAU3"]
    popt, pcov = curve_fit(lambda T,P,SIGMA2,A2,TAU2,A3,TAU3: logfunc3(T,P,a1,SIGMA2,tau1,t0,A2,TAU2,A3,TAU3),signal.wvf_x[max_in-buffer1:buffer2],np.log(signal.wvf[max_in-buffer1:buffer2]),p0 = initial2, bounds = bounds2, method = "trf")
    perr2 = np.sqrt(np.diag(pcov))

    # sigma2 = popt[0];a2 = popt[1];tau2 = popt[2];a3 = popt[3];tau3 = popt[4]
    p = popt[0];sigma2 = popt[1];a2 = popt[2];tau2 = popt[3];a3 = popt[4];tau3 = popt[5]
    param = [p,a1,sigma2,tau1,t0,a2,tau2,a3,tau3]

    output_file.write("%.2E \t\u00B1\t %.2E\n"%(p,perr1[0]))
    output_file.write("%.2E \t\u00B1\t %.2E\n"%(t0,perr1[4]))
    output_file.write("%.2E \t\u00B1\t %.2E\n"%(sigma2,perr2[0]))
    output_file.write("%.2E \t\u00B1\t %.2E\n"%(a1,perr1[1]))
    output_file.write("%.2E \t\u00B1\t %.2E\n"%(tau1,perr1[3]))
    output_file.write("%.2E \t\u00B1\t %.2E\n"%(a2,perr2[1]))
    output_file.write("%.2E \t\u00B1\t %.2E\n"%(tau2,0))
    output_file.write("%.2E \t\u00B1\t %.2E\n"%(a3,perr2[2]))
    output_file.write("%.2E \t\u00B1\t %.2E\n"%(tau3,0))

    print("\n--- SECOND FIT VALUES (SLOW) ---")
    for i in range(len(initial2)):
        print("%s: \t%.2E \u00B1 %.2E"%(variables2[i],popt[i],perr2[i]))
    print("--------------------------------\n")

    plt.plot(signal.wvf_x[max_in-buffer1:],func3(signal.wvf_x[max_in-buffer1:],*param),c="tab:orange")
    plt.plot(signal.wvf_x[max_in-buffer1:buffer2],func3(signal.wvf_x[max_in-buffer1:buffer2],*param),c="tab:green",label="FITTED WAVEFORM")
    # plt.plot(np.arange(0,20e-6,4e-9),func3(np.arange(0,20e-6,4e-9),*param),c="tab:green")
    plt.plot(signal.wvf_x,signal.wvf,zorder=0,c="tab:blue",label=labels[counter])
    plt.xlabel("Time in s"); plt.ylabel("Amp in a.u.")
    plt.axhline(np.mean(signal.wvf[:max_in-buffer1]),c="grey",ls=":")
    plt.axvline(signal.wvf_x[buffer2],ls = "--",c = "k")
    plt.semilogy()
    plt.legend();plt.show()
    output_file.close()
    print("Generated output file: %s"%(file+"_SiPM_FIT.txt"))
    counter = counter+1
