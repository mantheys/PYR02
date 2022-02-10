def q_vs_amp(dfin,run,ch,charge):
    
    np_q = dfin[charge].to_numpy()
    np_amp = dfin["Amp"].to_numpy()

    plt.scatter(np_amp, np_q, s = 25)
    plt.title("Q vs. Amp (RUN %i CH %i)"%(run,ch)) 
    plt.xlabel("Amp (ADC counts)", fontsize=16)
    plt.ylabel(charge, fontsize=16)
    plt.show()