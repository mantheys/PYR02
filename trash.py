#from ROOT import TF2, TH1D, TF1, TFile, TCanvas, gROOT

def plot_hist(df1,colors,ch,var,wall_var,wall,axes=None,nbins=250,xlabel="",ylabel="Counts"):
    #fig = plt.figure(figsize=(9,6), dpi= 100, facecolor='w', edgecolor='k')
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.set_title(var+" "+wall_var)
    plt.hist(df1[(df1["ch"]==ch) ][var],nbins,axes,histtype ='bar', ec='black',color=colors[0])
    plt.hist(df1[(df1["ch"]==ch) & (df1[wall_var]<wall)][var],nbins,axes,histtype ='bar', ec='black',color=colors[1])
    plt.hist(df1[(df1["ch"]==ch) & (df1[wall_var]>wall)][var],nbins,axes,histtype ='bar', ec='black',color=colors[2])
    plt.minorticks_on()
    plt.xlabel(xlabel, fontsize=16)  
    plt.ylabel(ylabel, fontsize=16)

def plot_multi(df1,colors,ch,var,axes=None,nbins=100,leyenda=[],log=False):
    fig = plt.figure(figsize=(9,6), dpi= 100, facecolor='w', edgecolor='k')
    for i in range(len(var)):
        plt.hist(df1[(df1["ch"]==ch) ][var[i]],nbins,axes,histtype ='step',color=colors[i])

    plt.xlabel("Charge (pC)", fontsize=16)  
    plt.ylabel("Counts", fontsize=16)            
    plt.minorticks_on()
    plt.legend(leyenda)
    if log==True:
        plt.yscale('log')

"""
#___ASK FOR PLOT TO BE DISPLAYED___
q1 = input("\nShow plot? (y/n): ")
if q1 == "y" or q1 =="yes" or q1 == "Y":
    plt.show()   
"""

"""
run_list = []
run = 1
while run != 0:
    run = int(input("Input run number (0 to stop): "))
    if run != 0:
        run_list.append(run)
"""