#___IMPORT MODULES___
from tkinter import Y
import numpy as np
import itertools as iter
import matplotlib.pyplot as plt

from ir02_lib.lib_ntuple_to_pd import ntuple_to_pd
from ir02_lib.lib_terminal_input import terminal_input
from ir02_lib.lib_charge_map import charge_map
from ir02_lib.lib_q_vs_amp import q_vs_amp

####################################################################################################################################
# PROGRAMM STARTS WITH INITIALISED VARIABLES AND CONFIGURATION (IN FUTURE VERSIONS MIGHT BE INCLUDED WITH AN EXTERNAL CONFIG FILE) #
####################################################################################################################################

#___INITIALISE VARIABLES TO BE USED IN THE SCRIPT___
fit = []
rates = []

#___PATH TO NTUPLE___
path = "/Users/sergio/DATA/DIC_21/AnalysisROOT/"
#path = "/pc/choozdsk01/palomare/SiPM/SC_Fuente_Alpha_Dic/AnalysisROOT/"

#___PARTIAL AND TOTAL CHARGE RANGES FOR EACH CHANNEL___
#chargeconfig = [          [] #SiPM 1                    [] #SiPM 2                     [] PMT                  [] SC            ]
#chargeconfig = [["QPeakRange2", "QFixRange4"],["QPeakRange2", "QFixRange4"],["QPeakRange1", "QFixRange4"],["QPeakRange1", "QPeak"]]
chargeconfig = [["QPeakRange1", "QFixRange2"],["QPeakRange1", "QFixRange2"],["QPeakRange1", "QPeak"],["QPeakRange1", "QFixRange3"]]

####################################################################################################################################
# PROGRAMM WELCOMES THE USER AND ASKS FOR RUN AND CH NUMBERS TO CREATE THE DATA FRAME FROM NTUPLES                                 #
####################################################################################################################################

#___HELLO WORLD!___
print("\nWELCOME! THIS SCRIPT CAN BE USED TO ITERACTIVELY SELECT EVENTS IN A CHARGE MAP\n")

#___TERMINAL INPUT OF RUN AND CHANNEL ARRAYS___
run_list, ch_list, month = terminal_input()

#___NTUPLE TO PANDAS DATA FRAME TO BE USED BY PYTHON___
raw_data = ntuple_to_pd(path,run_list,ch_list)

####################################################################################################################################
# MAIN ANALYSIS AND COMPUTATIONS                                                                                                   #
####################################################################################################################################

automatic = input("Use automatic selection? (yes / no): ")

#___EXECUTE THE FUNCTIONS THAT YOU WANT ITERATING OVER RUN NUMBER AND CHANNEL___
for run, ch in iter.product(range(0,np.size(run_list)),range(0,np.size(ch_list))):
    
    charge = chargeconfig[ch_list[ch]]
    zoom = False

    rate, nevents, dfout = charge_map(raw_data[run][ch],run_list[run],ch_list[ch],month,charge,automatic,zoom)
    #charge_map(list[run][ch],run_list[run],ch_list[ch],charge,automatic,zoom)
    q_vs_amp(dfout,run_list[run],ch_list[ch],month,charge[0])

    print("\n___RUN %i CH %i___"%(run_list[run],ch_list[ch]))
    print("\nSelected %i of %i events"%(dfout["evt"].size,raw_data[run][ch]["evt"].size))
    print("Evaluated rate = %f Hz\n"%rate)
    rates.append(rate)
    
####################################################################################################################################
# COMBINED RESULTS AND GRPAHICAL OUTPUTS                                                                                           #
####################################################################################################################################

print("Mean rate = %.2f\n"%np.mean(rates))
print("STD rate = %.2f\n"%np.std(rates))
