from ir02_lib import deconvolution_noise, my_wvf

path = '../DATA/JAN_22/AnalysisROOT/'

raw_signal = 'run79_ScintProfFirstSignalBin_SC.root'
pure_signal = 'run56_ScintProfFirstSignalBin_SC.root'
noise = 'run81_ScintProfFirstSignalBin_SC_Noise.root'

deconvolution_noise(path, raw_signal, pure_signal, noise, 2490, 450, 1000, 0.02)
