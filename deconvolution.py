from ir02_lib.lib import deconvolution_noise, my_wvf

path = '../DATA/FEB_22/AnalysisROOT/'

raw_signal = 'ScintProf_Single_SC_RUN14_CH6.root'
pure_signal = 'ScintProf_Single_SC_RUN03_CH6.root'
noise = 'ScintProf_Single_SC_RUN03_CH6.root'

deconvolution_noise(path, raw_signal, pure_signal, noise, 2490, 450, 1000, 0.02)
