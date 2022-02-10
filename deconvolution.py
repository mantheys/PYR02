import ROOT
import numpy as np
import matplotlib.pyplot as plt

from scipy.fftpack import fft, ifft
from ir02_lib import deconvolution_noise, my_wvf

path = 'DATA/JAN_22/AnalysisROOT/'

raw_signal = 'run79_ScintProfFirstSignalBin_SC.root'
pure_signal = 'run56_ScintProfFirstSignalBin_SC.root'
noise = 'run81_ScintProfFirstSignalBin_SC_Noise.root'

deconvolution_noise(path, raw_signal, pure_signal, noise, 2490, 450, 1000, 0.02)

"""
inSiPMName = path+raw_signal
inSiPM = ROOT.TFile.Open(inSiPMName ,"READ")
listkeys_SiPM = inSiPM.GetListOfKeys()
#Load laser name
inLaserName = path+pure_signal
inLaser = ROOT.TFile.Open(inLaserName ,"READ")
listkeys_Laser = inLaser.GetListOfKeys()
# Load Noise name
inNoiseName = path+noise
inNoise = ROOT.TFile.Open(inNoiseName ,"READ")
listkeys_Noise = inNoise.GetListOfKeys()

sipm = my_wvf("vector",inSiPMName,listkeys_SiPM[0].GetName())
laser = my_wvf("vector",inLaserName,listkeys_Laser[0].GetName())
noise = my_wvf("vector",inNoiseName,listkeys_Noise[0].GetName())

smooth=0.6;
laser.apply_smooth(smooth)
sipm.apply_smooth(smooth)

# Plot signals and apply smooth
# fig=plt.figure(figsize=(6,4), dpi= 150, facecolor='w', edgecolor='k')

shift_sipm = np.roll(np.array(sipm.wvf),np.argmax(laser.wvf)-np.argmax(sipm.wvf))
x_time = np.linspace(0,20,5002)

plt.plot(laser.wvf)
plt.plot(shift_sipm)
plt.legend(["Laser Raw", "Signal Raw"]); plt.xlabel("Bin number (1bin = 4ns)")
plt.show()

plt.plot(laser.wvf/max(laser.wvf))
plt.plot(shift_sipm/max(shift_sipm))
plt.legend(["Laser Raw", "Signal Raw"]); plt.xlabel("Bin number (1bin = 4ns)")
plt.show()

sipm_fft = fft(shift_sipm)
laser_fft = fft(np.array(laser.wvf))
noise_fft = fft(np.array(noise.wvf))

sr = 2.5e8
N = len(sipm_fft)
n = np.arange(N)
T = N/sr
freq = n/T


wiener = abs(sipm_fft)**2/(abs(sipm_fft)**2+abs(noise_fft)**2)
wiener_laser = abs(laser_fft)**2/(abs(laser_fft)**2+abs(noise_fft)**2)

#Calculate deconvolved
#deconvolved = sipm.wvf_F/laser.wvf_F
#deconvolved = sipm.wvf_F/laser.wvf_F*wiener_laser
deconvolved = sipm_fft/laser_fft*wiener*wiener_laser
decon = ifft(deconvolved)

# Plot freqs
plt.plot(freq,abs(deconvolved))
plt.plot(freq,abs(sipm_fft))
plt.plot(freq,abs(laser_fft))
plt.grid(which="both"); plt.semilogy(); plt.semilogx();
plt.legend(["Deconvolved", "SiPM", "Laser"])
plt.xlabel("Freq (Hz)")
plt.show()

plt.plot(x_time,decon)
plt.show()

smooth=0.6;
laser.apply_smooth(smooth)
sipm.apply_smooth(smooth)
plt.plot(laser.wvf/max(laser.wvf))
plt.plot(sipm.wvf/max(sipm.wvf))
#plt.semilogy()
plt.legend(["Laser Smooth", "Signal Smooth"]); plt.xlabel("Bin number (1bin = 4ns)")
plt.show()

# Plot freqs
plt.plot(sipm.wvf_F_x,abs(deconvolved))
plt.plot(sipm.wvf_F_x,abs(sipm.wvf_F))
plt.plot(sipm.wvf_F_x,abs(laser.wvf_F))
plt.grid(which="both"); plt.semilogy(); plt.semilogx();
plt.legend(["Deconvolved", "SiPM", "Laser"])
plt.xlabel("Freq (Hz)")
plt.show()

"""
#deconvolution_noise(path, raw_signal, pure_signal, noise, 2400, 200, 1000, 0, "Today")