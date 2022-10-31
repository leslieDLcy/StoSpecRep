from collections import namedtuple
import pywt
import numpy as np
import matplotlib.pyplot as plt
from .console import console



class CWTx():
    """ We assume we are using 'morl' wavelet """


    def __init__(self, signal, fs, t_axis):
        self.fs = fs
        self.dt = 1 / fs
        self.signal = signal
        self.t_axis = t_axis

    goto_scales = np.logspace(start=2, stop=10, num=30,  base=2, endpoint=True)


    def scale2freq(self, scales):
        """ 
        a direct func to get frequency values in Hertz
        Due to the base function has weird units;
        """

        freq_hz = pywt.scale2frequency(wavelet='morl', scale=scales) / self.dt
        return freq_hz


    def freq2scale(self, freqs):
        """ input frequency output scales """

        freqs = freqs * self.dt
        return pywt.frequency2scale('morl', freqs)



    def freqhelper(self, scales):
        ''' when given a proposed scales, show the range of frequency value'''

        freq_range = self.scale2freq(scales)
        FreqRange = namedtuple('FreqRange', 'low high')
        return FreqRange(freq_range[-1], freq_range[0])


    def check_scales(self, a, b, num):
        checked_scales = np.logspace(start=a, stop=b, num=num,  base=2, endpoint=True)
        console.print(f'scales_range({checked_scales[0]}, {checked_scales[-1]}) ==> {self.freqhelper(checked_scales)}')



    def propose_scales(self, a, b, num):
        self._proposed_scales = np.logspace(start=a, stop=b, num=num,  base=2, endpoint=True)



    def computeEPSD(self, chosen_scales):
        coef, self._freqs = pywt.cwt(
                        data=self.signal, 
                        scales=chosen_scales, 
                        wavelet='morl', 
                        sampling_period=self.dt)
        self._pwr_coef = np.square(np.abs(coef)) * 2 * self.dt



    def plot_waveletEPSD(self, option='3d'):
            """ Plot the computed EPSD by wavelet transform """
            
            if option == '2d':
                fig, ax =plt.subplots()
                im = ax.pcolormesh(self.t_axis, self._freqs, self._pwr_coef, cmap='BuPu', shading='gouraud')
                ax.set_ylim([0, 20])
                ax.set_xlabel("time")
                ax.set_ylabel("freq")
                plt.colorbar(im)
            elif option == '3d':
                fig = plt.figure(figsize=(8,8))
                ax = plt.axes(projection='3d')
                X, Y = np.meshgrid(self.t_axis, self._freqs)
                Z = self._pwr_coef
                ax.plot_surface(X, Y, Z, cmap='coolwarm')
                ax.set_xlabel('time (s)')
                ax.set_ylabel('frequency (Hz)')
                ax.set_zlabel('PSD')
