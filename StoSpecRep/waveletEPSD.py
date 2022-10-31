from collections import namedtuple
import pywt
import numpy as np
from rich.console import Console
console = Console()



class CWTx():
    """ We assume we are using 'morl' wavelet """


    def __init__(self, signal, fs):
        self.fs = fs
        self.dt = 1 / fs
        self.signal = signal

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


    def propose_scales(self, a, b, num):
        proposed_scales = np.logspace(start=a, stop=b, num=num,  base=2, endpoint=True)
        console.print(f'scales_range({proposed_scales[0]}, {proposed_scales[-1]})==>{self.freqhelper(proposed_scales)}')
