from collections import namedtuple
import seaborn as sns
import pywt
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import ticker
from functools import partial
from .console import console
from .SpecRepMethod import SRM_formula



class CWTx():
    """ We assume we are using 'morl' wavelet """


    def __init__(self, signal, fs, t_axis):
        self.fs = fs
        self.dt = 1 / fs
        self.signal = signal
        self.t_axis = t_axis

    @property
    def goto_scales(self,):
        _goto_scales = np.logspace(start=2, stop=10, num=30,  base=2, endpoint=True)
        console.print(f'(a={2},b={10}) ==> scales_range({_goto_scales[0]}, {_goto_scales[-1]}) ==> {self.freqhelper(_goto_scales)}')
        return _goto_scales


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
        console.print(f'(a={a},b={b}) ==> scales_range({checked_scales[0]}, {checked_scales[-1]}) ==> {self.freqhelper(checked_scales)}')



    def propose_scales(self, a, b, num):
        self._proposed_scales = np.logspace(start=a, stop=b, num=num,  base=2, endpoint=True)
        console.print("You've proposed scales:")
        self.check_scales(a, b, num)


    def computeEPSD(self, externaldata=None):
        ''' Compute EPSD by wavelet using `self._proposed_scales` 
        
        By default, in computes the EPSD from `self.signal`, which is 
        desired in most cases. Given a signal, construct an object and compute the EPSD.
        But from some special uses, we may think the `self.signal` as target and try to compare 
        it with that of some external data, eg. reconstructions.
        
        '''

        if externaldata is not None:
            console.print("Yo! Computing EPSD of external data with the proposed scales")
            
            coef, freqs = pywt.cwt(
                data=externaldata, 
                scales=self._proposed_scales, 
                wavelet='morl', 
                sampling_period=self.dt,
                method='fft',
                axis=1)
            # we store the EPSD across ensemble as an instance variable
            self._EPSD_pwr_coef_all = np.square(np.abs(coef)) * 2 * self.dt
            EPSD_pwr_coef_ensemnle_mean = np.mean(self._EPSD_pwr_coef_all, axis=1)
            return EPSD_pwr_coef_ensemnle_mean, freqs, self.t_axis
        else:
            coef, self._freqs = pywt.cwt(
                        data=self.signal, 
                        scales=self._proposed_scales, 
                        wavelet='morl', 
                        sampling_period=self.dt)
            self._pwr_coef = np.square(np.abs(coef)) * 2 * self.dt
            
            console.print("Yo! Computing EPSD with the proposed scales")
            console.print(f"Swt shape: {self._pwr_coef.shape}")



    def plot_wavelet2dEPSD(self, external_EPSDbundle=None):
        """ Plot the computed EPSD by wavelet transform in 2D """
            
        if external_EPSDbundle is None:
            fig, ax =plt.subplots()
            im = ax.pcolormesh(
                    self.t_axis, self._freqs, self._pwr_coef, 
                    cmap='BuPu', 
                    shading='gouraud',
                    rasterized=True)
            ax.set_xlabel("Time")
            ax.set_ylabel("Frequency (Hz)")
            plt.colorbar(im)
        else:
            fig, ax =plt.subplots()
            im = ax.pcolormesh(
                    external_EPSDbundle[2], external_EPSDbundle[1], external_EPSDbundle[0], 
                    cmap='BuPu', 
                    shading='gouraud',
                    rasterized=True)
            ax.set_xlabel("Time")
            ax.set_ylabel("Frequency (Hz)")
            plt.colorbar(im)






    def plot_3dEPSD(self, external_EPSDbundle=None, *, x_low, x_high, y_low, y_high):
        """ Plot the computed EPSD by wavelet transform of a certain range 
        
        Parameters
        ----------
        x_low: int
            the lower limit of the time axis;
        y_low: int
            the lower limit of the frequency axis;=
        """

        if external_EPSDbundle is None:
            fig = plt.figure(figsize=(8, 8))
            ax = plt.axes(projection='3d')
            X, Y = np.meshgrid(self.t_axis, self._freqs)
            Z = self._pwr_coef
            # a workaround to set up the range
            Z = np.where((X > x_low) & (X < x_high), Z, None)
            Z = np.where((Y > y_low) & (Y < y_high), Z, None)
            ax.plot_surface(X, Y, Z, cmap='coolwarm')
            ax.set_xlabel('Time (s)')
            ax.set_ylabel('Frequency (Hz)')
            ax.set_zlabel('PSD')
            ax.set_xlim3d(left=x_low, right=x_high)
            ax.set_ylim(bottom=y_low, top=y_high)
        else:
            print('You have to feed in an EPSD bundle')
            fig = plt.figure(figsize=(8, 8))
            ax = plt.axes(projection='3d')
            X, Y = np.meshgrid(external_EPSDbundle[2], external_EPSDbundle[1])
            Z = external_EPSDbundle[0]
            # a workaround to set up the range
            Z = np.where((X > x_low) & (X < x_high), Z, None)
            Z = np.where((Y > y_low) & (Y < y_high), Z, None)
            ax.plot_surface(X, Y, Z, cmap='coolwarm')
            ax.set_xlabel('Time (s)')
            ax.set_ylabel('Frequency (Hz)')
            ax.set_zlabel('PSD')
            ax.set_xlim3d(left=x_low, right=x_high)
            ax.set_ylim(bottom=y_low, top=y_high)



    ##### SRM part #####
    def g_a_SRMsimu(self, external_EPSD=None):
        ''' Draw simulations from the estimated EPSD by SRM
        '''

        if external_EPSD is None:
            trial_simulation = SRM_formula(
                Stw=self._pwr_coef, 
                f_vec=self._freqs, 
                t_vec=self.t_axis)
            return trial_simulation
        else:
            trial_simulation = SRM_formula(*external_EPSD)
            return trial_simulation



    def g_ensemble_simus(self, ensemble_size, external_EPSD=None):
        """ draw an ensemble of sample realizations """

        if external_EPSD is None:
            ensemble_list = [self.g_a_SRMsimu() for i in range(ensemble_size)]
            return np.vstack(ensemble_list)
        else:
            single_realization_func = partial(self.g_a_SRMsimu, external_EPSD=external_EPSD)
            ensemble_list = [single_realization_func() for i in range(ensemble_size)]
            return np.vstack(ensemble_list)




    def EPSD_uncertainty_givenranges(self, external_EPSDbundle):
        """ Shown the EPSD uncertain over ensemble size, given predefined ranges

        Parameters
        """

        fig, ax = plt.subplots()

        CriterionArray = np.array([[4.0, 5.0], [5.0, 6.0], [6.0, 7.0], [8.0, 10.0]])
        cordinates_ofinterest = []


        # for each range, give back a pair of coordinates   
        for i in CriterionArray:
            lower_lmt, upper_lmt = i
            cordinates_one = np.argwhere((upper_lmt > external_EPSDbundle[0]) & (external_EPSDbundle[0] > lower_lmt))
            if cordinates_one.size != 0:
                cordinates_ofinterest.append(cordinates_one[-1])
            else:
                continue


        # with an array of coordinates
        cordinates_ofinterest = np.vstack(cordinates_ofinterest)

        # new implementations above break at here
        _a_containter = []
        _indentifier = []

        for item in cordinates_ofinterest:
            x, y = item
            _a_containter.append(self._EPSD_pwr_coef_all[x, :, y])
            _indentifier.append(np.full(shape=self._EPSD_pwr_coef_all[x, :, y].shape, 
                                        fill_value=f"f={external_EPSDbundle[1][x]:.1f},t={external_EPSDbundle[2][y]:.1f}"))

        for_the_data = np.concatenate(_a_containter, axis=None)
        for_the_identifier = np.concatenate(_indentifier, axis=None)
        df = pd.DataFrame({'Swt': for_the_data, 'Location': for_the_identifier})

        ax = sns.kdeplot(
           data=df, x="Swt", hue="Location",
           legend=True,
           fill=True, 
           common_norm=False, 
           palette="crest",
           alpha=.5, 
           linewidth=0,
           cut=0,
        )
        ax.yaxis.set_major_locator(ticker.MaxNLocator(nbins=5))
        # ax.set_xlim([0, 8])
        # ax.legend(fontsize=9, frameon=False, title=None)



        ''' only for gettung the label in the plot '''

        the_labels= []
        for row in cordinates_ofinterest:
            x, y = row
            the_labels.append(r"$f=$" + f"{external_EPSDbundle[1][x]:.1f}" + ', ' + r"$t=$" + f"{external_EPSDbundle[2][y]:.1f}")

        the_labels.reverse()
        ax.legend(labels=the_labels, loc='upper left')


        ax.set_xlabel(r'$S(f, t)$')
        return {"cordinates_ofinterest":cordinates_ofinterest}
