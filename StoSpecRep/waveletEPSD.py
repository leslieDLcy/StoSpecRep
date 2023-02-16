from collections import namedtuple
import seaborn as sns
import pywt
import numpy as np
import pandas as pd
import itertools
import matplotlib.pyplot as plt
from matplotlib import ticker
from functools import partial
from .console import console
from .SpecRepMethod import SRM_formula
from cycler import cycler
import matplotlib as mpl
import matplotlib.lines as mlines
from itertools import cycle



WaveletBundle = namedtuple('WaveletBundle', ['EPSD_ensemble_mean', 'freqs', 't_axis', 'EPSD_ensemble_all'])

class CWTx():
    """ We assume using 'morl' wavelet """


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
        """ Check freq ranges given scales """

        checked_scales = np.logspace(start=a, stop=b, num=num,  base=2, endpoint=True)
        console.print(f'(a={a},b={b}) ==> scales_range({checked_scales[0]}, {checked_scales[-1]}) ==> {self.freqhelper(checked_scales)}')



    def propose_scales(self, a, b, num):
        """ propose scales for later computation """            

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
            console.print("Yo! Computing EPSD from external data with the proposed scales")
            
            coef, freqs = pywt.cwt(
                data=externaldata, 
                scales=self._proposed_scales, 
                wavelet='morl', 
                sampling_period=self.dt,
                method='fft',
                axis=1)
                
            # we store the EPSD across ensemble as an instance variable
            EPSD_pwr_coef_all = np.square(np.abs(coef)) * 2 * self.dt
            EPSD_pwr_coef_ensemble_mean = np.mean(EPSD_pwr_coef_all, axis=1)
            return WaveletBundle(EPSD_pwr_coef_ensemble_mean, freqs, self.t_axis, EPSD_pwr_coef_all)
        else:
            print('Computing based on the recording')
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
            print('Plotting an external EPSD bundle')
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
            ax.zaxis.set_rotate_label(False)  # disable automatic rotation
            ax.set_zlabel('PSD', rotation=270)
            ax.set_xlim3d(left=x_low, right=x_high)
            ax.set_ylim(bottom=y_low, top=y_high)
        else:
            print('Plotting an exteranl EPSD bundle')
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
            ax.zaxis.set_rotate_label(False)  # disable automatic rotation
            ax.set_zlabel('PSD', rotation=270)
            # ax.set_zlabel('PSD')
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
            trial_simulation = SRM_formula(external_EPSD[0], external_EPSD[1], external_EPSD[2])
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



    # @mpl.rc_context({'axes.prop_cycle': cycler(marker=['o', 'x', '.'])})
    def EPSD_UncertaintyAlongTime(self, external_EPSDbundle, fixedFreqIndex=250, time_range=[6, 7, 8]):
        """ Shown the EPSD uncertain over ensemble size, along the time axis

        Note
        ----
        Since working with evolutionary spectrum, we can fix a `f` and see
        the distribution shape change along the time axis;

        Parameters
        ----------
        external_EPSDbundle : array
            External EPSD arrays from ensemble
        freq : float
            the fixed frequency value
        time_range : array
            the times to show in seconds
        """

        fig, ax = plt.subplots()

        # designate a `f` value;  # f=2.0hz
        fixedfreq = external_EPSDbundle[1][fixedFreqIndex]
        print('the frequency working with', fixedfreq)

        # time in seconds
        console.log(f"select times: {time_range}")
        time_range = np.array(time_range)

        Sft_by_time = {}
        ground_truth= [] 
        the_labels = []

        # simply get Sft values for fixed 'f' for t in every second
        for time in time_range:
            Sft_by_time[f't{time}'] = external_EPSDbundle.EPSD_ensemble_all[fixedFreqIndex, :, int(time * self.fs)]
            ground_truth.append(self._pwr_coef[fixedFreqIndex, int(time * self.fs)])
            the_labels.append(r'$t$' + f' ={time:.1f}')

        selected_Swt = pd.DataFrame.from_dict(Sft_by_time)

        ax = sns.kdeplot(
           data = selected_Swt, 
           bw_adjust=2,
           legend=True,
           linewidth=1.5,
           color='blue',
           ax=ax)

        legend = ax.get_legend()
        handles_sns = legend.legendHandles

        # set the same color for gt vlines
        line_styles = plt.rcParams['axes.prop_cycle'].by_key()['linestyle']
        line_styles.reverse()
        line_styles_cycles = cycle(line_styles) 

        for gt in ground_truth:
            ax.axvline(x=gt, ymin=0, ymax=1, linestyle=next(line_styles_cycles), linewidth=1, color='purple')

        # # set up the legend marker for gt lines
        target_line = mlines.Line2D([], [], color='purple', linestyle=next(line_styles_cycles),
                                    label='target')
        # # append the lebel for the target artist
        the_labels.append('target')

        # ''' only for getting the legend in the plot '''

        # ax.legend(handles=[*handles_sns, target_line], labels=the_labels, loc=0)

        ax.text(3, 0.3, r'$f=$' + f' {fixedfreq:.1f}Hz', fontsize=12)

        # the_labels = ['a', 'b', 'c', 'd']
        ax.legend(labels=the_labels, handlelength=3)
        ax.set_xlabel(r'$S(f, t)$')
        ax.grid(axis='both')
        ax.yaxis.set_major_locator(ticker.MaxNLocator(nbins=5))



    def EPSD_uncertainty_givenranges(self, external_EPSDbundle, style='fillin'):
        """ Shown the EPSD uncertain over ensemble size, given predefined ranges of EPS

        Parameters
        ----------
        external_EPSDbundle : array
            External EPSD arrays from ensemble
        """

        fig, ax = plt.subplots()

        # CriterionArray = np.array([[4.0, 5.0], [5.0, 6.0], [6.0, 7.0], [8.0, 10.0]])
        
        CriterionArray = np.array([[2.0, 4.0], [6.0, 8.0], [8.0, 10.0]])
        cordinates_ofinterest = []

        # for each range, give back a pair of coordinates   
        for i in CriterionArray:
            lower_lmt, upper_lmt = i
            cordinates_one = np.argwhere((upper_lmt > external_EPSDbundle[0]) & (external_EPSDbundle[0] > lower_lmt))
            if cordinates_one.size != 0:
                rnd = np.random.choice(len(cordinates_one), replace=False)
                cordinates_ofinterest.append(cordinates_one[rnd])
            else:
                continue

        # with an array of coordinates
        cordinates_ofinterest = np.vstack(cordinates_ofinterest)

        # new implementations above break at here
        _a_containter = []
        _indentifier = []
        ground_truth = []

        for item in cordinates_ofinterest:
            x, y = item
            _a_containter.append(external_EPSDbundle.EPSD_ensemble_all[x, :, y])
            _indentifier.append(np.full(shape=external_EPSDbundle.EPSD_ensemble_all[x, :, y].shape, 
                                        fill_value=f"f={external_EPSDbundle[1][x]:.1f},t={external_EPSDbundle[2][y]:.1f}"))
            ground_truth.append(self._pwr_coef[x,y])

        for_the_data = np.concatenate(_a_containter, axis=None)
        for_the_identifier = np.concatenate(_indentifier, axis=None)
        df = pd.DataFrame({'Swt': for_the_data, 'Location': for_the_identifier})

        if style == 'fillin':

            # fill-in
            ax = sns.kdeplot(
               data=df, 
               x="Swt", 
               hue="Location",
               bw_adjust=2,
               legend=True,
               fill=True, 
               common_norm=False, 
               palette="crest",
               alpha=.5, 
               linewidth=0,
            )
        elif style == 'nofillin':
            # no fill-in
            ax = sns.kdeplot(
            data=df, 
            x="Swt", 
            hue="Location",
            bw_adjust=2,
            legend=True,
            )

        ax.yaxis.set_major_locator(ticker.MaxNLocator(nbins=5))

        ''' only for getting the label in the plot '''
        the_labels= []
        for row in cordinates_ofinterest:
            x, y = row
            the_labels.append(r"$f=$" + f"{external_EPSDbundle[1][x]:.1f}" + ', ' + r"$t=$" + f"{external_EPSDbundle[2][y]:.1f}")

        the_labels.reverse()
        ax.legend(labels=the_labels, loc=0)
        ax.set_xlabel(r'$S(f, t)$')


        # set up the color of the vertical lines (ground truth)
        palette = itertools.cycle(sns.color_palette("crest"))
        
        colors = []
        for i in range(len(ground_truth)):
            colors.append(next(palette))

        colors.reverse()

        for gt, color in zip(ground_truth, colors):
            ax.axvline(x=gt, ymin=0, ymax=1, color=color, linestyle='--', linewidth=2)

        return {"cordinates_ofinterest":cordinates_ofinterest}


















