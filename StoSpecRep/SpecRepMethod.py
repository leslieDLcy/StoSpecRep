# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# ==============================================================================

""" implementation for Spectral Representation method by Shinozuka """

import numpy as np
import matplotlib.pyplot as plt
import itertools
import pandas as pd
from scipy import interpolate
from scipy.interpolate import griddata
from .utils import EPSD_show


def SRM_formula(Stw, f_vec, t_vec):
    """ Redo the SRM based on a computed wavelet EPSD 

    Hint
    ----
    Swt    --> obj._pwr_coef.shape
    f_vec  --> obj._freqs
    t_vec  --> obj.t_axis
    """

    N1 = len(f_vec)
    w_vec = 2 * np.pi * f_vec
    phi_n = np.random.uniform(0, 2 * np.pi, N1)
    delta_f = np.abs(np.append(arr=np.diff(f_vec),  values=np.diff(f_vec)[-1]))
    delta_w = 2 * np.pi * delta_f

    # main course
    sum = 0
    for i in range(1, N1):
        sum = sum + np.sqrt(4 * Stw[i] * delta_w[i]) * \
            np.cos(w_vec[i] * t_vec + phi_n[i])

    # compensate for the two-sided into one-sided
    simulation = sum / 4
    return simulation


class SRM:
    """ From a given spectra (a stationary Sww or a nonstationary Swt, generate simulations. 

    Hint:
    ----
    Note that when given an estimated Swt, we need to interpolate the Swt first.
    An object is characterised by an EPSD (named as Swt) and resulting sample realizations;
    """

    def __init__(self, wu=100, N1=1024, fs=200, duration=16):
        self.wu = wu  # cutoff frequency
        self.N1 = N1  # total number N
        self.fs = fs
        self.duration = duration
        self.t_axis_4simu = np.arange(
            0, self.duration, 1 / self.fs)  # dt = 1 / Fs
        self.w_axis_4simu = np.arange(0, self.wu, self.wu/self.N1)

    def __str__(self):
        return '\n'.join([
            f'*** Generation started by {SRM.__name__} ***',
            f'cutoff frequency wu (rad/s): {self.wu:.2f}',
            f'the expected shape of Swt (w, t) axes: ({self.w_axis_4simu.shape}, {self.t_axis_4simu.shape})',
            f'the lower limit of sampling frequency (Hz): {np.ceil(1 / (2 * np.pi / (2 * self.wu)))}',
        ])

    def SpecRepsentation0(self, Sww):
        ''' The SRM fomula which generates sample simulation from given Swt/Sww.
        
        Note
        ----
        Each call will result in a different realization.
        For now, this func received a spectra as an argument,
        which may be obtained from 'getSww_from_a_model' func.
        This func is also used in simulations for nonstationary Swt.
        '''

        # create the t axis
        n = np.arange(self.N1)
        delta_w = self.wu / self.N1
        w_n = n * delta_w
        A_n = np.sqrt(2 * Sww * delta_w)
        phi_n = np.random.uniform(0, 2 * np.pi, self.N1)

        # Note that A0=0 or S(w0)=0
        sum = 0
        for i in range(1, self.N1):
            sum = sum + A_n[i] * np.cos(w_n[i] * self.t_axis_4simu + phi_n[i])
        simulation = sum * np.sqrt(2)
        simulation = simulation / (2 * 2)
        return simulation

    def _SpecRepsentation0(self, Sww, plot=True):
        ''' backup version of function `self.SpecRepsentation0`
        
        Note
        ----
        For now, this func received a spectra as argument,
        which may be obtained from 'getSww_from_a_model' func.
        '''

        # create the t axis
        n = np.arange(self.N1)
        delta_w = self.wu / self.N1
        w_n = n * delta_w
        A_n = np.sqrt(2 * Sww * delta_w)
        phi_n = np.random.uniform(0, 2 * np.pi, self.N1)

        # Note that A0=0 or S(w0)=0
        sum = 0
        for i in range(1, self.N1):
            sum = sum + A_n[i] * np.cos(w_n[i] * self.t_axis_4simu + phi_n[i])
        simulation = sum * np.sqrt(2)

        if plot:
            plt.plot(self.t_axis_4simu, simulation)
            plt.title('A simulation by SRM')
            plt.xlabel('time [s]')
            plt.ylabel('Amp')
        return simulation

    @staticmethod
    def getSww_from_a_model(model, w_axis):
        return model(w_axis)

    def _interpo_spectra(self, Pxx, freqs, t_bins):
        """ Given Swt estimated by STFT (one-sided), ie. Pxx, freqs, t_bins, im = ax2.specgram(...)

        Note
        -----
        Follow the example procedures: Points, values, and then new coordinates
        """
        points = list(itertools.product(t_bins, freqs))
        sample_df = pd.DataFrame()
        sample_df['X'] = [xy[0] for xy in points]
        sample_df['Y'] = [xy[1] for xy in points]

        up_reversed = np.flipud(Pxx)
        values = np.ravel(up_reversed, order='F')

        x_min, x_max = 0.0, sample_df['X'].max()
        y_min, y_max = 0.0, sample_df['Y'].max()

        new_x_coord = np.linspace(x_min, x_max, self.duration * self.fs)
        new_y_coord = np.linspace(y_min, y_max, self.N1)
        xx, yy = np.meshgrid(new_x_coord, new_y_coord)

        grid_z0 = griddata(points, values, xi=(xx, yy), method='nearest')
        grid_z1 = np.flipud(grid_z0)

        self._interpolated_bundle = (grid_z1, new_y_coord, new_x_coord)
        return grid_z1

    def get_interpolated_spectra(self, format):
        EPSD_show(*self._interpolated_bundle, format=format,
                  title_name='interpolated spectra')

    def nonsta_simulation(self, Pxx, freqs, t_bins):
        """ For given estimated spectra, do interpolation first and then do simulation;

        Paramters
        ---------
        Pxx, freqs, t_bins: 
        Given the estimated spectra by STFT

        **kwargs
            Used to control if showing the estimated spectra in 2d or 3d
        """
        interpolated_Swt = self._interpo_spectra(Pxx, freqs, t_bins)
        x = self.SpecRepsentation0(interpolated_Swt)
        return x

    ''' !!! the high-level func called in the main file '''

    def nonsta_esmb_simus(self, Pxx, freqs, t_bins, ensemble_num):
        """ high-level func that generates an ensemble of simulations.

        Steps
        -----
        Interpolate an estimated Swt to the expected shape of a PSDF;
        Then simulate realizations from the interpolated Swt;

        Hint
        ----
        Similar to the func above, with only difference that 
        we generate multiple simulations here
        """

        interpolated_Swt = self._interpo_spectra(Pxx, freqs, t_bins)

        simus = []
        for i in range(ensemble_num):
            x = self.SpecRepsentation0(interpolated_Swt)
            simus.append(x)
        simus = np.vstack(simus)
        return simus

    def SpecRepsentation3(self, Sww, t_bins, plot='y'):
        ''' For now, this func received a spectra as argument,
        which may be obtained from 'getSww_from_a_model' func
        '''

        # create the t axis
        n = np.arange(self.N1)
        delta_w = self.wu / self.N1
        w_n = n * delta_w
        A_n = np.sqrt(2 * Sww * delta_w)
        phi_n = np.random.uniform(0, 2 * np.pi, self.N1)

        # Note that A0=0 or S(w0)=0
        sum = 0
        for i in range(1, len(Sww)):
            f = interpolate.interp1d(t_bins,
                                     A_n[i],
                                     kind='next',
                                     bounds_error=False,
                                     fill_value='extrapolate')

            A_new = f(self.t_axis_4simu)
            sum = sum + A_new * np.cos(w_n[i] * self.t_axis_4simu + phi_n[i])

        simulation = sum * np.sqrt(2)
        # print('sampling frequency:', Fs)
        t_upper_limit = 2 * np.pi / (2 * self.wu)
        print("the lower limit of sampling frequency:",
              np.ceil(1 / t_upper_limit))
        print("the length of the simulation", simulation.shape)
        if plot == 'y':
            plt.plot(self.t_axis_4simu, simulation)
        return simulation

    def SpecRepsentation2(self, Sww, t_bins, freqs, plot='y'):
        """
        Idea
        ----

        From an estimated evolutionary spectra, do simulation.
        a estimated EPSD with shape, e.g. (65, 41);

        Based on 'wu', 'N1' -> the freq axis
        Based on 'Fs', 'duration' -> the time axis
        ergo, we will interpolate a spectra with shape (N1, duration * Fs)
        """

        # create the t axis
        n = np.arange(self.N1)
        delta_w = self.wu / self.N1
        w_n = n * delta_w
        A_n = np.sqrt(2 * Sww * delta_w)
        phi_n = np.random.uniform(0, 2 * np.pi, self.N1)

        points = np.meshgrid(t_bins, freqs)
        values = Sww
        grid_x, grid_y = np.mgrid[self.t_axis_4simu, self.w_axis_4simu]
        grid_z0 = griddata(points, values, (grid_x, grid_y), method='nearest')

        return self.SpecRepsentation(grid_z0)
