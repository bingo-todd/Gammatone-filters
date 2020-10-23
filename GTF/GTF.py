import ctypes
import os
import numpy as np
import matplotlib.pyplot as plt
from numpy.ctypeslib import ndpointer


class GTF:
    """Python interface for all-pole gammatone filters wrote using C
        a python version filters is also included"""

    _package_dir = os.path.dirname(os.path.abspath(__file__))
    _lib_fpath = os.path.join(_package_dir, 'libGTF.so')
    _c_fpath = os.path.join(_package_dir, 'GTF.c')

    def __init__(self, fs, cf_low=None, cf_high=None, freq_low=None,
                 freq_high=None, n_band=1):
        """
        Args:
            fs: sample frequency
            cf_low,cf_high: the lowest and highest center frequency
            freq_low,freq_high: the lowest and highest cutoff frequency
            n_band: frequency bands number
        """
        self.bw_factor = 1.019
        # args check
        if cf_low is None:
            if freq_low is not None:
                # set freq_low as the lowest cutoff frequency
                cf_low = ((2*freq_low+self.bw_factor*24.7)
                          / (2-self.bw_factor*24.7*4.37/1000))
            else:
                raise Exception('neither cf_low or freq_low is specified')
        if cf_high is None:
            if freq_high is not None:
                # make freq_high as the highest cutoff frequency
                cf_high = ((2*freq_high-self.bw_factor*24.7)
                           / (2+self.bw_factor*24.7*4.37/1000))
            else:
                raise Exception('neither cf_high or freq_high is specified')

        # center frequencies
        cfs = self.divide_freq_space(freq_low=cf_low,
                                     freq_high=cf_high,
                                     n_band=n_band)
        # bandwidths
        bws = self.cal_bw(cfs)
        if (cfs is None) or (bws is None):
            raise Exception('cfs and bws uninitlaized')

        self.cf_low = cf_low
        self.cf_high = cf_high
        self.fs = np.int32(fs)
        self.cfs = cfs
        self.bws = bws
        self.n_band = n_band

        # load shared library
        self._load_lib()

    def _load_lib(self):
        """Load c library
        """
        if not os.path.exists(self._lib_fpath):
            if not os.path.exists(self._c_fpath):
                raise OSError(0, 'missing file', 'libGTF.so and GTF.c')
            else:
                print('missing libGTF.so, compile from C code')
                os.system('gcc -fPIC -shared {} -o {}'.format(self._c_fpath,
                                                              self._lib_fpath))
                print('compile success')

        _cmodel = ctypes.cdll.LoadLibrary(self._lib_fpath)
        self._gt_filter = _cmodel.GTF

        # double* GTF(double*x, int x_len, int fs,
        #             double*cfs, double*bws, int n_band,
        #             int is_env_aligned, int is_fine_aligned,
        #             int delay_common, int is_gain_norm)
        self._gt_filter.argtypes = (
            ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'),
            ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'),
            ctypes.c_int, ctypes.c_int,
            ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'),
            ndpointer(ctypes.c_double, flags='C_CONTIGUOUS'),
            ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_int,
            ctypes.c_int)
        self._gt_filter.restype = None
        # ctypes.POINTER(ctypes.c_double)

        self._free_mem = _cmodel.free_mem
        self._free_mem.argtypes = (ctypes.POINTER(ctypes.c_double),)
        self._free_mem.restype = None

    def divide_freq_space(self, freq_low, freq_high, n_band,
                          divide_type='ERB'):
        """Divide frequency range (freq_low~freq_high) equally in erb
        scale(default)
        Args:
            freq_low: low bound of frequency range
            freq_high: high bound of frequency range
            n_band: segments number frequency range to be divided
            divide_type: default to ERB
        """
        if divide_type == 'ERB':
            if n_band == 1:
                return np.asarray(freq_low, dtype=np.float).reshape(1, )
            low_erb = self.Hz2ERBscal(freq_low)
            high_erb = self.Hz2ERBscal(freq_high)
            erb_elem = (high_erb-low_erb)/(n_band-1)
            f = self.ERBscal2Hz(low_erb+erb_elem*np.arange(n_band))
        else:
            raise Exception('unsupport Divide type')
        return f

    def Hz2ERBscal(self, freq):
        """convert Hz to ERB scale"""
        return 21.4*np.log10(4.37*freq/1e3+1)

    def ERBscal2Hz(self, erb_num):
        """convert ERB scale to Hz"""
        return (10**(erb_num/21.4)-1)/4.37*1e3

    def cal_ERB(self, cf):
        """calculate the ERB(Hz) of given center frequency based on equation
        given by Glasberg and Moore
        Args
            cf: center frequency Hz, single value or numpy array
        """
        return 24.7*(4.37*cf/1000+1.0)

    def cal_bw(self, cf):
        """calculate the 3-dB bandwidth
        Args
            cf: center frequency Hz, single value or numpy array
        """
        erb = self.cal_ERB(cf)
        return self.bw_factor*erb

    def filter_py(self, x, align_env=False, align_tfs=False,
                  delay_common=None):
        """Filters in Python
        Args:
        x: signal with shape of [x_len,n_band], if x only has single dimension,
        n_band will be added as 1
        align_env: aligned peaks of Gammatone filter impulse response
        align_tfs: aligned the time fine structure
        delay_common: if aligned, give the same delay to all channels,
        default, aligned to the maximum delay
        Returns:
        fitler result with the shape of [n_band,x_len,n_band]
        """
        # constant variables
        tpt = 2*np.pi*(1.0/self.fs)

        x = x.copy()
        if len(x.shape) > 2:
            raise Exception('two many dimensions for x')
        # ensure x is 2-D array
        x_len = x.shape[0]
        if len(x.shape) == 1:
            x = np.reshape(x, [-1, 1])
        n_chann = x.shape[1]

        if align_env is True:
            delays = np.round(3.0/(2.0*np.pi*self.bws)*self.fs)/self.fs
            if delay_common is not None:
                delay_common = np.max(delays)
                delays = np.int(delays-delay_common)
        else:
            delays = np.zeros(self.n_band)

        # IIR and FIR filters outputs
        out_a = np.zeros((5, n_chann), dtype=np.complex)
        coefs_a = np.zeros(5)
        out_b = 0
        coefs_b = np.zeros(4)

        norm_factors = np.zeros(self.n_band)
        y = np.zeros((self.n_band, x_len, n_chann), dtype=np.float)
        for band_i in range(self.n_band):
            bw = self.bws[band_i]
            cf = self.cfs[band_i]
            k = np.exp(-tpt*bw)
            # filter coefs
            coefs_a = [1, 4*k, -6*k**2, 4*k**3, -k**4]
            coefs_b = [1, 1, 4*k, k**2]
            #
            norm_factors[band_i] = (1-k)**4/(1+4*k+k**2)*2
            delay_len_band = np.int(delays[band_i]*self.fs)
            if align_tfs:
                phi_init = -3.0*cf/self.bws[band_i]
            else:
                phi_init = 0
            for sample_i in range(x_len):
                freq_shiftor = np.exp(-1j*(tpt*cf*sample_i))
                # IIR part
                out_a[0, :] = x[sample_i, :]*freq_shiftor*np.exp(1j*phi_init)
                for order_i in range(1, 5):
                    out_a[0, :] = (out_a[0, :]
                                   + coefs_a[order_i]*out_a[order_i, :])
                # FIR part
                out_b = 0
                for order_i in range(1, 4):
                    out_b = out_b+coefs_b[order_i]*out_a[order_i, :]

                if sample_i > delay_len_band:
                    y[band_i, sample_i-delay_len_band, :] = (
                                    norm_factors[band_i]
                                    * np.real(out_b
                                              * np.conjugate(freq_shiftor)))
                # update IIR output
                for order_i in range(4, 0, -1):
                    out_a[order_i, :] = out_a[order_i-1, :]
        return np.squeeze(y)

    def filter(self, x, is_env_aligned=False, is_fine_aligned=False,
               delay_common=-1, is_gain_norm=True):
        """filter input signal using c library
        Args:
            x: signal with shape of [x_len,n_chann], if x only has single
               dimension, n_chann will be added as 1
            is_env_aligned: whether to align the envelope of ir,
                default to False
            is_fine_aligned: whether to align the fine-structure if envelope is
                             aligned, default to False
            delay_common: set delay(s) of filter ir if it is aligned, default
                          to -1, which aligne all filter to the maximum peak
                          position
            is_gain_norm: whether to normlize gains at center frequency
        Returns:
            fitler result with the shape of [n_band,x_len,n_chann]
        """
        # check inputs
        if not isinstance(x, np.ndarray):
            raise Exception()

        x = x.copy()

        if len(x.shape) > 2:
            raise Exception('two many dimensions for x')

        # ensure x is 2-D array
        x_len = x.shape[0]
        if len(x.shape) == 1:
            x = np.reshape(x, [-1, 1])
        n_chann = x.shape[1]

        # convert bool values to int values(0,1) as the inputs of c function
        x_band_all = []
        for chann_i in range(n_chann):
            # data in slice are not stored in continue memmory
            x_band = np.zeros((self.n_band, x_len), dtype=np.float)
            self._gt_filter(x_band,
                            x[:, chann_i].astype(np.float),
                            x_len,
                            self.fs, self.cfs, self.bws, self.n_band,
                            is_env_aligned, is_fine_aligned,
                            delay_common, is_gain_norm)
            x_band_all.append(x_band[:, :, np.newaxis])
        x_filtered = np.concatenate(x_band_all, axis=2)
        return x_filtered

    def plot_delay_gain_cfs(self):
        """ Plot delay and center-frequency gain of gammatone filter
        before alignment and gain-normalization
        """
        # k = np.exp(-2*np.pi/self.fs*self.bws)
        Q = np.divide(self.cfs, self.bws)  # quality of filter

        temp1 = 8*Q*(1-4*Q**2)
        temp2 = np.multiply(Q**2, Q*(16*Q**2-24))
        phi_delay = (np.arctan(np.divide(temp1, temp2+2))
                     - np.arctan(np.divide(temp1, temp2+1)))
        delays = phi_delay/(2*np.pi)*1e3

        correct_factor = np.sqrt(((temp2+1)**2+temp1**2)
                                 / ((temp2+2)**2+temp1**2))
        gains = 10*np.log10(3/(2*np.pi*self.bws))*correct_factor

        fig = plt.figure(figsize=[8, 3])
        ax = fig.subplots(1, 2)
        # gains at cfs
        ax[0].plot(self.cfs/1000, gains, linewidth=2)
        ax[0].set_xlabel('Center frequency(kHz)')
        ax[0].set_ylabel('Gain(dB)')
        # delay at cfs
        ax[1].plot(self.cfs/1000, delays, linewidth=2)
        ax[1].set_yscale('log')
        ax[1].set_xlabel('Center frequency(kHz)')
        ax[1].set_ylabel('Delay(ms)')
        plt.tight_layout()

        return fig

    def plot_filter_spectrum(self, cf=4e3):
        order = 4
        fs = self.fs
        bw = self.cal_ERB(cf)
        freq_bins = np.arange(1, fs/2)  # frequency resolution: 1Hz
        # n_freq_bin = freq_bins.shape[0]
        gain_func = 6 / (((2*np.pi*bw)**order)
                         * (np.divide(1, 1+1j*(freq_bins-cf)/bw)**order
                            + np.divide(1, 1+1j*(freq_bins+cf)/bw)**order))

        amp_spectrum = np.abs(gain_func)

        phase_spectrum = np.angle(gain_func)
        cf_bin_index = np.int16(cf)
        # unwrap based on phase at cf
        phase_spectrum[:cf_bin_index] = np.flip(
                                         np.unwrap(
                                          np.flip(
                                            phase_spectrum[:cf_bin_index])))
        phase_spectrum[cf_bin_index:] = np.unwrap(
                                            phase_spectrum[cf_bin_index:])
        # delays = np.divide(phase_spectrum,freq_bins)

        linewidth = 2
        # Amplitude-phase spectrum
        fig = plt.figure()
        ax1 = fig.subplots(1, 1)
        color = 'tab:red'
        ax1.semilogy(freq_bins/1000, amp_spectrum, color=color,
                     linewidth=linewidth, label='amp')
        ax1.set_ylabel('dB', color=color)
        ax1.set_xlabel('Frequency(kHz)')
        ax1.legend(loc='upper left')
        ax1.tick_params(axis='y', labelcolor=color)
        ax1.set_title('cf={}Hz'.format(cf))

        ax2 = ax1.twinx()
        color = 'tab: blue'
        ax2.plot(freq_bins/1000, phase_spectrum, color=color,
                 linewidth=linewidth, label='phase')
        ax2.legend(loc='upper right')
        ax2.plot([cf/1000, cf/1000], [-8, 8], '-.', color='black')
        ax2.plot([0, fs/2/1000], [0, 0], '-.', color='black')
        ax2.set_ylabel('phase(rad)', color=color)
        ax2.tick_params(axis='y', labelcolor=color)

        return fig

    def plot_ir_spec(self, ir, fs=None, cfs=None, fig=None):
        """plot the waveform and spectrum of given impulse response
        Args:
            ir: impulse response
            fs: sample frequency,use self.fs as default
            fig: handle of matplotlib figure, if not given, not figure will
                be created
            title: title for ir waveform sub-panel
        """

        if fs is None:
            fs = self.fs
        if cfs is None:
            cfs = self.cfs

        ir_len = ir.shape[1]
        N_fft = ir_len
        N_fft_half = np.int(N_fft/2)

        index = np.flip(np.argsort(cfs))
        cfs = cfs[index]
        ir = ir[index, :]

        spec = np.abs(np.fft.fft(ir, N_fft, axis=1))[:, :N_fft_half]

        time_ticks = np.arange(ir_len)/self.fs
        freq_ticks = np.arange(N_fft_half)/N_fft*self.fs
        x_lim_max = 0.08
        linewidth = 2
        if fig is None:
            fig = plt.figure(figsize=[8, 3])
        ax = fig.subplots(1, 2)
        ax[0].plot(time_ticks, ir.T, linewidth=linewidth)
        ax[0].set_xlabel('Time(s)')
        ax[0].set_title('irs')
        ax[0].set_xlim([0, x_lim_max])

        ax[1].plot(freq_ticks, spec.T, linewidth=linewidth)
        ax[1].set_xlim([self.cf_low/8.0, self.cf_high*1.5])
        ax[1].set_xlabel('Frequency(Hz)')
        ax[1].set_title('spectrum')
        plt.tight_layout()
        return fig

    def get_ir(self, ir_duration=1, is_env_aligned=False,
               is_fine_aligned=False, delay_common=-1, is_gain_norm=True):
        """Get impulse responses of gammatone filter bank
        Args:
            ir_duration: time length of impulse response (s)
            is_env_aligned: whether to align the envelope of ir,
                            default to False
            is_fine_aligned: whether to align the fine-structure if envelope is
                             aligned, default to False
            delay_common: set delay(s) of filter ir if it is aligned, default
                          to -1, which aligne all filter to the maximum peak
                          position
            is_gain_norm: whether to normlize gains at center frequency
        Returns:
            filter bank impulse response, [n_band,ir_len]
        """
        n_sample = np.int(self.fs*ir_duration)
        # impulse stimuli
        x = np.zeros((n_sample, 1))
        x[100] = 1  # spike
        irs = self.filter(x, is_env_aligned, is_fine_aligned, delay_common,
                          is_gain_norm)
        return irs

    def get_ir_equation(self, t=None):
        """
        """
        if t is None:
            t = np.arange(self.fs)/self.fs

        n_sample = t.shape[0]
        ir = np.zeros((self.n_band, n_sample))

        order = 4

        part1 = t**(order-1)
        for band_i in range(self.n_band):
            part2 = np.multiply(np.cos(2*np.pi*self.cfs[band_i]*t),
                                np.exp(-2*np.pi*self.bws[band_i]*t))
            ir[band_i] = np.multiply(part1, part2)

        return ir
