import time
import ctypes
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from numpy.ctypeslib import ndpointer

class CParamsType:
    """Data type convertor for c function arguments"""
    def from_param(self,param): # called by ctypes
        typename = type(param).__name__
        if hasattr(self,'from_'+typename):
            return getattr(self,'from_'+typename)(param)
        elif isinstance(param,ctypes.Array):
            return param
        else:
            raise TypeError("Can't convert {}".format(typename))

    def from_bool(self,param):
        val = int(param)
        return val

    def from_array(self,param):
        if param.typecode != 'd':
            raise TypeError('must be an array of doubles')
        ptr,_ = param.buffer_info()
        return ctypes.cast(ptr,ctypes.POINTER(ctypes.c_double))

    def from_list(self,param):
        val = ((ctypes.c_double*len(param)))(*param)
        return val

    from_tuple = from_list #

    def from_ndarray(self,param):
        return param.ctypes.data_as(ctypes.POINTER(ctypes.c_double))


class gtf:
    """Python interface for all-pole gammatone filters wrote using C
        a python version filters is also included"""

    _package_dir = os.path.dirname(os.path.abspath(__file__))
    _lib_path = os.path.join(_package_dir,'libgtf.so')

    def __init__(self,fs,cf_low=None,cf_high=None,
                 freq_low=None,freq_high=None,
                 n_band=1):
        """
        Args:
            fs: sample frequency
            cf_low,cf_high: the lowest and highest center frequency
            freq_low,freq_high: the lowest and highest cutoff frequency
            n_band: frequency bands number
        """
        # args check
        if cf_low is None:
            if freq_low is not None:
                # set freq_low as the lowest cutoff frequency
                cf_low = (2*freq_low+24.7)/(2-24.7*4.37/1000)
            else:
                raise Exception('neither cf_low or freq_low is specified')
        if cf_high is None:
            if freq_high is not None:
                # make freq_high as the highest cutoff frequency
                cf_high = (2*freq_high-24.7)/(2+24.7*4.37/1000)
            else:
                raise Exception('neither cf_high or freq_high is specified')

        # center frequencies
        cfs = self.divide_freq_space(freq_low=cf_low,
                                     freq_high=cf_high,
                                     N=n_band)
        # bandwidths
        bws = self.cal_ERB(cfs)
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
        if not os.path.exists(self._lib_path):
            c_file_path = os.path.join(self._package_dir,'gtf.c')
            if not os.path.exists(c_file_path):
                raise OSError(0,'missing file','libgtf.so and gtf.c')
            else:
                print('missing libgtf.so, compile from C code')
                os.system('gcc -fPIC -shared {} -o {}'.format(c_file_path,
                                                                self._lib_path))
                print('compile success')

        _cmodel = ctypes.cdll.LoadLibrary(self._lib_path)
        self._gt_filter = _cmodel.gtf

        CParams = CParamsType()
        #double* gtf(double*x, int x_len,
        #            int fs, double*cfs, double*bs, int N_band,
        #            int is_env_aligned, int is_fine_aligned, int delay_common,
        #            int is_gain_norm)
        self._gt_filter.argtypes = (CParams,ctypes.c_int,
                                    ctypes.c_int,CParams,CParams,ctypes.c_int,
                                    CParams,CParams,ctypes.c_int,
                                    CParams)
        self._gt_filter.restype = ctypes.POINTER(ctypes.c_double)

        self._free_mem = _cmodel.free_mem
        self._free_mem.argtypes = (ctypes.POINTER(ctypes.c_double),)
        self._free_mem.restype = None


    def divide_freq_space(self,freq_low,freq_high,N,divide_type='ERB'):
        """Divide frequency range (freq_low~freq_high) equally in erb
        scale(default)
        Args:
            freq_low: low bound of frequency range
            freq_high: high bound of frequency range
            N: segments number frequency range to be divided
            divide_type: default to ERB
        """
        if divide_type is 'ERB':
            if N == 1:
                return np.asarray(freq_low,dtype=np.float).reshape(1,)
            low_erb = self.Hz2ERbwscal(freq_low)
            high_erb = self.Hz2ERbwscal(freq_high)
            erb_elem = (high_erb-low_erb)/(N-1)
            f = self.ERbwscal2Hz(low_erb+erb_elem*np.arange(N))
        else:
            raise Exception('unsupport Divide type')
        return f


    def Hz2ERbwscal(self,freq):
        """convert Hz to ERB scale"""
        return 21.4*np.log10(4.37*freq/1e3+1)


    def ERbwscal2Hz(self,erb_num):
        """convert ERB scale to Hz"""
        return (10**(erb_num/21.4)-1)/4.37*1e3


    def cal_ERB(self,cf):
        """calculate the ERB(Hz) of given center frequency based on equation
        given by Glasberg and Moore
        Args
            cf: center frequency Hz
        """
        return 24.7*(4.37*cf/1000+1.0)


    def filter_py(self,x,is_aligned=False,delay_common=None):
        """Filters in Python
        Args:
        x: signal with shape of [x_len,n_band], if x only has single dimension,
        n_band will be added as 1
        is_aligned: aligned peaks of Gammatone filter impulse response
        delay_common: if aligned, give the same delay to all channels,
        default, aligned to the maximum delay
        Returns:
        fitler result with the shape of [n_band,x_len,n_band]
        """
        #constant variables
        tpt = 2*np.pi*(1.0/self.fs)

        if is_aligned is True:
            delays = np.round(3.0/(2.0*np.pi*self.bws)*self.fs)/self.fs
            if delay_common is not None:
                delay_common = np.max(delays)
                delays = np.int(delays-delay_common)
        else:
            delays = np.zeros(self.n_band)

        x_len,n_chann = x.shape
        x_filtered = np.zeros((self.n_band,x_len,n_chann))

        # IIR and FIR filters outputs
        out_a = np.zeros((5,n_chann),dtype=np.complex); coefs_a = np.zeros(5);
        out_b = 0; coefs_b = np.zeros(4)

        # padd 0 to x, for the convinientce of delaying manipulation
        x_padd = np.concatenate((x,np.zeros((np.max(np.int16(delays*self.fs)),
                                            n_chann))),
                                axis=0)
        norm_factors = np.zeros(self.n_band)
        y = np.zeros((self.n_band,x_len,n_chann),dtype=np.float)
        for band_i in range(self.n_band):
            bw = self.bws[band_i]
            cf = self.cfs[band_i]
            k = np.exp(-tpt*bw)
            # filter coefs
            coefs_a[0]=1; coefs_a[1]=4*k; coefs_a[2]=-6*k**2;
            coefs_a[3]=4*k**3; coefs_a[4]=-k**4;
            coefs_b[0]=1; coefs_b[1]=1; coefs_b[2]=4*k; coefs_b[3]=k**2;
            #
            norm_factors[band_i] = (1-k)**4/(1+4*k+k**2)*2
            delay_len_band = np.int(delays[band_i]*self.fs)
            phi_0 = -2*np.pi*cf*np.int(self.fs*3.0/(2.0*np.pi*self.bws[band_i]))/self.fs
            #3.0/(2.0*np.pi*self.bws[band_i])
            # phi_0 = 0
            for sample_i in range(x_len+delay_len_band):
                # 2*pi*cf*n/fs -2*pi*cf*delay
                freq_shiftor = np.exp(-1j*(tpt*cf*sample_i))
                # IIR part
                out_a[0,:] = x_padd[sample_i,:]*freq_shiftor*np.exp(1j*phi_0)
                for order_i in range(1,5):
                    out_a[0,:] = out_a[0,:] + coefs_a[order_i]*out_a[order_i,:]
                # FIR part
                out_b = 0
                for order_i in range(1,4):
                    out_b = out_b+coefs_b[order_i]*out_a[order_i,:]
                    #
                if sample_i > delay_len_band:
                    y[band_i,sample_i-delay_len_band,:] = norm_factors[band_i]*\
                                    np.real(out_b*np.conjugate(freq_shiftor))
                # update IIR output
                for order_i in range(4,0,-1):
                    out_a[order_i,:] = out_a[order_i-1,:]
        return np.squeeze(y)


    def filter_c(self, x, is_env_aligned=False, is_fine_aligned=False,
                 delay_common=-1, is_gain_norm=True):
        """filter input signal using c library
        Args:
            x: signal with shape of [x_len,n_band], if x only has single
               dimension, n_band will be added as 1
            is_env_aligned: whether to align the envelope of ir,default to False
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
        if not isinstance(x,np.ndarray): #
            raise Exception()

        if len(x.shape) > 2:
            raise Exception('two many dimensions for x')

        # ensure x is 2-D array
        x_len = x.shape[0]
        if len(x.shape) == 1:
            x.shape = (x_len,1)
        n_chann = x.shape[1]

        # convert bool values to int values(0,1) as the inputs of c function
        x_filtered = np.zeros((self.n_band,x_len,n_chann))
        for chann_i in range(n_chann):
            # data in slice are not stored in continue memmory
            x_chann = np.copy(x[:,chann_i])
            temp =  self._gt_filter(x_chann,x_chann.shape[0],
                                    self.fs,self.cfs,self.bws,self.n_band,
                                    is_env_aligned,is_fine_aligned,delay_common,
                                    is_gain_norm)
            x_filtered[:,:,chann_i] = np.array(np.fromiter(temp,dtype=np.float64,count=x_len*self.n_band)).reshape(self.n_band,x_len)
            self._free_mem(temp) # avoid memmory leak
        return np.squeeze(x_filtered)


    def filter_spectrum(self):

        order = 4
        freq_bins = np.arange(1,self.fs/2)
        n_freq_bin = freq_bins.shape[0]

        cf =4e3
        bw = self.cal_ERB(cf)
        gain_funcs = 6/((2*np.pi*bw)**order)*\
                            (np.divide(1,1+1j*(freq_bins-cf)/bw)**order+
                             np.divide(1,1+1j*(freq_bins+cf)/bw)**order)
        amp_spectrum = np.abs(gain_funcs)

        phase_spectrum1 = np.flip(np.unwrap(np.angle(np.flip(gain_funcs[:4001]))))
        phase_spectrum2 = np.unwrap(np.angle(gain_funcs[4000:]))
        phase_spectrum = np.concatenate((phase_spectrum1[:4000],phase_spectrum2))
        delays = np.divide(phase_spectrum,freq_bins)

        fig,ax1 = plt.subplots()
        linewidth = 2

        color = 'tab:red'
        ax1.semilogy(freq_bins/1000,amp_spectrum,color=color,linewidth=linewidth)
        ax1.set_ylabel('dB',color=color )
        ax1.set_xlabel('frequency(kHz)')
        ax1.tick_params(axis='y',labelcolor=color)
        ax1.set_title('cf=4kHz')

        ax2 = ax1.twinx()
        color='tab:blue'
        ax2.plot(freq_bins/1000,delays,color=color,linewidth=linewidth)
        # ax2.get_x
        ax2.plot([4,4],[-8,8],'-.',color='black')
        ax2.plot([0,self.fs/2/1000],[0,0],'-.',color='black')
        ax2.set_ylabel('delay(ms)',color=color )
        ax2.tick_params(axis='y',labelcolor=color)

        fig.savefig('images/filter_spectrum.png')


    def cal_delay_gain_cfs(self,is_plot=False,fig=None):
        """ Calculate delay and center-frequency gain of gammatone filter
        before alignment and gain-normalization
        Returns:
            [delays,gains]
            the delays and gains at center freuqency of each frequency band
            delays: in milisecond
        """
        k = np.exp(-2*np.pi/self.fs*self.bws);
        Q = np.divide(self.cfs,self.bws) # quality of filter

        temp1 = 8*Q*(1-4*Q**2)
        temp2 = np.multiply(Q**2,Q*(16*Q**2-24))
        phi_delay = np.arctan(np.divide(temp1,temp2+2))-np.arctan(np.divide(temp1,temp2+1));
        delays = phi_delay/(2*np.pi)*1e3

        correct_factor = np.sqrt(((temp2+1)**2+temp1**2)/((temp2+2)**2+temp1**2))
        gains = 10*np.log10(3/(2*np.pi*self.bws))*correct_factor

        if is_plot:
            if fig is None:
                fig = plt.figure(figsize=[8,3])#,dpi=100)
            else:
                fig.set_size_inches(8,3)

            axes = fig.subplots(1,2)
            axes[0].plot(self.cfs,delays); axes[0].set_yscale('log')
            axes[0].set_xlabel('Frquency(Hz)'); axes[0].set_ylabel('Delay(ms)');
            axes[1].plot(self.cfs,gains);
            axes[1].set_xlabel('Frquency(Hz)'); axes[1].set_ylabel('Gain(dB)');
            plt.tight_layout()

        return [delays,gains]


    def plot_ir_spec(self,ir,fs=None,cfs=None,fig=None,title='ir'):
        """plot the waveform and spectrum of given impulse response
        Args:
            ir: impulse response
            fs: sample frequency,use self.fs as default
            fig: handle of matplotlib figure, if not given, not figure will be created
            title: title for ir waveform sub-panel
        """

        if fs == None:
            fs = self.fs

        ir_len = ir.shape[1]
        N_fft = ir_len
        N_fft_half = np.int(N_fft/2)

        index = np.flip(np.argsort(cfs))
        cfs = cfs
        ir = ir[index,:]

        spec = np.abs(np.fft.fft(ir,N_fft,axis=1))[:,:N_fft_half]

        time_ticks = np.arange(ir_len)/self.fs
        freq_ticks = np.arange(N_fft_half)/N_fft*self.fs
        if fig is None:
            fig = plt.figure(figsize=[8,3])
        x_lim_max = 0.08
        linewidth = 2
        axes = fig.subplots(1,2)
        axes[0].plot(time_ticks,ir.T,linewidth=linewidth);
        axes[0].set_xlabel('Time(s)'); axes[0].set_title(title)
        axes[0].set_xlim([0,x_lim_max])
        axes[0].legend(['{:.0f}'.format(cf) for cf in cfs])

        axes[1].plot(freq_ticks,spec.T,linewidth=linewidth);
        axes[1].set_xlim([self.cf_low/8.0,self.cf_high*1.5])
        axes[1].set_xlabel('Frequency(Hz)');
        axes[1].set_title('spectrum')
        axes[1].legend(['{:.0f}'.format(cf) for cf in cfs])
        plt.tight_layout()

        return fig


    def get_ir(self,ir_duration=1,
               is_env_aligned=False,is_fine_aligned=False,delay_common=-1,
               is_gain_norm=True):
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
        N_sample = np.int(self.fs*ir_duration)
        # impulse stimuli
        x = np.zeros((N_sample,1))
        x[200] = 1# spike
        irs = self.filter_c(x,is_env_aligned,is_fine_aligned,delay_common,
                            is_gain_norm)
        return irs


    def get_ir_equation(self,t=None):
        """
        """
        if t is None:
            t = np.arange(self.fs)/self.fs

        N_sample = t.shape[0]
        ir = np.zeros((self.n_band,N_sample))

        order = 4

        part1 = t**(order-1)
        for band_i in range(self.n_band):
            part2 = np.multiply(np.cos(2*np.pi*self.cfs[band_i]*t),np.exp(-2*np.pi*self.bws[band_i]*t))
            ir[band_i] = np.multiply(part1,part2)


        delay_size = np.int16(20e-3*self.fs)
        ir = np.concatenate((np.zeros((self.n_band,delay_size)),ir[:,:-delay_size]),axis=1)
        return ir


def example():
    fs = 16e3
    gt_filter = gtf(fs,freq_low=80,freq_high=1e3,n_band=4)

    # delays and gains
    fig = plt.figure()
    gt_filter.cal_delay_gain_cfs(is_plot=True,fig=fig)
    fig.savefig('images/delay_gain.png')

    # impulse response direct from equation
    ir_equation = gt_filter.get_ir_equation()
    ir_equation = ir_equation/np.max(np.abs(ir_equation))

    # ir: filter impulse signal
    ir = gt_filter.get_ir(is_gain_norm=False)
    # print('ir',ir)
    ir = ir/np.max(np.abs(ir))

    # gain normalization
    ir_norm = gt_filter.get_ir()

    # ir envelope compensated
    ir_norm_env_aligned = gt_filter.get_ir(is_env_aligned=True)
    ir_norm_all_aligned = gt_filter.get_ir(is_env_aligned=True,
                                           is_fine_aligned=True)

    # plot ir
    if not os.path.exists('images'):
        os.mkdir('images')

    fig_ir_eq = gt_filter.plot_ir_spec(ir_equation,gt_filter.fs,gt_filter.cfs)
    fig_ir_eq.savefig('images/ir_equation.png')

    fig_ir = gt_filter.plot_ir_spec(ir,gt_filter.fs,gt_filter.cfs)
    fig_ir.savefig('images/ir.png')

    fig_norm = gt_filter.plot_ir_spec(ir_norm,gt_filter.fs,gt_filter.cfs)
    fig_norm.savefig('images/ir_norm.png')

    fig_ir_norm_env_aligned = gt_filter.plot_ir_spec(ir_norm_env_aligned,
                                                     gt_filter.fs,
                                                     gt_filter.cfs)
    fig_ir_norm_env_aligned.savefig('images/ir_norm_env_aligned.png')

    fig_ir_norm_all_aligned = gt_filter.plot_ir_spec(ir_norm_all_aligned,
                                                     gt_filter.fs,
                                                     gt_filter.cfs)
    fig_ir_norm_all_aligned.savefig('images/ir_norm_all_aligned.png')

    gt_filter.filter_spectrum()

def efficiency_check():
    fs = 16e3
    gt_filter = gtf(fs,freq_low=80,freq_high=5e3,n_band=16)

    ir_duration = 1
    N_sample = np.int(fs*ir_duration)
    x = np.zeros((N_sample,1))
    x[200] = 1

    t_start = time.time()
    ir_c = gt_filter.filter_c(x);
    t_comsum_c = time.time()-t_start

    t_start = time.time()
    ir_py = gt_filter.filter_py(x)
    t_comsum_py = time.time()-t_start

    print('time consumed(s) \n \t {:<10}:{:.2f} \n\t {:<10}:{:.2f}'.format('c',t_comsum_c,'python',t_comsum_py))

if __name__ == '__main__':
    if 'example' in sys.argv:
        example()
    if 'efficiency' in sys.argv:
        efficiency_check()

    if 'ir_python' in sys.argv:
        fs = 16e3
        gt_filter = gtf(fs,freq_low=80,freq_high=1e3,n_band=4)
        ir_duration = 1
        N_sample = np.int(fs*ir_duration)
        x = np.zeros((N_sample,1))
        x[200] = 1
        irs = gt_filter.filter_py(x,is_aligned=True)
        fig = gt_filter.plot_ir_spec(irs,gt_filter.fs,gt_filter.cfs)
        fig.savefig('test.png')