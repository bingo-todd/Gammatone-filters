import matplotlib.pyplot as plt
import numpy as np
import sys
import os
sys.path.append('..\\')
from gtf import gtf

def filter_spectrum():
    fs = 16e3

    gt_filter = gtf(fs,cf_low=4e3,cf_high=4e3,n_band=1)

    order = 4
    freq_bins = np.arange(1,fs/2)
    n_freq_bin = freq_bins.shape[0]

    cf = 4e3
    bw = gt_filter.cal_ERB(cf)
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
    ax2.plot([0,fs/2/1000],[0,0],'-.',color='black')
    ax2.set_ylabel('delay(ms)',color=color )
    ax2.tick_params(axis='y',labelcolor=color)

    fig_dir = '..\images\\filter_spectrum'
    if not os.path.exists(fig_dir):
        os.mkdir(fig_dir)
    fig.savefig(os.path.join(fig_dir,'filter_spectrum.png'))

if __name__ == "__main__":
    filter_spectrum()
