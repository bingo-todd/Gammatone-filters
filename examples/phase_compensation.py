import matplotlib.pyplot as plt
import numpy as np
import sys
import os
sys.path.append('..\\')
from gtf import gtf

def phase_compensation_plot(irs,fs):
    t = np.arange(irs.shape[1])/fs*1e3
    fig = plt.figure()
    ha = fig.subplots(1,1)
    ha.plot(t,np.flipud(irs).T,linewidth=3) # plot high-cf filter ir first
    ha.set_xlabel('time(ms)')
    return fig


def phase_compensation_test():
    """"""
    fs = 16e3
    gtf_filter = gtf(fs=fs,freq_low=80,freq_high=1e3,n_band=4)

    irs = gtf_filter.get_ir()
    irs_env_aligned = gtf_filter.get_ir(is_env_aligned=True,
                                        delay_common=0)
    irs_all_aligned = gtf_filter.get_ir(is_env_aligned=True,
                                        is_fine_aligned=True,
                                        delay_common=0)

    ir_len = 1000
    fig_dir = '..\images\phase_compensation'
    if not os.path.exists(fig_dir):
        os.makedirs(fig_dir)

    fig = phase_compensation_plot(irs[:,:ir_len],fs)
    fig.savefig(os.path.join(fig_dir,'irs.png'))

    fig = phase_compensation_plot(irs_env_aligned[:,:ir_len],fs)
    fig.savefig(os.path.join(fig_dir,'irs_env_aligned.png'))

    fig = phase_compensation_plot(irs_all_aligned[:,:ir_len],fs)
    fig.savefig(os.path.join(fig_dir,'irs_all_aligned.png'))


if __name__ == "__main__":
    phase_compensation_test()
