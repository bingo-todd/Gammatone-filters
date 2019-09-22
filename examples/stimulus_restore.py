import matplotlib.pyplot as plt
import numpy as np
import sys
import os
sys.path.append('..\\')
from gtf import gtf

def stimulus_restore_plot(irs,fs):
    t = np.arange(irs.shape[1])/fs*1e3
    fig = plt.figure()

    ha = fig.subplots(2,1)
    ha[0].plot(t,irs.T)
    ha[0].set_title('irs')
    ha[0].set_xlabel('time(ms)')

    ha[1].plot(t,np.sum(irs,axis=0))
    ha[1].set_title('irs_sum')
    ha[1].set_xlabel('time(ms)')

    return fig


def stimulus_restore_test():
    fs = 16e3
    gt_filter = gtf(fs,freq_low=80,freq_high=5e3,n_band=16)

    # gain normalization for all
    # not aligned
    ir = gt_filter.get_ir()
    # envelope aligned
    ir_env_aligned = gt_filter.get_ir(is_env_aligned=True)
    # envelope & fine structure aligned
    ir_all_aligned = gt_filter.get_ir(is_env_aligned=True,
                                           is_fine_aligned=True)

    fig_dir = '..\images\stimulus_restore'
    if not os.path.exists(fig_dir):
        os.makedirs(fig_dir)

    ir_len = 1000
    fig = stimulus_restore_plot(ir[:,:ir_len],fs)
    fig.savefig(os.path.join(fig_dir,'irs.png'))

    fig = stimulus_restore_plot(ir_env_aligned[:,:ir_len],fs)
    fig.savefig(os.path.join(fig_dir,'irs_env_aligned.png'))

    fig = stimulus_restore_plot(ir_all_aligned[:,:ir_len],fs)
    fig.savefig(os.path.join(fig_dir,'irs_all_aligned.png'))

if __name__ == "__main__":
    stimulus_restore_test()
