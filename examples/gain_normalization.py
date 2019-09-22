import matplotlib.pyplot as plt
import numpy as np
import time
import os
import sys
sys.path.append('..\\')
from gtf import gtf

def gain_norm_test():
    fs = 16e3
    gt_filter = gtf(fs,freq_low=80,freq_high=5e3,n_band=16)

    # ir: filter impulse signal
    irs = gt_filter.get_ir(is_gain_norm=False)

    # gain normalization
    irs_norm = gt_filter.get_ir()

    # plot ir
    fig_dir = '..\images\gain_normalization'
    if not os.path.exists(fig_dir):
        os.mkdir(fig_dir)

    fig = gt_filter.plot_ir_spec(irs,gt_filter.fs,gt_filter.cfs)
    fig.savefig(os.path.join(fig_dir,'irs.png'))

    fig = gt_filter.plot_ir_spec(irs_norm,gt_filter.fs,gt_filter.cfs)
    fig.savefig(os.path.join(fig_dir,'irs_norm.png'))

    # delays and gains
    fig = plt.figure()
    gt_filter.cal_delay_gain_cfs(is_plot=True,fig=fig)
    fig.savefig(os.path.join(fig_dir,'delay_gain.png'))

if __name__ == '__main__':
    gain_norm_test()
