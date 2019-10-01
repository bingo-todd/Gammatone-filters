import matplotlib.pyplot as plt
import numpy as np
import time
import os
import sys
sys.path.append('..\\')
from gtf import gtf

fig_dir = '../images/gain_normalization/'
if not os.path.exists(fig_dir):
    os.makedirs(fig_dir)


def gain_norm_test():
    fs = 16e3
    gt_filter = gtf(fs,freq_low=80,freq_high=5e3,n_band=16)

    # ir: filter impulse signal
    irs = gt_filter.get_ir(is_gain_norm=False)
    fig = gt_filter.plot_ir_spec(irs)
    savefig(fig,'irs.png')

    # gain normalization
    irs_norm = gt_filter.get_ir()
    fig = gt_filter.plot_ir_spec(irs_norm)
    savefig(fig,'irs_norm.png')

    # delays and gains
    fig = gt_filter.plot_delay_gain_cfs()
    savefig(fig,'delay_gain.png')

if __name__ == '__main__':
    gain_norm_test()
