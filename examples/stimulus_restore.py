import matplotlib.pyplot as plt
import numpy as np
import sys
import os
sys.path.append('..\\')
from GTF import GTF  # noqa

fig_dir = r'images/stimulus_restore'
if not os.path.exists(fig_dir):
    os.makedirs(fig_dir)


def savefig(fig, fig_name):
    fig_path = os.path.join(fig_dir, fig_name)
    print(fig_path)
    fig.savefig(fig_path)


def stimulus_restore_plot(irs, fs):
    t = np.arange(irs.shape[1])/fs*1e3
    fig = plt.figure()

    ha = fig.subplots(2, 1)
    ha[0].plot(t, irs.T)
    ha[0].set_title('irs')
    ha[0].set_xlabel('time(ms)')

    ha[1].plot(t, np.sum(irs, axis=0))
    ha[1].set_title('irs_sum')
    ha[1].set_xlabel('time(ms)')
    plt.tight_layout()
    return fig


def stimulus_restore_test():
    fs = 16e3
    gt_filter = GTF(fs, freq_low=80, freq_high=5e3, n_band=16)

    ir_len = np.int16(50e-3*fs)

    # gain normalization for all
    # not aligned
    ir = gt_filter.get_ir()[:, :, 0]
    fig = stimulus_restore_plot(ir[:, :ir_len], fs)
    savefig(fig, 'irs.png')

    # envelope aligned
    ir_align_enved = gt_filter.get_ir(align_env=True)[:, :, 0]
    fig = stimulus_restore_plot(ir_align_enved[:, :ir_len], fs)
    savefig(fig, 'irs_align_enved.png')

    # envelope & fine structure aligned
    ir_all_aligned = gt_filter.get_ir(align_env=True, align_fine=True)[:, :, 0]
    fig = stimulus_restore_plot(ir_all_aligned[:, :ir_len], fs)
    savefig(fig, 'irs_all_aligned.png')


if __name__ == "__main__":
    stimulus_restore_test()
