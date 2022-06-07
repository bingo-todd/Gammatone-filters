import os
import sys
sys.path.append('..\\')
from GTF import GTF  # noqa

fig_dir = '../images/gain_normalization/'
if not os.path.exists(fig_dir):
    os.makedirs(fig_dir)


def savefig(fig, fig_name):
    fig_path = f'{fig_dir}/{fig_name}'
    fig.savefig(fig_path)


def gain_norm_test():
    fs = 16e3
    gt_filter = GTF(fs, freq_low=80, freq_high=5e3, n_band=16)

    # ir: filter impulse signal
    irs = gt_filter.get_ir(norm_gain=False)[:, :, 0]
    fig, ax = gt_filter.plot_ir_spec(irs)
    savefig(fig, 'irs.png')

    # gain normalization
    irs_norm = gt_filter.get_ir()[:, :, 0]
    fig, ax = gt_filter.plot_ir_spec(irs_norm)
    savefig(fig, 'irs_norm.png')

    # delays and gains
    fig, ax = gt_filter.plot_delay_gain_cfs()
    savefig(fig, 'delay_gain.png')


if __name__ == '__main__':
    gain_norm_test()
