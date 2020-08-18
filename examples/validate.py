import numpy as np
import os
import matplotlib.pyplot as plt
from GTF import GTF as gtf_proposed
from gammatone import filters as gtf_reference

fig_dir = 'images/validate/'
if not os.path.exists(fig_dir):
    os.makedirs(fig_dir)


def savefig(fig, fig_name):
    fig_fpath = os.path.join(fig_dir, fig_name)
    fig.savefig(fig_fpath)


def validate():
    fs = 16e3
    # impulse
    x_len = np.int16(fs)
    x = np.zeros(x_len)
    x[1] = 1

    gtf_obj = gtf_proposed(fs, cf_low=100, cf_high=2000, n_band=4)
    irs = gtf_obj.filter(x)
    # fig1 = gtf_obj.plot_ir_spec(irs1[:, :1000])
    # savefig(fig1, 'proposed.png')

    coefs = gtf_reference.make_erb_filters(fs, gtf_obj.cfs)
    irs_ref = gtf_reference.erb_filterbank(x, coefs)
    # fig2 = gtf_obj.plot_ir_spec(irs2[:, :1000])
    # savefig(fig2, "reference.png")

    irs_eq = gtf_obj.get_ir_equation()

    fig, ax = plt.subplots(1, 2, sharex=True, sharey=True, tight_layout=True)
    ax[0].plot(irs[3]/np.max(irs[3]), label='todd')
    ax[0].plot(irs_eq[3]/np.max(irs_eq[3]), label='eq')
    ax[0].legend()
    ax[0].set_xlim([0, 200])

    ax[1].plot(irs_ref[3]/np.max(irs_ref[3]), label='detly')
    ax[1].plot(irs_eq[3]/np.max(irs_eq[3]), label='eq')
    ax[1].legend()
    ax[1].set_xlim([0, 200])
    savefig(fig, 'compare.png')


if __name__ == "__main__":
    validate()
