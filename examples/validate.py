import matplotlib.pyplot as plt
import numpy as np
import sys
import os
import time
sys.path.append('..\\')
from gtf import gtf as gtf_proposed
from gammatone import filters as gtf_reference

fig_dir = '../images/validate/'
if not os.path.exists(fig_dir):
    os.makedirs(fig_dir)

def savefig(fig,fig_name):
    fig_fpath = os.path.join(fig_dir,fig_name)
    fig.savefig(fig_fpath)

def validate():
    fs = 16e3
    x_len = np.int16(fs)

    x = np.zeros(x_len)
    x[200] = 1

    gtf_obj = gtf_proposed(fs,cf_low=100,cf_high=2000,n_band=4)
    irs1 = gtf_obj.filter(x)
    fig1 = gtf_obj.plot_ir_spec(irs1[:,:1000])
    savefig(fig1,'proposed.png')

    print(gtf_obj.cfs)
    
    coefs = gtf_reference.make_erb_filters(fs,gtf_obj.cfs)
    irs2 = gtf_reference.erb_filterbank(x,coefs)
    fig2 = gtf_obj.plot_ir_spec(irs2[:,:1000])
    savefig(fig2,"reference.png")

if __name__ == "__main__":
    validate()
