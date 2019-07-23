#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Alexandru - George Rusu (2019). All Rights Reserved.
"""
from lmso_algorithm import algorithm
from lmso_algorithm import common
from lmso_algorithm import echoproc
import matplotlib.pyplot as plt
import scipy.io as sio
import numpy as np


fs = 8000 # sample rate
ford = 512 # adpative filter length
sim_len = 24000 # Simulation length in sample count
test_signal = 0 # 0 for WGN as input signal, 1 for AR(1) process
beta = 0.9 # a1 coeff of the first order system
ENR = 20 # echo-to-noise-ratio
nmse_w = np.zeros([sim_len, 1]) # computed normalized misalignment
nmse_g = np.zeros([sim_len, 1]) # computed normalized misalignment
data = sio.loadmat('ir_1.mat') # measured echo path
echo_path = data['g1'].ravel()[0 : ford] # read only the first 512 coeffs
eps = abs(echo_path[0]) # initialization constant


if(test_signal):
    far_end_signal = common.ar1(sim_len, beta)
else:
    far_end_signal = common.wgn(sim_len)

echo_signal, sgmv = echoproc.generate_echo(far_end_signal, echo_path, ENR)


if __name__ == "__main__":
    lmsow_obj = algorithm.lmso(ford = ford, cst = eps, lam = 0.95, sgmv = sgmv, sgmw = 0.0, lniproc = 2, alpha = 1.0, delta = 1e-6)
    lmsog_obj = algorithm.lmso(ford = ford, cst = eps, lam = 0.95, sgmv = sgmv, sgmw = 0.0, lniproc = 2, alpha = 1.0, delta = 1e-6)
    
    for i in range(0, sim_len):
        lmsow_obj.x.cb_push(far_end_signal[i])
        lmsog_obj.x.cb_push(far_end_signal[i])
        ew = lmsow_obj.lmso_w(echo_signal[i])
        eg = lmsog_obj.lmso_g(echo_signal[i])
        nmse_w[i] = echoproc.normalized_misalignment(lmsow_obj.h.coeffs, echo_path)
        nmse_g[i] = echoproc.normalized_misalignment(lmsog_obj.h.coeffs, echo_path)
    
    t = np.linspace(0, len(far_end_signal)/fs, len(far_end_signal))
    fig = plt.figure()
    line_1, = plt.plot(t, 20 * np.log10(nmse_w), label = r'$LMSO-W$')
    line_2, = plt.plot(t, 20 * np.log10(nmse_g), label = r'$LMSO-G$')
    plt.legend()