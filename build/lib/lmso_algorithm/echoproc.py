#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Alexandru - George Rusu (2019). All Rights Reserved.
"""
from . import common
import scipy.signal as ss
import numpy as np


def shift(seq, n = 0):
    """
    Shifts the input sequence with n samples.
    :param n: shift positions.
    """
    a = n % len(seq)
    return np.concatenate((seq[-a:], seq[:-a]), axis = 0)


def normalized_misalignment(estimated_echo_path, real_echo_path):
    """
    Computes the normalized misalignment between the real echo path
    and the estimated echo path.
    :param esimated_echo_path: coeffs of the estimated echo path.
    :param real_echo_path: coeffs of the real echo path.
    """
    return np.linalg.norm((estimated_echo_path).T - real_echo_path) / np.linalg.norm(real_echo_path)


def generate_echo(input_signal, echo_path, ENR):
    """
    Filters the input signal through the echo path and sets the measurement
    noise according to the input echo to noise ratio (ENR).
    :param input_signal: the input signal that must be filtered.
    :param echo_path: fir filter coeffs.
    :param ENR: the echo to noise ratio in dB.
    """
    comfort_noise = np.random.randn(len(input_signal), 1).ravel()
    echo_signal = ss.lfilter(echo_path, 1, input_signal)
    echo_power = common.vrms(echo_signal)
    noise_power = common.vrms(comfort_noise)
    factor = echo_power / noise_power * 10**(-ENR/20)
    measurement_noise = factor * comfort_noise
    echo_signal += measurement_noise
    sgmv = np.var(measurement_noise)
    
    return echo_signal, sgmv