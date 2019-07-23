#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Alexandru - George Rusu (2019). All Rights Reserved.
"""
import scipy.signal as ss
import numpy as np


def time_dif(time, all = False, string = False):
    hours, rest = divmod(time, 3600)
    minutes, seconds = divmod(rest, 60)
    if(all == True):
        return hours, minutes, seconds, str(int(hours)).zfill(2) + ":" + str(int(minutes)).zfill(2) + ":" + str(round(seconds)).zfill(2)
    elif(string == True):
        return str(int(hours)).zfill(2) + ":" + str(int(minutes)).zfill(2) + ":" + str(round(seconds)).zfill(2)
    else:
        return hours, minutes, seconds


def vrms(seq):
    """
    Computes the root mean square (RMS) power.
    :param seq: the input sequence for which the VRMS is returned.
    """
    return ((seq.astype(float)**2).sum() / len(seq))**0.5


def wgn(length):
    """
    Generates a white Gaussian noise sequence.
    :param length: the length of the output sequence. 
    """
    return np.random.randn(length, 1).ravel()


def ar1(length, beta):
    """
    Generates an AR(1) process resulted by filtering a white Gaussian 
    noise through a first-order system.
    :param length: the length of the output sequence. 
    :param beta: correlation factor.
    """
    b = np.array([1.0])
    a = np.array([1.0, -beta])
    ar1_process = wgn(length)
    ar1_process = ar1_process * np.sqrt(1 - beta**2)
    ar1_process = ss.lfilter(b, a, ar1_process)
    
    return ar1_process
