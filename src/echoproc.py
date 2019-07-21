#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Alexandru - George Rusu (2019). All Rights Reserved.
"""
import numpy as np

def shift(seq, n = 0):
    a = n % len(seq)
    return np.concatenate((seq[-a:], seq[:-a]), axis = 0)


