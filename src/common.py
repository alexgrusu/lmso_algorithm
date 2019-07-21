#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Alexandru - George Rusu (2019). All Rights Reserved.
"""
def time_dif(time, all = False, string = False):
    hours, rest = divmod(time, 3600)
    minutes, seconds = divmod(rest, 60)
    if(all == True):
        return hours, minutes, seconds, str(int(hours)).zfill(2) + ":" + str(int(minutes)).zfill(2) + ":" + str(round(seconds)).zfill(2)
    elif(string == True):
        return str(int(hours)).zfill(2) + ":" + str(int(minutes)).zfill(2) + ":" + str(round(seconds)).zfill(2)
    else:
        return hours, minutes, seconds


def vrms(a):
    return ((a.astype(float)**2).sum() / len(a))**0.5

