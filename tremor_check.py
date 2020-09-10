#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 24 12:22:02 2020

@author: root
"""

import obspy
from obspy import read
import numpy as np
from numpy import genfromtxt
import matplotlib.pyplot as plt
from obspy import Stream
from obspy.signal.trigger import classic_sta_lta, recursive_sta_lta
from obspy.signal.trigger import plot_trigger, trigger_onset
from obspy import UTCDateTime
from matplotlib import pyplot
import scipy
from obspy import UTCDateTime
import datetime as dt
import matplotlib.dates as md




cat = genfromtxt("/Users/william/Documents/Fuego_catalogue/Fuego_trem_v3_inf_test4b.csv", delimiter=',',skip_header=1)

print("Inf Cat:")
for x in range(0,len(cat)):
    print(cat[x,0])
    print(UTCDateTime(cat[x,34])," to ", UTCDateTime(cat[x,37]), '{:.1f}'.format((cat[x,37] - cat[x,34])/60) )
    print(UTCDateTime(cat[x,35])," to ", UTCDateTime(cat[x,36]), '{:.1f}'.format((cat[x,36] - cat[x,35])/60 ))
    print("")





cat = genfromtxt("/Users/william/Documents/Fuego_catalogue/Fuego_trem_v3_seis_test4b.csv", delimiter=',',skip_header=1)
print("")
print("Seis Cat:")
for x in range(0,len(cat)):
    print(cat[x,0])
    print(UTCDateTime(cat[x,9])," to ", UTCDateTime(cat[x,12]), '{:.1f}'.format((cat[x,12] - cat[x,9])/60) )
    print(UTCDateTime(cat[x,10])," to ", UTCDateTime(cat[x,11]), '{:.1f}'.format((cat[x,11] - cat[x,10])/60 ))
    print("")














































































