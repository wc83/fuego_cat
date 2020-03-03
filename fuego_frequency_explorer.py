#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 10:17:43 2020

@author: root
"""

import os
from collections import OrderedDict
import numpy as np
import obspy
import scipy.signal as sgn
import matplotlib.pyplot as plt 
import matplotlib.mlab as mlab
from scipy import integrate
from obspy.clients.earthworm import Client
from obspy import UTCDateTime
from obspy.signal.trigger import trigger_onset
from obspy import Stream
from numpy import argmax


# STATION, CHANNEL (DDF --> 400 Hz), NETWWORK AND LOCATION CODES 
sta = 'FG8' # STATION 
cha = 'BHZ' # CHANNEL
net = 'GI'  # 
loc = '00'    # location, it depends mostly of which network you are in. 

# t1. and t2 are in hours:minutes:seconds
# Get data from (Liverpool Winston default) wave server between times t1 and t2 for all stations in stalist      
client = Client('138.253.113.19', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23
t1 = UTCDateTime(2018, 5, 21, 9, 40, 0) #the format is year:day_of_the_year:month
#t1=UTCDateTime(1521072000)
t2 = t1 + 120
sts = Stream()
sts = client.get_waveforms(net, sta, loc, cha, t1 , t2)
#print(st)

#2014-11-25 20:20:16 (UTC)

sts.detrend(type='linear')
sts.detrend(type='demean')

sts1 = sts.copy()
sts2 = sts.copy()
sts3 = sts.copy()
sts4 = sts.copy()
sts5 = sts.copy()
sts6 = sts.copy()
sts7 = sts.copy()


fb1 = sts1.filter(type='bandpass',freqmin=0.5, freqmax=3)
fb1.plot(color='k',starttime=t1, endtime=t2)

sts.filter(type='bandpass',freqmin=0.5, freqmax=10)
sts.plot(color='b')

fb2 = sts2.filter(type='bandpass',freqmin=12, freqmax=20)
fb2.plot(color='k',starttime=t1, endtime=t2)


#%%

t1 = UTCDateTime(2018, 5, 21, 0, 0, 0)
t2 = t1 + 1*60*60
st2 = client.get_waveforms(net, sta, loc, cha, t1 , t2)


st2.detrend(type='linear')
st2.detrend(type='demean')

st2a = st2.copy()
st2b = st2.copy()

st2a.filter(type='bandpass',freqmin=12, freqmax=20)
st2a.plot(type='dayplot',interval=30,starttime=t1, endtime=t2)

st2.plot(type='dayplot',interval=30,starttime=t1, endtime=t2)

st2b.filter(type='bandpass',freqmin=0.5, freqmax=3)
st2b.plot(type='dayplot',interval=30,starttime=t1, endtime=t2)















