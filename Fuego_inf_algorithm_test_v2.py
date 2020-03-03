#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 20 16:34:12 2020

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
from obspy.signal.trigger import classic_sta_lta, recursive_sta_lta
from obspy.signal.trigger import plot_trigger, trigger_onset
from obspy import Stream
from numpy import argmax




# STATION, CHANNEL (DDF --> 400 Hz), NETWWORK AND LOCATION CODES 
sta = 'VF03' # STATION 
cha = 'HDF' # CHANNEL
net = 'XZ'  # 
loc = ''    # location, it depends mostly of which network you are in. 
client = Client('138.253.112.23', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23

t1 = UTCDateTime(2018, 5, 22, 0, 0, 0) #the format is year:day_of_the_year:month
t2 = t1 + 60*60*24
st_day = Stream()
st_day = client.get_waveforms(net, sta, loc, cha, t1 , t2)
st_day.detrend(type='linear')
st_day.detrend(type='demean')
sr = st_day[0].stats.sampling_rate
st_day.filter(type='bandpass',freqmin=0.1, freqmax=(sr/2)-(sr/20))
st_day.plot(type='dayplot',starttime=t1, endtime=t2)

trace = st_day[0]

#%%
#window endpoints
start= t1 + 10 #time window start 
end= t2
trs = trace.slice(starttime = start  , endtime= end) 

trs.plot(type='relative',color='b')#, starttime=start , endtime=end)

nsta=int(3*sr)
nlta=int(15*sr)
stream=trs.data
cft=recursive_sta_lta(stream, nsta, nlta)
trig_on=3
trig_off=0.25
plot_trigger(trs, cft, trig_on, trig_off) 
on_off = trigger_onset(cft,trig_on,trig_off)




#%%

#catalogue = np.zeros(shape=(0,1))
num=0

shift2 = 1000
for x in range(0,len(on_off)):
    event = trace[on_off[x,0]-1000:on_off[x,1]+1000]    
#    start = t1 + on_off[x,0]/sr -20
#    end = t1 + on_off[x,1]/sr + 20
#    duration = end-start
#    
#    trace_fft = trace.slice(starttime = start  , endtime= end) 
    
#    window=end-start
#    tr_data=trace_fft.data
#    m=np.mean(tr_data)
#    tr_data = tr_data-m
#    famp = abs(np.fft.fft(tr_data))
#    
#    # double peak ratio
#    fps=int(len(famp)/(sr/0.5))
#    midp=int(len(famp)/(sr/2))
#    spe=int(len(famp)/(sr/4.5))
#    
#    f_peak_m = np.mean(famp[fps:midp])/(window*sr)
#    s_peak_m = np.mean(famp[midp:spe])/(window*sr)
#    peak_r = f_peak_m/s_peak_m
#    
#    # 2Hz drop ratio
#    
#    fpbs=int(len(famp)/(sr/0.5))
#    fpbe=int(len(famp)/(sr/1.8))
#    
#    spbs=int(len(famp)/(sr/2.2))
#    spbe=int(len(famp)/(sr/3.5))
#    
#    peaks_ma = np.mean(famp[fpbs:fpbe])/(window*sr)
#    peaks_mb = np.mean(famp[spbs:spbe])/(window*sr)
#    peaks_m = (peaks_ma + peaks_mb)/2
#    
#    trough_m = np.mean(famp[fpbe:spbs])/(window*sr)    
#    trough_r = peaks_m/trough_m
    
    mean_lev = np.mean(abs(event))
    max_lev = max(event)
    rat_amps = max_lev/mean_lev
    loc_max = event.argmax()
    len_ev = len(event)
    max_loc = (len_ev - loc_max)
    
    
    if len(event)/(sr*60) > 1:
        if max_loc < 2*nlta:
            if rat_amps < 40:
                print(UTCDateTime(start + on_off[x,0]/sr ),'trem with length: ',len(event)/600, 'mins')
                plt.figure()
                plt.plot(event,color='r')
                plt.title(UTCDateTime(start + on_off[x,0]/sr ))
            else:
                plt.figure()
                plt.plot(event,color='k')
                plt.title(rat_amps)
        else:
            print(UTCDateTime(start + on_off[x,0]/sr ),'trem with length: ',len(event)/600, 'mins')
            plt.figure()
            plt.plot(event,color='r')
            plt.title(UTCDateTime(start + on_off[x,0]/sr ))
        
    else:
        print(UTCDateTime(start + on_off[x,0]/sr ),'iso with length: ',len(event)/600, 'mins')
        plt.figure()
        plt.plot(event,color='b')
        plt.title(UTCDateTime(start + on_off[x,0]/sr ))
    

    
        
##%%
#print('number of detections =',num)
#print('number of rejections =',len(on_off)-num)
#
#for x in range(0,len(catalogue)):
#    
#    print(UTCDateTime(catalogue[x]))
#
#
#
#
#













