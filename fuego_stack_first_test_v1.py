#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 22 13:01:01 2020

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
sta1 = 'FG8' # STATION 
sta2 = 'FG12'
sta3 = 'FG13'

cha_inf = 'BDF' 
cha_seis = 'BHZ'

net = 'GI'  # 

loc_0 = '00'
loc_1 = '01'
loc_2 = '02'  
   
client = Client('138.253.112.23', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23

t1 = UTCDateTime(2019, 10, 1, 0, 50, 10) #the format is year:day_of_the_year:month
length = 60*20
t2 = t1 + length

#%% infrasound stations

st_8_inf = Stream()
st_8_inf = client.get_waveforms(net, sta1, loc_2, cha_inf, t1 , t2)
st_8_inf.detrend(type='linear')
st_8_inf.detrend(type='demean')
st_inf8 = st_8_inf.copy()
fb8 = st_inf8.filter(type='bandpass',freqmin=0.1, freqmax=22)
trace8i=fb8[0]

st_12_inf = Stream()
st_12_inf = client.get_waveforms(net, sta2, loc_2, cha_inf, t1 , t2)
st_12_inf.detrend(type='linear')
st_12_inf.detrend(type='demean')
st_inf12 = st_12_inf.copy()
fb12 = st_inf12.filter(type='bandpass',freqmin=0.1, freqmax=22)
trace12i=fb12[0]

st_13_inf = Stream()
st_13_inf = client.get_waveforms(net, sta3, loc_2, cha_inf, t1 , t2)
st_13_inf.detrend(type='linear')
st_13_inf.detrend(type='demean')
st_inf13 = st_13_inf.copy()
fb13 = st_inf13.filter(type='bandpass',freqmin=0.1, freqmax=22)
trace13i=fb13[0]

#%% seismic stations

sts8 = Stream()
sts8 = client.get_waveforms(net, sta1, loc_0, cha_seis, t1 , t2)
sts8.detrend(type='linear')
sts8.detrend(type='demean')
sts8l = sts8.copy()
sts8h = sts8.copy()
fb8h = sts8h.filter(type='bandpass',freqmin=12, freqmax=22)
fb8l = sts8l.filter(type='bandpass',freqmin=0.5, freqmax=5)
trace8ls=fb8l[0]
trace8hs=fb8h[0]

sts12 = Stream()
sts12 = client.get_waveforms(net, sta2, loc_0, cha_seis, t1 , t2)
sts12.detrend(type='linear')
sts12.detrend(type='demean')
sts12l = sts12.copy()
sts12h = sts12.copy()
fb12h = sts12h.filter(type='bandpass',freqmin=12, freqmax=22)
fb12l = sts12l.filter(type='bandpass',freqmin=0.5, freqmax=5)
trace12ls=fb12l[0]
trace12hs=fb12h[0]

sts13 = Stream()
sts13 = client.get_waveforms(net, sta3, loc_0, cha_seis, t1 , t2)
sts13.detrend(type='linear')
sts13.detrend(type='demean')
sts13l = sts13.copy()
sts13h = sts13.copy()
fb13h = sts13h.filter(type='bandpass',freqmin=10, freqmax=22)
fb13l = sts13l.filter(type='bandpass',freqmin=0.5, freqmax=5)
trace13ls=fb13l[0]
trace13hs=fb13h[0]

#%%

#trace8i.plot(color='r',starttime=t1, endtime=t2)
#trace12i.plot(color='r',starttime=t1, endtime=t2)
#trace13i.plot(color='r',starttime=t1, endtime=t2)
#
#trace8hs.plot(color='k',starttime=t1, endtime=t2)
#trace12hs.plot(color='k',starttime=t1, endtime=t2)
#trace13hs.plot(color='k',starttime=t1, endtime=t2)

#%% sta_lta to find time difference

start= t1 + 10 #time window start 
end= t2

trs8i = trace8i.slice(starttime = start  , endtime= end) 
trs8i.plot(color='r')
sr8i = trace8i.stats.sampling_rate
nsta=int(0.5*sr8i)
nlta=int(25*sr8i)
stream8i=trs8i.data
cft=recursive_sta_lta(stream8i, nsta, nlta)
trig_on=25
trig_off=1
plot_trigger(trs8i, cft, trig_on, trig_off) 
on_off8i = trigger_onset(cft,trig_on,trig_off)

trs12i = trace12i.slice(starttime = start  , endtime= end) 
trs12i.plot(color='r')
sr12i = trace12i.stats.sampling_rate
nsta=int(0.5*sr12i)
nlta=int(25*sr12i)
stream12i=trs12i.data
cft=recursive_sta_lta(stream12i, nsta, nlta)
trig_on=25
trig_off=1
plot_trigger(trs12i, cft, trig_on, trig_off) 
on_off12i = trigger_onset(cft,trig_on,trig_off)

trs13i = trace13i.slice(starttime = start  , endtime= end) 
trs13i.plot(color='r')
sr13i = trace13i.stats.sampling_rate
nsta=int(0.5*sr13i)
nlta=int(25*sr13i)
stream13i=trs13i.data
cft=recursive_sta_lta(stream13i, nsta, nlta)
trig_on=25
trig_off=1
plot_trigger(trs13i, cft, trig_on, trig_off) 
on_off13i = trigger_onset(cft,trig_on,trig_off)

trs8hs = trace8hs.slice(starttime = start  , endtime= end) 
trs8hs.plot(color='k')
sr8s = trace8hs.stats.sampling_rate
nsta=int(1*sr8s)
nlta=int(20*sr8s)
stream8hs=trs8hs.data
cft=recursive_sta_lta(stream8hs, nsta, nlta)
trig_on=8
trig_off=1
plot_trigger(trs8hs, cft, trig_on, trig_off) 
on_off8hs = trigger_onset(cft,trig_on,trig_off)

trs12hs = trace12hs.slice(starttime = start  , endtime= end) 
trs12hs.plot(color='k')
sr12s = trace12hs.stats.sampling_rate
nsta=int(2*sr12s)
nlta=int(25*sr12s)
stream12hs=trs12hs.data
cft=recursive_sta_lta(stream12hs, nsta, nlta)
trig_on=8
trig_off=1
plot_trigger(trs12hs, cft, trig_on, trig_off) 
on_off12hs = trigger_onset(cft,trig_on,trig_off)

trs13hs = trace13hs.slice(starttime = start  , endtime= end) 
trs13hs.plot(color='k')
sr13s = trace13hs.stats.sampling_rate
nsta=int(0.1*sr13s)
nlta=int(30*sr13s)
stream13hs=trs13hs.data
cft=recursive_sta_lta(stream13hs, nsta, nlta)
trig_on=8
trig_off=1
plot_trigger(trs13hs, cft, trig_on, trig_off) 
on_off13hs = trigger_onset(cft,trig_on,trig_off)

#%%

print('')
print('Time of FG8 inf detection =', start + on_off8i[0,0]/sr8i)
print('Time of FG12 inf detection =', start + on_off12i[0,0]/sr12i)
print('Time of FG13 inf detection =', start + on_off13i[0,0]/sr13i)

print('Time of FG8 seis detection =', start + on_off8hs[0,0]/sr8s)
print('Time of FG12 seis detection =', start + on_off12hs[0,0]/sr12s)
print('Time of FG13 seis detection =', start + on_off13hs[0,0]/sr13s)

fg8i_diff = on_off8i[0,0]/sr8i -  on_off8i[0,0]/sr8i
fg12i_diff = on_off8i[0,0]/sr8i - on_off12i[0,0]/sr12i
fg13i_diff = on_off8i[0,0]/sr8i -  on_off13i[0,0]/sr13i

fg8s_diff = on_off8i[0,0]/sr8i -  on_off8hs[0,0]/sr8s
fg12s_diff = on_off8i[0,0]/sr8i -  on_off12hs[0,0]/sr12s
fg13s_diff = on_off8i[0,0]/sr8i -  on_off13hs[0,0]/sr13s

print('')
print('Time diff to FG8i from FG8 inf detection =', fg8i_diff )
print('Time diff to FG8i from FG12 inf detection =',  fg12i_diff)
print('Time diff to FG8i from FG13 inf detection =', fg13i_diff )

print('Time diff to FG8i from FG8 seis detection =', fg8s_diff )
print('Time diff to FG8i from FG12 seis detection =', fg12s_diff )
print('Time diff to FG8i from FG13 seis detection =', fg13s_diff )



#%% stack waveforms

start1i = start - fg8i_diff
start2i = start - fg12i_diff
start3i = start - fg13i_diff

start1s = start - fg8s_diff
start2s = start - fg12s_diff
start3s = start - fg13s_diff

length_trim = length -30

trs8i2 = trace8i.slice(starttime = start1i, endtime= start1i + length_trim ) 
trs12i2 = trace12i.slice(starttime = start2i , endtime= start2i + length_trim) 
trs13i2 = trace13i.slice(starttime = start3i , endtime= start3i + length_trim) 

trs8s2 = trace8hs.slice(starttime = start1s , endtime= start1s + length_trim) 
trs12s2 = trace12hs.slice(starttime = start2s , endtime= start2s + length_trim) 
trs13s2 = trace13hs.slice(starttime = start3s , endtime= start3s + length_trim) 

trs8i2.normalize()
trs12i2.normalize()
trs13i2.normalize()

trs8s2.normalize()
trs12s2.normalize()
trs13s2.normalize()

inf_stack = trs8i2.data + trs12i2.data + trs13i2.data 
seis_stack = trs8s2.data + trs12s2.data# + trs13s2.data 

plt.figure()
plt.plot(inf_stack,'k')
plt.figure()
plt.plot(trs12i2.data,'g')
plt.plot(trs8i2.data,'r')
plt.plot(trs13i2.data,'b')

plt.figure()
plt.plot(seis_stack,'k')
#plt.plot(trs12s2.data,'g')
#plt.figure()
#plt.plot(trs8s2.data,'r')
#plt.plot(trs13s2.data,'b')


#%% sta_lta on stacks



#nsta=int(0.5*50)
#nlta=int(25*50)
#cft=recursive_sta_lta(inf_stack, nsta, nlta)
#trig_on=20
#trig_off=1
#plot_trigger(inf_stack, cft, trig_on, trig_off) 
#on_off_is = trigger_onset(cft,trig_on,trig_off)



# FIX ISSUE OF NOT BEING ABLE TO DO STA_LTA WITH OBSPY - NO STATS WITH NUMPY ARRAY





























































