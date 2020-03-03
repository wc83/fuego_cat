#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 16:39:25 2020

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


#%% Gas rich reference

# STATION, CHANNEL (DDF --> 400 Hz), NETWWORK AND LOCATION CODES 
sta = 'VF03' # STATION 
cha = 'HDF' # CHANNEL
net = 'XZ'  # 
loc = ''    # location, it depends mostly of which network you are in. 
client = Client('138.253.112.23', 16022) # ip, port - ip's 138.253.113.19 or 138.253.112.23

t1r = UTCDateTime(2018, 5, 22, 8, 4, 6) #the format is year:day_of_the_year:month
t2r = t1r + 20
gr_r = Stream()
gr_r = client.get_waveforms(net, sta, loc, cha, t1r , t2r)
gr_r.detrend(type='linear')
gr_r.detrend(type='demean')
sr_r = gr_r[0].stats.sampling_rate
gr_r.filter(type='bandpass',freqmin=0.1, freqmax=(sr_r/2)-(sr_r/20))
gr_r.plot(type='relative',color='r',starttime=t1r, endtime=t2r)

#%% get day data

sta = 'FG13' # STATION 
cha = 'BDF' # CHANNEL
net = 'GI'  # 
loc = '03'    # location, it depends mostly of which network you are in. 

t1 = UTCDateTime(2019, 5, 21, 0, 0, 0) #the format is year:day_of_the_year:month
t2 = t1 + 60*60*24
st_day = Stream()
st_day = client.get_waveforms(net, sta, loc, cha, t1 , t2)
st_day.detrend(type='linear')
st_day.detrend(type='demean')
sr = st_day[0].stats.sampling_rate
st_day.filter(type='bandpass',freqmin=0.2, freqmax=(sr/2)-(sr/20))
st_day.plot(type='dayplot',starttime=t1, endtime=t2)

trace_t = st_day[0]
trace_i = st_day[0]


#%% Scan for Tremor
#window endpoints
start= t1 + 10 #time window start 
end= t2
trs = trace_t.slice(starttime = start  , endtime= end) 

trs.plot(type='relative',color='b')#, starttime=start , endtime=end)

nsta=int(5*sr)
nlta=int(60*sr)
stream=abs(trs.data)
cft=recursive_sta_lta(stream, nsta, nlta)
trig_on=5
trig_off=0.1
plot_trigger(trs, cft, trig_on, trig_off) 
on_off = trigger_onset(cft,trig_on,trig_off)


catalogue_t1 = np.zeros(shape=(0,3))
catalogue_t = np.zeros(shape=(0,3))
cat_T = np.zeros(shape=(0,3))
catalogue_a = np.zeros(shape=(0,3))
num=0
num1=0
numa = 0

shift2 = 1000
for x in range(0,len(on_off)):
    evlen = len(trace_t[on_off[x,0]:on_off[x,1]])/sr
    event = trace_t[on_off[x,0]-1000:on_off[x,1]+1000] 
    pre_event = trace_t[int(on_off[x,0]-(sr*10)):on_off[x,0]] 
    

    
    if ((on_off[x,1] - on_off[x,0])/sr) > 120:
        t_start = t1 + on_off[x,0]/sr  
        t_end = t1 + on_off[x,1]/sr 
        event1 = trs.slice(starttime = t_start , endtime= t_end )
        
        
        tr = event1
        sf = tr.stats.sampling_rate
        dom,cf, bwid50 = freq_info_NqV(tr,t_start,t_end,sf)
        med_dat = np.median(abs(event1.data))
        
#        print(dom, cf)
#        event1.plot(color='r')
        
        if dom > 1 and cf < 10:
        
#            print('')
#            print('central f =', cf, 'Hz')
#            print('dominant f =', dom, 'Hz')
#            print('bandwith 50% =', bwid50, 'Hz')
#            print('median amp =', med_dat)
#            print(event1.stats.starttime,'to', event1.stats.endtime )
            event1.plot(color='r')
            
            catalogue_t1 = np.lib.pad(catalogue_t1, ((0,1),(0,0)), 'constant', constant_values=(0))
            catalogue_t1[num1][0]= start + on_off[x,0]/sr #start time
            catalogue_t1[num1][1]= start + on_off[x,0]/sr  + evlen # end time
            catalogue_t1[num1][2]= evlen #duration
            num1 +=1


if len(catalogue_t1) > 0:
    catalogue_t = np.lib.pad(catalogue_t, ((0,1),(0,0)), 'constant', constant_values=(0))
    catalogue_t[0][0]= catalogue_t1[0][0]
    catalogue_t[0][1]= catalogue_t1[0][1]
    catalogue_t[0][2]= catalogue_t1[0][2]


    for x in range(1,len(catalogue_t1)):
        catalogue_t = np.lib.pad(catalogue_t, ((0,1),(0,0)), 'constant', constant_values=(0))
        catalogue_t[x][1]= catalogue_t1[x][1]
        
        t_gap = catalogue_t1[x,0] - catalogue_t1[x-1,1]
    #    T_len = catalogue_t1[x-1,2] + catalogue_t1[x,2]
        
        if  t_gap > 120 : 
            catalogue_t[x][0]= catalogue_t1[x][0]
            catalogue_t[x][2]= catalogue_t[x][1] - catalogue_t[x][0] 
        else:
            print(x)
            catalogue_t[x][0]= catalogue_t[x-1][0]
            catalogue_t[x][2]= catalogue_t[x][1] - catalogue_t[x][0]  
        
    
    for x in range(0,len(catalogue_t)-1):
        if catalogue_t[x,0] != catalogue_t[x+1,0]:
            catalogue_a = np.lib.pad(catalogue_a, ((0,1),(0,0)), 'constant', constant_values=(0))
            catalogue_a[numa][0]= catalogue_t1[x][0]
            catalogue_a[numa][1]= catalogue_t1[x][2]
            catalogue_a[numa][2]= 0
            
            cat_T = np.lib.pad(cat_T, ((0,1),(0,0)), 'constant', constant_values=(0))
            cat_T[numa][0]= catalogue_t1[x][0]
            cat_T[numa][1]= catalogue_t1[x][1]
            cat_T[numa][2]= catalogue_t1[x][2]
            
            numa +=1

    catalogue_a = np.lib.pad(catalogue_a, ((0,1),(0,0)), 'constant', constant_values=(0))
    catalogue_a[numa][0]= catalogue_t1[-1][0]
    catalogue_a[numa][1]= catalogue_t1[-1][2]
    catalogue_a[numa][2]= 0
    
    cat_T = np.lib.pad(cat_T, ((0,1),(0,0)), 'constant', constant_values=(0))
    cat_T[numa][0]= catalogue_t1[-1][0]
    cat_T[numa][1]= catalogue_t1[-1][1]
    cat_T[numa][2]= catalogue_t1[-1][2]       

else:
    print("No Tremor Detected")
        

#%% Scan for Isolated events
#window endpoints
start= t1 + 10 #time window start 
end= t2
trs = trace_i.slice(starttime = start  , endtime= end) 

#trs.plot(type='relative',color='b')#, starttime=start , endtime=end)

nsta=int(1*sr)
nlta=int(20*sr)
stream=trs.data
cft=recursive_sta_lta(stream, nsta, nlta)
trig_on=10
trig_off=0.1
plot_trigger(trs, cft, trig_on, trig_off) 
on_off = trigger_onset(cft,trig_on,trig_off)




#%%

catalogue_i = np.zeros(shape=(0,3))
num_i=0

shift = 200
if len(on_off) > 0:
    for x in range(0,len(on_off)):
        elen = len(trace_i[on_off[x,0]:on_off[x,1]])/sr
        event = trace_i[on_off[x,0]-1000:on_off[x,1]+2000]    
        event40 = trace_i[on_off[x,0]:on_off[x,0]+int(20*sr)]   
        
                    
        if elen/60 < 2 :
            
            e_start = t1 + on_off[x,0]/sr  
            e_end = t1 + on_off[x,1]/sr 
            event2 = trs.slice(starttime = e_start - 5 , endtime= e_end + 25 )        
            
    #        print(UTCDateTime(start + on_off[x,0]/sr ),'iso')
           
            et = int((start + on_off[x,0]/sr).timestamp)
            if len(catalogue_t1) > 0:
                near,ind=find_nearest(catalogue_a[:,0], et )
            else:
                near = 0
            
            if abs(et-near) > 30:   
                top_gr,top,corellgr = corel(gr_r[0],event40,shift)
                            
                catalogue_i = np.lib.pad(catalogue_i, ((0,1),(0,0)), 'constant', constant_values=(0))
                catalogue_i[num_i][0]= start + on_off[x,0]/sr 
                catalogue_i[num_i][1]= elen
                
                if top_gr > 0.5:
                    catalogue_i[num_i][2]= 1
                    
                else:
                    catalogue_i[num_i][2]= 2
                    
                num_i += 1
                
                catalogue_a = np.lib.pad(catalogue_a, ((0,1),(0,0)), 'constant', constant_values=(0))
                catalogue_a[numa][0]= start + on_off[x,0]/sr 
                catalogue_a[numa][1]= elen
                
                if top_gr > 0.5:
                    catalogue_a[numa][2]= 1
                    event2.plot(color='b')
                else:
                    catalogue_a[numa][2]= 2
                    event2.plot(color='g')
                numa+=1
       
#%% sort combined catalogue in terms of time

if len(catalogue_a) > 0:
    catalogue = catalogue_a[np.argsort(catalogue_a[0:-1,0])]


#%% label events

# 0 = Tremor
# 1 = Isolated Gas rich event
# 2 = Isolated Ash rich event
# 3 = Large Burst During Tremor

    event_num = len(catalogue)
    catalogue_d = np.zeros(shape=(event_num,3))
    catalogue_d[:,0]=catalogue[:,0]
    catalogue_d[:,1]=catalogue[:,2]
    catalogue_d[:,2]=catalogue[:,1]
    
    
    
    for x in range(0,len(catalogue_d)):
        
        if catalogue_d[x,1] == 0:
            trem_end = catalogue_d[x,0] + catalogue_d[x,2]
            
    #        print(UTCDateTime(catalogue_d[x,0]), UTCDateTime(trem_end))
            
            for y in range(x+1,len(catalogue)):
                if catalogue_d[y,0] <  trem_end :
                    if catalogue_d[y,1] > 0:
    #                    print(UTCDateTime(catalogue_d[y,0]))
                        catalogue_d[y,1] = 3
                
    #%%
if len(catalogue_d)>0:       
    for x in range(0,len(catalogue_d)):
        
        print(UTCDateTime(catalogue_d[x,0]),int(catalogue_d[x,1]),catalogue_d[x,2])
   
#%%  
    
#    Test on other days
#    Check vs swarm
#    harmonic tremor categorise
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    








