#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 09:55:05 2020

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



#%%

Day_start = UTCDateTime(2018,5,20,0,0,0).timestamp

st1,st2,st3,st4,st5,st6,station_activity = get_all_Fuego_stations(Day_start)

#%% loop over stations

print(len(station_activity),"stations running")
if len(station_activity) > 0:
    for i in range(0,len(station_activity)):
        
        st_ID = int(station_activity[i,0])
        
        if st_ID == 1:
            st_day = st1
        if st_ID == 2:
            st_day = st2
        if st_ID == 3:
            st_day = st3
        if st_ID == 4:
            st_day = st4
        if st_ID == 5:
            st_day = st5
        if st_ID == 6:
            st_day = st6

        t1 = UTCDateTime(station_activity[i,1])
        t2 = UTCDateTime(station_activity[i,2])


        
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
#                    event1.plot(color='r')
                    
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
            print("No Tremor Detected by station", i)
                
        
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
#                            event2.plot(color='b')
                        else:
                            catalogue_a[numa][2]= 2
#                            event2.plot(color='g')
                        numa+=1
        
#%%
        if len(catalogue_a) > 0:
            catalogue_s = catalogue_a[np.argsort(catalogue_a[0:-1,0])]
        
        
        #%% label events
        
        # 0 = Tremor
        # 1 = Isolated Gas rich event
        # 2 = Isolated Ash rich event
        # 3 = Large Burst During Tremor
        
            event_num = len(catalogue_s)
            catalogue_d = np.zeros(shape=(event_num,3))
            catalogue_d[:,0]=catalogue_s[:,0]
            catalogue_d[:,1]=catalogue_s[:,2]
            catalogue_d[:,2]=catalogue_s[:,1]
            
            for x in range(0,len(catalogue_d)):
                if catalogue_d[x,1] == 0:
                    trem_end = catalogue_d[x,0] + catalogue_d[x,2]
                    
                    for y in range(x+1,len(catalogue_s)):
                        if catalogue_d[y,0] <  trem_end :
                            if catalogue_d[y,1] > 0:
                                catalogue_d[y,1] = 3
        
        if station_activity[i,0] == 1:
            catalogue_VF01 = catalogue_d
        if station_activity[i,0] == 2:
            catalogue_VF02 = catalogue_d
        if station_activity[i,0] == 3:
            catalogue_VF03 = catalogue_d
        if station_activity[i,0] == 4:
            catalogue_VF04 = catalogue_d
        if station_activity[i,0] == 5:
            catalogue_VF05 = catalogue_d
        if station_activity[i,0] == 6:
            catalogue_VF06 = catalogue_d
            
        #Combined station events    
        if i == 0:
            li = len(catalogue_i)
            lt = len(cat_T)
            
            cat_i = np.zeros(shape=(li,9))
            cat_t = np.zeros(shape=(lt,4))

            cat_i[:,0] = catalogue_i[:,0]
            cat_i[:,8] = catalogue_i[:,1]
            cat_i[:,1] = catalogue_i[:,2]
            
            stat_ind = st_ID +1
            cat_i[:,stat_ind]= 1
            print("STATION INDEX i=0 is:",stat_ind)
            
            
            
            cat_t[:,0] = cat_T[:,0]
            cat_t[:,1] = cat_T[:,1]
            cat_t[:,2] = cat_T[:,2]
    
            cat_t[:][3]= st_ID
            
        else:
             # add new isolated explosions if detected on new stations
             for e in range(0,len(catalogue_i)):
                ev = catalogue_i[e,0]
                near,ind=find_nearest(cat_i[:,0], ev )
                if abs(ev-near) > 30: 
                    cat_i = np.lib.pad(cat_i, ((0,1),(0,0)), 'constant', constant_values=(0))
                    cat_i[-1][0]=  catalogue_i[e,0]
                    cat_i[-1][8]=  catalogue_i[e,1]
                    cat_i[-1][1]=  catalogue_i[e,2]
                    stat_ind = st_ID +1
                    cat_i[-1][stat_ind]=1
                # Add labels of detecting stations
                else:
                    stat_ind = st_ID +1
                    cat_i[ind][stat_ind]=1
                    
             # add all tremor from all stations - to be organised later
             for e in range(0,len(cat_T)):
                cat_t = np.lib.pad(cat_t, ((0,1),(0,0)), 'constant', constant_values=(0))
                cat_t[-1][0]=  cat_T[e,0]
                cat_t[-1][1]=  cat_T[e,1]
                cat_t[-1][2]=  cat_T[e,2]
                cat_t[-1][3]=  st_ID
                            

#%% order detections in time
                                
catalogue_day_iso = cat_i[np.argsort(cat_i[0:-1,0])]
catalogue_day_trem_list = cat_t[np.argsort(cat_t[0:-1,0])]

if catalogue_day_iso[0][0] == 1:
    catalogue_day_iso = np.delete(catalogue_day_iso, 0, 0)    
if catalogue_day_trem_list[0][0] == 1:
    catalogue_day_trem_list = np.delete(catalogue_day_trem_list, 0, 0)


#%% combine tremor detections from multiple stations
    
    
catalogue_day_trem = np.zeros(shape=(0,8))

if len(catalogue_day_trem_list)>0: 
    # make first entry 
    catalogue_day_trem = np.lib.pad(catalogue_day_trem, ((0,1),(0,0)), 'constant', constant_values=(0))
    catalogue_day_trem[0][0]= catalogue_day_trem_list[0][0]
    # add station tick 
    st_det = catalogue_day_trem_list[0][3]
    st_det_row = int(st_det )
    catalogue_day_trem[0][st_det_row] = 1
    # get list of durations for 1st entry
    duration_l = np.zeros(shape=(1,1))
    duration_l[0][0] = catalogue_day_trem_list[0][2]

    # loop over all events and save new events, or add station ticks
    for e in range(1,len(catalogue_day_trem_list)):
        #if stating within 30 seconds - add tick
        if catalogue_day_trem_list[e][0] - catalogue_day_trem[-1][0] < 30:
            st_det = catalogue_day_trem_list[e][3]
            st_det_row = int(st_det )
            catalogue_day_trem[-1][st_det_row] = 1
            #add to list of durations
            duration_l = np.lib.pad(duration_l, ((0,1),(0,0)), 'constant', constant_values=(0))
            duration_l[-1][0] = catalogue_day_trem_list[e][2]
        # get duration for previous event, save new event   
        else:
            if len(duration_l)== 2:
                med_dur = np.min(duration_l)
            else:
                med_dur = np.median(duration_l)
            catalogue_day_trem[-1][7] = med_dur
            
            #reset duration list
            duration_l = np.zeros(shape=(1,1))
            duration_l[0][0] = catalogue_day_trem_list[e][2]
            
            catalogue_day_trem = np.lib.pad(catalogue_day_trem, ((0,1),(0,0)), 'constant', constant_values=(0))
            catalogue_day_trem[-1][0]= catalogue_day_trem_list[e][0]
            
            st_det = catalogue_day_trem_list[e][3]
            st_det_row = int(st_det )
            catalogue_day_trem[-1][st_det_row] = 1
            
            
if len(duration_l)== 2:
    med_dur = np.min(duration_l)
else:
    med_dur = np.median(duration_l)
catalogue_day_trem[-1][7] = med_dur            

#only start times confirmed by 2 or more
catalogue_day_trem_co = np.zeros(shape=(0,8))

if len(catalogue_day_trem)>0: 
    for x in range(0,len(catalogue_day_trem)):
        if sum(catalogue_day_trem[x,1:7]) > 1:
            catalogue_day_trem_co = np.lib.pad(catalogue_day_trem_co, ((0,1),(0,0)), 'constant', constant_values=(0))
            catalogue_day_trem_co[-1][:]= catalogue_day_trem[x][:]    
# add stations which catch end half of an event
        if len(catalogue_day_trem_co) > 0:
            if catalogue_day_trem_co[-1][0] < catalogue_day_trem[x,0] < catalogue_day_trem_co[-1][0] + catalogue_day_trem_co[-1][7]:
                # add station ticks
                if catalogue_day_trem[x][1] == 1:
                    catalogue_day_trem_co[-1][1] = 1
                if catalogue_day_trem[x][2] == 1:
                    catalogue_day_trem_co[-1][2] = 1
                if catalogue_day_trem[x][3] == 1:
                    catalogue_day_trem_co[-1][3] = 1
                if catalogue_day_trem[x][4] == 1:
                    catalogue_day_trem_co[-1][4] = 1
                if catalogue_day_trem[x][5] == 1:
                    catalogue_day_trem_co[-1][5] = 1
                if catalogue_day_trem[x][6] == 1:
                    catalogue_day_trem_co[-1][6] = 1
      
# remove tremor which is only half the event
for x in range(0,len(catalogue_day_trem_co)-1): 
    try:
        for y in range(x+1,len(catalogue_day_trem_co)):
            if catalogue_day_trem_co[y][0] < catalogue_day_trem_co[x][0] + catalogue_day_trem_co[x][7]:
                catalogue_day_trem_co = np.delete(catalogue_day_trem_co, y, 0)    
        
    except:
        do_nothing=0
                                
#%% print out catalogues

if len(catalogue_day_iso)>0:       
    for x in range(0,len(catalogue_day_iso)):
        
        print(UTCDateTime(catalogue_day_iso[x,0]),catalogue_day_iso[x,1],catalogue_day_iso[x,2],catalogue_day_iso[x,3],catalogue_day_iso[x,4],catalogue_day_iso[x,5],catalogue_day_iso[x,6],catalogue_day_iso[x,7],catalogue_day_iso[x,8])
        
print("")                
if len(catalogue_day_trem_co)>0:       
    for x in range(0,len(catalogue_day_trem_co)):
        
        print(UTCDateTime(catalogue_day_trem_co[x,0]), 'to' ,UTCDateTime(catalogue_day_trem_co[x,7] + catalogue_day_trem_co[x,0]), catalogue_day_trem_co[x,1],catalogue_day_trem_co[x,2],catalogue_day_trem_co[x,3],catalogue_day_trem_co[x,4],catalogue_day_trem_co[x,5],catalogue_day_trem_co[x,6])
                
                






















