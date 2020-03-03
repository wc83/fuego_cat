#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 08:15:32 2020

@author: root
"""

def get_all_Fuego_stations(timestamp):
    
    from obspy.clients.earthworm import Client
    from obspy import UTCDateTime
    from obspy import Stream
    import numpy as np

    
    station_activity = np.zeros(shape=(0,3))
    num = 0
    
    t1 = UTCDateTime(timestamp) #the format is year:day_of_the_year:month
    t2 = t1 + 24*60*60
#%%
    sta = 'VF01' # STATION VF01
    cha = 'HDF' # CHANNEL - inf
    net = 'XZ'  # Fuego volcano
    loc = ''    # location
    
    client = Client('138.253.112.23', 16022)
    
    t1_2s = UTCDateTime(1526774400)
    t2_2s = t1_2s + 2
    
    st_blank = Stream()
    st_blank = client.get_waveforms(net, sta, loc, cha, t1_2s , t2_2s)
    st_blank.detrend(type='linear')
    st_blank.detrend(type='demean')
    st_blank.filter(type='bandpass',freqmin=0.5, freqmax=5)
    
    
#%% VF01    
    try:  
        ID = 1
        sta = 'VF01' # STATION VF01
        cha = 'HDF' # CHANNEL - inf
        net = 'XZ'  # Fuego volcano
        loc = ''    # location, 
            
        st1 = Stream()
        st1 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        st1.detrend(type='linear')
        st1.detrend(type='demean')
        
        break_test=st1
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        sr = st1[0].stats.sampling_rate
        
        if sum(abs(st1[0].data)) > 10 and st1[0].stats.npts > 7200*sr:
            station_activity = np.lib.pad(station_activity, ((0,1),(0,0)), 'constant', constant_values=(0))
    
            station_activity[num][0] = ID
            station_activity[num][1] = st1[0].stats.starttime
            station_activity[num][2] = st1[0].stats.endtime
            num+=1
        
        else:
            st1 = st_blank
            
    except:
        st1 = st_blank
            
            
        
#%% VF02   
    try:  
        ID = 2
        sta = 'VF02' # STATION VF01
        cha = 'HDF' # CHANNEL - inf
        net = 'XZ'  # Fuego volcano
        loc = ''    # location, 
    
        st2 = Stream()
        st2 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        st2.detrend(type='linear')
        st2.detrend(type='demean')
        
        break_test=st2
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        sr = st2[0].stats.sampling_rate
        
        if sum(abs(st2[0].data)) > 10 and st2[0].stats.npts > 7200*sr:
            station_activity = np.lib.pad(station_activity, ((0,1),(0,0)), 'constant', constant_values=(0))
            station_activity[num][0] = ID
            station_activity[num][1] = st2[0].stats.starttime
            station_activity[num][2] = st2[0].stats.endtime
            num+=1
            
        
        else:
            st2 = st_blank
            
    except:
        st2 = st_blank        
        
        
#%% VF03   
    try:  
        ID = 3
        sta = 'VF03' # STATION VF01
        cha = 'HDF' # CHANNEL - inf
        net = 'XZ'  # Fuego volcano
        loc = ''    # location, 
        
        st3 = Stream()
        st3 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        st3.detrend(type='linear')
        st3.detrend(type='demean')
        
        break_test=st3
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        sr = st3[0].stats.sampling_rate
        
        if sum(abs(st3[0].data)) > 10 and st3[0].stats.npts > 7200*sr:
            station_activity = np.lib.pad(station_activity, ((0,1),(0,0)), 'constant', constant_values=(0))
            station_activity[num][0] = ID
            station_activity[num][1] = st3[0].stats.starttime
            station_activity[num][2] = st3[0].stats.endtime
            num+=1
            
        
        else:
            st3 = st_blank
            
    except:
        st3 = st_blank        
        
        
        
 #%% VF04  
    try:  
        ID = 4
        sta = 'VF04' # STATION VF01
        cha = 'HDF' # CHANNEL - inf
        net = 'XZ'  # Fuego volcano
        loc = ''    # location, 

        st4 = Stream()
        st4 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        st4.detrend(type='linear')
        st4.detrend(type='demean')
        
        break_test=st1
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        sr = st4[0].stats.sampling_rate
        
        if sum(abs(st4[0].data)) > 10 and st4[0].stats.npts > 7200*sr:
            station_activity = np.lib.pad(station_activity, ((0,1),(0,0)), 'constant', constant_values=(0))
            station_activity[num][0] = ID
            station_activity[num][1] = st4[0].stats.starttime
            station_activity[num][2] = st4[0].stats.endtime
            num+=1
            
        
        else:
            st4 = st_blank
            
    except:
        st4 = st_blank        
 
               
 #%% VF05   
    try:  
        ID = 5
        sta = 'VF05' # STATION VF01
        cha = 'HDF' # CHANNEL - inf
        net = 'XZ'  # Fuego volcano
        loc = ''    # location, 

        st5 = Stream()
        st5 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        st5.detrend(type='linear')
        st5.detrend(type='demean')
        
        break_test=st5
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        sr = st5[0].stats.sampling_rate
        
        
        if sum(abs(st5[0].data)) > 10 and st5[0].stats.npts > 7200*sr:
            station_activity = np.lib.pad(station_activity, ((0,1),(0,0)), 'constant', constant_values=(0))
            station_activity[num][0] = ID
            station_activity[num][1] = st5[0].stats.starttime
            station_activity[num][2] = st5[0].stats.endtime
            num+=1
            
        
        else:
            st5 = st_blank
            
    except:
        st5 = st_blank        
        
        
 #%% VF06  
    try:  
        ID = 6
        sta = 'VF06' # STATION VF01
        cha = 'HDF' # CHANNEL - inf
        net = 'XZ'  # Fuego volcano
        loc = ''    # location, 

        st6 = Stream()
        st6 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        st6.detrend(type='linear')
        st6.detrend(type='demean')
        
        break_test=st6
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        sr = st6[0].stats.sampling_rate
        
        
        if sum(abs(st6[0].data)) > 10 and st6[0].stats.npts > 7200*sr:
            station_activity = np.lib.pad(station_activity, ((0,1),(0,0)), 'constant', constant_values=(0))
            station_activity[num][0] = ID
            station_activity[num][1] = st6[0].stats.starttime
            station_activity[num][2] = st6[0].stats.endtime
            num+=1
            
        
        else:
            st6 = st_blank
            
    except:
        st6 = st_blank        
        
        
  




      
                       
  #%%      give all traces and list of active stations
        
        
    return(st1,st2,st3,st4,st5,st6,station_activity)    
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        