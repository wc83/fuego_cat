#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 15:16:55 2019

@author: root
"""

def get_activity_fuego(day):

    from obspy.clients.earthworm import Client
    from obspy import UTCDateTime
    from obspy import Stream
    
    year1=2018
    month1=3
    day1=14
    hour1=0
    minute1=0
    second1=0
    
    num=0
    nums=0
    numa=0
    
    t0 = UTCDateTime(year1, month1, day1, hour1, minute1, second1) #the format is year:day_of_the_year:month
    t1 = t0 + day*24*60*60
    t2 = t1 + 23*60*60 + 59*60 +59.999 #UTCDateTime(year2, month2, day2, hour2, minute2, second2) # notice we have here 10 minutes, but we can select our times. 
        
    print(UTCDateTime(t1))
    
    
    port = '138.253.113.19' # ip, port - ip's 138.253.113.19 or 138.253.112.23
#%% FG10 seismic    
    try:  
        
        sta = 'FG10' # STATION LB01
        cha = 'BHZ' # CHANNEL - Vertical
        net = 'GI'  # Santiaguito volcano
        loc = ''    # location, it depends mostly of which network you are in. 

        client = Client(port, 16022) 
        st1 = Stream()
        st1 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        
        # st is a stream, we can operate normally as in obspy
        st1.detrend(type='linear')
        st1.detrend(type='demean')
        break_test=st1
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        
        FG10s=1
        num += 1
        nums += 1
                  
    except: 
        
        FG10s=0
        
#%% FG11 Seismic
    try:  
        
        sta = 'FG11' # STATION LB01
        cha = 'BHZ' # CHANNEL - Vertical
        net = 'GI'  # Santiaguito volcano
        loc = ''    # location, it depends mostly of which network you are in. 

        client = Client(port, 16022) 
        st1 = Stream()
        st1 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        
        # st is a stream, we can operate normally as in obspy
        st1.detrend(type='linear')
        st1.detrend(type='demean')
        break_test=st1
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
       
        FG11s=1
        num += 1
        nums += 1
            
    except:  

        FG11s=0

#%% FG12 Seismic
    try:  
        
        sta = 'FG12' # STATION LB01
        cha = 'BHZ' # CHANNEL - Vertical
        net = 'GI'  # Santiaguito volcano
        loc = '00'    # location, it depends mostly of which network you are in. 

        client = Client(port, 16022) 
        st1 = Stream()
        st1 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        
        # st is a stream, we can operate normally as in obspy
        st1.detrend(type='linear')
        st1.detrend(type='demean')
        break_test=st1
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        
        FG12s=1
        num += 1
        nums += 1
            
    except: 
        
        FG12s=0
        
#%% FG13 Seismic
    try:  
        
        sta = 'FG13' # STATION LB01
        cha = 'BHZ' # CHANNEL - Vertical
        net = 'GI'  # Santiaguito volcano
        loc = '00'    # location, it depends mostly of which network you are in. 

        client = Client(port, 16022) 
        st1 = Stream()
        st1 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        
        # st is a stream, we can operate normally as in obspy
        st1.detrend(type='linear')
        st1.detrend(type='demean')
        break_test=st1
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        
        FG13s=1
        num += 1
        nums += 1
            
    except:
        

        FG13s=0
        
#%% FG14 Seismic
    try:  
        
        sta = 'FG14' # STATION LB01
        cha = 'BHZ' # CHANNEL - Vertical
        net = 'GI'  # Santiaguito volcano
        loc = '01'    # location, it depends mostly of which network you are in. 

        client = Client(port, 16022) 
        st1 = Stream()
        st1 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        
        # st is a stream, we can operate normally as in obspy
        st1.detrend(type='linear')
        st1.detrend(type='demean')
        break_test=st1
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        
        FG14s=1
        num += 1
        nums += 1
            
    except: 
        

        FG14s=0
        

#%% FG16 Seismic
    try:  
        
        sta = 'FG16' # STATION LB01
        cha = 'BHZ' # CHANNEL - Vertical
        net = 'GI'  # Santiaguito volcano
        loc = '00'    # location, it depends mostly of which network you are in. 

        client = Client(port, 16022) 
        st1 = Stream()
        st1 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        
        # st is a stream, we can operate normally as in obspy
        st1.detrend(type='linear')
        st1.detrend(type='demean')
        break_test=st1
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        
        FG16s=1
        num += 1
        nums += 1
          
    except: 
        

        FG16s=0


#%% FG8 Seismic
    try:  
        
        sta = 'FG8' # STATION LB01
        cha = 'BHZ' # CHANNEL - Vertical
        net = 'GI'  # Santiaguito volcano
        loc = '00'    # location, it depends mostly of which network you are in. 

        client = Client(port, 16022) 
        st1 = Stream()
        st1 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        
        # st is a stream, we can operate normally as in obspy
        st1.detrend(type='linear')
        st1.detrend(type='demean')
        break_test=st1
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        
        FG8s=1
        num += 1
        nums += 1
                   
    except:     

        FG8s=0
        
#%% FG3 Seismic
    try:  
        
        sta = 'FG3' # STATION LB01
        cha = 'SHZ' # CHANNEL - Vertical
        net = 'GI'  # Santiaguito volcano
        loc = '01'    # location, it depends mostly of which network you are in. 

        client = Client(port, 16022) 
        st1 = Stream()
        st1 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        
        # st is a stream, we can operate normally as in obspy
        st1.detrend(type='linear')
        st1.detrend(type='demean')
        break_test=st1
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        
        FG3s=1
        num += 1
        nums += 1

    except: 

        FG3s=0



#%%
        




#%% FG10 inf  
    try:  
        
        sta = 'FG10' # STATION LB01
        cha = 'BDF' # CHANNEL - Vertical
        net = 'GI'  # Santiaguito volcano
        loc = '01'    # location, it depends mostly of which network you are in. 

        client = Client(port, 16022) 
        st1 = Stream()
        st1 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        
        # st is a stream, we can operate normally as in obspy
        st1.detrend(type='linear')
        st1.detrend(type='demean')
        break_test=st1
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        
        FG10a=1
        num += 1
        numa += 1
                  
    except: 
        
        FG10a=0
        
#%% FG11 inf
    try:  
        
        sta = 'FG11' # STATION LB01
        cha = 'BDF' # CHANNEL - Vertical
        net = 'GI'  # Santiaguito volcano
        loc = '01'    # location, it depends mostly of which network you are in. 

        client = Client(port, 16022) 
        st1 = Stream()
        st1 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        
        # st is a stream, we can operate normally as in obspy
        st1.detrend(type='linear')
        st1.detrend(type='demean')
        break_test=st1
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
       
        FG11a=1
        num += 1
        numa += 1
            
    except:  

        FG11a=0

#%% FG12 inf
    try:  
        
        sta = 'FG12' # STATION LB01
        cha = 'BDF' # CHANNEL - Vertical
        net = 'GI'  # Santiaguito volcano
        loc = '01'    # location, it depends mostly of which network you are in. 

        client = Client(port, 16022) 
        st1 = Stream()
        st1 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        
        # st is a stream, we can operate normally as in obspy
        st1.detrend(type='linear')
        st1.detrend(type='demean')
        break_test=st1
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        
        FG12a=1
        num += 1
        numa += 1
            
    except: 
        
        FG12a=0
        
#%% FG13 inf
    try:  
        
        sta = 'FG13' # STATION LB01
        cha = 'BDF' # CHANNEL - Vertical
        net = 'GI'  # Santiaguito volcano
        loc = '01'    # location, it depends mostly of which network you are in. 

        client = Client(port, 16022) 
        st1 = Stream()
        st1 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        
        # st is a stream, we can operate normally as in obspy
        st1.detrend(type='linear')
        st1.detrend(type='demean')
        break_test=st1
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        
        FG13a=1
        num += 1
        numa += 1
            
    except:
        
        FG13a=0
        

#%% FG15 inf
    try:  
        
        sta = 'FG15' # STATION LB01
        cha = 'BDF' # CHANNEL - Vertical
        net = 'GI'  # Santiaguito volcano
        loc = '01'    # location, it depends mostly of which network you are in. 

        client = Client(port, 16022) 
        st1 = Stream()
        st1 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        
        # st is a stream, we can operate normally as in obspy
        st1.detrend(type='linear')
        st1.detrend(type='demean')
        break_test=st1
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        
        FG15a=1
        num += 1
        numa += 1
          
    except: 
        
        FG15a=0



#%% FG8 inf
    try:  
        
        sta = 'FG8' # STATION LB01
        cha = 'BDF' # CHANNEL - Vertical
        net = 'GI'  # Santiaguito volcano
        loc = '01'    # location, it depends mostly of which network you are in. 

        client = Client(port, 16022) 
        st1 = Stream()
        st1 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        
        # st is a stream, we can operate normally as in obspy
        st1.detrend(type='linear')
        st1.detrend(type='demean')
        break_test=st1
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        
        FG8a=1
        num += 1
        numa += 1
          
    except: 
        
        FG8a=0


#%% FV01 inf
    try:  
        
        sta = 'FV01' # STATION LB01
        cha = 'HDF' # CHANNEL - Vertical
        net = 'GI'  # Santiaguito volcano
        loc = ''    # location, it depends mostly of which network you are in. 

        client = Client(port, 16022) 
        st1 = Stream()
        st1 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        
        # st is a stream, we can operate normally as in obspy
        st1.detrend(type='linear')
        st1.detrend(type='demean')
        break_test=st1
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        
        FV01a=1
        num += 1
        numa += 1
          
    except: 
        
        FV01a=0

#%% FV02 inf
    try:  
        
        sta = 'FV02' # STATION LB01
        cha = 'HDF' # CHANNEL - Vertical
        net = 'GI'  # Santiaguito volcano
        loc = ''    # location, it depends mostly of which network you are in. 

        client = Client(port, 16022) 
        st1 = Stream()
        st1 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        
        # st is a stream, we can operate normally as in obspy
        st1.detrend(type='linear')
        st1.detrend(type='demean')
        break_test=st1
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        
        FV02a=1
        num += 1
        numa += 1
          
    except: 
        
        FV02a=0
        
#%% FV03 inf
    try:  
        
        sta = 'FV03' # STATION LB01
        cha = 'HDF' # CHANNEL - Vertical
        net = 'GI'  # Santiaguito volcano
        loc = ''    # location, it depends mostly of which network you are in. 

        client = Client(port, 16022) 
        st1 = Stream()
        st1 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        
        # st is a stream, we can operate normally as in obspy
        st1.detrend(type='linear')
        st1.detrend(type='demean')
        break_test=st1
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        
        FV03a=1
        num += 1
        numa += 1
          
    except: 
        
        FV03a=0

#%% FV04 inf
    try:  
        
        sta = 'FV04' # STATION LB01
        cha = 'HDF' # CHANNEL - Vertical
        net = 'GI'  # Santiaguito volcano
        loc = ''    # location, it depends mostly of which network you are in. 

        client = Client(port, 16022) 
        st1 = Stream()
        st1 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        
        # st is a stream, we can operate normally as in obspy
        st1.detrend(type='linear')
        st1.detrend(type='demean')
        break_test=st1
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        
        FV04a=1
        num += 1
        numa += 1
          
    except: 
        
        FV04a=0



#%% VF01 inf
    try:  
        
        sta = 'VF01' # STATION LB01
        cha = 'HDF' # CHANNEL - Vertical
        net = 'XZ'  # Santiaguito volcano
        loc = ''    # location, it depends mostly of which network you are in. 

        client = Client(port, 16022) 
        st1 = Stream()
        st1 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        
        # st is a stream, we can operate normally as in obspy
        st1.detrend(type='linear')
        st1.detrend(type='demean')
        break_test=st1
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        
        VF01a=1
        num += 1
        numa += 1
          
    except: 
        
        VF01a=0


#%% VF02 inf
    try:  
        
        sta = 'VF02' # STATION LB01
        cha = 'HDF' # CHANNEL - Vertical
        net = 'XZ'  # Santiaguito volcano
        loc = ''    # location, it depends mostly of which network you are in. 

        client = Client(port, 16022) 
        st1 = Stream()
        st1 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        
        # st is a stream, we can operate normally as in obspy
        st1.detrend(type='linear')
        st1.detrend(type='demean')
        break_test=st1
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        
        VF02a=1
        num += 1
        numa += 1
          
    except: 
        
        VF02a=0
        
        
#%% VF03 inf
    try:  
        
        sta = 'VF03' # STATION LB01
        cha = 'HDF' # CHANNEL - Vertical
        net = 'XZ'  # Santiaguito volcano
        loc = ''    # location, it depends mostly of which network you are in. 

        client = Client(port, 16022) 
        st1 = Stream()
        st1 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        
        # st is a stream, we can operate normally as in obspy
        st1.detrend(type='linear')
        st1.detrend(type='demean')
        break_test=st1
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        
        VF03a=1
        num += 1
        numa += 1
          
    except: 
        
        VF03a=0

#%% VF04 inf
    try:  
        
        sta = 'VF04' # STATION LB01
        cha = 'HDF' # CHANNEL - Vertical
        net = 'XZ'  # Santiaguito volcano
        loc = ''    # location, it depends mostly of which network you are in. 

        client = Client(port, 16022) 
        st1 = Stream()
        st1 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        
        # st is a stream, we can operate normally as in obspy
        st1.detrend(type='linear')
        st1.detrend(type='demean')
        break_test=st1
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        
        VF04a=1
        num += 1
        numa += 1
          
    except: 
        
        VF04a=0

#%% VF05 inf
    try:  
        
        sta = 'VF05' # STATION LB01
        cha = 'HDF' # CHANNEL - Vertical
        net = 'XZ'  # Santiaguito volcano
        loc = ''    # location, it depends mostly of which network you are in. 

        client = Client(port, 16022) 
        st1 = Stream()
        st1 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        
        # st is a stream, we can operate normally as in obspy
        st1.detrend(type='linear')
        st1.detrend(type='demean')
        break_test=st1
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        
        VF05a=1
        num += 1
        numa += 1
          
    except: 
        
        VF05a=0

#%% VF06 inf
    try:  
        
        sta = 'VF06' # STATION LB01
        cha = 'HDF' # CHANNEL - Vertical
        net = 'XZ'  # Santiaguito volcano
        loc = ''    # location, it depends mostly of which network you are in. 

        client = Client(port, 16022) 
        st1 = Stream()
        st1 = client.get_waveforms(net, sta, loc, cha, t1 , t2)
        
        # st is a stream, we can operate normally as in obspy
        st1.detrend(type='linear')
        st1.detrend(type='demean')
        break_test=st1
        break_test = break_test[0].filter("bandpass", freqmin=1,freqmax=10)
        
        VF06a=1
        num += 1
        numa += 1
          
    except: 
        
        VF06a=0








    #%% return all stations
    
    return(num,nums,numa,FG3s,FG8s,FG10s,FG11s,FG12s,FG13s,FG14s,FG16s,FG8a,FG10a,FG11a,FG12a,FG13a,FG15a,FV01a,FV02a,FV03a,FV04a,VF01a,VF02a,VF03a,VF04a,VF05a,VF06a) 