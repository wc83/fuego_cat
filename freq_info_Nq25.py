#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  6 16:03:40 2019

@author: william
"""

def freq_info_Nq25(tr,start,end):
    from numpy import argmax
    import numpy as np
        
    #frequency info
   
    window=end-start
    tr_data=tr.data
    m=np.mean(tr_data)
    tr_data = tr_data-m
    famp = abs(np.fft.fft(tr_data))
    
    
    # dominant f
    peak= argmax(abs(famp))/window
    if peak > 25:
        peak = 50-peak        
    
    # center f 
    hal =sum(famp)/4 #half of the first half of the complete fft (ie. less than nq)
    num=0
    misf_i=1000000000000
    for t in range(0,len(famp),1):
        n=famp[t]                  
        num += n
        misf = abs(num-hal)
        if misf < misf_i and t != 0:
            misf_i = misf
            cf=t/window
    
    # bandwidth 50%
    values = famp
    searchval = max(famp)/2
    ii = np.where(values >= searchval)[0]
    
    above=[]
    for p in range(0,len(ii)):
        if ii[p]/window < 25:
            above.append(ii[p]/window)        
    bwid50=max(above)-min(above)    
    
        # bandwidth 25%
    values2 = famp
    searchval2 = max(famp)/4
    iii = np.where(values2 >= searchval2)[0]
    
    above2=[]
    for x in range(0,len(iii)):
        if iii[x]/window < 25:
            above2.append(iii[x]/window)
    bwid25=max(above2)-min(above2)
        
    
    return(peak,cf, bwid50, bwid25)

























