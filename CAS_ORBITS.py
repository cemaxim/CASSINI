#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 29 10:47:47 2024

@author: cemaxim
"""

from SPICE_POS import SPICE_POS
from scipy.signal import argrelextrema
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

"""
start_time = 'Jan 1, 2004'
end_time = 'Jan 1, 2010'
time = [start_time, end_time]

ksm = pd.read_pickle("/home/cemaxim/pickles/MAG_KSM.pkl")
ksm = ksm[['BX', 'BY', 'BZ']][time[0]:time[1]]

df = SPICE_POS(ksm,time)
df.insert(1,'r',np.sqrt(df.x**2 + df.y**2 + df.z**2))
df.r = df.r.drop_duplicates(keep='first')
df.r = df.r.interpolate()
df = df.iloc[::15,:]

n=600
df['perikrones'] = df.iloc[argrelextrema(df.r.values,np.less_equal,order=n)[0]]['r']
perikrones = df.iloc[argrelextrema(df.r.values,np.less_equal,order=n)[0]]['r']
perikrones = perikrones.to_frame()
perikrones.insert(1,'perikrone',np.arange(0,len(perikrones.index)))

perikrones.to_pickle('/home/cemaxim/pickles/PERIKRONES.pkl')
"""
def ORBIT(n):
    perikrones = pd.read_pickle('/home/cemaxim/pickles/PERIKRONES.pkl')
    t1 = perikrones.index[n].strftime('%Y-%m-%d %X')
    t2 = perikrones.index[n+1].strftime('%Y-%m-%d %X')
  
    #orb = df[['x','y']][t1:t2]
    #orb.plot(y=['x','y'], use_index=True)
    return [t1,t2]

"""
fig = plt.figure(figsize=(10,10))
df.plot(y=['r'], use_index=True)
plt.scatter(df.index, df['perikrones'], c='r')


"""