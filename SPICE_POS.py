#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 28 16:57:00 2024

@author: cemaxim
"""
import pandas as pd
import glob
import spiceypy as spice
#from CAS_ORBITS import ORBIT


def SPICE_POS(df,time):
    spice.kclear()
    mk = glob.glob('/home/cemaxim/cosp_1000/extras/mk/*(copy).tm')
    spice.furnsh(mk)
    step = 4000
    
    etOne = spice.str2et(time[0])
    etTwo = spice.str2et(time[1])
    
    times = [x*(etTwo-etOne)/step + etOne for x in range(step)]
    positions, lightimes = spice.spkpos('CASSINI', times, 'CASSINI_KSM', 'NONE', 'SATURN')
    positions = positions.T/60268
    
    dfp = pd.DataFrame(positions.T, spice.et2datetime(times))
    dfp.columns = ['x', 'y', 'z']
    dfp.index.name = 'datetime'
    dfp.index = dfp.index.round('s')
    dfp.index = dfp.index.tz_convert(None)
    
    df = pd.merge_asof(df, dfp, left_index=True, right_index=True, direction='backward')
    return(df)


mag_df = pd.read_pickle('/home/cemaxim/pickles/MAG_DF.pkl')
t1 = '2006-09-25 11:58:31'
t2 = '2006-11-20 09:31:31'
time = [t1[0], t2[1]]
