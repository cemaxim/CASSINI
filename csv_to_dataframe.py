#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 28 14:53:44 2024

@author: cemaxim
"""
import glob
import pandas as pd

files = glob.glob('/media/padelamere/DATAPART1/CAPS/MAG/COMAG_4002/DATA/*/*KSM_1M.TAB')

li = []

for filename in files:
    df = pd.read_csv(filename, index_col=None,header=None,delim_whitespace=True,parse_dates=True, usecols=[0,1,2,3])
    df[0] = pd.to_datetime(df[0])
    df.columns =['datetime','BX','BY','BZ']
    #df.set_index('datetime',inplace=True)
    li.append(df)
  
df = pd.concat(li, axis=0, ignore_index=True)
df.set_index('datetime',inplace=True)
df.to_pickle('/home/cemaxim/pickles/MAG_KSM1.pkl')

ksm1 = pd.read_pickle("/home/cemaxim/pickles/MAG_KSM1.pkl")
