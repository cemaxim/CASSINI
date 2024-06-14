#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 23 09:40:26 2024

@author: cemaxim
"""
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.markers import MarkerStyle
import matplotlib.dates as mdates
import numpy as np
from scipy import signal
from SPICE_POS import SPICE_POS
from CAS_ORBITS import ORBIT
import pycwt as wavelet
import pywt
#start_time = 'Feb 1, 2007'
#end_time = 'Mar 15, 2007'
#time = [start_time, end_time]

n=24 #orbit number 0-119
m=36
t1 = ORBIT(n)
t2 = ORBIT(m)
time = [t1[0], t1[1]]

start = t1[0]
end = t2[1]
#dti = pd.date_range(start, end, freq = 'min')

"GET DATA"
mag_df = pd.read_pickle("/home/cemaxim/pickles/MAG_DF.pkl")
mag_df = mag_df[time[0] : time[1]]
#mag_df = mag_df.resample('1T',origin=t1[0]).fillna(method=None)  #fill data for missing datetimes with NaN

bc_df = pd.read_pickle('/home/cemaxim/pickles/BC_ID.pkl')
df = SPICE_POS(mag_df,time)
df = pd.merge_asof(df, bc_df, left_index=True, right_index=True, direction='forward')
df.insert(0,'R',np.sqrt(df.x**2 + df.y**2 + df.z**2))


def plot_orbit(id, data):#, start_time, end_time ,n):
    #t = [start_time, end_time]
    #c_data = id[t[0] : t[1]]
    #data = data[t[0] : t[1]]
    fig, ax = plt.subplots(figsize=(9,9))  
    for i in range(len(id.index)-1):
        t = data.index
        t1 = id.index[i]
        t2 = id.index[i+1]
        wh = (t>t1) & (t<t2)
        x = data.iloc[wh]['x'].to_numpy()
        y = data.iloc[wh]['y'].to_numpy()
        if id.iloc[i]['ID'] == -1:
            clr = '-r'
            #label = 'solar wind'
        if id.iloc[i]['ID'] == 2:
            clr = 'lime'
            #label = 'sheath'
        if id.iloc[i]['ID'] == 1:
            clr = '-b'
            #label = 'magnetosphere'
        #if id.iloc[i]['ID'] == -2:
            #clr = '-c'
        ax.plot(x, y, clr)#, 'o', markersize=10, markevery=[0,-1])
        ax.set_xlabel('x ($R_S$)')
        ax.set_ylabel('y ($R_S$)')
        ax.plot(0,0,markersize = 15, color = 'black', markerfacecoloralt='white',marker=MarkerStyle("o", fillstyle='left'))
        plt.gca().set_aspect('equal')
        ax.grid(True, linestyle=':')
        ax.set_title('Cassini Orbit '+str(n))#Position from\n'+str(time[0])+' to '+str(time[1]))
    
def B_vs_R(data, R_min, win, start_time, end_time):
    
    data = data[start_time : end_time]
    #data = pd.merge_asof(data, id, left_index=True, right_index=True, direction='forward')
    data = data.loc[(data.R > R_min) & (data.ID ==1)]
    data = data.resample('1T',origin=t1[0]).fillna(method=None)
    #is_nan = data['BX'].isna() 
    
    fig, axs = plt.subplots(3,sharey=True)
    fig.set_size_inches((12,8))
    fig.supxlabel('Radial Distance ($R_S$)')
    axs[0].plot(data.R.values, data.BR.values)
    axs[0].plot(data.R.values, data.BR.rolling(win,center=True).mean(),'-r')
    axs[0].set_ylabel('$B_R$ (nT)')
    axs[1].plot(data.R.values, data.BTHETA.values)
    axs[1].plot(data.R.values, data.BTHETA.rolling(win,center=True).mean(),'-r')
    axs[1].set_ylabel('$B_{\theta}$ (nT)')
    axs[2].plot(data.R.values, data.BPHI.values)
    axs[2].plot(data.R.values, data.BPHI.rolling(win,center=True).mean(),'-r')
    axs[2].set_ylabel('$B_{\phi}$ (nT)')
"""
def filter(data, R_lim, ID):
    data = data.loc[(data.R > R_lim) & (data.ID ==1)]
    return data
"""

def MFA(data, win, R_min, start, end): #mean field align
    data = data.loc[(data.R > R_min) & (data.ID ==1)]
    means = pd.DataFrame({'BX': data.BX.rolling(win,center=True).mean(),
                          'BY': data.BY.rolling(win,center=True).mean(),
                          'BZ': data.BZ.rolling(win,center=True).mean()})
    means = means.dropna()
    means = means[:][start : end]
    means_mags = np.sqrt(means.BX**2+means.BY**2+means.BZ**2)
    raw = data[['BX', 'BY', 'BZ']]
    raw = raw.reindex_like(means,None)
    raw = raw[:][start : end]
    raw_mags = np.sqrt(data.BX**2+data.BY**2+data.BZ**2)
    li = []
    for t in range(len(means.index)):
        z_hat = means.iloc[t]/means_mags.iloc[t]
        R_hat = raw.iloc[t]/raw_mags.iloc[t]
        y_hat = np.cross(z_hat,R_hat)/np.sqrt(np.sum(np.cross(z_hat,R_hat)**2))
        x_hat = np.cross(y_hat,z_hat)
        bx = np.dot(raw.iloc[t].to_numpy()-means.iloc[t].to_numpy(),x_hat)
        by = np.dot(raw.iloc[t].to_numpy()-means.iloc[t].to_numpy(),y_hat)
        bz = np.dot(raw.iloc[t].to_numpy()-means.iloc[t].to_numpy(),z_hat)   
        #temp = pd.DataFrame({'B_PERP1':[bx], 'B_PERP2':[by], 'B_PAR':[bz]})#, index = means.index[t])
        #li.append(temp)
        li.append({'B_PERP1':bx, 'B_PERP2':by, 'B_PAR':bz})
    mfa_data = pd.DataFrame(li)
    mfa_data = mfa_data.set_index(means.index)
    data_R = data[['R']][mfa_data.index[0] : mfa_data.index[-1]]
    mfa_data = pd.merge_asof(data_R, mfa_data, left_index=True, right_index=True, direction='forward')
    

    fig, axs = plt.subplots(4, sharex=True)
    fig.set_size_inches((12,8))
    fig.tight_layout()
    fig.supxlabel('TIME')
    axs[0].plot(mfa_data.index.values, mfa_data.B_PERP1.values)
    axs[0].set_ylabel('B_PERP1 (nT)')
    #axs[0].grid()
    axs[1].plot(mfa_data.index.values, mfa_data.B_PERP2.values)
    axs[1].set_ylabel('B_PERP2 (nT)')
    #axs[1].grid()
    axs[2].plot(mfa_data.index.values, mfa_data.B_PAR.values)
    axs[2].set_ylabel('B_PAR (nT)')
    #axs[2].grid()
    axs[3].plot(mfa_data.index.values, mfa_data.R.values)
    axs[3].set_ylabel('Radial Distance ($R_S$)')
   # axs[3].get_shared_x_axes().join(axs[0], axs[1], axs[2], axs[3])
    return mfa_data


def CWT(data, freq):
    data = data.asfreq(freq)
    #data = data.resample(n, origin = 'start').asfreq()
    #data = data.iloc[::n] #downsample
    et = (data.index.to_series() - data.index.min()).dt.total_seconds()
    t = et.values
    dt = t[1]-t[0]
        
    omega0=6
    Fs = 1/dt
    N = len(data)
    f = np.arange(0,N/2+1)*Fs/N
    widths = omega0*Fs/(2*np.pi*f[1:])
    
    cwtmatr = signal.cwt(data, signal.morlet2, widths, w = omega0)
    #cwtmatr, freqs = wavelet.cwt(data,widths,"cmor2.5-1.0", sampling_period = dt)
    
    ha = plt.subplot(111)
    ha.pcolormesh(t,f,abs(cwtmatr[:,:-1])**2, cmap='turbo')
    ha.set_yscale('log')
    ha.set_ylim(f[1],f[-1])
    ha.set_xlabel('Time')
    ha.set_ylabel('Frequency [Hz]')
    
    
    mother = wavelet.Morlet(6)
    s0 = dt
    dj = 1/12
    J = 7/dj
    #alpha, _, _ = wavelet.ar1(data)  # Lag-1 autocorrelation for red noise
    
    wave, scales, freqs, coi, fft, fftfreqs = wavelet.cwt(data.values, dt, dj, s0, wavelet=mother)
    
    power = np.abs(wave)**2
    fft_power = np.abs(fft)**2
    period = (1/freqs)/3600 #hours
    """
    signif, fft_theor = wavelet.significance(1.0, dt, scales, 0, alpha,
                                         significance_level=0.95,
                                         wavelet=mother)
    sig95 = np.ones([1, N]) * signif[:, None]
    sig95 = power / sig95
    """
    plt.figure()
    ha = plt.subplot(111)
    ha.pcolormesh(t,period,power, cmap='turbo')
    #ha.contour(t,period,sig95)
    ha.set_yscale('log')
    ha.set_ylim(period.min(),period.max())
    ha.set_xlabel('Time (s)')
    #ha.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d-%h'))
    #for label in ax.get_xticklabels(which='major'):
     #   label.set(rotation=30, horizontalalignment='right')
    #ha.set_xticklabels(t,rotation=45)
    ha.set_ylabel('Period (hr)')
    ha.fill(np.concatenate([t, t[-1:] + dt, t[-1:] + dt,
                           t[:1] - dt, t[:1] - dt]),
        np.concatenate([coi/3600, [1e-1], period[-1:],
                           period[-1:], [1e-1]]),
        'k', alpha=0.3, hatch='x')
    