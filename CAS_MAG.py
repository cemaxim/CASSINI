#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 23 09:40:26 2024

@author: cemaxim
"""
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.markers import MarkerStyle
import numpy as np
from SPICE_POS import SPICE_POS
from CAS_ORBITS import ORBIT
import pycwt as wavelet
import os
import datetime
import scipy.integrate as integrate
from scipy.signal import find_peaks

def get_data(orbit_n = None, start_time = None, end_time = None): #get data from start and end time of orbit(s)
    #start = ORBIT(orbit_i)[0]
    if orbit_n is None:
        start = start_time
        end = end_time
    else:
        start = ORBIT(orbit_n)[0]
        end = ORBIT(orbit_n)[1]
 
    time = [start, end]
    mag_df = pd.read_pickle("/home/cemaxim/pickles/MAG_DF.pkl")
    mag_df = mag_df[time[0] : time[1]]
    bc_df = pd.read_pickle('/home/cemaxim/pickles/BC_ID.pkl')
    bc_df = bc_df.drop(bc_df.loc[bc_df.ID ==-2].index)
    df = SPICE_POS(mag_df,time)
    df = pd.merge_asof(df, bc_df, left_index=True, right_index=True, direction='backward')
    df.insert(0,'R',np.sqrt(df.x**2 + df.y**2 + df.z**2))
    return time, df, bc_df


def plot_orbit(data, id, R_min, n, time):
    fig, axs = plt.subplots(figsize=(6,6))  
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
        axs.plot(x, y, clr, linestyle=':')
        axs.set_xlabel('x ($R_S$)')
        axs.set_ylabel('y ($R_S$)')
        axs.plot(0,0,markersize = 15, color = 'black', markerfacecoloralt='white',marker=MarkerStyle("o", fillstyle='left'))
        plt.gca().set_aspect('equal')
        axs.grid(True, linestyle=':')
        axs.set_title('Cassini Orbit '+str(n)+'\n['+str(time[0])+' to '+str(time[1])+']')        
    data = data.loc[(data.R > R_min)]
    axs.plot(data.x.values,data.y.values,'-b', label='$R \geq$'+str(R_min))
    fig.tight_layout()
    plt.legend()

def save_orbits(n):
    folder_path = '/home/cemaxim/CAS/ORBITS'
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)    
        time, data, id = get_data(n)
        plot_orbit(data,id,25,n,time)     
        plt.savefig('/home/cemaxim/CAS/ORBITS/'+str(n))

        
def B_vs_R(data, R_min, win, start_time, end_time):   
    data = data[start_time : end_time]
    data = data.loc[(data.R > R_min) & (data.ID ==1)] 
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


def MFA(data, win, R_min): #mean field align
    data = data.loc[(data.R > R_min) & (data.ID ==1)]
    means = pd.DataFrame({'BX': data.BX.rolling(win,center=True).mean(),
                          'BY': data.BY.rolling(win,center=True).mean(),
                          'BZ': data.BZ.rolling(win,center=True).mean()})
    means = means.dropna()
    means_mags = np.sqrt(means.BX**2+means.BY**2+means.BZ**2)
    raw = data[['BX', 'BY', 'BZ']]
    raw = raw.reindex_like(means,None)
    raw_mags = np.sqrt(raw.BX**2+raw.BY**2+raw.BZ**2)
    li = []
    for t in range(len(means.index)):
        z_hat = means.iloc[t]/means_mags.iloc[t]
        z_hat = z_hat.to_numpy()
        R_hat = np.array([1,0,0])#raw.iloc[t]/raw_mags.iloc[t]
        y_hat = np.cross(z_hat,R_hat)/np.sqrt(np.sum(np.cross(z_hat,R_hat)**2))
        x_hat = np.cross(y_hat,z_hat)
      # Perturbation Components:
        db = raw.iloc[t].to_numpy()-means.iloc[t].to_numpy()
        bx = np.dot(db, x_hat) 
        by = np.dot(db, y_hat)
        bz = np.dot(db, z_hat)   
        li.append({'B_PERP1':bx, 'B_PERP2':by, 'B_PAR':bz})
    mfa_data = pd.DataFrame(li)
    mfa_data = mfa_data.set_index(means.index)
    data_R = data[['R']][mfa_data.index[0] : mfa_data.index[-1]]
    mfa_data = pd.merge_asof(data_R, mfa_data, left_index=True, right_index=True, direction='forward')
    mfa_data.insert(4,'B_PERP',np.sqrt(mfa_data.B_PERP1**2+mfa_data.B_PERP2**2))
    return mfa_data#, means, raw, means_mags, raw_mags

def MFA_plot(mfa_data):#, means, raw, means_mags, raw_mags):
    #means_mags = np.sqrt(means_mags.BX**2+means_mags.BY**2+means_mags.BZ**2)
    fig, axs = plt.subplots(4, sharex=True)
    fig.set_size_inches((12,8))
   
    fig.supxlabel('TIME')
    fig.suptitle('Magnetic Field Perturbations (Mean-Field-Aligned)')
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
    
    #p_mag = np.sqrt(mfa_data.B_PERP1**2+mfa_data.B_PERP2**2+mfa_data.B_PAR**2)
    #d = raw_mags - means_mags
    #axs[4].plot(p_mag)
    #d = raw - means
    #d_mag = np.sqrt(d.BX**2+d.BY**2+d.BZ**2)
    #axs[5].plot(d_mag)
    fig.tight_layout()


def CWT(data, n, dt, Rs=None):#freq='15min'): #Continuous Wavelet Analysis... Downsample data w freq
    freq = str(dt)+'min'
    data = data.asfreq(freq)
    data = data.interpolate()
    et = (data.index.to_series() - data.index.min()).dt.total_seconds() #elapsed time in seconds dataframe 
    #t = et.values
    t = data.index
    #dt = t[1]-t[0]
    
    
    omega0=6
    Fs = 1/dt
    N = len(data)
    f = np.arange(0,N/2+1)*Fs/N
    widths = omega0*Fs/(2*np.pi*f[1:])
    """   
    cwtmatr = signal.cwt(data, signal.morlet2, widths, w = omega0)
    #cwtmatr, freqs = wavelet.cwt(data,widths,"cmor2.5-1.0", sampling_period = dt)
 
    ha = plt.subplot(111)
    ha.pcolormesh(t,f,abs(cwtmatr[:,:-1])**2, cmap='turbo')
    ha.set_yscale('log')
    ha.set_ylim(f[1],f[-1])
    ha.set_xlabel('Time')
    ha.set_ylabel('Frequency [Hz]')
    """
    mother = wavelet.Morlet(6)
    #alpha, _, _ = wavelet.ar1(data)
    
    wave, scales, freqs, coi, fft, fftfreqs = wavelet.cwt(data.values, dt, wavelet=mother, freqs = f)
    
    power = np.abs(wave)**2
    #fft_power = np.abs(fft)**2
    period = (1/freqs)/60 #hours
    coi = coi**-1
  # Remove Cone of Infuluence
    for i, col in enumerate(power.T):
            col_num = len(col) - i
            coi_start_index = min(range(len(freqs)),
                                  key=lambda i: abs(freqs[i] - coi[col_num]))
            power[:coi_start_index, col_num] = np.zeros(coi_start_index) 
            
    #power = power.astype('float')
    #power[power ==0]= 'nan'
    """
   # Significance Test
    signif, fft_theor = wavelet.significance(1.0, dt, scales, 0, alpha,
                                         significance_level=0.95,
                                         wavelet=mother)
    sig95 = np.ones([1, N]) * signif[:, None]
    sig95 = power / sig95
    """
   
    
    #plt.figure()
    fig, axs = plt.subplots(2,1, sharex=True, figsize=(9,6), layout='constrained')
    fig.supxlabel('Date-Time')
    fig.suptitle(str(data.name)+' Wavelet Power Spectrum (Orbit '+str(n)+')')
    axs[0].plot(t,data)
    #axs[0].title.set_text(str(data.name)+'('+str(freq)+' down-sampled)')
    axs[0].set_ylabel('Magnitude (nT)')
    #axs[1].title.set_text(str(data.name)+' Wavelet Power Spectrum')
    axs[1].pcolormesh(t,period,power, cmap='turbo')
    #axs[1].contour(t,period,sig95)
    axs[1].set_yscale('log')
    axs[1].set_ylim(period.min(),period.max())
    axs[1].set_ylabel('Period (hr)')
    axs[1].plot([t.min(),t.max()],[10,10],'--', color='white',linewidth=0.25)
    axs[1].plot([t.min(),t.max()],[10.7,10.7],'--', color='white',linewidth=1)
    axs[1].plot([t.min(),t.max()],[15,15],'--',color='white',linewidth=0.25)
    axs[1].plot([t.min(),t.max()],[30,30],'--',color='white',linewidth=0.25)
    plt.xticks(rotation=30, horizontalalignment='right')
    plt.grid('--', linewidth=0.25)
    plt.colorbar(axs[1].pcolormesh(t,period,power, cmap='turbo'), label='($nT^2$/Hz)')

    #fig.tight_layout()
    
   
    """
    axs[1].fill(np.concatenate([t, t[-1:] + dt, t[-1:] + dt,
                           t[:1] - dt, t[:1] - dt]),
        np.concatenate([coi/3600, [1e-1], period[-1:],
                           period[-1:], [1e-1]]),
        'k', alpha=0.5, hatch='x')

    if Rs is not None:
        def R(t):
            R = Rs[t].to_numpy()
            return R
       
     """
    psd = np.zeros(len(freqs))
    np.nan_to_num(power)
    for i in range(0, len(freqs)):
        psd[i] = (2/len(freqs)) * integrate.trapz(power[i, :], range(0, len(power[i,:])))
    
    peaks, _ = find_peaks(psd)
    #print(period[peaks])
    max_period = period[psd.argmax()]
    print(max_period)
    plt.figure()
    plt.loglog(period, psd)
    #plt.plot(period[peaks],psd[peaks])
    for i in range(1,5):
        plt.loglog([max_period/i,max_period/i],[np.min(psd),np.max(psd)],':', label='m='+str(i))
    
    #plt.loglog([max_period/2,max_period/2],[np.min(psd),np.max(psd)],':', label='m=2')
    #plt.loglog([max_period/3,max_period/3],[np.min(psd),np.max(psd)],':', label='m=3')
    #plt.loglog([max_period/4,max_period/4],[np.min(psd),np.max(psd)],':', label='m=4')

    plt.title('Power Spectral Density')
    plt.xlabel('Period (Hr)')
    plt.ylabel('PSD ($nT^2$/Hz)')
    plt.legend()
   

    
def plots(R_min, win, dt, n = None, start_time=None, end_time=None):
    time , data, id = get_data(n,start_time,end_time)
    #folder_path = '/home/cemaxim/CAS/CWTs/'
   # if not os.path.exists(folder_path):
     #   os.makedirs(folder_path)
    #plot_orbit(data, id, R_min, n, time)
   # plt.savefig(folder_path+str(n)+'_ORBIT')
    mfa_data = MFA(data, win, R_min)
    #MFA_plot(mfa_data)
    CWT(mfa_data.B_PERP, n, dt)
   # plt.savefig(folder_path+str(n)+'_WPS')
    #for i in range(1,4):
      #  CWT(mfa_data.iloc[:,i],'15min')

def model_wave(n,T1,T2,T3):
    data = get_data(n)[1]
    data = data.drop(data.loc[data.R < 25].index)
    time = [data.index[0],data.index[-1]]
    f1 = 1/(T1*60)
    f2 = 1/(T2*60)
    f3 = 1/(T3*60)
    
    t = pd.date_range(time[0],time[1],freq='min')
    et = ((t.to_series() - t.min()).dt.total_seconds()).values/60
    y = np.sin(2*np.pi*f1*et) + 0.75*np.sin(2*np.pi*f2*et) + 0.5*np.sin(2*np.pi*f3*et)
    noise = np.zeros(len(y))
    noise = np.random.normal(0,1)
    wave = pd.DataFrame(data=y,index=t)
    wave.plot()
    return wave