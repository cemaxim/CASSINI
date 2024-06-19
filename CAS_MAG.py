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
from scipy import signal
from SPICE_POS import SPICE_POS
from CAS_ORBITS import ORBIT
import pycwt as wavelet
import os

def get_data(orbit_i,orbit_f=None): #get data from start and end time of orbit(s)
    start = ORBIT(orbit_i)[0]
    if orbit_f is not None:
       end = ORBIT(orbit_f)[1]
    else:
        end = ORBIT(orbit_i)[1]
    time = [start, end]
    mag_df = pd.read_pickle("/home/cemaxim/pickles/MAG_DF.pkl")
    mag_df = mag_df[time[0] : time[1]]
    bc_df = pd.read_pickle('/home/cemaxim/pickles/BC_ID.pkl')
    df = SPICE_POS(mag_df,time)
    df = pd.merge_asof(df, bc_df, left_index=True, right_index=True, direction='forward')
    df.insert(0,'R',np.sqrt(df.x**2 + df.y**2 + df.z**2))
    return time, df, bc_df


def plot_orbit(data, id, R_min, n, time):
    fig, axs = plt.subplots(figsize=(9,9))  
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


def CWT(data, freq='15min'): #Continuous Wavelet Analysis... Downsample data w freq
    data = data.asfreq(freq)
    data = data.interpolate()
    et = (data.index.to_series() - data.index.min()).dt.total_seconds() #elapsed time in seconds dataframe 
    t = et.values
    dt = t[1]-t[0]
        
    omega0=6
    Fs = 1/dt
    N = len(data)
    f = np.arange(0,N/2+1)*Fs/N
    widths = omega0*Fs/(2*np.pi*f[1:])
    
    cwtmatr = signal.cwt(data, signal.morlet2, widths, w = omega0)
    #cwtmatr, freqs = wavelet.cwt(data,widths,"cmor2.5-1.0", sampling_period = dt)
    """
    ha = plt.subplot(111)
    ha.pcolormesh(t,f,abs(cwtmatr[:,:-1])**2, cmap='turbo')
    ha.set_yscale('log')
    ha.set_ylim(f[1],f[-1])
    ha.set_xlabel('Time')
    ha.set_ylabel('Frequency [Hz]')
    """
    mother = wavelet.Morlet(6)
    alpha, _, _ = wavelet.ar1(data)
    
    wave, scales, freqs, coi, fft, fftfreqs = wavelet.cwt(data.values, dt, wavelet=mother, freqs = f)
    
    power = np.abs(wave)**2
    #fft_power = np.abs(fft)**2
    period = (1/freqs)/3600 #hours
    coi = coi**-1
  # Remove Cone of Infuluence
    for i, col in enumerate(power.T):
            col_num = len(col) - i
            coi_start_index = min(range(len(freqs)),
                                  key=lambda i: abs(freqs[i] - coi[col_num]))
            power[:coi_start_index, col_num] = np.zeros(coi_start_index)        
   # Significance Test
    signif, fft_theor = wavelet.significance(1.0, dt, scales, 0, alpha,
                                         significance_level=0.95,
                                         wavelet=mother)
    sig95 = np.ones([1, N]) * signif[:, None]
    sig95 = power / sig95
    
    plt.figure()
    axs = plt.subplot(111)
    axs.set_title(str(data.name)+' Wavelet Power Spectrum')
    axs.pcolormesh(t,period,power, cmap='turbo')
    #axs.contour(t,period,power) sig95
    axs.set_yscale('log')
    axs.set_ylim(period.min(),period.max())
    axs.set_xlabel('Time (s)')
    axs.set_ylabel('Period (hr)')
    axs.plot([t.min(),t.max()],[10,10],'--', color='white',linewidth=1)
    axs.plot([t.min(),t.max()],[15,15],'--',color='white',linewidth=1)
    axs.plot([t.min(),t.max()],[30,30],'--',color='white',linewidth=1)
    
    axs.fill(np.concatenate([t, t[-1:] + dt, t[-1:] + dt,
                           t[:1] - dt, t[:1] - dt]),
        np.concatenate([coi/3600, [1e-1], period[-1:],
                           period[-1:], [1e-1]]),
        'k', alpha=0.5, hatch='x')
    plt.tight_layout()
    #plt.colorbar(ha.pcolormesh(t,period,power, cmap='turbo'))

    
def plots(n, R_min, win, freq):
    time , data, id = get_data(n)
    plot_orbit(data, id, R_min, n, time)
    mfa_data = MFA(data, win, R_min)
    MFA_plot(mfa_data)
    CWT(mfa_data.iloc[:,1],'15min')

"""
def save_plots(data,id,R_min,win,freq='15min'):
    folder_path = '/home/cemaxim/CAS/ORBIT_'+str(n)
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
    plot_orbit(data, id, R_min)
    plt.savefig(folder_path+'/ORBIT')
    mfa = MFA(data, win, R_min)
    MFA_plot(mfa)
    plt.savefig(folder_path+'/MFA_PERTURBS')
    CWT(mfa.B_PERP1, freq)
    plt.savefig(folder_path+'/WPS')
"""
