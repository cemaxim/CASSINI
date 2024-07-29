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
from SPICE_POS import SPICE_POS
from CAS_ORBITS import ORBIT
import pycwt as wavelet
import os
import scipy.integrate as integrate
from scipy.signal import find_peaks
from scipy import signal as sig
import matplotlib.colors as mcolors
from scipy import stats



def get_data(orbit_n = None, start_time = None, end_time = None): 
    """
    Read in magnetometer and position data of an orbit or between a start and end time

    Parameters
    ----------
    orbit_n : int, optional
        Cassini's orbit number [0-120]. The default is None.
    start_time : str, optional
        Start time [Y-m-d H:M:S]. The default is None.
    end_time : str, optional
        End time [Y-m-d H:M:S]. The default is None.

    Returns
    -------
    time : list
        Start and end datetimes as strings.
    df : pandas.core.frame.DataFrame
        Dataframe containing 1 minute averaged magnetometer data.  
            R:
                Radial Saturn-Cassini distance(Rs)
            BX, BY, BZ:
                X, Y, and Z components of magetic-field(nT) in KSM cartesian coordinate system: X points from Saturn to Sun, the X-Z
                plane contains Saturn's centered magetic dipole axis, and Y completes the right handed set.
            BR, BPHI, BTHETA:
                R, Phi, and Theta components of magnetic-field(nT) in KRTP spherical coordinate system: R points from Saturn to Cassini,
                Phi is parallel to Saturn's equator and Theta completes the right handed set.
           x, y, z:
               Cassini's Position in KSM coordinates
            ID:
                Boundary condition IDs. 1=magnetosphere; 2=sheath; -1=solar wind
    ids : TYPE
        Boundary conditions dataframe. 

    """
    if orbit_n is None:
        start = start_time
        end = end_time
    else:
        start = ORBIT(orbit_n)[0]
        end = ORBIT(orbit_n)[1]
 
    time = [start, end]
    mag_df = pd.read_pickle("/home/cemaxim/pickles/MAG_DF.pkl")
    mag_df = mag_df[time[0] : time[1]]
    ids = pd.read_pickle('/home/cemaxim/pickles/BC_ID.pkl')
    ids = ids.drop(ids.loc[ids.ID ==-2].index)
    data = SPICE_POS(mag_df,time)
    data = pd.merge_asof(data, ids, left_index=True, right_index=True, direction='backward')
    data.insert(0,'R',np.sqrt(data.x**2 + data.y**2 + data.z**2))
    return time, data


def bc_ids():
    ids = pd.read_pickle('/home/cemaxim/pickles/BC_ID.pkl')
    ids = ids.drop(ids.loc[ids.ID ==-2].index)
    return ids

def RPWS_data(time):
    rpws = pd.read_pickle('/home/cemaxim/CAS/Codes/pickles/RPWS.pkl')
    rpws = rpws[time[0] : time[1]]
    return rpws

def CAPS_data(time,R_min=25):
    CAPS = pd.read_pickle('/home/cemaxim/CAS/Codes/pickles/CAPS.pkl')
    CAPS = CAPS[time[0] : time[1]]
    CAPS = CAPS.loc[(CAPS.R>R_min)]
    return CAPS.density

def MIMI_data(time):
    MIMI = pd.read_pickle('/home/cemaxim/CAS/Codes/pickles/MIMI.pkl')
    MIMI = MIMI[time[0] : time[1]]
    return MIMI

def date_format(axs):
    locator = mdates.AutoDateLocator(minticks=5, maxticks=20)
    formatter = mdates.ConciseDateFormatter(locator)
    axs.xaxis.set_major_locator(locator)
    axs.xaxis.set_major_formatter(formatter)

def plot_orbit(data, ids, R_min, n, time, plot=True, axs=None):
    """
    Plots orbit of Cassini spacecraft around Saturn.

    Parameters
    ----------
    data : pandas.core.frame.DataFrame
        Dataframe with Cassini's radial distance from Saturn in Rs and datetime index.
    ids : pandas.core.frame.DataFrame
        Dataframe with boundary condition ids and datetime index.
    R_min : int
        Highlight orbit outside of specified radial distance.
    n : int
        Cassini's orbit number.
    time : list
        Start and end time of plot.

    Returns
    -------
    None.

    """
    if plot is True:
        fig, axs = plt.subplots(figsize=(6,6))  
    for i in range(len(ids.index)-1):
        t = data.index
        t1 = ids.index[i]
        t2 = ids.index[i+1]
        wh = (t>t1) & (t<t2)
        x = data.iloc[wh]['x'].to_numpy()
        y = data.iloc[wh]['y'].to_numpy()
        
        if ids.iloc[i]['ID'] == -1:
            clr = '-r'
            # label = 'solar wind'
        if ids.iloc[i]['ID'] == 2:
            clr = 'lime'
            # label = 'sheath'
        if ids.iloc[i]['ID'] == 1:
            clr = '-b'
            # label = 'magnetosphere'
        axs.plot(x, y, clr)#, linestyle=':')
    axs.set_xlabel('x ($R_S$)')
    axs.set_ylabel('y ($R_S$)')
    axs.plot(0,0,markersize = 15, color = 'black', markerfacecoloralt='white',marker=MarkerStyle("o", fillstyle='left'))
    plt.gca().set_aspect('equal')
    axs.grid(True, linestyle=':')
    axs.set_title('Cassini Orbit '+str(n)+'\n['+str(time[0])+' to '+str(time[1])+']')        
    # data = data.loc[(data.R > R_min)]
    # axs.plot(data.x.values,data.y.values,'-b', label='$R \geq$'+str(R_min))
    # fig.tight_layout()
    # plt.legend()

def save_orbits(n):
    folder_path = '/home/cemaxim/CAS/ORBITS'
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)    
        time, data = get_data(n)
        ids = bc_ids()
        plot_orbit(data,ids,25,n,time)     
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


def MFA(data, win=60, R_min=25):
    """
    Mean-field align magnetometer data.

    Parameters
    ----------
    data : pandas.core.frame.DataFrame
        Dataframe with R, BX, BY, and BZ columns and datetime index.
    win : int
        Rolling mean window size.
    R_min : int, optional
        Limit magnetometer data to outside of specified radius.

    Returns
    -------
    mfa_data : pandas.core.frame.DataFrame
        Dataframe with mean-field-aligned perturbation components and datetime index.

    """
    data = data.loc[(data.R > R_min)] #& (data.ID ==1)]
    means = pd.DataFrame({'BX': data.BX.rolling(win,center=True).mean(),
                          'BY': data.BY.rolling(win,center=True).mean(),
                          'BZ': data.BZ.rolling(win,center=True).mean()})
    means = means.dropna()
    means_mags = np.sqrt(means.BX**2+means.BY**2+means.BZ**2)
    raw = data[['BX', 'BY', 'BZ']]
    raw = raw.reindex_like(means,None)
    li = []
    for t in range(len(means.index)):
        z_hat = means.iloc[t]/means_mags.iloc[t]
        z_hat = z_hat.to_numpy()
        R_hat = np.array([1,0,0])#raw.iloc[t]/raw_mags.iloc[t]
        y_hat = np.cross(z_hat,R_hat)/np.sqrt(np.sum(np.cross(z_hat,R_hat)**2))
        x_hat = np.cross(y_hat,z_hat)
        " Perturbation Components:"
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
    return mfa_data

def MFA_plot(mfa_data):    
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
    fig.tight_layout()

def get_B_PERP(n):
    time, data = get_data(n)
    mfa = MFA(data)
    B_PERP = mfa.B_PERP
    return B_PERP

def get_B_PAR(n):
    time, data = get_data(n)
    mfa = MFA(data)
    B_PAR = mfa.B_PAR
    return B_PAR

def downsample(signal, dt):
    freq = str(dt)+'min'
    signal = signal.asfreq(freq)
    signal = signal.interpolate()
    return signal    

def CWT(signal, dt=15):
    """
    Performs a continuous wavelet transfrom anaylsis on a signal.

    Parameters
    ----------
    signal : pandas.core.series.Series
        Signal to analyze with datetime index.
    dt : int
        Time interval in minutes between samples in singal.

    Returns
    -------
    periods : numpy.ndarray
        Periods range.
    power : numpy.ndarray
        2D power spectrum array.

    """
    freq = str(dt)+'min'
    signal = signal.asfreq(freq,'bfill')
    signal = signal.interpolate()
    Fs = 1/dt
    N = len(signal)
    f = np.arange(1,N/2)*Fs/N
     
    mother = wavelet.Morlet(6)
    wave, scales, freqs, coi, fft, fftfreqs = wavelet.cwt(signal.values, dt, wavelet=mother, freqs = f) 
    power = np.abs(wave)**2
    periods = (1/freqs)/60 #hours
    coi = coi**-1
    
    "Remove cone of infuluence"
    for i, col in enumerate(power.T):
            col_num = len(col) - i
            coi_start_index = min(range(len(freqs)),
                                  key=lambda i: abs(freqs[i] - coi[col_num]))
            power[:coi_start_index, col_num] = np.zeros(coi_start_index)        
            power[power ==0]= 'nan'
    return signal, periods, power

def CWT_plot(signal, axs=None):#title='$\delta$$b_\perp$ Wavelet Power Spectrum',
    signal, periods, power = CWT(signal)
    t = signal.index
    if axs is None:
        fig, axs = plt.subplots(1, figsize=(9,6), layout='constrained') 
    # axs.set_title(title)
    axs.pcolormesh(t,periods,power, cmap = 'turbo')
    axs.set_xlabel('Datetime')
    axs.set_yscale('log')
    axs.set_ylim(periods.min(),periods.max())
    axs.set_ylabel('Period (hr)',fontsize=8)
    axs.plot([t.min(),t.max()],[10,10],'--', color='white',linewidth=0.25)
    axs.plot([t.min(),t.max()],[10.7,10.7],'--', color='white',linewidth=1)
    axs.plot([t.min(),t.max()],[15,15],'--',color='white',linewidth=0.25)
    axs.plot([t.min(),t.max()],[30,30],'--',color='white',linewidth=0.25)
    date_format(axs)
    axs.grid('--', linewidth=0.25)
    plt.colorbar(axs.pcolormesh(t,periods,power,cmap='turbo'))#,norm=mcolors.Normalize(vmin=10e-6, vmax=1)), label='($nT^2$/Hz)')

def pgram(signal, T_min, T_max):
    '''
    Generate Lomb-Scargle periodogram from signal.

    Parameters
    ----------
    signal : pandas.core.series.Series
        Signal to analyze with datetime index.
    T_min/T_max : int
        Min/max period for analysis.

    Returns
    -------
    None.

    '''
    signal, periods, power = CWT(signal)
    power = np.nan_to_num(power)
    t = signal.index.values
    wh = np.where((periods>T_min) & (periods<T_max))[0]
    p = power[wh]
    p_T = p.T
    mean = np.array([])
    T = np.linspace(1,54,5000)
    T = T*3.6e12
    w = 2*np.pi/T
    for i in range(0, len(p_T)):
        mean = np.append(mean, p_T[i].mean())
    pgram = sig.lombscargle(t,mean,w,normalize=True)
    fig,axs = plt.subplots(3, layout = 'constrained')
    axs[0].pcolormesh(t,periods[wh],p,vmax=0.05)
    axs[0].grid(axis='x',linestyle='--',alpha=0.5)
    axs[0].set_ylabel('Period(hr)')
    axs[1].plot(t,mean)
    axs[1].set_ylabel('Power ($nT^2$/Hz)')
    date_format(axs[1])
    axs[1].grid(axis='x',linestyle='--',alpha=0.5)
    axs[0].sharex(axs[1])
    axs[0].tick_params(labelbottom=False)
    axs[2].plot(T/3.6e12,pgram)
    axs[2].set_ylabel('PSD')
    axs[2].set_xlabel('Period (hr)')
    axs[2].minorticks_on()
    axs[2].set_title('Lomb-Scargle Periodogram')
    
def find_nearest(array, value):
        array = np.asarray(array)
        idx = (np.abs(array - value)).argmin()
        return array[idx]
    
def lin_corr(signal_a,signal_b):
    
    # signals = [signal_a,signal_b]
    # for i in signals:
    #     signal, periods, power = CWT(i)
    #     wh = np.where(periods == find_nearest(periods,1))[0]
    #     p = power[wh]
    #     t = signal.index
        
    signal_a, periods_a, power_a = CWT(signal_a)
    signal_b, periods_b, power_b = CWT(signal_b)
    wh_a = np.where(periods_a == find_nearest(periods_a,3))[0]
    wh_b = np.where(periods_b == find_nearest(periods_b,3))[0]
    p_a = np.nan_to_num(power_a[wh_a])
    p_b = np.nan_to_num(power_b[wh_b])
    t_a = signal_a.index
    t_b = signal_b.index
    fig, axs = plt.subplots(2,sharex=True)
    axs[0].plot(t_a,p_a[0])
    axs[0].set_yscale('log')

    axs[0].plot(t_b,p_b[0])
    axs[0].set_yscale('log')
    date_format(axs[1])
    coeff = stats.pearsonr(p_a[0,0:1000],p_b[0,0:1000])
    print (coeff)
    
   
     
def PSD(signal, plot=True, axs=None):
    signal, periods, power = CWT(signal)
    psd = np.zeros(len(periods))
    power = np.nan_to_num(power)
    for i in range(0, len(periods)):
        psd[i] = (2/len(periods)) * integrate.trapz(power[i, :], range(0, len(power[i,:])))
    peaks, _ = find_peaks(psd,prominence=1e-2)
    # prominences = peak_prominences(psd, peaks)[0]
    
    # print(period[peaks])
    # print(prominences)
    fund_period = 10.7    
    if plot is True:
        fig, axs = plt.subplots()
    if axs is not None or plot is True:
        axs.loglog(periods, psd)
        # axs.plot(period[peaks], psd[peaks], 'x')
        for i in range(1,5):
            axs.loglog([fund_period/i,fund_period/i],[np.min(psd),np.max(psd)],':', label='m='+str(i))
        axs.set_title('Power Spectral Density')
        axs.set_xlabel('Period (Hr)')
        axs.set_xlim(xmax=periods.max())
        axs.set_ylabel('PSD ($nT^2$/Hz)')
        axs.legend()
    return psd

def peaks_hist(axs=None):
    orbs = [24,25,26,27,28,29,30,33,36,37,112,113,114,119]
    a = np.array([])
    for i in enumerate(orbs):
        time, data = get_data(i[1])
        mfa = MFA(data,60,25)
        signal, periods, power = CWT(mfa.B_PERP, 15)
        psd = PSD(signal, plot=False, axs=axs[i[0]])
        peaks, _ = find_peaks(psd,prominence=1e-2)
        peaks = periods[peaks]
        peaks = peaks[peaks < 50]
        a = np.concatenate([a,peaks])
    # np.save('/home/cemaxim/CAS/Codes/peaks',a)
    # a = np.load('/home/cemaxim/CAS/Codes/peaks.npy')
    fig, axs = plt.subplots()
    bins = np.arange(0,50)
    axs.hist(a, bins)#, log=True)
    axs.minorticks_on()
    axs.set_xlabel('Period (Hr)')
    axs.set_ylabel('Counts')
    axs.yaxis.set_tick_params(which='minor',left=False)
    

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
    # noise = np.zeros(len(y))
    # noise = np.random.normal(0,1)
    wave = pd.DataFrame(data=y,index=t)
    wave.plot()
    return wave


def plots(R_min=25, win=60, dt=15, n = None, start_time=None, end_time=None):
    time , data = get_data(n,start_time,end_time)
    mfa = MFA(data, win, R_min)
    signal = mfa.B_PERP
    signal_ds = downsample(signal,dt)
    fig, ax = plt.subplots(2, sharex = True, layout='constrained', figsize = (9,6))
    title = '$\delta$$b_\perp$ '
    if n is not None:
        title += 'Orbit '+str(n)
    else:
        title += '['+str(time[0])+' to '+str(time[1])+']'
    ax[0].plot(signal_ds.index,signal_ds)
    ax[0].set_title(title)
    CWT_plot(signal,axs=ax[1])
    # plt.savefig('/home/cemaxim/CAS/CWTs/'+str(n)+'_CWT')
    # plot_orbit(data, ids, R_min, n ,time, plot=False, axs=ax[3])


def CAPS_CWT():
    time = ORBIT(114)
    caps = CAPS_data(time)
    B_PERP = get_B_PERP(114)
    fig, axs = plt.subplots(3,sharex=True,layout='constrained',figsize=(8,5))
    fontsize = 8
    axs[0].plot(caps*10e-3,linewidth=0.4)
    axs[0].set_yscale('log')
    axs[0].set_ylabel('Electron Density ($m^{-3}$)',fontsize=fontsize)
    CWT_plot(B_PERP,axs=axs[2])
    axs[1].plot(B_PERP,linewidth=0.4)
    axs[1].set_ylabel('$\delta$$b_\perp$ Magnitude (nT)',fontsize=fontsize)
    # axs[1].set_title('$\delta$$b_\perp$')
    # axs[0].set_title('CAPS Electron Density')

def pgrams():
    orbs = [24,25,26,27,28,29,30,33,36,37,112,113,114,119]
    for i in orbs:
        B_PERP = get_B_PERP(i)
        pgram(B_PERP,0,1)
        # plt.savefig('/home/cemaxim/CAS/CWTs/db_perp/'+str(i)+'_pgram')
        
def mag_vs_den():
    B_PERP = get_B_PERP(114)
    e_DENS = CAPS_data(ORBIT(114))
    x = pd.merge_asof(B_PERP,e_DENS,left_index=True, right_index=True, direction='backward')
    a = x.B_PERP.rolling(50,center=True).mean()
    a = a.dropna()
    b = x.density.rolling(50,center=True).mean()
    b = b.dropna()
    # peaks_a, _ = find_peaks(a,prominence=0.)
    # peaks_b, _ = find_peaks(b,prominence=2)
    fig, axs = plt.subplots(2,sharex=True)
    axs[0].plot(a)
    # axs[0].plot(a.iloc[peaks_a],'x')
    # axs[1].plot(x.density,alpha=0.5)
    # axs[1].plot(b.iloc[peaks_b],'x')
    axs[1].plot(b)
    axs[1].set_yscale('log')
    date_format(axs[1])
    corr = stats.pearsonr(a,b)
    print(corr)

def data_plots():
    CAPS = CAPS_data(ORBIT(114))
    MIMI = MIMI_data(ORBIT(114))
    fig,axs = plt.subplots(2,sharex=True)
    axs[0].plot(CAPS)
    axs[0].set_yscale('log')
    axs[1].plot(MIMI.A0)
    axs[1].set_yscale('log')
    date_format(axs[1])