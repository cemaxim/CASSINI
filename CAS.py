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
import scipy.integrate as integrate
from scipy.signal import find_peaks
from scipy import signal as sig
import matplotlib.colors as mcolors
from scipy.signal import hilbert



def get_data(orbit_n = None, start_time = None, end_time = None): 
    """
    Read in magnetometer and position data of an orbit or between a start and end time.

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

    """
    if orbit_n is None:
        start = start_time
        end = end_time
    else:
        start = ORBIT(orbit_n)[0]
        end = ORBIT(orbit_n)[1]
 
    time = [start, end]
    mag_df = pd.read_pickle("/home/cemaxim/CAS/Scrpts/pickles/MAG_2.pkl")
    mag_df = mag_df[time[0] : time[1]]
    ids = pd.read_pickle('/home/cemaxim/CAS/Scrpts/pickles/BC_ID.pkl')
    ids = ids.drop(ids.loc[ids.ID ==-2].index)
    data = SPICE_POS(mag_df,time)
    data = pd.merge_asof(data, ids, left_index=True, right_index=True, direction='backward')
    data.insert(0,'R',np.sqrt(data.x**2 + data.y**2 + data.z**2))
    return time, data


def bc_ids():
    ids = pd.read_pickle('/home/cemaxim/CAS/Scrpts/pickles/BC_ID.pkl')
    ids = ids.drop(ids.loc[ids.ID ==-2].index)
    return ids

def RPWS_data(n):
    time = ORBIT(n)
    rpws = pd.read_pickle('/home/cemaxim/CAS/Scrpts/pickles/RPWS.pkl')
    rpws = rpws[time[0] : time[1]]
    return rpws

def CAPS_data(n,R_min=25):
    time = ORBIT(n)
    CAPS = pd.read_pickle('/home/cemaxim/CAS/Scrpts/pickles/CAPS.pkl')
    CAPS = CAPS[time[0] : time[1]]
    CAPS = CAPS.loc[(CAPS.R>R_min)]
    CAPS = CAPS.rolling(10,center=True).mean()
    return CAPS.density

def MIMI_data(time):
    MIMI = pd.read_pickle('/home/cemaxim/CAS/Scrpts/pickles/MIMI.pkl')
    MIMI = MIMI[time[0] : time[1]]
    MIMI = MIMI.drop(MIMI.index[MIMI.E0 == -1e38])
    return MIMI

def date_format(axs):
    locator = mdates.AutoDateLocator(minticks=5, maxticks=20)
    formatter = mdates.ConciseDateFormatter(locator)
    axs.xaxis.set_major_locator(locator)
    axs.xaxis.set_major_formatter(formatter)

def plot_orbit(n, R_min, plot=True, axs=None):
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
    time, data = get_data(n)
    ids = bc_ids()
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
        axs.plot(x, y, clr, linestyle=':')
    axs.set_xlabel('x ($R_S$)')
    axs.set_ylabel('y ($R_S$)')
    axs.plot(0,0,markersize = 15, color = 'black', markerfacecoloralt='white',marker=MarkerStyle("o", fillstyle='left'))
    plt.gca().set_aspect('equal')
    axs.grid(True, linestyle=':')
    axs.set_title('Cassini Orbit '+str(n)+'\n['+str(time[0])+' to '+str(time[1])+']')        
    data = data.loc[(data.R > R_min)]
    axs.plot(data.x.values,data.y.values,'-b', label='$R \geq$'+str(R_min))
    # fig.tight_layout()
    # plt.legend()
"""
def save_orbits(n):
    folder_path = '/home/cemaxim/CAS/ORBITS'
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)    
        time, data = get_data(n)
        ids = bc_ids()
        plot_orbit(data,ids,25,n,time)     
        plt.savefig('/home/cemaxim/CAS/ORBITS/'+str(n))

        
def B_vs_R(da  a = signal_a.asfreq(freq='60s',method='bfill')
    b = signal_b.asfreq(freq='60s',method='bfill')
    # b = signal_b.asof(a.index)
    # a = signal_a
    # b = signal_b
    a = a.dropna()
    b = b.dropna()
    
    fs = 1/60#(a.index[1]-a.index[0]).total_seconds()
    window = 'hann'   
    
    nperseg = 256*m
    #noverlap = nperseg/2
    f, Pxy = sig.csd(a,b,fs=fs,window=window,nperseg=nperseg)#,noverlap=noverlap)
    f2, Cxy = sig.coherence(a,b,fs=fs,window=window,nperseg=nperseg)#,noverlap=noverlap)
    T = 1/f/3600
    T = T[1:]
    Pxy = Pxy[1:]  
    Cxy = Cxy[1:]ta, R_min, win, start_time, end_time):   
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
"""

def MFA(data, win=30, R_min=25):
    """
    Mean-field align magnetometer data (Khurana & Kivelson 1989).

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
    # mfa_data.insert(5,'B_PERP1',)
    return mfa_data

def MFA_CAPS(data, win=30):
    means = data.rolling(win,center=True).mean()
    means = means.dropna()
    raw = data
    raw = raw.reindex_like(means,None)
    li = []
    for t in range(len(means.index)):
 
        " Perturbation Components:"
        dN = raw.iloc[t]-means.iloc[t]
  
        li.append({'dN':dN})
    mfa_data = pd.DataFrame(li)
    mfa_data = mfa_data.set_index(means.index)
    return mfa_data.dN
    
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

def get_MAG(n):
    time,data = get_data(n)
    mag = MFA(data)
    return mag
    
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

def get_dN(n):
    data = CAPS_data(n)
    data = MFA_CAPS(data)
    return data

def downsample(signal, dt):
    freq = str(dt)+'min'
    signal = signal.asfreq(freq)
    signal = signal.interpolate()
    return signal    

def align(signal_a,signal_b):
    a = signal_a.asfreq(freq = '60s', method = 'bfill')
    b = signal_b.asof(a.index)
    a = a.dropna()
    b = b.dropna()
    return a,b

def split(signal, date):
    a = signal[signal.index[0] : date]
    b = signal[date : signal.index[-1]]
    return a, b

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
    title = '$\delta$$b_\parallel$'
    axs.set_title(title)
    n=mcolors.Normalize(vmin=0.001, vmax=0.18)
    axs.pcolormesh(t,periods,power, cmap = 'turbo',norm=n)
    axs.set_xlabel('Datetime')
    axs.set_yscale('log')
    axs.set_ylim(periods.min(),periods.max())
    axs.set_ylabel('Period (hr)',fontsize=10)
    axs.plot([t.min(),t.max()],[10,10],'--', color='white',linewidth=0.25)
    axs.plot([t.min(),t.max()],[10.7,10.7],'--', color='white',linewidth=1)
    axs.plot([t.min(),t.max()],[15,15],'--',color='white',linewidth=0.25)
    axs.plot([t.min(),t.max()],[30,30],'--',color='white',linewidth=0.25)
    date_format(axs)
    axs.grid('--', linewidth=0.25)
    plt.colorbar(axs.pcolormesh(t,periods,power,cmap='turbo',norm=n),label='($nT^2$/Hz)')
    
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
    T = np.linspace(1,120,1000)
    T = T*3.6e12
    w = 2*np.pi/T
    for i in range(0, len(p_T)):
        mean = np.append(mean, p_T[i].mean())
        # mean = np.nan_to_num(mean)
    pgram = sig.lombscargle(t,mean,w,normalize=True)
    peaks,_ = find_peaks(pgram, prominence=1e-2)
    fig,axs = plt.subplots(3, layout = 'constrained')
    axs[0].pcolormesh(t,periods[wh],p)#, norm = mcolors.Normalize(vmin=0.001, vmax=1e8))
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
    axs[2].plot((T/3.6e12)[peaks],pgram[peaks],'x')
    return((T/3.6e12)[peaks])
    
def pgram_hist():
    orbs = [24,25,26,27,28,29,30,33,36,37,112,113,114,119]
    peaks = np.array([])
    for i in enumerate(orbs):
        signal = get_B_PERP(i[1])
        peak = pgram(signal,2,4)
        peaks = np.concatenate([peaks,peak])
    fig, axs = plt.subplots()
    bins = np.arange(0,50)
    axs.hist(peaks, bins)#, log=True)
    axs.minorticks_on()
    axs.set_xlabel('Period (Hr)')
    axs.set_ylabel('Counts')
    axs.yaxis.set_tick_params(which='minor',left=False)
    return(peaks)
    
    
def find_nearest(array, value):
        array = np.asarray(array)
        idx = (np.abs(array - value)).argmin()
        return array[idx]
    # aa = np.array_split(a,n)
   # bb = np.array_split(b,n)
     
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
    

def model_wave(n,T1,T2,T3,phase):
    data = get_data(n)[1]
    data = data.drop(data.loc[data.R < 25].index)
    time = [data.index[0],data.index[-1]]
    
    f1 = 1/(T1*60)
    if T2 == 0:
        f2 = 0
    else:
        f2 = 1/(T2*60)
    if T3 == 0:
        f3 = 0
    else:
        f3 = 1/(T3*60)
    
    t = pd.date_range(time[0],time[1],freq='min')
    et = ((t.to_series() - t.min()).dt.total_seconds()).values/60
    y = np.sin(2*np.pi*f1*et+phase) + 0.75*np.sin(2*np.pi*f2*et+phase) + 0.5*np.sin(2*np.pi*f3*et+phase)
    # noise = np.zeros(len(y))
    # noise = np.random.normal(0,1).asfreq(freq='60s',method='bfill')
    wave = pd.DataFrame(data=y,index=t)
    wave.plot()
    return wave[0]


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


def cpsd(signal_a,signal_b,m,fs):
   
    a = signal_a
    b = signal_b
    a = a.dropna()
    b = b.dropna()
    
    window = 'hann'   
    
    nperseg = 256*m
    #noverlap = nperseg/2
    f, Pxy = sig.csd(a,b,fs=fs,window=window,nperseg=nperseg)#,noverlap=noverlap)
    f2, Cxy = sig.coherence(a,b,fs=fs,window=window,nperseg=nperseg)#,noverlap=noverlap)
    T = 1/f/3600

    """ 
    aa = np.array_split(a,n)
    bb = np.array_split(b,n)
    A = []
    B = []
    D = []
    for i in range(0,len(aa)):
        f, Pxy = sig.csd(aa[i],bb[i],fs=fs,window=window,nperseg=nperseg,noverlap=noverlap)
        Pxy = Pxy[1:]label='($nT^2$/Hz)
        f, Cxy = sig.coherence(aa[i],bb[i],fs=fs,window=window,nperseg=nperseg,noverlap=noverlap)
        P = np.tile(np.abs(np.real(Pxy)),(len(aa[i]),1)).T
        Cxy = Cxy[1:]
        C = np.tile(np.abs(Cxy),(len(aa[i]),1)).T
        B.append(C)
        A.append(P)
        ang = np.abs(np.angle(Pxy))
        ang = np.tile(ang,(len(aa[i]),1)).T
        D.append(ang)
    ANG = np.concatenate(D,axis=1)
    P = np.concatenate(A,axis=1)
    C = np.concatenate(B,axis=1)
    
    f, Pxy1 = sig.csd(aa[0],bb[0],fs=fs,window=window,nperseg=nperseg,noverlap=noverlap)
    f, Pxy2 = sig.csd(aa[1],bb[1],fs=fs,window=window,nperseg=nperseg,noverlap=noverlap)
    
    T = 1/f/3600
    T = T[1:]
    # Pxy1 = Pxy1[1:]  
    # Pxy2 = Pxy2[1:]
    # Cxy1 = Cxy1[1:]
    # P1 = np.tile(np.abs(np.real(Pxy1)),(len(aa[0]),1)).T
    # P2 = np.tile(np.abs(np.real(Pxy2)),(len(aa[1]),1)).T
    # # P = np.concatenate((P1,P2),axis=1)
    # # print(T.shape, t.shape, spect.shape)
    fig, ax = plt.subplots(3,sharex=True)
    ax[0].pcolormesh(a.index,T,P,cmap='binary')#,norm = 'log')
    ax[0].set_yscale('log')
    ax[0].set_title('CPSD')
    date_format(ax[0])
    ax[1].pcolormesh(a.index,T,C,cmap='binary')
    ax[1].set_yscale('log')
    ax[1].set_title('Coherence')
    ax[2].pcolormesh(a.index,T,ANG,cmap='binary')
    ax[2].set_title('Phase')
    ax[2].set_yscale('log')
    ax[2].set_ylim(1,40)
    ax[1].set_ylim(1,40)
    ax[0].set_ylim(1,40)
    
    """
    # peaks, _ = find_peaks(np.abs(Pxy),prominence=1e6)#,width=2)
    fig, axs = plt.subplots(5, figsize = (9,6),layout='constrained')
    axs[0].plot(a)#,linewidth=0.5)
    axs[0].sharex(axs[1])
    axs[0].tick_params(labelbottom=False)
    axs[0].set_ylabel('$\delta$$B_\parallel$ [nT]')
    axs[1].plot(b)#,linewidth=0.5)
    # date_format(axs[1])
    axs[1].set_ylabel('$\delta$N [$m^{-3}$]')
    axs[1].ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    axs[1].set_xlabel('Hours')
    axs[2].sharex(axs[4])
    axs[2].tick_params(labelbottom=False)
    axs[2].plot(T,np.abs(Pxy))
    axs[2].set_ylabel('CPSD [$nT^2$/Hz]')
    axs[2].set_ylim(bottom=10**-7)
    # axs[2].plot((1/f/3600)[peaks], np.abs(Pxy)[peaks],'x')
    axs[2].tick_params(axis='x',which='both',top=True,labeltop=True)
    axs[3].sharex(axs[4])

    axs[4].set_xscale('log')
    axs[3].tick_params(labelbottom=False)
    axs[3].plot(T, np.abs(Cxy))
    

    # axs[3].plot(T,np.zeros(len(Cxy))+0.75,linestyle='--')
    axs[3].set_ylabel('Coherence')
    axs[3].set_ylim(0,1.1)
    axs[2].set_yscale('log')
    axs[4].set_xlabel('Period [Hour]')
    axs[4].plot(T,np.abs(np.angle(Pxy))/np.pi)
    # axs[4].plot(T,np.angle(Pxy)/np.pi)

    axs[4].set_ylabel('Phase')
    axs[4].set_xlim(0.2,10)
    axs[4].minorticks_on()
    for i in range(2,5):
        axs[i].axvline(4,linestyle='--',alpha=0.5)
        axs[i].axvline(1,linestyle='--',alpha=0.5)

    axs[4].grid(axis='x',linestyle='--',alpha=0.5)
    axs[3].grid(axis='x',linestyle='--',alpha=0.5)
    axs[2].grid(axis='x',linestyle='--',alpha=0.5)
    print(nperseg/60)
    


def get_hht(signal):
    import pyhht
    from pyhht.visualization import plot_imfs
    decomposer = pyhht.EMD(signal)
    imfs = decomposer.decompose()
    
    t = []
    for n in range (0,len(signal)):
        dt = signal.index[n]-signal.index[0]
        t.append(dt)
    t = pd.Series(t).dt.total_seconds().values
    
    plot_imfs(signal, imfs,t)

    s = np.shape(imfs)
    fig, ax = plt.subplots(s[0]-1,1)
    fig1, ax1 = plt.subplots(s[0]-1,1)
    T = np.zeros((len(imfs), len(t)-1))
    hs = np.zeros((len(imfs), len(t)))
    for i in range (0,s[0]-1):
        ax[i].plot(t/3600,imfs[i])
        
        analytic_signal = hilbert(imfs[i])
        amplitude_envelope = np.abs(analytic_signal)
        instantaneous_phase = np.unwrap(np.angle(analytic_signal))
        instantaneous_frequency = np.diff(instantaneous_phase) / (2 * np.pi * (t[1] - t[0]))

        # Plot the instantaneous frequency
        p =  (1/instantaneous_frequency)/3600 #hours
        T[i,:] = p
        hs[i, :] = amplitude_envelope ** 2
        ax1[i].plot(t[1:]/3600, smooth(p,10))
        p = p[200:-200]
        wh = p > 0
        med = np.round(np.median(smooth(p[wh],10)),2)
        sd = np.round(np.std(smooth(p[wh],10)),2)
        ax[i].text(np.max(t/3600),0.0,str(med)+' '+str(sd))
        ax1[i].set_ylim([0.5,100])
        ax1[i].set_yscale('log')
    for j in range (0,s[0]-2):
        ax[j].set_xticklabels([])
    ax[i].set_xlabel('Time (hr)')
    ax[int(s[0]/2)-1].set_ylabel('delta B (nT)')
    ax1[i].set_xlabel('Time (hr)')
    ax1[int(s[0]/2)-1].set_ylabel('Period')

    return imfs

def smooth(y, box_pts):
        box = np.ones(box_pts)/box_pts
        y_smooth = np.convolve(y, box, mode='same')
        return y_smooth

def cpsd_data(n,win):
    
    time,mag = get_data(n)
    mag = mag.loc[(mag.R>25)]
    density = CAPS_data(n)
    mag, density = align(mag,density)
    mag = mag.rolling(win,center=True).mean()
    density = density.rolling(win,center=True).mean()
    dB = MFA(mag)
    dN = MFA_CAPS(density)
    # dB,dN = align(dB,dN)
    return dB, dN
    
def phase_hist():
    orbs = [24,25,26,27,28,29,30,33,36,37,112,113,114,119]
    x = np.array([])
    
    for n in orbs:
        signal_a,signal_b = cpsd_data(n,30)
        T,Cxy,phase = cpsd(signal_a.B_PAR,signal_b,12,1/60)
        wh = np.where(T<4 and T>2)
        wh = np.where(Cxy>0.5)
        y = phase[wh]
        x = np.concatenate([x,y])
    fig,axs = plt.subplots()
    bins = np.linspace(0,1,20)
    axs.hist(x,bins)
    print(x)

