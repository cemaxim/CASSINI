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
import matplotlib.colors as mcolors
import numpy as np
from SPICE_POS import SPICE_POS
from CAS_ORBITS import ORBIT
import pycwt as wavelet
import os
import scipy.integrate as integrate
from scipy.signal import find_peaks
from scipy.signal import peak_prominences

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
    return time, data, ids


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
        time, data, ids = get_data(n)
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


def MFA(data, win, R_min=None):
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


def CWT(data, dt, n=None, Rs=None, plot=True, axs=None):
    """
    Performs a continuous wavelet transfrom anaylsis on a signal.

    Parameters
    ----------
    data : pandas.core.series.Series
        Signal to analyze with datetime index.
    n : int
        Orbit number for plot title.
    dt : int
        Time interval in minutes between samples in singal.
    Rs : TYPE, optional
        DESCRIPTION. The default is None.

    Returns
    -------
    freqs : numpy.ndarray
        Frequencies range.
    power : numpy.ndarray
        Power Spectrum.

    """
    freq = str(dt)+'min'
    data = data.asfreq(freq)
    data = data.interpolate()
    t = data.index    
    
  
    Fs = 1/dt
    N = len(data)
    f = np.arange(1,N/2)*Fs/N
    
    # omega0=6
    # widths = omega0*Fs/(2*np.pi*f[1:])     
    # cwtmatr = signal.cwt(data, signal.morlet2, widths, w = omega0)
    # #cwtmatr, freqs = wavelet.cwt(data,widths,"cmor2.5-1.0", sampling_period = dt)
 
    # ha = plt.subplot(111)
    # ha.pcolormesh(t,f,abs(cwtmatr[:,:-1])**2, cmap='turbo')
    # ha.set_yscale('log')
    # ha.set_ylim(f[1],f[-1])
    # ha.set_xlabel('Time')
    # ha.set_ylabel('Frequency [Hz]')
    
    mother = wavelet.Morlet(6)
    wave, scales, freqs, coi, fft, fftfreqs = wavelet.cwt(data.values, dt, wavelet=mother, freqs = f) 
    power = np.abs(wave)**2
    period = (1/freqs)/60 #hours
    coi = coi**-1
    
    "Remove cone of infuluence"
    for i, col in enumerate(power.T):
            col_num = len(col) - i
            coi_start_index = min(range(len(freqs)),
                                  key=lambda i: abs(freqs[i] - coi[col_num]))
            power[:coi_start_index, col_num] = np.zeros(coi_start_index)        
            power[power ==0]= 'nan'
    # "Significance Test"
    # alpha, _, _ = wavelet.ar1(data)
    # signif, fft_theor = wavelet.significance(1.0, dt, scales, 0, alpha,
    #                                      significance_level=0.95,
    #                                      wavelet=mother)
    # sig95 = np.ones([1, N]) * signif[:, None]
    # sig95 = power / sig95
    
    if plot is True:
        fig, axs = plt.subplots(2,1, sharex=True, figsize=(9,6), layout='constrained')
        fig.supxlabel('Date-Time')
        fig.suptitle('$b_\perp$ Wavelet Power Spectrum (Orbit '+str(n)+')')
        plt.xticks(rotation=30, horizontalalignment='right')
    axs[0].set_title('['+str(dt)+'-minute sampling interval]', fontsize=10)
    axs[0].plot(t,data)
    axs[0].set_ylabel('Magnitude (nT)')
    # power = np.log(power)
    wh =np.logical_not(np.isnan(power))
    axs[1].pcolormesh(t,period,power, cmap = 'turbo')#, norm=mcolors.LogNorm(vmin=power[wh].min(), vmax=power[wh].max()), cmap='turbo')
    #axs[1].contour(t,period,sig95)
    axs[1].set_xlabel('Datetime')
    axs[1].set_yscale('log')
    axs[1].set_ylim(period.min(),period.max())
    axs[1].set_ylabel('Period (hr)')
    axs[1].plot([t.min(),t.max()],[10,10],'--', color='white',linewidth=0.25)
    axs[1].plot([t.min(),t.max()],[10.7,10.7],'--', color='white',linewidth=1)
    axs[1].plot([t.min(),t.max()],[15,15],'--',color='white',linewidth=0.25)
    axs[1].plot([t.min(),t.max()],[30,30],'--',color='white',linewidth=0.25)
    
    plt.grid('--', linewidth=0.25)
    plt.colorbar(axs[1].pcolormesh(t,period,power,  norm=mcolors.LogNorm(vmin=10e-6, vmax=power[wh].max()),cmap='turbo'), label='($nT^2$/Hz)')
    return freqs, power


def PSD(freqs, power, plot=True, axs=None):
    psd = np.zeros(len(freqs))
    power = np.nan_to_num(power)
    for i in range(0, len(freqs)):
        psd[i] = (2/len(freqs)) * integrate.trapz(power[i, :], range(0, len(power[i,:])))
    period = (1/freqs)/60
    peaks, _ = find_peaks(psd,prominence=1e-3)
    prominences = peak_prominences(psd, peaks)[0]
    
    print(period[peaks])
    print(prominences)
    # fund_period = period[psd.argmax()]
    fund_period = 10.7    
    if plot is True:
        fig, axs = plt.subplots()
    axs.loglog(period, psd)
    axs.plot(period[peaks], psd[peaks], 'x')
    for i in range(1,5):
        axs.loglog([fund_period/i,fund_period/i],[np.min(psd),np.max(psd)],':', label='m='+str(i))
    axs.set_title('Power Spectral Density')
    axs.set_xlabel('Period (Hr)')
    axs.set_xlim(xmax=period.max())
    axs.set_ylabel('PSD ($nT^2$/Hz)')
    axs.legend()
    
    
# def plots(axs, R_min=25, win=30, dt=15, n = None, start_time=None, end_time=None):
    
#     folder_path = '/home/cemaxim/CAS/CWTs/'
#     if not os.path.exists(folder_path):
#         os.makedirs(folder_path)
        
#     time , data, ids = get_data(n,start_time,end_time)
#     plot_orbit(data, ids, R_min, n, time)
#     plt.savefig(folder_path+str(n)+'_ORBIT')
#     mfa_data = MFA(data, win, R_min)
#     MFA_plot(mfa_data)
#     CWT(mfa_data.B_PERP, dt, axs, n)    
#     plt.savefig(folder_path+str(n)+'_CWT')
#     freqs , power = CWT(mfa_data.B_PERP, dt, axs, n, plot=False)
#     PSD(freqs,power,axs)
#     plt.savefig(folder_path+str(n)+'_PSD')


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
    time , data, ids = get_data(n,start_time,end_time)
    mfa = MFA(data, win, R_min)
    # plot_orbit(data, ids, R_min, n, time)
    fig, ax = plt.subplots(3, layout='constrained', figsize = (9,6))
    title = '$b_\perp$ Wavelet Power Spectrum '
    if n:
        title += 'Orbit '+str(n)
    else:
        title += '['+str(time[0])+' to '+str(time[1])+']'
    fig.suptitle(title)
    # fig.suptitle('$b_\perp$ Wavelet Power Spectrum (Orbit '+str(n)+')')
    freqs, power = CWT(mfa.B_PERP, 15, plot=False, axs=ax)
    ax[0].tick_params(labelbottom=False)
    ax[1].sharex(ax[0])
    locator = mdates.AutoDateLocator(minticks=5, maxticks=20)
    formatter = mdates.ConciseDateFormatter(locator)
    ax[1].xaxis.set_major_locator(locator)
    ax[1].xaxis.set_major_formatter(formatter)
    # ax[1].tick_params(labelrotation=30)

    PSD(freqs, power, plot=False, axs=ax[2])
    # plt.savefig('/home/cemaxim/CAS/CWTs/'+str(n)+'_CWT')
    # plot_orbit(data, ids, R_min, n ,time, plot=False, axs=ax[3])



