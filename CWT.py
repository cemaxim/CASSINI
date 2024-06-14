#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 11 15:01:30 2024

@author: cemaxim
"""
from scipy import signal
import pycwt as wavelet
from scipy import integrate
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


time = ['2006-04-28 20:23:31', '2006-05-22 09:59:31']
#t = np.linspace(0,20,100)
#bx = np.sin(2*np.pi*t*60)
df = pd.read_pickle('/home/cemaxim/pickles/MAG_DF.pkl')
df = df[time[0]:time[1]]
df = df.iloc[::15,:]
bx = df.BTHETA.values


def get_cwt(bx):
 
    t = np.arange(0,len(bx))
    
    omega0=6
    print(t)
    dt = t[1]-t[0]
    Fs = 1/dt
    N = len(bx)
    f = np.arange(0,N/2+1)*Fs/N
    """  
    m3 = 1/((9.92/3)*60.*60.)
    m2 = 1/((9.92/2)*60.*60.)
    m1 = 1/(9.92*60.*60.)
    ten_hr = m1
    fifteen_hr = 1/(15*60.*60.)
    thirty_hr = 1/(30*60.*60.) 
    two_day = 1/(48*60.*60.)
    sixty_hr = 1/(60*60.*60.)
    """
    widths = omega0*Fs/(2*np.pi*f[1:])  #units of time
    cwtmatr = signal.cwt(bx, signal.morlet2, widths, w = omega0)  #continuous wavelet transform matrix
    #cwtmatr, freqs = wavelet.cwt(bx,widths,"cmor2.5-1.0", sampling_period = dt)
    
    ha = plt.subplot(111)
    print(np.size(t), np.size(f),np.shape(cwtmatr))
    
    ha.pcolormesh(t,f,abs(cwtmatr[:,:-1])**2, cmap='turbo')
    #ha.contourf(t,f,abs(cwtmatr[:,:-1])**2, cmap='jet')
    """
    ha.plot([t.min(),t.max()],[ten_hr,ten_hr],color='white')
    ha.plot([t.min(),t.max()],[sixty_hr/3, sixty_hr/3],color='white')
    ha.plot([t.min(),t.max()],[thirty_hr,thirty_hr],color='white')
    ha.plot([t.min(),t.max()],[sixty_hr,sixty_hr],color='white')
    """
    ha.set_yscale('log')
    ha.set_ylim(f[1],f[-1])
    ha.set_xlabel('Time')
    ha.set_ylabel('Frequency [Hz]')
   
    """
   
    mother = wavelet.Morlet(6)
    s0 = 1*dt
    dj = 1/12
    J = 7/dj
    alpha, _, _ = wavelet.ar1(bx)  # Lag-1 autocorrelation for red noise
    
    wave, scales, freqs, coi, fft, fftfreqs = wavelet.cwt(s0 = s0, signal=bx, dt=dt, wavelet=mother)
    
    power = np.abs(wave)**2
    fft_power = np.abs(fft)**2
    period = (1/freqs)/60 #hours
    
    signif, fft_theor = wavelet.significance(1.0, dt, scales, 0, alpha,
                                         significance_level=0.95,
                                         wavelet=mother)
    sig95 = np.ones([1, N]) * signif[:, None]
    sig95 = power / sig95
    
    plt.figure()
    ha = plt.subplot(111)
    ha.pcolormesh(t,period,power, cmap='turbo')
    ha.contour(t,period,sig95)
    ha.plot([t.min(),t.max()],[10,10],color='white')
    ha.plot([t.min(),t.max()],[15,15],color='white')
    ha.plot([t.min(),t.max()],[30,30],color='white')
    #ha.pcolormesh(t,period/3600,coi,cmap='jet')
    ha.set_yscale('log')
    ha.set_ylim(period.min(),period.max())
    ha.set_xlabel('Time (m)')
    ha.set_ylabel('Period (hr)')
    
    ha.fill(np.concatenate([t, t[-1:] + dt, t[-1:] + dt,
                           t[:1] - dt, t[:1] - dt]),
        np.concatenate([coi/3600, [1e-1], period[-1:],
                           period[-1:], [1e-1]]),
        'k', alpha=0.3, hatch='x')
    
    plt.show()

    psd = np.zeros(len(widths))
    for i in range(len(widths)):  #sum over time shifts such that PSD is a function of tau (width)
        psd[i] = (2/N)*integrate.trapz(np.abs(cwtmatr[i, :])**2, range(0, len(t)))
    
    plt.figure()
    plt.loglog(f[1:], psd)
    
    plt.plot([m1,m1],[psd.min(),psd.max()])
    plt.plot([m2,m2],[psd.min(),psd.max()])
    plt.plot([m3,m3],[psd.min(),psd.max()])   
    plt.plot([sixty_hr,sixty_hr],[psd.min(),psd.max()])
    plt.plot([thirty_hr,thirty_hr],[psd.min(),psd.max()])
    plt.plot([3*sixty_hr,3*sixty_hr],[psd.min(),psd.max()])
    
    plt.xlabel('freq (Hz)')
    plt.ylabel('PSD')
    plt.title('$\omega_0$ = '+str(omega0))
    #plt.loglog([50,50],[np.min(psd),np.max(psd)],':')
    #plt.loglog([120,120],[np.min(psd),np.max(psd)],':')
    plt.show()
    
    """