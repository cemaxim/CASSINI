# -*- coding: utf-8 -*-
"""
Created on Wed Oct  1 13:48:55 2025

@author: maxim
"""
import numpy as np
import matplotlib.pyplot as plt
import random


def synth_wave():
    x = np.linspace(0,60,1000)
    Ts = np.linspace(40,60,20)/60
    PPO = np.sin(2*np.pi*(1/10.7)*x)
    
    ys = []
    for i in Ts:
        y = np.sin(2*np.pi*(1/i)*x+random.uniform(0,2*np.pi))
        ys.append(y)
    y = np.sum(ys,axis=0)+PPO
    plt.figure()
    plt.plot(x,y)
    
        
    
    

        