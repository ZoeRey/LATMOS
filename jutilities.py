#Version Zoe
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 28 14:57:49 2019
"""
import numpy as np

#calcul p1 et t1 en fonction de p0 et t0
def cal(p0, t0, a, h0, h1):
    g = 9.80665
    R = 287.00
    if a!=0.:
        t1 = t0 + a * (h1 - h0)
        p1= p0 * (t1 / t0) ** (-g / a / R)
    else:
        t1 = t0
        p1 = p0 * np.exp(-g / R / t0 * (h1 - h0))
    return t1, p1

#Calcul le profil de température et pression
def isa(altitude_m,season):
    if season=='standard':
        t0=288.15
        a=[-0.0065,0.,0.001]
        h=[11000,20100,32100]
    elif season=='summer':
        t0=287
        a=[-0.0053,-0.007,0.,0.0014]
        h=[4700,10000,23000,31800]
    elif season=='winter':
        t0=257.1
        a=[0.003,-0.0032,-0.0068,0.,-0.0006,0.001]
        h=[1000,3200,8500,15500,25000,30000]
    else:
        print('incorrect season')
        return
    
    p0 = 101325
    prevh=0.
    isa={}
    if altitude_m < 0 or altitude_m > h[len(h)-1]:
        print("altitude must be in [0," +str(h[len(h)-1])+"]")
        return
    for i in range(len(h)):
        if altitude_m <= h[i]:
            temperature, pressure = cal(p0, t0, a[i], prevh, altitude_m)
            break;
        else:
            t0, p0 = cal(p0, t0, a[i], prevh, h[i])
            prevh = h[i]
    
    isa['t']=temperature
    isa['p']=pressure
    return isa

#calcul le coefficient d'extinction
def calc_apha(longueur_onde,p,t):
    p0 = p[len(p)-1] #Attention, doit être p[0] et t[0] si l'altitude est rangée dans l'autre sens
    t0 = t[len(p)-1]   
    alpha=1.66*((550/longueur_onde)**(4.09))*(p/p0)*(t0/t)*(10**(-5))
    return alpha

#calcul le coefficient de retrodiffusion
def calc_beta(alpha):
    beta=3/(8*np.pi)*alpha
    return beta

#besoin        
def calc_molecular_profile(source_mol,altitude,longonde):
    rayleigh={'t':[],'p':[],'alpha_m':[],'beta_m':[],'mod':[],'tm':[],'lambda': longonde}    
    for alt in range(len(altitude)):
        tmp=isa(altitude[alt]*1000,'standard')
        rayleigh['t'].append(tmp['t'])
        rayleigh['p'].append(tmp['p'])
    rayleigh['alpha_m']=calc_apha(rayleigh['lambda'],rayleigh['p'],rayleigh['t'])
    rayleigh['beta_m']=np.array(calc_beta(rayleigh['alpha_m']))
    rayleigh['mod']=[np.trapz(rayleigh['alpha_m'][:i],x=altitude[:i]*1000) for i in range(len(altitude))]
    rayleigh['tm']=[np.exp(-2*rayleigh['mod'][i]) for i in range(len(altitude))]
    
    
    return rayleigh
