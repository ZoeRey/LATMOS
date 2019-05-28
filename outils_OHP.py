#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 09:55:35 2019

@author: rey
"""
import numpy as np
import statistics
from scipy.optimize import curve_fit
#import matplotlib.pyplot as plt

#%%Fonction pour approcher la fonction de recouvrement

def func(x,a,b,c):
    f=c+(a*np.exp(-b*x))
    return f

#%%Calcul constante première itération

def calc_cste(datamoy,rayleigh,indice3,indice4):
    cste=statistics.mean(datamoy[indice3:indice4]/(rayleigh['beta_m'][indice3:indice4]*rayleigh['tm'][indice3:indice4]))   
    if cste<0:
        cste=-1*cste
    print(cste)
    return cste

#%%Calcul fonction de recouvrement
    
def calc_geom(data,cste,rayleigh,altitude):
    geom=[]
    for i in range(len(altitude)):
            geom.append((data[i]/(cste*rayleigh['beta_m'][i]*rayleigh['tm'][i])))
    return geom

def calc_geom_fit_moy(geom,altitude):
    #geom_fit_start=geom[:14]
    geom_fit=[]
    x=altitude[:200]
    y=geom[:200]
    popt,pcov=curve_fit(func,x,y)
    geom_fit=list(func(x,popt[0],popt[1],popt[2]))
    geom_fit2=geom_fit
    for i in range(len(geom_fit2)):
        if geom_fit2[i]<0:
            geom_fit2[i]=geom[i]
        if geom_fit2[i]>1 :
            geom_fit2[i]=1.
    geom_fit2[0]=0.000001
    for i in range(len(altitude)-len(geom_fit2)):
        geom_fit2.append(1)
        
    
    return geom_fit2

#%%Calcul constante seconde itération

def calc_cste2(datamoy,rayleigh,indice1,indice2,trans_aer):
    cste=statistics.mean(datamoy[indice1:indice2]/(rayleigh['beta_m'][indice1:indice2]*rayleigh['tm'][indice1:indice2]*trans_aer*trans_aer))
    if cste<0:
        cste=-1*cste    
    print(cste)
    return cste

#%%Calcul transmission aérosol à zmax
    
def calc_transmission_aer(LRp,altitude,datacorr,rayleigh):
    idx_3500=np.where(altitude==3.495)
    id35=idx_3500[0][0]

    imax=datacorr.index(max(datacorr[10:id35]))
    diff1=100.
    for i in range(imax,id35):
        diff=abs((datacorr[imax]-rayleigh['beta_m'][imax])/2+rayleigh['beta_m'][imax]-datacorr[i])
        if(diff<diff1):
            diff1=diff
            itop=i
    """
    beta_p_zmax=(datacorr[imax]/(rayleigh['tm'][imax]*0.975*0.975)*altitude[itop]*1000)
    beta_m_zmax=np.trapz(rayleigh['beta_m'][:itop],x=altitude[:itop]*1000)
    beta=[]
    for i in range(len(altitude)):
        beta.append(datacorr[i]/rayleigh['tm'][i])
        
    plt.figure()
    plt.plot(altitude[:id35],beta[:id35])
    plt.plot(altitude[:id35],rayleigh['beta_m'][:id35])
    plt.grid()
    plt.show()
    EO_p_zmax=beta_p_zmax-beta_m_zmax
    T_p=np.exp(-1.*LRp*EO_p_zmax)
    """
    
    EO_p=LRp*altitude[itop]*(datacorr[imax]/rayleigh['tm'][imax]-rayleigh['beta_m'][imax])
    T_p=np.exp(-1*EO_p)
    
    print("imax=",imax)
    print("itop=",itop)
    print("altitude[imax]=",altitude[imax])
    print("altitude[itop]",altitude[itop])
    #print(beta_p_zmax)
    #print(beta_m_zmax)
    print("Transmission=",T_p)
    print("Epaisseur optique=",EO_p)
    
    return T_p
    

#%%Calcul transmission aérosol à zmax
    
def calc_transmission_aer2(LRp,altitude,datacorr,rayleigh):
    R=[]
    for i in range(len(altitude)):
        R.append(datacorr[i]/(rayleigh['beta_m'][i]*rayleigh['tm'][i]))
    idx_3500=np.where(altitude==4.5)
    id35=idx_3500[0][0]

    imax=datacorr.index(max(datacorr[:id35]))
    diff1=100.
    for i in range(imax,id35):
        diff=abs((datacorr[imax]-rayleigh['beta_m'][imax])/6+rayleigh['beta_m'][imax]-datacorr[i])
        if(diff<diff1):
            diff1=diff
            itop=i
    T2=R[itop]/R[imax]
    
    print(imax)
    print(itop)
    print(T2**0.5)
    AOD=-1.*np.log(T2)/2.
    print(AOD)
    return T2**0.5

#%%calcul coefficient rétrodiffusion aérosols
    
def calc_beta_aer(rayleigh,altitude,datamoy,LRp,LRm,cste):
    a=[np.trapz(rayleigh['beta_m'][:i],x=altitude[:i]) for i in range(len(altitude))]
    b=[]
    for i in range(len(altitude)):
        b.append(datamoy[i]*np.exp(-2.*(LRp-LRm)*a[i]))
    c=[np.trapz(b[:i],x=altitude[:i]) for i in range(len(altitude))]
    beta_p=[]
    for i in range(len(altitude)):
        beta_p.append((datamoy[i]*np.exp(-2.*(LRp-LRm)*a[i]))/(cste-2.*LRp*c[i])-rayleigh['beta_m'][i])
    return beta_p
"""
Test pour avoir le fit
    beta_p_fit=beta_p[:id8]
    coeff_beta808_p_fit=np.polyfit(altitude[id8:],beta808_p[id8:],1)
    for i in range(len(altitude[id8:])):
        beta808_p_fit.append(coeff_beta808_p_fit[0]*altitude[i]+coeff_beta808_p_fit[1])
"""    
#%% calcul coefficient rétrodiffusion aérosols stabilité ++

def calc_beta_aer_stable(izr,datamoy_zr,altitude,rayleigh,datamoy,LRp,LRm):
    
    a=[]
    b=[]
    beta_p=[]
    for i in range(0,izr):
        a.append(np.trapz(rayleigh['beta_m'][i:izr],x=altitude[i:izr]))
        b.append(datamoy[i]*np.exp(2.*(LRp-LRm)*a[i]))
    c=[np.trapz(b[i:izr],x=altitude[i:izr]) for i in range(0,izr)]
    for i in range(0,izr):
        beta_p.append((datamoy[i]*np.exp(2.*(LRp-LRm)*a[i]))/(datamoy_zr/rayleigh['beta_m'][izr]+2.*LRp*c[i])-rayleigh['beta_m'][i])


    for i in range(0,len(altitude)-izr):
        a.append(np.trapz(rayleigh['beta_m'][izr-1:izr-1+i],x=altitude[izr-1:izr-1+i]))
        b.append(datamoy[izr-1+i]*np.exp(-2.*(LRp-LRm)*a[izr-1+i]))
    c=c+[np.trapz(b[izr-1:i],x=altitude[izr-1:i]) for i in range(0,len(altitude)-izr)]
    for i in range(0,len(altitude)-izr):
        beta_p.append((datamoy[izr-1+i]*np.exp(-2.*(LRp-LRm)*a[izr-1+i]))/(datamoy_zr/rayleigh['beta_m'][izr]-2.*LRp*c[izr-1+i])-rayleigh['beta_m'][izr-1+i])

    return beta_p



















