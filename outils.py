#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 09:55:35 2019

@author: rey
"""
import numpy as np
import statistics
from scipy.optimize import curve_fit

#%%Fonction pour approcher la fonction de recouvrement

def func(x,a,b,c):
    f=c+(a*np.exp(-b*x))
    return f

#%%Calcul constante première itération

def calc_cste(data,datamoy,rayleigh,indice1,indice2,indice3,indice4):
    cste=[0]*3
    cste[0]=statistics.mean(data[0][indice1:indice2]/(rayleigh['beta_m'][indice1:indice2]*rayleigh['tm'][indice1:indice2]))
    cste[1]=statistics.mean(data[1][indice1:indice2]/(rayleigh['beta_m'][indice1:indice2]*rayleigh['tm'][indice1:indice2]))
    cste[2]=statistics.mean(datamoy[indice3:indice4]/(rayleigh['beta_m'][indice3:indice4]*rayleigh['tm'][indice3:indice4]))
    
    for i in range(0,3):
        if cste[i]<0 :
            cste[i]=-cste[i]
    
    return cste

#%%Calcul fonction de recouvrement
    
def calc_geom(data,cste,rayleigh,altitude):
    geom=[[0]*len(altitude) for i in range(0,2)]
    for n in range(0,2):
        for i in range(len(altitude)):
            geom[n][i]=(data[n][i]/(cste[n]*rayleigh['beta_m'][i]*rayleigh['tm'][i]))
    return geom

def calc_geom_fit_moy(geom,altitude,indice1,indice2):
    geom_fit=[[0]*len(altitude) for i in range(0,2)]
    for n in range(0,2):
        x=altitude[indice1:indice2]
        y=geom[n][indice1:indice2]
        popt,pcov=curve_fit(func,x,y)
        geom_fit[n]=func(x,popt[0],popt[1],popt[2])

    geom_fit_moy=[]
    for i in range(len(geom_fit[0])):
        geom_fit_moy.append((geom_fit[0][i]+geom_fit[1][i])/2.) 
        
        
    geom_fit_tot=geom_fit_moy
    for i in range(len(altitude)-len(geom_fit_moy)):
        geom_fit_tot.append(1)
        
    return geom_fit_tot

#%%Calcul constante seconde itération

def calc_cste2(data,datamoy,rayleigh,indice1,indice2,indice3,indice4,trans_aer):
    cste=[0]*3
    cste[0]=statistics.mean(data[0][indice1:indice2]/(rayleigh['beta_m'][indice1:indice2]*rayleigh['tm'][indice1:indice2]*trans_aer*trans_aer))
    cste[1]=statistics.mean(data[1][indice1:indice2]/(rayleigh['beta_m'][indice1:indice2]*rayleigh['tm'][indice1:indice2]*trans_aer*trans_aer))
    cste[2]=statistics.mean(datamoy[indice3:indice4]/(rayleigh['beta_m'][indice3:indice4]*rayleigh['tm'][indice3:indice4]*trans_aer*trans_aer))
    
    for i in range(0,3):
        if cste[i]<0 :
            cste[i]=-cste[i]
    
    return cste

#%%Calcul transmission aérosol à zmax
    
def calc_transmission_aer(LRp,LRm,altitude,datamoy,datacorr,rayleigh):
    idx_3500=np.where(altitude==3.495)
    id35=idx_3500[0][0]

    imax=datamoy.index(max(datamoy[:id35]))

    beta_tot_zmax=datacorr[imax]
    beta_m_zmax=rayleigh['beta_m'][imax]
    beta_p_zmax=(beta_tot_zmax-beta_m_zmax)*altitude[imax]*1000
    T_p=np.exp(-1.*LRp*beta_p_zmax)
    
    return T_p

#%%calcul coefficient rétrodiffusion aérosols
    
def calc_beta_aer(rayleigh,altitude,datamoy,LRp,LRm,cste):
    a=[np.trapz(rayleigh['beta_m'][:i],x=altitude[:i]*1000) for i in range(len(altitude))]
    b=[]
    for i in range(len(altitude)):
        b.append(datamoy[i]*np.exp(-2.*(LRp-LRm)*a[i]))
    c=[np.trapz(b[:i],x=altitude[:i]*1000) for i in range(len(altitude))]
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
        a.append(np.trapz(rayleigh['beta_m'][i:izr],x=altitude[i:izr]*1000))
        b.append(datamoy[i]*np.exp(2.*(LRp-LRm)*a[i]))
    c=[np.trapz(b[i:izr],x=altitude[i:izr]*1000) for i in range(0,izr)]
    for i in range(0,izr):
        beta_p.append((datamoy[i]*np.exp(2.*(LRp-LRm)*a[i]))/(datamoy_zr/rayleigh['beta_m'][izr]+2.*LRp*c[i])-rayleigh['beta_m'][i])


    for i in range(0,len(altitude)-izr):
        a.append(np.trapz(rayleigh['beta_m'][izr-1:izr-1+i],x=altitude[izr-1:izr-1+i]*1000))
        b.append(datamoy[izr-1+i]*np.exp(-2.*(LRp-LRm)*a[izr-1+i]))
    c=c+[np.trapz(b[izr-1:i],x=altitude[izr-1:i]*1000) for i in range(0,len(altitude)-izr)]
    for i in range(0,len(altitude)-izr):
        beta_p.append((datamoy[izr-1+i]*np.exp(-2.*(LRp-LRm)*a[izr-1+i]))/(datamoy_zr/rayleigh['beta_m'][izr]-2.*LRp*c[izr-1+i])-rayleigh['beta_m'][izr-1+i])

    return beta_p



















