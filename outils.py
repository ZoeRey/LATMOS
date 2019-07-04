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

def calc_cste(data,datamoy,rayleigh,indice3,indice4):
    cste=statistics.mean(datamoy[indice3:indice4]/(rayleigh['beta_m'][indice3:indice4]*rayleigh['tm'][indice3:indice4]))
    
    return cste

#%%Calcul fonction de recouvrement
    
def calc_geom(data,cste,rayleigh,altitude):
    geom=[[0]*len(data[0]) for i in range(0,2)]
    for n in range(0,2):
        for i in range(len(data[0])):
            geom[n][i]=(data[n][i]/(cste[n]*rayleigh['beta_m'][i]*rayleigh['tm'][i]))
    for n in range(0,2):
        for i in range(len(altitude)-len(data[0])):
            geom[n].append(1)    
    return geom

def calc_geom_fit_moy(geom,altitude,indice2):
    geom_fit=[[0]*len(altitude) for i in range(0,2)]
    for n in range(0,2):
        x=altitude[:indice2]
        y=geom[n][:indice2]
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

def calc_cste2(data,datamoy,rayleigh,indice1,indice4,trans_aer):
    cste=statistics.mean(datamoy[indice1:indice4]/(rayleigh['beta_m'][indice1:indice4]*rayleigh['tm'][indice1:indice4]*trans_aer*trans_aer))
    return cste

#%%Calcul transmission aérosol à zmax
    
def calc_transmission_aer(LRp,LRm,altitude,datacorr,rayleigh):
    idx_3500=np.where(altitude==3.495)
    id35=idx_3500[0][0]

    imax=datacorr.index(max(datacorr[4:id35]))
    diff1=10000.
    for i in range(imax,id35):
        diff=abs((datacorr[imax]-rayleigh['beta_m'][imax])/2+rayleigh['beta_m'][imax]-datacorr[i])
        if(diff<diff1):
            diff1=diff
            itop=i

    EO_p=LRp*altitude[itop]*(datacorr[imax]/rayleigh['tm'][imax]-rayleigh['beta_m'][imax])
    T_p=np.exp(-1*EO_p)

    print("imax=",imax)
    print("altitude(imax)=",altitude[imax])
    print("altitude(itop)=",altitude[itop])
    print("AOD à zmax",EO_p)
    print("Transmission à zmax=",T_p)
    print("Rapport diffusion=",1.+(datacorr[imax]/rayleigh['beta_m'][imax]))
    print("\n")
    
    return T_p,imax

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

#%%
def calc_beta_aer_stable2(izr,datamoy_zr,altitude,rayleigh,datamoy,LRp,LRm,cste):
    """
    a=[]
    b=[]
    beta_p=[]
    for i in range(0,izr):
        a.append(np.trapz(rayleigh['beta_m'][i:izr],x=altitude[i:izr]))
        b.append(datamoy[i]*np.exp(2.*(LRp-LRm)*a[i]))
    c=[np.trapz(b[i:izr],x=altitude[i:izr]) for i in range(0,izr)]
    for i in range(0,izr):
        beta_p.append((datamoy[i]*np.exp(2.*(LRp-LRm)*a[i]))/(cste*rayleigh['tm'][izr]+2.*LRp*c[i])-rayleigh['beta_m'][i])
    return beta_p
    """
    
    altitude2=altitude[:izr]
    altitude3=altitude2[::-1]
    ray=rayleigh['beta_m'][:izr]
    ray2=ray[::-1]
    datamoy2=datamoy[:izr]
    datamoy3=datamoy2[::-1]

    a=[]
    b=[]
    beta_p=[]
    for i in range(len(ray)):
        a.append(np.trapz(ray2[:i],x=altitude3[:i]))
        b.append(datamoy3[i]*np.exp(-2.*(LRp-LRm)*a[i]))
    #print(a)
    c=[np.trapz(b[:i],x=altitude3[:i]) for i in range(len(ray))]
    for i in range(len(ray)):
        beta_p.append((datamoy3[i]*np.exp(-2.*(LRp-LRm)*a[i]))/(cste*rayleigh['tm'][izr]-2.*LRp*c[i])-ray2[i])
    beta_p=beta_p[::-1]
    return beta_p
    
        
    
    


#%%calcul alpha et transmission
    
def calc_alpha_AOD_trans(beta,LR,altitude):
    alpha=[]
    trans=[]
    for i in range(len(beta)):
        alpha.append(beta[i]*LR)
    AOD=[np.trapz(alpha[:i],x=altitude[:i]) for i in range(len(altitude))]
    for j in range(len(AOD)):
        trans.append(np.exp(-1*AOD[j]))
        
    return alpha,AOD,trans
    


