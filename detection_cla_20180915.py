#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  7 15:38:42 2019

@author: rey
"""
import numpy as np
import os.path
import matplotlib.pyplot as plt
import scipy.signal

#%%Lecture des fichiers de la journée

annee=["2018"]
dates=["20180915"]

if annee[0]=="2018" :
    trajet='/net/nfs/prj1/QUALAIR/microLIDAR/database/MILAN/ascii/2018/'+dates[0]+'/mli2MOY_'+dates[0]+'.000000.JUS'
else :
    trajet='/net/nfs/prj1/QUALAIR/microLIDAR/database/MILAN/ascii/2019/'+dates[0]+'/mli2MOY_'+dates[0]+'.000000.JUS'

file=np.loadtxt(trajet,skiprows=36)
altitude=file[:,0]
altitude=altitude[::-1]  
#zone aveugelement jusqu'à 75m
#altitude=altitude[6:]

num1=["0","1","2"]
num2=["0","1","2","3","4","5"]
num3=["0"]
num4=["0","1","2"]
num=["0","1","2","3","4","5","6","7","8","9"]
hours=[]
minutes=[]
secondes=[]


data532paral=[[0]*len(altitude) for k in range(0,143)]
i=0
for a in num1:
    for b in num:
        for c in num2:
            for d in num3:
                for e in num4:
                    for f in num:
                        trajet='/net/nfs/prj1/QUALAIR/microLIDAR/database/MILAN/ascii/'+annee[0]+'/'+dates[0]+'/'+'mli2MOY_'+dates[0]+'.'+a+b+c+d+e+f+'.JUS'
                        if os.path.exists(trajet):
                            #print(i)
                            hours.append(a+b)
                            minutes.append(c+d)
                            secondes.append(e+f)
                            file=np.loadtxt(trajet,skiprows=36)
                            col1=file[:,1]
                            col1reverse=col1[::-1]                    
                            data532paral[i]=col1reverse[:]
                            i=i+1

#%%
time=[]    
for i in range(len(data532paral)):
    time.append(float(hours[i])+((100./60.)*float(minutes[i])+float(secondes[i]))/100.)

#%%read fonction recouvrement 532 + data corrigé

trajet_geom='/net/nfs/home/rey/Documents/code_Zoe/MILAN/fct_geom532_JUS.txt'
geom532=list(np.loadtxt(trajet_geom))
geom532[0]=1.

if annee[0]=="2018" :
    for i in range(len(altitude)-len(geom532)):
        geom532.append(1.)
        

#%% recherche du max

iymax=[]
"""
for i in range(len(time_moy)):
    for j in range(len(altitude)):
        data532moy[i][j]=data532moy[i][j]/geom532[j]
"""        
        
plt.figure()
plt.plot(data532paral[55][:100],altitude[:100])
plt.grid()
plt.show()

for i in range(54,90): #de 9h à 16h
    if i<=60: #de 9h à 10h
        y=list(data532paral[i][:130])
        iymax.append(y.index(max(y)))
    if i>60 and i <78: #de 10h à 13h
        y=list(data532paral[i][:130])
        iymax.append(y.index(max(y)))
    if i>=78: #de 13h à 16h
        y=list(data532paral[i][:130])
        iymax.append(y.index(max(y)))

#%%filtre médian
    
time_hcla=time[54:90]
data532_medianfilter=[[] for i in range(len(time_hcla))]

for i in range(len(iymax)):
    if i<=5:
        data532_medianfilter[i]=data532_medianfilter[i]+list(scipy.signal.medfilt(data532paral[i+16][10:iymax[i]+50],3))
    else:
        data532_medianfilter[i]=data532_medianfilter[i]+list(scipy.signal.medfilt(data532paral[i+16][10:iymax[i]+80],3))
"""
#quelle heure tracer ? k=
k=8
plt.figure()
plt.plot(data532moy[k*2][iymax[k*2-16]-6:iymax[k*2-16]+60],altitude[iymax[k*2-16]-6:iymax[k*2-16]+60])
plt.plot(data532_medianfilter[k*2-16],altitude[iymax[k*2-16]-6:iymax[k*2-16]+60])
plt.grid()
plt.show()
"""
#%% Approximation polynôme de degré 5 et calcul dérivé 2nd de ce polynome

datafit=[[] for i in range(len(time_hcla))]
datafitderiv2=[[] for i in range(len(time_hcla))]
datafitderiv1=[[] for i in range(len(time_hcla))]

for i in range(len(time_hcla)):
    print(i)
    if i<=5:
        x=altitude[10:iymax[i]+50]
    else :
        x=altitude[10:iymax[i]+80]
    y=data532_medianfilter[i]
    coeff_poly=np.polyfit(x,y,5)

    for j in range(len(x)):
        datafit[i].append(coeff_poly[0]*x[j]*x[j]*x[j]*x[j]*x[j]+coeff_poly[1]*x[j]*x[j]*x[j]*x[j]+coeff_poly[2]*x[j]*x[j]*x[j]+coeff_poly[3]*x[j]*x[j]+coeff_poly[4]*x[j]+coeff_poly[5])
        datafitderiv1[i].append(5*coeff_poly[0]*x[j]*x[j]*x[j]*x[j]+4*coeff_poly[1]*x[j]*x[j]*x[j]+3*coeff_poly[2]*x[j]*x[j]+2*coeff_poly[3]*x[j]+coeff_poly[4])
        datafitderiv2[i].append(20*coeff_poly[0]*x[j]*x[j]*x[j]+12*coeff_poly[1]*x[j]*x[j]+6*coeff_poly[2]*x[j]+2*coeff_poly[3])
        
    plt.figure()
    plt.plot(datafit[i],x,color='black',linewidth=2.5,linestyle='--',label="approximation")
    #plt.plot(data532paral[i+49][iymax[i]-4:iymax4[i]+12],altitude[iymax[i]-4:iymax4[i]+12],label="S 532//")
    plt.plot(y,x)
    plt.xlabel("signal rétrodiffusé, 532// nm")
    plt.ylabel("altitude (km)")
    plt.legend()
    plt.grid()
    plt.show()
    

#%% Recherche point d'inflexion
#il peut y avoir plusieurs points d'inflexion

pt_inflexion=[[] for i in range(len(time_hcla))]

for i in range(len(time_hcla)):
    if i<=5:
        x=altitude[:40]
    else :
        x=altitude[20:iymax[i]+80]
    for k in range(len(x)-1):
        if(datafitderiv2[i][k+1]>=0 and datafitderiv2[i][k]<=0) or (datafitderiv2[i][k+1]<=0 and datafitderiv2[i][k]>=0):
            pt_inflexion[i].append(x[k])
            
#%%Vérif s'il y a des endroits sans points d'inflexion
                
for i in range(len(pt_inflexion)):
    if len(pt_inflexion[i])==0:
        pt_inflexion[i]=pt_inflexion[i-1]       

           
#%%
hcla=[]
#poser la question !!!!
print("A "+str(time[90])+" le sommet de la couche limite était plutôt à :")
print(pt_inflexion[len(pt_inflexion)-1])
hcla_15h=input("Choisissez: ")
hcla_15h=int(hcla_15h)
hcla.append(pt_inflexion[len(pt_inflexion)-1][hcla_15h])
for i in range(len(pt_inflexion)-1,0,-1):
    diff1=100.
    for j in range(len(pt_inflexion[i])):
        diff=abs(pt_inflexion[i][j]-hcla[abs(i-len(pt_inflexion)+1)])
        if diff<diff1:
            diff1=diff
            jgood=j
    hcla.append(pt_inflexion[i][jgood])  
hcla=hcla[::-1]

plt.figure()
plt.plot(time_hcla,hcla)
plt.xlabel("time")
plt.ylabel("hcla")
plt.grid()
plt.show()


#%% test en partant du début

hcla2=[]
#poser la question !!!!
print("A "+str(time[54])+", le sommet de la couche limite était plutôt à :")
print(pt_inflexion[0])
hcla_8h=input("Choisissez: ")
hcla_8h=int(hcla_8h)
hcla2.append(pt_inflexion[0][hcla_8h])
for i in range(1,len(pt_inflexion)):
    diff1=100.
    for j in range(len(pt_inflexion[i])):
        diff=abs(pt_inflexion[i][j]-hcla2[i-1])
        if diff<diff1:
            diff1=diff
            jgood=j
    hcla2.append(pt_inflexion[i][jgood])  



plt.figure()
plt.plot(time_hcla,hcla2)
plt.xlabel("time (h)")
plt.ylabel("hcla (km)")
plt.grid()
plt.show()
