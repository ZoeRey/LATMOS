#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  7 15:38:42 2019

@author: rey
"""

import numpy as np
import os.path
import matplotlib.pyplot as plt


#%%Lecture des des données moyennées sur 10 min de interntemp, greentemp et IRtemp,
#Indiquer journée voulue : 20190430 ou 20190429

annee=["2019"]
dates=["20190429"]

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

interntemp=[]
greentemp=[]
IRtemp=[]

i=0
for a in num1:
    for b in num:
        for c in num2:
            for d in num3:
                for e in num4:
                    for f in num:
                        trajet='/net/nfs/prj1/QUALAIR/microLIDAR/database/MILAN/ascii/'+annee[0]+'/'+dates[0]+'/'+'mli2MOY_'+dates[0]+'.'+a+b+c+d+e+f+'.JUS'
                        if os.path.exists(trajet):
                            with open(trajet,'r') as fich:
                                for t in range(30):
                                    ligne=fich.readline()
                                mot_ligne=ligne.split()
                            interntemp.append(float(mot_ligne[1]))
                            with open(trajet,'r') as fich:
                                for u in range(31):
                                    ligne=fich.readline()
                                mot_ligne=ligne.split()
                            greentemp.append(float(mot_ligne[1]))
                            with open(trajet,'r') as fich:
                                for v in range(32):
                                    ligne=fich.readline()
                                mot_ligne=ligne.split()
                            IRtemp.append(float(mot_ligne[1]))
                            #print(i)
                            hours.append(a+b)
                            minutes.append(c+d)
                            secondes.append(e+f)

                            
 #%%
time=[]    
for i in range(len(hours)):
    time.append(float(hours[i])+((100./60.)*float(minutes[i])+float(secondes[i]))/100.)

                           
#%%Plot températures en fonction du temps
    
plt.plot(time,interntemp)
plt.xlabel("temps (UTC)")
plt.ylabel("inter_temp (C)")
plt.title(dates[0])
plt.legend()
plt.show()

plt.plot(time,greentemp)
plt.xlabel("temps (UTC)")
plt.ylabel("Green_temp (C)")
plt.title(dates[0])
plt.legend()
plt.show()

plt.plot(time,IRtemp)
plt.xlabel("temps (UTC)")
plt.ylabel("IR_temp (C)")
plt.title(dates[0])
plt.legend()
plt.show()
