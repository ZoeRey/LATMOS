#Version Zoe
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 29 16:03:36 2019

Calcul du parametre geometrique et de la constante systeme pour le lidar Milan.

"""
#import statistics

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import pyplot

#from scipy.optimize import curve_fit


#%% Read file

annee=["2019"]
dates=["20190430"]
hour=["064918","065019","065119","065219","065319","065420","065520","065620","065720","065820","065920","070000"]
#hour=["053312"]

trajet='/net/nfs/prj1/QUALAIR/microLIDAR/database/MILAN/ascii/'+annee[0]+'/'+dates[0]+'/mli2_'+dates[0]+'.'+hour[0]+'.JUS'
file=np.loadtxt(trajet,skiprows=36)
altitude=file[:,0]
altitude=altitude[::-1]   

data532=[[0]*len(altitude) for i in range(len(hour))]

for i in range(len(hour)):
    trajet='/net/nfs/prj1/QUALAIR/microLIDAR/database/MILAN/ascii/'+annee[0]+'/'+dates[0]+'/mli2_'+dates[0]+'.'+hour[i]+'.JUS'
    file=np.loadtxt(trajet,skiprows=36) 
    voie532=file[:,1]
    data532[i]=voie532[::-1]
"""
plt.figure()
plt.title(dates[0]+","+hour[0])
plt.plot(data532[0][:600],altitude[:600])
plt.ylabel("altitude (km)")
plt.xlabel("voie1")
plt.legend()
plt.show()
"""
#%%Plot signal

plt.figure(figsize=(10, 8))
plt.subplot(221)
plt.title(dates[0]+","+hour[0])
plt.plot(data532[0][:600],altitude[:600])
plt.ylabel("altitude (km)")
plt.xlabel("voie1")

plt.subplot(222)
plt.title(dates[0]+","+hour[1])
plt.plot(data532[1][:600],altitude[:600])
plt.ylabel("altitude (km)")
plt.xlabel("voie1")

plt.subplot(223)
plt.title(dates[0]+","+hour[2])
plt.plot(data532[2][:600],altitude[:600])
plt.ylabel("altitude (km)")
plt.xlabel("voie1")

plt.subplot(224)
plt.title(dates[0]+","+hour[3])
plt.plot(data532[3][:600],altitude[:600])
plt.ylabel("altitude (km)")
plt.xlabel("voie1")

pyplot.tight_layout(pad=0.8, w_pad=1, h_pad=3.0)
plt.legend()
plt.show()


plt.figure(figsize=(10, 8))
plt.subplot(221)
plt.title(dates[0]+","+hour[4])
plt.plot(data532[4][:600],altitude[:600])
plt.ylabel("altitude (km)")
plt.xlabel("voie1")

plt.subplot(222)
plt.title(dates[0]+","+hour[5])
plt.plot(data532[5][:600],altitude[:600])
plt.ylabel("altitude (km)")
plt.xlabel("voie1")

plt.subplot(223)
plt.title(dates[0]+","+hour[6])
plt.plot(data532[6][:600],altitude[:600])
plt.ylabel("altitude (km)")
plt.xlabel("voie1")

plt.subplot(224)
plt.title(dates[0]+","+hour[7])
plt.plot(data532[7][:600],altitude[:600])
plt.ylabel("altitude (km)")
plt.xlabel("voie1")

pyplot.tight_layout(pad=0.8, w_pad=1, h_pad=3.0)
plt.legend()
plt.show()


plt.figure(figsize=(10, 8))
plt.subplot(221)
plt.title(dates[0]+","+hour[8])
plt.plot(data532[8][:600],altitude[:600])
plt.ylabel("altitude (km)")
plt.xlabel("voie1")

plt.subplot(222)
plt.title(dates[0]+","+hour[9])
plt.plot(data532[9][:600],altitude[:600])
plt.ylabel("altitude (km)")
plt.xlabel("voie1")

plt.subplot(223)
plt.title(dates[0]+","+hour[10])
plt.plot(data532[10][:600],altitude[:600])
plt.ylabel("altitude (km)")
plt.xlabel("voie1")

plt.subplot(224)
plt.title(dates[0]+","+hour[11])
plt.plot(data532[11][:600],altitude[:600])
plt.ylabel("altitude (km)")
plt.xlabel("voie1")

pyplot.tight_layout(pad=0.8, w_pad=1, h_pad=3.0)
plt.legend()
plt.show()


