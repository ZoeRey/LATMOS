#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 12 09:56:26 2019

@author: rey
"""

#Version Zoe
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 29 16:03:36 2019

Calcul du parametre geometrique et de la constante systeme

"""
import statistics
import jutilities as juti
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from matplotlib import pyplot

def func(x,a,b,c):
    f=c+(a*np.exp(-b*x))
    return f

#%% Read file


dates=["010000"]
trajet='/net/nfs/prj1/QUALAIR/microLIDAR/database/MILAN/ascii/2019/20190215/mli2MOY_20190215.'
file1=trajet+dates[0]+".JUS"
data=np.loadtxt(file1,skiprows=36)
altitude=data[:,0]
voie1=data[:,3]
var1=data[:,5]


#%% calculate mol 
source_mol='ISA'
rayleigh=juti.calc_molecular_profile(source_mol,altitude)

#%% plot figure

idx_3000m=np.where(altitude==3.)
id3=idx_3000m[0][0]
idx_9000m=np.where(altitude==9.0)
id9=idx_9000m[0][0]
idx_0m=np.where(altitude==0.)
id0=idx_0m[0][0]

#cste=1.
#cste=statistics.mean(voie1[id9:id3]/(rayleigh['beta_m'][id9:id3]*rayleigh['tm'][id9:id3]))


geom=[]
for i in range(len(altitude)):
    geom.append(voie1[i]/(rayleigh['beta_m'][i]*rayleigh['tm'][i]))

plt.figure()
plt.suptitle('02/15 1h00')
plt.subplot(211)
plt.plot(altitude,voie1,linewidth=1.)
plt.plot(altitude,rayleigh['beta_m']*rayleigh['tm'],linestyle='--',color='k')
plt.yscale('Log')
plt.xlim(0,3)
plt.xlabel("altitude")
plt.ylabel("log(voie)")

plt.subplot(212)
plt.plot(altitude,geom,linewidth=1.)
plt.xlim(0,3)
plt.ylim(0.1,1.5)
plt.xlabel("altitude")
plt.ylabel("facteur géométrique")
#plt.yscale('Log')


x=altitude[id3:id0]
y=geom[id3:id0]
popt,pcov=curve_fit(func,x,y)
geom_fit=func(x,popt[0],popt[1],popt[2])
#plt.plot(x,func(x,*popt))
plt.plot(x,geom[id3:id0],label='geom')
plt.plot(x,geom_fit,label='geom-fit')
plt.xlabel('altitude (km)')
plt.legend()


#plt.subplot(313)
#betaT2=voie1[id4:id3]*cste*cste_fit
#plt.plot(x,betaT2)
pyplot.gcf().subplots_adjust(left = 0.2, bottom = 0.2, right = 0.9, top = 0.9, wspace = 0.5, hspace = 0.5)
plt.show()

betaT2=voie1[id3:id0]/geom[id3:id0]
plt.plot(x,betaT2, label="signal corrigé")
plt.plot(x,voie1[id3:id0],label="signal brut")
plt.yscale('Log')
plt.legend()

#plt.plot(altitude,cste_fit)
#plt.plot(altitude,cste_fit,'o')



#%% save result
#geom=func(data[0]['altitude'],popt[0],popt[1])


