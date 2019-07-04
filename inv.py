#Version Zoe
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 29 16:03:36 2019

Calcul du parametre geometrique et de la constante systeme pour le lidar Milan.

"""
import statistics
import jutilities as juti
import matplotlib.pyplot as plt
import numpy as np
import outils
#from scipy.optimize import curve_fit
import scipy.signal 


#%% Read file

#les deux premiers fichiers sont utilisés pour calculer une moyenne de la fonction de recouvremenent (couche limite très peu développée)
#les autres sont utilisés pour le calcul de la constante instrumentale
#annee=["2019","2018","2018","2018","2018","2018","2018"]
annee=["2019","2018","2019","2019","2019","2019","2019"]
dates=["20190115","20190115","20190115","20190115","20190115","20190115","20190115"]
hour=["110000","111003","112006","113009","114012","115015","120000"]


trajet='/net/nfs/prj1/QUALAIR/microLIDAR/database/MILAN/ascii/'+annee[2]+'/'+dates[2]+'/mli2MOY_'+dates[2]+'.'+hour[2]+'.JUS'
file=np.loadtxt(trajet,skiprows=36)
altitude=file[:,0]
altitude=altitude[::-1]   

data808=[[0]*len(altitude) for i in range(len(dates))]
data532=[[0]*len(altitude) for i in range(len(dates))]
data808_dp=[[0]*len(altitude) for i in range(len(dates))]
data532_dp=[[0]*len(altitude) for i in range(len(dates))]

for i in range(len(dates)):
    trajet='/net/nfs/prj1/QUALAIR/microLIDAR/database/MILAN/ascii/'+annee[i]+'/'+dates[i]+'/mli2MOY_'+dates[i]+'.'+hour[i]+'.JUS'
    file=np.loadtxt(trajet,skiprows=36) 
    voie808=file[:,3]
    voie532=file[:,1]
    voie808_dp=file[:,4]
    voie532_dp=file[:,2]
    data808[i]=voie808[::-1]
    data532[i]=voie532[::-1]
    data808_dp[i]=voie808_dp[::-1]
    data532_dp[i]=voie532_dp[::-1]

data532moy=[0]*len(altitude)
data808moy=[0]*len(altitude)
data532moy_dp=[0]*len(altitude)
data808moy_dp=[0]*len(altitude)

for i in range(len(altitude)):
    for n in range(len(dates)):
        data808moy[i]=data808moy[i]+data808[n][i]
        data532moy[i]=data532moy[i]+data532[n][i]
        data808moy_dp[i]=data808moy_dp[i]+data808_dp[n][i]
        data532moy_dp[i]=data532moy_dp[i]+data532_dp[n][i]
    data808moy[i]=data808moy[i]/(len(dates))
    data532moy[i]=data532moy[i]/(len(dates)-2)
    data808moy_dp[i]=data808moy_dp[i]/(len(dates))
    data532moy_dp[i]=data532moy_dp[i]/(len(dates))

#%% dépo

depo532=[]
depo808=[]
for i in range(len(data532moy)):
    depo532.append(data532moy_dp[i]/data532moy[i])
    depo808.append(data808moy_dp[i]/data808moy[i])
    
plt.figure()
#plt.title('Dépolarisation')
plt.plot(depo532[:200],altitude[:200])
plt.xlabel('ratio de dépolarisation 532nm')
plt.ylabel('altitude(km)')
plt.xlim(0,0.1)
plt.show()

plt.figure()
#plt.title('Dépolarisation')
plt.plot(depo808[:200],altitude[:200])
plt.xlabel('ratio de dépolarisation 808np')
plt.ylabel('altitude(km)')
plt.xlim(0,0.1)
plt.show()


#%% calcul profil moleculaire et rayleigh
source_mol='ISA'
rayleigh808=juti.calc_molecular_profile(source_mol,altitude,808)
rayleigh532=juti.calc_molecular_profile(source_mol,altitude,532)

#%% limite de bruit étendue
z0=0.08

data1=[]
data2=[]
data3=[]
data4=[]
for i in range(len(altitude)):
    data1.append(data532moy[i]/(altitude[i]-z0)**2)
    data2.append(data532moy_dp[i]/(altitude[i]-z0)**2)
    data3.append(data808moy[i]/(altitude[i]-z0)**2)
    data4.append(data808moy_dp[i]/(altitude[i]-z0)**2)

bruit532=[]
bruit532_dp=[]
bruit808=[]
bruit808_dp=[]

ectype1=np.std(data1[1300:])
ectype2=np.std(data2[1300:])
ectype3=np.std(data3[1300:])
ectype4=np.std(data4[1300:])


for i in range(len(altitude)):
    bruit532.append(2.*ectype1*(altitude[i]-z0)*(altitude[i]-z0))
    bruit532_dp.append(2.*ectype2*(altitude[i]-z0)*(altitude[i]-z0))
    bruit808.append(2.*ectype3*(altitude[i]-z0)*(altitude[i]-z0))
    bruit808_dp.append(2.*ectype4*(altitude[i]-z0)*(altitude[i]-z0))
    
#%% calcul constante
idx_3000m=np.where(altitude==6)
id1=idx_3000m[0][0]
idx_4000m=np.where(altitude==8.01)
id2=idx_4000m[0][0]

idx_8000m=np.where(altitude==6)
id8=idx_8000m[0][0]
idx_9000m=np.where(altitude==7.005)
id9=idx_9000m[0][0]

idx_3000m=np.where(altitude==3)
id3=idx_3000m[0][0]
idx_4000m=np.where(altitude==4.005)
id4=idx_4000m[0][0]

print("808")
cste808=outils.calc_cste(data808,data808moy,rayleigh808,id8,id9)
print("532")
cste532=outils.calc_cste(data532,data532moy,rayleigh532,id8,id9)
print("808_dp")
cste808_dp=outils.calc_cste(data808,data808moy_dp,rayleigh808,id3,id4)
print("532_dp")
cste532_dp=outils.calc_cste(data532,data532moy_dp,rayleigh532,id3,id4)

plt.figure()
plt.title('Vérification superpostition courbes, '+dates[2]+' '+hour[2]+', 808')
plt.plot(data808moy,altitude,linewidth=1.,label="PR2")
plt.plot(bruit808,altitude,label="limite de bruit")
plt.plot(cste808*rayleigh808['beta_m']*rayleigh808['tm'],altitude,linestyle='--',color='k',label="Rayleigh*constante")
plt.ylim(0,10)
plt.xlim(0,100)
plt.ylabel("altitude (km)")
plt.xlabel("voie 3")
plt.grid()
plt.legend()
plt.show()

plt.figure()
plt.title('Vérification superpostition courbes, '+dates[2]+' '+hour[2]+', 808dp')
plt.plot(data808moy_dp,altitude,linewidth=1.,label="PR2")
plt.plot(bruit808_dp,altitude,label="limite de bruit")
plt.plot(cste808_dp*rayleigh808['beta_m']*rayleigh808['tm'],altitude,linestyle='--',color='k',label="Rayleigh*constante")
plt.ylim(0,10)
plt.xlim(-5,5)
plt.ylabel("altitude (km)")
plt.xlabel("voie 4")
plt.grid()
plt.legend()
plt.show()

plt.figure()
plt.title('Vérification superpostition courbes, '+dates[2]+' '+hour[2]+', 532')
plt.plot(data532moy,altitude,linewidth=1.,label="PR2")
plt.plot(bruit532,altitude,label="limite de bruit")
plt.plot(cste532*rayleigh532['beta_m']*rayleigh532['tm'],altitude,linestyle='--',color='k',label="Rayleigh*constante")
plt.ylim(0,4)
plt.xlim(0,50)
plt.ylabel("altitude (km)")
plt.xlabel("voie 1")
plt.grid()
plt.legend()
plt.show()

plt.figure()
plt.title('Vérification superpostition courbes, '+dates[2]+' '+hour[2]+', 532_dp')
plt.plot(data532moy_dp,altitude,linewidth=1.,label="PR2")
plt.plot(bruit532_dp,altitude,label="limite de bruit")
plt.plot(cste532_dp*rayleigh532['beta_m']*rayleigh532['tm'],altitude,linestyle='--',color='k',label="Rayleigh*constante")
plt.ylim(0,10)
plt.xlim(-5,5)
plt.ylabel("altitude (km)")
plt.xlabel("voie 2")
plt.grid()
plt.legend()
plt.show()


#%%Calcul fonction de recouvrement

trajet_geom='/net/nfs/prj1/QUALAIR/microLIDAR/database/MILAN/fgeom/milan_20190121_fgeom_v6.txt'
file=np.loadtxt(trajet_geom,skiprows=13)
fgeom532=file[:,1]
fgeom808=file[:,3]

geom532_fit=list(fgeom532[::-1])
geom808_fit=list(fgeom808[::-1])

for i in range(len(altitude)-len(geom532_fit)):
    geom532_fit.append(1)
    geom808_fit.append(1)

"""
geom808=outils.calc_geom(data808,cste808,rayleigh808,altitude)
geom532=outils.calc_geom(data532,cste532,rayleigh532,altitude)

geom808_fit=outils.calc_geom_fit_moy(geom808,altitude,id4)
geom532_fit=outils.calc_geom_fit_moy(geom532,altitude,id4)
"""

plt.figure()
plt.xlabel("altitude (km)")
plt.ylabel("fonction de recouvrement")
plt.plot(altitude,geom808_fit,label="808nm")
plt.plot(altitude,geom532_fit,label="532nm")
plt.xlim(0,2)
plt.legend()
plt.show()


#%%Calcul data corrigé de constante

data808cor=[]
data532cor=[]
data532cor_dp=[]
data808cor_dp=[]
for i in range(len(altitude)):
    data808cor.append(data808moy[i]/cste808)
    data808cor_dp.append(data808moy_dp[i]/cste808_dp)
    data532cor.append(data532moy[i]/cste532)
    data532cor_dp.append(data532moy_dp[i]/cste532_dp)

#%%

"""
plt.figure()
plt.title('Vérification superpostition courbes, '+dates[2]+' '+hour[2]+', 532')
plt.plot(data532cor[:50],altitude[:50],linewidth=1.,label="PR2")
#plt.plot(bruit532,altitude,label="limite de bruit")
plt.plot(rayleigh532['beta_m'][:50]*rayleigh532['tm'][:50],altitude[:50],linestyle='--',color='k',label="Rayleigh*constante")
plt.ylabel("altitude (km)")
plt.xlabel("voie 1")
plt.xscale('log')
plt.grid()
plt.legend()
plt.show()
"""

plt.figure()
plt.title('Vérification superpostition courbes, '+dates[2]+' '+hour[2]+', 532')
plt.plot(data532cor[:80],altitude[:80],linewidth=1.,label="PR2")
#plt.plot(bruit808,altitude,label="limite de bruit")
plt.plot(rayleigh532['beta_m'][5:80]*rayleigh532['tm'][5:80],altitude[5:80],linestyle='--',color='k',label="Rayleigh*constante")
plt.ylabel("altitude (km)")
plt.xlabel("voie 3")
plt.grid()
plt.legend()
plt.show()

#%%
"""
approx=[]
x=altitude[8:10]
y=data532cor[8:10]
coeff_poly=np.polyfit(x,y,1)
for i in range(0,20):
    approx.append(coeff_poly[0]*altitude[i]+coeff_poly[1])


plt.figure()
plt.title('Vérification superpostition courbes, '+dates[2]+' '+hour[2]+', 532')
plt.plot(data532cor[3:80],altitude[3:80],linewidth=1.,label="PR2")
plt.plot(approx,altitude[0:20],label="apro")
plt.ylabel("altitude (km)")
plt.xlabel("voie 3")
plt.grid()
plt.legend()
plt.show()

c=[]
for i in range(0,9):
    #print(altitude[i])
    print(data532cor[i]/approx[i])
    c.append(data532cor[i]/approx[i])

datacor2=[]
for i in range(len(c)):
    datacor2.append(data532cor[i]/c[i])
    
for i in range(len(data532cor)-len(datacor2)):
    #print(i)
    datacor2.append(data532cor[i+9])

plt.figure()
plt.title('Vérification superpostition courbes, '+dates[2]+' '+hour[2]+', 532')
plt.plot(datacor2[:200],altitude[:200],linewidth=1.,label="PR2")
plt.ylabel("altitude (km)")
plt.xlabel("voie 3")
plt.grid()
plt.legend()
plt.show()

k=[]
a=geom532_fit[:9]
for i in range(len(c)):
    k.append(a[i]*c[i])
k=k[::-1]
for i in range(len(k)):
    print(k[i])
"""
#%%Estimation transmission aérosols à zmax

LRp=40.
LRm=(8*np.pi)/3.

T_p532,imax532=outils.calc_transmission_aer(LRp,LRm,altitude,data532cor,rayleigh532)
T_p808,imax808=outils.calc_transmission_aer(LRp,LRm,altitude,data808cor,rayleigh808)
T_p532_dp,imax532=outils.calc_transmission_aer(LRp,LRm,altitude,data532cor_dp,rayleigh532)
T_p808_dp,imax808=outils.calc_transmission_aer(LRp,LRm,altitude,data808cor_dp,rayleigh808)


#%%Nouvelle estimation de la constante instrumentale

cste808=outils.calc_cste2(data808,data808moy,rayleigh808,id8,id9,T_p808)
cste532=outils.calc_cste2(data532,data532moy,rayleigh532,id8,id9,T_p532)
cste808_dp=outils.calc_cste2(data808,data808moy_dp,rayleigh808,id8,id9,T_p808_dp)
cste532_dp=outils.calc_cste2(data532,data532moy_dp,rayleigh532,id8,id9,T_p532_dp)


#%%Zone aveuglement

i_avg=15

altitude1=altitude
altitude=altitude[i_avg:]
data808moy1=data808moy
data532moy=data532moy[i_avg:]
data808moy=data808moy[i_avg:]
data532moy_dp=data532moy_dp[i_avg:]
data808moy_dp=data808moy_dp[i_avg:]

rayleigh532['beta_m']=rayleigh532['beta_m'][i_avg:]
rayleigh532['tm']=rayleigh532['tm'][i_avg:]
rayleigh808['beta_m']=rayleigh808['beta_m'][i_avg:]
rayleigh808['tm']=rayleigh808['tm'][i_avg:]
geom532_fit=geom532_fit[i_avg:]
geom808_fit=geom808_fit[i_avg:]
                        

#%%Nouveau calcul de data corrigé

data808cor=[]
data532cor=[]
data532cor_dp=[]
data808cor_dp=[]
for i in range(len(altitude)):
    data808cor.append(data808moy[i]/cste808/geom808_fit[i])
    data808cor_dp.append(data808moy_dp[i]/cste808_dp/geom808_fit[i])
    data532cor.append(data532moy[i]/cste532/geom532_fit[i])
    data532cor_dp.append(data532moy_dp[i]/cste532_dp/geom532_fit[i])
    
    
plt.figure()
plt.title('Vérification superpostition courbes, date+heure, 808')
plt.plot(data808cor[:300],altitude[:300],linewidth=1.,label="PR2 corrigé")
plt.plot(rayleigh808['beta_m'][:300]*rayleigh808['tm'][:300],altitude[:300],linestyle='--',color='k',label="Rayleigh")
plt.ylabel("altitude (km)")
plt.xlabel("voie3")
plt.grid()
plt.legend()
plt.show()

plt.figure()
plt.title('Vérification superpostition courbes, date+heure, 532dp')
plt.plot(data532cor_dp[:300],altitude[:300],linewidth=1.,label="PR2 corrigé")
plt.plot(rayleigh532['beta_m'][:300]*rayleigh532['tm'][:300],altitude[:300],linestyle='--',color='k',label="Rayleigh")
plt.ylabel("altitude (km)")
plt.xlabel("voie2")
plt.grid()
plt.legend()
plt.show()

plt.figure()
plt.title('Vérification superpostition courbes, date+heure, 532')
plt.plot(data532cor[:300],altitude[:300],linewidth=1.,label="PR2 corrigé")
plt.plot(rayleigh532['beta_m'][:300]*rayleigh532['tm'][:300],altitude[:300],linestyle='--',color='k',label="Rayleigh")
plt.ylabel("altitude (km)")
plt.xlabel("voie1")

plt.grid()
plt.legend()
plt.show()

#%%Calcul inversion


izrkm=np.where(altitude==5.01)
izr=izrkm[0][0]

data808moy_zr=statistics.mean(data808moy[izr-10:izr+10])
data532moy_zr=statistics.mean(data532moy[izr-10:izr+10])
data808moy_dp_zr=statistics.mean(data808moy_dp[izr-10:izr+10])
data532moy_dp_zr=statistics.mean(data532moy_dp[izr-10:izr+10])

beta808_aer=outils.calc_beta_aer_stable2(izr,data808moy_zr,altitude,rayleigh808,data808moy,LRp,LRm,cste808[2])
beta532_aer=outils.calc_beta_aer_stable2(izr,data532moy_zr,altitude,rayleigh532,data532moy,LRp,LRm,cste532[2])

beta808_aer_dp=outils.calc_beta_aer_stable2(izr,data808moy_dp_zr,altitude,rayleigh808,data808moy,LRp,LRm,cste808_dp[2])
beta532_aer_dp=outils.calc_beta_aer_stable2(izr,data532moy_dp_zr,altitude,rayleigh532,data532moy,LRp,LRm,cste532_dp[2])

#%%Calcul alpha et transmission

alpha532_aer,AOD532,trans532_aer=outils.calc_alpha_AOD_trans(beta532_aer,LRp,altitude[:izr])
alpha808_aer,AOD808,trans808_aer=outils.calc_alpha_AOD_trans(beta808_aer,LRp,altitude[:izr])

alpha532_aer_dp,AOD532_dp,trans532_aer_dp=outils.calc_alpha_AOD_trans(beta532_aer_dp,LRp,altitude[:izr])
alpha808_aer_dp,AOD808_dp,trans808_aer_dp=outils.calc_alpha_AOD_trans(beta808_aer_dp,LRp,altitude[:izr])

#%%plot transmission
"""
plt.figure()
plt.title("transmission")
plt.plot(trans808_aer[:izr],altitude[:izr])
plt.grid()
plt.xlabel("transmission 808// aer")
plt.ylabel("altitude (km)")
plt.legend()
plt.show()
"""
#%%plt profil rétroddifusion

plt.figure()
plt.title("Profil de rétrodiffusion "+dates[2]+' '+hour[2])
plt.plot(beta808_aer[7:izr],altitude[7:izr])
#plt.plot(rayleigh808['beta_m'][6:izr],altitude[6:izr],label="beta808_mol")
plt.grid()
plt.xlabel("beta 808// (km⁻¹.sr⁻¹)")
plt.ylabel("altitude (km)")
plt.legend()
plt.show()


plt.figure()
plt.title("Profil de rétrodiffusion "+dates[2]+' '+hour[2])
plt.plot(beta532_aer[7:izr],altitude[7:izr])
#plt.plot(rayleigh808['beta_m'][6:izr],altitude[6:izr],label="beta808_mol")
plt.grid()
plt.xlabel("beta 532// (km⁻¹.sr⁻¹)")
plt.ylabel("altitude (km)")
plt.legend()
plt.show()
"""
plt.figure()
plt.title("Profil de rétrodiffusion "+dates[2]+' '+hour[2])
plt.plot(beta532_aer_dp[7:izr],altitude[7:izr])
#plt.plot(rayleigh808['beta_m'][6:izr],altitude[6:izr],label="beta808_mol")
plt.grid()
plt.xlabel("beta 532_dp (km⁻¹.sr⁻¹)")
plt.ylabel("altitude (km)")
plt.legend()
plt.show()

plt.figure()
plt.title("Profil de rétrodiffusion "+dates[2]+' '+hour[2])
plt.plot(beta808_aer_dp[7:izr],altitude[7:izr])
#plt.plot(rayleigh808['beta_m'][6:izr],altitude[6:izr],label="beta808_mol")
plt.grid()
plt.xlabel("beta 808_dp (km⁻¹.sr⁻¹)")
plt.ylabel("altitude (km)")
plt.legend()
plt.show()
"""
#%% plot épaiseur optique
plt.figure()
plt.title("Epaisseur optique")
plt.plot(AOD808[10:izr],altitude[10:izr])
plt.grid()
plt.xlabel("AOD 808//")
plt.ylabel("altitude (km)")
plt.legend()
plt.show()

#%%plot PR2 brut
plt.figure()
plt.title("PR2")
plt.plot(data808moy[:izr],altitude[:izr])
plt.grid()
plt.xlabel("PR2 808// aer")
plt.ylabel("altitude (km)")
plt.legend()
plt.show()

#%% plot indice de couleur

chi=[]
for i in range(len(beta532_aer)):
    chi.append(beta808_aer[i]/beta532_aer[i])


plt.figure()
plt.title("Indice de couleur "+dates[2]+' '+hour[2])
plt.plot(chi[:izr],altitude[:izr])
plt.grid()
plt.xlabel("indice de couleur 808/532")
plt.ylabel("altitude (km)")
#plt.xlim(0,10)
plt.legend()
plt.show()

#%%Coeff angstrom

A=[]
for i in range(len(beta532_aer)):
    A.append(np.log(AOD532[i]/AOD808[i])/np.log(0.808/0.532))
    
plt.figure()
plt.title("Coeff Angtröm "+dates[2]+' '+hour[2])
plt.plot(A[:izr],altitude[:izr])
plt.grid()
plt.xlabel("coefficient d'Angström")
plt.ylabel("altitude (km)")
plt.legend()
plt.show()

print("A=",A[len(A)-1])

#%%Rapport de diffusion

R532=[]
R808=[]
for i in range(len(beta532_aer)):
    R532.append(1.+beta532_aer[i]/rayleigh532['beta_m'][i])
    R808.append(1.+beta808_aer[i]/rayleigh808['beta_m'][i])
    

R808_mdfilter=scipy.signal.medfilt(R808,15)
l1=[]
for i in range(0,400):
    l1.append(1)



plt.figure()
plt.title("Rapport de diffusion "+dates[2]+' '+hour[2])
plt.plot(R532[:izr],altitude[:izr])
plt.grid()
plt.xlabel("Rapport de diffusion 532// nm")
plt.ylabel("altitude (km)")
plt.legend()
plt.show()

plt.figure()
plt.title("Rapport de diffusion "+dates[2]+' '+hour[2])
plt.plot(R808[:izr],altitude[:izr])
plt.plot(l1[:280],altitude[:280],color='black',linestyle='--')
#plt.grid()
plt.xlabel("Rapport de diffusion 808// nm")
plt.ylabel("altitude (km)")
plt.legend()
plt.show()


plt.figure()
plt.title("Rapport de diffusion "+dates[2]+' '+hour[2])
plt.plot(R808_mdfilter[:],altitude[:len(R808_mdfilter)])
plt.plot(l1[:izr],altitude[:izr],color='black',linestyle='--')
plt.ylim(0,3)
#plt.grid()
plt.xlabel("Rapport de diffusion 808// nm")
plt.ylabel("altitude (km)")
plt.legend()
plt.show()

#%%Ratio de dépolarisation




