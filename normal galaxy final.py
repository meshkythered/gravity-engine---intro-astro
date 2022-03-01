# -*- coding: utf-8 -*-
"""
Created on Tue Jan 11 18:31:59 2022

@author: ASUS
"""

import numpy as np
import matplotlib.pyplot as plt
import time
import copy
import pandas as pd
start_time = time.time()


n=500
T0= 1 #temprature
day=24*3600
year=365.25*day
ly=3*10**8*year
cutoff=200*ly
m=6.6*10**20*(10**13)*2 
k=4*10**-4
X=Y=Z=2000*ly*10
G=6.67*10**-11




def g_force (d):
    if d>cutoff:
        return (G*m**2/d**2)
    else:
        return (G*m**2/cutoff**2)

def g_sum_a (i,x,y,z):
    fx=0
    fy=0
    fz=0
    n=len(x)
    for j in range (0,n):
        if i!=j:
            dx= x[j]-x[i]
            dy= y[j]-y[i]
            dz= z[j]-z[i]
            d= np.sqrt (dx**2+dy**2+dz**2)
            f=g_force (d)
            fx+=f*dx/d
            fy+=f*dy/d
            fz+=f*dz/d
    return (fx/m,fy/m,fz/m)





def circle_galaxy(n, xc,yc,zc, sig, omega, spin):  #spin will be in forms like: [1,0,0] meaning its turning in yz plane
    n=100
    print (n)
    '''
    r=np.random.uniform (0,sig,n)
    theta= np.random.uniform (0,np.pi,n)
    phi= np.random.uniform (0,2*np.pi,n)
    zf=r*np.cos(theta)
    xf=r*np.sin(theta)*np.cos(phi)
    yf=r*np.sin(theta)*np.sin(phi)
    '''
    xf=[]
    yf=[]
    zf=[]
    i=0
    while i<n:
        xn=np.random.uniform (-sig,sig)
        yn=np.random.uniform (-sig,sig)
        zn=np.random.uniform (-sig,sig)
        if (xn**2+yn**2+zn**2*25)<sig**2:
            xf.append(xn)
            yf.append(yn)
            zf.append(zn)
            i+=1
    
    vxf=np.array(yf)*-omega
    vyf=np.array(xf)*omega
    vzf=np.zeros (n)
    xf=np.array(xf)+X/2
    yf=np.array(yf)+Y/2
    zf=np.array(zf)+Z/2
    axf=np.zeros (n)
    ayf=np.zeros (n)
    azf=np.zeros (n)
    for i in range (0,n):
        axf[i],ayf[i],azf[i]=g_sum_a(i,xf,yf,zf)
    
    #xfp,yfp,zfp,vxp,vyp,vzp=time_step(xf, yf, zf, vxf, vyf, vzf, 3*10**9*year, 10**7*year, 2)
    xfp,yfp,zfp,vxp,vyp,vzp=xf, yf, zf, vxf, vyf, vzf
    xf=xf-X/2
    yf=yf-Y/2
    zf=zf-Z/2
    if spin[0]!=0:
        s=spin[0]
        return (zfp*s+xc,xfp*s+yc,yfp*s+zc,vzp*s,vxp*s,vyp*s)
    elif spin[1]!=0:
        s=spin[1]
        return (yfp*s+xc,zfp*s+yc,xfp*s+zc,vyp*s,vzp*s,vxp*s)
    else:
        s=spin[2]
        return (xfp*s+xc,yfp*s+yc,zfp*s+zc,vxp*s,vyp*s,vzp*s)
    
    
        
        
    
    
# forming cluster starting conditions
'''
vmax=np.sqrt(3/4*k*T0/m)
for dx in np.arange(X/20,X/2+1,X/20):
    for dy in np.arange(Y/11,Y/11*10+1,Y/11):
        x.append(dx)
        y.append(dy)
#x=np.random.rand (n)*X
#y=np.random.rand (n)*Y
vx=vmax*(.5-np.random.rand (n))*2
vy=vmax*(.5-np.random.rand (n))*2
vz=vmax*(.5-np.random.rand (n))*2
T00= 2*m* np.average (vx**2+vy**2)/k
ax=np.zeros (n)
ay=np.zeros (n)
az=np.zeros (n)'''




def time_step (x,y,z,
               vx,vy,vz,
               T,h,show):
    n=len(x)
    Time=[]
    ax=np.zeros (n)
    ay=np.zeros (n)
    az=np.zeros (n)
    for i in range (0,n):
        ax[i],ay[i],az[i]=g_sum_a(i,x,y,z)
    for t in range (1, int(T/h)+1):
        if show!=0:
            plt.figure()
            plt.ylim(0,Y)
            plt.xlim (0,X)
            plt.scatter(y[0:int(n/2)],x[0:int(n/2)],c=(1,1,0))
            plt.scatter(y[int(n/2):n],x[int(n/2):n],c=(1,0,1))
            plt.scatter(y[0:int(n/2)],x[0:int(n/2)])
            plt.title('2 galaxy cluster t=%s * 10^6 year'%(t-1))
            plt.xlabel('x')
            plt.ylabel('y')
            if show==2:
                plt.savefig('2xy %s'%(t-1))
            plt.figure()
            plt.ylim(0,Y)
            plt.xlim (0,Z)
            plt.scatter(y,z)
            plt.title('2 galaxy cluster t=%s * 10^6 year'%(t-1))
            plt.xlabel('z')
            plt.ylabel('y')
            if show==2:
                plt.savefig('2zy %s'%(t-1))
            
        Time.append((t-1)*h)
        
        for i in range (0,n):
            dx=(vx[i]*h+ax[i]*h**2/2)
            x[i]=(x[i]+dx)%X
            y[i]=(y[i]+vy[i]*h+ay[i]*h**2/2)%Y
            z[i]=(z[i]+vz[i]*h+az[i]*h**2/2)%Z
            apx,apy,apz=g_sum_a(i,x,y,z)
            vx[i]+=(apx+ax[i])/2*h
            vy[i]+=(apy+ay[i])/2*h
            vz[i]+=(apz+az[i])/2*h
            ax[i]=apx
            ay[i]=apy
            az[i]=apz
        
    return (x,y,z,
            vx,vy,vz)

x1,y1,z1,vx1,vy1,vz1=circle_galaxy(n, X/3,Y+Y/15,Z/2, 2000*ly, 0, [0,0,1])
x2,y2,z2,vx2,vy2,vz2=circle_galaxy(n, 2*X/3,Y-Y/15,Z/2, 2000*ly, 0, [0,0,1])

vmax=2152
vx1=vx1 -vmax/5
vx2=vx2+ vmax/5
xf= np.append(x1,x2)
yf=np.append(y1,y2)
zf=np.append(z1,z2)
vxf=np.append(vx1,vx2)
vyf=np.append(vy1,vy2)
vzf=np.append(vz1,vz2)

xfp,yfp,zfp,vxp,vyp,vzp=time_step(xf, yf, zf, vxf, vyf, vzf, 10**10*year, 10**7*year*2,2)
print("--- %s seconds ---" % (time.time() - start_time)) #runtime
