
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
cutoff=30*ly
m=6.6*10**20*(10**13)
k=4*10**-4
X=Y=Z=2000*ly*10
G=6.67*10**-11




def g_force (d,cut):
    if d>cut:
        return (G/d**2)
    else:
        return (G/cut**2)

def g_sum_a (i,x,y,z,m_arr):
    fx=0
    fy=0
    fz=0
    n=len(x)
    for j in range (0,n):
        if i!=j:
            if m_arr[j]!=1:
                dx= x[j]-x[i]
                dy= y[j]-y[i]
                dz= z[j]-z[i]
                d= np.sqrt (dx**2+dy**2+dz**2)
                f=g_force (d, (m_arr[i]+m_arr[j])*cutoff)
                f=f*m_arr[i]*m_arr[j]*m**2
                fx+=f*dx/d
                fy+=f*dy/d
                fz+=f*dz/d
    return (fx/m,fy/m,fz/m)





def circle_galaxy(n, xc,yc,zc, sig, omega, spin):  #spin will be in forms like: [1,0,0] meaning its turning in yz plane
    n=50
    m_arr=np.zeros(n)+1
    m_arr[0]=sig**3*omega**2/m/G
    r=np.random.uniform (100*cutoff,sig*.9,n)
    theta= np.random.uniform (0,np.pi,n)
    phi= np.random.uniform (0,2*np.pi,n)
    zf=r*np.cos(theta)
    xf=r*np.sin(theta)*np.cos(phi)
    yf=r*np.sin(theta)*np.sin(phi)
    xf[0]=0
    yf[0]=0
    zf[0]=0
    '''
    xf=[]
    yf=[]
    zf=[]
    i=0
    while i<n:
        xn=np.random.uniform (-sig,sig)
        yn=np.random.uniform (-sig,sig)
        zn=np.random.uniform (-sig,sig)
        if (xn**2+yn**2+zn**2)<sig**2:
            xf.append(xn)
            yf.append(yn)
            zf.append(zn)
            i+=1
    '''
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
        axf[i],ayf[i],azf[i]=g_sum_a(i,xf,yf,zf,m_arr)
    
    xfp,yfp,zfp,vxp,vyp,vzp=time_step(xf, yf, zf, vxf, vyf, vzf, axf, ayf, azf,m_arr, 1000*10**6*year*2, 10**6*year, 1)
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
               ax,ay,az,
               m_arr,
               T,h,show):
    Time=[]
    n=len(x)
    for t in range (1, int(T/h)+1):
        if show!=0:
            if (t-1)%(int(T/h)/10)==0:
                plt.figure()
                plt.ylim(0,Y)
                plt.xlim (0,X)
                plt.scatter(x, y)
                plt.title('t=%s light year'%(t-1))
                if show==2:
                    plt.savefig('%s'%(t-1))
            
        Time.append((t-1)*h)
        
        for i in range (0,n):
            x[i]=(x[i]+vx[i]*h+ax[i]*h**2/2)%X
            y[i]=(y[i]+vy[i]*h+ay[i]*h**2/2)%Y
            z[i]=(z[i]+vz[i]*h+az[i]*h**2/2)%Z
            apx,apy,apz=g_sum_a(i,x,y,z,m_arr)
            vx[i]+=(apx+ax[i])/2*h
            vy[i]+=(apy+ay[i])/2*h
            vz[i]+=(apz+az[i])/2*h
            ax[i]=apx
            ay[i]=apy
            az[i]=apz
        
    return (x,y,z,
            vx,vy,vz)

x,y,z,vx,vy,vz=circle_galaxy(n, X/2,Y/2,Z/2, 2000*ly, 1/(250*10**6*year), [0,0,1])

print("--- %s seconds ---" % (time.time() - start_time)) #runtime
