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
T0= 10**6 #temprature
day=24*3600
year=365.25*day
ly=3*10**8*year
cutoff=2*ly
m=6.6*10**20
mi=9.27*10**-23 #mass of iron
k=1.4*10**-23
X=Y=Z=ly
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

def p_sum_a (i,x,y,z, p,ns,dx):
    fx=0
    fy=0
    fz=0
    n=len(x)
    for i in range (0,n):
        xx=int(x[i]/dx)
        yy=int(y[i]/dx)
        zz=int(z[i]/dx)
        aa=int(X/dx)
        fx+=p[(xx-1)%aa][yy][zz]-p[(xx+1)%aa][yy][zz]
        fy+=p[xx][(yy-1)%aa][zz]-p[xx][(yy+1)%aa][zz]
        fz+=p[xx][(yy-1)%aa][zz]-p[xx][(yy+1)%aa][zz]
        nss=ns[xx][yy][zz]
    return (fx/m/nss,fy/m/nss,fz/m/nss)

def Temp (vx,vy,vz):
    return(m/3*np.sum(vx**2+vy**2+vz**2))

def pressure (x,y,z,vx,vy,vz, dx):
    p=np.zeros ((int(X/dx),int(Y/dx),int(Z/dx)))
    ns=np.zeros ((int(X/dx),int(Y/dx),int(Z/dx)))
    for xx in range (0,int(X/dx)):
        for yy in range (0,int(Y/dx)):
            for zz in range (0,int(Z/dx)):
                n=0
                for i in range (0,n):
                    if x[i]>xx*dx and x[i]<xx*dx+dx:
                        if y[i]>yy*dx and y[i]<yy*dx+dx:
                            if z[i]>zz*dx and z[i]<zz*dx+dx:
                                p[xx][yy][zz]+=m*8.3/k/dx**3/3*(vx[i]**2+vy[i]**2+vz[i]**2)
                                n+=1
                    ns[xx][yy][zz]=n
    return(p,ns)



def supernova(n, xc,yc,zc, sig, T):  #spin will be in forms like: [1,0,0] meaning its turning in yz plane
    n=200
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
        if (xn**2+yn**2+zn**2)<sig**2:
            xf.append(xn)
            yf.append(yn)
            zf.append(zn)
            i+=1
    
    vxf=np.random.uniform(-(2*3.14*k*T0/mi/3)**.5, (2*3.14*k*T0/mi/3)**.5, n)
    vyf=np.random.uniform(-(2*3.14*k*T0/mi/3)**.5, (2*3.14*k*T0/mi/3)**.5, n)
    vzf=np.random.uniform(-(2*3.14*k*T0/mi/3)**.5, (2*3.14*k*T0/mi/3)**.5, n)
    xf=np.array(xf)+X/2
    yf=np.array(yf)+Y/2
    zf=np.array(zf)+Z/2
    axf=np.zeros (n)
    ayf=np.zeros (n)
    azf=np.zeros (n)
    for i in range (0,n):
        axf[i],ayf[i],azf[i]=g_sum_a(i,xf,yf,zf)
    
    xfp,yfp,zfp,vxp,vyp,vzp=time_step(xf, yf, zf, vxf, vyf, vzf, 10**10*5*year, 10**8*year, 2)
    #xfp,yfp,zfp,vxp,vyp,vzp=xf, yf, zf, vxf, vyf, vzf
    xf=xf-X/2
    yf=yf-Y/2
    zf=zf-Z/2
    return (zfp,xfp,yfp,vzp,vxp,vyp)
    
        
        
    
    
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
    dex=int (ly/10)
    p,ns=pressure (x,y,z,vx,vy,vz, dex)
    for i in range (0,n):
        apx,apy,apz=g_sum_a(i,x,y,z)
        ll,oo,pp=p_sum_a (i,x,y,z, p,ns,dex)
        apx+=ll
        apy+=oo
        apz+=pp
        ax[i],ay[i],az[i]=apx,apy,apz
    for t in range (1, int(T/h)+1):
        if show!=0:
            plt.figure()
            plt.ylim(0,Y)
            plt.xlim (0,X)
            plt.scatter(x, y)
            plt.title('t=%s year'%(t-1))
            
            if show==2:
                plt.savefig('%s'%(t-1))
            
        Time.append((t-1)*h)
        
        for i in range (0,n):
            dx=(vx[i]*h+ax[i]*h**2/2)
            x[i]=(x[i]+dx)%X
            y[i]=(y[i]+vy[i]*h+ay[i]*h**2/2)%Y
            z[i]=(z[i]+vz[i]*h+az[i]*h**2/2)%Z
            apx,apy,apz=g_sum_a(i,x,y,z)
            ll,oo,pp=p_sum_a (i,x,y,z, p,ns,dx)
            apx+=ll
            apy+=oo
            apz+=pp
            vx[i]+=(apx+ax[i])/2*h
            vy[i]+=(apy+ay[i])/2*h
            vz[i]+=(apz+az[i])/2*h
            ax[i]=apx
            ay[i]=apy
            az[i]=apz
        
    return (x,y,z,
            vx,vy,vz)

x1,y1,z1,vx1,vy1,vz1=supernova(n, X/3,Y,Z/2, ly/4,10**9)
#x2,y2,z2,vx2,vy2,vz2=circle_galaxy(n, 2*X/3,Y,Z/2, 2000*ly, 1/(250*10**6*year), [0,0,1])
'''
vx1=vx1+ 1000*ly*10*5
vx2=vx2- 1000*ly*10*5
xf= np.append(x1,x2)
yf=np.append(y1,y2)
zf=np.append(z1,z2)
vxf=np.append(vx1,vx2)
vyf=np.append(vy1,vy2)
vzf=np.append(vz1,vz2)
'''
#xfp,yfp,zfp,vxp,vyp,vzp=time_step(xf, yf, zf, vxf, vyf, vzf, 1000*10**6*year, 10**6*year, 1)
print("--- %s seconds ---" % (time.time() - start_time)) #runtime
