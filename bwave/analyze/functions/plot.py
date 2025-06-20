#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  3 15:54:47 2020

@author: lloverasdan
"""

import numpy as np
from wrf import getvar, smooth2d, dbz, interplevel
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
from scipy.integrate import simps

CP = 1004
TR = 270.

#-------------------- PLAN VIEWS --------------------#

### Surface temperature and pressure
def surface(ncfile, x1, x2, y1, y2, dx, time=0, roll=0, save_file=None, width=9., height=6., filetype='PDF', \
            hline=None, hlinemin=0, hlinemax=0, hliney=0, vline=None, vlinemin=0, vlinemax=0, vlinex=0, title='Surface'):
    
    matplotlib.use(filetype)
    p = getvar(ncfile,'pressure',timeidx=time)[0,:,:]
    p_surf = np.asarray(p)
    p_surf = np.roll(p_surf,roll,axis=1)
    #p_surf = smooth2d(p_surf,150)
    p_surf = p_surf[y1:y2,x1:x2]
    tk = getvar(ncfile,'tk',timeidx=time)
    t_surf = np.asarray(tk)[0,:,:]-273.15
    t_surf = np.roll(t_surf,roll,axis=1)
    #t_surf = smooth2d(t_surf,150)
    t_surf = t_surf[y1:y2,x1:x2]
    fig = plt.figure(figsize=(width,height))
    ax = fig.add_subplot(1,1,1)
    ax.set_aspect('equal')
    plt.contourf(t_surf,levels=np.arange(-25,30,5),cmap=get_cmap('coolwarm'),extend='both')
    #plt.contourf(t_surf,levels=np.arange(-30,32.5,2.5),cmap=get_cmap('bwr'),extend='both')
    cbar = plt.colorbar(fraction=0.024, pad=0.04)
    cbar.ax.tick_params(labelsize=14) 
    cs = ax.contour(p_surf,levels=np.arange(800,1200,4),colors='k',linewidths=1.0)
    if np.amax(p_surf) - np.amin(p_surf) > 20:
        ax.clabel(cs,fmt='%1.0f',inline=1,levels=np.arange(800,1200,8),fontsize=10,colors='k')
    if hline:
        plt.hlines(hliney,hlinemin,hlinemax,colors='orange',linewidths=2.0)
    if vline:
        plt.vlines(vlinex,vlinemin,vlinemax,colors='orange',linewidths=2.0)        
    yax1 = np.arange(y1*dx,y2*dx,1000).astype(int)
    yax2 = (yax1 - y1*dx)/float(dx)
    yax2 = yax2.astype(int)
    xax1 = np.arange(x1*dx,x2*dx,1000).astype(int)
    xax2 = (xax1 - x1*dx)/float(dx)
    xax2 = xax2.astype(int)
    plt.tick_params(axis='both',labelsize=14)
    plt.xlabel('West-East (km)',fontsize=18,labelpad=6)
    plt.xticks(xax2,xax1)
    plt.ylabel('South-North (km)',fontsize=18,labelpad=6)
    plt.yticks(yax2,yax1)
    plt.title(title,fontsize=18)
    if save_file is not None:
        plt.savefig(save_file,dpi=300,bbox_inches='tight')
        
    return fig

### Surface wind
def surf_wind(ncfile, x1, x2, y1, y2, dx, time=0, roll=0, save_file=None, width=9., height=6., filetype='PDF', \
            hline=None, hlinemin=0, hlinemax=0, hliney=0, vline=None, vlinemin=0, vlinemax=0, vlinex=0, title='Surface'):
    
    matplotlib.use(filetype)
    p = getvar(ncfile,'pressure',timeidx=time)[0,:,:]
    p_surf = np.asarray(p)
    p_surf = np.roll(p_surf,roll,axis=1)
    #p_surf = smooth2d(p_surf,150)
    p_surf = p_surf[y1:y2,x1:x2]
    spd = getvar(ncfile,'wspd10',timeidx=time)
    spd = np.roll(spd,roll,axis=1)
    spd = spd[y1:y2,x1:x2]
    fig = plt.figure(figsize=(width,height))
    ax = fig.add_subplot(1,1,1)
    ax.set_aspect('equal')
    plt.contourf(spd,levels=np.arange(2,30,2),cmap=get_cmap('rainbow'),extend='max')
    cbar = plt.colorbar(fraction=0.024, pad=0.04)
    cbar.ax.tick_params(labelsize=14) 
    cs = ax.contour(p_surf,levels=np.arange(800,1200,4),colors='k',linewidths=1.0)
    if np.amax(p_surf) - np.amin(p_surf) > 20:
        ax.clabel(cs,fmt='%1.0f',inline=1,levels=np.arange(800,1200,8),fontsize=10,colors='k')
    if hline:
        plt.hlines(hliney,hlinemin,hlinemax,colors='orange',linewidths=2.0)
    if vline:
        plt.vlines(vlinex,vlinemin,vlinemax,colors='orange',linewidths=2.0)        
    yax1 = np.arange(y1*dx,y2*dx,1000).astype(int)
    yax2 = (yax1 - y1*dx)/float(dx)
    yax2 = yax2.astype(int)
    xax1 = np.arange(x1*dx,x2*dx,1000).astype(int)
    xax2 = (xax1 - x1*dx)/float(dx)
    xax2 = xax2.astype(int)
    plt.tick_params(axis='both',labelsize=14)
    plt.xlabel('West-East (km)',fontsize=18,labelpad=6)
    plt.xticks(xax2,xax1)
    plt.ylabel('South-North (km)',fontsize=18,labelpad=6)
    plt.yticks(yax2,yax1)
    plt.title(title,fontsize=18)
    if save_file is not None:
        plt.savefig(save_file,dpi=300,bbox_inches='tight')
        
    return fig

### Cloud top temperature
def cloud_tops(ncfile, x1, x2, y1, y2, dx, time=0, roll=0, save_file=None, width=9., height=6., filetype='PDF',\
               hline=None, hlinemin=0, hlinemax=0, hliney=0, vline=None, vlinemin=0, vlinemax=0, vlinex=0):

    matplotlib.use(filetype)
    p = getvar(ncfile,'slp',timeidx=time)
    p_surf = np.asarray(p)
    p_surf = np.roll(p_surf,roll,axis=1)
    #p_surf = smooth2d(p_surf,150)
    p_surf = p_surf[y1:y2,x1:x2]
    ctt = getvar(ncfile,'ctt',timeidx=time,fill_nocloud=True)
    ctt = np.asarray(ctt)
    ctt = np.roll(ctt,roll,axis=1)[y1:y2,x1:x2]
    fig = plt.figure(figsize=(width,height))
    ax = fig.add_subplot(1,1,1)
    ax.set_aspect('equal')
    ax.set_facecolor('lightblue')
    plt.contourf(ctt,levels=np.arange(-80,-10,10),cmap='Greys',extend='min')
    cbar = plt.colorbar(fraction=0.024, pad=0.04)
    cbar.ax.tick_params(labelsize=14)
    cs = ax.contour(p_surf,levels=np.arange(800,1200,10),colors='k',linewidths=2.0,alpha=0.5)
    if np.amax(p_surf) - np.amin(p_surf) > 20:
        ax.clabel(cs,fmt='%1.0f',inline=1,levels=np.arange(800,1200,20),fontsize=14,colors='k')
    if hline:
        plt.hlines(hliney,hlinemin,hlinemax,colors='orange',linewidths=2.0)
    if vline:
        plt.vlines(vlinex,vlinemin,vlinemax,colors='orange',linewidths=2.0)        
    yax1 = np.arange(y1*dx,y2*dx,1000).astype(int)
    yax2 = (yax1 - y1*dx)/float(dx)
    yax2 = yax2.astype(int)
    xax1 = np.arange(x1*dx,x2*dx,1000).astype(int)
    xax2 = (xax1 - x1*dx)/float(dx)
    xax2 = xax2.astype(int)
    plt.tick_params(axis='both',labelsize=14)
    plt.xlabel('West-East (km)',fontsize=18,labelpad=6)
    plt.xticks(xax2,xax1)
    plt.ylabel('South-North (km)',fontsize=18,labelpad=6)
    plt.yticks(yax2,yax1)
    if save_file is not None:
        plt.savefig(save_file,dpi=300,bbox_inches='tight')

    return fig

### Precipitation rate
def precip(ncfile, x1, x2, y1, y2, dx, time=0, timeres=6, roll=0, save_file=None, width=9., height=6., filetype='PDF', \
           coarse=True, coarse_val=4, rmin=1, rmax=16,title=None):

    matplotlib.use(filetype)
    rain1 = getvar(ncfile,'RAINC',timeidx=time)
    rainn1 = getvar(ncfile,'RAINNC',timeidx=time)
    precip1 = np.asarray(rain1) + np.asarray(rainn1)
    rain2 = getvar(ncfile,'RAINC',timeidx=time-1)
    rainn2 = getvar(ncfile,'RAINNC',timeidx=time-1)
    precip2 = np.asarray(rain2) + np.asarray(rainn2)
    precip_rate = precip1 #(precip1 - precip2)/timeres
    precip_rate = np.roll(precip_rate,roll,axis=1)
    precip_rate = precip_rate[y1:y2,x1:x2]
    if coarse:
        temp = precip_rate.reshape((precip_rate.shape[0] // coarse_val, coarse_val,\
                                    precip_rate.shape[1] // coarse_val, coarse_val))
        precip_rate = np.mean(temp, axis=(1,3))
        yax1 = np.arange(y1*dx,y2*dx,400).astype(int)
        yax2 = (yax1 - y1*dx)/(float(dx)*coarse_val)
        yax2 = yax2.astype(int)
        xax1 = np.arange(x1*dx,x2*dx,400).astype(int)
        xax2 = (xax1 - x1*dx)/(float(dx)*coarse_val)
        xax2 = xax2.astype(int)
    else:
        yax1 = np.arange(y1*dx,y2*dx,1000).astype(int)
        yax2 = (yax1 - y1*dx)/float(dx)
        yax2 = yax2.astype(int)
        xax1 = np.arange(x1*dx,x2*dx,1000).astype(int)
        xax2 = (xax1 - x1*dx)/float(dx)
        xax2 = xax2.astype(int)
    fig = plt.figure(figsize=(width,height))
    ax = fig.add_subplot(1,1,1)
    ax.set_aspect('equal')
    plt.contourf(precip_rate, levels=np.arange(rmin,rmax,1),cmap=get_cmap('jet'),extend='max')
    cbar = plt.colorbar(fraction=0.024, pad=0.04)
    cbar.ax.tick_params(labelsize=14)
    plt.tick_params(axis='both',labelsize=14)
    plt.xlabel('West-East (km)',fontsize=18,labelpad=6)
    plt.xticks(xax2,xax1)
    plt.ylabel('South-North (km)',fontsize=18,labelpad=6)
    plt.yticks(yax2,yax1)
    if save_file is not None:
        plt.savefig(save_file,dpi=300,bbox_inches='tight')
        
    plt.title(title,fontsize=18)    

    return fig

### Synthetic radar
def radar(ncfile, x1, x2, y1, y2, dx, time=0, roll=0, save_file=None, width=9., height=6., filetype='PDF',title=None):

    matplotlib.use(filetype)
    p = np.asarray(getvar(ncfile,'p',timeidx=time))
    tk = np.asarray(getvar(ncfile,'tk',timeidx=time))
    qv = np.asarray(getvar(ncfile,'QVAPOR',timeidx=time))
    qr = np.asarray(getvar(ncfile,'QRAIN',timeidx=time))
    qs = np.asarray(getvar(ncfile,'QSNOW',timeidx=time))
    qg = np.asarray(getvar(ncfile,'QGRAUP',timeidx=time))
    ps = getvar(ncfile,'slp',timeidx=time)
    p_surf = np.asarray(ps)
    p_surf = np.roll(p_surf,roll,axis=1)
    #p_surf = smooth2d(p_surf,150)
    p_surf = p_surf[y1:y2,x1:x2]
    ref = np.max(dbz(p,tk,qv,qr,qs),0)
    ref = np.roll(ref,roll,axis=1)
    ref = ref[y1:y2,x1:x2]    
    yax1 = np.arange(y1*dx,y2*dx,400).astype(int)
    yax2 = (yax1 - y1*dx)/float(dx)
    yax2 = yax2.astype(int)
    xax1 = np.arange(x1*dx,x2*dx,400).astype(int)
    xax2 = (xax1 - x1*dx)/float(dx)
    xax2 = xax2.astype(int)
    fig = plt.figure(figsize=(width,height))
    ax = fig.add_subplot(1,1,1)
    ax.set_aspect('equal')
    plt.contourf(ref,levels=[2*n for n in range(21)],cmap=get_cmap("jet"),extend='max')
    cbar = plt.colorbar(fraction=0.024, pad=0.04)
    cbar.ax.tick_params(labelsize=14)
    cs = ax.contour(p_surf,levels=np.arange(800,1200,4),colors='k',linewidths=2.0)
    if np.amax(p_surf) - np.amin(p_surf) > 8:
        ax.clabel(cs,fmt='%1.0f',inline=1,levels=np.arange(800,1200,8),fontsize=14,colors='k')
    plt.tick_params(axis='both',labelsize=14)
    plt.xlabel('West-East (km)',fontsize=20,labelpad=6)
    plt.xticks(xax2,xax1,fontsize=20)
    plt.ylabel('South-North (km)',fontsize=20,labelpad=6)
    plt.yticks(yax2,yax1,fontsize=20)
    plt.grid()
    if save_file is not None:
        plt.savefig(save_file,dpi=300,bbox_inches='tight')

    plt.title(title,fontsize=18)     
        
    return fig

### 500 hPa height
def hgt500(ncfile, x1, x2, y1, y2, dx, time=0, roll=0, save_file=None, width=10., height=10., filetype='PDF',title='500 hPa'):
    
    matplotlib.use(filetype)
    p = np.asarray(getvar(ncfile,'pressure',timeidx=time))
    ps = getvar(ncfile,'slp',timeidx=time)
    p_surf = np.asarray(ps)
    p_surf = np.roll(p_surf,roll,axis=1)
    #p_surf = smooth2d(p_surf,150)
    p_surf = p_surf[y1:y2,x1:x2]
    ht = getvar(ncfile,'z',units='dm',timeidx=time)
    ht_500 = interplevel(ht, p, 500.0)
    ht_500 = np.roll(ht_500,roll,axis=1)
    ht_500 = ht_500[y1:y2,x1:x2]
    fig = plt.figure(figsize=(width,height))
    ax = fig.add_subplot(1,1,1)
    ax.set_aspect('equal')
    plt.contourf(ht_500,levels=np.arange(480,590,5),cmap=get_cmap('plasma'),extend='both')
    cbar = plt.colorbar(fraction=0.024, pad=0.04)
    cbar.ax.tick_params(labelsize=14) 
    cs = ax.contour(p_surf,levels=np.arange(800,1200,10),colors='k',linewidths=2.0)
    if np.amax(p_surf) - np.amin(p_surf) > 20:
        ax.clabel(cs,fmt='%1.0f',inline=1,levels=np.arange(800,1200,20),fontsize=14,colors='k')
    yax1 = np.arange(y1*dx,y2*dx,1000).astype(int)
    yax2 = (yax1 - y1*dx)/float(dx)
    yax2 = yax2.astype(int)
    xax1 = np.arange(x1*dx,x2*dx,1000).astype(int)
    xax2 = (xax1 - x1*dx)/float(dx)
    xax2 = xax2.astype(int)
    plt.tick_params(axis='both',labelsize=14)
    plt.xlabel('West-East (km)',fontsize=18,labelpad=6)
    plt.xticks(xax2,xax1)
    plt.ylabel('South-North (km)',fontsize=18,labelpad=6)
    plt.yticks(yax2,yax1)
    plt.title(title,fontsize=18)
    if save_file is not None:
        plt.savefig(save_file,dpi=300,bbox_inches='tight')
        
    return fig

### Integrated vapor transport
def ivt(ncfile, x1, x2, y1, y2, dx, time=0, roll=0, save_file=None, width=9., height=6., filetype='PDF'):
    
    matplotlib.use(filetype)
    p = np.asarray(getvar(ncfile,'p',timeidx=time))
    ps = getvar(ncfile,'slp',timeidx=time)
    p_surf = np.asarray(ps)
    u = getvar(ncfile,'ua',timeidx=time)
    v = getvar(ncfile,'va',timeidx=time)
    qv = np.asarray(getvar(ncfile,'QVAPOR',timeidx=time))
    ht = getvar(ncfile,'z',units='dm',timeidx=time)
    ht_700 = interplevel(ht,p/100.,700.0)
    ht_700 = np.roll(ht_700,roll,axis=1)
    ht_700 = ht_700[y1:y2,x1:x2]
    ht_700 = smooth2d(ht_700,50)
    u_ivt = -simps(u*qv,p,axis=0)/9.81
    v_ivt = -simps(v*qv,p,axis=0)/9.81 
    ivt = np.sqrt(u_ivt**2 + v_ivt**2)
    u_ivt = np.roll(u_ivt,roll,axis=1)
    u_ivt = u_ivt[y1:y2,x1:x2]
    v_ivt = np.roll(v_ivt,roll,axis=1)
    v_ivt = v_ivt[y1:y2,x1:x2]
    ivt = np.roll(ivt,roll,axis=1)
    ivt = ivt[y1:y2,x1:x2]
    fig = plt.figure(figsize=(width,height))
    ax = fig.add_subplot(1,1,1)
    ax.set_aspect('equal')
    [X,Y] = np.meshgrid(np.arange(0,int(x2-x1)),np.arange(0,int(y2-y1)))
    plt.contourf(X,Y,ivt,levels=np.arange(100,700,100),cmap=get_cmap('GnBu'),extend='max')
    cbar = plt.colorbar(fraction=0.024, pad=0.04)
    cbar.ax.tick_params(labelsize=14)
    cs = ax.contour(X,Y,ht_700,levels=np.arange(100,500,10),colors='k',linewidths=1.0)
    ax.clabel(cs,fmt='%1.0f',inline=1,levels=np.arange(100,500,20),fontsize=14,colors='k')
    ax.quiver(X[100::100,100::100],Y[100::100,100::100],u_ivt[100::100,100::100],v_ivt[100::100,100::100],\
              units='xy',scale=3.,pivot='tip',headwidth=3.,headlength=4.,headaxislength=3.5)
    yax1 = np.arange(y1*dx,y2*dx,1000).astype(int)
    yax2 = (yax1 - y1*dx)/float(dx)
    yax2 = yax2.astype(int)
    xax1 = np.arange(x1*dx,x2*dx,1000).astype(int)
    xax2 = (xax1 - x1*dx)/float(dx)
    xax2 = xax2.astype(int)
    plt.tick_params(axis='both',labelsize=14)
    plt.xlabel('West-East (km)',fontsize=18,labelpad=6)
    plt.xticks(xax2,xax1)
    plt.ylabel('South-North (km)',fontsize=18,labelpad=6)
    plt.yticks(yax2,yax1)
    if save_file is not None:
        plt.savefig(save_file,dpi=300,bbox_inches='tight')
        
    return fig

### PV on theta surface
def pv_th(ncfile, x1, x2, y1, y2, dx, th_lvl=310, time=0, roll=0, save_file=None, width=9., height=6., filetype='PDF',title='PV'):
    
    matplotlib.use(filetype)
    p = np.asarray(getvar(ncfile,'pressure',timeidx=time))
    ps = getvar(ncfile,'slp',timeidx=time)
    p_surf = np.asarray(ps)
    p_surf = np.roll(p_surf,roll,axis=1)
    p_surf = smooth2d(p_surf,150)
    p_surf = p_surf[y1:y2,x1:x2]
    th = getvar(ncfile,'th',timeidx=time)
    pv = getvar(ncfile,'pvo',timeidx=time)
    pv_th = interplevel(pv, th, th_lvl)
    pv_th = np.roll(pv_th,roll,axis=1)
    pv_th = pv_th[y1:y2,x1:x2]
    fig = plt.figure(figsize=(width,height))
    ax = fig.add_subplot(1,1,1)
    ax.set_aspect('equal')
    plt.contourf(pv_th,levels=np.arange(0,5.5,0.5),cmap=get_cmap('coolwarm_r'),extend='both')
    cbar = plt.colorbar(fraction=0.024, pad=0.04)
    cbar.ax.tick_params(labelsize=14) 
    cs = ax.contour(p_surf,levels=np.arange(952,1060,8),colors='k',linewidths=2.0)
    if np.amax(p_surf) - np.amin(p_surf) > 20:
        ax.clabel(cs,fmt='%1.0f',inline=1,levels=np.arange(960,1060,16),fontsize=14,colors='k')
    yax1 = np.arange(y1*dx,y2*dx,1000).astype(int)
    yax2 = (yax1 - y1*dx)/float(dx)
    yax2 = yax2.astype(int)
    xax1 = np.arange(x1*dx,x2*dx,2000).astype(int)
    xax2 = (xax1 - x1*dx)/float(dx)
    xax2 = xax2.astype(int)
    plt.tick_params(axis='both',labelsize=14)
    plt.xlabel('West-East (km)',fontsize=18,labelpad=6)
    plt.xticks(xax2,xax1)
    plt.ylabel('South-North (km)',fontsize=18,labelpad=6)
    plt.yticks(yax2,yax1)
    plt.title(title,fontsize=18)
    if save_file is not None:
        plt.savefig(save_file,dpi=300,bbox_inches='tight')
        
    return fig

### Pressure on theta surface
def pres_th(ncfile, x1, x2, y1, y2, dx, th_lvl=300, time=0, roll=0, save_file=None, width=9., height=6., filetype='PDF'):
    
    matplotlib.use(filetype)
    p = np.asarray(getvar(ncfile,'pressure',timeidx=time))
    ps = getvar(ncfile,'slp',timeidx=time)
    p_surf = np.asarray(ps)
    p_surf = np.roll(p_surf,roll,axis=1)
    p_surf = smooth2d(p_surf,150)
    p_surf = p_surf[y1:y2,x1:x2]
    th = getvar(ncfile,'th',timeidx=time)
    p_th = interplevel(p, th, th_lvl)
    p_th = np.roll(p_th,roll,axis=1)
    p_th = p_th[y1:y2,x1:x2]
    fig = plt.figure(figsize=(width,height))
    ax = fig.add_subplot(1,1,1)
    ax.set_aspect('equal')
    plt.contourf(p_th,levels=np.arange(250,700,50),cmap=get_cmap('coolwarm'),extend='both')
    cbar = plt.colorbar(fraction=0.024, pad=0.04)
    cbar.ax.tick_params(labelsize=14) 
    cs = ax.contour(p_surf,levels=np.arange(800,1200,10),colors='k',linewidths=2.0)
    if np.amax(p_surf) - np.amin(p_surf) > 20:
        ax.clabel(cs,fmt='%1.0f',inline=1,levels=np.arange(800,1200,20),fontsize=14,colors='k')
    yax1 = np.arange(y1*dx,y2*dx,1000).astype(int)
    yax2 = (yax1 - y1*dx)/float(dx)
    yax2 = yax2.astype(int)
    xax1 = np.arange(x1*dx,x2*dx,1000).astype(int)
    xax2 = (xax1 - x1*dx)/float(dx)
    xax2 = xax2.astype(int)
    plt.tick_params(axis='both',labelsize=14)
    plt.xlabel('West-East (km)',fontsize=18,labelpad=6)
    plt.xticks(xax2,xax1)
    plt.ylabel('South-North (km)',fontsize=18,labelpad=6)
    plt.yticks(yax2,yax1)
    if save_file is not None:
        plt.savefig(save_file,dpi=300,bbox_inches='tight')
        
    return fig

### Pressure on PV surface
def pres_pv(ncfile, x1, x2, y1, y2, dx, pv_lvl=1.5, time=0, roll=0, save_file=None, width=9., height=6., filetype='PDF'):
    
    matplotlib.use(filetype)
    p = np.asarray(getvar(ncfile,'pressure',timeidx=time))
    ps = getvar(ncfile,'slp',timeidx=time)
    p_surf = np.asarray(ps)
    p_surf = np.roll(p_surf,roll,axis=1)
    p_surf = smooth2d(p_surf,150)
    p_surf = p_surf[y1:y2,x1:x2]
    pv = getvar(ncfile,'pvo',timeidx=time)
    p_pv = interplevel(p, pv, pv_lvl)
    p_pv = np.roll(p_pv,roll,axis=1)
    p_pv = p_pv[y1:y2,x1:x2]
    fig = plt.figure(figsize=(width,height))
    ax = fig.add_subplot(1,1,1)
    ax.set_aspect('equal')
    plt.contourf(p_pv,levels=np.arange(180,400,20),cmap=get_cmap('coolwarm_r'),extend='both')
    cbar = plt.colorbar(fraction=0.024, pad=0.04)
    cbar.ax.tick_params(labelsize=14) 
    cs = ax.contour(p_surf,levels=np.arange(800,1200,10),colors='k',linewidths=2.0)
    if np.amax(p_surf) - np.amin(p_surf) > 20:
        ax.clabel(cs,fmt='%1.0f',inline=1,levels=np.arange(800,1200,20),fontsize=14,colors='k')
    yax1 = np.arange(y1*dx,y2*dx,1000).astype(int)
    yax2 = (yax1 - y1*dx)/float(dx)
    yax2 = yax2.astype(int)
    xax1 = np.arange(x1*dx,x2*dx,1000).astype(int)
    xax2 = (xax1 - x1*dx)/float(dx)
    xax2 = xax2.astype(int)
    plt.tick_params(axis='both',labelsize=14)
    plt.xlabel('West-East (km)',fontsize=18,labelpad=6)
    plt.xticks(xax2,xax1)
    plt.ylabel('South-North (km)',fontsize=18,labelpad=6)
    plt.yticks(yax2,yax1)
    if save_file is not None:
        plt.savefig(save_file,dpi=300,bbox_inches='tight')
        
    return fig

### Theta on PV surface
def th_pv(ncfile, x1, x2, y1, y2, dx, pv_lvl=1.5, time=0, roll=0, save_file=None, width=9., height=6., filetype='PDF'):
    
    matplotlib.use(filetype)
    p = np.asarray(getvar(ncfile,'pressure',timeidx=time))
    ps = getvar(ncfile,'slp',timeidx=time)
    p_surf = np.asarray(ps)
    p_surf = np.roll(p_surf,roll,axis=1)
    p_surf = smooth2d(p_surf,150)
    p_surf = p_surf[y1:y2,x1:x2]
    pv = getvar(ncfile,'pvo',timeidx=time)
    th = getvar(ncfile,'th',timeidx=time)
    th_pv = interplevel(th, pv, pv_lvl)
    th_pv = np.roll(th_pv,roll,axis=1)
    th_pv = th_pv[y1:y2,x1:x2]
    fig = plt.figure(figsize=(width,height))
    ax = fig.add_subplot(1,1,1)
    ax.set_aspect('equal')
    plt.contourf(th_pv,levels=np.arange(270,325,5),cmap=get_cmap('coolwarm'),extend='both')
    cbar = plt.colorbar(fraction=0.024, pad=0.04)
    cbar.ax.tick_params(labelsize=14) 
    cs = ax.contour(p_surf,levels=np.arange(800,1200,10),colors='k',linewidths=2.0)
    if np.amax(p_surf) - np.amin(p_surf) > 20:
        ax.clabel(cs,fmt='%1.0f',inline=1,levels=np.arange(800,1200,20),fontsize=14,colors='k')
    yax1 = np.arange(y1*dx,y2*dx,1000).astype(int)
    yax2 = (yax1 - y1*dx)/float(dx)
    yax2 = yax2.astype(int)
    xax1 = np.arange(x1*dx,x2*dx,1000).astype(int)
    xax2 = (xax1 - x1*dx)/float(dx)
    xax2 = xax2.astype(int)
    plt.tick_params(axis='both',labelsize=14)
    plt.xlabel('West-East (km)',fontsize=18,labelpad=6)
    plt.xticks(xax2,xax1)
    plt.ylabel('South-North (km)',fontsize=18,labelpad=6)
    plt.yticks(yax2,yax1)
    if save_file is not None:
        plt.savefig(save_file,dpi=300,bbox_inches='tight')
        
    return fig

### Vertical velocity
def w_plan(ncfile, x1, x2, y1, y2, dx, z_lvl=4, time=0, roll=0, save_file=None, width=9., height=6., filetype='PDF', \
           levels=None):

    matplotlib.use(filetype)
    p = np.asarray(getvar(ncfile,'pressure',timeidx=time))
    ps = getvar(ncfile,'slp',timeidx=time)
    p_surf = np.asarray(ps)
    p_surf = np.roll(p_surf,roll,axis=1)
    p_surf = smooth2d(p_surf,150)
    p_surf = p_surf[y1:y2,x1:x2]
    ht = getvar(ncfile,'z',units='km',timeidx=time)
    w = getvar(ncfile,'wa',timeidx=time)
    w_z = interplevel(w, ht, z_lvl)
    w_z = np.roll(w_z,roll,axis=1)
    w_z = w_z[y1:y2,x1:x2]
    yax1 = np.arange(y1*dx,y2*dx,1000).astype(int)
    yax2 = (yax1 - y1*dx)/float(dx)
    yax2 = yax2.astype(int)
    xax1 = np.arange(x1*dx,x2*dx,1000).astype(int)
    xax2 = (xax1 - x1*dx)/float(dx)
    xax2 = xax2.astype(int)
    fig = plt.figure(figsize=(width,height))
    ax = fig.add_subplot(1,1,1)
    ax.set_aspect('equal')
    p = matplotlib.colors.ListedColormap(['blue','dodgerblue','deepskyblue','lightblue','white','navajowhite',\
                                          'orange','red','firebrick'])
    p.set_over('maroon')
    p.set_under('navy')
    if levels is not None:
        lvls = levels
    else:
        lvls = [-1.0,-0.5,-0.1,-0.05,-0.01,0.01,0.05,0.1,0.5,1.0]
    nrm = matplotlib.colors.BoundaryNorm(lvls,p.N)
    plt.contourf(w_z,levels=lvls,cmap=p,norm=nrm,extend='both')
    cbar = plt.colorbar(fraction=0.024, pad=0.04)
    cbar.ax.tick_params(labelsize=14)
    cs = ax.contour(p_surf,levels=np.arange(800,1200,10),colors='k',linewidths=2.0)
    if np.amax(p_surf) - np.amin(p_surf) > 20:
        ax.clabel(cs,fmt='%1.0f',inline=1,levels=np.arange(800,1200,20),fontsize=14,colors='k')
    plt.tick_params(axis='both',labelsize=14)
    plt.xlabel('West-East (km)',fontsize=18,labelpad=6)
    plt.xticks(xax2,xax1)
    plt.ylabel('South-North (km)',fontsize=18,labelpad=6)
    plt.yticks(yax2,yax1)
    if save_file is not None:
        plt.savefig(save_file,dpi=300,bbox_inches='tight')

    return fig

### Column-Maximum DTE
def dte_plan(base_ncfile, pert_ncfile, x1, x2, y1, y2, dx, time=0, roll=0, save_file=None, width=9., height=6.,\
             filetype='PDF',title='DTE'):
    
    matplotlib.use(filetype)
    u_base = np.asarray(getvar(base_ncfile,'ua',timeidx=time))
    v_base = np.asarray(getvar(base_ncfile,'va',timeidx=time))
    tk_base = np.asarray(getvar(base_ncfile,'tk',timeidx=time))
    u_pert = np.asarray(getvar(pert_ncfile,'ua',timeidx=time))
    v_pert = np.asarray(getvar(pert_ncfile,'va',timeidx=time))
    tk_pert = np.asarray(getvar(pert_ncfile,'tk',timeidx=time))
    du = u_base - u_pert
    dv = v_base - v_pert
    dt = tk_base - tk_pert
    dte = 0.5*(du**2 + dv**2 + CP/TR*dt**2)
    dte_max = np.amax(dte,axis=0)
    dte_max = np.roll(dte_max,roll,axis=1)[y1:y2,x1:x2]
    ps = np.asarray(getvar(base_ncfile,'slp',timeidx=time))
    #ps = getvar(pert_ncfile,'slp',timeidx=int(time-8))
    p_surf = np.asarray(ps)
    p_surf = np.roll(p_surf,roll,axis=1)
    #p_surf = smooth2d(p_surf,150)
    p_surf = p_surf[y1:y2,x1:x2]
    fig = plt.figure(figsize=(width,height),dpi=300)
    ax = fig.add_subplot(1,1,1)
    ax.set_aspect('equal')
    p = matplotlib.colors.LinearSegmentedColormap.from_list("",["lightblue","cyan","limegreen","gold","magenta"])
    p.set_over("purple")
    lvls = [0.5,2,8,32,128,512]
    #lvls = [0.001,0.01,0.1,1,10,100]
    nm = matplotlib.colors.BoundaryNorm(lvls,p.N)
    plt.contourf(dte_max,levels=lvls,cmap=p,norm=nm,extend='max')
    cbar = plt.colorbar(fraction=0.024, pad=0.04)
    cbar.ax.tick_params(labelsize=14) 
    cs = ax.contour(p_surf,levels=np.arange(800,1200,10),colors='k',linewidths=2.0)
    if np.amax(p_surf) - np.amin(p_surf) > 20:
        ax.clabel(cs,fmt='%1.0f',inline=1,levels=np.arange(800,1200,20),fontsize=14,colors='k')
    yax1 = np.arange(y1*dx,y2*dx,1000).astype(int)
    yax2 = (yax1 - y1*dx)/float(dx)
    yax2 = yax2.astype(int)
    xax1 = np.arange(x1*dx,x2*dx,1000).astype(int)
    xax2 = (xax1 - x1*dx)/float(dx)
    xax2 = xax2.astype(int)
    plt.tick_params(axis='both',labelsize=14)
    plt.xlabel('West-East (km)',fontsize=18,labelpad=6)
    plt.xticks(xax2,xax1)
    plt.ylabel('South-North (km)',fontsize=18,labelpad=6)
    plt.yticks(yax2,yax1)
    plt.title(title,fontsize=18)
    if save_file is not None:
        plt.savefig(save_file,dpi=300,bbox_inches='tight')
        
    return fig

#-------------------- CROSS SECTIONS --------------------#

### Jet
def jet(ncfile, y1, y2, z1, z2, x, dx, time=0, save_file=None, width=9., height=6., filetype='PDF'):

    matplotlib.use(filetype)
    u = np.asarray(getvar(ncfile,'ua',timeidx=time))[z1:z2,y1:y2,x]
    th = np.asarray(getvar(ncfile,'th',timeidx=time))[z1:z2,y1:y2,x]
    qv = np.asarray(getvar(ncfile,'QVAPOR',timeidx=time))[z1:z2,y1:y2,x]
    pv = np.asarray(getvar(ncfile,'pvo',timeidx=time))[z1:z2,y1:y2,x]
    z = np.asarray(getvar(ncfile,'z',timeidx=time))[z1:z2,y1:y2,x]/1000.
    zgrid = np.arange(z1,z2)
    ygrid = np.arange(y1*dx,y2*dx,dx)
    y,zmesh = np.meshgrid(ygrid,zgrid)
    fig = plt.figure(figsize=(width,height))
    ax = fig.add_subplot(1,1,1)
    plt.contourf(y,z,qv*1000,levels=np.arange(2,16,2),extend='max')
    plt.set_cmap('GnBu')
    plt.colorbar(fraction=0.024, pad=0.04)
    cs = ax.contour(y,z,th,levels=np.arange(200,560,5),colors='lightgrey')
    ax.clabel(cs,fmt='%1.0f',inline=1,levels=[280,320,360],fontsize=8,colors='black')
    cs2 = ax.contour(y,z,u,levels=np.arange(-100,105,5),colors='red')
    ax.clabel(cs2,fmt='%1.0f',inline=1,fontsize=8)
    cs3 = ax.contour(y,z,pv,levels=[1.5],colors='k',linewidths=2.5)
    plt.tick_params(axis='both',labelsize=14)
    plt.xlabel('South-North (km)',fontsize=12)
    plt.ylabel('Height (km)',fontsize=12)
    #plt.grid()
    if save_file is not None:
        plt.savefig(save_file,dpi=300,bbox_inches='tight')

    return fig
    

### Vertical velocity
def w_cross_x(ncfile, x1, x2, z1, z2, y, dx, time=0, roll=0, save_file=None, levels=None, \
              width=9., height=6., filetype='PDF'):

    matplotlib.use(filetype)
    w = np.roll(np.asarray(getvar(ncfile,'wa',timeidx=time)),roll,axis=2)[z1:z2,y,x1:x2]
    th = np.roll(np.asarray(getvar(ncfile,'th',timeidx=time)),roll,axis=2)[z1:z2,y,x1:x2]
    z = np.roll(np.asarray(getvar(ncfile,'z',timeidx=time)),roll,axis=2)[z1:z2,y,x1:x2]/1000.
    zgrid = np.arange(z1,z2)
    xgrid = np.arange(x1*dx,x2*dx,dx)
    x,zmesh = np.meshgrid(xgrid,zgrid)    
    fig = plt.figure(figsize=(width,height))
    ax = fig.add_subplot(1,1,1)
    p = matplotlib.colors.ListedColormap(['blue','dodgerblue','deepskyblue','lightblue','white','navajowhite',\
                                          'orange','red','firebrick'])
    p.set_over('maroon')
    p.set_under('navy')
    if levels is not None:
        lvls = levels
    else:
        lvls = [-1.0,-0.5,-0.1,-0.05,-0.01,0.01,0.05,0.1,0.5,1.0]
    nrm = matplotlib.colors.BoundaryNorm(lvls,p.N)
    plt.contourf(x,z,w,levels=lvls,cmap=p,norm=nrm,extend='both')
    cbar = plt.colorbar(fraction=0.024, pad=0.04)
    cbar.ax.tick_params(labelsize=14)
    cs = ax.contour(x,z,th,levels=np.arange(100,500,4),colors='k',linewidths=2.0)
    ax.clabel(cs,fmt='%1.0f',inline=1,levels=np.arange(100,500,20),fontsize=14,colors='k')
    plt.tick_params(axis='both',labelsize=14) 
    plt.xlabel('West-East (km)',fontsize=18,labelpad=6)
    plt.ylabel('Height (km)',fontsize=18,labelpad=6)
    if save_file is not None:
        plt.savefig(save_file,dpi=300,bbox_inches='tight')

    return fig

def w_cross_y(ncfile, y1, y2, z1, z2, x, dx, time=0, roll=0, save_file=None, levels=None, \
              width=9., height=6., filetype='PDF'):

    matplotlib.use(filetype)
    w = np.roll(np.asarray(getvar(ncfile,'wa',timeidx=time)),roll,axis=2)[z1:z2,y1:y2,x]
    th = np.roll(np.asarray(getvar(ncfile,'th',timeidx=time)),roll,axis=2)[z1:z2,y1:y2,x]
    z = np.roll(np.asarray(getvar(ncfile,'z',timeidx=time)),roll,axis=2)[z1:z2,y1:y2,x]/1000.
    zgrid = np.arange(z1,z2)
    ygrid = np.arange(y1*dx,y2*dx,dx)
    y,zmesh = np.meshgrid(ygrid,zgrid)    
    fig = plt.figure(figsize=(width,height))
    ax = fig.add_subplot(1,1,1)
    p = matplotlib.colors.ListedColormap(['blue','dodgerblue','deepskyblue','lightblue','white','navajowhite',\
                                          'orange','red','firebrick'])
    p.set_over('maroon')
    p.set_under('navy')
    if levels is not None:
        lvls = levels
    else:
        lvls = [-1.0,-0.5,-0.1,-0.05,-0.01,0.01,0.05,0.1,0.5,1.0]
    nrm = matplotlib.colors.BoundaryNorm(lvls,p.N)
    plt.contourf(y,z,w,levels=lvls,cmap=p,norm=nrm,extend='both')
    cbar = plt.colorbar(fraction=0.024, pad=0.04)
    cbar.ax.tick_params(labelsize=14)
    cs = ax.contour(y,z,th,levels=np.arange(100,500,4),colors='k',linewidths=2.0)
    ax.clabel(cs,fmt='%1.0f',inline=1,levels=np.arange(100,500,20),fontsize=14,colors='k')
    plt.tick_params(axis='both',labelsize=14) 
    plt.xlabel('South-North (km)',fontsize=18,labelpad=6)
    plt.ylabel('Height (km)',fontsize=18,labelpad=6)
    if save_file is not None:
        plt.savefig(save_file,dpi=300,bbox_inches='tight')

    return fig
                    
### Moisture
def moist_cross_x(ncfile, x1, x2, z1, z2, y, dx, time=0, roll=0, save_file=None, \
              width=9., height=6., filetype='PDF'):

    matplotlib.use(filetype)
    ice = np.roll(np.asarray(getvar(ncfile,'QICE',timeidx=time)),roll,axis=2)[z1:z2,y,x1:x2]*1000.
    cloud = np.roll(np.asarray(getvar(ncfile,'QCLOUD',timeidx=time)),roll,axis=2)[z1:z2,y,x1:x2]*1000.
    rain = np.roll(np.asarray(getvar(ncfile,'QRAIN',timeidx=time)),roll,axis=2)[z1:z2,y,x1:x2]*1000.
    snow = np.roll(np.asarray(getvar(ncfile,'QSNOW',timeidx=time)),roll,axis=2)[z1:z2,y,x1:x2]*1000.                 
    z = np.roll(np.asarray(getvar(ncfile,'z',timeidx=time)),roll,axis=2)[z1:z2,y,x1:x2]/1000.
    zgrid = np.arange(z1,z2)
    xgrid = np.arange(x1*dx,x2*dx,dx)
    x,zmesh = np.meshgrid(xgrid,zgrid)    
    fig = plt.figure(figsize=(width,height))
    ax = fig.add_subplot(1,1,1)
    plt.tick_params(axis='both',labelsize=14)
    p1 = matplotlib.colors.ListedColormap(['lightblue','deepskyblue','dodgerblue','blue','darkblue'])
    p1.set_over('midnightblue')
    p2 = matplotlib.colors.ListedColormap(['palegreen','mediumspringgreen','mediumseagreen','forestgreen','darkgreen'])
    p2.set_over('darkslategray')
    lvls1=[0.001,0.004,0.016,0.064,0.256,1.024]
    lvls2=[0.001,0.005,0.025,0.125,0.625,3.125]
    nrm1 = matplotlib.colors.BoundaryNorm(lvls1,p1.N)
    nrm2 = matplotlib.colors.BoundaryNorm(lvls2,p2.N)
    plot1 = plt.contourf(x,z,ice,levels=lvls1,cmap=p1,norm=nrm1,extend='max') 
    plot2 = plt.contourf(x,z,cloud,levels=lvls2,cmap=p2,norm=nrm2,extend='max')
    cs1 = ax.contour(x,z,rain,levels=[0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0],colors='orange',linewidths=2.0)
    cs2 = ax.contour(x,z,snow,levels=[0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0],colors='pink',linewidths=2.0)
    plt.xlabel('West-East (km)',fontsize=18,labelpad=6)
    plt.ylabel('Height (km)',fontsize=18,labelpad=6)
    plt.subplots_adjust(bottom=0.1,right=0.8,top=0.9,left=0.2)
    cax1 = plt.axes([0.82,0.18,0.02,0.6])
    cax2 = plt.axes([0.10,0.18,0.02,0.6])
    cb1 = plt.colorbar(plot1,cax=cax1)
    cb2 = plt.colorbar(plot2,cax=cax2)
    cax2.yaxis.set_ticks_position('left')
    if save_file is not None:
        plt.savefig(save_file,dpi=300,bbox_inches='tight')

    return fig                    

def moist_cross_y(ncfile, y1, y2, z1, z2, x, dx, time=0, roll=0, save_file=None, width=9., height=6., filetype='PDF'):

    matplotlib.use(filetype)
    ice = np.roll(np.asarray(getvar(ncfile,'QICE',timeidx=time)),roll,axis=2)[z1:z2,y1:y2,x]*1000.
    cloud = np.roll(np.asarray(getvar(ncfile,'QCLOUD',timeidx=time)),roll,axis=2)[z1:z2,y1:y2,x]*1000.
    rain = np.roll(np.asarray(getvar(ncfile,'QRAIN',timeidx=time)),roll,axis=2)[z1:z2,y1:y2,x]*1000.
    snow = np.roll(np.asarray(getvar(ncfile,'QSNOW',timeidx=time)),roll,axis=2)[z1:z2,y1:y2,x]*1000.                 
    z = np.roll(np.asarray(getvar(ncfile,'z',timeidx=time)),roll,axis=2)[z1:z2,y,x1:x2]/1000.
    zgrid = np.arange(z1,z2)
    ygrid = np.arange(y1*dx,y2*dx,dx)
    y,zmesh = np.meshgrid(ygrid,zgrid)    
    fig = plt.figure(figsize=(width,height))
    ax = fig.add_subplot(1,1,1)
    plt.tick_params(axis='both',labelsize=14)
    p1 = matplotlib.colors.ListedColormap(['lightblue','deepskyblue','dodgerblue','blue','darkblue'])
    p1.set_over('midnightblue')
    p2 = matplotlib.colors.ListedColormap(['palegreen','mediumspringgreen','mediumseagreen','forestgreen','darkgreen'])
    p2.set_over('darkslategray')
    lvls1=[0.001,0.004,0.016,0.064,0.256,1.024]
    lvls2=[0.001,0.005,0.025,0.125,0.625,3.125]
    nrm1 = matplotlib.colors.BoundaryNorm(lvls1,p1.N)
    nrm2 = matplotlib.colors.BoundaryNorm(lvls2,p2.N)
    plot1 = plt.contourf(y,z,ice,levels=lvls1,cmap=p1,norm=nrm1,extend='max') 
    plot2 = plt.contourf(y,z,cloud,levels=lvls2,cmap=p2,norm=nrm2,extend='max')
    cs1 = ax.contour(y,z,rain,levels=[0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0],colors='orange',linewidths=2.0)
    cs2 = ax.contour(y,z,snow,levels=[0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0],colors='pink',linewidths=2.0)
    plt.xlabel('South-North (km)',fontsize=18,labelpad=6)
    plt.ylabel('Height (km)',fontsize=18,labelpad=6)
    plt.subplots_adjust(bottom=0.1,right=0.8,top=0.9,left=0.2)
    cax1 = plt.axes([0.82,0.18,0.02,0.6])
    cax2 = plt.axes([0.10,0.18,0.02,0.6])
    cb1 = plt.colorbar(plot1,cax=cax1)
    cb2 = plt.colorbar(plot2,cax=cax2)
    cax2.yaxis.set_ticks_position('left')
    if save_file is not None:
        plt.savefig(save_file,dpi=300,bbox_inches='tight')

    return fig

### PV
def pv_cross(ncfile, y1, y2, z1, z2, x, dx, time=0, save_file=None, width=9., height=6., filetype='PDF',roll=0):

    matplotlib.use(filetype)
    pv = np.roll(np.asarray(getvar(ncfile,'pvo',timeidx=time)),roll,axis=2)[z1:z2,y1:y2,x]
    z = np.roll(np.asarray(getvar(ncfile,'z',timeidx=time)),roll,axis=2)[z1:z2,y1:y2,x]/1000.
    zgrid = np.arange(z1,z2)
    ygrid = np.arange(y1*dx,y2*dx,dx)
    y,zmesh = np.meshgrid(ygrid,zgrid)
    fig = plt.figure(figsize=(width,height))
    ax = fig.add_subplot(1,1,1)
    plt.contourf(y,z,pv,levels=np.arange(0,10.5,0.5),cmap=get_cmap('tab20'),extend='both')
    plt.colorbar(fraction=0.024, pad=0.04)
    plt.tick_params(axis='both',labelsize=14)
    plt.xlabel('South-North (km)',fontsize=12)
    plt.ylabel('Height (km)',fontsize=12)
    if save_file is not None:
        plt.savefig(save_file,dpi=300,bbox_inches='tight')

    return fig

#-------------------- LINE PLOTS --------------------#

### Static stability
def stability(ncfile, x1, x2, z1, z2, y, dx, time=0, save_file=None, width=9., height=6., filetype='PDF'):
    
    matplotlib.use(filetype)
    th = np.mean(np.asarray(getvar(ncfile,'th',timeidx=time))[z1:z2,y,x1:x2],axis=1)
    z = np.mean(np.asarray(getvar(ncfile,'z',timeidx=time))[z1:z2,y,x1:x2],axis=1)
    nz = int(z2-z1)
    n = np.zeros(nz)
    for k in range(1,nz-1):
        dtdz = (th[k+1] - th[k-1])/(z[k+1] - z[k-1])
        n[k] = np.sqrt(9.81/th[k]*(dtdz))
    dtdz_0 = (th[1] - th[0])/(z[1] - z[0])
    n[0] = np.sqrt(9.81/th[0]*(dtdz_0))
    dtdz_top = (th[nz-1] - th[nz-2])/(z[nz-1] - z[nz-2])
    n[nz-1] = np.sqrt(9.81/th[nz-1]*(dtdz_top))
    fig = plt.figure(figsize=(width,height))
    ax = fig.add_subplot(1,1,1)
    plt.plot(n,z/1000.,color='blue',linewidth=3.0)
    plt.xlim(0.005,0.03)
    plt.grid()
    plt.tick_params(axis='both',labelsize=14) 
    plt.xlabel('N (s$^{-1}$)',fontsize=18,labelpad=6)
    plt.ylabel('Height (km)',fontsize=18,labelpad=6)
    if save_file is not None:
        plt.savefig(save_file,dpi=300,bbox_inches='tight')

    return fig