#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  3 15:54:47 2020

@author: lloverasdan
"""

import numpy as np
from wrf import getvar, smooth2d, interplevel
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap

def surf_anim(ncfile, x1, x2, y1, y2, dx, time, roll=0, width=12., height=6., tmin=-30, tmax=35):
    
    matplotlib.use('Agg')
    p = getvar(ncfile,'pressure',timeidx=time)
    p_surf = np.asarray(p)[0,:,:]
    p_surf = np.roll(p_surf,roll,axis=1)
    p_surf = smooth2d(p_surf,150)
    p_surf = p_surf[y1:y2,x1:x2]
    tk = getvar(ncfile,'tk',timeidx=time)
    t_surf = np.asarray(tk)[0,:,:]-273.15
    t_surf = np.roll(t_surf,roll,axis=1)
    t_surf = smooth2d(t_surf,150)
    t_surf = t_surf[y1:y2,x1:x2]
    fig = plt.figure(figsize=(width,height),dpi=300)
    ax = fig.add_subplot(1,1,1)
    plt.contourf(t_surf,levels=np.arange(tmin,tmax,5),cmap=get_cmap('bwr'),extend='both')
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
    hour = time*6
    plt.title('Surface Pressure and Temperature at Hour %i' %hour)
    fig.canvas.draw()
    image = np.frombuffer(fig.canvas.tostring_rgb(), dtype='uint8')
    image = image.reshape(fig.canvas.get_width_height()[::-1] + (3,))
    
    return image

def cloud_anim(ncfile, x1, x2, y1, y2, dx, time, roll=0, width=12., height=6.):
               
    matplotlib.use('Agg')
    p = getvar(ncfile,'pressure',timeidx=time)
    p_surf = np.asarray(p)[0,:,:]
    p_surf = np.roll(p_surf,roll,axis=1)
    p_surf = smooth2d(p_surf,150)
    p_surf = p_surf[y1:y2,x1:x2]
    ctt = getvar(ncfile,'ctt',timeidx=time,fill_nocloud=True)
    ctt = np.asarray(ctt)
    ctt = np.roll(ctt,roll,axis=1)[y1:y2,x1:x2]               
    fig = plt.figure(figsize=(width,height),dpi=300)
    ax = fig.add_subplot(1,1,1)
    ax.set_facecolor('lightblue')
    plt.contourf(ctt,levels=np.arange(-80,10,10),cmap='Greys',extend='min')
    cbar = plt.colorbar(fraction=0.024, pad=0.04)
    cbar.ax.tick_params(labelsize=14)
    cs = ax.contour(p_surf,levels=np.arange(800,1200,10),colors='k',linewidths=2.0,alpha=0.5)
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
    hour = time*6
    plt.title('Cloud Top Temperature at Hour %i' %hour)
    fig.canvas.draw()
    image = np.frombuffer(fig.canvas.tostring_rgb(), dtype='uint8')
    image = image.reshape(fig.canvas.get_width_height()[::-1] + (3,))
    
    return image

def hgt500_anim(ncfile, x1, x2, y1, y2, dx, time, roll=0, width=12., height=6.):
    
    matplotlib.use('Agg')
    p = getvar(ncfile,'pressure',timeidx=time)
    p_surf = np.asarray(p)[0,:,:]
    p_surf = np.roll(p_surf,roll,axis=1)
    p_surf = smooth2d(p_surf,150)
    p_surf = p_surf[y1:y2,x1:x2]
    ht = getvar(ncfile,'z',units='dm',timeidx=time)
    ht_500 = interplevel(ht, p, 500.0)
    ht_500 = np.roll(ht_500,roll,axis=1)
    ht_500 = ht_500[y1:y2,x1:x2]
    fig = plt.figure(figsize=(width,height),dpi=300)
    ax = fig.add_subplot(1,1,1)
    plt.contourf(ht_500,levels=np.arange(480,590,10),cmap=get_cmap('plasma'),extend='both')
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
    hour = time*6
    plt.title('500 hPa Geopotential Height at Hour %i' %hour)
    fig.canvas.draw()
    image = np.frombuffer(fig.canvas.tostring_rgb(), dtype='uint8')
    image = image.reshape(fig.canvas.get_width_height()[::-1] + (3,))
    
    return image

def pvth_anim(ncfile, x1, x2, y1, y2, dx, time, th_lvl=300., roll=0, width=12., height=6.):
    
    matplotlib.use('Agg')
    p = getvar(ncfile,'pressure',timeidx=time)
    p_surf = np.asarray(p)[0,:,:]
    p_surf = np.roll(p_surf,roll,axis=1)
    p_surf = smooth2d(p_surf,150)
    p_surf = p_surf[y1:y2,x1:x2]
    th = getvar(ncfile,'th',timeidx=time)
    pv = getvar(ncfile,'pvo',timeidx=time)
    pv_th = interplevel(pv, th, th_lvl)
    pv_th = np.roll(pv_th,roll,axis=1)
    pv_th = pv_th[y1:y2,x1:x2]
    fig = plt.figure(figsize=(width,height),dpi=300)
    ax = fig.add_subplot(1,1,1)
    plt.contourf(pv_th,levels=np.arange(0.0,5.5,0.5),cmap=get_cmap('coolwarm_r'),extend='both')
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
    hour = time*6
    plt.title('300 K PV at Hour %i' %hour)
    fig.canvas.draw()
    image = np.frombuffer(fig.canvas.tostring_rgb(), dtype='uint8')
    image = image.reshape(fig.canvas.get_width_height()[::-1] + (3,))
    
    return image 

def pth_anim(ncfile, x1, x2, y1, y2, dx, time, th_lvl=300., roll=0, width=12., height=6.):
    
    matplotlib.use('Agg')
    p = getvar(ncfile,'pressure',timeidx=time)
    p_surf = np.asarray(p)[0,:,:]
    p_surf = np.roll(p_surf,roll,axis=1)
    p_surf = smooth2d(p_surf,150)
    p_surf = p_surf[y1:y2,x1:x2]
    th = getvar(ncfile,'th',timeidx=time)
    p_th = interplevel(p, th, th_lvl)
    p_th = np.roll(p_th,roll,axis=1)
    p_th = p_th[y1:y2,x1:x2]
    fig = plt.figure(figsize=(width,height))
    ax = fig.add_subplot(1,1,1)
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
    hour = time*6
    plt.title('300 K Pressure at Hour %i' %hour)
    fig.canvas.draw()
    image = np.frombuffer(fig.canvas.tostring_rgb(), dtype='uint8')
    image = image.reshape(fig.canvas.get_width_height()[::-1] + (3,))
    
    return image

def ppv_anim(ncfile, x1, x2, y1, y2, dx, time, pv_lvl=1.5, roll=0, width=12., height=6.):
    
    matplotlib.use('Agg')
    p = getvar(ncfile,'pressure',timeidx=time)
    p_surf = np.asarray(p)[0,:,:]
    p_surf = np.roll(p_surf,roll,axis=1)
    p_surf = smooth2d(p_surf,150)
    p_surf = p_surf[y1:y2,x1:x2]
    pv = getvar(ncfile,'pvo',timeidx=time)
    p_pv = interplevel(p, pv, pv_lvl)
    p_pv = np.roll(p_pv,roll,axis=1)
    p_pv = p_pv[y1:y2,x1:x2]
    fig = plt.figure(figsize=(width,height))
    ax = fig.add_subplot(1,1,1)
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
    hour = time*6
    plt.title('1.5 PVU Pressure at Hour %i' %hour)
    fig.canvas.draw()
    image = np.frombuffer(fig.canvas.tostring_rgb(), dtype='uint8')
    image = image.reshape(fig.canvas.get_width_height()[::-1] + (3,))
    
    return image

def thpv_anim(ncfile, x1, x2, y1, y2, dx, time, pv_lvl=1.5, roll=0, width=12., height=6.):
    
    matplotlib.use('Agg')
    p = getvar(ncfile,'pressure',timeidx=time)
    p_surf = np.asarray(p)[0,:,:]
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
    hour = time*6
    plt.title('1.5 PVU Theta at Hour %i' %hour)
    fig.canvas.draw()
    image = np.frombuffer(fig.canvas.tostring_rgb(), dtype='uint8')
    image = image.reshape(fig.canvas.get_width_height()[::-1] + (3,))
    
    return image