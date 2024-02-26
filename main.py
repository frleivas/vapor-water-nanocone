#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  1 15:52:03 2020

@author: frleivas
"""

import pandas as pd
import numpy as np
import math as mt
import matplotlib.pyplot as plt
import matplotlib.style as style 
from random import seed
from random import random
seed(8)

style.use('seaborn-poster')

dr = 0.142                             #Atoms Distance
phi = 60*mt.pi/180                     #Graphene atoms angle
apx = 19.2*mt.pi/180                   #Apex angle nanocone

Lg = 5                                 #Graphene sheet size

no = 2197                              #Oxygen Atoms
n = int((Lg+0.1)/(2*dr*mt.sin(phi)))
m = int((Lg+0.1)/(dr*mt.cos(phi)+2*dr))
N = 4*m*n                              #Number of atoms in the sheet

# Builds graphene sheet by lines


x = []
y = []
dy = (2*dr*mt.cos(phi)+2*dr)
drsin = dr*mt.sin(phi)
drcos = dr*mt.cos(phi)
for j in range (0, m):
    for i in range (0,n):
        x.append(2*drsin*i)
        y.append(j*dy)
    for i in range (0,n):
        x.append(x[i] + drsin)
        y.append(j*dy + drcos)
    for i in range (0,n):
        x.append(x[i] + drsin)
        y.append(j*dy + (drcos+dr))
    for i in range (0,n):
        x.append(x[i])
        y.append(j*dy + (2*drcos+dr))

#Creating coordinates of the two graphene walls
xs = []; ys = []; zs = []                   #1st
xsh = []; ysh = []; zsh = []                #2nt

h = 4.1*mt.sin(apx/2)/mt.cos(apx/2)         #Size of the hole on the graphene wall in constact with nanocone
for i in range (0,N):
        xs.append(x[i])
        ys.append(y[i])
        zs.append(-0.3)
        if (mt.sqrt((x[i]-Lg/2)*(x[i]-Lg/2) + (y[i]-Lg/2)*(y[i]-Lg/2)) > h ):
            xsh.append(x[i])
            ysh.append(y[i])
            zsh.append(Lg) 
            
# Creating molecules water coordinates - Equally Spaced

theta0 = 104.52                           #Angle H-O-H
rOH = 0.09572                             #Raio R-O, nanometros

dist = (Lg)/mt.pow(no+1, 1./3.);          #Distance of the wall to not create atoms over Graphene sheet
distx = (Lg-1)/mt.pow(no+1, 1./3.);  

#Oxygen positions
f=0 ; x_O = []; y_O = [] ; z_O = [] 
for i in range (0,(int)(mt.pow(no+1, 1./3.))): 
    for j in range (0,(int)(mt.pow(no+1, 1./3.))):
        for k in range (0,(int)(mt.pow(no+1, 1./3.))):   
            x_O.append(i*distx+random()/10)
            y_O.append(j*dist+random()/10)
            z_O.append(k*dist+random()/10)
            f+=1
            if (f==no): break
            
#Hydrogen1 positions	  
x_H = []; y_H = [] ; z_H = []; x_H2 = []; y_H2 = [] ; z_H2 = []
for i in range (0,no):		
    a = 2*random()*mt.pi; b = random()*mt.pi;
    x_H.append(x_O[i] + rOH*mt.cos(a)*mt.sin(b))
    y_H.append(y_O[i] + rOH*mt.sin(a)*mt.sin(b))
    z_H.append(z_O[i] + rOH*mt.cos(b))
#Hydrogen2 positions		
    c = random()*mt.pi; s = (theta0)*(mt.pi/180); d = mt.pi/2 - a;		
    x_H2.append(x_O[i] + mt.cos(d)*rOH*mt.sin(s)*mt.cos(c) + mt.sin(d)*mt.cos(b)*rOH*mt.sin(s)*mt.sin(c) + mt.sin(d)*mt.sin(b)*rOH*mt.cos(s)) 
    y_H2.append(y_O[i] - mt.sin(d)*rOH*mt.sin(s)*mt.cos(c) + mt.cos(d)*mt.cos(b)*rOH*mt.sin(s)*mt.sin(c) + mt.sin(b)*mt.cos(d)*rOH*mt.cos(s)) 
    z_H2.append(z_O[i] - mt.sin(b)*rOH*mt.sin(s)*mt.sin(c) + mt.cos(b)*rOH*mt.cos(s))
	
#PBC on x and z directions
for i in range (0,no):                 
		if (y_O[i]<=0):   y_O[i]=Lg+y_O[i];
		if (z_O[i]<=0):   z_O[i]=Lg+z_O[i];
		if (y_O[i]>=Lg):  y_O[i]=y_O[i]-Lg;
		if (z_O[i]>=Lg):  z_O[i]=z_O[i]-Lg;	
		if (y_H[i]<=0):   y_H[i]=Lg+y_H[i];
		if (z_H[i]<=0):   z_H[i]=Lg+z_H[i];
		if (y_H[i]>=Lg):  y_H[i]=y_H[i]-Lg;
		if (z_H[i]>=Lg):  z_H[i]=z_H[i]-Lg;
		if (y_H2[i]<=0):   y_H2[i]=Lg+y_H2[i];
		if (z_H2[i]<=0):   z_H2[i]=Lg+z_H2[i];
		if (y_H2[i]>=Lg):  y_H2[i]=y_H2[i]-Lg;
		if (z_H2[i]>=Lg):  z_H2[i]=z_H2[i]-Lg;
		        

#Coordenadas of the nanocone
f=pd.read_table("nanocone-p.xyz",  sep = " ")

#Plotting Figure

# Plotting Water molecules
fig = plt.figure()
ax = plt.axes(projection='3d')
for i in range (0,no):         
    if (random()<0.1):                   #Not plor all the atoms to facilitate visualization 
            ax.scatter(x_O[i], y_O[i],z_O[i],  c='red', s=30);
            ax.scatter(x_H[i], y_H[i],z_H[i],  c='blue',s=15)
            ax.scatter(x_H2[i], y_H2[i],z_H2[i], c='blue', s=15)

# Plotting Graphene Walls
for i in range (0,len(xs)):
    for j in range (i,len(xs)):
        xi = xs[i] ; yi = ys[i] ; zi = zs[i]; xj = xs[j] ; yj = ys[j] ; zj = zs[j] 
        if(abs(mt.sqrt((xi-xj)*(xi-xj)+(yi-yj)*(yi-yj)+(zi-zj)*(zi-zj))-dr) < 0.01): 
            ax.plot3D([zi, zj], [xi, xj], [yi, yj], 'gray')
 
for i in range (0,len(xsh)):
    for j in range (i,len(xsh)):
        xi = xsh[i] ; yi = ysh[i] ; zi = zsh[i]; xj = xsh[j] ; yj = ysh[j] ; zj = zsh[j] 
        if(abs(mt.sqrt((xi-xj)*(xi-xj)+(yi-yj)*(yi-yj)+(zi-zj)*(zi-zj))-dr) < 0.01): 
            ax.plot3D([zi, zj], [xi, xj], [yi, yj], 'gray')
            
# Plotting nanocone,adjusting the scale and changing the position to the middle of the carbon sheet
Nf = len(f) 
for i in range (0,Nf):
    for j in range (0,Nf):
        xi = f.x[i]/10 +2.5 ; yi = f.y[i]/10 +2.5 ; zi = f.z[i]/10 + 10; xj = f.x[j]/10 + 2.5 ; yj = f.y[j]/10 + 2.5 ; zj = f.z[j]/10 + 10
        if(abs(mt.sqrt((xi-xj)*(xi-xj)+(yi-yj)*(yi-yj)+(zi-zj)*(zi-zj))-dr) < 0.01):
            ax.plot3D([zi, zj], [xi, xj], [yi, yj], 'cornflowerblue')
