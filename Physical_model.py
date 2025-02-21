'''
Created on Aug 25, 2022

Conversion of physical model from matlab to python

@author: garmin
'''

#import necessary packages
from numpy import *
from matplotlib.pyplot import plot, xlim, ylim, show, xlabel, ylabel,subplot, figure, tight_layout
from pandas import *
from math import *
from math import pi
from time import *
from pylab import *
from sf_coast import *

#create the time grid
sperd=60*60*24   #seconds per day
spery=365*sperd  #seconds per year
nyr=1            #number of years to simulate. This can be chagned by the user
T=spery*nyr      #length of time in seconds to simulate
dt=spery/720      #time step
t=arange(0,T+dt,dt) #create a time array 
Nt=T/dt+1         #number of time grid points
Nt=int(Nt)

#create the spatial grid
Lx=5e5           #length in meters
dx=1e3           #grid spacing
Nx=(Lx/dx)+1      #number of grid points
Nx=int(Nx)
xx=arange(0,Lx+dx,dx)  #create a space array

#wind forcing
U0=0.1     #velocity (m/s) subject to change

#initial conditions of nutrient variables NO3(Nnut) and PO4 (Pnut)
Pnut=2*ones(Nx)  #creating a ones array the size of the number of spatial grid points
Pnut_i=Pnut[0]      #initial value of phosphorus at the coastal boundary
Rup_Po=10*Pnut_i/spery #baseline P uptake rate (will be changed with CFM-Phyto)
Nnut=Pnut*15        #assume NO3 concentration higher than PO4. Change based on field observations
Nnut_i=Pnut_i*15     #intial value of NO3 available at the coastal boundary
Rup_No=10*Nnut_i/spery  #baseline N uptake rate (will be changed with CFM-Phyto)

#initial condition of biomass variables-to be replaced with CFM later
Pbio=0.01*Pnut
Pbio_i=0.01*Pnut_i    #phytoplankton P at coast (boundary condition)
Nbio=0.01*Nnut
Nbio_i=0.01*Nnut_i    #phytoplankton N at coast (boundary condition)
print(size(Nbio))

#initial biological parameters
Kp=0.1   #half-saturation constant for Pnut
Kn=1     #half-saturation constant for Nnut
mu= 1/sperd  #growth rate per sec
phi=0.5     #fraction of uptake remineralized locally
Snp=16     #redfield N:P ratio of phytoplankton
m2=mu*0.2  #quadratic mortality

#Period, amplitude, frequency
Period=1   #period of oscillation in forcing (velocity) (yr)
w=(2*pi)/Period  #frequency of oscillation
A0=0.5   #amplitude of oscillation
nn=0    #year counter

#model equations
it_Nx=arange(0,Nx-1,1)
it_Nt=arange(0,Nt+1,1)
f=[]
Ua=[]
    #vary the circulation rates
for y in t:
    f=A0*(sin(w*y/spery))
    #fn.append(f)
    U0_array=full_like(f,U0)
    Ua=U0*f
    U=U0+Ua

for n in it_Nt:
    #calculate the biological rates-to be replaced by CFM
    RgrowN=mu*Nbio*(Nnut/(Nnut+Kn))
    RmortN=m2*Nbio**2
    RbioN=RgrowN-RmortN
    RnutN=-RgrowN+phi*RmortN
    RbioP=RbioN/Snp
    RnutP=RnutN/Snp
    
    #update the distribution: Advection scheme
    for i in it_Nx:
        Pnut[i+1]=((dt/dx)*U*Pnut[i]+Pnut[i+1]+RnutP[i]*dt)/(1+dt/dx*U)
        Nnut[i+1]=((dt/dx)*U*Nnut[i]+Nnut[i+1]+RnutN[i]*dt)/(1+dt/dx*U)
        
        Pbio[i+1]=((dt/dx)*U*Pbio[i]+Pbio[i+1]+RbioP[i]*dt)/(1+dt/dx*U)
        Nbio[i+1]=((dt/dx)*U*Nbio[i]+Nbio[i+1]+RbioN[i]*dt)/(1+dt/dx*U)

#plotting
ax=figure(1)
subplot(2,2,1)
x=arange(0,Lx+dx,dx)
x=x*1e-3
plot(x,Pnut,marker='o',color='orange')
xlabel('horizontal distance (km)')
ylabel('PO4 (uM)')
subplot(2,2,2)
plot(x,Nnut,marker='o',color='green')
xlabel('horizontal distance (km)')
ylabel('NO3 (uM)')
subplot(2,2,3)
plot(x,Pbio,marker='o',color='red')
xlabel('horizontal distance (km)')
ylabel('Phyto P (uM)')
subplot(2,2,4)
plot(x,Nbio,marker='o',color='blue')
xlabel('horizontal distance (km)')
ylabel('Phyto N (uM)') 
tight_layout()
sf('Nutrient_Concentrations_nocfm',500)


show()














