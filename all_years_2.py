'''
Updated on July 8, 2024

Incorporating cellular equations into coastal model
adding allocation
the excess N problem fixed 

@author: garmin
'''
#import necessary packages
from numpy import *
from matplotlib.pyplot import plot, xlim, ylim, show, xlabel, ylabel,subplot, figure, tight_layout,stackplot,title,legend,scatter,xticks, yticks
from math import *
from math import pi
from time import *
from pylab import *
from Climitation_max_growth_rate_with_temp import *
from Nlimitation_max_growth_rate_with_temp import *
from allocation_fxn_with_t import *
from allocation_fxn_excessN_with_t import *
from matplotlib.lines import *
from sf_coast import *

def coast_sim(Lx,U0,Te1,Te2,Te_step):
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
    Lx=Lx           #length in meters
    dx=1e3           #grid spacing
    Nx=(Lx/dx)+1      #number of grid points
    Nx=int(Nx)
    xx=arange(0,Lx+dx,dx)  #create a space array
    
    #wind forcing
    U0=U0     #velocity (m/s) subject to change
    
    #initial conditions of nutrient variables NO3(Nnut) and PO4 (Pnut)
    Nnut=2*ones(Nx)*15        #assume NO3 concentration higher than PO4. Change based on field observations
    Nnut_i=Nnut[0]     #intial value of NO3 available at the coastal boundary
    
    #initial condition of biomass variables-to be replaced with CFM later
    Nbio=0.01*Nnut
    Nbio_i=0.01*Nnut_i    #phytoplankton N at coast (boundary condition)
    Pbio=Nbio*(1/16)
    Cbio=Nbio*(106/16)
    Cbio_i=Nbio_i*(106/16)
    
    #temperature inputs
    Ea=70000
    R=8.3
    Tref=293
    A=Ea/R

    Te=arange(Te1,Te2,Te_step)
    Arr=exp(-A*((1/Te[0])-(1/Tref)))
    Arr_a=zeros(size(Nnut))
    Arr_a[0]=Arr

    Qn=Nbio/Cbio
    groN=(Nlim_mut(Qn,Arr))/sperd

    Qc=1.00*10**(-12)/12    #molC/cell
    aNO3= 3e-20/Qc      #from Monod paper supplementary material- just a starting point-can be changed
    aPO4=2.31e-15/Qc
    aNO3=aNO3*2
    
    phi=0.5     #fraction of uptake remineralized locally
    
    
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
       
    for n in it_Nt:    
        f=A0*(sin(w*n/spery))
        Ua=U0*f
        U=U0+Ua
    
    #starting to incorporate CFM principles here
        Qn=Nbio/Cbio
        #Qfe=Febio/Cbio
            
        groN=(Nlim_mut(Qn,Arr))/sperd
        groC=(Clim_mut(Arr))/sperd
        #groFe=(Felim_mu(Qfe))/sperd
        gro=groN
        gro[groN>groC]=groC
        mc=gro*0.2
        
        Vn=(aNO3*Nnut)#-(Qn_ex/sperd)
        #Vfe=afe*fenut
         
        RgrowN=Vn*Cbio
        #Rgrowfe=Vfe*Cbio
        RgrowC=gro*Cbio
        
        RmortC=mc*Cbio**2
        RmortN=(mc*Cbio**2)*Qn
        #Rmortfe=(mc*Cbio**2)*Qfe
        
        RbioN=RgrowN-RmortN
        RnutN=-RgrowN+phi*RmortN#+Qn_ex
        #Rbiofe=Rgrowfe-Rmortfe
        #Rnutfe=-Rgrowfe+phi*Rmortfe
        RbioC=RgrowC-RmortC
    
        
        #update the distribution: Advection scheme
        for i in it_Nx:
            Nnut[i+1]=((dt/dx)*U*Nnut[i]+Nnut[i+1]+RnutN[i]*dt)/(1+dt/dx*U)
            #Pnut[i+1]=((dt/dx)*U*Pnut[i]+Pnut[i+1]+RnutP[i]*dt)/(1+dt/dx*U)
                  
            Nbio[i+1]=((dt/dx)*U*Nbio[i]+Nbio[i+1]+RbioN[i]*dt)/(1+dt/dx*U)
            #Pbio[i+1]=((dt/dx)*U*Pbio[i]+Pbio[i+1]+RbioP[i]*dt)/(1+dt/dx*U)
            Cbio[i+1]=((dt/dx)*U*Cbio[i]+Cbio[i+1]+RbioC[i]*dt)/(1+dt/dx*U)
            
            Arr=exp(-A*((1/Te[i])-(1/Tref)))
            
            Arr_a[i+1]=Arr
    
        Qn_max=allo_ext(Arr)
         
        Qn_ex=Qn-Qn_max
         
        Qn[Qn_ex>0]=Qn_max
         
        Qn_ex[Qn_ex<0]=0
        exc=Qn_ex*Cbio
        Nnut=Nnut+exc
    
        Nbio=Nbio-exc
        
         
    return Nbio,Nnut,gro,Qn,Te,Cbio,Lx,dx,Arr


neu_Nbio,neu_Nnut,neu_gro,neu_Qn,neu_Te,neu_Cbio,neu_Lx,neu_dx,neu_Arr=coast_sim(8e5,0.1,283,286.0037,0.0037)
el_Nbio,el_Nnut,el_gro,el_Qn,el_Te,el_Cbio,el_Lx,el_dx,el_Arr=coast_sim(8e5,0.05,285,287.0025,0.0025)
la_Nbio,la_Nnut,la_gro,la_Qn,la_Te,la_Cbio,la_Lx,la_dx,la_Arr=coast_sim(8e5,0.15,281,283.0025,0.0025)


def mac_all_sim(gro,Qn,loop_size,Arr):
     
    loop=arange(0,loop_size,1)
      
    Qc_oth= []
    Qc_chl = []
    Qn_chl = []
    Qc_pro_pho = []
    Qc_pro_bio=[]
    Qc_pro= []
    Qn_pro = []
    Qc_RNA = []
    Qn_RNA = []
    Qn_DNA= []
    Qc_nuc = []
    Qc_thy = []
    Qn_sto = []
    Qc_nsto=[]
    Qc_csto= []
     
     
     
    #Qc_oth, Qc_chl, Qn_chl, Qc_pro_pho, Qc_pro_bio,Qc_pro,Qn_pro, Qc_RNA, Qn_RNA, Qn_DNA,Qc_nuc, Qc_thy, Qn_sto, Qc_nsto, Qc_csto
    for i in loop:
        al_data=allo_t(gro[i],Qn[i],Arr)
        Qc_oth.append(al_data[0])
        Qc_chl.append(al_data[1])
        Qn_chl.append(al_data[2])
        Qc_pro_pho.append(al_data[3])
        Qc_pro_bio.append(al_data[4])
        Qc_pro.append(al_data[5])
        Qn_pro.append(al_data[6])
        Qc_RNA.append(al_data[7])
        Qn_RNA.append(al_data[8])
        Qn_DNA.append(al_data[9])
        Qc_nuc.append(al_data[10])
        Qc_thy.append(al_data[11])
        Qn_sto.append(al_data[12])
        Qc_nsto.append(al_data[13])
        Qc_csto.append(al_data[14])
      
    #carbon allocation   
    Qc_oth_plot=100*array([Qc_oth])
    Qc_chl_plot=100*array([Qc_chl])
    Qc_pro_pho_plot=100*array([Qc_pro_pho])
    Qc_pro_bio_plot=100*array([Qc_pro_bio])
    Qc_pro_plot=100*array([Qc_pro])
    Qc_pro_oth_plot=Qc_pro_plot-Qc_pro_bio_plot-Qc_pro_pho_plot
    Qc_RNA_plot=100*array([Qc_RNA])
    Qc_nuc_plot=100*array([Qc_nuc])
    Qc_DNA_plot=Qc_nuc_plot-Qc_RNA_plot
    Qc_thy_plot=100*array([Qc_thy])
    Qc_nsto_plot=100*array([Qc_nsto])
    Qc_csto_plot=100*array([Qc_csto])
      
    Photo_plot=Qc_chl_plot+Qc_pro_pho_plot+Qc_thy_plot
    Bio_plot=Qc_pro_bio_plot+Qc_RNA_plot
    cStorage_plot=Qc_csto_plot
    Essential_plot=Qc_oth_plot+Qc_DNA_plot+Qc_pro_oth_plot
    nstorage_plot=Qc_nsto_plot
       
     
    #nitrogen allocation
    Qn_chl_plot=array(Qn_chl)
    Qn_pro_plot=array(Qn_pro)
    Qn_RNA_plot=array(Qn_RNA)
    Qn_DNA_plot=array(Qn_DNA)
    Qn_sto_plot=array(Qn_sto)
    Qn_pro_pho_plot=(1/3.82)*array(Qc_pro_pho)
    Qn_pro_bio_plot=(1/3.82)*array(Qc_pro_bio)
    Qn_pro_oth_plot=Qn_pro_plot-Qn_pro_bio_plot-Qn_pro_pho_plot
     
    n_pho=Qn_pro_pho_plot+Qn_chl_plot
    n_bio=Qn_RNA_plot+Qn_pro_bio_plot
    n_ess=Qn_DNA_plot+Qn_pro_oth_plot
    n_sto=Qn_sto_plot
     
    return Qc_chl_plot,Photo_plot,Bio_plot,cStorage_plot,Essential_plot,nstorage_plot,n_pho,n_bio,n_ess,n_sto
  
  
neu_Qc_chl_plot,neu_Photo_plot,neu_Bio_plot,neu_cStorage_plot,neu_Essential_plot,neu_nstorage_plot,neu_n_pho,neu_n_bio,neu_n_ess,neu_n_sto=mac_all_sim(neu_gro,neu_Qn,801,neu_Arr)
el_Qc_chl_plot,el_Photo_plot,el_Bio_plot,el_cStorage_plot,el_Essential_plot,el_nstorage_plot,el_n_pho,el_n_bio,el_n_ess,el_n_sto=mac_all_sim(el_gro,el_Qn,801,el_Arr)
la_Qc_chl_plot,la_Photo_plot,la_Bio_plot,la_cStorage_plot,la_Essential_plot,la_nstorage_plot,la_n_pho,la_n_bio,la_n_ess,la_n_sto=mac_all_sim(la_gro,la_Qn,801,la_Arr)
   
ne_leg=Line2D([], [], color='#DDDDDD',linewidth=3, marker='o', label='Neutral')
el_leg=Line2D([], [], color='#FFCCCC', linewidth=3,marker='v',label='El Ni$\~n$o')
la_leg=Line2D([], [], color='#77AADD', linewidth=3,marker='p', label='La Ni$\~n$a')
han=[ne_leg,el_leg,la_leg]


#savetxt("fo.csv", neu_Nbio, delimiter=",")

# 
# # # #plotting
# ax=figure(1,figsize=(8,6.5))
x1=arange(0,neu_Lx+neu_dx,neu_dx)
x1=x1*1e-3
# x2=arange(0,el_Lx+el_dx,el_dx)
# x2=x2*1e-3
# x3=arange(0,la_Lx+la_dx,la_dx)
# x3=x3*1e-3
# subplot(2,1,1)
# plot(x1,neu_Nnut,marker='o',color='#DDDDDD',label='Neutral')
# plot(x2,el_Nnut,marker='v',color='#FFCCCC',label='El Ni$\~n$o')
# plot(x3,la_Nnut,marker='p',color='#77AADD',label='La Ni$\~n$a')
# ylim(0,40)
# legend(handles=han,loc='upper center', bbox_to_anchor=(0.5,-0.2),ncol=3,frameon=False,fontsize='large')
# ylabel('Water Column N (\u00B5M)',fontsize=14)
#    
# subplot(2,1,2)
# plot(x1,neu_Nbio,marker='o',color='#DDDDDD')
# plot(x2,el_Nbio,marker='v',color='#FFCCCC')
# plot(x3,la_Nbio,marker='p',color='#77AADD')
# ylabel('Biomass N (\u00B5M)',fontsize=14) 
# ylim(0,4)
# xlabel('Distance from coast (km)',fontsize=16)
# tight_layout()
# sf('wat_bio_01',500)
#  
#   
# figure(2,figsize=(8,6.5))
# plot(x1,neu_gro,marker='o',color='#DDDDDD')
# plot(x2,el_gro,marker='v',color='#FFCCCC')
# plot(x3,la_gro,marker='p',color='#77AADD')
# legend(handles=han,loc='upper right', frameon=False,fontsize='large')
# ylim(0,1.4e-5)
# xlim(0)
# ylabel('\u00B5 (s$^{-1}$)',fontsize=16)
# xlabel('Distance from coast (km)',fontsize=16)
# sf('growth',500)
# # #  
# bio_pat=Line2D([], [], color='#44AA99',linewidth=5, label='Bio')
# pho_pat=Line2D([], [], color='#CC6677', linewidth=5,label='Photo')
# ess_pat=Line2D([], [], color='#AA4499', linewidth=5, label='Ess')
# csto_pat=Line2D([], [], color='#DDCC77', linewidth=5, label='C. sto')
# hand=[ess_pat,pho_pat,bio_pat,csto_pat]
# 
# figure(3,figsize=(8,6.5))
# plot(x3,la_Qn,marker='p',color='#77AADD')
# plot(x1,neu_Qn,marker='o',color='#DDDDDD')
# plot(x2,el_Qn,marker='v',color='#FFCCCC')
# legend(handles=han,loc='upper right', frameon=False)
# ylim(0,0.30)
# xlim(0)
# xlabel('Distance from coast (km)')
# ylabel('N:C')
# #legend(handles=hand,frameon=False,loc='lower center',ncol=4)
# sf('Qn',500)
# # 
# # # #  
stackc=['#AA4499','#882255','#CC6677','#44AA99','#DDCC77']
figure(4,figsize=(8,6.5))
stackplot(x1,neu_Essential_plot,neu_nstorage_plot,neu_Photo_plot,neu_Bio_plot,neu_cStorage_plot,colors=stackc)
#xlabel('Distance from coast (km)',fontsize=16)
ylabel('C Allocation (%)',fontsize=16)
xticks(fontsize=14)
yticks(fontsize=14)
title('Neutral',fontsize=18)
# sf('Neu_C',500)
# 
# 
# 
# figure(5,figsize=(8,6.5))
# stackplot(x2,el_Essential_plot,el_nstorage_plot,el_Photo_plot,el_Bio_plot,el_cStorage_plot,colors=stackc)
# #xlabel('Distance from coast (km)',fontsize=16)
# title('El Ni$\~n$o',fontsize=18)
# xticks(fontsize=14)
# yticks(fontsize=14)
# sf('El_C',500)
#   
# figure(6,figsize=(8,6.5))
# stackplot(x3,la_Essential_plot,la_nstorage_plot,la_Photo_plot,la_Bio_plot,la_cStorage_plot,colors=stackc)
# #xlabel('Distance from coast (km)',fontsize=16)
# xticks(fontsize=14)
# yticks(fontsize=14)
# title('La Ni$\~n$a',fontsize=18)
# sf('La_C',500)
# # 
# stackn=['#AA4499','#CC6677','#44AA99','#882255']
# 
# figure(7,figsize=(8,6.5))
# stackplot(x1,neu_n_ess,neu_n_pho,neu_n_bio,neu_n_sto,colors=stackn)
# xlabel('Distance from coast (km)',fontsize=16)
# ylabel('N Allocation (mol N mol C$^{-1}$)',fontsize=16)
# xticks(fontsize=14)
# yticks(fontsize=14)
# ylim(0,0.30)
# xlim(0)
# #title('Neutral',fontsize=18)
# sf('neu_NC',500)
#  
# figure(8,figsize=(8,6.5))
# stackplot(x2,el_n_ess,el_n_pho,el_n_bio,el_n_sto,colors=stackn)
# xlabel('Distance from coast (km)',fontsize=16)
# xticks(fontsize=14)
# yticks(fontsize=14)
# ylim(0,0.30)
# xlim(0)
# #title('El Ni$\~n$o',fontsize=18)
# sf('el_NC',500)
#    
# figure(9,figsize=(8,6.5))
# stackplot(x3,la_n_ess,la_n_pho,la_n_bio,la_n_sto,colors=stackn)
# xlabel('Distance from coast (km)',fontsize=16)
# xticks(fontsize=14)
# yticks(fontsize=14)
# ylim(0,0.30)
# xlim(0)
# #title('La Ni$\~n$a',fontsize=18)
# sf('la_NC',500)
# 
# 
show()


print(neu_Bio_plot)
print(el_Bio_plot)
print(la_Bio_plot)







