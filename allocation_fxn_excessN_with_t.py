'''
Created on May 29, 2024

making a function that calculates the allocation 

@author: garmin

'''

from numpy import *
from matplotlib.pyplot import plot, xlim, ylim, show, xlabel, ylabel,subplot, figure, tight_layout
from math import *
from pylab import *
from Climitation_max_growth_rate_with_temp import *

def allo_ext(Arr):
    
    mu=Clim_mut(Arr)
    I=1000

    E=0.7741553538213819
    
    #m=3.93e-1 #3.79146798299876E-19                                #3.79146798299876E-19
    m=0.3930994004773114
    vIm=276.923478810005 #2.77e2 #0.00320513285659728                          #0.00320513285659728        #maximum photosynthetic rate called Pmax in OG CFM code
    AI= 0.00863364097132997                           #0.00863364097132997         #called OT in CFM
    vI=vIm*(1-exp(-AI*I))
    vI=vIm 
    Achl=(1+E)/vI
    Bchl=m/vI
     
    #proportionality constants
    #Apho= 1.60e1 # 3.56099164557551          #1.60e1  
    Apho=15.988852488634041     
    #Abio=  2.71e-1 # 4.34728279914354E-10   # 2.71e-1  
    Abio=0.2711013856688124
    Abio=Abio/Arr
    #AP_RNA= 4.23e-3#6212.59249917364     # 4.23e-3
    AP_RNA=0.004234953828227311
    AP_RNA=AP_RNA/Arr
    A_pchl_pho= 0.0281633095303638 
    
    #stoichiometric constants
    NC_Chl=4/55
    NC_pro=1/3.82
    
    #================================
    #Stoichiometric parameters
    #================================
    CG=0.563                   #GC% not CG but I started with CG so I stick with it; it does not matter as "AT GC".   [http://www.ncbi.nlm.nih.gov/genome/13522 (accessed 06/18/2016)]

    AT=1-CG
    AU=1-CG
    #Values per P in mol/mol from "05 Review nucleic acid composition.xlsx"
    #RNA  
    C_CG_RNA = 19/2
    N_CG_RNA = 8/2
    
    C_AU_RNA = 19/2
    N_AU_RNA = 7/2
  
    #DNA  
    C_CG_DNA = 19/2
    N_CG_DNA = 8/2
    
    C_AT_DNA = 20/2
    N_AT_DNA = 7/2

   
    YdnaC_N = (C_CG_DNA*CG + C_AT_DNA*AT)/(N_CG_DNA*CG + N_AT_DNA*AT)
    YrnaC_N = (C_CG_RNA*CG + C_AU_RNA*AU)/(N_CG_RNA*CG + N_AU_RNA*AU)
   
    YnucacidN_P = N_CG_RNA*CG + N_AU_RNA*AU #same for RNA and DNA
    YnucacidP_N = 1/YnucacidN_P
    
    CP_RNA=YrnaC_N/YnucacidP_N
        
    NP_RNA=1/YnucacidP_N #3.8
    NC_DNA=1/YdnaC_N
    CP_Plip=40
    Qc=1.00*10**(-12)/12
    #macromolecular constants
    Qp_RNA_min=0.00022282635312159368#2.23e-4             #got this from inomura 2020 supplemental
    Qc_DNA=0.0009414247529173031#9.41e-4                            #got this from inomura 2020 supplemental
    Qc_pro_oth=0.239947521088736#Nconst_protein/Qc*(1/NC_pro)
    
    Qc_oth=0.01821441041892576
    Qc_chl=(Achl*mu)+Bchl
    Qc_pro_pho=Apho*Qc_chl
    Qc_pro_bio=Abio*mu
    Qc_pro=Qc_pro_pho+Qc_pro_bio+Qc_pro_oth
    Qc_RNA=CP_RNA*(AP_RNA*mu*Qc_pro+Qp_RNA_min)
    Qc_nuc=Qc_RNA+Qc_DNA
    Qc_thy=CP_Plip*A_pchl_pho*Qc_chl
    
    Qn_RNA=NP_RNA*(AP_RNA*mu*Qc_pro+Qp_RNA_min)#/Qc
    Qn_DNA=Qc_DNA*NC_DNA#/Qc
    Qn_pro=Qc_pro*NC_pro#/Qc
    Qn_chl=Qc_chl*NC_Chl#/Qc
    
    Qn_sto=3.5e-2 #molN/molC maximum N storage from Inomura 2022

    
    Qn_cell=Qn_RNA+Qn_DNA+Qn_pro+Qn_chl#+Qn_sto
    
    #ex=Qn-Qn_cell
    
    return Qn_cell
    
    
    
