'''
Created on May 8, 2024

@author: garmin
'''
from numpy import *
from af001_energy_calculation import *

#CFM constants either from Kei's code emailed to me or supplemental table in inomura 2020
def Clim_mut(Arr):
#light and photosynthetic constants/equations
    I=1000
    E3=evalue()
    #E=E3.E
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
    CP_Plip=40

    #================================
    CG=0.563                   #GC% not CG but I started with CG so I stick with it; it does not matter as "AT GC".   [http://www.ncbi.nlm.nih.gov/genome/13522 (accessed 06/18/2016)]
    AU=1-CG
    #Values per P in mol/mol from "05 Review nucleic acid composition.xlsx"
    #RNA  
    C_CG_RNA = 19/2
    N_CG_RNA = 8/2
    
    C_AU_RNA = 19/2
    N_AU_RNA = 7/2
  
   
    YnucacidN_P = N_CG_RNA*CG + N_AU_RNA*AU #same for RNA and DNA
    YnucacidP_N = 1/YnucacidN_P
    
    YrnaC_N = (C_CG_RNA*CG + C_AU_RNA*AU)/(N_CG_RNA*CG + N_AU_RNA*AU)

    CP_RNA=YrnaC_N/YnucacidP_N

    #macromolecular constants
    Qp_RNA_min=0.00022282635312159368#2.23e-4             #got this from inomura 2020 supplemental
    Qc_DNA=0.0009414247529173031#9.41e-4                            #got this from inomura 2020 supplemental

    Qc_oth=0.01821441041892576#1.51786753491048E-15/Qc     #Cessential in CFM
    Qc_pro_oth=0.239947521088736#Nconst_protein/Qc*(1/NC_pro)
    Qc_other=Qc_pro_oth+Qc_DNA+Qc_oth
    
    #Climitation
    a=CP_RNA*AP_RNA*((Apho*Achl)+Abio)
    b=Achl*(1+Apho+CP_Plip*A_pchl_pho)+CP_RNA*AP_RNA*(Apho*Bchl+Qc_pro_oth)+Abio
    c=Bchl*(1+Apho+CP_Plip*A_pchl_pho)+CP_RNA*Qp_RNA_min+Qc_other-1
    
    #Nlimitation
    #a=NP_RNA*AP_RNA*((Apho*Achl)+Abio)
    #b=NC_Chl*Achl+NC_pro*(Apho*Achl+Abio)+NP_RNA*AP_RNA*(Apho*Bchl+Qc_pro_oth)
    #c=NC_Chl*Bchl+NC_pro*(Apho*Bchl+Qc_pro_oth)+NC_DNA*Qc_DNA+NP_RNA*Qp_RNA_min-Qn
    
    mu_max=(-b+(b**2-(4*a*c))**(1/2))/(2*a)
    
    return(mu_max)


Ea=70000
R=8.3
Tref=293
Te=(25+273)        #temperature of water column moving away from the coast
A=Ea/R
Arr=exp(-A*((1/Te)-(1/Tref)))


#print(Clim_mut(Arr))









