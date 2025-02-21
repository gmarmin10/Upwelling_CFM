'''
Created on May 8, 2024

@author: garmin
'''
from numpy import *
from af001_energy_calculation import *
import csv

#CFM constants either from Kei's code emailed to me or supplemental table in inomura 2020
def Felim_mut(Qfe,Arr):
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
    #stoichiometric constants
    NC_Chl=4/55
    NC_pro=1/4.49

    #================================
    #Stoichiometric parameters
    #================================
    CG=0.563                   #GC% not CG but I started with CG so I stick with it; it does not matter as "AT GC".   [http://www.ncbi.nlm.nih.gov/genome/13522 (accessed 06/18/2016)]
    YnucacidP_N=1/(3.5*(1-CG)+4*CG)               #(molP molN-1) P/N molar ratio of RNA (193-26) values (193-28) excel file "08 N to P ratio in DNA and RNA.xlsx"
    YdnaC_N=3.5*(1-CG)+2.5*CG       #(molC molN-1) C/N molar ratio of dna (based on "10 Amino acids and nucleic acids stoichiometry.xlsx)
    
    NP_RNA=1/YnucacidP_N #3.8
    NC_DNA=1/YdnaC_N
    
    #macromolecular constants
    Qp_RNA_min=0.00022282635312159368#2.23e-4             #got this from inomura 2020 supplemental
    Qc_DNA=0.0009414247529173031#9.41e-4                            #got this from inomura 2020 supplemental
    Qc_pro_oth=0.239947521088736#Nconst_protein/Qc*(1/NC_pro)

    Afe_pho=5.83e-3 #molFe(molCinChl)-1 found in inomura supplemental table from inomura 2022

    a=Qfe/(Afe_pho*Achl)
    b=Bchl/Achl

    mu_max=(a-b)

    return(mu_max)




