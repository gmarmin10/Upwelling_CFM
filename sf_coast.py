'''
Created on Aug 25, 2022
Specifically for the coastal model no cfm figures

'''
from pylab import * 
from matplotlib.pyplot import savefig

def sf(figName,Dpi):
    First_part="C:/Users/19046/Documents/General_Research/Coastal_Modeling/figures/eclipse_output/"
    #Second_part=savefolder+"\\"
    Figure_name=str(figName)
    Last_part=".tif"
    savefig(First_part+Figure_name+Last_part,dpi=Dpi)
