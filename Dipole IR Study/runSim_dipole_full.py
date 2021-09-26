from matplotlib import pyplot as plt
from LHeC_shieldingStudy import dipoleOptimised_full as lf
import pybdsim
import numpy as np
import pandas as pd
import os
import csv

os.chdir("/Users/connormonaghan/Documents/LHeC/Final_Sims")

pybdsim.Data.LoadROOTLibraries()
plt.style.use(['science','no-latex','ieee', 'grid'])

material    = 'Pb'
runKey      = 'Dipole_full'

thicknesses = [0.021,0.022,0.025,0.03,0.04,0.06,0.08]
ngenerate   = 10000
nruns       = 30

percentAtten = np.zeros((len(thicknesses),12))

buffer = np.zeros((len(thicknesses), nruns))

i = 0
for t in thicknesses:
    s = lf.shieldingStudy(material, ngenerate, nruns, t, runKey)

    s.genGMAD()

    value, err, val_range = s.runStudy()
    totalPhotons, err_num = s.getTotalPhotons()
    ZpCut, err_zp         = s.getZpCut()
    eNum, pNum, err_e, err_p = s.getTotalAper()

    buffer[i,:] = s.getBuffer()

    percentAtten[i,:] = [value, err, val_range[0], val_range[1], 
                        totalPhotons, err_num, 
                        eNum, err_e, pNum, err_p,
                        ZpCut, err_zp]
    i = i +1 

# Quick plot to see trend in data
fig = plt.figure(1)
plt.errorbar(percentAtten[:,0], thicknesses, xerr=percentAtten[:,1])

fig.savefig("Plot-Outputs/{}_{}.pdf".format(runKey, material), format='pdf')

df = pd.DataFrame(data=percentAtten, index=thicknesses, columns=["Frac. survival", "Frac. err", "Frac. min", "Frac. max",
                                                                "Tot. w/ cut", "Tot. w/ cut err",
                                                                "Tot. e Aper", "Tot. e Aper err", "Tot. p Aper", "Tot. p Aper err",
                                                                "Tot. w/o zp", "Tot. w/o zp err"])

df.to_csv('DATA/{}/{}_runData.csv'.format(runKey, material), float_format='%f')

print(df)

buff = pd.DataFrame(data=buffer)
buff.to_csv('DATA/{}/{}_runBuffer.csv'.format(runKey, material), float_format='%f')