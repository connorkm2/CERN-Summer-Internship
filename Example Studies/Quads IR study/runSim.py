from matplotlib import pyplot as plt
from LHeC_shieldingStudy import quadsHalfQuads_full as lf
import pybdsim
import numpy as np
import pandas as pd
import os
import csv

os.chdir("/Users/connormonaghan/Documents/LHeC/IRstudy_full")

pybdsim.Data.LoadROOTLibraries()
plt.style.use(['science','no-latex','ieee', 'grid'])

material    = 'Pb'
runKey      = 'quads_noExtra'

thicknesses = [0.021, 0.022, 0.025, 0.03, 0.04, 0.06, 0.08]
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

print(buffer)
# Quick plot to see trend in data
fig = plt.figure(1)
plt.errorbar(percentAtten[:,0], thicknesses, xerr=percentAtten[:,1])

fig.savefig("runSim-output/{}_{}_absorbed.pdf".format(material, runKey), format='pdf')

df = pd.DataFrame(data=percentAtten, index=thicknesses, columns=["Frac. survival", "Frac. err", "Frac. min", "Frac. max",
                                                                "Tot. w/ cut", "Tot. w/ cut err",
                                                                "Tot. e Aper", "Tot. e Aper err", "Tot. p Aper", "Tot. p Aper err",
                                                                "Tot. w/o zp", "Tot. w/o zp err"])

df.to_csv('DATA/run_{}_{}/{}_runData_{}.csv'.format(material, runKey, material, runKey), float_format='%f')

buff = pd.DataFrame(data=buffer)
buff.to_csv('DATA/run_{}_{}/{}_runBuffer_{}.csv'.format(material,runKey, material,runKey), float_format='%f')

#np.savetxt('DATA/{}/{}_absorbed.csv'.format(runKey, material), a, delimiter=',')

#f = open('DATA/{}/{}_cutInfo.csv'.format(runKey, material), 'w')
#
#w = csv.writer(f)
#
#w.writerow(['Thickness [m]', 'Cut', 'Num before', 'Num After'])
#i = 0
#for t in thicknesses:
#    w.writerow([t, 'none',          cutDataBefore[i, 0], cutDataAfter[i, 0]])
#    w.writerow([t, 'pos>=0.005',   cutDataBefore[i, 1], cutDataAfter[i, 1]])
#    w.writerow([t, 'zp>=0',         cutDataBefore[i, 2], cutDataAfter[i, 2]])
#    w.writerow([t, 'pos+zp',        cutDataBefore[i, 3], cutDataAfter[i, 3]])
#    i = i+1
#
