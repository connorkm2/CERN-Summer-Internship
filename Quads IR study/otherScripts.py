import numpy as np
import pybdsim
import csv
import os
import scipy.optimize as opt
from matplotlib import pyplot as plt

os.chdir("/Users/connormonaghan/Documents/LHeC/IRstudy")
pybdsim.Data.LoadROOTLibraries()
plt.style.use(['science','no-latex','ieee', 'grid'])

def genBDSIMcsv():
    material = "Pb"
    runKey   = "run2b"
    
    thicknesses = [0.021, 0.022, 0.025, 0.027, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08]
    percentAtten = np.zeros(len(thicknesses))
    
    i = 0
    for t in thicknesses:
        d = pybdsim.Data.Load("DATA/{}/{}-{}m.root".format(runKey,material,t))
        target  = pybdsim.Data.SamplerData(d,'target')
        drift   = pybdsim.Data.SamplerData(d,'DRIFT_0')
        numBefore = len(drift.data['energy'][(drift.data['partID']==22)&(drift.data['zp']>=0)&(drift.data['x']>=0.005)])
        numAfter  = len(target.data['energy'][(target.data['partID']==22)&(target.data['zp']>=0)&(target.data['x']>=0.005)])
    
        percentAtten[i] = 1-(numAfter/numBefore)
        i = i+1
    
    a = np.asarray([thicknesses,  percentAtten])
    np.savetxt('DATA/{}/{}_absorbed.csv'.format(runKey, material), a, delimiter=',')

def plotX():
    d = pybdsim.Data.Load("other-outputs/output-noEM.root")
    target = pybdsim.Data.SamplerData(d,'target')
    drift = pybdsim.Data.SamplerData(d,'DRIFT_0')

    X = np.asarray([target.data['x'][target.data['partID']==11], target.data['y'][target.data['partID']==11]])

    fig = plt.figure(1)
    plt.hist(drift.data['x'][drift.data['partID']==11], 50, (-0.0047, 0.0047))
    plt.xlabel("x [m]")
    plt.ylabel("Freq.")
    plt.title("Electrons")

    fig2 = plt.figure(2)
    plt.hist(drift.data['x'][drift.data['partID']==22], 50, (min(drift.data['x'][drift.data['partID']==22]), 0.06), label="Photons", histtype='step', edgecolor='green', ls="-")
    plt.xlim((-0.13,0.01))
    plt.xlabel("x [m]")
    plt.ylabel("Number of particles")
    plt.title("Photons")
    fig2.savefig("/Users/connormonaghan/Documents/LHeC/IRstudy/runSim-output/xPhotons.pdf")

    fig.savefig("/Users/connormonaghan/Documents/LHeC/IRstudy/runSim-output/xElectrons.pdf")

    fig = plt.figure(3)
    plt.hist(drift.data['x'][drift.data['partID']==11], 50, (-0.0047, 0.0047), label="Electrons", histtype='step', edgecolor='red', fill=False, ls="-", facecolor='none')
    plt.hist(drift.data['x'][drift.data['partID']==22], 50, (-0.0047, 0.0047), label="Photons", histtype='step', edgecolor='green', ls="-")
    plt.xlabel("x [m]")
    plt.ylabel("Number of particles")

    fig.savefig("/Users/connormonaghan/Documents/LHeC/IRstudy/runSim-output/ep.pdf")

    print(len(drift.data['x'][drift.data['partID']==11]))
    print(len(drift.data['x'][drift.data['partID']==22]))

    fig = plt.figure(4)
    plt.hist2d(drift.data['x'][drift.data['partID']==11], drift.data['y'][drift.data['partID']==11], 50)
    fig.savefig("/Users/connormonaghan/Documents/LHeC/IRstudy/runSim-output/2dElectrons.pdf")

    fig = plt.figure(5)
    plt.hist2d(drift.data['x'][drift.data['partID']==22], drift.data['y'][drift.data['partID']==22], 50)
    fig.savefig("/Users/connormonaghan/Documents/LHeC/IRstudy/runSim-output/2dPhotons.pdf")

# Write here the name of the scirpt you want to run