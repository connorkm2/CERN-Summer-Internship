from matplotlib import pyplot as plt
import pybdsim
import numpy as np
import os

os.chdir("/Users/connormonaghan/Documents/LHeC/materials")

pybdsim.Data.LoadROOTLibraries()

# Load a set of prebuilt styles for matplotlib.
# Source can be accessed: https://github.com/garrettj403/SciencePlots
plt.style.use(['science','no-latex','ieee', 'grid'])

# Setup info for storing and loading the data, make sure the folders exist before continuing (see below)
runKey      = 'run2'
beamEnergy  = 0.0002    # GeV
material    = 'steelMagnetite'       # The target material to be studied

# The range of thicknesses to generate data for, this will need to be adjusted, have a guess before running and then
# fine tune it as required 
thicknesses = [0.12, 0.13, 0.135, 0.15, 0.2]
ngenerate   = 100000        # num photons to generate

# setup for storage of data
percentAtten = np.zeros(len(thicknesses))

i = 0
for t in thicknesses:
    # Begin building the acelerator lattice
    a = pybdsim.Builder.Machine()

    a.AddIncludePre("bariteConcrete.gmad")
    a.AddDrift('d1', 0.1)
    #a.AddDrift('airCylinder', 0.5, apertureType="circularVacuum", vacuumMaterial="air", aper1=0.4)
    a.AddRCol('target', t, material=material, xsize=0, ysize=0)
    a.AddDrift('airCylinder_1', 0.5, apertureType="circularVacuum", aper1=0.4)
    a.AddSampler('all')

    b = pybdsim.Beam.Beam(particletype='gamma',
                        energy=beamEnergy, 
                        distrtype='gauss', 
                        sigmaX='10*um', 
                        sigmaXp=0, 
                        sigmaY='10*um',
                        sigmaYp=0,
                        sigmaE='0.01')
    a.AddBeam(b)

    # Set the physics options
    o = pybdsim.Options.Options()
    o.SetPhysicsList('em')
    o.SetStopSecondaries(True)
    a.AddOptions(o)

    # write the accelerator to the GMAD file
    a.Write('GMAD/target')

    # run the simulation with bdsim
    pybdsim.Run.Bdsim('GMAD/target.gmad', 
                      'DATA/{}/{}keV_{}-{}m'.format(runKey,int(beamEnergy*1e6),material,t), ngenerate=ngenerate)
    # load the data
    pybdsim.Run.Rebdsim("analysisConfig.txt", "DATA/{}/{}keV_{}-{}m.root".format(runKey,int(beamEnergy*1e6),material,t), "DATA/{}/out-{}keV_{}-{}m.root".format(runKey,int(beamEnergy*1e6),material,t))
    d = pybdsim.Data.Load("DATA/{}/out-{}keV_{}-{}m.root".format(runKey,int(beamEnergy*1e6),material,t))

    # airCylinder_1 = pybdsim.Data.SamplerData(d,'airCylinder_1')
    
    # calculate the percentage which has been attenuated
    percentAtten[i] = 1-(d.histogramspy['Event/SimpleHistograms/NPhotonsInX'].entries/ngenerate)
    i = i +1 

# Quick plot to see trend in data
# Use other plotting script to make prettier plots
fig = plt.figure(1)
plt.plot(percentAtten, thicknesses)
plt.xscale('log')
print(percentAtten)
fig.savefig("pyBDSIM-outputs/{}_{}keV_{}_absorbed.pdf".format(runKey,int(beamEnergy*1e6), material), format='pdf')

a = np.asarray([thicknesses,  percentAtten])
np.savetxt('DATA/{}/{}keV_{}_absorbed.csv'.format(runKey, int(beamEnergy*1e6), material), a, delimiter=',')
