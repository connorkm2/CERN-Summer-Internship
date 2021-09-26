from itertools import cycle
from matplotlib import pyplot as plt
from scipy import optimize as opt
import numpy as np
import pandas as pd
import pybdsim
import os

os.chdir("/Users/connormonaghan/Documents/LHeC/Final_Sims")

pybdsim.Data.LoadROOTLibraries()
plt.style.use(['science','no-latex','ieee', 'grid'])

materials = ["Cu", "Pb", "W", "steelMagnetite", "bariteConcrete"]
runKey = "Dipole_half"

d1 = np.genfromtxt("DATA/Dipole_half/Cu_runData.csv", delimiter=',', skip_header=1)
d2 = np.genfromtxt("DATA/Dipole_half/Pb_runData.csv", delimiter=',', skip_header=1)
d3 = np.genfromtxt("DATA/Dipole_half/W_runData.csv",  delimiter=',', skip_header=1)

fig2 = plt.figure()
fig2.set_figheight(3)
fig2.set_figwidth(2.5)

plt.errorbar(d1[:,1]*100, d1[:,0]*100, fmt=".r", xerr=d1[:,2]*100, capsize=1.2, elinewidth=0.7, label="Copper", markersize=5)
plt.errorbar(d2[:,1]*100, d2[:,0]*100, fmt=".k", xerr=d2[:,2]*100, capsize=1.2, elinewidth=0.7, label="Lead", markersize=5)
plt.errorbar(d3[:,1]*100, d3[:,0]*100, fmt=".c", xerr=d3[:,2]*100, capsize=1.2, elinewidth=0.7, label="Tungsten",  markersize=5)

plt.xlabel("% absorbed")
plt.ylabel("Material thickness [cm]")
plt.xlim([99, 100.01])

fig2.legend(loc='upper left',frameon=True, framealpha=1,fancybox=False)

fig2.tight_layout()
fig2.savefig("Dipole_half_1.pdf", format='pdf')

d4 = np.genfromtxt("DATA/Dipole_half/steelMagnetite_runData.csv", delimiter=',', skip_header=1)
d5 = np.genfromtxt("DATA/Dipole_half/bariteConcrete_runData.csv", delimiter=',', skip_header=1)

fig3 = plt.figure()
fig3.set_figheight(3)
fig3.set_figwidth(2.5)

plt.errorbar(d4[:,1]*100, d4[:,0]*100, fmt=".g", xerr=d4[:,2]*100, capsize=1.2, elinewidth=0.7, label="Steel magnetite", markersize=5)
plt.errorbar(d5[:,1]*100, d5[:,0]*100, fmt=".b", xerr=d5[:,2]*100, capsize=1.2, elinewidth=0.7, label="Barite",  markersize=5)

plt.xlabel("% absorbed")
plt.ylabel("Material thickness [cm]")
plt.xlim([99, 100.01])

fig3.legend(loc='upper left',frameon=True, framealpha=1,fancybox=False)

fig3.tight_layout()
fig3.savefig("Dipole_half_2.pdf", format='pdf')