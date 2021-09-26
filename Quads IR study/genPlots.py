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

fig1, axs1 = plt.subplots(len(materials), 1)
#lines = cycle(["-","-.","--",":"])
fig1.set_figheight(10)
fig1.set_figwidth(2)


for i in range(len(materials)):
    d = np.genfromtxt("DATA/{}/{}_runData.csv".format(runKey, materials[i]), delimiter=',')
    
    axs1[i].plot(d[1,:]*100, d[0,:]*100, ".k", label=materials[i], markersize=4)

    axs1[i].set_xticks(np.arange(99,99.98, 0.49))
    axs1[i].set_xlim([99, 99.98])
    axs1[i].set_xlabel("% absorbed")
    axs1[i].set_ylabel("Material thickness [cm]")
    axs1[i].set_title(materials[i])

fig1.savefig("figs.pdf", format='pdf')

fig, axs = plt.subplots(3, 3)
#lines = cycle(["-","-.","--",":"])
fig.set_figheight(8)
fig.set_figwidth(8)

d = np.genfromtxt("DATA/Dipole_half/Cu_runData.csv", delimiter=',')
axs[0,0].plot(d[1,:]*100, d[0,:]*100, ".k", label="Cu", markersize=4)
#axs[0,0].set_xticks(np.arange(99,99.98, 0.49))
axs[0,0].set_xlim([99, 99.98])
axs[0,0].set_ylim([5,11])
axs[0,0].set_ylabel("Material thickness [cm]")
axs[0,0].set_xlabel("% absorbed")
axs[0,0].set_title("Cu")

d = np.genfromtxt("DATA/Dipole_half/Pb_runData.csv", delimiter=',')
axs[0,1].plot(d[1,:]*100, d[0,:]*100, ".k", label="Pb", markersize=4)
#axs[0,1].set_xticks(np.arange(99,99.98, 0.49))
axs[0,1].set_xlim([99, 100.01])
axs[0,1].set_ylim([0,9])
axs[0,1].set_ylabel("Material thickness [cm]")
axs[0,1].set_xlabel("% absorbed")
axs[0,1].set_title("Pb")

d = np.genfromtxt("DATA/Dipole_half/W_runData.csv", delimiter=',')
axs[0,2].plot(d[1,:]*100, d[0,:]*100, ".k", label="W", markersize=4)
#axs[0,2].set_xticks(np.arange(99,99.98, 0.49))
axs[0,2].set_xlim([99, 100.01])
#axs[0,2].set_ylim([5,11])
axs[0,2].set_ylabel("Material thickness [cm]")
axs[0,2].set_xlabel("% absorbed")
axs[0,2].set_title("W")

d = np.genfromtxt("DATA/Dipole_half/steelMagnetite_runData.csv", delimiter=',')
axs[1,0].plot(d[1,:]*100, d[0,:]*100, ".k", label="steelMagnetite", markersize=4)
#axs[1,0].set_xticks(np.arange(99,99.98, 0.49))
axs[1,0].set_xlim([99, 99.98])
#axs[1,0].set_ylim([5,11])
axs[1,0].set_ylabel("Material thickness [cm]")
axs[1,0].set_xlabel("% absorbed")
axs[1,0].set_title("Steel Magnetite")

d = np.genfromtxt("DATA/Dipole_half/bariteConcrete_runData.csv", delimiter=',')
axs[1,1].plot(d[1,:]*100, d[0,:]*100, ".k", label="bariteConcrete", markersize=4)
#axs[1,1].set_xticks(np.arange(99,99.98, 0.49))
axs[1,1].set_xlim([99, 99.98])
#axs[1,1].set_ylim([5,11])
axs[1,1].set_ylabel("Material thickness [cm]")
axs[1,1].set_xlabel("% absorbed")
axs[1,1].set_title("Barite Concrete")

fig.tight_layout()
fig.savefig("figs1.pdf", format='pdf')

d1 = np.genfromtxt("DATA/Dipole_half/Cu_runData.csv", delimiter=',')
d2 = np.genfromtxt("DATA/Dipole_half/Pb_runData.csv", delimiter=',')
d3 = np.genfromtxt("DATA/Dipole_half/W_runData.csv",  delimiter=',')

fig2 = plt.figure()
fig2.set_figheight(3)
fig2.set_figwidth(2.5)

plt.plot(d1[1,:]*100, d1[0,:]*100, ".k", label="Copper", markersize=5)
plt.plot(d2[1,:]*100, d2[0,:]*100, ".r", label="Lead", markersize=5)
plt.plot(d3[1,:]*100, d3[0,:]*100, ".b", label="Tungsten",  markersize=5)

plt.xlabel("% absorbed")
plt.ylabel("Material thickness [cm]")
plt.xlim([99, 100.01])

fig2.legend(bbox_to_anchor=(0.9, 0.5), loc='center left',frameon=True, framealpha=1,fancybox=False)

fig2.tight_layout()
fig2.savefig("figs2.pdf", format='pdf')

d4 = np.genfromtxt("DATA/Dipole_half/steelMagnetite_runData.csv", delimiter=',')
d5 = np.genfromtxt("DATA/Dipole_half/bariteConcrete_runData.csv", delimiter=',')

fig3 = plt.figure()
fig3.set_figheight(3)
fig3.set_figwidth(2.5)

plt.plot(d4[1,:]*100, d4[0,:]*100, ".k", label="Steel Magnetite", markersize=5)
plt.plot(d5[1,:]*100, d5[0,:]*100, ".r", label="Barite", markersize=5)

plt.xlabel("% absorbed")
plt.ylabel("Material thickness [cm]")
plt.xlim([99, 100.01])

fig3.legend(bbox_to_anchor=(0.9, 0.5), loc='center left',frameon=True, framealpha=1,fancybox=False)

fig3.tight_layout()
fig3.savefig("figs3.pdf", format='pdf')