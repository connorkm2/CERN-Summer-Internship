import pybdsim
from matplotlib import pyplot as plt

plt.style.use(['science','no-latex','ieee', 'grid'])

d = pybdsim.Data.Load("/Users/connormonaghan/Documents/LHeC/Final_Sims/EnergyDistr/combi.root")

fig = pybdsim.Plot.Histogram1D(d.histogramspy['Event/SimpleHistograms/E_photons_dDRIFT50'], figsize=(4,3))
ax = plt.gca()
ax.set_xlabel("Photon energy [GeV]")
ax.set_ylabel("Number of photons")
ax.set_yscale("log")

fig.tight_layout()
fig.savefig("energy_END.pdf")

fig2 = pybdsim.Plot.Histogram1D(d.histogramspy['Event/SimpleHistograms/E_photons_uDRIFT30_1'], figsize=(4,3))
ax2 = plt.gca()
ax2.set_xlabel("Photon energy [GeV]")
ax2.set_ylabel("Number of photons")
ax2.set_yscale("log")

fig2.tight_layout()
fig2.savefig("energy_BEND.pdf")

fig3 = pybdsim.Plot.Histogram1D(d.histogramspy['Event/SimpleHistograms/E_photons_COL_END_0'], figsize=(4,3))
ax3 = plt.gca()
ax3.set_xlabel("Photon energy [GeV]")
ax3.set_ylabel("Number of photons")
ax3.set_yscale("log")

fig3.tight_layout()
fig3.savefig("energy_END_af.pdf")