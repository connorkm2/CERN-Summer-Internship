from matplotlib import pyplot as plt
import pybdsim
import numpy as np

plt.style.use(['science','no-latex','ieee', 'grid'])

d = pybdsim.Data.Load("/Users/connormonaghan/Documents/LHeC/IRstudy_full/DATA/run_Pb_quads/Pb-0.025m_quads.root")

DRIFT_1 = pybdsim.Data.SamplerData(d, 'dDRIFT50')
COL_0 = pybdsim.Data.SamplerData(d, 'COL_END_0')

fig1, axs1 = plt.subplots(2,1)
fig1.set_figheight(4)
fig1.set_figwidth(6)

axs1[0].hist(DRIFT_1.data['x'][DRIFT_1.data['partID']==22], 80, (-0.005, 0.4), histtype='step', edgecolor='green', ls="-")
axs1[0].axvline(x=0.00, color='b', linestyle='--', linewidth=0.5)
axs1[0].axvline(x=-0.005, color='b', linestyle='-', linewidth=0.5)
axs1[0].axvline(x=0.005, color='b', linestyle='-', linewidth=0.5)

axs1[0].axvline(x=0.106105, color='r', linestyle='--', linewidth=0.5)
axs1[0].axvline(x=0.106105-0.02, color='r', linestyle='-', linewidth=0.5)
axs1[0].axvline(x=0.106105+0.02, color='r', linestyle='-', linewidth=0.5)
axs1[0].set_yscale("log")

axs1[1].hist(COL_0.data['x'][COL_0.data['partID']==22], 80, (-0.005, 0.4), histtype='step', edgecolor='green', ls="-")
axs1[1].axvline(x=0.00, color='b', linestyle='--', linewidth=0.5)
axs1[1].axvline(x=-0.005, color='b', linestyle='-', linewidth=0.5)
axs1[1].axvline(x=0.005, color='b', linestyle='-', linewidth=0.5)

axs1[1].axvline(x=0.106105, color='r', linestyle='--', linewidth=0.5)
axs1[1].axvline(x=0.106105-0.02, color='r', linestyle='-', linewidth=0.5)
axs1[1].axvline(x=0.106105+0.02, color='r', linestyle='-', linewidth=0.5)
axs1[1].set_yscale("log")

axs1[0].set_ylabel('Number photons')
axs1[1].set_ylabel('Number photons')
axs1[1].set_xlabel('x [m]')

plt.subplots_adjust(hspace=.0)

fig1.savefig('fig1-update.pdf')
