# A much more detailed plotting script to produce a four panneled plot comparing bdsim data
# to the analytical data.

from pandas import Series as pdS
from matplotlib import pyplot as plt
from itertools import cycle
import matplotlib as mpl
import pandas as pd
import numpy as np

# Function to add BDSIM generated data from genGMAD.py.
# 'file' is a csv file which contains the attenuation data and the specified thickness of material
def plotBDSIMdata(ax, file, style, label):
    d = np.genfromtxt(file, delimiter=',')
    ax.plot((d[1,:])*100, d[0,:]*100, style, label=label, markersize=4)

# Load a set of prebuilt styles for matplotlib.
# Source can be accessed: https://github.com/garrettj403/SciencePlots
plt.style.use(['science','no-latex','ieee','grid'])

# Import the spreadsheet containing the information for various materials and energies
data = pd.read_excel('/Users/connormonaghan/Documents/LHeC/Materials/Atten-Coeffs.xlsx') 

print(data)
# Some setup for plotting the correct BDSIM data
runKey      = 'run1'
beamEnergy  = [0.0002, 0.0003] # GeV
material    = ['bariteConcrete','steel','Cu', 'W', 'Pb', 'U']
print(int(len(material)/2))

# Generate an array of percentages to investigate the lower range < 1%
#intensity_range     = np.linspace(2e-4, 0.01, 20)   # percentage
intensity_range     = np.linspace(0.01, 0.0002, 50)    # percentage 0.00972  0.00083
intensity_incident  = 1                             # percentage, i.e. 100%

# Apply cut to data
data = data.drop(columns=['n','Iron', 'Water'])     # Remove these two for now as they require large amounts of material.
dataCut = data.iloc[3:5, :]                           # [Select only 200 & 300 keV, Remove first column] selecting '3' will select only 200, 3:5 selects 200 & 300

# Empty array for the thickness of material to be stored in.
t = np.zeros((len(dataCut.columns), len(dataCut), len(intensity_range)))  # cm

# Loop over each material and the selected coeffs
for k in range(len(dataCut.columns)):
    # Select mass-attenuation coefficients for the coresponding energy, energy stored in column 'n' (MeV).
    # Density is stored in position 0 of the column.
    coeff   = pdS.to_numpy(dataCut.iloc[:, k])  # cm^2/g
    density = data.iloc[0, k]                   # g/cm^3
    ## CM 21/07: Identified error here in original code. The correct density was not being selected as 
    ##           I was using variable 'dataCut' which did not have the density information stored in it.

    for j in range(len(coeff)):
        for i in range(len(intensity_range)):
            # Calculate the thickness
            t[k,j,i] = np.log(intensity_incident/intensity_range[i])/(coeff[j]*density)

fig, axs = plt.subplots(2, 2)
lines = cycle(["-","-.","--",":"])
fig.set_figheight(4)
fig.set_figwidth(4)

axs[0,0].plot((1-intensity_range)*100, t[0,0,:], "b-", label=dataCut.columns[0])
axs[0,0].plot((1-intensity_range)*100, t[1,0,:], "g--", label=dataCut.columns[1])
axs[0,0].plot((1-intensity_range)*100, t[2,0,:], "r:", label=dataCut.columns[2])
plotBDSIMdata(axs[0,0], "/Users/connormonaghan/Documents/LHeC/Materials/DATA/{}/200keV_bariteConcrete_absorbed.csv".format(runKey), "bD", "BDSIM-BC")
plotBDSIMdata(axs[0,0], "/Users/connormonaghan/Documents/LHeC/Materials/DATA/run2/200keV_steelMagnetite_absorbed.csv", "go", "BDSIM-SM")
plotBDSIMdata(axs[0,0], "/Users/connormonaghan/Documents/LHeC/Materials/DATA/{}/200keV_Cu_absorbed.csv".format(runKey), "rs", "BDSIM-Cu")
#axs[0,0].set_xlabel("% absorbed")
#axs[0,0].set_ylabel("Material thickness [cm]")
axs[0,0].set_xticks(np.arange(99,99.98, 0.49))
axs[0,0].set_title("200keV")
axs[0,0].tick_params(axis="x", direction="inout")

axs[0,1].plot((1-intensity_range)*100, t[3,0,:], "c-", label=dataCut.columns[3])
axs[0,1].plot((1-intensity_range)*100, t[4,0,:], "k--", label=dataCut.columns[4])
axs[0,1].plot((1-intensity_range)*100, t[5,0,:], "m:", label=dataCut.columns[5])
plotBDSIMdata(axs[0,1], "/Users/connormonaghan/Documents/LHeC/Materials/DATA/{}/200keV_W_absorbed.csv".format(runKey), "cD", "BDSIM-W")
plotBDSIMdata(axs[0,1], "/Users/connormonaghan/Documents/LHeC/Materials/DATA/{}/200keV_Pb_absorbed.csv".format(runKey),"ko", "BDSIM-Pb")
plotBDSIMdata(axs[0,1], "/Users/connormonaghan/Documents/LHeC/Materials/DATA/run2/200keV_U_absorbed.csv","ms", "BDSIM-U")
#axs[0,1].set_xlabel("% absorbed")
#axs[0,1].set_ylabel("Material thickness [cm]")
axs[0,1].set_xticks(np.arange(99,99.98, 0.49))
axs[0,1].set_title("200keV")
axs[0,1].tick_params(axis="x", direction="inout")

axs[1,0].plot((1-intensity_range)*100, t[0,1,:], "b-", label=dataCut.columns[0])
axs[1,0].plot((1-intensity_range)*100, t[1,1,:], "g--", label=dataCut.columns[1])
axs[1,0].plot((1-intensity_range)*100, t[2,1,:], "r:", label=dataCut.columns[2])
plotBDSIMdata(axs[1,0], "/Users/connormonaghan/Documents/LHeC/Materials/DATA/{}/300keV_bariteConcrete_absorbed.csv".format(runKey), "bD", "BDSIM-BC")
plotBDSIMdata(axs[1,0], "/Users/connormonaghan/Documents/LHeC/Materials/DATA/run2/300keV_steelMagnetite_absorbed.csv", "go", "BDSIM-SM")
plotBDSIMdata(axs[1,0], "/Users/connormonaghan/Documents/LHeC/Materials/DATA/{}/300keV_Cu_absorbed.csv".format(runKey), "rs", "BDSIM-Cu")
#axs[0,0].set_xlabel("% absorbed")
#axs[0,0].set_ylabel("Material thickness [cm]")
axs[1,0].set_xticks(np.arange(99,99.98, 0.49))
axs[1,0].set_title("300keV")
axs[1,0].tick_params(axis="x", direction="inout")

axs[1,1].plot((1-intensity_range)*100, t[3,1,:], "c-", label=dataCut.columns[3])
axs[1,1].plot((1-intensity_range)*100, t[4,1,:], "k--", label=dataCut.columns[4])
axs[1,1].plot((1-intensity_range)*100, t[5,1,:], "m:", label=dataCut.columns[5])
plotBDSIMdata(axs[1,1], "/Users/connormonaghan/Documents/LHeC/Materials/DATA/{}/300keV_W_absorbed.csv".format(runKey), "cD", "BDSIM-W")
plotBDSIMdata(axs[1,1], "/Users/connormonaghan/Documents/LHeC/Materials/DATA/{}/300keV_Pb_absorbed.csv".format(runKey),"ko", "BDSIM-Pb")
plotBDSIMdata(axs[1,1], "/Users/connormonaghan/Documents/LHeC/Materials/DATA/run2/300keV_U_absorbed.csv","ms", "BDSIM-U")
#axs[0,1].set_xlabel("% absorbed")
#axs[0,1].set_ylabel("Material thickness [cm]")
axs[1,1].set_xticks(np.arange(99,99.98, 0.49))
axs[1,1].set_title("300keV")
axs[1,1].tick_params(axis="x", direction="inout")

handles = []
labels = []

axHandles, axLabel = axs[0,0].get_legend_handles_labels()
handles.extend(axHandles)
labels.extend(axLabel)
axHandles, axLabel = axs[0,1].get_legend_handles_labels()
handles.extend(axHandles)
labels.extend(axLabel)

fig.supylabel("Material thickness [cm]")
fig.supxlabel("% absorbed")

fig.legend(handles, labels, bbox_to_anchor=(1, 0.5), loc='center left',frameon=True, framealpha=1,fancybox=False)
fig.tight_layout()

# for e in range(len(beamEnergy)):
#     for m in range(int(len(material)/2)):
#         axs[e,m].plot((1-intensity_range)*100, np.divide(t[k,e,:], 100), label=dataCut.columns[k], linestyle=next(lines))
#         # for k in range(len(dataCut.columns)):
#         #     axs[e,m].plot((1-intensity_range)*100, np.divide(t[k,e,:], 100), label=dataCut.columns[k], linestyle=next(lines))
#         plotBDSIMdata(axs[e,m], '/Users/connormonaghan/Documents/LHeC/Materials/DATA/{}/{}keV_{}_absorbed.csv'.format(runKey, int(beamEnergy[e]*1e6), material[m]), 
#                     'x', 'BDSIM-{}'.format(material[m]))

fig.savefig("/Users/connormonaghan/Documents/LHeC/Materials/plotThickness-output/fig-updated.pdf", format='pdf')
