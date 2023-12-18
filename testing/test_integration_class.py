import sys
import pathlib as p
sys.path.append('.')
print(p.WindowsPath(sys.path[-1]).resolve())

from zscan1 import (Integration,Fitting,SILICA_BETA,N_COMPONENTS,INTEGRATION_STEPS)

import numpy as np

import matplotlib
matplotlib.rcParams.update(matplotlib.rcParamsDefault)
from matplotlib import pyplot as plt

z_range=40E-3 # meter
nop=201
lda=800e-9 # meter
silica_n2 = 2.8203E-20 - 3E-27/(lda) + 2E-33/(lda)**2
d0=300E-3 # meter
ra=0.5E-3 # meter
DPhi0 = 0.5
silicaCA_beamwaist = 43E-6 # meter

data_file = np.genfromtxt('./data/2021_07_22__10_10__silica_0-0_1600-0_2.txt',skip_header=23)


positions = np.array([z_range*zz/nop-z_range/2 for zz in range(nop)]) # express in m

tos = []
tcs = []
for d0 in [0,240E-3,300E-3]:
    silica_curves = Integration(
        SILICA_BETA,silica_n2,DPhi0,positions,d0,ra,lda,silicaCA_beamwaist,N_COMPONENTS,INTEGRATION_STEPS)

    calculation = Fitting(silica_curves,DPhi0,silicaCA_beamwaist,0,0,nop,data_file)
    result = calculation.automatic(z_range,"silica","CA",data_file)
    tos.append(silica_curves.open_sum)
    tcs.append(silica_curves.closed_sum)

for i,_ in enumerate([0,0.2,0.566]):
    plt.plot(positions*1000,tcs[i]/tos[i])
    plt.xlabel("z (mm)")
plt.show()