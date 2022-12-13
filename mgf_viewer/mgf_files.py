from pyteomics import mgf, pepxml, mass, pylab_aux
import os
import pylab
import numpy as np
from operator import add
from pathlib import Path
import matplotlib.pyplot as plt

BASE_LOCATION = Path(__file__).resolve().parent.parent
DATA_LOCATION = BASE_LOCATION.joinpath('data')
FILE_LOCATION = BASE_LOCATION.joinpath('mgf_viewer')

with mgf.read(str(DATA_LOCATION.joinpath('high.human.PXD004424/peaks.db.mgf'))) as spectra:
    spectrum = spectra[50]

print(f'Sequence: {spectrum["params"]["seq"]}\n'
      f'Charge: {spectrum["params"]["charge"]}\n'
      f'Mass: {spectrum["params"]["pepmass"][0]}\n'
      f'Single Charge Mass: {mass.calculate_mass(spectrum["params"]["seq"])}')

print(spectrum['m/z array'])
# print(spectrum['intensity array'])
array_1 = []
for i in range(len(spectrum['m/z array'])):
    array_1.append((spectrum['m/z array'][i]/mass.calculate_mass(spectrum["params"]["seq"])))

print(array_1)
vector_1 = np.array((array_1))
vector_2 = np.flip(np.array(array_1))

print(mass.calculate_mass(spectrum["params"]["seq"][0:5], ))

vector3 = vector_1[:, None]+vector_2
print(vector3.shape)
plt.matshow(np.log(vector3))
plt.show()
# print((mass.calculate_mass(spectrum["params"]["seq"]) - spectrum['m/z array'][-2]))
# print(spectrum["params"]["seq"][-1])

# pylab.figure()
# # pylab_aux.plot_spectrum(spectrum,title=spectrum['params']['title'])
# pylab_aux.annotate_spectrum(spectrum, peptide=spectrum['params']['seq'],title=spectrum['params']['seq'], maxcharge=spectrum['params']['charge'][0])
# pylab.show()
