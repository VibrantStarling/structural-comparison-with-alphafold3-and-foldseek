# -*- coding: utf-8 -*-
"""
Created on Fri Mar  7 12:05:43 2025

@author: telma

USES CODE FROM https://github.com/gerdos/PyRAMA/tree/master
"""




import Bio.PDB
from Bio.PDB.MMCIFParser import MMCIFParser
import math
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib.colors as mplcolors


cif_file = 'model.cif'



def ramachandran(cif_file):
    """Accepts a cif file name (string) and returns two lists of phi 
    and psi angles, respectively.
    """
    parser = Bio.PDB.MMCIFParser()
    structure = parser.get_structure('protein', cif_file)
    polypeptide = Bio.PDB.PPBuilder().build_peptides(structure[0])
    
    phi = []
    psi = []

    for strand in polypeptide:
        phipsi = strand.get_phi_psi_list()
        for point in phipsi:
            try:
                phi_point = point[0] * (180 / 3.14159)
                psi_point = point[1] * (180 / 3.14159)
                phi.append(phi_point)
                psi.append(psi_point)
            except TypeError:
                pass
    
    return phi, psi

def _cache_RAMA_PREF_VALUES():
    RAMA_PREF_VALUES = {}
    for key, val in RAMA_PREFERENCES.items():
        RAMA_PREF_VALUES[key] = np.full((360, 360), 0, dtype=np.float64)
        with open(val["file"]) as fn:
            for line in fn:
                if line.startswith("#"):
                    continue
                else:
                    x = int(float(line.split()[1]))
                    y = int(float(line.split()[0]))
                    RAMA_PREF_VALUES[key][x + 180][y + 180] \
                        = RAMA_PREF_VALUES[key][x + 179][y + 179] \
                        = RAMA_PREF_VALUES[key][x + 179][y + 180] \
                        = RAMA_PREF_VALUES[key][x + 180][y + 179] \
                        = float(line.split()[2]) 
    return RAMA_PREF_VALUES

phi, psi = ramachandran(cif_file)
cmap = mplcolors.ListedColormap(['#FFFFFF', '#B3E8FF', '#7FD9FF'])

RAMA_PREFERENCES = {
    "General": {
        "file": "general-ramachandran-data.tsv",
        "cmap": mplcolors.ListedColormap(['#FFFFFF', '#B3E8FF', '#7FD9FF']),
        "bounds": [0, 0.0005, 0.02, 1],
    }}

RAMA_PREF_VALUES = _cache_RAMA_PREF_VALUES()


for idx, (key, val) in enumerate(sorted(RAMA_PREFERENCES.items(), key=lambda x: x[0].lower())):
        fig,ax = plt.subplots(figsize=(8, 8))
        plt.title(key)
        plt.imshow(RAMA_PREF_VALUES[key], cmap=RAMA_PREFERENCES[key]["cmap"],
                   norm=mplcolors.BoundaryNorm(RAMA_PREFERENCES[key]["bounds"], RAMA_PREFERENCES[key]["cmap"].N),
                   extent=(-180, 180, 180, -180))
        plt.scatter(phi, psi, s=2, c='black')
        plt.xlabel('phi $\phi$, degrees')
        plt.ylabel('psi $\psi$, degrees')
        plt.grid(True, linestyle="--", alpha=0.5)
        plt.xlim([-180, 180])
        plt.ylim([-180, 180])
        plt.plot([-180, 180], [0, 0], color="black")
        plt.plot([0, 0], [-180, 180], color="black")
        plt.locator_params(axis='x', nbins=7)
        plt.grid()

plt.tight_layout()
# plt.savefig("asd.png", dpi=300)
plt.show()

