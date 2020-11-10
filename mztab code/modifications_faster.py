import operator
import os
import pyteomics.auxiliary
import pandas as pd
import numpy as np
import sys
import argparse, pathlib

parser = argparse.ArgumentParser()
parser.add_argument('bait', type=str)
parser.add_argument('unimod_path', type=str)
parser.add_argument('tol_py_csv', type=pathlib.Path)
args = parser.parse_args()

output_path = args.bait + '.csv'

#   load unimod file:
unimod = pd.read_csv(args.unimod_path)

#   load psm file
psms_modified = pd.read_csv(args.tol_py_csv)

mod_names, mod_masses = [], []

for mass_tol_neg, mass_tol_pos in zip(psms_modified['mass_tol_neg'], psms_modified['mass_tol_pos']):
    idx_start, idx_stop = unimod['mod_mass'].searchsorted([mass_tol_neg, mass_tol_pos])
    if idx_stop > idx_start:
        potential_mods = unimod[['mod_mass', 'mod_name']].iloc[np.arange(idx_start, idx_stop)]
        mod_names.append(' / '.join(potential_mods['mod_name'].astype(str)))
        mod_masses.append(' / '.join(potential_mods['mod_mass'].astype(str)))
    else:
        mod_names.append('No direct match found in Unimod')
        mod_masses.append('')

psms_modified['mod'] = mod_names
psms_modified['mod_mass'] = mod_masses
psms_modified.to_csv(output_path, sep = ",")