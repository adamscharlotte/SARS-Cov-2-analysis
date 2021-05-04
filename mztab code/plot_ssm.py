# /Users/adams/anaconda3/envs/ann_solo/bin/python3
# "/Users/adams/Documents/PhD/SARS-CoV-2/Data/Workspace/mztab/qx017160"
#

# import argparse
import os
import urllib.parse as urlparse

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from spectrum_utils import plot
from spectrum_utils.spectrum import PeptideFragmentAnnotation

from ann_solo import reader
from ann_solo import spectrum_match
from ann_solo.config import config
from ann_solo.spectrum import process_spectrum

def set_matching_peaks(library_spectrum, query_spectrum):
    peak_matches = spectrum_match.get_best_match(
        query_spectrum, [library_spectrum],
        config.fragment_mz_tolerance, config.allow_peak_shifts)[2]
    query_spectrum.annotation = np.full_like(query_spectrum.mz, None, object)
    for peak_match in peak_matches:
        library_annotation = library_spectrum.annotation[peak_match[1]]
        if library_annotation is not None and library_annotation.ion_type in 'by':
            query_spectrum.annotation[peak_match[0]] = library_annotation
        else:
            fragment_annotation = PeptideFragmentAnnotation(1, 1, 'z', 0)
            fragment_annotation.ion_type = 'unknown'
            query_spectrum.annotation[peak_match[0]] =\
                library_spectrum.annotation[peak_match[1]] =\
                fragment_annotation

# Read the mzTab file.
metadata = {}
with open("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Workspace/mztab/qx017160.mztab") as f_mztab:
    for line in f_mztab:
        line_split = line.strip().split('\t')
        if line_split[0] == 'MTD':
            metadata[line_split[1]] = line_split[2]
        else:
            break   # Metadata lines should be on top.

ssms = reader.read_mztab_ssms("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Workspace/mztab/qx017160.mztab")
# make sure the SSM ids are strings.
ssms.index = ssms.index.map(str)

# Recreate the search configuration.
settings = []
# Search settings.
for key in metadata:
    if 'software[1]-setting' in key:
        param = metadata[key][: metadata[key].find(' ')]
        value = metadata[key][metadata[key].rfind(' ') + 1:]
        if value != 'None':
            if value != 'False':
                settings.append('--{}'.format(param))
            if value not in ('False', 'True'):
                settings.append(value)

# File names.
settings.append('dummy_spectral_library_filename')
settings.append('dummy_query_filename')
settings.append('dummy_output_filename')
config.parse(' '.join(settings))

# Retrieve information on the requested query.
query_id = "qx017160.9348.9348.2"
query_usi = "mzspec:PXD018117:qx017160:scan:9348:[-28.031300]?GYGCSCDQLR"     #[UNIMOD:21]    [UNIMOD:21#g1]
txt = 'm/z'
query_spectrum_number = "controllerType=0 controllerNumber=1 scan=9348"
query_uri = urlparse.urlparse(urlparse.unquote(
    metadata['ms_run[1]-location']))
query_filename = os.path.abspath(os.path.join(
    query_uri.netloc, query_uri.path))
ssm = ssms.loc[query_id]
library_id = ssm['accession']
library_uri = urlparse.urlparse(urlparse.unquote(ssm['database']))
library_filename = os.path.abspath(os.path.join(
    library_uri.netloc, library_uri.path))
score = ssm['search_engine_score[1]']

# Read library and query spectrum.
with reader.SpectralLibraryReader(library_filename) as lib_reader:
	library_spectrum = lib_reader.get_spectrum(library_id, False)

query_spectrum = None
for spec in reader.read_mgf(query_filename):
	if spec.identifier == query_spectrum_number:
		query_spectrum = spec
		break

# verify that the query spectrum was found
if query_spectrum is None:
	raise ValueError('Could not find the specified query spectrum')

# Set the matching peaks in the query spectrum to correctly color them.
set_matching_peaks(library_spectrum, query_spectrum)
# Modify the colors to differentiate non-matching peaks.
plot.colors[None] = '#757575'

# Plot the match.
fig, ax = plt.subplots(figsize=(12, 6))
# Plot with annotations.
plot.mirror(query_spectrum, library_spectrum,
            {'color_ions': True, 'annotate_ions': True}, ax)

# ax.set_ylim(-1.1, 1.05)

ax.text(0.5, 1.06, f'{query_usi}',
        horizontalalignment='center', verticalalignment='bottom',
        fontsize='x-large', fontweight='bold', transform=ax.transAxes)
ax.text(0.5, 1.02,  f'Precursor ${txt}$ (top): {query_spectrum.precursor_mz:.4f}, '
                    f'Library ${txt}$ (bottom): {library_spectrum.precursor_mz:.4f}, '
                    f'Charge: {query_spectrum.precursor_charge}',
        horizontalalignment='center', verticalalignment='bottom',
        fontsize='large', transform=ax.transAxes)

plt.savefig("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Figures/Mirror/SNO/"f'{query_id}.png', dpi=300, bbox_inches='tight')
plt.close()
