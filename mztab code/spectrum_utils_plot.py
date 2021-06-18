# /Users/adams/anaconda3/envs/ann_solo/bin/python3
import matplotlib.pyplot as plt
import spectrum_utils.plot as sup
import spectrum_utils.spectrum as sus
from pyteomics import mgf

#	Mirror plot
spectra = []
for spectrum_dict in mgf.read('/Users/adams/Documents/PhD/SARS-CoV-2/Data/Workspace/mgf/qx017105.mgf'):
	if 'scan=15917' in spectrum_dict['params']['title']:
		identifier = spectrum_dict['params']['title']
		precursor_mz_spec = spectrum_dict['params']['pepmass'][0]
		precursor_charge = spectrum_dict['params']['charge'][0]
		mz = spectrum_dict['m/z array']
		intensity = spectrum_dict['intensity array']
		retention_time = float(spectrum_dict['params']['rtinseconds'])
		peptide = 'GQGVPINTNSSPDDQIGYYR'
		modifications = {10: 79.966331}			# -28.031300

# Create the MS/MS spectrum.
spectra.append(sus.MsmsSpectrum(identifier, precursor_mz_spec,
								precursor_charge, mz, intensity,
								retention_time=retention_time,
								peptide=peptide,
								modifications=modifications)
				.annotate_peptide_fragments(0.02, 'Da', ion_types='aby'))

for spectrum_dict in mgf.read('/Users/adams/Documents/PhD/SARS-CoV-2/Data/Workspace/mgf/qx017105.mgf'):
	if 'scan=15100' in spectrum_dict['params']['title']:
		identifier = spectrum_dict['params']['title']
		precursor_mz = spectrum_dict['params']['pepmass'][0]
		precursor_charge = spectrum_dict['params']['charge'][0]
		mz = spectrum_dict['m/z array']
		intensity = spectrum_dict['intensity array']
		retention_time = float(spectrum_dict['params']['rtinseconds'])
		peptide = 'GQGVPINTNSSPDDQIGYYR'

spectra.append(sus.MsmsSpectrum(identifier, precursor_mz,
								precursor_charge, mz, intensity,
								retention_time=retention_time,
								peptide=peptide)
				.annotate_peptide_fragments(0.02, 'Da', ion_types='aby'))

fig, ax = plt.subplots(figsize=(12, 6))
spectrum_top, spectrum_bottom = spectra
sup.mirror(spectrum_top, spectrum_bottom, ax=ax)

ax.text(0.5, 1.10, f'Top: mzspec:PXD018117:qx017105:scan:15917:GQGVPINTNSSPDDQIGYYR',
        horizontalalignment='center', verticalalignment='bottom',
        fontsize='x-large', fontweight='bold', transform=ax.transAxes)
ax.text(0.5, 1.06, f'Bottom: mzspec:PXD018117:qx017105:15100:scan:GQGVPINTNSSPDDQIGYYR',
        horizontalalignment='center', verticalalignment='bottom',
        fontsize='x-large', fontweight='bold', transform=ax.transAxes)
ax.text(0.5, 1.02, f'Precursor m/z: {precursor_mz_spec:.4f}, '
                   f'Charge: {precursor_charge:.4f}',
        horizontalalignment='center', verticalalignment='bottom',
        fontsize='large', transform=ax.transAxes)
plt.savefig(f'/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Figures/Mirror/Phospho/mzspec:PXD018117:qx017105:scan:15917:GQGVPINTNSSPDDQIGYYR-2.png', 
	dpi=300, bbox_inches='tight')
plt.close()