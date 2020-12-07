# /Users/adams/anaconda3/envs/ann_solo/bin/python3
import matplotlib.pyplot as plt
import spectrum_utils.plot as sup
import spectrum_utils.spectrum as sus
from pyteomics import mgf

# Read the spectrum from an MGF file using Pyteomics.
spectrum_dict = mgf.get_spectrum(                   #Read one spectrum
	'/Users/adams/Documents/PhD/SARS-CoV-2/Data/Workspace/mgf/qx017178.mgf',
	'controllerType=0 controllerNumber=1 scan=18575')
identifier = spectrum_dict['params']['title']
precursor_mz = spectrum_dict['params']['pepmass'][0]
precursor_charge = spectrum_dict['params']['charge'][0]
mz = spectrum_dict['m/z array']
intensity = spectrum_dict['intensity array']
retention_time = float(spectrum_dict['params']['rtinseconds'])
peptide = 'ALSSGGSITSPPLSPALPK'

# Create the MS/MS spectrum.
spectrum = sus.MsmsSpectrum(
	identifier, precursor_mz, precursor_charge, mz, intensity,
	retention_time=retention_time, peptide=peptide)#, modifications=modifications)

# Process the MS/MS spectrum.
fragment_tol_mass = 10
fragment_tol_mode = 'ppm'
spectrum = (spectrum.set_mz_range(min_mz=100, max_mz=1400)
			.remove_precursor_peak(fragment_tol_mass, fragment_tol_mode)
			.filter_intensity(min_intensity=0.05, max_num_peaks=50)
			.annotate_peptide_fragments(fragment_tol_mass, fragment_tol_mode,
										ion_types='aby'))

#   Plot the MS/MS spectrum.
fig, ax = plt.subplots(figsize=(12, 6))
sup.spectrum(spectrum, ax=ax)
plt.show()
plt.close()

#   Variable modification
peptide = 'ALSSGGSITSPPLSPALPK'
modifications = {9: 79.966331}  #7 is the index of the modified AA, 0-baseds

spectrum = sus.MsmsSpectrum(
	identifier, precursor_mz, precursor_charge, mz, intensity,
	retention_time=retention_time, peptide=peptide, modifications=modifications)

# Process the MS/MS spectrum.
fragment_tol_mass = 10
fragment_tol_mode = 'ppm'
spectrum = (spectrum.set_mz_range(min_mz=100, max_mz=1400)
			.remove_precursor_peak(fragment_tol_mass, fragment_tol_mode)
			.filter_intensity(min_intensity=0.05, max_num_peaks=50)
			.annotate_peptide_fragments(fragment_tol_mass, fragment_tol_mode,
										ion_types='aby'))

#   Plot the MS/MS spectrum
fig, ax = plt.subplots(figsize=(12, 6))
sup.spectrum(spectrum, ax=ax)
plt.show()
plt.close()



#########################################################################################################################################################

spectra = []
for spectrum_dict in mgf.read('/Users/adams/Documents/PhD/SARS-CoV-2/Data/Workspace/mgf/qx017106.mgf'):
	if 'scan=15952' in spectrum_dict['params']['title']:
		identifier = spectrum_dict['params']['title']
		precursor_mz = spectrum_dict['params']['pepmass'][0]
		precursor_charge = spectrum_dict['params']['charge'][0]
		mz = spectrum_dict['m/z array']
		intensity = spectrum_dict['intensity array']
		retention_time = float(spectrum_dict['params']['rtinseconds'])
		peptide = 'GQGVPINTNSSPDDQIGYYR'
		modifications = {10: 79.966331}

# Create the MS/MS spectrum.
spectra.append(sus.MsmsSpectrum(identifier, precursor_mz,
								precursor_charge, mz, intensity,
								retention_time=retention_time,
								peptide=peptide,
								modifications=modifications)
				.set_mz_range(min_mz=100, max_mz=1400)
				.filter_intensity(0.01, 50)
				.scale_intensity('root')
				.annotate_peptide_fragments(0.5, 'Da', ion_types='aby'))

spectra.append(sus.MsmsSpectrum(identifier, precursor_mz,
								precursor_charge, mz, intensity,
								retention_time=retention_time,
								peptide=peptide)
				.set_mz_range(min_mz=100, max_mz=1400)
				.filter_intensity(0.01, 50)
				.scale_intensity('root')
				.annotate_peptide_fragments(0.5, 'Da', ion_types='aby'))

fig, ax = plt.subplots(figsize=(12, 6))
spectrum_top, spectrum_bottom = spectra
sup.mirror(spectrum_top, spectrum_bottom, ax=ax)
plt.show()

