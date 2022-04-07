import pandas as pd
import sys
import getopt
import pyopenms as  pyms
import argparse, pathlib

parser = argparse.ArgumentParser()
parser.add_argument('mztab', type=pathlib.Path)
parser.add_argument('idxml', type=pathlib.Path)
args = parser.parse_args()

# path = "/Users/adams/Documents/master/master 2/Stage en thesis/Data/sars_cov_2/"
# bait = "qx017077"
# mztab = path + 'mztab/' + bait + '.mztab'
# output_path = path + 'peptideIndexer/test/' + bait + '.idxml'

output_path = args.idxml

run_name = 'unknown'
skiplines = 0
with open(args.mztab) as f_in:
    line = next(f_in)
    while line.split('\t', 1)[0] != 'PSH':
        if 'ms_run[1]-location' in line:
            run_name = line.split('\t')[2]
        line = next(f_in)
        skiplines += 1

# run_name = 'unknown'
# skiplines = 0
# with open(mztab) as f_in:
#     line = next(f_in)
#     while line.split('\t', 1)[0] != 'PSH':
#         if 'ms_run[1]-location' in line:
#             run_name = line.split('\t')[2]
#         line = next(f_in)
#         skiplines += 1

psms = pd.read_csv(args.mztab, sep='\t', header=skiplines, index_col='PSM_ID')
# psms = pd.read_csv(args.mztab, sep='\t', header=skiplines)

# psms = pd.read_csv(mztab, sep='\t', header=skiplines, index_col='PSM_ID')
peptide_ids = []
for _, psm in psms.iterrows():
    peptide_id = pyms.PeptideIdentification()
    peptide_id.setRT(psm['retention_time'])
    peptide_id.setMZ(psm['exp_mass_to_charge'])
    peptide_id.setScoreType('q-value')
    peptide_id.setHigherScoreBetter(False)
    peptide_id.setIdentifier(run_name)
    peptide_hit = pyms.PeptideHit()
    peptide_hit.setScore(psm['search_engine_score[2]'])
    peptide_hit.setRank(1)
    peptide_hit.setCharge(psm['charge'])
    peptide_hit.setSequence(pyms.AASequence.fromString(psm['sequence']))
    peptide_id.setHits([peptide_hit])
    peptide_ids.append(peptide_id)

protein_id = pyms.ProteinIdentification()
protein_id.setIdentifier(run_name)
pyms.IdXMLFile().store(output_path, [protein_id], peptide_ids)