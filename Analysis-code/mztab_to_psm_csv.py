import sys
import pandas as pd
import sys
import argparse, pathlib

parser = argparse.ArgumentParser()
parser.add_argument('mztab', type=pathlib.Path)
parser.add_argument('csv', type=pathlib.Path)
args = parser.parse_args()

output_path = args.csv

with open (args.mztab,'r') as f:
	outF = open(output_path, 'w')
	for line in f:
		if line.startswith('MTD'):
			continue
		if line.startswith('PRT'):
			continue
		if line.startswith('PRH'):
			continue
		else:
			outF.write(line)
	outF.close()