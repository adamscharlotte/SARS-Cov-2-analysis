import sys
#print(sys.executable)       #sys.executable contains full path of the currently running Python interpreter
import pandas as pd
import sys
import argparse, pathlib

parser = argparse.ArgumentParser()
parser.add_argument('mztab', type=pathlib.Path)
parser.add_argument('bait', type=str)
args = parser.parse_args()

path = '/Users/adams/Documents/PhD/SARS-CoV-2/Data/Workspace/'
output_path = path + 'csv/mztab/' + args.bait + '.csv'
# output_path = args.bait + '.csv'

with open (args.mztab,'r') as f:
    outF = open(output_path, 'w')
    for line in f:
        if line.startswith('MTD'):
            #mztabHead.append(line)
            continue
        else:
            outF.write(line)
    outF.close()
            #mztabList.append(line)

#data = pd.read_csv(output_path, sep = '\t')