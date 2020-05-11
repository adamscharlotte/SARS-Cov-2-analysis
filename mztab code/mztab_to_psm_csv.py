import sys
#print(sys.executable)       #sys.executable contains full path of the currently running Python interpreter
import pandas as pd
import sys
import argparse, pathlib

parser = argparse.ArgumentParser()
parser.add_argument('mztab', type=pathlib.Path)
parser.add_argument('bait', type=str)
args = parser.parse_args()

# path = "/Users/adams/Documents/master/master 2/Stage en thesis/Data/sars_cov_2/peptideIndexer/"
# name = "AGK--_-_Partition_49_of_200"
# split = name.split("_")
# bait = split[1]
# split = name.split("--")
# bait = split[0]
# input_path = path + 'mztab/' + bait + '.mztab'

# output_path = path + 'csv/' + args.bait + '.csv'
output_path = args.bait + '.csv'

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
            #mztabList.append(line)