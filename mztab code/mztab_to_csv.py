import sys
#print(sys.executable)       #sys.executable contains full path of the currently running Python interpreter
import pandas as pd
import sys
import argparse, pathlib

parser = argparse.ArgumentParser()
parser.add_argument('mztab', type=pathlib.Path)
parser.add_argument('bait', type=str)
args = parser.parse_args()

#path = "/Users/adams/Documents/master/master 2/Stage en thesis/Data/MassIVE/26.11.2019/"
#name = "AGK--_-_Partition_49_of_200"

#split = bait.split("--")
#bait = split[0]

#   CWL output
#output_path = path + 'csv/' + bait + '.csv'
output_path = args.bait + '.csv'

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