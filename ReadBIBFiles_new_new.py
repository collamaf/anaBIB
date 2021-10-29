import glob
import argparse
import sys
import numpy as np

def is_interactive():
    import __main__ as main
    return not hasattr(main, '__file__')


parser = argparse.ArgumentParser(description='Read data path')
parser.add_argument('--dataPath', type=str, help='input file path')
parser.add_argument('--outPath', type=str, help='out file path')


if is_interactive():
    sys.argv = ['-f']

args = parser.parse_args()
if args.dataPath:
    path=args.dataPath
else:
    path = "FLUKAresults/provanew" #/*"
    #path = "FLUKAresults/3TeVFiles_BH_WM/*"

if args.outPath:
    outfile=args.outPath
else:
    #outfile = "FLUKAresults/aa"
    outfile1 = "FLUKAresults/3TeV_point_BH_WM_NEW"
    outfile2 = "FLUKAresults/3TeV_point_BH_WM_NEW_electrons"


numMuonsTot=0
for filename in sorted(glob.glob(path)):
    with open(filename, 'r') as f:
        print("Aperto File", filename)
        for line in f:
            if  line.strip():
                if line[1] == '#':
                    if line[4:8] == 'numb':
                        numMuons=int(line.split()[4])
                        numMuonsTot+=numMuons
                        continue
                    else:
                        continue
print("Eventi TOT=", numMuonsTot)


with open(outfile1,"w") as fOut1, open(outfile2,"w") as fOut2:
    for filename in sorted(glob.glob(path)):
        with open(filename, 'r') as f:
            print("Aperto File", filename)
            for line in f:
                if  line.strip():
                    if line[1:5] == '## f':
                        ele=line.split()
                        fOut2.write("{:d}\t{:e}\t{:e}\t{:e}\n"
                                   .format(int(ele[4]), float(ele[5]), float(ele[6]), float(ele[7])))
                    if line[0:2] == '  ':
                        ele=line.split()
                        fOut1.write("{:d}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\n"
                                   .format(int(ele[0]), float(ele[1]), float(ele[2]), float(ele[3]), float(ele[4]), float(ele[5])*2e12/numMuonsTot, float(ele[6]),
                                           float(ele[7]), float(ele[8]), float(ele[9]), float(ele[10]), float(ele[11]), float(ele[12]),
                                           float(ele[13]), float(ele[14]), float(ele[15]), float(ele[16]), float(ele[17]), float(ele[18]), float(ele[19])))
                    else:
                        continue
