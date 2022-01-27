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
    path = "3TeV_base_pointa001_DUMP"

if args.outPath:
    outfile=args.outPath
else:
    outfile = "prova"


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

count=0
with open(outfile,"w") as fOut1, open(outfile+'_ele',"w") as fOut2:
    for filename in sorted(glob.glob(path)):
        with open(filename, 'r') as f:
            print("Aperto File", filename)
            for line in f:
                if  line[8:12] == 'muon' and count<1:
                    nmd=line.split()
                    count=2
                if  line.strip():
                    if line[1:5] == '## f':
                        ele=line.split()
                        if float(ele[12])==float(ele[13])==float(ele[14])==0.00000000:
                            ele[14]=-2000.
                        fOut2.write("{:d}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\n"
                                   .format(int(ele[4]), float(ele[5]), float(ele[6]), float(ele[7]),float(ele[8]), float(ele[9]), float(ele[10]),float(ele[11]), float(ele[12]),
                                           float(ele[13]),float(ele[14]), float(nmd[7])*2e12/numMuonsTot))
                    if line[0:2] == '  ':
                        ele=line.split()
                        if float(ele[5])<1.:
                            fOut1.write("{:d}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\n"
                                   .format(int(ele[0]), float(ele[1]), float(ele[2]), float(ele[3]), float(ele[4]), float(ele[5])*2e12/numMuonsTot, float(ele[6]),
                                           float(ele[7]), float(ele[8]), float(ele[9]), float(ele[10]), float(ele[11]), float(ele[12]),
                                           float(ele[13]), float(ele[14]), float(ele[15]), float(ele[16]), float(ele[17]), float(ele[18]), float(ele[19])))
                    else:
                        continue
