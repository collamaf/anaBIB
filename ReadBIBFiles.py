#!/usr/bin/env python
# coding: utf-8

# In[6]:


import glob
import argparse
parser = argparse.ArgumentParser(description='Read data path')
parser.add_argument('--dataPath', type=str, help='path')
#path = "FLUKAresults/3TeV/prova3tev*
args = parser.parse_args()

#path = "FLUKAresults/1.5TeV/*"
if args.dataPath:
    path=args.dataPath
else:
    path = "FLUKAresults/1.5TeV/*"


# In[3]:


numMuonsTot=0
for filename in sorted(glob.glob(path)):
    with open(filename, 'r') as f:
        print("Aperto File", filename)
        for line in f:
            if  line.strip():
#        line=f.readline()
#            print(line[0]+"\n")
                if line[1] == '#':
                    if line[4:8] == 'numb':
                        numMuons=int(line.split()[4])
                        numMuonsTot+=numMuons
                        continue
                    else:
                        continue
print("Eventi TOT=", numMuonsTot)


# In[4]:


with open("OUT1p5","w") as fOut:
    for filename in sorted(glob.glob(path)):
        with open(filename, 'r') as f:
            print("Aperto File", filename)
            for line in f:
                if  line.strip():
                    if line[1] == '#':
                        continue
                    else:
                        ele=line.split()
                        fOut.write("{:d}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{}\t{:e}\t{:e}\t{:e}\t{:d}\t{}\t{:d}\n"
                                   .format(int(ele[0]), float(ele[1]), float(ele[2]), float(ele[3]), float(ele[4]), float(ele[5])*2e12/numMuonsTot, float(ele[6]), 
                                           float(ele[7]), float(ele[8]), float(ele[9]), ele[10], float(ele[11]), float(ele[12]), 
                                           float(ele[13]), int(ele[14]), ele[15], int(ele[16]) ))


# In[ ]:





# In[ ]:




