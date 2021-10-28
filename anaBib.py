#!/usr/bin/env python
# coding: utf-8

# # BIB ANALYSIS

# ## Imports

# In[1]:


import math
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import pylab as pl
import math
import random
import collections, numpy
import glob
import argparse
import sys
from matplotlib.pyplot import pie, axis, show
from particle import PDGID
from particle import Particle
#plt.rcParams['figure.dpi'] = 300
#plt.rcParams['savefig.dpi'] = 300

matplotlib.rcParams.update({'font.size': 14})

def is_interactive():
    import __main__ as main
    return not hasattr(main, '__file__')

parser = argparse.ArgumentParser(description='Read data path')
parser.add_argument('--runName', type=str, help='run name')
parser.add_argument('--fileList', nargs='+', help='input file or files')
parser.add_argument('--labelList', nargs='+', help='file label o labels')

if is_interactive():
    sys.argv = ['-f']
    
args = parser.parse_args()
if args.fileList:
    inputFilesList=args.fileList
    labelList=args.labelList
else:
    inputFilesList=["DigFiles/NozzleModNorm", "DigFiles/NozzleMod"]
    labelList=["Norm", "Mod"]
if args.runName:
    runName=args.runName+"_"
else:
    runName="BIB_TEST_"

print("Leggo Files: ", inputFilesList)


# ## Initial Flags and Variables

# In[2]:


flagApplyPaperEnCut=False
flagApplyZCut=True

listChargedHadrons=[321, 311, 211, 2212, 3122, 3112, -321, -311, -211, -2212, -3122, -3112]

nbins=50
nbinsH=200
nbinsZ=100
binwidth=0.8 #Unit?
binWidthZ=1 # [m]

colsToRead=["PDGcode", "KinE", "PX","PY","PZ","Weight","PosX","PosY","PosZ", "Time", "Elem","PosXmu","PosYmu","PosZmu","ind1","Elem2","ind2"]

# Energy Cuts (to reproduce MAP published results)
enCutPh=0.2e-3
enCutNeu=0.1e-3
enCutElPos=0.2e-3
enCutChHad=1e-3
enCutMu=1e-3

zCut=2500. #in cm

#inputFilesList=["DigFiles/Part1.5TeV.dump", "DigFiles/Part3TeV.dump"]
#labelList=["1p5TeV", "3TeV"]
#inputFilesList=["DigFiles/NozzleModNorm", "DigFiles/NozzleMod"]
#labelList=["Norm", "Mod"]

#inputFilesList=TEST
#labelList=["Norm", "Mod"]


# In[3]:


datasetList=[]


# ## Utility Functions

# In[44]:


def plot_arrays(array1, array2=None, label1="", label2="", title="", array3=None, label3=""):
    plt.figure(1, figsize= (14,11))
    ax1=plt.subplot(211)
    ax1.hist(array1, bins=100, color='r', label =label1, histtype='step') 
    if (array2 is not None):
        ax1.hist(array2.flatten(), bins=100, color='g',  label =label2, histtype='step')
    if (array3 is not None):
        ax1.hist(array3.flatten(), bins=100, color='b',  label =label3, histtype='step', linestyle='-.')
    ax1.set(xlabel='HU', ylabel='[#]', title=title)
    ax1.legend()
    plt.yscale('log')
    plt.show()

def plot1D(ax,x,plotTitle="", label="",xlabel="x",ylabel="y",log=True,col="r", weights=None,bins=None,rng=None, numPart=None, ls="solid"):
    ax.set_xlabel(xlabel,fontsize='14')
    ax.set_ylabel(ylabel,fontsize='14')
    ax.set_title(plotTitle)
   
    if log:
        ax.set_ylim(auto=True)
    if numPart:
        ax.hist(x,log=log,histtype='step', label=label+str.format(r' N={:.2e} $\bar x$={:.2e}',numPart if numPart else 0,x.mean()), bins=bins, range=rng, weights=weights, linestyle=ls)
    else:
        ax.hist(x,log=log,histtype='step', label=label+str.format(r' $\bar x$={:.2e}',x.mean()), bins=bins, range=rng, weights=weights, linestyle=ls)


#    ax.hist(x,log=log,histtype='step', label=label+str.format(r' N={:.2e} $\bar x$={:.2e}',numPart if numPart else 0,x.mean()), color=col, bins=bins, range=rng, weights=weights)

#    ax.hist(x,log=log,histtype='step', label=str.format(r'N={:.2e} $\bar x$={:.2e}',numPart if numPart else 0,x.mean()), color=col, bins=bins, range=rng, weights=weights)


    ax.legend(loc= "best")
    ax.grid(True)
    return ax

def plot1Dold(ax,data,label="",xlabel="x",ylabel="y",log=False,col="r", weights=None,bins=None,rng=None, numPart=None):
    ax.set_xlabel(xlabel,fontsize='14')
    ax.set_ylabel(ylabel,fontsize='14')
    if numPart:
        if numPart==0:
            ax.set_title('NO PARTS ')
        else:
#            ax.set_title('N of particles '+ str.format('{:.2e}',numPart)+"\n Mean "+ str.format('{:.2e} MeV',x.mean()*1e3),fontsize=12)
            ax.set_title('N of particles '+ str.format('{:.2e}',numPart),fontsize=12)

    if log:
        ax.set_ylim(auto=True)
    ax.hist(data,log=log,histtype='step', label=label+ str.format(r'  $\bar x$= {:.2e}',data.mean()), color=col, bins=bins, range=rng, weights=weights)
    ax.legend(loc= "best")
    ax.grid(True)
    return ax

def getMomentumOld(part, absFlag=False, bFlag=False):
    if bFlag:
        if abs:
            return numpy.sqrt(datasetB[abs(datasetB["PDGcode"])==part]["PX"]**2+datasetB[abs(datasetB["PDGcode"])==part]["PY"]**2+datasetB[abs(datasetB["PDGcode"])==part]["PZ"]**2)
        else:
            return numpy.sqrt(datasetB[datasetB["PDGcode"]==part]["PX"]**2+datasetB[datasetB["PDGcode"]==part]["PY"]**2+datasetB[datasetB["PDGcode"]==part]["PZ"]**2)  
    else:
        if abs:
            return numpy.sqrt(dataset[abs(dataset["PDGcode"])==part]["PX"]**2+dataset[abs(dataset["PDGcode"])==part]["PY"]**2+dataset[abs(dataset["PDGcode"])==part]["PZ"]**2)
        else:
            return numpy.sqrt(dataset[dataset["PDGcode"]==part]["PX"]**2+dataset[dataset["PDGcode"]==part]["PY"]**2+dataset[dataset["PDGcode"]==part]["PZ"]**2)

        
def getMomentum2(data, part, absFlag=False):
    if abs:
        return numpy.sqrt(data[abs(data["PDGcode"])==part]["PX"]**2+data[abs(data["PDGcode"])==part]["PY"]**2+data[abs(data["PDGcode"])==part]["PZ"]**2)
    else:
        return numpy.sqrt(data[data["PDGcode"]==part]["PX"]**2+data[data["PDGcode"]==part]["PY"]**2+data[data["PDGcode"]==part]["PZ"]**2)

def getInfo(dataset, part, info):
    if isinstance(part, int):
        return dataset[dataset["PDGcode"]==part][info]
    else:
        if isinstance(part, list):
            return dataset[dataset["PDGcode"].isin(part)][info]
    
def getMomentum(data, part, absFlag=False):
    if isinstance(part, int):
        return numpy.sqrt(data[data["PDGcode"]==part]["PX"]**2+data[data["PDGcode"]==part]["PY"]**2+data[data["PDGcode"]==part]["PZ"]**2)
    else:
        if isinstance(part, list):
            return numpy.sqrt(data[data["PDGcode"].isin(part)]["PX"]**2+data[data["PDGcode"].isin(part)]["PY"]**2+data[data["PDGcode"].isin(part)]["PZ"]**2)

def getParticleNumber(dataset, part):
    if isinstance(part, int):
        return sum(dataset[dataset["PDGcode"]==part]["Weight"])
    else:
        if isinstance(part, list):
            return sum(dataset[dataset["PDGcode"].isin(part)]["Weight"])
        
def drawPie2(var, figName, bFlag=False, title=""):
    fig=plt.figure(figsize=(5,5))
    plt.title(title)
    if bFlag:
        sums = dataset.groupby(dataset[var])["Weight"].sum()
    else:
        sums = datasetB.groupby(datasetB[var])["Weight"].sum()    
    axis('equal');
    cmap = plt.get_cmap('Spectral')
    colors = [cmap(i) for i in np.linspace(0, 1, 8)]
    pie(sums, autopct='%1.1f%%',labels=sums.index,colors=colors)
    figname=figName+str(title)
    pl.savefig(figname)

            
def drawPie(dataset, var, title=""):
    fig=plt.figure(figsize=(5,5))
    plt.title(runName+title)
    sums = dataset.groupby(dataset[var])["Weight"].sum()    
    axis('equal');
    cmap = plt.get_cmap('Spectral')
    colors = [cmap(i) for i in np.linspace(0, 1, 8)]
    pie(sums, autopct='%1.1f%%',labels=sums.index,colors=colors)
#    figname=+str(title)
    pl.savefig(runName+title)
    

    
def plotDistribution(datasetList, variable, plotTitle="", xlabel="", ylabel="Arb. Units", nbins=nbins, log=True, figTitle="", xrange=None, ymax=None):
    fig, ax = plt.subplots(nrows=1, ncols=len(datasetList)+1, figsize=((len(datasetList)+1)*8,8), sharey=False)
    plt.suptitle(runName+plotTitle)
    
    for i, dataset in enumerate(datasetList):
        ax[i].set_xlabel(xlabel,fontsize='14')
        ax[i].set_ylabel(ylabel,fontsize='14')

        temp=ax[i].hist(getInfo(dataset, 22, variable),histtype='step', bins=nbins, weights=getInfo(dataset, 22, "Weight"), log=log, range=xrange, label="$\gamma$"+str.format(r' N={:.2e} $\bar x$={:.2e}',getParticleNumber(dataset,22),getInfo(dataset,22,variable).mean()))
        ax[i].hist(getInfo(dataset, [-11,11], variable),histtype='step', bins=nbins, weights=getInfo(dataset, [11,-11], "Weight"), log=log, range=xrange, label="e+e-"+str.format(r' N={:.2e} $\bar x$={:.2e}',getParticleNumber(dataset,[11,-11]),getInfo(dataset,[11,-11],variable).mean()))
        ax[i].hist(getInfo(dataset, listChargedHadrons, variable),histtype='step', bins=nbins, weights=getInfo(dataset, listChargedHadrons, "Weight"), log=log, range=xrange, label="Ch. Had"+str.format(r' N={:.2e} $\bar x$={:.2e}',getParticleNumber(dataset,listChargedHadrons),getInfo(dataset,listChargedHadrons,variable).mean()))
        ax[i].hist(getInfo(dataset, 2112, variable),histtype='step', bins=nbins, weights=getInfo(dataset, 2112, "Weight"), log=log, range=xrange, label="Neutrons"+str.format(r' N={:.2e} $\bar x$={:.2e}',getParticleNumber(dataset,2112),getInfo(dataset,2112,variable).mean()))
        ax[i].hist(getInfo(dataset, [-13,13], variable),histtype='step', bins=nbins, weights=getInfo(dataset, [13,-13], "Weight"), log=log, range=xrange, label="Mu+Mu-"+str.format(r' N={:.2e} $\bar x$={:.2e}',getParticleNumber(dataset,[13,-13]),getInfo(dataset,[13,-13],variable).mean()))

        if i==0:
            maxHeight=temp[0].max() #Calcolo l'altezza massima del bin in modo da forzare gli N grafici ad avere la stessa scala verticale facilitando il confronto
            #print(maxHeight)
        ax[i].set_title(labelList[i])
        ax[i].legend()
        if ymax!=None:
            ax[i].axis(ymin=1e1, ymax=ymax)
        else:
            ax[i].axis(ymin=1e1, ymax=maxHeight*2)

        ax[len(datasetList)].hist(getInfo(dataset, 22, variable),histtype='step', bins=nbins, weights=getInfo(dataset, 22, "Weight"), log=log, range=xrange, label=labelList[i]+" $\gamma$")
        ax[len(datasetList)].hist(getInfo(dataset, [-11,11], variable),histtype='step', bins=nbins, weights=getInfo(dataset, [11,-11], "Weight"), log=log, range=xrange, label=labelList[i]+" e+e-")
        ax[len(datasetList)].hist(getInfo(dataset, listChargedHadrons, variable),histtype='step', bins=nbins, weights=getInfo(dataset, listChargedHadrons, "Weight"), log=log, range=xrange, label=labelList[i]+" Ch. Had")
        ax[len(datasetList)].hist(getInfo(dataset, 2112, variable),histtype='step', bins=nbins, weights=getInfo(dataset, 2112, "Weight"), log=log, range=xrange, label=labelList[i]+" Neutrons")
        ax[len(datasetList)].hist(getInfo(dataset, [-13,13], variable),histtype='step', bins=nbins, weights=getInfo(dataset, [13,-13], "Weight"), log=log, range=xrange, label=labelList[i]+" Mu+Mu-")

        if ymax!=None:
            ax[len(datasetList)].axis(ymin=1e1, ymax=ymax)
        else:
            ax[len(datasetList)].axis(ymin=1e1, ymax=maxHeight*2)
        
    ax[len(datasetList)].set_title("Comparison")
    ax[len(datasetList)].legend()
    ax[len(datasetList)].set_xlabel(xlabel,fontsize='14')
    ax[len(datasetList)].set_ylabel(ylabel,fontsize='14')

    figname=runName+figTitle
    pl.savefig(figname)
        
    
def scatter_histo(x, y, ax, ax_histx, ax_histy, weights=None, xlabel="", ylabel="", xrange=[-750, 750], yrange=[-30,30]):
    ax_histx.tick_params(axis="x", labelbottom=False)
    ax_histy.tick_params(axis="y", labelleft=False)
    ax.hist2d(x, y,weights=weights,bins=100,norm=matplotlib.colors.LogNorm(), cmap='Blues')
    ax_histx.hist(x,log=True, weights=weights,histtype='step',bins=100,rwidth=binwidth,color='b')
    ax_histy.hist(y,log=True, weights=weights,histtype='step',bins=100,rwidth=binwidth,color='b', orientation='horizontal')
    plt.gca().invert_yaxis()
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    
def scatter_histo2(x, y, ax, ax_histx, ax_histy, weights=None, xlabel="", ylabel="", xrange=[-750, 750], yrange=[-30,30]):
    ax_histx.tick_params(axis="x", labelbottom=False)
    ax_histy.tick_params(axis="y", labelleft=False)
    ax.hist2d(x, y,weights=weights,range=[xrange, yrange],bins=100,norm=matplotlib.colors.LogNorm(), cmap='Blues')
    ax_histx.hist(x,log=True, weights=weights,histtype='step',range=(-750,750),bins=100,rwidth=binwidth,color='b')
    ax_histy.hist(y,log=True, weights=weights,histtype='step',range=(-30,30),bins=100,rwidth=binwidth,color='b', orientation='horizontal')
    plt.gca().invert_yaxis()
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    
def scatter_histoOld(x, y, ax, ax_histx, ax_histy, weights=None, xlabel="", ylabel=""):
    ax_histx.tick_params(axis="x", labelbottom=False)
    ax_histy.tick_params(axis="y", labelleft=False)
    ax.hist2d(x, y,weights=weights,range=[[-750, 750], [-30, 30]],bins=100,norm=matplotlib.colors.LogNorm(), cmap='Blues')
    ax_histx.hist(x,log=True, weights=weights,histtype='step',range=(-750,750),bins=100,rwidth=binwidth,color='b')
    ax_histy.hist(y,log=True, weights=weights,histtype='step',range=(-30,30),bins=100,rwidth=binwidth,color='b', orientation='horizontal')
    plt.gca().invert_yaxis()
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

left, width = 0.1, 0.65
bottom, height = 0.1, 0.65
spacing = 0.005
rect_scatter = [left, bottom, width, height]
rect_histx = [left, bottom + height + spacing, width, 0.2]
rect_histy = [left + width + spacing, bottom, 0.2, height]


# ## Read Datasets

# In[5]:


for fileNumber, fileName in enumerate(inputFilesList):
    print("Read file: ",fileNumber, fileName)
    temp=pd.read_csv(fileName, header=None, names=colsToRead, delim_whitespace=True)
    datasetList.append(temp)


# ------------------------------------------

# ## Let's have a look at muon decay z position

# In[6]:


nBinZ=[]
histoCumA=[]
histoCumAnorm=[]


# In[7]:


fig, axs = plt.subplots(nrows=4, ncols=1, figsize=(10,12), sharex=True)
axs2b=axs[3].twinx() #Same x, but different y scale

for datasetNumber, dataset in enumerate(datasetList):
    nBinZ.append(int(((dataset["PosZmu"]/100).max()-(dataset["PosZmu"]/100).min())/binWidthZ))
    plot1D(axs[0],dataset["PosZmu"]/100, weights=dataset["Weight"],bins=nBinZ[-1], plotTitle="Muon Decay Z", label=labelList[datasetNumber], xlabel='', ylabel='Arb. Units' )
    axs[1].hist(dataset["PosZmu"]/100, weights=dataset["Weight"], bins=nBinZ[-1], cumulative=1,histtype='step', label="cum "+labelList[datasetNumber], density=False)
    histoCumAnorm.append(axs[2].hist(dataset["PosZmu"]/100, weights=dataset["Weight"], bins=nBinZ[-1], cumulative=1,histtype='step', label="cum "+labelList[datasetNumber], density=True))
    plot1D(axs[3],dataset["PosZmu"]/100, weights=dataset["Weight"],bins=nBinZ[-1], plotTitle="Muon Decay Z with Cumulatives", label=labelList[datasetNumber], xlabel='$z_{\mu \,dec}$ (m)', ylabel='Arb. Units' )
    histoCumA=axs2b.hist(dataset["PosZmu"]/100, weights=dataset["Weight"], bins=nBinZ[-1], cumulative=1,histtype='step', label="cum "+labelList[datasetNumber], density=False, linestyle=':')


#axs[0].axis(ymin=100, ymax=1e9)
axs[0].grid(True, which="both", axis='y')
axs[0].locator_params(axis="x", nbins=20)
axs[2].grid(True, which="both")
axs[2].set_title("Cumulative Function Norm")
axs[2].legend(loc= "best", fontsize='x-small')

axs[1].legend(loc= "best", fontsize='x-small')
axs[1].set_title("Cumulative Function Not Norm")
axs[1].grid(True, which="both")

figname=runName+"MuDecZ"
pl.savefig(figname)


# ### If Needed, let's apply a z-cut

# In[8]:


for i, dataset in enumerate(datasetList):
    if flagApplyZCut:
        print("### Z CUT REQUESTED! ", zCut)
        datasetList[i]=datasetList[i][datasetList[i]["PosZmu"]<zCut] #in cm
    if flagApplyPaperEnCut:
        print("### PAPER CUT REQUESTED!")
        datasetList[i]=datasetList[i].drop(datasetList[i][(datasetList[i]["PDGcode"] == 22) & (datasetList[i]["KinE"] < enCutPh)].index)
        datasetList[i]=datasetList[i].drop(datasetList[i][((datasetList[i]["PDGcode"]==11)|(datasetList[i]["PDGcode"]==-11)) & (datasetList[i]["KinE"] < enCutElPos)].index)
        datasetList[i]=datasetList[i].drop(datasetList[i][((datasetList[i]["PDGcode"]).isin(listChargedHadrons)) & (dataset["KinE"] < enCutChHad)].index)
        datasetList[i]=datasetList[i].drop(datasetList[i][((datasetList[i]["PDGcode"]==13)|(datasetList[i]["PDGcode"]==-13)) & (datasetList[i]["KinE"] < enCutMu)].index)
        datasetList[i]=datasetList[i].drop(datasetList[i][(datasetList[i]["PDGcode"] == 2112) & (datasetList[i]["KinE"] < enCutNeu)].index)


# ## List of found particles' IDs

# In[9]:


foundParticlesList=[]
foundParticlesFreqList=[]
foundParticlesUnique=set([])
foundParticlesUniqueFreq=[]
foundParticlesUniqueEntries=[]


# ### Creiamo una lista con tutte e sole le particelle che sono comparse in almeno un dataset

# In[10]:


for dataset in datasetList:
    (lista, freq)=np.unique(dataset["PDGcode"], return_counts=True)
    foundParticlesList.append(lista)
    foundParticlesFreqList.append(freq)
    foundParticlesUnique=foundParticlesUnique.union(foundParticlesUnique,lista)
    
foundParticlesUnique=list(foundParticlesUnique)
foundParticlesUnique.sort()
print("List of all (and unique) found particles in all datasets\n", foundParticlesUnique)


# ### Recupero i nomi di queste particelle

# In[11]:


unknownParticle="??"

particleNamesList=[]
for particle in foundParticlesUnique:
    #print(PDGID(particle))
    try: 
        Particle.from_pdgid(particle)
        #print(Particle.from_pdgid(particle))
        particleNamesList.append(Particle.from_pdgid(particle).name)
    except:
        #print("non trovata")
        particleNamesList.append(unknownParticle)
print("List of all (and unique) found particles in all datasets\n", particleNamesList)


# ### Per ogni particella, calcoliamo la frequenza con cui è comparsa in ciascun dataset

# In[12]:


for i, dataset in enumerate(datasetList):
    foundParticlesUniqueEntries.append([]) # Aggiungo una lista nuova (e vuota) per ciascun dataset
    for particle in foundParticlesUnique: # Giro su tutte le particelle della lista univoca
        if not dataset[dataset["PDGcode"]==particle].empty: # Se in questo dataset quela particella è comparsa ne calcolo il peso
            foundParticlesUniqueEntries[i].append(sum(dataset[dataset["PDGcode"]==particle]["Weight"]))
        else:
            foundParticlesUniqueEntries[i].append(0)
        #print("Per particella {} ho aggiunto {}".format(particle,foundParticlesUniqueEntries[i][-1]))
 #   print("Dataset {} - List of particles entries: {}\n\n".format(i, np.around(foundParticlesUniqueEntries[i],2)))


# In[59]:


fig, axs = plt.subplots(nrows=len(datasetList)+1, ncols=1, figsize=(18,len(datasetList)*8))
fig.suptitle(runName+"Particles Frequencies")
for i, dataset in enumerate(datasetList):
    axs[i].bar(range(len(foundParticlesUnique)), foundParticlesUniqueEntries[i], align='center', log=True)
    axs[i].set_title(labelList[i]) 
    plt.sca(axs[i])
    plt.xticks(range(len(foundParticlesUnique)), particleNamesList, size='large')
    plt.xticks(rotation=90)
    plt.ylabel("Occurency", size='large')
    
    for j, v in enumerate(foundParticlesUniqueEntries[i]):
        if v!=0:
            axs[i].text(j , v, "{:.2e}".format(v))
            #axs[len(datasetList)].text(j , v, "{:.2e}".format(v), color="tab:blue")
            
    axs[len(datasetList)].bar(range(len(foundParticlesUnique)), foundParticlesUniqueEntries[i], align='center', log=True, label=labelList[i], alpha=0.5)
    axs[len(datasetList)].set_title('Comparison')
    plt.sca(axs[len(datasetList)])
    plt.xticks(range(len(foundParticlesUnique)), particleNamesList, size='large')
    plt.xticks(rotation=90)
    plt.ylabel("Occurency", size='large')
    plt.legend()
    
fig.subplots_adjust(top=0.93)
#fig.tight_layout()
plt.subplots_adjust(hspace = 0.3)
figname=runName+"ParticleDistribution"
pl.savefig(figname)


# ## Count Particle Numbers

# In[14]:


for i, dataset in enumerate(datasetList):
    print("DATASET ", labelList[i])
    print("N Photons\t{:.2e}".format(getParticleNumber(dataset,22)))
    print("N Positrons\t{:.2e}".format(getParticleNumber(dataset,-11)))
    print("N Electrons\t{:.2e}".format(getParticleNumber(dataset,11)))
    print("N ElePos\t{:.2e}".format(getParticleNumber(dataset,[-11,11])))
    print("N Protons\t{:.2e}".format(getParticleNumber(dataset,2212)))
    print("N Neutrons\t{:.2e}".format(getParticleNumber(dataset,2112)))
    print("N CharHad\t{:.2e}".format(getParticleNumber(dataset,listChargedHadrons)))
    print("N MuonPlus\t{:.2e}".format(getParticleNumber(dataset,-13)))
    print("N MuonMin\t{:.2e}".format(getParticleNumber(dataset,13)))
    print("N MuonPM\t{:.2e}".format(getParticleNumber(dataset,[13,-13])))
    print()


# ## Plot Energy Spectra

# ### All relevant particles

# In[27]:


fig, axs = plt.subplots(nrows=1, ncols=5, figsize=(22,5), sharey=False)
plt.suptitle(runName+"Energy Spectra")

for i, dataset in enumerate(datasetList):
    axs[0].hist(getInfo(dataset,22,"KinE"),histtype='step', bins=nbins, weights=getInfo(dataset, 22, "Weight"), range=[0,np.percentile(getInfo(datasetList[0],22,"KinE").to_numpy(),99.999)], log=True, label=labelList[i]+str.format(r' N={:.2e} $\bar x$={:.2e}',getParticleNumber(dataset,22),getInfo(dataset,22,"KinE").mean()))
    axs[1].hist(getInfo(dataset,2212,"KinE"),histtype='step', bins=nbins, weights=getInfo(dataset, 2212, "Weight"), range=[0,np.percentile(getInfo(datasetList[0],2212,"KinE").to_numpy(),99.999)], log=True, label=labelList[i]+str.format(r' N={:.2e} $\bar x$={:.2e}',getParticleNumber(dataset,2212),getInfo(dataset,2212,"KinE").mean()))
    axs[2].hist(getInfo(dataset,[11,-11],"KinE"),histtype='step', bins=nbins, weights=getInfo(dataset, [11,-11], "Weight"), range=[0,np.percentile(getInfo(datasetList[0],[11,-11],"KinE").to_numpy(),99.999)], log=True, label=labelList[i]+str.format(r' N={:.2e} $\bar x$={:.2e}',getParticleNumber(dataset,[11,-11]),getInfo(dataset,[11,11],"KinE").mean()))
    axs[3].hist(getInfo(dataset,listChargedHadrons,"KinE"),histtype='step', bins=nbins, weights=getInfo(dataset, listChargedHadrons, "Weight"), range=[0,np.percentile(getInfo(datasetList[0],listChargedHadrons,"KinE").to_numpy(),99.999)], log=True, label=labelList[i]+str.format(r' N={:.2e} $\bar x$={:.2e}',getParticleNumber(dataset,listChargedHadrons),getInfo(dataset,listChargedHadrons,"KinE").mean()))
    axs[4].hist(getInfo(dataset,[13,-13],"KinE"),histtype='step', bins=nbins, weights=getInfo(dataset, [13,-13], "Weight"), range=[0,np.percentile(getInfo(datasetList[0],[13,-13],"KinE").to_numpy(),99.999)], log=True, label=labelList[i]+str.format(r' N={:.2e} $\bar x$={:.2e}',getParticleNumber(dataset,[13,-13]),getInfo(dataset,[13,-13],"KinE").mean()))

for i in range(0,5):
    axs[i].set_ylim(auto=True)
    axs[i].legend(loc= "best")
    axs[i].grid(True)
    axs[i].set_xlabel("$E_{kin}$ [GeV]",fontsize='14')
    axs[i].set_ylabel("Arb. Units",fontsize='14')
axs[0].set_title("$\gamma$")
axs[1].set_title("n")
axs[2].set_title("$e^-~/~e^+$")
axs[3].set_title("ch. had")
axs[4].set_title("$\mu^-~/~\mu^+$")

fig.subplots_adjust(left = 0.05,right = 0.99,wspace = 0.5, hspace = 0.5, top=0.85, bottom= 0.2)

figname=runName+"EnergySpectra"
pl.savefig(figname)


# ### Photons and e+/e-

# In[31]:


fig, ax = plt.subplots(nrows=1, ncols=len(datasetList)+1, figsize=((len(datasetList)+1)*8,8), sharey=False)
plt.suptitle(runName+"E.M. Components")

for i, dataset in enumerate(datasetList):
    ax[i].set_xlabel("p [GeV/c]",fontsize='14')
    ax[i].set_ylabel("Arb. Units",fontsize='14')
    
    temp=ax[i].hist(getMomentum(dataset, 22),histtype='step', bins=nbinsH, weights=getInfo(dataset, 22, "Weight"), log=True, range=[0,np.percentile(getInfo(datasetList[0],22,"KinE").to_numpy(),99.999)], label="$\gamma$"+str.format(r' N={:.2e} $\bar x$={:.2e}',getParticleNumber(dataset,22),getInfo(dataset,22,"KinE").mean()))
    ax[i].hist(getMomentum(dataset, [11,-11]),histtype='step', bins=nbinsH, weights=getInfo(dataset, [11,-11], "Weight"), log=True, range=[0,np.percentile(getInfo(datasetList[0],22,"KinE").to_numpy(),99.999)], label="$e^-~/~e^+$"+str.format(r' N={:.2e} $\bar x$={:.2e}',getParticleNumber(dataset,[11,-11]),getInfo(dataset,[11,-11],"KinE").mean()))

    if i==0:
        maxHeight=temp[0].max() #Calcolo l'altezza massima del bin in modo da forzare gli N grafici ad avere la stessa scala verticale facilitando il confronto
        #print(maxHeight)
        
    ax[i].set_title(labelList[i])
    ax[i].legend()
    ax[i].axis(ymin=1e1, ymax=maxHeight*2)

    ax[len(datasetList)].hist(getMomentum(dataset, 22),histtype='step', bins=nbinsH, weights=getInfo(dataset, 22, "Weight"), log=True, range=[0,np.percentile(getInfo(datasetList[0],22,"KinE").to_numpy(),99.999)], label=labelList[i]+" $\gamma$"+str.format(r' N={:.2e} $\bar x$={:.2e}',getParticleNumber(dataset,22),getInfo(dataset,22,"KinE").mean()))
    ax[len(datasetList)].hist(getMomentum(dataset, [11,-11]),histtype='step', bins=nbinsH, weights=getInfo(dataset, [11,-11], "Weight"), log=True, range=[0,np.percentile(getInfo(datasetList[0],22,"KinE").to_numpy(),99.999)], label=labelList[i]+" $e^-~/~e^+$"+str.format(r' N={:.2e} $\bar x$={:.2e}',getParticleNumber(dataset,[11,-11]),getInfo(dataset,[11,-11],"KinE").mean()))
    ax[len(datasetList)].axis(ymin=1e1, ymax=maxHeight*2)

ax[len(datasetList)].set_title("Comparison")
ax[len(datasetList)].legend()
ax[len(datasetList)].set_xlabel("p [GeV/c]",fontsize='14')
ax[len(datasetList)].set_ylabel("Arb. Units",fontsize='14')

figname=runName+"EneEGamma"
pl.savefig(figname)


# ### Hadrons

# In[39]:


fig, ax = plt.subplots(nrows=1, ncols=len(datasetList)+1, figsize=((len(datasetList)+1)*8,8), sharey=False)
plt.suptitle(runName+"Hadrons")

for i, dataset in enumerate(datasetList):
    ax[i].set_xlabel("p [GeV/c]",fontsize='14')
    ax[i].set_ylabel("Arb. Units",fontsize='14')
    
    ax[i].hist(getMomentum(dataset, listChargedHadrons),histtype='step', bins=nbinsH, weights=getInfo(dataset, listChargedHadrons, "Weight"), log=True, range=[0,1], label="Ch. Had."+str.format(r' N={:.2e} $\bar x$={:.2e}',getParticleNumber(dataset,22),getInfo(dataset,22,"KinE").mean()))
    temp=ax[i].hist(getMomentum(dataset, 2112),histtype='step', bins=nbinsH, weights=getInfo(dataset, 2112, "Weight"), log=True, range=[0,1], label="Neutrons"+str.format(r' N={:.2e} $\bar x$={:.2e}',getParticleNumber(dataset,[11,-11]),getInfo(dataset,[11,-11],"KinE").mean()))

    if i==0:
        maxHeight=temp[0].max() #Calcolo l'altezza massima del bin in modo da forzare gli N grafici ad avere la stessa scala verticale facilitando il confronto
        #print(maxHeight)
    ax[i].set_title(labelList[i])
    ax[i].legend()
    ax[i].axis(ymin=1e1, ymax=maxHeight*2)
    
    ax[len(datasetList)].hist(getMomentum(dataset, listChargedHadrons),histtype='step', bins=nbinsH, weights=getInfo(dataset, listChargedHadrons, "Weight"), log=True, range=[0,1], label=labelList[i]+" Ch. Had."+str.format(r' N={:.2e} $\bar x$={:.2e}',getParticleNumber(dataset,22),getInfo(dataset,22,"KinE").mean()))
    ax[len(datasetList)].hist(getMomentum(dataset, 2112),histtype='step', bins=nbinsH, weights=getInfo(dataset, 2112, "Weight"), log=True, range=[0,1], label=labelList[i]+" Neutrons"+str.format(r' N={:.2e} $\bar x$={:.2e}',getParticleNumber(dataset,[11,-11]),getInfo(dataset,[11,-11],"KinE").mean()))
    ax[len(datasetList)].axis(ymin=1e1, ymax=maxHeight*2)


ax[len(datasetList)].set_title("Comparison")
ax[len(datasetList)].legend()
ax[len(datasetList)].set_xlabel("p [GeV/c]",fontsize='14')
ax[len(datasetList)].set_ylabel("Arb. Units",fontsize='14')

#plt.xlim((0,1))
figname=runName+"EneHadrons"
pl.savefig(figname)


# ## Plot Time Distributions

# In[40]:


plotDistribution(datasetList, "Time", "Time Distribution", "t [ns]", "Arb. Units", nbinsH, True, "Time", xrange=(-30,100))


# ## Pie Plot Charts

# In[45]:


for i, dataset in enumerate(datasetList):
    drawPie(dataset, "Elem","PieDet_"+labelList[i])


# In[46]:


for i, dataset in enumerate(datasetList):
    drawPie(dataset, "Elem2","PieFirstInt_"+labelList[i])


# ## Plot Muons' Decay Position

# ### Global

# In[48]:


fig=plt.figure(figsize=(6,5))
plt.suptitle(runName+"Muon Decay Z")
plt.gca().set_xlabel('$z_{\mu \,dec}$ [m]')
plt.gca().set_ylabel('Arb. Units')

for i, dataset in enumerate(datasetList):
    fig.gca().hist(dataset["PosZmu"]/100,histtype='step', bins=nbinsH, weights=dataset["Weight"], log=True, label=labelList[i])

plt.legend()
#plt.ylim((100, 1e9))
figname=runName+"MuDec"
pl.savefig(figname)


# ### Per Particle

# In[49]:


plotDistribution(datasetList, "PosZmu", "Muon Decay Z Per Particle", '$z_{\mu \,dec}$ [cm]', "Arb. Units", nbinsZ, True, "MuDecPart", ymax=1e7)


# In[50]:


fig, ax = plt.subplots(nrows=len(datasetList), ncols=1, figsize=(9,len(datasetList)*4))
plt.suptitle(runName+"Nozzle")
if len(datasetList)>1:
    for i, dataset in enumerate(datasetList):
        ax[i].set_title(labelList[i])
        ax[i].hist2d(dataset["PosZ"],dataset["PosX"],norm=matplotlib.colors.LogNorm(),bins=500, cmap='plasma')
        ax[i].axis(xmin=-800, xmax=800)
        ax[i].axis(ymin=-250, ymax=250)
        ax[i].set_ylabel('x [cm]',fontsize=14)
else:
    ax.set_title(labelList[i])
    ax.hist2d(dataset["PosZ"],dataset["PosX"],norm=matplotlib.colors.LogNorm(),bins=500, cmap='plasma')
    ax.axis(xmin=-800, xmax=800)
    ax.axis(ymin=-250, ymax=250)
    ax.set_ylabel('x [cm]',fontsize=14)

plt.xlabel('z (cm)',fontsize=14)
#plt.ylabel('x (cm)',fontsize=14)

figname=runName+"ZvsX_FLUKA"
pl.savefig(figname)


# ### Scatter Plots

# In[51]:


fig = plt.figure(figsize=(5,5))
fig.suptitle(runName+"photons")
ax = fig.add_axes(rect_scatter)
ax_histx = fig.add_axes(rect_histx, sharex=ax)
ax_histy = fig.add_axes(rect_histy, sharey=ax)

scatter_histo2(getInfo(datasetList[0], 22, "PosZ"), getInfo(datasetList[0], 22, "Time"), ax, ax_histx, ax_histy, weights=getInfo(datasetList[0], 22,"Weight"), xlabel="PosZ", ylabel="Time", xrange=[-400,400], yrange=[-20,20])
figname=runName+"decph"
pl.savefig(figname)


# In[ ]:




