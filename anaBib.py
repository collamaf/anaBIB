#!/usr/bin/env python
# coding: utf-8

# # BIB ANALYSIS

# ### Last Update: 16-2-2022 by collamaf

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
parser.add_argument('--ele', dest='ele', action='store_true')
parser.add_argument('--noele', dest='ele', action='store_false')
parser.add_argument('--allPlots', dest='allPlots', action='store_true')

parser.set_defaults(ele=True, allPlots=False)

if is_interactive():
    sys.argv = ['-f']
    
args = parser.parse_args()

flagReadEle=args.ele
flagAllPlots=args.allPlots


if args.fileList:
    inputFilesList=args.fileList
    labelList=args.labelList
else:
    #inputFilesList=["../Dump_new/MARSresults/MARS1e5TeVmupiu", "../Dump_new/MARSresults/MARS1e5TeVmumeno"]
    #inputFilesList=["local_data/PR_3TeV_real", "local_data/PR_3TeV_real_ok"]
    #inputFilesList=["local_data/NEW_1e5TeV_base_point", "../Dump_new/MARSresults/MARS1e5TeVmumeno"]
    inputFilesList=["local_data/CV_1e5TeV_base_SMALL", "local_data/CV_3TeV_base_SMALL"]
    #labelList=["MARS+", "MARS-"]
    #labelList=["FLUKA3TeVreal", "FLUKA3TeVreal"]
    #labelList=["FLUKA", "MARS"]
    labelList=["1.5TeV", "3TeV"]
if args.runName:
    runName=args.runName+"_"
else:
    #runName="1e5TeVMARS+vsMARS-_"
    #runName="FLUKA3TeVrealvsFLUKA3TeVreal_"
    #runName="1e5TeV_FLUKAvsMARS_"
    runName="1e5vs3bis_"

print("Leggo Files: ", inputFilesList, flagReadEle)


# ## Initial Flags and Variables

# In[2]:


flagApplyPaperEnCut=False
flagApplyZCut=False

listChargedHadrons=[321, 311, 211, 2212, 3122, 3112, -321, -311, -211, -2212, -3122, -3112]
lineStyles=['solid', 'dotted', 'dashed', 'dashdot' ]

nbins=50
nbinsH=200
nbinsZ=100
binwidth=0.8 #Unit?
binWidthZ=1 # [m]
nSlicesErrors=5

#colsToRead=["PDGcode", "KinE", "PX","PY","PZ","Weight","PosX","PosY","PosZ", "Time", "Elem","PosXmu","PosYmu","PosZmu","ind1","Elem2","ind2"]
colsToRead=["PDGcode", "KinE", "PX","PY","PZ","Weight","PosX","PosY","PosZ", "Time","PosXmu","PosYmu","PosZmu","EneEle","CX","CY","CZ", "PosXFI", "PosYFI", "PosZFI"]
colsToReadEle=["NumPart","PosXmu","PosYmu","PosZmu","EneEle","CX","CY","CZ", "PosXEle", "PosYEle", "PosZEle", "Weight"]

# Energy Cuts (to reproduce published results) - mod CC 11/1/22
enCutPh=0.1e-3
enCutNeu=1.0e-12
enCutElPos=0.1e-3
enCutChHad=0.1e-3
enCutMu=0.1e-3

zCut=2500. #in cm
timeCut=(-1,15) #in ns


# In[3]:


datasetList=[]
if flagReadEle:
    datasetEleList=[]


# ## Utility Functions

# In[83]:


def getInfo(dataset, part, info):
    if isinstance(part, int):
        return dataset[dataset["PDGcode"]==part][info]
    else:
        if isinstance(part, list):
            return dataset[dataset["PDGcode"].isin(part)][info]
        
def plot1D(ax,x,plotTitle="", label="",xlabel="x",ylabel="y",log=True,col=None, weights=None,bins=None,rng=None, numPart=None, ls="solid"):
    ax.set_xlabel(xlabel,fontsize='14')
    ax.set_ylabel(ylabel,fontsize='14')
    ax.set_title(plotTitle)
   
    if log:
        ax.set_ylim(auto=True)
    if numPart:
        ax.hist(x,log=log,histtype='step', label=label+str.format(r' N={:.2e} $\bar x$={:.2e}',numPart if numPart else 0,x.mean()), bins=bins, range=rng, weights=weights, linestyle=ls, color=col)
    else:
        ax.hist(x,log=log,histtype='step', label=label, bins=bins, range=rng, weights=weights, linestyle=ls, color=col)
    ax.legend(loc= "best")
    ax.grid(True)
    return ax
    
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
        
def getParticlesNumbersErrors(dataset, nSlices):
    errPartSlices=np.empty(len(foundParticlesUnique))
    df=dataset
    chunk=int(len(df)/nSlices)+1
    dfs={}
    numPartSlices=np.empty([nSlices,len(foundParticlesUnique)])
    for n in range((df.shape[0] // chunk + 1)):
        #print("fetta numero ", n)
        df_temp=df.iloc[n*chunk:(n+1)*chunk]
        df_temp=df_temp.reset_index(drop=True)
        dfs[n]=df_temp
        for iPart, particle in enumerate(foundParticlesUnique):
           # print(n,particle,sum(df_temp[df_temp["PDGcode"]==particle]["Weight"]))
            numPartSlices[n, iPart]=sum(df_temp[df_temp["PDGcode"]==particle]["Weight"])
    errPartSlices=numPartSlices.std(0)*nSlices
    #print(errPartSlices[iDataset])
    return errPartSlices
        
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
    pl.savefig(figname,transparent=False, facecolor='white')

            
def drawPie(dataset, var, title=""):
    fig=plt.figure(figsize=(5,5))
    plt.title(runName+title)
    sums = dataset.groupby(dataset[var])["Weight"].sum()    
    axis('equal');
    cmap = plt.get_cmap('Spectral')
    colors = [cmap(i) for i in np.linspace(0, 1, 8)]
    pie(sums, autopct='%1.1f%%',labels=sums.index,colors=colors)
#    figname=+str(title)
    pl.savefig(runName+title,transparent=False, facecolor='white')
    
def plotVariablePerEachRelevantParticle(datasetList, variable, plotTitle="", xlabel="", ylabel="Arb. Units", nbins=nbins, log=True, figTitle="", xrange=None, ymax=None,trange=None, alsoWithTime=False):
    ## This function plots a given variable for Gammas, e+e-, ch.had, n. and mu+mi-.
    ## A plot for each dataset is produced, plus one last plot with all datasets superimposed
    fig, ax = plt.subplots(nrows=1, ncols=len(datasetList)+1, figsize=((len(datasetList)+1)*5,5), sharey=False)

    if trange:
        plt.suptitle(plotTitle+str.format(' tmin={} [ns] tmax={} [ns]', trange[0],trange[1]))
    else:
        plt.suptitle(plotTitle)
    
    for i, dataset in enumerate(datasetList):
        ax[i].set_xlabel(xlabel,fontsize='14')
        ax[i].set_ylabel(ylabel,fontsize='14')
        if trange and alsoWithTime:
            temp=ax[i].hist(getInfo(dataset, 22, variable),histtype='step', bins=nbins, weights=getInfo(dataset, 22, "Weight"), log=log, range=xrange,color='r', alpha=0.5)
            ax[i].hist(getInfo(dataset, 11, variable),histtype='step', bins=nbins, weights=getInfo(dataset, 11, "Weight"), log=log, range=xrange,color='k', alpha=0.5)
            ax[i].hist(getInfo(dataset, -11, variable),histtype='step', bins=nbins, weights=getInfo(dataset, -11, "Weight"), log=log, range=xrange,color='y', alpha=0.5)
            #ax[i].hist(getInfo(dataset, listChargedHadrons, variable),histtype='step', bins=nbins, weights=getInfo(dataset, listChargedHadrons, "Weight"), log=log, range=xrange, label="Ch. Had")
            ax[i].hist(getInfo(dataset, 2112, variable),histtype='step', bins=nbins, weights=getInfo(dataset, 2112, "Weight"), log=log, range=xrange,color='blue', alpha=0.5)
            #ax[i].hist(getInfo(dataset, [-13,13], variable),histtype='step', bins=nbins, weights=getInfo(dataset, [13,-13], "Weight"), log=log, range=xrange, label="Mu+Mu-")

            ax[len(datasetList)].hist(getInfo(dataset, 22, variable),histtype='step', bins=nbins, weights=getInfo(dataset, 22, "Weight"), log=log, range=xrange,color='r', linestyle=lineStyles[i])
            ax[len(datasetList)].hist(getInfo(dataset, 11, variable),histtype='step', bins=nbins, weights=getInfo(dataset, 11, "Weight"), log=log, range=xrange,color='k', linestyle=lineStyles[i])
            ax[len(datasetList)].hist(getInfo(dataset, -11, variable),histtype='step', bins=nbins, weights=getInfo(dataset, -11, "Weight"), log=log, range=xrange,color='y', linestyle=lineStyles[i])
            #ax[len(datasetList)].hist(getInfo(dataset, listChargedHadrons, variable),histtype='step', bins=nbins, weights=getInfo(dataset, listChargedHadrons, "Weight"), log=log, range=xrange, label=labelList[i]+" Ch. Had")
            ax[len(datasetList)].hist(getInfo(dataset, 2112, variable),histtype='step', bins=nbins, weights=getInfo(dataset, 2112, "Weight"), log=log, range=xrange,color='blue', linestyle=lineStyles[i])
            #ax[len(datasetList)].hist(getInfo(dataset, [-13,13], variable),histtype='step', bins=nbins, weights=getInfo(dataset, [13,-13], "Weight"), log=log, range=xrange, label=labelList[i]+" Mu+Mu-")
            
            
            dataset=dataset[(dataset["Time"]>trange[0]) & (dataset["Time"]<trange[1])]
            
            temp=ax[i].hist(getInfo(dataset, 22, variable),histtype='step', bins=nbins, weights=getInfo(dataset, 22, "Weight"), log=log, range=xrange,color='r', label="$\gamma$")
            ax[i].hist(getInfo(dataset, 11, variable),histtype='step', bins=nbins, weights=getInfo(dataset, 11, "Weight"), log=log, range=xrange,color='k', label="e-")
            ax[i].hist(getInfo(dataset, -11, variable),histtype='step', bins=nbins, weights=getInfo(dataset, -11, "Weight"), log=log, range=xrange,color='y', label="e+")
            #ax[i].hist(getInfo(dataset, listChargedHadrons, variable),histtype='step', bins=nbins, weights=getInfo(dataset, listChargedHadrons, "Weight"), log=log, range=xrange, label="Ch. Had")
            ax[i].hist(getInfo(dataset, 2112, variable),histtype='step', bins=nbins, weights=getInfo(dataset, 2112, "Weight"), log=log, range=xrange,color='blue', label="n")
            #ax[i].hist(getInfo(dataset, [-13,13], variable),histtype='step', bins=nbins, weights=getInfo(dataset, [13,-13], "Weight"), log=log, range=xrange, label="Mu+Mu-")

            ax[len(datasetList)].hist(getInfo(dataset, 22, variable),histtype='step', bins=nbins, weights=getInfo(dataset, 22, "Weight"), log=log, range=xrange,color='r', linestyle=lineStyles[i], label=labelList[i]+" $\gamma$")
            ax[len(datasetList)].hist(getInfo(dataset, 11, variable),histtype='step', bins=nbins, weights=getInfo(dataset, 11, "Weight"), log=log, range=xrange,color='k', linestyle=lineStyles[i], label=labelList[i]+" e-")
            ax[len(datasetList)].hist(getInfo(dataset, -11, variable),histtype='step', bins=nbins, weights=getInfo(dataset, -11, "Weight"), log=log, range=xrange,color='y', linestyle=lineStyles[i], label=labelList[i]+" e-")
            #ax[len(datasetList)].hist(getInfo(dataset, listChargedHadrons, variable),histtype='step', bins=nbins, weights=getInfo(dataset, listChargedHadrons, "Weight"), log=log, range=xrange, label=labelList[i]+" Ch. Had")
            ax[len(datasetList)].hist(getInfo(dataset, 2112, variable),histtype='step', bins=nbins, weights=getInfo(dataset, 2112, "Weight"), log=log, range=xrange,color='blue', linestyle=lineStyles[i], label=labelList[i]+" n")
            #ax[len(datasetList)].hist(getInfo(dataset, [-13,13], variable),histtype='step', bins=nbins, weights=getInfo(dataset, [13,-13], "Weight"), log=log, range=xrange, label=labelList[i]+" Mu+Mu-")

#        temp=ax[i].hist(getInfo(dataset, 22, variable),histtype='step', bins=nbins, weights=getInfo(dataset, 22, "Weight"), log=log, range=xrange, label="$\gamma$"+str.format(r' N={:.2e} $\bar x$={:.2e}',getParticleNumber(dataset,22),getInfo(dataset,22,variable).mean()))
#        ax[i].hist(getInfo(dataset, [-11,11], variable),histtype='step', bins=nbins, weights=getInfo(dataset, [11,-11], "Weight"), log=log, range=xrange, label="e+e-"+str.format(r' N={:.2e} $\bar x$={:.2e}',getParticleNumber(dataset,[11,-11]),getInfo(dataset,[11,-11],variable).mean()))
#        ax[i].hist(getInfo(dataset, listChargedHadrons, variable),histtype='step', bins=nbins, weights=getInfo(dataset, listChargedHadrons, "Weight"), log=log, range=xrange, label="Ch. Had"+str.format(r' N={:.2e} $\bar x$={:.2e}',getParticleNumber(dataset,listChargedHadrons),getInfo(dataset,listChargedHadrons,variable).mean()))
#        ax[i].hist(getInfo(dataset, 2112, variable),histtype='step', bins=nbins, weights=getInfo(dataset, 2112, "Weight"), log=log, range=xrange, label="Neutrons"+str.format(r' N={:.2e} $\bar x$={:.2e}',getParticleNumber(dataset,2112),getInfo(dataset,2112,variable).mean()))
#        ax[i].hist(getInfo(dataset, [-13,13], variable),histtype='step', bins=nbins, weights=getInfo(dataset, [13,-13], "Weight"), log=log, range=xrange, label="Mu+Mu-"+str.format(r' N={:.2e} $\bar x$={:.2e}',getParticleNumber(dataset,[13,-13]),getInfo(dataset,[13,-13],variable).mean()))
        else:
            if trange:
                dataset=dataset[(dataset["Time"]>trange[0]) & (dataset["Time"]<trange[1])]
            temp=ax[i].hist(getInfo(dataset, 22, variable),histtype='step', bins=nbins, weights=getInfo(dataset, 22, "Weight"), log=log, range=xrange,color='r', label="$\gamma$")
            ax[i].hist(getInfo(dataset, 11, variable),histtype='step', bins=nbins, weights=getInfo(dataset, 11, "Weight"), log=log, range=xrange,color='k', label="e-")
            ax[i].hist(getInfo(dataset, -11, variable),histtype='step', bins=nbins, weights=getInfo(dataset, -11, "Weight"), log=log, range=xrange,color='y', label="e+")
            #ax[i].hist(getInfo(dataset, listChargedHadrons, variable),histtype='step', bins=nbins, weights=getInfo(dataset, listChargedHadrons, "Weight"), log=log, range=xrange, label="Ch. Had")
            ax[i].hist(getInfo(dataset, 2112, variable),histtype='step', bins=nbins, weights=getInfo(dataset, 2112, "Weight"), log=log, range=xrange,color='blue', label="n")
            #ax[i].hist(getInfo(dataset, [-13,13], variable),histtype='step', bins=nbins, weights=getInfo(dataset, [13,-13], "Weight"), log=log, range=xrange, label="Mu+Mu-")



            ax[len(datasetList)].hist(getInfo(dataset, 22, variable),histtype='step', bins=nbins, weights=getInfo(dataset, 22, "Weight"), log=log, range=xrange,color='r', linestyle=lineStyles[i], label=labelList[i]+" $\gamma$")
            ax[len(datasetList)].hist(getInfo(dataset, 11, variable),histtype='step', bins=nbins, weights=getInfo(dataset, 11, "Weight"), log=log, range=xrange,color='k', linestyle=lineStyles[i], label=labelList[i]+" e-")
            ax[len(datasetList)].hist(getInfo(dataset, -11, variable),histtype='step', bins=nbins, weights=getInfo(dataset, -11, "Weight"), log=log, range=xrange,color='y', linestyle=lineStyles[i], label=labelList[i]+" e-")
            #ax[len(datasetList)].hist(getInfo(dataset, listChargedHadrons, variable),histtype='step', bins=nbins, weights=getInfo(dataset, listChargedHadrons, "Weight"), log=log, range=xrange, label=labelList[i]+" Ch. Had")
            ax[len(datasetList)].hist(getInfo(dataset, 2112, variable),histtype='step', bins=nbins, weights=getInfo(dataset, 2112, "Weight"), log=log, range=xrange,color='blue', linestyle=lineStyles[i], label=labelList[i]+" n")
            #ax[len(datasetList)].hist(getInfo(dataset, [-13,13], variable),histtype='step', bins=nbins, weights=getInfo(dataset, [13,-13], "Weight"), log=log, range=xrange, label=labelList[i]+" Mu+Mu-")


        if i==0:
            maxHeight=temp[0].max() #Calcolo l'altezza massima del bin in modo da forzare gli N grafici ad avere la stessa scala verticale facilitando il confronto
            #print(maxHeight)
        ax[i].set_title(labelList[i])
        #box = ax[i].get_position()
        #ax[i].set_position([box.x0, box.y0 , box.width, box.height * 0.8])
        #ax[i].legend(loc='upper center', bbox_to_anchor=(0.5, 1.2), ncol=3)
        ax[i].legend()
        if ymax!=None:
            ax[i].axis(ymin=1e2, ymax=ymax)
        else:
            ax[i].axis(ymin=1e2, ymax=maxHeight*2)

        if ymax!=None:
            ax[len(datasetList)].axis(ymin=1e2, ymax=ymax)
        else:
            ax[len(datasetList)].axis(ymin=1e2, ymax=maxHeight*2)
        
    ax[len(datasetList)].set_title("Comparison")
    #ax[len(datasetList)].legend(fontsize="x-small")
    box = ax[len(datasetList)].get_position()
   # ax[len(datasetList)].set_position([box.x0, box.y0 , box.width, box.height * 0.8])
   # ax[len(datasetList)].legend(loc='upper center', bbox_to_anchor=(0.5, 1.25), ncol=2, fontsize="x-small")
    ax[len(datasetList)].set_xlabel(xlabel,fontsize='14')
    ax[len(datasetList)].set_ylabel(ylabel,fontsize='14')
    ax[len(datasetList)].legend(fontsize='9')
    fig.subplots_adjust(left = 0.05,right = 0.99,wspace = 0.5, hspace = 0.5, top=0.85, bottom= 0.2)

    fig.tight_layout()
    figname=runName+figTitle
    pl.savefig(figname,transparent=False, facecolor='white')

def plotMomenta(datasetList, particleList, particleLabel, title, xlabel="p [GeV/c]", ylabel="Arb. Units", nbins=nbins, log=True, figName="", xrange=None, ymax=None):
    fig, ax = plt.subplots(nrows=1, ncols=len(datasetList)+1, figsize=((len(datasetList)+1)*8,8), sharey=False)
    plt.suptitle(title)

    for i, dataset in enumerate(datasetList):
        ax[i].set_xlabel(xlabel,fontsize='14')
        ax[i].set_ylabel(ylabel,fontsize='14')

        for iPart, particle in enumerate(particleList):
            temp=ax[i].hist(getMomentum(dataset, particle),histtype='step', bins=nbins, 
                   weights=getInfo(dataset, particle, "Weight"), log=log, range=xrange, 
                   label=particleLabel[iPart]+str.format(r' N={:.2e} $\bar x$={:.2e}',getParticleNumber(dataset,particle),getInfo(dataset,particle,"KinE").mean()))
            
            ax[len(datasetList)].hist(getMomentum(dataset, particle),histtype='step', bins=nbins, weights=getInfo(dataset, particle, "Weight"), log=log, range=xrange, label=labelList[i]+" "+particleLabel[iPart])
        
            if iPart==0 and i==0:
                maxHeight=temp[0].max() #Calcolo l'altezza massima del bin in modo da forzare gli N grafici ad avere la stessa scala verticale facilitando il confronto
        print("{:.2e}".format(maxHeight))

        if ymax!=None:
            ax[i].axis(ymin=1e1, ymax=ymax)
        else:
            ax[i].axis(ymin=1e1, ymax=maxHeight*2)
            
        ax[i].set_title(labelList[i])
        ax[i].legend()
    
    if ymax!=None:
        ax[len(datasetList)].axis(ymin=1e1, ymax=ymax)
    else:
        ax[len(datasetList)].axis(ymin=1e1, ymax=maxHeight*2)
            
    ax[len(datasetList)].set_title("Comparison")
    ax[len(datasetList)].legend()
    ax[len(datasetList)].set_xlabel("p [GeV/c]",fontsize='14')
    ax[len(datasetList)].set_ylabel("Arb. Units",fontsize='14')

    figname=runName+figName
    pl.savefig(figname,transparent=False, facecolor='white')
    
def plotAllEnergySpectra(datasetList, nbins=nbins, logY=True, logX=False):
    fig, axs = plt.subplots(nrows=1, ncols=5, figsize=(22,5), sharey=False)
    plt.suptitle("Energy Spectra")

    for i, dataset in enumerate(datasetList):
        axs[0].hist(getInfo(dataset,22,"KinE"),histtype='step', bins=nbins, weights=getInfo(dataset, 22, "Weight"), range=[0,np.percentile(getInfo(datasetList[0],22,"KinE").to_numpy(),99.999)], log=logY, label=labelList[i]+str.format(r' N={:.2e} $\bar x$={:.2e}',getParticleNumber(dataset,22),getInfo(dataset,22,"KinE").mean()))
        axs[1].hist(getInfo(dataset,2112,"KinE"),histtype='step', bins=nbins, weights=getInfo(dataset, 2112, "Weight"), range=[0,np.percentile(getInfo(datasetList[0],2112,"KinE").to_numpy(),99.999)], log=logY, label=labelList[i]+str.format(r' N={:.2e} $\bar x$={:.2e}',getParticleNumber(dataset,2112),getInfo(dataset,2112,"KinE").mean()))
        axs[2].hist(getInfo(dataset,[11,-11],"KinE"),histtype='step', bins=nbins, weights=getInfo(dataset, [11,-11], "Weight"), range=[0,np.percentile(getInfo(datasetList[0],[11,-11],"KinE").to_numpy(),99.999)], log=logY, label=labelList[i]+str.format(r' N={:.2e} $\bar x$={:.2e}',getParticleNumber(dataset,[11,-11]),getInfo(dataset,[11,11],"KinE").mean()))
        axs[3].hist(getInfo(dataset,listChargedHadrons,"KinE"),histtype='step', bins=nbins, weights=getInfo(dataset, listChargedHadrons, "Weight"), range=[0,np.percentile(getInfo(datasetList[0],listChargedHadrons,"KinE").to_numpy(),99.999)], log=logY, label=labelList[i]+str.format(r' N={:.2e} $\bar x$={:.2e}',getParticleNumber(dataset,listChargedHadrons),getInfo(dataset,listChargedHadrons,"KinE").mean()))
        axs[4].hist(getInfo(dataset,[13,-13],"KinE"),histtype='step', bins=nbins, weights=getInfo(dataset, [13,-13], "Weight"), range=[0,np.percentile(getInfo(datasetList[0],[13,-13],"KinE").to_numpy(),99.999)], log=logY, label=labelList[i]+str.format(r' N={:.2e} $\bar x$={:.2e}',getParticleNumber(dataset,[13,-13]),getInfo(dataset,[13,-13],"KinE").mean()))

    for i in range(0,5):
        axs[i].set_ylim(auto=True)
        if not logX:
            axs[i].legend(loc= "best")
        axs[i].grid(True)
        axs[i].set_xlabel("$E_{kin}$ [GeV]",fontsize='14')
        axs[i].set_ylabel("Arb. Units",fontsize='14')
        if logX:
            axs[i].set_xscale("log")
        
    axs[0].set_title("$\gamma$")
    axs[1].set_title("n")
    axs[2].set_title("$e^-~/~e^+$")
    axs[3].set_title("ch. had")
    axs[4].set_title("$\mu^-~/~\mu^+$")

    fig.subplots_adjust(left = 0.05,right = 0.99,wspace = 0.5, hspace = 0.5, top=0.85, bottom= 0.2)

    figname=runName+"EnergySpectra"
    if logX:
        figname=figname+"logX"
    pl.savefig(figname,transparent=False, facecolor='white')
    
   
def plotLethargy(datasetList, nbins=nbins, logY=True, logX=False, yrange=None, trange=None, plotTitle="Energy Spectra", xrange=None):
    lineStyles=['solid', 'dotted', 'dashed', 'dashdot' ]
    fig, axs = plt.subplots(nrows=1, ncols=len(datasetList)+1, figsize=((len(datasetList)+1)*5,5), sharey=False)
#    plt.suptitle("Lethergy plots")
    if trange:
        plt.suptitle(plotTitle+str.format(' tmin={} [ns] tmax={} [ns]', trange[0],trange[1]))
    else:
        plt.suptitle(plotTitle)

    for i, dataset in enumerate(datasetList):
        if trange:
            dataset=dataset[(dataset["Time"]>trange[0]) & (dataset["Time"]<trange[1])]
        axs[i].hist(np.log10(getInfo(dataset,22,"KinE")),histtype='step', bins=nbins,range=xrange, weights=getInfo(dataset, 22, "Weight"), log=logY,color='r',linestyle=lineStyles[i], label='$\gamma$')
        axs[i].hist(np.log10(getInfo(dataset,11,"KinE")),histtype='step', bins=nbins,range=xrange, weights=getInfo(dataset, 11, "Weight"), log=logY,color='k',linestyle=lineStyles[i], label='$e^-$')
        axs[i].hist(np.log10(getInfo(dataset,-11,"KinE")),histtype='step', bins=nbins,range=xrange, weights=getInfo(dataset, -11, "Weight"), log=logY,color='y',linestyle=lineStyles[i], label='$e^+$')
        axs[i].hist(np.log10(getInfo(dataset,2112,"KinE")),histtype='step', bins=nbins,range=xrange, weights=getInfo(dataset, 2112, "Weight"), log=logY,color='b',linestyle=lineStyles[i], label="n")

        axs[i].set_title(labelList[i])
        axs[i].legend(loc= 'upper left')
        
        axs[len(datasetList)].hist(np.log10(getInfo(dataset,22,"KinE")),histtype='step', bins=nbins,range=xrange, weights=getInfo(dataset, 22, "Weight"), log=logY,color='r',linestyle=lineStyles[i], label='$\gamma$')
        axs[len(datasetList)].hist(np.log10(getInfo(dataset,11,"KinE")),histtype='step', bins=nbins,range=xrange, weights=getInfo(dataset, 11, "Weight"), log=logY,color='k',linestyle=lineStyles[i], label='$e^-$')
        axs[len(datasetList)].hist(np.log10(getInfo(dataset,-11,"KinE")),histtype='step', bins=nbins,range=xrange, weights=getInfo(dataset, -11, "Weight"), log=logY,color='y',linestyle=lineStyles[i], label='$e^+$')
        axs[len(datasetList)].hist(np.log10(getInfo(dataset,2112,"KinE")),histtype='step', bins=nbins,range=xrange, weights=getInfo(dataset, 2112, "Weight"), log=logY,color='b',linestyle=lineStyles[i], label="n")
        axs[len(datasetList)].set_title("Comparison")
        
        if yrange:
            axs[i].set_ylim(yrange)
            axs[len(datasetList)].set_ylim(yrange)
        if logX:
            axs[i].set_xscale("log")
        #axs[i].grid(True)
        axs[i].set_xlabel("$Log(E_{kin})$ [GeV]",fontsize=14)
        axs[i].set_ylabel('dN/dlog(E$_{kin}$)',fontsize=14)

    #axs[len(datasetList)].grid(True)
    axs[len(datasetList)].set_xlabel("$Log(E_{kin})$ [GeV]",fontsize=14)
    axs[len(datasetList)].set_ylabel('dN/dlog(E$_{kin}$)',fontsize=14)
    fig.subplots_adjust(left = 0.05,right = 0.99,wspace = 0.5, hspace = 0., top=0.85, bottom= 0.2)
    figname=runName+"Lethergy"
    if logX:
        figname=figname+"logX"
    pl.savefig(figname,transparent=False, facecolor='white')
    
def plotLethargyAlsoWithTimeCut(datasetList, nbins=nbins, logY=True, logX=False, yrange=None, trange=None, plotTitle="Energy Spectra", xrange=None):
    lineStyles=['solid', 'dotted', 'dashed', 'dashdot' ]
    fig, axs = plt.subplots(nrows=1, ncols=len(datasetList), figsize=((len(datasetList))*5,5), sharey=False)
#    plt.suptitle("Lethergy plots")

    plt.suptitle(plotTitle+str.format(' w/ and w/o Time Cut tmin={} [ns] tmax={} [ns]', trange[0],trange[1]))
    for i, dataset in enumerate(datasetList):
        axs[i].hist(np.log10(getInfo(dataset,22,"KinE")),histtype='step', bins=nbins,range=xrange, weights=getInfo(dataset, 22, "Weight"), log=logY,color='r',linestyle=lineStyles[0], label='$\gamma$', alpha=0.5)
        axs[i].hist(np.log10(getInfo(dataset,11,"KinE")),histtype='step', bins=nbins,range=xrange, weights=getInfo(dataset, 11, "Weight"), log=logY,color='k',linestyle=lineStyles[0], label='$e^-$', alpha=0.5)
        axs[i].hist(np.log10(getInfo(dataset,-11,"KinE")),histtype='step', bins=nbins,range=xrange, weights=getInfo(dataset, -11, "Weight"), log=logY,color='y',linestyle=lineStyles[0], label='$e^+$', alpha=0.5)
 #       axs[i].hist(np.log10(getInfo(dataset,13,"KinE")),histtype='step', bins=nbins,range=xrange, weights=getInfo(dataset, 13, "Weight"), log=logY,color='g',linestyle=lineStyles[0], label='$\mu^-$')
 #       axs[i].hist(np.log10(getInfo(dataset,-13,"KinE")),histtype='step', bins=nbins,range=xrange, weights=getInfo(dataset, -13, "Weight"), log=logY,color='m',linestyle=lineStyles[0], label='$\mu^+$')
        axs[i].hist(np.log10(getInfo(dataset,2112,"KinE")),histtype='step', bins=nbins,range=xrange, weights=getInfo(dataset, 2112, "Weight"), log=logY,color='b',linestyle=lineStyles[0], label="n", alpha=0.5)
        if trange:
            dataset=dataset[(dataset["Time"]>trange[0]) & (dataset["Time"]<trange[1])]
        axs[i].hist(np.log10(getInfo(dataset,22,"KinE")),histtype='step', bins=nbins,range=xrange, weights=getInfo(dataset, 22, "Weight"), log=logY,color='r',linestyle=lineStyles[0])
        axs[i].hist(np.log10(getInfo(dataset,11,"KinE")),histtype='step', bins=nbins,range=xrange, weights=getInfo(dataset, 11, "Weight"), log=logY,color='k',linestyle=lineStyles[0])
        axs[i].hist(np.log10(getInfo(dataset,-11,"KinE")),histtype='step', bins=nbins,range=xrange, weights=getInfo(dataset, -11, "Weight"), log=logY,color='y',linestyle=lineStyles[0])
 #       axs[i].hist(np.log10(getInfo(dataset,13,"KinE")),histtype='step', bins=nbins,range=xrange, weights=getInfo(dataset, 13, "Weight"), log=logY,color='g',linestyle=lineStyles[0])
 #       axs[i].hist(np.log10(getInfo(dataset,-13,"KinE")),histtype='step', bins=nbins,range=xrange, weights=getInfo(dataset, -13, "Weight"), log=logY,color='m',linestyle=lineStyles[0]) 
        axs[i].hist(np.log10(getInfo(dataset,2112,"KinE")),histtype='step', bins=nbins,range=xrange, weights=getInfo(dataset, 2112, "Weight"), log=logY,color='b',linestyle=lineStyles[0])

        axs[i].set_title(labelList[i])
        axs[i].legend(loc= 'upper left')
        
        if yrange:
            axs[i].set_ylim(yrange)
        if logX:
            axs[i].set_xscale("log")
        #axs[i].grid(True)
        axs[i].set_xlabel("$Log(E_{kin})$ [GeV]",fontsize=14)
        axs[i].set_ylabel('dN/dLog(E$_{kin}$)',fontsize=14)

    #axs[len(datasetList)].grid(True)
    fig.subplots_adjust(left = 0.05,right = 0.99,wspace = 0.5, hspace = 0., top=0.85, bottom= 0.2)
    figname=runName+"Lethergy"
    if logX:
        figname=figname+"logX"
    pl.savefig(figname,transparent=False, facecolor='white')
    

left, width = 0.1, 0.65
bottom, height = 0.1, 0.65
spacing = 0.005
rect_scatter = [left, bottom, width, height]
rect_histx = [left, bottom + height + spacing, width, 0.2]
rect_histy = [left + width + spacing, bottom, 0.2, height]

def plotSingleDistributionLongMomNozzle(datasetList, plotTitle="", xlabel="", ylabel="Arb. Units", nbins=nbins, log=True, figTitle="", xrange=None, ymax=None, secondaryFlag=True):
    fig, ax = plt.subplots(nrows=1, ncols=len(datasetList)+1, figsize=((len(datasetList)+1)*8,8), sharey=False)
    plt.suptitle(plotTitle)
    for i, dataset in enumerate(datasetList):
        if secondaryFlag:
            dataset=dataset[dataset["NumPart"]>0]
        else:
            dataset=dataset[dataset["NumPart"]==0]
        ax[i].set_xlabel(xlabel,fontsize='14')
        ax[i].set_ylabel(ylabel,fontsize='14')

        ax[i].hist(-dataset["CZ"][dataset["PosZEle"]<600]*dataset["EneEle"][dataset["PosZEle"]<600],histtype='step', bins=nbins, weights=dataset["Weight"][dataset["PosZEle"]<600], log=log, range=xrange, label=labelList[i])
        ax[i].set_title(labelList[i])
        ax[i].legend()
        ax[len(datasetList)].hist(-dataset["CZ"][dataset["PosZEle"]<600]*dataset["EneEle"][dataset["PosZEle"]<600],histtype='step', bins=nbins, weights=dataset["Weight"][dataset["PosZEle"]<600], log=log, range=xrange, label=labelList[i])
     
    ax[len(datasetList)].set_title("Comparison")
    ax[len(datasetList)].legend()
    ax[len(datasetList)].set_xlabel(xlabel,fontsize='14')
    ax[len(datasetList)].set_ylabel(ylabel,fontsize='14')
    figname=runName+figTitle
    pl.savefig(figname,transparent=False, facecolor='white')
    
def plotSingleDistributionTransMomNozzle(datasetList, plotTitle="", xlabel="", ylabel="Arb. Units", nbins=nbins, log=True, figTitle="", xrange=None, ymax=None, secondaryFlag=True):
    fig, ax = plt.subplots(nrows=1, ncols=len(datasetList)+1, figsize=((len(datasetList)+1)*8,8), sharey=False)
    plt.suptitle(plotTitle)
    for i, dataset in enumerate(datasetList):
        if secondaryFlag:
            dataset=dataset[dataset["NumPart"]>0]
        else:
            dataset=dataset[dataset["NumPart"]==0]
        ax[i].set_xlabel(xlabel,fontsize='14')
        ax[i].set_ylabel(ylabel,fontsize='14')

        ax[i].hist(np.sqrt( dataset["CX"][dataset["PosZEle"]<600]**2 + dataset["CY"][dataset["PosZEle"]<600]**2)*dataset["EneEle"][dataset["PosZEle"]<600]
                        ,histtype='step', bins=nbins, weights=dataset["Weight"][dataset["PosZEle"]<600], log=log, range=xrange, label=labelList[i])
        ax[i].set_title(labelList[i])
        ax[i].legend()
        ax[len(datasetList)].hist(np.sqrt( dataset["CX"][dataset["PosZEle"]<600]**2 + dataset["CY"][dataset["PosZEle"]<600]**2)*dataset["EneEle"][dataset["PosZEle"]<600],histtype='step', bins=nbins, weights=dataset["Weight"][dataset["PosZEle"]<600], log=log, range=xrange, label=labelList[i])
     
    ax[len(datasetList)].set_title("Comparison")
    ax[len(datasetList)].legend()
    ax[len(datasetList)].set_xlabel(xlabel,fontsize='14')
    ax[len(datasetList)].set_ylabel(ylabel,fontsize='14')
    figname=runName+figTitle
    pl.savefig(figname,transparent=False, facecolor='white')
    
        
def plotSingleDistributionBendingRadius(datasetList, plotTitle="", xlabel="", ylabel="Arb. Units", nbins=nbins, log=True, figTitle="", xrange=None, ymax=None, secondaryFlag=True):
    fig, ax = plt.subplots(nrows=1, ncols=len(datasetList)+1, figsize=((len(datasetList)+1)*8,8), sharey=False)
    plt.suptitle(plotTitle)
    qB=1.071e9
    for i, dataset in enumerate(datasetList):
        if secondaryFlag:
            dataset=dataset[dataset["NumPart"]>0]
        else:
            dataset=dataset[dataset["NumPart"]==0]
        ax[i].set_xlabel(xlabel,fontsize='14')
        ax[i].set_ylabel(ylabel,fontsize='14')

        ax[i].hist(np.sqrt( dataset["CX"][dataset["PosZEle"]<600]**2 + dataset["CY"][dataset["PosZEle"]<600]**2)*dataset["EneEle"][dataset["PosZEle"]<600]*1e9/qB,histtype='step', bins=nbins, weights=dataset["Weight"][dataset["PosZEle"]<600], log=log, range=xrange, label=labelList[i])


        ax[i].set_title(labelList[i])
        ax[i].legend()
        ax[len(datasetList)].hist(np.sqrt( dataset["CX"][dataset["PosZEle"]<600]**2 + dataset["CY"][dataset["PosZEle"]<600]**2)*dataset["EneEle"][dataset["PosZEle"]<600]*1e9/qB,histtype='step', bins=nbins, weights=dataset["Weight"][dataset["PosZEle"]<600], log=log, range=xrange, label=labelList[i])
     
    ax[len(datasetList)].set_title("Comparison")
    ax[len(datasetList)].legend()
    ax[len(datasetList)].set_xlabel(xlabel,fontsize='14')
    ax[len(datasetList)].set_ylabel(ylabel,fontsize='14')
    figname=runName+figTitle
    pl.savefig(figname,transparent=False, facecolor='white')
    
def plotSingleVariable2D(datasetList, variableX, variableY, plotTitle="", xlabel="", ylabel="Arb. Units", nbins=nbins, log=True, figTitle="", range=None,vmin=None,vmax=None, hasBIB=-1):
    ## This function plots two given variables, with no particle selection (it is therefore mostly used for datasetEle, where we have ony electrons)
    ## "hasBIB" flag allows to select particles that generated BIB (1), not generated BIB (0), and all particles (-1, default)
    ## A plot for each dataset is produced, plus one last plot with all datasets superimposed
    fig, ax = plt.subplots(nrows=1, ncols=len(datasetList)+1, figsize=((len(datasetList)+1)*9,8), sharey=False)
    plt.suptitle(plotTitle)
    
    for i, dataset in enumerate(datasetList):
        if hasBIB==1:
            dataset=dataset[dataset["NumPart"]>0]
        elif hasBIB==0:
            dataset=dataset[dataset["NumPart"]==0]
        ax[i].set_xlabel(xlabel,fontsize='14')
        ax[i].set_ylabel(ylabel,fontsize='14')
       # ax[i].hist(dataset[variable],histtype='step', bins=nbins, weights=dataset["Weight"], log=log, range=xrange)
        ax[i].hist2d(dataset[variableX], dataset[variableY], range=range,weights=dataset["Weight"], bins=nbins, norm=matplotlib.colors.LogNorm(vmin=vmin,vmax=vmax), cmap='Reds')
        PCM=ax[i].get_children()[0] #get the mappable, the 1st and the 2nd are the x and y axes
        cb=plt.colorbar(PCM, ax=ax[i]) 
        cb.set_label('Arb. Units')
        ax[i].set_title(labelList[i])
#        ax[len(datasetList)].hist(dataset[variable],histtype='step', bins=nbins, weights=dataset["Weight"], log=log, range=xrange, label=labelList[i])
        ax[len(datasetList)].hist2d(dataset[variableX] ,dataset[variableY], range=range,weights=dataset["Weight"], bins=nbins, norm=matplotlib.colors.LogNorm())
      #  PCMb=ax[len(datasetList)].get_children()[0] #get the mappable, the 1st and the 2nd are the x and y axes
      #  plt.colorbar(PCMb, ax=ax[len(datasetList)]) 

    ax[len(datasetList)].set_title("Comparison")
   # ax[len(datasetList)].legend()
    ax[len(datasetList)].set_xlabel(xlabel,fontsize='14')
    ax[len(datasetList)].set_ylabel(ylabel,fontsize='14')

    figname=runName+figTitle
    pl.savefig(figname,transparent=False, facecolor='white')   
    
def plotStackElePlotsBibNoBib(datasetList, variable, title, xlabel, figname, log=False, nbins=100):
    ## This function plots a given variable of the electronDataset, superimposing the distributions of BIB and NoBIB (and the total)
    ## A pair of plots for each dataset is produced, with the bottom one with the histogram stacked (without the total one)
    fig, ax = plt.subplots(nrows=2, ncols=len(datasetList), figsize=(16,12), sharex=False)
    fig.suptitle(title)
    for i, dataset in enumerate(datasetEleList):
        ax[0][i].set_title(labelList[i])
        ax[0][i].hist(dataset[variable], bins=nbins, label=["All"], histtype="step", log=log)
        ax[0][i].hist(dataset[dataset["NumPart"]==0][variable], bins=nbins, label=["NoBib"], histtype="step", log=log)
        ax[0][i].hist(dataset[dataset["NumPart"]>0][variable], bins=nbins, label=["Bib"], histtype="step", log=log)
        ax[0][i].set_xlabel(xlabel)
        ax[0][i].legend()

        ax[1][i].set_title("Stacked "+labelList[i])
        ax[1][i].hist([dataset[dataset["NumPart"]==0][variable],dataset[dataset["NumPart"]>0][variable]], bins=nbins, stacked=True, label=["NoBib","BiB"], log=log)

        ax[1][i].set_xlabel(xlabel)
        ax[1][i].legend()
    fig.tight_layout()
    plt.savefig(runName+runName+figname,transparent=False, facecolor='white')


        
def plotEleDistrWithCut(dataset, variable, cutVariable, cutCenter, cutRange, nbins=50, title="title", xlabel="x", figname="trash", log=False):
    plt.figure(figsize=(12,6))
    plt.title(variable + " cut on "+cutVariable +" around "+ str(cutCenter))
    plt.hist(dataset[variable], bins=nbins, histtype="step", log=log, label="NoCut")
    for cut in cutRange:
        plt.hist(dataset[abs(dataset[cutVariable]-cutCenter)<cut][variable], bins=nbins, histtype="step", log=log, label="Cut_"+str(cut))
    fig.tight_layout()
    plt.xlabel(xlabel)
    plt.legend()
    plt.savefig(runName+figname,transparent=False, facecolor='white')

    

def plotSingleVariable(datasetList, variable, plotTitle="", xlabel="", ylabel="Arb. Units", nbins=nbins, log=True, figTitle="", xrange=None, ymax=None, hasBIB=-1):
    ## This function plots a given variable, with no particle selection (it is therefore mostly used for datasetEle, where we have ony electrons)
    ## "hasBIB" flag allows to select particles that generated BIB (1), not generated BIB (0), and all particles (-1, default)
    ## A plot for each dataset is produced, plus one last plot with all datasets superimposed
    fig, ax = plt.subplots(nrows=1, ncols=len(datasetList)+1, figsize=((len(datasetList)+1)*8,8), sharey=False)
    plt.suptitle(plotTitle)
    for i, dataset in enumerate(datasetList):
        if hasBIB==1:
            dataset=dataset[dataset["NumPart"]>0]
        elif hasBIB==0:
            dataset=dataset[dataset["NumPart"]==0]
        ax[i].set_xlabel(xlabel,fontsize='14')
        ax[i].set_ylabel(ylabel,fontsize='14')
        ax[i].hist(dataset[variable],histtype='step', bins=nbins, weights=dataset["Weight"], log=log, range=xrange, label=labelList[i]+str.format(r' N={:.2e} $\bar x$={:.2e}',sum(dataset["Weight"]),dataset[variable].mean()))
        ax[i].set_title(labelList[i])
        ax[i].legend()
        ax[len(datasetList)].hist(dataset[variable],histtype='step', bins=nbins, weights=dataset["Weight"], log=log, range=xrange, label=labelList[i])
     
    ax[len(datasetList)].set_title("Comparison")
    ax[len(datasetList)].legend()
    ax[len(datasetList)].set_xlabel(xlabel,fontsize='14')
    ax[len(datasetList)].set_ylabel(ylabel,fontsize='14')

    figname=runName+figTitle
    pl.savefig(figname,transparent=False, facecolor='white')
    


# ## Read Datasets

# In[5]:


for fileNumber, fileName in enumerate(inputFilesList):
    print("Read file: ",fileNumber, fileName)
    temp=pd.read_csv(fileName, header=None, names=colsToRead, delim_whitespace=True)
    datasetList.append(temp)
    if flagReadEle:
        tempEle=pd.read_csv(fileName+"_ele", header=None, names=colsToReadEle, delim_whitespace=True)
        datasetEleList.append(tempEle)


# ------------------------------------------

# ## Let's have a look at muon decay z position

# In[6]:


nBinZ=[]
histoCumA=[]
histoCumAnorm=[]


# ### Z of first interaction of decay electrons (from datasetBIB)
# In this case, each interaction is counted as many time as the number of BIB particles he eventually generated

# In[7]:


fig, axs = plt.subplots(nrows=4, ncols=1, figsize=(10,12), sharex=True)
axs2b=axs[3].twinx() #Same x, but different y scale
binWidthZ=0.1
for datasetNumber, dataset in enumerate(datasetList):
    nBinZ.append(int(((dataset["PosZFI"]/100).max()-(dataset["PosZFI"]/100).min())/binWidthZ))
    plot1D(axs[0],dataset["PosZFI"]/100, weights=dataset["Weight"],bins=nBinZ[-1], plotTitle="PosZEle BIB", label=labelList[datasetNumber], xlabel='', ylabel='Arb. Units' )
    axs[1].hist(dataset["PosZFI"]/100, weights=dataset["Weight"], bins=nBinZ[-1], cumulative=1,histtype='step', label="cum "+labelList[datasetNumber], density=False)
    histoCumAnorm.append(axs[2].hist(dataset["PosZFI"]/100, weights=dataset["Weight"], bins=nBinZ[-1], cumulative=1,histtype='step', label="cum "+labelList[datasetNumber], density=True))
    plot1D(axs[3],dataset["PosZFI"]/100, weights=dataset["Weight"],bins=nBinZ[-1], plotTitle="Ele Decay Z with Cumulatives", label=labelList[datasetNumber], xlabel='$z_{e \,dec}$ (m)', ylabel='Arb. Units' )
    histoCumA=axs2b.hist(dataset["PosZFI"]/100, weights=dataset["Weight"], bins=nBinZ[-1], cumulative=1,histtype='step', label="cum "+labelList[datasetNumber], density=False, linestyle=':')


#axs[0].axis(ymin=100, ymax=1e9)
axs[0].set_xlim(-2,20)
axs[0].grid(True, which="both", axis='y')
axs[0].locator_params(axis="x", nbins=20)
axs[2].grid(True, which="both")
axs[2].set_title("Cumulative Function Norm")
axs[2].legend(loc= "best", fontsize='x-small')

axs[1].legend(loc= "best", fontsize='x-small')
axs[1].set_title("Cumulative Function Not Norm")
axs[1].grid(True, which="both")

figname=runName+"BIBEleDecZCumDaBIB"
pl.savefig(figname,transparent=False, facecolor='white')


# same plot as above but inverse cumulative

# In[8]:


if flagAllPlots:
    fig, axs = plt.subplots(nrows=4, ncols=1, figsize=(10,12), sharex=True)
    axs2b=axs[3].twinx() #Same x, but different y scale
    binWidthZ=0.1
    for datasetNumber, dataset in enumerate(datasetList):
        nBinZ.append(int(((dataset["PosZFI"]/100).max()-(dataset["PosZFI"]/100).min())/binWidthZ))
        plot1D(axs[0],dataset["PosZFI"]/100, weights=dataset["Weight"],bins=nBinZ[-1], plotTitle="PosZEle BIB", label=labelList[datasetNumber], xlabel='', ylabel='Arb. Units' )
        axs[1].hist(dataset["PosZFI"]/100, weights=dataset["Weight"], bins=nBinZ[-1], cumulative=-1,histtype='step', label="cum "+labelList[datasetNumber], density=False)
        histoCumAnorm.append(axs[2].hist(dataset["PosZFI"]/100, weights=dataset["Weight"], bins=nBinZ[-1], cumulative=-1,histtype='step', label="cum "+labelList[datasetNumber], density=True))
        plot1D(axs[3],dataset["PosZFI"]/100, weights=dataset["Weight"],bins=nBinZ[-1], plotTitle="Ele Decay Z with Cumulatives", label=labelList[datasetNumber], xlabel='$z_{e \,dec}$ (m)', ylabel='Arb. Units' )
        histoCumA=axs2b.hist(dataset["PosZFI"]/100, weights=dataset["Weight"], bins=nBinZ[-1], cumulative=1,histtype='step', label="cum "+labelList[datasetNumber], density=False, linestyle=':')


    #axs[0].axis(ymin=100, ymax=1e9)
    axs[0].set_xlim(-2,20)
    axs[0].grid(True, which="both", axis='y')
    axs[0].locator_params(axis="x", nbins=20)
    axs[2].grid(True, which="both")
    axs[2].set_title("Cumulative Function Norm")
    axs[2].legend(loc= "best", fontsize='x-small')

    axs[1].legend(loc= "best", fontsize='x-small')
    axs[1].set_title("Cumulative Function Not Norm")
    axs[1].grid(True, which="both")

    figname=runName+"BIBEleDecZInvCumDaBIB"
    pl.savefig(figname,transparent=False, facecolor='white')


# ### Z of first interaction of decay electrons (from datasetEle)
# In this case, instead, each interaction is counted once, independently from the number of BIB particles he eventually generated

# In[9]:


if flagAllPlots:
    fig, axs = plt.subplots(nrows=4, ncols=1, figsize=(10,12), sharex=True)
    axs2b=axs[3].twinx() #Same x, but different y scale
    binWidthZ=0.1
    for datasetNumber, dataset in enumerate(datasetEleList):
        nBinZ.append(int(((dataset["PosZEle"]/100).max()-(dataset["PosZEle"]/100).min())/binWidthZ))
        plot1D(axs[0],dataset[dataset["NumPart"]>0]["PosZEle"]/100, weights=dataset[dataset["NumPart"]>0]["Weight"],bins=nBinZ[-1], plotTitle="PosZEle BIB", label=labelList[datasetNumber], xlabel='', ylabel='Arb. Units' )
        axs[1].hist(dataset[dataset["NumPart"]>0]["PosZEle"]/100, weights=dataset[dataset["NumPart"]>0]["Weight"], bins=nBinZ[-1], cumulative=1,histtype='step', label="cum "+labelList[datasetNumber], density=False)
        histoCumAnorm.append(axs[2].hist(dataset[dataset["NumPart"]>0]["PosZEle"]/100, weights=dataset[dataset["NumPart"]>0]["Weight"], bins=nBinZ[-1], cumulative=1,histtype='step', label="cum "+labelList[datasetNumber], density=True))
        plot1D(axs[3],dataset[dataset["NumPart"]>0]["PosZEle"]/100, weights=dataset[dataset["NumPart"]>0]["Weight"],bins=nBinZ[-1], plotTitle="PosZEle BIB with Cumulatives", label=labelList[datasetNumber], xlabel='$z_{e }$ (m)', ylabel='Arb. Units' )
        histoCumA=axs2b.hist(dataset[dataset["NumPart"]>0]["PosZEle"]/100, weights=dataset[dataset["NumPart"]>0]["Weight"], bins=nBinZ[-1], cumulative=1,histtype='step', label="cum "+labelList[datasetNumber], density=False, linestyle=':')


    #axs[0].axis(ymin=100, ymax=1e9)
    axs[0].set_xlim(-2,20)
    axs[0].grid(True, which="both", axis='y')
    axs[0].locator_params(axis="x", nbins=20)
    axs[2].grid(True, which="both")
    axs[2].set_title("Cumulative Function Norm")
    axs[2].legend(loc= "best", fontsize='x-small')

    axs[1].legend(loc= "best", fontsize='x-small')
    axs[1].set_title("Cumulative Function Not Norm")
    axs[1].grid(True, which="both")

    figname=runName+"BIBEleDecZCum"
    pl.savefig(figname,transparent=False, facecolor='white')


# ### Z of decay position o muons that eventually generated BIB
# 

# In[10]:


fig, axs = plt.subplots(nrows=4, ncols=1, figsize=(10,12), sharex=True)
axs2b=axs[3].twinx() #Same x, but different y scale
binWidthZ=0.1

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
pl.savefig(figname,transparent=False, facecolor='white')


# ### Z of exiting position of BIB particles
# 

# In[11]:


fig, axs = plt.subplots(nrows=4, ncols=1, figsize=(10,12), sharex=True)
axs2b=axs[3].twinx() #Same x, but different y scale
binWidthZ=0.1
for datasetNumber, dataset in enumerate(datasetList):
    nBinZ.append(int(((dataset["PosZ"]/100).max()-(dataset["PosZ"]/100).min())/binWidthZ))
    plot1D(axs[0],dataset["PosZ"]/100, weights=dataset["Weight"],bins=nBinZ[-1], plotTitle="Pos ZZ", label=labelList[datasetNumber], xlabel='', ylabel='Arb. Units' )
    axs[1].hist(dataset["PosZ"]/100, weights=dataset["Weight"], bins=nBinZ[-1], cumulative=1,histtype='step', label="cum "+labelList[datasetNumber], density=False)
    histoCumAnorm.append(axs[2].hist(dataset["PosZ"]/100, weights=dataset["Weight"], bins=nBinZ[-1], cumulative=1,histtype='step', label="cum "+labelList[datasetNumber], density=True))
    plot1D(axs[3],dataset["PosZ"]/100, weights=dataset["Weight"],bins=nBinZ[-1], plotTitle="BIB Z exit with Cumulatives", label=labelList[datasetNumber], xlabel='$z$ (m)', ylabel='Arb. Units' )
    histoCumA=axs2b.hist(dataset["PosZ"]/100, weights=dataset["Weight"], bins=nBinZ[-1], cumulative=1,histtype='step', label="cum "+labelList[datasetNumber], density=False, linestyle=':')


#axs[0].axis(ymin=100, ymax=1e9)
#axs[0].set_xlim(0,20)
axs[0].grid(True, which="both", axis='y')
axs[0].locator_params(axis="x", nbins=20)
axs[2].grid(True, which="both")
axs[2].set_title("Cumulative Function Norm")
axs[2].legend(loc= "best", fontsize='x-small')

axs[1].legend(loc= "best", fontsize='x-small')
axs[1].set_title("Cumulative Function Not Norm")
axs[1].grid(True, which="both")

figname=runName+"BIBExitX"
pl.savefig(figname,transparent=False, facecolor='white')


# ### If Needed, let's apply a z-cut

# In[12]:


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

# In[13]:


foundParticlesList=[]
foundParticlesFreqList=[]
foundParticlesUnique=set([])
foundParticlesUniqueFreq=[]
foundParticlesUniqueEntries=[]
foundParticlesUniqueEntriesTimeCut=[]


# ### Let's start by creating a list with all and only the particles that appeared in at least a dataset
# 

# In[14]:


for dataset in datasetList:
    (lista, freq)=np.unique(dataset["PDGcode"], return_counts=True)
    foundParticlesList.append(lista)
    foundParticlesFreqList.append(freq)
    foundParticlesUnique=foundParticlesUnique.union(foundParticlesUnique,lista)
    
foundParticlesUnique=list(foundParticlesUnique)
foundParticlesUnique.sort()
print("List of all (and unique) found particles in all datasets\n", foundParticlesUnique)


# ### Retrieve names for these particles

# In[15]:


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
print("Their Names\n", particleNamesList)


# ### For each particle, let's evaluate the frequency with wich it appeared in each dataset
# 

# In[16]:


for i, dataset in enumerate(datasetList):
    foundParticlesUniqueEntries.append([]) # Aggiungo una lista nuova (e vuota) per ciascun dataset
    foundParticlesUniqueEntriesTimeCut.append([])
    for particle in foundParticlesUnique: # Giro su tutte le particelle della lista univoca
        if not dataset[dataset["PDGcode"]==particle].empty: # Se in questo dataset quela particella è comparsa ne calcolo il peso
            foundParticlesUniqueEntries[i].append(sum(dataset[dataset["PDGcode"]==particle]["Weight"]))
            foundParticlesUniqueEntriesTimeCut[i].append(sum(dataset[(dataset["Time"]>timeCut[0]) & (dataset["Time"]<timeCut[1]) & (dataset["PDGcode"]==particle)]["Weight"]))
        else:
            foundParticlesUniqueEntries[i].append(0)
            foundParticlesUniqueEntriesTimeCut[i].append(0)


# In[17]:


temp=foundParticlesUniqueEntries[0]
foundParticlesUnique=[x for _,x in sorted(zip(temp,foundParticlesUnique), reverse=True)]
particleNamesList=[x for _,x in sorted(zip(temp,particleNamesList), reverse=True)]
for i in range(len(datasetList)-1,-1, -1):
    foundParticlesUniqueEntries[i]=[x for _,x in sorted(zip(temp,foundParticlesUniqueEntries[i]), reverse=True)]
    foundParticlesUniqueEntriesTimeCut[i]=[x for _,x in sorted(zip(temp,foundParticlesUniqueEntriesTimeCut[i]), reverse=True)]


# In[18]:


print(foundParticlesUniqueEntries)
print()
print(foundParticlesUniqueEntriesTimeCut)


# ### Let's plot all particles' frequencies

# In[19]:


if flagAllPlots:
    fig, axs = plt.subplots(nrows=len(datasetList)+2, ncols=1, figsize=(18,(len(datasetList)+1)*8))
    fig.suptitle(runName+"Particles Frequencies")
    width=0.3
    for i, dataset in enumerate(datasetList):
        axs[i].bar(np.arange(len(foundParticlesUnique)), foundParticlesUniqueEntries[i], align='center', log=True, yerr=getParticlesNumbersErrors(dataset,nSlicesErrors), ecolor="blue", capsize=10)
    #    axs[i].bar(range(len(foundParticlesUnique)), getParticlesNumbersErrors(dataset,nSlicesErrors), align='center', log=True)

        axs[i].set_title(labelList[i]) 
        plt.sca(axs[i])
        plt.xticks(range(len(foundParticlesUnique)), particleNamesList, size='large')
        plt.xticks(rotation=90)
        plt.ylabel("Occurency", size='large')

        for j, v in enumerate(foundParticlesUniqueEntries[i]):
            if v!=0:
                axs[i].text(j , v, "{:.2e}".format(v))
                #axs[len(datasetList)].text(j , v, "{:.2e}".format(v), color="tab:blue")

        axs[len(datasetList)].bar(np.arange(len(foundParticlesUnique))+i*width, foundParticlesUniqueEntries[i],width=width, yerr=getParticlesNumbersErrors(dataset,nSlicesErrors), ecolor="blue" , align='center', log=True, label=labelList[i], alpha=1, capsize=10)
        axs[len(datasetList)].set_title('Comparison')
        plt.sca(axs[len(datasetList)])
        plt.xticks(np.arange(len(foundParticlesUnique))+i*width/2, particleNamesList, size='large')
        plt.xticks(rotation=90)
        plt.ylabel("Occurency", size='large')
        plt.legend()

    if len(datasetList)==2:
        temp=(np.array(foundParticlesUniqueEntries[1])-np.array(foundParticlesUniqueEntries[0]))/np.array(foundParticlesUniqueEntries[0])*100
        lastPlotData=np.array(foundParticlesUniqueEntries[1])-np.array(foundParticlesUniqueEntries[0])
        lastPlotData=np.where(np.isinf(temp), 0, temp)

        lastPlotLabel="Difference [%]"
        lastPlotTitle=lastPlotLabel +"{}-{}".format(labelList[1],labelList[0])
        tempA=np.array(foundParticlesUniqueEntries[1])
        tempAerr=getParticlesNumbersErrors(datasetList[1],nSlicesErrors)
        tempB=np.array(foundParticlesUniqueEntries[0])
        tempBerr=getParticlesNumbersErrors(datasetList[0],nSlicesErrors)
        tempC=tempB
        tempCerr=tempBerr
        tempTotErr=np.sqrt(np.power(1/tempC,2)*np.power(tempAerr,2) + np.power(1/tempC,2)* np.power(tempBerr,2) + np.power((tempA-tempB)/np.power(tempC,2),2)*np.power(tempCerr,2) )
        tempErr=np.where(np.isinf(tempTotErr), 0, tempTotErr)

    else:
        tempMatrix=np.asmatrix(foundParticlesUniqueEntries)
        lastPlotData=np.array(tempMatrix.std(0)/tempMatrix.mean(0)*100).flatten()
        tempErr=np.zeros(len(lastPlotData))

        lastPlotLabel="All Run RMS / Mean [%]"
        lastPlotTitle=lastPlotLabel
    axs[len(datasetList)+1].bar(np.arange(len(foundParticlesUnique)), lastPlotData, yerr=tempErr*100, align='center', log=False, label=labelList[i], capsize=10)
    axs[len(datasetList)+1].set_title(lastPlotTitle)
    axs[len(datasetList)+1].set_ylim(-200,200)
    axs[len(datasetList)+1].locator_params(axis="y", nbins=20)
    axs[len(datasetList)+1].grid(True, which="both")
    for j, v in enumerate(lastPlotData):
        if v!=0 and tempErr[j]*100<100:
            axs[len(datasetList)+1].text(j , v, "{:.1f}$\pm${:.1f}".format(v,tempErr[j]*100))


    plt.sca(axs[len(datasetList)+1])
    plt.xticks(range(len(foundParticlesUnique)), particleNamesList, size='large')
    plt.xticks(rotation=90)
    plt.ylabel(lastPlotLabel, size='large')


    fig.subplots_adjust(top=0.93)
    fig.tight_layout()
    #plt.subplots_adjust(hspace = 0.1)
    figname=runName+"ParticleDistribution"
    pl.savefig(figname, transparent=False, facecolor='white')


# ### Show only the 4 more frequent particles

# In[20]:


if flagAllPlots:
    fig, axs = plt.subplots(nrows=len(datasetList)+2, ncols=1, figsize=(8,(len(datasetList)+1)*8))
    fig.suptitle(runName+"Particles Frequencies")
    width=(1-0.2)/len(datasetList)
    for i, dataset in enumerate(datasetList):
        axs[i].bar(np.arange(4), foundParticlesUniqueEntries[i][:4], align='center', log=True, yerr=getParticlesNumbersErrors(dataset,nSlicesErrors)[:4], ecolor="blue", capsize=10)
    #    axs[i].bar(range(len(foundParticlesUnique)), getParticlesNumbersErrors(dataset,nSlicesErrors), align='center', log=True)

        axs[i].set_title(labelList[i]) 
        plt.sca(axs[i])
        plt.xticks(range(4), particleNamesList[:4], size='large')
        plt.xticks(rotation=90)
        plt.ylabel("Occurency", size='large')

        for j, v in enumerate(foundParticlesUniqueEntries[i][:4]):
            if v!=0:
                axs[i].text(j , v, "{:.2e}".format(v))
                #axs[len(datasetList)].text(j , v, "{:.2e}".format(v), color="tab:blue")

        axs[len(datasetList)].bar(np.arange(4)+i*width, foundParticlesUniqueEntries[i][:4],width=width, ecolor="blue" , align='center',yerr=getParticlesNumbersErrors(dataset,nSlicesErrors)[:4], log=True, label=labelList[i], alpha=1, capsize=10)
        axs[len(datasetList)].set_title('Comparison')
        plt.sca(axs[len(datasetList)])
        plt.xticks(np.arange(4)+i*width/2, particleNamesList[:4], size='large')
        plt.xticks(rotation=90)
        plt.ylabel("Occurency", size='large')
        plt.legend()

    if len(datasetList)==2:
        temp=(np.array(foundParticlesUniqueEntries[1][:4])-np.array(foundParticlesUniqueEntries[0][:4]))/np.array(foundParticlesUniqueEntries[0][:4])*100
        lastPlotData=np.array(foundParticlesUniqueEntries[1][:4])-np.array(foundParticlesUniqueEntries[0][:4])
        lastPlotData=np.where(np.isinf(temp), 0, temp)

        lastPlotLabel="Difference [%]"
        lastPlotTitle=lastPlotLabel +"{}-{}".format(labelList[1],labelList[0])
        tempA=np.array(foundParticlesUniqueEntries[1][:4])
        tempAerr=getParticlesNumbersErrors(datasetList[1],nSlicesErrors)[:4]
        tempB=np.array(foundParticlesUniqueEntries[0][:4])
        tempBerr=getParticlesNumbersErrors(datasetList[0],nSlicesErrors)[:4]
        tempC=tempB
        tempCerr=tempBerr
        tempTotErr=np.sqrt(np.power(1/tempC,2)*np.power(tempAerr,2) + np.power(1/tempC,2)* np.power(tempBerr,2) + np.power((tempA-tempB)/np.power(tempC,2),2)*np.power(tempCerr,2) )
        tempErr=np.where(np.isinf(tempTotErr), 0, tempTotErr)

    else:
        tempMatrix=np.asmatrix(foundParticlesUniqueEntries)
        lastPlotData=np.array(tempMatrix.std(0)/tempMatrix.mean(0)*100).flatten()[:4]
        tempErr=np.zeros(4)

        lastPlotLabel="All Run RMS / Mean [%]"
        lastPlotTitle=lastPlotLabel
    axs[len(datasetList)+1].bar(np.arange(4), lastPlotData, yerr=tempErr*100, align='center', log=False, label=labelList[i], capsize=10)
    axs[len(datasetList)+1].set_title(lastPlotTitle)
    axs[len(datasetList)+1].set_ylim(-200,200)
    axs[len(datasetList)+1].locator_params(axis="y", nbins=20)
    axs[len(datasetList)+1].grid(True, which="both")
    for j, v in enumerate(lastPlotData):
        if v!=0 and tempErr[j]*100<100:
            axs[len(datasetList)+1].text(j , v, "{:.1f}$\pm${:.1f}".format(v,tempErr[j]*100))


    plt.sca(axs[len(datasetList)+1])
    plt.xticks(range(4), particleNamesList[:4], size='large')
    plt.xticks(rotation=90)
    plt.ylabel(lastPlotLabel, size='large')


    fig.subplots_adjust(top=0.93)
    fig.tight_layout()
    #plt.subplots_adjust(hspace = 0.1)
    figname=runName+"ParticleDistribution"
    pl.savefig(figname, transparent=False, facecolor='white')


# ## With and without Time Cut together

# In[21]:


fig, axs = plt.subplots(nrows=len(datasetList)+2, ncols=1, figsize=(16,(len(datasetList)+1)*8))
fig.suptitle(runName+"Particles Frequencies w/ and w/o Time Cut  tmin={} [ns] tmax={} [ns]".format(timeCut[0],timeCut[1]))
width=(1-0.2)/len(datasetList)
for i, dataset in enumerate(datasetList):
    axs[i].bar(np.arange(4), foundParticlesUniqueEntries[i][:4], align='center', log=True, yerr=getParticlesNumbersErrors(dataset,nSlicesErrors)[:4], edgecolor="tab:blue", ecolor="tab:blue", color="white", capsize=10, alpha=1, label="NoTimeCut")
    axs[i].bar(np.arange(4), foundParticlesUniqueEntriesTimeCut[i][:4], align='center', log=True, yerr=getParticlesNumbersErrors(dataset[(dataset["Time"]>timeCut[0]) & (dataset["Time"]<timeCut[1])],nSlicesErrors)[:4], linestyle="--", edgecolor="tab:blue", ecolor="tab:blue", capsize=10, alpha=1, label="TimeCut")
    axs[i].legend()

#    axs[i].bar(range(len(foundParticlesUnique)), getParticlesNumbersErrors(dataset,nSlicesErrors), align='center', log=True)

    axs[i].set_title(labelList[i]) 
    plt.sca(axs[i])
    plt.xticks(range(4), particleNamesList[:4], size='large')
    plt.xticks(rotation=0)
    plt.ylabel("Occurency", size='large')
    
    for j, v in enumerate(foundParticlesUniqueEntriesTimeCut[i][:4]):
        if v!=0:
            axs[i].text(j , 0.8*v, "{:.2e}".format(v))
#            axs[len(datasetList)].text(j+0.2 if i%2 else j-0.2 , v*0.5, "{:.2e}".format(v), rotation=30)
            axs[len(datasetList)].text(j if i==0 else j+width if i==1 else  j+2*width , v*0.5, "{:.2e}".format(v), rotation=30)


            #axs[len(datasetList)].text(j , v, "{:.2e}".format(v), color="tab:blue")
            
    for j, v in enumerate(foundParticlesUniqueEntries[i][:4]):
        if v!=0:
            axs[i].text(j , 1.2*v, "{:.2e}".format(v))
 #           axs[len(datasetList)].text(j+0.2 if i%2 else j-0.2 , v*1.2, "{:.2e}".format(v), rotation=30)
            axs[len(datasetList)].text(j if i==0 else j+width if i==1 else  j+2*width , v*1.2, "{:.2e}".format(v), rotation=30)


            
    axs[len(datasetList)].bar(np.arange(4)+i*width, foundParticlesUniqueEntries[i][:4],width=width, edgecolor="tab:blue", ecolor="tab:blue", color="white", align='center',yerr=getParticlesNumbersErrors(dataset[(dataset["Time"]>timeCut[0]) & (dataset["Time"]<timeCut[1])],nSlicesErrors)[:4], log=True, label=labelList[i], alpha=1, capsize=10)
    axs[len(datasetList)].bar(np.arange(4)+i*width, foundParticlesUniqueEntriesTimeCut[i][:4],width=width, linestyle="--", edgecolor="tab:blue", ecolor="tab:blue", align='center',yerr=getParticlesNumbersErrors(dataset,nSlicesErrors)[:4], log=True, label=labelList[i]+"_TimeCut", alpha=1, capsize=10)


    #axs[len(datasetList)].set_title('Comparison')
    plt.sca(axs[len(datasetList)])
    plt.xticks(np.arange(4)+i*width/2, particleNamesList[:4], size='large')
    plt.xticks(rotation=0)
    plt.ylabel("Occurency", size='large')
    plt.legend()
    
if len(datasetList)==2:
    temp=(np.array(foundParticlesUniqueEntriesTimeCut[1][:4])-np.array(foundParticlesUniqueEntriesTimeCut[0][:4]))/np.array(foundParticlesUniqueEntriesTimeCut[0][:4])*100
    lastPlotData=np.array(foundParticlesUniqueEntriesTimeCut[1][:4])-np.array(foundParticlesUniqueEntriesTimeCut[0][:4])
    lastPlotData=np.where(np.isinf(temp), 0, temp)

    lastPlotLabel="Difference [%]"
    lastPlotTitle=lastPlotLabel +"{}-{}".format(labelList[1],labelList[0])
    tempA=np.array(foundParticlesUniqueEntriesTimeCut[1][:4])
    tempAerr=getParticlesNumbersErrors(datasetList[1],nSlicesErrors)[:4]
    tempB=np.array(foundParticlesUniqueEntriesTimeCut[0][:4])
    tempBerr=getParticlesNumbersErrors(datasetList[0],nSlicesErrors)[:4]
    tempC=tempB
    tempCerr=tempBerr
    tempTotErr=np.sqrt(np.power(1/tempC,2)*np.power(tempAerr,2) + np.power(1/tempC,2)* np.power(tempBerr,2) + np.power((tempA-tempB)/np.power(tempC,2),2)*np.power(tempCerr,2) )
    tempErr=np.where(np.isinf(tempTotErr), 0, tempTotErr)
    
else:
    tempMatrix=np.asmatrix(foundParticlesUniqueEntries)
    lastPlotData=np.array(tempMatrix.std(0)/tempMatrix.mean(0)*100).flatten()[:4]
    tempErr=np.zeros(4)

    lastPlotLabel="All Run RMS / Mean [%]"
    lastPlotTitle=lastPlotLabel
axs[len(datasetList)+1].bar(np.arange(4), lastPlotData, yerr=tempErr*100, align='center', log=False, label=labelList[i], capsize=10)
axs[len(datasetList)+1].set_title(lastPlotTitle)
axs[len(datasetList)+1].set_ylim(-200,200)
axs[len(datasetList)+1].locator_params(axis="y", nbins=20)
axs[len(datasetList)+1].grid(True, which="both")
for j, v in enumerate(lastPlotData):
    if v!=0 and tempErr[j]*100<100:
        axs[len(datasetList)+1].text(j , v, "{:.1f}$\pm${:.1f}".format(v,tempErr[j]*100))


plt.sca(axs[len(datasetList)+1])
plt.xticks(range(4), particleNamesList[:4], size='large')
plt.xticks(rotation=0)
plt.ylabel(lastPlotLabel, size='large')


fig.subplots_adjust(top=0.93)
fig.tight_layout()
#plt.subplots_adjust(hspace = 0.1)
figname=runName+"ParticleDistributionAlsoTimeCut"
pl.savefig(figname, transparent=False, facecolor='white')


# ### With Time Cut

# In[22]:


if flagAllPlots:
    fig, axs = plt.subplots(nrows=len(datasetList)+2, ncols=1, figsize=(12,(len(datasetList)+1)*8))
    fig.suptitle(runName+"Particles Frequencies Time Cut  tmin={} [ns] tmax={} [ns]".format(timeCut[0],timeCut[1]))
    width=(1-0.2)/len(datasetList)
    for i, dataset in enumerate(datasetList):
        axs[i].bar(np.arange(4), foundParticlesUniqueEntriesTimeCut[i][:4], align='center', log=True, yerr=getParticlesNumbersErrors(dataset[(dataset["Time"]>timeCut[0]) & (dataset["Time"]<timeCut[1])],nSlicesErrors)[:4], ecolor="blue", capsize=10)
    #    axs[i].bar(range(len(foundParticlesUnique)), getParticlesNumbersErrors(dataset,nSlicesErrors), align='center', log=True)

        axs[i].set_title(labelList[i]) 
        plt.sca(axs[i])
        plt.xticks(range(4), particleNamesList[:4], size='large')
        plt.xticks(rotation=90)
        plt.ylabel("Occurency", size='large')

        for j, v in enumerate(foundParticlesUniqueEntriesTimeCut[i][:4]):
            if v!=0:
                axs[i].text(j , v, "{:.2e}".format(v))
                axs[len(datasetList)].text(j+0.2 if i%2 else j-0.2 , v*1.2, "{:.2e}".format(v))
                #axs[len(datasetList)].text(j , v, "{:.2e}".format(v), color="tab:blue")

        axs[len(datasetList)].bar(np.arange(4)+i*width, foundParticlesUniqueEntriesTimeCut[i][:4],width=width, ecolor="blue" , align='center',yerr=getParticlesNumbersErrors(dataset[(dataset["Time"]>timeCut[0]) & (dataset["Time"]<timeCut[1])],nSlicesErrors)[:4], log=True, label=labelList[i], alpha=1, capsize=10)
        axs[len(datasetList)].set_title('Comparison')
        plt.sca(axs[len(datasetList)])
        plt.xticks(np.arange(4)+i*width/2, particleNamesList[:4], size='large')
        plt.xticks(rotation=90)
        plt.ylabel("Occurency", size='large')
        plt.legend()

    if len(datasetList)==2:
        temp=(np.array(foundParticlesUniqueEntriesTimeCut[1][:4])-np.array(foundParticlesUniqueEntriesTimeCut[0][:4]))/np.array(foundParticlesUniqueEntriesTimeCut[0][:4])*100
        lastPlotData=np.array(foundParticlesUniqueEntriesTimeCut[1][:4])-np.array(foundParticlesUniqueEntriesTimeCut[0][:4])
        lastPlotData=np.where(np.isinf(temp), 0, temp)

        lastPlotLabel="Difference [%]"
        lastPlotTitle=lastPlotLabel +"{}-{}".format(labelList[1],labelList[0])
        tempA=np.array(foundParticlesUniqueEntriesTimeCut[1][:4])
        tempAerr=getParticlesNumbersErrors(datasetList[1],nSlicesErrors)[:4]
        tempB=np.array(foundParticlesUniqueEntriesTimeCut[0][:4])
        tempBerr=getParticlesNumbersErrors(datasetList[0],nSlicesErrors)[:4]
        tempC=tempB
        tempCerr=tempBerr
        tempTotErr=np.sqrt(np.power(1/tempC,2)*np.power(tempAerr,2) + np.power(1/tempC,2)* np.power(tempBerr,2) + np.power((tempA-tempB)/np.power(tempC,2),2)*np.power(tempCerr,2) )
        tempErr=np.where(np.isinf(tempTotErr), 0, tempTotErr)

    else:
        tempMatrix=np.asmatrix(foundParticlesUniqueEntries)
        lastPlotData=np.array(tempMatrix.std(0)/tempMatrix.mean(0)*100).flatten()[:4]
        tempErr=np.zeros(4)

        lastPlotLabel="All Run RMS / Mean [%]"
        lastPlotTitle=lastPlotLabel
    axs[len(datasetList)+1].bar(np.arange(4), lastPlotData, yerr=tempErr*100, align='center', log=False, label=labelList[i], capsize=10)
    axs[len(datasetList)+1].set_title(lastPlotTitle)
    axs[len(datasetList)+1].set_ylim(-200,200)
    axs[len(datasetList)+1].locator_params(axis="y", nbins=20)
    axs[len(datasetList)+1].grid(True, which="both")
    for j, v in enumerate(lastPlotData):
        if v!=0 and tempErr[j]*100<100:
            axs[len(datasetList)+1].text(j , v, "{:.1f}$\pm${:.1f}".format(v,tempErr[j]*100))


    plt.sca(axs[len(datasetList)+1])
    plt.xticks(range(4), particleNamesList[:4], size='large')
    plt.xticks(rotation=90)
    plt.ylabel(lastPlotLabel, size='large')


    fig.subplots_adjust(top=0.93)
    fig.tight_layout()
    #plt.subplots_adjust(hspace = 0.1)
    figname=runName+"ParticleDistributionTimeCut"
    pl.savefig(figname, transparent=False, facecolor='white')


# ## Count Particle Numbers

# In[23]:


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


# ## Plot Energy Spectra - Lethargy Plots

# ### New version - No Time Cut

# In[24]:


plotLethargy(datasetList, nbins=200, logY=True, logX=False, yrange=(1e2,3e7))


# ### New version - With Time Cut

# In[25]:


plotLethargy(datasetList, nbins=200, logY=True, logX=False, yrange=(1e2,3e7), trange=timeCut, xrange=[-6,1])


# ### With and Without Time Cut

# In[26]:


plotLethargyAlsoWithTimeCut(datasetList, nbins=200, logY=True, logX=False, yrange=(1e2,3e7), trange=timeCut)


# In[27]:


if flagAllPlots:
    plotAllEnergySpectra(datasetList, nbins=nbins, logY=True, logX=False)
    plotAllEnergySpectra(datasetList, nbins=10000, logY=True, logX=True)


# In[28]:


### Photons and e+/e-
if flagAllPlots:
    xmax=np.percentile(getMomentum(datasetList[0],22,"KinE").to_numpy(),99.999) # Get rid of outliers
    plotMomenta(datasetList=datasetList, particleList=[22, [11,-11]], particleLabel=["$\gamma$", "$e^-~/~e^+$"], title="EneEM", nbins=nbinsH, figName="EneEGamma", xrange=[0,xmax])
    ### Hadrons
    plotMomenta(datasetList=datasetList, particleList=[2112, listChargedHadrons], particleLabel=["Neutrons","Ch. Had"], title="EneHad", nbins=nbinsH, figName="EneHadrons", xrange=[0,1])


# ## Plot Time Distributions

# In[29]:


plotVariablePerEachRelevantParticle(datasetList=datasetList, variable="Time", plotTitle="Time Distribution", xlabel="t [ns]", ylabel="Arb. Units", nbins=200, log=True, figTitle="Time", xrange=(-25,100))


# In[30]:


plotVariablePerEachRelevantParticle(datasetList=datasetList, variable="Time", plotTitle="Time DistributionZoom", xlabel="t [ns]", ylabel="Arb. Units", nbins=200, log=True, figTitle="TimeZoom", xrange=(-1,15))


# ## Plot Muons' Decay Position

# ### Global

# In[31]:


fig=plt.figure(figsize=(6,5))
plt.suptitle(runName+"Muon Decay Z")
plt.gca().set_xlabel('$z_{\mu \,dec}$ [m]')
plt.gca().set_ylabel('Arb. Units')

for i, dataset in enumerate(datasetList):
    fig.gca().hist(dataset["PosZmu"]/100,histtype='step', bins=nbinsH, weights=dataset["Weight"], log=True, label=labelList[i])

plt.legend()
#plt.ylim((100, 1e9))
figname=runName+"MuDec"
pl.savefig(figname,transparent=False, facecolor='white')


# ### Z Position Plots

# In[32]:


plotVariablePerEachRelevantParticle(datasetList=datasetList, variable="PosZ", plotTitle="Zexit", xlabel="z [cm]", ylabel="Arb. Units", nbins=200, log=True, figTitle="PosZ")

plotVariablePerEachRelevantParticle(datasetList=datasetList, variable="PosZ", plotTitle="ZexitZoom", xlabel="z [cm]", ylabel="Arb. Units", nbins=200, log=True, figTitle="PosZZoom", xrange=[-200,200])

plotVariablePerEachRelevantParticle(datasetList=datasetList, variable="PosZ", plotTitle="ZexitZoom", xlabel="z [cm]", ylabel="Arb. Units", nbins=200, log=True, figTitle="PosZZoomTimeCut", xrange=[-200,200],trange=[-1,15])

plotVariablePerEachRelevantParticle(datasetList=datasetList, variable="PosZ", plotTitle="ZexitZoom w/ and w/o Time Cut", xlabel="z [cm]", ylabel="Arb. Units", nbins=200, log=True, figTitle="PosZZoomAlsoTimeCut", xrange=[-200,200],trange=[-1,15], alsoWithTime=True)


# fig=plt.figure(figsize=(14,5))
# plt.suptitle("z$_{\mu}$ decay")
# plt.subplot(1, 3, 1)
# plt.gca().set_title(labelList[0])
# plt.hist(getInfo(datasetList[0],22,"PosZmu")/100,histtype='step',range=[-25,25], weights=getInfo(datasetList[0], 22, "Weight"), bins=100,color='r',linestyle=":", label= '$\gamma$ ')
# plt.hist(getInfo(datasetList[0],11,"PosZmu")/100,histtype='step',range=[-25,25], weights=getInfo(datasetList[0], 11, "Weight"), bins=100,color='k',linestyle=":", label= '$e^-$ ')
# plt.hist(getInfo(datasetList[0],-11,"PosZmu")/100,histtype='step',range=[-25,25], weights=getInfo(datasetList[0], -11, "Weight"), bins=100,color='y',linestyle=":", label= '$e^+$ ')
# plt.hist(getInfo(datasetList[0],2112,"PosZmu")/100,histtype='step', range=[-25,25], weights=getInfo(datasetList[0], 2112, "Weight"),bins=100,color='blue',linestyle=":", label= 'n ')
# plt.xlabel('z$_{\mu}$ (m)',fontsize=14)
# plt.ylabel('dN/dz$_{\mu}$',fontsize=14)
# plt.yscale('log')
# plt.legend(loc= 'upper right')
# plt.ylim((1e2,2e7))
# plt.subplot(1, 3, 2)
# plt.gca().set_title(labelList[1])
# plt.hist(getInfo(datasetList[1],22,"PosZmu")/100,histtype='step',range=[-25,25], weights=getInfo(datasetList[1], 22, "Weight"), bins=100,color='r', label= '$\gamma$ ')
# plt.hist(getInfo(datasetList[1],11,"PosZmu")/100,histtype='step',range=[-25,25], weights=getInfo(datasetList[1], 11, "Weight"), bins=100,color='k', label= '$e^-$ ')
# plt.hist(getInfo(datasetList[1],-11,"PosZmu")/100,histtype='step',range=[-25,25], weights=getInfo(datasetList[1], -11, "Weight"), bins=100,color='y', label= '$e^+$ ')
# plt.hist(getInfo(datasetList[1],2112,"PosZmu")/100,histtype='step', range=[-25,25], weights=getInfo(datasetList[1], 2112, "Weight"),bins=100,color='blue', label= 'n ')
# plt.xlabel('z$_{\mu}$ (m)',fontsize=14)
# plt.ylabel('dN/dz$_{\mu}$',fontsize=14)
# plt.yscale('log')
# plt.legend(loc= 'upper right')
# plt.ylim((1e2,2e7))
# plt.subplot(1, 3, 3)
# plt.gca().set_title(labelList[0]+" & "+labelList[1])
# plt.hist( getInfo(datasetList[0],22,"PosZmu")/100,histtype='step',range=[-25,25], weights=getInfo(datasetList[0], 22, "Weight"), bins=100,color='r',linestyle=":", label= '$\gamma$ ')
# plt.hist( getInfo(datasetList[0],11,"PosZmu")/100,histtype='step',range=[-25,25], weights=getInfo(datasetList[0], 11, "Weight"), bins=100,color='k',linestyle=":", label= '$e^-$ ')
# plt.hist( getInfo(datasetList[0],-11,"PosZmu")/100,histtype='step',range=[-25,25], weights=getInfo(datasetList[0], -11, "Weight"), bins=100,color='y',linestyle=":", label= '$e^+$ ')
# plt.hist( getInfo(datasetList[0],2112,"PosZmu")/100,histtype='step', range=[-25,25], weights=getInfo(datasetList[0], 2112, "Weight"),bins=100,color='blue',linestyle=":", label= 'n ')
# plt.hist( getInfo(datasetList[1],22,"PosZmu")/100,histtype='step',range=[-25,25], weights=getInfo(datasetList[1], 22, "Weight"), bins=100,color='r', label= '$\gamma$ ')
# plt.hist( getInfo(datasetList[1],11,"PosZmu")/100,histtype='step',range=[-25,25], weights=getInfo(datasetList[1], 11, "Weight"), bins=100,color='k', label= '$e^-$ ')
# plt.hist( getInfo(datasetList[1],-11,"PosZmu")/100,histtype='step',range=[-25,25], weights=getInfo(datasetList[1], -11, "Weight"), bins=100,color='y', label= '$e^+$ ')
# plt.hist( getInfo(datasetList[1],2112,"PosZmu")/100,histtype='step', range=[-25,25], weights=getInfo(datasetList[1], 2112, "Weight"),bins=100,color='blue', label= 'n ')
# plt.xlabel('z$_{\mu}$ (m)',fontsize=14)
# plt.ylabel('dN/dz$_{\mu}$',fontsize=14)
# plt.yscale('log')
# #plt.legend(loc= 'upper left')
# plt.ylim((1e2,2e7))
# fig.subplots_adjust(left = 0.05,right = 0.99,wspace = 0.5, hspace = 0.5,top=0.85, bottom= 0.2)
# figname=runName+"PosZmu"
# pl.savefig(figname)

# In[86]:


fig, axs = plt.subplots(nrows=len(datasetList), ncols=1, figsize=(10,13), sharex=False)
binWidthZ=0.5

for datasetNumber, dataset in enumerate(datasetList):
    nBinZ.append(int(((dataset["PosZmu"]/100).max()-(dataset["PosZmu"]/100).min())/binWidthZ))
    
    plot1D(axs[datasetNumber],getInfo(datasetList[0],22,"PosZmu")/100, weights=getInfo(datasetList[0],22,"Weight")/100,bins=nBinZ[-1], plotTitle="Muon Decay Z", label="$\gamma$", xlabel='', ylabel='dN/dz$_{\mu}$', col="r" )
    plot1D(axs[datasetNumber],getInfo(datasetList[0],11,"PosZmu")/100, weights=getInfo(datasetList[0],11,"Weight")/100,bins=nBinZ[-1], plotTitle="Muon Decay Z", label="e-", xlabel='', ylabel='dN/dz$_{\mu}$', col="k" )
    plot1D(axs[datasetNumber],getInfo(datasetList[0],-11,"PosZmu")/100, weights=getInfo(datasetList[0],-11,"Weight")/100,bins=nBinZ[-1], plotTitle="Muon Decay Z", label="e+", xlabel='', ylabel='dN/dz$_{\mu}$', col="y" )
    plot1D(axs[datasetNumber],getInfo(datasetList[0],2112,"PosZmu")/100, weights=getInfo(datasetList[0],2112,"Weight")/100,bins=nBinZ[-1], plotTitle="Muon Decay Z", label="n", xlabel='', ylabel='dN/dz$_{\mu}$', col="b" )
    axs[datasetNumber].twinx().hist(dataset["PosZmu"]/100, weights=dataset["Weight"], bins=nBinZ[-1], cumulative=1,histtype='step', label="cum "+labelList[datasetNumber], density=True, linestyle=':', color="black")
    axs[datasetNumber].legend(loc= "best")
    axs[datasetNumber].set_title(labelList[datasetNumber])
    axs[datasetNumber].set_xlabel("$Z_{\mu} [m]$")

axs[0].set_xlim(0,45)
#axs[0].axis(ymin=100, ymax=1e9)
#axs[0].grid(True, which="both", axis='y')
#axs[0].locator_params(axis="x", nbins=20)
#axs[2].grid(True, which="both")
#axs[2].set_title("Cumulative Function Norm")
#axs[2].legend(loc= "best", fontsize='x-small')

#axs[1].legend(loc= "best")
#axs[1].set_title("Cumulative Function Not Norm")
#axs[1].grid(True, which="both")

figname=runName+"MuDecZComponents"
pl.savefig(figname,transparent=False, facecolor='white')


# ### Z vs x

# In[34]:


fig, ax = plt.subplots(nrows=len(datasetList), ncols=1, figsize=(9,len(datasetList)*4))
plt.suptitle(runName+"Nozzle")
if len(datasetList)>1:
    for i, dataset in enumerate(datasetList):
        ax[i].set_title(labelList[i])
        ax[i].hist2d(dataset["PosZ"],dataset["PosX"],weights=dataset["Weight"],norm=matplotlib.colors.LogNorm(vmin=1,vmax=1e7),bins=500, range=([-800,800],[-250,250]), cmap='plasma')
        ax[i].axis(xmin=-800, xmax=800)
        ax[i].axis(ymin=-250, ymax=250)
        ax[i].set_ylabel('x [cm]',fontsize=14)
        PCM=ax[i].get_children()[0] #get the mappable, the 1st and the 2nd are the x and y axes
        cb=plt.colorbar(PCM, ax=ax[i]) 
        cb.set_label('Arb. Units')
else:
    ax.set_title(labelList[i])
    ax.hist2d(dataset["PosZ"],dataset["PosX"],weights=dataset["Weight"],norm=matplotlib.colors.LogNorm(vmin=1,vmax=1e7),bins=500, range=([-800,800],[-250,250]), cmap='plasma')
    ax.axis(xmin=-800, xmax=800)
    ax.axis(ymin=-250, ymax=250)
    ax.set_ylabel('x [cm]',fontsize=14)
    
plt.xlabel('z (cm)',fontsize=14)
#plt.ylabel('x (cm)',fontsize=14)

figname=runName+"ZvsX_FLUKA"
pl.savefig(figname,transparent=False, facecolor='white')


# In[35]:


fig, ax = plt.subplots(nrows=len(datasetList), ncols=1, figsize=(9,len(datasetList)*4))
plt.suptitle(runName+"Nozzle")
if len(datasetList)>1:
    for i, dataset in enumerate(datasetList):
        ax[i].set_title(labelList[i])
        ax[i].hist2d(dataset["PosZ"],dataset["PosX"],weights=dataset["Weight"],norm=matplotlib.colors.LogNorm(vmin=1,vmax=1e5),bins=500, range=([-80,80],[-25,25]), cmap='plasma')
        ax[i].axis(xmin=-80, xmax=80)
        ax[i].axis(ymin=-25, ymax=25)
        ax[i].set_ylabel('x [cm]',fontsize=14)
        PCM=ax[i].get_children()[0] #get the mappable, the 1st and the 2nd are the x and y axes
        cb=plt.colorbar(PCM, ax=ax[i]) 
        cb.set_label('Arb. Units')
else:
    ax.set_title(labelList[i])
    ax.hist2d(dataset["PosZ"],dataset["PosX"],weights=dataset["Weight"],norm=matplotlib.colors.LogNorm(vmin=1,vmax=1e5),bins=500, range=([-80,80],[-25,25]), cmap='plasma')
    ax.axis(xmin=-80, xmax=80)
    ax.axis(ymin=-25, ymax=25)
    ax.set_ylabel('x [cm]',fontsize=14)
    
plt.xlabel('z (cm)',fontsize=14)
#plt.ylabel('x (cm)',fontsize=14)

figname=runName+"ZvsX_FLUKAZoom"
pl.savefig(figname,transparent=False, facecolor='white')


# ### Theta vs E for BIB electrons

# In[36]:


if flagAllPlots:
    fig, ax = plt.subplots(nrows=1, ncols=len(datasetList)+1, figsize=((len(datasetList)+1)*9,8), sharey=False)
    plt.suptitle("")

    for i, dataset in enumerate(datasetEleList):
        ax[i].set_xlabel("Theta",fontsize='14')
        ax[i].set_ylabel("E [GeV]",fontsize='14')
    #    ax[i].set_ylabel(ylabel,fontsize='14')
       # ax[i].hist(dataset[variable],histtype='step', bins=nbins, weights=dataset["Weight"], log=log, range=xrange)
        ax[i].hist2d(np.sqrt(dataset[dataset["PosZmu"]<1400]["CX"]**2+dataset[dataset["PosZmu"]<1400]["CY"]**2), dataset[dataset["PosZmu"]<1400]["EneEle"], range=([0,0.0005],[0,1500]), weights=dataset[dataset["PosZmu"]<1400]["Weight"], bins=100, norm=matplotlib.colors.LogNorm(), cmap='Reds')
        PCM=ax[i].get_children()[0] #get the mappable, the 1st and the 2nd are the x and y axes
        cb=plt.colorbar(PCM, ax=ax[i]) 
        cb.set_label('Arb. Units')
        ax[i].set_title(labelList[i])
    #        ax[len(datasetList)].hist(dataset[variable],histtype='step', bins=nbins, weights=dataset["Weight"], log=log, range=xrange, label=labelList[i])
        ax[len(datasetList)].hist2d(np.sqrt(dataset[dataset["PosZmu"]<1400]["CX"]**2+dataset[dataset["PosZmu"]<1400]["CY"]**2), dataset[dataset["PosZmu"]<1400]["EneEle"], range=([0,0.0005],[0,1500]),weights=dataset[dataset["PosZmu"]<1400]["Weight"], bins=nbins, norm=matplotlib.colors.LogNorm())
      #  PCMb=ax[len(datasetList)].get_children()[0] #get the mappable, the 1st and the 2nd are the x and y axes
      #  plt.colorbar(PCMb, ax=ax[len(datasetList)]) 

    ax[len(datasetList)].set_title("Comparison")
    # ax[len(datasetList)].legend()
    ax[len(datasetList)].set_xlabel("Theta",fontsize='14')
    ax[len(datasetList)].set_ylabel("E [GeV]",fontsize='14')

    figname=runName+"ThetaVsE"
    pl.savefig(figname,transparent=False, facecolor='white')   


# In[37]:


if flagAllPlots:
    plotSingleVariable2D(datasetList=datasetList, variableX="PosZFI", variableY="PosZ", plotTitle="Parent Electron First Interaction vs Z Exit Bib",
    xlabel="Z$_{e}$ [cm]", ylabel="Z [cm]", nbins=nbins*2, log=True, vmin=1,vmax=1e7,
    figTitle="ZFIvsZBib", range=[[-200, 1000], [-200, 800]], hasBIB=-1)


# In[38]:


fig, ax = plt.subplots(nrows=1, ncols=len(datasetList)+1, figsize=((len(datasetList)+1)*9,8), sharey=False)
plt.suptitle("")

for i, dataset in enumerate(datasetList):
    ax[i].set_ylabel("x_e [m]",fontsize='14')
    ax[i].set_xlabel("z_e [m]",fontsize='14')
    ax[i].hist2d(dataset["PosZFI"]/100,dataset["PosXFI"]/100, range=([-5,100],[-5,0.3]), weights=dataset["Weight"], bins=100, norm=matplotlib.colors.LogNorm(vmin=1,vmax=1e7), cmap='Reds')
    
    PCM=ax[i].get_children()[0] #get the mappable, the 1st and the 2nd are the x and y axes
    cb=plt.colorbar(PCM, ax=ax[i]) 
    cb.set_label('Arb. Units')
    ax[i].set_title(labelList[i])
    ax[len(datasetList)].hist2d(dataset["PosZFI"]/100,dataset["PosXFI"]/100, range=([-5,100],[-5,0.3]),weights=dataset["Weight"], bins=100, norm=matplotlib.colors.LogNorm(vmin=1,vmax=1e7))
    
    ax[len(datasetList)].set_title("Comparison")
    # ax[len(datasetList)].legend()
    ax[len(datasetList)].set_xlabel("z_e [cm]",fontsize='14')
    ax[len(datasetList)].set_ylabel("x_e [cm]",fontsize='14')

    figname=runName+"z_evsx_e"
    pl.savefig(figname,transparent=False, facecolor='white')   


# In[39]:


fig, ax = plt.subplots(nrows=1, ncols=len(datasetList)+1, figsize=((len(datasetList)+1)*9,8), sharey=False)
plt.suptitle("")

for i, dataset in enumerate(datasetList):
    ax[i].set_ylabel("x_e [m]",fontsize='14')
    ax[i].set_xlabel("z_e [m]",fontsize='14')
    ax[i].hist2d(dataset["PosZFI"]/100,dataset["PosXFI"]/100, range=([-2,20],[-0.1,0.1]), weights=dataset["Weight"], bins=100, norm=matplotlib.colors.LogNorm(vmin=1,vmax=1e7), cmap='Reds')
    
    PCM=ax[i].get_children()[0] #get the mappable, the 1st and the 2nd are the x and y axes
    cb=plt.colorbar(PCM, ax=ax[i]) 
    cb.set_label('Arb. Units')
    ax[i].set_title(labelList[i])
    ax[len(datasetList)].hist2d(dataset["PosZFI"]/100,dataset["PosXFI"]/100, range=([-2,20],[-0.1,0.1]),weights=dataset["Weight"], bins=100, norm=matplotlib.colors.LogNorm(vmin=1,vmax=1e7))
    
    ax[len(datasetList)].set_title("Comparison")
    # ax[len(datasetList)].legend()
    ax[len(datasetList)].set_xlabel("z_e [cm]",fontsize='14')
    ax[len(datasetList)].set_ylabel("x_e [cm]",fontsize='14')

    figname=runName+"z_evsx_e_zoom"
    pl.savefig(figname,transparent=False, facecolor='white')  


# In[40]:


fig, ax = plt.subplots(nrows=1, ncols=len(datasetList)+1, figsize=((len(datasetList)+1)*9,8), sharey=False)
plt.suptitle("")

for i, dataset in enumerate(datasetList):
    ax[i].set_ylabel("y_e [m]",fontsize='14')
    ax[i].set_xlabel("z_e [m]",fontsize='14')
    ax[i].hist2d(dataset["PosZFI"]/100,dataset["PosYFI"]/100, range=([-5,100],[-0.1,0.1]), weights=dataset["Weight"], bins=100, norm=matplotlib.colors.LogNorm(vmin=1,vmax=1e7), cmap='Reds')
    
    PCM=ax[i].get_children()[0] #get the mappable, the 1st and the 2nd are the x and y axes
    cb=plt.colorbar(PCM, ax=ax[i]) 
    cb.set_label('Arb. Units')
    ax[i].set_title(labelList[i])
    ax[len(datasetList)].hist2d(dataset["PosZFI"]/100,dataset["PosYFI"]/100, range=([-5,100],[-0.1,0.1]),weights=dataset["Weight"], bins=100, norm=matplotlib.colors.LogNorm(vmin=1,vmax=1e7))
    
    ax[len(datasetList)].set_title("Comparison")
    # ax[len(datasetList)].legend()
    ax[len(datasetList)].set_xlabel("z_e [cm]",fontsize='14')
    ax[len(datasetList)].set_ylabel("y_e [cm]",fontsize='14')

    figname=runName+"z_evsy_e"
    pl.savefig(figname,transparent=False, facecolor='white')  


# In[41]:


fig, ax = plt.subplots(nrows=1, ncols=len(datasetList)+1, figsize=((len(datasetList)+1)*9,8), sharey=False)
plt.suptitle("")

for i, dataset in enumerate(datasetList):
    ax[i].set_ylabel("y_e [m]",fontsize='14')
    ax[i].set_xlabel("z_e [m]",fontsize='14')
    ax[i].hist2d(dataset["PosZFI"]/100,dataset["PosYFI"]/100, range=([-5,10],[-0.1,0.1]), weights=dataset["Weight"], bins=100, norm=matplotlib.colors.LogNorm(vmin=1,vmax=1e7), cmap='Reds')
    
    PCM=ax[i].get_children()[0] #get the mappable, the 1st and the 2nd are the x and y axes
    cb=plt.colorbar(PCM, ax=ax[i]) 
    cb.set_label('Arb. Units')
    ax[i].set_title(labelList[i])
    ax[len(datasetList)].hist2d(dataset["PosZFI"]/100,dataset["PosYFI"]/100, range=([-5,10],[-0.1,0.1]),weights=dataset["Weight"], bins=100, norm=matplotlib.colors.LogNorm(vmin=1,vmax=1e7))
    
    ax[len(datasetList)].set_title("Comparison")
    # ax[len(datasetList)].legend()
    ax[len(datasetList)].set_xlabel("z_e [cm]",fontsize='14')
    ax[len(datasetList)].set_ylabel("y_e [cm]",fontsize='14')

    figname=runName+"z_evsy_e_zoom"
    pl.savefig(figname,transparent=False, facecolor='white')  


# ## Parent Electron Plots

# In[42]:


if flagReadEle:
    print("Plots regarding parent electrons requested")
    plotSingleVariable(datasetList=datasetEleList, variable="EneEle", plotTitle="Energy of Decay Electron [BIB + NoBIB]", xlabel="E$_{e}$ [GeV]", ylabel="Arb. Units", nbins=100, log=False, figTitle="EnElGenitore", hasBIB=-1)
    plotSingleVariable(datasetList=datasetEleList, variable="EneEle", plotTitle="Energy of Decay Electron That DID Generate BIB", xlabel="E$_{e}$ [GeV]", ylabel="Arb. Units", nbins=100, log=False, figTitle="EnElGenitoreBib", hasBIB=1)
    plotSingleVariable(datasetList=datasetEleList, variable="EneEle", plotTitle="Energy of Decay Electron That DID NOT Generate BIB", xlabel="E$_{e}$ [GeV]", ylabel="Arb. Units", nbins=100, log=False, figTitle="EnElGenitoreNoBib", hasBIB=0)

    plotSingleVariable(datasetList=datasetEleList, variable="PosZEle", plotTitle="Z of First Interaction of Decay Electron [BIB + NoBIB]", xlabel="Z$_{e}$ [cm]", ylabel="Arb. Units", nbins=100, log=True, figTitle="ZElGenitore", xrange=[-200,1500], hasBIB=-1)
    plotSingleVariable(datasetList=datasetEleList, variable="PosZEle", plotTitle="Z of First Interaction of Decay Electron That DID Generate BIB", xlabel="Z$_{e}$ [cm]", ylabel="Arb. Units", nbins=100, log=True, figTitle="ZElGenitoreBib", xrange=[-200,1500], hasBIB=1)
    plotSingleVariable(datasetList=datasetEleList, variable="PosZEle", plotTitle="Z of First Interaction of Decay Electron That DID NOT Generate BIB", xlabel="Z$_{e}$ [cm]", ylabel="Arb. Units", nbins=100, log=True, figTitle="ZElGenitoreNoBib", xrange=[-200,1500], hasBIB=0)

    if flagAllPlots:
        plotSingleVariable(datasetList=datasetEleList, variable="CZ", plotTitle="CosZ of Decay Electron [BIB + NoBIB]", xlabel="cos(z) ", ylabel="Arb. Units", nbins=40, log=True, figTitle="CZElGenitore", hasBIB=-1)
        plotSingleVariable(datasetList=datasetEleList, variable="CZ", plotTitle="CosZ of Decay Electron That DID Generate BIB", xlabel="cos(z) ", ylabel="Arb. Units", nbins=40, log=True, figTitle="CZElGenitoreBib", hasBIB=1)
        plotSingleVariable(datasetList=datasetEleList, variable="CZ", plotTitle="CosZ of Decay Electron That DID NOT Generate BIB", xlabel="cos(z) ", ylabel="Arb. Units", nbins=40, log=True, figTitle="CZElGenitoreNoBib", hasBIB=0)

        plotSingleDistributionLongMomNozzle(datasetList=datasetEleList, plotTitle="Longitudinal Momentum In Nozzle Bib", xlabel="P [GeV] ", ylabel="Arb. Units", nbins=60, log=True, figTitle="LongMomNozzleBib",secondaryFlag=True)
        plotSingleDistributionLongMomNozzle(datasetList=datasetEleList, plotTitle="Longitudinal Momentum In Nozzle NoBib", xlabel="P [GeV] ", ylabel="Arb. Units", nbins=60, log=True, figTitle="LongMomNozzleNoBib",secondaryFlag=False)

        plotSingleDistributionTransMomNozzle(datasetList=datasetEleList, plotTitle="Transverse Momentum In Nozzle Bib", xlabel="P [GeV] ", ylabel="Arb. Units", nbins=60, log=True, figTitle="TransMomNozzleBib",secondaryFlag=True)
        plotSingleDistributionTransMomNozzle(datasetList=datasetEleList, plotTitle="Transverse Momentum In Nozzle NoBib", xlabel="P [GeV] ", ylabel="Arb. Units", nbins=60, log=True, figTitle="TransMomNozzleNoBib",secondaryFlag=False)

        plotSingleDistributionBendingRadius(datasetList=datasetEleList, plotTitle="Bending Radius in B Bib", xlabel="r [m] ", ylabel="Arb. Units", nbins=60, log=True, figTitle="BendRadBib",secondaryFlag=True)
        plotSingleDistributionBendingRadius(datasetList=datasetEleList, plotTitle="Bending Radius in B NoBib", xlabel="r [m] ", ylabel="Arb. Units", nbins=60, log=True, figTitle="BendRadNoBib",secondaryFlag=False)

    plotSingleVariable2D(datasetList=datasetEleList, variableX="EneEle", variableY="PosZmu", plotTitle="Decay Electron Ene vs Z of Muon Decay BIB", 
                             xlabel="E$_{e}$ [GeV]", ylabel="$z_{\mu}$ [cm]", nbins=nbins*2, log=True,vmin=1e2,vmax=5e4, 
                             figTitle="Eevszmu_ElGenitore_Bib", range=[[0, 1500], [0, 3000]], hasBIB=1)
    plotSingleVariable2D(datasetList=datasetEleList, variableX="EneEle", variableY="PosZmu", plotTitle="Decay Electron Ene vs Z of Muon Decay NoBIB", 
                             xlabel="E$_{e}$ [GeV]", ylabel="$z_{\mu}$ [cm]", nbins=nbins*2, log=True, vmin=1e2,vmax=5e4,
                             figTitle="Eevszmu_ElGenitore_NoBib", range=[[0, 1500], [0, 3000]], hasBIB=0)

    plotSingleVariable2D(datasetList=datasetEleList, variableX="EneEle", variableY="PosZEle", plotTitle="Decay Electron Ene vs Z of First Interaction BIB", 
                             xlabel="E$_{e}$ [GeV]", ylabel="$z_{e}$ [cm]", nbins=nbins*2, log=True,vmin=1e2,vmax=5e4, 
                             figTitle="Eevszele_ElGenitore_Bib", range=[[0, 1500], [-200, 800]], hasBIB=1)

    plotSingleVariable2D(datasetList=datasetEleList, variableX="EneEle", variableY="PosZEle", plotTitle="Decay Electron Ene vs Z of First Interaction NoBIB", 
                             xlabel="E$_{e}$ [GeV]", ylabel="$z_{e}$ [cm]", nbins=nbins*2, log=True,vmin=1e2,vmax=5e4,
                             figTitle="Eevszele_ElGenitore_NoBib", range=[[0, 1500], [-200, 800]], hasBIB=0)
    
    plotSingleVariable2D(datasetList=datasetEleList, variableX="PosZmu", variableY="PosZEle", plotTitle="Decay Electron Z vs Muon Decay Z [BIB + NoBIB]", 
                             xlabel="$z_{\mu}$ [cm]", ylabel="$z_{e}$ [cm]", nbins=nbins*2, log=True, vmin=1e2,vmax=5e4, 
                             figTitle="ZMuVsZEle_ElGenitore_Bib", range=[[0, 4000], [-200, 800]], hasBIB=-1)
    plotSingleVariable2D(datasetList=datasetEleList, variableX="PosZmu", variableY="PosZEle", plotTitle="Decay Electron Z vs Muon Decay Z BIB", 
                             xlabel="$z_{\mu}$ [cm]", ylabel="$z_{e}$ [cm]", nbins=nbins*2, log=True, vmin=1e2,vmax=5e4, 
                             figTitle="ZMuVsZEle_ElGenitore_Bib", range=[[0, 4000], [-200, 800]], hasBIB=1)
    plotSingleVariable2D(datasetList=datasetEleList, variableX="PosZmu", variableY="PosZEle", plotTitle="Decay Electron Z vs Muon Decay Z NoBIB", 
                             xlabel="$z_{\mu}$ [cm]", ylabel="$z_{e}$ [m]", nbins=nbins*2, log=True, vmin=1e2,vmax=5e4,
                             figTitle="ZMuVsZEle_ElGenitore_NoBib", range=[[0, 4000], [-200, 800]], hasBIB=0)
    
    if flagAllPlots:
        plotSingleVariable2D(datasetList=datasetEleList, variableX="EneEle", variableY="CY", plotTitle="Elettrone Genitore Ene vs CY Bib", 
                                 xlabel="E$_{e}$ [GeV]", ylabel="cos(y)", nbins=nbins*2, log=True, 
                                 figTitle="Eevscosz_ElGenitore_Bib", range=[[0, 1500], [-0.002, 0.002]], hasBIB=1)
        plotSingleVariable2D(datasetList=datasetEleList, variableX="EneEle", variableY="CY", plotTitle="Elettrone Genitore Ene vs CY NoBib", 
                                 xlabel="E$_{e}$ [GeV]", ylabel="cos(y)", nbins=nbins*2, log=True, 
                             figTitle="Eevscosz_ElGenitore_NoBib", range=[[0, 1500], [-0.002, 0.002]], hasBIB=0)
    
    if flagAllPlots:
        ## Stack Plots
        plotStackElePlotsBibNoBib(datasetList=datasetEleList, variable="EneEle", nbins=100, title="Parent Electron Energy Bib/NoBib", xlabel="E$_{e}$ [m]", figname="ParentEleEneStack")

        plotStackElePlotsBibNoBib(datasetList=datasetEleList, variable="PosZEle", nbins=200, title="Parent Electron First Interaction Point Bib/NoBib", xlabel="Z$_{e}$ [cm]", figname="ParentEleZStack")

        plotStackElePlotsBibNoBib(datasetList=datasetEleList, variable="CZ", nbins=50, title="Parent Electron Energy Bib/NoBib", xlabel="cos(z)", figname="ParentEleZStack", log=True)
        ## Energy Cut Plots
        plotEleDistrWithCut(dataset=datasetEleList[0], variable="PosXEle", cutVariable="EneEle", cutCenter=500, cutRange=[100, 50, 10],nbins=100, log=True, xlabel="PosXEle [cm]", figname="PosXEleCutOnEneEle")

        plotEleDistrWithCut(dataset=datasetEleList[0], variable="PosYEle", cutVariable="EneEle", cutCenter=500, cutRange=[100, 50, 10],nbins=100, log=True, xlabel="PosYEle [cm]", figname="PosYEleCutOnEneEle")

        plotEleDistrWithCut(dataset=datasetEleList[0], variable="PosZEle", cutVariable="EneEle", cutCenter=500, cutRange=[100, 50, 10],nbins=100, log=True, xlabel="PosZEle [cm]", figname="PosZEleCutOnEneEle")

        plotEleDistrWithCut(dataset=datasetEleList[0],variable="PosXmu", cutVariable="EneEle", cutCenter=500, cutRange=[100, 50, 10],nbins=100, log=True, xlabel="PosXmu [cm]", figname="PosXMuCutOnEneEle")

        plotEleDistrWithCut(dataset=datasetEleList[0],variable="PosYmu", cutVariable="EneEle", cutCenter=500, cutRange=[100, 50, 10],nbins=100, log=True, xlabel="PosYmu [cm]", figname="PosYMuCutOnEneEle")

        plotEleDistrWithCut(dataset=datasetEleList[0],variable="PosZmu", cutVariable="EneEle", cutCenter=500, cutRange=[100, 50, 10],nbins=100, log=True, xlabel="PosZmu [cm]", figname="PosZMuCutOnEneEle")

        plotEleDistrWithCut(dataset=datasetEleList[0],variable="CX", cutVariable="EneEle", cutCenter=500, cutRange=[500, 450,100, 50, 10],nbins=100, log=True, xlabel="cos(x)", figname="CXCutOnEneEle")

        plotEleDistrWithCut(dataset=datasetEleList[0],variable="CY", cutVariable="EneEle", cutCenter=500, cutRange=[500, 450,100, 50, 10],nbins=100, log=True, xlabel="cos(y)", figname="CYCutOnEneEle")

        plotEleDistrWithCut(dataset=datasetEleList[0],variable="CZ", cutVariable="EneEle", cutCenter=500, cutRange=[500, 450,100, 50, 10],nbins=100, log=True, xlabel="cos(z)", figname="CZCutOnEneEle")

        plotEleDistrWithCut(dataset=datasetEleList[0],variable="EneEle", cutVariable="EneEle", cutCenter=500, cutRange=[500,100, 50, 10],nbins=100, log=True, xlabel="E [GeV]", figname="EnEleCutOnEneEle")
else:
    print("Plots regarding parent electrons NOT requested")


# In[43]:


fig, ax = plt.subplots(nrows=1, ncols=len(datasetList)+1, figsize=((len(datasetList)+1)*9,8), sharey=False)
plt.suptitle("")

for i, dataset in enumerate(datasetEleList):
    ax[i].set_ylabel("x_e [m]",fontsize='14')
    ax[i].set_xlabel("z_e [m]",fontsize='14')
    ax[i].hist2d(dataset["PosZEle"]/100,dataset["PosXEle"]/100, range=([-5,100],[-5,0.3]), weights=dataset["Weight"], bins=100, norm=matplotlib.colors.LogNorm(vmin=1,vmax=1e7), cmap='Reds')
    
    PCM=ax[i].get_children()[0] #get the mappable, the 1st and the 2nd are the x and y axes
    cb=plt.colorbar(PCM, ax=ax[i]) 
    cb.set_label('Arb. Units')
    ax[i].set_title(labelList[i])
    ax[len(datasetList)].hist2d(dataset["PosZEle"]/100,dataset["PosXEle"]/100, range=([-5,100],[-5,0.3]),weights=dataset["Weight"], bins=100, norm=matplotlib.colors.LogNorm(vmin=1,vmax=1e7))
    
    ax[len(datasetList)].set_title("Comparison")
    # ax[len(datasetList)].legend()
    ax[len(datasetList)].set_xlabel("z_e [cm]",fontsize='14')
    ax[len(datasetList)].set_ylabel("x_e [cm]",fontsize='14')

    figname=runName+"Elez_evsx_e"
    pl.savefig(figname,transparent=False, facecolor='white')   


# In[44]:


fig, ax = plt.subplots(nrows=1, ncols=len(datasetList)+1, figsize=((len(datasetList)+1)*9,8), sharey=False)
plt.suptitle("")

for i, dataset in enumerate(datasetEleList):
    ax[i].set_ylabel("y_e [m]",fontsize='14')
    ax[i].set_xlabel("z_e [m]",fontsize='14')
    ax[i].hist2d(dataset["PosZEle"]/100,dataset["PosYEle"]/100, range=([-5,100],[-0.1,0.1]), weights=dataset["Weight"], bins=100, norm=matplotlib.colors.LogNorm(vmin=1,vmax=1e7), cmap='Reds')
    
    PCM=ax[i].get_children()[0] #get the mappable, the 1st and the 2nd are the x and y axes
    cb=plt.colorbar(PCM, ax=ax[i]) 
    cb.set_label('Arb. Units')
    ax[i].set_title(labelList[i])
    ax[len(datasetList)].hist2d(dataset["PosZEle"]/100,dataset["PosYEle"]/100, range=([-5,100],[-0.1,0.1]),weights=dataset["Weight"], bins=100, norm=matplotlib.colors.LogNorm(vmin=1,vmax=1e7))
    
    ax[len(datasetList)].set_title("Comparison")
    # ax[len(datasetList)].legend()
    ax[len(datasetList)].set_xlabel("z_e [cm]",fontsize='14')
    ax[len(datasetList)].set_ylabel("y_e [cm]",fontsize='14')

    figname=runName+"Elez_evsy_e"
    pl.savefig(figname,transparent=False, facecolor='white')  


# In[ ]:




