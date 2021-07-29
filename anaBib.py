#!/usr/bin/env python
# coding: utf-8

# # BIB ANALYSIS

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
from matplotlib.pyplot import pie, axis, show
from particle import PDGID
from particle import Particle
flag3TeV=False
flagCFR=True

if (not flagCFR):
    if (not flag3TeV):
        inputFile="DigFiles/Part1.5TeV.dump"
    else:
        inputFile="DigFiles/Part3TeV.dump"
else:
    inputFileA="DigFiles/Part1.5TeV.dump"
    labelA="1p5TeV"
#    inputFileB="DigFiles/Part1.5TeV.mineNotOrd"
#    inputFileB="DigFiles/Part3TeV.dump"
    inputFileB="DigFiles/Part1.5TeV_200.dump"
#    labelB="3TeV"
    labelB="1p5TeV_200"


        #inputFile=folder+"/part_new_sigma_25_nocut"
        
listChargedHadrons=[321, 311, 211, 2212, 3122, 3112, -321, -311, -211, -2212, -3122, -3112]
nbins=50
nbinsH=200
nbinsZ=100
binwidth=0.8
cols=["PDGcode", "KinE", "PX","PY","PZ","Weight","PosX","PosY","PosZ", "Time", "Elem","PosXmu","PosYmu","PosZmu","ind1","Elem2","ind2"]


# In[2]:


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
#    plt.show()

def plot1D(ax,x,plotTitle="", label="",xlabel="x",ylabel="y",log=True,col="r", weights=None,bins=None,rng=None, numPart=None):
    ax.set_xlabel(xlabel,fontsize='14')
    ax.set_ylabel(ylabel,fontsize='14')
    ax.set_title(plotTitle)
   
    if log:
        ax.set_ylim(auto=True)
    if numPart:
        ax.hist(x,log=log,histtype='step', label=label+str.format(r' N={:.2e} $\bar x$={:.2e}',numPart if numPart else 0,x.mean()), color=col, bins=bins, range=rng, weights=weights)
    else:
        ax.hist(x,log=log,histtype='step', label=label+str.format(r'$\bar x$={:.2e}',x.mean()), color=col, bins=bins, range=rng, weights=weights)


#    ax.hist(x,log=log,histtype='step', label=label+str.format(r' N={:.2e} $\bar x$={:.2e}',numPart if numPart else 0,x.mean()), color=col, bins=bins, range=rng, weights=weights)

#    ax.hist(x,log=log,histtype='step', label=str.format(r'N={:.2e} $\bar x$={:.2e}',numPart if numPart else 0,x.mean()), color=col, bins=bins, range=rng, weights=weights)


    ax.legend(loc= "best", fontsize='x-small')
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
    ax.hist(data,log=log,histtype='step', label=label+ str.format(r' $\bar x$= {:.2e}',data.mean()), color=col, bins=bins, range=rng, weights=weights)
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

    
def getMomentum(data, part, absFlag=False):
    if isinstance(part, int):
        return numpy.sqrt(data[data["PDGcode"]==part]["PX"]**2+data[data["PDGcode"]==part]["PY"]**2+data[data["PDGcode"]==part]["PZ"]**2)
    else:
        if isinstance(part, list):
            return numpy.sqrt(data[data["PDGcode"].isin(part)]["PX"]**2+data[data["PDGcode"].isin(part)]["PY"]**2+data[data["PDGcode"].isin(part)]["PZ"]**2)


        
def drawPie(var, figName, bFlag=False, title=""):
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

def scatter_histo(x, y, ax, ax_histx, ax_histy, weights=None):
    ax_histx.tick_params(axis="x", labelbottom=False)
    ax_histy.tick_params(axis="y", labelleft=False)
    ax.hist2d(x, y,weights=weights,range=[[-750, 750], [-30, 30]],bins=100,norm=matplotlib.colors.LogNorm(), cmap='Blues')
    ax_histx.hist(x,log=True, weights=weights,histtype='step',range=(-750,750),bins=100,rwidth=binwidth,color='b')
    ax_histy.hist(y,log=True, weights=weights,histtype='step',range=(-30,30),bins=100,rwidth=binwidth,color='b', orientation='horizontal')
    plt.gca().invert_yaxis()
    ax.set_xlabel('z (cm)',fontsize=14)
    ax.set_ylabel('t (nsec)',fontsize=14)

left, width = 0.1, 0.65
bottom, height = 0.1, 0.65
spacing = 0.005
rect_scatter = [left, bottom, width, height]
rect_histx = [left, bottom + height + spacing, width, 0.2]
rect_histy = [left + width + spacing, bottom, 0.2, height]


# ## Read Dataset

# In[3]:


if (not flagCFR):
    dataset = pd.read_csv(inputFile, header=None, names=cols, delim_whitespace=True)
else:
    dataset = pd.read_csv(inputFileA, header=None, names=cols, delim_whitespace=True)
    datasetB = pd.read_csv(inputFileB, header=None, names=cols, delim_whitespace=True)
#dataset.info()


# ## Let's have a look at muon decay z position

# In[4]:


bw=1
nBinZ=int(((dataset["PosZmu"]/100).max()-(dataset["PosZmu"]/100).min())/bw)
nBinZB=int(((datasetB["PosZmu"]/100).max()-(datasetB["PosZmu"]/100).min())/bw)

#fig=plt.figure(figsize=(10,10))
fig, axs = plt.subplots(nrows=3, ncols=1, figsize=(16,14))


plot1D(axs[0],dataset["PosZmu"]/100, weights=dataset["Weight"],bins=nBinZ,col='tab:blue', plotTitle="Muon Decay Z", label=labelA, xlabel='', ylabel='Arb. Units' )
plot1D(axs[0],datasetB["PosZmu"]/100, weights=datasetB["Weight"],bins=nBinZB,col='tab:orange', plotTitle="Muon Decay Z", label=labelB, xlabel='', ylabel='Arb. Units' )
axs[1].hist(dataset["PosZmu"]/100, bins=nBinZ, cumulative=1,histtype='step', label="cum "+labelA, density=True)
axs[1].hist(datasetB["PosZmu"]/100, bins=nBinZB, cumulative=1,histtype='step', label="cum "+labelB, density=True)

axs[0].axis(ymin=100, ymax=1e9)
axs[1].grid(True, which="both")
axs[1].set_title("Cumulative Function")

axs[1].legend(loc= "best", fontsize='x-small')


axs2b=axs[2].twinx()

plot1D(axs[2],dataset["PosZmu"]/100, weights=dataset["Weight"],bins=nBinZ,col='tab:blue', plotTitle="Muon Decay Z", label=labelA, xlabel='$z_{\mu \,dec}$ (m)', ylabel='Arb. Units' )
plot1D(axs[2],datasetB["PosZmu"]/100, weights=datasetB["Weight"],bins=nBinZB,col='tab:orange', plotTitle="Muon Decay Z", label=labelB, xlabel='$z_{\mu \,dec}$ (m)', ylabel='Arb. Units' )
histoCumA=axs2b.hist(dataset["PosZmu"]/100, bins=nBinZ, cumulative=1,histtype='step', label="cum "+labelA, density=True, linestyle=':')
histoCumB=axs2b.hist(datasetB["PosZmu"]/100, bins=nBinZB, cumulative=1,histtype='step', label="cum "+labelB, density=True, linestyle=':')

cut95valA=histoCumA[1][np.argmax(histoCumA[0]>0.95)]
cut95valB=histoCumB[1][np.argmax(histoCumB[0]>0.95)]

axs[1].plot([cut95valA,cut95valA],[0,0.95], linestyle=":", color="tab:blue")
axs[1].plot([cut95valB,cut95valB],[0,0.95], linestyle=":", color="tab:orange")

#axs[1].axvline(x=10, ymax=0.5)
print("95% of entries are within: ", cut95valA, cut95valB)
axs[2].axis(ymin=100, ymax=1e9)
axs2b.grid(True, which="both", linestyle=":", color="grey")

axs2b.legend(loc= "best", fontsize='x-small')

axs[1].locator_params(axis="x", nbins=20)
axs[1].locator_params(axis="y", nbins=20)

axs[1].text(cut95valA, 0.5, "95% val: "+str("{:.2f} m".format(cut95valA)), c="tab:blue")

axs[1].text(cut95valB, 0.7, "95% val: "+str("{:.2f} m".format(cut95valB)), c="tab:orange")




#plt.ylim((100, 1e9))
figname="BIB_MuDecPreCut"
pl.savefig(figname)


# ### If Needed, let's apply a z-cut

# In[5]:


dataset=dataset[dataset["PosZmu"]<2500] #in cm
datasetB=datasetB[datasetB["PosZmu"]<2500] #in cm


# ## List of found particles' IDs

# In[6]:


(lista, freq)=np.unique(dataset["PDGcode"], return_counts=True)
(listaB, freqB)=np.unique(datasetB["PDGcode"], return_counts=True)

interestingParticles=list(set(lista).union(set(listaB)))
interestingParticles.sort()


# In[7]:


entriesA=[]
entriesB=[]
for particle in interestingParticles:
    #print(particle)
    if not dataset[dataset["PDGcode"]==particle].empty:
        #tempA=dataset.PDGcode.value_counts()[particle]
        tempA=sum(dataset[dataset["PDGcode"]==particle]["Weight"])

        #print("A: ",tempA, sum(dataset[dataset["PDGcode"]==particle]["Weight"]))
    else:
        tempA=0
    if not datasetB[datasetB["PDGcode"]==particle].empty:
#        tempB=datasetB.PDGcode.value_counts()[particle]
        tempB=sum(datasetB[datasetB["PDGcode"]==particle]["Weight"])

        #print("B: ",tempB, sum(datasetB[datasetB["PDGcode"]==particle]["Weight"]))
    else:
        tempB=0
    entriesA.append(tempA)
    entriesB.append(tempB)
    #print()


# In[8]:


unknownParticle="??"

particleNamesList=[]
for particle in interestingParticles:
    #print(PDGID(particle))
    try: 
        Particle.from_pdgid(particle)
        #print(Particle.from_pdgid(particle))
        particleNamesList.append(Particle.from_pdgid(particle).name)
    except:
        #print("non trovata")
        particleNamesList.append(unknownParticle)


# In[9]:


print(interestingParticles)
print(particleNamesList)
print(entriesA)
print(entriesB)


# In[10]:


#print("Found the following particles (in descending multiplicity order)\n",particleNamesList,"\n", particleList)

#print("B: Found the following particles (in descending multiplicity order)\n",particleNamesListB,"\n", particleListB)

fig, axs = plt.subplots(nrows=3, ncols=1, figsize=(18,16))

axs[0].bar(range(len(interestingParticles)), entriesA, align='center', log=True, color="tab:blue")
axs[0].set_title('1.5TeV') 
plt.sca(axs[0])
plt.xticks(range(len(interestingParticles)), particleNamesList, size='large')
plt.xticks(rotation=90)
plt.ylabel("Occurency", size='large')
for i, v in enumerate(entriesA):
    if v!=0:
        axs[0].text(i , v, "{:.2e}".format(v), color="tab:blue")
        #fig = matplotlib.pyplot.gcf()
#fig.set_size_inches(18.5, 10.5)

axs[1].bar(range(len(interestingParticles)), entriesB, align='center', log=True, color="tab:orange")
axs[1].set_title(labelB)
plt.sca(axs[1])
#plt.xticks(range(len(dataset["PDGcode"].value_counts())), dataset["PDGcode"].value_counts().index.values, size='small')
plt.xticks(range(len(interestingParticles)), particleNamesList, size='large')
plt.xticks(rotation=90)
plt.ylabel("Occurency", size='large')
for i, v in enumerate(entriesB):
    if v!=0:
        axs[1].text(i , v, "{:.2e}".format(v), color="tab:orange")
#fig = matplotlib.pyplot.gcf()
#fig.set_size_inches(18.5, 10.5)

if True:
    axs[2].bar(range(len(interestingParticles)), entriesA, align='center', log=True, color="tab:blue", label=labelA)
    axs[2].bar(range(len(interestingParticles)), entriesB, alpha=0.5, align='center', log=True, color="tab:orange",label=labelB)
    axs[2].set_title('Comparison')
    plt.sca(axs[2])
    #plt.xticks(range(len(dataset["PDGcode"].value_counts())), dataset["PDGcode"].value_counts().index.values, size='small')
    plt.xticks(range(len(interestingParticles)), particleNamesList, size='large')
    plt.xticks(rotation=90)
    plt.ylabel("Occurency", size='large')
    plt.legend()
    
    for i, v in enumerate(entriesA):
        if v!=0:
            axs[2].text(i , v, "{:.2e}".format(v), color="tab:blue")
    for i, v in enumerate(entriesB):
        if v!=0:
            axs[2].text(i , v, "{:.2e}".format(v), color="tab:orange")
pl.savefig("BIB_ParticleDistribution")
#plt.show()


# In[11]:


if True:
    timeOffset=3600
    datasetB["Time"]=datasetB["Time"]+timeOffset


# In[12]:


data_ph=dataset[dataset["PDGcode"]==22]
data_pos=dataset[dataset["PDGcode"]==-11]
data_el=dataset[dataset["PDGcode"]==11]
data_elpos=dataset[(dataset["PDGcode"]==11)|(dataset["PDGcode"]==-11)]
data_pr=dataset[dataset["PDGcode"]==2212]
data_neu=dataset[dataset["PDGcode"]==2112]
data_mum=dataset[dataset["PDGcode"]==13]
data_mup=dataset[dataset["PDGcode"]==-13]
data_mupmum=dataset[(dataset["PDGcode"]==-13)|(dataset["PDGcode"]==13)]
data_chh=dataset[(dataset["PDGcode"]).isin(listChargedHadrons)]


# In[13]:


data_phB=datasetB[datasetB["PDGcode"]==22]
data_posB=datasetB[datasetB["PDGcode"]==-11]
data_elB=datasetB[datasetB["PDGcode"]==11]
data_elposB=datasetB[(datasetB["PDGcode"]==11)|(datasetB["PDGcode"]==-11)]
data_prB=datasetB[datasetB["PDGcode"]==2212]
data_neuB=datasetB[datasetB["PDGcode"]==2112]
data_mumB=datasetB[datasetB["PDGcode"]==13]
data_mupB=datasetB[datasetB["PDGcode"]==-13]
data_mupmumB=datasetB[(datasetB["PDGcode"]==-13)|(datasetB["PDGcode"]==13)]
data_chhB=datasetB[(datasetB["PDGcode"]).isin(listChargedHadrons)]


# In[14]:


n_ph= sum(data_ph["Weight"])
n_pos= sum(data_pos["Weight"])
n_el=sum(data_el["Weight"])
n_pr= sum(data_pr["Weight"])
n_chh= sum(data_chh["Weight"])
n_neu= sum(data_neu["Weight"])
n_mum= sum(data_mum["Weight"])
n_mup= sum(data_mup["Weight"])
n_elpos=n_el+n_pos
n_mu=n_mup+n_mum


# In[15]:


n_phB= sum(data_phB["Weight"])
n_posB= sum(data_posB["Weight"])
n_elB=sum(data_elB["Weight"])
n_prB= sum(data_prB["Weight"])
n_chhB= sum(data_chhB["Weight"])
n_neuB= sum(data_neuB["Weight"])
n_mumB= sum(data_mumB["Weight"])
n_mupB= sum(data_mupB["Weight"])
n_elposB=n_elB+n_posB
n_muB=n_mupB+n_mumB


# ## Print Particle Numbers

# In[16]:


print('N photons= ', "{:.2e}".format(n_ph), '\nn positrons= ', "{:.2e}".format(n_pos), '\nn electrons= ', "{:.2e}".format(n_el), '\nn protons= ', "{:.2e}".format(n_pr), '\nn neutrons= ', "{:.2e}".format(n_neu), '\nn ch hadr= ', "{:.2e}".format(n_chh), '\nn mum= ', "{:.2e}".format(n_mum), '\nn mup= ', "{:.2e}".format(n_mup))
print()
print('N photons= ', "{:.2e}".format(n_phB), '\nn positrons= ', "{:.2e}".format(n_posB), '\nn electrons= ', "{:.2e}".format(n_elB), '\nn protons= ', "{:.2e}".format(n_prB), '\nn neutrons= ', "{:.2e}".format(n_neu), '\nn ch hadr= ', "{:.2e}".format(n_chhB), '\nn mum= ', "{:.2e}".format(n_mumB), '\nn mup= ', "{:.2e}".format(n_mupB))


# ## Plot Energy Spectra

# ### All relevant particles

# In[17]:


fig, axs = plt.subplots(nrows=1, ncols=5, figsize=(22,5), sharey=False)
plot1D(axs[0],data_ph["KinE"], plotTitle="$\gamma$", label= labelA,xlabel="$E_{kin}$ [GeV]",ylabel="Arb. Units", bins= nbins, weights=data_ph["Weight"], numPart=n_ph)
plot1D(axs[1],data_neu["KinE"], plotTitle="n.", label= labelA,xlabel="$E_{kin}$ [GeV]",ylabel="Arb. Units", bins= nbins, weights=data_neu["Weight"], numPart=n_neu)
plot1D(axs[2],data_elpos["KinE"], plotTitle="$e^-~/~e^+$", label= labelA,xlabel="$E_{kin}$ [GeV]",ylabel="Arb. Units", bins= nbins, weights=data_elpos["Weight"], numPart=n_elpos)
plot1D(axs[3],data_chh["KinE"], plotTitle="ch. had.", label= labelA,xlabel="$E_{kin}$ [GeV]",ylabel="Arb. Units", bins= nbins, weights=data_chh["Weight"], numPart=n_chh)
plot1D(axs[4],data_mupmum["KinE"], plotTitle="$\mu^-/\mu^+$", label= labelA,xlabel="$E_{kin}$ [GeV]",ylabel="Arb. Units", bins= nbins, weights=data_mupmum["Weight"], numPart=n_mu)

plot1D(axs[0],data_phB["KinE"], plotTitle="$\gamma$", label= labelB,xlabel="$E_{kin}$ [GeV]",ylabel="Arb. Units", bins= nbins, weights=data_phB["Weight"], numPart=n_phB, col="b")
plot1D(axs[1],data_neuB["KinE"], plotTitle="n.", label= labelB,xlabel="$E_{kin}$ [GeV]",ylabel="Arb. Units", bins= nbins, weights=data_neuB["Weight"], numPart=n_neuB, col="b")
plot1D(axs[2],data_elposB["KinE"], plotTitle="$e^-~/~e^+$", label= labelB,xlabel="$E_{kin}$ [GeV]",ylabel="Arb. Units", bins= nbins, weights=data_elposB["Weight"], numPart=n_elposB, col="b")
plot1D(axs[3],data_chhB["KinE"], plotTitle="ch. had", label= labelB,xlabel="$E_{kin}$ [GeV]",ylabel="Arb. Units", bins= nbins, weights=data_chhB["Weight"], numPart=n_chhB, col="b")
plot1D(axs[4],data_mupmumB["KinE"], plotTitle="$\mu^-/\mu^+$", label= labelB,xlabel="$E_{kin}$ [GeV]",ylabel="Arb. Units", bins= nbins, weights=data_mupmumB["Weight"], numPart=n_muB, col="b")

fig.subplots_adjust(left = 0.05,right = 0.99,wspace = 0.5, hspace = 0.5, top=0.9, bottom= 0.2)

figname="BIB_EnergySpectra"
pl.savefig(figname)


# ### Photons and e+/e-

# In[18]:


fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(24,8), sharey=False)

plot1D(ax[0], getMomentum(dataset, 22), weights=data_ph["Weight"], bins=200, label="$\gamma$", xlabel='p (GeV/c)', ylabel='Arb. Units', numPart=n_ph)
plot1D(ax[0], getMomentum(dataset, [-11,11]), weights=data_elpos["Weight"], bins=200, plotTitle="E.M. Components "+labelA, label="$e^-~/~e^+$", xlabel='p (GeV/c)', ylabel='Arb. Units', col="black", numPart=n_elpos)

plot1D(ax[1], getMomentum(datasetB, 22), weights=data_phB["Weight"], bins=200, label="$\gamma$", xlabel='p (GeV/c)', ylabel='Arb. Units', numPart=n_phB, col="magenta")
plot1D(ax[1], getMomentum(datasetB, [-11,11]), weights=data_elposB["Weight"], bins=200, plotTitle="E.M. Components "+labelB, label="$e^-~/~e^+$", xlabel='p (GeV/c)', ylabel='Arb. Units', col="gray", numPart=n_elposB)

plot1D(ax[2], getMomentum(dataset, 22), weights=data_ph["Weight"], bins=200, label="$\gamma$ "+labelA, xlabel='p (GeV/c)', ylabel='Arb. Units', numPart=n_ph)
plot1D(ax[2], getMomentum(dataset, [-11,11]), weights=data_elpos["Weight"], bins=200, plotTitle="E.M. Components", label="$e^-~/~e^+$ 1.5TeV", xlabel='p (GeV/c)', ylabel='Arb. Units', col="black", numPart=n_elpos)
plot1D(ax[2], getMomentum(datasetB, 22), weights=data_phB["Weight"], bins=200, label="$\gamma$ "+labelB, xlabel='p (GeV/c)', ylabel='Arb. Units', numPart=n_phB, col="magenta")
plot1D(ax[2], getMomentum(datasetB, [-11,11]), weights=data_elposB["Weight"], bins=200, plotTitle="E.M. Components CFR", label="$e^-~/~e^+$ 1.5TeV", xlabel='p (GeV/c)', ylabel='Arb. Units', col="gray", numPart=n_elposB)

#plt.ylim((1e3,2e9))
#plt.xlim((0,0.2))
ax[0].axis(ymin=1e1, ymax=3e7, xmin=0, xmax=0.45)
ax[1].axis(ymin=1e1, ymax=3e7, xmin=0, xmax=0.45)
ax[2].axis(ymin=1e1, ymax=3e7, xmin=0, xmax=0.45)
plt.legend()
figname="BIB_EneEGamma"
pl.savefig(figname)


# ### Hadrons

# In[19]:


fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(24,8), sharey=False)

plot1D(ax[0], getMomentum(dataset, listChargedHadrons), weights=data_chh["Weight"], bins=nbins, rng=(0,1), label="ch. had", xlabel='p (GeV/c)', ylabel='Arb. Units', col="m", numPart=n_pr)
plot1D(ax[0], getMomentum(dataset, [2112]), weights=data_neu["Weight"], bins=nbinsH, rng=(0,1), plotTitle="Hadrons "+labelA, label="neutrons", xlabel='p (GeV/c)', ylabel='Arb. Units', col="tab:blue",numPart=n_neu)

plot1D(ax[1], getMomentum(datasetB, listChargedHadrons), weights=data_chhB["Weight"], bins=nbins, rng=(0,1), label="ch. had", xlabel='p (GeV/c)', ylabel='Arb. Units', col="hotpink", numPart=n_prB)
plot1D(ax[1], getMomentum(datasetB, [2112]), weights=data_neuB["Weight"], bins=nbinsH, rng=(0,1), plotTitle="Hadrons "+labelB, label="neutrons", xlabel='p (GeV/c)', ylabel='Arb. Units', col="tab:cyan",numPart=n_neuB)

plot1D(ax[2], getMomentum(dataset, listChargedHadrons), weights=data_chh["Weight"], bins=nbins, rng=(0,1), label="ch. had "+labelA, xlabel='p (GeV/c)', ylabel='Arb. Units', col="m", numPart=n_pr)
plot1D(ax[2], getMomentum(dataset, [2112]), weights=data_neu["Weight"], bins=nbinsH, rng=(0,1), label="n "+labelA, xlabel='p (GeV/c)', ylabel='Arb. Units', col="tab:blue",numPart=n_neu)
plot1D(ax[2], getMomentum(datasetB, listChargedHadrons), weights=data_chhB["Weight"], bins=nbins, rng=(0,1), label="ch. had "+labelB, xlabel='p (GeV/c)', ylabel='Arb. Units', col="hotpink", numPart=n_prB)
plot1D(ax[2], getMomentum(datasetB, [2112]), weights=data_neuB["Weight"], bins=nbinsH, rng=(0,1), plotTitle="Hadrons CFR", label="n "+labelB, xlabel='p (GeV/c)', ylabel='Arb. Units', col="tab:cyan",numPart=n_neuB)

ax[0].axis(ymin=1e1, ymax=3e7)
ax[1].axis(ymin=1e1, ymax=3e7)
ax[2].axis(ymin=1e1, ymax=3e7)

#plt.ylim((1e2,2e7))
plt.xlim((0,1))
plt.legend()
figname="BIB_EneHadrons"
pl.savefig(figname)


# ## Plot Time Distributions

# In[20]:


fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(24,8), sharey=False)

plot1D(ax[0], data_ph["Time"], weights=data_ph["Weight"],rng=(-30,100), bins=nbinsH, col="r", label="Photon", xlabel="T [ns]", ylabel="Arb. Units")
plot1D(ax[0], data_elpos["Time"], weights=data_elpos["Weight"],rng=(-30,100), bins=nbinsH, col="black", label="e+e-", xlabel="T [ns]", ylabel="Arb. Units")
plot1D(ax[0], data_chh["Time"], weights=data_chh["Weight"],rng=(-30,100), bins=nbinsH, col="m", label="ch. had", xlabel="T [ns]", ylabel="Arb. Units")
plot1D(ax[0], data_neu["Time"], weights=data_neu["Weight"],rng=(-30,100), bins=nbinsH, col="blue", label="neutron", xlabel="T [ns]", ylabel="Arb. Units")
plot1D(ax[0], data_mupmum["Time"], weights=data_mupmum["Weight"],rng=(-30,100), bins=nbinsH, col="chartreuse", label="mu+ mu-", xlabel="T [ns]", ylabel="Arb. Units", plotTitle="Time Distribution "+labelA)

plot1D(ax[1], data_phB["Time"], weights=data_phB["Weight"],rng=(-30,100), bins=nbinsH, col="salmon", label="Photon", xlabel="T [ns]", ylabel="Arb. Units")
plot1D(ax[1], data_elposB["Time"], weights=data_elposB["Weight"],rng=(-30,100), bins=nbinsH, col="dimgrey", label="e+e-", xlabel="T [ns]", ylabel="Arb. Units")
plot1D(ax[1], data_chhB["Time"], weights=data_chhB["Weight"],rng=(-30,100), bins=nbinsH, col="hotpink", label="ch. had", xlabel="T [ns]", ylabel="Arb. Units")
plot1D(ax[1], data_neuB["Time"], weights=data_neuB["Weight"],rng=(-30,100), bins=nbinsH, col="cyan", label="neutron", xlabel="T [ns]", ylabel="Arb. Units")
plot1D(ax[1], data_mupmumB["Time"], weights=data_mupmumB["Weight"],rng=(-30,100), bins=nbinsH, col="springgreen", label="mu+ mu-", xlabel="T [ns]", ylabel="Arb. Units", plotTitle="Time Distribution "+labelB)

plot1D(ax[2], data_ph["Time"], weights=data_ph["Weight"],rng=(-30,100), bins=nbinsH, col="r", label="Photon", xlabel="T [ns]", ylabel="Arb. Units")
plot1D(ax[2], data_elpos["Time"], weights=data_elpos["Weight"],rng=(-30,100), bins=nbinsH, col="black", label="e+e-", xlabel="T [ns]", ylabel="Arb. Units")
plot1D(ax[2], data_chh["Time"], weights=data_chh["Weight"],rng=(-30,100), bins=nbinsH, col="m", label="ch. had", xlabel="T [ns]", ylabel="Arb. Units")
plot1D(ax[2], data_neu["Time"], weights=data_neu["Weight"],rng=(-30,100), bins=nbinsH, col="blue", label="neutron", xlabel="T [ns]", ylabel="Arb. Units")
plot1D(ax[2], data_mupmum["Time"], weights=data_mupmum["Weight"],rng=(-30,100), bins=nbinsH, col="chartreuse", label="mu+ mu-", xlabel="T [ns]", ylabel="Arb. Units", plotTitle="Time Distribution "+labelA)

plot1D(ax[2], data_phB["Time"], weights=data_phB["Weight"],rng=(-30,100), bins=nbinsH, col="salmon", label="Photon", xlabel="T [ns]", ylabel="Arb. Units")
plot1D(ax[2], data_elposB["Time"], weights=data_elposB["Weight"],rng=(-30,100), bins=nbinsH, col="dimgrey", label="e+e-", xlabel="T [ns]", ylabel="Arb. Units")
plot1D(ax[2], data_chhB["Time"], weights=data_chhB["Weight"],rng=(-30,100), bins=nbinsH, col="hotpink", label="ch. had", xlabel="T [ns]", ylabel="Arb. Units")
plot1D(ax[2], data_neuB["Time"], weights=data_neuB["Weight"],rng=(-30,100), bins=nbinsH, col="cyan", label="neutron", xlabel="T [ns]", ylabel="Arb. Units")
plot1D(ax[2], data_mupmumB["Time"], weights=data_mupmumB["Weight"],rng=(-30,100), bins=nbinsH, col="springgreen", label="mu+ mu-", xlabel="T [ns]", ylabel="Arb. Units", plotTitle="Time Distribution CFR")



ax[0].axis(ymin=1e1, ymax=3e7)
ax[1].axis(ymin=1e1, ymax=3e7)
ax[2].axis(ymin=1e1, ymax=3e7)



#plt.xlabel('t (nsec)',fontsize=14)
#plt.ylabel('Arb. Units',fontsize=14)
#plt.ylim((100,2e9))
# plt.xlim((-30,100))
plt.legend()
figname="BIB_Time"
pl.savefig(figname)


# ## Plot Pie Charts

# In[21]:


drawPie("Elem", "BIB_PieDet", title=labelA)
drawPie("Elem", "BIB_PieDet", bFlag=True, title=labelB)


# In[22]:


drawPie("Elem2", "BIB_PieFirstInt", title=labelA)
drawPie("Elem2", "BIB_PieFirstInt", title=labelB, bFlag=True)


# ## Plot Muons' Decay Position

# ### Global

# In[23]:


fig=plt.figure(figsize=(6,5))
plot1D(fig.gca(),dataset["PosZmu"]/100, weights=dataset["Weight"],bins=nbinsH,col='m', plotTitle="Muon Decay Z", label=labelA, xlabel='$z_{\mu \,dec}$ (m)', ylabel='Arb. Units' )

plot1D(fig.gca(),datasetB["PosZmu"]/100, weights=datasetB["Weight"],bins=nbinsH,col='hotpink', plotTitle="Muon Decay Z", label=labelB, xlabel='$z_{\mu \,dec}$ (m)', ylabel='Arb. Units' )


plt.ylim((100, 1e9))
figname="BIB_MuDec"
pl.savefig(figname)


# ### Per Particle

# In[24]:


fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(21,5), sharey=False)


plot1D(ax[0], dataset[dataset["PDGcode"]==22]["PosZmu"]/100,log=True, weights=dataset[dataset["PDGcode"]==22]["Weight"],bins=nbinsZ,rng=(-1,25),col='r', xlabel= "$z_{\mu \,dec}$ (m)", ylabel="Arb. Units", label='$\gamma$ '+labelA)
plot1D(ax[0], dataset[abs(dataset["PDGcode"])==11]["PosZmu"]/100,log=True, weights=dataset[abs(dataset["PDGcode"])==11]["Weight"],bins=nbinsZ,rng=(-1,25),col='black', xlabel= "$z_{\mu \,dec}$ (m)", ylabel="Arb. Units", label= 'e+e- '+labelA)
plot1D(ax[0], dataset[abs(dataset["PDGcode"])==2212]["PosZmu"]/100,log=True, weights=dataset[abs(dataset["PDGcode"])==2212]["Weight"],bins=nbinsZ,rng=(-1,25),col='m', xlabel= "$z_{\mu \,dec}$ (m)", ylabel="Arb. Units", label= 'ch. had '+labelA)
plot1D(ax[0], dataset[dataset["PDGcode"]==2112]["PosZmu"]/100,log=True, weights=dataset[dataset["PDGcode"]==2112]["Weight"],bins=nbinsZ,rng=(-1,25),col='blue', xlabel= "$z_{\mu \,dec}$ (m)", ylabel="Arb. Units", label= 'n '+labelA)
plot1D(ax[0], dataset[abs(dataset["PDGcode"])==13]["PosZmu"]/100,log=True, weights=dataset[abs(dataset["PDGcode"])==13]["Weight"],bins=nbinsZ,rng=(-1,25),col='chartreuse', xlabel= "$z_{\mu \,dec}$ (m)", ylabel="Arb. Units", plotTitle="Mu decay point", label= 'mu+mu- '+labelA)

plot1D(ax[1], datasetB[datasetB["PDGcode"]==22]["PosZmu"]/100,log=True, weights=datasetB[datasetB["PDGcode"]==22]["Weight"],bins=nbinsZ,rng=(-1,25),col='salmon', xlabel= "$z_{\mu \,dec}$ (m)", ylabel="Arb. Units", label='$\gamma$ '+labelB)
plot1D(ax[1], datasetB[abs(datasetB["PDGcode"])==11]["PosZmu"]/100,log=True, weights=datasetB[abs(datasetB["PDGcode"])==11]["Weight"],bins=nbinsZ,rng=(-1,25),col='dimgrey', xlabel= "$z_{\mu \,dec}$ (m)", ylabel="Arb. Units", label= 'e+e- '+labelB)
plot1D(ax[1], datasetB[abs(datasetB["PDGcode"])==2212]["PosZmu"]/100,log=True, weights=datasetB[abs(datasetB["PDGcode"])==2212]["Weight"],bins=nbinsZ,rng=(-1,25),col='hotpink', xlabel= "$z_{\mu \,dec}$ (m)", ylabel="Arb. Units", label= 'ch. had '+labelB)
plot1D(ax[1], datasetB[datasetB["PDGcode"]==2112]["PosZmu"]/100,log=True, weights=datasetB[datasetB["PDGcode"]==2112]["Weight"],bins=nbinsZ,rng=(-1,25),col='cyan', xlabel= "$z_{\mu \,dec}$ (m)", ylabel="Arb. Units", label= 'n '+labelB)
plot1D(ax[1], datasetB[abs(datasetB["PDGcode"])==13]["PosZmu"]/100,log=True, weights=datasetB[abs(datasetB["PDGcode"])==13]["Weight"],bins=nbinsZ,rng=(-1,25),col='springgreen', xlabel= "$z_{\mu \,dec}$ (m)", ylabel="Arb. Units", plotTitle="Mu decay point", label= 'mu+mu- '+labelB)


plot1D(ax[2], dataset[dataset["PDGcode"]==22]["PosZmu"]/100,log=True, weights=dataset[dataset["PDGcode"]==22]["Weight"],bins=nbinsZ,rng=(-1,25),col='r', xlabel= "$z_{\mu \,dec}$ (m)", ylabel="Arb. Units", label='$\gamma$ '+labelA)
plot1D(ax[2], dataset[abs(dataset["PDGcode"])==11]["PosZmu"]/100,log=True, weights=dataset[abs(dataset["PDGcode"])==11]["Weight"],bins=nbinsZ,rng=(-1,25),col='black', xlabel= "$z_{\mu \,dec}$ (m)", ylabel="Arb. Units", label= 'e+e- '+labelA)
plot1D(ax[2], dataset[abs(dataset["PDGcode"])==2212]["PosZmu"]/100,log=True, weights=dataset[abs(dataset["PDGcode"])==2212]["Weight"],bins=nbinsZ,rng=(-1,25),col='m', xlabel= "$z_{\mu \,dec}$ (m)", ylabel="Arb. Units", label= 'ch. had '+labelA)
plot1D(ax[2], dataset[dataset["PDGcode"]==2112]["PosZmu"]/100,log=True, weights=dataset[dataset["PDGcode"]==2112]["Weight"],bins=nbinsZ,rng=(-1,25),col='blue', xlabel= "$z_{\mu \,dec}$ (m)", ylabel="Arb. Units", label= 'n '+labelA)
plot1D(ax[2], dataset[abs(dataset["PDGcode"])==13]["PosZmu"]/100,log=True, weights=dataset[abs(dataset["PDGcode"])==13]["Weight"],bins=nbinsZ,rng=(-1,25),col='chartreuse', xlabel= "$z_{\mu \,dec}$ (m)", ylabel="Arb. Units", plotTitle="Mu decay point", label= 'mu+mu- '+labelA)

plot1D(ax[2], datasetB[datasetB["PDGcode"]==22]["PosZmu"]/100,log=True, weights=datasetB[datasetB["PDGcode"]==22]["Weight"],bins=nbinsZ,rng=(-1,25),col='salmon', xlabel= "$z_{\mu \,dec}$ (m)", ylabel="Arb. Units", label='$\gamma$ '+labelB)
plot1D(ax[2], datasetB[abs(datasetB["PDGcode"])==11]["PosZmu"]/100,log=True, weights=datasetB[abs(datasetB["PDGcode"])==11]["Weight"],bins=nbinsZ,rng=(-1,25),col='dimgrey', xlabel= "$z_{\mu \,dec}$ (m)", ylabel="Arb. Units", label= 'e+e- '+labelB)
plot1D(ax[2], datasetB[abs(datasetB["PDGcode"])==2212]["PosZmu"]/100,log=True, weights=datasetB[abs(datasetB["PDGcode"])==2212]["Weight"],bins=nbinsZ,rng=(-1,25),col='hotpink', xlabel= "$z_{\mu \,dec}$ (m)", ylabel="Arb. Units", label= 'ch. had '+labelB)
plot1D(ax[2], datasetB[datasetB["PDGcode"]==2112]["PosZmu"]/100,log=True, weights=datasetB[datasetB["PDGcode"]==2112]["Weight"],bins=nbinsZ,rng=(-1,25),col='cyan', xlabel= "$z_{\mu \,dec}$ (m)", ylabel="Arb. Units", label= 'n '+labelB)
plot1D(ax[2], datasetB[abs(datasetB["PDGcode"])==13]["PosZmu"]/100,log=True, weights=datasetB[abs(datasetB["PDGcode"])==13]["Weight"],bins=nbinsZ,rng=(-1,25),col='springgreen', xlabel= "$z_{\mu \,dec}$ (m)", ylabel="Arb. Units", plotTitle="Mu decay point", label= 'mu+mu- '+labelB)




plt.xlabel('$z_{\mu \,dec}$ (m)',fontsize=14)
plt.ylabel('Arb. Units',fontsize=14)
plt.ylim((10.,2e8))
plt.legend()
figname="BIB_Mudec2"
pl.savefig(figname)


# In[25]:


fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(9,8), sharey=False)

ax[0].set_title('FLUKA '+ labelA,fontsize=22)
ax[0].hist2d(dataset["PosZ"],dataset["PosX"],norm=matplotlib.colors.LogNorm(),bins=500, cmap='plasma')

ax[1].set_title('FLUKA '+labelB,fontsize=22)

ax[1].hist2d(datasetB["PosZ"],datasetB["PosX"],norm=matplotlib.colors.LogNorm(),bins=500, cmap='plasma')

plt.xlabel('z (cm)',fontsize=14)
plt.ylabel('x (cm)',fontsize=14)
ax[0].axis(ymin=-250, ymax=250)
ax[1].axis(ymin=-250, ymax=250)

#ax[1].ylim((-250,250))

#plt.axes().set_aspect('equal')
#cb=plt.colorbar(shrink=0.8)
#cb.set_label('Arb. Units')
#plt.tight_layout(pad=.5)
figname="BIB_ZvsX_FLUKA"
pl.savefig(figname)


# ### Scatter Plots

# In[26]:


fig = plt.figure(figsize=(8, 8))
fig.suptitle("photons")
ax = fig.add_axes(rect_scatter)
ax_histx = fig.add_axes(rect_histx, sharex=ax)
ax_histy = fig.add_axes(rect_histy, sharey=ax)

scatter_histo(dataset[abs(dataset["PDGcode"])==22]["PosZ"], dataset[abs(dataset["PDGcode"])==22]["Time"], ax, ax_histx, ax_histy, weights=dataset[abs(dataset["PDGcode"])==22]["Weight"])
figname="BIB_decph"
pl.savefig(figname)


# In[27]:


fig = plt.figure(figsize=(8, 8))
fig.suptitle("photons")
ax = fig.add_axes(rect_scatter)
ax_histx = fig.add_axes(rect_histx, sharex=ax)
ax_histy = fig.add_axes(rect_histy, sharey=ax)

scatter_histo(datasetB[abs(datasetB["PDGcode"])==22]["PosZ"], datasetB[abs(datasetB["PDGcode"])==22]["Time"], ax, ax_histx, ax_histy, weights=datasetB[abs(datasetB["PDGcode"])==22]["Weight"])
figname="BIB_decphB"
pl.savefig(figname)


# In[28]:


fig = plt.figure(figsize=(8, 8))
fig.suptitle("$e^+$ $e^-$")
ax = fig.add_axes(rect_scatter)
ax_histx = fig.add_axes(rect_histx, sharex=ax)
ax_histy = fig.add_axes(rect_histy, sharey=ax)

scatter_histo(dataset[abs(dataset["PDGcode"])==11]["PosZ"], dataset[abs(dataset["PDGcode"])==11]["Time"], ax, ax_histx, ax_histy, weights=dataset[abs(dataset["PDGcode"])==11]["Weight"])
figname="BIB_decelpos"
pl.savefig(figname)


# In[29]:


fig = plt.figure(figsize=(8, 8))
fig.suptitle("$e^+$ $e^-$")
ax = fig.add_axes(rect_scatter)
ax_histx = fig.add_axes(rect_histx, sharex=ax)
ax_histy = fig.add_axes(rect_histy, sharey=ax)

scatter_histo(datasetB[abs(datasetB["PDGcode"])==11]["PosZ"], datasetB[abs(datasetB["PDGcode"])==11]["Time"], ax, ax_histx, ax_histy, weights=datasetB[abs(datasetB["PDGcode"])==11]["Weight"])
figname="BIB_decelposB"
pl.savefig(figname)


# In[30]:


fig = plt.figure(figsize=(8, 8))
fig.suptitle("charged hadrons")
ax = fig.add_axes(rect_scatter)
ax_histx = fig.add_axes(rect_histx, sharex=ax)
ax_histy = fig.add_axes(rect_histy, sharey=ax)

scatter_histo(dataset[abs(dataset["PDGcode"])==2112]["PosZ"], dataset[abs(dataset["PDGcode"])==2112]["Time"], ax, ax_histx, ax_histy, weights=dataset[abs(dataset["PDGcode"])==2112]["Weight"])
figname="BIB_decneu"
pl.savefig(figname)


# In[31]:


fig = plt.figure(figsize=(8, 8))
fig.suptitle("charged hadrons")
ax = fig.add_axes(rect_scatter)
ax_histx = fig.add_axes(rect_histx, sharex=ax)
ax_histy = fig.add_axes(rect_histy, sharey=ax)

scatter_histo(datasetB[abs(datasetB["PDGcode"])==2112]["PosZ"], datasetB[abs(datasetB["PDGcode"])==2112]["Time"], ax, ax_histx, ax_histy, weights=datasetB[abs(datasetB["PDGcode"])==2112]["Weight"])
figname="BIB_decneuB"
pl.savefig(figname)


# In[ ]:





# In[ ]:





# In[ ]:




