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

flag3TeV=True

if (not flag3TeV):
    inputFile="DigFiles/1.5TeV"
else:
    inputFile="DigFiles/Part3TeV.dump"
#inputFile=folder+"/part_new_sigma_25_nocut"

nbins=50
nbinsH=200
nbinsZ=100
binwidth=0.8
cols=["PDGcode", "KinE", "PX","PY","PZ","Weight","PosX","PosY","PosZ", "Time", "Elem","PosXmu","PosYmu","PosZmu","ind1","Elem2","ind2"]


# ## Read Dataset

# In[2]:


dataset = pd.read_csv(inputFile, header=None, names=cols, delim_whitespace=True)
dataset.info()


# In[3]:


n_ph= sum(dataset[dataset["PDGcode"]==22]["Weight"])
n_pos= sum(dataset[dataset["PDGcode"]==-11]["Weight"])
n_el=sum(dataset[dataset["PDGcode"]==11]["Weight"])
n_pr= sum(dataset[abs(dataset["PDGcode"])==2212]["Weight"])
n_neu= sum(dataset[dataset["PDGcode"]==2112]["Weight"])
n_mum= sum(dataset[dataset["PDGcode"]==13]["Weight"])
n_mup= sum(dataset[dataset["PDGcode"]==-13]["Weight"])
n_elpos=n_el+n_pos
n_mu=n_mup+n_mum


# ## Print Particle Numbers

# In[4]:


print('N photons= ', "{:.2e}".format(n_ph), '\nn positrons= ', "{:.2e}".format(n_pos), '\nn electrons= ', "{:.2e}".format(n_el), '\nn protons= ', "{:.2e}".format(n_pr), '\nn neutrons= ', "{:.2e}".format(n_neu), '\nn mum= ', "{:.2e}".format(n_mum), '\nn mup= ', "{:.2e}".format(n_mup))
   
   


# In[5]:


def plot1D(ax,x,label="",xlabel="x",ylabel="y",log=False,col="r", weights=None,bins=None,rng=None, numPart=None):
    ax.set_xlabel(xlabel,fontsize='14')
    ax.set_ylabel(ylabel,fontsize='14')
    if numPart:
        if numPart==0:
            ax.set_title('NO PARTS ')
        else:
            ax.set_title('N of particles '+ str.format('{:.2e}',numPart),fontsize=12)
    if log:
        ax.set_ylim(auto=True)
    ax.hist(x,log=log,histtype='step', label=label, color=col, bins=bins, range=rng, weights=weights)
    ax.legend(loc= "best")
    ax.grid(True)
    return ax

def getMomentum(part, absFlag=False):
    if abs:
        return numpy.sqrt(dataset[abs(dataset["PDGcode"])==part]["PX"]**2+dataset[abs(dataset["PDGcode"])==part]["PY"]**2+dataset[abs(dataset["PDGcode"])==part]["PZ"]**2)
    else:
        return numpy.sqrt(dataset[dataset["PDGcode"]==part]["PX"]**2+dataset[dataset["PDGcode"]==part]["PY"]**2+dataset[dataset["PDGcode"]==part]["PZ"]**2)

def drawPie(var, figName):
    fig=plt.figure(figsize=(5,5))
    sums = dataset.groupby(dataset[var])["Weight"].sum()
    axis('equal');
    cmap = plt.get_cmap('Spectral')
    colors = [cmap(i) for i in np.linspace(0, 1, 8)]
    pie(sums, autopct='%1.1f%%',labels=sums.index,colors=colors)
    figname=figName
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


# ## Plot Energy Spectra

# ### All relevant particles

# In[25]:


fig, axs = plt.subplots(nrows=1, ncols=5, figsize=(15,3), sharey=False)
plot1D(axs[0],dataset[dataset["PDGcode"]==22]["KinE"], "Photons","$E_{kin}$ [GeV]","Arb. Units",log=True, bins= nbins, weights=dataset[dataset["PDGcode"]==22]["Weight"], numPart=n_ph)
plot1D(axs[1],dataset[dataset["PDGcode"]==2112]["KinE"], "Neutrons","$E_{kin}$ [GeV]","Arb. Units",log=True, bins= nbins, weights=dataset[dataset["PDGcode"]==2112]["Weight"], numPart=n_neu)
plot1D(axs[2],dataset[abs(dataset["PDGcode"])==11]["KinE"], "$e^-~/~e^+$","$E_{kin}$ [GeV]","Arb. Units",log=True, bins= nbins, weights=dataset[abs(dataset["PDGcode"])==11]["Weight"], numPart=n_elpos)
plot1D(axs[3],dataset[abs(dataset["PDGcode"])==2212]["KinE"], "Ch. Had","$E_{kin}$ [GeV]","Arb. Units",log=True, bins= nbins, weights=dataset[abs(dataset["PDGcode"])==2212]["Weight"], numPart=n_pr)
plot1D(axs[4],dataset[abs(dataset["PDGcode"])==13]["KinE"], "$\mu^-/\mu^+$","$E_{kin}$ [GeV]","Arb. Units",log=True, bins= nbins, weights=dataset[abs(dataset["PDGcode"])==13]["Weight"], numPart=n_mu)

fig.subplots_adjust(left = 0.05,right = 0.99,wspace = 0.5, hspace = 0.5, top=0.9, bottom= 0.2)

figname="BIB_EnergySpectra"
pl.savefig(figname)


# ### Photons and e+/e-

# In[26]:


fig=plt.figure(figsize=(8,8))
#plot1D(fig.gca(), getMomentum(22), log=True, weights=dataset[dataset["PDGcode"]==22]["Weight"], bins=200, rng=(0,0.2), label="$\gamma$", xlabel='p (GeV/c)', ylabel='Arb. Units')
#plot1D(fig.gca(), getMomentum(11, absFlag=True), log=True, weights=dataset[abs(dataset["PDGcode"])==11]["Weight"], bins=200, rng=(0,0.2), label="$e^-~/~e^+$", xlabel='p (GeV/c)', ylabel='Arb. Units', col="black")

plot1D(fig.gca(), getMomentum(22), log=True, weights=dataset[dataset["PDGcode"]==22]["Weight"], bins=200, label="$\gamma$", xlabel='p (GeV/c)', ylabel='Arb. Units')
plot1D(fig.gca(), getMomentum(11, absFlag=True), log=True, weights=dataset[abs(dataset["PDGcode"])==11]["Weight"], bins=200, label="$e^-~/~e^+$", xlabel='p (GeV/c)', ylabel='Arb. Units', col="black")

#plt.ylim((1e3,2e9))
#plt.xlim((0,0.2))
plt.legend()
figname="BIB_EneEGamma"
pl.savefig(figname)


# ### Hadrons

# In[27]:


fig=plt.figure(figsize=(8,8))
plot1D(fig.gca(), getMomentum(2212, absFlag=True), log=True, weights=dataset[abs(dataset["PDGcode"])==2212]["Weight"], bins=nbinsH, rng=(0,1), label="ch. had", xlabel='p (GeV/c)', ylabel='Arb. Units', col="m")
plot1D(fig.gca(), getMomentum(2112), log=True, weights=dataset[dataset["PDGcode"]==2112]["Weight"], bins=nbinsH, rng=(0,1), label="neutrons", xlabel='p (GeV/c)', ylabel='Arb. Units', col="blue")

plt.ylim((1e2,2e7))
plt.xlim((0,1))
plt.legend()
figname="BIB_EneHadrons"
pl.savefig(figname)


# ## Plot Time Distributions

# In[9]:


fig=plt.figure(figsize=(8,8))
plot1D(fig.gca(), dataset[dataset["PDGcode"]==22]["Time"], log=True, weights=dataset[dataset["PDGcode"]==22]["Weight"],rng=(-30,100), bins=nbinsH, col="r", label="Photon")
plot1D(fig.gca(), dataset[abs(dataset["PDGcode"])==11]["Time"], log=True, weights=dataset[abs(dataset["PDGcode"])==11]["Weight"],rng=(-30,100), bins=nbinsH, col="black", label="e+e-")
plot1D(fig.gca(), dataset[abs(dataset["PDGcode"])==2212]["Time"], log=True, weights=dataset[abs(dataset["PDGcode"])==2212]["Weight"],rng=(-30,100), bins=nbinsH, col="m", label="ch. had")
plot1D(fig.gca(), dataset[dataset["PDGcode"]==2112]["Time"], log=True, weights=dataset[dataset["PDGcode"]==2112]["Weight"],rng=(-30,100), bins=nbinsH, col="blue", label="neutron")
plot1D(fig.gca(), dataset[abs(dataset["PDGcode"])==13]["Time"], log=True, weights=dataset[abs(dataset["PDGcode"])==13]["Weight"],rng=(-30,100), bins=nbinsH, col="chartreuse", label="mu+ mu-")

plt.xlabel('t (nsec)',fontsize=14)
plt.ylabel('Arb. Units',fontsize=14)
plt.ylim((100,2e9))
# plt.xlim((-30,100))
plt.legend()
figname="BIB_Time"
pl.savefig(figname)


# ## Plot Pie Charts

# In[11]:


drawPie("Elem", "BIB_PieDet")


# In[12]:


drawPie("Elem2", "BIB_PieFirstInt")


# ## Plot Muons' Decay Position

# ### Global

# In[13]:


fig=plt.figure(figsize=(6,5))
plot1D(fig.gca(),dataset["PosZmu"]/100,log=True, weights=dataset["Weight"],bins=nbinsH,col='m', label="Muon Decay Z", xlabel='$z_{\mu \,dec}$ (m)', ylabel='Arb. Units' )

plt.ylim((100, 1e9))
figname="BIB_MuDec"
pl.savefig(figname)


# ### Per Particle

# In[14]:


fig=plt.figure(figsize=(6,5))
plot1D(fig.gca(), dataset[dataset["PDGcode"]==22]["PosZmu"]/100,log=True, weights=dataset[dataset["PDGcode"]==22]["Weight"],bins=nbinsZ,rng=(-1,25),col='r', label= 'photon')
plot1D(fig.gca(), dataset[abs(dataset["PDGcode"])==11]["PosZmu"]/100,log=True, weights=dataset[abs(dataset["PDGcode"])==11]["Weight"],bins=nbinsZ,rng=(-1,25),col='black', label= 'e+e-')
plot1D(fig.gca(), dataset[abs(dataset["PDGcode"])==2212]["PosZmu"]/100,log=True, weights=dataset[abs(dataset["PDGcode"])==2212]["Weight"],bins=nbinsZ,rng=(-1,25),col='m', label= 'ch. had')
plot1D(fig.gca(), dataset[dataset["PDGcode"]==2112]["PosZmu"]/100,log=True, weights=dataset[dataset["PDGcode"]==2112]["Weight"],bins=nbinsZ,rng=(-1,25),col='blue', label= 'neutron')
plot1D(fig.gca(), dataset[abs(dataset["PDGcode"])==13]["PosZmu"]/100,log=True, weights=dataset[abs(dataset["PDGcode"])==13]["Weight"],bins=nbinsZ,rng=(-1,25),col='chartreuse', label= 'mu+mu-')

plt.xlabel('$z_{\mu \,dec}$ (m)',fontsize=14)
plt.ylabel('Arb. Units',fontsize=14)
plt.ylim((10.,2e8))
plt.legend()
figname="BIB_Mudec2"
pl.savefig(figname)


# In[15]:


fig=plt.figure(figsize=(9,3))
plt.gca().set_title('FLUKA',fontsize=22)
plt.hist2d(dataset["PosZ"],dataset["PosX"],norm=matplotlib.colors.LogNorm(),bins=500, cmap='plasma')
plt.xlabel('z (cm)',fontsize=14)
plt.ylabel('x (cm)',fontsize=14)
plt.ylim((-250,250))
#plt.axes().set_aspect('equal')
cb=plt.colorbar(shrink=0.8)
cb.set_label('Arb. Units')
plt.tight_layout(pad=.5)
figname="BIB_ZvsX_FLUKA"
pl.savefig(figname)


# ### Scatter Plots

# In[22]:


fig = plt.figure(figsize=(8, 8))
fig.suptitle("photons")
ax = fig.add_axes(rect_scatter)
ax_histx = fig.add_axes(rect_histx, sharex=ax)
ax_histy = fig.add_axes(rect_histy, sharey=ax)

scatter_histo(dataset[abs(dataset["PDGcode"])==22]["PosZ"], dataset[abs(dataset["PDGcode"])==22]["Time"], ax, ax_histx, ax_histy, weights=dataset[abs(dataset["PDGcode"])==22]["Weight"])
figname="BIB_decph"
pl.savefig(figname)


# In[24]:


fig = plt.figure(figsize=(8, 8))
fig.suptitle("$e^+$ $e^-$")
ax = fig.add_axes(rect_scatter)
ax_histx = fig.add_axes(rect_histx, sharex=ax)
ax_histy = fig.add_axes(rect_histy, sharey=ax)

scatter_histo(dataset[abs(dataset["PDGcode"])==11]["PosZ"], dataset[abs(dataset["PDGcode"])==11]["Time"], ax, ax_histx, ax_histy, weights=dataset[abs(dataset["PDGcode"])==11]["Weight"])
figname="BIB_decelpos"
pl.savefig(figname)


# In[25]:


fig = plt.figure(figsize=(8, 8))
fig.suptitle("charged hadrons")
ax = fig.add_axes(rect_scatter)
ax_histx = fig.add_axes(rect_histx, sharex=ax)
ax_histy = fig.add_axes(rect_histy, sharey=ax)

scatter_histo(dataset[abs(dataset["PDGcode"])==2112]["PosZ"], dataset[abs(dataset["PDGcode"])==2112]["Time"], ax, ax_histx, ax_histy, weights=dataset[abs(dataset["PDGcode"])==2112]["Weight"])
figname="BIB_decneu"
pl.savefig(figname)


# In[ ]:




