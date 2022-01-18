# Estimate read lengths

* DNA fgrom a single drone was purified to perform 3 libraries.
* Using SMRTBell template Prep Kit 1.0 (PacBio), a DNA and END damage repair step was performed on 15µg of unsheared sample. Then blunt hairpin adapters were ligated to the libraries. The libraries were treated with an exonuclease cocktail to digest unligated DNA fragments.
* A size selection step using a 7kb (Library 1 "PacificBee") or 9kb (libraries 2 "Abeille1" and 3 "Abeille2") 
* SMRTbell libraries were sequenced on 36 SMRTcells on RSII instrument

## Extract read lengths for the 36 SMRTcells:

```
#!/usr/bin/env python3
# -*- coding: Utf-8 -*-

import pandas as pd
import numpy as np
import os
import subprocess
import glob

fastqList = glob.glob('/home/gencel/vignal/work/PacificBeeAll/Run*/Analyse_H5toFastq*/*.gz')
for i in fastqList:
    path = i.split('/')
    run_smrt_lib = path[6].split('_')
    run = run_smrt_lib[1]
    smrt_lib = run_smrt_lib[2].split('.')
    smrt = smrt_lib[1]
    lib = path[8].split('.')[0]
    cmd = "zcat " + i + " | awk '{if(NR%4==2) print length($1)}' > /home/gencel/vignal/work/PacificBeeAll/ReadLength/ReadLength_" + run + "_" + smrt + "_" + lib + ".txt"
    os.system(cmd)
```

* Produces 36 files: single columns with read lengths
* Copied in /Users/avignal/Documents/Stats/2016_PacificBee/ReadLength

## Plot read lengths

```
import csv
import pandas as pd
import numpy as np
import os
import subprocess
import glob
from matplotlib import pyplot as plt
from IPython.display import display
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset

#Get list of read lengths files
readLengthsList = glob.glob('/Users/avignal/Documents/Stats/2016_PacificBee/ReadLength/*.txt')

#Import the data into a dataframe
for i in readLengthsList:
    path = i.split('/')
    run_smrt_lib1 = path[7].split('.')[0]
    run_smrt_lib2 = run_smrt_lib1.split("_")
    run_smrt_lib3 = run_smrt_lib2[1] + "\n" + run_smrt_lib2[2] + "\n" + run_smrt_lib2[3]

#Draw grid
xbins = 6
ybins = 6
zbins = 100
fig, ax = plt.subplots(ybins, xbins, figsize=(10,12), sharex=True, sharey=True)
fig.suptitle('PacBio read distribution')
# no space between subplots
fig.subplots_adjust(hspace=0, wspace=0)
# make some room for the axes' labels
fig.subplots_adjust(bottom=0.13, right=0.97, top=0.93)

# don't draw x and y-axes' labels (keep ticks)
for xi in range(ybins):
    for yi in range(xbins):
        plt.setp(ax[xi,yi].get_xticklabels(), visible=False)
        plt.setp(ax[xi,yi].get_yticklabels(), visible=False)

#Import the data and plot histograms       
count = 0
for yi in range(ybins):
    for xi in range(xbins):
        path = readLengthsList[count].split('/')
        run_smrt_lib1 = path[7].split('.')[0]
        run_smrt_lib2 = run_smrt_lib1.split("_")
        run_smrt_lib3 = run_smrt_lib2[1] + "\n" + run_smrt_lib2[2] + "\n" + run_smrt_lib2[3]     
        data=pd.read_csv(readLengthsList[count], names=['Length'])
        lengthsList = data['Length']
        ax[xi,yi].hist(lengthsList,bins=zbins)
        ax[xi,yi].text(0.1, 0.9, run_smrt_lib3, verticalalignment='top', transform=ax[xi,yi].transAxes)
        count = count + 1
#add axes names
axtot = fig.add_subplot(1,1,1)
axtot.axes.patch.set_alpha(0.)
axtot.set_xticks([])
axtot.set_yticks([])

#Zoomed x axis
# tie the zoomed in x-axis to the bottom plot, third from the left
ax1 = ax[ybins-1,2]
# zoom factor is such that about half the space is used
zoom_factor = float(xbins) / 1.5
# create the zoomed in axis based on the original axes
axins = zoomed_inset_axes(parent_axes=ax1, zoom=zoom_factor, loc=8,
                          bbox_to_anchor=(0.5, 0.),
                          bbox_transform=ax1.transAxes,
                          axes_kwargs=dict(sharex=ax1, sharey=ax1),
                          borderpad=-1.5, # in fraction of font size
                          )
# dotted lines to show that this is the x-axis
pp, p1, p2 = mark_inset(parent_axes=ax1, inset_axes=axins,
                        loc1=3, loc2=4.,
                        linestyle="dotted")

# only want to draw some of the lines
pp.set_visible(False)

# hide the plotting area
#axins.axesPatch.set_alpha(0.)
axins.axes.patch.set_alpha(0.)

# we want to draw the bottom spine only
axins.set_frame_on(True)
axins.spines['top'].set_visible(False)
axins.spines['left'].set_visible(False)
axins.spines['right'].set_visible(False)

# don't draw the y axis ticks or labels
axins.set_yticks([])
axins.set_yticklabels([])

# only draw the bottom (x) axis
axins.xaxis.set_ticks_position('bottom')
axins.xaxis.set_label_position('bottom')
major_ticks = np.arange(0,75000,10000)
axins.set_xticks(major_ticks)
#major_ticks = np.arange(0,df_chrom['EndMb'].max(),1)

# The inset axis label
axins.set_xlabel('read length (bp)')
```

