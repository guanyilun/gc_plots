"""This script plots planck and act temperature maps side by side as a 3x2 grid

2021 Mar 26
- add colorbar and axis coordinates
"""

import argparse, os, os.path as op
import numpy as np, glob
import matplotlib as mpl
mpl.rcParams['xtick.direction']='in'
mpl.rcParams['ytick.direction']='in'
mpl.rcParams["xtick.major.size"]=2
mpl.rcParams["ytick.major.size"]=2
mpl.rcParams["xtick.minor.size"]=1
mpl.rcParams["ytick.minor.size"]=1

import matplotlib.pyplot as plt
from matplotlib import colors
from pixell import enmap, colorize

# setup common plotstyle
import plotstyle
# load common variables and utility functions
from common import *

parser = argparse.ArgumentParser()
parser.add_argument("-o","--odir", default='plots')
parser.add_argument("--oname", default="test.pdf")
parser.add_argument("--use", default='coadd')
parser.add_argument("--area", default='half')
parser.add_argument("--figsize", default="(9,6.5)")
args = parser.parse_args()
if not op.exists(args.odir): os.makedirs(args.odir)
if args.figsize: figsize = eval(args.figsize)
else: figsize = None
box = boxes[args.area]

# Build figure
# load a temperorary file to get wcs right
imap = load_map(filedb['f090']['planck'], box, fcode='f090')/1e9
fig, axes = plt.subplots(3, 2, figsize=figsize, gridspec_kw={'width_ratios':[0.93,1]},
                         subplot_kw={'projection': imap.wcs})

plot_opts = {
    'origin': 'lower',
    'cmap': 'planck_half',
}

for i, freq in zip(range(3), ['f090','f150','f220']):
    for j in range(2):
        # plotstyle
        if   i == 0: norm = colors.LogNorm(vmin=10**-0.5, vmax=10**1.5)
        elif i == 1: norm = colors.LogNorm(vmin=10**-0.5, vmax=10**1.5)
        else:        norm = colors.LogNorm(vmin=10**0.5,  vmax=10**2)
        # setup axes projection
        ax = axes[i,j]
        if i == 2 and j == 1:
            xticks, yticks = True, False
        elif i == 2 and j == 0:
            xticks, yticks = True, True
        elif j == 0:
            xticks, yticks = False, True
        else:
            xticks, yticks = False, False
        if i == 2:
            ax.set_xlabel('$l$')
        if j == 0:
            ax.set_ylabel('$b$')            
        plotstyle.setup_axis(ax, fmt="d.d", xticks=xticks, yticks=yticks, nticks=[10,5])
        # load data
        if j == 0: use = 'planck'
        else: use = args.use
        imap = load_map(filedb[freq][use], box, fcode=freq)/1e9
        im = ax.imshow(imap[0], norm=norm, **plot_opts)
        if j == 1:
            cb = plotstyle.add_colorbar(fig, axes[i,1])
            # only add colorbar label to the middle panel            
            if i == 1:
                fig.colorbar(im, cax=cb).set_label(texify('Total Intensity [MJy/sr]'), fontsize=14)
            else:
                fig.colorbar(im, cax=cb)

for i in range(3):  # row
    for j in range(2):  # col
        ax = axes[i,j]
        if i == 2 and j == 0: pass
        else: 
            ax.axes.yaxis.set_ticklabels([])
            ax.axes.xaxis.set_ticklabels([])
        if j == 0:
            # ax.yaxis.set_ticks_position('left')
            ax.coords[1].set_ticks_position('left')
            # ax.spines['right'].set_visible(False)
        if j == 1:
            # ax.yaxis.set_ticks_position('right')
            ax.coords[1].set_ticks_position('right')            
            # ax.spines['left'].set_visible(False)
        # if i != 2:
            # ax.yaxis.set_ticks_position('right')

axes[1,0].axes.yaxis.set_ticklabels([])
axes[-1,1].axes.xaxis.set_ticklabels([])

axes[-1,0].set_xlabel('$l$')
axes[-1,0].set_ylabel('$b$')

# setup labels: Planck, ACT+Planck
axes[0,0].text(0.4, 1.05, texify('Planck'),
               verticalalignment='bottom', horizontalalignment='left',
               transform=axes[0,0].transAxes, fontsize=14)
axes[0,1].text(0.35, 1.05, texify('ACT+Planck'),
               verticalalignment='bottom', horizontalalignment='left',
               transform=axes[0,1].transAxes, fontsize=14)

# setup labels: f090, f150, f220

# axes[0,0].text(-0.12, 0.40, 'f090', rotation=90,
#                verticalalignment='bottom', horizontalalignment='left',
#                transform=axes[0,0].transAxes, fontsize=12)
# axes[1,0].text(-0.12, 0.40, 'f150', rotation=90,
#                verticalalignment='bottom', horizontalalignment='left',
#                transform=axes[1,0].transAxes, fontsize=12)
# axes[2,0].text(-0.12, 0.40, 'f220', rotation=90,
#                verticalalignment='bottom', horizontalalignment='left',
#                transform=axes[2,0].transAxes, fontsize=12)
plt.tight_layout(w_pad=0.01, h_pad=0)
props = dict(alpha=1, facecolor='white')
plt.text(0.44, 0.91,  texify('f090'), transform=fig.transFigure, fontsize=12,
         usetex=True, horizontalalignment='left', bbox=props)
plt.text(0.44, 0.595, texify('f150'), transform=fig.transFigure, fontsize=12,
         usetex=True, horizontalalignment='left', bbox=props)
plt.text(0.44, 0.29,  texify('f220'), transform=fig.transFigure, fontsize=12,
         usetex=True, horizontalalignment='left', bbox=props)
ofile = op.join(args.odir, args.oname)
print("Writing:", ofile)
plt.savefig(ofile, bbox_inches='tight')
