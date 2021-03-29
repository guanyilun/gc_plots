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

def process_map(imap):
    return np.log10(imap[0])

parser = argparse.ArgumentParser()
parser.add_argument("-o","--odir", default='plots')
parser.add_argument("--oname", default="test.pdf")
parser.add_argument("--use", default='coadd')
parser.add_argument("--area", default='half')
args = parser.parse_args()
if not op.exists(args.odir): os.makedirs(args.odir)

# Build figure
fig, axes = plt.subplots(3, 2, figsize=(9,6.5), gridspec_kw={'width_ratios':[0.93,1]})
box = boxes[args.area]

plot_opts = {
    'origin': 'lower',
    'cmap': 'planck_half',
    'extent': box2extent(box)/np.pi*180
}

for i, freq in zip(range(3), ['f090','f150','f220']):
    # plotstyle
    if   i == 0: norm = colors.LogNorm(vmin=10**-0.5, vmax=10**1.5)
    elif i == 1: norm = colors.LogNorm(vmin=10**-0.5, vmax=10**1.5)
    else:        norm = colors.LogNorm(vmin=10**0.5,  vmax=10**2)
    # load data
    imap = load_map(filedb[freq]['planck'], box, fcode=freq)/1e9
    im_p = axes[i,0].imshow(imap[0], norm=norm, **plot_opts)
    imap = load_map(filedb[freq][args.use], box, fcode=freq)/1e9
    im_c = axes[i,1].imshow(imap[0], norm=norm, **plot_opts)
    cb   = plotstyle.add_colorbar(fig, axes[i,1])
    # only add colorbar label to the middle panel
    if i == 1: fig.colorbar(im_c, cax=cb).set_label('Total Intensity [MJy/sr]')
    else: fig.colorbar(im_c, cax=cb)

for i in range(3):  # row
    for j in range(2):  # col
        ax = axes[i,j]
        if i == 2 and j == 0: pass
        else: 
            ax.axes.yaxis.set_ticklabels([])
            ax.axes.xaxis.set_ticklabels([])
        if j == 0:
            ax.yaxis.set_ticks_position('left')
            # ax.spines['right'].set_visible(False)
        if j == 1:
            ax.yaxis.set_ticks_position('right')
            # ax.spines['left'].set_visible(False)
        # if i != 2:
            # ax.yaxis.set_ticks_position('right')

axes[1,0].axes.yaxis.set_ticklabels([])
axes[-1,1].axes.xaxis.set_ticklabels([])

axes[-1,0].set_xlabel('l [deg]')
axes[-1,0].set_ylabel('b [deg]')

# setup labels: Planck, ACT+Planck
axes[0,0].text(0.4, 1.05, 'Planck',
               verticalalignment='bottom', horizontalalignment='left',
               transform=axes[0,0].transAxes, fontsize=12)
axes[0,1].text(0.35, 1.05, 'ACT+Planck',
               verticalalignment='bottom', horizontalalignment='left',
               transform=axes[0,1].transAxes, fontsize=12)

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
plt.tight_layout(w_pad=-0.3, h_pad=-0.3)
props = dict(alpha=1, facecolor='white')
plt.text(0.47, 0.9, 'f090', transform=fig.transFigure, fontsize=12,
         usetex=True, horizontalalignment='left', bbox=props)
plt.text(0.47, 0.61, 'f150', transform=fig.transFigure, fontsize=12,
         usetex=True, horizontalalignment='left', bbox=props)
plt.text(0.47, 0.32, 'f220', transform=fig.transFigure, fontsize=12,
         usetex=True, horizontalalignment='left', bbox=props)
ofile = op.join(args.odir, args.oname)
print("Writing:", ofile)
plt.savefig(ofile, bbox_inches='tight')
