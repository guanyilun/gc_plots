"""Improved version of plot_qu_in_3bands with better layout and
colorbar.

From original script: this script plots Q/U maps in 3 bands

"""

import argparse, os, os.path as op
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.ticker as ticker
from mpl_toolkits.axes_grid1 import make_axes_locatable
from pixell import utils as u
import plotstyle
from common import *

parser = argparse.ArgumentParser()
parser.add_argument("-o","--odir",default="plots")
parser.add_argument('--smooth', type=float, default=0)
parser.add_argument("--oname",default="QU_3bands.pdf")
parser.add_argument("--min-f090",type=float, default=-0.3)
parser.add_argument("--max-f090",type=float, default= 0.3)
parser.add_argument("--min-f150",type=float, default=-0.3)
parser.add_argument("--max-f150",type=float, default= 0.3)
parser.add_argument("--min-f220",type=float, default=-0.3)
parser.add_argument("--max-f220",type=float, default= 0.3)
parser.add_argument("--area", default="half")
parser.add_argument("--figsize", default="(9,7)")
parser.add_argument("--cmap", default="planck")
parser.add_argument("--axis", action='store_true')
args = parser.parse_args()
if not op.exists(args.odir): os.makedirs(args.odir)

# load datafile
def process_map(m, component='Q'):
    if   component == 'Q': i = 1
    elif component == 'U': i = 2
    else: raise ValueError("Only Q/U are allowed")
    if args.smooth > 0:
        omap = enmap.smooth_gauss(m[i], args.smooth*u.fwhm*u.arcmin)
    else: omap = m[i]
    return omap/1e9
box = boxes[args.area]

if args.figsize: figsize = eval(args.figsize)
else: figsize = None
# initialize plot
fig, axes = plt.subplots(3,2,figsize=figsize)
plot_opts = {
    'origin': 'lower',
    'cmap': args.cmap,
    'extent': box2extent(box)/np.pi*180
}
# row 0: f090
plot_opts.update({'vmin': args.min_f090, 'vmax': args.max_f090})
imap = load_map(filedb['f090']['coadd'], box, fcode='f090')
qmap = process_map(imap, 'Q')
im_f090_q = axes[0,0].imshow(qmap, **plot_opts)
umap = process_map(imap, 'U')
im_f090_u = axes[0,1].imshow(umap, **plot_opts)
# add colorbar <-- not a good approach
# divider_f090 = make_axes_locatable(axes[0,1])
# cb_f090 = divider_f090.append_axes('right', size="5%", pad=0.1)
# fig.colorbar(im_f090_q, cax=cb_f090, orientation='vertical')

# row 1: f150
plot_opts.update({'vmin': args.min_f150, 'vmax': args.max_f150})
imap = load_map(filedb['f150']['coadd'], box, fcode='f150')
qmap = process_map(imap, 'Q')
im_f150_q = axes[1,0].imshow(qmap, **plot_opts)
umap = process_map(imap, 'U')
im_f150_u = axes[1,1].imshow(umap, **plot_opts)

# row 2: f220
plot_opts.update({'vmin': args.min_f220, 'vmax': args.max_f220})
imap = load_map(filedb['f220']['coadd'], box, fcode='f220')
qmap = process_map(imap, 'Q')
im_f220_q = axes[2,0].imshow(qmap, **plot_opts)
umap = process_map(imap, 'U')
im_f220_u = axes[2,1].imshow(umap, **plot_opts)

# axes settings
for i in range(3):
    for j in range(2):
        ax = axes[i,j]
        if not args.axis:
            ax.set_axis_off()
        else:
            if j>0:
                ax.axes.yaxis.set_ticklabels([])
            if i<2:
                ax.axes.xaxis.set_ticklabels([])
            if i==2:
                ax.xaxis.set_major_locator(ticker.FixedLocator([1.5,1,0.5,0,-0.5,-1,-1.5]))
                ax.set_xlabel('l [deg]')
            if j==0:
                ax.yaxis.set_major_locator(ticker.FixedLocator([-0.5,0,0.5]))
                ax.set_ylabel('b [deg]')

# labels
axes[0,0].text(0.5, 1.05, 'Q',
               verticalalignment='bottom', horizontalalignment='left',
               transform=axes[0,0].transAxes, fontsize=12)
axes[0,1].text(0.5, 1.05, 'U',
               verticalalignment='bottom', horizontalalignment='left',
               transform=axes[0,1].transAxes, fontsize=12)

# setup labels: f090, f150, f220
for i, label in zip(range(3), ['f090','f150','f220']):
    axes[i,0].text(-0.05, 0.35, label, rotation=90,
                   verticalalignment='bottom', horizontalalignment='left',
                   transform=axes[i,0].transAxes, fontsize=12)

# add colorbar
fig.subplots_adjust(right=0.83,wspace=0.02, hspace=0.02)
# left, bottom, width, height
# cb_f090 = fig.add_axes([0.83, 0.6, 0.02, 0.3])
# cb_f150 = fig.add_axes([0.83, 0.3, 0.02, 0.3])
# cb_f220 = fig.add_axes([0.83, 0.03, 0.02, 0.3])
# fig.colorbar(im_f090_q, cax=cb_f090, orientation='vertical')
# fig.colorbar(im_f150_q, cax=cb_f150, orientation='vertical')
# fig.colorbar(im_f220_q, cax=cb_f220, orientation='vertical')
cb_f090 = fig.add_axes([0.85, 0.12, 0.02, 0.76])
fig.colorbar(im_f090_q, cax=cb_f090, orientation='vertical').set_label("MJy/sr")
# fig.colorbar(im, ax=axes, shrink=0.8).set_label(label="mJy/sr")

ofile = op.join(args.odir, args.oname)
print("Writing:", ofile)
plt.savefig(ofile, bbox_inches='tight')
