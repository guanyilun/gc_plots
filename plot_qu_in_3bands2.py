"""Improved version of plot_qu_in_3bands with better layout and
colorbar.

From original script: this script plots Q/U maps in 3 bands

"""

import argparse, os, os.path as op
import numpy as np
import matplotlib as mpl
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
parser.add_argument("--min-f090",type=float, default=-0.2)
parser.add_argument("--max-f090",type=float, default= 0.2)
parser.add_argument("--min-f150",type=float, default=-0.2)
parser.add_argument("--max-f150",type=float, default= 0.2)
parser.add_argument("--min-f220",type=float, default=-1)
parser.add_argument("--max-f220",type=float, default= 1)
parser.add_argument("--area", default="half")
parser.add_argument("--figsize", default="(9,7)")
parser.add_argument("--cmap", default="planck")
parser.add_argument("--axis", action='store_true')
parser.add_argument("--IAU", action='store_true')
parser.add_argument("--sep-colorbar", action='store_true')
args = parser.parse_args()
if not op.exists(args.odir): os.makedirs(args.odir)
if args.axis:
    mpl.rcParams['xtick.direction']='in'
    mpl.rcParams['ytick.direction']='in'
    mpl.rcParams["xtick.major.size"]=2
    mpl.rcParams["ytick.major.size"]=2
    mpl.rcParams["xtick.minor.size"]=1
    mpl.rcParams["ytick.minor.size"]=1
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
# load temp file to get wcs
imap = load_map(filedb['f090']['coadd'], box, fcode='f090', IAU=args.IAU)
spkw = {'projection': imap.wcs}
if not args.sep_colorbar:
    fig, axes = plt.subplots(3,2,figsize=figsize, subplot_kw=spkw)
else:
    fig, axes = plt.subplots(3,2,figsize=figsize, subplot_kw=spkw,
                             gridspec_kw={'width_ratios': [0.942, 1]})    
plot_opts = {
    'origin': 'lower',
    'cmap': args.cmap,
    # 'extent': box2extent(box)/np.pi*180
}
# row 0: f090
plotstyle.setup_axis(axes[0,0], nticks=[10,5], xticks=False, yticks=True)
plotstyle.setup_axis(axes[0,1], nticks=[10,5], xticks=False, yticks=False)
plot_opts.update({'vmin': args.min_f090, 'vmax': args.max_f090})
imap = load_map(filedb['f090']['coadd'], box, fcode='f090', IAU=args.IAU)
qmap = process_map(imap, 'Q')
im_f090_q = axes[0,0].imshow(qmap, **plot_opts)
umap = process_map(imap, 'U')
im_f090_u = axes[0,1].imshow(umap, **plot_opts)
# add colorbar
if args.sep_colorbar:
    cb_f090 = plotstyle.add_colorbar(fig, axes[0,1], size="3%")
    fig.colorbar(im_f090_q, cax=cb_f090, orientation='vertical')
    cb_f090.yaxis.set_ticks([-0.15,0,0.15]);    

# row 1: f150
plotstyle.setup_axis(axes[1,0], nticks=[10,5], xticks=False, yticks=True)
plotstyle.setup_axis(axes[1,1], nticks=[10,5], xticks=False, yticks=False)
plot_opts.update({'vmin': args.min_f150, 'vmax': args.max_f150})
imap = load_map(filedb['f150']['coadd'], box, fcode='f150', IAU=args.IAU)
qmap = process_map(imap, 'Q')
im_f150_q = axes[1,0].imshow(qmap, **plot_opts)
umap = process_map(imap, 'U')
im_f150_u = axes[1,1].imshow(umap, **plot_opts)
if args.sep_colorbar:
    cb_f150 = plotstyle.add_colorbar(fig, axes[1,1], size="3%")
    fig.colorbar(im_f150_q, cax=cb_f150, orientation='vertical').set_label(texify("MJy/sr"), fontsize=14)
    cb_f150.yaxis.set_ticks([-0.15,0,0.15]);    

# row 2: f220
plotstyle.setup_axis(axes[2,0], nticks=[10,5], xticks=True, yticks=True)
plotstyle.setup_axis(axes[2,1], nticks=[10,5], xticks=True, yticks=False)
plot_opts.update({'vmin': args.min_f220, 'vmax': args.max_f220})
imap = load_map(filedb['f220']['coadd'], box, fcode='f220', IAU=args.IAU)
qmap = process_map(imap, 'Q')
im_f220_q = axes[2,0].imshow(qmap, **plot_opts)
umap = process_map(imap, 'U')
im_f220_u = axes[2,1].imshow(umap, **plot_opts)
# add colorbar
if args.sep_colorbar:
    cb_f220 = plotstyle.add_colorbar(fig, axes[2,1], size="3%")
    fig.colorbar(im_f220_q, cax=cb_f220, orientation='vertical')
    # cb_f150.yaxis.set_ticks([-0.15,0,0.15]);        

# axes settings
for i in range(3):
    for j in range(2):
        ax = axes[i,j]
        if not args.axis:
            ax.set_axis_off()
        else:
            if j == 0: ax.yaxis.set_ticks_position('left')
            if j == 1:
                ax.yaxis.set_ticks_position('right')
                ax.axes.yaxis.set_ticklabels([])
            # if (j!=0) or (i !=2):
            #     ax.axes.yaxis.set_ticklabels([])
            #     ax.axes.xaxis.set_ticklabels([])
            # else:
            if i == 2 and j == 0:
                ax.set_xlabel('$l$')
                ax.set_ylabel('$b$')
            else:
                ax.set_xlabel(' ')
                ax.set_ylabel(' ')
            # if i==2:
                # ax.xaxis.set_major_locator(ticker.FixedLocator([1.5,1,0.5,0,-0.5,-1,-1.5]))
            # if j==0:
                # ax.yaxis.set_major_locator(ticker.FixedLocator([-0.5,0,0.5]))
# labels
axes[0,0].text(0.5, 1.05, texify('Q'),
               verticalalignment='bottom', horizontalalignment='left',
               transform=axes[0,0].transAxes, fontsize=14)
axes[0,1].text(0.5, 1.05, texify('U'),
               verticalalignment='bottom', horizontalalignment='left',
               transform=axes[0,1].transAxes, fontsize=14)

# setup labels: f090, f150, f220
# for i, label in zip(range(3), ['f090','f150','f220']):
#     axes[i,0].text(-0.05, 0.35, label, rotation=90,
#                    verticalalignment='bottom', horizontalalignment='left',
#                    transform=axes[i,0].transAxes, fontsize=12)
props = dict(alpha=1, facecolor='white')
plt.text(0.48, 0.82, texify('f090'), transform=fig.transFigure, fontsize=12,
         usetex=True, horizontalalignment='left', bbox=props)
plt.text(0.48, 0.56, texify('f150'), transform=fig.transFigure, fontsize=12,
         usetex=True, horizontalalignment='left', bbox=props)
plt.text(0.48, 0.31, texify('f220'), transform=fig.transFigure, fontsize=12,
         usetex=True, horizontalalignment='left', bbox=props)

# add colorbar
if not args.sep_colorbar:
    fig.subplots_adjust(right=0.85,wspace=0,hspace=0.)
    cb_f090 = fig.add_axes([0.86, 0.11, 0.02, 0.77])
    fig.colorbar(im_f090_q, cax=cb_f090, orientation='vertical').set_label("MJy/sr")
else:
    fig.subplots_adjust(wspace=0,hspace=0.)
# fig.colorbar(im, ax=axes, shrink=0.8).set_label(label="mJy/sr")

ofile = op.join(args.odir, args.oname)
print("Writing:", ofile)
plt.savefig(ofile, bbox_inches='tight')
