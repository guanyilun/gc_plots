"""similar to plot source2 but focused on molecular clouds,
the major difference is that for molecular clouds we are more
interested in the f150 and f220, instead of f090

"""

import argparse, os, os.path
import numpy as np
import lib
from common import *
import plotstyle
from pixell import enplot, utils as u, enmap
from matplotlib import pyplot as plt
from scipy.optimize import leastsq
from mpl_toolkits.axes_grid1 import make_axes_locatable

parser = argparse.ArgumentParser()
parser.add_argument("-o","--odir", default="plots")
parser.add_argument("--oname", required=True)
parser.add_argument("--smooth", type=float, default=2)
parser.add_argument("--tonly", action='store_true', default=False)
parser.add_argument("--method", type=int, default=0)
parser.add_argument("-l", help="in deg", type=float, default=None)
parser.add_argument("-b", help="in deg", type=float, default=None)
parser.add_argument("--area", default=None)
parser.add_argument("--margin", help="margin in deg", type=float, default=0.1)
parser.add_argument("--cmap", default='planck')
parser.add_argument("--title", default=None)
parser.add_argument("--tmin", default="0,0,0")
parser.add_argument("--tmax", default="1,1,1")
parser.add_argument("--pmin", default="0,0,0")
parser.add_argument("--pmax", default="1,1,1")
parser.add_argument("--figsize", default="(8,7)")
parser.add_argument("--colorbar", action='store_true')
parser.add_argument("--use", default='coadd')
parser.add_argument("--freqs", default="f090,f150,f220")
args = parser.parse_args()
if not op.exists(args.odir): os.makedirs(args.odir)
tmin  = args.tmin.split(",")
tmax  = args.tmax.split(",")
pmin  = args.pmin.split(",")
pmax  = args.pmax.split(",")
freqs = args.freqs.split(',')

# define box
if args.area:
    box = boxes[args.area]
else:
    xmin, xmax = args.l-args.margin, args.l+args.margin
    ymin, ymax = args.b-args.margin, args.b+args.margin
    box = np.array([[ymin, xmax], [ymax, xmin]]) / 180*np.pi
# set figure size if needed
if args.figsize is not None: figsize=eval(args.figsize)
else: figsize=None

# reload maps in the actual box of interests
imap_f090 = load_map(filedb['f090'][args.use], fcode='f090', box=box)/1e9  # MJy/sr
imap_f150 = load_map(filedb['f150'][args.use], fcode='f150', box=box)/1e9
imap_f220 = load_map(filedb['f220'][args.use], fcode='f220', box=box)/1e9    

# short names for I and P maps
I_f090 = imap_f090[0]
I_f150 = imap_f150[0]
I_f220 = imap_f220[0]
P_f090 = np.sum(imap_f090[1:]**2, axis=0)**0.5
P_f150 = np.sum(imap_f150[1:]**2, axis=0)**0.5
P_f220 = np.sum(imap_f220[1:]**2, axis=0)**0.5
# make figure
popts = {
    'origin': 'lower',
    'cmap': args.cmap,
    'extent': box2extent(box)/np.pi*180
}
maps = ['T','P']
if args.tonly:
    nrow, ncol = 1, len(freqs)
else:
    nrow, ncol = 2, len(freqs)

fig, axes = plt.subplots(nrow, ncol, figsize=figsize, sharex=True, sharey=True)
data = {
    'f090': [I_f090, P_f090],
    'f150': [I_f150, P_f150],
    'f220': [I_f220, P_f220],    
}
# data = [[I_f090, I_f150, I_f220], [P_f090, P_f150, P_f220]]
mins = [tmin, pmin]
maxs = [tmax, pmax]

for i in range(nrow):
    for j in range(ncol):
        if args.tonly:
            if len(freqs) == 1:
                ax = axes
            else:
                ax = axes[j]
        else:
            if len(freqs) == 1:
                ax = axes[i]
            else:
                ax = axes[i,j]
        im = ax.imshow(data[freqs[j]][i], vmin=mins[i][j], vmax=maxs[i][j], **popts)
        if args.colorbar:
            divider = make_axes_locatable(ax)
            cb = divider.append_axes('right', size="5%", pad=0.1)
            fig.colorbar(im, cax=cb, orientation='vertical').set_label("MJy/sr")
for i in range(nrow):
    for j in range(ncol):
        if args.tonly:
            if len(freqs) == 1:
                ax = axes
            else:
                ax = axes[j]
        else:
            if len(freqs) == 1:
                ax = axes[i]
            else:
                ax = axes[i,j]        
        if i == 0:
            ax.text(0.4, 1.02, freqs[j],
                   verticalalignment='bottom', horizontalalignment='left',
                   transform=ax.transAxes, fontsize=14)
        if j == 0:
            # ax.text(-0.2, 0.5, maps[i], rotation=90,
            #        verticalalignment='bottom', horizontalalignment='left',
            #        transform=ax.transAxes, fontsize=14)
            ax.set_ylabel('b [deg]')
        if i == nrow-1:
            ax.set_xlabel('l [deg]')
if args.title is not None: plt.suptitle(args.title, y = 1.02, fontsize=16)
plt.tight_layout(h_pad=0, w_pad=0.1)
ofile = op.join(args.odir, args.oname)
print("Writing:", ofile)
plt.savefig(ofile, bbox_inches='tight')
