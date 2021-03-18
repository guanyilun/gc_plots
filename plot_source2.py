"""hopefully a more consice version of plot_source.py, particularly
in terms of dust removal. I will give the option to do different
dust removal approach here

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
parser.add_argument("--oname", default="mouse.pdf")
parser.add_argument("--smooth", type=float,default=2)
parser.add_argument("--tonly", action='store_true', default=False)
parser.add_argument("--method", type=int, default=0)
parser.add_argument("--dust-removal", action='store_true')
parser.add_argument("--dust-factor-f090", type=float, default=1)
parser.add_argument("--dust-factor-f150", type=float, default=1)
parser.add_argument("--auto-adj", action='store_true')
parser.add_argument("-l", help="in deg", type=float, default=None)
parser.add_argument("-b", help="in deg", type=float, default=None)
parser.add_argument("--margin", help="margin in deg", type=float, default=0.1)
parser.add_argument("--cmap", default='planck')
parser.add_argument("--title", default=None)
parser.add_argument("--tmin", default="0,0,0")
parser.add_argument("--tmax", default="1,1,1")
parser.add_argument("--pmin", default="0,0,0")
parser.add_argument("--pmax", default="1,1,1")
parser.add_argument("--dust-area", help='an area used for dust factor estimation', default=None)
parser.add_argument("--figsize", default="(8,7)")
parser.add_argument("--colorbar", action='store_true')

args = parser.parse_args()
if not op.exists(args.odir): os.makedirs(args.odir)
tmin = args.tmin.split(",")
tmax = args.tmax.split(",")
pmin = args.pmin.split(",")
pmax = args.pmax.split(",")

# define box
xmin, xmax = args.l-args.margin, args.l+args.margin
ymin, ymax = args.b-args.margin, args.b+args.margin
box = np.array([[ymin, xmax], [ymax, xmin]]) / 180*np.pi
# set figure size if needed
if args.figsize is not None: figsize=eval(args.figsize)
else: figsize=None

# estimate dust removing factor, the idea is to set the factor as x
# and try to find the best x such that f090-x*f220 gives the minimum dust
if args.dust_removal:
    # load map for dust area only
    if args.dust_area: dust_box = boxes[args.dust_area]
    else: dust_box = box
    imap_f090 = load_map(filedb['f090']['coadd'], fcode='f090', box=dust_box)/1e9  # MJy/sr
    imap_f150 = load_map(filedb['f150']['coadd'], fcode='f150', box=dust_box)/1e9
    imap_f220 = load_map(filedb['f220']['coadd'], fcode='f220', box=dust_box)/1e9    
    # form fiducial dust template at each frequency
    beta = 1.53
    s_f090 = (100/217)**(beta+2)
    s_f150 = (143/217)**(beta+2)
    dust_f090 = lib.beam_match(imap_f220, 'f090', 'f220')*s_f090
    dust_f150 = lib.beam_match(imap_f220, 'f150', 'f220')*s_f150
    
    if args.method == 0:
        # method 0: simply subtract dust template
        x_f090 = x_f150 = 1
    elif args.method == 1:
        # method 2: use user specified values
        x_f090 = args.dust_factor_f090
        x_f150 = args.dust_factor_f150        
    elif args.method == 2:
        # method 1: use least square
        # optionally use a specific area to find dust index
        I_f090_ = imap_f090.submap
        par, _ = leastsq(lambda x: np.ravel(imap_f090[0]-x*dust_f090[0]), x0=1)
        x_f090 = par[0]
        par, _ = leastsq(lambda x: np.ravel(imap_f150[0]-x*dust_f150[0]), x0=1)
        x_f150 = par[0]
    elif args.method == 3:
        # method 3: force median of f090-x*f220 to 0
        # or in other words, we minimize abs(median(f090-x*f220))
        from scipy.optimize import minimize
        res = minimize(lambda x: np.abs(np.median(imap_f090[0]-x*dust_f090[0])), x0=1, method="Nelder-Mead")
        x_f090 = res.x[0]
        res = minimize(lambda x: np.abs(np.median(imap_f150[0]-x*dust_f150[0])), x0=1, method="Nelder-Mead")
        x_f150 = res.x[0]
else:
    x_f090 = f150 = 0
    
print("dust factor in f090:", x_f090)
print("dust factor in f150:", x_f150)

# reload maps in the actual box of interests
imap_f090 = load_map(filedb['f090']['coadd'], fcode='f090', box=box)/1e9  # MJy/sr
imap_f150 = load_map(filedb['f150']['coadd'], fcode='f150', box=box)/1e9
imap_f220 = load_map(filedb['f220']['coadd'], fcode='f220', box=box)/1e9    
# form fiducial dust template at each frequency
beta = 1.53
s_f090 = (100/217)**(beta+2)
s_f150 = (143/217)**(beta+2)
dust_f090 = lib.beam_match(imap_f220, 'f090', 'f220')*s_f090
dust_f150 = lib.beam_match(imap_f220, 'f150', 'f220')*s_f150
imap_f090 -= x_f090 * dust_f090
imap_f150 -= x_f150 * dust_f150
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
if args.dust_removal: names = ['f090 - f220','f150 - f220']
else: names = ['f090', 'f150', 'f220']
maps = ['T','P']
if args.tonly:
    nrow, ncol = 1, len(names)
else:
    nrow, ncol = 2, len(names)

fig, axes = plt.subplots(nrow, ncol, figsize=figsize, sharex=True, sharey=True)
data = [[I_f090, I_f150, I_f220], [P_f090, P_f150, P_f220]]
mins = [tmin, pmin]
maxs = [tmax, pmax]
for i in range(nrow):
    for j in range(ncol):
        if args.tonly: ax = axes[j]
        else: ax = axes[i,j]
        im = ax.imshow(data[i][j], vmin=mins[i][j], vmax=maxs[i][j], **popts)
        if args.colorbar:
            divider = make_axes_locatable(ax)
            cb = divider.append_axes('right', size="5%", pad=0.1)
            fig.colorbar(im, cax=cb, orientation='vertical').set_label("MJy/sr")
for i in range(nrow):
    for j in range(ncol):
        if args.tonly: ax = axes[j]
        else: ax = axes[i,j]        
        if i == 0:
            ax.text(0.3, 1.02, names[j],
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
plt.tight_layout(h_pad=0.1, w_pad=0.1)
ofile = op.join(args.odir, args.oname)
print("Writing:", ofile)
plt.savefig(ofile, bbox_inches='tight')
