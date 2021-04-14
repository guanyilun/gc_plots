"""This script produces 3x1 subplots of magnetic field in each
frequency band

"""

import argparse, os, os.path as op
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import colors
from pixell import utils as u
from common import *
import lib, plotstyle

def process_map(imap):
    return np.log10(imap[0])

parser = argparse.ArgumentParser()
parser.add_argument("-o","--odir",default="plots")
parser.add_argument('--smooth', help='fwhm of smoothing kernel in arcmin', type=float, default=3.5)
parser.add_argument('--len', help='length as a fraction of input data shape', type=float, default=0.25)
parser.add_argument('-L', help='absolute length in degree', type=float)
parser.add_argument('--area', default='half')
parser.add_argument('--oname',default='mag_3bands.pdf')
parser.add_argument('--cmap', default='planck')
parser.add_argument('--axis', action='store_true')
parser.add_argument('--figsize', default=None)
args = parser.parse_args()
if not op.exists(args.odir): os.makedirs(args.odir)
if args.figsize: figsize = eval(args.figsize)
else: figsize = None
if args.L: L = args.L * 120  # assuming 0.5 arcmin pixel size
else: L = None
# define box of interests
box = boxes[args.area]
# initialize figure
# load tmp file for wcs
imap = load_map(filedb['f090']['coadd'], box, fcode='f090')/1e9
fig, axes = plt.subplots(3,1,figsize=figsize, sharex=True, subplot_kw={'projection': imap.wcs})

plot_opts = {
    'origin': 'lower',
    'cmap': args.cmap,
    'interpolation': 'nearest',
}
def process_map(imap, fill=1e-5):
    imap[0,imap[0]<=0] = fill
    return imap[0]

for i, fcode in zip(range(3), ['f090','f150','f220']):
    # load data
    imap = load_map(filedb[fcode]['coadd'], box, fcode=fcode)/1e9
    # plotstyle
    if   i == 0: norm = colors.LogNorm(vmin=10**-0.5, vmax=10**1.5)
    elif i == 1: norm = colors.LogNorm(vmin=10**-0.5, vmax=10**1.5)
    else:        norm = colors.LogNorm(vmin=10**0.5,  vmax=10**2)
    tmap = imap[0]
    # print(np.sum(tmap<=0))
    # continue
    if args.smooth > 0: pmap = enmap.smooth_gauss(imap, sigma=args.smooth/2.355*u.arcmin)[1:]
    else:               pmap = imap[1:]
    # compute texture
    theta = lib.Bangle(pmap[0], pmap[1], toIAU=True)
    # if L is not None, it will be used, otherwise length (fraction) will be used
    texture = lib.LIC_texture(theta, length=args.len, L=L)
    if i != 2:
        plotstyle.setup_axis(axes[i], xticks=False, yticks=True, nticks=[10,5])
    else:
        plotstyle.setup_axis(axes[i], xticks=True, yticks=True, nticks=[10,5])
    im = axes[i].imshow(tmap, norm=norm, **plot_opts)
    cbar = plotstyle.add_colorbar(fig, axes[i], size="3%")
    if i == 1: fig.colorbar(im, cax=cbar).set_label(texify("Total Intensity [MJy/sr]"), fontsize=14)
    else:      fig.colorbar(im, cax=cbar)
    if not args.axis: axes[i].axis('off')
    else:
        if i == 2:
            axes[i].set_xlabel('$l$')
            axes[i].set_ylabel('$b$')
        else:
            axes[i].set_ylabel('$b$')
    axes[i].imshow(texture, origin='lower', cmap='binary', alpha=0.7, interpolation='nearest')
    props = dict(alpha=1, facecolor='white')
    axes[i].text(0.03, 0.85, texify(fcode), bbox=props,
                 verticalalignment='bottom', horizontalalignment='left',
                 transform=axes[i].transAxes, fontsize=14)

plt.tight_layout(h_pad=0.01)
ofile = op.join(args.odir, args.oname)
print("Writing:", ofile)
plt.savefig(ofile, bbox_inches='tight')
