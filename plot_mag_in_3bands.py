"""This script produces 3x1 subplots of magnetic field in each
frequency band

"""

import argparse, os, os.path as op
import numpy as np
from matplotlib import pyplot as plt
from pixell import utils as u
from common import *
import lib, plotstyle

def process_map(imap):
    return np.log10(imap[0])

parser = argparse.ArgumentParser()
parser.add_argument("-o","--odir",default="plots")
parser.add_argument('--smooth', help='fwhm of smoothing kernel in arcmin', type=float, default=3.5)
parser.add_argument('--len', help='length as a fraction of input data shape', type=float, default=0.25)
parser.add_argument('--area', default='half')
parser.add_argument('--oname',default='mag_3bands.pdf')
parser.add_argument('--figsize', default=None)
args = parser.parse_args()
if not op.exists(args.odir): os.makedirs(args.odir)
if args.figsize: figsize = eval(args.figsize)
else: figsize = None

# define box of interests
box = boxes[args.area]
# initialize figure
fig, axes = plt.subplots(3,1,figsize=figsize)

for i, fcode in zip(range(3), ['f090','f150','f220']):
    # load data
    imap = load_map(filedb[fcode]['coadd'], box, fcode=fcode)
    tmap = process_map(imap)
    if args.smooth > 0: pmap = enmap.smooth_gauss(imap, sigma=args.smooth/2.355*u.arcmin)[1:]
    else:               pmap = imap[1:]
    # compute texture
    theta = lib.Bangle(pmap[0], pmap[1], toIAU=True)
    texture = lib.LIC_texture(theta, length=args.len)
    plot_opts = {
        'origin': 'lower',
        'cmap': 'planck',
        'vmin': 7,
        'vmax': 10.5,
        'extent': box2extent(box)/np.pi*180,
    }
    if fcode == 'f220':
        plot_opts.update({'vmin':8,'vmax':11})
    axes[i].imshow(tmap, **plot_opts)
    axes[i].axis('off')
    ax = fig.add_subplot(3,1,i+1)
    ax.imshow(texture, origin='lower', cmap='binary', alpha=0.5)
    ax.axis('off')
    axes[i].text(-0.03, 0.4, fcode, rotation=90,
                 verticalalignment='bottom', horizontalalignment='left',
                 transform=axes[i].transAxes, fontsize=12)

plt.tight_layout(h_pad=0.2)
ofile = op.join(args.odir, args.oname)
print("Writing:", ofile)
plt.savefig(ofile, bbox_inches='tight')
