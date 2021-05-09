"""This will make a two panel plot with upper panel being
the total intensity and lower panel being the polarization
intensity for comparison. It actually does no computing but
to load saved maps computed separately using plot_multifreq
or plot_multifreq2 (for polarization).

"""

import argparse, os, os.path as op
import numpy as np
import matplotlib.pyplot as plt
from pixell import enmap, enplot, utils as u
from matplotlib import colors

from common import *
import lib
import plotstyle

# parser defined in common
parser.add_argument("-T", default=None)
parser.add_argument("-P", default=None)
parser.add_argument("--axis", action='store_true')
parser.add_argument("--figsize", default="(8,4)")
parser.add_argument("--texture", help='file to store texture', default='mf_texture.npy')
parser.add_argument('--force', action='store_true')
args = parser.parse_args()
if not op.exists(args.odir): os.makedirs(args.odir)
if args.figsize: figsize=eval(args.figsize)
else: figsize=None
# define box of interests
box = boxes[args.area]
# load a map for wcs only
imap = load_map(filedb['f090']['coadd'], box=box, fcode='f090')
# load two maps
tmap = np.load(args.T)

# start plotting
popts = {
    'origin': 'lower',
}

# plots:
# -> upper panel: temperature
fig, ax = plt.subplots(1, 1, figsize=figsize, subplot_kw={'projection':imap.wcs})
if not args.axis:
    ax.axis('off')
    plt.tight_layout(h_pad=0.5)
else:
    ax.tick_params(axis='x', colors='white', which='both', labelcolor='black')
    ax.tick_params(axis='y', colors='white', which='both', labelcolor='black')
    ax.set_aspect('equal')
    for side in ['left','right','top','bottom']:
        ax.spines[side].set_visible(True)
        ax.spines[side].set_color('white')
    plotstyle.setup_axis(ax, nticks=[10,5])
    ax.set_ylabel("$b$")
    ax.set_xlabel('$l$')
    plt.tight_layout(h_pad=0.1)

# polarization angle plot
# reload imap to get the original resolution
# seed = enmap.rand_gauss(imap[0].shape, imap.wcs)
# seed = enmap.smooth_gauss(seed, 0.5*u.arcmin*u.fwhm)
seed = None
imap = enmap.smooth_gauss(imap, 5*u.arcmin*u.fwhm)
P    = np.sum(imap[1:]**2,axis=0)**0.5
if not op.exists(args.texture) or args.force:
    theta = lib.Bangle(imap[1], imap[2], toIAU=True)
    # no need to add for LIC pi/2
    texture = lib.LIC_texture(theta, length=0.1, seed=seed, contrast=True)
    np.save(args.texture, texture)
else:
    texture = np.load(args.texture)
# boost contrast
curve = lambda x: 1/(1+np.exp(-(x-0.5)))
texture = curve(texture)
alpha = np.min([np.ones_like(texture), 1.2*(P/P.max())**0.7],axis=0)
textures = np.stack([np.ones_like(texture)*alpha]*3+[0.6*texture], axis=2)
ax.imshow(tmap, **popts)
ax.imshow(textures, origin='lower')


ofile = op.join(args.odir, args.oname)
print("Writing:", ofile)
plt.savefig(ofile, bbox_inches='tight')
