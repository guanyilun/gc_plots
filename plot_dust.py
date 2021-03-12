"""This script focuses on a dust dominated region and plot the spectral index"""

import os, os.path as op
import numpy as np
import matplotlib.pyplot as plt
from common import *
import plotstyle
import lib

# parser defined in common
parser.add_argument("--beam-match", action='store_true')
args = parser.parse_args()
if not op.exists(args.odir): os.makedirs(args.odir)

# define area of interest
box = boxes[args.area]

# load data from f150 and f220, since we are interested in dusts
imap_f150 = load_map(filedb['f150']['coadd'], fcode='f090')
imap_f220 = load_map(filedb['f220']['coadd'], fcode='f220')
# match the beam
if args.beam_match: imap_f220 = lib.beam_match(imap_f220, 'f150', 'f220')
# calculate spectral index
bmap = lib.calc_beta(imap_f150, imap_f220, fcenters['f150'], fcenters['f220'])
bmap = enmap.submap(bmap, box)
# initialize figure
fig = plt.figure()
ax = fig.add_subplot(111)
opts = {
    'origin': 'lower',
    'extent': box2extent(box)/np.pi*180,
    'cmap': args.cmap,
    'vmin': args.min,
    'vmax': args.max
}
im = ax.imshow(bmap, **opts)
ax.set_xlabel('l [deg]')
ax.set_ylabel('b [deg]')
plt.colorbar(im).set_label(r"Spectral index $\beta$")
ofile = op.join(args.odir, args.oname)
print("Writing:", ofile)
plt.savefig(ofile, bbox_inches='tight')
plt.clf()
import seaborn as sns
sns.distplot(np.ravel(bmap), hist=True, kde=True, 
             bins=np.linspace(1.2,2.2,21), color = 'k', 
             hist_kws={'edgecolor':'black'},
             kde_kws={'linewidth': 4})
plt.xlabel(r'Spectral Index $\beta$')
plt.ylabel('Count')
plt.savefig(ofile.replace('.p','_hist.p'), bbox_inches='tight')
