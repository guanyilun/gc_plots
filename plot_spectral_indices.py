"""This script aims to plot a map of spectral indices
in terms of difference to synchtrotron radiation

"""

import argparse, os, os.path as op
import numpy as np

import lib
from common import *
import plotstyle
from matplotlib import pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument('-o', '--odir', default='plots')
parser.add_argument('-v', action='store_true', default=False)
parser.add_argument("--oname", default="spec.pdf")
parser.add_argument("--area", default='half')
parser.add_argument("--cmap", default='planck')
parser.add_argument("--freq", default="f090")
args = parser.parse_args()

if not op.exists(args.odir): os.makedirs(args.odir)

# first load maps from different frequencies
# define a box of interests
box = boxes[args.area]

imap_f090 = load_map(filedb['f090']['coadd'], fcode='f090', box=box)
imap_f150 = load_map(filedb['f150']['coadd'], fcode='f150', box=box)
imap_f220 = load_map(filedb['f220']['coadd'], fcode='f220', box=box)

# calculate the ratio between f090 and f150
if args.freq == 'f090':
    map1 = imap_f090
    map2 = imap_f150
    f1 = 100
    f2 = 143
    beta = -3.1
    label= r"$\beta - \beta_s$"
elif args.freq == 'f150':
    map1 = imap_f150
    map2 = imap_f220
    f1 = 143
    f2 = 217
    beta = 1.59
    label= r"$\beta - \beta_d$"

else: raise ValueError("Unknown spectral index")
# fiducial synchrotron spectral index

# suppose synchrontron radiation has a power law behavior
#
#   T(\nu) = T(\nu_0) (\nu / \nu_0)^(\beta_s+2),
#
# then we have,
#
#   T(90GHz) / T(150GHz) ~ (90 / 150)^(\beta_s+2).
#
# With a fiducial value of \beta_s, which we take
# from Choi and Page (2015).
#
#   beta = -3.1 (sync)
#   beta = 1.59 (dust)
#
# we can plot a map of spectral index with
#
#   \beta_s = log(T(90GHz)/T(150GHz)) / log(90/150) - 2
#
# and we want to plot the spectral index normalized to \beta_s = 1
# hence, (using Planck frequency)
omap = np.log(map1[0]/map2[0])/np.log(f1/f2) - 2
omap -= beta

opts = {
    'origin': 'lower',
    'cmap': args.cmap,
    'extent': box2extent(box)/np.pi * 180,
}
plt.imshow(omap, **opts)
plt.colorbar(shrink=0.55).set_label(label)
plt.xlabel('l [deg]')
plt.ylabel('b [deg]')
ofile = op.join(args.odir, args.oname)
print("Writing:", ofile)
plt.savefig(ofile, bbox_inches='tight')
