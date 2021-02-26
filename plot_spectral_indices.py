"""This script aims to plot a map of spectral indices"""

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
args = parser.parse_args()

if not op.exists(args.odir): os.makedirs(args.odir)

# first load maps from different frequencies
# define a box of interests
box = np.array([[-1,2],[1,-2]]) / 180*np.pi

imap_f090 = load_map(filedb['f090']['coadd'], fcode='f090', box=box)
imap_f150 = load_map(filedb['f150']['coadd'], fcode='f150', box=box)
imap_f220 = load_map(filedb['f220']['coadd'], fcode='f220', box=box)

# calculate the ratio between f090 and f150
rmap = imap_f090 / imap_f150

# fiducial synchrotron spectral index

# suppose synchrontron radiation has a power law behavior
#
#   T(\nu) = T(\nu_0) (\nu / \nu_0)^\beta_s,
# 
# then we have,
# 
#   T(90GHz) / T(150GHz) ~ (90 / 150)^\beta_s. 
#
# With a fiducial value of \beta_s, which we take
# from Choi and Page (2015).
beta = -3.1
# we can plot a map of spectral index with
#
#   \beta_s = log(T(90GHz)/T(150GHz)) / log(90/150)
#
# and we want to plot the spectral index normalized to \beta_s = 1
# hence,
omap = np.log(imap_f090[0]/imap_f150[0])/np.log(100/143)
omap /= beta

opts = {
    'origin': 'lower',
    'cmap': 'planck',
    'extent': box2extent(box)/np.pi * 180,
}
plt.imshow(omap, **opts)
plt.colorbar(shrink=0.55).set_label(r'$\beta/\beta_s$')
plt.xlabel('l [deg]')
plt.ylabel('b [deg]')
ofile = op.join(args.odir, args.oname)
print("Writing:", ofile)
plt.savefig(ofile, bbox_inches='tight')
