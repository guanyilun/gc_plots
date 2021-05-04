"""This script plots the noise level in a region"""

import argparse, os, os.path as op
import numpy as np
from common import *
import lib
import plotstyle
from matplotlib import pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument('-o', '--odir', default='plots')
parser.add_argument('-v', action='store_true', default=False)
parser.add_argument('--area', default='full')
parser.add_argument('--min', type=float, default=None)
parser.add_argument('--max', type=float, default=None)
parser.add_argument('--pol', help='use polarization intensity', action='store_true')
parser.add_argument('--method', help='method used to calculate pol uncertainty', type=int, default=1)
parser.add_argument('--log', help='plot log10(snr) instead', action='store_true')
parser.add_argument('--smooth', help='smoothening fwhm in arcmin', type=float, default=0)
parser.add_argument('--save', help='optionally save the snr fits file in this directory', default=None)
parser.add_argument('--downgrade', type=int, default=1)
parser.add_argument('--oname', default='noise_unnamed.pdf')
parser.add_argument("--freq", default='f220')
parser.add_argument("--title", default=None)
args = parser.parse_args()
if not op.exists(args.odir): os.makedirs(args.odir)

# box of interests
box = boxes[args.area]

# load maps
fcode = args.freq
title = args.title

# load map and inverse variance map
imap = load_map(filedb[fcode]['coadd'], fcode=fcode, mJy=True, box=box)/1e9
ivar = load_ivar(filedb[fcode]['coadd_ivar'], fcode=fcode, mJy=True, box=box)*1e18
# calculate in uK unit, should give the same result apart from
# the cib monopole part
# imap = load_map(filedb[fcode]['coadd'], fcode=fcode, mJy=False)
# ivar = load_ivar(filedb[fcode]['coadd_ivar'], fcode=fcode, mJy=False)

# optionally apply smoothing
if args.smooth > 0:
    imap = enmap.smooth_gauss(imap, args.smooth*u.fwhm*u.arcmin)
    ivar = enmap.smooth_gauss(ivar, args.smooth*u.fwhm*u.arcmin)
    s    = sfactor(args.freq, args.smooth)
else:
    s    = 1
# optionally apply downgrade
if args.downgrade > 1:
    imap = imap.downgrade(args.downgrade)
    ivar = ivar.downgrade(args.downgrade)
# decide whether to look at I or P
if not args.pol:
    noise = ivar[0]**-0.5/s
    label = "Noise (I)"
    if title: label = f'{title}: ' + label
else:
    noise = lib.P_error(imap, ivar*s**2)
    label = "Noise (P)"    
    if title: label = f'{title}: ' + label    

# start plotting
opts = {
    'cmap': 'gray',
    'vmin': args.min,
    'vmax': args.max
}
fig = plt.figure()
ax = plt.subplot(111, projection=imap.wcs)
plotstyle.setup_axis(ax, nticks=[5,5], fmt=None)
im = ax.imshow(noise, **opts)
plt.xlabel('$l$')
plt.ylabel('$b$')

cax = plotstyle.add_colorbar_hpad(ax)
fig.colorbar(im, cax=cax, orientation='horizontal').set_label(
    texify(f"{label} {fcode} [MJy/sr]"), fontsize=12)
cax.xaxis.set_ticks_position('top')
cax.xaxis.set_label_position('top')
# save figure: make sure the output file has proper name
ofile = op.join(args.odir, args.oname)
print("Writing:", ofile)
plt.savefig(ofile, bbox_inches='tight')
# optionally save the snr map as a fits file
if args.save is not None:
    if not args.pol: ofile=f'snr_{fcode}.fits'
    else: ofile = f'snr_{fcode}_pol.fits'
    ofile = op.join(args.save, ofile)
    print("Writing:", ofile)
    enmap.write_map(ofile, snr)
