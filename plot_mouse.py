"""Plot mouse PWN"""

import argparse, os, os.path
import numpy as np

import lib
from common import *
import plotstyle
from pixell import enplot, utils as u, enmap
from matplotlib import pyplot as plt
from scipy.optimize import leastsq

# parameters
fwhms = {
    'f090': 2.05,
    'f150': 1.40,
    'f220': 0.98
}

parser = argparse.ArgumentParser()
parser.add_argument("-o","--odir", default="plots")
parser.add_argument("--oname", default="mouse.pdf")
parser.add_argument("--smooth", type=float,default=2)
parser.add_argument("--dust-removal-f090", type=float, default=0.11)
parser.add_argument("--dust-removal-f150", type=float, default=0.147)
parser.add_argument("--cmap",default='planck')
args = parser.parse_args()
if not op.exists(args.odir): os.makedirs(args.odir)

box = np.array([[-0.9,-0.65],[-0.7,-0.8]]) / 180*np.pi

imap_f090 = load_map(filedb['f090']['coadd'])
imap_f150 = load_map(filedb['f150']['coadd'])
imap_f220 = load_map(filedb['f220']['coadd'])

# form dust template at each frequency
fmap = enmap.fft(imap_f220)
l = imap_f220.modlmap()
bmap_f090 = np.exp(-0.5*l**2*(fwhms['f090']*u.fwhm*u.arcmin)**2)
bmap_f150 = np.exp(-0.5*l**2*(fwhms['f150']*u.fwhm*u.arcmin)**2)
bmap_f220 = np.exp(-0.5*l**2*(fwhms['f220']*u.fwhm*u.arcmin)**2)
dust_f090 = enmap.ifft(fmap * (bmap_f090 / np.maximum(bmap_f220, 1e-3))).real
dust_f150 = enmap.ifft(fmap * (bmap_f150 / np.maximum(bmap_f220, 1e-3))).real

def preprocess_T(imap, dust_removal=None, box=None, auto_adj=True):
    omap = enmap.submap(imap[0], box=box)
    if dust_removal is not None:
        dust = enmap.submap(dust_removal[0], box=box)
        if auto_adj:
            par, _ = leastsq(lambda x: np.ravel(omap-x*dust), x0=1)
            factor = par[0]
        else: factor = 1
        print("factor:", factor)
        omap -= dust * factor
    return omap
def preprocess_P(imap, dust_removal=None, smooth=None, box=None):
    omap = np.sum(imap[1:]**2,axis=0)**0.5
    if dust_removal is not None:
        omap -= np.sum(dust_removal[1:]**2,axis=0)**0.5
    if smooth is not None:
        omap = enmap.smooth_gauss(omap, smooth*u.fwhm*u.arcmin)
    return enmap.submap(omap, box=box)

# preprocess
I_f090 = preprocess_T(imap_f090, dust_removal=args.dust_removal_f090*dust_f090, box=box)
I_f150 = preprocess_T(imap_f150, dust_removal=args.dust_removal_f150*dust_f150, box=box)
# I_f220 = preprocess_T(imap_f220, box=box)
P_f090 = preprocess_P(imap_f090, dust_removal=args.dust_removal_f090*dust_f090,
                      smooth=args.smooth, box=box)
P_f150 = preprocess_P(imap_f150, dust_removal=args.dust_removal_f150*dust_f150,
                      smooth=args.smooth, box=box)
# P_f220 = preprocess_P(imap_f220, smooth=args.smooth, box=box)

popts = {
    'origin': 'lower',
    'cmap': args.cmap
}
names = ['f090 - f220','f150 - f220']
maps = ['T','P']
fig, axes = plt.subplots(2,2, figsize=(7,9))
# axes[0,0] = plt.imshow(diff_func(0.106).upgrade(10), vmin=1000, vmax=6000)
# axes[0,0].imshow(I_f090, vmin=2500, vmax=7000, **popts)
axes[0,0].imshow(I_f090, vmin=-500, vmax=3000, **popts)
# axes[0,1].imshow(I_f150, vmin=1000, vmax=4000, **popts)
axes[0,1].imshow(I_f150, vmin=-500, vmax=2500, **popts)
# axes[0,2].imshow(I_f220, vmin=1000, vmax=2.9e4, **popts)
axes[1,0].imshow(P_f090, vmin=0, vmax=300, **popts)
axes[1,1].imshow(P_f150, vmin=0, vmax=250, **popts)
# axes[1,2].imshow(P_f220, vmin=0, vmax=900, **popts)
for i in range(2):
    for j in range(2):
        ax = axes[i,j]
        if i == 0:
            ax.text(0.3, 1.02, names[j],
                   verticalalignment='bottom', horizontalalignment='left',
                   transform=ax.transAxes, fontsize=14)
        if j == 0:
            ax.text(-0.1, 0.5, maps[i], rotation=90,
                   verticalalignment='bottom', horizontalalignment='left',
                   transform=ax.transAxes, fontsize=14)
        axes[i,j].axis('off')
plt.suptitle("Mouse", y = 1.02, fontsize=16)
plt.tight_layout(h_pad=0.1,w_pad=0.1)
ofile = op.join(args.odir, args.oname)
print("Writing:", ofile)
plt.savefig(ofile, bbox_inches='tight')
