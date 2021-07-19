"""This script aims to investigate the amount of residue leakages in the maps"""

import argparse, os, os.path as op
import matplotlib as mpl
import plotstyle
from matplotlib import pyplot as plt
from common import *
parser = argparse.ArgumentParser()
parser.add_argument('--odir', default="plots")
parser.add_argument("--area", default='full')
parser.add_argument("--smooth", type=float, default=0)
parser.add_argument("--corr", action='store_true')
parser.add_argument("--oname", default='test_leakage.pdf')
parser.add_argument("--freq", default='f150')
args = parser.parse_args()

if not op.exists(args.odir): os.makedirs(args.odir)

box = boxes[args.area]

# load pa5 f150 maps
if not args.corr:
    imap_pa5 = load_map(op.join(map_dir, 'act', f"map_pa5_{args.freq}_night*.fits"), fcode=args.freq, box=box)/1e9
    imap_pa6 = load_map(op.join(map_dir, 'act', f"map_pa6_{args.freq}_night*.fits"), fcode=args.freq, box=box)/1e9
else:
    imap_pa5 = load_map(op.join(map_dir, 'act_leakage_corr',
                                f"map_fpol_fq_fu_pa5_{args.freq}*_corr.fits"), fcode=args.freq, box=box)/1e9
    imap_pa6 = load_map(op.join(map_dir, 'act_leakage_corr',
                                f"map_fpol_fq_fu_pa6_{args.freq}*_corr.fits"), fcode=args.freq, box=box)/1e9    
# smooth
if args.smooth:
    imap_pa5 = enmap.smooth_gauss(imap_pa5, args.smooth*u.arcmin*u.fwhm)
    imap_pa6 = enmap.smooth_gauss(imap_pa6, args.smooth*u.arcmin*u.fwhm)

# calculate polarized intensity
P_pa5 = np.sum(imap_pa5[1:]**2, axis=0)**0.5
P_pa6 = np.sum(imap_pa6[1:]**2, axis=0)**0.5

diff = np.abs(P_pa6 - P_pa5)
# true sky
imap_true = load_map(filedb[args.freq]['coadd'], box=box, fcode=args.freq)/1e9
P_true = np.sum(imap_true[1:]**2, axis=0)**0.5

fig = plt.figure()
ax = fig.add_subplot(111, projection=imap_pa5.wcs)
im = ax.imshow(diff/P_true, cmap='jet', vmax=1)
# im = ax.imshow(diff/imap_true[0], cmap='magma', vmax=0.02)
ax.set_xlabel('$l$')
ax.set_ylabel('$b$')
cax = plotstyle.add_colorbar(fig, ax, size="2%")
cbar = fig.colorbar(im, cax=cax)

plt.savefig(op.join(args.odir, args.oname), bbox_inches='tight')
