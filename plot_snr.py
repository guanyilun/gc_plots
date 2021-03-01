"""This script plots the signal to noise ratio"""

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
parser.add_argument('--min', type=float, default=0)
parser.add_argument('--max', type=float, default=30)
parser.add_argument('--pol', help='use polarization intensity', action='store_true')
parser.add_argument('--log', help='plot log10(snr) instead', action='store_true')
parser.add_argument('--smooth', type=float, default=0)
parser.add_argument('--save-mask', action='store_true')
args = parser.parse_args()
if not op.exists(args.odir): os.makedirs(args.odir)

# box of interests
box = boxes[args.area]

# load maps
for fcode in ['f090', 'f150','f220']:
    # load map and inverse variance map
    # imap = load_map(filedb[fcode]['coadd'], fcode=fcode, mJy=True)
    # ivar = load_ivar(filedb[fcode]['coadd'], fcode=fcode, mJy=True)
    imap = load_map(filedb[fcode]['coadd'], fcode=fcode, mJy=False)
    ivar = load_ivar(filedb[fcode]['coadd_ivar'], fcode=fcode, mJy=False)
    # ugly hack: e.g. f220 ivar is actually div which has shape
    # 3x3xpixells what we really need is the diagonal components
    if len(ivar.shape)==4:
        ivar_correct = enmap.zeros(imap.shape, imap.wcs)
        ivar_correct[0] = ivar[0,0]
        ivar_correct[1] = ivar[1,1]
        ivar_correct[2] = ivar[2,2]
        del ivar
        ivar = ivar_correct
    if args.smooth > 0:
        imap = enmap.smooth_gauss(imap, args.smooth*u.fwhm*u.arcmin)
        ivar = enmap.smooth_gauss(ivar, args.smooth*u.fwhm*u.arcmin)
    if not args.pol:
        snr = imap[0] * ivar[0]**0.5
    else:
        perr = lib.P_error(imap, ivar, method=1)
        snr = np.sum(imap[1:]**2,axis=0)**0.5/perr
    if args.pol: ofile = op.join(args.odir, f'snr_{fcode}_pol.pdf')
    else: ofile = op.join(args.odir, f'snr_{fcode}.pdf')
    # plot
    opts = {
        'origin': 'lower',
        'cmap': 'magma',
        'extent': box2extent(box)/np.pi*180,
        'vmin': args.min,
        'vmax': args.max
    }
    if args.log: plt.imshow(np.log10(snr), **opts)
    else: plt.imshow(snr, **opts)
    plt.xlabel('l [deg]')
    plt.ylabel('b [deg]')
    if not args.pol:
        plt.colorbar(shrink=0.55).set_label('Total intensity S/N')
    else:
        plt.colorbar(shrink=0.55).set_label('Polarization intensity S/N')
    print("Writing:", ofile)
    plt.savefig(ofile, bbox_inches='tight')
    plt.clf()
