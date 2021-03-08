import argparse, os, os.path as op
import numpy as np
import matplotlib.pyplot as plt
from pixell import enmap, enplot, utils as u
from astropy.visualization import make_lupton_rgb
from matplotlib import colors

from common import *
import lib
import plotstyle

# parameters
fwhms = {
    'f090': 2.05,
    'f150': 1.40,
    'f220': 0.98
}

parser = argparse.ArgumentParser()
parser.add_argument("-o","--odir", default="plots")
parser.add_argument("--oname", default="multifreq.pdf")
parser.add_argument("-Q", help="softening", type=float, default=5)
parser.add_argument("-s", help="stretch", type=float, default=0.2)
parser.add_argument("--beam-match", help="whether to match beam to f090",
                    action="store_true", default=False)
parser.add_argument("--min", type=float, default=1e3)
parser.add_argument("--max", type=float, default=5e4)
parser.add_argument("--box", help="box in deg with format like [[ymin,xmin],[ymax,xmax]]",
                    default=None)
parser.add_argument("--area", default='full')
parser.add_argument("--norm", help="normalization method", type=int, default=1)
parser.add_argument("--pol", help="plot polarization intensity instead", action='store_true')
parser.add_argument("--smooth", help="optionally apply a smoothing kernel", type=float, default=0)
parser.add_argument("--downgrade", help="downgrade the map", type=int, default=1)
parser.add_argument("--snr", help="snr mask", type=float, default=None)
parser.add_argument("--mask-method", help="snr mask method", type=int, default=1)
parser.add_argument("--mask-alpha", help='show masked region with given alpha', type=float, default=1)
parser.add_argument("--save", help="save omap as data file", default=None)
args = parser.parse_args()
if not op.exists(args.odir): os.makedirs(args.odir)

if args.box is not None:
    box = np.array(eval(args.box)) / 180*np.pi
else:
    box = boxes[args.area]  # full view by default
    
imap_f090 = load_map(filedb['f090']['coadd'], box, fcode='f090')
imap_f150 = load_map(filedb['f150']['coadd'], box, fcode='f150')
imap_f220 = load_map(filedb['f220']['coadd'], box, fcode='f220')

rmap_f090 = imap_f090
if not args.beam_match:
    rmap_f150 = imap_f150
    rmap_f220 = imap_f220
else:
    # convert to the same beam as f090
    l = imap_f220.modlmap()
    bmap_f090 = np.exp(-0.5*l**2*(fwhms['f090']*u.fwhm*u.arcmin)**2)
    bmap_f150 = np.exp(-0.5*l**2*(fwhms['f150']*u.fwhm*u.arcmin)**2)
    bmap_f220 = np.exp(-0.5*l**2*(fwhms['f220']*u.fwhm*u.arcmin)**2)
    rmap_f150 = enmap.ifft(enmap.fft(imap_f150) * (bmap_f090 / np.maximum(bmap_f150, 1e-3))).real
    rmap_f220 = enmap.ifft(enmap.fft(imap_f220) * (bmap_f090 / np.maximum(bmap_f220, 1e-3))).real

# decide whether to plot total intensity (imap[0]) or the polarization
# intensity: imap[1]**2+imap[2]**2)**0.5
if not args.pol: fun = lambda x: x[0]
else: fun = lambda x: np.sum(x[1:]**2, axis=0)**0.5

rmap_f090 = fun(rmap_f090)
rmap_f150 = fun(rmap_f150)
rmap_f220 = fun(rmap_f220)

# optionally apply a filter
if args.smooth > 0:
    rmap_f090 = enmap.smooth_gauss(rmap_f090, args.smooth*u.fwhm*u.arcmin)
    rmap_f150 = enmap.smooth_gauss(rmap_f150, args.smooth*u.fwhm*u.arcmin)
    rmap_f220 = enmap.smooth_gauss(rmap_f220, args.smooth*u.fwhm*u.arcmin)

if args.norm == 1:
    # normalization method 1:
    # -> calibrate such that each map has the same variance
    # load mask first to test calculating std in the masked region
    if args.snr is not None:
        # load snr and mask a given ratio
        mask_f090 = load_snr('f090', pol=args.pol, box=box) < args.snr
        mask_f150 = load_snr('f150', pol=args.pol, box=box) < args.snr
        mask_f220 = load_snr('f220', pol=args.pol, box=box) < args.snr
        s_f090 = 1
        s_f150 = np.std(rmap_f150[~mask_f150])/np.std(rmap_f090[~mask_f090])
        s_f220 = np.std(rmap_f220[~mask_f220])/np.std(rmap_f090[~mask_f090])
    else:
        s_f090 = 1
        s_f150 = np.std(rmap_f150)/np.std(rmap_f090)
        s_f220 = np.std(rmap_f220)/np.std(rmap_f090)
elif args.norm == 2:
    # normalization method 2:
    # -> calibrate such that synchrotron appears white
    beta = -3.1
    s_f090 = 1
    s_f150 = (143/100)**beta
    s_f220 = (217/100)**beta
elif args.norm == 3:
    # normalization method 2:
    # -> calibrate such that dust appears white
    beta = 1.59
    s_f090 = 1
    s_f150 = (143/100)**beta
    s_f220 = (217/100)**beta

print(s_f090,s_f150,s_f220)
print("f090:", np.percentile(rmap_f090/s_f090, [25,50,75]))
print("f150:", np.percentile(rmap_f150/s_f150, [25,50,75]))
print("f220:", np.percentile(rmap_f220/s_f220, [25,50,75]))

if args.downgrade > 1:
    rmap_f090 = rmap_f090.downgrade(args.downgrade)
    rmap_f150 = rmap_f150.downgrade(args.downgrade)
    rmap_f220 = rmap_f220.downgrade(args.downgrade)

# make rgb image from these maps
omap = make_lupton_rgb(
    colors.Normalize(vmin=args.min, vmax=args.max)(rmap_f090/s_f090),
    colors.Normalize(vmin=args.min, vmax=args.max)(rmap_f150/s_f150),
    colors.Normalize(vmin=args.min, vmax=args.max)(rmap_f220/s_f220),
    stretch=args.s,
    Q=args.Q,
)

# optionally apply a mask
if args.snr is not None:
    # load snr and mask a given ratio
    mask_f090 = load_snr('f090', pol=args.pol, box=box, downgrade=args.downgrade) < args.snr
    mask_f150 = load_snr('f150', pol=args.pol, box=box, downgrade=args.downgrade) < args.snr
    mask_f220 = load_snr('f220', pol=args.pol, box=box, downgrade=args.downgrade) < args.snr
    # apply the mask by masking maps through alpha
    # first add an alpha channel to the output image
    alpha = np.ones(omap.shape[:-1]+(1,), dtype=float)
    omap = np.concatenate([omap/255, alpha], axis=-1)
    # change the alpha of masked region to args.mask_alpha
    if args.mask_method == 1:
        mask = np.logical_and.reduce([mask_f090, mask_f150, mask_f220])
        omap[mask,3] = args.mask_alpha
    elif args.mask_method == 2:
        omap[mask_f090,0] = 0
        omap[mask_f150,1] = 0
        omap[mask_f220,2] = 0

if args.save: np.save(args.save, omap)

# start plotting
popts = {
    'origin': 'lower',
}
plt.figure()
plt.imshow(omap, **popts)
plt.axis('off')
plt.tight_layout()
ofile = op.join(args.odir, args.oname)
print("Writing:", ofile)
plt.savefig(ofile, bbox_inches='tight')
