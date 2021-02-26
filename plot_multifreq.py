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
args = parser.parse_args()
if not op.exists(args.odir): os.makedirs(args.odir)

if args.box is not None:
    box = np.array(eval(args.box)) / 180*np.pi
else:
    # box = np.array([[-0.27,0.92],[0.235,-0.73]]) / 180*np.pi  # to compare with gismo
    # box = np.array([[-1,2],[1,-2]]) / 180*np.pi  # a narrower view
    box = None  # full view
    

imap_f090 = load_map(filedb['f090']['coadd'], box, mJy=False)
imap_f150 = load_map(filedb['f150']['coadd'], box, mJy=False)
imap_f220 = load_map(filedb['f220']['coadd'], box, mJy=False)

# convert to the same beam as f090
l = imap_f220.modlmap()
bmap_f090 = np.exp(-0.5*l**2*(fwhms['f090']*u.fwhm*u.arcmin)**2)
bmap_f150 = np.exp(-0.5*l**2*(fwhms['f150']*u.fwhm*u.arcmin)**2)
bmap_f220 = np.exp(-0.5*l**2*(fwhms['f220']*u.fwhm*u.arcmin)**2)

rmap_f090 = imap_f090
if not args.beam_match:
    rmap_f150 = imap_f150
    rmap_f220 = imap_f220
else:
    rmap_f150 = enmap.ifft(enmap.fft(imap_f150) * (bmap_f090 / np.maximum(bmap_f150, 1e-3))).real
    rmap_f220 = enmap.ifft(enmap.fft(imap_f220) * (bmap_f090 / np.maximum(bmap_f220, 1e-3))).real

# calibrate such that each map has the same variance
s_f090 = np.std(rmap_f090)
s_f150 = np.std(rmap_f150)
s_f220 = np.std(rmap_f220)
print(s_f090,s_f150,s_f220)
omap = make_lupton_rgb(
    colors.Normalize(vmin=args.min, vmax=args.max)(rmap_f090[0]),
    colors.Normalize(vmin=args.min, vmax=args.max)(rmap_f150[0]*s_f090/s_f150),
    colors.Normalize(vmin=args.min, vmax=args.max)(rmap_f220[0]*s_f090/s_f220),
    stretch=args.s,
    Q=args.Q,
)
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
