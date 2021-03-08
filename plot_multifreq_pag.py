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
parser.add_argument("--axis", action='store_true')
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

if args.norm == 1:
    # normalization method 1:
    # -> calibrate such that each map has the same variance
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
print("f090:", np.percentile(rmap_f090[0]/s_f090, [25,50,75]))
print("f150:", np.percentile(rmap_f150[0]/s_f150, [25,50,75]))
print("f220:", np.percentile(rmap_f220[0]/s_f220, [25,50,75]))
omap = make_lupton_rgb(
    colors.Normalize(vmin=args.min, vmax=args.max)(rmap_f090[0]/s_f090),
    colors.Normalize(vmin=args.min, vmax=args.max)(rmap_f150[0]/s_f150),
    colors.Normalize(vmin=args.min, vmax=args.max)(rmap_f220[0]/s_f220),
    stretch=args.s,
    Q=args.Q,
)
popts = {
    'origin': 'lower',
}

# Open objects database
import pandas as pd
df_extended = pd.read_csv('extended_sources_from_googledoc.csv')

plt.figure()
ax = plt.subplot(projection=imap_f090.wcs)
plt.imshow(omap, **popts)
if args.axis:
    ax.coords[0].set_format_unit('deg')
    ax.coords[0].set_ticks(number=10)
    ax.coords[1].set_format_unit('deg')
    ax.coords[1].set_ticks(number=10)
    ax.coords[0].display_minor_ticks(True)
    ax.coords[1].display_minor_ticks(True)
    ax.coords[0].set_minor_frequency(10)
    ax.coords[1].set_minor_frequency(10)
    ax.spines['bottom'].set_color('red')
    ax.spines['left'].set_color('red')
    ax.coords[0].set_major_formatter('d.dd')
    ax.coords[1].set_major_formatter('d.dd')
#    plt.grid(color='white', ls='solid')
# plot markings
for j in range(len(df_extended)):
    l, b = df_extended.l.values[j], df_extended.b.values[j]
    # size = df_extended['size'].values[j]
    size=300
    plt.scatter([l], [b], 
                transform=ax.get_transform('world'), 
                edgecolor='white', facecolor='none', lw=0.8, s=size, 
                alpha=0.5, label=df_extended.SourceName.values[j])
    #if df_extended['label_loc'].values[j] == 'top':
    #    deltab = +0.08
    #else:
    #    deltab = 0.08
    dx,dy = df_extended.labelx.values[j], df_extended.labely.values[j]
    plt.text(l+dx, b+dy, df_extended.SourceName.values[j],
             transform=ax.get_transform('world'), color='white',
             fontsize=5, ha='center')
# end plot markings
if not args.axis:
    plt.axis('off')
plt.tight_layout()
ofile = op.join(args.odir, args.oname)
print("Writing:", ofile)
plt.savefig(ofile, bbox_inches='tight')
plt.show()
