"""LIC version

"""


from common import *
import matplotlib.pyplot as plt
import lib
import plotstyle
from matplotlib import colors

# parser defined in common
parser.add_argument("--back", default="external/meerkat/MeerKAT_radio_bubbles.fits")
parser.add_argument("--title", default="Sgr A^*")
parser.add_argument("--figsize", default=None)
parser.add_argument("--smooth", type=float, default=None)
parser.add_argument("--freq", default='f090')
parser.add_argument("--mask", action='store_true')
parser.add_argument("--min2", type=float, default=None)
parser.add_argument("--max2", type=float, default=None)
# parser.add_argument("-L", default=10)
args = parser.parse_args()
if not op.exists(args.odir): os.makedirs(args.odir)
box = boxes[args.area]

# load coadd map
imap = load_map(filedb[args.freq]['coadd'], fcode=args.freq, box=box) / 1e9
ivar = load_ivar(filedb[args.freq]['coadd_ivar'], fcode=args.freq, box=box) * 1e18

# load external data
irmap = enmap.read_map(args.back)
irmap = enmap.submap(irmap, box=box)
imap  = imap.project(irmap.shape, irmap.wcs)
# smooth if necessary
if args.smooth:
    imap = enmap.smooth_gauss(imap, args.smooth*u.fwhm*u.arcmin)
    ivar = enmap.smooth_gauss(ivar, args.smooth*u.fwhm*u.arcmin)
    s    = sfactor(args.freq, args.smooth)
else:
    imap = imap
    ivar = ivar
    s    = 1

# setup two panel plot
if args.figsize: figsize=eval(args.figsize)
else: figsize=None
Y, X = imap.posmap()/np.pi*180
opts = {
    'cmap': args.cmap,
    'norm': colors.LogNorm(vmin=1e-5, vmax=1e-2),    
}
fig = plt.figure(figsize=figsize)

ax = plt.subplot(111, projection=irmap.wcs)
plotstyle.setup_axis(ax, nticks=[5,5], fmt=None)
im = ax.imshow(irmap, **opts)
# polarization angle plot
# reload imap to get the original resolution
if not op.exists('texture.npy'):
    theta = lib.Bangle(imap[1], imap[2], toIAU=True)
    # no need to add for LIC pi/2
    texture = lib.LIC_texture(theta, L=100)
    np.save("texture.npy", texture)
else:
    texture = np.load("texture.npy")

ax.imshow(texture, origin='lower', cmap='binary', alpha=0.8)
ax.set_xlabel('$l$')
ax.set_ylabel('$b$')
ax.set_aspect('equal')
# colorbar
cax = plotstyle.add_colorbar(fig, ax)
fig.colorbar(im, cax=cax).set_label(texify("Total Intensity [MJy/sr]"), fontsize=12)

plt.subplots_adjust(hspace=0)
if args.title: plt.suptitle(texify(args.title), fontsize=16)
ofile = op.join(args.odir, args.oname)
print("Writing:", ofile)
plt.savefig(ofile, bbox_inches='tight')
