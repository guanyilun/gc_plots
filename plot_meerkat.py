"""In this script, I aim to make a MeerKAT image
and compare with the coadd multicolor image

"""

from common import *
import matplotlib.pyplot as plt
import lib
import plotstyle
from matplotlib import colors

# parser defined in common
parser.add_argument("--back", default="external/meerkat/MeerKAT_radio_bubbles.fits")
parser.add_argument("--title", default=None)
parser.add_argument("--figsize", default=None)
parser.add_argument("--smooth", type=float, default=None)
parser.add_argument("--colorbar", action='store_true')
parser.add_argument("--log", action='store_true')
parser.add_argument("--axis", action='store_true')
args = parser.parse_args()
if not op.exists(args.odir): os.makedirs(args.odir)
box = boxes[args.area]

# load external data
irmap = enmap.read_map(args.back)
irmap = enmap.submap(irmap, box=box)

# setup two panel plot
if args.figsize: figsize=eval(args.figsize)
else: figsize=None
fig = plt.figure(figsize=figsize)

ax = plt.subplot(111, projection=irmap.wcs)
plotstyle.setup_axis(ax, nticks=[10,5])
if args.log:
    opts = {
        'cmap': 'planck_half',
        'norm': colors.LogNorm(vmin=1e-5, vmax=1e-2)
    }
else:
    opts = {
        'cmap': 'planck_half',
        'vmin': args.min,
        'vmax': args.max
    }
# plotstyle.setup_axis(ax, nticks=[10,5], fmt=None)
irmap[irmap<0] = 1e-6

im = ax.imshow(irmap, **opts)
if args.axis:
    ax.set_xlabel('$l$')
    ax.set_ylabel('$b$')
    ax.set_aspect('equal')
else:
    ax.axis('off')
plt.tight_layout()
# colorbar
if args.colorbar:
    cax = plotstyle.add_colorbar(fig, ax)
    fig.colorbar(im, cax=cax).set_label(texify("Total Intensity [MJy/sr]"), fontsize=12)
if args.title: plt.suptitle(texify(args.title), fontsize=16)
ofile = op.join(args.odir, args.oname)
print("Writing:", ofile)
plt.savefig(ofile, bbox_inches='tight')
