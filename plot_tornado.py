"""This script makes the tornado plot, it will be based of
plot_source2 script, with some special tuning for tornado
objects.

"""

import plotstyle
from common import *
import matplotlib.pyplot as plt

# parser defined in common
parser.add_argument("--tmin", type=float, default=None)
parser.add_argument("--tmax", type=float, default=None)
parser.add_argument("--pmin", type=float, default=None)
parser.add_argument("--pmax", type=float, default=None)
parser.add_argument("--tcmap", default="planck_half")
parser.add_argument("--pcmap", default="magma")
parser.add_argument("--freq", default="f090")
parser.add_argument("--title", default='Tornado')
parser.add_argument("--figsize", default=None)
args = parser.parse_args()
if not op.exists(args.odir): os.makedirs(args.odir)
if args.figsize: figsize = eval(args.figsize)
else: figsize = None

# define area of interests
box = boxes[args.area]

# load data
imap = load_map(filedb[args.freq]['coadd'], box=box, fcode=args.freq) / 1e9

# 2 panel plot
fig, axes = plt.subplots(1,2,figsize=figsize,subplot_kw={'projection': imap.wcs})
# fig = plt.figure(figsize=figsize)

# left: total intensity
opts = {
    'vmin': args.tmin,
    'vmax': args.tmax,
    'cmap': args.tcmap,
    'interpolation': 'nearest'
}
ax = plotstyle.setup_axis(axes[0], nticks=[5,5])
im = ax.imshow(imap[0], **opts)
ax.coords[0].set_axislabel(r"$l$")
ax.coords[1].set_axislabel(r"$b$")
# by default, AxesDivider makes new axes class based on the parent axes,
# for projection based axes this causes problem, so I had to force a
# non wcs-based Axis here
cax = plotstyle.add_colorbar(fig, ax)
fig.colorbar(im, cax=cax).set_label(texify("Total intensity [MJy/sr]"))

# right: polarization intensity
opts = {
    'vmin': args.pmin,
    'vmax': args.pmax,
    'cmap': args.pcmap,
    'interpolation': 'nearest'
}
ax = plotstyle.setup_axis(axes[1], nticks=[5,5], yticks=False)
P  = np.sum(imap[1:]**2, axis=0)**0.5
im = ax.imshow(P, **opts)
ax.tick_params(axis='x', colors='white', which='both', labelcolor='black')
ax.tick_params(axis='y', colors='white', which='both', labelcolor='black')
ax.set_xlabel('$l$')
for side in ['left','right','top','bottom']:
    ax.spines[side].set_visible(True)
    ax.spines[side].set_color('white')
cax = plotstyle.add_colorbar(fig, ax)
fig.colorbar(im, cax=cax).set_label(texify("Polarization intensity [MJy/sr]"))

# title
if args.title: plt.suptitle(texify(args.title), fontsize=16)
# IO
ofile = op.join(args.odir, args.oname)
print("Writing:", ofile)
plt.savefig(ofile, bbox_inches='tight')

