"""make a contour plot of noise level in each map. """

from common import *
import plotstyle
from pixell import utils as u
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

# parser defined in common
parser.add_argument("--smooth", type=float, default=0)
parser.add_argument("--figsize", default=None)
args = parser.parse_args()

# define area of interests
box = boxes[args.area]
# load data
freqs = ['f090','f150','f220']
# temperature
ivars = [load_ivar(filedb[f]['coadd_ivar'], box=box, fcode=f, mJy=False) for f in freqs]

if args.figsize: figsize=eval(args.figsize)
else: figsize=None
# start plotting
fig, axes = plt.subplots(3, 1, figsize=figsize, sharex=True)
opts = {
    'origin': 'lower',
    'cmap': args.cmap,
    'extent': box2extent(box)/np.pi*180,
    'vmin': args.min,
    'vmax': args.max
}

for i in range(3):
    ax = axes[i]
    nlev = ivars[i][0]**-0.5
    nlev[np.isinf(nlev)] = 0
    print(f"{freqs[i]}: {np.median(nlev[nlev!=0]):.2f} uK arcmin")
    im = ax.imshow(nlev/2, **opts)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='3%', pad=0.1)
    if i==1: fig.colorbar(im, cax=cax, orientation='vertical').set_label("Noise level [$\mu$K arcmin]")
    else: fig.colorbar(im, cax=cax, orientation='vertical')
    
axes[-2].set_ylabel('b [deg]')
axes[-1].set_xlabel('l [deg]')
plt.tight_layout(h_pad=-0.8)
ofile = op.join(args.odir, args.oname)
print("Writing:", ofile)
plt.savefig(ofile, bbox_inches='tight')
