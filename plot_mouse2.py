"""Plot mouse PWN: v2

2 panel plot: f090 only -- brainstorm below

f090 total intensity, with magnetic field orientation plotted on top
f090 polarization fraction, with contours showing polarization fraction
 
"""

from common import *
import matplotlib.pyplot as plt
import plotstyle, lib

# parser defined in common
parser.add_argument("--tmin", type=float, default=None)
parser.add_argument("--tmax", type=float, default=None)
parser.add_argument("--pmin", type=float, default=None)
parser.add_argument("--pmax", type=float, default=None)
parser.add_argument("--smooth", type=float, default=None)
parser.add_argument("--freq", default='f090')
parser.add_argument("--title", default="Mouse PWN")
parser.add_argument("--figsize", default=None)
args = parser.parse_args()
if not op.exists(args.odir): os.makedirs(argsn.odir)
if args.figsize: figsize = eval(args.figsize)
else: figsize = None
# define box of interests
box = boxes[args.area]

# load map
imap = load_map(filedb[args.freq]['coadd'], fcode=args.freq, box=box) / 1e9 # MJy/sr
ivar = load_ivar(filedb[args.freq]['coadd_ivar'], fcode=args.freq, box=box) * 1e18 # MJy/sr
# optionally smooth the map
if args.smooth:
    imap_sm = enmap.smooth_gauss(imap, args.smooth*u.fwhm*u.arcmin)
    ivar_sm = enmap.smooth_gauss(ivar, args.smooth*u.fwhm*u.arcmin)
    s       = sfactor(args.freq, args.smooth)
else:
    s       = 1

# initialize figure
fig, axes = plt.subplots(1, 2, figsize=figsize, subplot_kw={'projection': imap.wcs})

##############
# left panel #
##############
ax   = plotstyle.setup_axis(axes[0], nticks=[5,5])
opts = {
    'cmap': 'planck_half',
    'vmin': args.tmin,
    'vmax': args.tmax,
    'interpolation': 'nearest'
}
# background: total intensity
im = ax.imshow(imap_sm[0], **opts)
cax = plotstyle.add_colorbar(fig, ax)
fig.colorbar(im, cax=cax).set_label(texify("Total intensity [MJy/sr]"))
ax.set_xlabel(r"$l$")
ax.set_ylabel(r"$b$")
# foreground: magnetic field orientation

theta = lib.Bangle(imap_sm[1], imap_sm[2], toIAU=True)
theta+= np.pi/2  # this gets the B-field angle corrected
# x- and y-components of magnetic field
Bx = np.cos(theta)
By = np.sin(theta)
# calculate polarization error as masking
P     = np.sum(imap[1:]**2,axis=0)**0.5
P_err = lib.P_error(imap, ivar*s**2)
P_snr = P / P_err
mask  = P_snr < 3
import ipdb; ipdb.set_trace()
cmap_ = plt.get_cmap('binary')  # actual cmap doesn't matter
color = cmap_(np.ones_like(P))
color[ mask,-1] = 0.3
color[~mask,-1] = 1
color = color.reshape(color.shape[0]*color.shape[1],4)
ax.quiver(Bx,By,pivot='middle', headlength=0, headaxislength=0, color=color)

# right panel
ax = plotstyle.setup_axis(axes[1], nticks=[5,5], xticks=True, yticks=False)
ax.tick_params(axis='x', colors='white', which='both', labelcolor='black')
ax.tick_params(axis='y', colors='white', which='both', labelcolor='black')
for side in ['left','right','top','bottom']:
    ax.spines[side].set_visible(True)
    ax.spines[side].set_color('white')
ax.set_xlabel('$l$')
opts = {
    'cmap': 'magma',
    'vmin': args.pmin,
    'vmax': args.pmax,
    'interpolation': 'nearest'
}
P   = np.sum(imap_sm[1:]**2, axis=0)**0.5
im  = ax.imshow(P, **opts)
cax = plotstyle.add_colorbar(fig, ax)
fig.colorbar(im, cax=cax).set_label(texify("Polarization intensity [MJy/sr]"))
if args.title: plt.suptitle(texify(args.title), fontsize=16)
fig.subplots_adjust(hspace=0)
ofile = op.join(args.odir, args.oname)
print("Writing:", ofile)
plt.savefig(ofile, bbox_inches='tight')

