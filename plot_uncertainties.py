"""This script is similar to the plot angular dispersion script but it focuses 
on the uncertainties of Q/U instead of the polarization angle. 

"""

import matplotlib.pyplot as plt
from common import *
import plotstyle
from scipy.signal import convolve2d

# parser defined in common
parser.add_argument("--freq", default='f090')
parser.add_argument("--smooth", type=float, default=None)
parser.add_argument("--figsize", default=None)
parser.add_argument("--underlay", default='dQ')
parser.add_argument("--overlay", default='Q')
parser.add_argument("--use", default='coadd')
parser.add_argument("--mjy", action='store_true')
parser.add_argument("--nker",type=int, default=4)
args = parser.parse_args()
if not op.exists(args.odir): os.makedirs(args.odir)
if args.figsize: figsize = eval(args.figsize)
else: figsize = None
fig = plt.figure(figsize=figsize)
ax  = fig.add_subplot(111)

# load map and ivar
box = boxes[args.area]
imap = load_map(filedb[args.freq][args.use], fcode=args.freq, box=box, mJy=args.mjy, cib_monopole=False)
ivar = load_ivar(filedb[args.freq][f'{args.use}_ivar'], fcode=args.freq, box=box, mJy=args.mjy)
if args.mjy:
    imap /= 1e9
    ivar *= 1e18
    
############
# underlay #
############
kernel = np.ones((args.nker,args.nker), dtype=float)
if args.underlay == 'dQ':
    z     = imap[1]
    if args.mjy: label = r'$\log_{10}\sigma_Q$ [MJy/sr]'
    else: label = r'$\log_{10}\sigma_Q$ [$\mu$K arcmin]'
    z_i   = convolve2d(z, kernel, mode='same')
    z_i2  = convolve2d(z**2, kernel, mode='same')
    n     = np.sum(kernel)
    ddof  = 1.5  # most unbiased based on numerical test
    zstd  = ((z_i2/n-(z_i/n)**2) * (n/(n-ddof)))**0.5
    if args.smooth:
        s = sfactor(args.freq, args.smooth)
        zstd /= s
    zstd  = np.log10(zstd)
elif args.underlay == 'dU':
    z     = imap[2]
    if args.mjy: label = r'$\log_{10}\sigma_U$ [MJy/sr]'
    else: label = r'$\log_{10}\sigma_U$ [$\mu$K arcmin]'
    z_i   = convolve2d(z, kernel, mode='same')
    z_i2  = convolve2d(z**2, kernel, mode='same')
    n     = np.sum(kernel)
    ddof  = 1.5  # most unbiased based on numerical test
    zstd  = ((z_i2/n-(z_i/n)**2) * (n/(n-ddof)))**0.5
    if args.smooth:
        s = sfactor(args.freq, args.smooth)
        zstd /= s
    zstd  = np.log10(zstd)    

opts = {
    'origin': 'lower',
    'extent': box2extent(box)/np.pi*180,
    'cmap': args.cmap,
    'vmin': args.min,
    'vmax': args.max
}
im = ax.imshow(zstd, **opts)
ax.set_xlabel('l [deg]')
ax.set_ylabel('b [deg]')
plt.colorbar(im, shrink=0.7).set_label(label)

###########
# overlay #
###########
# we only smooth the overlay
# optionally smooth the map and ivar
if args.smooth:
    imap = enmap.smooth_gauss(imap, args.smooth*u.fwhm*u.arcmin)
    ivar = enmap.smooth_gauss(ivar, args.smooth*u.fwhm*u.arcmin)
    s    = sfactor(args.freq, args.smooth)
    ivar*= s**2

if args.overlay == 'Q':
    zstd   = ivar[1]**-0.5    
    # zstd   = (2*ivar[0])**-0.5
    zstd   = np.log10(zstd)
    levels = np.linspace(np.min(zstd),np.max(zstd),10)
    # levels = None # np.linspace(-3,-1,10)
    # [0,0.017,0.0175,0.0180,0.0185,0.0190,0.020]
    # label  = 'Q [MJy/sr]'
elif args.overlay == 'U':
    zstd   = ivar[2]**-0.5
    zstd   = np.log10(zstd)
    levels = np.linspace(np.min(zstd),np.max(zstd),10)
Y, X = imap.posmap()/np.pi*180
cs   = ax.contour(X, Y, zstd, levels=levels, alpha=1, cmap='planck_half')
plt.clabel(cs, levels, fontsize=6, colors='w', fmt='%.2f')

########
# save #
########
ofile = op.join(args.odir, args.oname)
print("Writing:", ofile)
plt.savefig(ofile, bbox_inches='tight')
