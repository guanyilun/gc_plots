"""This script aims to calculate angular dispersion defined as

S^2 = \sum_N (psi(x) - psi(x+d))^2  where psi_i are taken in a fixed annulus d
"""

import matplotlib.pyplot as plt
from common import *
import plotstyle
from scipy.signal import convolve2d

# parser defined in common
parser.add_argument("--freq", default='f090')
parser.add_argument("--delta", help='delta in annulus in arcmin', type=float, default=2)
parser.add_argument("--ddelta", help='d delta of in arcmin', type=float, default=0.5)
parser.add_argument("--smooth", type=float, default=None)
parser.add_argument("--figsize", default=None)
args = parser.parse_args()
if not op.exists(args.odir): os.makedirs(args.odir)
# load map and ivar
box = boxes[args.area]
imap = load_map(filedb[args.freq]['coadd'], fcode=args.freq, box=box)
ivar = load_ivar(filedb[args.freq][f'coadd_ivar'], fcode=args.freq, box=box)

# build an annulus map
# rmap = imap.modrmap()
# delta, ddelta = args.delta*u.arcmin, args.ddelta*u.arcmin
# annulus = np.abs(rmap-delta) <= ddelta/2
# print(np.sum(annulus))
# annulus = annulus.astype(float)
# S^2 = <psi_i^2>+<psi>^2-2 psi_i psi

# calculate B-field orientations
psi    = lib.Bangle(imap[1], imap[2], toIAU=True)/np.pi*180
psi2   = psi**2
# psi_i  = enmap.ifft(enmap.fft(enmap.fftshift(annulus))*enmap.fft(psi)).real
# psi_i2 = enmap.ifft(enmap.fft(enmap.fftshift(annulus))*enmap.fft(psi2)).real
# calculate dispersion
# n = np.sum(annulus)
# S2 = (psi_i2+psi2-2*psi_i*psi)/np.sum(annulus)
# S2 = (psi_i2/n-(psi_i/n)**2)

# method 2: manual convolve
# k = np.array([
#     [0,0,0,0,0],
#     [0,1,1,1,0],
#     [0,1,1,1,0],
#     [0,1,1,1,0],
#     [0,0,0,0,0],
# ])
k = np.ones((4,4))
psi_i  = convolve2d(psi, k, mode='same')  # edges are not trustworthy
psi_i2 = convolve2d(psi**2, k, mode='same')
n = np.sum(k)
ddof = 1.5  # most unbiased based on numerical test, with mean~1, std~0.18
S2 = (psi_i2/n-(psi_i/n)**2) * (n/(n-ddof))
opts = {
    'origin': 'lower',
    'extent': box2extent(box)/np.pi*180,
    'cmap': args.cmap,
    'vmin': args.min,
    'vmax': args.max
}
if args.figsize: figsize = eval(args.figsize)
else: figsize = None
fig = plt.figure(figsize=figsize)
ax  = fig.add_subplot(111)
im  = ax.imshow(np.sqrt(S2), **opts)
# plt.imshow(enmap.ifft(enmap.fft(enmap.fftshift(annulus))*enmap.fft(psi**2)).real, **opts)
ax.set_xlabel('l [deg]')
ax.set_ylabel('b [deg]')
plt.colorbar(im, shrink=0.7).set_label(r'$\sigma_\psi$ [deg]')
# add a contour plot
# optionally smooth the map and ivar
if args.smooth:
    imap = enmap.smooth_gauss(imap, args.smooth*u.fwhm*u.arcmin)
    ivar = enmap.smooth_gauss(ivar, args.smooth*u.fwhm*u.arcmin)
    s    = sfactor(args.freq, args.smooth)
    ivar*= s**2
# not accurate but close
P_err   = lib.Pangle_error(imap, ivar, deg=True)
label   = r'$\delta \psi$ [deg]'
level   = [0,5,10,15,30,45,60]
Y, X    = imap.posmap()/np.pi*180
cs      = ax.contour(X, Y, P_err, levels=level, alpha=1, cmap='gray')
plt.clabel(cs, level, fontsize=6, colors='k', fmt='%d')
# revert x axis
ofile = op.join(args.odir, args.oname)
print("Writing:", ofile)
plt.savefig(ofile, bbox_inches='tight')
