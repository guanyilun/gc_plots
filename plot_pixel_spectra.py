from common import *
import os, os.path as op, numpy as np
import matplotlib.pyplot as plt
import plotstyle
# import matplotlib.ticker as ticker
from matplotlib import colors

def parse(expr):
    if 'np' in expr:
        return eval(expr)
    if ':' in expr:
        start, end, step = expr.split(':')
        return np.arange(float(start), float(end), float(step))
    else:
        return np.array([float(v) for v in expr.split(',')])
    
# the rest of the parser definition are in common.py
parser.add_argument('-l', default="0")
parser.add_argument('-b', default="0")
parser.add_argument('--color', default='x')
parser.add_argument('--use', default='coadd')
parser.add_argument("--title", default=None)
args = parser.parse_args()
if not op.exists(args.odir): os.makedirs(args.odir) 

# load maps and pixel variance
maps = []
ivars = []
freqs = ['f090','f150','f220']
if args.use == 'planck': freqs += ['f350']
for f in freqs:
    maps.append(load_map(filedb[f][args.use], fcode=f)/1e9)
    ivars.append(load_ivar(filedb[f][f'{args.use}_ivar'], fcode=f)*1e18)
# locate pixel and the values in the map
xs = parse(args.l) / 180*np.pi
ys = parse(args.b) / 180*np.pi
# broadcast
if len(xs) == 1: xs = np.ones_like(ys)*xs[0]
if len(ys) == 1: ys = np.ones_like(xs)*ys[0] 
# plot
if args.color == 'x':
    norm = colors.Normalize(vmin=0,vmax=xs.max())
else:
    norm = colors.Normalize(vmin=0,vmax=ys.max())
fig, ax = plt.subplots(1,1)
cmap = plt.get_cmap(args.cmap)
for x, y in zip(xs, ys):
    vals = np.array([m.at([y, x])[0] for m in maps])
    errs = np.array([m.at([y, x])[0] for m in ivars])
    if args.color == 'x': c = x
    else: c = y
    fc = [fcenters[f] for f in freqs]
    ax.errorbar(fc, vals, yerr=errs**-0.5, color=cmap(norm(np.abs(c))), alpha=0.5)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylabel('Total Intensity [MJy/sr]')
ax.set_xlabel('frequency [GHz]')
if args.title: plt.title(args.title)

ofile = op.join(args.odir, args.oname)
print("Writing:", ofile)
plt.savefig(ofile)
