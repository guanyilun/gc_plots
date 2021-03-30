"""This script is a simplified version of plot_multifreq_pag
which simply loads the pregenerated image and add annotations

"""

import argparse, os, os.path as op
import numpy as np
import matplotlib.pyplot as plt
from pixell import enmap, enplot, utils as u
from matplotlib import colors
# Open objects database
import pandas as pd

from common import *
import lib
import plotstyle

# parser defined in common
parser.add_argument("-T", default=None)
parser.add_argument("--axis", action='store_true')
parser.add_argument("--figsize", default="(8,4)")
args = parser.parse_args()
if not op.exists(args.odir): os.makedirs(args.odir)
if args.figsize: figsize=eval(args.figsize)
else: figsize=None
# define box of interests
box = boxes[args.area]
# load a map to get the right wcs 
imap = load_map(filedb['f090']['coadd'], box=box, fcode='f090')
# load two maps
tmap = np.load(args.T)

# start plotting
popts = {
    'origin': 'lower',
    # 'extent': box2extent(box)/np.pi * 180
}
fig = plt.figure(figsize=figsize)
ax  = plt.subplot(projection=imap.wcs)
# fig, ax = plt.subplots(1, 1, figsize=figsize, projection=imap.wcs)
ax.imshow(tmap, **popts)

if not args.axis:
    ax.axis('off')
else:
    ax.tick_params(axis='x', colors='white', which='both', labelcolor='black')
    ax.tick_params(axis='y', colors='white', which='both', labelcolor='black')
    for side in ['left','right','top','bottom']:
        ax.spines[side].set_color('white')
    ax.set_aspect('equal')
    ax.set_xlabel('l')
    ax.set_ylabel('b')
    ax.coords[0].set_format_unit('deg')
    ax.coords[0].set_ticks(number=10)
    ax.coords[1].set_format_unit('deg')
    ax.coords[1].set_ticks(number=5)
    ax.coords[0].display_minor_ticks(True)
    ax.coords[1].display_minor_ticks(True)
    ax.coords[0].set_minor_frequency(10)
    ax.coords[1].set_minor_frequency(5)
    ax.coords[0].set_major_formatter('d.dd')
    ax.coords[1].set_major_formatter('d.dd')    

# add annotation
df_extended = pd.read_csv('extended_sources_from_googledoc.csv')
for j in range(len(df_extended)):
    l, b = df_extended.l.values[j], df_extended.b.values[j]
    # size = df_extended['size'].values[j]
    size=300
    plt.scatter([l], [b], 
                transform=ax.get_transform('world'), 
                edgecolor='white', facecolor='none', lw=0.8, s=size, 
                alpha=0.5, label=df_extended.SourceName.values[j])

    dx,dy = df_extended.labelx.values[j], df_extended.labely.values[j]
    plt.text(l+dx, b+dy, df_extended.SourceName.values[j],
             transform=ax.get_transform('world'), color='white',
             fontsize=5, ha='center')
    
ofile = op.join(args.odir, args.oname)
print("Writing:", ofile)
plt.savefig(ofile, bbox_inches='tight')
