"""This script aims to plot the passbands from Naess et al. (2020),
and show the location of molecular lines on top of them


"""
import argparse, os, os.path as op
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import ticker
import plotstyle
from common import texify

parser = argparse.ArgumentParser()
parser.add_argument("-o","--odir",default="plots")
parser.add_argument("--oname", default="passband.pdf")
args = parser.parse_args()
if not op.exists(args.odir): os.makedirs(args.odir)

# load passbands 
passbands = np.loadtxt("external/act_planck_dr5.01_s08s18_bandpasses.txt")
# the order of passbands is given by
# pa1_f150, pa2_f150, pa3_f090, pa3_f150, pa4_f150, pa4_f220, pa5_f090, pa5_f150, pa6_f090, pa6_f150, ar1_f150 ,ar2_f220, planck_f090, planck_f150 and planck_f220. 

# columns look up table
pb = {
    'f090': {
        'act': [7, 9],
        'plk': [-3],
    },
    'f150': {
        'act': [5, 8, 10],
        'plk': [-2],
    },
    'f220': {
        'act': [6],
        'plk': [-1],
    }
} 

# initialize figure
fig = plt.figure(figsize=(9,3))
ax = fig.add_subplot(111)

freqs = passbands[:,0]
act_band_f090 = np.mean(passbands[:, pb['f090']['act']], axis=1)
plk_band_f090 = np.mean(passbands[:, pb['f090']['plk']], axis=1)
act_band_f090 /= act_band_f090.max()
plk_band_f090 /= plk_band_f090.max()

act_band_f150 = np.mean(passbands[:, pb['f150']['act']], axis=1)
plk_band_f150 = np.mean(passbands[:, pb['f150']['plk']], axis=1)
act_band_f150 /= act_band_f150.max()
plk_band_f150 /= plk_band_f150.max()

act_band_f220 = np.mean(passbands[:, pb['f220']['act']], axis=1)
plk_band_f220 = np.mean(passbands[:, pb['f220']['plk']], axis=1)
act_band_f220 /= act_band_f220.max()
plk_band_f220 /= plk_band_f220.max()

# act_band = act_band_f090 + act_band_f150 + act_band_f220
# plk_band = plk_band_f090 + plk_band_f150 + plk_band_f220
l_act, = ax.plot(freqs[::5], act_band_f090[::5], label=texify('AdvACT'), color='r')
ax.plot(freqs[::5], act_band_f150[::5], color='r')
ax.plot(freqs[::5], act_band_f220[::5], color='r')
l_plk, = ax.plot(freqs[::5], plk_band_f090[::5], label=texify('Planck'), color='k')
ax.plot(freqs[::5], plk_band_f150[::5], color='k')
ax.plot(freqs[::5], plk_band_f220[::5], color='k')
ax.set_ylim(bottom=1e-5)

locator = ticker.MaxNLocator(prune='both', nbins=10)
ax.xaxis.set_major_locator(locator)
plt.xlabel(texify('Frequency [GHz]'))
plt.ylabel(texify('Response'))
ax.grid(linestyle='--')
plt.legend(handles=[l_act, l_plk])

# start to show molecular lines
lines = {
    'CO(1-0)': 115.3,
    'CO(2-1)': 230.5,
    'HCN': 88.6,
}
for l, f in lines.items():
    ax.axvline(f, c='g', linestyle='--')
    ax.text(f, 1.1, l, ha='center', color='g')

ofile = op.join(args.odir, args.oname)
print("Writing:", ofile)
plt.savefig(ofile, bbox_inches='tight')
