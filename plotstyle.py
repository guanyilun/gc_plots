"""Everything about plotting"""
import matplotlib as mpl

mpl.rcParams['font.size']=8
mpl.rcParams['figure.dpi']=140
mpl.rcParams['xtick.direction']='in'
mpl.rcParams['ytick.direction']='in'
mpl.rcParams['xtick.top']=True
mpl.rcParams['ytick.right']=True

from pixell import colorize
colorize.mpl_register('planck')
