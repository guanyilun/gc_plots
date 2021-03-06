"""shared io related variables"""

from pixell import enmap
from pixell.utils import arcmin
from pixell import utils as u
import os.path as op, numpy as np
import glob, lib

# Input files
map_dir = "/projects/ACT/yilung/work/galactic_center"
filedb = {
    'f090': {
        'coadd':  op.join(map_dir, 'coadd_corr', 'gc_f090_map.fits'),
        'coadd_ivar':  op.join(map_dir, 'coadd_corr', 'gc_f090_ivar.fits'),
        'planck': op.join(map_dir, 'planck', 'planck_npipe_100_split*_map.fits'),
        'planck_ivar': op.join(map_dir, 'planck', 'planck_npipe_100_split*_ivar.fits'),
        'act': op.join(map_dir, 'act_coadd', 'gc_f090_map.fits'),
        'act_ivar': op.join(map_dir, 'act_coadd', 'gc_f090_div.fits')
    },
    'f150': {
        'coadd':  op.join(map_dir, 'coadd_corr', 'gc_f150_map.fits'),
        'coadd_ivar':  op.join(map_dir, 'coadd_corr', 'gc_f150_ivar.fits'),
        'planck': op.join(map_dir, 'planck', 'planck_npipe_143_split*_map.fits'),
        'planck_ivar': op.join(map_dir, 'planck', 'planck_npipe_143_split*_ivar.fits'),
        'act': op.join(map_dir, 'act_coadd', 'gc_f150_map.fits'),
        'act_ivar': op.join(map_dir, 'act_coadd', 'gc_f150_div.fits')        
    },
    'f220': {
        'coadd':  op.join(map_dir, 'coadd_corr', 'gc_f220_map.fits'),
        'coadd_ivar':  op.join(map_dir, 'coadd_corr', 'gc_f220_ivar.fits'),
        'planck': op.join(map_dir, 'planck', 'planck_npipe_217_split*_map.fits'),
        'planck_ivar': op.join(map_dir, 'planck', 'planck_npipe_217_split*_ivar.fits'),
        'act': op.join(map_dir, 'act_coadd', 'gc_f220_map.fits'),
        'act_ivar': op.join(map_dir, 'act_coadd', 'gc_f220_div.fits')
    },
    'f350': {
        'planck': op.join(map_dir, 'planck', 'planck_npipe_353_split*_map.fits'),
        'planck_ivar': op.join(map_dir, 'planck', 'planck_npipe_353_split*_ivar.fits')
    },
}

def load_map(path, box=None, mJy=True, fcode=None, cib_monopole=True, IAU=False):
    files = glob.glob(path)
    if len(files) == 1:
        imap = enmap.read_map(files[0])
    elif len(files) == 2:
        imap1 = enmap.read_map(files[0])
        imap2 = enmap.read_map(files[1])
        imap  = (imap1 + imap2) / 2
    else: raise ValueError("Unknown format")
    # optionally convert everything to mJy beam^-1
    # numbers taken from Table 6 from Planck IX (2013)
    if mJy:
        if fcode == 'f090':
            imap *= 244.1*1e3
        elif fcode == 'f150':
            imap *= 371.74*1e3
        elif fcode == 'f220':
            imap *= 483.69*1e3
        elif fcode == 'f350':
            imap *= 287.5*1e3
        else: raise ValueError("Unknown fcode")
        # since the numbers below are given in mJy sr^-1, I will only
        # correct cib monopole when data has been converted to mJy
        # see Planck 2018 III table 12
        if cib_monopole:
            if fcode == 'f090':
                imap[0] -= 0.003 * 1e9  # mJy sr^-1
            elif fcode == 'f150':
                imap[0] -= 0.0079 * 1e9
            elif fcode == 'f220':
                imap[0] -= 0.033 * 1e9
            elif fcode == 'f350':
                imap[0] -= 0.13 * 1e9
            else: raise ValueError("Unknown fcode")
    if IAU: imap[2] *= -1
    if box is not None: return imap.submap(box)
    else: return imap

def load_ivar(path, box=None, mJy=True, fcode=None):
    files = glob.glob(path)
    if len(files) == 1:
        ivar = enmap.read_map(files[0])
    elif len(files) == 2:
        ivar1 = enmap.read_map(files[0])
        ivar2 = enmap.read_map(files[1])
        ivar  = (1/4*ivar1**-1 + 1/4*ivar2**-1)**-1
    else: raise ValueError("Unknown format")
    # ugly hack: e.g. f220 ivar is actually div which has shape
    # 3x3xpixells what we really need is the diagonal components
    if len(ivar.shape)==4:
        ivar_correct = enmap.zeros(ivar.shape[1:], ivar.wcs)
        ivar_correct[0] = ivar[0,0]
        ivar_correct[1] = ivar[1,1]
        ivar_correct[2] = ivar[2,2]
        del ivar
        ivar = ivar_correct
    if mJy:
        if fcode == 'f090':
            ivar /= (244.1*1e3)**2
        elif fcode == 'f150':
            ivar /= (371.74*1e3)**2
        elif fcode == 'f220':
            ivar /= (483.69*1e3)**2
        elif fcode == 'f350':
            ivar /= (287.5*1e3)**2
        else: raise ValueError("Unknown fcode")
    if box is not None: return ivar.submap(box)
    else: return ivar

def box2extent(box, pad=0.25*arcmin):
    """convert pixell box to plt extent"""
    return np.array([box[0][1]+pad, box[1][1]-pad, box[0][0]-pad, box[1][0]+pad])

def load_snr(fcode=None, box=None, pol=False, downgrade=1, use='coadd'):
    # load map and ivar
    imap = load_map(filedb[fcode][use], fcode=fcode, box=box)
    ivar = load_ivar(filedb[fcode][f"{use}_ivar"], fcode=fcode, box=box)
    if not pol:  # total intensity
        snr = imap[0] * ivar[0]**0.5
    else:  # polarization intensity
        P = np.sum(imap[1:]**2,axis=0)**0.5
        P_err = lib.P_error(imap, ivar)
        snr = P / P_err
    return snr
    # optionally downgrade
    if downgrade > 1:
        snr = snr.downgrade(downgrade)*downgrade
    if box is not None: return snr.submap(box)
    else: return snr

def sfactor(fcode, fwhm):
    """find signal to noise improvement factor after smoothing, fwhm
    should be in arcmin"""
    return np.sqrt(fwhms[fcode]**2+fwhm**2)/fwhms[fcode]

# common boxes: [[ymin, xmin], [ymax, xmax]]
boxes = {}
boxes['full']  = np.array([[-2,4],[2,-4]]) / 180*np.pi
boxes['half']  = np.array([[-1,2],[1,-2]]) / 180*np.pi
boxes['pilot'] = np.array([[-0.15,1.5],[0.15,-1.5]]) / 180*np.pi
boxes['quat']  = np.array([[-0.5,1],[0.5,-1]]) / 180*np.pi
boxes['trim']  = np.array([[-1,3.8],[1,-3.8]]) / 180*np.pi
boxes['trim2'] = np.array([[-1.2,3.6],[1.2,-3.6]]) / 180*np.pi
boxes['trim3'] = np.array([[-1,2],[0.7,-2.7]]) / 180*np.pi
boxes['gismo'] = np.array([[-0.27,0.92],[0.235,-0.73]]) / 180*np.pi
boxes['saga']  = np.array([[-0.17,0.08],[0.10,-0.20]]) / 180*np.pi
boxes['saga2'] = np.array([[-0.17,0.1],[0.10,-0.20]]) / 180*np.pi
boxes['saga3'] = np.array([[-0.3,0.25],[0.15,-0.20]]) / 180*np.pi  # much larger
boxes['gcra']  = np.array([[-0.3,0.21],[0.15,0.1]]) / 180*np.pi
boxes['mouse'] = np.array([[-0.9,-0.65],[-0.7,-0.8]]) / 180*np.pi
boxes['tndo']  = np.array([[-0.25, -2.15],[0.05, -2.45]]) / 180*np.pi
boxes['dust1'] = np.array([[-0.8,0.4],[-0.5,0.7]]) / 180*np.pi
boxes['dust2'] = np.array([[-0.125,-0.1],[0,-0.2]]) / 180*np.pi
boxes['dust3'] = np.array([[0,3.5],[1,2.5]]) / 180*np.pi
boxes['tlp']   = np.array([[-0.13,0.18],[-0.02,0.025]]) / 180*np.pi  # enlarged
boxes['tlp2']  = np.array([[-0.15,0.2],[0,0.005]]) / 180*np.pi  # enlarged
boxes['brick'] = np.array([[-0.025,0.30],[0.075,0.2]]) / 180*np.pi
boxes['brick2']= np.array([[-0.1,0.40],[0.2,0.1]]) / 180*np.pi
boxes['dust4'] = np.array([[-0.82, 0.38],[-0.76,0.34]]) / 180*np.pi
# boxes['l1.3']  = np.array([[-0.25,1.3], [0.25,1.1]]) / 180*np.pi
boxes['l1.3']  = np.array([[-0.25,1.4], [0.25,1.1]]) / 180*np.pi
# boxes['gcra']  = np.array([[-0.5,0.28],[0.5,0]]) / 180*np.pi
# boxes['tlp']   = np.array([[-0.12,0.18],[-0.04,0.025]]) / 180*np.pi
# boxes['dust1'] = np.array([[0,3.7],[0.2,3.5]]) / 180*np.pi
# boxes['dust1'] = np.array([[0,3.4],[0.6,3]]) / 180*np.pi
# boxes['dust4'] = np.array([[-0.9, 0.7],[-0.35,0.2]]) / 180*np.pi
# boxes['dust4'] = np.array([[-0.85, 0.4],[-0.73,0.32]]) / 180*np.pi

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-o", "--odir", default='plots')
parser.add_argument("-v", "--verbose", action='store_true')
parser.add_argument("--area", default='full')
parser.add_argument("--min", type=float, default=None)
parser.add_argument("--max", type=float, default=None)
parser.add_argument("--cmap", default='planck')
parser.add_argument("--oname", default=None)

# beam parameters
fwhms = {
    'f090': 2.05,
    'f150': 1.40,
    'f220': 0.98
}

fcenters = {
    'f090': 100,
    'f150': 143,
    'f220': 217,
    'f350': 353
}

def texify(text):
    text = text.replace(" ", "~")
    return r"${\rm %s}$" % text
