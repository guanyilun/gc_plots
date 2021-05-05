"""This script aims to prepare the data release to have the right headers.
etc. The header information can be found here:
https://lambda.gsfc.nasa.gov/common/fits_header.cfm?fitsfile=%2Fdata%2Fsuborbital%2FACT%2FACT%5Fdr5%2Fmaps%2Fact%5Fdr5%2E01%5Fs08s18%5FAA%5Ff090%5Fdaynight%5Fivar%2Efits
"""

from astropy.io import fits
import argparse, os, os.path as op
import numpy as np
from common import *

parser = argparse.ArgumentParser()
parser.add_argument("--odir", help='data release directory')
parser.add_argument("--tag", help='data release tag', default='sr-mw1')
args = parser.parse_args()
if not op.exists(args.odir): os.makedirs(args.odir)

# loop over maps and modify headers for each
freqs = ['f090', 'f150', 'f220']
for f in freqs:
    # load map and ivar
    imap_ifile = filedb[f]['coadd']
    ivar_ifile = filedb[f]['coadd_ivar']
    # define output filename
    tag = args.tag.replace('-','_')
    imap_ofile = op.join(args.odir, f'act_planck_{tag}_s19_{f}_map.fits')
    ivar_ofile = op.join(args.odir, f'act_planck_{tag}_s19_{f}_ivar.fits')
    # 1. convert to IAU: flip sign U -> -U
    imap = enmap.read_map(imap_ifile)
    imap[2] *= -1
    print("Writing:", imap_ofile)
    enmap.write_map(imap_ofile, imap)
    # 2. ivar drop QU covariance if it exists
    ivar = enmap.read_map(ivar_ifile)
    if len(ivar.shape) == 3:
        print("Writing:", ivar_ofile)
        enmap.write_map(ivar_ofile, ivar)
    else:
        ivar = np.stack([ivar[0,0],ivar[1,1],ivar[2,2]], axis=0)
        print("Writing:", ivar_ofile)
        enmap.write_map(ivar_ofile, ivar)
    # 3. update headers
    for name, fname in zip(['map', 'ivar'], [imap_ofile, ivar_ofile]):
        hdul = fits.open(fname)
        hdr = hdul[0].header
        del hdr["WCSAXES"]
        if 'RADESYS' in hdr: del hdr['RADESYS']
        if 'LONPOLE' in hdr: del hdr['LONPOLE']
        if 'LATPOLE' in hdr: del hdr['LATPOLE']
        hdr.set('CTYPE1', 'GLON-CAR', 'Galactic longitude, plate caree projection')
        hdr.set('CTYPE2', 'GLAT-CAR', 'Galactic latitude, plate caree projection')
        # hdr.set('LONPOLE',  0.0, "[deg] Native longitude of galactic pole")
        # hdr.set('LATPOLE', 90.0, "[deg] Native latitude of galactic pole")
        hdr.insert(8,  ('CRPIX3', 1.0, 'Pixel coordinate of reference point, Stokes T'))
        hdr.insert(11, ('CDELT3', 1.0, 'Stokes component increment'))
        hdr.insert(14, ('CUNIT3', '', ''))
        hdr.insert(17, ('CTYPE3', 'STOKES', 'I,Q,U Stokes components'))
        hdr.insert(20, ('CRVAL3', 1.0, 'Coordinate value at reference point, Stokes T'))
        hdr.append(('COORDSYS', 'GALACTIC', 'Galactic coordinate system'))
        hdr.append(('COMMENT', '', ''), end=True)
        hdr.append(('COMMENT', '*** ACT dataset keys ***', ''), end=True)
        hdr.append(('COMMENT', '', ''), end=True)
        hdr.append(('TELESCOP', 'act+planck', ''), end=True)
        hdr.append(('INSTRUME', 'actpol+planckHFI', ''), end=True)
        hdr.append(('RELEASE', args.tag, 'Data release tag'), end=True)
        hdr.append(('SEASON', 's19', 'Observation season'), end=True)
        hdr.append(('PATCH', 'galactic_18h', 'Survey patch'), end=True)
        hdr.append(('FREQ', f, 'Frequency tag'), end=True)
        hdr.append(('ACTTAGS', 'night,prelim', 'ACT data quality class'), end=True)
        hdr.append(('COMMENT', '', ''), end=True)
        hdr.append(('COMMENT', '*** ACT product keys ***', ''), end=True)
        hdr.append(('COMMENT', '', ''), end=True)
        hdr.append(('POLCCONV', 'IAU', 'Polarization convention'), end=True)
        if name == 'map':
            hdr.append(('BUNIT', 'uK', 'Physical (pixel) units'), end=True)
            hdr.append(('EXTNAME', 'FREQ-MAP', 'Extension name'), end=True)
        if name == 'ivar':
            hdr.append(('BUNIT', 'uK^-2', 'Physical (pixel) units'), end=True)
            hdr.append(('EXTNAME', 'IVAR-MAP', 'Extension name'), end=True)
        hdr.append(('FILENAME', op.basename(fname), 'FITS filename'), end=True)
        print("Writing headers:", fname)
        hdul.writeto(fname, overwrite=True, checksum=True)
