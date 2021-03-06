import numpy as np
from pixell import enmap
from common import *

# from https://github.com/astrosica/misc-functions/blob/master/functions_lic.py
def Bangle(Q_data,U_data,toIAU=False):
    '''
    Computes the plane-of-sky magnetic field angle. Stokes maps must be in IAU convention!

    Input
    Q_data : Stokes Q map                                                       [float]
    U_data : Stokes U map                                                       [float]
    toIAU  : if True, converts Stokes Q and U maps from COSMO to IAU convention [default=False]
    Output
    Bangle : magnetic field angles between 0 and pi                             [radians]
    '''

    # copy to ensure we're not rotating original data
    Q = np.copy(Q_data)
    U = np.copy(U_data)

    # flip both Q and U to convert from electric field to magnetic field orientation
    Q *= -1.
    U *= -1.

    if toIAU==True:
        # if converting from COSMOS to IAU convention, flip sign of Stokes U
        U *= -1.

    # compute magnetic field angle on the domain [0,pi)
    Bangle = np.mod(0.5*np.arctan2(U,Q), np.pi)

    return Bangle

def LIC_texture(Bangle,length=0.25,L=None, **kwargs):
    '''
    Computes the LIC texture to be overplotted on an image using LicPy.
    Input
    Bangle : plane-of-sky magnetic field angles on the domain 0 and pi [radians]
    length : fraction of image length to compute LIC across            [default=0.25]
    Output
    texture : LIC texture                                              [float]
    '''
    import lic
    # LIC measures angles from the horizontal while IAU polarization angles are 
    # measured from the vertical; this translates the magnetic field angle accordingly
    Bangle += (np.pi/2.)

    # x- and y-components of magnetic field
    b_x = np.sin(Bangle)
    b_y = np.cos(Bangle)

    # length scale; typically 25% of image size but can be adjusted
    if not L:
        L_z = np.shape(Bangle)
        L   = int(length*L_z[0])
    else:
        L   = int(L)

    # run lic to get texture
    texture = lic.lic(b_x,b_y,length=L, **kwargs)

    return texture

# unit conversion
def fPlanck_KCMB_MJysr(data,freq):
    '''
    Converts Planck Stokes maps from K_CMB to MJy/sr.
    See Table 6 from Planck IX (2013).
    '''
    if freq==100:
        fac = 244.1   # MJy/sr / K_CMB
    elif freq==143:
        fac = 371.74  # MJy/sr / K_CMB
    elif freq==217:
        fac = 483.690 # MJy/sr / K_CMB
    elif freq==353:
        fac = 287.450 # MJy/sr / K_CMB
    elif freq==545:
        fac = 58.04   # MJy/sr / K_CMB
    elif freq==857:
        fac = 2.27    # MJy/sr / K_CMB

    data_MJysr = data*fac

    return data_MJysr

def Pangle_error(data, ivar, method=1, deg=True):
    """Polarization angle uncertainties, inputs are imap and ivar.
    parameter `method` denotes the approximation """
    Q, U = data[1], data[2]
    P = np.sqrt(Q**2+U**2)
    if method == 1:
        QQ = ivar[1]**-1
        UU = ivar[2]**-1
        QU = 0
    elif method == 2:
        QQ = 2*ivar[0]**-1
        UU = 2*ivar[0]**-1
        QU = 0
    if deg: factor = 28.65  # 0.5 (in radian) converted to deg
    else:   factor = 0.5    # in radian
    err = factor*np.sqrt(U**2*QQ + Q**2*UU - 2*Q*U*QU) / (P**2)
    return err

def P_error(data, ivar, method=1):
    Q, U = data[1], data[2]
    P = np.sqrt(Q**2+U**2)
    if method == 1:
        QQ = ivar[1]**-1
        UU = ivar[2]**-1
        QU = 0
    elif method == 2:
        QQ = 2*ivar[0]**-1
        UU = 2*ivar[0]**-1
        QU = 0
    else: raise ValueError(f"Unknown method: {method}")
    Perr = np.sqrt((Q**2*QQ + U**2*UU + 2*Q*U*QU)/(P**2))
    return Perr

def calc_beta(map1, map2, f1, f2):
    return np.log(map1[0]/map2[0])/np.log(f1/f2) - 2

def beam_match(imap, f1, f2):
    """f1, f2 are fcodes instead of the actual frequency centers. It
    assumes that the first one (f1) has larger beam, so f2 will be matched
    to it.
    """
    from common import fwhms
    l = imap.modlmap()
    bmap_f1 = np.exp(-0.5*l**2*(fwhms[f1]*u.fwhm*u.arcmin)**2)
    bmap_f2 = np.exp(-0.5*l**2*(fwhms[f2]*u.fwhm*u.arcmin)**2)
    rmap = enmap.ifft(enmap.fft(imap) * (bmap_f1 / np.maximum(bmap_f2, 1e-3))).real
    return rmap
