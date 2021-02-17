import numpy as np
import lic

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

def LIC_texture(Bangle,length=0.25):
        '''
        Computes the LIC texture to be overplotted on an image using LicPy.
        Input
        Bangle : plane-of-sky magnetic field angles on the domain 0 and pi [radians]
        length : fraction of image length to compute LIC across            [default=0.25]
        Output
        texture : LIC texture                                              [float]
        '''

        # LIC measures angles from the horizontal while IAU polarization angles are 
        # measured from the vertical; this translates the magnetic field angle accordingly
        Bangle += (np.pi/2.)

        # x- and y-components of magnetic field
        b_x = np.sin(Bangle)
        b_y = np.cos(Bangle)

        # length scale; typically 25% of image size but can be adjusted
        L_z = np.shape(Bangle)
        L   = int(length*L_z[0])

        # run lic to get texture
        texture = lic.lic(b_x,b_y,length=L)

        return texture
