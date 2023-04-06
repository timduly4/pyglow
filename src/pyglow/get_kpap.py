from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from datetime import datetime
import numpy as np

from . import generate_kpap

# Fetch geophysical indices (Global variable to make it fast):
GEOPHYSICAL_INDICES = generate_kpap.fetch()

def cubic_interp1d(x0, x, y):
    """
    Interpolate a 1-D function using cubic splines.
      x0 : a float or an 1d-array
      x : (N,) array_like
          A 1-D array of real/complex values.
      y : (N,) array_like
          A 1-D array of real values. The length of y along the
          interpolation axis must be equal to the length of x.

    Implement a trick to generate at first step the cholesky matrice L of
    the tridiagonal matrice A (thus L is a bidiagonal matrice that
    can be solved in two distinct loops).

    additional ref: www.math.uh.edu/~jingqiu/math4364/spline.pdf
    """
    x = np.asfarray(x)
    y = np.asfarray(y)

    # remove non finite values
    # indexes = np.isfinite(x)
    # x = x[indexes]
    # y = y[indexes]

    # check if sorted
    if np.any(np.diff(x) < 0):
        indexes = np.argsort(x)
        x = x[indexes]
        y = y[indexes]

    size = len(x)

    xdiff = np.diff(x)
    ydiff = np.diff(y)

    # allocate buffer matrices
    Li = np.empty(size)
    Li_1 = np.empty(size-1)
    z = np.empty(size)

    # fill diagonals Li and Li-1 and solve [L][y] = [B]
    Li[0] = (2*xdiff[0])**.5
    Li_1[0] = 0.0
    B0 = 0.0 # natural boundary
    z[0] = B0 / Li[0]

    for i in range(1, size-1, 1):
        Li_1[i] = xdiff[i-1] / Li[i-1]
        Li[i] = (2*(xdiff[i-1]+xdiff[i]) - Li_1[i-1] * Li_1[i-1])**.5
        Bi = 6*(ydiff[i]/xdiff[i] - ydiff[i-1]/xdiff[i-1])
        z[i] = (Bi - Li_1[i-1]*z[i-1])/Li[i]

    i = size - 1
    Li_1[i-1] = xdiff[-1] / Li[i-1]
    Li[i] = (2*xdiff[-1] - Li_1[i-1] * Li_1[i-1])**.5
    Bi = 0.0 # natural boundary
    z[i] = (Bi - Li_1[i-1]*z[i-1])/Li[i]

    # solve [L.T][x] = [y]
    i = size-1
    z[i] = z[i] / Li[i]
    for i in range(size-2, -1, -1):
        z[i] = (z[i] - Li_1[i-1]*z[i+1])/Li[i]

    # find index
    index = x.searchsorted(x0)

    np.clip(index, 1, size-1, index)

    xi1, xi0 = x[index], x[index-1]
    yi1, yi0 = y[index], y[index-1]
    zi1, zi0 = z[index], z[index-1]
    hi1 = xi1 - xi0

    # calculate cubic
    f0 = zi0/(6*hi1)*(xi1-x0)**3 + \
         zi1/(6*hi1)*(x0-xi0)**3 + \
         (yi1/hi1 - zi1*hi1/6)*(x0-xi0) + \
         (yi0/hi1 - zi0*hi1/6)*(xi1-x0)
    return f0

def get_kpap(dn):
    """
    Function: get_kpap(dn)
    ---------------------
    returns the geophysical indices for datetime object dn

    Inputs:
    --------
        dn : datetime object of the requested time

    Outputs:
    --------
        kp, ap, f107, f107a, f107p, daily_kp, daily_ap, dst, ae, ap1

        ap1 == the interpolated ap index to avoid DWM discontinuities

    History:
    --------
        7/21/12 Created, Timothy Duly (duly2@illinois.edu)

    """
    dn_floor = datetime(dn.year, dn.month, dn.day)
    day_index = (dn_floor - generate_kpap.EPOCH).days
    hour_index = int(np.floor(dn.hour/3.))


    kp = GEOPHYSICAL_INDICES[hour_index, day_index]
    ap = GEOPHYSICAL_INDICES[hour_index+8, day_index]

    apall = GEOPHYSICAL_INDICES[:, day_index - 1][8:16]
    apall = np.append(apall,GEOPHYSICAL_INDICES[:, day_index][8:16])
    apall = np.append(apall,GEOPHYSICAL_INDICES[:, day_index + 1][8:16])

    aphours = np.array([1.5, 4.5, 7.5, 10.5, 13.5, 16.5, 19.5, 22.5])
    aphours = np.arange(-22.5,48,3)
    mestime = dn.hour+dn.minute/60.+dn.second/3600.
    aphoursn= np.sort(np.append(np.linspace(np.min(aphours),np.max(aphours),10000),mestime))

    aphours = aphours[~np.isnan(apall)]
    apall = apall[~np.isnan(apall)]


    ap1 = cubic_interp1d([mestime,mestime], aphours, apall)[0]

    f107 = GEOPHYSICAL_INDICES[16, day_index]
    f107a = GEOPHYSICAL_INDICES[17, day_index]
    f107p = GEOPHYSICAL_INDICES[16, day_index-1]

    daily_kp = GEOPHYSICAL_INDICES[18, day_index]
    daily_ap = GEOPHYSICAL_INDICES[19, day_index]

    dst = GEOPHYSICAL_INDICES[20+dn.hour, day_index]
    ae = GEOPHYSICAL_INDICES[44+dn.hour, day_index]

    return kp, ap, f107, f107a, f107p, daily_kp, daily_ap, dst, ae, ap1
