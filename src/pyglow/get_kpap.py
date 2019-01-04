from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from datetime import datetime
import numpy as np

from . import generate_kpap

# Fetch geophysical indices (Global variable to make it fast):
GEOPHYSICAL_INDICES = generate_kpap.fetch()


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
        kp, ap, f107, f107a, f107p, daily_kp, daily_ap, dst, ae

    History:
    --------
        7/21/12 Created, Timothy Duly (duly2@illinois.edu)

    """
    dn_floor = datetime(dn.year, dn.month, dn.day)
    day_index = (dn_floor - generate_kpap.EPOCH).days
    hour_index = int(np.floor(dn.hour/3.))

    kp = GEOPHYSICAL_INDICES[hour_index, day_index]
    ap = GEOPHYSICAL_INDICES[hour_index+8, day_index]
    f107 = GEOPHYSICAL_INDICES[16, day_index]
    f107a = GEOPHYSICAL_INDICES[17, day_index]
    f107p = GEOPHYSICAL_INDICES[16, day_index-1]

    daily_kp = GEOPHYSICAL_INDICES[18, day_index]
    daily_ap = GEOPHYSICAL_INDICES[19, day_index]

    dst = GEOPHYSICAL_INDICES[20+dn.hour, day_index]
    ae = GEOPHYSICAL_INDICES[44+dn.hour, day_index]

    return kp, ap, f107, f107a, f107p, daily_kp, daily_ap, dst, ae
