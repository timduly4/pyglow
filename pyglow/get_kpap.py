from datetime import datetime
import numpy as np

import generate_kpap

# Fetch geophysical indices:
# (Global variable to make it fast)
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
        kp, ap, f107, f107a, daily_kp, daily_ap, dst, ae

    History:
    --------
        7/21/12 Created, Timothy Duly (duly2@illinois.edu)

    """
    dn_floor = datetime(dn.year, dn.month, dn.day)
    day_index = (dn_floor - generate_kpap.EPOCH).days
    hour_index = int(np.floor(dn.hour/3.))


    kp = GEOPHYSICAL_INDICES[hour_index, day_index]
    ap = GEOPHYSICAL_INDICES[hour_index+8, day_index]
    f107  = GEOPHYSICAL_INDICES[16, day_index]
    f107a = GEOPHYSICAL_INDICES[17, day_index]
    
    daily_kp = GEOPHYSICAL_INDICES[18, day_index]
    daily_ap = GEOPHYSICAL_INDICES[19, day_index]
    
    dst = GEOPHYSICAL_INDICES[20+dn.hour, day_index]
    ae = GEOPHYSICAL_INDICES[44+dn.hour, day_index]

    return kp, ap, f107, f107a, daily_kp, daily_ap, dst, ae

def test_get_kpap():
    # Test it out:
    print("getting first set:" )
    dn = datetime(2008,2,3,15,4)
    kp, ap, f107, f107a, daily_kp, daily_ap, dst, ae = get_kpap(dn)
    print(kp, ap, f107, f107a, daily_kp, daily_ap, dst, ae)

    print("getting second set:")
    kp, ap, f107, f107a, daily_kp, daily_ap, dst, ae = get_kpap(datetime(1932,1,1))
    print(kp, ap, f107, f107a, daily_kp, daily_ap, dst, ae)


if __name__ == '__main__':
    test_get_kpap()
    
    
