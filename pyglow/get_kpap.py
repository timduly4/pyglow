from generate_kpap import geophysical_indices, epoch
from datetime import timedelta
from datetime import datetime
import numpy as np

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
        kp, ap, f107, f107a, daily_kp, daily_ap, dst

    History:
    --------
        7/21/12 Created, Timothy Duly (duly2@illinois.edu)

    """
    day_index = ( datetime(dn.year, dn.month, dn.day) - epoch).days
    hour_index = np.floor(dn.hour/3.)
    kp = geophysical_indices[hour_index, day_index]
    ap = geophysical_indices[hour_index+8, day_index]
    f107  = geophysical_indices[16, day_index]
    f107a = geophysical_indices[17, day_index]
    
    daily_kp = geophysical_indices[18, day_index]
    daily_ap = geophysical_indices[19, day_index]
    
    dst = geophysical_indices[20+dn.hour, day_index]

    return kp, ap, f107, f107a, daily_kp, daily_ap, dst



def test_get_kpap():
    # Test it out:
    print "getting first set:" 
    dn = datetime(2008,2,3,15,4)
    kp, ap, f107, f107a, daily_kp, daily_ap, dst = get_kpap(dn)
    print kp, ap, f107, f107a, daily_kp, daily_ap, dst

    print "getting second set:"
    kp, ap, f107, f107a, daily_kp, daily_ap, dst = get_kpap(datetime(1932,1,1))
    print kp, ap, f107, f107a, daily_kp, daily_ap, dst

if __name__ == '__main__':
    test_get_kpap()

      
    
    
