from __future__ import absolute_import
import numpy as np
from datetime import datetime
from datetime import timedelta
from scipy.stats import nanmean
from .get_kpap import get_kpap

def get_apmsis(dn):
    """
    Function: get_apmsis(dn)
    ---------------------
    returns an array of calculated ap indices suitable for MSIS.
    MSIS requires an array of ap values, described in nrlmsise00.f.
    This Python function formulates the various ap values for MSIS. From the
    fortran subroutine, we see that

        AP - MAGNETIC INDEX(DAILY) OR WHEN SW(9)=-1. :
           - ARRAY CONTAINING:
             (1) DAILY AP
             (2) 3 HR AP INDEX FOR CURRENT TIME
             (3) 3 HR AP INDEX FOR 3 HRS BEFORE CURRENT TIME
             (4) 3 HR AP INDEX FOR 6 HRS BEFORE CURRENT TIME
             (5) 3 HR AP INDEX FOR 9 HRS BEFORE CURRENT TIME
             (6) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 12 TO 33 HRS
                 PRIOR   TO CURRENT TIME
             (7) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 36 TO 57 HRS
                 PRIOR  TO CURRENT TIME

    Inputs:
    --------
        dn : datetime object of the requested time

    Outputs:
    --------
        out : a 1x7 array of the caclulated ap indices

    History:
    --------
        7/21/12 Created, Timothy Duly (duly2@illinois.edu)

    """
    out = float('nan')*np.zeros(7)

    # (1) DAILY AP
    tmp1, ap, tmp2, tmp3, tmp4, daily_ap, tmp5 = get_kpap(dn)
    out[0] = daily_ap

    # (2) 3 HR AP INDEX FOR CURRENT TIME
    out[1] = ap

    # (3) 3 HR AP INDEX FOR 3 HRS BEFORE CURRENT TIME
    tmp1, ap, tmp2, tmp3, tmp4, tmp5, tmp6 = get_kpap(dn+timedelta(hours=-3))
    out[2] = ap

    # (4) 3 HR AP INDEX FOR 6 HRS BEFORE CURRENT TIME
    tmp1, ap, tmp2, tmp3, tmp4, tmp5, tmp6 = get_kpap(dn+timedelta(hours=-6))
    out[3] = ap

    # (5) 3 HR AP INDEX FOR 9 HRS BEFORE CURRENT TIME
    tmp1, ap, tmp2, tmp3, tmp4, tmp5, tmp6 = get_kpap(dn+timedelta(hours=-9))
    out[4] = ap

    # (6) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 12 TO 33 HRS
    #     PRIOR   TO CURRENT TIME

    temp = np.zeros(8)

    tmp1, temp[0], tmp2, tmp3, tmp4, tmp5, tmp6 = get_kpap(dn+timedelta(hours=-12))
    tmp1, temp[1], tmp2, tmp3, tmp4, tmp5, tmp6 = get_kpap(dn+timedelta(hours=-15))
    tmp1, temp[2], tmp2, tmp3, tmp4, tmp5, tmp6 = get_kpap(dn+timedelta(hours=-18))
    tmp1, temp[3], tmp2, tmp3, tmp4, tmp5, tmp6 = get_kpap(dn+timedelta(hours=-21))
    tmp1, temp[4], tmp2, tmp3, tmp4, tmp5, tmp6 = get_kpap(dn+timedelta(hours=-24))
    tmp1, temp[5], tmp2, tmp3, tmp4, tmp5, tmp6 = get_kpap(dn+timedelta(hours=-27))
    tmp1, temp[6], tmp2, tmp3, tmp4, tmp5, tmp6 = get_kpap(dn+timedelta(hours=-30))
    tmp1, temp[7], tmp2, tmp3, tmp4, tmp5, tmp6 = get_kpap(dn+timedelta(hours=-33))

    out[5] = nanmean(temp)

    # (7) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 36 TO 57 HRS
    #     PRIOR  TO CURRENT TIME

    temp = np.zeros(8)

    tmp1, temp[0], tmp2, tmp3, tmp4, tmp5, tmp6 = get_kpap(dn+timedelta(hours=-36))
    tmp1, temp[1], tmp2, tmp3, tmp4, tmp5, tmp6 = get_kpap(dn+timedelta(hours=-39))
    tmp1, temp[2], tmp2, tmp3, tmp4, tmp5, tmp6 = get_kpap(dn+timedelta(hours=-42))
    tmp1, temp[3], tmp2, tmp3, tmp4, tmp5, tmp6 = get_kpap(dn+timedelta(hours=-45))
    tmp1, temp[4], tmp2, tmp3, tmp4, tmp5, tmp6 = get_kpap(dn+timedelta(hours=-48))
    tmp1, temp[5], tmp2, tmp3, tmp4, tmp5, tmp6 = get_kpap(dn+timedelta(hours=-51))
    tmp1, temp[6], tmp2, tmp3, tmp4, tmp5, tmp6 = get_kpap(dn+timedelta(hours=-54))
    tmp1, temp[7], tmp2, tmp3, tmp4, tmp5, tmp6 = get_kpap(dn+timedelta(hours=-57))

    out[6] = nanmean(temp)

    return out




def test_get_apmsis():
    out = get_apmsis(datetime(2000,3,23,0))
    print("ap indices for msis are:\n",out)

if __name__ == '__main__':
    test_get_apmsis()


