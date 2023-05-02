from __future__ import print_function
from __future__ import absolute_import
import numpy as np
from datetime import timedelta

from .get_kpap import get_kpap


def get_apmsis(dn):
    """
    Function: get_apmsis(dn)
    ---------------------
    returns an array of calculated ap indices suitable for MSIS.
    MSIS requires an array of ap values, described in nrlmsise00.f.
    This Python function formulates the various ap values for MSIS. From the
    fortran subroutine, we see that:

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
        out : a 1x7 array of the calculated ap indices

    History:
    --------
        7/21/12 Created, Timothy Duly (duly2@illinois.edu)

    """
    out = float('nan')*np.zeros(8)

    # (1) DAILY AP
    _, ap, _, _, _, _, daily_ap, _, _, ap1 = get_kpap(dn)
    out[0] = daily_ap

    # (2) 3 HR AP INDEX FOR CURRENT TIME
    out[1] = ap

    # (3) 3 HR AP INDEX FOR 3 HRS BEFORE CURRENT TIME
    _, ap, _, _, _, _, _, _, _, _ = get_kpap(dn+timedelta(hours=-3))
    out[2] = ap

    # (4) 3 HR AP INDEX FOR 6 HRS BEFORE CURRENT TIME
    _, ap, _, _, _, _, _, _, _, _ = get_kpap(dn+timedelta(hours=-6))
    out[3] = ap

    # (5) 3 HR AP INDEX FOR 9 HRS BEFORE CURRENT TIME
    _, ap, _, _, _, _, _, _, _, _ = get_kpap(dn+timedelta(hours=-9))
    out[4] = ap

    # (6) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 12 TO 33 HRS
    #     PRIOR   TO CURRENT TIME

    temp = np.zeros(8)

    _, temp[0], _, _, _, _, _, _, _, _ = get_kpap(dn+timedelta(hours=-12))
    _, temp[1], _, _, _, _, _, _, _, _ = get_kpap(dn+timedelta(hours=-15))
    _, temp[2], _, _, _, _, _, _, _, _ = get_kpap(dn+timedelta(hours=-18))
    _, temp[3], _, _, _, _, _, _, _, _ = get_kpap(dn+timedelta(hours=-21))
    _, temp[4], _, _, _, _, _, _, _, _ = get_kpap(dn+timedelta(hours=-24))
    _, temp[5], _, _, _, _, _, _, _, _ = get_kpap(dn+timedelta(hours=-27))
    _, temp[6], _, _, _, _, _, _, _, _ = get_kpap(dn+timedelta(hours=-30))
    _, temp[7], _, _, _, _, _, _, _, _ = get_kpap(dn+timedelta(hours=-33))

    out[5] = np.nan if all(np.isnan(temp)) else np.nanmean(temp)

    # (7) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 36 TO 57 HRS
    #     PRIOR  TO CURRENT TIME

    temp = np.zeros(8)

    _, temp[0], _, _, _, _, _, _, _, _ = get_kpap(dn+timedelta(hours=-36))
    _, temp[1], _, _, _, _, _, _, _, _ = get_kpap(dn+timedelta(hours=-39))
    _, temp[2], _, _, _, _, _, _, _, _ = get_kpap(dn+timedelta(hours=-42))
    _, temp[3], _, _, _, _, _, _, _, _ = get_kpap(dn+timedelta(hours=-45))
    _, temp[4], _, _, _, _, _, _, _, _ = get_kpap(dn+timedelta(hours=-48))
    _, temp[5], _, _, _, _, _, _, _, _ = get_kpap(dn+timedelta(hours=-51))
    _, temp[6], _, _, _, _, _, _, _, _ = get_kpap(dn+timedelta(hours=-54))
    _, temp[7], _, _, _, _, _, _, _, _ = get_kpap(dn+timedelta(hours=-57))

    out[6] = np.nan if all(np.isnan(temp)) else np.nanmean(temp)

    # (8) CUBIC INTERPOLATED AP INDEX AT THE TIME dn

    out[7] = ap1
    return out
