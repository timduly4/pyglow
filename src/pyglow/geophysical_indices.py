import numpy as np

from .get_kpap import get_kpap
from .get_apmsis import get_apmsis
from .constants import nan


class Indice(object):

    def __init__(self, dn):
        """
        :param dn: datetime of indice data
        """

        # Store datetime associated with indice:
        self.dn = dn
        # Assign to member variables:
        self.kp = nan
        self.ap = nan
        self.f107 = nan
        self.f107a = nan
        self.f107p = nan
        self.kp_daily = nan
        self.ap_daily = nan
        self.dst = nan
        self.ae = nan
        self.ap1 = nan

        # AP values for MSIS:
        self.apmsis = [nan, ] * 7

    def run(self):
        """
        Calculates the geophysical indices
        """

        # Geophysical indices:
        kp, ap, f107, f107a, f107p, kp_daily, ap_daily, dst, ae, ap1 = \
            get_kpap(self.dn)

        # Assign to member variables:
        self.kp = kp
        self.ap = ap
        self.f107 = f107
        self.f107a = f107a
        self.f107p = f107p
        self.kp_daily = kp_daily
        self.ap_daily = ap_daily
        self.dst = dst
        self.ae = ae
        self.ap1 = ap1

        # AP values for MSIS:
        self.apmsis = get_apmsis(self.dn)

    def all_nan(self):
        """
        Returns a boolean indicating if all indices are NaN
        """

        all_nan = True
        all_nan &= np.isnan(self.kp)
        all_nan &= np.isnan(self.ap)
        all_nan &= np.isnan(self.f107)
        all_nan &= np.isnan(self.f107a)
        all_nan &= np.isnan(self.f107p)
        all_nan &= np.isnan(self.kp_daily)
        all_nan &= np.isnan(self.ap_daily)
        all_nan &= np.isnan(self.dst)
        all_nan &= np.isnan(self.ae)
        all_nan &= all(np.isnan(_) for _ in self.apmsis)

        return all_nan
