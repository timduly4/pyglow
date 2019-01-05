import os
import numpy as np

from .constants import DIR_FILE, nan
from hwm93py import gws5 as hwm93
from hwm07py import hwmqt as hwm07
from hwm14py import hwm14


class HWM(object):

    def __init__(self):
        """ Constructor for HWM representation """

        self.u = nan
        self.v = nan
        self.hwm_version = None

        # Data path:
        self.data_path_stub = DIR_FILE
        self.testing_data_stub = False
        # Override if using local source (typically for testing):
        if 'src' in self.data_path_stub:
            self.data_path_stub = "src/pyglow/models/dl_models"
            self.testing_data_stub = True

    def run(self, location_time, version,
            f107=None, f107a=None, ap=None, ap_daily=None):
        """
        Wrapper to call various HWM models

        :param location_time: Instance of LocationTime
        :param version: Version of HWM to run
        :param f107: f107 indice (used for 93, 07)
        :param f107a: f107a indice (used for 93, 07)
        :param ap: ap indice (used for 07, 14)
        :param ap_daily: ap_daily indice (used for 93)
        """

        # HWM93:
        if version == 1993:
            if not f107 or not f107a or not ap_daily:
                raise ValueError(
                    "Must supply f107, f107a, and ap_daily for HWM93"
                )
            self._run_hwm93(location_time, f107, f107a, ap_daily)

        # HWM07:
        elif version == 2007:
            if not f107 or not f107a or not ap:
                raise ValueError(
                    "Must supply f107, f107a, and ap for HWM07"
                )
            self._run_hwm07(location_time, f107, f107a, ap)

        # HWM14:
        elif version == 2014:
            if not ap:
                raise ValueError(
                    "Must supply ap for HWM14"
                )
            self._run_hwm14(location_time, ap)

        # Unknown version:
        else:
            raise ValueError(
                "Invalid version of {} for HWM.\n".format(version) +
                "Either 2014, 2007, or 1993 is valid."
            )

        return self

    def _run_hwm93(self, location_time, f107, f107a, ap_daily):
        """
        HWM 1993 Climatological model.

        :param location_time: Instance of LocationTime
        :param f107: f107 indice
        :param f107a: f107a indice
        :param ap_daily: ap_daily indice
        """

        # Call HWM93 wrapper:
        w = hwm93(
            location_time.iyd,
            location_time.utc_sec,
            location_time.alt,
            location_time.lat,
            np.mod(location_time.lon, 360),
            location_time.slt_hour,
            f107a,
            f107,
            ap_daily,
        )
        self.v = w[0]
        self.u = w[1]
        self.hwm_version = '93'

        return self

    def _run_hwm07(self, location_time, f107, f107a, ap):
        """
        HWM 2007 Climatological model.

        :param location_time: Instance of LocationTime
        :param f107: f107 indice
        :param f107a: f107a indice
        :param ap: ap indice
        """

        # Grab current directory:
        my_pwd = os.getcwd()

        # Figure out HWM07 data folder:
        if self.testing_data_stub:
            folder = "hwm07"
        else:
            folder = "hwm07_data"
        hwm07_data_path = os.path.join(
            self.data_path_stub,
            folder,
        )

        # Change directory to HWM07 data path:
        os.chdir(hwm07_data_path)

        # Call HWM07 wrapper:
        w = hwm07(
            location_time.iyd,
            location_time.utc_sec,
            location_time.alt,
            location_time.lat,
            np.mod(location_time.lon, 360),
            location_time.slt_hour,
            f107a,
            f107,
            [nan, ap],
        )

        # Change back to original directory:
        os.chdir(my_pwd)

        # Assign outputs:
        self.v = w[0]
        self.u = w[1]
        self.hwm_version = '07'

        return self

    def _run_hwm14(self, location_time, ap):
        """
        HWM 2014 Climatological model.

        :param location_time: Instance of LocationTime
        :param ap: ap indice
        """

        # Grab current directory:
        my_pwd = os.getcwd()

        # Figure out HWM14 data folder:
        if self.testing_data_stub:
            folder = "hwm14"
        else:
            folder = "hwm14_data"
        hwm14_data_path = os.path.join(
            self.data_path_stub,
            folder,
        )

        # Change directory to HWM14 data path:
        os.chdir(hwm14_data_path)

        # Call HWM14 wrapper:
        v, u = hwm14(
            location_time.iyd,
            location_time.utc_sec,
            location_time.alt,
            location_time.lat,
            np.mod(location_time.lon, 360),
            nan,
            nan,
            nan,
            [nan, ap],
        )

        # Change back to original directory:
        os.chdir(my_pwd)

        # Assign outputs:
        self.v = v
        self.u = u
        self.hwm_version = '14'

        return self
