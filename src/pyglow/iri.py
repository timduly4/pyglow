import numpy as np
import os

from iri12py import iri_sub as iri12
from iri16py import iri_sub as iri16
from iri16py import read_ig_rz, readapf107
from .constants import DIR_FILE, nan

IONS = ['O+', 'H+', 'HE+', 'O2+', 'NO+']

# Global variable indicating if IRI 2016 has been initialized with the contents
# of the ionosphere global index (ig_rz.dat) and Ap/F10.7 index (apf107.dat)
# files. IRI 2016 initialization is required only once per session.
__INIT_IRI16 = False


class IRI(object):

    def __init__(self):
        """
        Constructor for IRI representation
        """

        # IRI member variables:
        self.ne = nan
        self.ni = {ion: nan for ion in IONS}

        # Temperatures:
        self.Ti = nan
        self.Te = nan
        self.Tn = nan

        # Ionospheric parameters:
        self.NmF2 = nan
        self.hmF2 = nan

    @staticmethod
    def init_iri16():
        """
        If required (depending on the global variable *__INIT_IRI16*),
        initialize IRI 2016. Return `True` if the model was
        initialized and `False` otherwise.
        """
        if not globals()['__INIT_IRI16']:
            read_ig_rz()
            readapf107()
            globals()['__INIT_IRI16'] = True
            return True
        else:
            return False

    def run(
        self,
        location_time,
        version=2016,
        NmF2=None,
        hmF2=None,
        compute_Ne=True,
        compute_Te_Ti=True,
        compute_Ni=True,
        f107=None,
        f107a=None,
        f1_layer=True,
        bil2000=False,
    ):
        """
        Run IRI model at point time/location and update the object state
        accordingly. If *NmF2* (in [cm^{-3}}]) or *hmF2* (in [km]) are
        specified, input them to the model (see documentation for IRI_SUB)).
        Override the model with *version* --- valid options are currently 2016
        or 2012. Output debugging information if *debug* is true. The toggles
        *compute_Ne*, *compute_Te_Ti*, and *compute_Ni* control, respectively,
        whether electron density, electron and ion temperatures, and ion
        density are computed (restricting the model to only what is required
        can reduce run time) or set to `NaN`.

        :param location_time: Instance of LocationTime
        :param version: Version of IRI to run
        :param NmF2: User-specified NmF2 [cm^-3]
        :param hmF2: User-specified hmF2 [km]
        :param compute_Ne: Switch to compute Ne
        :param compute_Te_Ti: Switch to compute Te and Ti
        :param compute_Ni: Switch to compute Ni
        :param f107: User specified F107
        :param f107a: User specified F107A
        :param f1_layer: If True (default) include F1-layer (JF switches 19,20 = True)
        :param bil2000: If True, use Bil-2000 model for bottomside (JF switch 4). Default IRI is False
        """

        if version == 2016:
            iri_data_stub = 'iri16_data/'
            iri = iri16
            init_iri = IRI.init_iri16
        elif version == 2012:
            iri_data_stub = 'iri12_data/'
            iri = iri12

            def init_iri(): return False
        else:
            raise ValueError(
                "Invalid version of {} for IRI.\n".format(version) +
                "Either 2016 (default) or 2012 is valid."
            )

        jf = np.ones((50,))  # JF switches
        # Standard IRI model flags
        #             | FORTRAN Index
        #             |
        #             V
        jf[3] = 0  # 4 B0,B1 other model-31
        jf[4] = 0  # 5  foF2 - URSI
        jf[5] = 0  # 6  Ni - RBV-10 & TTS-03
        jf[20] = 0  # 21 ion drift not computed
        jf[22] = 0  # 23 Te_topside (TBT-2011)
        jf[27] = 0  # 28 spreadF prob not computed
        jf[28] = 0  # 29 (29,30) => NeQuick
        jf[29] = 0  # 30
        # (Brian found a case that stalled IRI when on):
        jf[32] = 0  # 33 Auroral boundary model off
        jf[34] = 0  # 35 no foE storm update
        # Not standard, but outputs same as values as standard so not an issue
        jf[21] = 0  # 22 ion densities in m^-3 (not %)
        jf[33] = 0  # 34 turn messages off

        if not compute_Ne:
            jf[0] = 0

        if not compute_Te_Ti:
            jf[1] = 0

        if not compute_Ni:
            jf[2] = 0

        oarr = np.zeros((100,))

        if NmF2 is not None:
            # User specified F2 peak density
            jf[7] = 0
            oarr[0] = NmF2 * 100.**3  # IRI expects [m^-3]

        if hmF2 is not None:
            # User specified F2 peak height
            jf[8] = 0
            oarr[1] = hmF2

        if f107 and f107a:
            # User supplied indices

            # Set jf(25) switch to false (in Fortran) which is jf[24] in Python
            jf[24] = 0

            # Set jf(32) switch to false (in Fortran) which is jf[31] in Python
            jf[31] = 0

            # Store user indice for F10.7 in oarr:
            oarr[40] = f107

            # Store user index for F10.7 81 day average in oarr:
            oarr[45] = f107a

            # Reference:
            # https://github.com/timduly4/pyglow/issues/34#issuecomment-340645358
        
        if not f1_layer:
            # Set jf(19) and jf(20) to False (in Fortran index)
            jf[18] = 0
            jf[19] = 0
            
        if bil2000:
            # Set jf(4) to True (in Fortran index)
            jf[3] = True

        # Get current directory:
        my_pwd = os.getcwd()

        # We need to change directories into where the IRI data are located in
        # order to run IRI:
        iri_data_path = os.path.join(
            DIR_FILE,
            iri_data_stub,
        )

        # Override if using local source (typically for testing):
        if 'src' in iri_data_path:
            iri_data_path = "src/pyglow/models/dl_models/iri{}".format(
                version-2000,
            )

        # Change directories to data path:
        os.chdir(iri_data_path)

        # Initalize IRI:
        init_iri()

        # Call f2py wrapper:
        outf = iri(
            jf,
            0,
            location_time.lat,
            location_time.lon,
            location_time.year,
            -location_time.doy,
            location_time.utc_sec/3600.+25.,
            location_time.alt,
            location_time.alt+1,
            1,
            oarr,
        )
        os.chdir(my_pwd)

        # Densities [cm^-3]:
        # Electron density (m^-3):
        self.ne = outf[0, 0] / 100.**3
        if not compute_Ne:
            self.ne = nan
        # O+ Density (%, or m^-3 with JF(22) = 0)
        self.ni['O+'] = outf[4, 0] / 100.**3

        # H+ Density (%, or m^-3 with JF(22) = 0)
        self.ni['H+'] = outf[5, 0] / 100.**3

        # HE+ Density (%, or m^-3 with JF(22) = 0):
        self.ni['HE+'] = outf[6, 0] / 100.**3

        # O2+ Density (%, or m^-3 with JF(22) = 0)
        self.ni['O2+'] = outf[7, 0] / 100.**3

        # NO+ Density (%, or m^-3 with JF(22) = 0)
        self.ni['NO+'] = outf[8, 0] / 100.**3

        if not compute_Ni:
            for ion in IONS:
                self.ni[ion] = nan

        # Temperatures:
        if compute_Te_Ti:
            self.Te = outf[3, 0]  # Electron temperature from IRI (K)
            self.Ti = outf[2, 0]  # Ion temperature from IRI (K)
        else:
            self.Te = float('NaN')
            self.Ti = float('NaN')
        self.Tn = outf[1, 0]  # Neutral temperature from IRI (K)

        # Ionospheric parameters:
        self.NmF2 = oarr[0] / 100.**3  # [items/cm^3]
        self.hmF2 = oarr[1]

        return self
