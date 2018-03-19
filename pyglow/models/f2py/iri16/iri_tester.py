from __future__ import division
from __future__ import print_function

import numpy as np
from datetime import datetime, timedelta

import iri16py
reload(iri16py)
from iri16py import iri_sub as iri

lat = 40.
lon = -80.
dn = datetime(2014, 3, 23, 14, 30)
alt = 250.

doy = dn.timetuple().tm_yday
utc_sec = dn.hour*3600. + dn.minute*60

ni = {}
# ---------------------------------------------------

jf = np.ones((50,))
jf[4]  = 0 # 5  foF2 - URSI (what does that mean?)
jf[5]  = 0 # 6  Ni - RBV-10 & TTS-03
jf[20] = 0 # 21 ion drift not computed
jf[21] = 0 # 22 ion densities in m^-3
jf[22] = 0 # 23 Te_topside (TBT-2011)
jf[28] = 0 # 29 (29,30) => NeQuick
jf[29] = 0 # 30
jf[33] = 0 # 34 messages [on|off]
jf[32] = 0 #  33    Auroral boundary model on/off
           # Brian found a case that stalled IRI

[outf,oarr] = iri(jf,0,\
        lat,\
        lon,\
        int(dn.year),\
        -doy,\
        (utc_sec/3600.+25.),\
        alt,\
        alt+1,\
        1)


Te        = outf[3,0] # electron temperature from IRI (K)
Ti        = outf[2,0] # ion temperature from IRI (K)
Tn_iri    = outf[1,0] # neutral temperature from IRI (K)


ne        = outf[0,0] # electron density (m^-3)
ni['O+']  = outf[4,0] # O+ Density (%, or m^-3 with JF(22) = 0)
ni['H+']  = outf[5,0] # H+ Density (%, or m^-3 with JF(22) = 0)
ni['HE+'] = outf[6,0] # HE+ Density (%, or m^-3 with JF(22) = 0)
ni['O2+'] = outf[7,0] # O2+ Density (%, or m^-3 with JF(22) = 0)
ni['NO+'] = outf[8,0] # NO+ Density (%, or m^-3 with JF(22) = 0)

print("Te =     =", Te)
print("Ti =     =", Ti)
print("Tn_iri = =", Tn_iri)
print("ne =     =", ne)
print("ni['O+'] =", ni['O+'])
print("ni['H+'] =", ni['H+'])
print("ni['HE+']=", ni['HE+'])
print("ni['O2+']=", ni['O2+'])
print("ni['NO+']=", ni['NO+'])
