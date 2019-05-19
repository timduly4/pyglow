#!/usr/bin/env python3
""" Testing out the different 'version' keywords """

from datetime import datetime
from pyglow import Point

dn = datetime(2010, 3, 23, 15, 30)
lat = 40.
lon = -80.
alt = 250.

pt = Point(dn, lat, lon, alt)

pt.run_hwm(version=1993)
pt.run_hwm(version=2007)
pt.run_hwm(version=2014)
pt.run_hwm()

pt.run_msis()
pt.run_msis(version=2000)

pt.run_igrf()
pt.run_igrf(version=11)
pt.run_igrf(version=12)

pt.run_iri()
pt.run_iri(version=2016)
pt.run_iri(version=2012)

try:
    pt.run_iri(version=2020)
except ValueError as e:
    print("Caught an exception: `{}`".format(e))
