#! /usr/bin/env python3

from datetime import datetime

import pyglow
from pyglow import Point

print(
    "pyglow version = {}".format(
        pyglow.__version__,
    )
)

# Make sure we're using a pyglow version
# that implements the user defined indices:
assert(float(pyglow.__version__) > 1.0)

lat = 30.
lon = -10.
alt = 250.
dn = datetime(2000, 1, 1)

print("Testing user defined indices in IRI:")
print("------- ---- ------- ------- -- ----")
# Using default F10.7:
pt = Point(dn, lat, lon, alt)
pt.run_iri()
print("Using default F10.7:")
print("f10.7 = {:3.2f}".format(pt.f107))
print("f10.7a = {:3.2f}".format(pt.f107a))
print("ne = {:3.2f}".format(pt.ne))
print("")

# Using user-defined F10.7:
pt = Point(dn, lat, lon, alt, user_ind=True)
pt.f107 = 200  # We need to define our F10.7, now
pt.f107a = 200  # We need to define our F10.7a, now
pt.run_iri()
print("Using user defined F10.7:")
print("f10.7 = {:3.2f}".format(pt.f107))
print("f10.7a = {:3.2f}".format(pt.f107a))
print("ne = {:3.2f}".format(pt.ne))
print("")

# Using user-defined F10.7:
pt = Point(dn, lat, lon, alt, user_ind=True)
pt.f107 = 125.6  # We need to define our F10.7, now
pt.f107a = 160.55  # We need to define our F10.7a, now
pt.run_iri()
print("Using user defined F10.7:")
print("f10.7 = {:3.2f}".format(pt.f107))
print("f10.7a = {:3.2f}".format(pt.f107a))
print("ne = {:3.2f}".format(pt.ne))
print("")
