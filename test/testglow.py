from pyglow.pyglow import Point
from datetime import datetime

dn = datetime(2011, 3, 23, 9, 30)
lat = 0.
lon = -80.
alt = 250.

pt = Point(dn, lat, lon, alt)

print "Before running any models:"
print pt

pt.run_igrf()
pt.run_hwm93()
pt.run_msis()
pt.run_iri()

print "After running models:"
print pt