
#!/usr/bin/env python3
""" Profile plot to compare IRI 2012 and IRI 2016. """
import matplotlib.pyplot as plt
from datetime import datetime
import numpy as np

import pyglow

# Inputs:
lat = 40.
lon = -80.
alt = 250.
alts = np.linspace(100., 500., 101)
dn = datetime(2015, 3, 23, 15, 30)

ne_2012 = []
ne_2016 = []

# Calculate for both IRI model year 2012 and 2016:
for alt in alts:
    print("Computing alt=%3.1f km..." % (alt))
    pt = pyglow.Point(dn, lat, lon, alt)

    pt.run_iri()  # default year is 2016
    ne_2016.append(pt.ne)

    pt.run_iri(version=2012)  # Can revert back to 2012 model, if necessary.
    ne_2012.append(pt.ne)

# Plot
plt.figure(1)
plt.clf()
plt.semilogx(ne_2016, alts, 'bo-', label='IRI Model Year: 2016')
plt.semilogx(ne_2012, alts, 'r.--', label='IRI Model Year: 2012')
plt.grid()
plt.xlabel(r'$n_e$ [cm$^{-3}$]')
plt.ylabel('Altitude [km]')
plt.title(r'%s UT, lat=%3.1f$^\circ$, lon=%3.1f$^\circ$' %
          (dn.strftime('%Y-%m-%d %H:%M:%S'), lat, lon))
plt.legend(loc=0)
plt.draw()
plt.show()
