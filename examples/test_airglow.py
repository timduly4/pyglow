#!/usr/bin/env python3

'''
Profile 7774 and 6300-nm emissions.
2016-09-17
'''
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import numpy as np

import pyglow

# Inputs:
lat = 40.
lon = -80.
alt = 250.
alts = np.linspace(100., 500., 101)
dn = datetime(2015, 3, 23, 15, 30)

# Calculate and save airglow values:
ag6300 = []
ag7774 = []
for alt in alts:
    print("Computing alt={:3.1f} km...".format(alt))
    pt = pyglow.Point(dn, lat, lon, alt)
    pt.run_airglow()
    ag6300.append(pt.ag6300)
    ag7774.append(pt.ag7774)

# Plot:
plt.figure(1, figsize=(8,8));
plt.clf()
plt.plot(
    ag6300,
    alts,
    'ro-',
    label='630.0-nm airglow VER',
)
plt.plot(
    ag7774,
    alts,
    'go--',
    label='777.4-nm airglow VER',
)
plt.grid()
plt.xlabel('Volume Emission Rate (VER) [ph/cm$^3$/s]')
plt.ylabel('Altitude [km]')
plt.title(
    'Testing: Point.run_airglow() \n' \
            '{} UT, Lat = {:3.1f}$^\circ$, Lon = {:3.1f}$^\circ$'.\
            format(dn.strftime('%Y-%m-%d %H:%M:%S'), lat, lon)
)
plt.legend(loc=1, fontsize=10)
plt.draw()
plt.show()
