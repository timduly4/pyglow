#!/usr/bin/env python3
""" Example Python script using pyglow and its climatological models to plot
profile of airglow emission and electron density """

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime, timedelta

import pyglow

matplotlib.rcParams.update({'font.size': 16})

# Setting lat, lon, and a range of altitudes
lat = 18.37  # Arecibo
lon = -66.62
alts = np.linspace(85, 1000, 500)

# Set time (i.e., dn)
dn_lt = datetime(2012, 3, 22, 0, 0)
tz = np.ceil(lon/15.)
dn_ut = dn_lt - timedelta(hours=tz)

# Airglow calculation using pyglow:
ag, ne = [], []
for alt in alts:
    pt = pyglow.Point(dn_ut, lat, lon, alt)

    pt.run_iri()
    pt.run_msis()
    pt.run_airglow()

    ag.append(pt.ag6300)
    ne.append(pt.ne)

# Plot data:
plt.figure(1, figsize=(7, 8))
plt.clf()
plt.semilogy(ag, alts, '-k', lw=4, label='$V_{630.0}$')
plt.grid()
plt.ylim([alts[0], alts[-1]])
the_yticks = np.arange(100, 1001, 100)
plt.yticks(the_yticks, the_yticks)
plt.xlabel('Volume Emission Rate [ph/cm$^3$/$s$]')
plt.ylabel('Altitude [km]')
plt.legend(loc='lower right')

ax2 = plt.gca().twiny()
ax2.semilogx(ne, alts, '--b', lw=4, label='$n_e$')
plt.xlabel('Electron density [items/cm$^3$]', color='blue')
ax2.legend(loc='upper right')

plt.draw()
plt.subplots_adjust(left=0.15)

# Save png:
png = '/tmp/test_ne_airglow_profile.png'
print("Saving: {}".format(png))
plt.savefig(png, dpi=200)

plt.show()
