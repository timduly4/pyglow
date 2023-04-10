from datetime import datetime as dt
import matplotlib.pyplot as plt
from src import pyglow
import numpy as np
import os

'''
This script shows the HWM/DWM test.
'''

AP0 = []
AP1 = []
TIME = []
lat = 0
lon = -145

indd = 0

for lat in np.arange(45,75.1, 5):

    fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(10, 5), sharex=True)

    fig.subplots_adjust(hspace=0, wspace=0)
    fig.suptitle('Ap index effect and DWM test at Latitude = %2.1f'%lat, fontsize=16)
    for i in range(2):
        ax[-1, i].set_xlabel('UT (hours of 2023/02/17)')
        for j in range(2):

            if i == 0:
                if j < 1:
                    ax[j, 0].set_ylabel('Winds\n(m/s)', weight='bold')
                else:
                    ax[j, 0].set_ylabel('Ap index', weight='bold')
    AP = []
    AP1 = []
    FFF107 = []
    FFF107a = []

    apn = 0
    for dwm in [None, 'off','smooth']:
        UUU, VVV = [], []
        TIME = []
        for hour in np.arange(0, 24):
            for minute in np.arange(0, 60, 30):
                dn = dt(2023, 2, 17, hour, minute, 30)
                TIME.append(hour + minute / 60. + 30 / 3600.)

                pt = pyglow.Point(dn, lat, lon, 250)
                pt.run_hwm(version=2014, dwm=dwm)

                Vhwm = pt.v
                Uhwm = pt.u

                UUU.append(Uhwm)
                VVV.append(Vhwm)
                if apn == 0:
                    AP.append(pt.ap)
                    AP1.append(pt.ap1)
                    FFF107.append(pt.f107)
                    FFF107a.append(pt.f107a)
        if dwm == 'off':
            linestyle = '--'
            linewidth = 0.5
        elif (dwm == 'on') or (dwm == None):
            linestyle = '-'
            linewidth = 0.5
        else:
            linestyle = '-'
            linewidth = 1.0

        if dwm == None:
            label = 'None'
        else:
            label = dwm.capitalize()
        ax[0, 0].plot(TIME, UUU, color='blue', linestyle=linestyle, linewidth=linewidth, label=label)
        ax[0, 1].plot(TIME, VVV, color='red', linestyle=linestyle, linewidth=linewidth, label=label)

        apn = apn + 1
    ax[1, 0].plot(TIME, AP, color='black', linestyle='-', label='Ap')
    ax[1, 1].plot(TIME, AP1, color='black', linestyle='-', label='Ap modified')

    ax[1, 1].text(1.035, 0.5, 'Ap', rotation=90, horizontalalignment='center',
                  verticalalignment='center', transform=ax[1, 1].transAxes)

    ax[0, 0].set_title('Zonal')
    ax[0, 1].set_title('Meridional')

    for i in range(2):
        for j in range(2):
            ax[j, i].grid()
            ax[j, i].legend(loc='best')
            if j < 1:
                ax[j, i].set_ylim(-130, 130)
            else:
                ax[j, i].set_ylim(0, 21)
            ax[j, i].locator_params(axis='y', nbins=7)
            plt.setp(ax[j, 1].get_yticklabels()[:], visible=False)

    plt.savefig(os.getcwd() + '/DWM_Test_%03d.png' % indd,
                bbox_inches='tight', dpi=300)
    indd = indd + 1
    print('DWM Test at %2.2f' % lat)
    plt.close()
