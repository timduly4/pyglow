from __future__ import division
from __future__ import print_function
from future import standard_library
standard_library.install_aliases()
from builtins import str
from builtins import range
from datetime import datetime, timedelta, date
from pickle import load, dump
import numpy as np
import os
from numpy import nanmean
import sys
import glob


"""
Description:
------------
This library parses the geophysical indices from the yearly files located in
./kpap/ and ./dst/ .  The files are obtained by running pyglow.update_indices()

After parsing the files, a variable is created, "geophysical_indices", whose
dimensions are 20 x (number of days), with:

    geophysical_indices[0,:] = kp value from 00:00 - 03:00 UTC
    geophysical_indices[1,:] = kp value from 03:00 - 06:00
    geophysical_indices[2,:] = kp value from 06:00 - 09:00
    geophysical_indices[3,:] = kp value from 09:00 - 12:00
    geophysical_indices[4,:] = kp value from 12:00 - 15:00
    geophysical_indices[5,:] = kp value from 15:00 - 18:00
    geophysical_indices[6,:] = kp value from 18:00 - 21:00
    geophysical_indices[7,:] = kp value from 21:00 - 24:00

    geophysical_indices[8,:]  = ap value from 00:00 - 03:00
    geophysical_indices[9,:]  = ap value from 03:00 - 06:00
    geophysical_indices[10,:] = ap value from 06:00 - 09:00
    geophysical_indices[11,:] = ap value from 09:00 - 12:00
    geophysical_indices[12,:] = ap value from 12:00 - 15:00
    geophysical_indices[13,:] = ap value from 15:00 - 18:00
    geophysical_indices[14,:] = ap value from 18:00 - 21:00
    geophysical_indices[15,:] = ap value from 21:00 - 24:00

    geophysical_indices[16,:] = f107 value
    geophysical_indices[17,:] = f107a value

    geophysical_indices[18,:] = daily_kp value
    geophysical_indices[19,:] = daily_ap value

    geophysical_indices[20,:] = dst value from 00:00 - 01:00 UTC
    geophysical_indices[21,:] = dst value from 01:00 - 02:00
    ...
    geophysical_indices[43,:] = dst value from 23:00 - 24:00

    geophysical_indices[44,:] = ae value from 00:00 - 01:00 UTC
    geophysical_indices[45,:] = ae value from 01:00 - 02:00
    ...
    geophysical_indices[67,:] = ae value from 23:00 - 24:00

This script is meant to be run in conjunction with "get_kpap.py" to grab the
geophysical indices from a specific datetime object.

-----------
Side Notes:
-----------
An easy way to grab the yearly data, in Python of course:
(make sure the folder "kpap" is created in the current directory):

import os
for y in range(1932, END_YEAR):
    os.system(
        'wget
        ftp://ftp.ngdc.noaa.gov/STP/GEOMAGNETIC_DATA/INDICES/KP_AP/%4i ./kpap/'
        % y
    )

--------
History:
--------
    7/21/12 : Created, Timothy Duly (duly2@illinois.edu)
    1/07/15 : Changed END_YEAR so that it doesn't crash during
              January. Brian Harding (bhardin2@illinois.edu)
    1/08/15 : Added DST index. Brian Harding (bhardin2@illinois.edu)
    4/13/16 : Added AE index. Daniel Fisher (dfisher2@illinois.edu)
   10/03/16 : Use file cache to optimize initial spin-up.
              Mark D. Butala (butala@illinois.edu)
"""

# TODO
# [ ] Is dividing by 10 correct for kp?
# [ ] daily kp -- should it be the average of sum? right now it's the sum...
# [ ] fix path issue... how to make it callable anywhere?

# Load all data up to and including this year, if it exists
END_YEAR = date.today().year + 1
EPOCH = datetime(1932, 1, 1)
PYGLOW_PATH = os.path.dirname(__file__)
if not PYGLOW_PATH:  # Corner case: we are in the installation directory
    PYGLOW_PATH = './'

# File to store table of index file modification times:
MTIME_TABLE_FNAME = os.path.join(PYGLOW_PATH, 'mtime_table.pkl')

# File to store cached geophysical_indices array:
GEOPHYSICAL_INDICES_FNAME = os.path.join(
    PYGLOW_PATH,
    'geophysical_indices.npy',
)


def get_mtime_table():
    """

    """
    mtime_table = {}

    # kpap files:
    for y in range(1932, END_YEAR):
        fn = PYGLOW_PATH + "/kpap/%4i" % y
        if os.path.isfile(fn):
            mtime_table[fn] = os.path.getmtime(fn)

    # dst files:
    oldies = ['1957_1969', '1970_1989', '1990_2004']
    dst_path = '%s/dst/' % PYGLOW_PATH
    files = glob.glob('%s??????' % dst_path)  # files like 201407
    old_files = [
        '%s%s' % (dst_path, old) for old in oldies
    ]  # older files listed above
    files.extend(old_files)  # a list of every dst file
    for fn in files:
        mtime_table[fn] = os.path.getmtime(fn)

    # ae files:
    ae_path = '%s/ae/' % PYGLOW_PATH
    files = glob.glob('%s*' % ae_path)  # find files like 1975
    for fn in files:
        mtime_table[fn] = os.path.getmtime(fn)

    return mtime_table


def update_required():
    """
    Determines if new data are available to index into geophysical_indices.npy
    """
    if not os.path.isfile(MTIME_TABLE_FNAME):
        return True
    else:
        with open(MTIME_TABLE_FNAME, 'rb') as fid:
            mtime_table = load(fid)
        new_mtime_table = get_mtime_table()
        for key, val in new_mtime_table.items():
            try:
                if mtime_table[key] != new_mtime_table[key]:
                    return True
            except (IndexError, KeyError):
                # Differing keys in mtime_tables,
                # so we need to update:
                return True
        # All current file modification times matched those associated
        # with cached geophysical_indices
        return False


def generate_kpap():
    """

    """

    print('[generate_kpap.py] Generating geophysical indices table...')

    # Files have changed --- update geophysical indices"""
    mtime_table = get_mtime_table()

    """ Part 1: Parsing the raw data files """
    # Create empty dictionaries:
    ap = {}
    kp = {}
    f107 = {}
    daily_kp = {}
    daily_ap = {}
    f107a = {}
    dst = {}
    ae = {}

    for y in range(1932, END_YEAR):
        fn = PYGLOW_PATH + "/kpap/%4i" % y
        if os.path.isfile(fn):  # If the file has been downloaded
            f = open(fn)
            for x in f.readlines():
                # parse the line for the year, only last 2 digits:
                year = int(x[0:2])

                if year < 30:  # need to change this is 2030....
                    year = year + 2000
                else:
                    year = year + 1900

                month = int(x[2:4])  # parse the line for month
                day = int(x[4:6])  # ... and the days

                # Parse the values for kp, ap, daily_kp, and daily_ap:
                kp[datetime(year, month, day, 0)] = int(x[12:14])/10.
                kp[datetime(year, month, day, 3)] = int(x[14:16])/10.
                kp[datetime(year, month, day, 6)] = int(x[16:18])/10.
                kp[datetime(year, month, day, 9)] = int(x[18:20])/10.
                kp[datetime(year, month, day, 12)] = int(x[20:22])/10.
                kp[datetime(year, month, day, 15)] = int(x[22:24])/10.
                kp[datetime(year, month, day, 18)] = int(x[24:26])/10.
                kp[datetime(year, month, day, 21)] = int(x[26:28])/10.

                ap[datetime(year, month, day, 0)] = int(x[31:34])
                ap[datetime(year, month, day, 3)] = int(x[34:37])
                ap[datetime(year, month, day, 6)] = int(x[37:40])
                ap[datetime(year, month, day, 9)] = int(x[40:43])
                ap[datetime(year, month, day, 12)] = int(x[43:46])
                ap[datetime(year, month, day, 15)] = int(x[46:49])
                ap[datetime(year, month, day, 18)] = int(x[49:52])
                ap[datetime(year, month, day, 21)] = int(x[52:55])

                daily_kp[datetime(year, month, day)] = int(x[28:31])
                daily_ap[datetime(year, month, day)] = int(x[55:58])

                try:
                    temp = float(x[65:70])  # f107
                except ValueError:
                    temp = float('NaN')  # If the string is empty, just use NaN

                if temp == 0.:  # Replace 0's of f107 with NaN
                    temp = float('NaN')

                f107[datetime(year, month, day)] = temp

            f.close()

    # Caculate f107a:
    for dn, value in f107.items():
        f107_81values = []
        for dday in range(-40, 41):  # 81 day sliding window
            dt = timedelta(dday)
            try:
                f107_81values.append(f107[dn+dt])
            except KeyError:
                f107_81values.append(float('NaN'))
        f107a[dn] = np.nan if all(np.isnan(f107_81values)) else \
            np.nanmean(f107_81values)

    # Read in DST values.
    # (1) Read old files that were shipped with pyglow.
    # (2) Search the dst/ folder for all new files, and read them.
    oldies = ['1957_1969', '1970_1989', '1990_2004']
    dst_path = '%s/dst/' % PYGLOW_PATH
    files = glob.glob('%s??????' % dst_path)  # files like 201407
    # Older files listed above:
    old_files = ['%s%s' % (dst_path, old) for old in oldies]
    files.extend(old_files)  # a list of every dst file
    for fn in files:
        with open(fn, 'r') as f:
            s = f.readlines()
        for x in s:
            if len(x) <= 1:
                break  # reached last line. Done with this file.
            yr23 = x[3:5]  # 3rd and 4th digits of year
            month = int(x[5:7])
            day = int(x[8:10])
            yr12 = x[14:16]  # 1st and 2nd digits of year
            base = int(x[16:20])  # "Base value, unit 100 nT"
            year = int('%02s%02s' % (yr12, yr23))
            dst_per_hour = np.zeros(24)
            for i in range(24):
                dsthr = base*100 + int(x[20+4*i:24+4*i])
                if dsthr == 9999:
                    dsthr = np.nan
                dst_per_hour[i] = dsthr
            dst[datetime(year, month, day)] = dst_per_hour

    # Read in AE values.
    ae_path = '%s/ae/' % PYGLOW_PATH
    files = glob.glob('%s*' % ae_path)  # find files like 1975
    for fn in files:
        with open(fn, 'r') as f:
            s = f.readlines()
        for x in s:
            if len(x) <= 1:
                break  # reached last line. Done with this file.
            yr23 = int(x[0:2])  # 3rd and 4th digits of year
            month = int(x[2:4])
            day = int(x[4:6])
            year = int(fn.split("/")[-1][:4])
            # in case file contains an extra day:
            if yr23 != int(str(year)[2:4]):
                break
            hour = int(x[6:8])
            if hour == 0:
                ae_per_hour = np.zeros(24)
            aehr = int(x[8:14])
            if aehr == 99999:
                aehr = np.nan
            ae_per_hour[hour] = aehr
            ae[datetime(year, month, day)] = ae_per_hour

    """ Part 2: placing indices in 'geophysical_indices' array """
    # time to put the values into a numpy array, geophysical_indices
    end_day = datetime.today()
    total_days = (end_day-EPOCH).days+1
    geophysical_indices = np.zeros((68, total_days))*float('nan')

    i = 0
    while i < total_days:  # Try every day. Some will be NaN.
        dn = EPOCH + timedelta(i)  # Increment a day

        try:  # This will fail if kp/ap data doesn't exist on this day:
            geophysical_indices[0, i] = \
                kp[datetime(dn.year, dn.month, dn.day, 0)]
            geophysical_indices[1, i] = \
                kp[datetime(dn.year, dn.month, dn.day, 3)]
            geophysical_indices[2, i] = \
                kp[datetime(dn.year, dn.month, dn.day, 6)]
            geophysical_indices[3, i] = \
                kp[datetime(dn.year, dn.month, dn.day, 9)]
            geophysical_indices[4, i] = \
                kp[datetime(dn.year, dn.month, dn.day, 12)]
            geophysical_indices[5, i] = \
                kp[datetime(dn.year, dn.month, dn.day, 15)]
            geophysical_indices[6, i] = \
                kp[datetime(dn.year, dn.month, dn.day, 18)]
            geophysical_indices[7, i] = \
                kp[datetime(dn.year, dn.month, dn.day, 21)]

            geophysical_indices[8, i] = \
                ap[datetime(dn.year, dn.month, dn.day, 0)]
            geophysical_indices[9, i] = \
                ap[datetime(dn.year, dn.month, dn.day, 3)]
            geophysical_indices[10, i] = \
                ap[datetime(dn.year, dn.month, dn.day, 6)]
            geophysical_indices[11, i] = \
                ap[datetime(dn.year, dn.month, dn.day, 9)]
            geophysical_indices[12, i] = \
                ap[datetime(dn.year, dn.month, dn.day, 12)]
            geophysical_indices[13, i] = \
                ap[datetime(dn.year, dn.month, dn.day, 15)]
            geophysical_indices[14, i] = \
                ap[datetime(dn.year, dn.month, dn.day, 18)]
            geophysical_indices[15, i] = \
                ap[datetime(dn.year, dn.month, dn.day, 21)]

            geophysical_indices[18, i] = \
                daily_kp[datetime(dn.year, dn.month, dn.day)]
            geophysical_indices[19, i] = \
                daily_ap[datetime(dn.year, dn.month, dn.day)]

        except KeyError:
            pass

        try:  # This will fail if no f10.7 data are available on this day
            geophysical_indices[16, i] = \
                f107[datetime(dn.year, dn.month, dn.day)]
            geophysical_indices[17, i] = \
                f107a[datetime(dn.year, dn.month, dn.day)]
        except KeyError:
            pass

        try:  # This will fail if no dst data are available on this day
            geophysical_indices[20:44, i] = dst[dn]
        except KeyError:
            pass

        try:  # This will fail if no ae data are available on this day
            geophysical_indices[44:68, i] = ae[dn]
        except KeyError:
            pass

        i = i + 1

    # Update file cache:
    with open(MTIME_TABLE_FNAME, 'wb') as fid:
        dump(mtime_table, fid, -1)
    np.save(GEOPHYSICAL_INDICES_FNAME, geophysical_indices)
    print("[generate_kpap.py] Generated: {}".format(GEOPHYSICAL_INDICES_FNAME))

    return geophysical_indices


def fetch():
    """
    Main interface to retrieve geophysical indices

    :return geophysical_indices: data structure containing geophysical indices

    """

    # Do we need to update / generate a file?
    update = update_required()

    if update:
        # Fetch indices (it will also save a file):
        geophysical_indices = generate_kpap()
    else:
        # Update not required, load cached indices:
        geophysical_indices = np.load(GEOPHYSICAL_INDICES_FNAME)

    return geophysical_indices
