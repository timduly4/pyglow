
#def generate_kpap():
"""
Description:
------------
This script parses the geophysical indices from the yearly files located in
./kpap/ and ./dst/ .  The files are obtained by running pyglow.update_indices()

After parsing the files, a variable is created, "geophysical_indices", whose
dimensions are 20 x (number of days), with

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
An easy way to grab the yearly data, in python of course:
(make sure the folder "kpap" is created in the current directory):

import os
for y in range(1932,end_year):
    os.system('wget ftp://ftp.ngdc.noaa.gov/STP/GEOMAGNETIC_DATA/INDICES/KP_AP/%4i ./kpap/' % y)

--------
History:
--------
    7/21/12 : Created, Timothy Duly (duly2@illinois.edu)
    1/07/15 : Changed end_year so that it doesn't crash during
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


from datetime import datetime, timedelta, date
from cPickle import load, dump
import numpy as np
import os
from numpy import nanmean
import sys
import pyglow
import glob


# Load all data up to and including this year, if it exists
end_year = date.today().year + 1

epoch = datetime(1932,1,1)

pyglow_path = '/'.join(pyglow.__file__.split("/")[:-1])

"""File to store table of index file modification times"""
mtime_table_fname = os.path.join(pyglow_path, 'mtime_table.pkl')

"""File to store cached geophysical_indices array"""
geophysical_indices_fname = os.path.join(pyglow_path, 'geophysical_indices.npy')


def get_mtime_table():
    mtime_table = {}
    # kpap files
    for y in range(1932,end_year):
        fn = pyglow_path + "/kpap/%4i" % y
        if os.path.isfile(fn):
            mtime_table[fn] = os.path.getmtime(fn)
    # dst files
    oldies = ['1957_1969','1970_1989','1990_2004']
    dst_path = '%s/dst/' % pyglow_path
    files = glob.glob('%s??????' % dst_path) # files like 201407
    old_files = ['%s%s' % (dst_path, old) for old in oldies] # older files listed above
    files.extend(old_files) # a list of every dst file
    for fn in files:
        mtime_table[fn] = os.path.getmtime(fn)
    # ae files
    ae_path = '%s/ae/' % pyglow_path
    files = glob.glob('%s*' % ae_path) # find files like 1975
    for fn in files:
        mtime_table[fn] = os.path.getmtime(fn)
    return mtime_table


def update_required():
    if not os.path.isfile(mtime_table_fname):
        return True
    with open(mtime_table_fname) as fid:
        mtime_table = load(fid)
    new_mtime_table = get_mtime_table()
    for key, val in new_mtime_table.iteritems():
        try:
            if mtime_table[key] != new_mtime_table[key]:
                return True
        except IndexError:
            return True
    # all current file modification times matched those associated
    # with cached geophysical_indices
    return False


if update_required():
    print('Updating geophysical indices table')

    """Files have changed --- update geophysical indices"""
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

    for y in range(1932,end_year):
        #f = open(os.getcwd() + "/kpap/%4i" % y)
        #f = open('/Users/duly/Dropbox/research/f2py/model_atmosphere/kpap/%4i' % y)
        # pyglow_path = '/'.join(pyglow.__file__.split("/")[:-1])
        fn = pyglow_path + "/kpap/%4i" % y
        if os.path.isfile(fn): # If the file has been downloaded
            f = open(fn)
            for x in f.readlines():
                year = int(x[0:2]) # parse the line for the year, only last 2 digits

                if year < 30: # need to change this is 2030....
                    year = year + 2000
                else:
                    year = year + 1900

                month = int(x[2:4]) # parse the line for month
                day = int(x[4:6])   # ... and the days

                # parse the values for kp, ap, daily_kp, and daily_ap:
                kp[ datetime(year,month,day,0)  ] = int(x[12:14])/10.
                kp[ datetime(year,month,day,3)  ] = int(x[14:16])/10.
                kp[ datetime(year,month,day,6)  ] = int(x[16:18])/10.
                kp[ datetime(year,month,day,9)  ] = int(x[18:20])/10.
                kp[ datetime(year,month,day,12) ] = int(x[20:22])/10.
                kp[ datetime(year,month,day,15) ] = int(x[22:24])/10.
                kp[ datetime(year,month,day,18) ] = int(x[24:26])/10.
                kp[ datetime(year,month,day,21) ] = int(x[26:28])/10.

                ap[ datetime(year,month,day,0)  ] = int(x[31:34])
                ap[ datetime(year,month,day,3)  ] = int(x[34:37])
                ap[ datetime(year,month,day,6)  ] = int(x[37:40])
                ap[ datetime(year,month,day,9)  ] = int(x[40:43])
                ap[ datetime(year,month,day,12) ] = int(x[43:46])
                ap[ datetime(year,month,day,15) ] = int(x[46:49])
                ap[ datetime(year,month,day,18) ] = int(x[49:52])
                ap[ datetime(year,month,day,21) ] = int(x[52:55])

                daily_kp[ datetime(year,month,day) ] = int(x[28:31])
                daily_ap[ datetime(year,month,day) ] = int(x[55:58])

                try:
                    temp = float(x[65:70]) # f107
                except:
                    temp = float('NaN') # if the string is empty, just use NaN

                if temp == 0.: # replace 0's of f107 with NaN
                    temp = float('NaN')

                f107[ datetime(year,month,day) ]  = temp


            f.close()

    # Caculate f107a:
    for dn, value in f107.items():
        f107_81values = []
        for dday in range(-40,41): # 81 day sliding window
            dt = timedelta(dday)
            try:
                f107_81values.append(f107[dn+dt])
            except:
                f107_81values.append(float('NaN'))
        f107a[dn] = nanmean(f107_81values)

    # Read in DST values.
    # (1) Read old files that were shipped with pyglow.
    # (2) Search the dst/ folder for all new files, and read them.
    oldies = ['1957_1969','1970_1989','1990_2004']
    # pyglow_path = '/'.join(pyglow.__file__.split("/")[:-1])
    dst_path = '%s/dst/' % pyglow_path
    files = glob.glob('%s??????' % dst_path) # files like 201407
    old_files = ['%s%s' % (dst_path, old) for old in oldies] # older files listed above
    files.extend(old_files) # a list of every dst file
    for fn in files:
        with open(fn,'r') as f:
            s = f.readlines()
        for x in s:
            if len(x) <= 1:
                break # reached last line. Done with this file.
            yr23 = x[3:5] # 3rd and 4th digits of year
            month = int(x[5:7])
            day   = int(x[8:10])
            yr12  = x[14:16] # 1st and 2nd digits of year
            base  = int(x[16:20]) # "Base value, unit 100 nT"
            year = int('%02s%02s' % (yr12,yr23))
            dst_per_hour = np.zeros(24)
            for i in range(24):
                dsthr = base*100 + int(x[20+4*i:24+4*i])
                if dsthr==9999:
                    dsthr = np.nan
                dst_per_hour[i] = dsthr
            dst[datetime(year,month,day)] = dst_per_hour


    # Read in AE values.
    # pyglow_path = '/'.join(pyglow.__file__.split("/")[:-1])
    ae_path = '%s/ae/' % pyglow_path
    files = glob.glob('%s*' % ae_path) # find files like 1975
    for fn in files:
        with open(fn,'r') as f:
            s = f.readlines()
        for x in s:
            if len(x) <= 1:
                break # reached last line. Done with this file.
            yr23 = int(x[0:2]) # 3rd and 4th digits of year
            month = int(x[2:4])
            day   = int(x[4:6])
            year = int(fn.split("/")[-1][:4])
            # in case file contains an extra day:
            if yr23 != int(str(year)[2:4]):
                break
            hour = int(x[6:8])
            if hour == 0:
                ae_per_hour = np.zeros(24)
            aehr = int(x[8:14])
            if aehr==99999:
                aehr = np.nan
            ae_per_hour[hour] = aehr
            ae[datetime(year,month,day)] = ae_per_hour



    """ Part 2: placing indices in 'geophysical_indices' array """

    # time to put the values into a numpy array, geophysical_indices
    end_day = datetime.today()
    total_days = (end_day-epoch).days+1
    geophysical_indices = np.zeros((68,total_days))*float('nan')

    i = 0
    while i < total_days: # Try every day. Some will be nan.
        dn = epoch + timedelta(i) # increment a day

        try: # This will fail if kp/ap data doesn't exist on this day
            geophysical_indices[0,i]  = kp[ datetime(dn.year, dn.month, dn.day, 0)]
            geophysical_indices[1,i]  = kp[ datetime(dn.year, dn.month, dn.day, 3)]
            geophysical_indices[2,i]  = kp[ datetime(dn.year, dn.month, dn.day, 6)]
            geophysical_indices[3,i]  = kp[ datetime(dn.year, dn.month, dn.day, 9)]
            geophysical_indices[4,i]  = kp[ datetime(dn.year, dn.month, dn.day, 12)]
            geophysical_indices[5,i]  = kp[ datetime(dn.year, dn.month, dn.day, 15)]
            geophysical_indices[6,i]  = kp[ datetime(dn.year, dn.month, dn.day, 18)]
            geophysical_indices[7,i]  = kp[ datetime(dn.year, dn.month, dn.day, 21)]

            geophysical_indices[8,i]  = ap[ datetime(dn.year, dn.month, dn.day, 0)]
            geophysical_indices[9,i]  = ap[ datetime(dn.year, dn.month, dn.day, 3)]
            geophysical_indices[10,i] = ap[ datetime(dn.year, dn.month, dn.day, 6)]
            geophysical_indices[11,i] = ap[ datetime(dn.year, dn.month, dn.day, 9)]
            geophysical_indices[12,i] = ap[ datetime(dn.year, dn.month, dn.day, 12)]
            geophysical_indices[13,i] = ap[ datetime(dn.year, dn.month, dn.day, 15)]
            geophysical_indices[14,i] = ap[ datetime(dn.year, dn.month, dn.day, 18)]
            geophysical_indices[15,i] = ap[ datetime(dn.year, dn.month, dn.day, 21)]

            geophysical_indices[18,i] = daily_kp[ datetime(dn.year, dn.month, dn.day)]
            geophysical_indices[19,i] = daily_ap[ datetime(dn.year, dn.month, dn.day)]
        except KeyError:
            pass

        try: # This will fail if no f10.7 data are available on this day
            geophysical_indices[16,i] = f107[ datetime(dn.year, dn.month, dn.day)]
            geophysical_indices[17,i] = f107a[ datetime(dn.year, dn.month, dn.day)]
        except KeyError:
            pass

        try: # This will fail if no dst data are available on this day
            geophysical_indices[20:44,i] = dst[dn]
        except KeyError:
            pass

        try: # This will fail if no ae data are available on this day
            geophysical_indices[44:68,i] = ae[dn]
        except KeyError:
            pass

        i = i + 1

    """Update file cache"""
    with open(mtime_table_fname, 'w') as fid:
        dump(mtime_table, fid, -1)
    np.save(geophysical_indices_fname, geophysical_indices)
else:
    """Update not required --- load cached indices"""
    geophysical_indices = np.load(geophysical_indices_fname)


# print "saving geophysical indices"
# np.savez('kpap.npz', \
#     geophysical_indices = geophysical_indices, \
#     epoch = epoch, \
#     )
# need this for later: ?
####np.save('geophysical_indices',geophysical_indices)





# if __name__ == '__main__':
#     generate_kpap()


'''
# check to make sure we have everything
dn_start = datetime(1932,1,1)
dn_end = datetime(2010,1,1)

no_days = (dn_end-dn_start).days

kk = 0

while kk < no_days:
    dn = dn_start + timedelta(kk)
    for jj in [0, 3, 6, 9, 12, 15, 18, 21]:
        print dn+timedelta(jj/24.),kp[dn+timedelta(jj/24.)],ap[dn+timedelta(jj/24.)],f107[dn],f107a[dn],daily_kp[dn],daily_ap[dn]
    kk+=1
'''
