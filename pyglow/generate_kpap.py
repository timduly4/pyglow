
#def generate_kpap():
"""
Description:
------------
This script parses the geophysical indices from the yearly files located in
./kpap/ .  The files are obtained from
 ftp://ftp.ngdc.noaa.gov/STP/GEOMAGNETIC_DATA/INDICES/KP_AP/

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
"""

# TODO
# [ ] Is dividing by 10 correct for kp?
# [ ] daily kp -- should it be the average of sum? right now it's the sum...
# [ ] fix path issue... how to make it callable anywhere?


from datetime import datetime, timedelta, date
import numpy as np
import os
from scipy.stats import nanmean
import sys
import pyglow


""" Part 1: Parsing the raw data files """

end_year = date.today().year + 1

# Create empty dictionaries:
ap = {}    
kp = {}
f107 = {}
daily_kp = {}
daily_ap = {}
f107a = {}

for y in range(1932,end_year):
    #f = open(os.getcwd() + "/kpap/%4i" % y)
    #f = open('/Users/duly/Dropbox/research/f2py/model_atmosphere/kpap/%4i' % y)
    pyglow_path = '/'.join(pyglow.__file__.split("/")[:-1])
    f = open(pyglow_path + "/kpap/%4i" % y)
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
        
""" Part 2: placing indices in 'geophysical_indices' array """

# time to put the values into an numpy array:
geophysical_indices = np.zeros((20,len(f107)))*float('nan')

i = 0
epoch = datetime(1932,1,1)
while i < len(f107): # assume that len(f107) = number of days parsed
    dn = epoch + timedelta(i) # increment a day
    
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

    geophysical_indices[16,i] = f107[ datetime(dn.year, dn.month, dn.day)]
    geophysical_indices[17,i] = f107a[ datetime(dn.year, dn.month, dn.day)]

    geophysical_indices[18,i] = daily_kp[ datetime(dn.year, dn.month, dn.day)]
    geophysical_indices[19,i] = daily_ap[ datetime(dn.year, dn.month, dn.day)]

    i = i + 1

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

