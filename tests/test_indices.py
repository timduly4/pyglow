#!/usr/bin/env python

from datetime import datetime, timedelta

import pyglow

dns = [datetime(2010, 3, 23, 15, 30), 
        datetime(2015, 3, 23, 15, 30)]


for dn in dns:
    kp, ap, f107, f107a, daily_kp, daily_ap, dst, ae = \
            pyglow.get_kpap.get_kpap(dn)
    
    print("\n{}".format(dn))
    print('-'*len(str(dn)))
    
    print("kp = {}".format(kp))
    print("ap = {}".format(ap))
    print("f107 = {}".format(f107))
    print("f107a = {:3.2f}".format(f107a))
    print("daily_kp = {}".format(daily_kp))
    print("daily_ap = {}".format(daily_ap))
    print("dst = {}".format(dst))
    print("ae = {}".format(ae))

