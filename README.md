![alt text](https://remote2.csl.illinois.edu/~duly/pyglow/logo.png "pyglow")

1. Download the climatological models and wrap them with f2py

    $ cd ./pyglow/models/
    $ make all

    If successful, there should be a *.so file in each of the
    ./models/dl_models/<model>/ directories

2. Install the python package

    $ python ./setup.py install (--user)

   On a mac, the folder "pyglow" and *.so files from
   ./models/dl_models/<model>/ should go in

   /Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages

3.  Test

    In Python, run:

from pyglow.pyglow import Point

dn = datetime(2010,3,1, 9,30)
lat = 0.
lon = -80.
alt = 250.

pt = Point(dn, lat, lon, alt)

pt.run_igrf()
pt.run_hwm93()
pt.run_msis()
pt.run_iri()

print pt
