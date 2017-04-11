![alt text](https://raw.github.com/timduly4/pyglow/master/logo.png "pyglow")
[_(airglow viewed aboard the ISS)_](http://en.wikipedia.org/wiki/File:Cupola_above_the_darkened_Earth.jpg)

# Overview

pyglow is a Python module that wraps several upper atmosphere climatoglogical
models written in FORTRAN, such as the Horizontal Wind Model (HWM), the
International Geomagnetic Reference Field (IGRF), the International Reference
Ionosphere (IRI), and the Mass Spectrometer and Incoherent Scatter Radar (MSIS).

It includes the following upper atmospheric models:

  * HWM 2014
  * HWM 2007
  * HWM 1993
  * IGRF 12
  * IGRF 11
  * IRI 2016
  * IRI 2012
  * MSIS 2000

pyglow also provides access to the the following geophysical indices:
  * AP
  * Kp
  * F10.7
  * DST
  * AE

pyglow offers access to these models & indices in a convenient, high-level object-oriented interface within Python.

# Installation

### I'm Feeling Lucky

First, checkout the repository:

```
    $ git clone git://github.com/timduly4/pyglow.git;
```

Change directories into the repository folder and run the installation script:
```
    $ cd pyglow/
    $ ./pyglow_install.sh
```

### Individual installation steps

If you have troubles, follow the individual installation steps:

(1) Download the package:
```
    $ git clone git://github.com/timduly4/pyglow.git
    $ cd pyglow/
```

(2) Download the climatological models and wrap them with f2py:
```
    $ cd ./pyglow/models/
    $ make all
```
  * If successful, there should be a `*.so` file in each of the
    `./models/dl_models/<model>/` directories:

    ```
    $ find . -name "*.so"
    ./dl_models/hwm07/hwm07py.so
    ./dl_models/hwm93/hwm93py.so
    ./dl_models/hwm14/hwm14py.so
    ./dl_models/igrf11/igrf11py.so
    ./dl_models/igrf12/igrf12py.so
    ./dl_models/iri12/iri12py.so
    ./dl_models/iri16/iri16py.so
    ./dl_models/msis/msis00py.so
    ```

(3) Install the python package
```
    $ cd ../../   # get back to root directory
    $ python ./setup.py install --user
```
  * On a mac, the folder `pyglow` and `*.so` files from
    `./models/dl_models/<model>/` should be in
    `/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages`
  * If you are denied permission, I recommend adding `--user` flag in command

(4) Download the geophysical indices

```
	$ cd ~/
	$ python -c "from pyglow import pyglow; pyglow.update_indices()"
```


# Testing / Example

* In Python, run:

```
from pyglow.pyglow import Point
import datetime as dt

dn = dt.datetime(2011, 3, 23, 9, 30)
lat = 0.0
lon = -80.0
alt = 250.0

pt = Point(dn, lat, lon, alt)

print "Before running any models:"
print pt

pt.run_igrf()
pt.run_hwm93()
pt.run_msis()
pt.run_iri()
pt.run_airglow()

print "After running models:"
pt.flatten()
test_line = "IGRF {:d}:\n(Bx, By, Bz) = ({:.4g}, {:.4g}, {:.4g})\n".format(pt.igrf_version, pt.By, pt.Bx, pt.Bz)
test_line = "{:s}{:>15s}{:.4g}\n".format(test_line, "B = ", pt.B)
test_line = "{:s}{:>15s}{:.4g}\n\n".format(test_line, "dip = ", pt.dip)
test_line = "{:s}{:>15s}{:d}\n".format(test_line, "HWM Version = ", pt.hwm_version)
test_line = "{:s}{:>15s}({:.2f}, {:.2f})\n\n".format(test_line, "(u, v) = ", pt.u, pt.v)
test_line = "{:s}{:>15s}{:d}\n".format(test_line, "IRI Version = ", pt.iri_version)
test_line = "{:s}{:>15s}({:.2f}, {:.4g})\n\n".format(test_line, "(hmF2, NmF2) = ", pt.hmF2, pt.NmF2)
test_line = "{:s}{:>15s}{:d}\n".format(test_line, "MSIS Version = ", pt.msis_version)
test_line = "{:s}{:>15s}{:.4g}\n\n".format(test_line, "Tn = ", pt.Tn_msis)
test_line = "{:s}{:>15s}{:.4g}\n".format(test_line, "Airglow at 6300 = ", pt.ag6300)
test_line = "{:s}{:>15s}{:.4g}".format(test_line, "Airglow at 7774 = ", pt.ag7774)
print test_line
```

* The output should be as follows:

```
Before running any models:
                    date and time = 2011-03-23 09:30:00
lat, lon, alt min, max, step (km) = 0.00, -80.00, 250.00, 250.00, 10.00

Geophysical Indices:
---------------------
                               kp = 2.70
                               ap = 12.00
                             f107 = 104.00
                            f107a = 110.47
                              Dst = nan
                               ae = nan

Model Data Available For:
-----------
                             IGRF = False
                              HWM = False
                              IRI = False
                             MSIS = False
                          Airglow = False

# After running each model, an updated status will print, ending with:

                    date and time = 2011-03-23 09:30:00
lat, lon, alt min, max, step (km) = 0.00, -80.00, 250.00, 250.00, 10.00

Geophysical Indices:
---------------------
                               kp = 2.70
                               ap = 12.00
                             f107 = 104.00
                            f107a = 110.47
                              Dst = nan
                               ae = nan

Model Data Available For:
-----------
                             IGRF = 12
                              HWM = 1993
                              IRI = 2016
                             MSIS = 2000
                          Airglow = True

After running models:

IGRF 12:
(Bx, By, Bz) = (2.448e-05, -6.433e-07, -9.925e-06)
           B = 2.643e-05
         dip = 22.06

 HWM Version = 1993
      (u, v) = (14.07, 12.18)

 IRI Version = 2016
(hmF2, NmF2) = (266.28, 7.869e+04)

MSIS Version = 2000
          Tn = 779.9

Airglow at 6300 = 4.032
Airglow at 7774 = 0.00749
```


# Hints

### General
1. Use tab completion in ipython to view the full set of member data and
   variables available in the Point class.
  * For example, in the test code, run `pt.<TAB><TAB>` and class information
    will be listed.

### Updating geophysical indices with `update_indices()`
1. You'll need to download the geophysical indices as they become available.
   The `update_indices()` function is available in pyglow that enables you do
   this:

```
from pyglow.pyglow import update_indices
update_indices([2012, 2013]) # grabs indices for 2012 and 2013
update_indices() # grabs all indices starting from 1932 to the current year
```

  * Only need to run this function when you would like to update the indices.


# Uninstallation 

1. The install directory for pyglow is outputted when you run the
   `python ./setup.py install` command.  For example, for macs this is usually
    `/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages`.
2.  Simply remove the `*.so` climatological models in this directory, as well
    as the `pyglow` and `pyglow_trash` folders.

