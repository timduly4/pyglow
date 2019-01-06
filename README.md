![alt text](https://raw.github.com/timduly4/pyglow/master/logo.png "pyglow")
[_(airglow viewed aboard the ISS)_](http://en.wikipedia.org/wiki/File:Cupola_above_the_darkened_Earth.jpg)

# Overview

`pyglow` is a Python module that wraps several upper atmosphere climatological models written in FORTRAN, such as the Horizontal Wind Model (HWM), the International Geomagnetic Reference Field (IGRF), the International Reference Ionosphere (IRI), and the Mass Spectrometer and Incoherent Scatter Radar (MSIS).

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

`pyglow` offers access to these models & indices in a convenient, high-level object-oriented interface within Python.

# Prerequisites

`pyglow` requires the following packages for installation:

1. `gfortran` (`$ sudo apt-get install gfortran`)
2. `f2py` (`$ pip install numpy --upgrade`)

# Installation

### I'm Feeling Lucky:

First, checkout the repository:

```
    $ git clone git://github.com/timduly4/pyglow.git pyglow
```

Change directories into the repository folder, compile the f2py bindings, then install:
```
    $ cd pyglow/
    $ make -C src/pyglow/models source
    $ python setup.py install --user
```

### Individual installation steps:

If you have troubles, follow the individual installation steps:

(1) Download the package:
```
    $ git clone git://github.com/timduly4/pyglow.git
    $ cd pyglow/
```

(2) Download the climatological models and wrap them with f2py:
```
    $ cd ./src/pyglow/models/
    $ make all
```
  * If successful, there should be a `*.so` file in each of the `./models/dl_models/<model>/` directories:

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
    $ cd ../../../   # get back to root directory
    $ python setup.py install --user --prefix=
```
  * On a mac, the folder `pyglow` and `*.so` files from `./models/dl_models/<model>/` should be in `/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages`
  * If you are denied permission, I recommend adding `--user` flag in command

(4) Download the geophysical indices

```
	$ cd ~/
	$ python -c "import pyglow; pyglow.update_indices()"
```


# Testing / Examples

See unit tests in `./test`.  For example, run the unittest suite with:

`$ python -m unittest test.test_suite_pyglow`)

(However, be sure that the f2py modules have been compiled via `$ make -C src/pyglow/models source`, first.)

See example scripts located in `./examples` for example calls to `pyglow`.

# Hints

### General
1. Use tab completion in ipython to view the full set of member data and variables available in the Point class.
  * For example, in the test code, run `pt.<TAB><TAB>` and class information will be listed.

### Updating geophysical indices with `pyglow.update_indices()`
1. You'll need to download the geophysical indices as they become available.  The `update_indices()` function is available in pyglow that enables you do this:

```
# Grabs indices for 2017 and 2018:
~ $ python -c "import pyglow; pyglow.update_indices([2017, 2018])"

# Grabs all indices starting from 1932 to the current year:
~ $ python -c "import pyglow; pyglow.update_indices()"
```

  * Note: you only need to run this function when you would like to update the indices.


# Uninstallation

The install directory for pyglow can be outputted via `python -c "import pyglow; print(pyglow.__file__)"`.  For example:
```
~ $ python -c "import pyglow; print(pyglow.__file__)"
/Users/duly/Library/Python/2.7/lib/python/site-packages/pyglow/__init__.pyc
```
This tells you the installation location, and then you can remove the package with:
```
~ $ rm -rf /Users/duly/Library/Python/2.7/lib/python/site-packages/pyglow
```
