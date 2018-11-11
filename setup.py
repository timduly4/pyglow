#!/usr/bin/env python

import os
from io import open  # For Python 2 compatibility
from numpy.distutils.core import setup, Extension

DL_MODELS = 'src/pyglow/models/dl_models'


def reencode(dosfile, target='utf-8'):
    """
    Remove invalid unicode characters that appear in the comments
    """
    with open(dosfile, 'r', encoding='cp1251', errors='ignore') as f:
        content = f.read()
    with open(dosfile, 'w', encoding=target) as f:
        f.write(content)

    return


# IGRF 11:
igrf11 = Extension(
    name='igrf11py',
    sources=[os.path.join(DL_MODELS, 'igrf11', fname) for fname in [
        'igrf11_modified.f',
        'sig_file_patched.pyf',
    ]],
)

# IGRF 12:
igrf12 = Extension(
    name='igrf12py',
    sources=[os.path.join(DL_MODELS, 'igrf12', fname) for fname in [
        'igrf12_modified.f',
        'sig_file_patched.pyf',
    ]],
)

# HWM 93:
hwm93 = Extension(
    name='hwm93py',
    sources=[os.path.join(DL_MODELS, 'hwm93', fname) for fname in [
        'hwm93_modified.f',
        'sig_file_patched.pyf',
    ]],
    extra_f77_compile_args=['-std=legacy'],
)

# Rencode HWM07 sources:
hwm07_sources = [os.path.join(DL_MODELS, 'hwm07', fname) for fname in [
    'hwm07e_modified.f90',
    'apexcord.f90',
]]
for source in hwm07_sources:
    reencode(source)

# Use the makefile to generate a signature
os.system('make -Cpyglow/models/dl_models/hwm07 sig')

# HWM07:
hwm07 = Extension(
    name='hwm07py',
    sources=hwm07_sources + [os.path.join(DL_MODELS, 'hwm07', 'sig_file.pyf')],
    # f2py_options=['only: hwmqt :'],  # where is the right place to put this?
)

# HWM14:
hwm14 = Extension(
    name='hwm14py',
    sources=[os.path.join(DL_MODELS, 'hwm14', fname) for fname in [
        'hwm14.f90',
        'sig_file.pyf',
    ]],
    extra_f77_compile_args=['-std=legacy'],
)

# IRI 12:
iri12 = Extension(
    name='iri12py',
    sources=[os.path.join(DL_MODELS, 'iri12', fname) for fname in [
        'cira.for',
        'igrf.for',
        'iridreg_modified.for',
        'irifun.for',
        'irisub.for',
        'iritec.for',
        'iriflip.for',
        'sig_file_patched.pyf',
    ]],
    extra_f77_compile_args=[
        '-std=legacy',
        '-w',
        '-O2',
        '-fbacktrace',
        '-fno-automatic',
        '-fPIC',
    ],
)

# IRI16:
iri16 = Extension(
    name='iri16py',
    sources=[os.path.join(DL_MODELS, 'iri16', fname) for fname in [
        'cira.for',
        'igrf.for',
        'iridreg_modified.for',
        'irifun.for',
        'irisub.for',
        'iritec.for',
        'iriflip_modified.for',
        'cosd_sind.for',
        'sig_file_patched.pyf',
    ]],
    extra_f77_compile_args=[
        '-std=legacy',
        '-w',
        '-O2',
        '-fbacktrace',
        '-fno-automatic',
        '-fPIC',
    ],
)

# MSIS00:
msis00 = Extension(
    name='msis00py',
    sources=[os.path.join(DL_MODELS, 'msis', fname) for fname in [
        'nrlmsise00_sub_patched.for',
        'sig_file_patched.pyf'
    ]],
    extra_f77_compile_args=['-std=legacy'],
)

# Distutils setup:
setup(
    name='pyglow',
    url='https://github.com/timduly4/pyglow',
    author='Timothy M. Duly',
    author_email='timduly4@gmail.com',
    packages=['src/pyglow'],
    ext_modules=[
        igrf11,
        igrf12,
        hwm93,
        hwm07,
        hwm14,
        iri12,
        iri16,
        msis00,
    ],
    data_files=[
        ('pyglow_trash', ['src/pyglow/models/Makefile']),
        ('pyglow_trash', ['src/pyglow/models/get_models.py']),
        ('pyglow_trash', ['src/pyglow/models/dl_models/hwm07/dummy.txt']),
        ('pyglow_trash', ['src/pyglow/models/dl_models/hwm93/dummy.txt']),
        ('pyglow_trash', ['src/pyglow/models/dl_models/igrf11/dummy.txt']),
        ('pyglow_trash', ['src/pyglow/models/dl_models/igrf12/dummy.txt']),
        ('pyglow_trash', ['src/pyglow/models/dl_models/iri12/dummy.txt']),
        ('pyglow_trash', ['src/pyglow/models/dl_models/iri16/dummy.txt']),
        ('pyglow_trash', ['src/pyglow/models/dl_models/msis/dummy.txt']),
        ('pyglow_trash', ['src/pyglow/models/dl_models/hwm14/Makefile']),
        ('pyglow_trash', ['src/pyglow/models/f2py/hwm07/hwm07e.patch']),
        ('pyglow_trash', ['src/pyglow/models/f2py/hwm07/Makefile']),
        ('pyglow_trash', ['src/pyglow/models/f2py/hwm93/hwm93.patch']),
        ('pyglow_trash', ['src/pyglow/models/f2py/hwm93/Makefile']),
        ('pyglow_trash', ['src/pyglow/models/f2py/hwm93/sig.patch']),
        ('pyglow_trash', ['src/pyglow/models/f2py/igrf11/igrf11.patch']),
        ('pyglow_trash', ['src/pyglow/models/f2py/igrf11/Makefile']),
        ('pyglow_trash', ['src/pyglow/models/f2py/igrf11/sig.patch']),
        ('pyglow_trash', ['src/pyglow/models/f2py/igrf12/igrf12.patch']),
        ('pyglow_trash', ['src/pyglow/models/f2py/igrf12/Makefile']),
        ('pyglow_trash', ['src/pyglow/models/f2py/igrf12/sig.patch']),
        ('pyglow_trash', [
            'src/pyglow/models/f2py/iri12/delete_iriflip_comments.py'
        ]),
        ('pyglow_trash', ['src/pyglow/models/f2py/iri12/Makefile']),
        ('pyglow_trash', ['src/pyglow/models/f2py/iri12/sig.patch']),
        ('pyglow_trash', ['src/pyglow/models/f2py/iri12/iridreg.patch']),
        ('pyglow_trash', ['src/pyglow/models/f2py/msis/Makefile']),
        ('pyglow_trash', ['src/pyglow/models/f2py/msis/nrlmsise00_sub.patch']),
        ('pyglow_trash', ['src/pyglow/models/f2py/msis/sig.patch']),
        ('pyglow/hwm07_data/', [
            'src/pyglow/models/dl_models/hwm07/apexgrid.dat',
            'src/pyglow/models/dl_models/hwm07/dwm07b_104i.dat',
            'src/pyglow/models/dl_models/hwm07/hwm071308e.dat',
        ]),
        (
            'pyglow/hwm14_data/',
            [
                'src/pyglow/models/dl_models/hwm14/gd2qd.dat',
                'src/pyglow/models/dl_models/hwm14/dwm07b_104i.dat',
                'src/pyglow/models/dl_models/hwm14/hwm14-beta.bin',
                'src/pyglow/models/dl_models/hwm14/hwm14.f90',
            ],
        ),
        (
            'pyglow/iri12_data/',
            [
                'src/pyglow/models/dl_models/iri12/apf107.dat',
                'src/pyglow/models/dl_models/iri12/ccir11.asc',
                'src/pyglow/models/dl_models/iri12/ccir12.asc',
                'src/pyglow/models/dl_models/iri12/ccir13.asc',
                'src/pyglow/models/dl_models/iri12/ccir14.asc',
                'src/pyglow/models/dl_models/iri12/ccir15.asc',
                'src/pyglow/models/dl_models/iri12/ccir16.asc',
                'src/pyglow/models/dl_models/iri12/ccir17.asc',
                'src/pyglow/models/dl_models/iri12/ccir18.asc',
                'src/pyglow/models/dl_models/iri12/ccir19.asc',
                'src/pyglow/models/dl_models/iri12/ccir20.asc',
                'src/pyglow/models/dl_models/iri12/ccir21.asc',
                'src/pyglow/models/dl_models/iri12/ccir22.asc',
                'src/pyglow/models/dl_models/iri12/dgrf1945.dat',
                'src/pyglow/models/dl_models/iri12/dgrf1950.dat',
                'src/pyglow/models/dl_models/iri12/dgrf1955.dat',
                'src/pyglow/models/dl_models/iri12/dgrf1960.dat',
                'src/pyglow/models/dl_models/iri12/dgrf1965.dat',
                'src/pyglow/models/dl_models/iri12/dgrf1970.dat',
                'src/pyglow/models/dl_models/iri12/dgrf1975.dat',
                'src/pyglow/models/dl_models/iri12/dgrf1980.dat',
                'src/pyglow/models/dl_models/iri12/dgrf1985.dat',
                'src/pyglow/models/dl_models/iri12/dgrf1990.dat',
                'src/pyglow/models/dl_models/iri12/dgrf1995.dat',
                'src/pyglow/models/dl_models/iri12/dgrf2000.dat',
                'src/pyglow/models/dl_models/iri12/dgrf2005.dat',
                'src/pyglow/models/dl_models/iri12/ig_rz_IPS.dat',
                'src/pyglow/models/dl_models/iri12/ig_rz_SEC.dat',
                'src/pyglow/models/dl_models/iri12/ig_rz.dat',
                'src/pyglow/models/dl_models/iri12/igrf2010.dat',
                'src/pyglow/models/dl_models/iri12/igrf2010s.dat',
                'src/pyglow/models/dl_models/iri12/ursi11.asc',
                'src/pyglow/models/dl_models/iri12/ursi12.asc',
                'src/pyglow/models/dl_models/iri12/ursi13.asc',
                'src/pyglow/models/dl_models/iri12/ursi14.asc',
                'src/pyglow/models/dl_models/iri12/ursi15.asc',
                'src/pyglow/models/dl_models/iri12/ursi16.asc',
                'src/pyglow/models/dl_models/iri12/ursi17.asc',
                'src/pyglow/models/dl_models/iri12/ursi18.asc',
                'src/pyglow/models/dl_models/iri12/ursi19.asc',
                'src/pyglow/models/dl_models/iri12/ursi20.asc',
                'src/pyglow/models/dl_models/iri12/ursi21.asc',
                'src/pyglow/models/dl_models/iri12/ursi22.asc',
            ],
        ),
        (
            'pyglow/iri16_data/',
            [
                'src/pyglow/models/dl_models/iri16/apf107.dat',
                'src/pyglow/models/dl_models/iri16/ccir11.asc',
                'src/pyglow/models/dl_models/iri16/ccir12.asc',
                'src/pyglow/models/dl_models/iri16/ccir13.asc',
                'src/pyglow/models/dl_models/iri16/ccir14.asc',
                'src/pyglow/models/dl_models/iri16/ccir15.asc',
                'src/pyglow/models/dl_models/iri16/ccir16.asc',
                'src/pyglow/models/dl_models/iri16/ccir17.asc',
                'src/pyglow/models/dl_models/iri16/ccir18.asc',
                'src/pyglow/models/dl_models/iri16/ccir19.asc',
                'src/pyglow/models/dl_models/iri16/ccir20.asc',
                'src/pyglow/models/dl_models/iri16/ccir21.asc',
                'src/pyglow/models/dl_models/iri16/ccir22.asc',
                'src/pyglow/models/dl_models/iri16/dgrf1945.dat',
                'src/pyglow/models/dl_models/iri16/dgrf1950.dat',
                'src/pyglow/models/dl_models/iri16/dgrf1955.dat',
                'src/pyglow/models/dl_models/iri16/dgrf1960.dat',
                'src/pyglow/models/dl_models/iri16/dgrf1965.dat',
                'src/pyglow/models/dl_models/iri16/dgrf1970.dat',
                'src/pyglow/models/dl_models/iri16/dgrf1975.dat',
                'src/pyglow/models/dl_models/iri16/dgrf1980.dat',
                'src/pyglow/models/dl_models/iri16/dgrf1985.dat',
                'src/pyglow/models/dl_models/iri16/dgrf1990.dat',
                'src/pyglow/models/dl_models/iri16/dgrf1995.dat',
                'src/pyglow/models/dl_models/iri16/dgrf2000.dat',
                'src/pyglow/models/dl_models/iri16/dgrf2005.dat',
                'src/pyglow/models/dl_models/iri16/dgrf2010.dat',
                'src/pyglow/models/dl_models/iri16/ig_rz.dat',
                'src/pyglow/models/dl_models/iri16/igrf2015.dat',
                'src/pyglow/models/dl_models/iri16/igrf2015s.dat',
                'src/pyglow/models/dl_models/iri16/mcsat11.dat',
                'src/pyglow/models/dl_models/iri16/mcsat12.dat',
                'src/pyglow/models/dl_models/iri16/mcsat13.dat',
                'src/pyglow/models/dl_models/iri16/mcsat14.dat',
                'src/pyglow/models/dl_models/iri16/mcsat15.dat',
                'src/pyglow/models/dl_models/iri16/mcsat16.dat',
                'src/pyglow/models/dl_models/iri16/mcsat17.dat',
                'src/pyglow/models/dl_models/iri16/mcsat18.dat',
                'src/pyglow/models/dl_models/iri16/mcsat19.dat',
                'src/pyglow/models/dl_models/iri16/mcsat20.dat',
                'src/pyglow/models/dl_models/iri16/mcsat21.dat',
                'src/pyglow/models/dl_models/iri16/mcsat22.dat',
                'src/pyglow/models/dl_models/iri16/ursi11.asc',
                'src/pyglow/models/dl_models/iri16/ursi12.asc',
                'src/pyglow/models/dl_models/iri16/ursi13.asc',
                'src/pyglow/models/dl_models/iri16/ursi14.asc',
                'src/pyglow/models/dl_models/iri16/ursi15.asc',
                'src/pyglow/models/dl_models/iri16/ursi16.asc',
                'src/pyglow/models/dl_models/iri16/ursi17.asc',
                'src/pyglow/models/dl_models/iri16/ursi18.asc',
                'src/pyglow/models/dl_models/iri16/ursi19.asc',
                'src/pyglow/models/dl_models/iri16/ursi20.asc',
                'src/pyglow/models/dl_models/iri16/ursi21.asc',
                'src/pyglow/models/dl_models/iri16/ursi22.asc',
            ],
        ),
        (
            'pyglow/kpap/',
            [
                'src/pyglow/kpap/1932',
                'src/pyglow/kpap/1933',
                'src/pyglow/kpap/1934',
                'src/pyglow/kpap/1935',
                'src/pyglow/kpap/1936',
                'src/pyglow/kpap/1937',
                'src/pyglow/kpap/1938',
                'src/pyglow/kpap/1939',
                'src/pyglow/kpap/1940',
                'src/pyglow/kpap/1941',
                'src/pyglow/kpap/1942',
                'src/pyglow/kpap/1943',
                'src/pyglow/kpap/1944',
                'src/pyglow/kpap/1945',
                'src/pyglow/kpap/1946',
                'src/pyglow/kpap/1947',
                'src/pyglow/kpap/1948',
                'src/pyglow/kpap/1949',
                'src/pyglow/kpap/1950',
                'src/pyglow/kpap/1951',
                'src/pyglow/kpap/1952',
                'src/pyglow/kpap/1953',
                'src/pyglow/kpap/1954',
                'src/pyglow/kpap/1955',
                'src/pyglow/kpap/1956',
                'src/pyglow/kpap/1957',
                'src/pyglow/kpap/1958',
                'src/pyglow/kpap/1959',
                'src/pyglow/kpap/1960',
                'src/pyglow/kpap/1961',
                'src/pyglow/kpap/1962',
                'src/pyglow/kpap/1963',
                'src/pyglow/kpap/1964',
                'src/pyglow/kpap/1965',
                'src/pyglow/kpap/1966',
                'src/pyglow/kpap/1967',
                'src/pyglow/kpap/1968',
                'src/pyglow/kpap/1969',
                'src/pyglow/kpap/1970',
                'src/pyglow/kpap/1971',
                'src/pyglow/kpap/1972',
                'src/pyglow/kpap/1973',
                'src/pyglow/kpap/1974',
                'src/pyglow/kpap/1975',
                'src/pyglow/kpap/1976',
                'src/pyglow/kpap/1977',
                'src/pyglow/kpap/1978',
                'src/pyglow/kpap/1979',
                'src/pyglow/kpap/1980',
                'src/pyglow/kpap/1981',
                'src/pyglow/kpap/1982',
                'src/pyglow/kpap/1983',
                'src/pyglow/kpap/1984',
                'src/pyglow/kpap/1985',
                'src/pyglow/kpap/1986',
                'src/pyglow/kpap/1987',
                'src/pyglow/kpap/1988',
                'src/pyglow/kpap/1989',
                'src/pyglow/kpap/1990',
                'src/pyglow/kpap/1991',
                'src/pyglow/kpap/1992',
                'src/pyglow/kpap/1993',
                'src/pyglow/kpap/1994',
                'src/pyglow/kpap/1995',
                'src/pyglow/kpap/1996',
                'src/pyglow/kpap/1997',
                'src/pyglow/kpap/1998',
                'src/pyglow/kpap/1999',
                'src/pyglow/kpap/2000',
                'src/pyglow/kpap/2001',
                'src/pyglow/kpap/2002',
                'src/pyglow/kpap/2003',
                'src/pyglow/kpap/2004',
                'src/pyglow/kpap/2005',
                'src/pyglow/kpap/2006',
                'src/pyglow/kpap/2007',
                'src/pyglow/kpap/2008',
                'src/pyglow/kpap/2009',
                'src/pyglow/kpap/2010',
                'src/pyglow/kpap/2011',
                'src/pyglow/kpap/2012',
                'src/pyglow/kpap/2013',
            ],
        ),
        (
            'pyglow/dst/',
            [
                'src/pyglow/dst/1957_1969',
                'src/pyglow/dst/1970_1989',
                'src/pyglow/dst/1990_2004',
            ],
        ),
        (
            'pyglow/ae/',
            [
                'src/pyglow/ae/1975',
                'src/pyglow/ae/1978',
                'src/pyglow/ae/1979',
                'src/pyglow/ae/1980',
                'src/pyglow/ae/1981',
                'src/pyglow/ae/1982',
                'src/pyglow/ae/1983',
                'src/pyglow/ae/1984',
                'src/pyglow/ae/1985',
                'src/pyglow/ae/1986',
                'src/pyglow/ae/1987',
                'src/pyglow/ae/1988',
                'src/pyglow/ae/1989',
                'src/pyglow/ae/1990',
                'src/pyglow/ae/1991',
                'src/pyglow/ae/1992',
                'src/pyglow/ae/1993',
                'src/pyglow/ae/1994',
                'src/pyglow/ae/1995',
                'src/pyglow/ae/1996',
                'src/pyglow/ae/1997',
                'src/pyglow/ae/1998',
                'src/pyglow/ae/1999',
            ],
        ),
    ],
)

print("... All done!")
