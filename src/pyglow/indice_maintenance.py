import contextlib
import dateutil.parser
from datetime import timedelta
import os
import urllib.request
import shutil

from .constants import DIR_FILE, VERSION
from .geophysical_indices import Indice


def _update_kpap():
    '''
    Update the Kp and Ap indices used in pyglow. The files will be downloaded
    from GFZ to your pyglow installation directory.
    '''

    # Pyglow directory:
    pyglow_dir = os.path.join(DIR_FILE, "kpap/")
    
    src = "https://kp.gfz-potsdam.de/app/files/Kp_ap_Ap_SN_F107_since_1932.txt"
    dest = pyglow_dir + "/Kp_ap_Ap_SN_F107_since_1932.txt"

    urllib.request.urlretrieve(src, dest)

    return


def download_dst(year, month, des):
    """
    Helper function to earch for the appropriate location
    and download the DST index file from WDC Kyoto for the
    given month and year. Save it to the specified file "des".
    Return True if successful, False if not.
    """
    # There are three possible sources of data. Search for
    # them in the following order:
    # 1) Final
    # 2) Provisional
    # 3) Realtime
    year_month = '%i%02i' % (year, month)
    wgdc_fn = 'dst%s%02i.for.request' % (str(year)[2:], month)
    src_final = 'http://wdc.kugi.kyoto-u.ac.jp/dst_final/%s/%s' % \
        (year_month, wgdc_fn)
    src_provisional = \
        'http://wdc.kugi.kyoto-u.ac.jp/dst_provisional/%s/%s' % \
        (year_month, wgdc_fn)
    src_realtime = 'http://wdc.kugi.kyoto-u.ac.jp/dst_realtime/%s/%s' % \
        (year_month, wgdc_fn)

    success = False
    for src in [src_final, src_provisional, src_realtime]:
        try:
            with contextlib.closing(urllib.request.urlopen(src)) as r:
                contents = r.read().decode('utf8')
                # If that succeeded, then the file exists
                print(
                    "\nDownloading\n{src}\nto\n{des}".format(
                        src=src,
                        des=des
                    )
                )
                with open(des, 'w') as f:
                    f.write(contents)
                success = True
                break
        except urllib.error.HTTPError:
            pass

    return success


def _update_dst(year):
    """
    Update the Dst index files used in pyglow.
    The files will be downloaded from WDC Kyoto
    to your pyglow installation directory.

    _update_dst(years=None)

    :param years: (optional) a list of years to download.
            If this input is not provided, the full
            range of years starting from 2005 to the
            current year will be downloaded. Pre-2005
            files are shipped with pyglow.
    """

    # Pyglow directory:
    pyglow_dir = os.path.join(DIR_FILE, "dst/")

    # Loop through each month:
    for month in range(1, 13):
        des = '%s%i%02i' % (pyglow_dir, year, month)
        download_dst(year, month, des)

    return


def download_ae(year, month, des):
    '''
    Helper function to earch for the appropriate location
    and download the AE index file from WDC Kyoto for the
    given month and year. Save it to the specified file "des".
    Return True if successful, False if not.
    '''
    # There are three possible sources of data. Search for
    # them in the following order:
    # 1) Final
    # 2) Provisional
    # 3) Realtime
    year_month = '%i%02i' % (year, month)
    wgdc_fn = 'ae%s%02i.for.request' % (str(year)[2:], month)
    src_provisional = \
        'http://wdc.kugi.kyoto-u.ac.jp/ae_provisional/%s/%s' % \
        (year_month, wgdc_fn)
    src_realtime = 'http://wdc.kugi.kyoto-u.ac.jp/ae_realtime/%s/%s' % \
        (year_month, wgdc_fn)

    success = False
    for src in [src_provisional, src_realtime]:
        try:
            with contextlib.closing(urllib.request.urlopen(src)) as r:
                contents = r.readlines()
                # If that succeeded, then the file exists
                print(
                    "\nDownloading\n{src}\nto\n{des}".format(
                        src=src,
                        des=des,
                    )
                )
                with open(des, 'w') as f:
                    # this shrinks the filesize to hourly
                    for c in contents:
                        c = c.decode('utf8')
                        f.write(
                            "%s%s%s\n" % (c[12:18], c[19:21], c[394:400])
                        )
                success = True
                break
        except urllib.error.HTTPError:
            pass
    return success


def _update_ae(year):
    '''
    Update the AE index files used in pyglow. The files will be downloaded from
    WDC Kyoto to your pyglow installation directory.

    :param year: Year to download
    '''

    # Pyglow directory:
    pyglow_dir = os.path.join(DIR_FILE, "ae/")

    # Loop through each month:
    for month in range(1, 13):
        des = '%s%i%02i' % (pyglow_dir, year, month)
        download_ae(year, month, des)

    return


def update_indices(year0, year1=None):
    '''
    Update all geophysical indices (e.g., KP, DST, AE).

    :param years: (optional) a list of years to download.
            If this input is not provided, default
            values will be used.
    '''

    # Display version:
    txt = "pyglow v{}".format(VERSION)
    print(txt)
    print("-"*len(txt))

    if year1:
        years = range(year0, year1)
    else:
        years = [year0]

    _update_kpap()
    # Update each of the indices for each year:
    for year in years:
        _update_dst(year)
        _update_ae(year)

    return


def check_stored_indices(date0, date1):
    """
    Helper function to determine which dates do not have indices

    :param date0: String start date
    :param date1: String end date
    """

    # Parse input dates:
    dn0 = dateutil.parser.parse(date0)
    dn1 = dateutil.parser.parse(date1)

    # Find date range:
    dns = [dn0 + timedelta(days=kk) for kk in range((dn1-dn0).days)]

    print("Checking: input date range:")
    print("  {}".format(dn0.strftime("%Y-%m-%d")))
    print("  to")
    print("  {}".format(dn1.strftime("%Y-%m-%d")))

    have_all = True
    for dn in dns:

        # Instantiate indice class:
        indice = Indice(dn)

        # Find the indices:
        indice.run()

        # Are all the indices NaN?
        if indice.all_nan():
            status = "--- FAIL ---"
            failed = True
            have_all = False
        else:
            failed = False
            status = "OK"

        # Report:
        if failed:
            print("{}: {}".format(dn.strftime("%Y-%m-%d"), status))

    # Report only if there were no issues:
    if have_all:
        print(
            ">> We have all of the geophysical indices files between these "
            "dates."
        )

    return
