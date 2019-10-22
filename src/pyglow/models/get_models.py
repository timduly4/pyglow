from __future__ import print_function

from future import standard_library
standard_library.install_aliases()
import os  # noqa E402
import urllib.request  # noqa E402
import urllib.error  # noqa E402
import urllib.parse  # noqa E402
import shutil  # noqa E402

GITHUB_BASE = "https://github.com/timduly4/pyglow/raw"

INDICE_URLS = [
    'http://irimodel.org/indices/apf107.dat',
    'http://irimodel.org/indices/ig_rz.dat',
]

FN_MSIS = 'nrlmsise00_sub.for'
FN_HWM93 = 'hwm93.f'

# Timeout to download a file:
TIMEOUT = 5  # [sec]

# Note on  MSIS00: the model does not seem to be available from any URL.
# We will now include the model as part of the pyglow code base.

# Note on HWM93: the CCMC FTP server is no longer active and an alternative
# source is unknown. Hosting the file in the repository for now.


# IRI 2012:
iri12 = {
    'folder': 'iri12',
    'name': 'iri12',
    'url': 'http://spdf.gsfc.nasa.gov/pub/models/iri/iri2012/00_iri2012.tar',
    'url_backup': [
        "{}/{}".format(GITHUB_BASE, 'master/static/00_iri2012.tar'),
    ],
    'filename': '00_iri2012.tar',
    'tar': True,
    'zip': False,
    'tar_folder': 'IRI12',
}

# IGRF 11:
igrf11 = {
    'folder': 'igrf11',
    'name': 'igrf11',
    'url': 'http://www.ngdc.noaa.gov/IAGA/vmod/igrf11.f',
    'url_backup': [
        "{}/{}".format(GITHUB_BASE, 'master/static/igrf11.f'),
        "{}/{}".format(GITHUB_BASE, 'igrf-backup/static/igrf11.f'),
    ],
    'filename': 'igrf11.f',
    'tar': False,
    'zip': False,
}

# IGRF 12:
igrf12 = {
    'folder': 'igrf12',
    'name': 'igrf12',
    'url': 'http://www.ngdc.noaa.gov/IAGA/vmod/igrf12.f',
    'url_backup': [
        "{}/{}".format(GITHUB_BASE, 'master/static/igrf12.f'),
        "{}/{}".format(GITHUB_BASE, 'igrf-backup/static/igrf12.f'),
    ],
    'filename': 'igrf12.f',
    'tar': False,
    'zip': False,
}

# HWM 07:
hwm07 = {
    'folder': 'hwm07',
    'name': 'hwm07',
    'url': "{}/{}".format(GITHUB_BASE, 'master/static/HWM07_all_files.zip'),
    'url_backup': [
        "http://remote2.ece.illinois.edu/~airglow/models/hwm07/"
        "HWM07_all_files.zip",
        "{}/{}".format(GITHUB_BASE, 'master/static/HWM07_all_files.zip'),
    ],
    'filename': 'HWM07_all_file.zip',
    'tar': False,
    'zip': True,
    'zip_folder': 'HWM07',
}

# HWM 93:
hwm93 = {
    'folder': 'hwm93',
    'name': 'hwm93',
    'url': "{}/{}/{}".format(GITHUB_BASE, 'master/static', FN_HWM93),
    'url_backup': [
        "{}/{}/{}".format(GITHUB_BASE, 'msis-fix/static', FN_HWM93),
    ],
    'tar': False,
    'zip': False,
    'filename': FN_HWM93,
}

# MSIS:
msis = {
    'folder': 'msis',
    'name': 'msis',
    'url': "{}/{}/{}".format(GITHUB_BASE, 'master/static', FN_MSIS),
    'url_backup': [
        "{}/{}/{}".format(GITHUB_BASE, 'msis-fix/static', FN_MSIS),
    ],
    'tar': False,
    'zip': False,
    'filename': FN_MSIS,
}

# files need to download or not:
is_downloaded = False

# Download each model:
for model in [igrf11, igrf12, hwm07, hwm93, iri12, msis]:

    is_downloaded = False

    # local file length:
    local_length = 0
    # remote file length:
    remote_length = 0

    # Parse model name:
    model_name = model['name']

    # data file name:
    fn = "./dl_models/{}/{}".format(
        model['folder'],
        model['filename'],
    )

    # avoid repeating to download:
    if os.path.isfile(fn):
        local_length = os.stat(fn).st_size

    # Download the tar or zip file:
    print(
        "Downloading files for {} at {}".format(
            model['name'],
            model['url'],
        )
    )
    modelfile = None
    modelfile_error = True
    try:
        modelfile = urllib.request.urlopen(model['url'], timeout=TIMEOUT)

        # Make sure we have a legit file (can be spoofed due to govt shutdown):
        if modelfile.getheader('Content-Type') == 'text/html':
            raise ValueError("Received an HTML file")

        modelfile_error = False

        # avoid repeating to download:
        remote_length = int(modelfile.getheader('Content-Length'))
        if local_length == remote_length:
            print("{} was already downloaded.".format(fn))
            is_downloaded = True

    except (urllib.error.HTTPError, urllib.error.URLError, ValueError) as e:

        # Try backups:
        print(
            "{} did not work, attempting backup urls...".format(model['url'])
        )
        if model['url_backup']:
            for url_backup in model['url_backup']:
                try:
                    print("Trying: {}".format(url_backup))
                    modelfile = urllib.request.urlopen(
                        url_backup,
                        timeout=TIMEOUT,
                    )

                    modelfile_error = False

                    # avoid repeating to download:
                    remote_length = int(modelfile.getheader('Content-Length'))

                except urllib.error.HTTPError as e:
                    pass
                if not modelfile_error:
                    break
            if local_length == remote_length:
                print("{} was already downloaded.".format(fn))
                is_downloaded = True
        else:
            raise e

    # Ensure that we downloaded the model:
    if modelfile_error:
        raise ValueError("Unable to download: {}".format(model['name']))
    
    if not is_downloaded:
        # Write data file:
        output = open(fn, 'wb')
        output.write(modelfile.read())
        output.close()
        print("Wrote: {}".format(fn))

    # Tar file?
    if model['tar']:
        print("time to tar...")
        cmd = 'tar -xvf ./dl_models/{folder}/{filename} '\
              '-C ./dl_models/{folder}/'.format(
                folder=model['folder'],
                filename=model['filename'],
              )
        print(cmd)
        os.system(cmd)

    # Zip file?
    if model['zip']:
        # unzip:
        cmd = 'unzip ./dl_models/{folder}/{filename}'.format(
            folder=model['folder'],
            filename=model['filename'],
        )
        print(cmd)
        os.system(cmd)

        # Move contents:
        os.system(
            'mv ./{zip_folder}/* ./dl_models/{folder}'.format(
                zip_folder=model['zip_folder'],
                folder=model['folder'],
            )
        )

        # Remove folder:
        os.system('rm -rf ./{}'.format(model['zip_folder']))

# IRI16, multiple files
model_folder = 'iri16'
model_urls = [
    'http://irimodel.org/IRI-2016/00_iri.tar',
    'http://irimodel.org/COMMON_FILES/00_ccir-ursi.tar',
]
for model_url in model_urls:

    is_downloaded = False

    # local file length:
    local_length = 0
    # remote file length:
    remote_length = 0

    tar_file = model_url.split('/')[-1]
    # data file name:
    fn = "./dl_models/{}/{}".format(model_folder, tar_file)

    # avoid repeating to download:
    if os.path.isfile(fn):
        local_length = os.stat(fn).st_size
    
    # Download the tar or zip file:
    print(
        "Downloading files for {} at {}".format(
            tar_file,
            model_url,
        )
    )

    # Open url for downloading tar file:
    model_file = urllib.request.urlopen(model_url)

    # avoid repeating to download:
    remote_length = int(model_file.getheader('Content-Length'))
    if local_length == remote_length:
        print("{} was already downloaded.".format(fn))
        is_downloaded = True

    if not is_downloaded:
        # Download tar file:
        output = open(fn, 'wb')
        output.write(model_file.read())
        output.close()

    # Untar:
    cmd = 'tar -xvf ./dl_models/{model_folder}/{tar_file} '\
          '-C ./dl_models/{model_folder}/'.format(
            model_folder=model_folder,
            tar_file=tar_file,
          )
    print(cmd)
    os.system(cmd)

# Download indice files:
for indice_url in INDICE_URLS:

    is_downloaded = False

    # local file length:
    local_length = 0
    # remote file length:
    remote_length = 0

    # Local filename:
    dat_file = indice_url.split('/')[-1]

    # Download filename:
    fn = "./dl_models/{}/{}".format(model_folder, dat_file)
        
    # avoid repeating to download:
    if os.path.isfile(fn):
        local_length = os.stat(fn).st_size

    # Make URL request:
    model_file = urllib.request.urlopen(indice_url)

    # avoid repeating to download:
    remote_length = int(model_file.getheader('Content-Length'))
    if local_length == remote_length:
        print("{} was already downloaded.".format(fn))
        is_downloaded = True
    
    if not is_downloaded:
        # Write file to disk:
        print("Downloading: {} to {}".format(indice_url, fn))
        output = open(fn, 'wb')
        output.write(model_file.read())
        output.close()
