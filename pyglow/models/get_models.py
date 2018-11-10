from __future__ import print_function

from future import standard_library
standard_library.install_aliases()
import os
import urllib.request, urllib.error, urllib.parse
import shutil

GITHUB_BASE = "https://github.com/timduly4/pyglow/raw"

INDICE_URLS = [
    'http://irimodel.org/indices/apf107.dat',
    'http://irimodel.org/indices/ig_rz.dat',
]

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
        "{}/{}".format(GITHUB_BASE, 'url-backup/static/00_iri2012.tar'),
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
    'filename': 'igrf11.f',
    'tar': False,
    'zip': False,
}

# IGRF 12:
igrf12 = {
    'folder': 'igrf12',
    'name': 'igrf12',
    'url': 'http://www.ngdc.noaa.gov/IAGA/vmod/igrf12.f',
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

# Download each model:
for model in [igrf11, igrf12, hwm07, iri12]:

    # Parse model name:
    model_name = model['name']

    # Download the tar or zip file:
    print("Downloading files for %s ..." % (model['name']))
    modelfile = None
    try:
        modelfile = urllib.request.urlopen(model['url'])
    except urllib.error.HTTPError as e:

        # Try backups:
        print(
            "{} did not work, attempting backup urls...".format(model['url'])
        )
        if model['url_backup']:
            for url_backup in model['url_backup']:
                try:
                    print("Trying: {}".format(url_backup))
                    modelfile = urllib.request.urlopen(url_backup)
                except urllib.error.HTTPError as e:
                    pass
                if modelfile:
                    break
        else:
            raise e

    # Ensure that we downloaded the model:
    if not modelfile:
        raise ValueError("Unable to download: {}".format(model['name']))

    # Write data file:
    output = open(
        "./dl_models/{}/{}".format(model['folder'], model['filename']),
        'wb',
    )
    output.write(modelfile.read())
    output.close()

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
    'http://irimodel.org/IRI-2016/00_iri2016.tar',
    'http://irimodel.org/COMMON_FILES/00_ccir-ursi.tar',
]
for model_url in model_urls:

    # Download tar file:
    tar_file = model_url.split('/')[-1]
    model_file = urllib.request.urlopen(model_url)
    output = open("./dl_models/{}/{}".format(model_folder, tar_file), 'wb')
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

    # Local filename:
    dat_file = indice_url.split('/')[-1]

    # Make URL request:
    model_file = urllib.request.urlopen(indice_url)

    # Download filename:
    fn = "./dl_models/{}/{}".format(model_folder, dat_file)

    # Write file to disk:
    print("Downloading: {} to {}".format(indice_url, fn))
    output = open(fn, 'wb')
    output.write(model_file.read())
    output.close()
