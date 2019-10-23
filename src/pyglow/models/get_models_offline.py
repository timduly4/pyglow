from __future__ import print_function

from future import standard_library
standard_library.install_aliases()
import os  # noqa E402
import urllib.request  # noqa E402
import urllib.error  # noqa E402
import urllib.parse  # noqa E402
import shutil  # noqa E402

URL_BASE = "http://127.0.0.1:8080"


INDICE_URLS = [
    '{}/apf107.dat'.format(URL_BASE),
    '{}/ig_rz.dat'.format(URL_BASE),
]

FN_MSIS = 'nrlmsise00_sub.for'
FN_HWM93 = 'hwm93.f'

# Timeout to download a file:
TIMEOUT = 5  # [sec]

# IRI 2012:
iri12 = {
    'folder': 'iri12',
    'name': 'iri12',
    'url': '{}/00_iri2012.tar'.format(URL_BASE),
    'filename': '00_iri2012.tar',
    'tar': True,
    'zip': False,
    'tar_folder': 'IRI12',
}

# IGRF 11:
igrf11 = {
    'folder': 'igrf11',
    'name': 'igrf11',
    'url': '{}/igrf11.f'.format(URL_BASE),
    'filename': 'igrf11.f',
    'tar': False,
    'zip': False,
}

# IGRF 12:
igrf12 = {
    'folder': 'igrf12',
    'name': 'igrf12',
    'url': '{}/igrf12.f'.format(URL_BASE),
    'filename': 'igrf12.f',
    'tar': False,
    'zip': False,
}

# HWM 07:
hwm07 = {
    'folder': 'hwm07',
    'name': 'hwm07',
    'url': "{}/HWM07_all_files.zip".format(URL_BASE),
    'filename': 'HWM07_all_file.zip',
    'tar': False,
    'zip': True,
    'zip_folder': 'HWM07',
}

# HWM 93:
hwm93 = {
    'folder': 'hwm93',
    'name': 'hwm93',
    'url': "{}/{}".format(URL_BASE, FN_HWM93),
    'tar': False,
    'zip': False,
    'filename': FN_HWM93,
}

# MSIS:
msis = {
    'folder': 'msis',
    'name': 'msis',
    'url': "{}/{}".format(URL_BASE, FN_MSIS),
    'tar': False,
    'zip': False,
    'filename': FN_MSIS,
}

# Download each model:
for model in [igrf11, igrf12, hwm07, hwm93, iri12, msis]:

    # Parse model name:
    model_name = model['name']

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

    except (urllib.error.HTTPError, urllib.error.URLError, ValueError) as e:

        raise e

    # Ensure that we downloaded the model:
    if modelfile_error:
        raise ValueError("Unable to download: {}".format(model['name']))

    # Write data file:
    fn = "./dl_models/{}/{}".format(
        model['folder'],
        model['filename'],
    )
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
    '{}/00_iri.tar'.format(URL_BASE),
    '{}/00_ccir-ursi.tar'.format(URL_BASE),
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
