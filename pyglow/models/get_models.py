
import os
import urllib2


# IRI 2012:
iri12 = {
    'folder' : 'iri12',
    'name' : 'iri12',
    #'url' : 'ftp://nssdcftp.gsfc.nasa.gov/models/ionospheric/iri/iri2012/00_iri2012.tar',
    'url' : 'http://spdf.gsfc.nasa.gov/pub/models/iri/iri2012/00_iri2012.tar',
    'filename' : '00_iri2012.tar',
    'tar' : True,
    'zip' : False,
    'tar_folder' : 'IRI12',
}

"""
The MSIS00 model does not seem to be available from any URL. As a
stopgap, we will now include the model as part of the pyglow code
base.
"""
# msis = {\
#         'folder'     : 'msis',\
#         'name'       : 'msis00',\
#         #'url'        : 'ftp://nssdcftp.gsfc.nasa.gov/models/atmospheric/msis/nrlmsise00/nrlmsise00_sub.for',\
#         'url'        : 'ftp://hanna.ccmc.gsfc.nasa.gov/pub/modelweb/atmospheric/msis/nrlmsise00/nrlmsise00_sub.for',\
#         'filename'   : 'nrlmsise00_sub.for',\
#         'tar'        : False,\
#         'zip'        : False,\
#         }

igrf11 = {
    'folder' : 'igrf11',
    'name' : 'igrf11',
    'url' : 'http://www.ngdc.noaa.gov/IAGA/vmod/igrf11.f',
    'filename' : 'igrf11.f',
    'tar' : False,
    'zip' : False,
}

igrf12 = {
    'folder' : 'igrf12',
    'name' : 'igrf12',
    'url' : 'http://www.ngdc.noaa.gov/IAGA/vmod/igrf12.f',
    'filename' : 'igrf12.f',
    'tar' : False,
    'zip' : False,
}

hwm07 = {
    'folder' : 'hwm07',
    'name' : 'hwm07',
    #'url' : 'http://nssdcftp.gsfc.nasa.gov/models/atmospheric/hwm07/HWM07_all_files.zip',
    'url' : 'http://remote2.ece.illinois.edu/~airglow/models/hwm07/HWM07_all_files.zip',
    'filename' : 'HWM07_all_file.zip',
    'tar' : False,
    'zip' : True,
    'zip_folder' : 'HWM07',
}

"""
Similar issue to MSIS --- the CCMC FTP server is no longer active
and an alternative source is unknown. Hosting the file in the
repository for now.
"""
# hwm93 = {\
#         'folder'     : 'hwm93',\
#         'name'       : 'hwm93',\
#         #'url'        : 'http://nssdcftp.gsfc.nasa.gov/models/atmospheric/hwm93/hwm93.txt',\
#         'url'        : 'ftp://hanna.ccmc.gsfc.nasa.gov/pub/modelweb/atmospheric/hwm93/hwm93.txt',\
#         'filename'   : 'hwm93.f',\
#         'tar'        : False,\
#         'zip'        : False,\
#         }

for model in [igrf11, igrf12, hwm07, iri12]:
    print "Downloading files for %s ..." % (model['name'])
    modelfile = urllib2.urlopen(model['url'])
    output = open(
        "./dl_models/{}/{}".format(model['folder'], model['filename']),
        'wb',
    )
    output.write(modelfile.read())
    output.close()

    # Tar file?
    if model['tar']:
        print "time to tar..."
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


# ---------------------------------------------------------
# IRI16 has multiple files. For now, download separately:
# ---------------------------------------------------------
model_folder = 'iri16'

# tar files:
# - - - - -
model_urls = [
    'http://irimodel.org/IRI-2016/00_iri2016.tar',
    'http://irimodel.org/COMMON_FILES/00_ccir-ursi.tar',
]

for model_url in model_urls:
    tar_file = model_url.split('/')[-1]
    model_file = urllib2.urlopen(model_url)
    output = open("./dl_models/{}/{}".format(model_folder, tar_file), 'wb')
    output.write(model_file.read())
    output.close()
    cmd = 'tar -xvf ./dl_models/{model_folder}/{tar_file} '\
            '-C ./dl_models/{model_folder}/'.format(
                model_folder=model_folder,
                tar_file=tar_file,
            )
    print(cmd)
    os.system(cmd)

indice_urls = [
    'http://irimodel.org/indices/apf107.dat',
    'http://irimodel.org/indices/ig_rz.dat',
]

# dat files (indices)
# - - - - -
for indice_url in indice_urls:
    dat_file = indice_url.split('/')[-1]
    model_file = urllib2.urlopen(indice_url)
    output = open(
        "./dl_models/{}/{}".format(model_folder, dat_file),
        'wb',
    )
    output.write(model_file.read())
    output.close()


