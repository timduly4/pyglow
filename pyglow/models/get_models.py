
import os
import urllib2

iri = {\
        'folder'     : 'iri',\
        'name'       : 'iri12',\
        #'url'        : 'ftp://nssdcftp.gsfc.nasa.gov/models/ionospheric/iri/iri2012/00_iri2012.tar', \
        'url'        : 'http://spdf.gsfc.nasa.gov/pub/models/iri/iri2012/00_iri2012.tar',\
        'filename'   : '00_iri2012.tar',\
        'tar'        : True,\
        'zip'        : False,\
        'tar_folder' : 'IRI12',\
        }

msis = {\
        'folder'     : 'msis',\
        'name'       : 'msis00',\
        #'url'        : 'ftp://nssdcftp.gsfc.nasa.gov/models/atmospheric/msis/nrlmsise00/nrlmsise00_sub.for',\
        'url'        : 'ftp://hanna.ccmc.gsfc.nasa.gov/pub/modelweb/atmospheric/msis/nrlmsise00/nrlmsise00_sub.for',\
        'filename'   : 'nrlmsise00_sub.for',\
        'tar'        : False,\
        'zip'        : False,\
        }

igrf = {\
        'folder'     : 'igrf',\
        'name'       : 'igrf11',\
        'url'        : 'http://www.ngdc.noaa.gov/IAGA/vmod/igrf11.f',\
        'filename'   : 'igrf11.f',\
        'tar'        : False,\
        'zip'        : False,\
        }

hwm07 = {\
        'folder'     : 'hwm07',\
        'name'       : 'hwm07',\
        #'url'        : 'http://nssdcftp.gsfc.nasa.gov/models/atmospheric/hwm07/HWM07_all_files.zip',\
        'url'        : 'http://remote2.ece.illinois.edu/~airglow/models/hwm07/HWM07_all_files.zip',\
        'filename'   : 'HWM07_all_file.zip',\
        'tar'        : False,\
        'zip'        : True,\
        'zip_folder' : 'HWM07',\
        }

hwm93 = {\
        'folder'     : 'hwm93',\
        'name'       : 'hwm93',\
        #'url'        : 'http://nssdcftp.gsfc.nasa.gov/models/atmospheric/hwm93/hwm93.txt',\
        'url'        : 'ftp://hanna.ccmc.gsfc.nasa.gov/pub/modelweb/atmospheric/hwm93/hwm93.txt',\
        'filename'   : 'hwm93.f',\
        'tar'        : False,\
        'zip'        : False,\
        }

for model in [msis, igrf, hwm07, hwm93, iri]:
    print "Downloading files for %s ..." % (model['name'])
    modelfile = urllib2.urlopen(model['url'])
    output = open("./dl_models/%s/%s" % (model['folder'], model['filename']), 'wb')
    output.write(modelfile.read())
    output.close()
    if model['tar']:
        print "time to tar..."
        cmd = 'tar -xvf ./dl_models/%s/%s -C ./dl_models/%s/' % (model['folder'], model['filename'], model['folder'])
        print cmd
        os.system(cmd)
    if model['zip']:
        # unzip:
        cmd = 'unzip ./dl_models/%s/%s' % (model['folder'],model['filename'])
        print cmd
        os.system(cmd)
        # move contents:
        os.system('mv ./%s/* ./dl_models/%s' % (model['zip_folder'],model['folder']))
        # remove folder:
        os.system('rm -rf ./%s' % (model['zip_folder']))







