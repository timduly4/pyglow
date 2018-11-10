
all: source
	cd ./dl_models/igrf11/; make compile

	cd ./dl_models/igrf12/; make compile

	cd ./dl_models/msis/; make compile;

	cd ./dl_models/hwm93/; make compile;

	cd ./dl_models/iri12/; make compile;

	cd ./dl_models/iri16/; make compile;

	cd ./dl_models/hwm07/; make compile;

	cd ./dl_models/hwm14/; make compile;

download:
	python get_models.py

source: download
	cp ./f2py/igrf11/* ./dl_models/igrf11/
	cd ./dl_models/igrf11/; make source

	cp ./f2py/igrf12/* ./dl_models/igrf12/
	cd ./dl_models/igrf12/; make source

	cp -v ./f2py/msis/* ./dl_models/msis/
	cd ./dl_models/msis/; make source;

	cp ./f2py/hwm93/* ./dl_models/hwm93/
	cd ./dl_models/hwm93/; make source;

	cp ./f2py/iri12/* ./dl_models/iri12/
	cd ./dl_models/iri12/; make source;

	cp ./f2py/iri16/* ./dl_models/iri16/
	cd ./dl_models/iri16/; make source;

	cp ./f2py/hwm07/* ./dl_models/hwm07/
	cd ./dl_models/hwm07/; make source;

	cd ./dl_models/hwm14/; make source;

clean:
	rm -rf *~
	make -Cdl_models/hwm93 clean
	rm -rf ./dl_models/hwm07/*
	rm -rf ./dl_models/hwm14/*.pyf
	rm -rf ./dl_models/hwm14/out*
	rm -rf ./dl_models/hwm14/*.so
	rm -rf ./dl_models/igrf11/*
	rm -rf ./dl_models/igrf12/*
	rm -rf ./dl_models/iri12/*
	rm -rf ./dl_models/iri16/*
	make -Cdl_models/msis clean
	touch ./dl_models/hwm07/dummy.txt
	touch ./dl_models/hwm93/dummy.txt
	touch ./dl_models/igrf11/dummy.txt
	touch ./dl_models/igrf12/dummy.txt
	touch ./dl_models/iri12/dummy.txt
	touch ./dl_models/iri16/dummy.txt
	touch ./dl_models/msis/dummy.txt
