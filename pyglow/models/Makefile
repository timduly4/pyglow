



all:
	python get_models.py

	cp ./f2py/igrf/* ./dl_models/igrf/
	cd ./dl_models/igrf/; make compile

	cp ./f2py/msis/* ./dl_models/msis/
	cd ./dl_models/msis/; make compile;

	cp ./f2py/hwm93/* ./dl_models/hwm93/
	cd ./dl_models/hwm93/; make compile;

	cp ./f2py/iri/* ./dl_models/iri/
	cd ./dl_models/iri/; make compile;

	cp ./f2py/hwm07/* ./dl_models/hwm07/
	cd ./dl_models/hwm07/; make compile;
	
	cd ./dl_models/hwm14/; make compile;

clean:
	rm -rf *~
	rm -rf ./dl_models/hwm07/*
	rm -rf ./dl_models/hwm93/*
	rm -rf ./dl_models/igrf/*
	rm -rf ./dl_models/iri/*
	rm -rf ./dl_models/msis/*
	touch ./dl_models/hwm07/dummy.txt
	touch ./dl_models/hwm93/dummy.txt
	touch ./dl_models/igrf/dummy.txt
	touch ./dl_models/iri/dummy.txt
	touch ./dl_models/msis/dummy.txt


