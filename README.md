# lc-data-validator
Downloads and validates light curve data. Currently set to work with Fermi-LAT LCR data.

### To download all of the LCR data do:
`python lcr_import_data_by_source.py`
This will download and format the json files into ASCII .txt files

### For accessing the usage options do:
```
python lcr_import_data_by_source.py -h
usage: lcr_import_data_by_source.py [-h] [-l] [-n] [-s] [-r] [sources [sources ...]]

Import source data from the LCR website

positional arguments:
  sources               Input source names (in quotes if the name has spaces) to import/update data for specific
                        sources

optional arguments:
  -h, --help            show this help message and exit
  -l, --list            Only list the available sources, importing no source data
  -n, --include-nearby  If using source names, add -n to also import/update all of their nearby sources
  -s, --save-list       Download source list as a .json after requesting it from the site
  -r, --read-list       Use a previously downloaded .json file as source list
```

