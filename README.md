# lc-data-validator
Downloads and validates light curve data. Currently set to work with Fermi-LAT LCR data. See slides 17 and 18 at:
https://drive.google.com/drive/folders/1EQETWQcRyXy8gwP3rj7Kq6SL9WS3_-Mm

### To download all of the LCR data do:
`python lcr_import_data_by_source.py`
This will download and format the json files into ASCII .txt files

### For accessing the usage options do, e.g.:
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


### To use the data filtering:

`python lcr_data_validator.py [SOURCE_NAME] [FLUX_TYPE] [CADENCE] [PARAMETER] [COMPARATOR] [VALUE]`

These six parameters should be input in that specific order. Example:

`python lcr_data_validator.py "4FGL J2359.2-3134" free weekly photon_index_interval lt 2.5`

Some constraints:

`SOURCE_NAME`: if the source name has spaces, like above, put it in "double quotes". It also accepts underscored names like 4FGL_J2359.2-3134, which does not need quotes.

`FLUX_TYPE`: Accepted values: fixed, free, all. ("all" means this filter will be applied in both fixed and free tables)

`CADENCE`: Accepted values: daily, weekly, monthly, all. ("all" means this filter will be applied in daily, weekly and monthly tables)

`PARAMETER`: Accepted values: ts, flux, flux_error, photon_index, photon_index_interval (we will change these names in the future)

`COMPARATOR`: Accepted values: lessthan, lt, morethan, mt (lt serves to abbreviate lessthan and mt abbreviates morethan).
            If you put lessthan/lt it means "I want to KEEP all lines whose parameter is LESS than this value".
	    If you put morethan/mt it means "I want to KEEP all lines whose parameter is MORE than this value".

`VALUE`: Accepts any float value, like 1 or 2.001 or 3.4e-06 . Use "." not "," as decimal point.

Optional parameters:
1) `-d` or `--directory`: The code assumes that all your sources' data is in ./sources . If that's not the case, you can specify the directory of the sources via `-d` . For example:

`python lcr_data_validator.py "4FGL J2359.2-3134" free weekly photon_index_interval lt 2.5 -d /home/my_user/source_data`

(this means that you have a folder `/home/my_user/source_data/4FGL_J2359.2-3134/` that contains the data for this source)

2) `-r` or `--reset`: The code creates history files that help remember what has been already filtered. If you want to undo these, you can reset the history with this parameter.
The data that was removed will return to the original table and the graphs will be reset as well.

`python lcr_data_validator.py "4FGL J2359.2-3134" fixed daily -r`

(this will only rease the history for fixed_flux and for the daily table for this source) (rememver that flux_type must go before cadence)

`python lcr_data_validator.py "4FGL J2359.2-3134" -r`

(this will erase all histories for this source)

RESULTS:
- The data removed by these filters will be moved to their corresponding removed_[cadence]_table.txt with a reason for removal "removed_by_filter"
- The code creates and updates history files for the flux types and cadences that record the filters already appplied.
- The code creates or re-creates the two plots where the lines will be displayed according to the history files for this source.
- To reset a history, use the argument --reset (or just -r). This clears the relevant history, resets the relevants tables, and resets the plots.
