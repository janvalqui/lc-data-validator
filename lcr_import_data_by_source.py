# Contact: Janeth Valverde <valverde@llr.in2p3.fr>

import requests
import json
import sys
import os
import time
from argparse import ArgumentParser
from auxiliary_classes import print_list_by_batches, angular_distance, source_spaced, source_underscored
from auxiliary_classes import FLUX_TYPES, CADENCES

URL = 'https://fermi.gsfc.nasa.gov/ssc/data/access/lat/LightCurveRepository/queryDB.php'
VARIABILITY_THRESHOLD = 21.67
SOURCE_FILE_NAME = './source_list.json'
FILTERED_SOURCE_FILE_NAME = './filtered_source_list.json'
LC_FILES_PATH = './sources'
PLACEHOLDER = 'NaN'

FLUX_TABLE_TITLES = ['ts', 'flux', 'flux_error', 'photon_index', 'photon_index_interval', 'fit_tolerance', 'fit_convergence']
                    # + ['dlogl', 'EG', 'GAL', 'bin_id']
UPPER_LIMIT_TABLE_TITLES = ['ts', 'flux_upper_limits', 'fit_tolerance', 'fit_convergence']

def request_get_with_timer(*args, **kwargs):
    start_time = time.time()
    result = requests.get(*args, **kwargs)
    end_time = time.time()
    print(f'({(end_time - start_time):.4g}s)', flush=True)
    return result

def fetch_light_curve(source_name, cadence = 'monthly', flux_type='fixed', ts_min = '1'):
    light_curve_params = {
        'typeOfRequest': 'lightCurveData',
        'cadence': cadence,
        'flux_type': 'photon',
        'index_type': flux_type,
        'ts_min': ts_min,
        'source_name': source_name,
    }
    print(f'Fetching {flux_type.upper()} {cadence.upper()} data for {source_name}', end=' ', flush=True)
    req = request_get_with_timer(URL, params=light_curve_params)
    return json.loads(req.text)

def format_light_curve_data(data, placeholder):
    timestamp_list = sorted(list(map(lambda entry: entry[0], data['ts'])))
    upper_limit_list = sorted(list(map(lambda entry: entry[0], data['flux_upper_limits'])))
    flux_list = [timestamp for timestamp in timestamp_list if timestamp not in upper_limit_list]

    flux_dict = { time:{ key:placeholder for key in FLUX_TABLE_TITLES } for time in flux_list }
    for key in FLUX_TABLE_TITLES:
        for entry in data[key]:
            if not isinstance(entry, list):
                break
            if entry[0] in flux_list:
                if len(entry) == 2:
                    flux_dict[entry[0]][key] = entry[1]
                if len(entry) == 3:
                    difference = abs((entry[2]-entry[1])/2)
                    flux_dict[entry[0]][key] = f'{difference:.3g}'

    upper_limit_dict = { time:{ key:placeholder for key in FLUX_TABLE_TITLES } for time in upper_limit_list }
    for key in UPPER_LIMIT_TABLE_TITLES:
        for entry in data[key]:
            if entry[0] in upper_limit_list:
                upper_limit_dict[entry[0]][key] = entry[1]
    
    flux_table = [['timestamp']]
    flux_table[0].extend(FLUX_TABLE_TITLES)
    for time in flux_list:
        flux_row = [time]
        flux_row.extend(flux_dict[time][key] for key in FLUX_TABLE_TITLES)
        flux_table.append(flux_row)

    upper_limit_table = [['timestamp']]
    upper_limit_table[0].extend(UPPER_LIMIT_TABLE_TITLES)
    for time in upper_limit_list:
        upper_limit_row = [time]
        upper_limit_row.extend(upper_limit_dict[time][key] for key in UPPER_LIMIT_TABLE_TITLES)
        upper_limit_table.append(upper_limit_row)  

    return flux_table, upper_limit_table
    
def fetch_source_list(threshold = None, save_list = False):
    print(f'Fetching the list of sources from {URL}', end=' ', flush=True)
    req = request_get_with_timer(URL, params={ 'typeOfRequest': 'SourceList', 'catalog': '4FGL' })
    source_list = json.loads(req.text)

    if save_list:
        print(f'Saving full source list containing {len(source_list)} sources to file {SOURCE_FILE_NAME}')
        with open(SOURCE_FILE_NAME, 'w') as file:
            file.write(json.dumps(source_list, indent=2))
    
    if threshold is None:
        return source_list

    filtered_list = [source for source in source_list if float(source['Variability_Index']) > threshold]

    if save_list:
        print(f'Saving (filtered) valid source list containing {len(filtered_list)} sources to file {FILTERED_SOURCE_FILE_NAME}')
        with open(FILTERED_SOURCE_FILE_NAME, 'w') as file:
            file.write(json.dumps(filtered_list, indent=2))

    return filtered_list

def read_source_list():
    print(f'Reading the list of sources from {FILTERED_SOURCE_FILE_NAME}')
    with open(FILTERED_SOURCE_FILE_NAME, 'r') as file:
        filtered_list = json.load(file)
    return filtered_list

def source_name_map(source_list):
    return [source['Source_Name'] for source in source_list]

def find_source_in_list(source_name, source_list):
    try:
        return next(source for source in source_list if source['Source_Name'] == source_name)
    except StopIteration:
        print(f'ERROR: Cound not find source "{source_name}" in list of {len(source_list)} sources.')
        sys.exit(1)

def is_neighbour(source_1, source_2):
    # First, verify they're not exactly the same source
    if source_1['Source_Name'] == source_2['Source_Name']:
        return False
    ascension_1, ascension_2 = float(source_1['RAJ2000']), float(source_2['RAJ2000'])
    declination_1, declination_2 = float(source_1['DEJ2000']), float(source_2['DEJ2000'])
    # Then, discard declinations too distant since it's much easier to filter
    if abs(declination_1 - declination_2) > 12:
        return False
    # Finally, evaluate angular distance and compare
    return angular_distance(ascension_1, declination_1, ascension_2, declination_2) <= 12

def select_nearby_sources(main_source_name, source_list=None):
    if source_list is None:
        source_list = read_source_list()
    main_source = find_source_in_list(main_source_name, source_list)
    return [source for source in source_list if is_neighbour(source, main_source)]

def import_data_for_single_source(source_name, counter):
    start_time = time.time()
    source_data_path = f"{LC_FILES_PATH}/{source_underscored(source_name)}"

    os.makedirs(f"{source_data_path}/fixed_flux_tables", exist_ok=True)
    os.makedirs(f"{source_data_path}/fixed_upper_limit_tables", exist_ok=True)
    os.makedirs(f"{source_data_path}/free_flux_tables", exist_ok=True)
    os.makedirs(f"{source_data_path}/free_upper_limit_tables", exist_ok=True)

    print(f"{counter}: Fetching light curve data for {source_name}")

    for flux_type in FLUX_TYPES:
        for cadence in CADENCES:
            json_values = fetch_light_curve(source_name, cadence = cadence, flux_type = flux_type, ts_min='1')
            flux_table, upper_limit_table = format_light_curve_data(json_values, placeholder=PLACEHOLDER)

            filename = f"{source_data_path}/{flux_type}_flux_tables/{cadence}_table.txt"
            with open(filename, 'w') as file:
                print(f"{counter}: Printing {flux_type.upper()} light curve {cadence.upper()} flux data to file {filename}")
                for index in range(len(flux_table)):
                    file.write(f"{', '.join(str(value) for value in flux_table[index])}\n")

            filename = f"{source_data_path}/{flux_type}_upper_limit_tables/{cadence}_table.txt"
            with open(filename, 'w') as file:
                print(f"{counter}: Printing {flux_type.upper()} light curve {cadence.upper()} flux upper limit data to file {filename}")
                for index in range(len(upper_limit_table)):
                    file.write(f"{', '.join(str(value) for value in upper_limit_table[index])}\n")
    end_time = time.time()
    print(f'{counter}: Finished printing data for {source_name} in {end_time - start_time:.3g}s\n')
    return

def get_arguments():
    parser = ArgumentParser(description='Import source data from the LCR website')
    parser.add_argument("sources", nargs='*', help='Input source names (in quotes if the name has spaces) to import/update data for specific sources')
    parser.add_argument("-l", "--list", default=False, action='store_true', help='Only list the available sources, importing no source data')
    parser.add_argument("-n", "--include-nearby", default=False, action='store_true', help='If using source names, add -n to also import/update all of their nearby sources')
    parser.add_argument("-s", "--save-list", default=False, action='store_true', help='Download source list as a .json after requesting it from the site')
    parser.add_argument("-r", "--read-list", default=False, action='store_true', help='Use a previously downloaded .json file as source list')

    return parser.parse_args()

if __name__ == "__main__":
    args = get_arguments()

    start_time = time.time()
    if args.read_list:
        valid_source_list = read_source_list()
    else:
        valid_source_list = fetch_source_list(threshold = VARIABILITY_THRESHOLD, save_list=args.save_list)
    
    valid_source_names = source_name_map(valid_source_list)
    print(f"Valid sources found: {len(valid_source_names)}")

    if args.list:
        print_list_by_batches(valid_source_names, batch_size=10)
        sys.exit(0)

    if len(args.sources) > 0:
        sources_to_import = [source_spaced(source) for source in args.sources if source_spaced(source) in valid_source_names]
        if args.include_nearby:
            nearby_sources_to_import = []
            for source_name in sources_to_import:
                print(f"Source {source_name} will be imported along with all its nearby sources.")
                nearby_source_list = select_nearby_sources(source_name, valid_source_list)
                nearby_source_names = source_name_map(nearby_source_list)
                print(f'The following {len(nearby_source_names)} sources are neighbouring sources to {source_name}:\n')
                print_list_by_batches(nearby_source_names, batch_size=10)

                # Remove eventually since this is only debug info
                print(f'Their coordinates are:')
                main_source = find_source_in_list(source_name, valid_source_list)
                print(f"Original >> NAME: {main_source['Source_Name']}, RAJ:{main_source['RAJ2000']}, DEC:{main_source['DEJ2000']}")
                for index, source in enumerate(nearby_source_list):
                    print(f"{index + 1} >> NAME: {source['Source_Name']}, RAJ:{source['RAJ2000']}, DEC:{source['DEJ2000']}")
                # End of debug info

                nearby_sources_to_import.extend(nearby_source_names)

            for source_name in nearby_sources_to_import:
                if source_name not in sources_to_import:
                    sources_to_import.append(source_name)

        invalid_sources = [source_spaced(source) for source in args.sources if source_spaced(source) not in valid_source_names]
        if len(invalid_sources) > 0:
            print(f"WARNING: {invalid_sources} not found among the list of valid sources. These sources will not be imported.")

    else:
        sources_to_import = valid_source_names[:]
    
    print(f"Importing a total of {len(sources_to_import)} sources:")
    print_list_by_batches(sources_to_import, batch_size=10)

    os.makedirs(LC_FILES_PATH, exist_ok=True)

    counter = 0
    time_counter = 0
    for source_name in sources_to_import:
        counter += 1
        start_time_for_source = time.time()
        import_data_for_single_source(source_name, counter)
        end_time_for_source = time.time()
        time_counter += end_time_for_source - start_time_for_source

    print(f'Finished printing files for {counter} sources')
    end_time = time.time()
    total_minutes, total_seconds = divmod(round(end_time - start_time), 60)
    total_hours, total_minutes = divmod(total_minutes, 60)
    print(f'Start time: {time.ctime(start_time)}\nEnd time: {time.ctime(end_time)}')
    average = time_counter/counter if counter > 0 else 0
    print(f'Total time: {total_hours:d}h {total_minutes:02d}min {total_seconds:02d}s. Average time per source: {average:.3g}s')
