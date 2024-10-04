# Contact: Janeth Valverde <valverde@llr.in2p3.fr>

import re, os
from argparse import ArgumentParser
from auxiliary_classes import source_underscored, CADENCES, FLUX_TYPES
from lcr_plot_generator import outliers_plot

SEPARATOR = ', '
LESS_THAN_COMPARATORS = ['lessthan', 'lt']
MORE_THAN_COMPARATORS = ['morethan', 'mt']
COMPARATORS = LESS_THAN_COMPARATORS + MORE_THAN_COMPARATORS

PARAMETERS=['ts', 'flux', 'flux_error', 'photon_index', 'photon_index_interval']
REMOVAL_REASON = 'removed_by_filter'

# FORMATTING AND STRING METHODS

def flux_type_directory(directory, source, flux_type):
    return f'{directory}/{source}/{flux_type}_flux_tables'

def history_file(directory, source, flux_type, cadence):
    return f'{directory}/{source}/{flux_type}_flux_tables/filter_{cadence}_history.log'

def old_history_file(flux_type_directory, cadence):
    return f'{flux_type_directory}/filter_{cadence}_history.old.log'

def validated_table_file(directory, source, flux_type, cadence):
    return f'{directory}/{source}/{flux_type}_flux_tables/validated_{cadence}_table.txt'

def history_line(*values):
    return SEPARATOR.join((str(value) for value in values)) + '\n'

def removed_file(file):
    return re.sub('validated_', 'removed_', file)

def timestamp_from_line(line):
    return int(line.split(',')[0])

# LOGICAL AND FORMULAIC METHODS

def is_valid_line(value_list, parameter_list, filter_list):
    if len(value_list) != len(parameter_list):
        raise ValueError(f'ERROR: List of values does not match label count on: {", ".join(value_list)} vs {", ".join(parameter_list)}')
    parameter_values = { parameter_list[index]:float(value_list[index]) for index in range(len(parameter_list)) }
    for filter_condition in filter_list:
        parameter, comparator, value = (filter_condition[key] for key in ['parameter', 'comparator', 'value'])
        value = float(value)
        if comparator in LESS_THAN_COMPARATORS:
            if parameter_values[parameter] > value:
                return False
        elif comparator in MORE_THAN_COMPARATORS:
            if parameter_values[parameter] < value:
                return False
        else:
            raise ValueError(f'ERROR: Comparator outside of acceptable values: "{comparator}"')
    return True

def print_results_lines(results_lines): # unused
    (print(line) for line in results_lines)
    return

def silently_remove_file(file):
    try:
        os.remove(file)
    except FileNotFoundError:
        pass

# FUNCTIONAL METHODS

def parse_history(history_location):
    with open(history_location, 'r') as history:
        filter_list = []
        for line in history:
            parameter, comparator, value = line.strip().split(SEPARATOR, 2)
            filter_list.append({ 'parameter': parameter, 'comparator': comparator, 'value': value })
    
    return filter_list

def filter_data(history_location, table_location):
    filter_list = parse_history(history_location)
    headers, valid_lines, removable_lines = None, [], []

    with open(table_location, 'r+') as lcr_data:
        for index, line in enumerate(lcr_data):
            if index == 0:
                parameter_list = [value for value in line.strip().split(', ')]
                headers = line
            elif is_valid_line(line.strip().split(', '), parameter_list, filter_list):
                valid_lines.append(line)
            else:
                removable_lines.append(line.replace('\n', f', {REMOVAL_REASON}\n'))
        
        if len(removable_lines) > 0:
            lcr_data.seek(0)
            lcr_data.write(headers)
            for line in valid_lines:
                lcr_data.write(line)
            lcr_data.truncate()

    lines_moved = len(removable_lines)
    if lines_moved > 0:
        with open(removed_file(table_location), 'r+') as removed_lcr_data:
            for index, line in enumerate(removed_lcr_data):
                if index == 0:
                    headers = line
                else:
                    removable_lines.append(line)
        
            removable_lines.sort(key=timestamp_from_line)
                
            removed_lcr_data.seek(0)
            removed_lcr_data.write(headers)
            for line in removable_lines:
                removed_lcr_data.write(line)
            removed_lcr_data.truncate()

    return lines_moved

def clear_all_filters(table_location):
    headers, restorable_lines, non_restorable_lines = None, [], []

    with open(removed_file(table_location), 'r+') as removed_lcr_data:
        for line in removed_lcr_data:
            if re.search(f', {REMOVAL_REASON}', line) is not None:
                restorable_lines.append(line.replace(f', {REMOVAL_REASON}\n', '\n'))
            else:
                non_restorable_lines.append(line)
            
        if len(restorable_lines) > 0:
            removed_lcr_data.seek(0)
            for line in non_restorable_lines:
                removed_lcr_data.write(line)
            removed_lcr_data.truncate()

    lines_moved = len(restorable_lines)
    if lines_moved > 0:
        with open(table_location, 'r+') as lcr_data:
            for index, line in enumerate(lcr_data):
                if index == 0:
                    headers = line
                else:
                    restorable_lines.append(line)
            
            restorable_lines.sort(key=timestamp_from_line)
            
            lcr_data.seek(0)
            lcr_data.write(headers)
            for line in restorable_lines:
                lcr_data.write(line)
            lcr_data.truncate()
    
    return lines_moved

def redraw_validation_plots(directory, source):
    full_filter_list = {}
    for flux_type in FLUX_TYPES:
        full_filter_list[flux_type] = { parameter:[] for parameter in PARAMETERS }
        for cadence in CADENCES:
            history_location = history_file(directory, source, flux_type, cadence)
            if os.path.isfile(history_location):
                filter_list = parse_history(history_location)
                for filter_entry in filter_list:
                    full_filter_list[flux_type][filter_entry['parameter']].append({ 'cadence': cadence, 'value': filter_entry['value'] })
    print(full_filter_list)
    outliers_plot(source, parameter_filters=full_filter_list)


# MAIN FUNCTION

def parse_arguments():
    parser = ArgumentParser(description='Set lower and upper limits for parameters and filter data only within those limits')
    parser.add_argument('source', help='Input source name (in quotes if the name has spaces)')
    parser.add_argument('flux_types', nargs='?', choices=[*FLUX_TYPES, 'all'], default='all', help='flux_types: Choose to filter data for fixed and/or free flux tables')
    parser.add_argument('cadences', nargs='?', choices=[*CADENCES, 'all'], default='all', help='cadences: Choose to filter data for the daily, weekly, and/or monthly tables')

    parser.add_argument('parameter', nargs='?', choices=PARAMETERS, help='The parameter to use in the filtering')
    parser.add_argument('comparator', nargs='?', choices=COMPARATORS, help='Choose to keep values either less than or more than [value]')
    parser.add_argument('value', nargs='?', type=float, help='The value to compare the parameters to')

    parser.add_argument('-d', '--directory', default = './sources', help='The directory containing all the sources. Default: "./sources"')
    parser.add_argument('-r', '--reset', action='store_true', help='Delete existing filter history and start from a fresh one')
    parser.add_argument('-c', '--clear', action='store_true', help='Delete all history files, do not create new ones. Delete validation graphs.')

    args = parser.parse_args()
    args.source = source_underscored(args.source)
    args.flux_types = FLUX_TYPES if args.flux_types == 'all' else [args.flux_types]
    args.cadences = CADENCES if args.cadences == 'all' else [args.cadences]
    
    if args.comparator in LESS_THAN_COMPARATORS:
        args.comparator = LESS_THAN_COMPARATORS[0]
    elif args.comparator in MORE_THAN_COMPARATORS:
        args.comparator = MORE_THAN_COMPARATORS[0]
    
    args.new_filter = args.parameter is not None and args.comparator is not None and args.value is not None
    return args

if __name__ == '__main__':
    args = parse_arguments()

    if args.clear:
        for flux_type in FLUX_TYPES:
            for cadence in CADENCES:
                clear_all_filters(validated_table_file(args.directory, args.source, flux_type, cadence))
                silently_remove_file(history_file(args.directory, args.source, flux_type, cadence))
        silently_remove_file(f'{args.directory}/{args.source}/{args.source}_outliers_fixed_flux.png')
        silently_remove_file(f'{args.directory}/{args.source}/{args.source}_outliers_free_flux.png')
        quit()

    writing_mode = 'w' if args.reset else 'a'
    for flux_type in args.flux_types:
        for cadence in args.cadences:
            history_location = history_file(args.directory, args.source, flux_type, cadence)
            with open(history_location, writing_mode) as file:
                if args.new_filter:
                    file.write(history_line(args.parameter, args.comparator, args.value))

    if args.reset:
        print(f'LINES RESTORED FOR {args.source}')
        for flux_type in args.flux_types:
            results_line = f'{flux_type.upper()}_FLUX ||'
            for cadence in args.cadences:
                lines_moved = clear_all_filters(validated_table_file(args.directory, args.source, flux_type, cadence))
                results_line += f' {cadence.upper()}: {lines_moved} |'
            print(results_line)
            
    if args.new_filter or not args.reset:
        print(f'LINES REMOVED FOR {args.source}')
        for flux_type in args.flux_types:
            results_line = f'{flux_type.upper()}_FLUX >> '
            for cadence in args.cadences:
                history_location = history_file(args.directory, args.source, flux_type, cadence)
                lines_moved = filter_data(history_location, validated_table_file(args.directory, args.source, flux_type, cadence))
                results_line += f' {cadence.upper()}: {lines_moved} |'
            print(results_line)

    redraw_validation_plots(args.directory, args.source)
