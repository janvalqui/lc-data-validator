# Contact: Janeth Valverde <valverde@llr.in2p3.fr>

import os, sys, time
from argparse import ArgumentParser
import matplotlib.pyplot as pyplot
from numpy import absolute, log10, isnan
from numpy import array, amax, amin, genfromtxt, linspace, polyfit, polyval
from numpy import count_nonzero as count, in1d as intersection, std as standard_deviation
import auxiliary_classes as aux
from auxiliary_classes import Limits, Logger, SourceData, FLUX_TYPES, CADENCES
from lcr_import_data_by_source import select_nearby_sources, read_source_list, find_source_in_list, FLUX_TABLE_TITLES

# CONSTANTS

LC_FILES_PATH = './sources'
SOURCE_NAMES = sorted(next(os.walk(LC_FILES_PATH))[1]) if os.path.exists(LC_FILES_PATH) else []

DEFAULT_ALLOWED_DEVIATIONS = { 'free': { 'daily': 4.0, 'weekly': 4.0, 'monthly': 5.0 }, 'fixed': { 'daily': 4.0, 'weekly': 4.0, 'monthly': 5.0 } }

NUMBER_OF_BINS = 14 # Rough estimate, check if a better one can be defined

DAY_BY_CADENCE = { 'daily': 3, 'weekly': 7, 'monthly': 30 }
COLOR_BY_CADENCE = { 'daily': 'C3', 'weekly': 'C1', 'monthly': 'C0' }
LLS_COLOR_BY_CADENCE = { 'daily': 'maroon', 'weekly': 'saddlebrown', 'monthly': 'darkblue' }
ZORDER_BY_CADENCE = { 'daily': 1, 'weekly': 2, 'monthly': 3 }

# AUXILIARY FUNCTIONS

def table_path(source=None, flux_type=None, cadence=None):
    return f'{LC_FILES_PATH}/{source}/{flux_type}_flux_tables/{cadence}_table.txt'

def validated_table_path(source=None, flux_type=None, cadence=None):
    return f'{LC_FILES_PATH}/{source}/{flux_type}_flux_tables/validated_{cadence}_table.txt'

def removed_table_path(source=None, flux_type=None, cadence=None):
    return f'{LC_FILES_PATH}/{source}/{flux_type}_flux_tables/removed_{cadence}_table.txt'

def correlation_plots_path(source):
    return f'{LC_FILES_PATH}/{source}/correlations'
    
def common_timestamp_filter(data_a, data_b):
    timestamps_a, timestamps_b = data_a[:,0], data_b[:,0]
    # These masks filter only common values on the two timestamp arrays so that fixed_timestamps[fixed_timestamp_mask] == free_timestamps[free_timestamp_mask]
    timestamp_a_mask = intersection(timestamps_a, timestamps_b, assume_unique=True)
    timestamp_b_mask = intersection(timestamps_b, timestamps_a, assume_unique=True)
    return data_a[timestamp_a_mask], data_b[timestamp_b_mask]


def simple_ts_and_flux_filter(data, flux_type, logger):
    logger.log(f'Filtering {len(data)} common timestamp rows')
    ts_data, flux_data, flux_error_data, spectral_index_error_data = data[:,1], data[:,2], data[:,3], data[:,5]

    ts_mask = ~isnan(ts_data)
    if any(isnan(ts_data)):
        logger.log(f'Reducing {len(data)} rows to {count(ts_mask)} rows due to {count(isnan(ts_data))} rows with null TS')

    ts_and_flux_mask = ts_mask & ~isnan(flux_data)
    if count(ts_mask) > count(ts_and_flux_mask):
        logger.log(f'Reducing {count(ts_mask)} rows to {count(ts_and_flux_mask)} rows due to {count(ts_mask) - count(ts_and_flux_mask)} rows with null Flux')

    ts_flux_and_flux_error_mask = ts_and_flux_mask & (flux_error_data != 0)
    if count(ts_and_flux_mask) > count(ts_flux_and_flux_error_mask):
        logger.log(f'WARNING, ignoring {count(ts_and_flux_mask) - count(ts_flux_and_flux_error_mask)} rows with Flux Error == 0')

    positive_ts_flux_and_flux_error_mask = ts_flux_and_flux_error_mask & (ts_data > 0)
    if count(ts_flux_and_flux_error_mask) > count(positive_ts_flux_and_flux_error_mask):
        logger.log(f'WARNING, ignoring {count(ts_flux_and_flux_error_mask) - count(positive_ts_flux_and_flux_error_mask)} rows with TS < 0')

    if flux_type == 'free':
        spectral_index_ts_and_flux_mask = positive_ts_flux_and_flux_error_mask & (spectral_index_error_data > 0.001)
        if count(positive_ts_flux_and_flux_error_mask) > count(spectral_index_ts_and_flux_mask):
            logger.log(f'WARNING, ignoring {count(positive_ts_flux_and_flux_error_mask) - count(spectral_index_ts_and_flux_mask)} rows with Spectral Index Error <= 0.001')  
        logger.log(f'Returning {count(spectral_index_ts_and_flux_mask)} rows of valid data')
        return SourceData(data[spectral_index_ts_and_flux_mask])
    else:
        logger.log(f'Returning {count(positive_ts_flux_and_flux_error_mask)} rows of valid data')
        return SourceData(data[positive_ts_flux_and_flux_error_mask])

def common_timestamp_ts_and_flux_filter(fixed_data, free_data, logger):
    fixed_timestamp_count, free_timestamp_count = len(fixed_data[:,0]), len(free_data[:,0])
    common_fixed_data, common_free_data = common_timestamp_filter(fixed_data, free_data)
    common_data_count = len(common_fixed_data[:,0])

    logger.log(f'Reduced {fixed_timestamp_count} fixed / {free_timestamp_count} free rows to {common_data_count} common timestamp rows')

    fixed_ts, free_ts = common_fixed_data[:,1], common_free_data[:,1]
    ts_mask = ~isnan(fixed_ts) & ~isnan(free_ts)
    if any(~ts_mask):
        logger.log(f'Further reducing {common_data_count} rows to {count(ts_mask)} rows due to {count(~ts_mask)} rows with null TS')

    fixed_flux, free_flux =  common_fixed_data[ts_mask][:,2], common_free_data[ts_mask][:,2]
    flux_mask = ~isnan(fixed_flux) & ~isnan(free_flux)
    if any(~flux_mask):
        logger.log(f'Further reducing {count(ts_mask)} rows to {count(flux_mask)} rows due to {count(~flux_mask)} rows with null Flux')
    logger.log(f'{count(flux_mask)} rows of both fixed and free data after final filtering.')

    return SourceData(common_fixed_data[ts_mask][flux_mask]), SourceData(common_free_data[ts_mask][flux_mask])

def least_squares_line_of_relative_flux(data):
    relative_flux_log = 2*log10(data.flux/data.flux_error)
    log_ts = log10(data.ts)
    coefficients = polyfit(log_ts, relative_flux_log, 1)
    return polyval(coefficients, log_ts)

def vertical_distances_to_relative_flux_lsl(data):
    relative_flux_log = 2*log10(data.flux/data.flux_error)
    return relative_flux_log - least_squares_line_of_relative_flux(data)

def calculate_outlier_mask(data, max_std_dev, logger):
    vertical_distances = vertical_distances_to_relative_flux_lsl(data)
    # Demonstrable(?) that the sum of the signed distances to the LLS line should be about zero, besides rounding errors
    # average = sum(vertical_distances)/vertical_distances.size
    std_dev = standard_deviation(vertical_distances)
    outlier_mask = absolute(vertical_distances) > std_dev * max_std_dev
    maximum_sigmas = amax(absolute(vertical_distances)/std_dev)
    maximum_allowed_sigmas =  amax(absolute(vertical_distances[~outlier_mask])/std_dev)

    logger.log('Calculating outliers')
    logger.log(f'MAX={amax(vertical_distances)}, MIN={amin(vertical_distances)}, STD={std_dev}, '
               f'MAXSIG={round(maximum_sigmas,2)}, ALLSIG={round(maximum_allowed_sigmas,2)}')
    logger.log(f'Found {count(outlier_mask)} outliers out of {vertical_distances.size} over {max_std_dev} std deviations away')
    if count(outlier_mask > 0):
        logger.log(f'{vertical_distances[outlier_mask]}')
    logger.log()
    return outlier_mask

# PLOT SETUP AND MAIN PLOT FUNCTIONS

### FLUX FLUX PLOT

def flux_flux_plot_spectral_setup():
    pyplot.ylabel('Flux - spectral index free')
    pyplot.xlabel('Flux - spectral index fixed')
    pyplot.yscale('log')
    pyplot.xscale('log')
    return

def flux_flux_plot(source, suffix):
    Logger(source).log('***CREATING FLUXFLUX PLOT***\n')
    flux_flux_plot_spectral_setup()

    for cadence in CADENCES:
        logger = Logger(source=source, cadence=cadence)

        fixed_table_path, free_table_path = (table_path(source, flux_type, cadence) for flux_type in FLUX_TYPES)
        logger.log(f'Checking files: {fixed_table_path}, {free_table_path}')
        with open(fixed_table_path, 'r')as fixed_table_data, open(free_table_path, 'r') as free_table_data:
            fixed_light_curve_temp = genfromtxt(fixed_table_data, skip_header=1, delimiter=',')
            free_light_curve_temp = genfromtxt(free_table_data, skip_header=1, delimiter=',')
            if len(fixed_light_curve_temp) == 0 or len(free_light_curve_temp) == 0:
                logger.log('Skipping plotting data for empty table.')
                continue
            fixed_data, free_data = common_timestamp_ts_and_flux_filter(fixed_light_curve_temp, free_light_curve_temp, logger)

            fixed_q_mask = fixed_data.quality != 0
            logger.log(f'Number of fixed bins with quality==nan are {count(isnan(fixed_data.quality))}')
            logger.log(f'Number of fixed bins with quality!=0 are {count(fixed_q_mask)}')
            # Probably a better way to represent the values below is to do them line by line
            logger.log(f'MJD: {fixed_data.timestamp[fixed_q_mask]}.\nTS: {fixed_data.ts[fixed_q_mask]}.\n'
                    f'Flux: {fixed_data.flux[fixed_q_mask]}.\nFluxErr: {fixed_data.flux_error[fixed_q_mask]}')

            free_q_mask = free_data.quality != 0
            logger.log(f'Number of free bins with quality==nan are {count(isnan(free_data.quality))}')
            logger.log(f'Number of free bins with quality!=0 are {count(free_q_mask)}')
            # Probably a better way to represent the values below is to do them line by line
            logger.log(f'MJD: {fixed_data.timestamp[free_q_mask]}.\nTS: {free_data.ts[free_q_mask]}.\n'
                    f'Flux: {free_data.flux[free_q_mask]}.\nFluxErr: {free_data.flux_error[free_q_mask]}')
            logger.log()

            common_errorbar_args = { 'color': COLOR_BY_CADENCE[cadence], 'zorder': ZORDER_BY_CADENCE[cadence], 'lw': 0.2 }

            pyplot.errorbar(
                fixed_data.flux, free_data.flux,
                xerr=fixed_data.flux_error, yerr=free_data.flux_error,
                label=f'LCR, {cadence}', fmt='.',
                **common_errorbar_args,
            )
            if any(fixed_q_mask):
                pyplot.errorbar(
                    fixed_data.flux[fixed_q_mask], free_data.flux[fixed_q_mask],
                    xerr=fixed_data.flux_error[fixed_q_mask], yerr=free_data.flux_error[fixed_q_mask],
                    label=f'LCR, {cadence}, $q\\neq 0$', fmt='o', fillstyle='none',
                    **common_errorbar_args,
                )
            if any(free_q_mask):
                pyplot.errorbar(
                    fixed_data.flux[free_q_mask], free_data.flux[free_q_mask],
                    xerr=fixed_data.flux_error[free_q_mask], yerr=free_data.flux_error[free_q_mask],
                    fmt='o', fillstyle='none',
                    **common_errorbar_args,
                )
    
    pyplot.legend(title=aux.source_spaced(source))
    file_name = f'{LC_FILES_PATH}/{source}/{source}_flux_flux{suffix}.png'
    pyplot.gcf().savefig(file_name, dpi=300)
    Logger(source=source).log(f'Saved file {file_name}\n')
    pyplot.clf()
    return

### OUTLIERS PLOT

def outliers_plot_setup(flux_type):
    rows, height = (4, 12) if flux_type == 'free' else (3, 9)

    _fig, axs = pyplot.subplots(nrows=rows, ncols=2, figsize=(11, height))
    pyplot.subplots_adjust(wspace=.25, hspace=.3)
    axs[0,0].set(ylabel='$(\\frac{F}{\sigma_F})^2$')
    axs[0,0].set(xlabel='$TS$')
    axs[0,0].grid(c='0.5',linewidth=0.2)
    axs[0,0].set_xscale('log')
    axs[0,0].set_yscale('log')

    axs[0,1].set(ylabel='occurrence / bin')
    axs[0,1].set(xlabel='Distance to MinSq line')
    axs[0,1].grid(c='0.5',linewidth=0.2)
    axs[0,1].set_xscale('linear') # No need for symlog since we're excluding every non-positive TS
    axs[0,1].set_yscale('log')
    
    axs[1,0].set(ylabel='occurrence / bin')
    axs[1,0].set(xlabel='$TS$')
    axs[1,0].grid(c='0.5',linewidth=0.2)
    axs[1,0].set_xscale('log') # No need for symlog since we're excluding every non-positive TS
    axs[1,0].set_yscale('log')

    axs[1,1].axis('off')

    axs[2,0].set(ylabel='occurrence / bin')
    axs[2,0].set(xlabel='$\log_{10}(\sigma_{F})$')
    axs[2,0].grid(c='0.5',linewidth=0.2)
    axs[2,0].set_yscale('log')

    axs[2,1].set(ylabel='occurrence / bin')
    axs[2,1].set(xlabel='$\log_{10}(F)$')
    axs[2,1].grid(c='0.5',linewidth=0.2)
    axs[2,1].set_yscale('log')

    if flux_type == 'free':
        axs[3,0].set(ylabel='occurrence / bin')
        axs[3,0].set(xlabel='$\sigma_{\\rm{Spectral\;index}}$')
        axs[3,0].grid(c='0.5',linewidth=0.2)
        axs[3,0].set_yscale('log')

        axs[3,1].set(ylabel='occurrence / bin')
        axs[3,1].set(xlabel='$\\rm{Spectral\;index}$')
        axs[3,1].grid(c='0.5',linewidth=0.2)
        axs[3,1].set_yscale('log')
    return axs

def outliers_plot(source, allowed_deviations=DEFAULT_ALLOWED_DEVIATIONS, suffix='', parameter_filters=None):
    Logger(source).log('***CREATING OUTLIERS PLOT***\n')
    Logger(source).log(f'Max standard deviations allowed for valid data is {allowed_deviations}\n')

    for flux_type in FLUX_TYPES:
        axs = outliers_plot_setup(flux_type)
        pyplot.suptitle(f'{aux.source_spaced(source)} {flux_type}_flux')

        ts_limits, flux_limits, flux_error_limits = Limits(), Limits(), Limits()
        for cadence in CADENCES:
            logger = Logger(source, flux_type, cadence)
            logger.log(f"Opening file {table_path(source, flux_type, cadence)}")
            light_curve_temp = genfromtxt(table_path(source, flux_type, cadence), skip_header=1, delimiter=',')
            if len(light_curve_temp) == 0:
                logger.log("Skipping plotting lines for empty table.")
                continue
            data = simple_ts_and_flux_filter(light_curve_temp, flux_type=flux_type, logger=logger)
            outlier_mask = calculate_outlier_mask(data, allowed_deviations[flux_type][cadence], logger)

            axs[0,0].plot(data.ts, 10**least_squares_line_of_relative_flux(data), LLS_COLOR_BY_CADENCE[cadence], label=f'MinSq, {cadence}', zorder=4, lw=1)
            axs[0,0].plot(
                data.ts[~outlier_mask],((data.flux/data.flux_error)**2)[~outlier_mask], f'.{COLOR_BY_CADENCE[cadence]}',
                label=f'LCR, {cadence}', zorder=1, lw=1
            )
            axs[0,0].plot(
                data.ts[outlier_mask], ((data.flux/data.flux_error)**2)[outlier_mask], LLS_COLOR_BY_CADENCE[cadence],
                label=f'Outliers, {cadence} ({count(outlier_mask)}), dev={allowed_deviations[flux_type][cadence]}$\sigma$', fillstyle='none', zorder=3, lw=1, marker='.', linestyle='none'
            )
            axs[0,0].plot(
                data.ts[data.quality != 0], ((data.flux/data.flux_error)**2)[data.quality != 0], f'o{COLOR_BY_CADENCE[cadence]}',
                label=f'LCR q>0, {cadence}', fillstyle='none', zorder=2, mew=1
            )

            common_hist_args = { 'color': COLOR_BY_CADENCE[cadence], 'histtype': 'step', 'zorder': 1 }
            any_quality_args = {**common_hist_args, 'linestyle': '--', 'lw': 1.2,  'label': f'LCR, {cadence}' }
            nonzero_quality_args = {**common_hist_args, 'lw': 1.8, 'label': f'LCR q>0, {cadence}' }
            outlier_quality_args = {**common_hist_args, 'color': LLS_COLOR_BY_CADENCE[cadence], 'lw': 1.5, 'linestyle': 'dashdot', 'label': f'Outliers, {cadence}' }

            vertical_distances = vertical_distances_to_relative_flux_lsl(data)
            sigma = standard_deviation(vertical_distances)
            vertical_distance_bins = linspace(amin(vertical_distances), amax(vertical_distances), num=NUMBER_OF_BINS)
            axs[0,1].hist(vertical_distances, bins=vertical_distance_bins, **{**any_quality_args, **{'label': f'{cadence}, $\sigma$={sigma:.3g}'}})
            axs[0,1].hist(vertical_distances[data.quality != 0], bins=vertical_distance_bins, **{**nonzero_quality_args, **{'label': f'{cadence}, q>0'}})

            ts_bins = linspace(amin(data.ts), amax(data.ts), num=NUMBER_OF_BINS)
            flux_error_bins = linspace(amin(log10(data.flux_error)), amax(log10(data.flux_error)), num=NUMBER_OF_BINS)
            flux_bins = linspace(amin(log10(data.flux)), amax(log10(data.flux)), num=NUMBER_OF_BINS)
            axs[1,0].hist(data.ts, bins=ts_bins, **any_quality_args)[2]  # No need for symlog since it's only positive TS. Should be log-binned?
            axs[1,0].hist(data.ts[data.quality != 0], bins=ts_bins, **nonzero_quality_args)[2]  # No need for symlog since it's only positive TS. Should be log-binned?
            axs[2,0].hist(log10(data.flux_error), bins=flux_error_bins, **any_quality_args)[2]
            axs[2,0].hist(log10(data.flux_error[data.quality != 0]), bins=flux_error_bins, **nonzero_quality_args)[2]
            axs[2,1].hist(log10(data.flux), bins=flux_bins, **any_quality_args)[2]
            axs[2,1].hist(log10(data.flux[data.quality != 0]), bins=flux_bins, **nonzero_quality_args)[2]

            if flux_type == 'free':
                spectral_index_error_bins = linspace(amin(data.spectral_index_error), amax(data.spectral_index_error), num=NUMBER_OF_BINS)
                spectral_index_bins = linspace(amin(data.spectral_index), amax(data.spectral_index), num=NUMBER_OF_BINS)
                axs[3,0].hist(data.spectral_index_error, bins=spectral_index_error_bins, **any_quality_args)[2]
                axs[3,0].hist(data.spectral_index_error[data.quality != 0], bins=spectral_index_error_bins, **nonzero_quality_args)[2]
                axs[3,1].hist(data.spectral_index, bins=spectral_index_bins, **any_quality_args)[2]
                axs[3,1].hist(data.spectral_index[data.quality != 0], bins=spectral_index_bins, **nonzero_quality_args)[2]

            outlier_data = SourceData(data.source_data[outlier_mask])

            axs[1,0].hist(outlier_data.ts, bins=ts_bins, **outlier_quality_args)[2]
            axs[2,0].hist(log10(outlier_data.flux_error), bins=flux_error_bins, **outlier_quality_args)[2]
            axs[2,1].hist(log10(outlier_data.flux), bins=flux_bins, **outlier_quality_args)[2]
            if flux_type == 'free':
                axs[3,0].hist(outlier_data.spectral_index_error, bins=spectral_index_error_bins, **outlier_quality_args)[2]
                axs[3,1].hist(outlier_data.spectral_index, bins=spectral_index_bins, **outlier_quality_args)[2]

            ts_limits.update(data.ts)
            flux_limits.update(log10(data.flux))
            flux_error_limits.update(log10(data.flux_error))

        if parameter_filters is not None and parameter_filters[flux_type] is not None:
            parameter_filters_by_name = parameter_filters[flux_type]
            for parameter_name in ['ts', 'flux', 'flux_error', 'photon_index', 'photon_index_interval']:
                filters = parameter_filters_by_name[parameter_name]
                for filter_entry in filters:
                    x_value = float(filter_entry['value'])
                    common_args = { 'color': COLOR_BY_CADENCE[filter_entry['cadence']], 'lw': 1, 'ls': '--' }
                    print(f'Adding value in {parameter_name}, {filter_entry["cadence"]} of value {filter_entry["value"]}')
                    if parameter_name == 'ts':
                        axs[0,0].axvline(x_value, **common_args)
                        axs[1,0].axvline(x_value, **common_args)
                    elif parameter_name == 'flux':
                        print(log10(x_value))
                        axs[2,1].axvline(x=log10(x_value), **common_args)
                    elif parameter_name == 'flux_error':
                        axs[2,0].axvline(x=log10(x_value), **common_args)
                    elif parameter_name == 'photon_index' and flux_type == 'free':
                        axs[3,1].axvline(x=x_value, **common_args)
                    elif parameter_name == 'photon_index_interval' and flux_type == 'free':
                        axs[3,0].axvline(x=x_value, **common_args)
        
        axs[1,0].set_xlim(ts_limits.as_tuple(with_margins=True, in_log_scale=True))
        axs[2,0].set_xlim(flux_error_limits.as_tuple(with_margins=True))
        axs[2,1].set_xlim(flux_limits.as_tuple(with_margins=True))
        axs[0,0].legend(loc='lower right', fontsize = 'xx-small')
        axs[0,1].legend(loc='upper left', fontsize = 'xx-small')
        axs[1,0].legend(loc='upper right', fontsize = 'xx-small')
        axs[2,0].legend(loc='upper right', fontsize = 'xx-small')
        axs[2,1].legend(loc='upper right', fontsize = 'xx-small')
        if flux_type == 'free':
            axs[3,0].legend(loc='upper right', fontsize = 'xx-small')
            axs[3,1].legend(loc='upper right', fontsize = 'xx-small')
    
        file_name = f'{LC_FILES_PATH}/{source}/{source}_outliers_{flux_type}_flux{suffix}.png'
        pyplot.gcf().savefig(file_name, dpi=300)
        Logger(source=source, flux_type=flux_type).log(f'Saved file {file_name}\n')
        # pyplot.show()
        pyplot.clf()
    return

def sun_and_moon_distance(source, placeholder):
    Logger(source).log('***ADDING SUN AND MOON DISTANCES TO THE TABLES***\n')
    
    filtered_source_list = read_source_list()
    source_content = find_source_in_list(aux.source_spaced(source), filtered_source_list)

    full_start_time = time.time()
    for flux_type in FLUX_TYPES:
        for cadence in CADENCES:
            logger = Logger(source, flux_type, cadence)
            existing_table_path = validated_table_path(source, flux_type, cadence)

            with open(existing_table_path, 'r+') as table_data:
                lines = table_data.readlines()
                table_data.seek(0)
                replace = False
                start_time = time.time()
                mid_time_before = start_time
                sun_counter, moon_counter = 0, 0
                sun_timer_counter, moon_timer_counter = 0, 0
                source_in_ecliptic_band = aux.is_source_within_ecliptic_band(source_content)

                for index, line in enumerate(lines):
                    if index == 0:
                        if line.strip().rsplit(', ', 1)[1] == 'distance_to_moon':
                            replace = True
                            logger.log(f'REPLACING sun and moon rows for table {existing_table_path}')
                            table_data.write(line)
                        else:
                            logger.log(f'Adding sun and moon rows for table {existing_table_path}')
                            table_data.write(f'{line.strip()}, distance_to_sun, distance_to_moon\n')
                    else:
                        if source_in_ecliptic_band:
                            timestamp = int(line.split(', ', 1)[0])
                            days = DAY_BY_CADENCE[cadence]
                            solar_timerange = linspace(timestamp - int(86400*days/2), timestamp + int(86400*days/2), days+1)
                            lunar_timerange = linspace(timestamp - int(86400*days/2), timestamp + int(86400*days/2), 13*days+1)

                            sun_timer_start = time.time()
                            min_distance_to_sun = aux.minimal_distance_to_body(source_content, 'Sun', timerange=solar_timerange)
                            sun_timer_counter += time.time() - sun_timer_start

                            moon_timer_start = time.time()
                            min_distance_to_moon = aux.minimal_distance_to_body(source_content, 'Moon', timerange=lunar_timerange)
                            moon_timer_counter += time.time() - moon_timer_start

                            mid_time_after=time.time()
                            
                            if (mid_time_after - mid_time_before > 3) or (index == len(lines) - 1):
                                logger.log(f'Lines completed so far: {index + 1}, Total time step: {(mid_time_after - mid_time_before):.4g}s, Total time elapsed: {(mid_time_after - start_time):.4g}s')
                                mid_time_before = mid_time_after

                            new_values = (f', {min_distance_to_sun:.3g}' if min_distance_to_sun < 12 else f', {placeholder}')
                            new_values += (f', {min_distance_to_moon:.3g}' if min_distance_to_moon < 12 else f', {placeholder}')

                            sun_counter += (1 if min_distance_to_sun < 12 else 0) 
                            moon_counter += (1 if min_distance_to_moon < 12 else 0)
                        else:
                            new_values = (f', {placeholder}, {placeholder}')

                        if replace: 
                            table_data.write(f'{line.strip().rsplit(", ", 2)[0]}{new_values}\n')
                        else:
                            table_data.write(f'{line.strip()}{new_values}\n')

                if not source_in_ecliptic_band:
                    logger.log(f'Source outside the 18° ecliptic band. AUTOFILLED values for sun and moon distances with {placeholder}')

                end_time=time.time()
                table_data.truncate()
                operation_count = len(lines) - 1 if len(lines) > 1 else 1
                logger.log(f'AVG sun calculation time: {(sun_timer_counter/operation_count):.3g}s. AVG moon calculation time: {(moon_timer_counter/operation_count):.3g}s.')
                logger.log(f'Total sun calculation time: {(sun_timer_counter):.3g}s. Total moon calculation time: {(moon_timer_counter):.3g}s. Total time: {(end_time - start_time):.4g}s')
                logger.log(f'Sun was close to the source {sun_counter} times. Moon was close the source {moon_counter} times.')
            
            logger.log(f'Finished updating table for {existing_table_path}')
    full_end_time = time.time()
    Logger(source=source).log(f'Total execution time for source: {int((full_end_time - full_start_time) // 60)}min {((full_end_time - full_start_time) % 60):.3g}s')
    return

# TODO: Refactor this method and the simple_ts_and_flux_filter and calculate_outlier_mask to work together
def validate_tables(source, placeholder, allowed_deviations=DEFAULT_ALLOWED_DEVIATIONS):
    Logger(source).log('***VALIDATING SOURCE DATA TABLES***\n')
    Logger(source).log(f'Max standard deviations allowed for valid data is {allowed_deviations}\n')
    for flux_type in FLUX_TYPES:
        for cadence in CADENCES:
            logger = Logger(source, flux_type, cadence)
            original_table_data = table_path(source, flux_type, cadence)
            logger.log(f'Validating table {original_table_data}')
            full_light_curve_data = genfromtxt(original_table_data, skip_header=1, delimiter=',')
            if len(full_light_curve_data) > 0:
                data = simple_ts_and_flux_filter(full_light_curve_data, flux_type=flux_type, logger=logger)
                outlier_mask = calculate_outlier_mask(data, allowed_deviations[flux_type][cadence], logger)
                valid_timestamps = data.source_data[(data.quality == 0) & ~outlier_mask][:, 0]
                outlier_timestamps = data.source_data[outlier_mask][:,0]
            else:
                logger.log("No data filtered for empty table.")

            validated_table_name = validated_table_path(source, flux_type, cadence)
            removed_data_table_name = removed_table_path(source, flux_type, cadence)
            logger.log(f'Saving validated table data in : {validated_table_name}, removed table data in {removed_data_table_name}')
            with open(original_table_data, 'r') as original, open(validated_table_name, 'w') as validate_data, open(removed_data_table_name, 'w') as removed_data:
                for index, line in enumerate(original):
                    stripped_line = line.rstrip('\n')
                    if index == 0:
                        validate_data.write(line)
                        removed_data.write(f'{stripped_line}, reasons_for_removal\n')
                    else:    
                        entries = array(stripped_line.split(', '), dtype=str)
                        timestamp, ts, flux, flux_error, _, spectral_index_error, _, quality = entries
                        if float(timestamp) in valid_timestamps:
                            validate_data.write(line)
                        elif timestamp is not None:
                            reasons_for_removal = []
                            if ts == placeholder:
                                reasons_for_removal.append('missing_ts')
                            elif float(ts) < 0:
                                reasons_for_removal.append('negative_ts')
                            if flux == placeholder:
                                reasons_for_removal.append('missing_flux')
                            if flux_error == placeholder:
                                reasons_for_removal.append('missing_flux_error')
                            elif float(flux_error) == 0:
                                reasons_for_removal.append('zero_flux_error')
                            if float(quality) != 0:
                                reasons_for_removal.append('non_zero_quality')
                            if float(timestamp) in outlier_timestamps:
                                reasons_for_removal.append('outlier_in_relative_flux')
                            if float(spectral_index_error) <= 0.001 and flux_type == 'free':
                                reasons_for_removal.append('spectral_index_error_too_small')
                            if len(reasons_for_removal) == 0:
                                reasons_for_removal = ['unidentified'] # Should not appear unless something has been coded wrong
                            removed_data.write(f'{stripped_line}, {"||".join(reasons_for_removal)}\n')

def flux_flux_plot_setup(source_a, source_b):
    pyplot.ylabel(f'Flux - {source_b}')
    pyplot.xlabel(f'Flux - {source_a}')
    pyplot.yscale('log')
    pyplot.xscale('log')
    return

def correlation_plot(source, suffix):
    source_dir, source_name = aux.source_underscored(source), aux.source_spaced(source)
    source_logger = Logger(source_name)
    source_logger.log('***CREATING COMPARATIVE NEARBY SOURCES PLOT***')
    nearby_source_list = select_nearby_sources(source_name)
    nearby_source_names = [aux.source_underscored(src["Source_Name"]) for src in nearby_source_list]
    source_logger.log(f'{len(nearby_source_list)} nearby sources found: {nearby_source_names}')
    missing_sources = [src for src in nearby_source_names if src not in SOURCE_NAMES]
    if len(missing_sources) > 0:
        source_logger.log(f'WARNING: Folders for the following sources {missing_sources} could not be found among downloaded source data. These sources will be ignored.')
        nearby_source_names = [src for src in nearby_source_names if src in SOURCE_NAMES]
    
    os.makedirs(correlation_plots_path(source_dir), exist_ok=True)
    correlation_table_name = f'{correlation_plots_path(source_dir)}/correlations.txt'
    with open(correlation_table_name, 'w') as correlation_table:
        correlation_table.write(', '.join(['source', 'nearby_source', 'pcc_daily', 'pcc_weekly', 'pcc_monthly']) + '\n')
        for nearby_source in nearby_source_names:
            pearson = { cadence:0 for cadence in CADENCES }
            flux_flux_plot_setup(source_name, aux.source_spaced(nearby_source))
            logger = Logger(source_name, flux_type='fixed')
            for cadence in CADENCES:
                source_table_path, nearby_source_table_path = (validated_table_path(src, 'fixed', cadence) for src in [source_dir, nearby_source])
                logger.log(f'Checking files: {source_table_path}, {nearby_source_table_path}')
                with open(source_table_path, 'r')as source_table_data, open(nearby_source_table_path, 'r') as nearby_source_table_data:
                    source_light_curve_temp = genfromtxt(source_table_data, skip_header=1, delimiter=',')
                    nearby_source_light_curve_temp = genfromtxt(nearby_source_table_data, skip_header=1, delimiter=',')

                if len(source_light_curve_temp) == 0 or len(nearby_source_light_curve_temp) == 0:
                    logger.log(f'Ignoring correlation data for empty table')
                    continue
                source_data, nearby_source_data = (SourceData(data) for data in common_timestamp_filter(source_light_curve_temp, nearby_source_light_curve_temp))

                pearson[cadence] = aux.pearson_correlation_coefficient(source_data.flux, nearby_source_data.flux)
                pyplot.errorbar(
                    source_data.flux, nearby_source_data.flux,
                    xerr=source_data.flux_error, yerr=nearby_source_data.flux_error,
                    label=f'LCR, {cadence}, pcc={pearson[cadence]:.3g}',
                    fmt='.', color=COLOR_BY_CADENCE[cadence], zorder=ZORDER_BY_CADENCE[cadence], lw=0.2,
                )
                pyplot.legend(loc='lower right', fontsize = 'xx-small')

            file_name = f'{correlation_plots_path(source_dir)}/{source_dir}_vs_{nearby_source}_correlation_plot{suffix}.png'
            pyplot.gcf().savefig(file_name, dpi=300)
            source_logger.log(f'Saved file {file_name}\n')
            pyplot.clf()

            correlation_table_row_data = [source_name, aux.source_spaced(nearby_source),
                                          f'{pearson["daily"]:.3g}', f'{pearson["weekly"]:.3g}', f'{pearson["monthly"]:.3g}']
            correlation_table.write(', '.join(correlation_table_row_data) + '\n')

    source_logger.log(f'Finished filling table {correlation_table_name}.')

def get_arguments():
    parser = ArgumentParser(description='Process source data to generate plots, validate data, add sun and moon distances, or compare nearby source data')
    parser.add_argument('sources', nargs='*', help='Input source names (in quotes if the name has spaces) to run the selected method only for specific sources')
    parser.add_argument('-m', '--method', default='all', choices=['all', 'tables', 'outliers', 'fluxflux', 'correlation', 'sunmoon'],
                        help='Choose a method to run, or use "all" to run all of them (in the order listed). Default: "all"')
    parser.add_argument('-e', '--free-dev', type=float, nargs=3, metavar=('DAILY', 'WEEKLY', 'MONTHLY'),
                        help='For "tables" and "outliers" methods, add the amount of std deviations allowed for FREE_flux tables, for DAILY, WEEKLY and MONTHLY data '
                        '(in that order) that will determine outliers. Default: 4(daily) 4(weekly) 5(monthly)')
    parser.add_argument('-x', '--fixed-dev', type=float, nargs=3, metavar=('DAILY', 'WEEKLY', 'MONTHLY'),
                        help='For "tables" and "outliers" methods, add the amount of std deviations allowed for FIXED_flux tables, for daily, weekly and monthly data '
                        '(in that order) that will determine outliers. Default: 4(daily) 4(weekly) 5(monthly)')
    parser.add_argument('--placeholder', default='NaN', help='For the methods "tables" and "sunmoon", choose the value that will represent an "empty cell".'
                        'Default: NaN')
    parser.add_argument('--suffix', default='', help='Add suffix for the plots that will be created e.g. --suffix=test will create plots named '
                        '"some_plot_test.png". Default: ""')
    parser.add_argument('-l', '--list', default=False, action='store_true', help='Does not run any method. Lists all available source folders.')

    parsed_args = parser.parse_args()
    parsed_args.suffix = '' if parsed_args.suffix == '' else f'_{parsed_args.suffix}'
    parsed_args.methods = ['tables', 'outliers', 'fluxflux', 'correlation', 'sunmoon'] if parsed_args.method == 'all' else [parsed_args.method]

    parsed_args.allowed_deviations = DEFAULT_ALLOWED_DEVIATIONS.copy()
    if parsed_args.free_dev is not None:
        parsed_args.allowed_deviations['free'] = {'daily': parsed_args.free_dev[0], 'weekly': parsed_args.free_dev[1], 'monthly': parsed_args.free_dev[2]}
    if parsed_args.fixed_dev is not None:
        parsed_args.allowed_deviations['fixed'] = {'daily': parsed_args.fixed_dev[0], 'weekly': parsed_args.fixed_dev[1], 'monthly': parsed_args.fixed_dev[2]}

    return parsed_args

# MAIN FUNCTION

if __name__ == '__main__':
    args = get_arguments()

    if args.list:
        aux.print_list_by_batches(SOURCE_NAMES, 10)
        sys.exit(0)

    if len(args.sources) > 0:
        argument_sources = [aux.source_underscored(source) for source in args.sources]
        selected_sources = [source for source in argument_sources if source in SOURCE_NAMES]

        invalid_sources = [source for source in argument_sources if source not in SOURCE_NAMES]
        if len(invalid_sources) > 0:
            print(f'WARNING: Folders for {invalid_sources} not found among downloaded source data. These sources will be ignored.')
    else:
        selected_sources = SOURCE_NAMES
    
    for method_name in args.methods:
        if method_name == 'tables':
            for source in selected_sources:
                validate_tables(source, placeholder=args.placeholder, allowed_deviations=args.allowed_deviations)
        elif method_name == 'outliers':
            for source in selected_sources:
                outliers_plot(source, allowed_deviations=args.allowed_deviations, suffix=args.suffix)
        elif method_name == 'fluxflux':
            for source in selected_sources:
                flux_flux_plot(source, args.suffix)
        elif  method_name == 'correlation':
            for source in selected_sources:
                correlation_plot(source, args.suffix)
        elif method_name == 'sunmoon':
            for source in selected_sources:
                sun_and_moon_distance(source, placeholder=args.placeholder)
        else:
            print(f'Incorrect method name: {method_name}')
