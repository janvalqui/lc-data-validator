# Contact: Janeth Valverde <valverde@llr.in2p3.fr>

from math import inf as infinite
from numpy import amin, amax, genfromtxt, ones
from astropy.coordinates import SkyCoord, get_body
from astropy.time import Time
from astropy import units
from math import sqrt, acos, cos, sin, pi

CADENCES = ['daily', 'weekly', 'monthly']
FLUX_TYPES = ['fixed', 'free']

def source_spaced(source):
    return source.replace('_', ' ')

def source_underscored(source):
    return source.replace(' ', '_')

def print_list_by_batches(target_list, batch_size):
    list_size = len(target_list)
    for index in range(list_size // batch_size):
        print(target_list[batch_size*index : batch_size*(index + 1)])
    if list_size % batch_size > 0:
        print(target_list[(list_size - list_size % batch_size) : list_size])

def mean(array):
    return sum(array)/len(array)

def std_dev(array):
    return sqrt(sum((array-mean(array)*ones(len(array)))**2)/len(array))

def pearson_correlation_coefficient(array_a, array_b):
    return sum((array_a-mean(array_a))*(array_b-mean(array_b)))/(len(array_a)*std_dev(array_a)*std_dev(array_b))

def convert_met_to_mjd(met_time):
    return 51910 + met_time/86400

def convert_eq_to_ecl(right_ascension, declination):
    ecl_coordinates = SkyCoord(right_ascension, declination, unit=units.deg).transform_to('barycentrictrueecliptic')
    return ecl_coordinates.lon.deg, ecl_coordinates.lat.deg

def angular_distance(ascension_1, declination_1, ascension_2, declination_2):
    radian_distance = acos(sin(pi*declination_2/180)*sin(pi*declination_1/180) + cos(pi*declination_2/180)*cos(pi*declination_1/180)*cos(pi*(ascension_2 - ascension_1)/180))
    return 180*radian_distance/pi

def is_source_within_ecliptic_band(source_data, band_width=18):
    source_ra, source_dec = float(source_data['RAJ2000']), float(source_data['DEJ2000'])

    # Since the moon can be at most ~5° away, any source 18° away from the ecliptic can never have the sun nor moon in the ROI.
    _longitude_from_ecliptic, latitude_from_ecliptic = convert_eq_to_ecl(source_ra, source_dec)
    return latitude_from_ecliptic <= band_width

def minimal_distance_to_body(source_data, body, timerange):
    # Avoid doing any calculations if we know source is not within 18° of the ecliptic
    if not is_source_within_ecliptic_band(source_data, 18):
        return infinite
    
    source_ra, source_dec = float(source_data['RAJ2000']), float(source_data['DEJ2000'])
    
    distance_from_body = lambda position : angular_distance(source_ra, source_dec, position.ra.deg, position.dec.deg)
    
    mid_point_index = len(timerange) // 2 # Works for range >= 2

    mid_mjd_timerange = [convert_met_to_mjd(timestamp) for timestamp in timerange[mid_point_index-1:mid_point_index+1]]
    mid_range_body_positions = get_body(body, Time(mid_mjd_timerange, format='mjd'))
    # mid_range_body_positions = [SkyCoord(194,0,1,unit=units.deg), SkyCoord(194,0,1,unit=units.deg)]

    lower_range_starting_distance, upper_range_starting_distance = distance_from_body(mid_range_body_positions[0]), distance_from_body(mid_range_body_positions[1])
    if lower_range_starting_distance > upper_range_starting_distance:
        truncated_timerange = timerange[mid_point_index:len(timerange)]
        truncated_index_range = range(0,len(truncated_timerange))
    else:
        truncated_timerange = timerange[0:(mid_point_index + len(timerange) % 2)] # "+ len(timerange) % 2" will include the middle point if length is odd
        truncated_index_range = range(len(truncated_timerange)-1, -1, -1) # evaluating in reverse to start from the original mid-point

    mjd_timerange = [convert_met_to_mjd(timestamp) for timestamp in truncated_timerange]
    body_positions = get_body(body, Time(mjd_timerange, format='mjd'))
    # body_positions = [SkyCoord(194,0,1,unit=units.deg)]*len(truncated_timerange)

    minimum = upper_range_starting_distance
    for index in truncated_index_range:
        position = body_positions[index]
        distance_at_timestamp = distance_from_body(position)
        if distance_at_timestamp > minimum:
            return minimum
        else:
            minimum = distance_at_timestamp

    return minimum

class Logger:
    source = None
    flux_type = None
    cadence = None

    def __init__(self, source=None, flux_type=None, cadence=None):
        self.source = source_spaced(source) if source is not None else None
        self.flux_type = flux_type
        self.cadence = cadence
    
    def log(self, message=None):
        if message is None:
            print()
        else:
            prefix = ', '.join(str(value) for value in [self.source, self.flux_type, self.cadence] if value is not None)
            print(f'{prefix}{": " if len(prefix) > 0 else ""}{message}')

class SourceData:
    source_data = None

    def __init__(self, source_data):
        self.source_data = source_data

    @property
    def timestamp(self):
        return self.source_data[:, 0]
    @property
    def ts(self):
        return self.source_data[:, 1]
    @property
    def flux(self):
        return self.source_data[:, 2]
    @property
    def flux_error(self):
        return self.source_data[:, 3]
    @property
    def spectral_index(self):
        return self.source_data[:, 4]
    @property
    def spectral_index_error(self):
        return self.source_data[:, 5]
    @property
    def fit_tolerance(self):
        return self.source_data[:, 6]
    @property
    def quality(self):
        return self.source_data[:, 7]

class Limits:
    lower_limit=None
    upper_limit=None

    def __init__(self, value_array=None):
        if value_array is not None:
            self.lower_limit, self.upper_limit = amin(value_array), amax(value_array)
        else:
            self.lower_limit, self.upper_limit = (infinite, -infinite)

    def __repr__(self):  
        return f"Limits Object (lower={self.lower_limit}, upper={self.upper_limit})"
    
    def update(self, value_array):
        self.lower_limit = min(self.lower_limit, amin(value_array))
        self.upper_limit = max(self.upper_limit, amax(value_array))
    
    def as_tuple(self, with_margins=False, in_log_scale=False):
        if self.lower_limit == infinite:
            return (-1,1)
        if not with_margins:
            return (self.lower_limit, self.upper_limit)
        if not in_log_scale:
            return self.lower_limit - 0.1*(self.upper_limit - self.lower_limit), self.upper_limit + 0.1*(self.upper_limit - self.lower_limit)
        # Only use in_log_scale=True for strictly positive ranges:
        return (self.lower_limit/((self.upper_limit/self.lower_limit)**(0.1))), (self.upper_limit*((self.upper_limit/self.lower_limit)**(0.1)))

######################################################################### 

def load_source_data(filename):
    with open(filename) as file:
        light_curve_raw = genfromtxt(file, skip_header=1, delimiter=',')
    return SourceData(light_curve_raw)

class OutliersPlot:
    axs = None

    def __init__(self, axs):
        self.axs = axs
    
    @property
    def relative_flux(self):
        return self.axs[0,0]
    @property
    def ts_dist(self):
        return self.axs[0,1]
    @property
    def flux_error_log_dist(self):
        return self.axs[1,0]
    @property
    def flux_log_dist(self):
        return self.axs[1,1]
    @property
    def spectral_index_error_dist(self):
        return self.axs[2,0]
    @property
    def spectral_index_dist(self):
        return self.axs[2,1]

