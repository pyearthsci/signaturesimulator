import netCDF4 as nc
import datetime as dt
import numpy as np
import warnings
#import pandas as pd


class StateVector:
    """Class to hold state vector data for optical and microwave canopy RT models.
    """
    def __init__(self):

        self.date_utc = None
        self.lai = None
        self.can_height = None
        self.leafChl = None
        self.leafWater = None
        self.soil_moisture = None
        self.soilAlbedo = None


def get_state_csv(fname):
    """Function that returns StateVector instance for a given file
    :param fname: Path to file from which to extract data, files must be csv with columns:
     (date, lai, canopy height, soil moisture)
    :type fname: str
    :return: StateVector instance
    """
    state_dat = np.loadtxt(fname, delimiter=',', dtype='str')
    state_inst = StateVector()
    state_inst.date_utc = [dt.datetime.strptime(dd, "%Y/%m/%d %H:%M") for dd in state_dat[:,0]]
    state_inst.lai = [float(lai) for lai in state_dat[:,1]]  # (m2 m-2)
    state_inst.can_height = [float(can_height) for can_height in state_dat[:,2]]  # (m)
    state_inst.soil_moisture = [float(soil_m) for soil_m in state_dat[:,3]]  # (m-3 m-3)
    return state_inst


def get_jules_state(nc_file, pft_idx=5):
    """Function that returns a stateVector instance for a given time.

    :param date_utc: datetime object of when to extract JULES output.
    :type date_utc: object
    :param nc_file: JULES output file from which to extract data.
    :type nc_file: str
    :return: Instance of stateVector class.
    :rtype: instance
    """
    nc_dat = nc.Dataset(nc_file, 'r')
    state_inst = StateVector()
    state_inst.date_utc = nc.num2date(nc_dat.variables['time'][:], nc_dat.variables['time'].units).tolist()
    state_inst.lai = nc_dat.variables['lai'][:, pft_idx, 0, 0].tolist()  # (m2 m-2)
    state_inst.can_height = nc_dat.variables['canht'][:, pft_idx, 0, 0].tolist()  # (m)
    state_inst.soil_moisture = (nc_dat.variables['smcl'][:, 0, 0, 0]/100).tolist() # (m-3 m-3)
    nc_dat.close()
    return state_inst


def get_jules_state_old(date_utc, nc_file):
    """Function that returns a stateVector instance for a given time.

    :param date_utc: datetime object of when to extract JULES output.
    :type date_utc: object
    :param nc_file: JULES output file from which to extract data.
    :type nc_file: str
    :return: Instance of stateVector class.
    :rtype: instance
    """
    nc_dat = nc.Dataset(nc_file, 'r')
    t_idx = nc.date2index(date_utc, nc_dat.variables['time'], select='nearest')
    state_inst = StateVector()
    state_inst.date_utc = nc.num2date(nc_dat.variables['time'][t_idx], nc_dat.variables['time'].units)
    state_inst.lai = nc_dat.variables['croplai'][t_idx, 0, 0, 0]  # (m2 m-2)
    state_inst.can_height = nc_dat.variables['cropcanht'][t_idx, 0, 0, 0]  # (m)
    state_inst.soil_moisture = nc_dat.variables['smcl'][t_idx, 0, 0, 0]  # (m-3 m-3)
    nc_dat.close()
    return state_inst


def find_nearest_date_idx(items, pivot):
    nearest = min(items, key=lambda x: abs(x - pivot))
    timedelta = abs(nearest - pivot)
    if timedelta.days > 31:
        warnings.warn('Your closest land state in time is specified over a month before your specified satellite ' \
                      'observation')
    idx = np.where(np.array(items)==nearest)[0][0]
    return idx, timedelta


#def get_date_list(year, month=1, days=365):
#    start_date = dt.datetime(year, month, 1, 12, 0)
#    date_list = pd.date_range(start_date, periods=days).tolist()
#    return date_list


def read(file_format='jules', file_str=None, year=None):
    """Reads output data to a dictionary of state vectors indexed by time.

    .. note:: This function requires sub-functions capable of reading specified file format.

    :param file_format: format of output to read.
    :type file_str: str
    :param file_str: location of file.
    :type file_str: str
    :param year: year of data to extract, if equal to None whole time series extracted
    :type year: int
    :return: state dictionary.
    :rtype: dict
    """
    if file_format == 'jules':
        state_dict = read_jules(file_str, year)
    else:
        state_dict = {}
    return state_dict


def read_jules(nc_file=None, year=None):
    """Reads jules output from netCDF file and writes it to a dictionary indexed by date.

    :param nc_file: location of nc_file.
    :type nc_file: str
    :param year: year of data to extract, if equal to None whole time series extracted.
    :type year: int
    :return: state dictionary.
    :rtype: dict
    """
    if nc_file is None:
        nc_file = 'jules/output/wallerfing_jules_1989_2012.nc'
        print("%s") % nc_file
    nc_dat = nc.Dataset(nc_file, 'r')
    state_dict = {}
    time = nc_dat.variables['time']
    if year is not None:
        strt_idx = nc.date2index(dt.datetime(year, 1, 1), time)
        end_idx = nc.date2index(dt.datetime(year, 12, 31), time)
        t_idx = np.arange(strt_idx, end_idx + 1)
        times = nc.num2date(time[strt_idx:end_idx], time.units)
    else:
        times = nc.num2date(time[:], time.units)
        t_idx = np.arange(len(times))
    for t in enumerate(times):
        state_dict[t[1]] = StateVector()
        state_dict[t[1]].lai = nc_dat.variables['croplai'][t_idx[t[0]], 0, 0, 0]  # (m2 m-2)
        state_dict[t[1]].can_height = nc_dat.variables['cropcanht'][t_idx[t[0]], 0, 0, 0]  # (m)
        state_dict[t[1]].soil_moisture = nc_dat.variables['smcl'][t_idx[t[0]], 0, 0, 0]  # (m-3 m-3)
        # state_dict[t[1]].soil_temp = nc_dat.variables['t_soil'][t_idx[t[0]], 0, 0, 0]  # (K)
        # figure out how to add others, which soil albedo to output?
    return state_dict
