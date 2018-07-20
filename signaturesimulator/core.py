# core module imports
import os
# third party imports
import numpy as np
import datetime as dt
import matplotlib.pyplot as plt
# signature simulator specific imports
import opticalcanopy_rt as opcan_rt
import activemicro_rt as actmicro_rt
import atmospheric_rt as atmo_rt
import satellite_geometry as satgeo
import state_vector as sv
import f90nml

# Define some example factory functions:


def single_reflectance(vza=5.5, vaa=286.3, sza=26.8, saa=157.0, lai=3.0, canopy_ht=1.0, soil_m=0.3):
    """
    Factory function computing a single reflectance given viewing geometries of sensor and land surface state variables
    :param vza: view zenith angle (deg)
    :param vaa: view azimuth angle (deg)
    :param sza: solar zenith angle (deg)
    :param saa: solar azimuth angle (deg)
    :param lai: leaf area index (m2 m-2)
    :param canopy_ht: canopy height (m)
    :param soil_m: soil moisture (m3 m-3)
    :return: simulated Sentinel 2 reflectance for given input
    """
    s = Simulator()
    s.get_geom = s.geom_default(vza=vza, vaa=vaa, sza=sza, saa=saa)
    s.get_land_state = s.state_default(lai=lai, canopy_ht=canopy_ht, soil_m=soil_m)
    s.run_rt = s.passive_optical
    s.run()
    return s.spectra.refl


def multi_backscat(geom_csv='/data/geometries/s1_example_const.csv',
                   state_csv='/data/state_variables/state_example.csv'):
    """
    Factory function computing multiple values of backscatter given viewing geometries of sensor and land surface state
    variables specified in two csv files
    :param geom_csv: location of csv containing viewing geometries and times
    :param state_csv: location of csv containing land surface state variables and times
    :return: simulated Sentinel 1 backscatter in the hv polarisation
    """
    s = Simulator()
    s.get_geom = s.geom_csv(geom_csv)
    s.get_land_state = s.state_csv(state_csv)
    s.run_rt = s.active_microwave
    s.run()
    return s.backscat.hv


class Simulator(object):
    """Class to simulate satellite observations.
    """
    def __init__(self, site_nml=None):
        """Calculate class attributes.
        """
        # find file path
        self.dir_path = os.path.dirname(os.path.realpath(__file__))
        # setup site parameter dictionary
        self.site_param_dic = {'lat':12.88, 'lon': 48.684, 'clay': 0.23, 'sand': 0.27, 'bulk': 1.65, 'freq': 5.405,
                               's': 0.015, 'lai_coeff': 0.1, 'omega': 0.1, 'mode': 'fast', 'rsl1': 0.2, 'sm_coeff': 0.5,
                               'cab': 75.0, 'cw': 0.01}
        # update parameters with values specifed in site.nml
        if site_nml is not None:
            self.set_site_params(site_nml)
        else:
            self.set_site_params(self.dir_path + '/site.nml')
        # setup land surface state variables
        self.get_land_state = self.state_default()
        # setup geometry for observations
        self.get_geom = self.geom_default()
        # setup radiative transfer model
        self.run_rt = self.passive_optical  # choose from self.active_microwave or self.passive_optical
        # setup atmospheric RT model for passive optical sensors
        self.run_atmos_rt = atmo_rt.fwd_atm
        # setup keys for plotting variables
        self.output_variables = []
        # paths to example geometry and land state csv files
        self.example_s1_geometries = self.dir_path + '/' + 'data/geometries/s1_example_const.csv'
        self.example_s2_geometries = self.dir_path + '/' + 'data/geometries/s2_example_const.csv'
        self.example_state_variables = self.dir_path + '/' + 'data/state_variables/state_example.csv'


    def state_default(self, date_utc=dt.datetime(2016, 6, 17, 9, 0), lai=3.0, canopy_ht=1.0, soil_m=0.3):
        """
        Defines a default land state given some input
        :param date_utc: date as a datetime object
        :type date_utc: object
        :param lai: Leaf area index at specified time (m2 m-2)
        :type lai: float
        :param canopy_ht: canopy height at specified time (m)
        :type canopy_ht: float
        :param soil_m: soil moisture at specified time (m3 m-3)
        :type soil_m: float
        :return: instance of StateVector class
        """
        state_inst = sv.StateVector()
        state_inst.date_utc = [date_utc]
        state_inst.lai = [lai]
        state_inst.can_height = [canopy_ht]
        state_inst.soil_moisture = [soil_m]
        return state_inst

    def state_csv(self, fname):
        """
        Extract data from a .csv file and create a state instance.
        :param fname: Filename and location of .csv to extract data from
        :type fname: str
        :return: instance of StateVector class
        """
        state_inst = sv.get_state_csv(fname)
        return state_inst

    def state_jules(self, land_cover='crop', pft_idx=5):
        """Extracts state variables from either a) a specified JULES netcdf file, or b) a default land cover profile
        from the list: ['crop', 'evergreen', 'broadleaf', 'grassland']
        :param land_cover: Either a filepath to JULES output netcdf file or a string from list:
        ['crop', 'evergreen', 'broadleaf', 'grassland']
        :type land_cover: str
        :return: instance of StateVector class
        """
        if land_cover in ['crop', 'crop_early', 'crop_late']:
            state_inst = sv.get_jules_state(self.dir_path+'/'+'../data/state_variables/'+land_cover+'.nc')
        else:
            state_inst = sv.get_jules_state(land_cover, pft_idx)
        return state_inst

    def geom_default(self, date_utc=dt.datetime(2016, 6, 17, 10, 25), vza=5.5, vaa=286.3, sza=26.8, saa=157.0):
        """
        Defines geometry for satellite observation
        :param date_utc: date of observation as a datetime object
        :type date_utc: object
        :param vza: view zenith angle (degrees)
        :type vza: float
        :param vaa: view azimuth angle (degrees)
        :type vaa: float
        :param sza: solar zenith angle (degrees)
        :type sza: float
        :param saa: solar azimuth angle (degrees)
        :type saa: float
        :return: instance of SensorGeometry class
        """
        geom_inst = satgeo.SensorGeometry()
        geom_inst.date_utc = [date_utc]
        geom_inst.vza = [vza]
        geom_inst.vaa = [vaa]
        geom_inst.sza = [sza]
        geom_inst.saa = [saa]
        return geom_inst

    def geom_csv(self, fname):
        """
        Extracts dates and viewing geomteries from a .csv and creates an instance of the SensorGeometry class
        :param fname: Filename and location of .cvs file containing dates and viewing angles for sensor
        :type fname: str
        :return: instance of SensorGeometry class
        """
        geom_inst = satgeo.get_geom_csv(fname)
        return geom_inst

    def geom_pyoribtal(self, start_date=dt.datetime(2016, 1, 1), num_days=365, altitude=0.0,
                       sat_mission='S2', tle_file=None):
        """
        Uses pyorbital to calculate sensor geometries for a specified time period and location on the earth
        :param start_date: start date for period in which to get sensor geometries
        :type start_date: datetime object
        :param num_days: number of days in period to get sensor geometries
        :type num_days: int
        :param lon: location longitude
        :type lon: float
        :param lat: location latitude
        :type lat: float
        :param altitude: location altitude
        :type altitude: float
        :param sat_mission: satellite mission
        :type sat_mission: str
        :param tle_file: tle file location with parameters for given sensor satellite mission
        :type tle_file: str
        :return: instance of SensorGeometry class
        """
        mission_dict = {'S2': 'Sentinel-2a', 'S1': 'Sentinel-1b'}
        if tle_file is None:
            tle_file = self.dir_path + "/data/tle/norad_resource_tle.txt"
        geom_inst = satgeo.get_satellite_geometry(start_date, num_days, self.site_param_dic['lon'],
                                                  self.site_param_dic['lat'], alt=altitude,
                                                  mission=mission_dict[sat_mission], tleFile=tle_file)
        return geom_inst

    def set_site_params(self, site_nml):
        """
        Function that sets values of site_param_dic from values specified in a fortran namelist file
        :param site_nml: path to fortran nml file
        :return: None
        """
        nml_dic = f90nml.read(site_nml)
        for key in nml_dic['site_params'].keys():
            if type(nml_dic['site_params'][key]) == float:
                self.site_param_dic[key] = nml_dic['site_params'][key]
            elif type(nml_dic['site_params'][key]) == list:
                self.site_param_dic[key] = nml_dic['site_params'][key][0]

    def passive_optical(self, state, geom,):
        """
        Function that runs semi discrete rt model for passive optical sensors and returns instance of the spectra class
        :param state: instance of the StateVector class
        :type state: object
        :param geom: instance of the SensorGeometry class
        :type geom: object
        :return: instance of the spectra class
        """
        spectra = opcan_rt.passive_optical_rt(state, geom, mode=self.site_param_dic['mode'],
                                              rsl1=self.site_param_dic['rsl1'],
                                              sm_coeff=self.site_param_dic['sm_coeff'], cab=self.site_param_dic['cab'],
                                              cw=self.site_param_dic['cw'])
        return spectra

    def active_microwave(self, state, geom,):
        """
        Function that runs sense model for active microwave sensors and returns instance of the backscatter class
        :param state: instance of the StateVector class
        :type state: object
        :param geom: instance of the SensorGeometry class
        :type geom: object
        :return: instance of the backscatter class
        """
        backscat = actmicro_rt.active_microwave_rt(state, geom, freq=self.site_param_dic['freq'],
                                                   s=self.site_param_dic['s'],
                                                   lai_coeff=self.site_param_dic['lai_coeff'],
                                                   omega=self.site_param_dic['omega'], clay=self.site_param_dic['clay'],
                                                   sand=self.site_param_dic['sand'])
        return backscat

    def run(self):
        """
        Class method that runs the simulator for specifed viewing geometries and land surface state variables
        """
        if self.run_rt == self.passive_optical:
            self.spectra = self.run_rt(self.get_land_state, self.get_geom)
            self.spectra = self.run_atmos_rt(self.spectra)
            self.output_variables += self.spectra.__dict__.keys()
            self.output_variables.remove('ftype')
            self.output_variables.remove('wavl')
        elif self.run_rt == self.active_microwave:
            self.backscat = self.run_rt(self.get_land_state, self.get_geom)
            self.output_variables += self.backscat.__dict__.keys()

    def plot(self, plot_key, band_idx=None, S2=1):
        """
        Function to plot output of simulator after completing a run
        :param plot_key: plot key of variable to plot (must be in self.plot_keys list)
        :type plot_key: str
        :param band_idx: index to reflectance band to plot if plot key "refl" is specified
        :type band_idx: int
        :return: requested plot
        """
        try:
            if plot_key not in self.output_variables:
                return "Plot key is not in specified plot keys list for simulator output"
            elif plot_key == "refl" and band_idx == None:
                return "For plot key refl you must also specify a band index using, band_idx=..."
            elif plot_key == "refl":
                plt.plot(self.spectra.date_sat_ob, np.array(self.spectra.__dict__["refl"])[:, band_idx], 'o')
                plt.xlabel('Date')
                if S2 == 1:
                    S2_band_labels = ['1', '2', '3', '4', '5', '6', '7', '8', '8a', '9', '10', '11', '12']
                    plt.ylabel('Band '+S2_band_labels[band_idx]+' reflectance')
                elif S2 == 0:
                    plt.ylabel('Band reflectance')
                #plt.show()
            elif self.run_rt == self.active_microwave:
                plt.plot(self.backscat.date_sat_ob, self.backscat.__dict__[plot_key], 'o')
                plt.xlabel('Date')
                plt.ylabel(plot_key)
                #plt.show()
            elif self.run_rt == self.passive_optical:
                plt.plot(self.spectra.date_sat_ob, self.spectra.__dict__[plot_key], 'o')
                plt.xlabel('Date')
                plt.ylabel(plot_key)
                #plt.show()
            else:
                return "Something went wrong, please check that you have run the simulator and are specifying a valid" \
                       "plot key"
        except AttributeError:
            print('Please run simulator before plotting, using e.g. simulator_instance.run()')