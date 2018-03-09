# third party imports
import numpy as np
# signature simulator specific imports
import state_vector as sv
import backscatter as bs
# import sense code
from sense import model as sense_mod
from sense import soil as sense_soil
from sense import canopy as sense_canopy


def active_microwave_rt(state, geom):
    """
    Function that simulates SAR backscattering, given surface biogeophysical variables and viewing geometries
    :param state: instance of StateVector class
    :type state: object
    :param geom: instance of SensorGeometry class
    :type geom: object
    :return: instance of BackScatter class
    """
    backscat = bs.BackScatter()
    backscat.date_sat_ob = []
    backscat.date_land_ob = []
    backscat.soil_moisture = []
    backscat.lai = []
    backscat.can_height = []
    backscat.hh = []
    backscat.hv = []
    backscat.vv = []
    for date_utc in enumerate(geom.date_utc):
        idx, timedelt = sv.find_nearest_date_idx(state.date_utc, date_utc[1])
        SAR = run_sense(state.soil_moisture[idx], state.lai[idx], state.can_height[idx], geom.vza[date_utc[0]])
        backscat.date_sat_ob.append(date_utc[1])
        backscat.date_land_ob.append(state.date_utc[idx])
        backscat.soil_moisture.append(state.soil_moisture[idx])
        backscat.lai.append(state.lai[idx])
        backscat.can_height.append(state.can_height[idx])
        backscat.hh.append(10 * np.log10(SAR.__dict__['stot']['hh']))
        backscat.hv.append(10 * np.log10(SAR.__dict__['stot']['hv']))
        backscat.vv.append(10 * np.log10(SAR.__dict__['stot']['vv']))
    return backscat



def run_sense(soil_m, lai, can_height, vza, freq=5.405, stype='turbid_rayleigh', surf='Oh92', s=0.015, lai_coef=0.1,
                        omega=0.1):
    """
    Function that runs the SenSE SAR ScattEring model given some inputs
    :param soil_m: soil moisture (m3 m-3)
    :type soil_m: float
    :param lai: leaf area index (m2 m-2)
    :type lai: float
    :param can_height: canopy height (m)
    :type can_height: float
    :param vza: view zenith angle of sensor (degrees)
    :type vza: float
    :param freq: frequency (GHz)
    :type freq: float
    :param stype: canopy scattering model
    :type stype: str
    :param surf: surface scattering model
    :type surf: str
    :param s: surface rms height (m)
    :type s: float
    :param lai_coef: coefficient of lai for which to calculate extinction and volume scattering coefficients
    :type lai_coef: float
    :param omega: coefficient for calculation of volume scattering coefficients
    :type omega: float
    :return: instance of sense SingleScatRT class
    """
    """Function that runs sense model, given surface biogeophysical variables and viewing geometries
    """
    models = {'surface': surf, 'canopy': stype}
    SAR = sense_mod.SingleScatRT(
          surface=sense_soil.Soil(mv=soil_m, f=freq, s=s, clay=0.23, sand=0.27),
          # surface=sense_soil.Soil(eps=self.eps, f=self.freq, s=self.s),
          canopy=sense_canopy.OneLayer(ke_h=lai_coef * lai,  # extinction coefficient
                                       ke_v=lai_coef * lai,
                                       d=can_height,
                                       ks_h=omega * (lai_coef * lai),  # volume scattering coeff
                                       ks_v=omega * (lai_coef * lai)),
          models=models,
          theta=np.deg2rad(vza),
          freq=freq)
    SAR.sigma0()
    return SAR