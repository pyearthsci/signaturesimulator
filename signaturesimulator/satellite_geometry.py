from pyorbital.orbital import Orbital
import datetime as dt
import pyorbital.astronomy as ast
import numpy as np


class SensorGeometry:
    """Class to hold sun-sensor geometry information.
    """
    def __init__(self):
        self.date_utc = None
        # view zenith angle
        self.vza = None
        # view azimuth angle
        self.vaa = None
        # solar zenith angle
        self.sza = None
        # solar azimuth angle
        self.saa = None

    def print_geom(self):
        """Prints currently specified class attributes.
        """
        print(self.date_utc, self.vza, self.vaa, self.sza, self.saa)


def get_geom_csv(fname='/data/geometries/s2_example.csv'):
    """Function that returns StateVector instance for a given file
    :param fname: Path to file from which to extract data, files must be csv with columns:
     (date, lai, canopy height, soil moisture)
    :type fname: str
    :return: StateVector instance
    """
    geom_dat = np.loadtxt(fname, delimiter=',', dtype='str')
    geom_inst = SensorGeometry()
    geom_inst.date_utc = [dt.datetime.strptime(dd, "%Y/%m/%d %H:%M") for dd in geom_dat[:,0]]
    geom_inst.vza = [float(vza) for vza in geom_dat[:,1]]
    geom_inst.vaa = [float(vaa) for vaa in geom_dat[:,2]]
    geom_inst.sza = [float(sza) for sza in geom_dat[:,3]]
    geom_inst.saa = [float(saa) for saa in geom_dat[:,4]]
    return geom_inst


def get_satellite_geometry(startDateUTC, lengthDays, lon, lat, alt=0.0, mission="Sentinel-2a",
                           tleFile="/data/tle/norad_resource_tle.txt"):
    """Calculate approximate geometry for Sentinel overpasses.
    Approximate because it assumes maximum satellite elevation
    is the time at which target is imaged.

    :param startDateUTC: a datetime object specifying when to start prediction.
    :type startDateUTC: object
    :param lengthDays: number of days over which to perform calculations.
    :type lengthDays: int
    :param lat: latitude of target.
    :type lat: float
    :param lon: longitude of target.
    :type lon: float
    :param alt: altitude of target (in km).
    :type alt: float
    :param mission: mission name as in TLE file.
    :type mission: str
    :param tleFile: TLE file.
    :type tleFile: str
    :return: a python list containing instances of the sensorGeometry class arranged in date order.
    :rtype: list
    """
    orb = Orbital(mission, tleFile)
    passes = orb.get_next_passes(startDateUTC, 24 * lengthDays, lon, lat, alt)

    geom_inst = SensorGeometry()
    geom_inst.date_utc = []
    geom_inst.vza = []
    geom_inst.vaa = []
    geom_inst.sza = []
    geom_inst.saa = []

    for p in passes:
        look = orb.get_observer_look(p[2], lon, lat, alt)
        vza = 90 - look[1]
        vaa = look[0]
        sza = ast.sun_zenith_angle(p[2], lon, lat)
        saa = np.rad2deg(ast.get_alt_az(p[2], lon, lat)[1])

        if mission == 'Sentinel-1b':
            if 75 < vaa < 105 and 20. < vza < 45.:  # vaa usually [0, 180], testing observation times
                geom_inst.date_utc.append(p[2])
                geom_inst.vza.append(vza)
                geom_inst.vaa.append(vaa)
                geom_inst.sza.append(sza)
                geom_inst.saa.append(saa)
        elif mission == 'Sentinel-1a':
            if 75 < vaa < 105 and 20. < vza < 45.:  # vaa usually [0, 180], testing observation times
                geom_inst.date_utc.append(p[2])
                geom_inst.vza.append(vza)
                geom_inst.vaa.append(vaa)
                geom_inst.sza.append(sza)
                geom_inst.saa.append(saa)
        elif mission == 'Sentinel-2a':
            if sza < 90 and vza < 10.3:
                geom_inst.date_utc.append(p[2])
                geom_inst.vza.append(vza)
                geom_inst.vaa.append(vaa)
                geom_inst.sza.append(sza)
                geom_inst.saa.append(saa)

    return geom_inst


def get_sentinel_geometry(startDateUTC, lengthDays, lat, lon, alt=0.0, mission="Sentinel-2a",
                          tleFile="/data/tle/norad_resource_tle.txt"):
    """Calculate approximate geometry for Sentinel overpasses.
    Approximate because it assumes maximum satellite elevation
    is the time at which target is imaged.

    :param startDateUTC: a datetime object specifying when to start prediction.
    :type startDateUTC: object
    :param lengthDays: number of days over which to perform calculations.
    :type lengthDays: int
    :param lat: latitude of target.
    :type lat: float
    :param lon: longitude of target.
    :type lon: float
    :param alt: altitude of target (in km).
    :type alt: float
    :param mission: mission name as in TLE file.
    :type mission: str
    :param tleFile: TLE file.
    :type tleFile: str
    :return: a python list containing instances of the sensorGeometry class arranged in date order.
    :rtype: list
    """
    orb = Orbital(mission, tleFile)
    passes = orb.get_next_passes(startDateUTC, 24 * lengthDays, lon, lat, alt)

    geomList = []

    for p in passes:
        look = orb.get_observer_look(p[2], lon, lat, alt)
        vza = 90 - look[1]
        vaa = look[0]
        sza = ast.sun_zenith_angle(p[2], lon, lat)
        saa = np.rad2deg(ast.get_alt_az(p[2], lon, lat)[1])

        if mission == 'Sentinel-1b':
            if 75 < vaa < 105 and 20. < vza < 45.:  # vaa usually [0, 180], testing observation times
                thisGeom = SensorGeometry()
                thisGeom.date_utc = p[2]
                thisGeom.vza = vza
                thisGeom.vaa = vaa
                thisGeom.sza = sza
                thisGeom.saa = saa
                geomList.append(thisGeom)
        elif mission == 'Sentinel-1a':
            if 75 < vaa < 105 and 20. < vza < 45.:  # vaa usually [0, 180], testing observation times
                thisGeom = SensorGeometry()
                thisGeom.date_utc = p[2]
                thisGeom.vza = vza
                thisGeom.vaa = vaa
                thisGeom.sza = sza
                thisGeom.saa = saa
                geomList.append(thisGeom)
        elif mission == 'Sentinel-2a':
            if sza < 90 and vza < 10.3:
                thisGeom = SensorGeometry()
                thisGeom.date_utc = p[2]
                thisGeom.vza = vza
                thisGeom.vaa = vaa
                thisGeom.sza = sza
                thisGeom.saa = saa
                geomList.append(thisGeom)

    return geomList


if __name__ == "__main__":
    """Simple test with S2a for Wallerfing
    """

    startDate = dt.datetime(2017, 7, 1, 0, 0, 0, 0, None)
    lon, lat = 12.880, 48.684
    alt = 0.
    days = 30

    geomList = get_sentinel_geometry(startDate, days, lat, lon, mission="Sentinel-1a", alt=alt)
    for g in geomList:
        g.print_geom()
