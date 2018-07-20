# python builtins:
import os
import sys
import subprocess
from copy import copy
from tempfile import mkstemp

# third party imports:
import netCDF4
import numpy as np

# sentienl synergy specific imports:
import spectra as sp
import satellite_geometry as satgeo
import state_vector as sv


def passive_optical_rt(state, geom, mode='fast', rsl1=0.2, sm_coeff=0.5, cab=75.0, cw=0.01):
    """Function that simulates reflectances given surface biogeophysical variables and viewing geometries.

    :param state: Instance of the StateVector class.
    :type state: instance
    :param geom: Instance of the SensorGeometry class.
    :type geom: instance
    :param mode: Run semiDiscrete in either fast ('fast') or slow ('slow') mode [optional].
    :type resln: str
    :param rsl1: weight of the first soil vector
    :type rsl1: float
    :param sm_coeff: weighting of soil moisture impact, bound between (0,1)
    :type sm_coeff: float
    :param cab: leaf chlorophyl concentration
    :type cab: float
    :param cw: equivelant leaf water thickness
    :type cw: float
    :return: Instance of the spectra class.
    :rtype: instance
    """
    spect = sp.spectra()
    spect.date_sat_ob = []
    spect.date_land_ob = []
    spect.soil_moisture = []
    spect.lai = []
    spect.can_height = []
    spect.refl = []

    for date_utc in enumerate(geom.date_utc):
        idx, timedelt = sv.find_nearest_date_idx(state.date_utc, date_utc[1])
        reflectance = run_semidiscrete(state.soil_moisture[idx], state.lai[idx], state.can_height[idx],
                                       geom.vza[date_utc[0]], geom.vaa[date_utc[0]], geom.sza[date_utc[0]],
                                       geom.saa[date_utc[0]], mode=mode, rsl1=rsl1, sm_coeff=sm_coeff, cab=cab, cw=cw)
        spect.date_sat_ob.append(date_utc[1])
        spect.date_land_ob.append(state.date_utc[idx])
        spect.soil_moisture.append(state.soil_moisture[idx])
        spect.lai.append(state.lai[idx])
        spect.can_height.append(state.can_height[idx])
        spect.refl.append(reflectance)
    return spect


def run_semidiscrete(soil_m, lai, can_height, vza, vaa, sza, saa, mode='fast', rsl1=0.2, sm_coeff=0.5, cab=75.0,
                     cw=0.01):
    """A python wrapper to the SemiDiscrete optical canopy RT model of Nadine Gobron.

    :param soil_m: soil moisture at specified time (m3 m-3)
    :type soil_m: float
    :param lai: Leaf area index at specified time (m2 m-2)
    :type lai: float
    :param can_height: canopy height at specified time (m)
    :type can_height: float
    :param vza: view zenith angle (degrees)
    :type vza: float
    :param vaa: view azimuth angle (degrees)
    :type vaa: float
    :param sza: solar zenith angle (degrees)
    :type sza: float
    :param saa: solar azimuth angle (degrees)
    :type saa: float
    :param mode: Run semiDiscrete in either fast ('fast') or slow ('slow') mode [optional].
    :type resln: str
    :param rsl1: weight of the first soil vector
    :type rsl1: float
    :param sm_coeff: weighting of soil moisture impact, bound between (0,1)
    :type sm_coeff: float
    :param cab: leaf chlorophyl concentration
    :type cab: float
    :param cw: equivelant leaf water thickness
    :type cw: float
    :return: Reflectance values for specified input.
    :rtype: array
    """

    # generate a tmp file and write geometry into it
    dir_path = os.path.dirname(os.path.realpath(__file__))
    fd, temp_path = mkstemp(prefix='/tmp/senSyntmp__', text=True)
    tmpFile = os.fdopen(fd, 'w')
    print >> tmpFile, "%s %s %s %s" % (vza, vaa, sza, saa)
    # print "%s %s %s %s" %(geom.vza,geom.vaa,geom.sza,geom.saa)
    tmpFile.close()

    # set up nadim command line
    if mode == 'fast':
        cmd = dir_path+"/semidiscrete_srf/semiD -srf "+dir_path+"/data/srf/s2a.srf -fast"
    else:
        cmd = dir_path+"/semidiscrete_srf/semiD -srf "+dir_path+"/data/srf/s2a.srf"
    if lai != None:
        cmd = cmd + " -LAI %f -hc %f -rsl1 %f -cab %f -cw %f" % (lai, can_height, rsl1 * (1. - sm_coeff * soil_m), cab,
                                                                 cw)
        # Think about soil moisture implementation here
    cmd = cmd + " < %s" % temp_path
    # print cmd

    # run process
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=sys.stderr, shell=True)
    out = p.stdout.readlines()
    p.wait()

    reflectance = np.array([float(val) for val in out[0].split()])
    return reflectance


def canopy_rt_optical(state, geom, resln=1.0):
    """A python wrapper to the SemiDiscrete optical
    canopy RT model of Nadine Gobron. Runs the
    model for the the whole of its valid spectra
    range at a resolution set by resln.

    :param state: Instance of the stateVector class.
    :type state: instance
    :param geom: Instance of the sensorGeomety class.
    :type geom: instance
    :param resln: the spectral resolution in nm [optional].
    :type resln: float
    :return: Instance of the spectra class.
    :rtype: instance
    """

    spect = sp.spectra()
    spect.wavl = np.arange(400, 2500 + resln, resln)

    # generate a tmp file and write
    # wavelengths and geometry into it
    fd, temp_path = mkstemp(prefix='/tmp/senSyntmp__', text=True)
    tmpFile = os.fdopen(fd, 'w')
    print >> tmpFile, "1 %d" % len(spect.wavl),
    for w in spect.wavl:
        print >> tmpFile, " %f" % w,
    print >> tmpFile, "\n%s %s %s %s" %(geom.vza,geom.vaa,geom.sza,geom.saa)
    tmpFile.close()

    # set up nadim command line
    cmd = "semiDiscrete"
    if state.lai != None:
        cmd = cmd + " -LAI %f -hc %f -rsl1 %f" %(state.lai, state.can_height, 0.2*(1.-0.5*(state.soil_moisture/100.)))
        # Think about soil moisture implementation here
    cmd = cmd + " < %s" % temp_path

    # run process
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=sys.stderr, shell=True)
    out = p.stdout.readlines()
    p.wait()

    # clean up
    os.remove(temp_path)

    # process RT model output
    rtOut = out[1].split()
    for i in xrange(4, len(rtOut)):
        reflTmp = np.append(spect.refl, float(rtOut[i]))
        spect.refl = copy(reflTmp)

    return spect


def canopy_rt_optical_fast(state, geom, resln=1.0, mode='fast'):
    """A python wrapper to the SemiDiscrete optical
    canopy RT model of Nadine Gobron. Runs the
    model for the the whole of its valid spectra
    range at a resolution set by resln.

    :param state: Instance of the StateVector class.
    :type state: instance
    :param geom: Instance of the SensorGeometry class.
    :type geom: instance
    :param mode: Run semiDiscrete in either fast ('fast') or slow ('slow') mode [optional].
    :type resln: str
    :return: Instance of the spectra class.
    :rtype: instance
    """
    spect = sp.spectra()
    spect.wavl = np.arange(400, 2500 + resln, resln)

    # generate a tmp file and write geometry into it
    fd, temp_path = mkstemp(prefix='/tmp/senSyntmp__', text=True)
    tmpFile = os.fdopen(fd, 'w')
    print >> tmpFile, "%s %s %s %s" %(geom.vza,geom.vaa,geom.sza,geom.saa)
    # print "%s %s %s %s" %(geom.vza,geom.vaa,geom.sza,geom.saa)
    tmpFile.close()

    # set up nadim command line
    if mode == 'fast':
        cmd = "semiD -srf ../data/srf/s2a.srf -fast"
    else:
        cmd = "semiD -srf ../data/srf/s2a.srf"
    if state.lai != None:
        cmd = cmd + " -LAI %f -hc %f -rsl1 %f" %(state.lai, state.can_height, 0.2*(1.-0.5*(state.soil_moisture/100.)))
        # Think about soil moisture implementation here
        # include -cab (leaf chlorophyl concentration, default:75.0) and -cw (equivelant leaf water thickness, default: 0.01)??
        # update 0.2 to tunable rsl1 (weight of the first soil vector, default: 0.2)
    cmd = cmd + " < %s" % temp_path
    # print cmd

    # run process
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=sys.stderr, shell=True)
    out = p.stdout.readlines()
    p.wait()
    spect.refl = np.array([float(val) for val in out[0].split()])
    return spect


if __name__ == "__main__":
    """This example opens a test output file from the JULES
    model, reads in LAI data, runs these through the optical RT
    model at given geometry and convoles the resulting spectra
    with Sentinel 2 spectra response functions.
    """

    from matplotlib import pyplot as plt

    # read example LAI data:
    allLai = netCDF4.Dataset('../testData/crop_germ_gl4.day.nc').variables['lai'][-365:, 5, 0, 0]

    # main classes:
    state = sv.stateVector()
    geom = satgeo.sensorGeometry()
    geom.sza = 30.

    # container for output:
    allSpect = []
    allBRF = []

    # orbit revist:
    revist = 10

    # plotting variables:
    xpnts = []
    lgnd = []

    # loop over all states and call the
    # canopy RT model and convolve the
    # retruned spectra to S2 bands
    for (n, L) in enumerate(allLai[::revist]):
        state.lai = L
        spect = canopy_rt_optical(state, geom)
        allSpect.append(sp.sentinel2(spect))
        allBRF.append(allSpect[n].refl)
        xpnts.append(n * revist)

    allBRF = np.array(allBRF)

    # sort out legend
    for i in xrange(len(allBRF[0, :])):
        lgnd.append('S2 band %d' % (i + 1))

    # do plots:
    lineObjs = plt.plot(xpnts, allBRF)
    plt.ylabel('reflectance (BRF)')
    plt.xlabel('Day of year')
    plt.xlim([0, 364])
    for i in xrange(7, len(lineObjs)):
        lineObjs[i].set_dashes([3, 1])
    plt.legend(iter(lineObjs), lgnd)
    plt.show()
    # plt.savefig('s2Sim_test1.png')
