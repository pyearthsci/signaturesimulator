signaturesimulator
==================

This Python package simulates satellite data from passive optical and active microwave sensors based on land surface
biogeophysical state variables (leaf area index, canopy height and soil moisture) and viewing geometries and times.
Here we provide documentation on this package. Below we illustrate a few basic examples of installing the package and
its use.

Installing the package
----------------------

To install the signaturesimulator Python package using pip simply run the following from the command line::

    pip install signaturesimulator

Example of calculating reflectance for a single acquisition
-----------------------------------------------------------

From an interactive console we first import the signaturesimulator package and setup an instance of the Simulator
class::

    import signaturesimulator as ss
    sim = ss.Simulator()

We can now set the viewing geometries and time of the satellite observation, to do this we also need to specify a
Python datetime object so we will also import the datetime module::

    import datetime as dt
    sim.get_geom = sim.geom_default(date_utc=dt.datetime(2016,6,17,10,25), vza=5.5, vaa=286.3, sza=26.8, saa=157.0)

Where vza and vaa are the view zenith and azimuth angles and sza and saa are the solar zenith and azimuth angles. We
must now specify the land surface state variables for which to calculate reflectance::

    sim.get_land_state = sim.state_default(date_utc=dt.datetime(2016,6,17,9,0), lai=3.0, canopy_ht=1.0, soil_m=0.3)

Where lai is leaf are index (m2 m-2), canpopy_ht is vegetation canopy height (m) and soil_m is the top layer soil
moisture (m3 m-3). We now specify the radiative transfer model we want to use and then run the simulator (currently
there is the choice between active microwave and passive optical with these models being setup for the Sentinel 1 and 2
missions respectively). Here we are interested in reflectance, so set the following::

    sim.run_rt = sim.passive_optical
    sim.run()

.. note::  Currently the passive optical model will only work if the signaturesimulator is installed on a linux machine,
    due to the packaged Semi-Discrete binaries. We are working to bring out a version compatiable with Windows and Mac
    OSx. Please see the next example if you are not on a linux machine. The active microwave model is usable on any
    architecture.


Once the simulator has finished its run it will write the output to a new spectra class. We can then look at the
reflectance values for different bands (ordered by wavelength) using the following command::

    sim.spectra.refl

Which for this example outputs::

    [array([ 0.018664,  0.020378,  0.026745,  0.02097 ,  0.035889,  0.213055,
             0.332095,  0.338523,  0.343047,  0.337906,  0.230565,  0.192891,
             0.086688])]

Next we show how to run the simulator for an active microwave sensor over a time series of specified geometries and
land surface state variables.

Example of calculating backscatter for multiple aquisitions
-----------------------------------------------------------

Using the Simulator instance specified in the previous example we show how we can run the simulator for multiple
aquisitions. First we specify the dates and view geometries of our observations using a csv file::

    sim.get_geom = sim.geom_csv(fname=sim.example_s1_geometries)

Here :code:`sim.example_s1_geometries` points to an example csv file included with the signaturesimulator package found
here::

$PYTHONPATH/signaturesimulator/data/geometries/s1_example_const.csv

Any geometry csv file must be of the following format::

    # date, vza, vaa, sza, saa
    2016/01/03 05:23,34.3773946043,100.545700717,105.298744327,107.406448412
    2016/01/08 05:31,23.4284120953,102.103838392,103.928256857,108.076934788
    ...

We next specify the dates and values for our land surface state variables using a csv file::

    sim.get_land_state = sim.state_csv(fname=sim.example_state_variables)

Again we are using an example csv file packaged with the simulator, for the land state csv files must follow the
format::

    # date, lai (m2 m-2), canopy_height (m), soil_moisture (m3 m-3)
    2016/01/01 12:00,0.001,0.001,0.339560508728
    2016/01/02 12:00,0.001,0.001,0.341959953308
    ...

We can now set the radiative transfer model we wish to use and run the simulator::

    sim.run_rt = sim.active_microwave
    sim.run()

Once the simulator has finished its run it will write the output to a new backscat class instance. We can see the output
variables available to us using::

    sim.output_variables

In this case this returns::

    ['lai', 'date_sat_ob', 'soil_moisture', 'hv', 'can_height', 'date_land_ob', 'hh', 'vv']

To inspect the output in the hv polarisation we use::

    sim.backscat.hv

Which for this example returns::

    [-14.839751698875441, -14.612031628695206, ...,  -14.647495031040052, -14.470503894767003]

To plot the backscatter in the hv polaristation we can use the following command::

    sim.plot('hv')

Which will return the plot:

.. image:: s1_hv.png

We can plot any of the output variables using the plot method of the Simulator class, for LAI::

    sim.plot('lai')

Returning:

.. image:: s1_lai.png

Here we can see the effect that leaf area index is having on the simulated observations of backscatter.

Source Code
-----------

www.github.com/pyearthsci/signaturesimulator

Support
-------

If you are having issues, please let us know.
Contact: e.pinnington@reading.ac.uk

License
-------

Details of licensing information. TBC.