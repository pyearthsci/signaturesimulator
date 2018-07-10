# !/home/db903833//dataLand01/enthoughtDistros/epd-7.2-2-rh5-x86_64/bin/python
# !/usr/bin/env python

# core python modules:
import subprocess
import os
# 3rd party modules:
import numpy as np
import sys
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import netCDF4 as nc
import f90nml
import glob


def crop_run_jules(year=2016, sow_date=110, run_id='crop_jules', output_dir='default'):
    """
    Factory function running the JULES land surface model and updating the value of sow date for crops.
    :param year: year to run JULES for, here only 2015 and 2016 available (str)
    :param sow_date: crop sow date (int)
    :param run_id: id for model run (str)
    :param output_dir: directory to store output (str)
    :return: 'run completed'
    """
    dir_path = os.path.dirname(os.path.realpath(__file__))
    example_nml = dir_path + '/' + 'data/jules/example_nml'  # example nml directory bundled with simulator
    if output_dir == 'default':
        output_dir = dir_path + '/data/jules/output'  # if default use example directory to store model output
    jules = Jules(nml_dir= example_nml, print_out = 0)
    jules.nml_dic['timesteps']['jules_time']['main_run_start'] = str(year) + '-01-01 00:00:00'
    jules.nml_dic['timesteps']['jules_time']['main_run_end'] = str(year) + '-12-31 21:00:00'
    jules.nml_dic['ancillaries']['jules_crop_props']['const_val'][2] = sow_date  # set crop sow date
    jules.nml_dic['output']['jules_output']['run_id'] = run_id  # set run id
    jules.nml_dic['output']['jules_output']['output_dir'] = output_dir  # set directory to store model output
    jules.run_jules()
    return 'run completed'


class Jules():
    """Class to read in nml files and run JULES.

    :param nml_dir: location of JULES nml file directory.
    :type nml_dir: str
    :param jules_exe: location of JULES executable.
    :type jules_exe: str

    .. note:: You must have JULES installed on local system with a version of 4.8 or higher.

    """
    def __init__(self, nml_dir, jules_exe='/home/if910917/jules/models/jules4.8/build/bin/jules.exe', print_out=0):
        self.nml_dir = nml_dir
        self.jules = jules_exe
        self.nml_dic = self.read_nml(nml_dir)
        self.dir_path = os.path.dirname(os.path.realpath(__file__))
        self.example_nml = self.dir_path + '/' + 'data/jules/example_nml'
        if print_out == 0:
            self.run_jules = self.run_jules_pipe
        elif print_out == 1:
            self.run_jules = self.run_jules_print

    def read_nml(self, nml_dir):
        """
        Reading all nml files in specified directory and stores them in a dictionary.
        :param nml_dir: directory to read nml files from
        :return: dictionary of nml files
        """
        nml_dic = {}
        for f_nml in glob.glob(nml_dir + '/*.nml'):
            nml_dic[f_nml.split('/')[-1][:-4]] = f90nml.read(f_nml)
        return nml_dic

    def write_nml(self):
        """
        Function to write dictionary of stored nml data to nml files in current working dir.
        :return: n/a
        """
        for key in self.nml_dic.keys():
            self.nml_dic[key].write(key + '.nml', force=True)

    def run_jules_pipe(self):
        """Write all NML files to disk. Run JULES in a subprocess. Check output for fatal errors.

        :return: stdout and stderr output from JULES model run.
        :rtype: tuple
        """

        # write all the nml files here so the
        # user doesn't have to remember to...

        # run JULES
        cwd = os.getcwd()
        os.chdir(self.nml_dir)
        self.write_nml()
        cmd = []
        cmd.append(self.jules)

        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out = p.stdout.readlines()
        err = p.stderr.readlines()
        p.wait()
        os.chdir(cwd)

        # catch "fatal" errors
        for line in out:
            if len(line.split()) == 0: continue
            if line.split()[0] == "[FATAL":
                print >> sys.stderr, "*** runJules: caught fatal error in JULES run:"
                print >> sys.stderr, line,
                sys.exit()

        # catch anything in stderr
        if len(err) > 0:
            for line in err:
                print >> sys.stderr, "*** runJules: caught output on stderr in JULES run:"
                print >> sys.stderr, line,
                # sys.exit()
        return (out, err)

    def run_jules_print(self):
        """Write all NML files to disk. Run JULES in a subprocess. Check output for fatal errors.

        :return: stdout and stderr output from JULES model run.
        :rtype: str
        """

        # write all the nml files here so the
        # user doesn't have to remember to...
        cwd = os.getcwd()
        os.chdir(self.nml_dir)
        self.write_nml()

        # run JULES
        cmd = []
        cmd.append(self.jules)
        proc = subprocess.Popen(cmd, shell=False)
        proc.communicate()
        os.chdir(cwd)
        return 'Done', 'Done'