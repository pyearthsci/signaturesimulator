"""A setuptools based setup module.
See:
https://packaging.python.org/en/latest/distributing.html
https://github.com/pypa/sampleproject
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path
import subprocess
from setuptools.command.install import install

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

class MyInstall(install):
    def run(self):
        try:
            # note cwd - this makes the current directory
            # the one with the Makefile.
            print path.join(here,'signaturesimulator/semidiscrete_srf/')
            subprocess.call(['make'], cwd=path.join(here, 'signaturesimulator/semidiscrete_srf/'))
        except Exception as e:
            print e
            print "Error compiling semiD.c.   Try running 'make'."
            exit(1)
        else:
            install.run(self)

# Arguments marked as "Required" below must be included for upload to PyPI.
# Fields marked as "Optional" may be commented out.

setup(
    # This is the name of your project. The first time you publish this
    # package, this name will be registered for you. It will determine how
    # users can install this project, e.g.:
    #
    # $ pip install sampleproject
    #
    # And where it will live on PyPI: https://pypi.org/project/sampleproject/
    #
    # There are some restrictions on what makes a valid project name
    # specification here:
    # https://packaging.python.org/specifications/core-metadata/#name
    name='signaturesimulator',  # Required

    # Versions should comply with PEP 440:
    # https://www.python.org/dev/peps/pep-0440/
    #
    # For a discussion on single-sourcing the version across setup.py and the
    # project code, see
    # https://packaging.python.org/en/latest/single_source_version.html
    version='0.2.4',  # Required

    # This is a one-line description or tagline of what your project does. This
    # corresponds to the "Summary" metadata field:
    # https://packaging.python.org/specifications/core-metadata/#summary
    description='A Python package that simulates satellite data from passive optical and active microwave sensors '
                'based on land surface biogeophysical state variables (leaf area index, canopy height and soil '
                'moisture) and viewing geometries and times.',  # Required

    # This is an optional longer description of your project that represents
    # the body of text which users will see when they visit PyPI.
    #
    # Often, this is the same as your README, so you can just read it in from
    # that file directly (as we have already done above)
    #
    # This field corresponds to the "Description" metadata field:
    # https://packaging.python.org/specifications/core-metadata/#description-optional
    long_description=long_description,  # Optional

    # This should be a valid link to your project's main homepage.
    #
    # This field corresponds to the "Home-Page" metadata field:
    # https://packaging.python.org/specifications/core-metadata/#home-page-optional
    url='https://github.com/pyearthsci/signaturesimulator',  # Optional

    # This should be your name or the name of the organization which owns the
    # project.
    author='T. Quaife and E. Pinnington',  # Optional

    # This should be a valid email address corresponding to the author listed
    # above.
    author_email='e.pinnington@reading.ac.uk',  # Optional

    # This field adds keywords for your project which will appear on the
    # project page. What does your project relate to?
    #
    # Note that this is a string of words separated by whitespace, not a list.
    keywords='satellite earth science land surface passive optical active microwave sentinel',  # Optional

    # You can just specify package directories manually here if your project is
    # simple. Or you can use find_packages().
    #
    # Alternatively, if you just want to distribute a single Python file, use
    # the `py_modules` argument instead as follows, which will expect a file
    # called `my_module.py` to exist:
    #
    #   py_modules=["my_module"],
    #
    packages=['signaturesimulator'],  # Required

    # This field lists other packages that your project depends on to run.
    # Any package you put here will be installed by pip when your project is
    # installed, so they must be valid existing projects.
    #
    # For an analysis of "install_requires" vs pip's requirements files see:
    # https://packaging.python.org/en/latest/requirements.html
    install_requires=[#'pandas',
                      'numpy',
                      'scipy',
                      'pyorbital',
                      'mock',
                      'matplotlib',
                      'netCDF4',
                      'f90nml'],  # Optional

    # If there are data files included in your packages that need to be
    # installed, specify them here.
    #
    # If using Python 2.6 or earlier, then these have to be included in
    # MANIFEST.in as well.
    package_data={  # Optional
        'signaturesimulator': ['site.nml', 'data/tle/*', 'data/srf/*', 'data/state_variables/*', 'data/geometries/*',
                               'sense/*.py', 'sense/dielectric/*', 'sense/surface/*', 'sense/vegetation/*',
                               'semidiscrete_srf/*', 'semidiscrete_srf/man/*', 'data/jules/*',
                               'data/jules/example_nml/*', 'data/jules/jules_data/*', 'data/jules/output/*']
    },
    cmdclass={'install': MyInstall},
)
