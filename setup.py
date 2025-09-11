from setuptools import setup

def readme():
    with open('README.md') as f:
        return f.read()

setup(
    name='InstruPy',
    version='0.2',
    description='Instrument module',
    author='BAERI',
    author_email='vravindra@baeri.org',
    packages=['instrupy'],
    scripts=[ 
    ],
    install_requires=['shapely', 'numpy', 'pandas', 'scipy', 'lowtran==2.4.1', 'sphinx', 'sphinx_rtd_theme==2.0.0', 'netCDF4','metpy','astropy','deepdiff'] 
)
