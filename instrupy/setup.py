from setuptools import setup

setup(
    name='InstruPy',
    version='0.1.0',
    description='Instrument module',
    author='BAERI',
    author_email='vinay.ravindra@nasa.gov',
    packages=['instrupy'],
    scripts=[
        'bin/InstruPy_get_coverage_specs.py',
        'bin/instrument_module.py'
    ],
    install_requires=['isodate', 'numpy', 'pandas', 'scipy']
)
