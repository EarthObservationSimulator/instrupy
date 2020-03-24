from setuptools import setup

def readme():
    with open('README.md') as f:
        return f.read()

setup(
    name='InstruPy',
    version='0.2',
    description='Instrument module',
    author='BAERI',
    author_email='vinay.ravindra@nasa.gov',
    packages=['instrupy'],
    scripts=[ # TODO: remove this? Does not seem to serve any purpose. 
        'bin/InstruPy_get_coverage_specs.py',
        'bin/instrument_module.py'
    ],
    install_requires=['numpy', 'pandas', 'scipy', 'lowtran'] # TODO: Specify lowtran as dependency?
)
