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
    scripts=[ 
    ],
    install_requires=['numpy', 'pandas', 'scipy', 'lowtran', 'nose', 'sphinx', 'sphinx_rtd_theme'] 
)
