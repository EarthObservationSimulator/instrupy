``tests`` --- Tests
********************

Testing is done using the :mod:`unittest` framework available in the standard distribution of Python. To run using :mod:`nose` module, run from the command line:

.. code-block:: shell

        python -m nose


The tests are classified into three types:

1.  Low-level tests: To test the smallest unit functions and/or classes as available in the API.
2.  Mid-level tests: Used to tests two to three interacting functions and/or classes *together*.
3.  High-level tests: End-to-end tests.


Directory structure:
::

    tests/
    ├── low_level/
    |    ├── test_util
    |    ├── test_public_library
    |    ├── ...
    ├── mid_level/
    |    ├── ...
    |    ├── ...
    ├── high_level/
    |    ├── ...
    |    ├── ...
