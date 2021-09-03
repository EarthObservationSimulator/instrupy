Tests
*******

Testing is done using the :mod:`unittest` module. Use the below command to run all the tests:

.. code-block:: shell

        make runtest

Testing strategies
^^^^^^^^^^^^^^^^^^^^

1. Test the expected format of the output files.
2. Use of random input values and calculated output values to validate the test.
       
        i.  Sometimes special cases are run, and some input values which do not make a difference to the 
            output are made random. Eg: In case of testing the sign of look-angles, only the sign of the angle at which
            the sensor is tilted (with respect to the nadir-pointing frame) is important, and the magnitude can be used as an random input. 
        ii. In cases where the random input values do influence the output, the expected output is calculated 
            (as much as possible) from methods, code other than used by the :code:`InstruPy` package. Sometimes
            the chosen validation methods are approximate in which case an *approximately equal to* assertion tests
            are used.
3. Using known inputs, and outputs from external sources (eg: literature).
4. Use results from a previous run (corresponding to an older version of the software) as truth data.
5. Run tests after making any revisions to the code. This helps to check that the revisions do not have unintended effects on the results.

Test Modules
============

:class:`test_base`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: test_base
   :members:

:class:`test_basic_sensor_model`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: test_basic_sensor_model
   :members:

:class:`test_passive_optical_scanner_model`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: test_passive_optical_scanner_model
   :members:

:class:`test_synthetic_aperture_radar_model`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: test_synthetic_aperture_radar_model
   :members:

:class:`test_radiometer_model`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: test_radiometer_model
   :members:

:class:`test_util`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. automodule:: test_util
   :members: