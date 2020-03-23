Input satellite states and access description
**********************************************

.. _poi_file_description:

Point-Of-interest (POI) file description
============================================

First row contains the following header elements, and the following rows contain the data corresponding to these headers.
                                 
Description of the header elements:

* :code:`POI`, Index of the point-of-interest.
* :code:`Lat [deg]`, :code:`Lon [deg]`, latitude, longitude of the point-of-interest.  


Example 
-------

.. csv-table:: POI example CSV table
   :widths: 5,10,10
   :header-rows: 1
    
    POI,lat[deg],lon[deg]
    0,0.000000000,90.000000000
    1,0.000000000,-90.000000000
    2,0.000000000,76.153846154
    3,60.000000000,76.153846154

.. note:: Make sure the header titles are as specified above, and the delimiter are commas.
 

.. _access_info_CSV_file_description:

Access Info CSV file description
==================================

The first three rows contain general information. Fourth row contains the epoch in Julian Date UT1. The field name is :code:`Epoch[JDUT1]`. See following `link <https://aa.usno.navy.mil/data/docs/JulianDate.php>`_ for
reference about Julian Dates. 

First four rows convey general information:
For example:

.. code-block:: python

    Satellite states are in Earth-Centered-Inertial equatorial plane.,,,,,,,,,,
    Epoch[JDUT1] is 2458562.344745371,,,,,,,,,,
    All time is referenced to the Epoch.
    Mission Duration [Days] is 1

First and third rows are same always as above.  Second row tells the epoch (in Julian Day UT1) to which the time-series 
n the file is referenced. Fourth row tells in entire Mission Duration (for which the access is calculated) in # Days. 
The fifth row  contains the following header elements, and the following rows contain the data corresponding to these headers.

Description of the header elements:

* :code:`AccessFrom[Days]`,  The time at which access starts in [days], referenced to epoch in row-4.
* :code:`Duration[s]`, Duration of access in [s].
* :code:`POI` indicating index of ground-point.
* :code:`EventNum` indicating index of event.
* :code:`Time[Days]`, Time in [Days] at which the alongside satellite-state is recorded. Referenced to the epoch specified in row-4.
* :code:`X[km]`, :code:`Y[km]` :code:`Z[km]`, cartesian spatial coordinates of satellite in Earth Centered Inertial frame with equatorial plane.
* :code:`VX[km/s]`, :code:`VY[km/s]`, :code:`VZ[km/s]`, velocity of spacecraft in Earth Centered Inertial frame with equatorial plane.

Example 
-------

.. csv-table:: Access Info example CSV table
   :header-rows: 5
   :widths: 10,10,10,10,10,10,10,10,10,10,10
    
    Satellite states are in Earth-Centered-Inertial equatorial plane.,,,,,,,,,,
    Epoch[JDUT1] is 2458562.344745371,,,,,,,,,,
    All time is referenced to the Epoch.,,,,,,,,,,
    Mission Duration [Days] is 1,,,,,,,,,,
    EventNum,POI,AccessFrom[Days],Duration[s],Time[Days],X[km],Y[km],Z[km],VX[km/s],VY[km/s],VZ[km/s]
    0,2,0.977120574,0.014081597,0.977120656,-702.601494297,1520.830748185,6879.150844801,6.544122557,-3.386908975,1.417249979
    1,6,0.019225555,0.013075769,0.019225625,1651.089340020,194.042016986,6882.343505376,6.288180699,-3.847571465,-1.399650576
    2,7,0.152406946,0.013075769,0.152407016,-728.062641367,1555.830560106,6868.665172092,6.489713921,-3.466789262,1.473203059

.. note:: Make sure the header titles are as specified above, and the delimiter are commas. 