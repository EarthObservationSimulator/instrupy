Basic Sensor Description
************************

Input JSON format specifications description
=========================================================

.. csv-table:: Input parameter description 
    :header: Parameter, Data type, Units, Description
    :widths: 10,10,5,40

    name, string, ,Full name of the instrument 
    acronym, string, ,Acronym or initialism or abbreviation.
    mass, float, kilograms, Total mass of this entity.
    volume, float, :math:`m^3`, Total volume of this entity.
    power, float, Watts, Nominal operating power.
    orientation, :ref:`orientation_json_string`, ,Orientation of the instrument with respect to Nadir-frame. 
    fieldOfView, :ref:`fieldOfView_json_string`, ,Field of view specification of instrument. 
    dataRate, float, Mega-bits-per-s,Rate of data recorded during nominal operations.
    bitsPerPixel, integer, ,Bits encoded per pixel of image.

.. _basic_sensor_csv_output:

Basic Sensor Level-0 CSV output file description
=================================================

Description of the header elements:

.. csv-table:: Level-0 output data-metrics description
    :widths: 8,4,4,20
    :header: Metric/Aux data,Data Type,Units,Description 
     
    :code:`Access From [JDUT1]`, float , Julian Date UT1, Access from time
    :code:`Access Duration [s]`, float, seconds, Duration of access
    :code:`POI index`, integer, , Index of point of interest
    :code:`Coverage [T/F]` , string ,   , Indicates if observation was  possible during the access event  (True/ False). 
    :code:`Incidence angle [deg]`  , float,  degrees  , Incidence angle at target point calculated assuming spherical Earth.
    :code:`Look angle [deg]`, float,  degrees , Look angle at target point calculated assuming spherical Earth.
    :code:`Observation Range [km]` , float, kilometers, Distance from satellite to ground-point during the observation acquisition. 
    :code:`Solar Zenith [deg]`, float, degrees, Solar Zenith during observation

Example 
-------

.. csv-table:: Basic Sensor typical data metrics example CSV output file
    :header: Access From [JDUT1],Access Duration [s],POI index,Observation Range [km],Look angle [deg],Incidence angle [deg],Solar Zenith [deg], Coverage [T/F]
    :widths: 10,10,10,10,10,10,10,10
    
    2458562.8678740812,0.014081597,3,2188.572929282959,62.646876593956996,80.39127235404261,63.67762735547993,True
    2458563.001060143,0.013075769,4,4820.959603738683,61.38895768323014,77.03910683668607,,True
    2458562.524382066,0.014081597,5,9254.51259835866,43.477067498737576,49.80135036934889,83.010563889425,True
    2458562.657561033,0.013075769,6,12506.798214191318,20.515922575438836,22.895084061309824,,True
    
.. _basic_sensor_data_metrics_calc:

Typical observation data metrics calculation
=============================================
See :ref:`satellite_to_target_viewing_geometry` for the calculation of the viewing geometry parameters.





