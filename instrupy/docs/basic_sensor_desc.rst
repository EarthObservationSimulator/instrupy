Basic Sensor Description
************************

Types of Operating Modes
==========================
There are no special types of operating modes (only a single mode type). 

Input JSON format specifications description
=========================================================

.. csv-table:: Input parameter description 
    :header: Parameter, Data type, Units, Description
    :widths: 10,10,5,40

    @type, string, , Must be *Basic Sensor*
    @id, string, , Unique identifier for the instrument.
    name, string, ,Full name of the instrument 
    acronym, string, ,Acronym or initialism or abbreviation.
    mass, float, kilograms, Total mass of this entity.
    volume, float, :math:`m^3`, Total volume of this entity.
    power, float, Watts, Nominal operating power.
    orientation, :ref:`orientation_json_obj`, ,Orientation of the instrument with respect to Nadir-frame. 
    fieldOfView, :ref:`fieldOfView_json_obj`, ,Field of view specification of instrument. 
    dataRate, float, Mega-bits-per-s,Rate of data recorded during nominal operations.
    bitsPerPixel, integer, ,Bits encoded per pixel of image.
    maneuverability, :ref:`maneuverability_json_object`, ,Payload maneuverability (see :ref:`manuv_desc`)

.. _basic_sensor_data_metrics_calc:

Output observation data metrics calculation
=============================================

.. csv-table:: Observation data metrics table
    :widths: 8,4,4,20
    :header: Metric,Data Type,Units,Description 
     
    Coverage [T/F], string ,, Indicates if observation was possible during the access event *(Always True)*. 
    Incidence angle [deg], float,  degrees, Incidence-angle at target point calculated assuming spherical Earth.
    Look angle [deg], float,  degrees, Look-angle to the target point calculated assuming spherical Earth.
    Observation Range [km], float, kilometers, Distance from satellite to ground-point during the observation. 
    Solar Zenith [deg], float, degrees, Solar-zenith-angle during observation

Viewing geometry
-----------------

See :ref:`satellite_to_target_viewing_geometry` for the calculation of the viewing geometry parameters.





