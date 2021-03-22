Basic Sensor Model Description
********************************

Input JSON format specifications description
=========================================================

.. csv-table:: Input parameter description 
    :header: Parameter, Data type, Units, Description
    :widths: 10,10,5,40

    @type, string, ,Must be *Basic Sensor*
    @id, string, , Unique identifier for the instrument.
    name, string, ,Full name of the instrument 
    mass, float, kilograms, Total mass of this entity.
    volume, float, :math:`m^3`, Total volume of this entity.
    power, float, Watts, Nominal operating power.
    orientation, :ref:`orientation_json_obj`, ,Orientation of the instrument.
    fieldOfView, :ref:`fieldOfViewGeometry_json_obj`, , Field of view spherical geometry specification of the instrument.
    maneuver, :ref:`maneuver_json_object`, , Maneuver specifications (see :ref:`maneuv_desc`).
    pointingOption, :ref:`pointing_opt_json_obj`, , List of orientations to which the instrument axis can be maneuvered.
    syntheticDataConfig, :ref:`syntheticDataConfig_json_obj`, , Synthetic data configuration.
    dataRate, float, Mega-bits-per-s, Rate of data recorded during nominal operations.
    bitsPerPixel, integer, ,Bits encoded per pixel of image.
    numberDetectorRows, integer, ,Number of detector rows (along the Y-axis of the SENOR_BODY_FIXED frame). 
    numberDetectorCols, integer, ,Number of detector columns (along the X-axis of the SENOR_BODY_FIXED frame). 
    
.. _basic_sensor_data_metrics_calc:

Output observation data metrics calculation
=============================================

The following data metrics at a target location (also referred to as a grid-point, point-of-interest) on the surface of Earth can be calculated using the basic sensor model.

.. csv-table:: Observation data metrics table
    :widths: 8,4,4,20
    :header: Metric,Data Type,Units,Description 
     
    Coverage [T/F], string,, Indicates if observation was possible during the access event *(Always True)*. 
    Incidence angle [deg], float,  degrees, Incidence-angle at target point calculated assuming spherical Earth.
    Look angle [deg], float,  degrees, Look-angle to the target point calculated assuming spherical Earth. Positive sign => look is in positive half-space made by the orbit-plane (i.e. orbit plane normal vector) and vice-versa.
    Observation Range [km], float, kilometers, Distance from satellite to ground-point during the observation. 
    Solar Zenith [deg], float, degrees, Solar-zenith-angle during observation

Viewing geometry
------------------

See :ref:`satellite_to_target_viewing_geometry` for the calculation of the observational data-metrics of the basic sensor model.





