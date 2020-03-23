""" 
.. module:: basic_sensor

:synopsis: *Module to handle sensor type with the most basic attributes.*

"""

import json
import numpy
import copy
import pandas, csv
from .util import Entity, Orientation, FieldOfView, MathUtilityFunctions, Constants

class BasicSensor(Entity):
    """A basic sensor class. 
      
        :ivar name: Full name of the instrument.
        :vartype name: str
        
        :ivar acronym: Acronym, initialism, or abbreviation.
        :vartype acronym: str

        :ivar mass: Total mass (kg) of this entity.
        :vartype mass: float

        :ivar volume: Total volume (m3) of this entity.
        :vartype volume: float
        
        :ivar power: Nominal operating power (W).
        :vartype power: float

        :ivar orientation: Orientation of the instrument with respect to Nadir-frame.
        :vartype orientation: :class:`instrupy.util.Orientation`

        :ivar fieldOfView: Field of view specification of instrument. 
        :vartype fieldOfView: :class:`instrupy.util.FieldOfView`   
       
        :ivar dataRate: Rate of data recorded (Mbps) during nominal operations.
        :vartype dataRate: float  

        :ivar bitsPerPixel: Bits encoded per pixel of image
        :vartype bitsPerPixel: int
   
    """

    def __init__(self, name=None, acronym=None, mass=None,
            volume=None, power=None,  orientation=None,
            fieldOfView=None, dataRate=None, bitsPerPixel = None, _id=None):
        """Initialize a Basic Sensor object.

        """
        self.name = str(name) if name is not None else None
        self.acronym = str(acronym) if acronym is not None else self.name
        self.mass = float(mass) if mass is not None else None
        self.volume = float(volume) if volume is not None else None
        self.power = float(power) if power is not None else None
        self.orientation = copy.deepcopy(orientation) if orientation is not None else None
        self.fieldOfView = copy.deepcopy(fieldOfView) if fieldOfView is not None else None
        self.dataRate = float(dataRate) if dataRate is not None else None
        self.bitsPerPixel = int(bitsPerPixel) if bitsPerPixel is not None else None            
        super(BasicSensor,self).__init__(_id, "Basic Sensor")

    @staticmethod
    def from_dict(d):
        """Parses an instrument from a normalized JSON dictionary."""
        return BasicSensor(
                name = d.get("name", None),
                acronym = d.get("acronym", None),
                mass = d.get("mass", None),
                volume = d.get("volume", None),
                power = d.get("power", None),
                orientation = Orientation.from_json(d.get("orientation", None)),
                fieldOfView = FieldOfView.from_json(d.get("fieldOfView", None)),
                dataRate = d.get("dataRate", None),
                bitsPerPixel = d.get("bitsPerPixel", None),
                _id = d.get("@id", None)
            )

    def calc_typ_data_metrics_over_one_access_interval(self, SpacecraftOrbitState, AccessInfo):
        ''' Calculate typical observation data metrics during the given access period.

            :param SpacecraftOrbitState: Spacecraft position at the middle (or as close as possible to the middle) of the access interval.

                               Dictionary keys are: 
                               
                               * :code:`Time[JDUT1]` (:class:`float`), Time in Julian Day UT1.
                               * :code:`x[km]` (:class:`float`), :code:`y[km]` (:class:`float`), :code:`z[km]` (:class:`float`), cartesian spatial coordinates of satellite in Earth Centered Inertial frame with equatorial plane.
                               * :code:`vx[km/s]` (:class:`float`), :code:`vy[km/s]` (:class:`float`), :code:`vz[km/s]` (:class:`float`), velocity of spacecraft in Earth Centered Inertial frame with equatorial plane.
            :paramtype SpacecraftOrbitState: dict

            
            :param AccessInfo: Access information.

                               Dictionary keys are: 
                                
                               * :code:`Access From [JDUT1]` (:class:`float`) Start absolute time of Access in Julian Day UT1.
                               * :code:`Duration [s]` (:class:`float`): Access duration in [s].
                               * :code:`Lat [deg]` (:class:`float`), :code:`Lon [deg]` (:class:`float`), indicating the corresponding ground-point accessed (latitude, longitude) in degrees.
            :paramtype AccessInfo: dict

            :returns: Typical calculated observation data metrics. It is assumed that the satellite takes *one* observation data sample per access event. Below metrics are 
                      calculated using satellite position data at or near the middle of the access event.
                      
                      Dictionary keys are: 
                    
                      * :code:`Coverage [T/F]` (:class:`bool`) indicating if observation was possible during the access event.
                      * :code:`Incidence angle [deg]` (:class:`float`) Incidence angle at target point calculated assuming spherical Earth.
                      * :code:`Look angle [deg]` (:class:`float`) Look angle at target point calculated assuming spherical Earth.
                      * :code:`Observation Range [km]` (:class:`float`) Distance from satellite to ground-point during the observation acquisition.
                      * :code:`Solar Zenith [deg]` (:class:`float`) Solar Zenith during observation

                      .. note:: Metrics are calculated at the instant corresponding to the supplied Satellite orbit state in the parameter :code:`SpacecraftOrbitState`.
            :rtype: dict


            .. note:: It is assumed that the instrument captures **one** *observation/ image* over the entire access interval. 
                      While this is true for stripmap radar imagers, push-broom imagers and whisk-broom imagers, it is not true for all
                      different types of instruments/ imaging modes.

            .. note:: We differentiate between **access** and **coverage**. **Access** is when the target location
                      falls under the sensor field-of-view. **Coverage** is when the target location falls under sensor field-of-view *and* can be observed.
                        
        '''
        # Observation time in Julian Day UT1
        tObs_JDUT1 = SpacecraftOrbitState["Time[JDUT1]"]

        # Calculate Target position in ECI frame
        TargetPosition_km = MathUtilityFunctions.geo2eci([AccessInfo["Lat [deg]"], AccessInfo["Lon [deg]"], 0.0], tObs_JDUT1)

        # Spacecraft position in Cartesian coordinates
        SpacecraftPosition_km = numpy.array([SpacecraftOrbitState["x[km]"], SpacecraftOrbitState["y[km]"], SpacecraftOrbitState["z[km]"]])  

        alt_km = numpy.linalg.norm(SpacecraftPosition_km) - Constants.radiusOfEarthInKM

        #  Calculate actual range between spacecraft and POI (Target)
        range_vector_km = TargetPosition_km - SpacecraftPosition_km

        range_km = numpy.linalg.norm(range_vector_km)

        # Calculate look angle
        look_angle = numpy.arccos(numpy.dot(MathUtilityFunctions.normalize(range_vector_km), -1*MathUtilityFunctions.normalize(SpacecraftPosition_km)))
        look_angle_deg = numpy.rad2deg(look_angle)
        
        # Look angle to corresponding incidence angle conversion for spherical Earth
        incidence_angle = numpy.arcsin(numpy.sin(look_angle)*(Constants.radiusOfEarthInKM +alt_km)/Constants.radiusOfEarthInKM)
        incidence_angle_deg =  numpy.rad2deg(incidence_angle)

        # Solar zenith angle
        [solar_zenith, solar_distance] = MathUtilityFunctions.compute_sun_zenith(tObs_JDUT1, TargetPosition_km)
        if solar_zenith is not None:
            solar_zenith_deg =  numpy.rad2deg(solar_zenith)
        else:
            solar_zenith_deg = numpy.nan

        # For a basic-sensor, coverage is true for all access-events
        isCovered = True 
    
        obsv_metrics = {}
        obsv_metrics["Observation Range [km]"] = range_km
        obsv_metrics["Look angle [deg]"] = look_angle_deg
        obsv_metrics["Incidence angle [deg]"] = incidence_angle_deg
        obsv_metrics["Solar Zenith [deg]"] = solar_zenith_deg
        obsv_metrics["Coverage [T/F]"] = isCovered

        return obsv_metrics

