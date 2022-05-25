"""
:synopsis: *Module to handle sensor model which is used for reflectometry.*

"""
import copy
from instrupy.util import Entity, Orientation, SphericalGeometry, Maneuver, Antenna, ViewGeometry, MathUtilityFunctions, GeoUtilityFunctions, Constants, SyntheticDataConfiguration
import numpy as np
import uuid
from collections import namedtuple
import logging
logger = logging.getLogger(__name__)

class ReflectometerModel(Entity):
    """A reflectometer model class. 
      
        :ivar name: Full name of the instrument.
        :vartype name: str

        :ivar mass: Total mass (kg) of this instrument.
        :vartype mass: float

        :ivar volume: Total volume (m3) of this instrument.
        :vartype volume: float
        
        :ivar power: Nominal operating power (W) of this instrument.
        :vartype power: float

        :ivar orientation: Orientation of the instrument.
        :vartype orientation: :class:`instrupy.util.Orientation`

        :ivar fieldOfView: Field of view of instrument specification (SphericalGeometry and Orientation).
        :vartype fieldOfView: :class:`instrupy.util.ViewGeometry`

        :ivar sceneFieldOfView: Scene field of view specification (SphericalGeometry and Orientation).
        :vartype fieldOfView: :class:`instrupy.util.ViewGeometry`

        :ivar maneuver: Maneuver specification of the instrument. 
        :vartype maneuver: :class:`instrupy.util.Maneuver`  

        :ivar fieldOfRegard: Field of regard of the instrument taking into account the sceneFOV and manueverability of the satellite-sensor system. 
                             Note that this shall be a list in order to accommodate non-intersecting view-geometries.
                             TODO: Modify behavior to have FOR = sceneFOV when no maneuver is specified (hence fixed pointing). Currently FOR = None if manuever is not specified.
        :vartype fieldOfRegard: list, :class:`instrupy.util.ViewGeometry`  
       
        :ivar pointingOption: List of ``Orientation`` objects which specify the orientations into which the instrument-axis can be maneuvered. 
                               The orientations must be specified in the NADIR_POINTING frame.
        :vartype pointingOption: list, :class:`orbitpy.util.Orientation`

        :ivar antenna: Antenna specifications. Only circular shape and uniform aperture excitation profile is accepted.
        :vartype antenna: :class:`instrupy.util.Antenna`

        :ivar operatingFrequency: Operating antenna center frequency in (Hertz).
        :vartype operatingFrequency: float

        :ivar dataRate: Rate of data recorded (Mega bits per sec) during nominal operations.
        :vartype dataRate: float

        :ivar bitsPerPixel: Number of bits encoded per pixel of image.
        :vartype bitsPerPixel: int

        :ivar _id: Unique instrument identifier.
        :vartype _id: str or int
   
    """
    def __init__(self, name=None, mass=None, volume=None, power=None,  orientation=None,
                 fieldOfViewGeometry=None, sceneFieldOfViewGeometry=None, maneuver=None, pointingOption=None, 
                 antenna=None, operatingFrequency=None, dataRate=None, bitsPerPixel=None, 
                 _id=None):
        """ Initialization. All except the below two parameters have identical description as that of the corresponding class instance variables.

        :param fieldOfViewGeometry: Instrument field-of-view spherical geometry.
        :paramtype fieldOfViewGeometry: :class:`instrupy.util.SphericalGeometry`

        :param sceneFieldOfViewGeometry: Scene field-of-view spherical geometry.
        :paramtype sceneFieldOfViewGeometry: :class:`instrupy.util.SphericalGeometry`

        :return: None
        :rtype: None

        """
        self.name = str(name) if name is not None else None
        self.mass = float(mass) if mass is not None else None
        self.volume = float(volume) if volume is not None else None
        self.power = float(power) if power is not None else None
        self.orientation = copy.deepcopy(orientation) if orientation is not None and isinstance(orientation, Orientation) else None
        self.fieldOfView = ViewGeometry(orien=self.orientation, sph_geom=fieldOfViewGeometry) if self.orientation is not None and fieldOfViewGeometry is not None and isinstance(fieldOfViewGeometry, SphericalGeometry) else None
        self.sceneFieldOfView = ViewGeometry(orien=self.orientation, sph_geom=sceneFieldOfViewGeometry) if self.orientation is not None and sceneFieldOfViewGeometry is not None and isinstance(sceneFieldOfViewGeometry, SphericalGeometry) else None
        self.maneuver = copy.deepcopy(maneuver) if maneuver is not None and isinstance(maneuver, Maneuver) else None
        self.fieldOfRegard = self.maneuver.calc_field_of_regard(self.sceneFieldOfView.sph_geom) if (self.maneuver is not None and self.sceneFieldOfView is not None) else None
        self.pointingOption = None
        if isinstance(pointingOption, list):
            if all(isinstance(x, Orientation) for x in pointingOption):
                self.pointingOption = pointingOption 
        elif isinstance(pointingOption, Orientation):
            self.pointingOption = [pointingOption] # make into single-element list
        self.antenna = copy.deepcopy(antenna) if antenna is not None and isinstance(antenna, Antenna) else None
        self.operatingFrequency = float(operatingFrequency) if operatingFrequency is not None else None
        self.dataRate = float(dataRate) if dataRate is not None else None
        self.bitsPerPixel = int(bitsPerPixel) if bitsPerPixel is not None else None                

        self._id = _id if _id is not None else None

        super(ReflectometerModel,self).__init__(_id, "Reflectometer")

    @staticmethod
    def from_dict(d):
        """ Parses a ``ReflectometerModel`` object from a normalized JSON dictionary. 
            Refer to :ref:`reflectometer_model_desc` for description of the accepted key/value pairs.
        
        The following default values are assigned to the object instance parameters in case of 
        :class:`None` values or missing key/value pairs in the input dictionary.

        .. csv-table:: Default values
            :header: Parameter, Default Value
            :widths: 10,40

            orientation, Referenced and aligned to the SC_BODY_FIXED frame.
            fieldOfViewGeometry, CIRCULAR shape with 25 deg diameter.
            sceneFieldOfViewGeometry, fieldOfViewGeometry
            numberDetectorRows, 4
            numberDetectorRows, 4
            _id, random string
        
        :param d: Normalized JSON dictionary with the corresponding model specifications. 
        :paramtype d: dict

        :returns: ``ReflectometerModel`` object initialized with the input specifications.
        :rtype: :class:`instrupy.basic_sensor_model.ReflectometerModel`
            
        """
        orien_dict = d.get("orientation", {"referenceFrame": "SC_BODY_FIXED", "convention": "REF_FRAME_ALIGNED"}) #  default orientation = referenced and aligned to the SC_BODY_FIXED frame.
        
        # parse the pointing options as a list of Orientation objects.
        pnt_opt_dict = d.get("pointingOption", None)
        _pointing_option = None
        if pnt_opt_dict:
            # translate to a list of Orientation objects
            if isinstance(pnt_opt_dict, list):
                _pointing_option = [Orientation.from_dict(x) for x in pnt_opt_dict]
            else:
                _pointing_option = [Orientation.from_dict(pnt_opt_dict)]
        
        antenna_dict = d.get("antenna", None)
        if antenna_dict:
                antenna = Antenna.from_dict(antenna_dict)
                antenna_fov_sph_geom = antenna.get_spherical_geometry(d.get("operatingFrequency", None))
                [dia, _] = antenna_fov_sph_geom.get_fov_height_and_width() # only circular diameter sensor FOV specs is supported.
                instru_fov_geom_dict = {"shape": "CIRCULAR", "diameter":dia} 
        else:
            antenna = None
            instru_fov_geom_dict = None

        scene_fov_geom_dict = d.get("sceneFieldOfViewGeometry", instru_fov_geom_dict)  # default sceneFOV geometry is the instrument FOV geometry

        return ReflectometerModel(
                name = d.get("name", None),
                mass = d.get("mass", None),
                volume = d.get("volume", None),
                power = d.get("power", None),
                orientation = Orientation.from_dict(orien_dict),
                fieldOfViewGeometry =  SphericalGeometry.from_json(instru_fov_geom_dict),
                sceneFieldOfViewGeometry = SphericalGeometry.from_json(scene_fov_geom_dict),
                maneuver =  Maneuver.from_json(d.get("maneuver", None)),
                dataRate = d.get("dataRate", None),
                bitsPerPixel = d.get("bitsPerPixel", None),
                pointingOption = _pointing_option,
                antenna = antenna,
                operatingFrequency = d.get("operatingFrequency", None),
                _id = d.get("@id", str(uuid.uuid4()))
                )

    def to_dict(self):
        """ Translate the ``ReflectometerModel`` object to a Python dictionary such that it can be uniquely reconstructed back from the dictionary.

        :returns: ``ReflectometerModel`` specifications as python dictionary.
        :rtype: dict

        """
        fieldOfViewGeometry_dict = self.fieldOfView.sph_geom.to_dict() if self.fieldOfView is not None and isinstance(self.fieldOfView, ViewGeometry) else None
        sceneFieldOfViewGeometry_dict = self.sceneFieldOfView.sph_geom.to_dict() if self.sceneFieldOfView is not None and isinstance(self.sceneFieldOfView, ViewGeometry) else None
        orientation_dict = self.orientation.to_dict() if self.orientation is not None and isinstance(self.orientation, Orientation) else None
        maneuver_dict = self.maneuver.to_dict() if self.maneuver is not None and isinstance(self.maneuver, Maneuver) else None
        pointing_opt_dict = [Orientation.to_dict(x) for x in self.pointingOption] if self.pointingOption is not None else None
        antenna_dict = self.antenna.to_dict() if self.antenna is not None and isinstance(self.antenna, Antenna) else None
        return dict({
                "@type": "Reflectometer",
                "name":self.name,
                "mass":self.mass,
                "volume":self.volume,
                "power":self.power,
                "fieldOfViewGeometry":fieldOfViewGeometry_dict,
                "sceneFieldOfViewGeometry": sceneFieldOfViewGeometry_dict,
                "orientation":orientation_dict,
                "maneuver":maneuver_dict,
                "pointingOption": pointing_opt_dict,
                "dataRate":self.dataRate,
                "bitsPerPixel": self.bitsPerPixel,
                "antenna": antenna_dict,
                "operatingFrequency": self.operatingFrequency,
                "@id": self._id
                })

    def __repr__(self):
        return "ReflectometerModel.from_dict({})".format(self.to_dict())

    def calc_data_metrics(self, sc_orbit_state, target_coords, tx_spc_orbit_state):
        """ ======= TBD =======
        
        Calculate typical observation data metrics. Refer to :ref:`basic_sensor_model_desc` for description.            

        :param sc_orbit_state: Spacecraft state at the time of observation.

                            Dictionary keys are: 
                            
                            * :code:`time [JDUT1]` (:class:`float`), Time in Julian Day UT1. Corresponds to the time of observation. 
                            * :code:`x [km]` (:class:`float`), :code:`y [km]` (:class:`float`), :code:`z [km]` (:class:`float`), Cartesian coordinates of satellite in EARTH_CENTERED_INERTIAL frame at the time of observation.
                            * :code:`vx [km/s]` (:class:`float`), :code:`vy [km/s]` (:class:`float`), :code:`vz [km/s]` (:class:`float`), Velocity of spacecraft in EARTH_CENTERED_INERTIAL frame at the time of observation.
        
        :paramtype sc_orbit_state: dict

        
        :param target_coords: Location of the observation. Also sometimes the Point-Of-Interest (POI).

                            Dictionary keys are: 
                            
                            * :code:`lat [deg]` (:class:`float`), :code:`lon [deg]` (:class:`float`) indicating the corresponding ground-point accessed (latitude, longitude) in degrees.
        
        :paramtype target_coords: dict

        :returns: Calculated observation data metrics.
                    
                    Dictionary keys are: 
                
                    * :code:`incidence angle [deg]` (:class:`float`) Incidence angle in degrees at the target point calculated assuming spherical Earth.
                    * :code:`look angle [deg]` (:class:`float`) Look angle in degrees at the target point calculated assuming spherical Earth. Positive sign => look is in positive half-space made by the orbit-plane (i.e. orbit plane normal vector) and vice-versa.
                    * :code:`observation range [km]` (:class:`float`) Distance in kilometers from satellite to ground-point during the observation.
                    * :code:`solar zenith [deg]` (:class:`float`) Solar Zenith angle in degrees during observation.
                
                    .. todo:: Include AT, CT footprint size calculations.

        :rtype: dict                      

        """
        # Observation time in Julian Day UT1
        tObs_JDUT1 = sc_orbit_state["time [JDUT1]"]

        # Calculate Target cartesian position in EARTH_CENTERED_INERTIAL frame
        target_pos = GeoUtilityFunctions.geo2eci([target_coords["lat [deg]"], target_coords["lon [deg]"], 0.0], tObs_JDUT1)

        # Spacecraft position in Cartesian coordinates in the EARTH_CENTERED_INERTIAL frame
        sc_pos = np.array([sc_orbit_state["x [km]"], sc_orbit_state["y [km]"], sc_orbit_state["z [km]"]])  
        sc_vel = np.array([sc_orbit_state["vx [km/s]"], sc_orbit_state["vy [km/s]"], sc_orbit_state["vz [km/s]"]])  

        alt_km = np.linalg.norm(sc_pos) - Constants.radiusOfEarthInKM # altitude
