""" 
.. module:: util

:synopsis: *Utility classes and functions for the :class:`instrupy` package.*

"""
from __future__ import division 
import json
import numpy as np
import math
from enum import Enum
from numbers import Number
import scipy.constants
import string
import lowtran

class Entity(object): 
    """An entity is an abstract class to aggregate common functionality.
    
    :ivar _id: Unique identifier for this entity.
    :vartype _id: str
    
    :ivar _type: Class type description for this entity.
    :vartype _type: str

    """

    def __init__(self, _id=None, _type="Entity"):
        """Initialize an entity.

        :param _id: (default: None)
        :paramtype _id: str 

        :param _type: (default: "Entity")
        :paramtype _type: str

        """
        self._id = _id
        self._type = _type

    def to_dict(self):
        """Convert this entity to a JSON-formatted dictionary."""
        # extract and copy the python dictionary
        json_dict = dict(self.__dict__)

        def recursive_normalize(d):
            """Helper function to recursively remove null values and serialize
            unserializable objects from dictionary."""
            # if not a dictionary, return immediately
            if not isinstance(d, dict):
                if isinstance(d, Entity): return d.to_dict()
                else: return d
            # otherwise loop through each key/value pair
            for key, value in list(d.items()):
                # if value is None remove key
                if value is None:
                    del d[key]
                # else if non-seralizable object, manually serialize to json
                elif isinstance(value, Entity):
                    d[key] = value.to_dict()
                # else if list, recursively serialize each list element
                elif isinstance(value, list):
                    d[key] = map(recursive_normalize, value)
                # otherwise recursively call function
                else: recursive_normalize(value)
            return d
        recursive_normalize(json_dict)
        # translate special python to json keys: _id to @id, _type to @type
        if json_dict.get("_id"): json_dict["@id"] = json_dict.pop("_id")
        if json_dict.get("_type"): json_dict["@type"] = json_dict.pop("_type")
        return json_dict

    def to_json(self, file=None, *args, **kwargs):
        """Serializes this entity to a JSON-formatted string or file."""
        if file is None:
            # return json string
            return json.dumps(self.to_dict(), *args, **kwargs)
        else:
            # write json file
            return json.dump(self.to_dict(), file, *args, **kwargs)

    @classmethod
    def from_json(cls, json_doc):
        """Parses an entity from a JSON-formatted string, dictionary, or file."""
        # convert json string or file to dictionary (if necessary)
        if isinstance(json_doc, str):
            json_doc = json.loads(json_doc)
        elif hasattr(json_doc, 'read'):
            json_doc = json.load(json_doc)
        # if pre-formatted, return directly
        if (json_doc is None or isinstance(json_doc, Number) or isinstance(json_doc, Entity)):
            return json_doc
        # if list, recursively parse each element and return mapped list
        if isinstance(json_doc, list):
            return map(lambda e: cls.from_json(e), json_doc)
        # otherwise use class method to initialize from normalized dictionary
        return cls.from_dict(json_doc)

    @staticmethod
    def from_dict(d):
        """Parses an entity from a normalized JSON dictionary."""
        return Entity(_id = d.get("@id", None))

    def __eq__(self, other):
        """Overrides the default check if this entity is equal to another by
        comparing unique identifiers if available.
        """
        # if specified, perform comparison on unique identifier
        if self._id is not None:
            return self._id == other._id
        # otherwise return the default comparison using `is` operator
        else:
            return self is other

    def __ne__(self, other):
        """Overrides the default check if this entity is not equal to another
        by comparing unique identifiers if available. (n.b. required for Python 2)
        """
        return not self.__eq__(other)

    def __hash__(self):
        """Overrides the default hash function by using the unique identifiers
        if available."""
        # if specified, return the hash of the unique identifier
        if self._id is not None:
            return hash(self._id)
        # otherwise return default hash from superclass
        return super(Entity, self).__hash__()

class EnumEntity(str, Enum):
    """Enumeration of recognized types."""

    @classmethod
    def get(cls, key):
        """Attempts to parse a type from a string, otherwise returns None."""
        if isinstance(key, cls):
            return key
        elif isinstance(key, list):
            return list(map(lambda e: cls.get(e), key))
        else:
            try: return cls(key.upper())
            except: return None

class Constants(object):
    """ Enumeration of various constants. Unless indicated otherwise, the constants 
        are in S.I. units. 
    """
    radiusOfEarthInKM = 6378.137 # Nominal equatorial radius of Earth
    speedOfLight = scipy.constants.speed_of_light
    GMe = 3.986004418e14*1e-9 # product of Gravitaional constant and Mass of Earth [km^3 s^-2]
    Boltzmann = scipy.constants.physical_constants['Boltzmann constant'][0]
    angularSpeedofEarthInRADpS = 7292115e-11 # WGS-84 nominal mean angular velocity of Earth
    Planck = scipy.constants.physical_constants['Planck constant'][0]
    SunBlackBodyTemperature = 6000.0 # Sun black body temperature in Kelvin
    SolarRadius = 6.95700e8 # Solar radius


class ManueverType(EnumEntity):
    """Enumeration of recognized manuvever types"""
    FIXED = "FIXED"
    CONE = "CONE",
    ROLLONLY = "ROLLONLY",
    YAW180 = "YAW180"
    YAW180ROLL = "YAW180ROLL" 
class OrientationConvention(EnumEntity):
    """Enumeration of recognized instrument orientation conventions."""
    XYZ = "XYZ"
    NADIR = "NADIR"
    SIDE_LOOK = "SIDE_LOOK"

class Orientation(Entity):
    """ Class to handle instrument orientation.
        Orientation is maintained internally via three rotation angles in the following order, xRotation -> yRotation -> zRotation, (intrinsic) with respect to "**Nadir-frame**" (A GMAT defined term). 

        At zero rotation, the instrument is aligned to the satellite body-frame which in turn is aligned to the Nadir-frame nominally. 
        Thus the instrument orientation with respect to the satellite body-frame is the same as it’s orientation with respect to the 
        Nadir-frame in the scenario that the satellite body-frame is aligned to the Nadir-frame. It is also assumed that the instrument
        imaging axis is along the instrument z-axis.

        **Nadir-frame:**
        * :math:`\\bf X_{nadir}` axis: :math:`-({\\bf Z_{nadir}} \\times {\\bf V})`, where :math:`\\bf V` is the Velocity vector of satellite in EarthFixed frame)
        * :math:`\\bf Y_{nadir}` axis: :math:`({\\bf Z_{nadir}} \\times {\\bf X_{nadir}})`
        * :math:`\\bf Z_{nadir}` axis: Aligned to Nadir vector (position vector of satellite in EarthFixed frame)

        It is also assumed that the instrument imaging axis is along the instrument z-axis.
     
        :ivar euler_angle1: (deg) Rotation about X-axis
        :vartype euler_angle1: float

        :ivar euler_angle2: (deg) Rotation about Y-axis
        :vartype euler_angle2: float

        :ivar euler_angle3: (deg) Rotation about Z-axis
        :vartype euler_angle3: float

        :ivar euler_seq1: Axis number to be rotated the first time.
        :vartype euler_angle1: int

        :ivar euler_seq2: Axis number to be rotated the second time.
        :vartype euler_angle2: int

        :ivar euler_seq3: Axis number to be rotated the third time.
        :vartype euler_angle3: int

        .. todo:: Add support to rotations specified by other Euler sequences (currently only 1,2,3 is supported.).
    
    """
    def __init__(self, euler_angle1, euler_angle2, euler_angle3, euler_seq1 = int(1), euler_seq2 = int(2), euler_seq3 = int(3), _id = None):
        self.euler_angle1 = float(euler_angle1)%360 
        self.euler_angle2 = float(euler_angle2)%360 
        self.euler_angle3 = float(euler_angle3)%360
        self.euler_seq1 = euler_seq1
        self.euler_seq2 = euler_seq2
        self.euler_seq3 = euler_seq3
        super(Orientation, self).__init__(_id, "Orientation")
    
    @classmethod
    def from_sideLookAngle(cls, side_look_angle_deg, _id = None):
        ''' Translate side-look-angle user input parameter to internal instrument orientation specs.
        
        :param side_look_angle_deg: Side look angle, also commonly called as nadir/ off-nadir angle in degrees. A positive angle corresponds to anti-clockwise rotation applied around the instrument y-axis.
        :paramtype side_look_angle_deg: float

        :param _id: Unique identifier
        :paramtype _id: str

        :return: Corresponding `Orientation` object
        :rtype: :class:`instrupy.util.Orientation`

        '''
        return Orientation(0.0, side_look_angle_deg, 0.0, 1,2,3,_id)
    
    @classmethod
    def from_XYZ_rotations(cls,  xRotation, yRotation, zRotation, _id = None):
        ''' Translate XYZ user input parameter to internal instrument orientation specs. 
            Orientation is specified via three rotation angles, xRotation -> yRotation -> zRotation..

        .. note:: This orientation specification convention input to this function is same as the GMAT-OC orientation convention specification.   
        
        :ivar x_rot_deg: (deg) Rotation about X-axis
        :vartype x_rot_deg: float

        :ivar y_rot_deg: (deg) Rotation about Y-axis
        :vartype y_rot_deg: float

        :ivar z_rot_deg: (deg) Rotation about Z-axis
        :vartype z_rot_deg: float`

        '''
        return Orientation(xRotation, yRotation, zRotation, 1,2,3,_id)       
        
    @staticmethod
    def from_dict(d):
        """Parses orientation specifications from a normalized JSON dictionary.
        
        :return: Parsed python object 
        :rtype: :class:`instrupy.util.Orientation`
        """
        _orien = OrientationConvention.get(d.get("convention", None))
        if(_orien == "XYZ"):
            return Orientation.from_XYZ_rotations(xRotation = d.get("xRotation", None), yRotation = d.get("yRotation", None), zRotation = d.get("zRotation", None), _id = d.get("@id", None))
        elif(_orien == "SIDE_LOOK"):
            return Orientation.from_sideLookAngle(d.get("sideLookAngle", None), d.get("@id", None))
        elif(_orien == "NADIR"):
            return Orientation.from_sideLookAngle(0, d.get("@id", None))
        else:
            raise Exception("Invalid or no Orientation convention specification")
            
class SensorGeometry(EnumEntity):
    """Enumeration of recognized instrument sensor geometries."""
    CONICAL = "CONICAL",
    RECTANGULAR = "RECTANGULAR",
    CUSTOM = "CUSTOM"

class FieldOfView(Entity):
        """ Class to handle instrument field-of-view, field-of-regard specifications.
        The Sensor FOV is maintained internally via vector of cone and clock angles. This is the same definition as that of :class:`CustomSensor` class definition in the GMAT-OC code.
       
        :ivar _geometry: Geometry of the sensor field-of-view. Accepted values are "CONICAL", "RECTANGULAR" or "CUSTOM".
        :vartype _geometry: str

        :ivar _coneAngleVec_deg: (deg) Array of cone angles measured from +Z sensor axis. If (:math:`xP`, :math:`yP`, :math:`zP`) is a unit vector describing a FOV point, then the 
                                 cone angle for the point is :math:`\\pi/2 - \\sin^{-1}zP`.
        :vartype _coneAngleVec_deg: list, float

        :ivar _clockAngleVec_deg: (deg) Array of clock angles (right ascensions) measured anti-clockwise from the + X-axis.  if (:math:`xP`, :math:`yP`, :math:`zP`) is a unit vector
                                  describing a FOV point, then the clock angle for the point is :math:`atan2(yP,xP)`.
        :vartype _clockAngleVec_deg: list, float

        :ivar _AT_fov_deg: Along track FOV in degrees (only if CONICAL or RECTANGULAR geometry)
        :vartype _AT_fov_deg: float

        :ivar _CT_fov_deg: Cross track FOV in degrees (only if CONICAL or RECTANGULAR geometry)
        :vartype _CT_fov_deg: float

        :ivar yaw180_flag: Flag applies in case of field-of-regard. If true, it signifies that the field-of-regard includes the field-of-view of payload rotated about the yaw axis by 180 deg. 
        :vartype yaw180_flag: bool

        .. note:: :code:`coneAnglesVec_deg[0]` ties to :code:`clockAnglesVec_deg[0]`, and so on. Except for the case of *CONICAL* FOV, in which we 
                  have just one :code:`clockAnglesVec_deg[0] = 1/2 full_cone_angle_deg` and no corresponding clock angle. 
        """
        def __init__(self, geometry = None, coneAnglesVec_deg = None, clockAnglesVec_deg = None, 
                     AT_fov_deg = None, CT_fov_deg = None, yaw180_flag = None, _id = None):                       

            if(coneAnglesVec_deg is not None):
                if(isinstance(coneAnglesVec_deg, list)):
                    self._coneAngleVec_deg = list(map(float, coneAnglesVec_deg))
                    self._coneAngleVec_deg = [x%360 for x in self._coneAngleVec_deg]
                else:
                    self._coneAngleVec_deg = [float(coneAnglesVec_deg)%360]
            else:
                self._coneAngleVec_deg = None
            
            if(clockAnglesVec_deg is not None):
                if(isinstance(clockAnglesVec_deg, list)):
                    self._clockAngleVec_deg = list(map(float, clockAnglesVec_deg))
                    self._clockAngleVec_deg = [x%360 for x in self._clockAngleVec_deg]
                else:
                    self._clockAngleVec_deg = [float(clockAnglesVec_deg)%360]
            else:
                self._clockAngleVec_deg = None

            self._geometry = geometry if geometry is not None else None
            self._AT_fov_deg = AT_fov_deg if AT_fov_deg is not None else None
            self._CT_fov_deg = CT_fov_deg if CT_fov_deg is not None else None
            self._yaw180_flag = bool(yaw180_flag) if yaw180_flag is not None else None

            super(FieldOfView, self).__init__(_id, "FieldOfView")

        @classmethod
        def from_customFOV(cls, coneAnglesVec_deg = None, clockAnglesVec_deg = None, yaw180_flag = False, _id = None):
            """  Return corresponding :class:`instrupy.util.FieldOfView` object.

                .. note:: The number of values in :code:`customConeAnglesVector` and :code:`customClockAnglesVector` should be the same (or) the number of 
                          values in :code:`customConeAnglesVector` should be just one (corresponding to CONICAL Sensor) and no values in 
                          :code:`customClockAnglesVector`.
            
            """
            if(coneAnglesVec_deg):
                if(not isinstance(coneAnglesVec_deg, list)):
                    coneAnglesVec_deg = [coneAnglesVec_deg]
            else:
                raise Exception("No cone angle vector specified!")
            
            if(clockAnglesVec_deg):
                if(not isinstance(clockAnglesVec_deg, list)):
                    clockAnglesVec_deg = [clockAnglesVec_deg]

            if(any(coneAnglesVec_deg) < 0 or any(coneAnglesVec_deg) > 90):
                raise Exception("CUSTOM cone angles specification must be in the range 0 deg to 90 deg.")

            if(len(coneAnglesVec_deg) == 1 and (clockAnglesVec_deg is not None)):
                raise Exception("With only one cone angle specified, there should be no clock angles specified.")
                
            if(not(len(coneAnglesVec_deg) == 1 and (clockAnglesVec_deg is None))):
                if(len(coneAnglesVec_deg) != len(clockAnglesVec_deg)):
                    raise Exception("With more than one cone angle specified, the length of cone angle vector should be the same as length of the clock angle vector.")
                
            return FieldOfView("CUSTOM", coneAnglesVec_deg, clockAnglesVec_deg, None, _id)

        @classmethod
        def from_conicalFOV(cls, full_cone_angle_deg = None, yaw180_flag = False, _id = None):
            ''' Convert user-given conical sensor specifications to cone, clock angles and return corresponding :class:`instrupy.util.FieldOfView` object.
            
            :param full_cone_angle_deg: Full conical angle of the Cone Sensor in degrees.
            :paramtype full_cone_angle_deg: float

            :param _id: Unique identifier
            :paramtype _id: str

            :return: Corresponding `FieldOfView` object
            :rtype: :class:`instrupy.util.FieldOfView`

            '''
            if full_cone_angle_deg is None:
                raise Exception("Please specify full-cone-angle of the CONICAL fov instrument.")

            if(full_cone_angle_deg < 0 or full_cone_angle_deg > 180):
                raise Exception("Specified full-cone angle of CONICAL sensor must be within the range 0 deg to 180 deg")

            return FieldOfView("CONICAL", 0.5*full_cone_angle_deg, None, full_cone_angle_deg, full_cone_angle_deg, yaw180_flag, _id)

        @classmethod
        def from_rectangularFOV(cls, along_track_fov_deg = None, cross_track_fov_deg = None, yaw180_flag = False, _id = None):
            ''' Convert the along-track (full) fov and cross-track (full) fov specs to clock, cone angles and return corresponding :class:`instrupy.util.FieldOfView` object.
                Along-track fov **must** be less than cross-track fov.

            :param along_track_fov_deg: Along track field-of-view (full) in degrees
            :paramtype along_track_fov_deg: float

            :param cross_track_fov_deg: Cross track field-of-view (full) in degrees
            :paramtype cross_track_fov_deg: float
            
            :param _id: Unique identifier
            :paramtype _id: str

            :return: Corresponding `FieldOfView` object
            :rtype: :class:`instrupy.util.FieldOfView`

            .. seealso:: The :class:`CustomSensor` class definition in GMAT-OC code. 

            .. warning:: The along-track FOV and cross-track FOV specs are assigned assuming the instrument is in nominal orientation, i.e. the instrument is
                         aligned to the satellite body frame, and further the satellite is aligned to nadir-frame.
                         If the instrument is rotated about the satellite body frame (by specifying the non-zero orientation angles in the instrument json specs file), the actual along-track
                         and cross-track fovs maybe different.                         

            '''
            if(along_track_fov_deg is None or cross_track_fov_deg is None):
                raise Exception("Please specify the along-track and cross-track fovs for the rectangular fov instrument.")
            
            if(along_track_fov_deg < 0 or along_track_fov_deg > 180 or cross_track_fov_deg < 0 or cross_track_fov_deg > 180):
                raise Exception("Specified along-track and cross-track fovs of RECTANGULAR sensor must be within the range 0 deg to 180 deg")       
            
            alongTrackFieldOfView = np.deg2rad(along_track_fov_deg)
            crossTrackFieldOfView = np.deg2rad(cross_track_fov_deg)

            cosCone = np.cos(alongTrackFieldOfView/2.0)*np.cos(crossTrackFieldOfView/2.0)
            cone = np.arccos(cosCone)

            sinClock =  np.sin(alongTrackFieldOfView/2.0) / np.sin(cone)

            clock = np.arcsin(sinClock)

            cone_deg = np.rad2deg(cone)
            clock_deg = np.rad2deg(clock)

            coneAngles_deg = [cone_deg, cone_deg, cone_deg, cone_deg]

            clockAngles_deg = [clock_deg, 180.0 - clock_deg, 180.0 + clock_deg, -clock_deg]

            return FieldOfView("RECTANGULAR", coneAngles_deg, clockAngles_deg,along_track_fov_deg, cross_track_fov_deg, yaw180_flag, _id)

        @staticmethod
        def from_dict(d):
            """Parses field-of-view specifications from a normalized JSON dictionary. If manueverability is included, the 
               parsed object corresponds to the field-of-regard. 
    
               :param d: Dictionary with the instrument field-of-view and maneverability (optional) specifications.
               :paramtype d: dict

               :return: Field-of-view/ Field-of-regard object 
               :rtype: :class:`instrupy.util.FieldOfView`

            """          
            # calulate the field-of-view
            _geo = SensorGeometry.get(d.get("sensorGeometry", None))
            if(_geo == "CONICAL"):
                fldofview = FieldOfView.from_conicalFOV(d.get("fullConeAngle", None), False, d.get("_id", None))
            elif(_geo == "RECTANGULAR"):
                fldofview = FieldOfView.from_rectangularFOV(d.get("alongTrackFieldOfView", None), d.get("crossTrackFieldOfView", None), False, d.get("_id", None))
            elif(_geo == "CUSTOM"):
                fldofview = FieldOfView.from_customFOV(d.get("customConeAnglesVector", None), d.get("customClockAnglesVector", None), False, d.get("_id", None))  
            else:
                raise Exception("Invalid Sensor FOV specified")

            if(d.get("maneuverability", None)):
                # Calculate the field-of-regard
                manuv = d.get("maneuverability")
                
                mv_type = ManueverType.get(manuv["@type"])
                if(mv_type == 'FIXED' or mv_type == 'YAW180'):
                    pass
                elif(mv_type == 'CONE'):
                    mv_cone = 0.5 * float(manuv["fullConeAngle"])
                elif(mv_type == 'ROLLONLY' or mv_type=='YAW180ROLL'):
                    mv_ct_range = float(manuv["rollMax"]) - float(manuv["rollMin"])
                else:
                    raise Exception('Invalid manuver type.')                                 

                if(mv_type == 'FIXED'):
                    fr_geom = fldofview._geometry
                    fr_at = fldofview._AT_fov_deg
                    fr_ct = fldofview._CT_fov_deg
                
                elif(mv_type == 'YAW180'):
                    fr_geom = fldofview._geometry
                    fr_at = fldofview._AT_fov_deg
                    fr_ct = fldofview._CT_fov_deg

                elif(mv_type == 'CONE'):
                    if(fldofview._geometry == 'CONICAL'):
                        fr_geom = 'CONICAL'
                        fr_at =2*(mv_cone + fldofview._coneAngleVec_deg[0])
                        fr_ct = fr_at

                    elif(fldofview._geometry == 'RECTANGULAR'):
                        fr_geom = 'CONICAL'
                        diag_half_angle = np.rad2deg(np.arccos(np.cos(np.deg2rad(0.5*fldofview._AT_fov_deg))*np.cos(np.deg2rad(0.5*fldofview._CT_fov_deg))))
                        fr_at = 2*(mv_cone +  diag_half_angle)
                        fr_ct = fr_at

                    else:
                        raise Exception('Invalid FOV geometry')    

                elif(mv_type == 'ROLLONLY'  or mv_type=='YAW180ROLL'):
                    if(fldofview._geometry == 'CONICAL'):
                        print("Approximating FOR as rectangular shape")
                        fr_geom = 'RECTANGULAR'
                        fr_at = 2*(fldofview._coneAngleVec_deg[0])
                        fr_ct = 2*(0.5*mv_ct_range + fldofview._coneAngleVec_deg[0])

                    elif(fldofview._geometry == 'RECTANGULAR'):
                        fr_geom = 'RECTANGULAR'
                        fr_at = fldofview._AT_fov_deg
                        fr_ct = mv_ct_range + fldofview._CT_fov_deg
                    else:
                        raise Exception('Invalid FOV geometry')

                if(mv_type=='YAW180ROLL' or mv_type=='YAW180'):
                    fr_yaw180_flag = True
                else:
                    fr_yaw180_flag = False

                if(fr_geom == 'CONICAL'):
                    fldofreg = FieldOfView.from_conicalFOV(fr_at, fr_yaw180_flag, d.get("_id", None))
                elif(fr_geom == 'RECTANGULAR'):
                    # Get the cone and clock angles from the rectangular FOV specifications.
                    fldofreg = FieldOfView.from_rectangularFOV(fr_at, fr_ct, fr_yaw180_flag, d.get("_id", None))
                
                return fldofreg
            else:
                # no manueverability specified, return just the calculated field-of-view
                return fldofview
        
        def get_cone_clock_fov_specs(self):
            """ Function to the get the cone and clock angle vectors from the resepective FieldOfView object.

                :return: Cone, Clock angles in degrees
                :rtype: list, float

            """
            return [self._coneAngleVec_deg, self._clockAngleVec_deg]

        def get_rectangular_fov_specs_from_custom_fov_specs(self):
            """ Function to get the rectangular fov specifications (along-track, cross-track fovs in degrees), from the sensor initialized
                clock, cone angles.           

                :return: Along-track and cross-track fovs in degrees
                :rtype: :class:`numpy.array`, float
      
                .. todo:: Make sure selected clock angle is from first quadrant. 
            """
            # Check if the instance does correspond to an rectangular fov.
            # Length of cone angle vector and clock angle vector must be 4.
            if(len(self._coneAngleVec_deg) != 4):
                raise Exception("This FieldOfView instance does not correspond to a rectangular fov.")

            # Check that all elements in the cone angle vector are the same value.
            if(len(set(self._coneAngleVec_deg))!= 1):
                raise Exception("This FieldOfView instance does not correspond to a rectangular fov.") 

            # The elements of the clock angle vector fold the following relationship: [theta, 180-theta, 180+theta, 360-theta]
            # Check for this relationship
            if(self._clockAngleVec_deg is not None):
                if(not math.isclose(self._clockAngleVec_deg[3],(360-self._clockAngleVec_deg[0])) or not math.isclose(self._clockAngleVec_deg[1], (180 - self._clockAngleVec_deg[0])) or not math.isclose(self._clockAngleVec_deg[2], (180 + self._clockAngleVec_deg[0]))):
                    raise Exception("This FieldOfView instance does not correspond to a rectangular fov.") 
            else:
                raise Exception("This FieldOfView instance does not correspond to a rectangular fov.") 
            
            theta = np.deg2rad(self._coneAngleVec_deg[0])
            omega = np.deg2rad(self._clockAngleVec_deg[0])

            alpha = np.arcsin(np.sin(theta)*np.sin(omega))
            beta = np.arccos(np.cos(theta)/np.cos(alpha))

            alTrFov_deg = 2*np.rad2deg(alpha)
            crTrFov_deg = 2*np.rad2deg(beta)

            return np.array([alTrFov_deg, crTrFov_deg]) 

        def get_ATCT_fov(self):
            """ Get the along-track and cross-track FOVs. Valid only for CONICAL and 
                RECTANGULAR FOV geometry.

                :return: Along-track and cross-track fovs in degrees
                :rtype: list, float
            """
            return [self._AT_fov_deg, self._CT_fov_deg]


                
class MathUtilityFunctions:
    """ Class aggregating various mathematical computation functions used in the InstruPy package. """

    @staticmethod
    def compute_satellite_footprint_speed(r,v):
        """ *Compute satellite footprint (**at Nadir** on ground-plane) linear speed*

            :param r: [distance] postion vector of satellite in ECI equatorial frame
            :paramtype r: list, float

            :param v: [distance/s] velocity vector of satellite in ECI equatorial frame
            :paramtype v: list, float

            :return: speed of satellite footprint in [m/s]
            :rtype: float

        """        
        # Calculate angular velocity of satellite
        r = np.array(r)
        v = np.array(v)

        omega = np.cross(r,v)/ (np.linalg.norm(r)**2) # angular velocity vector in radians per second
        
        # compensate for Earths rotation
        omega = omega - np.array([0,0,Constants.angularSpeedofEarthInRADpS])

        # Find linear speed of satellite footprint on ground [m/s].   
        return np.linalg.norm(omega)*Constants.radiusOfEarthInKM*1e3

    @staticmethod
    def latlonalt_To_Cartesian(lat_deg, lon_deg, alt_km):
        """ *LLA to ECEF vector considering Earth as sphere with a radius equal to the equatorial radius.*

            :param lat_deg: [deg] geocentric latitude 
            :paramtype lat_deg: float

            :param lon_deg: [deg] geocentric longitude
            :paramtype lon_deg: float

            :param alt_km: [km] Altitude
            :paramtype alt_km: float

            :return: [km] Position vector in ECEF
            :rtype: float, list
        
        """
        lat = np.deg2rad(lat_deg)
        if lon_deg<0:
            lon_deg = 360 + lon_deg
        lon = np.deg2rad(lon_deg)
        R = Constants.radiusOfEarthInKM + alt_km
        position_vector_km = np.array( [(np.cos(lon) * np.cos(lat)) * R,
                                (np.sin(lon) * np.cos(lat)) * R,
                                    np.sin(lat) * R] )     

        return position_vector_km


    @staticmethod
    def latlonaltGeodetic_To_Cartesian(lat_deg, lon_deg, alt_km):
        """ *LLA to ECEF vector considering Earth as WGS-84 ellipsoid.*

            :param lat_deg: [deg] WGS-84 geodetic latitude 
            :paramtype: float

            :param lon_deg: [deg] WGS-84 geodetic longitude
            :paramtype lon_deg: float

            :param alt_km: [km] Altitude
            :paramtype alt_km: float

            :return: [km] Position vector in WGS-84 ECEF
            :rtype: float, list
        
        """
        
        a = 6378137
        b = 6356752.3142
        f = (a - b) / a
        e_sq = f * (2-f)

        lamb = np.deg2rad(lat_deg)
        if lon_deg<0:
            lon_deg = 360 + lon_deg
        phi = np.deg2rad(lon_deg)
        h = alt_km*1e3

        s = math.sin(lamb)
        N = a / math.sqrt(1 - e_sq * s * s)

        sin_lambda = math.sin(lamb)
        cos_lambda = math.cos(lamb)
        sin_phi = math.sin(phi)
        cos_phi = math.cos(phi)

        x = (h + N) * cos_lambda * cos_phi
        y = (h + N) * cos_lambda * sin_phi
        z = (h + (1 - e_sq) * N) * sin_lambda

        position_vector_km = np.array( [x*1e-3, y*1e-3, z*1e-3] )     

        return position_vector_km

    
    @staticmethod
    def geo2eci(gcoord, JDtime):
        """ *Convert geographic spherical coordinates to Earth-centered inertial coords.*    
             
        :param gcoord: geographic coordinates of point [latitude [deg] ,longitude [deg], altitude [km]]. Geographic coordinates assume the Earth is a perfect sphere, with radius 
                     equal to its equatorial radius.
        :paramtype  gcoord: list, float

        :param JDtime: Julian Day time.
        :paramtype JDtime: float

        :return: A 3-element array of ECI [X,Y,Z] coordinates in kilometers. The TOD epoch is the supplied JDtime.                           
        :rtype: float

        .. seealso:: 
            * :mod:`JD2GMST`
            * `IDL Astronomy Users Library <https://idlastro.gsfc.nasa.gov/ftp/pro/astro/geo2eci.pro>`_
        
        EXAMPLES:

        .. code-block:: python

               ECIcoord = geo2eci([0,0,0], 2452343.38982663)
               print(ECIcoord)
              -3902.9606       5044.5548       0.0000000
        
        (The above is the ECI coordinates of the intersection of the equator and
        Greenwich's meridian on 2002/03/09 21:21:21.021)             
        
        """
        lat = np.deg2rad(gcoord[0])
        lon = np.deg2rad(gcoord[1])
        alt = gcoord[2]
               
        gst = MathUtilityFunctions.JD2GMST(JDtime)
        
        angle_sid=gst*2.0*np.pi/24.0 # sidereal angle

        theta = lon + angle_sid # azimuth
        r = ( alt + Constants.radiusOfEarthInKM ) * np.cos(lat)
        X = r*np.cos(theta)
        Y = r*np.sin(theta)
        Z = ( alt + Constants.radiusOfEarthInKM )*np.sin(lat)
                
        return [X,Y,Z]
        
   
    @staticmethod
    def normalize(v):
        """ Normalize a input vector.

            :param v: Input vector
            :paramtype v: list, float

            :return: Normalized vector
            :rtype: :obj:`np.array`, float
        
        """
        v = np.array(v)
        norm = np.linalg.norm(v)
        if(norm == 0):
            raise Exception("Encountered division by zero in vector normalization function.")
        return v / norm

            
        
    
    @staticmethod
    def angle_between_vectors(vector1, vector2):
        """ Find angle between two input vectors in radians. Use the dot-product relationship to obtain the angle.
        
            :param vector1: Input vector 1
            :paramtype vector1: list, float

            :param vector2: Input vector 2
            :paramtype vector2: list, float

            :return: [rad] Angle between the vectors, calculated using dot-product relationship.
            :rtype: float

        """
        unit_vec1 = MathUtilityFunctions.normalize(vector1)
        unit_vec2 = MathUtilityFunctions.normalize(vector2)
        return np.arccos(np.dot(unit_vec1, unit_vec2))


    @staticmethod
    def JD2GMST(JD):
        """ Convert Julian Day to Greenwich Mean Sidereal Time.
            Reference `USNO NAVY <https://aa.usno.navy.mil/faq/docs/GAST.php>`_

            :param JD: Julian Date UT1
            :paramtype JD: float

            :return: [hrs] Greenwich Mean Sidereal Time at the corresponding JD
            :rtype: float

        """
        _x = round(JD) + 0.5
        if(_x > JD):
            JD0 = round(JD) - 0.5
        else:
            JD0 = round(JD) + 0.5
      
        D = JD - 2451545.0
        D0 = JD0 - 2451545.0
        H = (JD - JD0)*24.0
        T = D / 36525.0

        # Greenwich mean sidereal time in hours
        GMST = 6.697374558 + 0.06570982441908*D0 + 1.00273790935*H + 0.000026*T**2
        
        GMST = GMST % 24

        return GMST


    @staticmethod
    def find_closest_value_in_array(array, value):
        """ Find the value in an array closest to the supplied scalar value.

            :param array: Array under consideration
            :paramtype array: list, float

            :param value: Value under consideration
            :paramtype value: float

            :return: Value in array corresponding to one which is closest to supplied value, Index of the value in array
            :rtype: list, float

        """
        array = np.asarray(array)
        idx = (np.abs(array - value)).argmin()
        return [array[idx], idx]
    
    @staticmethod
    def SunVector_ECIeq(Time_JDUT1):
        """Find Sun Vector in Earth Centered Inertial (ECI) equatorial frame.
           Algorithm in David A.Vallado, Fundamental of Astrodynamics and Applications, 4th ed, Page 280.
           Also see accompanying software scripts with the book.

           :param Time_JDUT1: Time in Julian Day UT1
           :paramtype Time_JDUT1: float

           :return: Sun-vector in ECI equatorial frame (km)
           :rtype: list, float
           
        """
        T_UT1 = (Time_JDUT1 - 2451545.0) / 36525

        T_TDB = T_UT1

        lambda_M_deg = (280.460 + 36000.77*T_UT1)
        lambda_M_deg = lambda_M_deg % 360
        
        M = np.deg2rad((357.5277233 + 35999.05034*T_TDB))%(2*np.pi)

        if M <0:
            M = 2*np.pi + M

        lambda_ecliptic_deg = (lambda_M_deg + 1.914666471* np.sin(M) + 0.019994643*np.sin(2*M))
        lambda_ecliptic = np.deg2rad(lambda_ecliptic_deg % 360)

        eps = np.deg2rad(23.439291 - 0.0130042*T_TDB)

        
        # find the unit Sun vector in ECI frame
        unitSv_ECI = [np.cos(lambda_ecliptic), np.cos(eps)*np.sin(lambda_ecliptic), np.sin(eps)*np.sin(lambda_ecliptic)]
        
        
        # magnitude of distance to the Sun from Earth center (note: not from the spacecraft)
        r_Sun_km = (1.000140612 - 0.016708617 * np.cos(M) - 0.000139589 * np.cos(2*M)) * 149597870.700

        
        # complete vector of Sun from Earth Center (ECI equatorial frame)
        Sv_ECI_km = [x*r_Sun_km for x in unitSv_ECI]
        
        return Sv_ECI_km

    @staticmethod
    def checkLOSavailability(object1_pos, object2_pos, obstacle_radius):
        """ Determine if line-of-sight exists
             Algorithm from Page 198 Fundamental of Astrodynamics and Applications, David A.Vallado is used. The first algorithm
             described is used.

             :param object1_pos: Object 1 position vector
             :paramtype object1_pos: float

             :param object2_pos: Object2 position vector
             :paramtype object2_pos: float

             :param obstacle_radius: Radius of spherical obstacle
             :paramtype obstacle_radius: float

             :return: T/F flag indicating availability of line of sight from object1 to object2.
             :rtype: bool

             .. note: The frame of reference for describing the object positions must be centered at spherical obstacle.
        
        """        
        obj1_unitVec = MathUtilityFunctions.normalize(object1_pos)
        obj2_unitVec = MathUtilityFunctions.normalize(object2_pos)  

        # This condition tends to give a numerical error, so solve for it independently.
        eps = 1e-9
        x = np.dot(obj1_unitVec, obj2_unitVec)
        
        if((x > -1-eps) and (x < -1+eps)):
            return False
        else:
            if(x>1): 
                x=1
            theta  = np.arccos(x)
                
        obj1_r = np.linalg.norm(object1_pos)
        if(obj1_r - obstacle_radius > 1e-5):
            theta1 = np.arccos(obstacle_radius/obj1_r)
        elif(abs(obj1_r - obstacle_radius) < 1e-5):
            theta1 =  0.0
        else:
            return False # object1 is inside the obstacle

        obj2_r = np.linalg.norm(object2_pos)
        if(obj2_r - obstacle_radius > 1e-5):
            theta2 = np.arccos(obstacle_radius/obj2_r)
        elif(abs(obj2_r - obstacle_radius) < 1e-5):
            theta2 =  0.0
        else:
            return False # object2 is inside the obstacle
                
        if (theta1 + theta2 < theta):
            return False
        else:
            return True     

    @staticmethod
    def compute_sun_zenith(time_JDUT1, pos_km):
        """ Compute the Sun zenith angle at the given time and location. 

            :param time_JDUT1: Time in Julian Day UT1
            :paramtype time_JDUT1: float

            :param pos_km: Position vector [km] of point in ECI frame at which the Sun elevation angle is to be calculated.  
            :paramtype pos_km: list, float

            :returns: Sun zenith angle in [radians], Distance from Sun in [km]
            :rtype: list, float
            
        """
        # calculate Sun-vector in ECI frame
        sunVector_km = MathUtilityFunctions.SunVector_ECIeq(time_JDUT1)
        # verify point-of-interest is below in night region
        if(MathUtilityFunctions.checkLOSavailability(pos_km, sunVector_km, Constants.radiusOfEarthInKM) is False):
            return [None, None]

        # calculate sun to target vector
        Sun2Tar_pos_km = np.array(pos_km) - np.array(sunVector_km)
        # calculate solar incidence angle
        solar_zenith_angle_rad = MathUtilityFunctions.angle_between_vectors(pos_km, -1*np.array(Sun2Tar_pos_km))

        solar_distance = np.linalg.norm(Sun2Tar_pos_km)

        return [solar_zenith_angle_rad, solar_distance]

    @staticmethod
    def calculate_derived_satellite_coords(tObs_JDUT1, obs_position_km, obs_vel_vec_kmps, target_position_km):
        """ The satellite state supplied to the function may not be exactly at the middle of the access interval, and hence the state supplied would
            not correspond to time at which the satellite to point-of-interest line is exactly orthogonal to the ground-track. 

            :param tObs_JDUT1: Observation time in Julian Day UT1
            :paramtype tObs_JDUT1: float

            :param obs_position_km: Spacecraft position vector
            :paramtype obs_position_km: list, float

            :param obs_vel_vec_kmps: Observer velocity vector
            :paramtype obs_vel_vec_kmps: list, float

            :param target_position_km: Target (object being observed) position vector
            :paramtype target_position_km: list, float

            :returns: Derived observation time, position, range, altitude and incidence angle as a dictionary.
            :rtype: dict, float
            
        """
        obs_position_km = np.array(obs_position_km)
        obs_vel_vec_kmps = np.array(obs_vel_vec_kmps)
        target_position_km = np.array(target_position_km)

        #  Calculate range vector between spacecraft and POI (Target)
        range_vector_km = target_position_km - obs_position_km

        alt_km = np.linalg.norm(obs_position_km) - Constants.radiusOfEarthInKM
        look_angle = np.arccos(np.dot(MathUtilityFunctions.normalize(range_vector_km), -1*MathUtilityFunctions.normalize(obs_position_km)))
        incidence_angle_rad = np.arcsin(np.sin(look_angle)*(Constants.radiusOfEarthInKM + alt_km)/Constants.radiusOfEarthInKM)

        # 1. Get vector perpendicular to the cross-orbital plane, using the satellite velocity vector.
        crossOrbPlaneNorVec = np.array(MathUtilityFunctions.normalize(obs_vel_vec_kmps))
        # 2. Calculate projection of range-vector onto the cross-orbital-plane        
        r_vec_projCrossOrbPlane_km = range_vector_km - (np.dot(range_vector_km, crossOrbPlaneNorVec)*crossOrbPlaneNorVec)      
        # 3. Calculate the range-vector and viewing geometry to the "derived" observer position
        derived_range_vec_km = r_vec_projCrossOrbPlane_km
        derived_obs_pos_km = target_position_km - derived_range_vec_km
        derived_look_angle = np.arccos(np.dot(MathUtilityFunctions.normalize(derived_range_vec_km), -1*MathUtilityFunctions.normalize(derived_obs_pos_km)))

        derived_alt_km = np.linalg.norm(derived_obs_pos_km) - Constants.radiusOfEarthInKM
        derived_incidence_angle_rad = np.arcsin(np.sin(derived_look_angle)*(Constants.radiusOfEarthInKM + derived_alt_km)/Constants.radiusOfEarthInKM)

        travel_dis_km = np.linalg.norm(derived_obs_pos_km - obs_position_km) 
        travel_time_s = travel_dis_km/ np.linalg.norm(obs_vel_vec_kmps)
        derived_obsTime_JDUT1 = tObs_JDUT1 + travel_time_s
        
        #return {"derived_obsTime_JDUT1": tObs_JDUT1, "derived_obs_pos_km": obs_position_km, "derived_range_vec_km": range_vector_km, "derived_alt_km": alt_km, "derived_incidence_angle_rad": incidence_angle_rad}
        return {"derived_obsTime_JDUT1": derived_obsTime_JDUT1, "derived_obs_pos_km": derived_obs_pos_km.tolist(), "derived_range_vec_km": derived_range_vec_km.tolist(), "derived_alt_km": derived_alt_km, "derived_incidence_angle_rad": derived_incidence_angle_rad}
      
    @staticmethod
    def get_transmission_Obs2Space(wav_low_m, wav_high_m, obs_zenith_angle_rad):
        ''' Observer to space transmission loss

            :returns: transmissitivity in steps of 20cm-1
            :rtype: array, float
        '''        
        obs_alt_km = 0
        wav_low_nm = wav_low_m*1e9
        wav_high_nm = wav_high_m*1e9
        obs_zenith_angle_deg = np.rad2deg(obs_zenith_angle_rad)
        wav_step_percm = 40
        
        c1 = {'model': 6, # US standard
            'h1': obs_alt_km,
            'angle': obs_zenith_angle_deg, 
            'wlshort': wav_low_nm,
            'wllong': wav_high_nm,
            'wlstep': wav_step_percm,
            }
        
        TR = lowtran.transmittance(c1)
        TR = TR.where(TR['wavelength_nm']!=0, drop=True) # LowTran sometimes returns a entry with '0' wavelength (when the bandwidth is not "compatible" with the step-size)
        return TR 
    
    @staticmethod
    def get_eca(fov_deg, alt_km):

        RE = Constants.radiusOfEarthInKM 
        sinRho = RE/(RE + alt_km)
        hfov_deg = 0.5*fov_deg
        elev_deg = np.rad2deg(np.arccos(np.sin(np.deg2rad(hfov_deg))/sinRho))
        lambda_deg = 90 - hfov_deg - elev_deg # half-earth centric angle 
        eca_deg = lambda_deg*2 # total earth centric angle

        return eca_deg

class FileUtilityFunctions:

    @staticmethod
    def from_json(json_doc):
        """Parses a dictionary from a JSON-formatted string, dictionary, or file."""
        # convert json string or file to dictionary (if necessary)
        if isinstance(json_doc, str):
            json_doc = json.loads(json_doc)
        elif hasattr(json_doc, 'read'):
            json_doc = json.load(json_doc)
        # if pre-formatted, return directly
        if (json_doc is None or isinstance(json_doc, Number) or isinstance(json_doc, Entity)):
            return json_doc
        # if list, recursively parse each element and return mapped list
        if isinstance(json_doc, list):
            return map(lambda e: FileUtilityFunctions.from_json(e), json_doc)
        # otherwise use class method to initialize from normalized dictionary
        return json_doc 