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
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import shapely
import cartopy.geodesic
import math
import json

from math import radians, cos, sin, asin, sqrt
#import lowtran #TEMPORARY: PLEASE REMOVE if commented

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

class SyntheticDataConfiguration(Entity):
    def __init__(self, sourceFilePaths=None, environVar=None, interplMethod=None):
        self.sourceFilePaths = sourceFilePaths if sourceFilePaths is not None else list()
        self.environVar = environVar if environVar is not None else None 
        self.interplMethod = interplMethod if interplMethod is not None else None
    
    @staticmethod
    def from_dict(d):
        return SyntheticDataConfiguration(sourceFilePaths = d.get("sourceFilePaths", None), 
                                          environVar      = d.get("environVar", None),
                                          interplMethod    = d.get("interplMethod", None))

    def to_dict(self):
        """ Return data members of the object as python dictionary. 
        
            :return: SyntheticDataConfiguration object as python dictionary
            :rtype: dict
        """
        syndataconf_dict = {"sourceFilePaths": self.sourceFilePaths, "environVar": self.environVar, "interplMethod": self.interplMethod}
        return syndataconf_dict

class ReferenceFrame(EnumEntity):
    """ Enumeration of recognized reference frames.
    
        :cvar EARTH_CENTERED_INERTIAL: Earth Centered Inertial reference frame.
        
                                        This is an Earth equator inertial reference frame identical to EarthMJ2000Eq coordinate system used in GMAT.
                                        The nominal x-axis points along the line formed by the intersection of the Earth’s 
                                        mean equatorial plane and the mean ecliptic plane (at the J2000 epoch), in the direction
                                        of Aries. The z-axis is normal to the Earth’s mean equator at the J2000 epoch and the 
                                        y-axis completes the right-handed system. The mean planes of the ecliptic and equator, 
                                        at the J2000 epoch, are computed using IAU-1976/FK5 theory with 1980 update for nutation.

        :vartype EARTH_CENTERED_INERTIAL: str

        :cvar EARTH_FIXED: Earth Fixed reference frame.

                            The Earth Fixed reference frame is referenced to the Earth's equator and the prime meridian 
                            and is computed using IAU-1976/FK5 theory. This system is identical to the EarthFixed coordinate
                            system used in GMAT.

        :vartype EARTH_FIXED: str

        :cvar NADIR_POINTING: Nadir-pointing reference frame.

                            The axis of the Nadir-pointing reference frame are defined as follows:

                            * :math:`\\bf X_{np}` axis: :math:`-({\\bf Z_{np}} \\times {\\bf V})`, where :math:`\\bf V` is the Velocity vector of satellite in EARTH_FIXED frame)
                                    
                            * :math:`\\bf Y_{np}` axis: :math:`({\\bf Z_{np}} \\times {\\bf X_{np}})`
                                    
                            * :math:`\\bf Z_{np}` axis: Aligned to Nadir vector (i.e. the negative of the position vector of satellite in EARTH_FIXED frame)

                            .. figure:: nadirframe.png
                                :scale: 100 %
                                :align: center

                            .. todo:: Verify the claim about position vector and velocity vector in EARTH_FIXED frame.

        :vartype NADIR_POINTING: str

        :cvar SC_BODY_FIXED: Spacecraft Body Fixed reference frame. The axis of this coordinate system are aligned with the axis of the Spacecraft Bus. 
        :vartype SC_BODY_FIXED: str

    """
    EARTH_CENTERED_INERTIAL = "EARTH_CENTERED_INERTIAL"
    EARTH_FIXED = "EARTH_FIXED"
    NADIR_POINTING = "NADIR_POINTING"
    SC_BODY_FIXED = "SC_BODY_FIXED"

class Orientation(Entity):
    """ Class to store and handle orientation. Orientation is parameterized as intrinsic rotations specified by Euler angles and sequence 
        with respect to the user-specified reference frame. The definition of the Euler angle rotation is identical to the 
        one used in the orbitpy->propcov->extern->gmatutil->util->AttitudeUtil, AttitudeConversionUtility C++ classes.

        A Euler sequence = 123 implies the following rotation: R = R3.R2.R1, where Ri is the rotation matrix about the ith axis. 
        A positive angle corresponds to an anti-clockwise rotation about the respective axis.
        Each rotation matrix rotates the coordinate system (not the vector). 
        See:
 
         * https://mathworld.wolfram.com/RotationMatrix.html
     
        :ivar ref_frame: Reference frame. Default in "NADIR_POINTING".
        :vartype ref_frame: str

        :ivar euler_angle1: (deg) Rotation angle corresponding to the first rotation. Default is 0.
        :vartype euler_angle1: float

        :ivar euler_angle2: (deg) Rotation angle corresponding to the second rotation. Default is 0.
        :vartype euler_angle2: float

        :ivar euler_angle3: (deg) Rotation angle corresponding to the third rotation. Default is 0.
        :vartype euler_angle3: float

        :ivar euler_seq1: Axis-number corresponding to the first rotation. Default is 1.
        :vartype euler_angle1: int

        :ivar euler_seq2: Axis-number corresponding to the second rotation. Default is 2.
        :vartype euler_angle2: int

        :ivar euler_seq3: Axis-number corresponding to the third rotation. Default is 3.
        :vartype euler_angle3: int

        :ivar _id: Unique identifier.
        :vartype _id: str
    
    """    
    def __init__(self, ref_frame="NADIR_POINTING", euler_angle1=0, euler_angle2=0, euler_angle3=0, euler_seq1=int(1), euler_seq2=int(2), euler_seq3=int(3), _id=None):
        
        self.ref_frame = ReferenceFrame.get(ref_frame)
        self.euler_angle1 = float(euler_angle1)%360 
        self.euler_angle2 = float(euler_angle2)%360 
        self.euler_angle3 = float(euler_angle3)%360
        self.euler_seq1 = euler_seq1
        self.euler_seq2 = euler_seq2
        self.euler_seq3 = euler_seq3
        super(Orientation, self).__init__(_id, "Orientation")
    
    class Convention(EnumEntity):
        """ Enumeration of recognized orientation conventions. The rotations below can be specified with respect to 
            any of the reference frames given in :class:`instrupy.util.ReferenceFrame`.

        :cvar XYZ: Rotations about the X, Y and Z axis in the order 123.
        :vartype XYZ: str

        :cvar REF_FRAME_ALIGNED: Aligned with respective to the underlying reference frame. Identity rotation matrix. 
        :vartype REF_FRAME_ALIGNED: str

        :cvar SIDE_LOOK: Rotation about the Y axis only. 
        :vartype SIDE_LOOK: str

        :cvar EULER: Rotation according to the specified Euler angles and sequence.
        :vartype EULER: str

        """
        XYZ = "XYZ"
        REF_FRAME_ALIGNED = "REF_FRAME_ALIGNED"
        SIDE_LOOK = "SIDE_LOOK"
        EULER = "EULER"

    @classmethod
    def from_sideLookAngle(cls, ref_frame="NADIR_POINTING", side_look_angle=0, _id=None):
        """ Return :class:`Orientation` object constructed from the side-look angle. 
        
        :param ref_frame: Reference frame. Default in "NADIR_POINTING".
        :paramtype ref_frame: str

        :param side_look_angle: (deg) Side look angle. A positive angle corresponds to anti-clockwise rotation applied around the y-axis. Default is 0.
        :paramtype side_look_angle: float

        :param _id: Unique identifier.
        :paramtype _id: str

        :return: Corresponding `Orientation` object.
        :rtype: :class:`instrupy.util.Orientation`

        """
        return Orientation(ref_frame, 0.0, side_look_angle, 0.0, 1,2,3,_id)
    
    @classmethod
    def from_XYZ_rotations(cls, ref_frame="NADIR_POINTING", x_rot=0, y_rot=0, z_rot=0, _id = None):
        """ Return :class:`Orientation` object by the user-specified XYZ rotation angles with  
            the sequence=123.

        :param ref_frame: Reference frame. Default in "NADIR_POINTING".
        :paramtype ref_frame: str

        :ivar x_rot: (deg) Rotation about X-axis. Default is 0.
        :vartype x_rot: float

        :ivar y_rot: (deg) Rotation about Y-axis. Default is 0.
        :vartype y_rot: float

        :ivar z_rot: (deg) Rotation about Z-axis. Default is 0.
        :vartype z_rot: float`

        :param _id: Unique identifier.
        :paramtype _id: str

        """
        return Orientation(ref_frame, x_rot, y_rot, z_rot, 1,2,3,_id)       
        
    @staticmethod
    def from_dict(d):
        """Parses orientation specifications from a dictionary.
        
            :return: Parsed python object. 
            :rtype: :class:`instrupy.util.Orientation`
        """
        orien_conv = Orientation.Convention.get(d.get("convention", None))
        ref_frame = ReferenceFrame.get(d.get("referenceFrame", "NADIR_POINTING")).value # default reference frame is NADIR_POINTING
        if(orien_conv == "XYZ"):
            return Orientation.from_XYZ_rotations(ref_frame=ref_frame, x_rot=d.get("xRotation", None), y_rot=d.get("yRotation", None), z_rot=d.get("zRotation", None), _id = d.get("@id", None))
        elif(orien_conv == "SIDE_LOOK"):
            return Orientation.from_sideLookAngle(ref_frame=ref_frame, side_look_angle=d.get("sideLookAngle", None), _id=d.get("@id", None))
        elif(orien_conv == "REF_FRAME_ALIGNED"):
            return Orientation.from_sideLookAngle(ref_frame=ref_frame, side_look_angle=0, _id=d.get("@id", None))
        elif(orien_conv == "EULER"):
            return Orientation(ref_frame=ref_frame, euler_angle1=d.get("eulerAngle1", None), euler_angle2=d.get("eulerAngle2", None), 
                               euler_angle3=d.get("eulerAngle3", None), euler_seq1=d.get("eulerSeq1", None), euler_seq2=d.get("eulerSeq2", None),
                               euler_seq3=d.get("eulerSeq3", None), _id=d.get("@id", None))
        else:
            raise Exception("Invalid or no Orientation convention specification")
    
    def to_list(self):
        """ Return data members of the instance as list,
        
            :return: Orientation object data attributes as list.
            :rtype: list

            .. todo:: Remove this function
        """
        return [self.ref_frame, self.euler_seq1, self.euler_seq2, self.euler_seq3, self.euler_angle1, self.euler_angle2, self.euler_angle3]

    def to_dict(self):
        """ Return data members of the instance as python dictionary. 
        
            :return: Orientation object data attributes as python dictionary.
            :rtype: dict
        """
        orien_dict = {
                      "referenceFrame": self.ref_frame.value, 
                      "convention": "EULER", 
                      "eulerAngle1": self.euler_angle1,  
                      "eulerAngle2": self.euler_angle2, 
                      "eulerAngle3": self.euler_angle3,
                      "eulerSeq1": self.euler_seq1, 
                      "eulerSeq2": self.euler_seq2, 
                      "eulerSeq3": self.euler_seq3,
                      "@id": self._id
                     }
        return orien_dict
    
    def __repr__(self):
        if isinstance(self._id, str):
            return "Orientation(ref_frame='{}',euler_angle1={},euler_angle2={},euler_angle3={},euler_seq1={},euler_seq2={},euler_seq3={},_id='{}')".format(self.ref_frame, self.euler_angle1, self.euler_angle2, self.euler_angle3,
                                                                                            self.euler_seq1, self.euler_seq2, self.euler_seq3, self._id)
        else:
            return "Orientation(ref_frame='{}',euler_angle1={},euler_angle2={},euler_angle3={},euler_seq1={},euler_seq2={},euler_seq3={},_id={})".format(self.ref_frame, self.euler_angle1, self.euler_angle2, self.euler_angle3,
                                                                                          self.euler_seq1, self.euler_seq2, self.euler_seq3, self._id)
            
class SphericalGeometry(Entity):
        """ Class to handle spherical geometries (spherical polygons and circles) which are used to characterize the sensor 
            field-of-view (FOV) / scene FOV.

            The spherical polygon/circle is maintained internally via vector of cone and clock angles defined in the sensor frame. This is the same 
            definition as that in the orbitpy->propcov->lib->propcov-cpp CustomSensor C++ class. 

            .. figure:: cone_clock_angle.png
                :scale: 100 %
                :align: center
       
        :ivar shape: Shape of the spherical geometry. Accepted values are "CIRCULAR", "RECTANGULAR" or "CUSTOM".
        :vartype shape: str

        :ivar cone_angle_vec: (deg) Array of cone angles measured from +Z sensor axis (pointing axis). If (:math:`xP`, :math:`yP`, :math:`zP`) is a unit vector describing a point on unit sphere, then the 
                                 cone angle for the point is :math:`\\pi/2 - \\sin^{-1}zP`.
        :vartype cone_angle_vec: list, float

        :ivar clock_angle_vec: (deg) Array of clock angles (right ascensions) measured anti-clockwise from the +X-axis. If (:math:`xP`, :math:`yP`, :math:`zP`) is a unit vector
                                  describing a point on unit sphere, then the clock angle for the point is :math:`atan2(yP,xP)`.
        :vartype clock_angle_vec: list, float

        :ivar angle_height: (deg) Spherical geometry angular width (about sensor X axis) (only for CIRCULAR or RECTANGULAR shapes). Corresponds to along-track angular width if sensor frame is aligned to NADIR_POINTING frame.
        :vartype angle_height: float

        :ivar angle_width: (deg) Spherical geometry angular height (about sensor Y axis)  (only for CIRCULAR or RECTANGULAR shapes). Corresponds to cross-track angular width if sensor frame is aligned to NADIR_POINTING frame.
        :vartype angle_width: float

        :param _id: Unique identifier.
        :paramtype _id: str

        .. note:: :code:`cone_angle_vec[0]` ties to :code:`clock_angle_vec[0]`, and so on. Except for the case of *CIRCULAR* shape, in which we 
                  have only one cone angle (:code:`cone_angle_vec[0] = 1/2 diameter`) and no corresponding clock angle. 

        """
        class Shape(EnumEntity):
            """Enumeration of recognized SphericalGeometry shapes.
            
            :cvar CIRCULAR: Circular shape definition, characterized by the radius of the circle around the Z-axis.
            :vartype CIRCULAR: str

            :cvar RECTANGULAR: Rectangular polygon definition, characterized by angular width (about Y-axis) and angular height (about X-axis). 
            :vartype RECTANGULAR: str

            :cvar CUSTOM: Custom polygon definition, where an arbitrary number of cone, clock angles
                          denoting the vertices of the polygon can be specified. 
            :vartype CUSTOM: str
            
            """
            CIRCULAR = "CIRCULAR"
            RECTANGULAR = "RECTANGULAR"
            CUSTOM = "CUSTOM"
            
        def __init__(self, shape=None, cone_angle_vec=None, clock_angle_vec=None, _id=None):   
            if(cone_angle_vec is not None):
                if(isinstance(cone_angle_vec, list)):
                    self.cone_angle_vec = list(map(float, cone_angle_vec))
                    self.cone_angle_vec = [x%360 for x in self.cone_angle_vec]
                else:
                    self.cone_angle_vec = [float(cone_angle_vec)%360]
            else:
                self.cone_angle_vec = None
            
            if(clock_angle_vec is not None):
                if(isinstance(clock_angle_vec, list)):
                    self.clock_angle_vec = list(map(float, clock_angle_vec))
                    self.clock_angle_vec = [x%360 for x in self.clock_angle_vec]
                else:
                    self.clock_angle_vec = [float(clock_angle_vec)%360]
            else:
                self.clock_angle_vec = None

            self.shape = SphericalGeometry.Shape.get(shape) if shape is not None else None
            if(self.shape is SphericalGeometry.Shape.CIRCULAR):
                self.angle_height = 2 * self.cone_angle_vec[0]
                self.angle_width = self.angle_height
            elif(self.shape is SphericalGeometry.Shape.RECTANGULAR):
                [self.angle_height, self.angle_width] = SphericalGeometry.get_rect_poly_specs_from_cone_clock_angles(self.cone_angle_vec, self.clock_angle_vec)
            else:
                self.angle_height = None
                self.angle_width = None

            super(SphericalGeometry, self).__init__(_id, "SphericalGeometry")

        def to_dict(self):
            """ Return data members of the object as python dictionary. 

                :return: SphericalGeometry object as python dictionary
                :rtype: dict 
            """
            if self.shape==SphericalGeometry.Shape.CIRCULAR:
                sph_geom_dict = {"shape": "CIRCULAR", "diameter": self.angle_height, "@id": self._id}
            elif self.shape==SphericalGeometry.Shape.RECTANGULAR:
                sph_geom_dict = {"shape": "RECTANGULAR", "angleHeight": self.angle_height, "angleWidth": self.angle_width, "@id": self._id}
            elif self.shape==SphericalGeometry.Shape.CUSTOM:
                sph_geom_dict = {"shape": "CUSTOM", 
                            "customConeAnglesVector": "[" + ','.join(map(str, self.cone_angle_vec))  + "]", 
                            "customClockAnglesVector": "[" + ','.join(map(str, self.clock_angle_vec))  + "]",
                            "@id": self._id
                           }
            else:
                sph_geom_dict = None
            return sph_geom_dict

        @classmethod
        def from_custom_specs(cls, cone_angle_vec=None, clock_angle_vec=None, _id=None):
            """  Return corresponding :class:`instrupy.util.SphericalGeometry` object from input cone and clock angles.

                :param cone_angle_vec: (deg) Array of cone angles measured from +Z sensor axis. If (:math:`xP`, :math:`yP`, :math:`zP`) is a unit vector describing a point on unit sphere, then the 
                                 cone angle for the point is :math:`\\pi/2 - \\sin^{-1}zP`.
                :paramtype cone_angle_vec: list, float

                :param clock_angle_vec: (deg) Array of clock angles (right ascensions) measured anti-clockwise from the + X-axis. If (:math:`xP`, :math:`yP`, :math:`zP`) is a unit vector
                                        describing a point on the unit sphere, then the clock angle for the point is :math:`atan2(yP,xP)`.
                :paramtype clock_angle_vec: list, float

                :param _id: Unique identifier.
                :paramtype _id: str

                .. note:: :code:`cone_angle_vec[0]` ties to :code:`clock_angle_vec[0]`, and so on. Except for the case of *CIRCULAR* shaped FOV, in which we 
                    have only one cone angle (:code:`cone_angle_vec[0] = 1/2 diameter`) and no corresponding clock angle. 
            
            """
            if(cone_angle_vec):
                if(not isinstance(cone_angle_vec, list)):
                    cone_angle_vec = [cone_angle_vec]
            else:
                raise Exception("No cone angle vector specified!")
            
            if(clock_angle_vec):
                if(not isinstance(clock_angle_vec, list)):
                    clock_angle_vec = [clock_angle_vec]

            if(any(cone_angle_vec) < 0 or any(cone_angle_vec) > 90):
                raise Exception("CUSTOM cone angles specification must be in the range 0 deg to 90 deg.")

            if(len(cone_angle_vec) == 1 and (clock_angle_vec is not None)):
                raise Exception("With only one cone angle specified, there should be no clock angles specified.")
                
            if(not(len(cone_angle_vec) == 1 and (clock_angle_vec is None))):
                if(len(cone_angle_vec) != len(clock_angle_vec)):
                    raise Exception("With more than one cone angle specified, the length of cone angle vector should be the same as length of the clock angle vector.")
                
            return SphericalGeometry("CUSTOM", cone_angle_vec, clock_angle_vec, _id)

        @classmethod
        def from_circular_specs(cls, diameter=None, _id=None):
            """ Convert input circular specs to cone, clock angles and return corresponding :class:`instrupy.util.SphericalGeometry` object.
            
            :param diameter: (deg) Diameter of the circle.
            :paramtype diameter: float

            :param _id: Unique identifier
            :paramtype _id: str

            :return: Corresponding `SphericalGeometry` object
            :rtype: :class:`instrupy.util.SphericalGeometry`

            """
            if diameter is None:
                raise Exception("Please specify diameter of the CIRCULAR fov.")

            if(diameter < 0 or diameter > 180):
                raise Exception("Specified diameter of CIRCULAR fov must be within the range 0 deg to 180 deg")

            return SphericalGeometry("CIRCULAR", 0.5*diameter, None, _id)

        @classmethod
        def from_rectangular_specs(cls, angle_height=None, angle_width=None, _id=None):
            """ Convert the angle_height and angle_width rectangular specs to clock, cone angles and return corresponding :class:`instrupy.util.SphericalGeometry` object.

            :param angle_height: (deg) Angular height (about sensor X axis). Corresponds to along-track FOV if sensor is aligned to NADIR_POINTING frame.
            :paramtype angle_height: float

            :param angle_width: (deg) Angular width (about sensor Y axis). Corresponds to cross-track FOV if sensor is aligned to NADIR_POINTING frame.
            :paramtype angle_width: float
            
            :param _id: Unique identifier
            :paramtype _id: str

            :return: Corresponding `SphericalGeometry` object
            :rtype: :class:`instrupy.util.SphericalGeometry`                      

            """
            if(angle_height is None or angle_width is None):
                raise Exception("Please specify the angle_height and angle_width for the RECTANGULAR fov.")
            
            if(angle_height < 0 or angle_height > 180 or angle_width < 0 or angle_width > 180):
                raise Exception("Specified angle_height and angle_width of the RECTANGULAR fov must be within the range 0 deg to 180 deg")       
            
            angle_height = np.deg2rad(angle_height)
            angle_width = np.deg2rad(angle_width)

            cosCone = np.cos(angle_height/2.0)*np.cos(angle_width/2.0)
            cone = np.arccos(cosCone)

            sinClock =  np.sin(angle_height/2.0) / np.sin(cone)

            clock = np.arcsin(sinClock)

            cone = np.rad2deg(cone)
            clock = np.rad2deg(clock)

            cone_angle_vec = [cone, cone, cone, cone]

            clock_angle_vec = [clock, 180.0-clock, 180.0+clock, -clock]

            return SphericalGeometry("RECTANGULAR", cone_angle_vec, clock_angle_vec, _id)

        @staticmethod
        def from_dict(d):
            """Parses spherical geometry specifications from a normalized JSON dictionary.
    
               :param d: Dictionary with the spherical geometry specifications.
               :paramtype d: dict

               :return: Spherical geometry object
               :rtype: :class:`instrupy.util.SphericalGeometry`

            """          
            shape = SphericalGeometry.Shape.get(d.get("shape", None))

            if(shape == "CIRCULAR"):
                sph_geom_dict = SphericalGeometry.from_circular_specs(d.get("diameter", None), d.get("@id", None))
            elif(shape == "RECTANGULAR"):
                sph_geom_dict = SphericalGeometry.from_rectangular_specs(d.get("angleHeight", None), d.get("angleWidth", None),  d.get("@id", None))
            elif(shape == "CUSTOM"):
                sph_geom_dict = SphericalGeometry.from_custom_specs(d.get("customConeAnglesVector", None), d.get("customClockAnglesVector", None),  d.get("@id", None))  
            else:
                raise Exception("Invalid spherical geometry shape specified.")

            return sph_geom_dict
        
        def get_cone_clock_fov_specs(self):
            """ Function to the get the cone and clock angle vectors from the respective SphericalGeometry object.

                :return: Cone, Clock angles in degrees
                :rtype: list, float

            """
            return [self.cone_angle_vec, self.clock_angle_vec]

        @staticmethod
        def get_rect_poly_specs_from_cone_clock_angles(cone_angle_vec, clock_angle_vec):
            """ Function to get the rectangular specifications (angle_height and angle_width), from input clock, cone angle vectors.           

                :param cone_angle_vec: (deg) Array of cone angles measured from +Z sensor axis. If (:math:`xP`, :math:`yP`, :math:`zP`) is a unit vector describing a point on unit sphere, then the 
                                 cone angle for the point is :math:`\\pi/2 - \\sin^{-1}zP`. 
                :paramtype cone_angle_vec: list, float

                :param clock_angle_vec: (deg) Array of clock angles (right ascensions) measured anti-clockwise from the + X-axis. If (:math:`xP`, :math:`yP`, :math:`zP`) is a unit vector
                                        describing a point on unit sphere, then the clock angle for the point is :math:`atan2(yP,xP)`.
                :paramtype clock_angle_vec: list, float

                :return: angle_height and angle_width in degrees
                :rtype: list, float
      
                .. todo:: Make sure selected clock angle is from first quadrant. 
            """
            # Check if the instance does correspond to an rectangular shape.

            # Length of cone angle vector and clock angle vector must be 4.
            if(len(cone_angle_vec) != 4) or (len(clock_angle_vec) != 4):
                raise Exception("This SphericalGeometry instance does not correspond to a rectangular shape.")
            # Check that all elements in the cone angle vector are the same value.
            if(len(set(cone_angle_vec))!= 1):
                raise Exception("This SphericalGeometry instance does not correspond to a rectangular shape.")
            # The elements of the clock angle vector satisfy the following relationship: [theta, 180-theta, 180+theta, 360-theta]
            # in case of rectangular shape. Check for this relationship.
            if(not math.isclose(clock_angle_vec[3],(360-clock_angle_vec[0])) or not math.isclose(clock_angle_vec[1], (180 - clock_angle_vec[0])) or not math.isclose(clock_angle_vec[2], (180 + clock_angle_vec[0]))):
                raise Exception("This SphericalGeometry instance does not correspond to a rectangular shape.") 
            
            theta = np.deg2rad(cone_angle_vec[0])
            omega = np.deg2rad(clock_angle_vec[0])

            alpha = np.arcsin(np.sin(theta)*np.sin(omega))
            beta = np.arccos(np.cos(theta)/np.cos(alpha))

            angle_height = 2*np.rad2deg(alpha)
            angle_width = 2*np.rad2deg(beta)

            return [angle_height, angle_width]

        def get_fov_height_and_width(self):
            """ Get the angle_height and angle_width. Valid only for CIRCULAR and 
                RECTANGULAR shapes.

                :return: angle_height and angle_width in degrees
                :rtype: list, float
            """
            return [self.angle_height, self.angle_width]

        def __repr__(self):
            return "SphericalGeometry.from_dict({})".format(self.to_dict())

class Maneuver(Entity):
    """ Class handling the maneuverability of the satellite and/or sensor. 
        
    The maneuverability is specified with reference to the NADIR_POINTING_FRAME. The maneuver specifications 
    describe the angular-space, where the pointing axis of the sensor can be positioned. 
    
    This class includes a function which can be used to obtain the Field-Of-Regard specifications
    in terms of a proxy-sensor setup. The proxy sensor setup is characterized by orientation (wrt the NADIR_POINTING frame) of the proxy-sensor 
    and a spherical polygon/circle specification of the proxy-sensor's field-of-view. This proxy-sensor setup allows to calculate all coverage opportunities
    by the satellite-sensor pair, taking into account the satellite and/or sensor maneuverability.    

    :ivar manuever_type: Type of manuevers. Accepted values are "FIXED", "CIRCULAR", "SINGLE_ROLL_ONLY", "DOUBLE_ROLL_ONLY".
    :vartype manuever_type: str

    :ivar A_roll_min: (deg) Minimum roll angle of the 1st ROLL_ONLY region. 
    :vartype A_roll_min: float

    :ivar A_roll_max: (deg) Maximum roll angle of the 1st ROLL_ONLY region. 
    :vartype A_roll_max: float

    :ivar B_roll_min: (deg) Minimum roll angle of the 2nd ROLL_ONLY region. 
    :vartype B_roll_min: float

    :ivar B_roll_max: (deg) Maximum roll angle of the 2nd ROLL_ONLY region. 
    :vartype B_roll_max: float

    :ivar diameter: (deg) Diameter of the circular maneuver region.
    :vartype diameter: float

    :param _id: Unique identifier.
    :paramtype _id: str
    
    """

    class Category(EnumEntity):
        """ Enumeration of recognized maneuver categories. All maneuvers are with respect to the NADIR_POINTING_FRAME.

        :cvar FIXED: This option is equivalent to specifying no-maneuver. The resulting FOR is a proxy-sensor instrument with the same 
                     orientation and FOV as the original input sensor.
        :vartype FIXED: str

        :cvar CIRCULAR: This maneuver option indicates that the pointing axis can be maneuvered within a circular region (corresponding to a
                        specified angular diameter) *around* the z-axis (nadir-direction). The rotation about the pointing axis is unrestricted. 
                        The resulting FOR is characterized by a proxy-sensor is as follows:

                            * Reference orientation is the proxy-sensor aligned to the NADIR_POINTING frame.

                            * If input sensor FOV is CIRCULAR: proxy-sensor FOV is CIRCULAR with diameter = maneuver diameter + input FOV diameter

                            * If input sensor FOV is RECTANGULAR: proxy-sensor FOV  is CIRCULAR with diameter = maneuver diameter + diagonal angle of the input rectangular FOV

                                where diagonal angle of the RECTANGULAR FOV = 2 acos( cos(angle_width/2) . cos(angle_height/2) )

                        .. figure:: circular_maneuver.png
                            :scale: 75 %
                            :align: center

        :vartype CIRCULAR: str

        :cvar SINGLE_ROLL_ONLY: This maneuver option indicates that the pointing axis can be maneuvered about the roll axis (= y-axis of the NADIR_POINTING_FRAME) 
                                over a (single) range indicated by minimum and maximum roll angles. The resulting FOR characterized by a proxy-sensor is as follows:
                                
                                    * Reference orientation is the proxy-sensor at a roll-position (wrt to the NADIR_POINTING_FRAME) as follows:
                                        
                                        roll position = rollMin + 0.5 * (rollMax - rollMin)

                                    * If input sensor FOV is CIRCULAR: proxy-sensor FOV is rectangular with:
                                        
                                        angle width = (rollMax - rollMin) + input FOV diameter

                                        angle height = input sensor diameter

                                    * If input sensor FOV is RECTANGULAR: proxy-sensor FOV is rectangular with:
                                        
                                        angle width  = (rollMax - rollMin) + input FOV angle width

                                        angle height = input FOV angle height

                                .. figure:: single_rollonly_maneuver.png
                                    :scale: 75 %
                                    :align: center

        :vartype SINGLE_ROLL_ONLY: str

        :cvar DOUBLE_ROLL_ONLY: This maneuver option is similar to the SINGLE_ROLL_ONLY case, except that **two** 
                                (potentially non-overlapping) ranges of roll-angles (minimum and maximum angles).

                                .. figure:: double_rollonly_maneuver.png
                                    :scale: 75 %
                                    :align: center

        :vartype DOUBLE_ROLL_ONLY: str


        """
        FIXED = "FIXED"
        CIRCULAR = "CIRCULAR"
        SINGLE_ROLL_ONLY = "SINGLE_ROLL_ONLY"
        DOUBLE_ROLL_ONLY = "DOUBLE_ROLL_ONLY"

    def __init__(self, manuever_type=None, A_roll_min=None, A_roll_max=None, B_roll_min=None, B_roll_max=None, diameter=None, _id=None):

        self.manuever_type = Maneuver.Category.get(manuever_type) if manuever_type is not None else None
        self.A_roll_min =  float(A_roll_min) if A_roll_min is not None else None
        self.A_roll_max = float(A_roll_max) if A_roll_max is not None else None
        self.B_roll_min =  float(B_roll_min) if B_roll_min is not None else None
        self.B_roll_max = float(B_roll_max) if B_roll_max is not None else None
        self.diameter = float(diameter) if diameter is not None else None

        super(Maneuver, self).__init__(_id, "Maneuver")

    @staticmethod
    def from_dict(d):
        """Parses an maneuverability object from a normalized JSON dictionary."""
        return Maneuver(
                manuever_type = d.get("@type", None),
                A_roll_min = d.get("A_rollMin", None),
                A_roll_max = d.get("A_rollMax", None),
                B_roll_min = d.get("B_rollMin", None),
                B_roll_max = d.get("B_rollMax", None),
                diameter = d.get("diameter", None),
                )
    
    def __repr__(self):
            return "Maneuver.from_dict({})".format(self.to_dict())

    def to_dict(self):
        if self.manuever_type == Maneuver.Category.FIXED:
            specs_dict= dict({ "@type": "FIXED"})
        elif self.manuever_type == Maneuver.Category.CIRCULAR:
            specs_dict= dict({"@type": "CIRCULAR", "diameter": self.diameter})
        elif self.manuever_type == Maneuver.Category.SINGLE_ROLL_ONLY:
            specs_dict= dict({"@type": "SINGLE_ROLL_ONLY", "A_rollMin": self.A_roll_min, "A_rollMax": self.A_roll_max})
        elif self.manuever_type == Maneuver.Category.DOUBLE_ROLL_ONLY:
            specs_dict= dict({"@type": "DOUBLE_ROLL_ONLY", "A_rollMin": self.A_roll_min, "A_rollMax": self.A_roll_max, "B_rollMin": self.B_roll_min, "B_rollMax": self.B_roll_max})
        else:
            raise Exception("Invalid or no Manuever type specification.")
        
        return specs_dict
    
    def calc_field_of_regard(self, field_of_view):
        """ Calculate the field-of-regard (FOR) in terms of a *proxy sensor setup* for an input sensor FOV/ scene-FOV. 
            
            The FOR is characterized by (list of) pair of :code:`Orientation` and :code:`SphericalGeometry` objects. 
            The pair of orientation and spherical geometry objects form a *proxy-sensor setup*, which can be utilized to run coverage calculations
            and calculate all possible access opportunites by the satellite-sensor taking into account the satellite/sensor maneuverability.
            Note that only CIRCULAR or RECTANGULAR shaped sensor FOV are permitted as inputs.
            
        The returned :code:`Orientation` object characterizes the orientation of the proxy-sensor with respect to the 
        NADIR_POINTING_FRAME. The returned :code:`SphericalGeometry` object specifies the FOV of the proxy sensor. 

        In some scenarios where the FOR can have non-overlapping angular spaces (e.g. sidelooking
        SARs which can point on either "side", but cannot point at the nadir), we shall have as return a list of 
        the pairs of :code:`Orientation` and :code:`SphericalGeometry` objects, where each element of the list corresponds to a separate proxy sensor setup.
        All the proxy-sensor setups in the list must be considered to form the FOR.

        Note that always, the proxy-sensor FOV >= input sensor FOV. 

        :param field_of_view:  Field-of-view (or scene field-of-view) of the sensor. Must be either CIRCULAR or RECTANGULAR shape.
        :paramtype field_of_view: :class:`instrupy.util.SphericalGeometry`

        :return: Field-of-Regard characterized by a proxy sensor setup consiting of orientation with respect to the NADIR_POINTING_FRAME, and the coresponding spherical geometry specifications.
        :rtype: list, tuple, :class:`instrupy.util.Orientation`, :class:`instrupy.util.SphericalGeometry`

        """
        field_of_regard = None 
        mv_type = self.manuever_type

        # evaluate FOR for FIXED maneuver
        if(mv_type == 'FIXED'):
            field_of_regard = [(None, field_of_view)] # None is returned for orientation. The orientation in this 'maneuver' case has to be the original instrument orientation.

        # evaluate FOR for CIRCULAR maneuver
        if(mv_type == 'CIRCULAR'):

            if(field_of_view.shape == 'CIRCULAR'):
                proxy_fov_diameter = self.diameter + field_of_view.angle_height # field-of-regard diameter

            elif(field_of_view.shape == 'RECTANGULAR'):
                diag_half_angle = np.rad2deg(np.arccos(np.cos(np.deg2rad(0.5*field_of_view.angle_height))*np.cos(np.deg2rad(0.5*field_of_view.angle_width))))
                proxy_fov_diameter = self.diameter +  2*diag_half_angle

            else:
                raise Exception('Invalid input FOV geometry')    

            field_of_regard = [(    Orientation(ref_frame="NADIR_POINTING_FRAME"), 
                                    SphericalGeometry.from_dict({"shape": 'CIRCULAR', "diameter": proxy_fov_diameter})
                                )]                              

        def get_roll_only_mv_proxy_sen_specs(roll_min, roll_max):
            mv_angle_width_range = roll_max - roll_min # angular maneuver range
            proxy_sen_roll_angle = roll_min + 0.5*mv_angle_width_range # reference orientation

            if(field_of_view.shape == 'CIRCULAR'):
                print("Approximating FOR as rectangular shape")
                proxy_fov_angle_height = field_of_view.angle_height
                proxy_fov_angle_width =  mv_angle_width_range + field_of_view.angle_width

            elif(field_of_view.shape == 'RECTANGULAR'):
                proxy_fov_angle_height = field_of_view.angle_height
                proxy_fov_angle_width = mv_angle_width_range + field_of_view.angle_width

            else:
                raise Exception('Invalid input FOV geometry') 

            return [proxy_sen_roll_angle, proxy_fov_angle_height, proxy_fov_angle_width]

        # evaluate FOR for SINGLE_ROLL_ONLY maneuver
        if(mv_type == 'SINGLE_ROLL_ONLY'):

            [w, x, y] = get_roll_only_mv_proxy_sen_specs(self.A_roll_min, self.A_roll_max)

            field_of_regard = [(    Orientation.from_sideLookAngle(ref_frame="NADIR_POINTING_FRAME", side_look_angle=w), 
                                    SphericalGeometry.from_dict({"shape":'RECTANGULAR', "angleHeight":x, "angleWidth":y})
                               )]

        # evaluate FOR for DOUBLE_ROLL_ONLY maneuver
        if(mv_type == 'DOUBLE_ROLL_ONLY'):

            [w1, x1, y1] = get_roll_only_mv_proxy_sen_specs(self.A_roll_min, self.A_roll_max)
            [w2, x2, y2] = get_roll_only_mv_proxy_sen_specs(self.B_roll_min, self.B_roll_max)

            field_of_regard = [(    Orientation.from_sideLookAngle(ref_frame="NADIR_POINTING_FRAME", side_look_angle=w1), 
                                    SphericalGeometry.from_dict({"shape":'RECTANGULAR', "angleHeight":x1, "angleWidth":y1})
                               ),
                               (    Orientation.from_sideLookAngle(ref_frame="NADIR_POINTING_FRAME", side_look_angle=w2), 
                                    SphericalGeometry.from_dict({"shape":'RECTANGULAR', "angleHeight":x2, "angleWidth":y2})
                               )]

        
        return field_of_regard
                
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
    def eci2geo(ecicoord, JDtime):
        """ *Convert Earth-centered inertial coords to geographic spherical coordinates.*    
             https://idlastro.gsfc.nasa.gov/ftp/pro/astro/eci2geo.pro
       
        """
        [X,Y,Z] = ecicoord
        theta = np.arctan2(Y,X)
        
        gst = MathUtilityFunctions.JD2GMST(JDtime)        
        angle_sid=gst*2.0*np.pi/24.0 # sidereal angle
        lon = np.mod((theta - angle_sid),2*np.pi) 

        r=np.sqrt(X*X + Y*Y)
        lat=np.arctan2(Z,r)                         
        
        alt=r/np.cos(lat) - Constants.radiusOfEarthInKM

        lat = np.rad2deg(lat)
        lon = np.rad2deg(lon)
                
        return [lat,lon,alt]

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


