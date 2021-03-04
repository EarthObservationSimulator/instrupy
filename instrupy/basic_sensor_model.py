""" 
.. module:: basic_sensor_model

:synopsis: *Module to handle sensor model with the most basic attributes.*

"""
import copy
from instrupy.util import Entity, Orientation, SphericalGeometry, Maneuver, ViewGeometry, MathUtilityFunctions, GeoUtilityFunctions, Constants, SyntheticDataConfiguration
from netCDF4 import Dataset
import cartopy.crs as ccrs
import numpy as np
import scipy.interpolate
import metpy.interpolate
import astropy.time
import uuid
from collections import namedtuple
import logging
logger = logging.getLogger(__name__)

class BasicSensorModel(Entity):
    """A basic sensor model class. 
      
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

        :ivar fieldOfViewGeometry: Field of view spherical geometry specification of the instrument. 
        :vartype fieldOfViewGeometry: :class:`instrupy.util.SphericalGeometry`   

        :ivar fieldOfView: Field of view of instrument (spherical geometry and Orientation).
        :vartype fieldOfView: :class:`instrupy.util.ViewGeometry`

        :ivar maneuver: Maneuver specification of the instrument. TODO: Modify behavior to have FOR =FOV when no maneuver is specified (hence fixed pointing).
        :vartype maneuver: :class:`instrupy.util.Maneuver`  

        :ivar fieldOfRegard: Field of regard of the instrument taking into account the sensor FOV and manueverability of the satellite-sensor system. 
                             Note that this shall be a list in order to accommodate non-intersecting view geometries.
        :vartype fieldOfRegard: list, :class:`instrupy.util.ViewGeometry`  
       
        :ivar dataRate: Rate of data recorded (Mega bits per sec) during nominal operations.
        :vartype dataRate: float  

        :ivar syntheticDataConfig: Synthetic data configuration for the instrument.
        :vartype syntheticDataConfig: :class:`instrupy.util.SyntheticDataConfiguration`  

        :ivar bitsPerPixel: Number of bits encoded per pixel of image.
        :vartype bitsPerPixel: int        

        :ivar numberDetectorRows: Number of detector rows (along the Y-axis of the SENOR_BODY_FIXED frame).
        :vartype numberDetectorRows: int

        :ivar numberDetectorCols: Number of detector columns (along the X-axis of the SENOR_BODY_FIXED frame).
        :vartype numberDetectorCols: int

        :ivar _id: Unique instrument identifier.
        :vartype _id: str

        .. figure:: detector_config.png
            :scale: 75 %
            :align: center
   
    """

    def __init__(self, name=None, mass=None, volume=None, power=None,  orientation=None,
                 fieldOfViewGeometry=None, maneuver=None, dataRate=None, syntheticDataConfig=None, bitsPerPixel=None, 
                 numberDetectorRows=None, numberDetectorCols=None, _id=None):
        """Initialization

        """
        self.name = str(name) if name is not None else None
        self.mass = float(mass) if mass is not None else None
        self.volume = float(volume) if volume is not None else None
        self.power = float(power) if power is not None else None
        self.orientation = copy.deepcopy(orientation) if orientation is not None and isinstance(orientation, Orientation) else None
        self.fieldOfViewGeometry = copy.deepcopy(fieldOfViewGeometry) if fieldOfViewGeometry is not None and isinstance(fieldOfViewGeometry, SphericalGeometry) else None
        self.fieldOfView = ViewGeometry(orien=self.orientation, sph_geom=self.fieldOfViewGeometry) if self.orientation is not None and self.fieldOfViewGeometry is not None else None
        self.maneuver = copy.deepcopy(maneuver) if maneuver is not None and isinstance(maneuver, Maneuver) else None
        self.fieldOfRegard = self.maneuver.calc_field_of_regard(self.fieldOfViewGeometry) if (self.maneuver is not None and self.fieldOfViewGeometry is not None) else None
        self.dataRate = float(dataRate) if dataRate is not None else None
        self.syntheticDataConfig = copy.deepcopy(syntheticDataConfig) if syntheticDataConfig is not None and isinstance(syntheticDataConfig, SyntheticDataConfiguration) else None
        self.bitsPerPixel = int(bitsPerPixel) if bitsPerPixel is not None else None                
        self.numberDetectorRows = int(numberDetectorRows) if numberDetectorRows is not None else None
        self.numberDetectorCols = int(numberDetectorCols) if numberDetectorCols is not None else None
        self._id = _id if _id is not None else None

        super(BasicSensorModel,self).__init__(_id, "Basic Sensor")

    @staticmethod
    def from_dict(d):
        """ Parses an basic sensor model from a normalized JSON dictionary.
        
        The following default values are assigned to the object instance parameters in case of 
        :class:`None` values or missing key/value pairs in the input dictionary.

        .. csv-table:: Default values
            :header: Parameter, Default Value
            :widths: 10,40

            orientation, Referenced and aligned to the SC_BODY_FIXED frame.
            fieldOfView, CIRCULAR shape with 25 deg diameter.
            numberDetectorRows, 4
            numberDetectorRows, 4
            _id, random string
        
        :param d: Normalized JSON dictionary with the corresponding model specifications. 
        :paramtype d: dict

        :returns: BasicSensorModel object initialized with the input specifications.
        :rtype: :class:`instrupy.BasicSensorModel`
            
        """
        default_fov = dict({'shape': 'CIRCULAR', 'diameter':25}) # default fov is a 25 deg diameter circular
        default_orien = dict({"referenceFrame": "SC_BODY_FIXED", "convention": "REF_FRAME_ALIGNED"}) #  default orientation = referenced and aligned to the SC_BODY_FIXED frame.
        return BasicSensorModel(
                name = d.get("name", None),
                mass = d.get("mass", None),
                volume = d.get("volume", None),
                power = d.get("power", None),
                orientation = Orientation.from_json(d.get("orientation", default_orien)),
                fieldOfViewGeometry =  SphericalGeometry.from_json(d.get("fieldOfViewGeometry", default_fov)),
                maneuver =  Maneuver.from_json(d.get("maneuver", None)),
                dataRate = d.get("dataRate", None),
                bitsPerPixel = d.get("bitsPerPixel", None),
                syntheticDataConfig = SyntheticDataConfiguration.from_json(d.get("syntheticDataConfig", None)),
                numberDetectorRows = d.get("numberDetectorRows", 4),
                numberDetectorCols = d.get("numberDetectorCols", 4),
                _id = d.get("@id", str(uuid.uuid4()))
                )

    def to_dict(self):
        """ Translate the BasicSensorModel object to a Python dictionary such that it can be uniquely reconstructed back from the dictionary.

        :returns: BasicSensorModel specifications as python dictionary.
        :rtype: dict

        """
        fieldOfViewGeometry_dict = self.fieldOfViewGeometry.to_dict() if self.fieldOfViewGeometry is not None and isinstance(self.fieldOfViewGeometry, SphericalGeometry) else None
        orientation_dict = self.orientation.to_dict() if self.orientation is not None and isinstance(self.orientation, Orientation) else None
        maneuver_dict = self.maneuver.to_dict() if self.maneuver is not None and isinstance(self.maneuver, Maneuver) else None
        syntheticDataConfig_dict = self.syntheticDataConfig.to_dict() if self.syntheticDataConfig is not None and isinstance(self.syntheticDataConfig, SyntheticDataConfiguration) else None
        return dict({
                "@type": "Basic Sensor",
                "name":self.name,
                "mass":self.mass,
                "volume":self.volume,
                "power":self.power,
                "fieldOfViewGeometry":fieldOfViewGeometry_dict,
                "orientation":orientation_dict,
                "maneuver":maneuver_dict,
                "dataRate":self.dataRate,
                "bitsPerPixel": self.bitsPerPixel,
                "syntheticDataConfig": syntheticDataConfig_dict,
                "numberDetectorRows": self.numberDetectorRows,
                "numberDetectorCols": self.numberDetectorCols,
                "@id": self._id
                })

    def __repr__(self):
        return "BasicSensorModel.from_dict({})".format(self.to_dict())

    def __eq__(self, other):
        """ Simple equality check. Returns True if the class attributes are equal, else returns False. 
            The derived attributes ``fieldOfView`` and ``fieldOfRegard`` are not checked.
            Note that the ``_id`` data attribute could be different.
        """
        return (self.name == other.name and self.mass == other.mass and self.volume == other.volume and self.power==other.power and self.fieldOfViewGeometry == other.fieldOfViewGeometry and
                self.orientation == other.orientation and self.maneuver == other.maneuver and self.dataRate == other.dataRate and self.bitsPerPixel ==other.bitsPerPixel and 
                self.syntheticDataConfig == other.syntheticDataConfig and self.numberDetectorCols == other.numberDetectorCols and self.numberDetectorRows == other.numberDetectorRows)

    def calc_data_metrics(self, sc_orbit_state, target_coords):
        """ Calculate typical observation data metrics. 
            
            .. figure:: target_geom_3D.png
                :scale: 75 %
                :align: center
            
            .. figure:: target_geom_2D.png
                :scale: 75 %
                :align: center

            *   :math:`\\mathbf{R = T - S}`
            *   :math:`\\gamma = \\cos^{-1}(\\mathbf{\\dfrac{R}{|R|}} \\cdot \\mathbf{\\dfrac{-S}{|S|}})`
            *   :math:`\\theta_i = \\sin^{-1}(\\sin\\gamma  \\hspace{1mm}  \\dfrac{R_E + h}{R_E})`

            Assuming spherical Earth of radius :math:`R_E`

            where,

            * :math:`\\mathbf{S}`: Position-vector of the satellite in the EARTH_CENTERED_INERTIAL frame.
            * :math:`\\mathbf{T}`: Position-vector of the target ground-point in the EARTH_CENTERED_INERTIAL frame.
            * :math:`\\mathbf{R}`: Range vector from satellite to target ground point.
            * :math:`\\gamma`:  Look-angle to target ground point from satellite.
            * :math:`\\theta_i`: Incidence-angle at the target ground point.
            * :math:`R_E`: Nominal equatorial radius of Earth.
            * :math:`h`: Altitude of satellite.
            
            
            Please refer to the :class:`instrupy.util.GeoUtilityFunctions.compute_sun_zenith` function for description of the calculation of the Sun-zenith angle.

        :param sc_orbit_state: Spacecraft state at the time of observation.

                            Dictionary keys are: 
                            
                            * :code:`time [JDUT1]` (:class:`float`), Time in Julian Day UT1. Corresponds to the time of observation. 
                            * :code:`x [km]` (:class:`float`), :code:`y [km]` (:class:`float`), :code:`z [km]` (:class:`float`), Cartesian spatial coordinates of satellite in EARTH_CENTERED_INERTIAL frame at the time of observation.
                            * :code:`vx [km/s]` (:class:`float`), :code:`vy [km/s]` (:class:`float`), :code:`vz [km/s]` (:class:`float`), Velocity of spacecraft in EARTH_CENTERED_INERTIAL frame at the time of observation.
        :paramtype sc_orbit_state: dict

        
        :param target_coords: Location of the observation. Also sometimes called the Point-Of-Interest (POI).

                            Dictionary keys are: 
                            
                            * :code:`lat [deg]` (:class:`float`), :code:`lon [deg]` (:class:`float`) indicating the corresponding ground-point accessed (latitude, longitude) in degrees.
        :paramtype target_coords: dict

        :returns: Calculated observation data metrics.
                    
                    Dictionary keys are: 
                
                    * :code:`coverage [T/F]` (:class:`bool`) indicating if observation was possible during the access event. Always true for the :class:`BasicSensor` type.
                    * :code:`incidence angle [deg]` (:class:`float`) Incidence angle in degrees at target point calculated assuming spherical Earth.
                    * :code:`look angle [deg]` (:class:`float`) Look angle in degrees at target point calculated assuming spherical Earth. Positive sign => look is in positive half-space made by the orbit-plane (i.e. orbit plane normal vector) and vice-versa.
                    * :code:`observation range [km]` (:class:`float`) Distance in kilometers from satellite to ground-point during the observation.
                    * :code:`solar zenith [deg]` (:class:`float`) Solar Zenith angle in degrees during observation.

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

        #  Calculate the range vector between spacecraft and POI (Target)
        range_vector_km = target_pos - sc_pos

        range_km = np.linalg.norm(range_vector_km)

        # Calculate look angle
        look_angle = np.arccos(np.dot(MathUtilityFunctions.normalize(range_vector_km), -1*MathUtilityFunctions.normalize(sc_pos)))
        look_angle_deg = np.rad2deg(look_angle)
        
        # Look angle to corresponding incidence angle conversion for spherical Earth
        incidence_angle = np.arcsin(np.sin(look_angle)*(Constants.radiusOfEarthInKM + alt_km)/Constants.radiusOfEarthInKM)
        incidence_angle_deg =  np.rad2deg(incidence_angle)

        # Solar zenith angle
        [solar_zenith, solar_distance] = GeoUtilityFunctions.compute_sun_zenith(tObs_JDUT1, target_pos)
        if solar_zenith is not None:
            solar_zenith_deg =  np.rad2deg(solar_zenith)
        else:
            solar_zenith_deg = np.nan

        # assign sign to look-angle. 
        # positive sign => Positive sign => look is in positive half-space made by the orbit-plane (i.e. orbit plane normal vector) and vice-versa.
        orbit_normal = np.cross(sc_pos, sc_vel)
        sgn = np.sign(np.dot(range_vector_km, orbit_normal))
        if(sgn==0):
            sgn = 1

        # For a basic-sensor, coverage is true for all access-events
        isCovered = True 
    
        obsv_metrics = {}
        obsv_metrics["observation range [km]"] = range_km
        obsv_metrics["look angle [deg]"] = sgn*look_angle_deg
        obsv_metrics["incidence angle [deg]"] = incidence_angle_deg
        obsv_metrics["solar zenith [deg]"] = solar_zenith_deg
        obsv_metrics["coverage [T/F]"] = isCovered

        return obsv_metrics

    def synthesize_observation(self, time_JDUT1, pixel_center_pos):
        """ Compute the value of the geophysical variable (associated with the instrument) at the 
            input pixel (center) positions. 
            
            The interpolation (in space) of the source data is done using the interpolation method
            indicated in the synthetic data configuration of the instrument. The data file closest 
            to the input time of observation is chosen as the source data file.

        :param time_JDUT1: Observation time in Julian Date UT1.
        :paramtype time_JDUT1: float

        :param pixel_center_pos: Array of 'N' pixel center positions (longitudes, latitudes) of the pixels in the observation.
        :paramtype pixel_center_pos: array_like, shape (N, 2), float
        
        :returns: The interpolated data along with the corresponding pixel-center positions. The geophysical variable is also returned. Contents of the
                  namestuple are: 
                  
                  1. Pixel center positions (lon, lat) 
                  2. Interpolated data at the Pixel center positions
                  3. Geophysical variable
        :rtype: namedtuple, (array_like, shape (N, 2), float), (array_like, shape (N, 1), float), (str)

        """        
        source_fps = self.syntheticDataConfig.sourceFilePaths
        geophy_var = self.syntheticDataConfig.geophysicalVar
        interpl = self.syntheticDataConfig.get_interpolator()
        
        # get the source data file closest to the input time
        # collect all the datasets at the different forecast times
        forecast_time_hrs = [] # list of (relative wrt initial_time) forecast times
        for _src_fp in source_fps:
            dataset  = Dataset(_src_fp, "r", format="NETCDF4")
            _initial_abs_time = dataset[geophy_var].initial_time # start of the forecast (absolute time)
            _forecast_time = dataset[geophy_var].forecast_time # forecast time relative to the initial_time
            _forecast_time_units = dataset[geophy_var].forecast_time_units
            if(_forecast_time_units != 'hours'):
                raise Exception('Forecast time units must be of hours')            
            forecast_time_hrs.append(_forecast_time)
        
        # convert the initial_time in UTC (string format) to JDUT1 (example string is: 12/11/2020 (12:00))
        t = astropy.time.Time.strptime(_initial_abs_time, '%m/%d/%Y (%H:%M)') 
        initial_time_UTC = astropy.time.Time(t, scale='utc')
        initial_time_UT1 = initial_time_UTC.ut1 # convert from UTC to UT1
        initial_time_JDUT1 = initial_time_UT1.jd # convert to Julian Date format

        # required forecast time relative to the initial_time
        required_forecast_time_hrs = (time_JDUT1 - initial_time_JDUT1)*24 

        # find dataset closest to the required forecast time
        idx = (np.abs(forecast_time_hrs - required_forecast_time_hrs)).argmin()

        dataset  = Dataset(source_fps[idx], "r", format="NETCDF4")
        lons = dataset.variables['lon_0'][:]
        lats = dataset.variables['lat_0'][:]
        var_data = dataset.variables[geophy_var][:]

        #print(lons)
        #print(lats)
        #print(temperature)

        # interpolate the variable data to the pixel center positions   
        interpl_data = interpl(lons, lats, var_data, pixel_center_pos)

        #print(result)

        '''
        ax = plt.axes(projection=ccrs.PlateCarree())
        ax.coastlines()

        plt.contourf(lons, lats, temperature,60, transform=ccrs.PlateCarree(), cmap=get_cmap("jet"))
        plt.colorbar(ax=ax, shrink=.98)
        plt.show()
        '''
        dataset.close()

        # return the interpolated variable data
        syn_data = namedtuple("syn_data", ["pixel_center_pos", "interpl_data", "geophy_var"])

        return syn_data(
            pixel_center_pos,
            interpl_data,
            geophy_var
        )
    
    def get_id(self):
        """ Get the object identifier.

        :returns: Object identifier.
        :rtype: str

        """
        return self._id
    
    def get_field_of_view(self):
        return self.fieldOfView

    def get_field_of_regard(self):
        """ Get the instrument field of regard.

        :returns: Field of regard. 
        :rtype: list, :class:`instrupy.util.ViewGeometry`

        """
        return self.fieldOfRegard

    def get_orientation(self):
        return self.orientation
    
    def get_pixel_config(self):
        """ Get pixel-configuration.
        
        :returns: Pixel-configuration, i.e. the number of detector rows and columns.
        :rtype: namedtuple, (int, int)
        
        """        
        pixel_config = namedtuple("pixel_config", ["numberDetectorRows", "numberDetectorCols"])

        return pixel_config(self.numberDetectorRows, self.numberDetectorCols)

